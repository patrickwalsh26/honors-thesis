"""
Differential Privacy mechanisms for privacy-preserving phenotype matching.

Implements:
- Laplace mechanism for (epsilon, 0)-differential privacy
- Gaussian mechanism for (epsilon, delta)-differential privacy
- Exponential mechanism for private selection
- Privacy-preserving similarity calculator wrapper
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Optional, Callable, Any
import logging

logger = logging.getLogger(__name__)


class DPMechanism(ABC):
    """Abstract base class for differential privacy mechanisms."""

    def __init__(self, epsilon: float, delta: float = 0.0):
        """
        Initialize DP mechanism.

        Args:
            epsilon: Privacy budget (lower = more private)
            delta: Probability of privacy failure (0 for pure DP)
        """
        if epsilon <= 0:
            raise ValueError("Epsilon must be positive")
        if delta < 0 or delta >= 1:
            raise ValueError("Delta must be in [0, 1)")

        self.epsilon = epsilon
        self.delta = delta

    @abstractmethod
    def add_noise(self, value: float, sensitivity: float) -> float:
        """
        Add calibrated noise to a value.

        Args:
            value: The true value to privatize
            sensitivity: The sensitivity of the query (max change from one record)

        Returns:
            Noisy value satisfying differential privacy
        """
        pass

    def privatize_array(self, values: np.ndarray, sensitivity: float) -> np.ndarray:
        """
        Add noise to an array of values.

        Args:
            values: Array of true values
            sensitivity: Sensitivity per element

        Returns:
            Array of noisy values
        """
        return np.array([self.add_noise(v, sensitivity) for v in values])


class LaplaceMechanism(DPMechanism):
    """
    Laplace mechanism for (epsilon, 0)-differential privacy.

    Adds noise drawn from Laplace(0, sensitivity/epsilon) distribution.
    Provides pure differential privacy guarantee.
    """

    def __init__(self, epsilon: float, delta: float = 0.0):
        """
        Initialize Laplace mechanism.

        Args:
            epsilon: Privacy budget
            delta: Must be 0 for Laplace mechanism (pure DP)
        """
        if delta != 0.0:
            logger.warning("Laplace mechanism ignores delta; using pure DP")
        super().__init__(epsilon, 0.0)

    def add_noise(self, value: float, sensitivity: float) -> float:
        """
        Add Laplace noise.

        Noise scale b = sensitivity / epsilon
        """
        scale = sensitivity / self.epsilon
        noise = np.random.laplace(0, scale)
        return value + noise

    def privatize_similarity(self, similarity: float) -> float:
        """
        Add noise to a similarity score.

        Similarity scores are in [0, 1], so sensitivity = 1.0.
        Output is clipped to [0, 1].

        Args:
            similarity: True similarity score in [0, 1]

        Returns:
            Noisy similarity score clipped to [0, 1]
        """
        noisy = self.add_noise(similarity, 1.0)
        return np.clip(noisy, 0.0, 1.0)

    def privatize_count(self, count: int) -> float:
        """
        Add noise to a count query.

        Counting queries have sensitivity = 1 (one record changes count by 1).

        Args:
            count: True count

        Returns:
            Noisy count (may be negative, typically floor/round in practice)
        """
        noisy = self.add_noise(float(count), 1.0)
        return max(0.0, noisy)


class GaussianMechanism(DPMechanism):
    """
    Gaussian mechanism for (epsilon, delta)-differential privacy.

    Adds noise drawn from N(0, sigma^2) where:
    sigma = sensitivity * sqrt(2 * ln(1.25/delta)) / epsilon

    Provides approximate differential privacy with probability 1-delta.
    """

    def __init__(self, epsilon: float, delta: float = 1e-5):
        """
        Initialize Gaussian mechanism.

        Args:
            epsilon: Privacy budget
            delta: Probability of privacy failure (must be > 0)
        """
        if delta <= 0:
            raise ValueError("Gaussian mechanism requires delta > 0")
        super().__init__(epsilon, delta)

        # Pre-compute sigma factor
        self._sigma_factor = np.sqrt(2 * np.log(1.25 / delta))

    def add_noise(self, value: float, sensitivity: float) -> float:
        """
        Add Gaussian noise.

        sigma = sensitivity * sqrt(2 * ln(1.25/delta)) / epsilon
        """
        sigma = sensitivity * self._sigma_factor / self.epsilon
        noise = np.random.normal(0, sigma)
        return value + noise

    def privatize_similarity(self, similarity: float) -> float:
        """Add noise to similarity score, clipped to [0, 1]."""
        noisy = self.add_noise(similarity, 1.0)
        return np.clip(noisy, 0.0, 1.0)


class ExponentialMechanism(DPMechanism):
    """
    Exponential mechanism for private selection from discrete options.

    Selects option with probability proportional to:
    exp(epsilon * utility(option) / (2 * sensitivity))

    Useful for private top-k selection where we want to select
    high-utility options without revealing exact utilities.
    """

    def add_noise(self, value: float, sensitivity: float) -> float:
        """Not applicable for exponential mechanism."""
        raise NotImplementedError("Exponential mechanism uses select(), not add_noise()")

    def select(
        self,
        options: List[Any],
        utility_fn: Callable[[Any], float],
        sensitivity: float
    ) -> Any:
        """
        Select an option with privacy.

        Args:
            options: List of options to choose from
            utility_fn: Function that returns utility score for each option
            sensitivity: Max change in utility from adding/removing one record

        Returns:
            Selected option (with probability proportional to exp(epsilon*u/(2*s)))
        """
        if not options:
            raise ValueError("Options list cannot be empty")

        # Compute utilities
        utilities = np.array([utility_fn(opt) for opt in options])

        # Compute log probabilities (for numerical stability)
        log_probs = (self.epsilon * utilities) / (2 * sensitivity)

        # Normalize using log-sum-exp trick
        log_probs = log_probs - np.max(log_probs)
        probs = np.exp(log_probs)
        probs = probs / probs.sum()

        # Sample
        idx = np.random.choice(len(options), p=probs)
        return options[idx]

    def select_top_k(
        self,
        options: List[Any],
        utility_fn: Callable[[Any], float],
        sensitivity: float,
        k: int
    ) -> List[Any]:
        """
        Select top-k options with privacy using iterative exponential mechanism.

        Note: Each selection consumes epsilon privacy budget.
        Total budget for k selections is k * epsilon.

        Args:
            options: List of options to choose from
            utility_fn: Function that returns utility score for each option
            sensitivity: Max change in utility from one record
            k: Number of options to select

        Returns:
            List of k selected options
        """
        selected = []
        remaining = list(options)

        for _ in range(min(k, len(options))):
            choice = self.select(remaining, utility_fn, sensitivity)
            selected.append(choice)
            remaining.remove(choice)

        return selected


class ReportNoisyMax:
    """
    Report Noisy Max mechanism for private argmax.

    Alternative to exponential mechanism that's often more practical.
    Adds noise to each option's score and returns the argmax.

    Privacy guarantee: epsilon-DP when using Laplace noise.
    """

    def __init__(self, epsilon: float):
        """
        Initialize Report Noisy Max.

        Args:
            epsilon: Privacy budget
        """
        self.epsilon = epsilon
        self._laplace = LaplaceMechanism(epsilon)

    def select(
        self,
        scores: np.ndarray,
        sensitivity: float = 1.0
    ) -> Tuple[int, float]:
        """
        Select argmax with privacy.

        Args:
            scores: Array of scores
            sensitivity: Score sensitivity

        Returns:
            Tuple of (selected index, noisy score)
        """
        # Add noise to each score
        noisy_scores = self._laplace.privatize_array(scores, sensitivity)

        # Return argmax
        idx = np.argmax(noisy_scores)
        return int(idx), noisy_scores[idx]

    def select_top_k(
        self,
        scores: np.ndarray,
        sensitivity: float,
        k: int
    ) -> List[Tuple[int, float]]:
        """
        Select top-k indices with privacy.

        Args:
            scores: Array of scores
            sensitivity: Score sensitivity
            k: Number of top results

        Returns:
            List of (index, noisy_score) tuples sorted by noisy score
        """
        noisy_scores = self._laplace.privatize_array(scores, sensitivity)

        # Get top-k indices
        top_indices = np.argsort(noisy_scores)[::-1][:k]

        return [(int(idx), noisy_scores[idx]) for idx in top_indices]


class PrivacyAccountant:
    """
    Track cumulative privacy budget spent.

    Implements basic composition theorem:
    - Sequential composition: epsilons add
    - Advanced composition available for tighter bounds
    """

    def __init__(self, total_budget: float):
        """
        Initialize privacy accountant.

        Args:
            total_budget: Maximum epsilon to spend
        """
        self.total_budget = total_budget
        self._spent = 0.0
        self._queries = 0

    def spend(self, epsilon: float) -> bool:
        """
        Attempt to spend privacy budget.

        Args:
            epsilon: Budget to spend

        Returns:
            True if budget available and spent, False if insufficient
        """
        if self._spent + epsilon > self.total_budget:
            return False

        self._spent += epsilon
        self._queries += 1
        return True

    @property
    def remaining(self) -> float:
        """Remaining privacy budget."""
        return self.total_budget - self._spent

    @property
    def spent(self) -> float:
        """Total budget spent."""
        return self._spent

    @property
    def num_queries(self) -> int:
        """Number of queries processed."""
        return self._queries

    def get_report(self) -> Dict[str, float]:
        """Get privacy accounting summary."""
        return {
            "total_budget": self.total_budget,
            "spent": self._spent,
            "remaining": self.remaining,
            "num_queries": self._queries,
            "utilization": self._spent / self.total_budget
        }


class PrivateSimilarityCalculator:
    """
    Wrapper for privacy-preserving similarity computation.

    Composes a DP mechanism with an existing similarity calculator.
    Follows the same interface as PhenopacketSimilarityCalculator.
    """

    def __init__(
        self,
        base_calculator,
        mechanism: DPMechanism,
        sensitivity: float = 1.0,
        privacy_budget: Optional[float] = None
    ):
        """
        Initialize private similarity calculator.

        Args:
            base_calculator: Underlying PhenopacketSimilarityCalculator
            mechanism: DP mechanism to use
            sensitivity: Sensitivity of similarity function
            privacy_budget: Optional total budget limit
        """
        self.base_calculator = base_calculator
        self.mechanism = mechanism
        self.sensitivity = sensitivity

        # Privacy accounting
        self._accountant = PrivacyAccountant(privacy_budget) if privacy_budget else None

    def compute_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict
    ) -> float:
        """
        Compute similarity with DP noise.

        Args:
            phenopacket1: First phenopacket
            phenopacket2: Second phenopacket

        Returns:
            Noisy similarity score
        """
        # Check privacy budget
        if self._accountant and not self._accountant.spend(self.mechanism.epsilon):
            logger.warning("Privacy budget exhausted")
            return 0.0

        # Get true similarity
        true_sim = self.base_calculator.compute_similarity(phenopacket1, phenopacket2)

        # Add noise (handles clipping to [0, 1])
        if isinstance(self.mechanism, (LaplaceMechanism, GaussianMechanism)):
            return self.mechanism.privatize_similarity(true_sim)
        else:
            return np.clip(self.mechanism.add_noise(true_sim, self.sensitivity), 0.0, 1.0)

    def find_most_similar(
        self,
        query_phenopacket: Dict,
        candidate_phenopackets: List[Dict],
        top_k: int = 10,
        use_exponential: bool = False
    ) -> List[Tuple[int, float]]:
        """
        Find most similar phenopackets with privacy.

        Uses Report Noisy Max by default (more practical).
        Can use Exponential Mechanism if specified.

        Args:
            query_phenopacket: Query phenopacket
            candidate_phenopackets: Candidates to search
            top_k: Number of results
            use_exponential: Use exponential mechanism instead of noisy max

        Returns:
            List of (index, noisy_similarity) tuples
        """
        # Compute all true similarities
        true_sims = np.array([
            self.base_calculator.compute_similarity(query_phenopacket, cand)
            for cand in candidate_phenopackets
        ])

        if use_exponential:
            # Use exponential mechanism
            exp_mech = ExponentialMechanism(self.mechanism.epsilon)
            indices = list(range(len(candidate_phenopackets)))

            selected = exp_mech.select_top_k(
                options=indices,
                utility_fn=lambda i: true_sims[i],
                sensitivity=self.sensitivity,
                k=top_k
            )

            # Return with noisy scores
            results = []
            for idx in selected:
                noisy_sim = self.mechanism.privatize_similarity(true_sims[idx])
                results.append((idx, noisy_sim))
            return results

        else:
            # Use Report Noisy Max (more practical)
            rnm = ReportNoisyMax(self.mechanism.epsilon)
            results = rnm.select_top_k(true_sims, self.sensitivity, top_k)

            # Clip scores to [0, 1]
            return [(idx, np.clip(score, 0.0, 1.0)) for idx, score in results]

    def get_privacy_report(self) -> Dict[str, Any]:
        """Get privacy accounting summary."""
        report = {
            "mechanism": type(self.mechanism).__name__,
            "epsilon_per_query": self.mechanism.epsilon,
            "delta": self.mechanism.delta,
            "sensitivity": self.sensitivity
        }

        if self._accountant:
            report.update(self._accountant.get_report())

        return report


def create_dp_mechanism(
    mechanism_type: str = "laplace",
    epsilon: float = 1.0,
    delta: float = 1e-5
) -> DPMechanism:
    """
    Factory function to create DP mechanism from config.

    Args:
        mechanism_type: Type of mechanism ("laplace" or "gaussian")
        epsilon: Privacy budget
        delta: Privacy failure probability (for Gaussian)

    Returns:
        Configured DP mechanism
    """
    if mechanism_type.lower() == "laplace":
        return LaplaceMechanism(epsilon)
    elif mechanism_type.lower() == "gaussian":
        return GaussianMechanism(epsilon, delta)
    else:
        raise ValueError(f"Unknown mechanism type: {mechanism_type}")


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Test Laplace mechanism
    laplace = LaplaceMechanism(epsilon=1.0)

    # Add noise to similarity score
    true_sim = 0.75
    noisy_sims = [laplace.privatize_similarity(true_sim) for _ in range(100)]

    print(f"True similarity: {true_sim}")
    print(f"Noisy mean: {np.mean(noisy_sims):.4f}")
    print(f"Noisy std: {np.std(noisy_sims):.4f}")

    # Test Gaussian mechanism
    gaussian = GaussianMechanism(epsilon=1.0, delta=1e-5)
    noisy_sims_g = [gaussian.privatize_similarity(true_sim) for _ in range(100)]

    print(f"\nGaussian noisy mean: {np.mean(noisy_sims_g):.4f}")
    print(f"Gaussian noisy std: {np.std(noisy_sims_g):.4f}")

    # Test privacy accountant
    accountant = PrivacyAccountant(total_budget=5.0)
    for i in range(6):
        success = accountant.spend(1.0)
        print(f"Query {i+1}: {'allowed' if success else 'DENIED'}, remaining: {accountant.remaining}")
