"""
Unified privacy-preserving phenotype similarity calculator.

Composes all privacy mechanisms:
- Private Set Intersection (PSI) for secure overlap
- Differential Privacy (DP) for noisy outputs
- k-Anonymity for result filtering
- Rare term suppression

Provides the same interface as PhenopacketSimilarityCalculator.
"""

from typing import Dict, List, Tuple, Optional, Set, Any
import numpy as np
import logging
import copy

from .differential_privacy import (
    DPMechanism, LaplaceMechanism, GaussianMechanism,
    PrivateSimilarityCalculator, PrivacyAccountant, create_dp_mechanism
)
from .k_anonymity import (
    RareTermFilter, KAnonymityGuard, PrivacyConfig, load_privacy_config
)
from .psi import HybridPSI, PSIPhenopacketMatcher, create_psi_matcher

logger = logging.getLogger(__name__)


class PrivacyPreservingCalculator:
    """
    Unified privacy-preserving phenotype similarity calculator.

    Composes multiple privacy mechanisms in a configurable pipeline:
    1. Rare term filtering (suppress potentially identifying terms)
    2. PSI overlap computation (secure intersection)
    3. DP noise addition (differential privacy)
    4. k-anonymity enforcement (suppress small result sets)

    Follows the same interface as PhenopacketSimilarityCalculator.
    """

    def __init__(
        self,
        base_similarity=None,
        config: Optional[PrivacyConfig] = None,
        hpo_manager=None,
        use_psi: bool = True,
        use_dp: bool = True,
        use_k_anonymity: bool = True,
        use_rare_filter: bool = True
    ):
        """
        Initialize privacy-preserving calculator.

        Args:
            base_similarity: Base HPOSimilarity metric (optional if using PSI)
            config: Privacy configuration
            hpo_manager: HPOManager for ontology operations
            use_psi: Enable Private Set Intersection
            use_dp: Enable Differential Privacy
            use_k_anonymity: Enable k-anonymity enforcement
            use_rare_filter: Enable rare term filtering
        """
        self.base_similarity = base_similarity
        self.config = config or PrivacyConfig.default()
        self.hpo_manager = hpo_manager

        # Feature flags
        self._use_psi = use_psi
        self._use_dp = use_dp
        self._use_k_anonymity = use_k_anonymity
        self._use_rare_filter = use_rare_filter

        # Initialize components
        self._init_components()

        # Privacy accounting
        self._queries = 0
        self._total_epsilon_spent = 0.0

        # Caching
        self._cache: Dict[Tuple[str, str], float] = {}
        self.cache_similarities = True

    def _init_components(self):
        """Initialize privacy mechanism components."""
        # Rare term filter
        self.rare_filter = RareTermFilter(
            min_prevalence=self.config.min_prevalence,
            hpo_manager=self.hpo_manager
        )

        # k-anonymity guard
        self.k_guard = KAnonymityGuard(k=self.config.k)

        # Differential privacy mechanism
        self.dp_mechanism = create_dp_mechanism(
            mechanism_type=self.config.dp_mechanism,
            epsilon=self.config.epsilon,
            delta=self.config.delta
        )

        # PSI matcher
        if self._use_psi:
            self.psi = HybridPSI(
                hpo_manager=self.hpo_manager,
                use_ancestors=True,
                simulated=True  # Use simulated mode for thesis evaluation
            )
        else:
            self.psi = None

        # Privacy accountant
        self._accountant = PrivacyAccountant(total_budget=float('inf'))

    def set_corpus(self, phenopackets: List[Dict]):
        """
        Set the corpus for computing term prevalence.

        Should be called before filtering to compute which terms are rare.

        Args:
            phenopackets: List of phenopackets for prevalence computation
        """
        self.rare_filter.compute_prevalence_from_corpus(phenopackets)
        logger.info(f"Computed prevalence for {len(self.rare_filter.term_prevalence)} terms")

    def _extract_terms(self, phenopacket: Dict) -> Set[str]:
        """Extract HPO terms from phenopacket."""
        terms = set()
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                terms.add(feature["type"]["id"])
        return terms

    def _filter_phenopacket(self, phenopacket: Dict) -> Dict:
        """Apply rare term filtering to phenopacket."""
        if not self._use_rare_filter:
            return phenopacket

        return self.rare_filter.filter_phenopacket(
            phenopacket,
            strategy=self.config.rare_term_strategy
        )

    def compute_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict
    ) -> float:
        """
        Compute similarity with privacy protections.

        Pipeline:
        1. Filter rare terms from both phenopackets
        2. Compute similarity (via PSI or base metric)
        3. Add DP noise if enabled

        Args:
            phenopacket1: First phenopacket
            phenopacket2: Second phenopacket

        Returns:
            Privacy-preserving similarity score
        """
        self._queries += 1

        # Check cache
        id1 = phenopacket1.get("id")
        id2 = phenopacket2.get("id")

        if self.cache_similarities and id1 and id2:
            cache_key = tuple(sorted([id1, id2]))
            if cache_key in self._cache:
                return self._cache[cache_key]

        # Step 1: Filter rare terms
        pp1 = self._filter_phenopacket(phenopacket1)
        pp2 = self._filter_phenopacket(phenopacket2)

        # Step 2: Compute similarity
        if self._use_psi and self.psi:
            terms1 = self._extract_terms(pp1)
            terms2 = self._extract_terms(pp2)
            similarity = self.psi.compute_private_jaccard(terms1, terms2)
        elif self.base_similarity:
            terms1 = list(self._extract_terms(pp1))
            terms2 = list(self._extract_terms(pp2))
            similarity = self.base_similarity.compute_pairwise_similarity(terms1, terms2)
        else:
            # Fallback: simple Jaccard
            terms1 = self._extract_terms(pp1)
            terms2 = self._extract_terms(pp2)
            if not terms1 and not terms2:
                similarity = 1.0
            elif not terms1 or not terms2:
                similarity = 0.0
            else:
                similarity = len(terms1 & terms2) / len(terms1 | terms2)

        # Step 3: Add DP noise
        if self._use_dp:
            similarity = self.dp_mechanism.privatize_similarity(similarity)
            self._total_epsilon_spent += self.config.epsilon

        # Cache result
        if self.cache_similarities and id1 and id2:
            self._cache[cache_key] = similarity

        return similarity

    def compute_private_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict,
        use_psi: bool = None,
        add_dp_noise: bool = None,
        filter_rare: bool = None
    ) -> float:
        """
        Compute similarity with fine-grained control over privacy mechanisms.

        Args:
            phenopacket1: First phenopacket
            phenopacket2: Second phenopacket
            use_psi: Override PSI setting (None = use default)
            add_dp_noise: Override DP setting (None = use default)
            filter_rare: Override rare filter setting (None = use default)

        Returns:
            Similarity score with specified privacy protections
        """
        # Temporarily override settings
        orig_psi = self._use_psi
        orig_dp = self._use_dp
        orig_filter = self._use_rare_filter

        if use_psi is not None:
            self._use_psi = use_psi
        if add_dp_noise is not None:
            self._use_dp = add_dp_noise
        if filter_rare is not None:
            self._use_rare_filter = filter_rare

        try:
            return self.compute_similarity(phenopacket1, phenopacket2)
        finally:
            # Restore settings
            self._use_psi = orig_psi
            self._use_dp = orig_dp
            self._use_rare_filter = orig_filter

    def find_most_similar(
        self,
        query_phenopacket: Dict,
        candidate_phenopackets: List[Dict],
        top_k: int = 10
    ) -> Optional[List[Tuple[int, float]]]:
        """
        Find most similar phenopackets with privacy protections.

        Applies all privacy mechanisms:
        1. Rare term filtering
        2. PSI or base similarity
        3. DP noise for ranking
        4. k-anonymity enforcement

        Args:
            query_phenopacket: Query phenopacket
            candidate_phenopackets: List of candidates
            top_k: Number of top results

        Returns:
            List of (index, similarity) tuples, or None if k-anonymity violated
        """
        # Compute all similarities
        results = []
        for idx, candidate in enumerate(candidate_phenopackets):
            sim = self.compute_similarity(query_phenopacket, candidate)
            results.append((idx, sim))

        # Sort by similarity (descending)
        results.sort(key=lambda x: x[1], reverse=True)
        results = results[:top_k]

        # Apply k-anonymity
        if self._use_k_anonymity:
            return self.k_guard.filter_results(results)

        return results

    def find_similar_patients(
        self,
        query: Dict,
        candidates: List[Dict],
        top_k: int = 10,
        threshold: Optional[float] = None
    ) -> Optional[List[Tuple[int, float]]]:
        """
        Privacy-preserving patient similarity search.

        Implements two-step reveal ladder:
        1. Threshold matching (optional)
        2. Top-k retrieval with DP noise
        3. k-anonymity enforcement

        Args:
            query: Query phenopacket
            candidates: List of candidates
            top_k: Number of results
            threshold: Optional minimum similarity threshold

        Returns:
            List of (index, similarity) tuples, or None if suppressed
        """
        # Compute similarities
        results = []
        for idx, candidate in enumerate(candidates):
            sim = self.compute_similarity(query, candidate)

            # Apply threshold if specified
            if threshold is None or sim >= threshold:
                results.append((idx, sim))

        # Sort by similarity
        results.sort(key=lambda x: x[1], reverse=True)
        results = results[:top_k]

        # Apply k-anonymity
        if self._use_k_anonymity:
            return self.k_guard.filter_results(results)

        return results

    def compute_similarity_matrix(
        self,
        phenopackets: List[Dict]
    ) -> np.ndarray:
        """
        Compute all-vs-all similarity matrix.

        Note: This reveals the full similarity structure.
        Consider privacy implications before using.

        Args:
            phenopackets: List of phenopackets

        Returns:
            NxN similarity matrix
        """
        n = len(phenopackets)
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    matrix[i, j] = 1.0
                else:
                    sim = self.compute_similarity(phenopackets[i], phenopackets[j])
                    matrix[i, j] = sim
                    matrix[j, i] = sim  # Symmetric

        return matrix

    def get_privacy_report(self) -> Dict[str, Any]:
        """
        Get comprehensive privacy accounting report.

        Returns:
            Dictionary with privacy statistics
        """
        report = {
            "total_queries": self._queries,
            "epsilon_per_query": self.config.epsilon if self._use_dp else 0,
            "total_epsilon_spent": self._total_epsilon_spent,
            "mechanisms_enabled": {
                "psi": self._use_psi,
                "dp": self._use_dp,
                "k_anonymity": self._use_k_anonymity,
                "rare_filter": self._use_rare_filter
            },
            "config": {
                "k": self.config.k,
                "epsilon": self.config.epsilon,
                "delta": self.config.delta,
                "min_prevalence": self.config.min_prevalence
            }
        }

        # Add component statistics
        if self._use_rare_filter:
            report["rare_filter_stats"] = self.rare_filter.get_statistics()

        if self._use_k_anonymity:
            report["k_anonymity_stats"] = self.k_guard.get_statistics()

        return report

    def reset_privacy_accounting(self):
        """Reset all privacy counters."""
        self._queries = 0
        self._total_epsilon_spent = 0.0
        self.rare_filter.reset_statistics()
        self.k_guard.reset_statistics()
        self._cache.clear()


class TwoStepRevealLadder:
    """
    Implements the two-step reveal ladder for privacy-preserving matching.

    Step 1: Threshold check using PSI
    - Determines if similarity exceeds threshold
    - No specific terms or similarity values revealed

    Step 2: Detailed matching (if step 1 passes)
    - Reveals only overlapping terms
    - Supports formal data sharing agreement initiation
    """

    def __init__(
        self,
        calculator: PrivacyPreservingCalculator,
        initial_threshold: float = 0.3,
        detailed_threshold: float = 0.5
    ):
        """
        Initialize two-step reveal ladder.

        Args:
            calculator: Privacy-preserving calculator
            initial_threshold: Threshold for step 1 (existence check)
            detailed_threshold: Threshold for step 2 (detailed reveal)
        """
        self.calculator = calculator
        self.initial_threshold = initial_threshold
        self.detailed_threshold = detailed_threshold

    def step1_existence_check(
        self,
        query: Dict,
        candidates: List[Dict]
    ) -> List[int]:
        """
        Step 1: Check which candidates might be relevant.

        Only returns indices, no similarity values or terms.

        Args:
            query: Query phenopacket
            candidates: Candidate phenopackets

        Returns:
            List of candidate indices exceeding initial threshold
        """
        matching_indices = []

        for idx, candidate in enumerate(candidates):
            # Use PSI threshold matching
            if self.calculator.psi:
                query_terms = self.calculator._extract_terms(query)
                cand_terms = self.calculator._extract_terms(candidate)

                if self.calculator.psi.threshold_match(
                    query_terms, cand_terms, self.initial_threshold
                ):
                    matching_indices.append(idx)
            else:
                # Fallback to similarity computation
                sim = self.calculator.compute_similarity(query, candidate)
                if sim >= self.initial_threshold:
                    matching_indices.append(idx)

        return matching_indices

    def step2_detailed_match(
        self,
        query: Dict,
        candidate: Dict
    ) -> Optional[Dict]:
        """
        Step 2: Detailed matching with controlled disclosure.

        Only proceeds if similarity exceeds detailed threshold.

        Args:
            query: Query phenopacket
            candidate: Single candidate phenopacket

        Returns:
            Match details if above threshold, None otherwise
        """
        # Compute similarity
        similarity = self.calculator.compute_similarity(query, candidate)

        if similarity < self.detailed_threshold:
            return None

        # Extract overlapping terms (controlled disclosure)
        query_terms = self.calculator._extract_terms(query)
        cand_terms = self.calculator._extract_terms(candidate)
        overlapping = query_terms & cand_terms

        return {
            "candidate_id": candidate.get("id"),
            "similarity_score": similarity,
            "num_overlapping_terms": len(overlapping),
            "overlapping_terms": list(overlapping),
            "ready_for_data_sharing_agreement": True
        }

    def full_workflow(
        self,
        query: Dict,
        candidates: List[Dict]
    ) -> Dict:
        """
        Execute full two-step workflow.

        Args:
            query: Query phenopacket
            candidates: Candidate phenopackets

        Returns:
            Workflow results
        """
        # Step 1
        step1_matches = self.step1_existence_check(query, candidates)

        if not step1_matches:
            return {
                "step1_matches": 0,
                "step2_matches": [],
                "workflow_complete": True
            }

        # Step 2 for each match
        step2_results = []
        for idx in step1_matches:
            result = self.step2_detailed_match(query, candidates[idx])
            if result:
                result["candidate_index"] = idx
                step2_results.append(result)

        return {
            "step1_matches": len(step1_matches),
            "step2_matches": step2_results,
            "workflow_complete": True
        }


def create_privacy_calculator(
    config_path: Optional[str] = None,
    base_similarity=None,
    hpo_manager=None,
    **kwargs
) -> PrivacyPreservingCalculator:
    """
    Factory function to create privacy-preserving calculator.

    Args:
        config_path: Path to config.yaml
        base_similarity: Base HPOSimilarity metric
        hpo_manager: HPOManager instance
        **kwargs: Override specific privacy settings

    Returns:
        Configured PrivacyPreservingCalculator
    """
    # Load config
    config = load_privacy_config(config_path)

    # Override with kwargs
    for key, value in kwargs.items():
        if hasattr(config, key):
            setattr(config, key, value)

    return PrivacyPreservingCalculator(
        base_similarity=base_similarity,
        config=config,
        hpo_manager=hpo_manager
    )


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Create sample phenopackets
    sample_phenopackets = [
        {
            "id": f"patient_{i}",
            "phenotypicFeatures": [
                {"type": {"id": f"HP:000000{j}"}, "excluded": False}
                for j in range(1, 6 + i % 3)
            ]
        }
        for i in range(10)
    ]

    # Create privacy-preserving calculator
    calc = create_privacy_calculator(
        epsilon=1.0,
        k=3,
        min_prevalence=0.05
    )

    # Set corpus for prevalence computation
    calc.set_corpus(sample_phenopackets)

    # Test similarity computation
    sim = calc.compute_similarity(sample_phenopackets[0], sample_phenopackets[1])
    print(f"Privacy-preserving similarity: {sim:.4f}")

    # Test find_most_similar
    results = calc.find_most_similar(
        sample_phenopackets[0],
        sample_phenopackets[1:],
        top_k=5
    )
    print(f"\nTop 5 similar patients: {results}")

    # Get privacy report
    report = calc.get_privacy_report()
    print("\nPrivacy Report:")
    for key, value in report.items():
        print(f"  {key}: {value}")

    # Test two-step reveal ladder
    print("\nTwo-step reveal ladder:")
    ladder = TwoStepRevealLadder(calc, initial_threshold=0.1, detailed_threshold=0.2)
    workflow_result = ladder.full_workflow(sample_phenopackets[0], sample_phenopackets[1:])
    print(f"Step 1 matches: {workflow_result['step1_matches']}")
    print(f"Step 2 detailed matches: {len(workflow_result['step2_matches'])}")
