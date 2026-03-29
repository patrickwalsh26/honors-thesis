"""
k-Anonymity and rare term filtering for privacy-preserving phenotype matching.

Implements:
- RareTermFilter: Suppress or generalize rare HPO terms
- KAnonymityGuard: Enforce minimum cohort sizes
- PrivacyConfig: Load privacy parameters from config
"""

from pathlib import Path

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    yaml = None
from typing import Dict, List, Set, Tuple, Optional
from collections import Counter
from dataclasses import dataclass
import logging
import copy

logger = logging.getLogger(__name__)


@dataclass
class PrivacyConfig:
    """Configuration container for privacy parameters."""

    # k-anonymity
    k: int = 5

    # Rare term filtering
    min_prevalence: float = 0.001
    rare_term_strategy: str = "suppress"  # "suppress", "generalize"

    # Differential privacy
    epsilon: float = 1.0
    delta: float = 1e-5
    dp_mechanism: str = "laplace"

    # PSI
    psi_enabled: bool = True
    psi_security_parameter: int = 128

    @classmethod
    def from_yaml(cls, config_path: str) -> "PrivacyConfig":
        """
        Load configuration from YAML file.

        Args:
            config_path: Path to config.yaml

        Returns:
            PrivacyConfig instance
        """
        if not YAML_AVAILABLE:
            logger.warning("PyYAML not available, using default config")
            return cls.default()

        with open(config_path) as f:
            config = yaml.safe_load(f)

        privacy = config.get("privacy", {})
        dp = privacy.get("differential_privacy", {})
        k_anon = privacy.get("k_anonymity", {})
        rare = privacy.get("rare_term_filter", {})
        psi = privacy.get("psi", {})

        return cls(
            k=k_anon.get("k", 5),
            min_prevalence=rare.get("min_prevalence", 0.001),
            epsilon=dp.get("epsilon", 1.0),
            delta=dp.get("delta", 1e-5),
            dp_mechanism=dp.get("mechanism", "laplace"),
            psi_enabled=psi.get("enabled", True),
            psi_security_parameter=psi.get("security_parameter", 128)
        )

    @classmethod
    def default(cls) -> "PrivacyConfig":
        """Return default configuration."""
        return cls()


class RareTermFilter:
    """
    Filter rare HPO terms that could be re-identifying.

    Rare phenotypes (low prevalence) can act as quasi-identifiers,
    potentially revealing patient identity when combined with other
    information. This class implements strategies to mitigate this risk.

    Strategies:
    - suppress: Remove rare terms entirely
    - generalize: Replace rare terms with more general parent terms
    """

    def __init__(
        self,
        min_prevalence: float = 0.001,
        term_prevalence: Optional[Dict[str, float]] = None,
        hpo_manager=None
    ):
        """
        Initialize rare term filter.

        Args:
            min_prevalence: Minimum term prevalence threshold
            term_prevalence: Pre-computed term prevalence values
            hpo_manager: HPOManager instance for term generalization
        """
        self.min_prevalence = min_prevalence
        self.term_prevalence = term_prevalence or {}
        self.hpo_manager = hpo_manager

        # Statistics
        self._suppressed_count = 0
        self._generalized_count = 0
        self._total_filtered = 0

    def compute_prevalence_from_corpus(
        self,
        phenopackets: List[Dict]
    ) -> Dict[str, float]:
        """
        Compute term prevalence from a phenopacket corpus.

        Prevalence = fraction of patients with the term.

        Args:
            phenopackets: List of phenopacket dictionaries

        Returns:
            Dictionary mapping term IDs to prevalence values
        """
        term_counts = Counter()
        n_patients = len(phenopackets)

        if n_patients == 0:
            return {}

        for pp in phenopackets:
            patient_terms = set()
            for feature in pp.get("phenotypicFeatures", []):
                if not feature.get("excluded", False):
                    term_id = feature["type"]["id"]
                    patient_terms.add(term_id)

            for term in patient_terms:
                term_counts[term] += 1

        self.term_prevalence = {
            term: count / n_patients
            for term, count in term_counts.items()
        }

        logger.info(f"Computed prevalence for {len(self.term_prevalence)} terms")
        return self.term_prevalence

    def is_rare(self, term_id: str) -> bool:
        """
        Check if a term is below the prevalence threshold.

        Args:
            term_id: HPO term ID

        Returns:
            True if term is rare (below threshold)
        """
        prevalence = self.term_prevalence.get(term_id, 0.0)
        return prevalence < self.min_prevalence

    def get_rare_terms(self, terms: List[str]) -> List[str]:
        """
        Get list of rare terms from a term list.

        Args:
            terms: List of HPO term IDs

        Returns:
            Subset of terms that are rare
        """
        return [t for t in terms if self.is_rare(t)]

    def filter_terms(
        self,
        terms: List[str],
        strategy: str = "suppress"
    ) -> List[str]:
        """
        Filter rare terms using the specified strategy.

        Args:
            terms: List of HPO term IDs
            strategy: Filtering strategy ("suppress" or "generalize")

        Returns:
            Filtered list of terms
        """
        if strategy == "suppress":
            return self._suppress_rare_terms(terms)
        elif strategy == "generalize":
            return self._generalize_rare_terms(terms)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

    def _suppress_rare_terms(self, terms: List[str]) -> List[str]:
        """Remove rare terms entirely."""
        filtered = [t for t in terms if not self.is_rare(t)]
        suppressed = len(terms) - len(filtered)
        self._suppressed_count += suppressed
        self._total_filtered += suppressed
        return filtered

    def _generalize_rare_terms(self, terms: List[str]) -> List[str]:
        """Replace rare terms with non-rare ancestors."""
        if not self.hpo_manager:
            logger.warning("No HPO manager; falling back to suppression")
            return self._suppress_rare_terms(terms)

        result = []
        for term in terms:
            if self.is_rare(term):
                # Find first non-rare ancestor
                replacement = self._find_non_rare_ancestor(term)
                if replacement:
                    result.append(replacement)
                    self._generalized_count += 1
                else:
                    # No non-rare ancestor found; suppress
                    self._suppressed_count += 1
                self._total_filtered += 1
            else:
                result.append(term)

        # Remove duplicates while preserving order
        seen = set()
        unique_result = []
        for t in result:
            if t not in seen:
                seen.add(t)
                unique_result.append(t)

        return unique_result

    def _find_non_rare_ancestor(self, term_id: str) -> Optional[str]:
        """Find the most specific non-rare ancestor of a term."""
        if not self.hpo_manager:
            return None

        # Get ancestors ordered by depth (most specific first)
        ancestors = self.hpo_manager.get_ancestors(term_id, include_self=False)

        # Sort by IC (higher IC = more specific = preferred)
        def get_ic(t):
            return self.hpo_manager.ic_values.get(t, 0.0) if hasattr(self.hpo_manager, 'ic_values') else 0.0

        sorted_ancestors = sorted(ancestors, key=get_ic, reverse=True)

        for ancestor in sorted_ancestors:
            if not self.is_rare(ancestor):
                return ancestor

        return None

    def filter_phenopacket(
        self,
        phenopacket: Dict,
        strategy: str = "suppress"
    ) -> Dict:
        """
        Return a filtered copy of a phenopacket.

        Args:
            phenopacket: Original phenopacket dictionary
            strategy: Filtering strategy

        Returns:
            New phenopacket with rare terms handled
        """
        # Deep copy to avoid modifying original
        filtered = copy.deepcopy(phenopacket)

        # Get current terms
        original_terms = []
        for feature in filtered.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                original_terms.append(feature["type"]["id"])

        # Filter terms
        filtered_terms = set(self.filter_terms(original_terms, strategy))

        # Update phenotypic features
        new_features = []
        for feature in filtered.get("phenotypicFeatures", []):
            term_id = feature["type"]["id"]
            if feature.get("excluded", False):
                # Keep excluded/negated terms
                new_features.append(feature)
            elif term_id in filtered_terms:
                # Keep non-rare observed terms
                new_features.append(feature)
            # Rare observed terms are removed

        filtered["phenotypicFeatures"] = new_features

        return filtered

    def get_statistics(self) -> Dict[str, int]:
        """Get filtering statistics."""
        return {
            "suppressed": self._suppressed_count,
            "generalized": self._generalized_count,
            "total_filtered": self._total_filtered
        }

    def reset_statistics(self):
        """Reset filtering counters."""
        self._suppressed_count = 0
        self._generalized_count = 0
        self._total_filtered = 0


class KAnonymityGuard:
    """
    k-Anonymity enforcement for query results.

    Ensures that query results never identify fewer than k patients,
    preventing re-identification through small cohort sizes.

    If fewer than k matches are found, the result is either:
    - Suppressed entirely (returns None)
    - Padded with dummy entries (less private)
    """

    def __init__(self, k: int = 5):
        """
        Initialize k-anonymity guard.

        Args:
            k: Minimum number of records to return
        """
        if k < 1:
            raise ValueError("k must be at least 1")
        self.k = k

        # Statistics
        self._total_queries = 0
        self._suppressed_queries = 0
        self._passed_queries = 0

    def check_anonymity(self, result_count: int) -> bool:
        """
        Check if result set satisfies k-anonymity.

        Args:
            result_count: Number of matching records

        Returns:
            True if k-anonymity is satisfied
        """
        return result_count >= self.k

    def filter_results(
        self,
        results: List[Tuple[int, float]],
        suppress_if_insufficient: bool = True
    ) -> Optional[List[Tuple[int, float]]]:
        """
        Filter results to enforce k-anonymity.

        Args:
            results: List of (index, score) tuples
            suppress_if_insufficient: If True, return None when < k results

        Returns:
            Filtered results, or None if suppressed
        """
        self._total_queries += 1

        if len(results) < self.k:
            if suppress_if_insufficient:
                self._suppressed_queries += 1
                logger.debug(f"Query suppressed: {len(results)} < k={self.k}")
                return None
            else:
                # Pad with dummy results (index=-1 indicates dummy)
                while len(results) < self.k:
                    results.append((-1, 0.0))

        self._passed_queries += 1
        return results

    def filter_count(
        self,
        count: int,
        suppress_value: int = 0
    ) -> int:
        """
        Filter a count query for k-anonymity.

        If count < k, return suppress_value instead.

        Args:
            count: True count
            suppress_value: Value to return if suppressed

        Returns:
            Original count or suppress_value
        """
        self._total_queries += 1

        if count < self.k:
            self._suppressed_queries += 1
            return suppress_value

        self._passed_queries += 1
        return count

    def get_suppression_rate(self) -> float:
        """
        Get fraction of queries that were suppressed.

        Returns:
            Suppression rate in [0, 1]
        """
        if self._total_queries == 0:
            return 0.0
        return self._suppressed_queries / self._total_queries

    def get_statistics(self) -> Dict[str, any]:
        """Get query statistics."""
        return {
            "k": self.k,
            "total_queries": self._total_queries,
            "suppressed_queries": self._suppressed_queries,
            "passed_queries": self._passed_queries,
            "suppression_rate": self.get_suppression_rate()
        }

    def reset_statistics(self):
        """Reset query counters."""
        self._total_queries = 0
        self._suppressed_queries = 0
        self._passed_queries = 0


class EquivalenceClassAnalyzer:
    """
    Analyze phenopacket cohort for k-anonymity compliance.

    Groups phenopackets into equivalence classes based on quasi-identifiers
    and identifies records that violate k-anonymity.
    """

    def __init__(self, k: int = 5):
        """
        Initialize analyzer.

        Args:
            k: k-anonymity parameter
        """
        self.k = k

    def group_by_equivalence_class(
        self,
        phenopackets: List[Dict],
        quasi_identifiers: List[str]
    ) -> Dict[Tuple, List[int]]:
        """
        Group phenopackets into equivalence classes.

        Args:
            phenopackets: List of phenopacket dictionaries
            quasi_identifiers: List of quasi-identifier field names

        Returns:
            Dictionary mapping QI tuples to list of phenopacket indices
        """
        classes = {}

        for idx, pp in enumerate(phenopackets):
            qi_values = []

            for qi in quasi_identifiers:
                if qi == "disease":
                    diseases = pp.get("diseases", [])
                    val = diseases[0]["term"]["label"] if diseases else "Unknown"
                elif qi == "sex":
                    val = pp.get("subject", {}).get("sex", "UNKNOWN_SEX")
                elif qi == "n_features":
                    features = pp.get("phenotypicFeatures", [])
                    observed = [f for f in features if not f.get("excluded", False)]
                    val = len(observed)
                elif qi == "n_features_bucket":
                    features = pp.get("phenotypicFeatures", [])
                    observed = [f for f in features if not f.get("excluded", False)]
                    n = len(observed)
                    # Bucket into ranges: 1-5, 6-10, 11-15, 16+
                    if n <= 5:
                        val = "1-5"
                    elif n <= 10:
                        val = "6-10"
                    elif n <= 15:
                        val = "11-15"
                    else:
                        val = "16+"
                else:
                    val = "Unknown"

                qi_values.append(val)

            key = tuple(qi_values)
            if key not in classes:
                classes[key] = []
            classes[key].append(idx)

        return classes

    def find_unsafe_records(
        self,
        phenopackets: List[Dict],
        quasi_identifiers: List[str]
    ) -> List[int]:
        """
        Find records that violate k-anonymity.

        Args:
            phenopackets: List of phenopacket dictionaries
            quasi_identifiers: List of quasi-identifier field names

        Returns:
            List of indices of unsafe records
        """
        classes = self.group_by_equivalence_class(phenopackets, quasi_identifiers)

        unsafe = []
        for indices in classes.values():
            if len(indices) < self.k:
                unsafe.extend(indices)

        return unsafe

    def analyze_cohort(
        self,
        phenopackets: List[Dict],
        quasi_identifiers: List[str]
    ) -> Dict:
        """
        Comprehensive k-anonymity analysis of a cohort.

        Args:
            phenopackets: List of phenopacket dictionaries
            quasi_identifiers: List of quasi-identifier field names

        Returns:
            Analysis report dictionary
        """
        classes = self.group_by_equivalence_class(phenopackets, quasi_identifiers)

        class_sizes = [len(indices) for indices in classes.values()]
        unsafe_classes = sum(1 for size in class_sizes if size < self.k)
        unsafe_records = sum(size for size in class_sizes if size < self.k)

        return {
            "k": self.k,
            "quasi_identifiers": quasi_identifiers,
            "total_records": len(phenopackets),
            "num_equivalence_classes": len(classes),
            "min_class_size": min(class_sizes) if class_sizes else 0,
            "max_class_size": max(class_sizes) if class_sizes else 0,
            "mean_class_size": sum(class_sizes) / len(class_sizes) if class_sizes else 0,
            "unsafe_classes": unsafe_classes,
            "unsafe_records": unsafe_records,
            "compliance_rate": 1.0 - (unsafe_records / len(phenopackets)) if phenopackets else 1.0,
            "is_k_anonymous": unsafe_classes == 0
        }


def load_privacy_config(config_path: Optional[str] = None) -> PrivacyConfig:
    """
    Load privacy configuration.

    Args:
        config_path: Path to config.yaml, or None for default

    Returns:
        PrivacyConfig instance
    """
    if config_path and Path(config_path).exists():
        return PrivacyConfig.from_yaml(config_path)
    else:
        # Try default location
        default_path = Path("config.yaml")
        if default_path.exists():
            return PrivacyConfig.from_yaml(str(default_path))

    logger.info("Using default privacy configuration")
    return PrivacyConfig.default()


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Create sample phenopackets
    sample_phenopackets = [
        {
            "id": f"patient_{i}",
            "subject": {"sex": "MALE" if i % 2 == 0 else "FEMALE"},
            "phenotypicFeatures": [
                {"type": {"id": f"HP:000000{j}"}, "excluded": False}
                for j in range(1, 6)
            ],
            "diseases": [{"term": {"label": "Disease A" if i < 5 else "Disease B"}}]
        }
        for i in range(10)
    ]

    # Test RareTermFilter
    filter = RareTermFilter(min_prevalence=0.2)
    filter.compute_prevalence_from_corpus(sample_phenopackets)

    print("Term prevalence:")
    for term, prev in filter.term_prevalence.items():
        print(f"  {term}: {prev:.3f}")

    # Test KAnonymityGuard
    guard = KAnonymityGuard(k=3)

    results_small = [(i, 0.9 - i * 0.1) for i in range(2)]
    results_large = [(i, 0.9 - i * 0.1) for i in range(5)]

    filtered_small = guard.filter_results(results_small)
    filtered_large = guard.filter_results(results_large)

    print(f"\nSmall result set: {filtered_small}")
    print(f"Large result set: {filtered_large}")
    print(f"Suppression rate: {guard.get_suppression_rate():.2%}")

    # Test EquivalenceClassAnalyzer
    analyzer = EquivalenceClassAnalyzer(k=3)
    report = analyzer.analyze_cohort(
        sample_phenopackets,
        quasi_identifiers=["disease", "sex"]
    )

    print("\nk-Anonymity Analysis:")
    for key, value in report.items():
        print(f"  {key}: {value}")
