"""
Tests for k-anonymity and rare term filtering.
"""

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.privacy.k_anonymity import (
    RareTermFilter,
    KAnonymityGuard,
    EquivalenceClassAnalyzer,
    PrivacyConfig,
    load_privacy_config
)


class TestRareTermFilter:
    """Tests for RareTermFilter."""

    def test_compute_prevalence(self, sample_phenopackets):
        """Test prevalence computation from corpus."""
        filter = RareTermFilter(min_prevalence=0.1)
        prevalence = filter.compute_prevalence_from_corpus(sample_phenopackets)

        assert len(prevalence) > 0
        assert all(0 <= p <= 1 for p in prevalence.values())

    def test_is_rare(self, sample_phenopackets):
        """Test rare term detection."""
        filter = RareTermFilter(min_prevalence=0.5)
        filter.compute_prevalence_from_corpus(sample_phenopackets)

        # High threshold should mark most terms as rare
        rare_count = sum(1 for term in filter.term_prevalence if filter.is_rare(term))
        assert rare_count > 0

    def test_filter_terms_suppress(self, sample_phenopackets):
        """Test term suppression."""
        filter = RareTermFilter(min_prevalence=0.8)  # Very high threshold
        filter.compute_prevalence_from_corpus(sample_phenopackets)

        terms = ["HP:0001166", "HP:9999999"]  # One common, one unknown
        filtered = filter.filter_terms(terms, strategy="suppress")

        # Unknown term should be suppressed
        assert "HP:9999999" not in filtered

    def test_filter_phenopacket(self, sample_phenopackets):
        """Test phenopacket filtering."""
        filter = RareTermFilter(min_prevalence=0.8)
        filter.compute_prevalence_from_corpus(sample_phenopackets)

        original = sample_phenopackets[0]
        filtered = filter.filter_phenopacket(original, strategy="suppress")

        # Filtered should have same or fewer terms
        assert len(filtered["phenotypicFeatures"]) <= len(original["phenotypicFeatures"])

        # Original should be unchanged
        assert filtered is not original
        assert filtered["id"] == original["id"]

    def test_statistics(self, sample_phenopackets):
        """Test statistics tracking."""
        filter = RareTermFilter(min_prevalence=0.5)
        filter.compute_prevalence_from_corpus(sample_phenopackets)

        # Filter some phenopackets
        for pp in sample_phenopackets[:5]:
            filter.filter_phenopacket(pp)

        stats = filter.get_statistics()
        assert "total_filtered" in stats
        assert stats["total_filtered"] >= 0


class TestKAnonymityGuard:
    """Tests for KAnonymityGuard."""

    def test_initialization(self):
        """Test valid initialization."""
        guard = KAnonymityGuard(k=5)
        assert guard.k == 5

    def test_invalid_k(self):
        """Test that k < 1 raises error."""
        with pytest.raises(ValueError):
            KAnonymityGuard(k=0)

    def test_check_anonymity(self):
        """Test anonymity checking."""
        guard = KAnonymityGuard(k=5)

        assert not guard.check_anonymity(3)
        assert guard.check_anonymity(5)
        assert guard.check_anonymity(10)

    def test_filter_results_suppression(self):
        """Test result suppression when < k."""
        guard = KAnonymityGuard(k=5)

        small_results = [(i, 0.9 - i * 0.1) for i in range(3)]
        filtered = guard.filter_results(small_results, suppress_if_insufficient=True)

        assert filtered is None

    def test_filter_results_pass(self):
        """Test results pass through when >= k."""
        guard = KAnonymityGuard(k=5)

        large_results = [(i, 0.9 - i * 0.1) for i in range(10)]
        filtered = guard.filter_results(large_results)

        assert filtered is not None
        assert len(filtered) == 10

    def test_filter_count(self):
        """Test count filtering."""
        guard = KAnonymityGuard(k=5)

        assert guard.filter_count(3) == 0
        assert guard.filter_count(10) == 10

    def test_suppression_rate(self):
        """Test suppression rate calculation."""
        guard = KAnonymityGuard(k=5)

        # Run some queries
        guard.filter_results([(0, 0.5)] * 3)  # Suppressed
        guard.filter_results([(i, 0.5) for i in range(10)])  # Passed
        guard.filter_results([(0, 0.5)] * 2)  # Suppressed

        rate = guard.get_suppression_rate()
        assert rate == 2/3  # 2 suppressed out of 3


class TestEquivalenceClassAnalyzer:
    """Tests for EquivalenceClassAnalyzer."""

    def test_group_by_disease(self, sample_phenopackets):
        """Test grouping by disease."""
        analyzer = EquivalenceClassAnalyzer(k=5)

        classes = analyzer.group_by_equivalence_class(
            sample_phenopackets,
            quasi_identifiers=["disease"]
        )

        # Should have 4 diseases (from fixture)
        assert len(classes) == 4

    def test_find_unsafe_records(self, sample_phenopackets):
        """Test finding unsafe records."""
        analyzer = EquivalenceClassAnalyzer(k=15)  # High k

        unsafe = analyzer.find_unsafe_records(
            sample_phenopackets,
            quasi_identifiers=["disease"]
        )

        # With 10 per disease and k=15, all should be unsafe
        assert len(unsafe) == len(sample_phenopackets)

    def test_analyze_cohort(self, sample_phenopackets):
        """Test cohort analysis."""
        analyzer = EquivalenceClassAnalyzer(k=5)

        report = analyzer.analyze_cohort(
            sample_phenopackets,
            quasi_identifiers=["disease"]
        )

        assert "k" in report
        assert "total_records" in report
        assert "num_equivalence_classes" in report
        assert "is_k_anonymous" in report


class TestPrivacyConfig:
    """Tests for PrivacyConfig."""

    def test_default(self):
        """Test default configuration."""
        config = PrivacyConfig.default()

        assert config.k == 5
        assert config.epsilon == 1.0
        assert config.min_prevalence == 0.001

    def test_custom_values(self):
        """Test custom configuration."""
        config = PrivacyConfig(
            k=10,
            epsilon=0.5,
            min_prevalence=0.01
        )

        assert config.k == 10
        assert config.epsilon == 0.5
        assert config.min_prevalence == 0.01
