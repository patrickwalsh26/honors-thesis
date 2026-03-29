"""
Integration tests for privacy-preserving calculator.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.privacy.privacy_calculator import (
    PrivacyPreservingCalculator,
    TwoStepRevealLadder,
    create_privacy_calculator
)
from src.privacy.k_anonymity import PrivacyConfig
from src.similarity.hpo_similarity import CosineSimilarity


class TestPrivacyPreservingCalculator:
    """Tests for PrivacyPreservingCalculator."""

    def test_initialization(self, privacy_config, ic_values):
        """Test initialization with various configurations."""
        metric = CosineSimilarity(ic_values)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=privacy_config,
            use_psi=False,
            use_dp=True,
            use_k_anonymity=True,
            use_rare_filter=False
        )

        assert calc._use_dp
        assert calc._use_k_anonymity
        assert not calc._use_rare_filter

    def test_compute_similarity(self, small_phenopackets, ic_values):
        """Test basic similarity computation."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=1.0, k=1)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=False,
            use_dp=True,
            use_k_anonymity=False
        )

        sim = calc.compute_similarity(
            small_phenopackets[0],
            small_phenopackets[1]
        )

        assert 0.0 <= sim <= 1.0

    def test_similarity_with_psi(self, small_phenopackets):
        """Test similarity with PSI enabled."""
        config = PrivacyConfig(epsilon=float('inf'), k=1)

        calc = PrivacyPreservingCalculator(
            config=config,
            use_psi=True,
            use_dp=False,
            use_k_anonymity=False,
            use_rare_filter=False
        )

        # Patient A and B share HP:0002345
        sim_ab = calc.compute_similarity(
            small_phenopackets[0],
            small_phenopackets[1]
        )

        # Patient A and C share no terms
        sim_ac = calc.compute_similarity(
            small_phenopackets[0],
            small_phenopackets[2]
        )

        assert sim_ab > sim_ac

    def test_find_most_similar(self, sample_phenopackets, ic_values):
        """Test finding most similar phenopackets."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=1.0, k=3)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=False,
            use_dp=True,
            use_k_anonymity=True
        )

        results = calc.find_most_similar(
            sample_phenopackets[0],
            sample_phenopackets,
            top_k=5
        )

        # Results should be list of tuples or None
        if results is not None:
            assert len(results) <= 5
            assert all(isinstance(r, tuple) for r in results)

    def test_k_anonymity_suppression(self, small_phenopackets, ic_values):
        """Test that k-anonymity suppresses small result sets."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=float('inf'), k=10)  # High k

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=False,
            use_dp=False,
            use_k_anonymity=True
        )

        results = calc.find_most_similar(
            small_phenopackets[0],
            small_phenopackets,  # Only 3 patients
            top_k=5
        )

        # Should be suppressed
        assert results is None

    def test_rare_term_filtering(self, sample_phenopackets, ic_values):
        """Test rare term filtering."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(
            epsilon=float('inf'),
            k=1,
            min_prevalence=0.5  # High threshold
        )

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=False,
            use_dp=False,
            use_k_anonymity=False,
            use_rare_filter=True
        )

        # Set corpus
        calc.set_corpus(sample_phenopackets)

        # Compute similarity
        sim = calc.compute_similarity(
            sample_phenopackets[0],
            sample_phenopackets[1]
        )

        # Should still return valid similarity
        assert 0.0 <= sim <= 1.0

    def test_privacy_report(self, sample_phenopackets, ic_values):
        """Test privacy report generation."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=1.0, k=5)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=False,
            use_dp=True,
            use_k_anonymity=True
        )

        # Run some queries
        for pp in sample_phenopackets[:3]:
            calc.compute_similarity(pp, sample_phenopackets[0])

        report = calc.get_privacy_report()

        assert "total_queries" in report
        assert "epsilon_per_query" in report
        assert "mechanisms_enabled" in report

    def test_compute_private_similarity_overrides(self, small_phenopackets, ic_values):
        """Test fine-grained control over privacy mechanisms."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=1.0, k=1)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=False,
            use_dp=True,
            use_k_anonymity=False
        )

        # Compute with DP disabled
        sim_no_dp = calc.compute_private_similarity(
            small_phenopackets[0],
            small_phenopackets[1],
            add_dp_noise=False
        )

        # Multiple calls should give same result
        sim_no_dp_2 = calc.compute_private_similarity(
            small_phenopackets[0],
            small_phenopackets[1],
            add_dp_noise=False
        )

        # Without noise, should be deterministic (cache aside)
        # Actually with cache it should be exactly the same
        assert sim_no_dp == sim_no_dp_2


class TestTwoStepRevealLadder:
    """Tests for TwoStepRevealLadder."""

    def test_step1_existence_check(self, sample_phenopackets, ic_values):
        """Test step 1 existence check."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=float('inf'), k=1)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=True,
            use_dp=False,
            use_k_anonymity=False
        )

        ladder = TwoStepRevealLadder(
            calculator=calc,
            initial_threshold=0.1,
            detailed_threshold=0.2
        )

        # Get same-disease patients (should match)
        matches = ladder.step1_existence_check(
            sample_phenopackets[0],
            sample_phenopackets
        )

        assert len(matches) > 0

    def test_step2_detailed_match(self, small_phenopackets, ic_values):
        """Test step 2 detailed match."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=float('inf'), k=1)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=True,
            use_dp=False,
            use_k_anonymity=False
        )

        ladder = TwoStepRevealLadder(
            calculator=calc,
            initial_threshold=0.0,
            detailed_threshold=0.0
        )

        result = ladder.step2_detailed_match(
            small_phenopackets[0],
            small_phenopackets[1]
        )

        assert result is not None
        assert "similarity_score" in result
        assert "overlapping_terms" in result

    def test_full_workflow(self, sample_phenopackets, ic_values):
        """Test complete workflow."""
        metric = CosineSimilarity(ic_values)
        config = PrivacyConfig(epsilon=float('inf'), k=1)

        calc = PrivacyPreservingCalculator(
            base_similarity=metric,
            config=config,
            use_psi=True,
            use_dp=False,
            use_k_anonymity=False
        )

        ladder = TwoStepRevealLadder(
            calculator=calc,
            initial_threshold=0.1,
            detailed_threshold=0.2
        )

        result = ladder.full_workflow(
            sample_phenopackets[0],
            sample_phenopackets
        )

        assert "step1_matches" in result
        assert "step2_matches" in result
        assert "workflow_complete" in result


class TestCreatePrivacyCalculator:
    """Tests for factory function."""

    def test_default_creation(self, ic_values):
        """Test creating with defaults."""
        metric = CosineSimilarity(ic_values)

        calc = create_privacy_calculator(
            base_similarity=metric,
            epsilon=1.0,
            k=5
        )

        assert calc.config.epsilon == 1.0
        assert calc.config.k == 5

    def test_with_overrides(self, ic_values):
        """Test creating with parameter overrides."""
        metric = CosineSimilarity(ic_values)

        calc = create_privacy_calculator(
            base_similarity=metric,
            epsilon=0.5,
            k=10,
            min_prevalence=0.01
        )

        assert calc.config.epsilon == 0.5
        assert calc.config.k == 10
        assert calc.config.min_prevalence == 0.01
