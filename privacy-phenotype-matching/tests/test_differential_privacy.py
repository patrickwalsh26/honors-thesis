"""
Tests for differential privacy mechanisms.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.privacy.differential_privacy import (
    LaplaceMechanism,
    GaussianMechanism,
    ExponentialMechanism,
    ReportNoisyMax,
    PrivacyAccountant,
    create_dp_mechanism
)


class TestLaplaceMechanism:
    """Tests for Laplace mechanism."""

    def test_initialization(self):
        """Test valid initialization."""
        mech = LaplaceMechanism(epsilon=1.0)
        assert mech.epsilon == 1.0
        assert mech.delta == 0.0

    def test_invalid_epsilon(self):
        """Test that invalid epsilon raises error."""
        with pytest.raises(ValueError):
            LaplaceMechanism(epsilon=0.0)
        with pytest.raises(ValueError):
            LaplaceMechanism(epsilon=-1.0)

    def test_noise_scale(self):
        """Test that noise scales inversely with epsilon."""
        mech_high = LaplaceMechanism(epsilon=10.0)
        mech_low = LaplaceMechanism(epsilon=0.1)

        # Generate samples
        n_samples = 1000
        samples_high = [mech_high.add_noise(0, 1.0) for _ in range(n_samples)]
        samples_low = [mech_low.add_noise(0, 1.0) for _ in range(n_samples)]

        # Low epsilon should have higher variance
        assert np.var(samples_low) > np.var(samples_high) * 10

    def test_unbiased(self):
        """Test that noise is unbiased (mean ≈ true value)."""
        mech = LaplaceMechanism(epsilon=1.0)
        true_value = 0.5

        samples = [mech.add_noise(true_value, 1.0) for _ in range(1000)]
        mean = np.mean(samples)

        # Mean should be close to true value
        assert abs(mean - true_value) < 0.1

    def test_privatize_similarity_bounds(self):
        """Test that privatized similarity stays in [0, 1]."""
        mech = LaplaceMechanism(epsilon=0.5)

        for _ in range(100):
            priv = mech.privatize_similarity(0.5)
            assert 0.0 <= priv <= 1.0

    def test_privatize_count(self):
        """Test count privatization."""
        mech = LaplaceMechanism(epsilon=1.0)
        true_count = 10

        # Should be non-negative
        for _ in range(100):
            priv = mech.privatize_count(true_count)
            assert priv >= 0


class TestGaussianMechanism:
    """Tests for Gaussian mechanism."""

    def test_initialization(self):
        """Test valid initialization."""
        mech = GaussianMechanism(epsilon=1.0, delta=1e-5)
        assert mech.epsilon == 1.0
        assert mech.delta == 1e-5

    def test_requires_delta(self):
        """Test that Gaussian requires delta > 0."""
        with pytest.raises(ValueError):
            GaussianMechanism(epsilon=1.0, delta=0.0)

    def test_noise_distribution(self):
        """Test that noise is normally distributed."""
        mech = GaussianMechanism(epsilon=1.0, delta=1e-5)

        samples = [mech.add_noise(0, 1.0) for _ in range(1000)]

        # Check that samples look roughly Gaussian
        # Skewness should be close to 0
        from scipy import stats
        _, p_value = stats.normaltest(samples)
        # Allow some deviation since we have finite samples
        assert p_value > 0.01 or abs(stats.skew(samples)) < 0.5


class TestExponentialMechanism:
    """Tests for Exponential mechanism."""

    def test_select(self):
        """Test selection from options."""
        mech = ExponentialMechanism(epsilon=1.0)
        options = ["a", "b", "c"]
        utilities = {"a": 0.0, "b": 0.5, "c": 1.0}

        # Should favor high-utility option
        selections = [
            mech.select(options, lambda x: utilities[x], 1.0)
            for _ in range(100)
        ]

        # c should be selected most often
        c_count = selections.count("c")
        assert c_count > 30  # Should be favored

    def test_select_top_k(self):
        """Test top-k selection."""
        mech = ExponentialMechanism(epsilon=1.0)
        options = list(range(10))

        selected = mech.select_top_k(
            options,
            utility_fn=lambda x: x,  # Higher index = higher utility
            sensitivity=1.0,
            k=3
        )

        assert len(selected) == 3


class TestReportNoisyMax:
    """Tests for Report Noisy Max mechanism."""

    def test_select(self):
        """Test argmax selection."""
        rnm = ReportNoisyMax(epsilon=1.0)
        scores = np.array([0.1, 0.3, 0.9, 0.5])

        # Run multiple times
        selections = [rnm.select(scores)[0] for _ in range(100)]

        # Index 2 (score 0.9) should be selected most often
        assert selections.count(2) > 30

    def test_select_top_k(self):
        """Test top-k selection."""
        rnm = ReportNoisyMax(epsilon=1.0)
        scores = np.array([0.1, 0.2, 0.9, 0.8, 0.5])

        results = rnm.select_top_k(scores, 1.0, k=3)

        assert len(results) == 3
        # Indices should be unique
        indices = [idx for idx, _ in results]
        assert len(set(indices)) == 3


class TestPrivacyAccountant:
    """Tests for privacy accountant."""

    def test_budget_tracking(self):
        """Test budget tracking."""
        accountant = PrivacyAccountant(total_budget=5.0)

        assert accountant.spend(1.0)
        assert accountant.spent == 1.0
        assert accountant.remaining == 4.0

    def test_budget_exhaustion(self):
        """Test that budget is enforced."""
        accountant = PrivacyAccountant(total_budget=2.0)

        assert accountant.spend(1.0)
        assert accountant.spend(1.0)
        assert not accountant.spend(0.5)  # Should fail

    def test_query_counting(self):
        """Test query counting."""
        accountant = PrivacyAccountant(total_budget=10.0)

        for _ in range(5):
            accountant.spend(1.0)

        assert accountant.num_queries == 5


class TestCreateDPMechanism:
    """Tests for factory function."""

    def test_create_laplace(self):
        """Test creating Laplace mechanism."""
        mech = create_dp_mechanism("laplace", epsilon=1.0)
        assert isinstance(mech, LaplaceMechanism)

    def test_create_gaussian(self):
        """Test creating Gaussian mechanism."""
        mech = create_dp_mechanism("gaussian", epsilon=1.0, delta=1e-5)
        assert isinstance(mech, GaussianMechanism)

    def test_invalid_type(self):
        """Test that invalid type raises error."""
        with pytest.raises(ValueError):
            create_dp_mechanism("invalid", epsilon=1.0)
