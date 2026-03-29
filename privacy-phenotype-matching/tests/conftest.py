"""
Shared pytest fixtures for privacy-phenotype-matching tests.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture
def sample_phenopackets():
    """Generate sample phenopackets for testing."""
    phenopackets = []
    diseases = ["Marfan Syndrome", "Ehlers-Danlos", "Achondroplasia", "Progeria"]

    for i in range(40):  # 10 per disease
        disease_idx = i // 10
        disease = diseases[disease_idx]

        # Base terms for each disease
        base_terms = {
            "Marfan Syndrome": ["HP:0001166", "HP:0001519", "HP:0002616"],
            "Ehlers-Danlos": ["HP:0000974", "HP:0001382", "HP:0002616"],
            "Achondroplasia": ["HP:0003498", "HP:0002007", "HP:0000256"],
            "Progeria": ["HP:0001510", "HP:0002059", "HP:0001635"]
        }

        # Add some variation
        terms = base_terms[disease].copy()
        # Add 1-3 additional terms
        extra_terms = [f"HP:000100{j}" for j in range(i % 3 + 1)]
        terms.extend(extra_terms)

        phenopackets.append({
            "id": f"patient_{i}",
            "phenotypicFeatures": [
                {"type": {"id": term, "label": f"Term {term}"}, "excluded": False}
                for term in terms
            ],
            "diseases": [
                {"term": {"id": f"OMIM:{disease_idx}", "label": disease}}
            ],
            "subject": {"sex": "MALE" if i % 2 == 0 else "FEMALE"}
        })

    return phenopackets


@pytest.fixture
def small_phenopackets():
    """Small set of phenopackets for quick tests."""
    return [
        {
            "id": "patient_a",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0001234"}, "excluded": False},
                {"type": {"id": "HP:0002345"}, "excluded": False}
            ],
            "diseases": [{"term": {"id": "OMIM:123", "label": "Disease A"}}]
        },
        {
            "id": "patient_b",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0002345"}, "excluded": False},
                {"type": {"id": "HP:0003456"}, "excluded": False}
            ],
            "diseases": [{"term": {"id": "OMIM:123", "label": "Disease A"}}]
        },
        {
            "id": "patient_c",
            "phenotypicFeatures": [
                {"type": {"id": "HP:0004567"}, "excluded": False},
                {"type": {"id": "HP:0005678"}, "excluded": False}
            ],
            "diseases": [{"term": {"id": "OMIM:456", "label": "Disease B"}}]
        }
    ]


@pytest.fixture
def privacy_config():
    """Default privacy configuration for tests."""
    from src.privacy.k_anonymity import PrivacyConfig
    return PrivacyConfig(
        k=3,
        epsilon=1.0,
        delta=1e-5,
        min_prevalence=0.05
    )


@pytest.fixture
def ic_values(sample_phenopackets):
    """Compute IC values from sample phenopackets."""
    from src.similarity.hpo_similarity import compute_empirical_ic
    return compute_empirical_ic(sample_phenopackets)


@pytest.fixture
def base_calculator(ic_values):
    """Create base similarity calculator."""
    from src.similarity.hpo_similarity import (
        CosineSimilarity, PhenopacketSimilarityCalculator
    )
    metric = CosineSimilarity(ic_values)
    return PhenopacketSimilarityCalculator(metric)
