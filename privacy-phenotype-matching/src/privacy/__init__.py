"""
Privacy-preserving phenotype matching module.

This module implements multiple privacy mechanisms:
- Differential Privacy (DP) for noisy outputs
- Private Set Intersection (PSI) for secure overlap computation
- k-Anonymity for result suppression
- Rare term filtering for re-identification prevention
"""

from .differential_privacy import (
    DPMechanism,
    LaplaceMechanism,
    GaussianMechanism,
    ExponentialMechanism,
    ReportNoisyMax,
    PrivacyAccountant,
    PrivateSimilarityCalculator,
    create_dp_mechanism,
)

from .k_anonymity import (
    PrivacyConfig,
    RareTermFilter,
    KAnonymityGuard,
    EquivalenceClassAnalyzer,
    load_privacy_config,
)

from .psi import (
    PSIConfig,
    DiffieHellmanPSI,
    OPRFBasedPSI,
    HybridPSI,
    PSIPhenopacketMatcher,
    create_psi_matcher,
)

from .privacy_calculator import (
    PrivacyPreservingCalculator,
    TwoStepRevealLadder,
    create_privacy_calculator,
)

__all__ = [
    # Differential Privacy
    "DPMechanism",
    "LaplaceMechanism",
    "GaussianMechanism",
    "ExponentialMechanism",
    "ReportNoisyMax",
    "PrivacyAccountant",
    "PrivateSimilarityCalculator",
    "create_dp_mechanism",
    # k-Anonymity
    "PrivacyConfig",
    "RareTermFilter",
    "KAnonymityGuard",
    "EquivalenceClassAnalyzer",
    "load_privacy_config",
    # PSI
    "PSIConfig",
    "DiffieHellmanPSI",
    "OPRFBasedPSI",
    "HybridPSI",
    "PSIPhenopacketMatcher",
    "create_psi_matcher",
    # Unified Calculator
    "PrivacyPreservingCalculator",
    "TwoStepRevealLadder",
    "create_privacy_calculator",
]
