"""
Evaluation module for privacy-preserving phenotype matching.

Provides:
- Standard retrieval metrics (Precision, Recall, F1, nDCG)
- Privacy-utility frontier analysis
- Privacy leakage audits (membership/attribute inference)
"""

from .metrics import (
    RetrievalMetrics,
    PhenopacketEvaluator,
    compute_metric_confidence_interval,
)

from .privacy_utility import (
    PrivacyUtilityFrontier,
)

from .leakage_audit import (
    AttackResults,
    MembershipInferenceAttack,
    AttributeInferenceAttack,
    LeakageAuditReport,
)

__all__ = [
    # Metrics
    "RetrievalMetrics",
    "PhenopacketEvaluator",
    "compute_metric_confidence_interval",
    # Privacy-Utility
    "PrivacyUtilityFrontier",
    # Leakage Audits
    "AttackResults",
    "MembershipInferenceAttack",
    "AttributeInferenceAttack",
    "LeakageAuditReport",
]
