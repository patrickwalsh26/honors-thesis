"""
Matchmaker Exchange (MME) API adapter for privacy-preserving phenotype matching.

Implements GA4GH Matchmaker Exchange API v1.1 specification:
https://github.com/ga4gh/mme-apis

Features:
- MME request/response format conversion
- Privacy-preserving matching integration
- Controlled disclosure of shared phenotypes
"""

from typing import Dict, List, Tuple, Optional, Any, Set
from dataclasses import dataclass, field, asdict
from datetime import datetime
import json
import logging
import copy

logger = logging.getLogger(__name__)


@dataclass
class MMEContact:
    """Contact information for MME match follow-up."""
    name: str
    institution: str
    href: str  # Email or URL
    roles: List[str] = field(default_factory=lambda: ["clinician"])


@dataclass
class MMEFeature:
    """
    Phenotypic feature in MME format.

    Corresponds to HPO term annotation with observed status.
    """
    id: str  # HPO term ID
    label: Optional[str] = None
    observed: str = "yes"  # "yes", "no", or "unknown"

    def to_dict(self) -> Dict:
        result = {"id": self.id, "observed": self.observed}
        if self.label:
            result["label"] = self.label
        return result


@dataclass
class MMEGenomicFeature:
    """Genomic variant feature in MME format."""
    gene: Dict  # {"id": "HGNC:...", "label": "..."}
    variant: Optional[Dict] = None  # {"assembly": "GRCh38", ...}
    zygosity: Optional[int] = None
    type: Optional[Dict] = None

    def to_dict(self) -> Dict:
        result = {"gene": self.gene}
        if self.variant:
            result["variant"] = self.variant
        if self.zygosity is not None:
            result["zygosity"] = self.zygosity
        if self.type:
            result["type"] = self.type
        return result


@dataclass
class MMEDisorder:
    """Disease/disorder in MME format."""
    id: str  # OMIM, Orphanet, etc.
    label: Optional[str] = None

    def to_dict(self) -> Dict:
        result = {"id": self.id}
        if self.label:
            result["label"] = self.label
        return result


@dataclass
class MMEPatient:
    """Patient record in MME format."""
    id: str
    contact: MMEContact
    features: List[MMEFeature] = field(default_factory=list)
    genomicFeatures: List[MMEGenomicFeature] = field(default_factory=list)
    disorders: List[MMEDisorder] = field(default_factory=list)
    label: Optional[str] = None
    sex: Optional[str] = None
    ageOfOnset: Optional[Dict] = None

    def to_dict(self) -> Dict:
        result = {
            "id": self.id,
            "contact": asdict(self.contact),
            "features": [f.to_dict() for f in self.features],
            "genomicFeatures": [g.to_dict() for g in self.genomicFeatures],
            "disorders": [d.to_dict() for d in self.disorders]
        }
        if self.label:
            result["label"] = self.label
        if self.sex:
            result["sex"] = self.sex
        if self.ageOfOnset:
            result["ageOfOnset"] = self.ageOfOnset
        return result


@dataclass
class MMERequest:
    """MME match request."""
    patient: MMEPatient

    def to_dict(self) -> Dict:
        return {"patient": self.patient.to_dict()}

    @classmethod
    def from_dict(cls, data: Dict) -> "MMERequest":
        """Parse MME request from dictionary."""
        patient_data = data.get("patient", data)

        contact_data = patient_data.get("contact", {})
        contact = MMEContact(
            name=contact_data.get("name", "Unknown"),
            institution=contact_data.get("institution", "Unknown"),
            href=contact_data.get("href", ""),
            roles=contact_data.get("roles", ["clinician"])
        )

        features = [
            MMEFeature(
                id=f.get("id"),
                label=f.get("label"),
                observed=f.get("observed", "yes")
            )
            for f in patient_data.get("features", [])
        ]

        genomic_features = [
            MMEGenomicFeature(
                gene=g.get("gene", {}),
                variant=g.get("variant"),
                zygosity=g.get("zygosity"),
                type=g.get("type")
            )
            for g in patient_data.get("genomicFeatures", [])
        ]

        disorders = [
            MMEDisorder(
                id=d.get("id"),
                label=d.get("label")
            )
            for d in patient_data.get("disorders", [])
        ]

        patient = MMEPatient(
            id=patient_data.get("id", "unknown"),
            contact=contact,
            features=features,
            genomicFeatures=genomic_features,
            disorders=disorders,
            label=patient_data.get("label"),
            sex=patient_data.get("sex"),
            ageOfOnset=patient_data.get("ageOfOnset")
        )

        return cls(patient=patient)


@dataclass
class MMEMatchScore:
    """Score information for a match result."""
    patient: float  # 0.0 to 1.0


@dataclass
class MMEResult:
    """Single match result in MME format."""
    patient: MMEPatient
    score: MMEMatchScore

    def to_dict(self) -> Dict:
        return {
            "patient": self.patient.to_dict(),
            "score": {"patient": self.score.patient}
        }


@dataclass
class MMEResponse:
    """MME match response."""
    results: List[MMEResult]
    disclaimer: Optional[str] = None

    def to_dict(self) -> Dict:
        result = {
            "results": [r.to_dict() for r in self.results],
            "_meta": {
                "timestamp": datetime.now().isoformat(),
                "apiVersion": "1.1"
            }
        }
        if self.disclaimer:
            result["_disclaimer"] = self.disclaimer
        return result


class MMEAdapter:
    """
    Adapter for Matchmaker Exchange API.

    Integrates with privacy-preserving calculator to provide
    MME-compatible matching with privacy protections.
    """

    def __init__(
        self,
        calculator,
        database: List[Dict],
        contact_info: MMEContact,
        score_threshold: float = 0.0,
        max_results: int = 10,
        privacy_enabled: bool = True
    ):
        """
        Initialize MME adapter.

        Args:
            calculator: Privacy-preserving calculator
            database: List of phenopackets in database
            contact_info: Contact info for this node
            score_threshold: Minimum score to return
            max_results: Maximum results per query
            privacy_enabled: Enable privacy-preserving features
        """
        self.calculator = calculator
        self.database = database
        self.contact_info = contact_info
        self.score_threshold = score_threshold
        self.max_results = max_results
        self.privacy_enabled = privacy_enabled

        # Build phenopacket index for quick lookup
        self._id_to_idx = {
            pp.get("id"): idx
            for idx, pp in enumerate(database)
        }

    def phenopacket_to_mme_patient(
        self,
        phenopacket: Dict,
        contact: Optional[MMEContact] = None
    ) -> MMEPatient:
        """
        Convert phenopacket to MME patient format.

        Args:
            phenopacket: GA4GH Phenopacket
            contact: Contact info (uses default if None)

        Returns:
            MME patient object
        """
        features = []
        for f in phenopacket.get("phenotypicFeatures", []):
            features.append(MMEFeature(
                id=f["type"]["id"],
                label=f["type"].get("label"),
                observed="no" if f.get("excluded", False) else "yes"
            ))

        disorders = []
        for d in phenopacket.get("diseases", []):
            disorders.append(MMEDisorder(
                id=d["term"]["id"],
                label=d["term"].get("label")
            ))

        subject = phenopacket.get("subject", {})

        return MMEPatient(
            id=phenopacket.get("id", "unknown"),
            contact=contact or self.contact_info,
            features=features,
            disorders=disorders,
            sex=subject.get("sex"),
            label=subject.get("id")
        )

    def mme_patient_to_phenopacket(self, patient: MMEPatient) -> Dict:
        """
        Convert MME patient to phenopacket format.

        Args:
            patient: MME patient object

        Returns:
            GA4GH Phenopacket dictionary
        """
        phenopacket = {
            "id": patient.id,
            "phenotypicFeatures": [
                {
                    "type": {"id": f.id, "label": f.label or ""},
                    "excluded": f.observed == "no"
                }
                for f in patient.features
            ],
            "diseases": [
                {"term": {"id": d.id, "label": d.label or ""}}
                for d in patient.disorders
            ]
        }

        if patient.sex:
            phenopacket["subject"] = {"sex": patient.sex}

        return phenopacket

    def _get_shared_features(
        self,
        query: Dict,
        matched: Dict
    ) -> List[MMEFeature]:
        """
        Get overlapping features for privacy-aware disclosure.

        Only reveals features that both patients share.

        Args:
            query: Query phenopacket
            matched: Matched phenopacket

        Returns:
            List of shared MME features
        """
        query_terms = {
            f["type"]["id"]
            for f in query.get("phenotypicFeatures", [])
            if not f.get("excluded", False)
        }

        shared = []
        for f in matched.get("phenotypicFeatures", []):
            if not f.get("excluded", False):
                term_id = f["type"]["id"]
                if term_id in query_terms:
                    shared.append(MMEFeature(
                        id=term_id,
                        label=f["type"].get("label"),
                        observed="yes"
                    ))

        return shared

    def match(self, request: MMERequest) -> MMEResponse:
        """
        Process MME match request.

        Args:
            request: MME match request

        Returns:
            MME response with privacy-protected results
        """
        # Convert request to phenopacket
        query_pp = self.mme_patient_to_phenopacket(request.patient)

        # Find matches using privacy calculator
        matches = self.calculator.find_most_similar(
            query_pp,
            self.database,
            top_k=self.max_results
        )

        # Handle suppressed results (k-anonymity)
        if matches is None:
            return MMEResponse(
                results=[],
                disclaimer="Query suppressed for privacy protection"
            )

        # Build results
        results = []
        for idx, score in matches:
            if score < self.score_threshold:
                continue

            matched_pp = self.database[idx]

            # Create anonymized patient for response
            if self.privacy_enabled:
                # Only reveal shared features
                shared_features = self._get_shared_features(query_pp, matched_pp)
                anonymized_id = f"match_{idx}"

                result_patient = MMEPatient(
                    id=anonymized_id,
                    contact=self.contact_info,
                    features=shared_features,
                    disorders=[]  # Don't reveal diagnosis
                )
            else:
                # Full disclosure (not recommended)
                result_patient = self.phenopacket_to_mme_patient(matched_pp)

            results.append(MMEResult(
                patient=result_patient,
                score=MMEMatchScore(patient=round(score, 4))
            ))

        return MMEResponse(
            results=results,
            disclaimer="Results from privacy-preserving phenotype matching"
        )

    def match_json(self, request_json: str) -> str:
        """
        Process MME request from JSON string.

        Args:
            request_json: JSON-encoded MME request

        Returns:
            JSON-encoded MME response
        """
        request_data = json.loads(request_json)
        request = MMERequest.from_dict(request_data)
        response = self.match(request)
        return json.dumps(response.to_dict(), indent=2)


def convert_phenopacket_to_mme(
    phenopacket: Dict,
    contact: Dict
) -> Dict:
    """
    Standalone conversion of phenopacket to MME format.

    Args:
        phenopacket: GA4GH Phenopacket
        contact: Contact information dict

    Returns:
        MME patient dictionary
    """
    mme_contact = MMEContact(
        name=contact.get("name", "Unknown"),
        institution=contact.get("institution", "Unknown"),
        href=contact.get("href", ""),
        roles=contact.get("roles", ["clinician"])
    )

    features = []
    for f in phenopacket.get("phenotypicFeatures", []):
        features.append({
            "id": f["type"]["id"],
            "label": f["type"].get("label", ""),
            "observed": "no" if f.get("excluded", False) else "yes"
        })

    disorders = []
    for d in phenopacket.get("diseases", []):
        disorders.append({
            "id": d["term"]["id"],
            "label": d["term"].get("label", "")
        })

    return {
        "patient": {
            "id": phenopacket.get("id", "unknown"),
            "contact": asdict(mme_contact),
            "features": features,
            "disorders": disorders
        }
    }


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Example usage
    contact = MMEContact(
        name="Dr. Smith",
        institution="Stanford Medicine",
        href="mailto:drsmith@stanford.edu",
        roles=["clinician", "researcher"]
    )

    # Create sample phenopacket
    sample_pp = {
        "id": "patient_001",
        "phenotypicFeatures": [
            {"type": {"id": "HP:0001234", "label": "Feature A"}, "excluded": False},
            {"type": {"id": "HP:0002345", "label": "Feature B"}, "excluded": False}
        ],
        "diseases": [
            {"term": {"id": "OMIM:123456", "label": "Sample Disease"}}
        ],
        "subject": {"sex": "MALE"}
    }

    # Convert to MME format
    mme_dict = convert_phenopacket_to_mme(sample_pp, asdict(contact))
    print("Phenopacket to MME conversion:")
    print(json.dumps(mme_dict, indent=2))

    # Create MME request
    request = MMERequest.from_dict(mme_dict)
    print(f"\nParsed MME request for patient: {request.patient.id}")
    print(f"Number of features: {len(request.patient.features)}")
