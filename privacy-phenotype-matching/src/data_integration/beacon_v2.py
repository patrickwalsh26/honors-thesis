"""
GA4GH Beacon v2 adapter for privacy-preserving phenotype queries.

Implements Beacon v2 phenotype query endpoint with:
- Differential privacy for count queries
- k-anonymity threshold enforcement
- Range-bucketed responses

Reference: https://beacon-project.io/
"""

from typing import Dict, List, Set, Tuple, Optional, Any
from dataclasses import dataclass, field, asdict
from datetime import datetime
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


@dataclass
class BeaconMeta:
    """Beacon response metadata."""
    beaconId: str
    apiVersion: str = "v2.0"
    returnedSchemas: List[Dict] = field(default_factory=list)

    def to_dict(self) -> Dict:
        return {
            "beaconId": self.beaconId,
            "apiVersion": self.apiVersion,
            "returnedSchemas": self.returnedSchemas
        }


@dataclass
class BeaconOrganization:
    """Beacon organization information."""
    id: str
    name: str
    description: Optional[str] = None
    address: Optional[str] = None
    contactUrl: Optional[str] = None

    def to_dict(self) -> Dict:
        result = {"id": self.id, "name": self.name}
        if self.description:
            result["description"] = self.description
        if self.address:
            result["address"] = self.address
        if self.contactUrl:
            result["contactUrl"] = self.contactUrl
        return result


@dataclass
class BeaconFilter:
    """
    Beacon query filter.

    Supports phenotype (HPO) and disease (OMIM/Orphanet) filters.
    """
    id: str  # Term ID (e.g., "HP:0001234" or "OMIM:123456")
    scope: str = "phenotype"  # "phenotype", "disease", "gene"
    operator: str = "="  # "=", "!=", ">"
    includeDescendantTerms: bool = True

    def to_dict(self) -> Dict:
        return {
            "id": self.id,
            "scope": self.scope,
            "operator": self.operator,
            "includeDescendantTerms": self.includeDescendantTerms
        }


@dataclass
class BeaconQuery:
    """Beacon v2 query structure."""
    filters: List[BeaconFilter]
    requestedGranularity: str = "count"  # "boolean", "count", "record"
    pagination: Optional[Dict] = None

    @classmethod
    def from_hpo_terms(
        cls,
        terms: List[str],
        granularity: str = "count"
    ) -> "BeaconQuery":
        """
        Create query from list of HPO terms.

        Args:
            terms: List of HPO term IDs
            granularity: Response granularity

        Returns:
            BeaconQuery instance
        """
        filters = [
            BeaconFilter(id=term, scope="phenotype")
            for term in terms
        ]
        return cls(filters=filters, requestedGranularity=granularity)

    def to_dict(self) -> Dict:
        result = {
            "filters": [f.to_dict() for f in self.filters],
            "requestedGranularity": self.requestedGranularity
        }
        if self.pagination:
            result["pagination"] = self.pagination
        return result


@dataclass
class BeaconResultset:
    """Beacon query resultset."""
    id: str
    setType: str
    exists: bool
    resultsCount: Optional[int] = None
    results: Optional[List[Dict]] = None

    def to_dict(self) -> Dict:
        result = {
            "id": self.id,
            "setType": self.setType,
            "exists": self.exists
        }
        if self.resultsCount is not None:
            result["resultsCount"] = self.resultsCount
        if self.results is not None:
            result["results"] = self.results
        return result


@dataclass
class BeaconResponse:
    """Complete Beacon v2 response."""
    meta: BeaconMeta
    responseSummary: Dict
    response: Dict

    def to_dict(self) -> Dict:
        return {
            "meta": self.meta.to_dict(),
            "responseSummary": self.responseSummary,
            "response": self.response
        }


class BeaconAdapter:
    """
    Beacon v2 adapter with privacy protections.

    Implements privacy-preserving count queries with:
    - Differential privacy noise
    - k-anonymity suppression
    - Range-bucketed responses
    """

    def __init__(
        self,
        database: List[Dict],
        beacon_id: str = "org.stanford.privacy-phenotype-matching",
        organization: Optional[BeaconOrganization] = None,
        epsilon: float = 1.0,
        k_anonymity: int = 5,
        hpo_manager=None
    ):
        """
        Initialize Beacon adapter.

        Args:
            database: List of phenopackets
            beacon_id: Unique beacon identifier
            organization: Organization information
            epsilon: DP epsilon for count queries
            k_anonymity: Minimum count for disclosure
            hpo_manager: HPOManager for term expansion
        """
        self.database = database
        self.beacon_id = beacon_id
        self.organization = organization or BeaconOrganization(
            id="org.stanford",
            name="Stanford University",
            description="Privacy-Preserving Phenotype Matching Beacon"
        )
        self.epsilon = epsilon
        self.k_anonymity = k_anonymity
        self.hpo_manager = hpo_manager

        # Initialize DP mechanism
        from ..privacy.differential_privacy import LaplaceMechanism
        self._dp_mechanism = LaplaceMechanism(epsilon=epsilon)

        # Build term index
        self._term_index = self._build_term_index()

        # Statistics
        self._total_queries = 0
        self._suppressed_queries = 0

    def _build_term_index(self) -> Dict[str, Set[int]]:
        """Build index mapping terms to patient indices."""
        index = defaultdict(set)

        for idx, pp in enumerate(self.database):
            for feature in pp.get("phenotypicFeatures", []):
                if not feature.get("excluded", False):
                    term_id = feature["type"]["id"]
                    index[term_id].add(idx)

                    # Also index ancestors if HPO manager available
                    if self.hpo_manager:
                        ancestors = self.hpo_manager.get_ancestors(
                            term_id, include_self=False
                        )
                        for anc in ancestors:
                            index[anc].add(idx)

        return index

    def _get_matching_patients(
        self,
        filters: List[BeaconFilter],
        logic: str = "AND"
    ) -> Set[int]:
        """
        Find patients matching all filters.

        Args:
            filters: List of BeaconFilter objects
            logic: Filter combination logic ("AND" or "OR")

        Returns:
            Set of matching patient indices
        """
        if not filters:
            return set(range(len(self.database)))

        matching_sets = []

        for f in filters:
            term_matches = set()

            # Direct term matches
            term_matches.update(self._term_index.get(f.id, set()))

            # Include descendants if requested
            if f.includeDescendantTerms and self.hpo_manager:
                descendants = self.hpo_manager.get_descendants(
                    f.id, include_self=False
                )
                for desc in descendants:
                    term_matches.update(self._term_index.get(desc, set()))

            # Handle operator
            if f.operator == "!=":
                # Negate: all patients NOT matching
                term_matches = set(range(len(self.database))) - term_matches

            matching_sets.append(term_matches)

        # Combine sets based on logic
        if logic == "AND":
            if not matching_sets:
                return set()
            result = matching_sets[0]
            for s in matching_sets[1:]:
                result = result.intersection(s)
            return result
        else:  # OR
            result = set()
            for s in matching_sets:
                result = result.union(s)
            return result

    def _apply_dp_noise(self, count: int) -> int:
        """Apply differential privacy noise to count."""
        noisy = self._dp_mechanism.privatize_count(count)
        return max(0, int(round(noisy)))

    def _apply_k_anonymity(
        self,
        count: int
    ) -> Tuple[bool, Optional[int]]:
        """
        Apply k-anonymity to count.

        Returns:
            Tuple of (exists, count_or_none)
            If count < k, returns (True, None) for boolean-only response
        """
        if count == 0:
            return False, 0
        elif count < self.k_anonymity:
            # Suppress exact count, confirm existence
            return True, None
        else:
            return True, count

    def count_query(
        self,
        query: BeaconQuery,
        add_noise: bool = True
    ) -> BeaconResponse:
        """
        Execute count query with privacy protections.

        Args:
            query: Beacon query
            add_noise: Whether to add DP noise

        Returns:
            Beacon response
        """
        self._total_queries += 1

        # Find matching patients
        matching = self._get_matching_patients(query.filters)
        true_count = len(matching)

        # Apply DP noise
        if add_noise:
            noisy_count = self._apply_dp_noise(true_count)
        else:
            noisy_count = true_count

        # Apply k-anonymity
        exists, disclosed_count = self._apply_k_anonymity(noisy_count)

        # Determine response granularity
        if disclosed_count is None:
            # k-anonymity suppressed - return boolean only
            self._suppressed_queries += 1
            return_granularity = "boolean"
            result_count = None
        elif query.requestedGranularity == "boolean":
            return_granularity = "boolean"
            result_count = None
        else:
            return_granularity = "count"
            result_count = disclosed_count

        # Build response
        meta = BeaconMeta(beaconId=self.beacon_id)

        response_summary = {
            "exists": exists,
            "numTotalResults": result_count
        }

        response = {
            "resultSets": [
                BeaconResultset(
                    id="phenotype_matches",
                    setType="individual",
                    exists=exists,
                    resultsCount=result_count
                ).to_dict()
            ]
        }

        return BeaconResponse(
            meta=meta,
            responseSummary=response_summary,
            response=response
        )

    def range_query(
        self,
        query: BeaconQuery,
        ranges: Optional[List[Tuple[int, int]]] = None
    ) -> Dict:
        """
        Execute range-bucketed count query.

        More privacy-friendly than exact counts.

        Args:
            query: Beacon query
            ranges: Count ranges [(0, 0), (1, 5), (6, 20), (21, inf)]

        Returns:
            Range bucket result
        """
        if ranges is None:
            ranges = [
                (0, 0),
                (1, 5),
                (6, 20),
                (21, 100),
                (101, float('inf'))
            ]

        # Get noisy count
        response = self.count_query(query, add_noise=True)
        count = response.responseSummary.get("numTotalResults")

        if count is None:
            # Suppressed
            return {
                "range": "1+",
                "note": "Exact count suppressed for privacy"
            }

        # Find bucket
        for low, high in ranges:
            if low <= count <= high:
                if high == float('inf'):
                    return {"range": f"{low}+"}
                elif low == high:
                    return {"range": str(low)}
                else:
                    return {"range": f"{low}-{high}"}

        return {"range": "unknown"}

    def phenotype_exists(self, hpo_term: str) -> bool:
        """
        Check if any patient has a phenotype (boolean query).

        Args:
            hpo_term: HPO term ID

        Returns:
            True if at least one patient has the phenotype
        """
        query = BeaconQuery.from_hpo_terms([hpo_term], granularity="boolean")
        response = self.count_query(query)
        return response.responseSummary["exists"]

    def get_info(self) -> Dict:
        """
        Get Beacon info response.

        Returns:
            Beacon information dictionary
        """
        return {
            "meta": {
                "beaconId": self.beacon_id,
                "apiVersion": "v2.0"
            },
            "response": {
                "id": self.beacon_id,
                "name": "Privacy-Preserving Phenotype Beacon",
                "organization": self.organization.to_dict(),
                "description": "Beacon with differential privacy and k-anonymity",
                "createDateTime": datetime.now().isoformat(),
                "datasets": [{
                    "id": "phenopackets",
                    "name": "Phenopacket Database",
                    "description": "GA4GH Phenopackets v2.0",
                    "assemblyId": "N/A",
                    "createDateTime": datetime.now().isoformat(),
                    "info": {
                        "numIndividuals": len(self.database),
                        "privacyProtections": [
                            f"Differential Privacy (epsilon={self.epsilon})",
                            f"k-Anonymity (k={self.k_anonymity})"
                        ]
                    }
                }]
            }
        }

    def get_statistics(self) -> Dict:
        """Get query statistics."""
        return {
            "total_queries": self._total_queries,
            "suppressed_queries": self._suppressed_queries,
            "suppression_rate": (
                self._suppressed_queries / self._total_queries
                if self._total_queries > 0 else 0
            ),
            "database_size": len(self.database),
            "unique_terms_indexed": len(self._term_index),
            "epsilon": self.epsilon,
            "k_anonymity": self.k_anonymity
        }


def create_beacon(
    database: List[Dict],
    beacon_id: str = "org.phenotype.beacon",
    epsilon: float = 1.0,
    k_anonymity: int = 5,
    hpo_manager=None
) -> BeaconAdapter:
    """
    Factory function to create Beacon adapter.

    Args:
        database: Phenopacket database
        beacon_id: Beacon identifier
        epsilon: DP epsilon
        k_anonymity: k-anonymity threshold
        hpo_manager: HPO manager for term expansion

    Returns:
        Configured BeaconAdapter
    """
    return BeaconAdapter(
        database=database,
        beacon_id=beacon_id,
        epsilon=epsilon,
        k_anonymity=k_anonymity,
        hpo_manager=hpo_manager
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Example usage
    sample_database = [
        {
            "id": f"patient_{i}",
            "phenotypicFeatures": [
                {"type": {"id": f"HP:000000{j}"}, "excluded": False}
                for j in range(1, 5 + i % 3)
            ]
        }
        for i in range(100)
    ]

    beacon = create_beacon(
        database=sample_database,
        epsilon=1.0,
        k_anonymity=5
    )

    # Test info endpoint
    print("Beacon Info:")
    info = beacon.get_info()
    print(f"  ID: {info['response']['id']}")
    print(f"  Individuals: {info['response']['datasets'][0]['info']['numIndividuals']}")

    # Test count query
    query = BeaconQuery.from_hpo_terms(["HP:0000001", "HP:0000002"])
    response = beacon.count_query(query)
    print(f"\nCount query result:")
    print(f"  Exists: {response.responseSummary['exists']}")
    print(f"  Count: {response.responseSummary['numTotalResults']}")

    # Test range query
    range_result = beacon.range_query(query)
    print(f"\nRange query result: {range_result}")

    # Statistics
    stats = beacon.get_statistics()
    print(f"\nStatistics:")
    for k, v in stats.items():
        print(f"  {k}: {v}")
