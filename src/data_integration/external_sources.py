"""
Data integration utilities for external rare disease data sources.

This module provides interfaces to connect to and retrieve data from
various rare disease databases and registries.

Supported data sources:
- OMIM (Online Mendelian Inheritance in Man)
- Orphanet (European rare disease database)
- ClinVar (Clinical variant database)
- HPO Annotations (phenotype.hpoa)
- DECIPHER (Database of Chromosomal Imbalance and Phenotype in Humans)

Future data sources (pending access approval):
- MIMIC-IV (de-identified ICU data)
- STARR-OMOP-deid (Stanford de-identified EHR)
"""

import json
import logging
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from dataclasses import dataclass, asdict
from datetime import datetime
import csv
import gzip

logger = logging.getLogger(__name__)


@dataclass
class DiseaseAnnotation:
    """A disease-phenotype annotation from HPO."""
    disease_id: str
    disease_name: str
    hpo_id: str
    hpo_name: str
    frequency: Optional[str] = None  # e.g., "HP:0040280" (Obligate)
    aspect: str = "P"  # P=Phenotype, I=Inheritance, C=Clinical course
    source: str = "HPOA"


@dataclass
class ORPHADisease:
    """Disease information from Orphanet."""
    orpha_code: str
    name: str
    prevalence: Optional[str] = None
    inheritance: List[str] = None
    age_of_onset: List[str] = None
    hpo_terms: List[str] = None


class HPOAnnotationsLoader:
    """
    Load and query HPO phenotype annotations (phenotype.hpoa).

    This file contains disease-phenotype associations curated by HPO.
    """

    def __init__(self, hpoa_path: str = "data/hpo_ontology/phenotype.hpoa"):
        """
        Initialize loader.

        Args:
            hpoa_path: Path to phenotype.hpoa file
        """
        self.hpoa_path = Path(hpoa_path)
        self.annotations: List[DiseaseAnnotation] = []
        self.disease_to_phenotypes: Dict[str, List[str]] = {}
        self.phenotype_to_diseases: Dict[str, List[str]] = {}
        self._loaded = False

    def download_hpoa(self, force: bool = False) -> Path:
        """Download the HPO annotations file."""
        url = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"

        self.hpoa_path.parent.mkdir(parents=True, exist_ok=True)

        if self.hpoa_path.exists() and not force:
            logger.info(f"HPOA file already exists at {self.hpoa_path}")
            return self.hpoa_path

        logger.info(f"Downloading HPO annotations from {url}")
        try:
            urllib.request.urlretrieve(url, self.hpoa_path)
            logger.info(f"Downloaded HPOA to {self.hpoa_path}")
        except Exception as e:
            logger.error(f"Failed to download HPOA: {e}")
            raise

        return self.hpoa_path

    def load(self) -> None:
        """Load annotations from file."""
        if self._loaded:
            return

        if not self.hpoa_path.exists():
            self.download_hpoa()

        logger.info(f"Loading HPO annotations from {self.hpoa_path}")

        self.annotations = []
        self.disease_to_phenotypes = {}
        self.phenotype_to_diseases = {}

        with open(self.hpoa_path, 'r', encoding='utf-8') as f:
            # Skip header lines starting with #
            reader = csv.reader(f, delimiter='\t')
            header = None

            for row in reader:
                if not row or row[0].startswith('#'):
                    # Check for column header
                    if row and row[0].startswith('#database_id'):
                        header = [c.lstrip('#') for c in row]
                    continue

                if header is None:
                    # Assume standard HPOA format
                    header = [
                        'database_id', 'disease_name', 'qualifier', 'hpo_id',
                        'reference', 'evidence', 'onset', 'frequency',
                        'sex', 'modifier', 'aspect', 'biocuration'
                    ]

                if len(row) < len(header):
                    continue

                data = dict(zip(header, row))

                # Only include phenotype annotations (aspect = P)
                aspect = data.get('aspect', 'P')
                if aspect != 'P':
                    continue

                disease_id = data.get('database_id', '')
                disease_name = data.get('disease_name', '')
                hpo_id = data.get('hpo_id', '')
                frequency = data.get('frequency', '')

                if not disease_id or not hpo_id:
                    continue

                annotation = DiseaseAnnotation(
                    disease_id=disease_id,
                    disease_name=disease_name,
                    hpo_id=hpo_id,
                    hpo_name="",  # Would need HPO lookup
                    frequency=frequency if frequency else None,
                    aspect=aspect
                )
                self.annotations.append(annotation)

                # Index by disease
                if disease_id not in self.disease_to_phenotypes:
                    self.disease_to_phenotypes[disease_id] = []
                if hpo_id not in self.disease_to_phenotypes[disease_id]:
                    self.disease_to_phenotypes[disease_id].append(hpo_id)

                # Index by phenotype
                if hpo_id not in self.phenotype_to_diseases:
                    self.phenotype_to_diseases[hpo_id] = []
                if disease_id not in self.phenotype_to_diseases[hpo_id]:
                    self.phenotype_to_diseases[hpo_id].append(disease_id)

        self._loaded = True
        logger.info(f"Loaded {len(self.annotations)} annotations for {len(self.disease_to_phenotypes)} diseases")

    def get_disease_phenotypes(self, disease_id: str) -> List[str]:
        """Get HPO terms associated with a disease."""
        if not self._loaded:
            self.load()
        return self.disease_to_phenotypes.get(disease_id, [])

    def get_phenotype_diseases(self, hpo_id: str) -> List[str]:
        """Get diseases associated with an HPO term."""
        if not self._loaded:
            self.load()
        return self.phenotype_to_diseases.get(hpo_id, [])

    def get_all_rare_diseases(self, max_prevalence_term: str = None) -> List[str]:
        """
        Get list of all rare diseases.

        Args:
            max_prevalence_term: HPO frequency term for max prevalence filter
        """
        if not self._loaded:
            self.load()
        return list(self.disease_to_phenotypes.keys())


class OMIMClient:
    """
    Client for OMIM (Online Mendelian Inheritance in Man) API.

    Note: OMIM requires API key for access. Get one at https://omim.org/api
    """

    BASE_URL = "https://api.omim.org/api"

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize OMIM client.

        Args:
            api_key: OMIM API key (required for API access)
        """
        self.api_key = api_key

    def search_by_phenotype(self, hpo_ids: List[str]) -> List[Dict]:
        """
        Search OMIM for diseases matching phenotypes.

        Args:
            hpo_ids: List of HPO term IDs

        Returns:
            List of matching disease entries
        """
        if not self.api_key:
            logger.warning("OMIM API key not configured. Set api_key to use this feature.")
            return []

        # This would make API calls to OMIM
        # Placeholder for actual implementation
        logger.info(f"Would search OMIM for {len(hpo_ids)} HPO terms")
        return []

    def get_disease_info(self, omim_id: str) -> Optional[Dict]:
        """
        Get detailed information for an OMIM entry.

        Args:
            omim_id: OMIM ID (e.g., "154700")

        Returns:
            Disease information dictionary
        """
        if not self.api_key:
            return None

        # Placeholder for API implementation
        return None


class OrphanetClient:
    """
    Client for Orphanet rare disease database.

    Orphanet provides European reference data on rare diseases.
    Data can be downloaded from: https://www.orphadata.com/
    """

    def __init__(self, data_dir: str = "data/orphanet"):
        """
        Initialize Orphanet client.

        Args:
            data_dir: Directory for Orphanet data files
        """
        self.data_dir = Path(data_dir)
        self.diseases: Dict[str, ORPHADisease] = {}
        self._loaded = False

    def download_data(self) -> None:
        """
        Download Orphanet data files.

        Downloads:
        - Disease nomenclature (en_product1.json)
        - Disease-HPO associations (en_product4.xml)
        """
        self.data_dir.mkdir(parents=True, exist_ok=True)

        # Orphanet provides data in various formats
        # This is a placeholder - actual download requires registration
        logger.info("Orphanet data download would require Orphadata access")
        logger.info("Register at: https://www.orphadata.com/")

    def load_disease_phenotypes(self, xml_path: str = None) -> None:
        """Load disease-phenotype associations from Orphanet XML."""
        if xml_path is None:
            xml_path = self.data_dir / "en_product4_HPO.xml"

        if not Path(xml_path).exists():
            logger.warning(f"Orphanet HPO file not found: {xml_path}")
            return

        # XML parsing would go here
        logger.info(f"Would load Orphanet data from {xml_path}")

    def get_disease_by_orpha(self, orpha_code: str) -> Optional[ORPHADisease]:
        """Get disease information by ORPHA code."""
        return self.diseases.get(orpha_code)

    def search_diseases(self, query: str) -> List[ORPHADisease]:
        """Search diseases by name."""
        if not self._loaded:
            return []

        query_lower = query.lower()
        return [
            d for d in self.diseases.values()
            if query_lower in d.name.lower()
        ]


class ClinVarClient:
    """
    Client for ClinVar variant database.

    Useful for finding variants associated with rare diseases.
    """

    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    def __init__(self, api_key: Optional[str] = None):
        """
        Initialize ClinVar client.

        Args:
            api_key: NCBI API key (optional, but recommended)
        """
        self.api_key = api_key

    def search_variants_by_condition(
        self,
        condition: str,
        significance: str = "pathogenic"
    ) -> List[Dict]:
        """
        Search ClinVar for variants associated with a condition.

        Args:
            condition: Disease name or OMIM ID
            significance: Clinical significance filter

        Returns:
            List of variant entries
        """
        # Build search query
        query = f"{condition}[DIS] AND {significance}[CLINSIG]"

        # E-utilities API call would go here
        logger.info(f"Would search ClinVar for: {query}")
        return []


class DECIPHERClient:
    """
    Client for DECIPHER database.

    DECIPHER contains data on chromosomal abnormalities and their phenotypes.
    Access: https://www.deciphergenomics.org/
    """

    def __init__(self):
        """Initialize DECIPHER client."""
        pass

    def search_by_phenotypes(self, hpo_ids: List[str]) -> List[Dict]:
        """Search DECIPHER for cases matching phenotypes."""
        # Would require DECIPHER API access
        logger.info(f"Would search DECIPHER for {len(hpo_ids)} HPO terms")
        return []


class DataSourceAggregator:
    """
    Aggregate data from multiple rare disease sources.

    This class provides a unified interface to query multiple databases
    and combine results.
    """

    def __init__(
        self,
        hpoa_path: str = "data/hpo_ontology/phenotype.hpoa",
        omim_api_key: Optional[str] = None
    ):
        """
        Initialize aggregator.

        Args:
            hpoa_path: Path to HPO annotations file
            omim_api_key: OMIM API key
        """
        self.hpoa = HPOAnnotationsLoader(hpoa_path)
        self.omim = OMIMClient(omim_api_key)
        self.orphanet = OrphanetClient()

    def get_disease_phenotypes_all_sources(
        self,
        disease_id: str
    ) -> Dict[str, List[str]]:
        """
        Get phenotypes for a disease from all available sources.

        Args:
            disease_id: Disease identifier (OMIM or ORPHA code)

        Returns:
            Dictionary mapping source to list of HPO terms
        """
        results = {}

        # HPO annotations
        hpoa_terms = self.hpoa.get_disease_phenotypes(disease_id)
        if hpoa_terms:
            results["HPOA"] = hpoa_terms

        # Additional sources would be added here
        return results

    def find_candidate_diseases(
        self,
        patient_hpo_terms: List[str],
        min_overlap: int = 2
    ) -> List[Dict]:
        """
        Find candidate diseases for a patient based on phenotypes.

        Args:
            patient_hpo_terms: List of patient HPO terms
            min_overlap: Minimum number of matching terms

        Returns:
            List of candidate diseases with overlap scores
        """
        self.hpoa.load()

        candidates = []
        patient_set = set(patient_hpo_terms)

        for disease_id, disease_terms in self.hpoa.disease_to_phenotypes.items():
            disease_set = set(disease_terms)
            overlap = len(patient_set & disease_set)

            if overlap >= min_overlap:
                candidates.append({
                    "disease_id": disease_id,
                    "overlap_count": overlap,
                    "patient_terms_matched": list(patient_set & disease_set),
                    "jaccard": overlap / len(patient_set | disease_set)
                })

        # Sort by overlap
        candidates.sort(key=lambda x: x["overlap_count"], reverse=True)
        return candidates


# Utility functions for data conversion

def convert_omim_to_phenopacket(
    omim_data: Dict,
    patient_id: str = None
) -> Dict:
    """
    Convert OMIM disease data to Phenopacket format.

    Args:
        omim_data: OMIM API response data
        patient_id: Patient ID for the phenopacket

    Returns:
        Phenopacket dictionary
    """
    if patient_id is None:
        patient_id = f"patient_{datetime.now().strftime('%Y%m%d%H%M%S')}"

    # Extract phenotypes from OMIM data
    # This is a placeholder - actual implementation depends on OMIM response format
    phenopacket = {
        "id": f"phenopacket_{patient_id}",
        "subject": {
            "id": patient_id,
            "sex": "UNKNOWN_SEX"
        },
        "phenotypicFeatures": [],
        "diseases": [],
        "metaData": {
            "created": datetime.now().isoformat(),
            "created_by": "OMIMConverter",
            "resources": [
                {
                    "id": "hp",
                    "name": "Human Phenotype Ontology",
                    "namespace_prefix": "HP"
                }
            ],
            "phenopacket_schema_version": "2.0"
        }
    }

    return phenopacket


def export_to_matchmaker_exchange(phenopacket: Dict) -> Dict:
    """
    Convert Phenopacket to Matchmaker Exchange (MME) format.

    MME is the GA4GH standard for patient matching across institutions.

    Args:
        phenopacket: GA4GH Phenopacket dictionary

    Returns:
        MME-compatible patient record
    """
    mme_record = {
        "id": phenopacket.get("id", "unknown"),
        "contact": {
            "name": "Privacy-Preserving Matcher",
            "institution": "Stanford University",
            "href": "mailto:research@stanford.edu"
        },
        "features": [],
        "genomicFeatures": [],
        "disorders": []
    }

    # Convert phenotypic features
    for feature in phenopacket.get("phenotypicFeatures", []):
        if not feature.get("excluded", False):
            mme_record["features"].append({
                "id": feature["type"]["id"],
                "label": feature["type"].get("label", ""),
                "observed": "yes"
            })
        else:
            mme_record["features"].append({
                "id": feature["type"]["id"],
                "label": feature["type"].get("label", ""),
                "observed": "no"
            })

    # Convert diseases
    for disease in phenopacket.get("diseases", []):
        mme_record["disorders"].append({
            "id": disease["term"]["id"],
            "label": disease["term"].get("label", "")
        })

    return mme_record


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Test HPOA loader
    hpoa = HPOAnnotationsLoader()
    hpoa.load()

    print(f"Loaded {len(hpoa.disease_to_phenotypes)} diseases")

    # Test disease lookup
    marfan_id = "OMIM:154700"
    terms = hpoa.get_disease_phenotypes(marfan_id)
    print(f"\nMarfan syndrome ({marfan_id}): {len(terms)} HPO terms")

    if terms:
        print(f"First 5 terms: {terms[:5]}")

    # Test aggregator
    aggregator = DataSourceAggregator()
    candidates = aggregator.find_candidate_diseases(
        ["HP:0001519", "HP:0001166", "HP:0002616"],  # Marfan-like phenotypes
        min_overlap=2
    )
    print(f"\nTop 5 candidate diseases:")
    for c in candidates[:5]:
        print(f"  {c['disease_id']}: {c['overlap_count']} matching terms (Jaccard: {c['jaccard']:.3f})")
