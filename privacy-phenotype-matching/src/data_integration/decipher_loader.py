"""
DECIPHER Database Loader for Privacy-Preserving Phenotype Matching.

This module provides tools to download, parse, and convert DECIPHER open-access
data into GA4GH Phenopacket format for evaluation.

DECIPHER (Database of Genomic Variation and Phenotype in Humans Using Ensembl
Resources) is a freely accessible database that shares anonymized phenotype-linked
variant data from rare disease patients.

Key Statistics (as of 2024):
- >44,000 patient records openly shared
- >57,000 variants
- >181,000 phenotype annotations (HPO terms)
- Data from >250 projects in ~40 countries
- Cited in >3,000 publications

Data Access:
- Open Access: https://www.deciphergenomics.org/
- API Documentation: https://www.deciphergenomics.org/api-docs
- Data Downloads: https://www.deciphergenomics.org/about/downloads

References:
    Foreman J, et al. (2022) DECIPHER: Supporting the interpretation and sharing
    of rare disease phenotype-linked variant data to advance diagnosis and research.
    Human Mutation. DOI: 10.1002/humu.24340
"""

import json
import logging
import urllib.request
import urllib.error
import csv
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Any, Iterator
from dataclasses import dataclass, asdict, field
from datetime import datetime
from collections import defaultdict
import time

logger = logging.getLogger(__name__)


@dataclass
class DECIPHERPatient:
    """
    A patient record from DECIPHER database.

    Attributes:
        patient_id: DECIPHER patient identifier
        sex: Patient sex (MALE, FEMALE, UNKNOWN)
        phenotypes: List of HPO term IDs
        phenotype_labels: Optional mapping of HPO IDs to labels
        variants: List of variant information (if available)
        inheritance: Mode of inheritance (if known)
        genes: Associated genes
        pathogenicity: Variant pathogenicity classifications
        contribution_id: Contributing project identifier
    """
    patient_id: str
    sex: str = "UNKNOWN_SEX"
    phenotypes: List[str] = field(default_factory=list)
    phenotype_labels: Dict[str, str] = field(default_factory=dict)
    variants: List[Dict] = field(default_factory=list)
    inheritance: Optional[str] = None
    genes: List[str] = field(default_factory=list)
    pathogenicity: List[str] = field(default_factory=list)
    contribution_id: Optional[str] = None


class DECIPHERLoader:
    """
    Loader for DECIPHER open-access patient data.

    This class provides methods to:
    1. Download DECIPHER open data dumps
    2. Parse patient records with HPO annotations
    3. Convert to GA4GH Phenopacket format
    4. Filter and preprocess for privacy experiments

    Example usage:
        >>> loader = DECIPHERLoader(data_dir="data/decipher")
        >>> loader.download_open_data()
        >>> patients = loader.load_patients()
        >>> phenopackets = loader.to_phenopackets(patients)
    """

    # DECIPHER data endpoints
    DECIPHER_DOWNLOADS_BASE = "https://www.deciphergenomics.org/about/downloads"
    DECIPHER_API_BASE = "https://www.deciphergenomics.org/api"

    # Alternative: Direct data file URLs (when available)
    OPEN_DATA_URLS = {
        "patients": "https://www.deciphergenomics.org/files/downloads/patients.txt.gz",
        "phenotypes": "https://www.deciphergenomics.org/files/downloads/phenotypes.txt.gz",
        "variants": "https://www.deciphergenomics.org/files/downloads/variants.txt.gz",
    }

    def __init__(
        self,
        data_dir: str = "data/decipher",
        cache_expiry_days: int = 30
    ):
        """
        Initialize DECIPHER loader.

        Args:
            data_dir: Directory for storing downloaded data
            cache_expiry_days: Days before re-downloading cached data
        """
        self.data_dir = Path(data_dir)
        self.cache_expiry_days = cache_expiry_days
        self.data_dir.mkdir(parents=True, exist_ok=True)

        # Internal data storage
        self._patients: Dict[str, DECIPHERPatient] = {}
        self._phenotype_index: Dict[str, List[str]] = {}  # HPO -> patient IDs
        self._loaded = False

    def download_open_data(self, force: bool = False) -> Dict[str, Path]:
        """
        Download DECIPHER open-access data files.

        Note: DECIPHER provides data through their website. For bulk downloads,
        you may need to use their Data Downloads page or request API access.

        Args:
            force: Force re-download even if files exist

        Returns:
            Dictionary mapping data type to file path
        """
        downloaded = {}

        for data_type, url in self.OPEN_DATA_URLS.items():
            output_path = self.data_dir / f"{data_type}.txt.gz"

            if output_path.exists() and not force:
                # Check cache age
                age_days = (datetime.now().timestamp() - output_path.stat().st_mtime) / 86400
                if age_days < self.cache_expiry_days:
                    logger.info(f"Using cached {data_type} data ({age_days:.1f} days old)")
                    downloaded[data_type] = output_path
                    continue

            logger.info(f"Downloading DECIPHER {data_type} data from {url}")
            try:
                urllib.request.urlretrieve(url, output_path)
                downloaded[data_type] = output_path
                logger.info(f"Downloaded {data_type} to {output_path}")
            except urllib.error.HTTPError as e:
                logger.warning(f"Could not download {data_type}: {e}")
                logger.info("DECIPHER data may require registration. See: https://www.deciphergenomics.org/about/downloads")
            except Exception as e:
                logger.error(f"Error downloading {data_type}: {e}")

        return downloaded

    def load_from_tsv(
        self,
        patients_file: str,
        phenotypes_file: Optional[str] = None,
        min_phenotypes: int = 1,
        max_patients: Optional[int] = None
    ) -> List[DECIPHERPatient]:
        """
        Load DECIPHER data from TSV/CSV files.

        This method handles the tab-separated format commonly used in
        DECIPHER data exports.

        Args:
            patients_file: Path to patients TSV file
            phenotypes_file: Path to phenotypes TSV file (optional)
            min_phenotypes: Minimum phenotypes required per patient
            max_patients: Maximum number of patients to load

        Returns:
            List of DECIPHERPatient objects
        """
        patients = {}
        patients_path = Path(patients_file)

        # Determine if gzipped
        open_func = gzip.open if str(patients_file).endswith('.gz') else open
        mode = 'rt' if str(patients_file).endswith('.gz') else 'r'

        logger.info(f"Loading DECIPHER patients from {patients_path}")

        try:
            with open_func(patients_path, mode, encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')

                for row in reader:
                    patient_id = row.get('patient_id') or row.get('id') or row.get('Patient ID')
                    if not patient_id:
                        continue

                    # Parse sex
                    sex_raw = row.get('sex', '').upper()
                    if 'MALE' in sex_raw and 'FEMALE' not in sex_raw:
                        sex = 'MALE'
                    elif 'FEMALE' in sex_raw:
                        sex = 'FEMALE'
                    else:
                        sex = 'UNKNOWN_SEX'

                    # Parse phenotypes (may be in various columns)
                    phenotypes = []
                    for col in ['phenotypes', 'hpo_terms', 'HPO', 'phenotype_ids']:
                        if col in row and row[col]:
                            # Handle various delimiters
                            terms = row[col].replace(';', ',').replace('|', ',').split(',')
                            phenotypes.extend([t.strip() for t in terms if t.strip().startswith('HP:')])

                    # Parse genes
                    genes = []
                    for col in ['genes', 'gene_symbol', 'Gene']:
                        if col in row and row[col]:
                            genes.extend([g.strip() for g in row[col].split(',') if g.strip()])

                    patient = DECIPHERPatient(
                        patient_id=patient_id,
                        sex=sex,
                        phenotypes=phenotypes,
                        genes=genes,
                        contribution_id=row.get('project_id') or row.get('contribution')
                    )
                    patients[patient_id] = patient

                    if max_patients and len(patients) >= max_patients:
                        break

        except FileNotFoundError:
            logger.error(f"File not found: {patients_path}")
            return []
        except Exception as e:
            logger.error(f"Error loading patients: {e}")
            return []

        # Load phenotypes from separate file if provided
        if phenotypes_file and Path(phenotypes_file).exists():
            self._load_phenotypes_file(phenotypes_file, patients)

        # Filter by minimum phenotypes
        filtered = [p for p in patients.values() if len(p.phenotypes) >= min_phenotypes]
        logger.info(f"Loaded {len(filtered)} patients with >= {min_phenotypes} phenotypes")

        return filtered

    def _load_phenotypes_file(
        self,
        phenotypes_file: str,
        patients: Dict[str, DECIPHERPatient]
    ) -> None:
        """Load phenotype annotations from separate file."""
        phenotypes_path = Path(phenotypes_file)
        open_func = gzip.open if str(phenotypes_file).endswith('.gz') else open
        mode = 'rt' if str(phenotypes_file).endswith('.gz') else 'r'

        logger.info(f"Loading phenotypes from {phenotypes_path}")

        try:
            with open_func(phenotypes_path, mode, encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')

                for row in reader:
                    patient_id = row.get('patient_id') or row.get('id')
                    hpo_id = row.get('hpo_id') or row.get('phenotype_id') or row.get('HPO')
                    hpo_label = row.get('hpo_label') or row.get('phenotype_name') or row.get('label')

                    if patient_id in patients and hpo_id:
                        if hpo_id not in patients[patient_id].phenotypes:
                            patients[patient_id].phenotypes.append(hpo_id)
                        if hpo_label:
                            patients[patient_id].phenotype_labels[hpo_id] = hpo_label

        except Exception as e:
            logger.warning(f"Error loading phenotypes file: {e}")

    def load_from_json(
        self,
        json_file: str,
        min_phenotypes: int = 1,
        max_patients: Optional[int] = None
    ) -> List[DECIPHERPatient]:
        """
        Load DECIPHER data from JSON format.

        Handles JSON API responses or JSON data exports.

        Args:
            json_file: Path to JSON file
            min_phenotypes: Minimum phenotypes required
            max_patients: Maximum patients to load

        Returns:
            List of DECIPHERPatient objects
        """
        json_path = Path(json_file)

        logger.info(f"Loading DECIPHER data from {json_path}")

        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
        except Exception as e:
            logger.error(f"Error loading JSON: {e}")
            return []

        # Handle different JSON structures
        if isinstance(data, list):
            records = data
        elif isinstance(data, dict):
            # JSON:API format or wrapped data
            records = data.get('data', data.get('patients', data.get('records', [data])))
            if isinstance(records, dict):
                records = [records]
        else:
            logger.error(f"Unexpected JSON structure: {type(data)}")
            return []

        patients = []
        for record in records:
            try:
                patient = self._parse_json_record(record)
                if len(patient.phenotypes) >= min_phenotypes:
                    patients.append(patient)

                if max_patients and len(patients) >= max_patients:
                    break
            except Exception as e:
                logger.debug(f"Error parsing record: {e}")
                continue

        logger.info(f"Loaded {len(patients)} patients from JSON")
        return patients

    def _parse_json_record(self, record: Dict) -> DECIPHERPatient:
        """Parse a single JSON record into DECIPHERPatient."""
        # Handle JSON:API format
        if 'type' in record and 'attributes' in record:
            patient_id = record.get('id', '')
            attrs = record.get('attributes', {})
        else:
            patient_id = record.get('patient_id') or record.get('id', '')
            attrs = record

        # Extract phenotypes
        phenotypes = []
        phenotype_labels = {}

        # Various possible field names
        pheno_data = attrs.get('phenotypes') or attrs.get('phenotypicFeatures') or attrs.get('hpo_terms', [])

        if isinstance(pheno_data, list):
            for p in pheno_data:
                if isinstance(p, str):
                    if p.startswith('HP:'):
                        phenotypes.append(p)
                elif isinstance(p, dict):
                    hpo_id = p.get('id') or p.get('hpo_id') or p.get('type', {}).get('id', '')
                    if hpo_id.startswith('HP:'):
                        phenotypes.append(hpo_id)
                        label = p.get('label') or p.get('hpo_label') or p.get('type', {}).get('label', '')
                        if label:
                            phenotype_labels[hpo_id] = label
        elif isinstance(pheno_data, str):
            phenotypes = [t.strip() for t in pheno_data.split(',') if t.strip().startswith('HP:')]

        # Extract sex
        sex_raw = str(attrs.get('sex', '')).upper()
        if 'MALE' in sex_raw and 'FEMALE' not in sex_raw:
            sex = 'MALE'
        elif 'FEMALE' in sex_raw:
            sex = 'FEMALE'
        else:
            sex = 'UNKNOWN_SEX'

        # Extract genes
        genes = []
        gene_data = attrs.get('genes') or attrs.get('gene_symbols', [])
        if isinstance(gene_data, list):
            genes = [g if isinstance(g, str) else g.get('symbol', '') for g in gene_data]
        elif isinstance(gene_data, str):
            genes = [g.strip() for g in gene_data.split(',') if g.strip()]

        return DECIPHERPatient(
            patient_id=str(patient_id),
            sex=sex,
            phenotypes=phenotypes,
            phenotype_labels=phenotype_labels,
            genes=genes,
            contribution_id=attrs.get('project_id') or attrs.get('contribution_id')
        )

    def to_phenopackets(
        self,
        patients: List[DECIPHERPatient],
        include_genes: bool = False
    ) -> List[Dict]:
        """
        Convert DECIPHER patients to GA4GH Phenopacket format.

        Args:
            patients: List of DECIPHER patients
            include_genes: Include gene information (may reduce privacy)

        Returns:
            List of Phenopacket dictionaries
        """
        phenopackets = []

        for patient in patients:
            # Build phenotypic features
            phenotypic_features = []
            for hpo_id in patient.phenotypes:
                feature = {
                    "type": {
                        "id": hpo_id,
                        "label": patient.phenotype_labels.get(hpo_id, f"HPO term {hpo_id}")
                    },
                    "excluded": False
                }
                phenotypic_features.append(feature)

            # Build phenopacket
            phenopacket = {
                "id": f"decipher_{patient.patient_id}",
                "subject": {
                    "id": patient.patient_id,
                    "sex": patient.sex
                },
                "phenotypicFeatures": phenotypic_features,
                "diseases": [],  # DECIPHER often doesn't include diagnosis labels
                "metaData": {
                    "created": datetime.now().isoformat(),
                    "created_by": "DECIPHERLoader",
                    "resources": [
                        {
                            "id": "hp",
                            "name": "Human Phenotype Ontology",
                            "url": "http://purl.obolibrary.org/obo/hp.owl",
                            "namespace_prefix": "HP"
                        }
                    ],
                    "phenopacket_schema_version": "2.0",
                    "external_references": [
                        {
                            "id": f"DECIPHER:{patient.patient_id}",
                            "description": "DECIPHER database patient record"
                        }
                    ]
                }
            }

            # Optionally include genes
            if include_genes and patient.genes:
                phenopacket["genes"] = [
                    {"id": gene, "symbol": gene}
                    for gene in patient.genes
                ]

            phenopackets.append(phenopacket)

        logger.info(f"Converted {len(phenopackets)} patients to Phenopacket format")
        return phenopackets

    def save_phenopackets(
        self,
        phenopackets: List[Dict],
        output_path: str,
        pretty: bool = True
    ) -> None:
        """Save phenopackets to JSON file."""
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, 'w', encoding='utf-8') as f:
            if pretty:
                json.dump(phenopackets, f, indent=2)
            else:
                json.dump(phenopackets, f)

        logger.info(f"Saved {len(phenopackets)} phenopackets to {output_path}")

    def get_phenotype_statistics(
        self,
        patients: List[DECIPHERPatient]
    ) -> Dict[str, Any]:
        """
        Compute statistics about phenotype distributions.

        Args:
            patients: List of DECIPHER patients

        Returns:
            Dictionary with statistics
        """
        all_phenotypes = []
        phenotype_counts = defaultdict(int)
        patients_per_phenotype = defaultdict(set)

        for patient in patients:
            all_phenotypes.extend(patient.phenotypes)
            for hpo in patient.phenotypes:
                phenotype_counts[hpo] += 1
                patients_per_phenotype[hpo].add(patient.patient_id)

        # Compute statistics
        phenotypes_per_patient = [len(p.phenotypes) for p in patients]

        stats = {
            "total_patients": len(patients),
            "total_phenotype_annotations": len(all_phenotypes),
            "unique_phenotypes": len(set(all_phenotypes)),
            "phenotypes_per_patient": {
                "mean": sum(phenotypes_per_patient) / len(phenotypes_per_patient) if phenotypes_per_patient else 0,
                "min": min(phenotypes_per_patient) if phenotypes_per_patient else 0,
                "max": max(phenotypes_per_patient) if phenotypes_per_patient else 0,
            },
            "most_common_phenotypes": sorted(
                phenotype_counts.items(),
                key=lambda x: x[1],
                reverse=True
            )[:20],
            "rare_phenotypes_count": sum(1 for c in phenotype_counts.values() if c == 1),
            "sex_distribution": {
                "MALE": sum(1 for p in patients if p.sex == "MALE"),
                "FEMALE": sum(1 for p in patients if p.sex == "FEMALE"),
                "UNKNOWN": sum(1 for p in patients if p.sex == "UNKNOWN_SEX"),
            }
        }

        return stats

    def filter_patients(
        self,
        patients: List[DECIPHERPatient],
        min_phenotypes: int = 3,
        max_phenotypes: Optional[int] = None,
        required_phenotypes: Optional[List[str]] = None,
        excluded_phenotypes: Optional[List[str]] = None,
        genes: Optional[List[str]] = None
    ) -> List[DECIPHERPatient]:
        """
        Filter patients based on various criteria.

        Args:
            patients: Input patient list
            min_phenotypes: Minimum number of phenotypes
            max_phenotypes: Maximum number of phenotypes
            required_phenotypes: HPO terms that must be present
            excluded_phenotypes: HPO terms that must not be present
            genes: Filter to patients with these genes

        Returns:
            Filtered patient list
        """
        filtered = []

        for patient in patients:
            # Check phenotype count
            n_pheno = len(patient.phenotypes)
            if n_pheno < min_phenotypes:
                continue
            if max_phenotypes and n_pheno > max_phenotypes:
                continue

            # Check required phenotypes
            if required_phenotypes:
                if not all(hp in patient.phenotypes for hp in required_phenotypes):
                    continue

            # Check excluded phenotypes
            if excluded_phenotypes:
                if any(hp in patient.phenotypes for hp in excluded_phenotypes):
                    continue

            # Check genes
            if genes:
                if not any(g in patient.genes for g in genes):
                    continue

            filtered.append(patient)

        logger.info(f"Filtered to {len(filtered)} patients (from {len(patients)})")
        return filtered

    def create_evaluation_cohort(
        self,
        patients: List[DECIPHERPatient],
        n_patients: int = 200,
        balance_by_phenotype_count: bool = True,
        random_seed: int = 42
    ) -> List[DECIPHERPatient]:
        """
        Create a balanced cohort for evaluation experiments.

        Args:
            patients: Full patient list
            n_patients: Target cohort size
            balance_by_phenotype_count: Balance across phenotype count ranges
            random_seed: Random seed for reproducibility

        Returns:
            Sampled patient cohort
        """
        import random
        random.seed(random_seed)

        if len(patients) <= n_patients:
            logger.warning(f"Requested {n_patients} but only {len(patients)} available")
            return patients.copy()

        if balance_by_phenotype_count:
            # Group by phenotype count ranges
            bins = [(1, 3), (4, 6), (7, 10), (11, 20), (21, float('inf'))]
            binned = defaultdict(list)

            for patient in patients:
                n = len(patient.phenotypes)
                for low, high in bins:
                    if low <= n <= high:
                        binned[(low, high)].append(patient)
                        break

            # Sample proportionally from each bin
            cohort = []
            per_bin = n_patients // len([b for b in binned.values() if b])

            for bin_key, bin_patients in binned.items():
                if bin_patients:
                    sample_size = min(per_bin, len(bin_patients))
                    cohort.extend(random.sample(bin_patients, sample_size))

            # Fill remainder
            remaining = [p for p in patients if p not in cohort]
            if len(cohort) < n_patients and remaining:
                cohort.extend(random.sample(remaining, min(n_patients - len(cohort), len(remaining))))
        else:
            cohort = random.sample(patients, n_patients)

        logger.info(f"Created evaluation cohort of {len(cohort)} patients")
        return cohort


class DECIPHERSimulator:
    """
    Simulate DECIPHER-like data when real data is not available.

    This is useful for testing and development when DECIPHER access
    is pending or restricted.
    """

    # Realistic phenotype distributions based on DECIPHER statistics
    PHENOTYPE_COUNT_DISTRIBUTION = [
        (1, 3, 0.15),   # 15% have 1-3 phenotypes
        (4, 6, 0.30),   # 30% have 4-6 phenotypes
        (7, 10, 0.30),  # 30% have 7-10 phenotypes
        (11, 15, 0.15), # 15% have 11-15 phenotypes
        (16, 25, 0.10), # 10% have 16+ phenotypes
    ]

    def __init__(
        self,
        hpo_terms_file: Optional[str] = None,
        seed: int = 42
    ):
        """
        Initialize simulator.

        Args:
            hpo_terms_file: Path to file with HPO terms (one per line)
            seed: Random seed
        """
        import random
        random.seed(seed)
        self.random = random

        # Load or generate HPO term pool
        self.hpo_pool = self._load_hpo_pool(hpo_terms_file)

    def _load_hpo_pool(self, hpo_file: Optional[str]) -> List[str]:
        """Load HPO terms for simulation."""
        if hpo_file and Path(hpo_file).exists():
            with open(hpo_file, 'r') as f:
                return [line.strip() for line in f if line.strip().startswith('HP:')]

        # Default pool of common phenotypic abnormality terms
        return [
            # Neurological
            "HP:0001249", "HP:0001250", "HP:0001251", "HP:0001252",
            "HP:0002360", "HP:0002376", "HP:0000708", "HP:0000739",
            # Craniofacial
            "HP:0000252", "HP:0000256", "HP:0000268", "HP:0000316",
            "HP:0000494", "HP:0000545", "HP:0000639", "HP:0000347",
            # Cardiovascular
            "HP:0001627", "HP:0001629", "HP:0001631", "HP:0001634",
            "HP:0001639", "HP:0001657", "HP:0002616", "HP:0001678",
            # Skeletal
            "HP:0000924", "HP:0001382", "HP:0001385", "HP:0002650",
            "HP:0002808", "HP:0002857", "HP:0003498", "HP:0002007",
            # Growth
            "HP:0001508", "HP:0001510", "HP:0001519", "HP:0000158",
            # Skin
            "HP:0000974", "HP:0000978", "HP:0001030", "HP:0001075",
            # Internal organs
            "HP:0001744", "HP:0002240", "HP:0001903", "HP:0001873",
            # Developmental
            "HP:0001263", "HP:0012758", "HP:0011968", "HP:0000717",
        ]

    def generate_patients(
        self,
        n_patients: int = 500,
        add_noise: bool = True
    ) -> List[DECIPHERPatient]:
        """
        Generate simulated DECIPHER-like patients.

        Args:
            n_patients: Number of patients to generate
            add_noise: Add realistic noise/variation

        Returns:
            List of simulated patients
        """
        patients = []

        for i in range(n_patients):
            # Determine phenotype count based on distribution
            n_phenotypes = self._sample_phenotype_count()

            # Sample phenotypes
            phenotypes = self.random.sample(
                self.hpo_pool,
                min(n_phenotypes, len(self.hpo_pool))
            )

            # Determine sex
            sex = self.random.choice(["MALE", "FEMALE", "UNKNOWN_SEX"])

            patient = DECIPHERPatient(
                patient_id=f"SIM_{i:06d}",
                sex=sex,
                phenotypes=phenotypes,
                contribution_id="SIMULATED"
            )
            patients.append(patient)

        logger.info(f"Generated {len(patients)} simulated patients")
        return patients

    def _sample_phenotype_count(self) -> int:
        """Sample phenotype count from distribution."""
        r = self.random.random()
        cumulative = 0

        for low, high, prob in self.PHENOTYPE_COUNT_DISTRIBUTION:
            cumulative += prob
            if r <= cumulative:
                return self.random.randint(low, high)

        return self.random.randint(5, 10)


def create_decipher_evaluation_dataset(
    output_dir: str = "data/decipher",
    n_patients: int = 500,
    use_simulated: bool = True
) -> str:
    """
    Create a DECIPHER evaluation dataset.

    This convenience function either loads real DECIPHER data or
    generates simulated data for testing.

    Args:
        output_dir: Output directory
        n_patients: Number of patients
        use_simulated: Use simulated data if real data unavailable

    Returns:
        Path to created phenopackets file
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    loader = DECIPHERLoader(data_dir=str(output_path))

    # Try to load real data first
    patients = []

    # Check for existing data files
    for ext in ['.json', '.tsv', '.txt.gz', '.csv']:
        for pattern in ['patients', 'decipher', 'data']:
            potential_file = output_path / f"{pattern}{ext}"
            if potential_file.exists():
                logger.info(f"Found existing data file: {potential_file}")
                if ext == '.json':
                    patients = loader.load_from_json(str(potential_file))
                else:
                    patients = loader.load_from_tsv(str(potential_file))
                break
        if patients:
            break

    # Fall back to simulated data
    if not patients and use_simulated:
        logger.info("No real data found, generating simulated DECIPHER-like data")
        simulator = DECIPHERSimulator()
        patients = simulator.generate_patients(n_patients)

    if not patients:
        raise RuntimeError("No data available and simulation disabled")

    # Create evaluation cohort
    if len(patients) > n_patients:
        patients = loader.create_evaluation_cohort(patients, n_patients)

    # Convert to phenopackets
    phenopackets = loader.to_phenopackets(patients)

    # Save
    phenopackets_path = output_path / f"decipher_cohort_{len(phenopackets)}.json"
    loader.save_phenopackets(phenopackets, str(phenopackets_path))

    # Print statistics
    stats = loader.get_phenotype_statistics(patients)
    logger.info(f"Dataset statistics:")
    logger.info(f"  Total patients: {stats['total_patients']}")
    logger.info(f"  Unique phenotypes: {stats['unique_phenotypes']}")
    logger.info(f"  Mean phenotypes/patient: {stats['phenotypes_per_patient']['mean']:.1f}")

    return str(phenopackets_path)


if __name__ == "__main__":
    # Example usage and testing
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    print("\n" + "=" * 70)
    print("  DECIPHER Data Loader - Integration Test")
    print("=" * 70)

    # Test simulated data generation
    print("\n1. Testing simulated data generation...")
    simulator = DECIPHERSimulator(seed=42)
    patients = simulator.generate_patients(100)

    print(f"   Generated {len(patients)} simulated patients")
    print(f"   Sample patient: {patients[0].patient_id}")
    print(f"   Phenotypes: {patients[0].phenotypes[:5]}...")

    # Test conversion to phenopackets
    print("\n2. Testing phenopacket conversion...")
    loader = DECIPHERLoader(data_dir="data/decipher")
    phenopackets = loader.to_phenopackets(patients[:10])

    print(f"   Converted {len(phenopackets)} phenopackets")
    print(f"   Sample ID: {phenopackets[0]['id']}")
    print(f"   Features: {len(phenopackets[0]['phenotypicFeatures'])}")

    # Test statistics
    print("\n3. Computing statistics...")
    stats = loader.get_phenotype_statistics(patients)

    print(f"   Total patients: {stats['total_patients']}")
    print(f"   Unique phenotypes: {stats['unique_phenotypes']}")
    print(f"   Mean phenotypes/patient: {stats['phenotypes_per_patient']['mean']:.1f}")
    print(f"   Rare phenotypes (n=1): {stats['rare_phenotypes_count']}")

    # Test creating evaluation dataset
    print("\n4. Creating evaluation dataset...")
    output_path = create_decipher_evaluation_dataset(
        output_dir="data/decipher",
        n_patients=200,
        use_simulated=True
    )
    print(f"   Saved to: {output_path}")

    print("\n" + "=" * 70)
    print("  Integration test complete!")
    print("=" * 70 + "\n")
