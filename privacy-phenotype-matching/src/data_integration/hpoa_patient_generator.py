"""
Generate synthetic patient cohorts from real HPO disease-phenotype annotations.

This module creates realistic patient datasets based on the HPO phenotype.hpoa file,
which contains curated disease-phenotype associations from OMIM, Orphanet, and DECIPHER.

The generated patients have:
- Phenotypes drawn from real disease profiles
- Realistic noise (missing phenotypes, phenotypic expansion)
- Ground truth labels for evaluation (patients with same disease are relevant)

This enables evaluation of phenotype matching algorithms with realistic data
while avoiding the need for real patient data access.
"""

import json
import logging
import random
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass, field, asdict
from collections import defaultdict
from datetime import datetime
import csv

logger = logging.getLogger(__name__)


@dataclass
class DiseaseProfile:
    """A disease with its associated phenotypes from HPOA."""
    disease_id: str
    disease_name: str
    phenotypes: List[str] = field(default_factory=list)
    source: str = "HPOA"  # OMIM, ORPHA, DECIPHER

    @property
    def phenotype_count(self) -> int:
        return len(self.phenotypes)


@dataclass
class SyntheticPatient:
    """A synthetic patient generated from a disease profile."""
    patient_id: str
    underlying_disease: str
    disease_name: str
    phenotypes: List[str] = field(default_factory=list)
    sex: str = "UNKNOWN_SEX"
    noise_phenotypes: List[str] = field(default_factory=list)  # Added noise
    missing_phenotypes: List[str] = field(default_factory=list)  # Removed from disease

    def to_phenopacket(self) -> Dict:
        """Convert to GA4GH Phenopacket format."""
        return {
            "id": f"phenopacket_{self.patient_id}",
            "subject": {
                "id": self.patient_id,
                "sex": self.sex
            },
            "phenotypicFeatures": [
                {
                    "type": {
                        "id": hpo_id,
                        "label": ""
                    }
                }
                for hpo_id in self.phenotypes
            ],
            "diseases": [
                {
                    "term": {
                        "id": self.underlying_disease,
                        "label": self.disease_name
                    }
                }
            ],
            "metaData": {
                "created": datetime.now().isoformat(),
                "created_by": "HPOAPatientGenerator",
                "phenopacket_schema_version": "2.0",
                "resources": [
                    {
                        "id": "hp",
                        "name": "Human Phenotype Ontology",
                        "namespace_prefix": "HP",
                        "url": "http://purl.obolibrary.org/obo/hp.owl"
                    }
                ]
            }
        }


class HPOAPatientGenerator:
    """
    Generate synthetic patients from real HPO disease-phenotype annotations.

    This generator creates realistic patient cohorts by:
    1. Loading disease profiles from phenotype.hpoa
    2. Sampling diseases based on various criteria
    3. Generating patients with realistic phenotype noise
    4. Providing ground truth for evaluation
    """

    def __init__(
        self,
        hpoa_path: str = "data/hpo_annotations/phenotype.hpoa",
        seed: int = None
    ):
        """
        Initialize generator.

        Args:
            hpoa_path: Path to phenotype.hpoa file
            seed: Random seed for reproducibility
        """
        self.hpoa_path = Path(hpoa_path)
        self.disease_profiles: Dict[str, DiseaseProfile] = {}
        self.all_phenotypes: Set[str] = set()
        self.phenotype_frequencies: Dict[str, int] = defaultdict(int)
        self._loaded = False

        if seed is not None:
            random.seed(seed)

    def load_annotations(self) -> None:
        """Load disease-phenotype annotations from HPOA file."""
        if self._loaded:
            return

        if not self.hpoa_path.exists():
            raise FileNotFoundError(
                f"HPOA file not found: {self.hpoa_path}\n"
                f"Download from: http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
            )

        logger.info(f"Loading HPOA annotations from {self.hpoa_path}")

        disease_phenotypes = defaultdict(lambda: {"name": "", "phenotypes": set()})

        with open(self.hpoa_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')

            for row in reader:
                if not row or row[0].startswith('#'):
                    continue

                if len(row) < 11:
                    continue

                # Parse HPOA columns
                disease_id = row[0]      # database_id (e.g., OMIM:154700)
                disease_name = row[1]    # disease_name
                qualifier = row[2]       # qualifier (NOT if excluded)
                hpo_id = row[3]          # hpo_id
                aspect = row[10] if len(row) > 10 else "P"  # aspect

                # Skip excluded phenotypes (qualifier = NOT)
                if qualifier.upper() == "NOT":
                    continue

                # Only include phenotype annotations (aspect = P)
                if aspect != "P":
                    continue

                if not disease_id or not hpo_id:
                    continue

                # Update disease profile
                disease_phenotypes[disease_id]["name"] = disease_name
                disease_phenotypes[disease_id]["phenotypes"].add(hpo_id)

                # Track phenotype frequencies
                self.all_phenotypes.add(hpo_id)
                self.phenotype_frequencies[hpo_id] += 1

        # Convert to DiseaseProfile objects
        for disease_id, data in disease_phenotypes.items():
            # Determine source from ID prefix
            if disease_id.startswith("OMIM:"):
                source = "OMIM"
            elif disease_id.startswith("ORPHA:"):
                source = "Orphanet"
            elif disease_id.startswith("DECIPHER:"):
                source = "DECIPHER"
            else:
                source = "Other"

            self.disease_profiles[disease_id] = DiseaseProfile(
                disease_id=disease_id,
                disease_name=data["name"],
                phenotypes=list(data["phenotypes"]),
                source=source
            )

        self._loaded = True
        logger.info(
            f"Loaded {len(self.disease_profiles)} diseases with "
            f"{len(self.all_phenotypes)} unique phenotypes"
        )

    def get_statistics(self) -> Dict:
        """Get dataset statistics."""
        if not self._loaded:
            self.load_annotations()

        phenotypes_per_disease = [p.phenotype_count for p in self.disease_profiles.values()]

        # Count by source
        sources = defaultdict(int)
        for profile in self.disease_profiles.values():
            sources[profile.source] += 1

        return {
            "total_diseases": len(self.disease_profiles),
            "total_unique_phenotypes": len(self.all_phenotypes),
            "phenotypes_per_disease": {
                "mean": sum(phenotypes_per_disease) / len(phenotypes_per_disease),
                "min": min(phenotypes_per_disease),
                "max": max(phenotypes_per_disease),
                "median": sorted(phenotypes_per_disease)[len(phenotypes_per_disease) // 2]
            },
            "diseases_by_source": dict(sources),
            "most_common_phenotypes": sorted(
                self.phenotype_frequencies.items(),
                key=lambda x: x[1],
                reverse=True
            )[:10]
        }

    def filter_diseases(
        self,
        min_phenotypes: int = 3,
        max_phenotypes: int = 50,
        sources: List[str] = None
    ) -> List[DiseaseProfile]:
        """
        Filter diseases by criteria.

        Args:
            min_phenotypes: Minimum phenotypes required
            max_phenotypes: Maximum phenotypes allowed
            sources: List of sources to include (OMIM, Orphanet, DECIPHER)

        Returns:
            List of matching disease profiles
        """
        if not self._loaded:
            self.load_annotations()

        filtered = []
        for profile in self.disease_profiles.values():
            if profile.phenotype_count < min_phenotypes:
                continue
            if profile.phenotype_count > max_phenotypes:
                continue
            if sources and profile.source not in sources:
                continue
            filtered.append(profile)

        return filtered

    def generate_patient(
        self,
        disease: DiseaseProfile,
        patient_id: str,
        phenotype_recall: float = 0.7,
        noise_rate: float = 0.1,
        sex: str = None
    ) -> SyntheticPatient:
        """
        Generate a synthetic patient from a disease profile.

        Args:
            disease: Disease profile to base patient on
            patient_id: Unique patient identifier
            phenotype_recall: Fraction of disease phenotypes to include (0-1)
            noise_rate: Fraction of phenotypes to add as noise (0-1)
            sex: Patient sex (random if None)

        Returns:
            Synthetic patient with realistic phenotype profile
        """
        if sex is None:
            sex = random.choice(["MALE", "FEMALE", "UNKNOWN_SEX"])

        # Sample phenotypes from disease (simulating incomplete phenotyping)
        n_to_keep = max(1, int(len(disease.phenotypes) * phenotype_recall))
        kept_phenotypes = random.sample(
            disease.phenotypes,
            min(n_to_keep, len(disease.phenotypes))
        )
        missing = [p for p in disease.phenotypes if p not in kept_phenotypes]

        # Add noise phenotypes from other diseases
        n_noise = int(len(kept_phenotypes) * noise_rate)
        noise_candidates = list(self.all_phenotypes - set(disease.phenotypes))
        noise_phenotypes = random.sample(
            noise_candidates,
            min(n_noise, len(noise_candidates))
        ) if noise_candidates and n_noise > 0 else []

        # Combine phenotypes
        all_phenotypes = kept_phenotypes + noise_phenotypes
        random.shuffle(all_phenotypes)

        return SyntheticPatient(
            patient_id=patient_id,
            underlying_disease=disease.disease_id,
            disease_name=disease.disease_name,
            phenotypes=all_phenotypes,
            sex=sex,
            noise_phenotypes=noise_phenotypes,
            missing_phenotypes=missing
        )

    def generate_cohort(
        self,
        n_patients: int = 300,
        n_diseases: int = None,
        patients_per_disease: int = None,
        min_phenotypes: int = 5,
        max_phenotypes: int = 30,
        sources: List[str] = None,
        phenotype_recall: float = 0.7,
        noise_rate: float = 0.1,
        balance_diseases: bool = True
    ) -> Tuple[List[SyntheticPatient], Dict]:
        """
        Generate a cohort of synthetic patients.

        Args:
            n_patients: Total number of patients to generate
            n_diseases: Number of diseases to sample (default: auto)
            patients_per_disease: Patients per disease (overrides n_patients)
            min_phenotypes: Minimum phenotypes per disease
            max_phenotypes: Maximum phenotypes per disease
            sources: Disease sources to include
            phenotype_recall: Fraction of disease phenotypes to keep
            noise_rate: Fraction of noise phenotypes to add
            balance_diseases: Whether to balance patients across diseases

        Returns:
            Tuple of (patients list, metadata dict)
        """
        if not self._loaded:
            self.load_annotations()

        # Filter candidate diseases
        candidates = self.filter_diseases(
            min_phenotypes=min_phenotypes,
            max_phenotypes=max_phenotypes,
            sources=sources
        )

        if not candidates:
            raise ValueError("No diseases match the specified criteria")

        logger.info(f"Found {len(candidates)} candidate diseases")

        # Determine sampling strategy
        if patients_per_disease is not None:
            # Fixed number per disease
            if n_diseases is None:
                n_diseases = min(n_patients // patients_per_disease, len(candidates))
            selected_diseases = random.sample(candidates, min(n_diseases, len(candidates)))
            n_patients = n_diseases * patients_per_disease
        else:
            # Sample diseases proportionally
            if n_diseases is None:
                n_diseases = min(n_patients // 3, len(candidates))  # ~3 patients per disease
            selected_diseases = random.sample(candidates, min(n_diseases, len(candidates)))

            if balance_diseases:
                patients_per_disease = n_patients // n_diseases
            else:
                # Variable patients per disease
                patients_per_disease = None

        # Generate patients
        patients = []
        disease_counts = defaultdict(int)

        if patients_per_disease:
            # Balanced assignment
            for disease in selected_diseases:
                for i in range(patients_per_disease):
                    patient_id = f"patient_{len(patients):04d}"
                    patient = self.generate_patient(
                        disease=disease,
                        patient_id=patient_id,
                        phenotype_recall=phenotype_recall,
                        noise_rate=noise_rate
                    )
                    patients.append(patient)
                    disease_counts[disease.disease_id] += 1
        else:
            # Random assignment
            for i in range(n_patients):
                disease = random.choice(selected_diseases)
                patient_id = f"patient_{i:04d}"
                patient = self.generate_patient(
                    disease=disease,
                    patient_id=patient_id,
                    phenotype_recall=phenotype_recall,
                    noise_rate=noise_rate
                )
                patients.append(patient)
                disease_counts[disease.disease_id] += 1

        # Shuffle patients
        random.shuffle(patients)

        # Build metadata
        phenotypes_per_patient = [len(p.phenotypes) for p in patients]
        metadata = {
            "n_patients": len(patients),
            "n_diseases": len(disease_counts),
            "generation_params": {
                "min_phenotypes": min_phenotypes,
                "max_phenotypes": max_phenotypes,
                "phenotype_recall": phenotype_recall,
                "noise_rate": noise_rate,
                "sources": sources
            },
            "phenotypes_per_patient": {
                "mean": sum(phenotypes_per_patient) / len(phenotypes_per_patient),
                "min": min(phenotypes_per_patient),
                "max": max(phenotypes_per_patient)
            },
            "disease_distribution": dict(disease_counts),
            "timestamp": datetime.now().isoformat()
        }

        logger.info(
            f"Generated {len(patients)} patients from {len(disease_counts)} diseases"
        )

        return patients, metadata

    def create_ground_truth(
        self,
        patients: List[SyntheticPatient],
        same_disease_relevant: bool = True,
        min_phenotype_overlap: int = None
    ) -> Dict[str, List[str]]:
        """
        Create ground truth relevance labels for evaluation.

        Args:
            patients: List of patients
            same_disease_relevant: Consider same-disease patients as relevant
            min_phenotype_overlap: Alternative: minimum phenotype overlap for relevance

        Returns:
            Dict mapping patient_id to list of relevant patient_ids
        """
        ground_truth = {}

        if same_disease_relevant:
            # Group patients by disease
            disease_patients = defaultdict(list)
            for patient in patients:
                disease_patients[patient.underlying_disease].append(patient.patient_id)

            # Build relevance lists
            for patient in patients:
                relevant = [
                    pid for pid in disease_patients[patient.underlying_disease]
                    if pid != patient.patient_id
                ]
                ground_truth[patient.patient_id] = relevant

        elif min_phenotype_overlap is not None:
            # Use phenotype overlap
            patient_phenotypes = {
                p.patient_id: set(p.phenotypes) for p in patients
            }

            for patient in patients:
                relevant = []
                query_pheno = patient_phenotypes[patient.patient_id]

                for other in patients:
                    if other.patient_id == patient.patient_id:
                        continue
                    other_pheno = patient_phenotypes[other.patient_id]
                    overlap = len(query_pheno & other_pheno)

                    if overlap >= min_phenotype_overlap:
                        relevant.append(other.patient_id)

                ground_truth[patient.patient_id] = relevant

        return ground_truth

    def export_cohort(
        self,
        patients: List[SyntheticPatient],
        output_path: str,
        format: str = "phenopackets"
    ) -> Path:
        """
        Export cohort to file.

        Args:
            patients: List of patients
            output_path: Output file path
            format: Output format (phenopackets, json, csv)

        Returns:
            Path to output file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format == "phenopackets":
            data = [p.to_phenopacket() for p in patients]
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2)

        elif format == "json":
            data = [asdict(p) for p in patients]
            with open(output_path, 'w') as f:
                json.dump(data, f, indent=2)

        elif format == "csv":
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    "patient_id", "disease_id", "disease_name",
                    "sex", "phenotypes", "n_phenotypes"
                ])
                for p in patients:
                    writer.writerow([
                        p.patient_id,
                        p.underlying_disease,
                        p.disease_name,
                        p.sex,
                        "|".join(p.phenotypes),
                        len(p.phenotypes)
                    ])
        else:
            raise ValueError(f"Unknown format: {format}")

        logger.info(f"Exported {len(patients)} patients to {output_path}")
        return output_path


def create_evaluation_dataset(
    output_dir: str = "data/hpoa_evaluation",
    n_patients: int = 500,
    n_diseases: int = 100,
    seed: int = 42
) -> Tuple[List[Dict], Dict, Dict]:
    """
    Create a complete evaluation dataset from HPOA.

    Args:
        output_dir: Output directory
        n_patients: Number of patients to generate
        n_diseases: Number of diseases to sample
        seed: Random seed

    Returns:
        Tuple of (phenopackets, ground_truth, metadata)
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize generator
    generator = HPOAPatientGenerator(
        hpoa_path="data/hpo_annotations/phenotype.hpoa",
        seed=seed
    )

    # Generate cohort
    patients, metadata = generator.generate_cohort(
        n_patients=n_patients,
        n_diseases=n_diseases,
        min_phenotypes=5,
        max_phenotypes=30,
        phenotype_recall=0.7,
        noise_rate=0.15,
        balance_diseases=True
    )

    # Create ground truth
    ground_truth = generator.create_ground_truth(
        patients,
        same_disease_relevant=True
    )

    # Convert to phenopackets
    phenopackets = [p.to_phenopacket() for p in patients]

    # Export files
    generator.export_cohort(
        patients,
        output_dir / "cohort_phenopackets.json",
        format="phenopackets"
    )

    generator.export_cohort(
        patients,
        output_dir / "cohort_patients.json",
        format="json"
    )

    generator.export_cohort(
        patients,
        output_dir / "cohort_patients.csv",
        format="csv"
    )

    # Save ground truth
    with open(output_dir / "ground_truth.json", 'w') as f:
        json.dump(ground_truth, f, indent=2)

    # Save metadata
    with open(output_dir / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    # Save dataset statistics
    stats = generator.get_statistics()
    with open(output_dir / "hpoa_statistics.json", 'w') as f:
        json.dump(stats, f, indent=2)

    logger.info(f"Created evaluation dataset in {output_dir}")

    return phenopackets, ground_truth, metadata


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    print("=" * 70)
    print("  HPOA Patient Generator Demo")
    print("=" * 70)

    # Initialize generator
    generator = HPOAPatientGenerator(
        hpoa_path="data/hpo_annotations/phenotype.hpoa",
        seed=42
    )

    # Load and show statistics
    print("\nLoading HPOA annotations...")
    generator.load_annotations()

    stats = generator.get_statistics()
    print(f"\nDataset Statistics:")
    print(f"  Total diseases: {stats['total_diseases']}")
    print(f"  Unique phenotypes: {stats['total_unique_phenotypes']}")
    print(f"\nDiseases by source:")
    for source, count in stats['diseases_by_source'].items():
        print(f"  {source}: {count}")
    print(f"\nPhenotypes per disease:")
    for key, val in stats['phenotypes_per_disease'].items():
        print(f"  {key}: {val:.1f}" if isinstance(val, float) else f"  {key}: {val}")

    # Generate cohort
    print("\nGenerating evaluation cohort...")
    patients, metadata = generator.generate_cohort(
        n_patients=300,
        n_diseases=60,
        min_phenotypes=5,
        max_phenotypes=25,
        phenotype_recall=0.75,
        noise_rate=0.1
    )

    print(f"\nGenerated cohort:")
    print(f"  Patients: {metadata['n_patients']}")
    print(f"  Diseases: {metadata['n_diseases']}")
    print(f"  Phenotypes/patient: {metadata['phenotypes_per_patient']['mean']:.1f}")

    # Show example patients
    print("\nExample patients:")
    for p in patients[:3]:
        print(f"\n  {p.patient_id}:")
        print(f"    Disease: {p.disease_name[:50]}...")
        print(f"    Phenotypes: {len(p.phenotypes)} ({len(p.noise_phenotypes)} noise)")
        print(f"    Sample: {', '.join(p.phenotypes[:3])}...")

    # Create ground truth
    ground_truth = generator.create_ground_truth(patients)
    relevant_counts = [len(v) for v in ground_truth.values()]
    print(f"\nGround truth (same-disease relevance):")
    print(f"  Mean relevant per patient: {sum(relevant_counts)/len(relevant_counts):.1f}")
    print(f"  Min: {min(relevant_counts)}, Max: {max(relevant_counts)}")

    # Export
    print("\nExporting cohort...")
    generator.export_cohort(
        patients,
        "data/hpoa_evaluation/demo_cohort.json",
        format="phenopackets"
    )

    print("\nDone!")
