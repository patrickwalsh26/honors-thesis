"""
Generate synthetic GA4GH Phenopackets for testing privacy-preserving matching.

This module creates realistic synthetic patient phenotypes using HPO terms
with appropriate prevalence distributions for rare diseases.
"""

import json
import random
import uuid
from datetime import datetime, timedelta
from pathlib import Path
from typing import List, Dict, Optional, Set
from dataclasses import dataclass, asdict
import logging

logger = logging.getLogger(__name__)


@dataclass
class PhenotypicFeature:
    """A phenotypic feature using HPO."""
    type: Dict[str, str]  # {"id": "HP:XXXXXXX", "label": "Feature name"}
    excluded: bool = False
    onset: Optional[Dict[str, str]] = None  # Age of onset


@dataclass
class Individual:
    """An individual/patient."""
    id: str
    sex: str  # MALE, FEMALE, UNKNOWN_SEX
    date_of_birth: Optional[str] = None
    age: Optional[Dict[str, str]] = None


@dataclass
class Disease:
    """A disease diagnosis."""
    term: Dict[str, str]  # {"id": "OMIM:XXXXXX", "label": "Disease name"}
    onset: Optional[Dict[str, str]] = None


@dataclass
class Phenopacket:
    """GA4GH Phenopacket structure (simplified)."""
    id: str
    subject: Individual
    phenotypic_features: List[PhenotypicFeature]
    diseases: List[Disease]
    meta_data: Dict


class RareDiseaseProfile:
    """Profile for a rare disease with characteristic HPO terms."""

    def __init__(
        self,
        disease_id: str,
        disease_name: str,
        core_terms: List[str],
        common_terms: Optional[List[str]] = None,
        rare_terms: Optional[List[str]] = None,
        prevalence: float = 0.0001
    ):
        """
        Initialize a rare disease profile.

        Args:
            disease_id: Disease identifier (e.g., "OMIM:154700")
            disease_name: Human-readable disease name
            core_terms: HPO terms always present (pathognomonic features)
            common_terms: HPO terms frequently present (>50% of patients)
            rare_terms: HPO terms occasionally present (<30% of patients)
            prevalence: Disease prevalence in population
        """
        self.disease_id = disease_id
        self.disease_name = disease_name
        self.core_terms = core_terms
        self.common_terms = common_terms or []
        self.rare_terms = rare_terms or []
        self.prevalence = prevalence


# Predefined rare disease profiles
RARE_DISEASE_PROFILES = {
    "marfan": RareDiseaseProfile(
        disease_id="OMIM:154700",
        disease_name="Marfan syndrome",
        core_terms=["HP:0001519", "HP:0001166", "HP:0002616"],  # Dolichostenomelia, Arachnodactyly, Aortic root aneurysm
        common_terms=[
            "HP:0000268",  # Dolichocephaly
            "HP:0000545",  # Myopia
            "HP:0001065",  # Striae distensae
            "HP:0001382",  # Joint hypermobility
            "HP:0002650",  # Scoliosis
            "HP:0000974",  # Hyperextensible skin
        ],
        rare_terms=[
            "HP:0000175",  # Cleft palate
            "HP:0001634",  # Mitral valve prolapse
            "HP:0002107",  # Pneumothorax
            "HP:0100775",  # Dural ectasia
        ],
        prevalence=0.0002
    ),

    "ehlers_danlos": RareDiseaseProfile(
        disease_id="OMIM:130000",
        disease_name="Ehlers-Danlos syndrome, classic type",
        core_terms=["HP:0000974", "HP:0001382", "HP:0001075"],  # Hyperextensible skin, Joint hypermobility, Atrophic scarring
        common_terms=[
            "HP:0000978",  # Bruising susceptibility
            "HP:0001385",  # Hip dysplasia
            "HP:0002617",  # Vascular fragility
            "HP:0100490",  # Pes planus
            "HP:0001065",  # Striae distensae
        ],
        rare_terms=[
            "HP:0000268",  # Dolichocephaly
            "HP:0002650",  # Scoliosis
            "HP:0001252",  # Muscular hypotonia
            "HP:0001629",  # Ventricular septal defect
        ],
        prevalence=0.0001
    ),

    "achondroplasia": RareDiseaseProfile(
        disease_id="OMIM:100800",
        disease_name="Achondroplasia",
        core_terms=["HP:0003498", "HP:0002007", "HP:0003307"],  # Disproportionate short stature, Frontal bossing, Hyperlordosis
        common_terms=[
            "HP:0000256",  # Macrocephaly
            "HP:0000470",  # Short neck
            "HP:0001385",  # Hip dysplasia
            "HP:0002857",  # Genu varum
            "HP:0002938",  # Lumbar hyperlordosis
            "HP:0005280",  # Depressed nasal bridge
        ],
        rare_terms=[
            "HP:0000160",  # Narrow mouth
            "HP:0002007",  # Frontal bossing
            "HP:0002808",  # Kyphosis
            "HP:0100490",  # Pes planus
        ],
        prevalence=0.00004
    ),

    "progeria": RareDiseaseProfile(
        disease_id="OMIM:176670",
        disease_name="Hutchinson-Gilford progeria syndrome",
        core_terms=["HP:0000230", "HP:0001051", "HP:0001597"],  # Aged appearance, Hypotrichosis, Alopecia
        common_terms=[
            "HP:0000252",  # Microcephaly
            "HP:0000494",  # Downslanted palpebral fissures
            "HP:0000664",  # Synophrys
            "HP:0001597",  # Alopecia
            "HP:0005214",  # Abdominal wall muscle weakness
        ],
        rare_terms=[
            "HP:0000164",  # Abnormality of the dentition
            "HP:0001251",  # Ataxia
            "HP:0001647",  # Bicuspid aortic valve
        ],
        prevalence=0.0000004
    ),
}


class PhenopacketGenerator:
    """Generate synthetic GA4GH Phenopackets."""

    def __init__(
        self,
        hpo_terms_file: Optional[str] = None,
        min_terms: int = 3,
        max_terms: int = 15,
        noise_probability: float = 0.1,
        negation_probability: float = 0.05,
        random_seed: Optional[int] = None
    ):
        """
        Initialize the phenopacket generator.

        Args:
            hpo_terms_file: Path to HPO OBO file (optional)
            min_terms: Minimum number of phenotypic features per patient
            max_terms: Maximum number of phenotypic features per patient
            noise_probability: Probability of adding noise terms
            negation_probability: Probability of negating a term
            random_seed: Random seed for reproducibility
        """
        self.min_terms = min_terms
        self.max_terms = max_terms
        self.noise_probability = noise_probability
        self.negation_probability = negation_probability

        if random_seed is not None:
            random.seed(random_seed)

        # Common HPO noise terms (non-specific symptoms)
        self.noise_terms = [
            ("HP:0001337", "Tremor"),
            ("HP:0002014", "Diarrhea"),
            ("HP:0002315", "Headache"),
            ("HP:0025406", "Asthenia"),
            ("HP:0012378", "Fatigue"),
            ("HP:0002829", "Arthralgia"),
            ("HP:0001324", "Muscle weakness"),
            ("HP:0001251", "Ataxia"),
            ("HP:0002094", "Dyspnea"),
            ("HP:0001945", "Fever"),
        ]

    def generate_phenopacket(
        self,
        disease_profile: RareDiseaseProfile,
        patient_id: Optional[str] = None
    ) -> Phenopacket:
        """
        Generate a single synthetic phenopacket.

        Args:
            disease_profile: Disease profile to generate from
            patient_id: Optional patient ID (will generate if not provided)

        Returns:
            Generated Phenopacket
        """
        if patient_id is None:
            patient_id = f"patient_{uuid.uuid4().hex[:8]}"

        # Generate patient demographics
        sex = random.choice(["MALE", "FEMALE"])
        age_years = random.randint(1, 70)

        subject = Individual(
            id=patient_id,
            sex=sex,
            age={"iso8601duration": f"P{age_years}Y"}
        )

        # Generate phenotypic features
        phenotypic_features = []
        selected_terms: Set[str] = set()

        # Always include core terms
        for term_id in disease_profile.core_terms:
            term_label = self._get_term_label(term_id)
            phenotypic_features.append(
                PhenotypicFeature(
                    type={"id": term_id, "label": term_label},
                    excluded=random.random() < self.negation_probability
                )
            )
            selected_terms.add(term_id)

        # Add common terms with high probability
        for term_id in disease_profile.common_terms:
            if random.random() < 0.7:  # 70% chance
                term_label = self._get_term_label(term_id)
                phenotypic_features.append(
                    PhenotypicFeature(
                        type={"id": term_id, "label": term_label},
                        excluded=random.random() < self.negation_probability
                    )
                )
                selected_terms.add(term_id)

        # Add rare terms with low probability
        for term_id in disease_profile.rare_terms:
            if random.random() < 0.3:  # 30% chance
                term_label = self._get_term_label(term_id)
                phenotypic_features.append(
                    PhenotypicFeature(
                        type={"id": term_id, "label": term_label},
                        excluded=random.random() < self.negation_probability
                    )
                )
                selected_terms.add(term_id)

        # Add noise terms
        num_noise = int(random.random() < self.noise_probability) * random.randint(1, 3)
        for _ in range(num_noise):
            term_id, term_label = random.choice(self.noise_terms)
            if term_id not in selected_terms:
                phenotypic_features.append(
                    PhenotypicFeature(
                        type={"id": term_id, "label": term_label}
                    )
                )
                selected_terms.add(term_id)

        # Ensure we meet min/max requirements
        while len(phenotypic_features) < self.min_terms:
            term_id, term_label = random.choice(self.noise_terms)
            if term_id not in selected_terms:
                phenotypic_features.append(
                    PhenotypicFeature(type={"id": term_id, "label": term_label})
                )
                selected_terms.add(term_id)

        if len(phenotypic_features) > self.max_terms:
            phenotypic_features = random.sample(phenotypic_features, self.max_terms)

        # Create disease entry
        diseases = [
            Disease(
                term={
                    "id": disease_profile.disease_id,
                    "label": disease_profile.disease_name
                }
            )
        ]

        # Create metadata
        meta_data = {
            "created": datetime.now().isoformat(),
            "created_by": "PhenopacketGenerator",
            "resources": [
                {
                    "id": "hp",
                    "name": "Human Phenotype Ontology",
                    "url": "http://purl.obolibrary.org/obo/hp.owl",
                    "version": "2025-11-24",
                    "namespace_prefix": "HP"
                }
            ],
            "phenopacket_schema_version": "2.0"
        }

        return Phenopacket(
            id=f"phenopacket_{patient_id}",
            subject=subject,
            phenotypic_features=phenotypic_features,
            diseases=diseases,
            meta_data=meta_data
        )

    def generate_cohort(
        self,
        disease_key: str,
        n_patients: int,
        heterogeneity: float = 0.2
    ) -> List[Phenopacket]:
        """
        Generate a cohort of patients with a specific disease.

        Args:
            disease_key: Key for disease profile in RARE_DISEASE_PROFILES
            n_patients: Number of patients to generate
            heterogeneity: Amount of phenotypic heterogeneity (0-1)

        Returns:
            List of generated Phenopackets
        """
        if disease_key not in RARE_DISEASE_PROFILES:
            raise ValueError(f"Unknown disease: {disease_key}. Available: {list(RARE_DISEASE_PROFILES.keys())}")

        disease_profile = RARE_DISEASE_PROFILES[disease_key]

        cohort = []
        for i in range(n_patients):
            phenopacket = self.generate_phenopacket(
                disease_profile,
                patient_id=f"{disease_key}_{i:04d}"
            )
            cohort.append(phenopacket)

        logger.info(f"Generated cohort of {n_patients} patients with {disease_profile.disease_name}")
        return cohort

    def generate_mixed_cohort(
        self,
        n_patients: int,
        disease_distribution: Optional[Dict[str, float]] = None
    ) -> List[Phenopacket]:
        """
        Generate a mixed cohort with multiple rare diseases.

        Args:
            n_patients: Total number of patients
            disease_distribution: Dict mapping disease keys to proportions
                                If None, uses uniform distribution

        Returns:
            List of generated Phenopackets
        """
        if disease_distribution is None:
            # Uniform distribution
            diseases = list(RARE_DISEASE_PROFILES.keys())
            disease_distribution = {d: 1.0 / len(diseases) for d in diseases}

        # Normalize distribution
        total = sum(disease_distribution.values())
        disease_distribution = {k: v / total for k, v in disease_distribution.items()}

        cohort = []
        for disease_key, proportion in disease_distribution.items():
            n_disease_patients = int(n_patients * proportion)
            if n_disease_patients > 0:
                disease_cohort = self.generate_cohort(disease_key, n_disease_patients)
                cohort.extend(disease_cohort)

        # Shuffle to mix diseases
        random.shuffle(cohort)

        logger.info(f"Generated mixed cohort of {len(cohort)} patients")
        return cohort

    def save_cohort(
        self,
        cohort: List[Phenopacket],
        output_path: str,
        pretty: bool = True
    ):
        """
        Save cohort to JSON file.

        Args:
            cohort: List of Phenopackets
            output_path: Where to save the file
            pretty: Whether to use pretty printing
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Convert to dictionaries
        cohort_dict = [self._phenopacket_to_dict(p) for p in cohort]

        with open(output_path, 'w') as f:
            if pretty:
                json.dump(cohort_dict, f, indent=2)
            else:
                json.dump(cohort_dict, f)

        logger.info(f"Saved {len(cohort)} phenopackets to {output_path}")

    def _phenopacket_to_dict(self, phenopacket: Phenopacket) -> Dict:
        """Convert Phenopacket to dictionary."""
        return {
            "id": phenopacket.id,
            "subject": asdict(phenopacket.subject),
            "phenotypicFeatures": [asdict(f) for f in phenopacket.phenotypic_features],
            "diseases": [asdict(d) for d in phenopacket.diseases],
            "metaData": phenopacket.meta_data
        }

    def _get_term_label(self, term_id: str) -> str:
        """Get human-readable label for HPO term. Placeholder for now."""
        # In a real implementation, this would query the HPO ontology
        # For now, return a generic label
        return f"Phenotype {term_id}"


def main():
    """Command-line interface for generating phenopackets."""
    import argparse

    parser = argparse.ArgumentParser(description="Generate synthetic phenopackets")
    parser.add_argument("--disease", type=str, default="marfan",
                        choices=list(RARE_DISEASE_PROFILES.keys()),
                        help="Disease to generate")
    parser.add_argument("--cohort-size", type=int, default=100,
                        help="Number of patients to generate")
    parser.add_argument("--output", type=str, default="data/synthetic/cohort.json",
                        help="Output file path")
    parser.add_argument("--mixed", action="store_true",
                        help="Generate mixed cohort with all diseases")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility")

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    generator = PhenopacketGenerator(random_seed=args.seed)

    if args.mixed:
        cohort = generator.generate_mixed_cohort(args.cohort_size)
    else:
        cohort = generator.generate_cohort(args.disease, args.cohort_size)

    generator.save_cohort(cohort, args.output)

    # Print summary
    print(f"\nGenerated {len(cohort)} phenopackets")
    print(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()
