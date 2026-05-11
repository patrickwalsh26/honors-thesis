"""
Data integration modules for external rare disease data sources.

This package provides loaders and converters for various rare disease
databases, converting them to GA4GH Phenopacket format for use in
privacy-preserving phenotype matching experiments.

Supported Data Sources:
- DECIPHER: Database of Genomic Variation and Phenotype in Humans
- HPO Annotations: phenotype.hpoa disease-phenotype associations
- Orphanet: European rare disease reference database
- OMIM: Online Mendelian Inheritance in Man
- ClinVar: Clinical variant database

GA4GH Standard Protocols:
- Beacon v2: Count and boolean queries
- Matchmaker Exchange: Federated patient matching
"""

from .decipher_loader import (
    DECIPHERLoader,
    DECIPHERPatient,
    DECIPHERSimulator,
    create_decipher_evaluation_dataset,
)

from .external_sources import (
    HPOAnnotationsLoader,
    DiseaseAnnotation,
    OMIMClient,
    OrphanetClient,
    ClinVarClient,
    DataSourceAggregator,
    convert_omim_to_phenopacket,
    export_to_matchmaker_exchange,
)

from .hpoa_patient_generator import (
    HPOAPatientGenerator,
    DiseaseProfile,
    SyntheticPatient,
    create_evaluation_dataset,
)

from .beacon_v2 import BeaconAdapter
from .matchmaker_exchange import MMEAdapter

__all__ = [
    # DECIPHER
    "DECIPHERLoader",
    "DECIPHERPatient",
    "DECIPHERSimulator",
    "create_decipher_evaluation_dataset",
    # External sources
    "HPOAnnotationsLoader",
    "DiseaseAnnotation",
    "OMIMClient",
    "OrphanetClient",
    "ClinVarClient",
    "DataSourceAggregator",
    "convert_omim_to_phenopacket",
    "export_to_matchmaker_exchange",
    # HPOA Patient Generator (real disease data)
    "HPOAPatientGenerator",
    "DiseaseProfile",
    "SyntheticPatient",
    "create_evaluation_dataset",
    # GA4GH protocols
    "BeaconAdapter",
    "MMEAdapter",
]
