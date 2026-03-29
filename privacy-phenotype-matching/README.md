# Privacy-Preserving Phenotype Matching for Rare Disease Cohorts

**Stanford CS Honors Thesis 2025-2026**
**Author:** Patrick Walsh
**Advisor:** Professor Stephen Montgomery

## Overview

This project develops and evaluates a privacy-preserving discovery tool that enables peer institutions to identify "patients like mine" using standardized phenotypes without sharing raw medical records or protected health information (PHI). The system uses cryptographic protocols and differential privacy to enable secure phenotype matching for rare disease research.

## Key Features

- **Two-Step Reveal Ladder**: Thresholded similarity testing followed by minimal data disclosure
- **Privacy Mechanisms**:
  - Private Set Intersection (PSI) for secure overlap computation
  - Differential Privacy (DP) for aggregate statistics
  - k-anonymity guardrails
  - Rare-term suppression
- **Standards Compliance**: GA4GH Phenopackets, Matchmaker Exchange (MME), Beacon v2
- **Risk-Aware Encoding**: Down-weighting potentially re-identifying HPO terms

## Project Structure

```
privacy-phenotype-matching/
├── data/                      # Data files and ontologies
│   ├── hpo_ontology/         # Human Phenotype Ontology files
│   ├── synthetic/            # Generated synthetic phenopackets
│   └── term_prevalence/      # Published HPO term prevalence data
├── src/                      # Source code
│   ├── data_generation/      # Synthetic phenopacket generator
│   ├── similarity/           # Baseline similarity metrics
│   ├── privacy/              # Privacy-preserving protocols (PSI, DP)
│   ├── evaluation/           # Evaluation metrics and analysis
│   └── utils/                # Utilities (HPO parsing, GA4GH adapters)
├── experiments/              # Experimental notebooks and scripts
├── tests/                    # Unit tests
├── docs/                     # Documentation
└── reports/                  # Progress and final reports
```

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd privacy-phenotype-matching

# Create virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Quick Start

### 1. Generate Synthetic Phenopackets

```python
from src.data_generation.synthetic_phenopackets import PhenopacketGenerator

generator = PhenopacketGenerator()
phenopackets = generator.generate_cohort(n_patients=100, disease="RARE_DISEASE_001")
generator.save_cohort(phenopackets, "data/synthetic/cohort_001.json")
```

### 2. Compute Baseline Similarity

```python
from src.similarity.hpo_similarity import ResnikSimilarity

sim = ResnikSimilarity()
score = sim.compute_similarity(phenopacket_1, phenopacket_2)
```

### 3. Run Privacy-Preserving Matching

```python
from src.privacy.psi_matcher import PSIMatcher

matcher = PSIMatcher(k_anonymity=5)
matches = matcher.find_matches(query_phenopacket, database_phenopackets)
```

## Methodology

### Phenotype Representation
- GA4GH Phenopackets with HPO term vectors
- Information content (IC) weighting
- Optional genomic hints (candidate genes, variant classes)

### Similarity Metrics
- Resnik IC-weighted similarity (deterministic baseline)
- Optional learned HPO embeddings
- Calibrated for top-k retrieval

### Privacy Mechanisms
1. **Private Set Intersection (PSI)**: Compute overlap without revealing non-intersecting terms
2. **Differential Privacy**: Add calibrated noise to counts and rankings
3. **k-anonymity**: Suppress results if fewer than k matches
4. **Rare-term filtering**: Down-weight potentially identifying rare phenotypes

### Evaluation Metrics
- Retrieval quality: Recall@k, nDCG
- Privacy-utility frontier: Performance vs. privacy budget (ε)
- Leakage audits: Attribute and membership inference attacks

## Development Timeline

### Fall Quarter 2024
- [x] Project repository setup
- [ ] Synthetic Phenopacket generator
- [ ] Baseline similarity metrics (Resnik IC)
- [ ] PSI prototype implementation
- [ ] Initial privacy-utility experiments

### Winter Quarter 2025
- [ ] LSH-ANN private retrieval
- [ ] MME/Beacon adapters
- [ ] Containerized service
- [ ] Real data integration (if approved)

### Spring Quarter 2025
- [ ] Comprehensive evaluation
- [ ] Leakage audits
- [ ] Final thesis and open-source release

## Data Sources

### Current (Fall 2024)
- **Synthetic Data**: Generated from HPO term prevalence statistics
- **HPO Ontology**: Public Human Phenotype Ontology
- **Published Prevalence**: Aggregated statistics from literature

### Planned (Winter/Spring 2025)
- **MIMIC-IV**: De-identified ICU EHR data (application pending)
- **STARR-OMOP-deid**: Stanford de-identified EHR (application pending)

## Dependencies

Core libraries:
- `pronto` or `ontobio` - HPO ontology parsing
- `phenopackets` - GA4GH Phenopackets schema
- `numpy`, `scipy` - Numerical computation
- `scikit-learn` - Machine learning utilities
- `cryptography` - Cryptographic primitives for PSI

See `requirements.txt` for complete list.

## Testing

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=src tests/
```

## Citation

If you use this code, please cite:

```
Walsh, P. (2026). Privacy-Preserving Phenotype Matching for Rare Disease Cohorts.
Stanford University CS Honors Thesis.
```

## License


## Contact

Patrick Walsh - walshp26@stanford.edu
Montgomery Lab - Stanford University School of Medicine

## Acknowledgments

- Professor Stephen Montgomery (Advisor)
- Stanford CS Honors Program
- Montgomery Lab members
