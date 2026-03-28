# Quick Start Guide

This guide will help you get the project up and running quickly.

## Prerequisites

- Python 3.9 or higher
- pip package manager
- Virtual environment (recommended)

## Setup Steps

### 1. Create and Activate Virtual Environment

```bash
# Navigate to project directory
cd privacy-phenotype-matching

# Create virtual environment
python3 -m venv venv

# Activate virtual environment
# On macOS/Linux:
source venv/bin/activate

# On Windows:
# venv\Scripts\activate
```

### 2. Install Dependencies

```bash
# Upgrade pip
pip install --upgrade pip

# Install required packages
pip install -r requirements.txt

# Install package in development mode
pip install -e .
```

### 3. Download HPO Ontology Data

```bash
# Run the download script
python scripts/download_hpo.py
```

This will download:
- HPO ontology (hp.obo) - the complete Human Phenotype Ontology
- HPO annotations (phenotype.hpoa) - disease-phenotype associations

Files will be saved to `data/hpo_ontology/`.

### 4. Generate Synthetic Phenopackets

```bash
# Generate a test cohort
python -m src.data_generation.synthetic_phenopackets --cohort-size 100 --output data/synthetic/test_cohort.json
```

### 5. Run Tests

```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=src tests/
```

## Project Structure Overview

```
privacy-phenotype-matching/
├── src/                    # Source code
│   ├── data_generation/   # Synthetic data generation
│   ├── similarity/        # Similarity metrics
│   ├── privacy/           # Privacy-preserving protocols
│   ├── evaluation/        # Evaluation and metrics
│   └── utils/             # Utilities (HPO, phenopackets)
├── data/                  # Data files
├── experiments/           # Jupyter notebooks for experiments
├── tests/                 # Unit tests
└── scripts/               # Utility scripts
```

## Next Steps

1. **Explore HPO Data**: Open a Python shell and try:
   ```python
   from src.utils.hpo_utils import HPOManager

   manager = HPOManager()
   ontology = manager.load_ontology()

   # Search for terms
   results = manager.search_terms("seizure")
   for term in results:
       print(f"{term.id}: {term.name}")
   ```

2. **Generate Synthetic Data**: Create synthetic phenopackets for testing

3. **Compute Similarity**: Try baseline similarity metrics

4. **Experiment with Privacy**: Test PSI-based matching

## Troubleshooting

### pronto installation fails
```bash
pip install --upgrade pip setuptools wheel
pip install pronto
```

### HPO download fails
- Check internet connection
- Try manual download from: https://purl.obolibrary.org/obo/hp.obo
- Place file in `data/hpo_ontology/hp.obo`

### Import errors
- Make sure virtual environment is activated
- Verify installation: `pip list | grep pronto`
- Reinstall: `pip install -e .`

## Running Experiments

Jupyter notebooks for experiments are in `experiments/`. Start Jupyter:

```bash
jupyter notebook experiments/
```

## Getting Help

- Check documentation in `docs/`
- Review example notebooks in `experiments/`
- See test files in `tests/` for usage examples

## Development Workflow

1. Make changes to code
2. Run tests: `pytest tests/`
3. Format code: `black src/`
4. Check style: `flake8 src/`
5. Commit changes

Happy researching!
