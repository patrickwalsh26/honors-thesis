#!/usr/bin/env python3
"""
Script to download HPO ontology and annotation data.

Usage:
    python scripts/download_hpo.py
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.utils.hpo_utils import HPOManager, download_hpo_annotations
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

logger = logging.getLogger(__name__)


def main():
    """Download HPO data files."""
    logger.info("Starting HPO data download...")

    # Download ontology
    logger.info("=" * 60)
    logger.info("Downloading HPO ontology")
    logger.info("=" * 60)

    manager = HPOManager()
    obo_path = manager.download_hpo(force=False)
    logger.info(f"HPO ontology saved to: {obo_path}")

    # Load to verify
    logger.info("\nVerifying ontology can be loaded...")
    ontology = manager.load_ontology()
    logger.info(f"✓ Successfully loaded {len(ontology)} HPO terms")

    # Download annotations
    logger.info("\n" + "=" * 60)
    logger.info("Downloading HPO annotations")
    logger.info("=" * 60)

    annot_path = download_hpo_annotations()
    logger.info(f"HPO annotations saved to: {annot_path}")

    # Print some stats
    logger.info("\n" + "=" * 60)
    logger.info("HPO Ontology Statistics")
    logger.info("=" * 60)
    logger.info(f"Total terms: {len(ontology)}")

    # Count top-level categories
    root_term = manager.get_term("HP:0000001")  # All
    if root_term:
        logger.info(f"Root term: {root_term.name}")

    # Example terms
    logger.info("\nExample HPO terms:")
    example_ids = ["HP:0000118", "HP:0000707", "HP:0001507", "HP:0002664"]
    for term_id in example_ids:
        term = manager.get_term(term_id)
        if term:
            logger.info(f"  {term_id}: {term.name}")

    logger.info("\n✓ HPO data download complete!")
    logger.info("\nNext steps:")
    logger.info("  1. Generate synthetic phenopackets")
    logger.info("  2. Implement similarity metrics")
    logger.info("  3. Run baseline experiments")


if __name__ == "__main__":
    main()
