"""
Utilities for working with the Human Phenotype Ontology (HPO).

This module provides functions to download, parse, and work with HPO terms.
"""

import os
import urllib.request
from pathlib import Path
from typing import Dict, List, Set, Optional
import logging

try:
    import pronto
except ImportError:
    pronto = None
    logging.warning("pronto not installed. Install with: pip install pronto")

logger = logging.getLogger(__name__)


class HPOManager:
    """Manager for HPO ontology data."""

    def __init__(self, obo_path: Optional[str] = None):
        """
        Initialize HPO Manager.

        Args:
            obo_path: Path to HPO OBO file. If None, will use default location.
        """
        if obo_path is None:
            obo_path = "data/hpo_ontology/hp.obo"
        self.obo_path = Path(obo_path)
        self.ontology = None
        self._ic_cache = {}

    def download_hpo(self, force: bool = False) -> Path:
        """
        Download the latest HPO ontology file.

        Args:
            force: If True, download even if file exists.

        Returns:
            Path to downloaded OBO file.
        """
        url = "https://purl.obolibrary.org/obo/hp.obo"

        # Create directory if it doesn't exist
        self.obo_path.parent.mkdir(parents=True, exist_ok=True)

        if self.obo_path.exists() and not force:
            logger.info(f"HPO ontology already exists at {self.obo_path}")
            return self.obo_path

        logger.info(f"Downloading HPO ontology from {url}")
        try:
            urllib.request.urlretrieve(url, self.obo_path)
            logger.info(f"Successfully downloaded HPO to {self.obo_path}")
        except Exception as e:
            logger.error(f"Failed to download HPO: {e}")
            raise

        return self.obo_path

    def load_ontology(self) -> pronto.Ontology:
        """
        Load the HPO ontology using pronto.

        Returns:
            Loaded ontology object.
        """
        if self.ontology is not None:
            return self.ontology

        if pronto is None:
            raise ImportError("pronto is required. Install with: pip install pronto")

        if not self.obo_path.exists():
            logger.info("HPO file not found, downloading...")
            self.download_hpo()

        logger.info(f"Loading HPO ontology from {self.obo_path}")
        self.ontology = pronto.Ontology(str(self.obo_path))
        logger.info(f"Loaded {len(self.ontology)} HPO terms")

        return self.ontology

    def get_term(self, term_id: str) -> Optional[pronto.Term]:
        """
        Get an HPO term by ID.

        Args:
            term_id: HPO term ID (e.g., 'HP:0000001')

        Returns:
            Term object or None if not found.
        """
        if self.ontology is None:
            self.load_ontology()

        try:
            return self.ontology[term_id]
        except KeyError:
            logger.warning(f"Term {term_id} not found in ontology")
            return None

    def get_ancestors(self, term_id: str, include_self: bool = False) -> Set[str]:
        """
        Get all ancestor terms for a given HPO term.

        Args:
            term_id: HPO term ID
            include_self: Whether to include the term itself

        Returns:
            Set of ancestor term IDs.
        """
        term = self.get_term(term_id)
        if term is None:
            return set()

        ancestors = {str(ancestor.id) for ancestor in term.superclasses()}

        if include_self:
            ancestors.add(term_id)

        return ancestors

    def get_descendants(self, term_id: str, include_self: bool = False) -> Set[str]:
        """
        Get all descendant terms for a given HPO term.

        Args:
            term_id: HPO term ID
            include_self: Whether to include the term itself

        Returns:
            Set of descendant term IDs.
        """
        term = self.get_term(term_id)
        if term is None:
            return set()

        descendants = {str(desc.id) for desc in term.subclasses()}

        if include_self:
            descendants.add(term_id)

        return descendants

    def compute_information_content(
        self,
        term_frequencies: Optional[Dict[str, int]] = None
    ) -> Dict[str, float]:
        """
        Compute information content (IC) for all HPO terms.

        IC(t) = -log(P(t)) where P(t) is the probability of observing term t.

        Args:
            term_frequencies: Dictionary mapping term IDs to occurrence counts.
                            If None, uses uniform distribution over all terms.

        Returns:
            Dictionary mapping term IDs to IC values.
        """
        if self.ontology is None:
            self.load_ontology()

        ic = {}

        # If no frequencies provided, use uniform distribution
        if term_frequencies is None:
            num_terms = len(self.ontology)
            for term in self.ontology.terms():
                ic[str(term.id)] = -1 * (1 / num_terms)
            return ic

        # Compute total count
        total_count = sum(term_frequencies.values())

        # Compute IC for each term
        for term in self.ontology.terms():
            term_id = str(term.id)

            # Get count for this term and all its descendants
            count = term_frequencies.get(term_id, 0)
            for desc_id in self.get_descendants(term_id):
                count += term_frequencies.get(desc_id, 0)

            # Avoid log(0)
            if count == 0:
                ic[term_id] = float('inf')
            else:
                prob = count / total_count
                ic[term_id] = -1 * (prob) if prob > 0 else float('inf')

        self._ic_cache = ic
        return ic

    def get_common_ancestors(self, term_id1: str, term_id2: str) -> Set[str]:
        """
        Get common ancestors of two HPO terms.

        Args:
            term_id1: First HPO term ID
            term_id2: Second HPO term ID

        Returns:
            Set of common ancestor term IDs.
        """
        ancestors1 = self.get_ancestors(term_id1, include_self=True)
        ancestors2 = self.get_ancestors(term_id2, include_self=True)

        return ancestors1 & ancestors2

    def get_most_informative_common_ancestor(
        self,
        term_id1: str,
        term_id2: str,
        ic: Optional[Dict[str, float]] = None
    ) -> Optional[str]:
        """
        Get the most informative common ancestor (MICA) of two terms.

        Args:
            term_id1: First HPO term ID
            term_id2: Second HPO term ID
            ic: Information content dictionary. If None, uses cached IC.

        Returns:
            Term ID of MICA or None if no common ancestor.
        """
        common_ancestors = self.get_common_ancestors(term_id1, term_id2)

        if not common_ancestors:
            return None

        if ic is None:
            ic = self._ic_cache
            if not ic:
                logger.warning("No IC values available, computing with uniform distribution")
                ic = self.compute_information_content()

        # Return ancestor with maximum IC
        mica = max(common_ancestors, key=lambda t: ic.get(t, 0))
        return mica

    def get_term_name(self, term_id: str) -> str:
        """
        Get human-readable name for an HPO term.

        Args:
            term_id: HPO term ID

        Returns:
            Term name or the ID if not found.
        """
        term = self.get_term(term_id)
        return term.name if term else term_id

    def search_terms(self, query: str, limit: int = 10) -> List[pronto.Term]:
        """
        Search for HPO terms by name.

        Args:
            query: Search query
            limit: Maximum number of results

        Returns:
            List of matching terms.
        """
        if self.ontology is None:
            self.load_ontology()

        query_lower = query.lower()
        matches = []

        for term in self.ontology.terms():
            if query_lower in term.name.lower():
                matches.append(term)
                if len(matches) >= limit:
                    break

        return matches


def download_hpo_annotations(output_path: str = "data/hpo_ontology/phenotype.hpoa") -> Path:
    """
    Download HPO phenotype annotations (disease-phenotype associations).

    Args:
        output_path: Where to save the annotation file.

    Returns:
        Path to downloaded file.
    """
    url = "http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"
    output_path = Path(output_path)

    # Create directory if needed
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        logger.info(f"Annotations already exist at {output_path}")
        return output_path

    logger.info(f"Downloading HPO annotations from {url}")
    try:
        urllib.request.urlretrieve(url, output_path)
        logger.info(f"Successfully downloaded annotations to {output_path}")
    except Exception as e:
        logger.error(f"Failed to download annotations: {e}")
        raise

    return output_path


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    manager = HPOManager()
    manager.download_hpo()
    ontology = manager.load_ontology()

    print(f"Loaded {len(ontology)} HPO terms")

    # Example: get term info
    term = manager.get_term("HP:0000001")
    if term:
        print(f"Root term: {term.name}")

    # Example: search
    results = manager.search_terms("seizure", limit=5)
    print(f"\nSearch results for 'seizure':")
    for term in results:
        print(f"  {term.id}: {term.name}")
