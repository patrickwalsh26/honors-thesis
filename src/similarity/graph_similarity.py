"""
Graph-based similarity metrics for HPO phenotype matching.

This module implements advanced semantic similarity metrics that leverage
the HPO ontology's directed acyclic graph (DAG) structure.

These metrics are more sophisticated than set-based approaches and capture
the semantic relationships between phenotype terms.

Reference:
    - Resnik P. (1995). "Using Information Content to Evaluate Semantic Similarity"
    - Lin D. (1998). "An Information-Theoretic Definition of Similarity"
    - Jiang & Conrath (1997). "Semantic Similarity Based on Corpus Statistics"
"""

import numpy as np
from typing import Dict, List, Set, Optional, Tuple
from collections import defaultdict
from functools import lru_cache
import logging

logger = logging.getLogger(__name__)


class OntologyGraphSimilarity:
    """
    Advanced semantic similarity using HPO ontology graph structure.

    This class implements proper Resnik, Lin, and Jiang-Conrath similarity
    metrics using the full HPO ancestor graph.
    """

    def __init__(
        self,
        hpo_manager=None,
        ic_method: str = "intrinsic"
    ):
        """
        Initialize graph-based similarity calculator.

        Args:
            hpo_manager: HPOManager instance with loaded ontology
            ic_method: Method for computing IC ("intrinsic", "corpus", "hybrid")
        """
        self.hpo_manager = hpo_manager
        self.ic_method = ic_method
        self.ic_values: Dict[str, float] = {}
        self._ancestor_cache: Dict[str, Set[str]] = {}

        if hpo_manager is not None:
            self._initialize_from_ontology()

    def _initialize_from_ontology(self):
        """Pre-compute IC values and cache ancestors."""
        if self.hpo_manager is None or self.hpo_manager.ontology is None:
            logger.warning("No ontology loaded, using intrinsic IC")
            return

        logger.info("Computing intrinsic information content from ontology structure...")

        # Count descendants for each term (intrinsic IC)
        ontology = self.hpo_manager.ontology
        total_terms = len(list(ontology.terms()))

        for term in ontology.terms():
            term_id = str(term.id)

            # Cache ancestors
            ancestors = self.hpo_manager.get_ancestors(term_id, include_self=True)
            self._ancestor_cache[term_id] = ancestors

            # Compute intrinsic IC based on descendants
            descendants = self.hpo_manager.get_descendants(term_id, include_self=True)
            n_descendants = len(descendants)

            # IC = -log(P(term)) where P(term) = descendants / total
            if n_descendants > 0:
                prob = n_descendants / total_terms
                self.ic_values[term_id] = -np.log(prob)
            else:
                self.ic_values[term_id] = 0.0

        logger.info(f"Computed IC for {len(self.ic_values)} terms")

    def set_corpus_ic(self, phenopackets: List[Dict]):
        """
        Compute IC from a corpus of phenopackets.

        Args:
            phenopackets: List of phenopacket dictionaries
        """
        term_counts = defaultdict(int)
        total_patients = len(phenopackets)

        for pp in phenopackets:
            observed_terms = set()
            for feature in pp.get("phenotypicFeatures", []):
                if not feature.get("excluded", False):
                    term_id = feature["type"]["id"]
                    observed_terms.add(term_id)

                    # Also count ancestors (propagation)
                    if term_id in self._ancestor_cache:
                        observed_terms.update(self._ancestor_cache[term_id])

            for term in observed_terms:
                term_counts[term] += 1

        # Compute IC
        for term_id, count in term_counts.items():
            freq = count / total_patients
            self.ic_values[term_id] = -np.log(freq) if freq > 0 else 0.0

        logger.info(f"Updated IC from corpus ({total_patients} patients)")

    def get_ic(self, term_id: str) -> float:
        """Get information content for a term."""
        return self.ic_values.get(term_id, 0.0)

    def get_common_ancestors(self, term1: str, term2: str) -> Set[str]:
        """Get common ancestors of two terms."""
        anc1 = self._ancestor_cache.get(term1, {term1})
        anc2 = self._ancestor_cache.get(term2, {term2})
        return anc1 & anc2

    def get_mica(self, term1: str, term2: str) -> Tuple[Optional[str], float]:
        """
        Get Most Informative Common Ancestor (MICA).

        Returns:
            Tuple of (MICA term ID, IC value)
        """
        common_ancestors = self.get_common_ancestors(term1, term2)

        if not common_ancestors:
            return None, 0.0

        # Find ancestor with maximum IC
        mica = max(common_ancestors, key=lambda t: self.get_ic(t))
        return mica, self.get_ic(mica)

    def resnik_similarity(self, term1: str, term2: str) -> float:
        """
        Resnik semantic similarity.

        sim(t1, t2) = IC(MICA(t1, t2))
        """
        _, ic = self.get_mica(term1, term2)
        return ic

    def lin_similarity(self, term1: str, term2: str) -> float:
        """
        Lin semantic similarity.

        sim(t1, t2) = 2 * IC(MICA) / (IC(t1) + IC(t2))
        """
        _, mica_ic = self.get_mica(term1, term2)
        ic1 = self.get_ic(term1)
        ic2 = self.get_ic(term2)

        if ic1 + ic2 == 0:
            return 0.0

        return (2 * mica_ic) / (ic1 + ic2)

    def jiang_conrath_distance(self, term1: str, term2: str) -> float:
        """
        Jiang-Conrath semantic distance.

        d(t1, t2) = IC(t1) + IC(t2) - 2 * IC(MICA)

        Lower is more similar. Convert to similarity with 1/(1+d).
        """
        _, mica_ic = self.get_mica(term1, term2)
        ic1 = self.get_ic(term1)
        ic2 = self.get_ic(term2)

        distance = ic1 + ic2 - 2 * mica_ic
        return distance

    def jiang_conrath_similarity(self, term1: str, term2: str) -> float:
        """Convert JC distance to similarity."""
        d = self.jiang_conrath_distance(term1, term2)
        return 1.0 / (1.0 + d)

    def compute_termset_similarity(
        self,
        terms1: List[str],
        terms2: List[str],
        method: str = "resnik",
        aggregation: str = "bma"
    ) -> float:
        """
        Compute similarity between two sets of HPO terms.

        Args:
            terms1: First set of HPO term IDs
            terms2: Second set of HPO term IDs
            method: Similarity method ("resnik", "lin", "jc")
            aggregation: Aggregation method ("bma", "max", "average")

        Returns:
            Similarity score
        """
        if not terms1 or not terms2:
            return 0.0

        # Select similarity function
        sim_func = {
            "resnik": self.resnik_similarity,
            "lin": self.lin_similarity,
            "jc": self.jiang_conrath_similarity
        }.get(method, self.resnik_similarity)

        if aggregation == "bma":
            return self._best_match_average(terms1, terms2, sim_func)
        elif aggregation == "max":
            return self._max_similarity(terms1, terms2, sim_func)
        elif aggregation == "average":
            return self._average_similarity(terms1, terms2, sim_func)
        else:
            raise ValueError(f"Unknown aggregation: {aggregation}")

    def _best_match_average(self, terms1: List[str], terms2: List[str], sim_func) -> float:
        """Best-match average (BMA) aggregation."""
        # Forward direction
        scores_1to2 = []
        for t1 in terms1:
            best = max((sim_func(t1, t2) for t2 in terms2), default=0.0)
            scores_1to2.append(best)

        # Reverse direction
        scores_2to1 = []
        for t2 in terms2:
            best = max((sim_func(t1, t2) for t1 in terms1), default=0.0)
            scores_2to1.append(best)

        all_scores = scores_1to2 + scores_2to1
        return np.mean(all_scores) if all_scores else 0.0

    def _max_similarity(self, terms1: List[str], terms2: List[str], sim_func) -> float:
        """Maximum pairwise similarity."""
        max_sim = 0.0
        for t1 in terms1:
            for t2 in terms2:
                sim = sim_func(t1, t2)
                max_sim = max(max_sim, sim)
        return max_sim

    def _average_similarity(self, terms1: List[str], terms2: List[str], sim_func) -> float:
        """Average of all pairwise similarities."""
        scores = []
        for t1 in terms1:
            for t2 in terms2:
                scores.append(sim_func(t1, t2))
        return np.mean(scores) if scores else 0.0


class GraphAwarePhenopacketCalculator:
    """
    High-level interface for graph-aware phenopacket similarity.

    This calculator uses the full HPO ontology structure for more
    accurate semantic similarity computation.
    """

    def __init__(
        self,
        graph_similarity: OntologyGraphSimilarity,
        method: str = "resnik",
        aggregation: str = "bma"
    ):
        """
        Initialize calculator.

        Args:
            graph_similarity: OntologyGraphSimilarity instance
            method: Similarity method ("resnik", "lin", "jc")
            aggregation: Aggregation method ("bma", "max", "average")
        """
        self.graph_similarity = graph_similarity
        self.method = method
        self.aggregation = aggregation
        self._cache = {}

    def extract_terms(self, phenopacket: Dict) -> List[str]:
        """Extract HPO terms from a phenopacket."""
        terms = []
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                term_id = feature["type"]["id"]
                terms.append(term_id)
        return terms

    def compute_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict
    ) -> float:
        """Compute similarity between two phenopackets."""
        id1 = phenopacket1.get("id", "")
        id2 = phenopacket2.get("id", "")

        cache_key = tuple(sorted([id1, id2]))
        if cache_key in self._cache:
            return self._cache[cache_key]

        terms1 = self.extract_terms(phenopacket1)
        terms2 = self.extract_terms(phenopacket2)

        sim = self.graph_similarity.compute_termset_similarity(
            terms1, terms2,
            method=self.method,
            aggregation=self.aggregation
        )

        self._cache[cache_key] = sim
        return sim

    def find_most_similar(
        self,
        query: Dict,
        candidates: List[Dict],
        top_k: int = 10
    ) -> List[Tuple[int, float]]:
        """Find most similar phenopackets to a query."""
        similarities = []

        for idx, candidate in enumerate(candidates):
            sim = self.compute_similarity(query, candidate)
            similarities.append((idx, sim))

        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities[:top_k]

    def compute_similarity_matrix(self, phenopackets: List[Dict]) -> np.ndarray:
        """Compute all-vs-all similarity matrix."""
        n = len(phenopackets)
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    matrix[i, j] = 1.0
                else:
                    sim = self.compute_similarity(phenopackets[i], phenopackets[j])
                    matrix[i, j] = sim
                    matrix[j, i] = sim

        return matrix


class PhenotypicOverlapScore:
    """
    Compute phenotypic overlap score between patients and diseases.

    This metric is commonly used in rare disease diagnosis tools like
    Exomiser and LIRICAL.
    """

    def __init__(self, graph_similarity: OntologyGraphSimilarity):
        """
        Initialize overlap score calculator.

        Args:
            graph_similarity: OntologyGraphSimilarity instance
        """
        self.graph_similarity = graph_similarity

    def compute_patient_disease_score(
        self,
        patient_terms: List[str],
        disease_terms: List[str]
    ) -> Tuple[float, Dict]:
        """
        Compute phenotypic overlap score between patient and disease.

        Returns:
            Tuple of (score, details)
        """
        if not patient_terms or not disease_terms:
            return 0.0, {}

        # For each patient term, find best matching disease term
        term_matches = []
        for p_term in patient_terms:
            best_match = None
            best_score = 0.0

            for d_term in disease_terms:
                sim = self.graph_similarity.resnik_similarity(p_term, d_term)
                if sim > best_score:
                    best_score = sim
                    best_match = d_term

            term_matches.append({
                "patient_term": p_term,
                "matched_disease_term": best_match,
                "score": best_score
            })

        # Aggregate scores
        scores = [m["score"] for m in term_matches]

        # Geometric mean (as in LIRICAL)
        if scores and all(s > 0 for s in scores):
            geometric_mean = np.exp(np.mean(np.log(scores)))
        else:
            geometric_mean = 0.0

        # Also compute coverage
        matched_count = sum(1 for s in scores if s > 0)
        coverage = matched_count / len(patient_terms)

        return geometric_mean, {
            "term_matches": term_matches,
            "coverage": coverage,
            "arithmetic_mean": np.mean(scores),
            "geometric_mean": geometric_mean
        }


def create_graph_calculator(obo_path: str = None) -> GraphAwarePhenopacketCalculator:
    """
    Factory function to create a graph-aware calculator.

    Args:
        obo_path: Path to HPO OBO file

    Returns:
        Configured GraphAwarePhenopacketCalculator
    """
    import sys
    from pathlib import Path

    # Try to import HPOManager
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))
        from src.utils.hpo_utils import HPOManager

        manager = HPOManager(obo_path)
        manager.load_ontology()

        graph_sim = OntologyGraphSimilarity(manager)
        return GraphAwarePhenopacketCalculator(graph_sim)

    except Exception as e:
        logger.warning(f"Could not create graph calculator: {e}")
        logger.warning("Falling back to basic calculator without ontology")

        graph_sim = OntologyGraphSimilarity(None)
        return GraphAwarePhenopacketCalculator(graph_sim)


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("Creating graph-aware similarity calculator...")
    calc = create_graph_calculator()

    # Test with synthetic terms
    terms1 = ["HP:0001519", "HP:0001166", "HP:0002616"]  # Marfan terms
    terms2 = ["HP:0001519", "HP:0002650", "HP:0000545"]  # Overlapping Marfan terms
    terms3 = ["HP:0003498", "HP:0002007", "HP:0003307"]  # Achondroplasia terms

    print("\nSimilarity tests:")

    sim_same = calc.graph_similarity.compute_termset_similarity(terms1, terms2)
    print(f"Marfan vs Marfan (overlap): {sim_same:.4f}")

    sim_diff = calc.graph_similarity.compute_termset_similarity(terms1, terms3)
    print(f"Marfan vs Achondroplasia: {sim_diff:.4f}")
