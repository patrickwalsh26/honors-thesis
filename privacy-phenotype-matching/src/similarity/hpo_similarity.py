"""
Similarity metrics for comparing HPO phenotype profiles.

Implements various similarity measures including:
- Resnik similarity (IC-weighted)
- Lin similarity
- Jaccard similarity
- Cosine similarity
"""

import numpy as np
from typing import List, Dict, Set, Optional, Tuple
from collections import defaultdict
import logging
import json

logger = logging.getLogger(__name__)


class HPOSimilarity:
    """Base class for HPO similarity metrics."""

    def __init__(self, ic_values: Optional[Dict[str, float]] = None):
        """
        Initialize similarity calculator.

        Args:
            ic_values: Information content values for HPO terms
        """
        self.ic_values = ic_values or {}

    def compute_pairwise_similarity(
        self,
        terms1: List[str],
        terms2: List[str]
    ) -> float:
        """
        Compute similarity between two sets of HPO terms.

        Args:
            terms1: First set of HPO term IDs
            terms2: Second set of HPO term IDs

        Returns:
            Similarity score (higher is more similar)
        """
        raise NotImplementedError

    def load_phenopacket_terms(self, phenopacket: Dict) -> List[str]:
        """
        Extract HPO terms from a phenopacket.

        Args:
            phenopacket: Phenopacket dictionary

        Returns:
            List of HPO term IDs (excluding negated terms)
        """
        terms = []
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                term_id = feature["type"]["id"]
                terms.append(term_id)
        return terms


class JaccardSimilarity(HPOSimilarity):
    """Jaccard similarity: |A ∩ B| / |A ∪ B|"""

    def compute_pairwise_similarity(
        self,
        terms1: List[str],
        terms2: List[str]
    ) -> float:
        """Compute Jaccard similarity."""
        set1 = set(terms1)
        set2 = set(terms2)

        if not set1 and not set2:
            return 1.0

        if not set1 or not set2:
            return 0.0

        intersection = len(set1 & set2)
        union = len(set1 | set2)

        return intersection / union if union > 0 else 0.0


class CosineSimilarity(HPOSimilarity):
    """Cosine similarity with optional IC weighting."""

    def compute_pairwise_similarity(
        self,
        terms1: List[str],
        terms2: List[str]
    ) -> float:
        """Compute cosine similarity."""
        # Create term frequency vectors
        all_terms = list(set(terms1 + terms2))

        if not all_terms:
            return 1.0

        # Build vectors
        vec1 = np.zeros(len(all_terms))
        vec2 = np.zeros(len(all_terms))

        term_to_idx = {term: i for i, term in enumerate(all_terms)}

        for term in terms1:
            idx = term_to_idx[term]
            weight = self.ic_values.get(term, 1.0) if self.ic_values else 1.0
            vec1[idx] = weight

        for term in terms2:
            idx = term_to_idx[term]
            weight = self.ic_values.get(term, 1.0) if self.ic_values else 1.0
            vec2[idx] = weight

        # Compute cosine similarity
        dot_product = np.dot(vec1, vec2)
        norm1 = np.linalg.norm(vec1)
        norm2 = np.linalg.norm(vec2)

        if norm1 == 0 or norm2 == 0:
            return 0.0

        return dot_product / (norm1 * norm2)


class ResnikSimilarity(HPOSimilarity):
    """
    Resnik similarity using information content.

    For two terms, similarity is the IC of their most informative common ancestor (MICA).
    For two sets of terms, we use best-match average.
    """

    def __init__(
        self,
        ic_values: Optional[Dict[str, float]] = None,
        ancestors: Optional[Dict[str, Set[str]]] = None
    ):
        """
        Initialize Resnik similarity.

        Args:
            ic_values: Information content values for HPO terms
            ancestors: Dictionary mapping each term to its ancestors
        """
        super().__init__(ic_values)
        self.ancestors = ancestors or {}

    def get_common_ancestors(self, term1: str, term2: str) -> Set[str]:
        """Get common ancestors of two terms."""
        if not self.ancestors:
            # If no ancestor information, assume only exact matches have common ancestors
            return {term1} if term1 == term2 else set()

        anc1 = self.ancestors.get(term1, {term1})
        anc2 = self.ancestors.get(term2, {term2})

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
        max_ic = 0.0
        mica = None

        for ancestor in common_ancestors:
            ic = self.ic_values.get(ancestor, 0.0)
            if ic > max_ic:
                max_ic = ic
                mica = ancestor

        return mica, max_ic

    def compute_term_similarity(self, term1: str, term2: str) -> float:
        """Compute Resnik similarity between two terms."""
        _, ic = self.get_mica(term1, term2)
        return ic

    def compute_pairwise_similarity(
        self,
        terms1: List[str],
        terms2: List[str],
        method: str = "bma"  # best-match average
    ) -> float:
        """
        Compute Resnik similarity between two term sets.

        Args:
            terms1: First set of HPO terms
            terms2: Second set of HPO terms
            method: Aggregation method ("bma" for best-match average)

        Returns:
            Similarity score
        """
        if not terms1 or not terms2:
            return 0.0

        if method == "bma":
            return self._best_match_average(terms1, terms2)
        elif method == "average":
            return self._average_similarity(terms1, terms2)
        else:
            raise ValueError(f"Unknown method: {method}")

    def _best_match_average(self, terms1: List[str], terms2: List[str]) -> float:
        """
        Best-match average (BMA) similarity.

        For each term in set1, find best matching term in set2, and vice versa.
        Return average of these best matches.
        """
        # For each term in terms1, find best match in terms2
        scores_1to2 = []
        for t1 in terms1:
            best_score = max(
                (self.compute_term_similarity(t1, t2) for t2 in terms2),
                default=0.0
            )
            scores_1to2.append(best_score)

        # For each term in terms2, find best match in terms1
        scores_2to1 = []
        for t2 in terms2:
            best_score = max(
                (self.compute_term_similarity(t1, t2) for t1 in terms1),
                default=0.0
            )
            scores_2to1.append(best_score)

        # Average of all best matches
        all_scores = scores_1to2 + scores_2to1
        return np.mean(all_scores) if all_scores else 0.0

    def _average_similarity(self, terms1: List[str], terms2: List[str]) -> float:
        """Average of all pairwise similarities."""
        scores = []
        for t1 in terms1:
            for t2 in terms2:
                scores.append(self.compute_term_similarity(t1, t2))

        return np.mean(scores) if scores else 0.0


class SimplifiedResnikSimilarity(HPOSimilarity):
    """
    Simplified Resnik similarity without requiring full ontology.

    Uses only exact term matches and IC weighting.
    Useful for initial testing without loading full HPO structure.
    """

    def compute_pairwise_similarity(
        self,
        terms1: List[str],
        terms2: List[str]
    ) -> float:
        """Compute simplified Resnik similarity."""
        set1 = set(terms1)
        set2 = set(terms2)

        if not set1 or not set2:
            return 0.0

        # Find intersection
        common_terms = set1 & set2

        if not common_terms:
            return 0.0

        # Weight by information content
        if self.ic_values:
            # Sum IC of common terms
            ic_sum = sum(self.ic_values.get(term, 1.0) for term in common_terms)
            # Normalize by average IC of all terms
            all_ic_sum = sum(self.ic_values.get(term, 1.0) for term in set1 | set2)
            return ic_sum / all_ic_sum if all_ic_sum > 0 else 0.0
        else:
            # Fall back to Jaccard if no IC values
            return len(common_terms) / len(set1 | set2)


class PhenopacketSimilarityCalculator:
    """High-level interface for computing phenopacket similarities."""

    def __init__(
        self,
        similarity_metric: HPOSimilarity,
        cache_similarities: bool = True
    ):
        """
        Initialize calculator.

        Args:
            similarity_metric: The similarity metric to use
            cache_similarities: Whether to cache pairwise similarities
        """
        self.similarity_metric = similarity_metric
        self.cache_similarities = cache_similarities
        self._cache = {}

    def compute_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict
    ) -> float:
        """
        Compute similarity between two phenopackets.

        Args:
            phenopacket1: First phenopacket
            phenopacket2: Second phenopacket

        Returns:
            Similarity score
        """
        id1 = phenopacket1.get("id")
        id2 = phenopacket2.get("id")

        # Check cache
        if self.cache_similarities and id1 and id2:
            cache_key = tuple(sorted([id1, id2]))
            if cache_key in self._cache:
                return self._cache[cache_key]

        # Extract terms
        terms1 = self.similarity_metric.load_phenopacket_terms(phenopacket1)
        terms2 = self.similarity_metric.load_phenopacket_terms(phenopacket2)

        # Compute similarity
        similarity = self.similarity_metric.compute_pairwise_similarity(terms1, terms2)

        # Cache result
        if self.cache_similarities and id1 and id2:
            self._cache[cache_key] = similarity

        return similarity

    def compute_similarity_matrix(
        self,
        phenopackets: List[Dict]
    ) -> np.ndarray:
        """
        Compute all-vs-all similarity matrix.

        Args:
            phenopackets: List of phenopackets

        Returns:
            NxN similarity matrix
        """
        n = len(phenopackets)
        matrix = np.zeros((n, n))

        for i in range(n):
            for j in range(i, n):
                if i == j:
                    matrix[i, j] = 1.0
                else:
                    sim = self.compute_similarity(phenopackets[i], phenopackets[j])
                    matrix[i, j] = sim
                    matrix[j, i] = sim  # Symmetric

        return matrix

    def find_most_similar(
        self,
        query_phenopacket: Dict,
        candidate_phenopackets: List[Dict],
        top_k: int = 10
    ) -> List[Tuple[int, float]]:
        """
        Find most similar phenopackets to a query.

        Args:
            query_phenopacket: Query phenopacket
            candidate_phenopackets: Candidate phenopackets to search
            top_k: Number of top results to return

        Returns:
            List of (index, similarity_score) tuples, sorted by similarity
        """
        similarities = []

        for idx, candidate in enumerate(candidate_phenopackets):
            sim = self.compute_similarity(query_phenopacket, candidate)
            similarities.append((idx, sim))

        # Sort by similarity (descending)
        similarities.sort(key=lambda x: x[1], reverse=True)

        return similarities[:top_k]


def load_phenopackets(file_path: str) -> List[Dict]:
    """Load phenopackets from JSON file."""
    with open(file_path, 'r') as f:
        return json.load(f)


def compute_uniform_ic(terms: List[str]) -> Dict[str, float]:
    """
    Compute uniform IC values (all terms equally informative).

    Args:
        terms: List of all HPO terms

    Returns:
        Dictionary of IC values
    """
    return {term: 1.0 for term in set(terms)}


def compute_empirical_ic(phenopackets: List[Dict]) -> Dict[str, float]:
    """
    Compute empirical IC from phenopacket corpus.

    IC(term) = -log(P(term)) where P(term) is frequency in corpus.

    Args:
        phenopackets: List of phenopackets

    Returns:
        Dictionary of IC values
    """
    term_counts = defaultdict(int)
    total_patients = len(phenopackets)

    # Count term occurrences
    for pp in phenopackets:
        terms = set()
        for feature in pp.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                terms.add(feature["type"]["id"])

        for term in terms:
            term_counts[term] += 1

    # Compute IC
    ic_values = {}
    for term, count in term_counts.items():
        freq = count / total_patients
        ic_values[term] = -np.log(freq) if freq > 0 else 0.0

    return ic_values


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Load test data
    phenopackets = load_phenopackets("data/synthetic/test_marfan_10.json")

    # Compute empirical IC
    ic_values = compute_empirical_ic(phenopackets)
    print(f"Computed IC for {len(ic_values)} terms")

    # Test different similarity metrics
    metrics = {
        "Jaccard": JaccardSimilarity(),
        "Cosine": CosineSimilarity(ic_values),
        "Simplified Resnik": SimplifiedResnikSimilarity(ic_values)
    }

    for name, metric in metrics.items():
        calc = PhenopacketSimilarityCalculator(metric)

        # Compute similarity between first two phenopackets
        sim = calc.compute_similarity(phenopackets[0], phenopackets[1])
        print(f"{name} similarity: {sim:.4f}")

        # Find similar patients to first one
        matches = calc.find_most_similar(
            phenopackets[0],
            phenopackets[1:],
            top_k=3
        )
        print(f"{name} top matches: {[(idx, f'{score:.4f}') for idx, score in matches]}")
