"""
Locality-Sensitive Hashing for approximate nearest neighbor search.

Implements:
- MinHashLSH: Efficient Jaccard similarity approximation
- LSHPhenopacketCalculator: Two-stage retrieval (LSH + exact re-ranking)

Enables scalable phenotype matching for large databases
while maintaining compatibility with privacy mechanisms.
"""

import numpy as np
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
import logging
import hashlib

logger = logging.getLogger(__name__)


class MinHash:
    """
    MinHash signature generator for set similarity.

    MinHash enables estimation of Jaccard similarity:
    J(A, B) ≈ (# matching MinHash values) / num_perm
    """

    def __init__(
        self,
        num_perm: int = 128,
        seed: int = 42
    ):
        """
        Initialize MinHash generator.

        Args:
            num_perm: Number of permutation hash functions
            seed: Random seed for reproducibility
        """
        self.num_perm = num_perm
        self.seed = seed

        # Generate hash function parameters
        np.random.seed(seed)
        self._max_hash = (1 << 32) - 1
        self._mersenne_prime = (1 << 61) - 1

        # a*x + b mod prime parameters
        self._a = np.random.randint(1, self._mersenne_prime, size=num_perm, dtype=np.uint64)
        self._b = np.random.randint(0, self._mersenne_prime, size=num_perm, dtype=np.uint64)

    def _hash_element(self, element: str) -> int:
        """Hash an element to integer."""
        return int(hashlib.sha256(element.encode()).hexdigest()[:16], 16)

    def compute_signature(self, elements: Set[str]) -> np.ndarray:
        """
        Compute MinHash signature for a set.

        Args:
            elements: Set of elements (e.g., HPO term IDs)

        Returns:
            MinHash signature array of shape (num_perm,)
        """
        if not elements:
            return np.full(self.num_perm, self._max_hash, dtype=np.uint64)

        # Initialize signature with max values
        signature = np.full(self.num_perm, self._max_hash, dtype=np.uint64)

        for element in elements:
            h = self._hash_element(element)
            # Apply each permutation
            hashes = (self._a * h + self._b) % self._mersenne_prime
            signature = np.minimum(signature, hashes)

        return signature

    def estimate_jaccard(
        self,
        sig1: np.ndarray,
        sig2: np.ndarray
    ) -> float:
        """
        Estimate Jaccard similarity from MinHash signatures.

        Args:
            sig1: First MinHash signature
            sig2: Second MinHash signature

        Returns:
            Estimated Jaccard similarity
        """
        return np.mean(sig1 == sig2)


class MinHashLSH:
    """
    MinHash LSH for approximate Jaccard similarity search.

    Uses banding technique to find candidate pairs:
    - Divide signature into b bands of r rows each
    - Two sets are candidates if any band matches exactly
    - Probability of being candidate ≈ 1 - (1 - s^r)^b for Jaccard s
    """

    def __init__(
        self,
        num_perm: int = 128,
        threshold: float = 0.5,
        num_bands: Optional[int] = None,
        seed: int = 42
    ):
        """
        Initialize MinHash LSH index.

        Args:
            num_perm: Number of permutations
            threshold: Similarity threshold for candidate generation
            num_bands: Number of bands (auto-computed if None)
            seed: Random seed
        """
        self.num_perm = num_perm
        self.threshold = threshold
        self.seed = seed

        # Compute optimal banding
        if num_bands is None:
            self.num_bands, self.rows_per_band = self._optimal_banding(threshold)
        else:
            self.num_bands = num_bands
            self.rows_per_band = num_perm // num_bands

        # Adjust if num_perm not evenly divisible
        if self.num_bands * self.rows_per_band > num_perm:
            self.num_bands = num_perm // self.rows_per_band

        # MinHash generator
        self._minhash = MinHash(num_perm, seed)

        # LSH index: band_idx -> hash -> set of indices
        self._index: Dict[int, Dict[int, Set[int]]] = {
            b: defaultdict(set) for b in range(self.num_bands)
        }

        # Store signatures
        self._signatures: Dict[int, np.ndarray] = {}

        # Store original data references
        self._phenopackets: List[Dict] = []

    def _optimal_banding(self, threshold: float) -> Tuple[int, int]:
        """
        Find optimal number of bands for target threshold.

        Finds (b, r) where b*r ≤ num_perm and collision probability
        at threshold is maximized.

        Args:
            threshold: Target Jaccard similarity threshold

        Returns:
            Tuple of (num_bands, rows_per_band)
        """
        best_b, best_r = 1, self.num_perm

        for b in range(1, self.num_perm + 1):
            if self.num_perm % b != 0:
                continue

            r = self.num_perm // b

            # Probability of collision at threshold
            # P(candidate | Jaccard=s) = 1 - (1 - s^r)^b
            prob_at_threshold = 1 - (1 - threshold**r)**b

            # We want high probability at threshold
            if prob_at_threshold >= 0.7:
                best_b, best_r = b, r
                break

        logger.debug(f"LSH banding: {best_b} bands x {best_r} rows (threshold={threshold})")
        return best_b, best_r

    def _hash_band(self, signature: np.ndarray, band: int) -> int:
        """Hash a band of the signature."""
        start = band * self.rows_per_band
        end = start + self.rows_per_band
        band_values = signature[start:end]
        return hash(tuple(band_values))

    def _extract_terms(self, phenopacket: Dict) -> Set[str]:
        """Extract HPO terms from phenopacket."""
        terms = set()
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                terms.add(feature["type"]["id"])
        return terms

    def index_phenopacket(self, idx: int, phenopacket: Dict):
        """
        Add a phenopacket to the LSH index.

        Args:
            idx: Index/ID for the phenopacket
            phenopacket: Phenopacket dictionary
        """
        terms = self._extract_terms(phenopacket)
        signature = self._minhash.compute_signature(terms)
        self._signatures[idx] = signature

        # Add to each band's hash table
        for band in range(self.num_bands):
            band_hash = self._hash_band(signature, band)
            self._index[band][band_hash].add(idx)

    def index_all(self, phenopackets: List[Dict]):
        """
        Index all phenopackets.

        Args:
            phenopackets: List of phenopackets to index
        """
        self._phenopackets = phenopackets

        for idx, pp in enumerate(phenopackets):
            self.index_phenopacket(idx, pp)

        logger.info(f"Indexed {len(phenopackets)} phenopackets in LSH")

    def query_candidates(self, phenopacket: Dict) -> Set[int]:
        """
        Find candidate indices that might be similar.

        Args:
            phenopacket: Query phenopacket

        Returns:
            Set of candidate indices
        """
        terms = self._extract_terms(phenopacket)
        query_sig = self._minhash.compute_signature(terms)

        candidates = set()
        for band in range(self.num_bands):
            band_hash = self._hash_band(query_sig, band)
            if band_hash in self._index[band]:
                candidates.update(self._index[band][band_hash])

        return candidates

    def query(
        self,
        phenopacket: Dict,
        top_k: int = 10
    ) -> List[Tuple[int, float]]:
        """
        Find approximate nearest neighbors.

        Args:
            phenopacket: Query phenopacket
            top_k: Number of results to return

        Returns:
            List of (index, estimated_jaccard) tuples
        """
        terms = self._extract_terms(phenopacket)
        query_sig = self._minhash.compute_signature(terms)

        # Get candidates
        candidates = self.query_candidates(phenopacket)

        if not candidates:
            return []

        # Estimate Jaccard for each candidate
        results = []
        for idx in candidates:
            cand_sig = self._signatures[idx]
            jaccard_est = self._minhash.estimate_jaccard(query_sig, cand_sig)
            results.append((idx, jaccard_est))

        # Sort by estimated similarity
        results.sort(key=lambda x: x[1], reverse=True)

        return results[:top_k]

    def get_statistics(self) -> Dict[str, Any]:
        """Get index statistics."""
        band_sizes = []
        for band_idx in range(self.num_bands):
            band_sizes.extend(len(indices) for indices in self._index[band_idx].values())

        return {
            "num_indexed": len(self._signatures),
            "num_bands": self.num_bands,
            "rows_per_band": self.rows_per_band,
            "threshold": self.threshold,
            "avg_bucket_size": np.mean(band_sizes) if band_sizes else 0,
            "max_bucket_size": max(band_sizes) if band_sizes else 0
        }


class LSHPhenopacketCalculator:
    """
    LSH-accelerated phenopacket similarity calculator.

    Two-stage retrieval:
    1. LSH for fast candidate generation
    2. Exact similarity for re-ranking

    Follows the same interface as PhenopacketSimilarityCalculator.
    """

    def __init__(
        self,
        exact_calculator=None,
        num_perm: int = 128,
        threshold: float = 0.3,
        refine: bool = True
    ):
        """
        Initialize LSH calculator.

        Args:
            exact_calculator: Calculator for exact re-ranking
            num_perm: Number of MinHash permutations
            threshold: LSH candidate threshold
            refine: Whether to re-rank with exact similarity
        """
        self.exact_calculator = exact_calculator
        self.refine = refine

        self.lsh = MinHashLSH(
            num_perm=num_perm,
            threshold=threshold
        )

        self._indexed = False
        self._phenopackets: List[Dict] = []

    def index(self, phenopackets: List[Dict]):
        """
        Build LSH index for phenopackets.

        Args:
            phenopackets: List of phenopackets to index
        """
        self._phenopackets = phenopackets
        self.lsh.index_all(phenopackets)
        self._indexed = True

    def compute_similarity(
        self,
        phenopacket1: Dict,
        phenopacket2: Dict
    ) -> float:
        """
        Compute similarity between two phenopackets.

        Uses exact calculator if available, otherwise MinHash estimate.

        Args:
            phenopacket1: First phenopacket
            phenopacket2: Second phenopacket

        Returns:
            Similarity score
        """
        if self.exact_calculator:
            return self.exact_calculator.compute_similarity(phenopacket1, phenopacket2)

        # Fallback to MinHash estimate
        terms1 = self._extract_terms(phenopacket1)
        terms2 = self._extract_terms(phenopacket2)

        sig1 = self.lsh._minhash.compute_signature(terms1)
        sig2 = self.lsh._minhash.compute_signature(terms2)

        return self.lsh._minhash.estimate_jaccard(sig1, sig2)

    def _extract_terms(self, phenopacket: Dict) -> Set[str]:
        """Extract HPO terms from phenopacket."""
        terms = set()
        for feature in phenopacket.get("phenotypicFeatures", []):
            if not feature.get("excluded", False):
                terms.add(feature["type"]["id"])
        return terms

    def find_most_similar(
        self,
        query_phenopacket: Dict,
        candidate_phenopackets: Optional[List[Dict]] = None,
        top_k: int = 10
    ) -> List[Tuple[int, float]]:
        """
        Find most similar phenopackets using LSH.

        Args:
            query_phenopacket: Query phenopacket
            candidate_phenopackets: Candidates (uses indexed if None)
            top_k: Number of results

        Returns:
            List of (index, similarity) tuples
        """
        # Re-index if candidates changed
        if candidate_phenopackets is not None:
            if candidate_phenopackets != self._phenopackets:
                self.index(candidate_phenopackets)

        if not self._indexed:
            raise ValueError("No phenopackets indexed. Call index() first.")

        # Stage 1: LSH candidate generation
        # Get more candidates than needed for re-ranking
        lsh_results = self.lsh.query(query_phenopacket, top_k=top_k * 3)

        if not lsh_results:
            return []

        if not self.refine or not self.exact_calculator:
            # Return LSH estimates directly
            return lsh_results[:top_k]

        # Stage 2: Exact similarity re-ranking
        refined = []
        for idx, _ in lsh_results:
            exact_sim = self.exact_calculator.compute_similarity(
                query_phenopacket,
                self._phenopackets[idx]
            )
            refined.append((idx, exact_sim))

        # Sort by exact similarity
        refined.sort(key=lambda x: x[1], reverse=True)

        return refined[:top_k]

    def find_above_threshold(
        self,
        query_phenopacket: Dict,
        threshold: float
    ) -> List[Tuple[int, float]]:
        """
        Find all phenopackets above similarity threshold.

        Uses LSH for efficient filtering.

        Args:
            query_phenopacket: Query phenopacket
            threshold: Minimum similarity threshold

        Returns:
            List of (index, similarity) tuples above threshold
        """
        if not self._indexed:
            raise ValueError("No phenopackets indexed. Call index() first.")

        # Get LSH candidates (they're already likely above threshold)
        candidates = self.lsh.query_candidates(query_phenopacket)

        results = []
        for idx in candidates:
            if self.exact_calculator:
                sim = self.exact_calculator.compute_similarity(
                    query_phenopacket,
                    self._phenopackets[idx]
                )
            else:
                terms1 = self._extract_terms(query_phenopacket)
                terms2 = self._extract_terms(self._phenopackets[idx])
                sig1 = self.lsh._minhash.compute_signature(terms1)
                sig2 = self.lsh._minhash.compute_signature(terms2)
                sim = self.lsh._minhash.estimate_jaccard(sig1, sig2)

            if sim >= threshold:
                results.append((idx, sim))

        results.sort(key=lambda x: x[1], reverse=True)
        return results

    def get_statistics(self) -> Dict[str, Any]:
        """Get LSH index statistics."""
        stats = self.lsh.get_statistics()
        stats["has_exact_calculator"] = self.exact_calculator is not None
        stats["refine_enabled"] = self.refine
        return stats


def create_lsh_calculator(
    exact_calculator=None,
    num_perm: int = 128,
    threshold: float = 0.3
) -> LSHPhenopacketCalculator:
    """
    Factory function to create LSH calculator.

    Args:
        exact_calculator: Calculator for re-ranking
        num_perm: Number of MinHash permutations
        threshold: LSH threshold

    Returns:
        Configured LSHPhenopacketCalculator
    """
    return LSHPhenopacketCalculator(
        exact_calculator=exact_calculator,
        num_perm=num_perm,
        threshold=threshold
    )


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Create sample phenopackets
    sample_phenopackets = [
        {
            "id": f"patient_{i}",
            "phenotypicFeatures": [
                {"type": {"id": f"HP:000000{j}"}, "excluded": False}
                for j in range(1, 5 + i % 3)
            ]
        }
        for i in range(100)
    ]

    # Create LSH calculator
    calc = LSHPhenopacketCalculator(
        num_perm=64,
        threshold=0.3,
        refine=False  # Use LSH estimates only
    )

    # Index phenopackets
    calc.index(sample_phenopackets)

    # Query
    query = sample_phenopackets[0]
    results = calc.find_most_similar(query, top_k=5)

    print("Top 5 similar patients (LSH):")
    for idx, sim in results:
        print(f"  Patient {idx}: {sim:.4f}")

    # Statistics
    stats = calc.get_statistics()
    print(f"\nLSH Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
