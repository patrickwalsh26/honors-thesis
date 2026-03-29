"""
Evaluation metrics for phenotype matching systems.

Implements standard information retrieval metrics:
- Precision@k, Recall@k, F1@k
- nDCG@k (Normalized Discounted Cumulative Gain)
- Mean Average Precision (MAP)
- Mean Reciprocal Rank (MRR)

Also provides a high-level PhenopacketEvaluator for evaluating
matching systems on phenopacket datasets.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


class RetrievalMetrics:
    """
    Standard information retrieval metrics.

    All metrics assume:
    - retrieved: List of item indices in ranked order
    - relevant: Set of indices that are relevant to the query
    """

    @staticmethod
    def precision_at_k(
        retrieved: List[int],
        relevant: Set[int],
        k: int
    ) -> float:
        """
        Precision@k: Fraction of top-k that are relevant.

        P@k = |retrieved_k ∩ relevant| / k

        Args:
            retrieved: Ranked list of retrieved item indices
            relevant: Set of relevant item indices
            k: Cutoff rank

        Returns:
            Precision@k score in [0, 1]
        """
        if k == 0:
            return 0.0

        retrieved_k = set(retrieved[:k])
        hits = len(retrieved_k & relevant)

        return hits / k

    @staticmethod
    def recall_at_k(
        retrieved: List[int],
        relevant: Set[int],
        k: int
    ) -> float:
        """
        Recall@k: Fraction of relevant items in top-k.

        R@k = |retrieved_k ∩ relevant| / |relevant|

        Args:
            retrieved: Ranked list of retrieved item indices
            relevant: Set of relevant item indices
            k: Cutoff rank

        Returns:
            Recall@k score in [0, 1]
        """
        if not relevant:
            return 0.0

        retrieved_k = set(retrieved[:k])
        hits = len(retrieved_k & relevant)

        return hits / len(relevant)

    @staticmethod
    def f1_at_k(
        retrieved: List[int],
        relevant: Set[int],
        k: int
    ) -> float:
        """
        F1@k: Harmonic mean of Precision@k and Recall@k.

        F1@k = 2 * P@k * R@k / (P@k + R@k)

        Args:
            retrieved: Ranked list of retrieved item indices
            relevant: Set of relevant item indices
            k: Cutoff rank

        Returns:
            F1@k score in [0, 1]
        """
        p = RetrievalMetrics.precision_at_k(retrieved, relevant, k)
        r = RetrievalMetrics.recall_at_k(retrieved, relevant, k)

        if p + r == 0:
            return 0.0

        return 2 * p * r / (p + r)

    @staticmethod
    def dcg_at_k(
        retrieved: List[int],
        relevance_scores: Dict[int, float],
        k: int
    ) -> float:
        """
        Discounted Cumulative Gain at k.

        DCG@k = Σ(i=1 to k) rel_i / log2(i + 1)

        Args:
            retrieved: Ranked list of retrieved item indices
            relevance_scores: Dictionary mapping indices to relevance scores
            k: Cutoff rank

        Returns:
            DCG@k score
        """
        dcg = 0.0

        for i, idx in enumerate(retrieved[:k]):
            rel = relevance_scores.get(idx, 0.0)
            # Position is 1-indexed, so use i+1
            # Add 1 more to avoid log2(1)=0 division issue
            dcg += rel / np.log2(i + 2)

        return dcg

    @staticmethod
    def ndcg_at_k(
        retrieved: List[int],
        relevance_scores: Dict[int, float],
        k: int
    ) -> float:
        """
        Normalized Discounted Cumulative Gain at k.

        nDCG@k = DCG@k / IDCG@k

        where IDCG@k is DCG@k with ideal ranking.

        Args:
            retrieved: Ranked list of retrieved item indices
            relevance_scores: Dictionary mapping indices to relevance scores
            k: Cutoff rank

        Returns:
            nDCG@k score in [0, 1]
        """
        dcg = RetrievalMetrics.dcg_at_k(retrieved, relevance_scores, k)

        # Compute ideal DCG (items sorted by relevance)
        ideal_order = sorted(
            relevance_scores.keys(),
            key=lambda x: relevance_scores[x],
            reverse=True
        )
        idcg = RetrievalMetrics.dcg_at_k(ideal_order, relevance_scores, k)

        if idcg == 0:
            return 0.0

        return dcg / idcg

    @staticmethod
    def average_precision(
        retrieved: List[int],
        relevant: Set[int]
    ) -> float:
        """
        Average Precision for a single query.

        AP = (1/|relevant|) * Σ P@k * rel(k)

        where rel(k) is 1 if item at rank k is relevant, 0 otherwise.

        Args:
            retrieved: Ranked list of retrieved item indices
            relevant: Set of relevant item indices

        Returns:
            Average Precision score in [0, 1]
        """
        if not relevant:
            return 0.0

        hits = 0
        sum_precision = 0.0

        for i, idx in enumerate(retrieved):
            if idx in relevant:
                hits += 1
                precision_at_i = hits / (i + 1)
                sum_precision += precision_at_i

        return sum_precision / len(relevant)

    @staticmethod
    def mean_average_precision(
        all_retrieved: List[List[int]],
        all_relevant: List[Set[int]]
    ) -> float:
        """
        Mean Average Precision across multiple queries.

        MAP = (1/|queries|) * Σ AP(query)

        Args:
            all_retrieved: List of retrieved lists, one per query
            all_relevant: List of relevant sets, one per query

        Returns:
            MAP score in [0, 1]
        """
        if not all_retrieved:
            return 0.0

        aps = []
        for retrieved, relevant in zip(all_retrieved, all_relevant):
            if relevant:  # Only count queries with relevant items
                ap = RetrievalMetrics.average_precision(retrieved, relevant)
                aps.append(ap)

        if not aps:
            return 0.0

        return np.mean(aps)

    @staticmethod
    def reciprocal_rank(
        retrieved: List[int],
        relevant: Set[int]
    ) -> float:
        """
        Reciprocal Rank: 1 / (rank of first relevant item).

        Args:
            retrieved: Ranked list of retrieved item indices
            relevant: Set of relevant item indices

        Returns:
            RR score in (0, 1], or 0 if no relevant items found
        """
        for i, idx in enumerate(retrieved):
            if idx in relevant:
                return 1.0 / (i + 1)

        return 0.0

    @staticmethod
    def mean_reciprocal_rank(
        all_retrieved: List[List[int]],
        all_relevant: List[Set[int]]
    ) -> float:
        """
        Mean Reciprocal Rank across multiple queries.

        MRR = (1/|queries|) * Σ RR(query)

        Args:
            all_retrieved: List of retrieved lists, one per query
            all_relevant: List of relevant sets, one per query

        Returns:
            MRR score in [0, 1]
        """
        if not all_retrieved:
            return 0.0

        rrs = []
        for retrieved, relevant in zip(all_retrieved, all_relevant):
            rr = RetrievalMetrics.reciprocal_rank(retrieved, relevant)
            rrs.append(rr)

        return np.mean(rrs)

    @staticmethod
    def hit_rate_at_k(
        retrieved: List[int],
        relevant: Set[int],
        k: int
    ) -> float:
        """
        Hit Rate@k: Whether any relevant item appears in top-k.

        Args:
            retrieved: Ranked list of retrieved item indices
            relevant: Set of relevant item indices
            k: Cutoff rank

        Returns:
            1.0 if hit, 0.0 otherwise
        """
        retrieved_k = set(retrieved[:k])
        return 1.0 if retrieved_k & relevant else 0.0


class PhenopacketEvaluator:
    """
    High-level evaluator for phenotype matching systems.

    Provides convenient methods for evaluating matching systems
    on phenopacket datasets using standard retrieval metrics.
    """

    def __init__(
        self,
        calculator,
        k_values: List[int] = [1, 5, 10, 20],
        metrics: Optional[List[str]] = None
    ):
        """
        Initialize evaluator.

        Args:
            calculator: Any calculator with find_most_similar method
            k_values: List of k values for @k metrics
            metrics: List of metric names to compute
        """
        self.calculator = calculator
        self.k_values = k_values
        self.metrics = metrics or ["precision", "recall", "f1", "ndcg"]
        self._retrieval_metrics = RetrievalMetrics()

    def evaluate_retrieval(
        self,
        queries: List[Dict],
        candidates: List[Dict],
        ground_truth: Dict[str, Set[int]],
        verbose: bool = False
    ) -> Dict[str, Dict[int, float]]:
        """
        Evaluate retrieval performance.

        Args:
            queries: List of query phenopackets
            candidates: List of candidate phenopackets
            ground_truth: Dict mapping query_id to set of relevant candidate indices
            verbose: Print progress

        Returns:
            Dict of metric_name -> {k: mean_score}
        """
        # Initialize results storage
        results = {
            metric: {k: [] for k in self.k_values}
            for metric in self.metrics
        }

        max_k = max(self.k_values)

        for i, query in enumerate(queries):
            query_id = query.get("id")

            if verbose and i % 10 == 0:
                logger.info(f"Evaluating query {i+1}/{len(queries)}")

            relevant = ground_truth.get(query_id, set())

            # Skip queries with no relevant items
            if not relevant:
                continue

            # Get retrieval results
            matches = self.calculator.find_most_similar(
                query, candidates, top_k=max_k
            )

            # Handle suppressed results (k-anonymity)
            if matches is None:
                # All metrics are 0 for suppressed queries
                for metric in self.metrics:
                    for k in self.k_values:
                        results[metric][k].append(0.0)
                continue

            retrieved = [idx for idx, _ in matches]

            # Compute metrics at each k
            for k in self.k_values:
                if "precision" in self.metrics:
                    results["precision"][k].append(
                        self._retrieval_metrics.precision_at_k(retrieved, relevant, k)
                    )

                if "recall" in self.metrics:
                    results["recall"][k].append(
                        self._retrieval_metrics.recall_at_k(retrieved, relevant, k)
                    )

                if "f1" in self.metrics:
                    results["f1"][k].append(
                        self._retrieval_metrics.f1_at_k(retrieved, relevant, k)
                    )

                if "ndcg" in self.metrics:
                    # Use binary relevance for nDCG
                    rel_scores = {idx: 1.0 for idx in relevant}
                    results["ndcg"][k].append(
                        self._retrieval_metrics.ndcg_at_k(retrieved, rel_scores, k)
                    )

                if "hit_rate" in self.metrics:
                    results["hit_rate"][k].append(
                        self._retrieval_metrics.hit_rate_at_k(retrieved, relevant, k)
                    )

        # Compute means
        return {
            metric: {k: np.mean(scores) if scores else 0.0 for k, scores in k_results.items()}
            for metric, k_results in results.items()
        }

    def cross_disease_evaluation(
        self,
        phenopackets: List[Dict],
        exclude_self: bool = True,
        verbose: bool = False
    ) -> Dict[str, Dict[int, float]]:
        """
        Evaluate ability to retrieve same-disease patients.

        Ground truth: Same disease = relevant.

        Args:
            phenopackets: List of phenopackets with disease annotations
            exclude_self: Exclude self from results
            verbose: Print progress

        Returns:
            Dict of metric_name -> {k: mean_score}
        """
        # Group by disease
        by_disease = defaultdict(list)
        for i, pp in enumerate(phenopackets):
            diseases = pp.get("diseases", [])
            if diseases:
                disease = diseases[0]["term"].get("label", "Unknown")
            else:
                disease = "Unknown"
            by_disease[disease].append(i)

        # Build ground truth
        ground_truth = {}
        for i, pp in enumerate(phenopackets):
            query_id = pp.get("id")
            diseases = pp.get("diseases", [])
            if diseases:
                disease = diseases[0]["term"].get("label", "Unknown")
            else:
                disease = "Unknown"

            # Same disease patients are relevant
            relevant = set(by_disease[disease])
            if exclude_self:
                relevant = relevant - {i}

            ground_truth[query_id] = relevant

        return self.evaluate_retrieval(
            phenopackets, phenopackets, ground_truth, verbose
        )

    def evaluate_with_details(
        self,
        queries: List[Dict],
        candidates: List[Dict],
        ground_truth: Dict[str, Set[int]]
    ) -> pd.DataFrame:
        """
        Evaluate with per-query details.

        Args:
            queries: List of query phenopackets
            candidates: List of candidate phenopackets
            ground_truth: Dict mapping query_id to relevant indices

        Returns:
            DataFrame with per-query metrics
        """
        records = []
        max_k = max(self.k_values)

        for query in queries:
            query_id = query.get("id")
            relevant = ground_truth.get(query_id, set())

            if not relevant:
                continue

            matches = self.calculator.find_most_similar(
                query, candidates, top_k=max_k
            )

            if matches is None:
                # Suppressed query
                record = {
                    "query_id": query_id,
                    "num_relevant": len(relevant),
                    "suppressed": True
                }
                for k in self.k_values:
                    for metric in self.metrics:
                        record[f"{metric}@{k}"] = 0.0
                records.append(record)
                continue

            retrieved = [idx for idx, _ in matches]

            record = {
                "query_id": query_id,
                "num_relevant": len(relevant),
                "suppressed": False
            }

            for k in self.k_values:
                record[f"precision@{k}"] = self._retrieval_metrics.precision_at_k(
                    retrieved, relevant, k
                )
                record[f"recall@{k}"] = self._retrieval_metrics.recall_at_k(
                    retrieved, relevant, k
                )
                record[f"f1@{k}"] = self._retrieval_metrics.f1_at_k(
                    retrieved, relevant, k
                )
                rel_scores = {idx: 1.0 for idx in relevant}
                record[f"ndcg@{k}"] = self._retrieval_metrics.ndcg_at_k(
                    retrieved, rel_scores, k
                )

            records.append(record)

        return pd.DataFrame(records)

    def disease_stratified_evaluation(
        self,
        phenopackets: List[Dict]
    ) -> pd.DataFrame:
        """
        Evaluate performance stratified by disease.

        Args:
            phenopackets: List of phenopackets with disease annotations

        Returns:
            DataFrame with per-disease metrics
        """
        # Group by disease
        by_disease = defaultdict(list)
        for i, pp in enumerate(phenopackets):
            diseases = pp.get("diseases", [])
            if diseases:
                disease = diseases[0]["term"].get("label", "Unknown")
            else:
                disease = "Unknown"
            by_disease[disease].append(i)

        records = []

        for disease, indices in by_disease.items():
            if len(indices) < 2:
                continue

            # Get phenopackets for this disease
            disease_pps = [phenopackets[i] for i in indices]

            # Build ground truth (all same-disease patients)
            ground_truth = {
                pp.get("id"): set(indices) - {indices[j]}
                for j, pp in enumerate(disease_pps)
            }

            # Evaluate
            results = self.evaluate_retrieval(
                disease_pps, phenopackets, ground_truth
            )

            record = {"disease": disease, "n_patients": len(indices)}
            for metric, k_scores in results.items():
                for k, score in k_scores.items():
                    record[f"{metric}@{k}"] = score

            records.append(record)

        return pd.DataFrame(records)


def compute_metric_confidence_interval(
    scores: List[float],
    confidence: float = 0.95
) -> Tuple[float, float, float]:
    """
    Compute confidence interval for a metric.

    Uses bootstrap percentile method.

    Args:
        scores: List of scores
        confidence: Confidence level (default 95%)

    Returns:
        Tuple of (mean, lower_bound, upper_bound)
    """
    if not scores:
        return 0.0, 0.0, 0.0

    mean = np.mean(scores)

    # Bootstrap
    n_bootstrap = 1000
    bootstrap_means = []

    for _ in range(n_bootstrap):
        sample = np.random.choice(scores, size=len(scores), replace=True)
        bootstrap_means.append(np.mean(sample))

    alpha = 1 - confidence
    lower = np.percentile(bootstrap_means, 100 * alpha / 2)
    upper = np.percentile(bootstrap_means, 100 * (1 - alpha / 2))

    return mean, lower, upper


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    # Test individual metrics
    retrieved = [1, 3, 5, 7, 9, 2, 4, 6, 8, 10]
    relevant = {1, 2, 3, 4, 5}

    print("Testing RetrievalMetrics:")
    print(f"  Precision@5: {RetrievalMetrics.precision_at_k(retrieved, relevant, 5):.4f}")
    print(f"  Recall@5: {RetrievalMetrics.recall_at_k(retrieved, relevant, 5):.4f}")
    print(f"  F1@5: {RetrievalMetrics.f1_at_k(retrieved, relevant, 5):.4f}")

    rel_scores = {i: 1.0 for i in relevant}
    print(f"  nDCG@5: {RetrievalMetrics.ndcg_at_k(retrieved, rel_scores, 5):.4f}")
    print(f"  AP: {RetrievalMetrics.average_precision(retrieved, relevant):.4f}")
    print(f"  RR: {RetrievalMetrics.reciprocal_rank(retrieved, relevant):.4f}")

    # Test confidence interval
    scores = [0.8, 0.75, 0.82, 0.79, 0.81, 0.77, 0.83, 0.76, 0.80, 0.78]
    mean, lower, upper = compute_metric_confidence_interval(scores)
    print(f"\nConfidence interval: {mean:.4f} [{lower:.4f}, {upper:.4f}]")
