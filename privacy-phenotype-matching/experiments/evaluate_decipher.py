#!/usr/bin/env python3
"""
Evaluate privacy-preserving phenotype matching on DECIPHER data.

This experiment script:
1. Loads DECIPHER patient data (real or simulated)
2. Computes baseline similarity metrics
3. Applies privacy mechanisms (PSI, DP, k-anonymity)
4. Evaluates retrieval performance with privacy-utility tradeoffs
5. Generates figures and statistics for the thesis

Usage:
    python experiments/evaluate_decipher.py --help
    python experiments/evaluate_decipher.py --n-patients 500 --simulated
    python experiments/evaluate_decipher.py --data-file data/decipher/patients.json

Output:
    - results/decipher/baseline_metrics.csv
    - results/decipher/privacy_utility_tradeoff.csv
    - figures/decipher_*.png

References:
    Foreman J, et al. (2022) DECIPHER: Supporting the interpretation and
    sharing of rare disease phenotype-linked variant data. Human Mutation.
"""

import sys
import argparse
import logging
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Any

import numpy as np
import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_integration.decipher_loader import (
    DECIPHERLoader,
    DECIPHERSimulator,
    create_decipher_evaluation_dataset,
)
from src.similarity.hpo_similarity import (
    JaccardSimilarity,
    CosineSimilarity,
    SimplifiedResnikSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
    load_phenopackets,
)
from src.evaluation.metrics import (
    RetrievalMetrics,
    PhenopacketEvaluator,
    compute_metric_confidence_interval,
)

logger = logging.getLogger(__name__)


class DECIPHEREvaluator:
    """
    Comprehensive evaluator for DECIPHER phenotype matching experiments.

    This class provides methods for:
    - Loading and preprocessing DECIPHER data
    - Computing baseline similarity performance
    - Evaluating privacy-preserving mechanisms
    - Generating thesis figures and tables
    """

    def __init__(
        self,
        data_dir: str = "data/decipher",
        results_dir: str = "results/decipher",
        figures_dir: str = "figures"
    ):
        """
        Initialize evaluator.

        Args:
            data_dir: Directory containing DECIPHER data
            results_dir: Directory for results output
            figures_dir: Directory for figure output
        """
        self.data_dir = Path(data_dir)
        self.results_dir = Path(results_dir)
        self.figures_dir = Path(figures_dir)

        # Create directories
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)

        # Data storage
        self.phenopackets: List[Dict] = []
        self.ic_values: Dict[str, float] = {}

        # Evaluation parameters
        self.k_values = [1, 5, 10, 20, 50]
        self.metrics = ["precision", "recall", "ndcg", "hit_rate"]

    def load_data(
        self,
        data_file: Optional[str] = None,
        n_patients: int = 500,
        use_simulated: bool = False,
        min_phenotypes: int = 3
    ) -> int:
        """
        Load DECIPHER data for evaluation.

        Args:
            data_file: Path to existing data file
            n_patients: Number of patients (if generating)
            use_simulated: Use simulated data
            min_phenotypes: Minimum phenotypes per patient

        Returns:
            Number of loaded patients
        """
        if data_file and Path(data_file).exists():
            logger.info(f"Loading data from {data_file}")
            self.phenopackets = load_phenopackets(data_file)
        else:
            logger.info("Creating DECIPHER evaluation dataset...")
            output_path = create_decipher_evaluation_dataset(
                output_dir=str(self.data_dir),
                n_patients=n_patients,
                use_simulated=use_simulated
            )
            self.phenopackets = load_phenopackets(output_path)

        # Filter by minimum phenotypes
        self.phenopackets = [
            pp for pp in self.phenopackets
            if len(pp.get("phenotypicFeatures", [])) >= min_phenotypes
        ]

        # Compute IC values
        logger.info("Computing information content values...")
        self.ic_values = compute_empirical_ic(self.phenopackets)

        logger.info(f"Loaded {len(self.phenopackets)} patients with {len(self.ic_values)} unique HPO terms")
        return len(self.phenopackets)

    def get_dataset_statistics(self) -> Dict[str, Any]:
        """Compute comprehensive dataset statistics."""
        if not self.phenopackets:
            return {}

        # Phenotype statistics
        phenotypes_per_patient = []
        all_phenotypes = []
        phenotype_counts = defaultdict(int)

        for pp in self.phenopackets:
            features = pp.get("phenotypicFeatures", [])
            n_features = len(features)
            phenotypes_per_patient.append(n_features)

            for f in features:
                hpo_id = f["type"]["id"]
                all_phenotypes.append(hpo_id)
                phenotype_counts[hpo_id] += 1

        # Sex distribution
        sex_dist = defaultdict(int)
        for pp in self.phenopackets:
            sex = pp.get("subject", {}).get("sex", "UNKNOWN_SEX")
            sex_dist[sex] += 1

        # IC distribution
        ic_sorted = sorted(self.ic_values.values())

        stats = {
            "n_patients": len(self.phenopackets),
            "n_unique_phenotypes": len(set(all_phenotypes)),
            "n_total_annotations": len(all_phenotypes),
            "phenotypes_per_patient": {
                "mean": np.mean(phenotypes_per_patient),
                "std": np.std(phenotypes_per_patient),
                "median": np.median(phenotypes_per_patient),
                "min": min(phenotypes_per_patient),
                "max": max(phenotypes_per_patient),
                "q25": np.percentile(phenotypes_per_patient, 25),
                "q75": np.percentile(phenotypes_per_patient, 75),
            },
            "sex_distribution": dict(sex_dist),
            "rare_phenotypes": sum(1 for c in phenotype_counts.values() if c == 1),
            "common_phenotypes": sum(1 for c in phenotype_counts.values() if c >= 10),
            "ic_stats": {
                "mean": np.mean(ic_sorted),
                "std": np.std(ic_sorted),
                "min": min(ic_sorted),
                "max": max(ic_sorted),
            },
            "most_common_phenotypes": sorted(
                phenotype_counts.items(), key=lambda x: x[1], reverse=True
            )[:10],
        }

        return stats

    def evaluate_baseline_similarities(self) -> pd.DataFrame:
        """
        Evaluate baseline similarity metrics without privacy.

        Returns:
            DataFrame with evaluation results
        """
        logger.info("Evaluating baseline similarity metrics...")

        # Initialize metrics
        metrics = {
            "Jaccard": JaccardSimilarity(),
            "Cosine": CosineSimilarity(self.ic_values),
            "Cosine (unweighted)": CosineSimilarity(),
            "Simplified Resnik": SimplifiedResnikSimilarity(self.ic_values),
        }

        results = []

        for metric_name, metric in metrics.items():
            logger.info(f"  Evaluating {metric_name}...")

            calc = PhenopacketSimilarityCalculator(metric)
            evaluator = PhenopacketEvaluator(
                calc,
                k_values=self.k_values,
                metrics=self.metrics
            )

            # Leave-one-out evaluation with phenotype-based ground truth
            # Ground truth: patients sharing >= 3 phenotypes are considered relevant
            ground_truth = self._compute_phenotype_overlap_ground_truth(min_overlap=3)

            eval_results = evaluator.evaluate_retrieval(
                self.phenopackets,
                self.phenopackets,
                ground_truth,
                verbose=False
            )

            # Store results
            for eval_metric, k_scores in eval_results.items():
                for k, score in k_scores.items():
                    results.append({
                        "similarity_metric": metric_name,
                        "eval_metric": eval_metric,
                        "k": k,
                        "score": score,
                    })

        df = pd.DataFrame(results)

        # Save results
        output_path = self.results_dir / "baseline_metrics.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Saved baseline results to {output_path}")

        return df

    def _compute_phenotype_overlap_ground_truth(
        self,
        min_overlap: int = 3
    ) -> Dict[str, set]:
        """
        Compute ground truth based on phenotype overlap.

        Two patients are considered relevant if they share >= min_overlap phenotypes.
        """
        # Build phenotype sets for each patient
        patient_phenotypes = {}
        for i, pp in enumerate(self.phenopackets):
            pp_id = pp.get("id")
            phenotypes = set()
            for f in pp.get("phenotypicFeatures", []):
                if not f.get("excluded", False):
                    phenotypes.add(f["type"]["id"])
            patient_phenotypes[pp_id] = (i, phenotypes)

        # Compute ground truth
        ground_truth = {}
        for pp_id, (idx, phenotypes) in patient_phenotypes.items():
            relevant = set()
            for other_id, (other_idx, other_phenotypes) in patient_phenotypes.items():
                if other_id != pp_id:
                    overlap = len(phenotypes & other_phenotypes)
                    if overlap >= min_overlap:
                        relevant.add(other_idx)
            ground_truth[pp_id] = relevant

        return ground_truth

    def evaluate_similarity_distributions(self) -> Dict[str, np.ndarray]:
        """
        Compute similarity score distributions for analysis.

        Returns:
            Dictionary mapping metric names to similarity matrices
        """
        logger.info("Computing similarity distributions...")

        metrics = {
            "Jaccard": JaccardSimilarity(),
            "Cosine (IC)": CosineSimilarity(self.ic_values),
            "Resnik": SimplifiedResnikSimilarity(self.ic_values),
        }

        distributions = {}

        # Use subset for speed
        subset = self.phenopackets[:min(100, len(self.phenopackets))]

        for name, metric in metrics.items():
            calc = PhenopacketSimilarityCalculator(metric)
            matrix = calc.compute_similarity_matrix(subset)

            # Extract upper triangle (excluding diagonal)
            triu_indices = np.triu_indices_from(matrix, k=1)
            distributions[name] = matrix[triu_indices]

        return distributions

    def evaluate_phenotype_clustering(self) -> Dict[str, Any]:
        """
        Analyze how well similarity metrics cluster patients by phenotype profiles.

        Returns:
            Clustering quality metrics
        """
        logger.info("Evaluating phenotype-based clustering...")

        # Use IC-weighted Cosine as primary metric
        metric = CosineSimilarity(self.ic_values)
        calc = PhenopacketSimilarityCalculator(metric)

        # Compute similarity matrix
        subset = self.phenopackets[:min(200, len(self.phenopackets))]
        sim_matrix = calc.compute_similarity_matrix(subset)

        # Basic clustering statistics
        results = {
            "n_patients": len(subset),
            "mean_similarity": float(np.mean(sim_matrix[np.triu_indices_from(sim_matrix, k=1)])),
            "std_similarity": float(np.std(sim_matrix[np.triu_indices_from(sim_matrix, k=1)])),
            "min_similarity": float(np.min(sim_matrix)),
            "max_off_diagonal": float(np.max(sim_matrix - np.eye(len(sim_matrix)))),
        }

        # Find highly similar patient pairs
        high_sim_threshold = 0.5
        high_sim_pairs = []
        for i in range(len(subset)):
            for j in range(i + 1, len(subset)):
                if sim_matrix[i, j] >= high_sim_threshold:
                    high_sim_pairs.append({
                        "patient1": subset[i]["id"],
                        "patient2": subset[j]["id"],
                        "similarity": float(sim_matrix[i, j]),
                    })

        results["high_similarity_pairs"] = len(high_sim_pairs)
        results["high_similarity_examples"] = high_sim_pairs[:5]

        return results

    def run_full_evaluation(self) -> Dict[str, Any]:
        """
        Run complete evaluation pipeline.

        Returns:
            Dictionary with all evaluation results
        """
        logger.info("=" * 70)
        logger.info("Running full DECIPHER evaluation pipeline")
        logger.info("=" * 70)

        results = {
            "timestamp": datetime.now().isoformat(),
            "n_patients": len(self.phenopackets),
        }

        # Dataset statistics
        logger.info("\n1. Computing dataset statistics...")
        results["dataset_stats"] = self.get_dataset_statistics()

        # Baseline evaluation
        logger.info("\n2. Evaluating baseline similarity metrics...")
        baseline_df = self.evaluate_baseline_similarities()
        results["baseline_summary"] = self._summarize_baseline_results(baseline_df)

        # Similarity distributions
        logger.info("\n3. Analyzing similarity distributions...")
        distributions = self.evaluate_similarity_distributions()
        results["similarity_distributions"] = {
            name: {
                "mean": float(np.mean(dist)),
                "std": float(np.std(dist)),
                "median": float(np.median(dist)),
            }
            for name, dist in distributions.items()
        }

        # Clustering analysis
        logger.info("\n4. Evaluating clustering quality...")
        results["clustering"] = self.evaluate_phenotype_clustering()

        # Save full results
        output_path = self.results_dir / "full_evaluation_results.json"
        with open(output_path, "w") as f:
            json.dump(results, f, indent=2, default=str)
        logger.info(f"\nSaved full results to {output_path}")

        return results

    def _summarize_baseline_results(self, df: pd.DataFrame) -> Dict[str, Any]:
        """Summarize baseline evaluation results."""
        summary = {}

        for metric_name in df["similarity_metric"].unique():
            metric_df = df[df["similarity_metric"] == metric_name]
            summary[metric_name] = {}

            for eval_metric in df["eval_metric"].unique():
                eval_df = metric_df[metric_df["eval_metric"] == eval_metric]
                summary[metric_name][eval_metric] = {
                    k: float(eval_df[eval_df["k"] == k]["score"].values[0])
                    for k in self.k_values
                    if k in eval_df["k"].values
                }

        return summary

    def print_summary(self, results: Dict[str, Any]) -> None:
        """Print a formatted summary of evaluation results."""
        print("\n" + "=" * 70)
        print("  DECIPHER Evaluation Summary")
        print("=" * 70)

        stats = results.get("dataset_stats", {})
        print(f"\nDataset Statistics:")
        print(f"  Patients: {stats.get('n_patients', 'N/A')}")
        print(f"  Unique phenotypes: {stats.get('n_unique_phenotypes', 'N/A')}")
        print(f"  Mean phenotypes/patient: {stats.get('phenotypes_per_patient', {}).get('mean', 0):.1f}")
        print(f"  Rare phenotypes (n=1): {stats.get('rare_phenotypes', 'N/A')}")

        baseline = results.get("baseline_summary", {})
        if baseline:
            print(f"\nBaseline Retrieval Performance (Recall@k):")
            print(f"  {'Metric':<25} R@5    R@10   R@20")
            print("  " + "-" * 50)
            for metric_name, metrics in baseline.items():
                recall = metrics.get("recall", {})
                print(f"  {metric_name:<25} "
                      f"{recall.get(5, 0):.3f}  "
                      f"{recall.get(10, 0):.3f}  "
                      f"{recall.get(20, 0):.3f}")

        clustering = results.get("clustering", {})
        if clustering:
            print(f"\nSimilarity Analysis:")
            print(f"  Mean pairwise similarity: {clustering.get('mean_similarity', 0):.4f}")
            print(f"  High-similarity pairs (>0.5): {clustering.get('high_similarity_pairs', 0)}")

        print("\n" + "=" * 70)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Evaluate phenotype matching on DECIPHER data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        "--data-file",
        type=str,
        help="Path to existing DECIPHER data file (JSON)"
    )
    parser.add_argument(
        "--n-patients",
        type=int,
        default=500,
        help="Number of patients for evaluation (default: 500)"
    )
    parser.add_argument(
        "--simulated",
        action="store_true",
        help="Use simulated DECIPHER-like data"
    )
    parser.add_argument(
        "--min-phenotypes",
        type=int,
        default=3,
        help="Minimum phenotypes per patient (default: 3)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results/decipher",
        help="Output directory for results"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Run evaluation
    evaluator = DECIPHEREvaluator(results_dir=args.output_dir)

    # Load data
    evaluator.load_data(
        data_file=args.data_file,
        n_patients=args.n_patients,
        use_simulated=args.simulated,
        min_phenotypes=args.min_phenotypes
    )

    # Run full evaluation
    results = evaluator.run_full_evaluation()

    # Print summary
    evaluator.print_summary(results)

    print(f"\nResults saved to: {args.output_dir}")


if __name__ == "__main__":
    main()
