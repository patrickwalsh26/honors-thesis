#!/usr/bin/env python3
"""
Evaluate phenotype matching on real disease data from HPO annotations.

This script evaluates baseline similarity metrics and privacy-preserving
approaches using synthetic patients generated from real OMIM/Orphanet
disease-phenotype associations.

Results are suitable for thesis figures and tables.
"""

import json
import logging
import sys
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from typing import Dict, List, Tuple

import numpy as np

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_integration.hpoa_patient_generator import (
    HPOAPatientGenerator,
    SyntheticPatient,
)
from src.similarity.hpo_similarity import (
    JaccardSimilarity,
    CosineSimilarity,
    SimplifiedResnikSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class HPOAEvaluator:
    """Evaluate phenotype matching on HPOA-derived patient cohorts."""

    def __init__(
        self,
        data_dir: str = "data/hpoa_evaluation",
        results_dir: str = "results/hpoa"
    ):
        self.data_dir = Path(data_dir)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)

        self.phenopackets = None
        self.ground_truth = None
        self.metadata = None
        self.ic_values = None

    def load_data(self) -> None:
        """Load pre-generated evaluation data."""
        logger.info(f"Loading data from {self.data_dir}")

        # Load phenopackets
        with open(self.data_dir / "cohort_phenopackets.json") as f:
            self.phenopackets = json.load(f)

        # Load ground truth
        with open(self.data_dir / "ground_truth.json") as f:
            self.ground_truth = json.load(f)

        # Load metadata
        with open(self.data_dir / "metadata.json") as f:
            self.metadata = json.load(f)

        logger.info(
            f"Loaded {len(self.phenopackets)} patients from "
            f"{self.metadata['n_diseases']} diseases"
        )

        # Compute IC values from corpus
        logger.info("Computing information content from corpus...")
        self.ic_values = compute_empirical_ic(self.phenopackets)
        logger.info(f"Computed IC for {len(self.ic_values)} phenotypes")

    def generate_fresh_cohort(
        self,
        n_patients: int = 500,
        n_diseases: int = 100,
        seed: int = 42
    ) -> None:
        """Generate fresh cohort instead of loading."""
        logger.info("Generating fresh cohort from HPOA...")

        generator = HPOAPatientGenerator(
            hpoa_path="data/hpo_annotations/phenotype.hpoa",
            seed=seed
        )

        patients, self.metadata = generator.generate_cohort(
            n_patients=n_patients,
            n_diseases=n_diseases,
            min_phenotypes=5,
            max_phenotypes=30,
            phenotype_recall=0.75,
            noise_rate=0.1,
            balance_diseases=True
        )

        self.phenopackets = [p.to_phenopacket() for p in patients]
        self.ground_truth = generator.create_ground_truth(
            patients, same_disease_relevant=True
        )

        # Compute IC
        self.ic_values = compute_empirical_ic(self.phenopackets)

        logger.info(f"Generated {len(self.phenopackets)} patients")

    def compute_precision_at_k(
        self,
        ranked_results: List[int],
        relevant_ids: List[str],
        patient_ids: List[str],
        k: int
    ) -> float:
        """Compute Precision@K."""
        if k == 0 or not relevant_ids:
            return 0.0

        relevant_set = set(relevant_ids)
        top_k = ranked_results[:k]
        hits = sum(1 for idx in top_k if patient_ids[idx] in relevant_set)
        return hits / k

    def compute_recall_at_k(
        self,
        ranked_results: List[int],
        relevant_ids: List[str],
        patient_ids: List[str],
        k: int
    ) -> float:
        """Compute Recall@K."""
        if not relevant_ids:
            return 0.0

        relevant_set = set(relevant_ids)
        top_k = ranked_results[:k]
        hits = sum(1 for idx in top_k if patient_ids[idx] in relevant_set)
        return hits / len(relevant_set)

    def compute_ndcg_at_k(
        self,
        ranked_results: List[int],
        relevant_ids: List[str],
        patient_ids: List[str],
        k: int
    ) -> float:
        """Compute NDCG@K."""
        if not relevant_ids or k == 0:
            return 0.0

        relevant_set = set(relevant_ids)

        # DCG
        dcg = 0.0
        for i, idx in enumerate(ranked_results[:k]):
            if patient_ids[idx] in relevant_set:
                dcg += 1.0 / np.log2(i + 2)

        # Ideal DCG
        n_relevant = min(len(relevant_set), k)
        idcg = sum(1.0 / np.log2(i + 2) for i in range(n_relevant))

        return dcg / idcg if idcg > 0 else 0.0

    def compute_hit_rate_at_k(
        self,
        ranked_results: List[int],
        relevant_ids: List[str],
        patient_ids: List[str],
        k: int
    ) -> float:
        """Compute Hit Rate@K (at least one relevant in top-K)."""
        if not relevant_ids or k == 0:
            return 0.0

        relevant_set = set(relevant_ids)
        top_k = ranked_results[:k]
        return 1.0 if any(patient_ids[idx] in relevant_set for idx in top_k) else 0.0

    def compute_mrr(
        self,
        ranked_results: List[int],
        relevant_ids: List[str],
        patient_ids: List[str]
    ) -> float:
        """Compute Mean Reciprocal Rank."""
        if not relevant_ids:
            return 0.0

        relevant_set = set(relevant_ids)
        for i, idx in enumerate(ranked_results):
            if patient_ids[idx] in relevant_set:
                return 1.0 / (i + 1)
        return 0.0

    def evaluate_baseline_similarities(
        self,
        k_values: List[int] = [1, 5, 10, 20, 50]
    ) -> Dict:
        """Evaluate baseline similarity metrics."""
        logger.info("Evaluating baseline similarity metrics...")

        metrics = {
            "Jaccard": JaccardSimilarity(),
            "Cosine": CosineSimilarity(self.ic_values),
            "Cosine (unweighted)": CosineSimilarity(),
            "Simplified Resnik": SimplifiedResnikSimilarity(self.ic_values),
        }

        # Build patient ID mapping
        patient_ids = [pp["subject"]["id"] for pp in self.phenopackets]
        id_to_idx = {pid: i for i, pid in enumerate(patient_ids)}

        results = {}

        for metric_name, metric in metrics.items():
            logger.info(f"  Evaluating {metric_name}...")
            calc = PhenopacketSimilarityCalculator(metric)

            # Compute similarity matrix
            sim_matrix = calc.compute_similarity_matrix(self.phenopackets)

            # Evaluate retrieval for each patient
            metric_results = {
                "precision": {k: [] for k in k_values},
                "recall": {k: [] for k in k_values},
                "ndcg": {k: [] for k in k_values},
                "hit_rate": {k: [] for k in k_values},
                "mrr": [],
            }

            for i, pp in enumerate(self.phenopackets):
                patient_id = patient_ids[i]
                relevant = self.ground_truth.get(patient_id, [])

                if not relevant:
                    continue

                # Get rankings (excluding self)
                scores = sim_matrix[i].copy()
                scores[i] = -np.inf  # Exclude self
                ranked = np.argsort(scores)[::-1]

                # Compute metrics at each K
                for k in k_values:
                    metric_results["precision"][k].append(
                        self.compute_precision_at_k(ranked, relevant, patient_ids, k)
                    )
                    metric_results["recall"][k].append(
                        self.compute_recall_at_k(ranked, relevant, patient_ids, k)
                    )
                    metric_results["ndcg"][k].append(
                        self.compute_ndcg_at_k(ranked, relevant, patient_ids, k)
                    )
                    metric_results["hit_rate"][k].append(
                        self.compute_hit_rate_at_k(ranked, relevant, patient_ids, k)
                    )

                metric_results["mrr"].append(
                    self.compute_mrr(ranked, relevant, patient_ids)
                )

            # Average results
            results[metric_name] = {
                "precision": {k: np.mean(v) for k, v in metric_results["precision"].items()},
                "recall": {k: np.mean(v) for k, v in metric_results["recall"].items()},
                "ndcg": {k: np.mean(v) for k, v in metric_results["ndcg"].items()},
                "hit_rate": {k: np.mean(v) for k, v in metric_results["hit_rate"].items()},
                "mrr": np.mean(metric_results["mrr"]),
            }

        return results

    def evaluate_similarity_distributions(self) -> Dict:
        """Analyze similarity score distributions."""
        logger.info("Analyzing similarity distributions...")

        metrics = {
            "Jaccard": JaccardSimilarity(),
            "Cosine (IC)": CosineSimilarity(self.ic_values),
            "Resnik": SimplifiedResnikSimilarity(self.ic_values),
        }

        results = {}

        for name, metric in metrics.items():
            calc = PhenopacketSimilarityCalculator(metric)
            sim_matrix = calc.compute_similarity_matrix(self.phenopackets)

            # Get off-diagonal values
            off_diag = sim_matrix[np.triu_indices_from(sim_matrix, k=1)]

            results[name] = {
                "mean": float(np.mean(off_diag)),
                "std": float(np.std(off_diag)),
                "median": float(np.median(off_diag)),
                "min": float(np.min(off_diag)),
                "max": float(np.max(off_diag)),
                "q25": float(np.percentile(off_diag, 25)),
                "q75": float(np.percentile(off_diag, 75)),
            }

        return results

    def evaluate_by_disease_type(self) -> Dict:
        """Analyze performance by disease source (OMIM vs Orphanet)."""
        logger.info("Analyzing performance by disease source...")

        # Load patient details to get disease info
        with open(self.data_dir / "cohort_patients.json") as f:
            patients_data = json.load(f)

        # Group by source
        omim_patients = []
        orphanet_patients = []

        for p in patients_data:
            if p["underlying_disease"].startswith("OMIM:"):
                omim_patients.append(p["patient_id"])
            elif p["underlying_disease"].startswith("ORPHA:"):
                orphanet_patients.append(p["patient_id"])

        return {
            "omim_count": len(omim_patients),
            "orphanet_count": len(orphanet_patients),
            "omim_fraction": len(omim_patients) / len(patients_data),
            "orphanet_fraction": len(orphanet_patients) / len(patients_data),
        }

    def run_full_evaluation(self) -> Dict:
        """Run complete evaluation pipeline."""
        logger.info("=" * 60)
        logger.info("Running full HPOA evaluation pipeline")
        logger.info("=" * 60)

        # Load or generate data
        if self.phenopackets is None:
            if (self.data_dir / "cohort_phenopackets.json").exists():
                self.load_data()
            else:
                self.generate_fresh_cohort()

        # Dataset statistics
        phenotypes_per_patient = [
            len(pp["phenotypicFeatures"]) for pp in self.phenopackets
        ]

        # Get unique phenotypes
        all_phenotypes = set()
        for pp in self.phenopackets:
            for f in pp["phenotypicFeatures"]:
                all_phenotypes.add(f["type"]["id"])

        dataset_stats = {
            "n_patients": len(self.phenopackets),
            "n_diseases": self.metadata.get("n_diseases", 0),
            "n_unique_phenotypes": len(all_phenotypes),
            "phenotypes_per_patient": {
                "mean": np.mean(phenotypes_per_patient),
                "std": np.std(phenotypes_per_patient),
                "median": np.median(phenotypes_per_patient),
                "min": min(phenotypes_per_patient),
                "max": max(phenotypes_per_patient),
            },
            "generation_params": self.metadata.get("generation_params", {}),
        }

        # Baseline evaluation
        baseline_results = self.evaluate_baseline_similarities()

        # Similarity distributions
        distributions = self.evaluate_similarity_distributions()

        # Disease type analysis
        disease_analysis = self.evaluate_by_disease_type()

        # Compile results
        results = {
            "timestamp": datetime.now().isoformat(),
            "dataset_stats": dataset_stats,
            "baseline_summary": baseline_results,
            "similarity_distributions": distributions,
            "disease_analysis": disease_analysis,
        }

        # Save results
        output_file = self.results_dir / "full_evaluation_results.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)
        logger.info(f"Saved results to {output_file}")

        # Save baseline metrics as CSV
        self._save_baseline_csv(baseline_results)

        return results

    def _save_baseline_csv(self, results: Dict) -> None:
        """Save baseline results as CSV for easy plotting."""
        csv_path = self.results_dir / "baseline_metrics.csv"

        with open(csv_path, "w") as f:
            f.write("similarity_metric,eval_metric,k,score\n")

            for sim_name, sim_results in results.items():
                for eval_metric in ["precision", "recall", "ndcg", "hit_rate"]:
                    for k, score in sim_results[eval_metric].items():
                        f.write(f"{sim_name},{eval_metric},{k},{score}\n")

                # MRR doesn't have K
                f.write(f"{sim_name},mrr,0,{sim_results['mrr']}\n")

        logger.info(f"Saved baseline CSV to {csv_path}")

    def print_summary(self, results: Dict) -> None:
        """Print human-readable summary."""
        print("\n" + "=" * 70)
        print("  HPOA Evaluation Results - Real Disease Data")
        print("=" * 70)

        stats = results["dataset_stats"]
        print(f"\nDataset Statistics:")
        print(f"  Patients: {stats['n_patients']}")
        print(f"  Diseases: {stats['n_diseases']}")
        print(f"  Unique phenotypes: {stats['n_unique_phenotypes']}")
        print(f"  Phenotypes/patient: {stats['phenotypes_per_patient']['mean']:.1f} "
              f"(range: {stats['phenotypes_per_patient']['min']}-"
              f"{stats['phenotypes_per_patient']['max']})")

        print(f"\nBaseline Retrieval Performance:")
        print("-" * 70)
        print(f"{'Metric':<25} {'P@5':>8} {'P@10':>8} {'R@5':>8} {'R@10':>8} "
              f"{'NDCG@10':>8} {'MRR':>8}")
        print("-" * 70)

        for sim_name, sim_results in results["baseline_summary"].items():
            print(f"{sim_name:<25} "
                  f"{sim_results['precision'][5]:>8.3f} "
                  f"{sim_results['precision'][10]:>8.3f} "
                  f"{sim_results['recall'][5]:>8.3f} "
                  f"{sim_results['recall'][10]:>8.3f} "
                  f"{sim_results['ndcg'][10]:>8.3f} "
                  f"{sim_results['mrr']:>8.3f}")

        print("\n" + "=" * 70)

        # Best performing metric
        best_metric = max(
            results["baseline_summary"].items(),
            key=lambda x: x[1]["ndcg"][10]
        )
        print(f"\nBest performing metric: {best_metric[0]}")
        print(f"  NDCG@10: {best_metric[1]['ndcg'][10]:.3f}")
        print(f"  MRR: {best_metric[1]['mrr']:.3f}")


def main():
    """Run HPOA evaluation."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Evaluate phenotype matching on HPOA-derived data"
    )
    parser.add_argument(
        "--n-patients", type=int, default=500,
        help="Number of patients (default: 500)"
    )
    parser.add_argument(
        "--n-diseases", type=int, default=100,
        help="Number of diseases to sample (default: 100)"
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed (default: 42)"
    )
    parser.add_argument(
        "--regenerate", action="store_true",
        help="Regenerate cohort even if exists"
    )
    args = parser.parse_args()

    evaluator = HPOAEvaluator()

    # Load or generate
    if args.regenerate or not (evaluator.data_dir / "cohort_phenopackets.json").exists():
        evaluator.generate_fresh_cohort(
            n_patients=args.n_patients,
            n_diseases=args.n_diseases,
            seed=args.seed
        )
    else:
        evaluator.load_data()

    # Run evaluation
    results = evaluator.run_full_evaluation()

    # Print summary
    evaluator.print_summary(results)


if __name__ == "__main__":
    main()
