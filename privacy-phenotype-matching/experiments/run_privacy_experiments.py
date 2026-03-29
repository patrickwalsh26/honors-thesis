#!/usr/bin/env python3
"""
Comprehensive privacy-utility evaluation experiments.

Generates all experimental data for thesis results including:
- Baseline (non-private) performance
- DP epsilon sweep experiments
- k-anonymity parameter sweep
- Privacy leakage audits
- Combined mechanism evaluation

Usage:
    python experiments/run_privacy_experiments.py --data data/synthetic/mixed_cohort_200.json
    python experiments/run_privacy_experiments.py --experiments baseline dp k_anon leakage
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from datetime import datetime
import pandas as pd
import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.similarity.hpo_similarity import (
    CosineSimilarity, JaccardSimilarity, SimplifiedResnikSimilarity,
    PhenopacketSimilarityCalculator, compute_empirical_ic, load_phenopackets
)
from src.privacy.privacy_calculator import PrivacyPreservingCalculator
from src.privacy.k_anonymity import PrivacyConfig
from src.evaluation.metrics import PhenopacketEvaluator, RetrievalMetrics
from src.evaluation.privacy_utility import PrivacyUtilityFrontier

logger = logging.getLogger(__name__)


def setup_logging(output_dir: Path, verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO

    handlers = [
        logging.StreamHandler(),
        logging.FileHandler(output_dir / "experiment.log")
    ]

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def run_baseline_experiments(
    phenopackets: list,
    ic_values: dict,
    output_dir: Path
) -> pd.DataFrame:
    """
    Run baseline (non-private) experiments.

    Tests different similarity metrics without privacy mechanisms.
    """
    logger.info("Running baseline experiments...")

    metrics = {
        "jaccard": JaccardSimilarity(),
        "cosine_ic": CosineSimilarity(ic_values),
        "resnik_simplified": SimplifiedResnikSimilarity(ic_values)
    }

    k_values = [1, 5, 10, 20]
    results = []

    for metric_name, metric in metrics.items():
        logger.info(f"  Testing {metric_name}...")

        calc = PhenopacketSimilarityCalculator(metric)
        evaluator = PhenopacketEvaluator(calc, k_values=k_values)
        eval_results = evaluator.cross_disease_evaluation(phenopackets)

        for k in k_values:
            results.append({
                "metric": metric_name,
                "privacy": "none",
                "k": k,
                "recall": eval_results["recall"][k],
                "precision": eval_results["precision"][k],
                "f1": eval_results["f1"][k],
                "ndcg": eval_results["ndcg"][k]
            })

    df = pd.DataFrame(results)
    df.to_csv(output_dir / "baseline_results.csv", index=False)

    logger.info(f"Baseline results saved to {output_dir / 'baseline_results.csv'}")
    return df


def run_dp_experiments(
    phenopackets: list,
    ic_values: dict,
    output_dir: Path,
    epsilon_values: list = None,
    n_trials: int = 10
) -> pd.DataFrame:
    """
    Run differential privacy experiments.

    Sweeps epsilon parameter and measures utility degradation.
    """
    logger.info("Running DP experiments...")

    if epsilon_values is None:
        epsilon_values = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]

    k_values = [5, 10, 20]
    base_metric = CosineSimilarity(ic_values)
    base_calc = PhenopacketSimilarityCalculator(base_metric)

    # Compute baseline for comparison
    evaluator = PhenopacketEvaluator(base_calc, k_values=k_values)
    baseline = evaluator.cross_disease_evaluation(phenopackets)

    results = []

    for epsilon in epsilon_values:
        logger.info(f"  Testing epsilon = {epsilon}...")

        for trial in range(n_trials):
            # Create private calculator
            config = PrivacyConfig(
                epsilon=epsilon,
                k=1,  # No k-anonymity suppression
                min_prevalence=0.0  # No rare term filtering
            )

            private_calc = PrivacyPreservingCalculator(
                base_similarity=base_metric,
                config=config,
                use_psi=False,
                use_dp=True,
                use_k_anonymity=False,
                use_rare_filter=False
            )

            # Evaluate
            private_eval = PhenopacketEvaluator(private_calc, k_values=k_values)
            metrics = private_eval.cross_disease_evaluation(phenopackets)

            for k in k_values:
                baseline_recall = baseline["recall"][k]
                results.append({
                    "epsilon": epsilon,
                    "trial": trial,
                    "k": k,
                    "recall": metrics["recall"][k],
                    "precision": metrics["precision"][k],
                    "f1": metrics["f1"][k],
                    "ndcg": metrics["ndcg"][k],
                    "baseline_recall": baseline_recall,
                    "relative_recall": metrics["recall"][k] / baseline_recall if baseline_recall > 0 else 0
                })

    df = pd.DataFrame(results)
    df.to_csv(output_dir / "dp_sweep_results.csv", index=False)

    logger.info(f"DP results saved to {output_dir / 'dp_sweep_results.csv'}")
    return df


def run_k_anonymity_experiments(
    phenopackets: list,
    ic_values: dict,
    output_dir: Path,
    k_values: list = None
) -> pd.DataFrame:
    """
    Run k-anonymity experiments.

    Sweeps k parameter and measures query suppression rates.
    """
    logger.info("Running k-anonymity experiments...")

    if k_values is None:
        k_values = [2, 3, 5, 10, 20, 50]

    retrieval_k_values = [5, 10, 20]
    base_metric = CosineSimilarity(ic_values)
    base_calc = PhenopacketSimilarityCalculator(base_metric)

    # Compute baseline
    evaluator = PhenopacketEvaluator(base_calc, k_values=retrieval_k_values)
    baseline = evaluator.cross_disease_evaluation(phenopackets)

    # Build ground truth for recall computation
    from collections import defaultdict
    by_disease = defaultdict(list)
    for i, pp in enumerate(phenopackets):
        diseases = pp.get("diseases", [])
        disease = diseases[0]["term"]["label"] if diseases else "Unknown"
        by_disease[disease].append(i)

    ground_truth = {}
    for i, pp in enumerate(phenopackets):
        diseases = pp.get("diseases", [])
        disease = diseases[0]["term"]["label"] if diseases else "Unknown"
        ground_truth[pp.get("id")] = {j for j in by_disease[disease] if j != i}

    results = []

    for k_anon in k_values:
        logger.info(f"  Testing k-anonymity k = {k_anon}...")

        config = PrivacyConfig(
            k=k_anon,
            epsilon=float('inf')  # No DP noise
        )

        private_calc = PrivacyPreservingCalculator(
            base_similarity=base_metric,
            config=config,
            use_psi=False,
            use_dp=False,
            use_k_anonymity=True,
            use_rare_filter=False
        )

        # Count successful vs suppressed queries
        for retrieval_k in retrieval_k_values:
            successful = 0
            suppressed = 0
            recall_sum = 0.0

            for pp in phenopackets:
                result = private_calc.find_most_similar(
                    pp, phenopackets, top_k=retrieval_k
                )

                relevant = ground_truth.get(pp.get("id"), set())

                if result is None:
                    suppressed += 1
                else:
                    successful += 1
                    retrieved = [idx for idx, _ in result]
                    recall = RetrievalMetrics.recall_at_k(retrieved, relevant, retrieval_k)
                    recall_sum += recall

            n_queries = len(phenopackets)
            baseline_recall = baseline["recall"][retrieval_k]

            results.append({
                "k_anonymity": k_anon,
                "retrieval_k": retrieval_k,
                "success_rate": successful / n_queries,
                "suppression_rate": suppressed / n_queries,
                "successful_queries": successful,
                "suppressed_queries": suppressed,
                "conditional_recall": recall_sum / successful if successful > 0 else 0,
                "baseline_recall": baseline_recall,
                "expected_recall": (successful / n_queries) * (recall_sum / successful if successful > 0 else 0)
            })

        # Reset statistics
        private_calc.k_guard.reset_statistics()

    df = pd.DataFrame(results)
    df.to_csv(output_dir / "k_anonymity_results.csv", index=False)

    logger.info(f"k-anonymity results saved to {output_dir / 'k_anonymity_results.csv'}")
    return df


def run_leakage_audit(
    phenopackets: list,
    ic_values: dict,
    output_dir: Path
) -> dict:
    """
    Run privacy leakage audits.

    Evaluates membership and attribute inference attacks.
    """
    logger.info("Running leakage audits...")

    try:
        from src.evaluation.leakage_audit import LeakageAuditReport
    except ImportError as e:
        logger.warning(f"Leakage audit skipped: {e}")
        return {"error": str(e)}

    base_metric = CosineSimilarity(ic_values)
    base_calc = PhenopacketSimilarityCalculator(base_metric)

    # Private calculator
    config = PrivacyConfig(epsilon=1.0, k=5)
    private_calc = PrivacyPreservingCalculator(
        base_similarity=base_metric,
        config=config,
        use_psi=False,
        use_dp=True,
        use_k_anonymity=True,
        use_rare_filter=False
    )

    audit = LeakageAuditReport(base_calc, private_calc, phenopackets)

    try:
        results = audit.run_full_audit(verbose=True)
    except Exception as e:
        logger.error(f"Leakage audit failed: {e}")
        results = {"error": str(e)}

    with open(output_dir / "leakage_audit.json", "w") as f:
        json.dump(results, f, indent=2, default=str)

    logger.info(f"Leakage audit results saved to {output_dir / 'leakage_audit.json'}")
    return results


def run_combined_experiments(
    phenopackets: list,
    ic_values: dict,
    output_dir: Path,
    n_trials: int = 5
) -> pd.DataFrame:
    """
    Run combined DP + k-anonymity experiments.
    """
    logger.info("Running combined mechanism experiments...")

    epsilon_values = [0.5, 1.0, 2.0, 5.0]
    k_anon_values = [3, 5, 10]
    retrieval_k_values = [5, 10, 20]

    base_metric = CosineSimilarity(ic_values)
    base_calc = PhenopacketSimilarityCalculator(base_metric)

    # Baseline
    evaluator = PhenopacketEvaluator(base_calc, k_values=retrieval_k_values)
    baseline = evaluator.cross_disease_evaluation(phenopackets)

    results = []
    total_combos = len(epsilon_values) * len(k_anon_values) * n_trials
    combo_num = 0

    for epsilon in epsilon_values:
        for k_anon in k_anon_values:
            for trial in range(n_trials):
                combo_num += 1
                if combo_num % 10 == 0:
                    logger.info(f"  Progress: {combo_num}/{total_combos}")

                config = PrivacyConfig(epsilon=epsilon, k=k_anon)

                private_calc = PrivacyPreservingCalculator(
                    base_similarity=base_metric,
                    config=config,
                    use_psi=False,
                    use_dp=True,
                    use_k_anonymity=True,
                    use_rare_filter=False
                )

                private_eval = PhenopacketEvaluator(private_calc, k_values=retrieval_k_values)
                metrics = private_eval.cross_disease_evaluation(phenopackets)

                k_stats = private_calc.k_guard.get_statistics()

                for k in retrieval_k_values:
                    results.append({
                        "epsilon": epsilon,
                        "k_anonymity": k_anon,
                        "trial": trial,
                        "retrieval_k": k,
                        "recall": metrics["recall"][k],
                        "precision": metrics["precision"][k],
                        "ndcg": metrics["ndcg"][k],
                        "baseline_recall": baseline["recall"][k],
                        "relative_recall": metrics["recall"][k] / baseline["recall"][k] if baseline["recall"][k] > 0 else 0,
                        "suppression_rate": k_stats.get("suppression_rate", 0)
                    })

    df = pd.DataFrame(results)
    df.to_csv(output_dir / "combined_results.csv", index=False)

    logger.info(f"Combined results saved to {output_dir / 'combined_results.csv'}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Run privacy-utility experiments for phenotype matching"
    )
    parser.add_argument(
        "--data",
        default="data/synthetic/mixed_cohort_200.json",
        help="Path to phenopacket data"
    )
    parser.add_argument(
        "--output",
        default="experiments/results",
        help="Output directory for results"
    )
    parser.add_argument(
        "--experiments",
        nargs="+",
        default=["baseline", "dp", "k_anon", "combined"],
        choices=["baseline", "dp", "k_anon", "leakage", "combined", "all"],
        help="Experiments to run"
    )
    parser.add_argument(
        "--n-trials",
        type=int,
        default=10,
        help="Number of trials for stochastic experiments"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Setup
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_logging(output_dir, args.verbose)

    # Expand "all" to all experiments
    if "all" in args.experiments:
        experiments = ["baseline", "dp", "k_anon", "leakage", "combined"]
    else:
        experiments = args.experiments

    # Load data
    logger.info(f"Loading phenopackets from {args.data}...")
    try:
        phenopackets = load_phenopackets(args.data)
    except FileNotFoundError:
        logger.error(f"Data file not found: {args.data}")
        sys.exit(1)

    logger.info(f"Loaded {len(phenopackets)} phenopackets")

    # Compute IC values
    ic_values = compute_empirical_ic(phenopackets)
    logger.info(f"Computed IC for {len(ic_values)} terms")

    # Run experiments
    results = {}

    if "baseline" in experiments:
        results["baseline"] = run_baseline_experiments(phenopackets, ic_values, output_dir)

    if "dp" in experiments:
        results["dp"] = run_dp_experiments(
            phenopackets, ic_values, output_dir, n_trials=args.n_trials
        )

    if "k_anon" in experiments:
        results["k_anon"] = run_k_anonymity_experiments(phenopackets, ic_values, output_dir)

    if "leakage" in experiments:
        results["leakage"] = run_leakage_audit(phenopackets, ic_values, output_dir)

    if "combined" in experiments:
        results["combined"] = run_combined_experiments(
            phenopackets, ic_values, output_dir, n_trials=args.n_trials // 2
        )

    # Save summary
    summary = {
        "timestamp": datetime.now().isoformat(),
        "data_path": args.data,
        "n_phenopackets": len(phenopackets),
        "n_unique_terms": len(ic_values),
        "experiments_run": experiments,
        "n_trials": args.n_trials
    }

    with open(output_dir / "experiment_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    logger.info(f"All results saved to {output_dir}")
    logger.info("Experiment complete!")


if __name__ == "__main__":
    main()
