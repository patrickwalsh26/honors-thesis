"""
Privacy-utility frontier analysis for phenotype matching.

Implements:
- Parameter sweeps (epsilon, k-anonymity)
- Privacy-utility tradeoff visualization
- Comparative analysis across privacy mechanisms
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
import logging
from pathlib import Path

from .metrics import PhenopacketEvaluator, RetrievalMetrics

logger = logging.getLogger(__name__)


class PrivacyUtilityFrontier:
    """
    Analyze and visualize privacy-utility tradeoffs.

    Sweeps privacy parameters and measures impact on retrieval utility.
    """

    def __init__(
        self,
        base_calculator,
        phenopackets: List[Dict],
        hpo_manager=None,
        k_values: List[int] = [5, 10, 20]
    ):
        """
        Initialize frontier analyzer.

        Args:
            base_calculator: Base PhenopacketSimilarityCalculator (no privacy)
            phenopackets: List of phenopackets for evaluation
            hpo_manager: HPOManager for privacy calculator
            k_values: List of k values for @k metrics
        """
        self.base_calculator = base_calculator
        self.phenopackets = phenopackets
        self.hpo_manager = hpo_manager
        self.k_values = k_values

        # Pre-compute ground truth (same disease = relevant)
        self._ground_truth = self._build_ground_truth()

        # Store results
        self._results: List[Dict] = []

    def _build_ground_truth(self) -> Dict[str, Set[int]]:
        """Build same-disease ground truth."""
        by_disease = defaultdict(list)
        for i, pp in enumerate(self.phenopackets):
            diseases = pp.get("diseases", [])
            if diseases:
                disease = diseases[0]["term"].get("label", "Unknown")
            else:
                disease = "Unknown"
            by_disease[disease].append(i)

        ground_truth = {}
        for i, pp in enumerate(self.phenopackets):
            query_id = pp.get("id")
            diseases = pp.get("diseases", [])
            if diseases:
                disease = diseases[0]["term"].get("label", "Unknown")
            else:
                disease = "Unknown"

            # Same disease (excluding self)
            ground_truth[query_id] = {j for j in by_disease[disease] if j != i}

        return ground_truth

    def compute_baseline(self) -> Dict[str, Dict[int, float]]:
        """
        Compute baseline (non-private) performance.

        Returns:
            Dict of metric_name -> {k: score}
        """
        evaluator = PhenopacketEvaluator(self.base_calculator, self.k_values)
        return evaluator.evaluate_retrieval(
            self.phenopackets,
            self.phenopackets,
            self._ground_truth
        )

    def sweep_epsilon(
        self,
        epsilon_values: List[float],
        n_trials: int = 10,
        mechanism: str = "laplace",
        use_psi: bool = False,
        use_k_anonymity: bool = False,
        k_anonymity_k: int = 5
    ) -> pd.DataFrame:
        """
        Sweep differential privacy epsilon parameter.

        Args:
            epsilon_values: List of epsilon values to test
            n_trials: Number of trials per epsilon (for variance estimation)
            mechanism: DP mechanism ("laplace" or "gaussian")
            use_psi: Enable PSI alongside DP
            use_k_anonymity: Enable k-anonymity alongside DP
            k_anonymity_k: k value for k-anonymity

        Returns:
            DataFrame with results
        """
        from ..privacy.privacy_calculator import PrivacyPreservingCalculator
        from ..privacy.k_anonymity import PrivacyConfig

        # Compute baseline first
        baseline = self.compute_baseline()

        results = []

        for epsilon in epsilon_values:
            logger.info(f"Testing epsilon = {epsilon}")

            for trial in range(n_trials):
                # Create privacy calculator
                config = PrivacyConfig(
                    epsilon=epsilon,
                    dp_mechanism=mechanism,
                    k=k_anonymity_k if use_k_anonymity else 1
                )

                # Get base similarity metric from calculator
                base_metric = None
                if hasattr(self.base_calculator, 'similarity_metric'):
                    base_metric = self.base_calculator.similarity_metric

                private_calc = PrivacyPreservingCalculator(
                    base_similarity=base_metric,
                    config=config,
                    hpo_manager=self.hpo_manager,
                    use_psi=use_psi,
                    use_dp=True,
                    use_k_anonymity=use_k_anonymity,
                    use_rare_filter=False
                )

                # Evaluate
                evaluator = PhenopacketEvaluator(private_calc, self.k_values)
                metrics = evaluator.evaluate_retrieval(
                    self.phenopackets,
                    self.phenopackets,
                    self._ground_truth
                )

                # Record results
                for k in self.k_values:
                    result = {
                        "epsilon": epsilon,
                        "trial": trial,
                        "k": k,
                        "mechanism": mechanism,
                        "use_psi": use_psi,
                        "use_k_anonymity": use_k_anonymity,
                        "k_anonymity_k": k_anonymity_k if use_k_anonymity else None
                    }

                    for metric_name in ["precision", "recall", "f1", "ndcg"]:
                        if metric_name in metrics:
                            result[metric_name] = metrics[metric_name][k]
                            # Relative to baseline
                            baseline_val = baseline.get(metric_name, {}).get(k, 0)
                            result[f"relative_{metric_name}"] = (
                                metrics[metric_name][k] / baseline_val
                                if baseline_val > 0 else 0
                            )

                    results.append(result)

        df = pd.DataFrame(results)
        self._results.extend(results)

        return df

    def sweep_k_anonymity(
        self,
        k_values: List[int],
        retrieval_k: int = 10,
        use_dp: bool = False,
        epsilon: float = 1.0
    ) -> pd.DataFrame:
        """
        Sweep k-anonymity parameter.

        Args:
            k_values: List of k-anonymity k values to test
            retrieval_k: k for retrieval evaluation
            use_dp: Enable DP alongside k-anonymity
            epsilon: Epsilon for DP (if enabled)

        Returns:
            DataFrame with results
        """
        from ..privacy.privacy_calculator import PrivacyPreservingCalculator
        from ..privacy.k_anonymity import PrivacyConfig

        # Compute baseline
        baseline = self.compute_baseline()

        results = []

        for k in k_values:
            logger.info(f"Testing k-anonymity k = {k}")

            config = PrivacyConfig(
                k=k,
                epsilon=epsilon if use_dp else float('inf')
            )

            base_metric = None
            if hasattr(self.base_calculator, 'similarity_metric'):
                base_metric = self.base_calculator.similarity_metric

            private_calc = PrivacyPreservingCalculator(
                base_similarity=base_metric,
                config=config,
                hpo_manager=self.hpo_manager,
                use_psi=False,
                use_dp=use_dp,
                use_k_anonymity=True,
                use_rare_filter=False
            )

            # Count successful vs suppressed queries
            successful = 0
            suppressed = 0
            utility_sum = {k_ret: 0.0 for k_ret in self.k_values}
            recall_sum = {k_ret: 0.0 for k_ret in self.k_values}

            for query in self.phenopackets:
                result = private_calc.find_most_similar(
                    query, self.phenopackets, top_k=max(self.k_values)
                )

                query_id = query.get("id")
                relevant = self._ground_truth.get(query_id, set())

                if result is None:
                    suppressed += 1
                else:
                    successful += 1
                    retrieved = [idx for idx, _ in result]

                    for k_ret in self.k_values:
                        recall = RetrievalMetrics.recall_at_k(
                            retrieved, relevant, k_ret
                        )
                        recall_sum[k_ret] += recall

            n_queries = len(self.phenopackets)
            success_rate = successful / n_queries if n_queries > 0 else 0

            for k_ret in self.k_values:
                baseline_recall = baseline.get("recall", {}).get(k_ret, 0)
                conditional_recall = (
                    recall_sum[k_ret] / successful if successful > 0 else 0
                )

                results.append({
                    "k_anonymity": k,
                    "retrieval_k": k_ret,
                    "use_dp": use_dp,
                    "epsilon": epsilon if use_dp else None,
                    "success_rate": success_rate,
                    "suppression_rate": suppressed / n_queries if n_queries > 0 else 0,
                    "successful_queries": successful,
                    "suppressed_queries": suppressed,
                    "conditional_recall": conditional_recall,
                    "baseline_recall": baseline_recall,
                    "relative_recall": conditional_recall / baseline_recall if baseline_recall > 0 else 0,
                    "expected_recall": success_rate * conditional_recall
                })

        df = pd.DataFrame(results)
        self._results.extend(results)

        return df

    def sweep_rare_term_threshold(
        self,
        threshold_values: List[float],
        retrieval_k: int = 10
    ) -> pd.DataFrame:
        """
        Sweep rare term prevalence threshold.

        Args:
            threshold_values: List of min_prevalence values to test
            retrieval_k: k for retrieval evaluation

        Returns:
            DataFrame with results
        """
        from ..privacy.privacy_calculator import PrivacyPreservingCalculator
        from ..privacy.k_anonymity import PrivacyConfig

        baseline = self.compute_baseline()

        results = []

        for threshold in threshold_values:
            logger.info(f"Testing rare term threshold = {threshold}")

            config = PrivacyConfig(
                min_prevalence=threshold,
                k=1  # No k-anonymity suppression
            )

            base_metric = None
            if hasattr(self.base_calculator, 'similarity_metric'):
                base_metric = self.base_calculator.similarity_metric

            private_calc = PrivacyPreservingCalculator(
                base_similarity=base_metric,
                config=config,
                hpo_manager=self.hpo_manager,
                use_psi=False,
                use_dp=False,
                use_k_anonymity=False,
                use_rare_filter=True
            )

            # Set corpus for prevalence computation
            private_calc.set_corpus(self.phenopackets)

            # Evaluate
            evaluator = PhenopacketEvaluator(private_calc, self.k_values)
            metrics = evaluator.evaluate_retrieval(
                self.phenopackets,
                self.phenopackets,
                self._ground_truth
            )

            # Get filtering statistics
            filter_stats = private_calc.rare_filter.get_statistics()

            for k in self.k_values:
                baseline_recall = baseline.get("recall", {}).get(k, 0)

                results.append({
                    "min_prevalence": threshold,
                    "retrieval_k": k,
                    "recall": metrics.get("recall", {}).get(k, 0),
                    "precision": metrics.get("precision", {}).get(k, 0),
                    "baseline_recall": baseline_recall,
                    "relative_recall": (
                        metrics.get("recall", {}).get(k, 0) / baseline_recall
                        if baseline_recall > 0 else 0
                    ),
                    "terms_filtered": filter_stats.get("total_filtered", 0)
                })

        return pd.DataFrame(results)

    def combined_sweep(
        self,
        epsilon_values: List[float],
        k_anonymity_values: List[int],
        n_trials: int = 5
    ) -> pd.DataFrame:
        """
        Combined sweep of epsilon and k-anonymity.

        Args:
            epsilon_values: Epsilon values to test
            k_anonymity_values: k-anonymity k values to test
            n_trials: Trials per combination

        Returns:
            DataFrame with results
        """
        from ..privacy.privacy_calculator import PrivacyPreservingCalculator
        from ..privacy.k_anonymity import PrivacyConfig

        baseline = self.compute_baseline()
        results = []

        total_combos = len(epsilon_values) * len(k_anonymity_values) * n_trials
        combo_num = 0

        for epsilon in epsilon_values:
            for k_anon in k_anonymity_values:
                for trial in range(n_trials):
                    combo_num += 1
                    if combo_num % 10 == 0:
                        logger.info(f"Progress: {combo_num}/{total_combos}")

                    config = PrivacyConfig(
                        epsilon=epsilon,
                        k=k_anon
                    )

                    base_metric = None
                    if hasattr(self.base_calculator, 'similarity_metric'):
                        base_metric = self.base_calculator.similarity_metric

                    private_calc = PrivacyPreservingCalculator(
                        base_similarity=base_metric,
                        config=config,
                        hpo_manager=self.hpo_manager,
                        use_psi=False,
                        use_dp=True,
                        use_k_anonymity=True,
                        use_rare_filter=False
                    )

                    evaluator = PhenopacketEvaluator(private_calc, self.k_values)
                    metrics = evaluator.evaluate_retrieval(
                        self.phenopackets,
                        self.phenopackets,
                        self._ground_truth
                    )

                    # Get suppression stats
                    k_stats = private_calc.k_guard.get_statistics()

                    for k in self.k_values:
                        baseline_recall = baseline.get("recall", {}).get(k, 0)

                        results.append({
                            "epsilon": epsilon,
                            "k_anonymity": k_anon,
                            "trial": trial,
                            "retrieval_k": k,
                            "recall": metrics.get("recall", {}).get(k, 0),
                            "precision": metrics.get("precision", {}).get(k, 0),
                            "ndcg": metrics.get("ndcg", {}).get(k, 0),
                            "baseline_recall": baseline_recall,
                            "relative_recall": (
                                metrics.get("recall", {}).get(k, 0) / baseline_recall
                                if baseline_recall > 0 else 0
                            ),
                            "suppression_rate": k_stats.get("suppression_rate", 0)
                        })

        return pd.DataFrame(results)

    def plot_frontier(
        self,
        results: pd.DataFrame,
        output_path: str,
        x_param: str = "epsilon",
        y_metric: str = "recall"
    ):
        """
        Generate privacy-utility frontier plot.

        Args:
            results: DataFrame from sweep methods
            output_path: Path to save figure
            x_param: Privacy parameter for x-axis
            y_metric: Utility metric for y-axis
        """
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
        except ImportError:
            logger.warning("matplotlib/seaborn not available; skipping plot")
            return

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Plot 1: Absolute utility
        ax1 = axes[0]
        for k in results["k"].unique() if "k" in results.columns else results["retrieval_k"].unique():
            k_col = "k" if "k" in results.columns else "retrieval_k"
            subset = results[results[k_col] == k]

            if "trial" in subset.columns:
                grouped = subset.groupby(x_param)[y_metric].agg(["mean", "std"])
                ax1.errorbar(
                    grouped.index,
                    grouped["mean"],
                    yerr=grouped["std"],
                    marker='o',
                    capsize=3,
                    label=f"k={k}"
                )
            else:
                ax1.plot(
                    subset[x_param],
                    subset[y_metric],
                    marker='o',
                    label=f"k={k}"
                )

        ax1.set_xlabel(f"Privacy Parameter ({x_param})")
        ax1.set_ylabel(f"{y_metric.capitalize()}")
        if x_param == "epsilon":
            ax1.set_xscale("log")
        ax1.legend()
        ax1.set_title(f"{y_metric.capitalize()} vs {x_param.capitalize()}")
        ax1.grid(True, alpha=0.3)

        # Plot 2: Relative utility
        ax2 = axes[1]
        rel_col = f"relative_{y_metric}"

        if rel_col in results.columns:
            for k in results["k"].unique() if "k" in results.columns else results["retrieval_k"].unique():
                k_col = "k" if "k" in results.columns else "retrieval_k"
                subset = results[results[k_col] == k]

                if "trial" in subset.columns:
                    grouped = subset.groupby(x_param)[rel_col].agg(["mean", "std"])
                    ax2.errorbar(
                        grouped.index,
                        grouped["mean"],
                        yerr=grouped["std"],
                        marker='s',
                        capsize=3,
                        label=f"k={k}"
                    )
                else:
                    ax2.plot(
                        subset[x_param],
                        subset[rel_col],
                        marker='s',
                        label=f"k={k}"
                    )

            ax2.axhline(y=1.0, linestyle='--', color='gray', alpha=0.7, label='Baseline')

        ax2.set_xlabel(f"Privacy Parameter ({x_param})")
        ax2.set_ylabel(f"Relative {y_metric.capitalize()} (vs Baseline)")
        if x_param == "epsilon":
            ax2.set_xscale("log")
        ax2.legend()
        ax2.set_title("Relative Utility Degradation")
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        logger.info(f"Saved frontier plot to {output_path}")

    def generate_report(
        self,
        output_dir: str
    ) -> Dict[str, Any]:
        """
        Generate comprehensive analysis report.

        Args:
            output_dir: Directory to save outputs

        Returns:
            Summary dictionary
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # Run sweeps if not already done
        logger.info("Running epsilon sweep...")
        epsilon_results = self.sweep_epsilon(
            [0.1, 0.5, 1.0, 2.0, 5.0, 10.0],
            n_trials=5
        )
        epsilon_results.to_csv(output_path / "epsilon_sweep.csv", index=False)

        logger.info("Running k-anonymity sweep...")
        k_anon_results = self.sweep_k_anonymity([2, 3, 5, 10, 20])
        k_anon_results.to_csv(output_path / "k_anonymity_sweep.csv", index=False)

        # Generate plots
        self.plot_frontier(
            epsilon_results,
            str(output_path / "epsilon_frontier.png"),
            x_param="epsilon"
        )

        # Summary statistics
        baseline = self.compute_baseline()

        summary = {
            "n_phenopackets": len(self.phenopackets),
            "baseline_metrics": baseline,
            "epsilon_sweep": {
                "best_epsilon": epsilon_results.groupby("epsilon")["recall"].mean().idxmax(),
                "mean_relative_recall_at_eps1": epsilon_results[
                    epsilon_results["epsilon"] == 1.0
                ]["relative_recall"].mean()
            },
            "k_anonymity_sweep": {
                "suppression_rates": k_anon_results.groupby("k_anonymity")[
                    "suppression_rate"
                ].mean().to_dict()
            }
        }

        return summary


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # This module requires actual phenopackets and calculator to run
    print("PrivacyUtilityFrontier module loaded successfully")
    print("Use with PhenopacketSimilarityCalculator and phenopacket data")
