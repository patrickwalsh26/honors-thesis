#!/usr/bin/env python3
"""
Generate publication-quality figures for thesis.

Creates visualizations for:
- Privacy-utility tradeoff curves
- Leakage audit comparisons
- Parameter sensitivity analysis
- Combined mechanism heatmaps

Usage:
    python experiments/generate_figures.py --results experiments/results
"""

import argparse
import json
import logging
from pathlib import Path
import sys

import pandas as pd
import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("Warning: matplotlib/seaborn not available")

logger = logging.getLogger(__name__)

# Publication-quality settings
FIGURE_PARAMS = {
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'figure.figsize': (6, 4),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'grid.alpha': 0.3,
}


def setup_plotting():
    """Configure matplotlib for publication-quality figures."""
    if PLOTTING_AVAILABLE:
        plt.rcParams.update(FIGURE_PARAMS)
        sns.set_style("whitegrid")


def plot_epsilon_frontier(
    dp_results: pd.DataFrame,
    output_dir: Path
):
    """
    Generate DP epsilon vs utility frontier plots.

    Creates:
    - Recall vs epsilon curves
    - Relative utility degradation plot
    """
    if not PLOTTING_AVAILABLE:
        logger.warning("Plotting not available")
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Color palette
    colors = sns.color_palette("husl", n_colors=len(dp_results["k"].unique()))

    # Plot 1: Absolute Recall
    ax1 = axes[0]
    for i, k in enumerate(sorted(dp_results["k"].unique())):
        subset = dp_results[dp_results["k"] == k]
        grouped = subset.groupby("epsilon")["recall"].agg(["mean", "std"])

        ax1.errorbar(
            grouped.index,
            grouped["mean"],
            yerr=grouped["std"],
            marker='o',
            capsize=3,
            color=colors[i],
            label=f"Recall@{k}"
        )

    ax1.set_xlabel("Privacy Budget (ε)")
    ax1.set_ylabel("Recall")
    ax1.set_xscale("log")
    ax1.legend(loc="lower right")
    ax1.set_title("Retrieval Recall vs Privacy Budget")

    # Plot 2: Relative Utility
    ax2 = axes[1]
    for i, k in enumerate(sorted(dp_results["k"].unique())):
        subset = dp_results[dp_results["k"] == k]
        grouped = subset.groupby("epsilon")["relative_recall"].agg(["mean", "std"])

        ax2.errorbar(
            grouped.index,
            grouped["mean"],
            yerr=grouped["std"],
            marker='s',
            capsize=3,
            color=colors[i],
            label=f"k={k}"
        )

    ax2.axhline(y=1.0, linestyle='--', color='gray', alpha=0.7, linewidth=2)
    ax2.set_xlabel("Privacy Budget (ε)")
    ax2.set_ylabel("Relative Recall (vs Non-Private)")
    ax2.set_xscale("log")
    ax2.legend(loc="lower right")
    ax2.set_title("Utility Degradation")

    plt.tight_layout()
    plt.savefig(output_dir / "epsilon_frontier.pdf")
    plt.savefig(output_dir / "epsilon_frontier.png", dpi=150)
    plt.close()

    logger.info(f"Saved epsilon frontier plot to {output_dir}")


def plot_k_anonymity_tradeoff(
    k_anon_results: pd.DataFrame,
    output_dir: Path
):
    """
    Generate k-anonymity tradeoff plots.

    Shows:
    - Query success rate vs k
    - Conditional recall vs k
    """
    if not PLOTTING_AVAILABLE:
        return

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Filter to k=10 for cleaner visualization
    subset = k_anon_results[k_anon_results["retrieval_k"] == 10]

    # Plot 1: Success Rate
    ax1 = axes[0]
    ax1.bar(
        range(len(subset)),
        subset["success_rate"],
        tick_label=[str(k) for k in subset["k_anonymity"]],
        color=sns.color_palette("Blues", n_colors=len(subset))
    )
    ax1.set_xlabel("k-Anonymity Parameter")
    ax1.set_ylabel("Query Success Rate")
    ax1.set_title("Query Success Rate vs k-Anonymity")
    ax1.set_ylim(0, 1.05)

    # Add suppression rate as text
    for i, (_, row) in enumerate(subset.iterrows()):
        ax1.text(
            i, row["success_rate"] + 0.02,
            f"{row['suppression_rate']:.0%}\nsuppressed",
            ha='center', va='bottom', fontsize=8
        )

    # Plot 2: Conditional Recall
    ax2 = axes[1]
    x = range(len(subset))
    width = 0.35

    ax2.bar(
        [i - width/2 for i in x],
        subset["conditional_recall"],
        width,
        label="Conditional Recall",
        color=sns.color_palette("Greens")[3]
    )
    ax2.bar(
        [i + width/2 for i in x],
        subset["baseline_recall"],
        width,
        label="Baseline Recall",
        color=sns.color_palette("Grays")[3],
        alpha=0.7
    )

    ax2.set_xlabel("k-Anonymity Parameter")
    ax2.set_ylabel("Recall@10")
    ax2.set_xticks(x)
    ax2.set_xticklabels([str(k) for k in subset["k_anonymity"]])
    ax2.legend()
    ax2.set_title("Conditional Recall (Successful Queries)")

    plt.tight_layout()
    plt.savefig(output_dir / "k_anonymity_tradeoff.pdf")
    plt.savefig(output_dir / "k_anonymity_tradeoff.png", dpi=150)
    plt.close()

    logger.info(f"Saved k-anonymity plot to {output_dir}")


def plot_leakage_comparison(
    leakage_results: dict,
    output_dir: Path
):
    """
    Generate leakage audit comparison plots.

    Compares attack success between baseline and private systems.
    """
    if not PLOTTING_AVAILABLE:
        return

    if "error" in leakage_results:
        logger.warning("Skipping leakage plot due to audit errors")
        return

    mia = leakage_results.get("membership_inference", {})

    if "baseline" not in mia or "private" not in mia:
        logger.warning("Missing MIA results for plotting")
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    metrics = ["accuracy", "precision", "recall", "auc"]
    x = np.arange(len(metrics))
    width = 0.35

    baseline_vals = [mia["baseline"].get(m, 0) for m in metrics]
    private_vals = [mia["private"].get(m, 0) for m in metrics]

    bars1 = ax.bar(
        x - width/2,
        baseline_vals,
        width,
        label="No Privacy",
        color=sns.color_palette("Reds")[3]
    )
    bars2 = ax.bar(
        x + width/2,
        private_vals,
        width,
        label="With Privacy",
        color=sns.color_palette("Blues")[3]
    )

    ax.axhline(y=0.5, linestyle='--', color='gray', alpha=0.7, label="Random Guess")

    ax.set_xlabel("Attack Metric")
    ax.set_ylabel("Score")
    ax.set_xticks(x)
    ax.set_xticklabels([m.capitalize() for m in metrics])
    ax.legend()
    ax.set_title("Membership Inference Attack Success")
    ax.set_ylim(0, 1.0)

    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width()/2,
                height + 0.02,
                f'{height:.2f}',
                ha='center', va='bottom', fontsize=8
            )

    plt.tight_layout()
    plt.savefig(output_dir / "leakage_comparison.pdf")
    plt.savefig(output_dir / "leakage_comparison.png", dpi=150)
    plt.close()

    logger.info(f"Saved leakage comparison plot to {output_dir}")


def plot_combined_heatmap(
    combined_results: pd.DataFrame,
    output_dir: Path
):
    """
    Generate combined mechanism heatmap.

    Shows utility across epsilon x k-anonymity grid.
    """
    if not PLOTTING_AVAILABLE:
        return

    # Filter to retrieval_k=10
    subset = combined_results[combined_results["retrieval_k"] == 10]

    # Pivot for heatmap
    pivot = subset.groupby(["epsilon", "k_anonymity"])["relative_recall"].mean().unstack()

    fig, ax = plt.subplots(figsize=(8, 6))

    sns.heatmap(
        pivot,
        annot=True,
        fmt=".2f",
        cmap="RdYlGn",
        vmin=0,
        vmax=1,
        ax=ax,
        cbar_kws={"label": "Relative Recall"}
    )

    ax.set_xlabel("k-Anonymity Parameter")
    ax.set_ylabel("Privacy Budget (ε)")
    ax.set_title("Combined Privacy Mechanisms: Utility Heatmap")

    plt.tight_layout()
    plt.savefig(output_dir / "combined_heatmap.pdf")
    plt.savefig(output_dir / "combined_heatmap.png", dpi=150)
    plt.close()

    logger.info(f"Saved combined heatmap to {output_dir}")


def plot_baseline_comparison(
    baseline_results: pd.DataFrame,
    output_dir: Path
):
    """
    Generate baseline metric comparison plot.
    """
    if not PLOTTING_AVAILABLE:
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    metrics = baseline_results["metric"].unique()
    k_values = sorted(baseline_results["k"].unique())

    x = np.arange(len(k_values))
    width = 0.25

    colors = sns.color_palette("husl", n_colors=len(metrics))

    for i, metric in enumerate(metrics):
        subset = baseline_results[baseline_results["metric"] == metric]
        recalls = [subset[subset["k"] == k]["recall"].values[0] for k in k_values]

        ax.bar(
            x + i * width,
            recalls,
            width,
            label=metric.replace("_", " ").title(),
            color=colors[i]
        )

    ax.set_xlabel("k (Top-k Retrieval)")
    ax.set_ylabel("Recall@k")
    ax.set_xticks(x + width)
    ax.set_xticklabels([f"k={k}" for k in k_values])
    ax.legend()
    ax.set_title("Baseline Similarity Metrics Comparison")

    plt.tight_layout()
    plt.savefig(output_dir / "baseline_comparison.pdf")
    plt.savefig(output_dir / "baseline_comparison.png", dpi=150)
    plt.close()

    logger.info(f"Saved baseline comparison plot to {output_dir}")


def generate_summary_table(
    results_dir: Path,
    output_dir: Path
):
    """
    Generate summary statistics table.
    """
    summary = []

    # Baseline best
    baseline_path = results_dir / "baseline_results.csv"
    if baseline_path.exists():
        df = pd.read_csv(baseline_path)
        best_metric = df.loc[df["recall"].idxmax()]
        summary.append({
            "Experiment": "Baseline (Best)",
            "Metric": best_metric["metric"],
            "Recall@10": f"{best_metric['recall']:.4f}",
            "Notes": "No privacy mechanisms"
        })

    # DP at epsilon=1
    dp_path = results_dir / "dp_sweep_results.csv"
    if dp_path.exists():
        df = pd.read_csv(dp_path)
        eps1 = df[(df["epsilon"] == 1.0) & (df["k"] == 10)]
        summary.append({
            "Experiment": "DP (ε=1.0)",
            "Metric": "Cosine IC",
            "Recall@10": f"{eps1['recall'].mean():.4f} ± {eps1['recall'].std():.4f}",
            "Notes": f"Relative: {eps1['relative_recall'].mean():.2%}"
        })

    # k-anonymity at k=5
    k_anon_path = results_dir / "k_anonymity_results.csv"
    if k_anon_path.exists():
        df = pd.read_csv(k_anon_path)
        k5 = df[(df["k_anonymity"] == 5) & (df["retrieval_k"] == 10)]
        if len(k5) > 0:
            summary.append({
                "Experiment": "k-Anonymity (k=5)",
                "Metric": "Cosine IC",
                "Recall@10": f"{k5['conditional_recall'].values[0]:.4f}",
                "Notes": f"Suppression: {k5['suppression_rate'].values[0]:.1%}"
            })

    # Create summary dataframe
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_dir / "summary_table.csv", index=False)

    # Also save as markdown
    with open(output_dir / "summary_table.md", "w") as f:
        f.write("# Experimental Results Summary\n\n")
        f.write(summary_df.to_markdown(index=False))

    logger.info(f"Saved summary table to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate thesis figures from experiment results"
    )
    parser.add_argument(
        "--results",
        default="experiments/results",
        help="Directory containing experiment results"
    )
    parser.add_argument(
        "--output",
        default="figures",
        help="Output directory for figures"
    )
    parser.add_argument(
        "--format",
        nargs="+",
        default=["pdf", "png"],
        choices=["pdf", "png", "svg"],
        help="Output formats"
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    if not PLOTTING_AVAILABLE:
        logger.error("matplotlib/seaborn required for figure generation")
        logger.error("Install with: pip install matplotlib seaborn")
        sys.exit(1)

    results_dir = Path(args.results)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    setup_plotting()

    # Generate figures based on available results
    if (results_dir / "dp_sweep_results.csv").exists():
        dp_results = pd.read_csv(results_dir / "dp_sweep_results.csv")
        plot_epsilon_frontier(dp_results, output_dir)

    if (results_dir / "k_anonymity_results.csv").exists():
        k_anon_results = pd.read_csv(results_dir / "k_anonymity_results.csv")
        plot_k_anonymity_tradeoff(k_anon_results, output_dir)

    if (results_dir / "leakage_audit.json").exists():
        with open(results_dir / "leakage_audit.json") as f:
            leakage_results = json.load(f)
        plot_leakage_comparison(leakage_results, output_dir)

    if (results_dir / "combined_results.csv").exists():
        combined_results = pd.read_csv(results_dir / "combined_results.csv")
        plot_combined_heatmap(combined_results, output_dir)

    if (results_dir / "baseline_results.csv").exists():
        baseline_results = pd.read_csv(results_dir / "baseline_results.csv")
        plot_baseline_comparison(baseline_results, output_dir)

    # Generate summary table
    generate_summary_table(results_dir, output_dir)

    logger.info(f"All figures saved to {output_dir}")


if __name__ == "__main__":
    main()
