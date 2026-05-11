#!/usr/bin/env python3
"""
Rank-based vs. score-based differential privacy on the real Phenopacket Store cohort.

Motivated by the synthetic-to-real gap (§5.1.5, §4.6): per-score Laplace noise
collapses to random under ε ≤ 10 on real patients because same-disease score
gaps are compressed into 0.2–0.3 while noise scale is 1/ε. We compare three
top-k release mechanisms at fixed ε:

  LAPLACE     (baseline): add Lap(1/ε) to each similarity, return top-k by
              noisy score. Releases the noisy score vector; ε-DP w.r.t. any
              single record (sensitivity 1).
  EXP-SCORE   Exponential mechanism over candidates with utility = similarity
              score, sensitivity 1. Budget ε split across the k iterative
              selections (ε/k per round). Same "compressed signal" pathology
              as Laplace because the utility range is the score range.
  EXP-RANK    Exponential mechanism over candidates with utility = -rank,
              sensitivity 1 (one record change shifts any other patient's rank
              by at most 1). Rank gaps are large by construction
              (rank 1 vs rank 1000 differ by 999) regardless of score
              distribution, so the same ε buys far stronger top-k separation.

Outputs:
    results/rank_based_dp/sweep.csv
    figures/fig14_rank_based_dp.{png,pdf}
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.data_integration.phenopacket_store_loader import CohortSpec, load_phenopacket_store
from src.similarity.hpo_similarity import (
    CosineSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("rank_dp")

REPO = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO / "results" / "rank_based_dp"
FIG_DIR = REPO / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Mechanisms
# ---------------------------------------------------------------------------

def topk_laplace(true_sims: np.ndarray, epsilon: float, k: int, rng: np.random.Generator) -> List[int]:
    """Laplace per-score, sort, take top-k. Noise scale 1/ε; sensitivity 1."""
    if not np.isfinite(epsilon):
        return np.argsort(true_sims)[::-1][:k].tolist()
    noise = rng.laplace(0.0, 1.0 / epsilon, size=true_sims.shape)
    noisy = true_sims + noise
    return np.argsort(noisy)[::-1][:k].tolist()


def topk_exp_score(true_sims: np.ndarray, epsilon: float, k: int, rng: np.random.Generator) -> List[int]:
    """Iterative exponential mechanism with score utility; budget ε split as ε/k.

    Sensitivity of similarity (∈ [0, 1]) is 1.
    """
    if not np.isfinite(epsilon):
        return np.argsort(true_sims)[::-1][:k].tolist()
    eps_round = epsilon / k
    remaining = list(range(len(true_sims)))
    selected: List[int] = []
    for _ in range(min(k, len(remaining))):
        utilities = true_sims[remaining]
        # Probability ∝ exp(ε_round · u / (2 · Δu)); Δu = 1.
        logits = (eps_round / 2.0) * utilities
        logits = logits - logits.max()  # numerical stability
        probs = np.exp(logits)
        probs /= probs.sum()
        idx_in_remaining = int(rng.choice(len(remaining), p=probs))
        selected.append(remaining.pop(idx_in_remaining))
    return selected


def topk_exp_rank(true_sims: np.ndarray, epsilon: float, k: int, rng: np.random.Generator) -> List[int]:
    """Iterative exponential mechanism with rank utility; budget ε split as ε/k.

    Sensitivity of the rank function is 1: adding/removing one record shifts
    any other record's rank by at most 1. Utility = -rank (higher = more
    similar). Rank gaps are O(n) by construction, so the mechanism recovers
    sharp top-k selection even when similarity scores are compressed.
    """
    if not np.isfinite(epsilon):
        return np.argsort(true_sims)[::-1][:k].tolist()
    # Compute ranks once: rank[i] is the position of patient i (1-indexed),
    # with rank 1 = most similar.
    order = np.argsort(true_sims)[::-1]
    rank = np.empty_like(order)
    rank[order] = np.arange(1, len(order) + 1)
    utilities_all = -rank.astype(float)  # higher = more similar
    eps_round = epsilon / k
    remaining = list(range(len(true_sims)))
    selected: List[int] = []
    for _ in range(min(k, len(remaining))):
        utilities = utilities_all[remaining]
        logits = (eps_round / 2.0) * utilities  # Δu = 1
        logits = logits - logits.max()
        probs = np.exp(logits)
        probs /= probs.sum()
        idx_in_remaining = int(rng.choice(len(remaining), p=probs))
        selected.append(remaining.pop(idx_in_remaining))
    return selected


MECHANISMS = {
    "Laplace (per-score)": topk_laplace,
    "Exp. mech. (score)":  topk_exp_score,
    "Exp. mech. (rank)":   topk_exp_rank,
}


# ---------------------------------------------------------------------------
# Retrieval evaluation
# ---------------------------------------------------------------------------

def evaluate_mechanism(
    sim_matrix: np.ndarray,
    ground_truth: List[str],
    mechanism,
    epsilon: float,
    k: int,
    rng: np.random.Generator,
) -> dict:
    """Leave-one-out retrieval using the supplied mechanism."""
    n = sim_matrix.shape[0]
    mrrs, ndcgs, recall5, recall10, p1 = [], [], [], [], []
    for q in range(n):
        # candidates = everything except the query itself
        true_sims = np.delete(sim_matrix[q], q)
        candidates = [i for i in range(n) if i != q]
        ranked = mechanism(true_sims, epsilon=epsilon, k=10, rng=rng)
        ranked_global = [candidates[i] for i in ranked]
        q_disease = ground_truth[q]
        rel = np.array([1 if ground_truth[ix] == q_disease else 0 for ix in ranked_global])
        relevant_total = sum(1 for j in candidates if ground_truth[j] == q_disease)
        if relevant_total == 0:
            continue
        # MRR
        if rel.any():
            mrrs.append(1.0 / (np.argmax(rel) + 1))
        else:
            mrrs.append(0.0)
        # nDCG@10
        gains = rel / np.log2(np.arange(2, len(rel) + 2))
        ideal = np.zeros_like(rel)
        ideal[: min(relevant_total, len(rel))] = 1
        idcg = (ideal / np.log2(np.arange(2, len(rel) + 2))).sum()
        ndcgs.append(gains.sum() / idcg if idcg else 0.0)
        recall5.append(rel[:5].sum() / min(relevant_total, 5))
        recall10.append(rel[:10].sum() / min(relevant_total, 10))
        p1.append(int(rel[0]) if len(rel) else 0)
    return {
        "MRR": float(np.mean(mrrs)),
        "nDCG@10": float(np.mean(ndcgs)),
        "Recall@5": float(np.mean(recall5)),
        "Recall@10": float(np.mean(recall10)),
        "P@1": float(np.mean(p1)),
        "n_queries": len(mrrs),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n-diseases", type=int, default=100)
    ap.add_argument("--max-per-disease", type=int, default=15)
    ap.add_argument("--min-phenotypes", type=int, default=3)
    ap.add_argument("--epsilons", nargs="+", type=float,
                    default=[0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, float("inf")])
    ap.add_argument("--top-k", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    spec = CohortSpec(
        min_phenotypes=args.min_phenotypes,
        min_cohort_size=2,
        max_per_disease=args.max_per_disease,
        top_n_diseases=args.top_n_diseases,
        seed=args.seed,
    )
    pps, gt_map, stats = load_phenopacket_store(spec=spec)
    log.info("Cohort: %d patients, %d diseases", stats["n_after_filter"], stats["n_diseases"])

    # Precompute the n×n similarity matrix once.
    ic = compute_empirical_ic(pps)
    metric = CosineSimilarity(ic_values=ic)
    calc = PhenopacketSimilarityCalculator(metric, cache_similarities=True)
    log.info("Computing %dx%d similarity matrix…", len(pps), len(pps))
    sim_matrix = calc.compute_similarity_matrix(pps)
    log.info("Similarity matrix computed (mean=%.3f, max=%.3f, nonzero=%d%%)",
             float(sim_matrix.mean()),
             float(sim_matrix.max()),
             int(100 * (sim_matrix > 0).mean()))

    ground_truth = [gt_map[pp["id"]] for pp in pps]
    rng = np.random.default_rng(args.seed)

    rows = []
    for mech_name, mech in MECHANISMS.items():
        for eps in args.epsilons:
            log.info("  %-22s ε=%s", mech_name, "∞" if not np.isfinite(eps) else f"{eps:g}")
            metrics = evaluate_mechanism(sim_matrix, ground_truth, mech,
                                         epsilon=eps, k=args.top_k, rng=rng)
            rows.append({"mechanism": mech_name, "epsilon": eps, **metrics})
            log.info("    %s", {k: round(v, 4) for k, v in metrics.items() if isinstance(v, float)})

    df = pd.DataFrame(rows)
    out = RESULTS_DIR / "sweep.csv"
    df.to_csv(out, index=False)
    log.info("Wrote %s", out)

    # ---- Plot ----
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = {"Laplace (per-score)": "C0", "Exp. mech. (score)": "C1", "Exp. mech. (rank)": "C2"}
    markers = {"Laplace (per-score)": "o", "Exp. mech. (score)": "s", "Exp. mech. (rank)": "^"}
    for mech in MECHANISMS:
        sub = df[df["mechanism"] == mech].sort_values("epsilon")
        finite = sub[np.isfinite(sub["epsilon"])]
        inf_row = sub[~np.isfinite(sub["epsilon"])]
        c, m = colors[mech], markers[mech]
        ax1.plot(finite["epsilon"], finite["nDCG@10"], m + "-", color=c, label=mech)
        ax2.plot(finite["epsilon"], finite["MRR"], m + "-", color=c, label=mech)
        if not inf_row.empty:
            x_inf = finite["epsilon"].max() * 3
            ax1.scatter([x_inf], inf_row["nDCG@10"], marker=m, facecolors="none", edgecolors=c, s=70)
            ax2.scatter([x_inf], inf_row["MRR"], marker=m, facecolors="none", edgecolors=c, s=70)
    for ax in (ax1, ax2):
        ax.set_xscale("log")
        ax.set_xlabel("Privacy budget ε (log scale)")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="lower right", fontsize=9)
    ax1.set_ylabel("nDCG@10")
    ax2.set_ylabel("MRR")
    ax1.set_ylim(0, 1.02)
    ax2.set_ylim(0, 1.02)
    fig.suptitle(
        "Rank-based DP recovers retrieval utility on a real-patient cohort\n"
        f"({stats['n_after_filter']} Phenopacket Store patients across {stats['n_diseases']} OMIM diseases)",
        y=1.02,
    )
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig14_rank_based_dp.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)
    log.info("Wrote fig14_rank_based_dp.{png,pdf}")


if __name__ == "__main__":
    main()
