#!/usr/bin/env python3
"""
Retrieval benchmark on the Monarch Phenopacket Store (real published case-report
phenopackets with confirmed OMIM diagnoses).

We evaluate two non-private retrieval baselines on the same cohort:

  COSINE-IC   IC-weighted cosine similarity (our system's baseline metric)
  RESNIK-BMA  Resnik semantic similarity with best-match-average aggregation
              over the full HPO DAG (Phenomizer-style; Köhler et al. 2009)

and then sweep the Laplace differential-privacy budget ε ∈ {0.5, 1, 2, 5, ∞}
for the cosine-IC pipeline to characterise the privacy cost on a real cohort.

Outputs:
    results/phenopacket_store/baselines.csv
    results/phenopacket_store/dp_sweep.csv
    figures/fig12_phenopacket_benchmark.{png,pdf}
"""

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.data_integration.phenopacket_store_loader import CohortSpec, load_phenopacket_store
from src.privacy.differential_privacy import LaplaceMechanism, PrivateSimilarityCalculator
from src.similarity.graph_similarity import GraphAwarePhenopacketCalculator, OntologyGraphSimilarity
from src.similarity.hpo_similarity import (
    CosineSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
)
from src.utils.hpo_utils import HPOManager

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger("ppstore_benchmark")

REPO = Path(__file__).resolve().parent.parent
OBO_PATH = REPO / "data" / "hpo_ontology" / "hp.obo"
RESULTS_DIR = REPO / "results" / "phenopacket_store"
FIG_DIR = REPO / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Retrieval metrics
# ---------------------------------------------------------------------------

def retrieval_metrics(calc, phenopackets, ground_truth, top_k=10):
    """Compute mean MRR, nDCG@k, Recall@k over leave-one-out queries.

    For each query patient, the relevant set is all other patients sharing
    the same OMIM disease ID. The candidate database is the cohort minus
    the query.
    """
    n = len(phenopackets)
    mrrs, ndcgs, recalls_5, recalls_10, prec_at_1 = [], [], [], [], []
    for q_idx, q in enumerate(phenopackets):
        q_disease = ground_truth[q["id"]]
        candidates = [p for j, p in enumerate(phenopackets) if j != q_idx]
        cand_disease = [ground_truth[p["id"]] for p in candidates]
        relevant_total = sum(1 for d in cand_disease if d == q_disease)
        if relevant_total == 0:
            continue
        ranked = calc.find_most_similar(q, candidates, top_k=top_k)
        rel = np.array([1 if cand_disease[i] == q_disease else 0 for i, _ in ranked])
        # MRR
        first_hit = np.argmax(rel) if rel.any() else -1
        mrrs.append(1.0 / (first_hit + 1) if rel.any() else 0.0)
        # nDCG@k
        gains = rel / np.log2(np.arange(2, len(rel) + 2))
        ideal_rel = np.zeros_like(rel)
        ideal_rel[: min(relevant_total, len(rel))] = 1
        ideal = ideal_rel / np.log2(np.arange(2, len(rel) + 2))
        ndcgs.append(gains.sum() / ideal.sum() if ideal.sum() else 0.0)
        # Recall@5, @10
        recalls_5.append(rel[:5].sum() / min(relevant_total, 5))
        recalls_10.append(rel[:10].sum() / min(relevant_total, 10))
        prec_at_1.append(int(rel[0]) if len(rel) else 0)
    return {
        "n_queries": len(mrrs),
        "MRR": float(np.mean(mrrs)),
        "nDCG@10": float(np.mean(ndcgs)),
        "Recall@5": float(np.mean(recalls_5)),
        "Recall@10": float(np.mean(recalls_10)),
        "P@1": float(np.mean(prec_at_1)),
    }


# ---------------------------------------------------------------------------
# Calculators
# ---------------------------------------------------------------------------

def build_cosine_calc(phenopackets):
    ic = compute_empirical_ic(phenopackets)
    metric = CosineSimilarity(ic_values=ic)
    return PhenopacketSimilarityCalculator(metric, cache_similarities=True), ic


def build_resnik_calc(phenopackets, obo_path):
    log.info("Loading HPO ontology from %s", obo_path)
    mgr = HPOManager(str(obo_path))
    mgr.load_ontology()  # eager: OntologyGraphSimilarity inspects mgr.ontology in __init__
    graph = OntologyGraphSimilarity(hpo_manager=mgr)
    log.info("Computing corpus IC with ancestor propagation (%d phenopackets)", len(phenopackets))
    graph.set_corpus_ic(phenopackets)
    return GraphAwarePhenopacketCalculator(graph, method="resnik", aggregation="bma")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n-diseases", type=int, default=50)
    ap.add_argument("--max-per-disease", type=int, default=10)
    ap.add_argument("--min-phenotypes", type=int, default=3)
    ap.add_argument("--epsilons", nargs="+", type=float, default=[0.5, 1.0, 2.0, 5.0, float("inf")])
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--skip-resnik", action="store_true",
                    help="Skip the (slower) Resnik+BMA baseline.")
    args = ap.parse_args()

    spec = CohortSpec(
        min_phenotypes=args.min_phenotypes,
        min_cohort_size=2,
        max_per_disease=args.max_per_disease,
        top_n_diseases=args.top_n_diseases,
        seed=args.seed,
    )
    phenopackets, ground_truth, stats = load_phenopacket_store(spec=spec)
    log.info("Cohort: %d patients, %d diseases, mean %.1f patients/disease",
             stats["n_after_filter"], stats["n_diseases"], stats["mean_patients_per_disease"])

    rows = []

    # ---- Cosine-IC (non-private + DP sweep) -------------------------------
    cosine_calc, _ic = build_cosine_calc(phenopackets)

    log.info("Cosine-IC: non-private retrieval")
    m = retrieval_metrics(cosine_calc, phenopackets, ground_truth)
    rows.append({"system": "Cosine-IC", "epsilon": np.inf, **m})
    log.info("  %s", {k: round(v, 4) for k, v in m.items() if isinstance(v, float)})

    for eps in args.epsilons:
        if not np.isfinite(eps):
            continue  # ε=∞ already covered above
        log.info("Cosine-IC: ε = %.2f", eps)
        priv = PrivateSimilarityCalculator(
            base_calculator=cosine_calc,
            mechanism=LaplaceMechanism(epsilon=eps),
            sensitivity=1.0,
        )
        m = retrieval_metrics(priv, phenopackets, ground_truth)
        rows.append({"system": "Cosine-IC", "epsilon": eps, **m})
        log.info("  %s", {k: round(v, 4) for k, v in m.items() if isinstance(v, float)})

    # ---- Resnik+BMA (Phenomizer-style, non-private only) ------------------
    if not args.skip_resnik:
        resnik_calc = build_resnik_calc(phenopackets, OBO_PATH)
        log.info("Resnik+BMA (Phenomizer-style): non-private retrieval")
        m = retrieval_metrics(resnik_calc, phenopackets, ground_truth)
        rows.append({"system": "Resnik-BMA (Phenomizer-style)", "epsilon": np.inf, **m})
        log.info("  %s", {k: round(v, 4) for k, v in m.items() if isinstance(v, float)})

    df = pd.DataFrame(rows)
    out = RESULTS_DIR / "dp_sweep.csv"
    df.to_csv(out, index=False)
    log.info("Wrote %s", out)

    # ---- Plot -------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sweep = df[df["system"] == "Cosine-IC"].sort_values("epsilon")
    finite = sweep[np.isfinite(sweep["epsilon"])]
    inf_row = sweep[~np.isfinite(sweep["epsilon"])]
    ax.plot(finite["epsilon"], finite["nDCG@10"], "o-", label="Cosine-IC, nDCG@10")
    ax.plot(finite["epsilon"], finite["MRR"], "s--", label="Cosine-IC, MRR")
    if not inf_row.empty:
        x_inf = finite["epsilon"].max() * 3
        ax.scatter([x_inf], inf_row["nDCG@10"], marker="o", facecolors="none", edgecolors="C0", s=80)
        ax.scatter([x_inf], inf_row["MRR"], marker="s", facecolors="none", edgecolors="C1", s=80)
        ax.annotate("ε=∞", (x_inf, inf_row["nDCG@10"].iloc[0]), textcoords="offset points", xytext=(5, 5))
    if not args.skip_resnik:
        rb = df[df["system"].str.startswith("Resnik")].iloc[0]
        ax.axhline(rb["nDCG@10"], ls=":", color="C3",
                   label=f"Resnik-BMA non-private, nDCG@10 = {rb['nDCG@10']:.3f}")
    ax.set_xscale("log")
    ax.set_xlabel("Privacy budget ε (log scale)")
    ax.set_ylabel("Retrieval performance")
    ax.set_ylim(0, 1.02)
    ax.set_title("Retrieval on Phenopacket Store (real published cases)\n"
                 f"{stats['n_after_filter']} patients across {stats['n_diseases']} OMIM diseases")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fig12_phenopacket_benchmark.{ext}", dpi=180)
    plt.close(fig)
    log.info("Wrote fig12_phenopacket_benchmark.{png,pdf}")
    log.info("Done.")


if __name__ == "__main__":
    main()
