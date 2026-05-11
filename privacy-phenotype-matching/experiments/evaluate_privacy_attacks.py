#!/usr/bin/env python3
"""
Membership inference + k-anonymity ablation against the HPOA evaluation cohort.

Outputs:
    results/privacy_attacks/mia_epsilon_sweep.csv
    results/privacy_attacks/kanon_ablation.csv
    figures/fig10_mia_pareto.{png,pdf}
    figures/fig11_kanon_ablation.{png,pdf}
"""

import argparse
import json
import logging
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.privacy.differential_privacy import (
    LaplaceMechanism,
    PrivateSimilarityCalculator,
)
from src.similarity.hpo_similarity import (
    CosineSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(message)s",
)
log = logging.getLogger("privacy_attacks")

REPO = Path(__file__).resolve().parent.parent
COHORT_PATH = REPO / "data" / "hpoa_evaluation" / "cohort_phenopackets.json"
RESULTS_DIR = REPO / "results" / "privacy_attacks"
FIG_DIR = REPO / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Setup helpers
# ---------------------------------------------------------------------------

def load_cohort():
    log.info("Loading cohort from %s", COHORT_PATH)
    with open(COHORT_PATH) as f:
        return json.load(f)


def build_base_calculator(phenopackets):
    ic_values = compute_empirical_ic(phenopackets)
    metric = CosineSimilarity(ic_values=ic_values)
    return PhenopacketSimilarityCalculator(metric, cache_similarities=True), ic_values


def patient_terms(pp):
    return {f["type"]["id"] for f in pp.get("phenotypicFeatures", []) if not f.get("excluded", False)}


# ---------------------------------------------------------------------------
# Membership-inference attacks
# ---------------------------------------------------------------------------
# Two attackers, both following the standard MI literature:
#
#   THRESHOLD: simplest possible — query the system, take the max similarity
#              over the database, threshold to predict membership. ROC AUC
#              is computed by sweeping the threshold (Yeom et al. 2018).
#
#   SHADOW:    shadow-model attack (Shokri et al. 2017): train shadow models
#              on random subsets, extract statistical features of returned
#              similarity distributions, learn a binary classifier.

def query_features(calc, target, database, top_k=10):
    matches = calc.find_most_similar(target, database, top_k=top_k)
    if not matches:
        return np.zeros(10)
    sims = np.array([s for _, s in matches])
    return np.array(
        [
            sims.max(),
            sims.mean(),
            sims.min(),
            sims.std() if len(sims) > 1 else 0.0,
            np.median(sims),
            sims[0] - sims[-1] if len(sims) > 1 else 0.0,
            (sims > 0.5).sum(),
            (sims > 0.3).sum(),
            np.percentile(sims, 75) if len(sims) >= 4 else sims.max(),
            np.percentile(sims, 25) if len(sims) >= 4 else sims.min(),
        ]
    )


def threshold_attack_auc(calc, phenopackets, member_idx, nonmember_idx, top_k=5, rng=None):
    """Yeom-style threshold attack: max-similarity is the attack score.

    Members are queried against a database that contains them; non-members are
    queried against the same database. AUC measures how separable the two
    distributions are.
    """
    member_db = [phenopackets[i] for i in member_idx]
    scores_in = []
    for i in member_idx:
        m = calc.find_most_similar(phenopackets[i], member_db, top_k=top_k)
        scores_in.append(max((s for _, s in m), default=0.0))
    scores_out = []
    for i in nonmember_idx:
        m = calc.find_most_similar(phenopackets[i], member_db, top_k=top_k)
        scores_out.append(max((s for _, s in m), default=0.0))
    y = np.concatenate([np.ones(len(scores_in)), np.zeros(len(scores_out))])
    s = np.concatenate([scores_in, scores_out])
    return float(roc_auc_score(y, s))


def shadow_attack_auc(
    calc,
    phenopackets,
    n_shadows=3,
    shadow_frac=0.5,
    n_per_class=80,
    rng=None,
):
    rng = rng or np.random.default_rng(0)
    n = len(phenopackets)
    X_train, y_train = [], []
    for _ in range(n_shadows):
        size = int(n * shadow_frac)
        in_idx = set(rng.choice(n, size=size, replace=False).tolist())
        shadow_db = [phenopackets[i] for i in in_idx]
        in_pick = list(in_idx)[:n_per_class]
        out_pick = [i for i in range(n) if i not in in_idx][:n_per_class]
        for i in in_pick:
            X_train.append(query_features(calc, phenopackets[i], shadow_db))
            y_train.append(1)
        for i in out_pick:
            X_train.append(query_features(calc, phenopackets[i], shadow_db))
            y_train.append(0)
    X_train, y_train = np.array(X_train), np.array(y_train)
    clf = RandomForestClassifier(n_estimators=100, max_depth=10, random_state=0, n_jobs=-1)
    clf.fit(X_train, y_train)

    # Held-out test on a fresh split
    size = int(n * shadow_frac)
    in_idx = set(rng.choice(n, size=size, replace=False).tolist())
    test_db = [phenopackets[i] for i in in_idx]
    test_in = list(in_idx)[:n_per_class]
    test_out = [i for i in range(n) if i not in in_idx][:n_per_class]
    X_test, y_test = [], []
    for i in test_in:
        X_test.append(query_features(calc, phenopackets[i], test_db))
        y_test.append(1)
    for i in test_out:
        X_test.append(query_features(calc, phenopackets[i], test_db))
        y_test.append(0)
    X_test, y_test = np.array(X_test), np.array(y_test)
    return float(roc_auc_score(y_test, clf.predict_proba(X_test)[:, 1]))


def run_mia_sweep(phenopackets, epsilons, seed=0):
    base_calc, _ = build_base_calculator(phenopackets)
    rng = np.random.default_rng(seed)
    n = len(phenopackets)
    # Fixed in/out split for threshold attack so ε values are comparable
    member_idx = rng.choice(n, size=n // 2, replace=False).tolist()
    member_set = set(member_idx)
    nonmember_idx = [i for i in range(n) if i not in member_set]

    rows = []
    for eps in epsilons:
        label = "inf" if not np.isfinite(eps) else f"{eps:g}"
        log.info("MIA sweep — ε = %s", label)
        if not np.isfinite(eps):
            calc = base_calc
        else:
            calc = PrivateSimilarityCalculator(
                base_calculator=base_calc,
                mechanism=LaplaceMechanism(epsilon=eps),
                sensitivity=1.0,
            )
        thr_auc = threshold_attack_auc(calc, phenopackets, member_idx, nonmember_idx, rng=rng)
        sh_auc = shadow_attack_auc(calc, phenopackets, rng=rng)
        rows.append(
            {
                "epsilon": eps if np.isfinite(eps) else np.inf,
                "epsilon_label": label,
                "threshold_attack_auc": thr_auc,
                "shadow_attack_auc": sh_auc,
                "threshold_advantage": max(thr_auc, 1 - thr_auc) - 0.5,
                "shadow_advantage": max(sh_auc, 1 - sh_auc) - 0.5,
            }
        )
        log.info("  threshold AUC = %.3f  shadow AUC = %.3f", thr_auc, sh_auc)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# k-anonymity singling-out ablation
# ---------------------------------------------------------------------------
# We pick rare HPO terms (held by ≤ K_MAX patients), then for each k, ask:
# given a query containing only that rare term, would the k-anon gate
# permit the singling-out result to be released?

def kanon_ablation(phenopackets, k_values, rare_max=50):
    """k-anonymity ablation against a single-rare-term adversary.

    Threat: the adversary knows their target patient possesses some rare HPO
    term t (e.g., from a clinical note or public registry). They query the
    system for patients with t. Without a k-anon gate, the system returns the
    set S_t of all patients holding t; with a gate at threshold k, the system
    returns S_t only if |S_t| ≥ k, otherwise emits ⊥.

    Re-identification probability for a successful query (uniform-random pick
    inside the returned cohort): 1 / |S_t|. Suppressed queries contribute 0.

    We aggregate over every rare term in the cohort (1 ≤ |S_t| ≤ rare_max).
    """
    term_to_patients = defaultdict(set)
    for i, pp in enumerate(phenopackets):
        for t in patient_terms(pp):
            term_to_patients[t].add(i)
    counts = np.array(
        sorted(
            len(s) for s in term_to_patients.values() if 1 <= len(s) <= rare_max
        )
    )
    log.info(
        "Considering %d rare terms with cohort size in [1, %d]; min=%d, median=%d, max=%d",
        len(counts), rare_max, counts.min(), int(np.median(counts)), counts.max(),
    )

    rows = []
    n = len(counts)
    # Pre-gate: every rare-term query released → re-id prob = 1/|S_t|.
    pre_reid = float(np.mean(1.0 / counts))
    pre_singling = float(np.mean(counts == 1))

    for k in k_values:
        released_mask = counts >= k
        suppressed = int((~released_mask).sum())
        # Post-gate re-id probability: 0 for suppressed; 1/|S_t| otherwise.
        per_query_reid = np.where(released_mask, 1.0 / counts, 0.0)
        post_reid = float(np.mean(per_query_reid))
        post_singling = float(np.mean(released_mask & (counts == 1)))
        mean_released_cohort = float(counts[released_mask].mean()) if released_mask.any() else 0.0
        rows.append(
            {
                "k": k,
                "n_rare_terms": n,
                "suppressed_queries": suppressed,
                "suppression_rate": suppressed / n,
                "reidentification_prob_pre_gate": pre_reid,
                "reidentification_prob_post_gate": post_reid,
                "singling_out_rate_pre_gate": pre_singling,
                "singling_out_rate_post_gate": post_singling,
                "mean_released_cohort_size": mean_released_cohort,
            }
        )
        log.info(
            "  k=%d: suppress=%.3f, re-id prob (post-gate)=%.4f, mean released cohort=%.2f",
            k,
            suppressed / n,
            post_reid,
            mean_released_cohort,
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_mia_pareto(df, out_base):
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    # Finite epsilons on log scale; inf as separate column
    finite = df[np.isfinite(df["epsilon"])].sort_values("epsilon")
    inf_row = df[~np.isfinite(df["epsilon"])]
    ax.plot(finite["epsilon"], finite["threshold_attack_auc"], "o-", label="Threshold attack (Yeom et al. 2018)")
    ax.plot(finite["epsilon"], finite["shadow_attack_auc"], "s-", label="Shadow-model attack (Shokri et al. 2017)")
    ax.axhline(0.5, ls=":", color="gray", label="Random guessing")
    if not inf_row.empty:
        eps_max = finite["epsilon"].max()
        x_inf = eps_max * 3
        ax.scatter([x_inf], inf_row["threshold_attack_auc"], marker="o", facecolors="none", edgecolors="C0", s=80, label="No DP (ε=∞)")
        ax.scatter([x_inf], inf_row["shadow_attack_auc"], marker="s", facecolors="none", edgecolors="C1", s=80)
        ax.annotate("ε=∞", (x_inf, max(inf_row["threshold_attack_auc"].iloc[0], inf_row["shadow_attack_auc"].iloc[0])), textcoords="offset points", xytext=(5, 5))
    ax.set_xscale("log")
    ax.set_xlabel("Privacy budget ε (log scale)")
    ax.set_ylabel("Attack ROC AUC")
    ax.set_ylim(0.45, 1.0)
    ax.set_title("Membership-inference attack vs. DP budget\n(Laplace mechanism, IC-weighted cosine similarity)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=180)
    plt.close(fig)


def plot_kanon_ablation(df, out_base):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    ax1.plot(df["k"], df["suppression_rate"], "o-", color="C2")
    ax1.set_xlabel("k-anonymity threshold k")
    ax1.set_ylabel("Suppression rate")
    ax1.set_title("(a) Fraction of rare-term queries blocked")
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 1.02)

    pre = df["reidentification_prob_pre_gate"].iloc[0]
    ax2.axhline(pre, ls="--", color="C3", label=f"No k-anon gate (baseline = {pre:.3f})")
    ax2.plot(df["k"], df["reidentification_prob_post_gate"], "o-", color="C0", label="With k-anon gate")
    ax2.set_xlabel("k-anonymity threshold k")
    ax2.set_ylabel("Expected re-identification probability")
    ax2.set_title("(b) Re-identification probability vs. k")
    ax2.legend(loc="best")
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(0, max(pre * 1.1, df["reidentification_prob_post_gate"].max() * 1.1))
    fig.suptitle("k-anonymity ablation: rare HPO-term singling-out attack", y=1.02)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(out_base.with_suffix(f".{ext}"), dpi=180, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--epsilons", nargs="+", type=float, default=[0.1, 0.5, 1.0, 2.0, 5.0, 10.0, float("inf")])
    p.add_argument("--k-values", nargs="+", type=int, default=[1, 2, 5, 10, 20, 50])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--skip-mia", action="store_true")
    p.add_argument("--skip-kanon", action="store_true")
    args = p.parse_args()
    np.random.seed(args.seed)

    pps = load_cohort()
    log.info("Loaded %d phenopackets", len(pps))

    if not args.skip_mia:
        mia_df = run_mia_sweep(pps, args.epsilons, seed=args.seed)
        out = RESULTS_DIR / "mia_epsilon_sweep.csv"
        mia_df.to_csv(out, index=False)
        log.info("Wrote %s", out)
        plot_mia_pareto(mia_df, FIG_DIR / "fig10_mia_pareto")
        log.info("Wrote MIA Pareto figure (fig10)")

    if not args.skip_kanon:
        kdf = kanon_ablation(pps, args.k_values)
        out = RESULTS_DIR / "kanon_ablation.csv"
        kdf.to_csv(out, index=False)
        log.info("Wrote %s", out)
        plot_kanon_ablation(kdf, FIG_DIR / "fig11_kanon_ablation")
        log.info("Wrote k-anon ablation figure (fig11)")

    log.info("Done.")


if __name__ == "__main__":
    main()
