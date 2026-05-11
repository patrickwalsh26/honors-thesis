"""
Pilot system for privacy-preserving phenotype matching.

Single-page Streamlit app that lets a clinician select observed HPO phenotypes
for a query patient, choose privacy mechanism parameters (differential privacy
budget ε, k-anonymity threshold, rare-term filtering prevalence), and obtain
top-k similar patients from a 1,500-patient reference cohort drawn from the
Monarch Phenopacket Store (Danis et al. 2025).

Deployment:
  - Local:           streamlit run app/streamlit_app.py
  - Streamlit Cloud: point a new app at this file in the public GitHub repo.

The app is self-contained: it loads app/data/pilot_cohort.json at startup
and never touches the heavy crypto / ontology stack.
"""

from __future__ import annotations

import json
import math
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import streamlit as st

APP_DIR = Path(__file__).resolve().parent
DATA_PATH = APP_DIR / "data" / "pilot_cohort.json"
REPO_URL = "https://github.com/patrickwalsh26/honors-thesis"

st.set_page_config(
    page_title="Privacy-Preserving Phenotype Matching — Pilot",
    page_icon=":lock:",
    layout="wide",
)


# ---------------------------------------------------------------------------
# Data loading and IC precomputation (cached once per session)
# ---------------------------------------------------------------------------

@st.cache_resource(show_spinner="Loading reference cohort…")
def load_cohort() -> dict:
    with open(DATA_PATH) as f:
        return json.load(f)


@st.cache_resource(show_spinner="Computing information-content priors…")
def compute_ic(_cohort: dict) -> Tuple[Dict[str, float], Counter]:
    """Corpus IC: -log2(prevalence). One vote per patient per term."""
    n = len(_cohort["patients"])
    counts: Counter = Counter()
    for pat in _cohort["patients"]:
        counts.update(set(pat["phenotypes"]))
    ic = {t: -math.log2(c / n) for t, c in counts.items()}
    return ic, counts


# ---------------------------------------------------------------------------
# Privacy mechanism implementations (minimal, app-local)
# ---------------------------------------------------------------------------

def rare_term_filter(
    terms: List[str],
    counts: Counter,
    n_total: int,
    prevalence_threshold: float,
) -> Tuple[List[str], List[str]]:
    """Suppress terms with prevalence < threshold. Returns (kept, filtered)."""
    kept, dropped = [], []
    cutoff = prevalence_threshold * n_total
    for t in terms:
        if counts.get(t, 0) >= cutoff:
            kept.append(t)
        else:
            dropped.append(t)
    return kept, dropped


def cosine_ic(query: List[str], target: List[str], ic: Dict[str, float]) -> float:
    """IC-weighted cosine over the shared term set."""
    qset, tset = set(query), set(target)
    if not qset or not tset:
        return 0.0
    inter = qset & tset
    num = sum(ic.get(t, 0.0) ** 2 for t in inter)
    qnorm = math.sqrt(sum(ic.get(t, 0.0) ** 2 for t in qset))
    tnorm = math.sqrt(sum(ic.get(t, 0.0) ** 2 for t in tset))
    if qnorm == 0 or tnorm == 0:
        return 0.0
    return num / (qnorm * tnorm)


def laplace_noise(epsilon: float, sensitivity: float = 1.0) -> float:
    """ε-DP Laplace noise; sensitivity defaults to 1 (similarity in [0,1])."""
    if not math.isfinite(epsilon):
        return 0.0
    return float(np.random.laplace(0.0, sensitivity / max(epsilon, 1e-9)))


# ---------------------------------------------------------------------------
# Retrieval pipeline (composing the mechanisms)
# ---------------------------------------------------------------------------

def run_query(
    query_terms: List[str],
    cohort: dict,
    ic: Dict[str, float],
    counts: Counter,
    *,
    epsilon: float,
    k_anon: int,
    rare_threshold: float,
    top_k: int = 5,
) -> dict:
    n_total = len(cohort["patients"])

    # 1) rare-term filtering on the query
    kept_terms, dropped_terms = rare_term_filter(query_terms, counts, n_total, rare_threshold)

    # 2) score all patients
    scores: List[Tuple[int, float, float]] = []
    for i, pat in enumerate(cohort["patients"]):
        true_sim = cosine_ic(kept_terms, pat["phenotypes"], ic)
        noisy_sim = float(np.clip(true_sim + laplace_noise(epsilon, 1.0), 0.0, 1.0))
        scores.append((i, true_sim, noisy_sim))

    # 3) sort by noisy similarity, take top-k
    ranked = sorted(scores, key=lambda x: x[2], reverse=True)[:top_k]
    n_matches_above_zero = sum(1 for _, _, ns in ranked if ns > 0.0)

    # 4) k-anonymity gate
    suppressed = n_matches_above_zero < k_anon
    return {
        "kept_terms": kept_terms,
        "dropped_terms": dropped_terms,
        "ranked": ranked,
        "suppressed": suppressed,
        "n_matches_above_zero": n_matches_above_zero,
        "epsilon": epsilon,
        "k_anon": k_anon,
    }


# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

def label_for(hpo_id: str, cohort: dict) -> str:
    lbl = cohort["hpo_labels"].get(hpo_id, "")
    return f"{hpo_id} — {lbl}" if lbl else hpo_id


def main():
    cohort = load_cohort()
    ic, counts = compute_ic(cohort)

    # --- header ---
    st.title("Privacy-Preserving Phenotype Matching")
    st.caption(
        "Pilot for Patrick Walsh's Stanford CS Honors Thesis (2026). "
        f"[Source code]({REPO_URL}). "
        f"Reference cohort: {cohort['n_patients']} patients across "
        f"{cohort['n_diseases']} OMIM diseases, drawn from the Monarch Phenopacket Store "
        "(Danis et al. 2025)."
    )
    st.info(
        "**Disclaimer.** This is a research prototype. Reference patients are "
        "drawn from a public corpus of published case reports; no identifiable "
        "patient data is processed here. Do not use for clinical decision-making.",
        icon=":material/info:",
    )

    # --- privacy controls ---
    with st.sidebar:
        st.header("Privacy controls")
        epsilon = st.select_slider(
            "Differential-privacy budget ε",
            options=[0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, math.inf],
            value=5.0,
            format_func=lambda v: "∞ (no DP)" if not math.isfinite(v) else f"{v:g}",
            help="Smaller ε is stronger privacy. Real cohorts (this one) need ε ≥ 20 for usable retrieval; "
                 "synthetic cohorts permit ε ≈ 1. See thesis §4.6.",
        )
        k_anon = st.select_slider(
            "k-anonymity threshold k",
            options=[1, 2, 5, 10, 20],
            value=5,
            help="Results are suppressed if fewer than k patients match. k=1 disables the gate.",
        )
        rare_threshold = st.slider(
            "Rare-term filter (min prevalence)",
            min_value=0.0,
            max_value=0.10,
            value=0.01,
            step=0.005,
            help="Query terms appearing in fewer than this fraction of cohort patients are suppressed before scoring.",
        )
        top_k = st.select_slider("Top-k results", options=[3, 5, 10], value=5)

        # session privacy budget tracker — render lazily so it reflects the
        # increment that happens after the search button is clicked.
        if "spent" not in st.session_state:
            st.session_state.spent = 0.0
            st.session_state.queries = 0
        st.markdown("---")
        st.subheader("Session privacy budget")
        budget_metric = st.empty()
        budget_caption = st.empty()
        if st.button("Reset session"):
            st.session_state.spent = 0.0
            st.session_state.queries = 0
            st.rerun()

    # --- query input ---
    st.subheader("Query phenotypes")
    all_hpo = sorted(cohort["hpo_labels"].keys())
    label_map = {h: label_for(h, cohort) for h in all_hpo}

    preset_label = "(None — use the selector below)"
    presets = {
        preset_label: [],
        "Example: classic intellectual-disability + seizure presentation": [
            "HP:0001249", "HP:0001250", "HP:0001263", "HP:0000252",
        ],
        "Example: connective-tissue presentation": [
            "HP:0001382", "HP:0001634", "HP:0001763",
        ],
        "Example: short stature + craniofacial": [
            "HP:0004322", "HP:0000316", "HP:0000218",
        ],
    }
    preset_choice = st.selectbox("Start from a preset (optional)", list(presets.keys()))
    default_terms = presets[preset_choice]

    query_terms = st.multiselect(
        "Observed HPO terms",
        options=all_hpo,
        default=default_terms,
        format_func=lambda h: label_map[h],
        help="Select 3–15 phenotypic features observed in your query patient.",
    )

    run_clicked = st.button("Search reference cohort", type="primary", disabled=not query_terms)
    if not run_clicked:
        # still render the budget so user sees current totals
        budget_metric.metric("Cumulative ε spent", f"{st.session_state.spent:.2f}")
        budget_caption.caption(f"{st.session_state.queries} query(s) this session")
        st.stop()

    if math.isfinite(epsilon):
        st.session_state.spent += epsilon
    st.session_state.queries += 1
    budget_metric.metric("Cumulative ε spent", f"{st.session_state.spent:.2f}")
    budget_caption.caption(f"{st.session_state.queries} query(s) this session")

    result = run_query(
        query_terms, cohort, ic, counts,
        epsilon=epsilon, k_anon=k_anon, rare_threshold=rare_threshold, top_k=top_k,
    )

    # --- results ---
    st.subheader("Top matches")
    if result["suppressed"]:
        st.warning(
            f"**Suppressed by k-anonymity** (k = {k_anon}). Only "
            f"{result['n_matches_above_zero']} patient(s) returned a non-zero similarity, "
            "fewer than the configured threshold. No results released. Reduce k or relax "
            "the rare-term filter to see matches."
        )
    else:
        rows = []
        for rank, (idx, true_sim, noisy_sim) in enumerate(result["ranked"], start=1):
            pat = cohort["patients"][idx]
            disease_lbl = cohort["disease_labels"].get(pat["disease_id"], "")
            rows.append(
                {
                    "Rank": rank,
                    "Patient": pat["anon_id"],
                    "Disease (true label)": f"{pat['disease_id']} — {disease_lbl}",
                    "Shared phenotypes (count)": len(set(result["kept_terms"]) & set(pat["phenotypes"])),
                    "DP-noised similarity": round(noisy_sim, 3),
                    "True similarity (debug)": round(true_sim, 3),
                }
            )
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        st.caption(
            "Disease labels are shown for evaluation purposes only. In a real federated "
            "deployment the server would return only ranked anonymous IDs; the querying "
            "clinician would request follow-up via a separate authenticated channel."
        )

    # --- privacy trace ---
    with st.expander("Privacy trace (mechanism-by-mechanism)"):
        st.markdown(f"**1. Rare-term filter (prevalence ≥ {rare_threshold:.3f})**")
        if result["dropped_terms"]:
            st.write(
                "Filtered as rare (quasi-identifier risk): "
                + ", ".join(label_map.get(t, t) for t in result["dropped_terms"])
            )
        else:
            st.write("No query terms were below the rare-term threshold.")
        st.write(f"Terms passed to similarity scoring: {len(result['kept_terms'])}")

        eps_str = "∞ (no DP)" if not math.isfinite(epsilon) else f"{epsilon:g}"
        st.markdown(f"**2. Differential privacy (Laplace, ε = {eps_str})**")
        if math.isfinite(epsilon):
            noise_scale = 1.0 / epsilon
            st.write(
                f"Laplace noise with scale 1/ε = {noise_scale:.3f} added to each "
                f"similarity score; clipped to [0, 1]. Privacy invariant I2 (Methods §3.1.2)."
            )
        else:
            st.write("DP disabled — true similarity scores released.")

        st.markdown(f"**3. k-anonymity gate (k = {k_anon})**")
        st.write(
            f"Result released only if at least {k_anon} candidate patient(s) returned a "
            f"non-zero similarity. Privacy invariant I3 (Methods §3.1.2)."
        )
        st.write(
            f"Match cardinality: {result['n_matches_above_zero']} "
            f"{'(gate triggered: results suppressed)' if result['suppressed'] else '(gate passed)'}."
        )

    # --- mechanism explanations ---
    with st.expander("What is this system doing? (background)"):
        st.markdown(
            """
The pilot composes three privacy mechanisms from the thesis:

- **Rare-term filtering** suppresses quasi-identifying phenotypes (terms held by very few patients) before any scoring happens. This defends against re-identification via rare-phenotype quasi-identifiers (Sweeney, 2002; El Emam et al., 2011).
- **Differential privacy (Laplace mechanism)** adds calibrated noise to each similarity score so a single record's inclusion or exclusion cannot be detected from the output. We measured the empirical MI defense in §4.5.1 of the thesis.
- **k-anonymity gate** suppresses the response entirely if fewer than *k* patients match the query, preventing singling-out of small cohorts.

A real federated deployment would additionally run **private set intersection** between the query institution and the database institution so that no phenotype terms ever leave their respective sides in cleartext. The pilot omits the cryptographic step for usability and runs the equivalent local computation; the privacy-utility tradeoff for the released scores is unchanged.

See the thesis for the formal threat model (§3.1.2), empirical privacy attacks (§4.5), real-cohort retrieval benchmark (§4.6), and the deployment recommendation (§5.3.1).
"""
        )


if __name__ == "__main__":
    main()
