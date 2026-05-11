"""
Loader for the Monarch Phenopacket Store release archive.

Phenopacket Store (https://github.com/monarch-initiative/phenopacket-store)
is a curated, gene-organized collection of GA4GH Phenopackets v2 derived
from published case reports. Each phenopacket has a confirmed OMIM disease
diagnosis and PMID provenance.

This loader filters the corpus to retrieval-evaluable patients (at least
``min_phenotypes`` observed features and a disease cohort of at least
``min_cohort_size`` patients), optionally subsamples per disease for
balanced evaluation, and returns a list of GA4GH phenopacket dicts plus a
ground-truth mapping ``phenopacket_id -> disease_id``.
"""

from __future__ import annotations

import glob
import json
import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

log = logging.getLogger(__name__)

DEFAULT_ROOT = "data/phenopacket_store"


@dataclass
class CohortSpec:
    """Filtering / subsampling parameters for a phenopacket-store cohort."""

    min_phenotypes: int = 3
    min_cohort_size: int = 2
    max_per_disease: Optional[int] = None  # cap per-disease count for balance
    top_n_diseases: Optional[int] = None   # keep only top-N most-populated
    seed: int = 0


def _resolve_release_dir(root: str) -> str:
    """Find the versioned subdirectory inside ``root``."""
    candidates = sorted(
        d for d in glob.glob(os.path.join(root, "*"))
        if os.path.isdir(d) and "." in os.path.basename(d)
    )
    if not candidates:
        raise FileNotFoundError(
            f"No versioned phenopacket-store release found under {root}. "
            f"Run scripts/download_phenopacket_store.sh first."
        )
    return candidates[-1]


def _extract_disease_id(pp: Dict) -> Optional[str]:
    for itp in pp.get("interpretations", []):
        d = itp.get("diagnosis", {}).get("disease", {}).get("id")
        if d:
            return d
    for d in pp.get("diseases", []):
        did = d.get("term", {}).get("id")
        if did:
            return did
    return None


def _observed_terms(pp: Dict) -> List[str]:
    return [
        f["type"]["id"]
        for f in pp.get("phenotypicFeatures", [])
        if not f.get("excluded", False)
    ]


def load_phenopacket_store(
    root: str = DEFAULT_ROOT,
    spec: Optional[CohortSpec] = None,
) -> Tuple[List[Dict], Dict[str, str], Dict[str, int]]:
    """Load and filter the corpus.

    Returns
    -------
    phenopackets : list of phenopacket dicts (passing all filters)
    ground_truth : mapping phenopacket id -> disease id (OMIM:NNN)
    stats        : counts useful for reporting (input, after-filter, etc.)
    """
    spec = spec or CohortSpec()
    rng = np.random.default_rng(spec.seed)
    release_dir = _resolve_release_dir(root)
    paths = glob.glob(os.path.join(release_dir, "*", "*.json"))
    log.info("Reading %d phenopackets from %s", len(paths), release_dir)

    by_disease: Dict[str, List[Tuple[Dict, str]]] = defaultdict(list)
    n_no_disease = 0
    n_too_few_pheno = 0
    for p in paths:
        with open(p) as fh:
            pp = json.load(fh)
        did = _extract_disease_id(pp)
        if not did:
            n_no_disease += 1
            continue
        if len(_observed_terms(pp)) < spec.min_phenotypes:
            n_too_few_pheno += 1
            continue
        by_disease[did].append((pp, p))

    # Drop diseases with too few patients.
    by_disease = {
        d: members for d, members in by_disease.items() if len(members) >= spec.min_cohort_size
    }
    if spec.top_n_diseases:
        keep = sorted(by_disease.items(), key=lambda kv: -len(kv[1]))[: spec.top_n_diseases]
        by_disease = dict(keep)
    if spec.max_per_disease:
        for d, members in list(by_disease.items()):
            if len(members) > spec.max_per_disease:
                idx = rng.choice(len(members), size=spec.max_per_disease, replace=False)
                by_disease[d] = [members[i] for i in sorted(idx.tolist())]

    phenopackets: List[Dict] = []
    ground_truth: Dict[str, str] = {}
    for d, members in by_disease.items():
        for pp, _path in members:
            phenopackets.append(pp)
            ground_truth[pp["id"]] = d

    stats = {
        "n_input": len(paths),
        "n_no_disease": n_no_disease,
        "n_too_few_phenotypes": n_too_few_pheno,
        "n_after_filter": len(phenopackets),
        "n_diseases": len(by_disease),
        "mean_patients_per_disease": (
            float(np.mean([len(m) for m in by_disease.values()])) if by_disease else 0.0
        ),
    }
    log.info("Filtered cohort: %s", stats)
    return phenopackets, ground_truth, stats
