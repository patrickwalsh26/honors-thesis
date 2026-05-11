#!/usr/bin/env python3
"""
Build a self-contained pilot-app cohort JSON.

Reads the unpacked Phenopacket Store archive (~9.6k phenopackets) under
data/phenopacket_store/, applies the same filter used in
experiments/evaluate_phenopacket_store.py (top-N OMIM diseases, capped per
disease, min phenotypes), strips patient identifiers down to anonymous indices,
and writes a single JSON the Streamlit pilot app can load at startup without
needing the full corpus or the project's heavy dependencies.

Output: app/data/pilot_cohort.json   (≈600 KB, ~1500 patients, 100 diseases)
"""

import argparse
import json
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.data_integration.phenopacket_store_loader import CohortSpec, load_phenopacket_store

logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
log = logging.getLogger("build_pilot_cohort")

REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUT = REPO / "app" / "data" / "pilot_cohort.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n-diseases", type=int, default=100)
    ap.add_argument("--max-per-disease", type=int, default=15)
    ap.add_argument("--min-phenotypes", type=int, default=3)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    spec = CohortSpec(
        min_phenotypes=args.min_phenotypes,
        min_cohort_size=2,
        max_per_disease=args.max_per_disease,
        top_n_diseases=args.top_n_diseases,
        seed=0,
    )
    pps, gt, stats = load_phenopacket_store(spec=spec)
    log.info("Filtered cohort: %s", stats)

    # Read HPO labels once from the canonical .obo to enrich the pilot UI.
    labels = {}
    obo = REPO / "data" / "hpo_ontology" / "hp.obo"
    if obo.exists():
        cur = None
        for line in obo.read_text().splitlines():
            if line == "[Term]":
                cur = None
            elif line.startswith("id: HP:"):
                cur = line.split("id: ", 1)[1].strip()
            elif line.startswith("name: ") and cur:
                labels[cur] = line.split("name: ", 1)[1].strip()
        log.info("Loaded %d HPO labels", len(labels))

    # Build a denormalised, identifier-free payload.
    out_patients = []
    for i, pp in enumerate(pps):
        terms = [
            f["type"]["id"]
            for f in pp.get("phenotypicFeatures", [])
            if not f.get("excluded", False)
        ]
        out_patients.append(
            {
                "anon_id": f"patient_{i:04d}",
                "disease_id": gt[pp["id"]],
                "phenotypes": terms,
                "sex": pp.get("subject", {}).get("sex", "UNKNOWN"),
            }
        )

    # Disease label table: take any patient's disease entry that has a label.
    disease_labels = {}
    for pp in pps:
        for d in pp.get("diseases", []):
            did = d.get("term", {}).get("id")
            lbl = d.get("term", {}).get("label")
            if did and lbl and did not in disease_labels:
                disease_labels[did] = lbl
        for itp in pp.get("interpretations", []):
            d = itp.get("diagnosis", {}).get("disease", {})
            if d.get("id") and d.get("label") and d["id"] not in disease_labels:
                disease_labels[d["id"]] = d["label"]

    payload = {
        "schema_version": 1,
        "n_patients": len(out_patients),
        "n_diseases": stats["n_diseases"],
        "mean_patients_per_disease": stats["mean_patients_per_disease"],
        "patients": out_patients,
        "disease_labels": disease_labels,
        "hpo_labels": labels,
        "source": "Monarch Phenopacket Store (Danis et al. 2025); filtered",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload))
    log.info("Wrote %s (%d patients, %d diseases, %d KB)",
             args.out, payload["n_patients"], payload["n_diseases"],
             args.out.stat().st_size // 1024)


if __name__ == "__main__":
    main()
