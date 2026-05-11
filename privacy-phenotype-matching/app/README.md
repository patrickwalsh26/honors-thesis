# Privacy-Preserving Phenotype Matching — Pilot App

Single-page Streamlit app demonstrating the thesis system. A clinician picks
observed HPO phenotypes, tunes privacy parameters (ε, k, rare-term filter), and
queries a reference cohort of 1,500 real published-case-report patients from the
Monarch Phenopacket Store (Danis et al. 2025).

## Run locally

```bash
pip install -r app/requirements.txt
streamlit run app/streamlit_app.py
```

Browser opens at `http://localhost:8501`. The cohort JSON
(`app/data/pilot_cohort.json`, ~1.3 MB) is loaded once and cached.

## Regenerate the cohort

```bash
# Requires the full corpus download (scripts/download_phenopacket_store.sh).
python scripts/build_pilot_cohort.py
```

## Deploy to Streamlit Community Cloud

1. Push to `main` (the app reads only files committed to the repo).
2. <https://share.streamlit.io> → **New app** → connect the `patrickwalsh26/honors-thesis` repository.
3. Set:
   - **Main file path:** `privacy-phenotype-matching/app/streamlit_app.py`
   - **Python version:** 3.11
   - **Branch:** `main`
4. **Advanced settings → Python dependencies:** point at `privacy-phenotype-matching/app/requirements.txt` (this is critical — the repo-root `requirements.txt` is the heavy research stack and is not needed for the app).
5. Deploy. First boot ~30 s; subsequent reloads are cached.

## What the app demonstrates

- IC-weighted cosine retrieval against the reference cohort
- Differential privacy (Laplace mechanism) on similarity scores, ε ∈ {0.5, 1, 2, 5, 10, 20, 50, ∞}
- k-anonymity gate on result release, k ∈ {1, 2, 5, 10, 20}
- Rare-term prevalence filter applied to the query before scoring
- Session privacy budget accounting (cumulative ε across queries)

What it does **not** demonstrate (out of scope for the v1 pilot):

- Live private-set-intersection between two institutions — the app runs the
  retrieval locally because PSI semantic security would not be visible to the
  user. The released score has the same privacy properties as the federated
  protocol would produce.
- EHR FHIR import, authenticated user identity, per-institution budget accounting.

See thesis §3.7 (Reference Implementation) and §4.7 (Pilot System) for the full
write-up.
