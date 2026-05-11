# STARR Dashboard Query Plan — HPO Term Prevalence for Empirical IC

**Purpose:** Replace circular HPOA-derived information content (IC) with priors calibrated against the STARR-OMOP-deid clinical population (~3M patients at Stanford). Output feeds the Resnik/Lin similarity computations and rare-term filtering thresholds in `src/similarity/hpo_similarity.py` and `src/privacy/k_anonymity.py`.

**Access constraint:** STARR Cohort Discovery dashboard only — aggregate counts, no patient-level export. This matches the Beacon v2 / Matchmaker Exchange threat model the thesis already evaluates, so the resulting numbers also serve as a real-world demonstration of that threat model.

---

## 1. What we need from STARR

For each HPO term `t` in our evaluation cohort, we need:

- `starr_count(t)` — patients in STARR-OMOP-deid with **any** encounter coded under the ICD-10 codes that crosswalk to `t`.
- `starr_denominator` — total distinct patients in STARR-OMOP-deid (record this **once**).

From these we compute empirical prevalence `p(t) = starr_count(t) / starr_denominator` and corpus IC `IC(t) = -log2(p(t))`. Terms with `starr_count(t) = 0` get IC = log2(starr_denominator) (effectively unseen → maximally informative, capped).

---

## 2. STARR Cohort Discovery setup (one-time)

1. Log in to the STARR Tools portal → open **Cohort Discovery** (or "STARR-OMOP-deid Self-Service" — whichever your dashboard surface labels it).
2. Confirm the data source is `starr_omop_deid_latest` or the most current de-identified release. Record the release tag in the paste-back CSV (`starr_release` column).
3. Create a baseline "All patients" cohort by leaving the criteria empty. Record the total patient count — this is the **denominator**. Put it in `data/starr/term_counts.csv` once.
4. For each term in §4 below, build a cohort whose only criterion is **"Condition occurrence — ICD-10-CM code in [list]"**. Use the codes in column `icd10_codes` exactly as written; the dashboard's value-set picker accepts comma-separated lists or wildcards.
5. Record the resulting patient count. **Do not** export rows; you only need the count.

Time estimate: 60–80 queries × ~30 sec each = 30–45 minutes once you're set up.

---

## 3. Two-pass strategy

**Pass A (today, ~45 min):** the top 60 most-frequent HPO terms from the evaluation cohort (this file's §4). These dominate the IC signal because they appear in the most retrieval comparisons.

**Pass B (optional, day 2):** mid-frequency terms (cohort frequency 5–11), another ~80 terms. If you run out of time, skip — Pass A alone is sufficient for the thesis result.

Long-tail rare terms (cohort frequency ≤ 4) do **not** need STARR queries: their HPOA-derived IC is already near the cap, and any STARR count will also be very low. They're noise-floor terms either way.

---

## 4. Term list — Pass A (top 60 by cohort frequency)

The full paste-back template is at `data/starr/term_counts.csv`. Below is the query list. Columns:

- `hpo_id` — HPO term to be calibrated
- `label` — for readability only, not used downstream
- `cohort_count` — how often the term appears in our 500-patient eval cohort (sanity check, not for STARR)
- `icd10_codes` — what to paste into the STARR Cohort Discovery condition filter
- `mapping_confidence` — `high` (1:1-ish), `medium` (lossy but defensible), `low` (skip if short on time; will be excluded from primary IC table)

| hpo_id | label | cohort_count | icd10_codes | mapping_confidence |
|---|---|---|---|---|
| HP:0001249 | Intellectual disability | 70 | F70, F71, F72, F73, F78, F79 | high |
| HP:0001263 | Global developmental delay | 61 | F88, F89, R62.50 | high |
| HP:0001250 | Seizure | 52 | R56.9, G40.* | high |
| HP:0002650 | Scoliosis | 46 | M41.* | high |
| HP:0001252 | Hypotonia | 38 | P94.2, R29.898 | medium |
| HP:0000639 | Nystagmus | 33 | H55.0* | high |
| HP:0000252 | Microcephaly | 29 | Q02 | high |
| HP:0004322 | Short stature | 28 | E34.3, R62.52 | high |
| HP:0002240 | Hepatomegaly | 26 | R16.0 | high |
| HP:0001260 | Dysarthria | 26 | R47.1 | high |
| HP:0001508 | Failure to thrive | 26 | R62.51, P92.6 | high |
| HP:0001347 | Hyperreflexia | 23 | R29.2 | medium |
| HP:0000347 | Micrognathia | 21 | Q67.4, K07.0 | medium |
| HP:0000407 | Sensorineural hearing impairment | 21 | H90.3, H90.5 | high |
| HP:0000505 | Visual impairment | 20 | H54.* | high |
| HP:0001344 | Absent speech | 20 | F80.0, R47.01 | medium |
| HP:0000750 | Delayed speech and language development | 20 | F80.1, F80.2, F80.9 | high |
| HP:0001324 | Muscle weakness | 19 | M62.81 | high |
| HP:0001945 | Fever | 19 | R50.9 | high |
| HP:0002094 | Dyspnea | 19 | R06.00, R06.02 | high |
| HP:0000729 | Autistic behavior | 19 | F84.0, F84.9 | high |
| HP:0001744 | Splenomegaly | 17 | R16.1 | high |
| HP:0000218 | High palate | 17 | Q38.5 | medium |
| HP:0001257 | Spasticity | 17 | R25.2 | medium |
| HP:0001251 | Ataxia | 17 | R27.0, R27.8, R27.9 | high |
| HP:0001288 | Gait disturbance | 17 | R26.0, R26.2, R26.89, R26.9 | high |
| HP:0003487 | Babinski sign | 17 | R29.2 | low |
| HP:0001321 | Cerebellar hypoplasia | 17 | Q04.3 | medium |
| HP:0001511 | Intrauterine growth retardation | 16 | P05.* | high |
| HP:0002059 | Cerebral atrophy | 16 | G31.9 | medium |
| HP:0002013 | Vomiting | 15 | R11.10, R11.11, R11.2 | high |
| HP:0000316 | Hypertelorism | 15 | Q75.2 | high |
| HP:0003510 | Severe short stature | 15 | E34.3 | medium |
| HP:0000175 | Cleft palate | 15 | Q35.* | high |
| HP:0000648 | Optic atrophy | 15 | H47.2* | high |
| HP:0001382 | Joint hypermobility | 15 | M35.7 | high |
| HP:0001903 | Anemia | 15 | D64.9 | high |
| HP:0002643 | Neonatal respiratory distress | 15 | P22.* | high |
| HP:0002808 | Kyphosis | 15 | M40.* | high |
| HP:0001332 | Dystonia | 14 | G24.* | high |
| HP:0001761 | Pes cavus | 14 | Q66.7 | high |
| HP:0001270 | Motor delay | 14 | R62.0 | medium |
| HP:0010864 | Severe intellectual disability | 14 | F72 | high |
| HP:0002342 | Moderate intellectual disability | 14 | F71 | high |
| HP:0000739 | Anxiety | 14 | F41.* | high |
| HP:0002014 | Diarrhea | 14 | R19.7, K59.1 | high |
| HP:0002119 | Ventriculomegaly | 13 | Q03.*, G91.9 | medium |
| HP:0002079 | Hypoplasia of the corpus callosum | 13 | Q04.0 | high |
| HP:0001387 | Joint stiffness | 13 | M25.60–M25.69 | high |
| HP:0002515 | Waddling gait | 13 | R26.0 | low |
| HP:0001272 | Cerebellar atrophy | 13 | G31.9 | medium |
| HP:0000518 | Cataract | 13 | H25.*, H26.* | high |
| HP:0002716 | Lymphadenopathy | 13 | R59.* | high |
| HP:0012733 | Macule | 13 | — | skip (low) |
| HP:0000083 | Renal insufficiency | 13 | N18.* | high |
| HP:0002007 | Frontal bossing | 13 | — | skip (low) |
| HP:0001873 | Thrombocytopenia | 12 | D69.6 | high |
| HP:0007663 | Reduced visual acuity | 12 | H54.7, H53.8 | medium |
| HP:0000486 | Strabismus | 12 | H50.* | high |
| HP:0005518 | Increased mean corpuscular volume | 12 | — | skip (lab value, not condition) |

**Mapping caveats** (these go into Methods §4.2 / the appendix):

- HPO terms are finer-grained than ICD-10. Some HPO terms (e.g., severity-specific intellectual disability subtypes) map to a single ICD-10 code; others (e.g., "Hypotonia" vs. neonatal-specific "Floppy infant") collapse multiple HPO concepts.
- ICD-10 codes appear in claim/encounter records; their prevalence reflects clinical recognition, not biological prevalence. Some highly informative HPO terms (e.g., "Babinski sign") are physical-exam findings rarely coded as primary conditions → STARR counts will under-estimate true prevalence. We report and discuss this directionally.
- We mark `skip` for terms with no defensible ICD-10 mapping (radiology signs, lab values, very general physical findings). These will retain their HPOA-derived IC in the final analysis — we'll note this explicitly in Results.

---

## 5. Paste-back CSV

`data/starr/` is gitignored (intentionally — it's where any populated STARR-derived data lives, and we don't ever want to risk committing clinical data). Templates live in `scripts/starr/`:

1. Copy `scripts/starr/term_counts_template.csv` → `data/starr/term_counts.csv`. Fill in:
   - `starr_count` — patient count from the dashboard for each row
   - `starr_release` — STARR-OMOP-deid release tag (same value for every row; e.g., `starr_omop_deid_2026q1`)
2. Copy `scripts/starr/starr_meta_template.csv` → `data/starr/starr_meta.csv`. Fill in:
   - `starr_denominator` — total patient count for the "All patients" baseline cohort
   - `starr_release`, `query_date`

When you're done, paste the CSV contents back in chat (or commit them locally — they're gitignored, so they won't push). I'll run the IC recomputation and rerun figures 3, 5, 6, 7 against STARR-calibrated priors.

---

## 6. What I will do once your CSV is back

1. `experiments/recompute_ic_from_starr.py` — load `data/starr/term_counts.csv`, fold in HPO ancestor closure (a parent term's count is the union of its descendants' counts in STARR), produce `data/starr/empirical_ic.json`.
2. Patch `src/similarity/hpo_similarity.py` to accept an external IC table (drop-in alternative to `compute_empirical_ic`).
3. Re-run `experiments/evaluate_hpoa.py` and `experiments/evaluate_decipher.py` with both HPOA-IC and STARR-IC, emit a side-by-side comparison table for Results §4.3.
4. Re-run privacy-utility sweeps with the new IC. Update Figs 3, 5, 6, 7.
5. Add a Methods sub-section (§4.2.3 "STARR-OMOP-deid IC calibration") and a Results paragraph quantifying the IC shift (we expect HPOA to over-weight rare Mendelian terms and under-weight common findings; STARR will compress IC for common phenotypes like Anemia / Seizure).
