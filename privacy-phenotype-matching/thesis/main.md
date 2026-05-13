# Privacy-Preserving Phenotype Matching for Rare Disease Cohort Discovery

---

**Patrick Walsh**

Honors Thesis

Department of Computer Science
Stanford University

April 2026

---

**Thesis Advisor:** Professor Stephen B. Montgomery, Ph.D.
Department of Pathology, Genetics, and Biomedical Data Science
Stanford University School of Medicine

---

## Abstract

Rare diseases affect 300 million people worldwide; the diagnostic odyssey averages 4.8–7 years. Federated patient matching across institutions can shorten it, but sharing phenotype data is risky — rare phenotype combinations act as quasi-identifiers. We present a privacy-preserving phenotype matching framework that composes Private Set Intersection, differential privacy, and k-anonymity with rare-term filtering on GA4GH Phenopackets, compatible with Matchmaker Exchange and Beacon v2. A two-party semi-honest threat model with three adversary goals (membership inference, attribute inference, singling-out re-identification) anchors per-step disclosure analysis and three composition-level privacy invariants.

Evaluation uses two cohorts: a synthetic cohort sampled from 12,974 OMIM/Orphanet disease profiles, and the Monarch Phenopacket Store (Danis et al., 2025) — 9,588 real published case-report phenopackets with confirmed OMIM diagnoses. Non-private IC-weighted cosine retrieval on the real cohort achieves MRR = 0.87 and nDCG@10 = 0.69, within the Phenomizer/LIRICAL band. Shadow-model membership-inference attack AUC drops from 0.98 to 0.50 (random) at ε ≤ 1; k-anonymity at k = 10 cuts re-identification probability against the rare-term adversary from 0.42 to 0.005.

The central deployment finding: the safe Laplace-DP budget on real patients is 20–50× larger than synthetic-cohort experiments imply, because real similarity-score distributions are compressed and per-score Laplace noise dominates them. We diagnose this gap and validate the principled fix — the iterative exponential mechanism on a rank utility, whose noise tracks rank gaps rather than score magnitudes — empirically: 90% of non-private nDCG@10 retained at ε = 5 versus 13% for Laplace at matched ε-DP, a 10× budget-efficiency improvement. The open-source implementation, a reproducibility pipeline that regenerates every figure in four minutes, and an interactive pilot system at <https://honors-thesis-54tubqjkgwqjm4zyegglxw.streamlit.app/> support immediate adoption.

**Keywords:** rare diseases, phenotype matching, differential privacy, private set intersection, k-anonymity, Human Phenotype Ontology, GA4GH Phenopackets, membership inference

---

## Acknowledgments

I am deeply grateful to my thesis advisor, **Professor Stephen B. Montgomery**, for his guidance, mentorship, and support throughout this project. His work on rare-disease genomics and the application of transcriptomics to diagnosis shaped the framing of this thesis.

I thank the members of the Montgomery Lab for valuable discussions and feedback, and the broader rare-disease research community — particularly the developers of the Human Phenotype Ontology, GA4GH Phenopackets, the Monarch Initiative, and Matchmaker Exchange — whose open standards and curated data make work of this kind possible.

---

## Table of Contents

1. [Introduction](#introduction)
2. [Literature Review](#literature-review)
3. [Methods](#methods)
4. [Results](#results)
5. [Discussion](#discussion)
6. [Future Work](#future-work)
7. [Conclusion](#conclusion)
8. [Bibliography](#bibliography)
9. [Appendices](#appendices)

---

## List of Figures

- **Figure 1:** System architecture for privacy-preserving phenotype matching
- **Figure 2:** Human Phenotype Ontology structure (example subgraph)
- **Figure 3:** Privacy-preserving pipeline workflow
- **Figure 4:** Similarity score distributions across metrics
- **Figure 5:** Precision-Recall curves by similarity metric
- **Figure 6:** nDCG@k vs. ε (synthetic-cohort differential privacy tradeoff)
- **Figure 7:** Suppression rate vs. k (k-anonymity tradeoff)
- **Figure 8:** Privacy-utility Pareto frontier (synthetic cohort)
- **Figure 9:** Membership inference attack ROC curves
- **Figure 10:** Empirical MI attack AUC vs. ε (threshold and shadow-model attackers)
- **Figure 11:** k-anonymity ablation: suppression and re-identification probability
- **Figure 12:** Real-cohort retrieval (Phenopacket Store): privacy-utility curve and Phenomizer-style baseline
- **Figure 13:** Pilot Streamlit application: clinician-facing demo of the privacy-preserving phenotype matching pipeline
- **Figure 14:** Rank-based DP recovers retrieval utility on the real cohort: comparison of Laplace, score-utility exponential mechanism, and rank-utility exponential mechanism at matched ε-DP

---

## List of Tables

- **Table 1:** HPO Annotations Corpus Summary
- **Table 2:** Evaluation Cohort Characteristics
- **Table 3:** Per-step disclosure to a semi-honest adversary (Methods §3.1.2)
- **Tables 1–4 (Results):** Cohort statistics and example phenotypes
- **Tables 5–8:** Baseline retrieval performance and similarity distributions
- **Tables 9–12:** Synthetic-cohort privacy-utility tradeoffs (DP, k-anonymity, rare-term, composition)
- **Table 13:** Membership-inference attack ROC AUC vs. ε
- **Table 14:** k-anonymity ablation against rare-term singling-out
- **Table 15:** Retrieval on the real Phenopacket Store cohort
- **Table 16:** Cosine-IC retrieval under Laplace DP on the real cohort
- **Table 17:** Three ε-DP top-k release mechanisms (Laplace, score-exp, rank-exp)
- **Table 18:** Computational overhead

---

<!-- SECTION: INTRODUCTION -->
<!-- Include: sections/introduction.md -->

---

<!-- SECTION: LITERATURE REVIEW -->
<!-- Include: sections/literature_review.md -->

---

<!-- SECTION: METHODS -->
<!-- Include: sections/methods.md -->

---

<!-- SECTION: RESULTS -->
<!-- Include: sections/results.md -->

---

<!-- SECTION: DISCUSSION -->
<!-- Include: sections/discussion.md -->

---

<!-- SECTION: FUTURE WORK -->
<!-- Include: sections/future_work.md -->

---

<!-- SECTION: CONCLUSION -->
<!-- Include: sections/conclusion.md -->

---

<!-- SECTION: BIBLIOGRAPHY -->
<!-- Include: sections/bibliography.md -->

---

<!-- SECTION: APPENDICES -->
<!-- Include: sections/appendices.md -->

---

*This thesis was prepared in partial fulfillment of the requirements for the degree of Bachelor of Science with Honors in Computer Science at Stanford University.*

*April 2026*
