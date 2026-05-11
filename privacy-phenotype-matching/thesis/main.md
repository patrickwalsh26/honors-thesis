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

Rare diseases affect roughly 300 million people worldwide, and the diagnostic odyssey averages 4.8–7 years. Federated patient matching across institutions can accelerate diagnosis, but sharing phenotype data is risky: rare phenotype combinations act as quasi-identifiers. We present a privacy-preserving phenotype matching framework that composes Private Set Intersection (PSI) for cryptographic protection of phenotype sets, differential privacy (DP) for bounded score leakage, and k-anonymity with rare-term filtering for quasi-identifier defense. The system operates on GA4GH Phenopackets over the Human Phenotype Ontology and is compatible with Matchmaker Exchange and Beacon v2.

We formalize a two-party semi-honest threat model with three concrete adversary goals (membership inference, attribute inference, singling-out re-identification) and prove the corresponding privacy invariants for our composition. We evaluate the system on two cohorts: a synthetic cohort sampled from 12,974 OMIM/Orphanet disease profiles, and the Monarch Phenopacket Store (Danis et al., 2025) — 9,588 real published case-report phenopackets with confirmed OMIM diagnoses. On the real cohort, non-private IC-weighted cosine retrieval achieves MRR = 0.87 and nDCG@10 = 0.69, placing the system within the Phenomizer/LIRICAL band. Shadow-model membership-inference attack ROC AUC drops from 0.98 (no DP) to 0.50 (random) at ε ≤ 1; k-anonymity at k = 10 reduces re-identification probability against the rare-term adversary from 0.42 to 0.005.

The most consequential finding concerns deployment: the safe Laplace-DP budget on real patients is 20–50× larger than synthetic-cohort experiments suggest, because real similarity-score distributions are compressed and per-score Laplace noise dominates them. We propose the iterative exponential mechanism on rank utility as the principled response and validate it empirically: it recovers 90% of non-private nDCG@10 at ε = 5 versus 13% for Laplace at matched ε-DP — a 10× budget efficiency improvement. Our open-source implementation provides a reusable benchmark, a deployment configuration calibrated against the real-cohort findings, and an interactive pilot system at <https://honors-thesis-54tubqjkgwqjm4zyegglxw.streamlit.app/>.

**Keywords:** rare diseases, phenotype matching, differential privacy, private set intersection, k-anonymity, Human Phenotype Ontology, GA4GH Phenopackets, membership inference

---

## Acknowledgments

I am deeply grateful to my thesis advisor, **Professor Stephen B. Montgomery**, for his guidance, mentorship, and support throughout this project. His pioneering work on rare disease genomics and the application of transcriptomics to diagnosis has been a constant source of inspiration.

I thank the members of the Montgomery Lab for valuable discussions and feedback. I am also grateful to the broader rare disease research community, including the developers of the Human Phenotype Ontology, GA4GH Phenopackets, and Matchmaker Exchange, whose open standards made this work possible.

This research was supported by [funding sources]. Computational resources were provided by [computing resources].

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

- **Table 1:** Rare disease statistics
- **Table 2:** HPO annotations dataset summary
- **Table 3:** Evaluation cohort characteristics
- **Table 4:** Baseline retrieval performance
- **Table 5:** Privacy mechanism comparison
- **Table 6:** Computational overhead analysis

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
