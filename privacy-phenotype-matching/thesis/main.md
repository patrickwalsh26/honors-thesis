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

Rare diseases affect approximately 300 million people worldwide, yet patients endure an average diagnostic odyssey of 4.8-7 years due to the scarcity of clinical expertise and fragmented medical knowledge. Federated patient matching—connecting individuals with similar phenotypic presentations across institutions—offers a promising avenue for accelerating diagnosis. However, sharing sensitive phenotype data raises significant privacy concerns, as rare phenotype combinations can serve as quasi-identifiers enabling patient re-identification.

We present a privacy-preserving phenotype matching framework that enables federated rare disease cohort discovery while protecting patient confidentiality. Our system integrates three complementary privacy mechanisms: Private Set Intersection (PSI) for secure phenotype overlap computation, differential privacy (DP) for quantifiable information leakage bounds, and k-anonymity for rare term filtering. The framework operates on GA4GH Phenopackets with Human Phenotype Ontology (HPO) terms, ensuring compatibility with established genomic data sharing standards.

We evaluate our approach using synthetic patients generated from 12,974 real disease profiles in the HPO annotation corpus, encompassing phenotype associations from OMIM and Orphanet. Baseline retrieval experiments achieve strong performance (nDCG@10 = 99.7%, MRR = 1.0) using information content-weighted Resnik similarity. We characterize privacy-utility tradeoffs across mechanism configurations and empirically measure privacy protection through membership and attribute inference attack simulations.

Our open-source implementation demonstrates that privacy-preserving phenotype matching is both technically feasible and practically useful, offering a path toward broader federated collaboration in rare disease research without compromising patient confidentiality.

**Keywords:** rare diseases, phenotype matching, privacy-preserving computation, private set intersection, differential privacy, Human Phenotype Ontology, GA4GH Phenopackets

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
- **Figure 6:** nDCG@k vs. ε (differential privacy tradeoff)
- **Figure 7:** Suppression rate vs. k (k-anonymity tradeoff)
- **Figure 8:** Privacy-utility Pareto frontier
- **Figure 9:** Membership inference attack ROC curves

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
