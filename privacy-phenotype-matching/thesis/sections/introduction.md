# Introduction

## 1.1 The Rare-Disease Diagnostic Odyssey

Rare diseases — conditions affecting fewer than 200,000 individuals each in the United States — collectively reach 300 million patients worldwide across roughly 7,000 entities (Nguengang Wakap et al., 2020; NORD, 2023). The path to diagnosis is long: an average of 4.8–7 years, 7.3 physicians consulted, multiple misdiagnoses (EURORDIS, 2009; Global Genes, 2023). Eighty percent of cases are genetic (Ferreira, 2019); whole-exome and whole-genome sequencing now solve 25–40% of previously undiagnosed cases (Boycott et al., 2017), and RNA-sequencing adds an additional 7.5–17% on top of that yield (Frésard et al., 2019). The remaining undiagnosed majority — disproportionately patients with atypical or unique phenotypic constellations — is the population this thesis targets.

## 1.2 Federated Patient Matching

A clinician confronting an undiagnosed patient with an unusual phenotype profile gains diagnostic power by locating other patients with similar profiles, particularly those already diagnosed. The Matchmaker Exchange network (Philippakis et al., 2015) operationalises this insight across seven major rare-disease databases, has facilitated over a hundred disease-gene discoveries since 2015, and is built on two enabling standards: the Human Phenotype Ontology for vocabulary (Köhler et al., 2021) and GA4GH Phenopackets for record exchange (Jacobsen et al., 2022). The clinical case for federated matching is settled. The privacy case is not.

## 1.3 The Privacy Bottleneck

Phenotype data carries a privacy signature distinct from typical clinical records. Rare phenotype combinations function as quasi-identifiers — a patient with "Seizures," "Macrocephaly," and "Supernumerary nipple" may be uniquely identifiable in any database — and even Boolean variant-presence responses in GA4GH Beacons admit membership-inference attacks with around 5,000 queries (Shringarpure & Bustamante, 2015; El Emam et al., 2011). Institutional review boards consequently restrict cross-border data sharing, even where patients themselves would consent (Ramoni et al., 2017). The gap between the clinical value of federated matching and the institutional risk of phenotype exchange is the gap this thesis closes with rigorous, quantifiable privacy guarantees.

## 1.4 Three Privacy Primitives

Three established privacy-preserving computation techniques map naturally onto phenotype matching. **Private Set Intersection** (Meadows, 1986; Pinkas et al., 2018) computes shared phenotypes between two patients without revealing non-matching terms. **Differential privacy** (Dwork et al., 2006) adds calibrated noise so that the inclusion or exclusion of any single record is statistically undetectable in released outputs. **k-anonymity** (Sweeney, 2002) suppresses responses that would single out cohorts smaller than $k$ patients. No single primitive defends against every threat; their composition does. Chapter 3 makes this precise.

## 1.5 Contributions

This thesis closes the clinical-utility / institutional-risk gap above with five contributions, ordered by novelty:

1. **A formal threat model with three composition-level privacy invariants.** A two-party semi-honest model with auxiliary information; three adversary goals (membership inference, attribute inference, singling-out re-identification); per-step disclosure analysis; and proofs that the composition of DH-PSI, Laplace/Gaussian/Exponential DP, and a k-anonymity gate satisfies (I1) PSI semantic security, (I2) $(\varepsilon, \delta)$-DP score release, and (I3) k-anonymity of result release (§3.1.2).

2. **The synthetic-to-real DP gap and its measured fix.** We are aware of no prior work that benchmarks rare-disease privacy mechanisms on real published patients. Doing so reveals a 20–50× discrepancy in safe ε between synthetic and real cohorts, traced mechanically to compressed similarity-score distributions. The principled response — iterative exponential mechanism on a rank utility — recovers 90% of non-private nDCG@10 at ε = 5 versus 13% for Laplace under identical ε-DP guarantee, a 10× budget-efficiency advantage (§4.7).

3. **Empirically validated retrieval on real published patients.** Non-private IC-weighted cosine retrieval on 1,500 Phenopacket Store patients (Danis et al., 2025) achieves MRR = 0.87 and nDCG@10 = 0.69, placing the system within the Phenomizer/LIRICAL band (§4.6). A full Resnik+BMA Phenomizer-style baseline runs on the same corpus and is essentially tied with the simpler metric.

4. **Empirical privacy measurement.** Yeom-threshold and Shokri-shadow membership-inference attacks against the DP score-release oracle reduce attack ROC AUC from 0.98 to 0.50 (random) at ε ≤ 1; a singling-out attack against the k-anonymity gate falls from 0.42 to 0.005 re-identification probability between k = 1 and k = 10 (§4.5).

5. **A reproducible, deployable artefact.** GA4GH-compatible open-source implementation with a Docker- and Make-driven reproduction pipeline (every figure regenerated in roughly four minutes from a clean checkout) and an interactive pilot system deployed at <https://honors-thesis-54tubqjkgwqjm4zyegglxw.streamlit.app/> (§3.8).

## 1.6 Organisation

Chapter 2 surveys the relevant standards, similarity measures, federated systems, and privacy literature, naming four gaps the contributions above close. Chapter 3 specifies the threat model, data representation, similarity metrics, privacy mechanisms, and pilot system. Chapter 4 reports the synthetic- and real-cohort retrieval, attack, and rank-DP results. Chapter 5 interprets the findings — most notably the synthetic-to-real generalisation gap — and revises the deployment recommendation accordingly. Chapter 6 lists follow-on work; Chapter 7 concludes.
