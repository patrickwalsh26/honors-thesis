# Conclusion

This thesis closes the gap between the clinical value of federated rare-disease matching and the institutional barriers to phenotype sharing. The framework composes Private Set Intersection, differential privacy, and k-anonymity within a pipeline standards-compatible with GA4GH Phenopackets, Beacon v2, and Matchmaker Exchange; the central technical contribution is a measured fix — the rank-utility exponential mechanism — that recovers 90% of non-private retrieval utility at ε = 5 on real patients, where standard Laplace DP requires ε ≥ 50 for equivalent performance.

## Summary of Contributions

**1. Formal threat model and modular privacy framework.** We specified a two-party semi-honest threat model with auxiliary information, three concrete adversary goals, and per-step disclosure analysis. The composition of Diffie-Hellman PSI, differentially-private score release, and a k-anonymity gate is proven against three privacy invariants (PSI semantic security under DDH, $(\varepsilon,\delta)$-DP of released similarities, k-anonymity of result release).

**2. Empirically validated retrieval on real published patients.** On 1,500 Phenopacket Store patients across 100 OMIM diseases, non-private IC-weighted cosine retrieval achieves MRR = 0.87 and nDCG@10 = 0.69 — within the Phenomizer/LIRICAL band. The full Resnik+BMA Phenomizer-style baseline runs on the same corpus and is essentially tied with our simpler metric.

**3. Empirically measured privacy defense.** Shadow-model membership-inference attack ROC AUC collapses from 0.98 (no DP) to 0.50 (random) at ε ≤ 1; k-anonymity at k = 10 cuts re-identification probability against the rare-term singling-out adversary from 0.42 to 0.005.

**4. The synthetic-to-real privacy budget gap, and a principled fix.** The safe Laplace-DP budget on real patients is 20–50× larger than synthetic-cohort experiments suggest. We diagnose the cause as the compression of same-disease similarity-score gaps in real cohorts and propose the iterative exponential mechanism on rank utility as the principled response. Empirically (§4.7), the rank-utility mechanism recovers 90% of non-private nDCG@10 at ε = 5 on the real cohort versus 13% for Laplace at matched ε-DP — a 10× budget efficiency advantage.

**5. GA4GH standards compatibility and open-source release.** Native Phenopackets v2.0, HPO, Beacon v2, and MME adapters. The complete implementation (~4,500 lines of Python) ships with a Dockerfile and Makefile that regenerate every figure and table in ~4 minutes from a clean checkout.

## Broader Impact

The synthetic-to-real DP gap is not specific to phenotype matching. Any privacy-preserving retrieval task whose underlying signal is compressed — drug-response similarity for clinical-trial recruitment, gene-expression outlier matching for diagnosis, epidemiological-cluster discovery on coded encounter data — is vulnerable to the same Laplace pathology and benefits from the same rank-utility fix. The framing this thesis offers is therefore generalisable: in retrieval-style privacy problems, the utility function for the underlying DP mechanism is a more consequential design decision than the mechanism family.

## Closing Remarks

Rare-disease diagnosis is fundamentally a collective inference problem; the 300 million patients living with these conditions have been ill-served by privacy regimes that treat data-sharing prohibition as the default and individual-institutional compliance as the means. The contribution of this thesis is to demonstrate, with measured numbers on real published patients, that a quantifiably private alternative exists and is deployable today. Larger cohorts, EHR-derived data, multi-omic extension, malicious-secure cryptographic protocols, and prospective clinical validation remain — but the foundational privacy-utility question now has a defensible empirical answer.

---

*"Alone we can do so little; together we can do so much."*
— Helen Keller
