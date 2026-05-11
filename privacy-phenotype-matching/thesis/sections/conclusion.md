# Conclusion

Rare diseases collectively affect 300 million people worldwide, yet individual patients often wait years for diagnosis due to the scarcity of clinical expertise and fragmented medical knowledge. Federated patient matching—connecting individuals with similar phenotypic presentations across institutions—offers a promising path to accelerate diagnosis. However, the sensitive nature of phenotype data has limited adoption, as institutions hesitate to share patient information across organizational boundaries.

This thesis demonstrates that privacy-preserving phenotype matching is both technically feasible and practically useful. We developed a modular framework integrating three complementary privacy mechanisms—Private Set Intersection, differential privacy, and k-anonymity—within a pipeline compatible with GA4GH Phenopackets, Beacon v2, and Matchmaker Exchange standards.

## Summary of Contributions

**1. Formal threat model and modular privacy framework.** We specified a two-party semi-honest threat model with auxiliary information, three concrete adversary goals, and per-step disclosure analysis. The composition of Diffie-Hellman PSI, differentially-private score release, and a k-anonymity gate is proven against three privacy invariants (PSI semantic security under DDH, $(\varepsilon,\delta)$-DP of released similarities, k-anonymity of result release).

**2. Empirically validated retrieval on real published patients.** On 1,500 Phenopacket Store patients across 100 OMIM diseases, non-private IC-weighted cosine retrieval achieves MRR = 0.87 and nDCG@10 = 0.69 — within the Phenomizer/LIRICAL band. The full Resnik+BMA Phenomizer-style baseline runs on the same corpus and is essentially tied with our simpler metric.

**3. Empirically measured privacy defense.** Shadow-model membership-inference attack ROC AUC collapses from 0.98 (no DP) to 0.50 (random) at ε ≤ 1; k-anonymity at k = 10 cuts re-identification probability against the rare-term singling-out adversary from 0.42 to 0.005.

**4. The synthetic-to-real privacy budget gap.** The safe Laplace-DP budget on real patients is 20–50× larger than synthetic-cohort experiments suggest. Per-score Laplace noise is the wrong choice for real cohorts because real similarity-score distributions are compressed; rank-based mechanisms (Report-Noisy-Max, Exponential) are the principled response. Privacy-utility claims grounded only in synthetic experiments should be treated as upper bounds on real-world performance.

**5. GA4GH standards compatibility and open-source release.** Native Phenopackets v2.0, HPO, Beacon v2, and MME adapters. The complete implementation (~4,500 lines of Python) ships with a Dockerfile and Makefile that regenerate every figure and table in ~4 minutes from a clean checkout.

## Broader Impact

This work contributes to the broader goal of enabling collaborative rare disease research while respecting patient privacy. By demonstrating that meaningful matching is achievable with quantifiable privacy guarantees, we hope to encourage participation by institutions currently excluded from federated networks due to data sharing concerns.

The technical approach generalizes beyond phenotype matching. The composition of PSI, DP, and k-anonymity could apply to other biomedical matching problems—drug response similarity, clinical trial recruitment, and epidemiological surveillance—where sensitive data must be compared across institutional boundaries.

More broadly, this thesis illustrates how privacy-enhancing technologies can expand the frontier of what is possible in biomedical research. Rather than treating privacy as an obstacle to be minimized, we approach it as a design constraint that, when satisfied, enables collaboration that would otherwise be impossible.

## Closing Remarks

The 300 million people living with rare diseases deserve the benefits of collaborative research. Privacy-preserving computation offers a path forward—protecting individual confidentiality while unlocking the collective knowledge distributed across institutions worldwide. This thesis represents one step on that path.

Much work remains. Rank-based DP mechanisms (Report-Noisy-Max, Exponential) need to be evaluated against the same real cohort to determine whether they recover utility at the ε values infeasible for Laplace. Larger and more heterogeneous cohorts — including EHR-derived data which would further stress the synthetic-to-real gap — would refine the privacy-utility curves. Multi-omic extension, malicious-secure cryptographic protocols, and integration with clinical workflows are all essential next steps for translation to practice. We hope this work provides a foundation for these efforts and contributes to the ultimate goal: ending the diagnostic odyssey for rare disease patients everywhere.

---

*"Alone we can do so little; together we can do so much."*
— Helen Keller
