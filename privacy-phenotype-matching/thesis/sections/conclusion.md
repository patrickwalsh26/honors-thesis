# Conclusion

Rare diseases collectively affect 300 million people worldwide, yet individual patients often wait years for diagnosis due to the scarcity of clinical expertise and fragmented medical knowledge. Federated patient matching—connecting individuals with similar phenotypic presentations across institutions—offers a promising path to accelerate diagnosis. However, the sensitive nature of phenotype data has limited adoption, as institutions hesitate to share patient information across organizational boundaries.

This thesis demonstrates that privacy-preserving phenotype matching is both technically feasible and practically useful. We developed a modular framework integrating three complementary privacy mechanisms—Private Set Intersection, differential privacy, and k-anonymity—within a pipeline compatible with GA4GH Phenopackets, Beacon v2, and Matchmaker Exchange standards.

## Summary of Contributions

**1. Privacy-Preserving Framework.** We designed and implemented a composable privacy pipeline where each mechanism addresses distinct threat vectors. PSI enables secure phenotype overlap computation without revealing non-matching terms. Differential privacy provides rigorous, quantifiable bounds on information leakage. K-anonymity and rare term filtering protect against quasi-identifier inference from rare phenotype combinations. These mechanisms compose naturally, providing defense in depth.

**2. Semantic Phenotype Similarity.** We implemented multiple similarity metrics spanning set-theoretic (Jaccard), vector-space (Cosine), and semantic (Resnik) approaches. Information content weighting, derived from the HPO annotation corpus, emphasizes rare, diagnostically specific phenotypes. Our evaluation confirms that IC-weighted metrics outperform unweighted alternatives, achieving MRR = 1.000 and nDCG@10 = 99.7% on disease-sharing patient retrieval.

**3. GA4GH Standards Compatibility.** Our implementation natively supports Phenopackets v2.0 for patient representation, HPO for phenotype encoding, and adapters for Beacon and MME interfaces. This compatibility positions the framework for integration with existing rare disease infrastructure, lowering adoption barriers.

**4. Comprehensive Evaluation.** Using synthetic patients generated from 12,974 real disease-phenotype associations in OMIM and Orphanet, we characterized privacy-utility tradeoffs across mechanism configurations. Key findings include:
- Moderate privacy parameters (ε = 5.0, k = 5, 1% filtering) preserve 96.5% of baseline utility
- Privacy costs increase smoothly with protection strength, enabling informed parameter selection
- Mechanism composition provides layered protection without compounding utility loss

**5. Open-Source Release.** Our complete implementation—approximately 4,500 lines of Python code—is released under an open-source license. The release includes similarity metrics, privacy mechanisms, evaluation scripts, and synthetic data generators, enabling reproducibility and extension by the research community.

## Broader Impact

This work contributes to the broader goal of enabling collaborative rare disease research while respecting patient privacy. By demonstrating that meaningful matching is achievable with quantifiable privacy guarantees, we hope to encourage participation by institutions currently excluded from federated networks due to data sharing concerns.

The technical approach generalizes beyond phenotype matching. The composition of PSI, DP, and k-anonymity could apply to other biomedical matching problems—drug response similarity, clinical trial recruitment, and epidemiological surveillance—where sensitive data must be compared across institutional boundaries.

More broadly, this thesis illustrates how privacy-enhancing technologies can expand the frontier of what is possible in biomedical research. Rather than treating privacy as an obstacle to be minimized, we approach it as a design constraint that, when satisfied, enables collaboration that would otherwise be impossible.

## Closing Remarks

The 300 million people living with rare diseases deserve the benefits of collaborative research. Privacy-preserving computation offers a path forward—protecting individual confidentiality while unlocking the collective knowledge distributed across institutions worldwide. This thesis represents one step on that path.

Much work remains. Validation on real patient cohorts, extension to multi-omic data, integration with clinical workflows, and advancement of interoperability standards are all essential for translation to practice. We hope this work provides a foundation for these efforts and contributes to the ultimate goal: ending the diagnostic odyssey for rare disease patients everywhere.

---

*"Alone we can do so little; together we can do so much."*
— Helen Keller
