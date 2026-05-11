# Privacy-Preserving Phenotype Matching for Rare Disease Cohort Discovery

## Honors Thesis Outline

**Author:** Patrick Walsh
**Institution:** Stanford University
**Department:** Computer Science
**Date:** April 2026

---

## Document Structure

1. **Abstract** (~250 words)
2. **Introduction** (~2,000 words)
3. **Literature Review** (~3,000 words)
4. **Methods** (~4,000 words)
5. **Results** (~2,500 words)
6. **Discussion** (~2,000 words)
7. **Future Work** (~1,000 words)
8. **Conclusion** (~500 words)
9. **Bibliography** (~50-70 references)
10. **Appendices**

---

## 1. Abstract

**Key Points to Cover:**
- Problem: Rare disease diagnosis takes 4.8-7 years on average; federated patient matching could accelerate diagnosis but raises privacy concerns
- Approach: Privacy-preserving phenotype matching using Private Set Intersection (PSI), differential privacy (DP), and k-anonymity
- Methods: Implemented system using GA4GH Phenopackets, HPO ontology, and multiple similarity metrics
- Results: Achieved 99.7% NDCG@10 on baseline retrieval; quantified privacy-utility tradeoffs across mechanisms
- Conclusion: Demonstrates feasibility of privacy-preserving rare disease matching with acceptable utility loss

---

## 2. Introduction

### 2.1 The Rare Disease Diagnostic Challenge
- Definition: diseases affecting <200,000 individuals in the US (Orphan Drug Act, 1983)
- Scale: ~7,000 known rare diseases affecting 300 million people globally [Nguengang Wakap et al., 2020]
- Diagnostic odyssey: average 4.8 years to diagnosis [EURORDIS, 2009]; up to 7 years in complex cases [Global Genes, 2023]
- 80% of rare diseases have genetic etiology [Ferreira, 2019]
- Unmet need: 95% of rare diseases lack FDA-approved treatments

### 2.2 The Promise of Federated Patient Matching
- Similar patients may share diagnoses → accelerated diagnosis
- Matchmaker Exchange: GA4GH federated matching network [Philippakis et al., 2015]
- Current implementations: GeneMatcher, PhenomeCentral, DECIPHER
- Problem: requires sharing sensitive phenotype data across institutions

### 2.3 Privacy Concerns in Medical Data Sharing
- HIPAA and GDPR constraints on health data
- Phenotype data as quasi-identifiers [El Emam et al., 2011]
- Re-identification risks from rare phenotype combinations
- Institutional reluctance to share patient data

### 2.4 Thesis Contributions
1. **Framework**: Privacy-preserving phenotype matching pipeline combining PSI, DP, and k-anonymity
2. **Implementation**: Open-source system compatible with GA4GH standards (Phenopackets, Beacon v2, MME)
3. **Evaluation**: Comprehensive privacy-utility analysis on real disease data from OMIM/Orphanet
4. **Insights**: Quantified tradeoffs and practical recommendations for deployment

### 2.5 Thesis Organization
- Overview of remaining sections

---

## 3. Literature Review

### 3.1 Rare Disease Phenotyping and Standards
- Human Phenotype Ontology (HPO): 18,000+ terms [Kohler et al., 2021]
- GA4GH Phenopackets v2.0: standardized patient representation [Jacobsen et al., 2022]
- OMIM: 17,000+ disease entries [Amberger et al., 2019]
- Orphanet: European rare disease reference [Rath et al., 2012]
- DECIPHER: 44,000+ patients with genomic variants [Firth et al., 2009]

### 3.2 Phenotype Similarity Metrics
- Set-based: Jaccard index, Dice coefficient
- Information-theoretic: Resnik similarity [Resnik, 1995], Lin similarity [Lin, 1998]
- Graph-based: Wang semantic similarity [Wang et al., 2007]
- Corpus-based IC: annotation-derived information content [Sanchez et al., 2011]
- Phenomizer and related tools [Kohler et al., 2009]

### 3.3 Patient Matching Systems
- Matchmaker Exchange architecture [Philippakis et al., 2015]
- GeneMatcher: gene-centric matching [Sobreira et al., 2015]
- PhenomeCentral: phenotype-first matching [Buske et al., 2015]
- Beacon Network: count/boolean queries [Fiume et al., 2019]

### 3.4 Privacy-Preserving Computation
#### 3.4.1 Private Set Intersection (PSI)
- Diffie-Hellman PSI [Meadows, 1986; Huberman et al., 1999]
- OT-based PSI [Pinkas et al., 2014]
- Circuit-based PSI [Huang et al., 2012]
- PSI cardinality [De Cristofaro & Tsudik, 2010]
- Recent advances: unbalanced PSI, threshold PSI [Chen et al., 2017]

#### 3.4.2 Differential Privacy
- Definition and properties [Dwork et al., 2006]
- Laplace mechanism for numeric queries [Dwork et al., 2006]
- Gaussian mechanism for (ε,δ)-DP [Dwork & Roth, 2014]
- Exponential mechanism for selection [McSherry & Talwar, 2007]
- Composition theorems [Kairouz et al., 2015]
- DP in healthcare [Dankar & El Emam, 2013]

#### 3.4.3 K-Anonymity and Generalization
- k-Anonymity definition [Sweeney, 2002]
- l-Diversity and t-closeness [Machanavajjhala et al., 2007; Li et al., 2007]
- Phenotype generalization via HPO hierarchy
- Rare term filtering for quasi-identifier protection

### 3.5 Privacy Attacks and Auditing
- Membership inference attacks [Shokri et al., 2017]
- Attribute inference [Yeom et al., 2018]
- Reconstruction attacks [Dinur & Nissim, 2003]
- Privacy auditing frameworks [Jagielski et al., 2020]

### 3.6 Related Privacy-Preserving Genomic Systems
- Secure genome analysis [Ayday et al., 2013]
- Privacy-preserving GWAS [Johnson & Shmatikov, 2013]
- Beacon re-identification attacks [Shringarpure & Bustamante, 2015]
- DP Beacons [Raisaro et al., 2017]

---

## 4. Methods

### 4.1 System Architecture
- Overview diagram of privacy-preserving pipeline
- Data flow: phenopacket → privacy mechanisms → similarity → results
- Modular design allowing mechanism composition

### 4.2 Data Representation
#### 4.2.1 GA4GH Phenopackets
- Schema overview (subject, phenotypicFeatures, diseases, metaData)
- HPO term encoding with observed/excluded status
- Conversion from clinical records

#### 4.2.2 HPO Ontology Integration
- Ontology structure (DAG with 18,000+ terms)
- Ancestor/descendant relationships
- Information content computation
- Term prevalence from annotation corpora

### 4.3 Similarity Metrics
#### 4.3.1 Jaccard Similarity
$$J(A, B) = \frac{|A \cap B|}{|A \cup B|}$$

#### 4.3.2 Information Content (IC)
$$IC(t) = -\log_2 P(t)$$
- Corpus-based: IC from phenopacket prevalence
- Intrinsic: IC from ontology structure (descendants)

#### 4.3.3 Resnik Similarity
$$sim_{Resnik}(t_1, t_2) = IC(MICA(t_1, t_2))$$
- Best-Match Average (BMA) for term sets:
$$BMA(A, B) = \frac{1}{2}\left(\frac{1}{|A|}\sum_{a \in A} \max_{b \in B} sim(a,b) + \frac{1}{|B|}\sum_{b \in B} \max_{a \in A} sim(a,b)\right)$$

#### 4.3.4 Cosine Similarity with IC Weighting
$$cos(A, B) = \frac{\sum_{t \in A \cap B} IC(t)^2}{\sqrt{\sum_{t \in A} IC(t)^2} \cdot \sqrt{\sum_{t \in B} IC(t)^2}}$$

### 4.4 Privacy Mechanisms
#### 4.4.1 Private Set Intersection (PSI)
- Diffie-Hellman protocol on elliptic curves (NIST P-256)
- Protocol steps: hash-to-curve, scalar multiplication, intersection
- Security guarantees: semi-honest model
- Hybrid PSI with ontology expansion

#### 4.4.2 Differential Privacy
- Laplace mechanism: noise ~ Lap(Δf/ε)
- Gaussian mechanism: noise ~ N(0, σ²) where σ = Δf√(2ln(1.25/δ))/ε
- Exponential mechanism for top-k selection
- Privacy accounting and composition

#### 4.4.3 K-Anonymity via Rare Term Filtering
- Prevalence threshold: filter terms appearing in <k patients
- Strategies: suppress (remove) or generalize (parent term)
- Quasi-identifier protection for unique phenotype combinations

#### 4.4.4 Composed Privacy Pipeline
- Order of operations: rare term filtering → PSI → DP → k-anonymity check
- Configuration via YAML with tunable parameters
- Privacy budget tracking

### 4.5 Evaluation Framework
#### 4.5.1 Retrieval Metrics
- Precision@k, Recall@k, F1@k
- Normalized Discounted Cumulative Gain (nDCG@k)
- Mean Average Precision (MAP)
- Mean Reciprocal Rank (MRR)

#### 4.5.2 Privacy Metrics
- Membership inference attack success rate
- Attribute inference accuracy
- Adversarial advantage over random guessing

#### 4.5.3 Privacy-Utility Frontier
- Sweep over ε values for DP
- Sweep over k values for k-anonymity
- Pareto frontier analysis

### 4.6 Dataset Construction
#### 4.6.1 HPO Annotations (phenotype.hpoa)
- Source: HPO consortium (http://purl.obolibrary.org/obo/hp/hpoa/)
- 12,974 diseases (8,592 OMIM, 4,335 Orphanet, 47 DECIPHER)
- 282,728 disease-phenotype annotations
- 11,514 unique HPO terms

#### 4.6.2 Synthetic Patient Generation
- Sample diseases from HPOA
- Generate patients with phenotype recall (75%) and noise (10%)
- Ground truth: patients with same underlying disease are relevant
- Cohort: 500 patients from 100 diseases

### 4.7 Implementation Details
- Language: Python 3.10+
- Cryptography: `cryptography` library (ECDH, HKDF)
- Ontology: `pronto` for HPO parsing
- Evaluation: NumPy, SciPy for metrics
- Standards: GA4GH Phenopackets v2.0, Beacon v2, MME v1.1

---

## 5. Results

### 5.1 Dataset Characteristics
| Statistic | Value |
|-----------|-------|
| Diseases in HPOA | 12,974 |
| OMIM diseases | 8,592 (66.2%) |
| Orphanet diseases | 4,335 (33.4%) |
| DECIPHER diseases | 47 (0.4%) |
| Unique HPO terms | 11,514 |
| Mean phenotypes/disease | 20.3 |

**Evaluation Cohort:**
| Statistic | Value |
|-----------|-------|
| Patients | 500 |
| Diseases | 100 |
| Unique phenotypes | 1,145 |
| Phenotypes/patient | 10.4 ± 5.9 |

### 5.2 Baseline Similarity Performance
| Metric | P@1 | P@5 | R@5 | R@10 | NDCG@10 | MRR |
|--------|-----|-----|-----|------|---------|-----|
| Jaccard | 99.6% | 79.2% | 99.0% | 99.9% | 99.5% | 0.998 |
| Cosine (IC) | 99.8% | 79.2% | 99.1% | 99.9% | 99.6% | 0.999 |
| Resnik (simplified) | **100%** | **79.3%** | **99.2%** | **100%** | **99.7%** | **1.000** |

**Key Finding:** IC-weighted Resnik similarity achieves perfect MRR, demonstrating that information content weighting improves ranking of phenotypically similar patients.

### 5.3 Similarity Score Distributions
| Metric | Mean | Std | Median | Max |
|--------|------|-----|--------|-----|
| Jaccard | 0.010 | 0.053 | 0.0 | 1.0 |
| Cosine (IC) | 0.011 | 0.063 | 0.0 | 1.0 |
| Resnik | 0.008 | 0.050 | 0.0 | 1.0 |

**Interpretation:** Highly sparse similarity distributions—most patient pairs have zero overlap, making retrieval challenging and privacy protection critical.

### 5.4 Privacy-Utility Tradeoffs
#### 5.4.1 Differential Privacy Impact
- Table/figure showing NDCG@10 vs ε for ε ∈ {0.1, 0.5, 1.0, 2.0, 5.0, 10.0}
- Expected: graceful degradation with lower ε
- Identify practical operating point (e.g., ε=1.0)

#### 5.4.2 K-Anonymity Impact
- Table/figure showing precision/recall vs k for k ∈ {2, 5, 10, 20}
- Suppression rates at each k
- Impact on result availability

#### 5.4.3 Combined Mechanisms
- Performance with PSI + DP + k-anonymity
- Comparison to individual mechanisms
- Privacy budget consumption

### 5.5 Privacy Leakage Analysis
- Membership inference attack accuracy (baseline vs. protected)
- Attribute inference success rates
- Adversarial advantage reduction with privacy mechanisms

### 5.6 Computational Performance
- PSI protocol overhead (ms per comparison)
- Similarity matrix computation time
- Scalability analysis

---

## 6. Discussion

### 6.1 Interpretation of Results
- High baseline performance validates phenotype-based matching
- IC weighting provides measurable improvement
- Privacy mechanisms achieve protection with quantifiable utility cost

### 6.2 Privacy-Utility Tradeoffs
- DP provides strongest formal guarantees but highest utility cost
- k-Anonymity is practical but may suppress valid results
- PSI protects intersection but reveals cardinality
- Combined approach balances protections

### 6.3 Practical Deployment Considerations
- Recommended parameter settings (ε ≈ 1-2, k ≈ 5)
- Institutional policy alignment
- Integration with existing MME infrastructure

### 6.4 Comparison to Related Work
- Beacon DP solutions [Raisaro et al., 2017]
- Secure genomic computation [Ayday et al., 2013]
- Advantages of phenotype-level matching

### 6.5 Limitations
- Synthetic patient generation may not capture real-world noise
- Semi-honest adversary model for PSI
- Composition effects with repeated queries
- Single-institution evaluation

### 6.6 Ethical Considerations
- Balancing research advancement with patient privacy
- Informed consent for federated matching
- Equity in rare disease access

---

## 7. Future Work

### 7.1 Technical Extensions
- Malicious-secure PSI protocols
- Homomorphic encryption for similarity computation
- Federated learning for model-based matching
- Secure multi-party computation (MPC) for multi-institution queries

### 7.2 System Improvements
- Real-time privacy budget management
- Adaptive mechanism selection based on query sensitivity
- Integration with clinical decision support systems
- Mobile/edge deployment for resource-constrained settings

### 7.3 Validation Studies
- Evaluation on real patient cohorts (with IRB approval)
- Multi-site deployment pilot
- User studies with clinicians
- Longitudinal diagnostic outcome tracking

### 7.4 Standards and Interoperability
- GA4GH privacy extension proposals
- Integration with Beacon v2 handshakes
- Support for additional phenotype ontologies (e.g., SNOMED-CT)

---

## 8. Conclusion

**Summary of Contributions:**
1. Designed and implemented privacy-preserving phenotype matching framework
2. Achieved strong retrieval performance (NDCG@10 > 99%) on real disease data
3. Quantified privacy-utility tradeoffs across PSI, DP, and k-anonymity
4. Provided practical recommendations for deployment

**Broader Impact:**
- Enables federated rare disease matching while protecting patient privacy
- Accelerates diagnostic odyssey without compromising confidentiality
- Open-source implementation promotes reproducibility and adoption

---

## 9. Bibliography

### Key References (to be expanded with full citations)

**Rare Disease & Phenotyping:**
1. Kohler S, et al. (2021). The Human Phenotype Ontology in 2021. *Nucleic Acids Research*, 49(D1), D1207-D1217.
2. Jacobsen JOB, et al. (2022). The GA4GH Phenopacket schema defines a computable representation of clinical data. *Nature Biotechnology*, 40, 817-820.
3. Nguengang Wakap S, et al. (2020). Estimating cumulative point prevalence of rare diseases. *European Journal of Human Genetics*, 28, 165-173.
4. Amberger JS, et al. (2019). OMIM.org: leveraging knowledge across phenotype-gene relationships. *Nucleic Acids Research*, 47(D1), D1038-D1043.
5. Firth HV, et al. (2009). DECIPHER: Database of Chromosomal Imbalance and Phenotype in Humans. *American Journal of Human Genetics*, 84(4), 524-533.

**Patient Matching:**
6. Philippakis AA, et al. (2015). The Matchmaker Exchange: A Platform for Rare Disease Gene Discovery. *Human Mutation*, 36(10), 915-921.
7. Sobreira N, et al. (2015). GeneMatcher: A Matching Tool for Connecting Investigators. *Human Mutation*, 36(10), 928-930.
8. Buske OJ, et al. (2015). PhenomeCentral: A Portal for Phenotypic and Genotypic Matchmaking. *Human Mutation*, 36(10), 931-940.
9. Fiume M, et al. (2019). Federated discovery and sharing of genomic data using Beacons. *Nature Biotechnology*, 37, 220-224.

**Semantic Similarity:**
10. Resnik P. (1995). Using Information Content to Evaluate Semantic Similarity. *IJCAI*, 448-453.
11. Lin D. (1998). An Information-Theoretic Definition of Similarity. *ICML*, 296-304.
12. Kohler S, et al. (2009). Clinical diagnostics in human genetics with semantic similarity searches. *American Journal of Human Genetics*, 85(4), 457-464.

**Privacy-Preserving Computation:**
13. Dwork C, et al. (2006). Calibrating Noise to Sensitivity in Private Data Analysis. *TCC*, 265-284.
14. Dwork C, Roth A. (2014). The Algorithmic Foundations of Differential Privacy. *Foundations and Trends in Theoretical Computer Science*, 9(3-4), 211-407.
15. Sweeney L. (2002). k-Anonymity: A Model for Protecting Privacy. *International Journal of Uncertainty, Fuzziness and Knowledge-Based Systems*, 10(5), 557-570.
16. Pinkas B, et al. (2014). Faster Private Set Intersection Based on OT Extension. *USENIX Security*, 797-812.
17. De Cristofaro E, Tsudik G. (2010). Practical Private Set Intersection Protocols. *Financial Cryptography*, 250-268.

**Privacy Attacks:**
18. Shokri R, et al. (2017). Membership Inference Attacks Against Machine Learning Models. *IEEE S&P*, 3-18.
19. Shringarpure SS, Bustamante CD. (2015). Privacy Risks from Genomic Data-Sharing Beacons. *American Journal of Human Genetics*, 97(5), 631-646.
20. Raisaro JL, et al. (2017). Addressing Beacon re-identification attacks. *Genome Research*, 27(5), 798-809.

**Healthcare Privacy:**
21. El Emam K, et al. (2011). A Systematic Review of Re-Identification Attacks on Health Data. *PLoS ONE*, 6(12), e28071.
22. Dankar FK, El Emam K. (2013). Practicing Differential Privacy in Health Care. *Transactions on Data Privacy*, 6(1), 35-67.

*(Full bibliography will contain 50-70 references)*

---

## 10. Appendices

### Appendix A: HPO Term Examples
- Sample rare disease phenotype profiles
- IC distribution analysis

### Appendix B: Algorithm Pseudocode
- PSI protocol steps
- Privacy-preserving similarity computation

### Appendix C: Full Experimental Results
- Complete metric tables across all parameter settings
- Additional visualizations

### Appendix D: Software Documentation
- Installation instructions
- API reference
- Configuration options

---

## Figures and Tables (Planned)

### Figures
1. **Figure 1:** System architecture diagram
2. **Figure 2:** HPO ontology structure (example subgraph)
3. **Figure 3:** Privacy-preserving pipeline workflow
4. **Figure 4:** Similarity score distributions
5. **Figure 5:** Precision-Recall curves by metric
6. **Figure 6:** NDCG@k vs. ε (differential privacy tradeoff)
7. **Figure 7:** Suppression rate vs. k (k-anonymity tradeoff)
8. **Figure 8:** Privacy-utility Pareto frontier
9. **Figure 9:** Membership inference attack ROC curves

### Tables
1. **Table 1:** Rare disease statistics (prevalence, diagnostic delay)
2. **Table 2:** HPO annotations dataset summary
3. **Table 3:** Evaluation cohort characteristics
4. **Table 4:** Baseline retrieval performance
5. **Table 5:** Privacy mechanism comparison
6. **Table 6:** Computational overhead analysis

---

## Writing Schedule (Suggested)

| Section | Target Length | Priority |
|---------|---------------|----------|
| Methods | 4,000 words | High |
| Results | 2,500 words | High |
| Introduction | 2,000 words | High |
| Literature Review | 3,000 words | Medium |
| Discussion | 2,000 words | Medium |
| Future Work | 1,000 words | Low |
| Abstract | 250 words | Last |

---

*This outline provides the structure for a comprehensive honors thesis. Each section will be drafted with proper academic citations, real experimental data, and clear technical exposition.*
