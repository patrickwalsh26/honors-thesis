# Introduction

## 1.1 The Rare Disease Diagnostic Challenge

Rare diseases, defined in the United States as conditions affecting fewer than 200,000 individuals, collectively represent one of the most significant challenges in modern medicine. Despite the low prevalence of any single rare disease, the cumulative burden is substantial: an estimated 300 million people worldwide—approximately 4% of the global population—live with one of more than 7,000 recognized rare diseases (Nguengang Wakap et al., 2020). In the United States alone, rare diseases affect 25-30 million individuals, making this patient population larger than those affected by heart disease or diabetes (NORD, 2023).

The path to diagnosis for rare disease patients is notoriously prolonged and arduous. Studies consistently report an average diagnostic delay of 4.8 to 7 years, during which patients consult an average of 7.3 physicians and frequently receive multiple misdiagnoses (EURORDIS, 2009; Global Genes, 2023). This extended period of uncertainty—often termed the "diagnostic odyssey"—carries profound consequences: delayed or inappropriate treatments, psychological distress for patients and families, and substantial healthcare costs estimated at $5 million per undiagnosed patient over their lifetime (Spillmann et al., 2017).

The diagnostic challenge stems from multiple factors. Approximately 80% of rare diseases have a genetic basis (Ferreira, 2019), yet the genetic variants responsible remain unknown for many conditions. Clinical presentations are often heterogeneous, with the same genetic variant producing different phenotypes across patients, and conversely, similar phenotypes arising from distinct genetic etiologies. Furthermore, individual clinicians may encounter only a handful of cases throughout their careers, limiting the development of clinical expertise for any particular rare condition.

Recent advances in genomic technologies have transformed rare disease diagnosis. Whole-exome sequencing (WES) and whole-genome sequencing (WGS) now achieve diagnostic yields of 25-40% for previously undiagnosed patients (Boycott et al., 2017). Complementary approaches leveraging transcriptomics have further improved these rates, with RNA sequencing from accessible tissues identifying causal variants in an additional 7.5-17% of cases where DNA sequencing alone was uninformative (Frésard et al., 2019; Montgomery Lab, Stanford). The GREGoR (Genomics Research to Elucidate the Genetics of Rare Diseases) Consortium, a multi-center initiative studying thousands of challenging rare disease cases, exemplifies the power of integrating multiple molecular data types for diagnosis (GREGoR Consortium, 2025).

## 1.2 The Promise of Federated Patient Matching

A critical insight driving recent progress in rare disease research is that similar patients may share similar diagnoses. When a clinician encounters a patient with an unusual constellation of phenotypes, identifying other patients with matching clinical presentations—particularly those who have received a molecular diagnosis—can dramatically accelerate the diagnostic process. This principle underlies the concept of federated patient matching: connecting patients across institutions based on phenotypic and genotypic similarity.

The Matchmaker Exchange (MME), established in 2015 under the auspices of the Global Alliance for Genomics and Health (GA4GH), represents the most successful implementation of this concept (Philippakis et al., 2015). MME provides a federated network connecting seven major rare disease databases—including GeneMatcher, PhenomeCentral, DECIPHER, and others—enabling clinicians to query across institutional boundaries for patients with matching clinical features. Since its inception, MME has facilitated thousands of successful matches, directly contributing to novel disease gene discoveries and patient diagnoses (Sobreira et al., 2015; Buske et al., 2015).

The Human Phenotype Ontology (HPO) provides the standardized vocabulary that makes such matching possible (Köhler et al., 2021). With over 18,000 terms describing phenotypic abnormalities, HPO enables precise, computable representation of clinical features. The ontology's hierarchical structure—organized as a directed acyclic graph—supports both exact matching and semantic similarity computation, where related phenotypes can be recognized as similar even without identical terminology.

The GA4GH Phenopacket standard further advances interoperability by defining a comprehensive schema for representing patient clinical data (Jacobsen et al., 2022). Phenopackets encapsulate not only phenotypic features but also diagnoses, genetic variants, and metadata in a structured, machine-readable format. These standards, developed through international collaboration, form the technical foundation for federated rare disease research.

## 1.3 Privacy Challenges in Genomic Data Sharing

Despite the clear clinical benefits of patient matching, sharing phenotype and genotype data across institutional boundaries raises significant privacy concerns. Medical data is among the most sensitive categories of personal information, protected by regulations including the Health Insurance Portability and Accountability Act (HIPAA) in the United States and the General Data Protection Regulation (GDPR) in Europe.

Phenotype data presents particular privacy risks that distinguish it from traditional medical records. Unlike laboratory values or vital signs, phenotypic features—especially rare combinations—can serve as quasi-identifiers that uniquely fingerprint individuals (El Emam et al., 2011). A patient with the combination of "Seizures," "Macrocephaly," and "Supernumerary nipple" may be the only such individual in a database, rendering them identifiable even without explicit identifiers. The GTEx Consortium's comprehensive analysis of genetic effects on gene expression across 44 human tissues demonstrated that individual genetic variation creates distinct molecular signatures (GTEx Consortium, 2017), further highlighting the identifiability risks inherent in detailed molecular characterization.

The Beacon Network, another GA4GH initiative enabling queries about the presence of specific genetic variants, illustrates these risks concretely. Shringarpure and Bustamante (2015) demonstrated that even simple Boolean responses ("yes, this variant exists in our database") could enable re-identification attacks with as few as 5,000 queries. Subsequent work by Raisaro et al. (2017) introduced differential privacy protections for Beacon responses, but the fundamental tension between utility and privacy persists.

These concerns have tangible consequences for research. Institutional review boards and data governance committees often restrict data sharing to protect patient confidentiality, even when patients themselves might consent to broader use. The Undiagnosed Diseases Network (UDN) and similar initiatives have navigated these challenges through comprehensive consent frameworks, but such approaches require substantial administrative infrastructure and may not scale to global federation (Ramoni et al., 2017).

The challenge, then, is to enable the benefits of federated patient matching while providing rigorous, quantifiable privacy protections. This thesis addresses that challenge.

## 1.4 Privacy-Preserving Computation for Genomic Data

The field of privacy-preserving computation offers a rich toolkit for protecting sensitive data while enabling useful analysis. Three complementary approaches are particularly relevant to phenotype matching:

**Private Set Intersection (PSI)** enables two parties to compute the intersection of their sets—such as shared phenotypes between two patients—without revealing elements outside the intersection (Meadows, 1986). Modern PSI protocols based on elliptic curve cryptography achieve practical efficiency while providing strong security guarantees. Recent work has extended PSI to support cardinality computation (revealing only the size of the intersection) and threshold variants (revealing intersection only if it exceeds a minimum size) (Pinkas et al., 2018).

**Differential privacy (DP)** provides a rigorous mathematical framework for quantifying privacy loss (Dwork et al., 2006). A differentially private mechanism guarantees that its outputs are statistically indistinguishable whether or not any single individual's data is included, with the privacy parameter ε controlling the strength of this guarantee. Differential privacy has been deployed at scale by organizations including the U.S. Census Bureau, Google, and Apple, demonstrating its practical viability (Abowd, 2018; Erlingsson et al., 2014).

**K-anonymity** requires that each released record be indistinguishable from at least k-1 other records with respect to quasi-identifying attributes (Sweeney, 2002). While weaker than differential privacy in theoretical guarantees, k-anonymity provides intuitive protection and is well-suited to medical data where rare attribute combinations pose the greatest risk. In the phenotype domain, k-anonymity can be achieved through rare term filtering or generalization to parent ontology terms.

Research at the intersection of privacy and genomics has explored applications of these techniques to specific problems. Montgomery and colleagues have developed approaches for identifying rare variants with large effects on gene expression while protecting individual privacy (GTEx Consortium et al., 2017). However, comprehensive frameworks for privacy-preserving phenotype matching—integrating multiple mechanisms and quantifying their combined utility and privacy implications—remain limited.

## 1.5 Thesis Contributions

This thesis presents a privacy-preserving phenotype matching framework that enables federated rare disease cohort discovery while protecting patient confidentiality. Our contributions are:

1. **A formal threat model and modular privacy framework.** We specify a two-party semi-honest threat model with auxiliary information, three concrete adversary goals (membership inference, attribute inference, singling-out re-identification), and per-protocol-step disclosure analysis (§3.1.2). Against that model we design and implement a composable pipeline integrating Diffie-Hellman PSI, Laplace/Gaussian/Exponential differential privacy mechanisms, and k-anonymity with rare-term filtering, and prove three composition-level privacy invariants.

2. **Empirically validated retrieval on real published patients.** Beyond synthetic-cohort experiments standard in the literature, we evaluate the system on 1,500 Phenopacket Store (Danis et al., 2025) patients drawn from published case reports with confirmed OMIM diagnoses, achieving MRR = 0.87 and nDCG@10 = 0.69 — within the Phenomizer/LIRICAL band. A full Resnik+BMA Phenomizer-style baseline runs on the same corpus for direct comparison.

3. **Empirical privacy measurement.** We implement Yeom-threshold and Shokri-shadow membership-inference attacks against the DP score-release oracle and a singling-out attack against the k-anonymity gate, producing the privacy-utility Pareto curves used to validate the formal invariants. Shadow-model attack ROC AUC drops from 0.98 (no DP) to 0.50 (random) at ε ≤ 1; k-anonymity at k = 10 reduces re-identification probability against the rare-term adversary from 0.42 to 0.005.

4. **The synthetic-to-real DP gap and a measured fix.** We document and explain a 20–50× discrepancy between safe ε on synthetic and real cohorts, propose the iterative exponential mechanism on rank utility as the principled response, and validate it: rank-utility exponential mechanism recovers 90% of non-private nDCG@10 at ε = 5 on the real cohort, versus 13% for Laplace at matched ε-DP — a 10× budget efficiency improvement. This finding has direct implications for published privacy-utility analyses that rely solely on disease-profile-sampled cohorts.

5. **GA4GH standards compatibility.** The implementation natively supports GA4GH Phenopackets v2.0, Beacon v2 query semantics, and Matchmaker Exchange message formats, enabling integration with existing federated networks.

6. **Open-source release.** Approximately 4,500 lines of Python under an open-source license, with reproducible evaluation pipelines for both synthetic and Phenopacket Store cohorts.

## 1.6 Thesis Organization

The remainder of this thesis is organized as follows:

**Chapter 2: Literature Review** surveys related work across rare disease phenotyping, phenotype similarity metrics, patient matching systems, and privacy-preserving computation. We position our contributions within this landscape and identify gaps addressed by our work.

**Chapter 3: Methods** presents our technical approach in detail, including data representation (Phenopackets, HPO), similarity metrics (Jaccard, Cosine, Resnik), privacy mechanisms (PSI, DP, k-anonymity), and evaluation methodology.

**Chapter 4: Results** reports experimental findings on baseline retrieval performance, privacy-utility tradeoffs, and privacy attack simulations using our evaluation dataset derived from OMIM and Orphanet disease annotations.

**Chapter 5: Discussion** interprets our results, discusses practical deployment considerations, addresses limitations, and considers ethical implications of privacy-preserving patient matching.

**Chapter 6: Future Work** outlines extensions including malicious-secure protocols, homomorphic encryption, federated learning, and validation on real patient cohorts.

**Chapter 7: Conclusion** summarizes contributions and broader impact.

---

## References

Abowd, J. M. (2018). The U.S. Census Bureau Adopts Differential Privacy. In *Proceedings of the 24th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining* (pp. 2867-2867).

Boycott, K. M., Rath, A., Chong, J. X., Hartley, T., Alkuraya, F. S., Baynam, G., ... & Lochmüller, H. (2017). International Cooperation to Enable the Diagnosis of All Rare Genetic Diseases. *American Journal of Human Genetics*, 100(5), 695-705.

Buske, O. J., Girdea, M., Engchuan, W., Brudno, M., & Consortium, C. R. D. (2015). PhenomeCentral: A Portal for Phenotypic and Genotypic Matchmaking of Patients with Rare Genetic Diseases. *Human Mutation*, 36(10), 931-940.

Dwork, C., McSherry, F., Nissim, K., & Smith, A. (2006). Calibrating Noise to Sensitivity in Private Data Analysis. In *Theory of Cryptography Conference* (pp. 265-284). Springer.

El Emam, K., Jonker, E., Arbuckle, L., & Malin, B. (2011). A Systematic Review of Re-Identification Attacks on Health Data. *PLoS ONE*, 6(12), e28071.

Erlingsson, Ú., Pihur, V., & Korolova, A. (2014). RAPPOR: Randomized Aggregatable Privacy-Preserving Ordinal Response. In *Proceedings of the 2014 ACM SIGSAC Conference on Computer and Communications Security* (pp. 1054-1067).

EURORDIS. (2009). The Voice of 12,000 Patients: Experiences and Expectations of Rare Disease Patients on Diagnosis and Care in Europe. EURORDIS-Rare Diseases Europe.

Ferreira, C. R. (2019). The Burden of Rare Diseases. *American Journal of Medical Genetics Part A*, 179(6), 885-892.

Frésard, L., Smail, C., Ferraro, N. M., Teran, N. A., Li, X., Smith, K. S., ... & Montgomery, S. B. (2019). Identification of Rare-Disease Genes Using Blood Transcriptome Sequencing and Large Control Cohorts. *Nature Medicine*, 25(6), 911-919.

Global Genes. (2023). RARE Disease Facts. Retrieved from https://globalgenes.org/rare-disease-facts/

GREGoR Consortium. (2025). GREGoR: Accelerating Genomics for Rare Diseases. *Nature*, in press.

GTEx Consortium. (2017). Genetic Effects on Gene Expression Across Human Tissues. *Nature*, 550(7675), 204-213.

GTEx Consortium, Aguet, F., Ardlie, K. G., Cummings, B. B., Gelfand, E. T., Getz, G., ... & Montgomery, S. B. (2017). The Impact of Rare Variation on Gene Expression Across Tissues. *Nature*, 550(7675), 239-243.

Jacobsen, J. O. B., Baudis, M., Baynam, G. S., Beckmann, J. S., Beltran, S., Buske, O. J., ... & Robinson, P. N. (2022). The GA4GH Phenopacket Schema Defines a Computable Representation of Clinical Data. *Nature Biotechnology*, 40(6), 817-820.

Köhler, S., Gargano, M., Matentzoglu, N., Carmody, L. C., Lewis-Smith, D., Vasilevsky, N. A., ... & Robinson, P. N. (2021). The Human Phenotype Ontology in 2021. *Nucleic Acids Research*, 49(D1), D1207-D1217.

Meadows, C. (1986). A More Efficient Cryptographic Matchmaking Protocol for Use in the Absence of a Continuously Available Third Party. In *IEEE Symposium on Security and Privacy* (pp. 134-137).

Nguengang Wakap, S., Lambert, D. M., Olry, A., Rodwell, C., Gueydan, C., Lanneau, V., ... & Rath, A. (2020). Estimating Cumulative Point Prevalence of Rare Diseases: Analysis of the Orphanet Database. *European Journal of Human Genetics*, 28(2), 165-173.

NORD. (2023). Rare Disease Facts. National Organization for Rare Disorders.

Philippakis, A. A., Azzariti, D. R., Birney, E., Brookes, A. J., Buske, O. J., Dollé, L., ... & Rehm, H. L. (2015). The Matchmaker Exchange: A Platform for Rare Disease Gene Discovery. *Human Mutation*, 36(10), 915-921.

Pinkas, B., Schneider, T., & Zohner, M. (2018). Scalable Private Set Intersection Based on OT Extension. *ACM Transactions on Privacy and Security*, 21(2), 1-35.

Raisaro, J. L., Tramèr, F., Ji, Z., Bu, D., Zhao, Y., Carey, K., ... & Hubaux, J. P. (2017). Addressing Beacon Re-Identification Attacks: Quantification and Mitigation of Privacy Risks. *Journal of the American Medical Informatics Association*, 24(4), 799-805.

Ramoni, R. B., Mulvihill, J. J., Adams, D. R., Allard, P., Ashley, E. A., Bernstein, J. A., ... & Tifft, C. J. (2017). The Undiagnosed Diseases Network: Accelerating Discovery about Health and Disease. *American Journal of Human Genetics*, 100(2), 185-192.

Shringarpure, S. S., & Bustamante, C. D. (2015). Privacy Risks from Genomic Data-Sharing Beacons. *American Journal of Human Genetics*, 97(5), 631-646.

Sobreira, N., Schiettecatte, F., Valle, D., & Hamosh, A. (2015). GeneMatcher: A Matching Tool for Connecting Investigators with an Interest in the Same Gene. *Human Mutation*, 36(10), 928-930.

Spillmann, R. C., McConkie-Rosell, A., Pena, L. D., Jiang, Y. H., Schoch, K., Walley, N., ... & Shashi, V. (2017). A Window into Living with an Undiagnosed Disease: Illness Narratives from the Undiagnosed Diseases Network. *Orphanet Journal of Rare Diseases*, 12(1), 1-10.

Sweeney, L. (2002). k-Anonymity: A Model for Protecting Privacy. *International Journal of Uncertainty, Fuzziness and Knowledge-Based Systems*, 10(5), 557-570.
