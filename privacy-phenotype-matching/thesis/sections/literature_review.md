# Literature Review

This chapter surveys the research landscape underlying privacy-preserving phenotype matching, spanning five interconnected domains: rare disease phenotyping standards, semantic similarity measures, federated patient matching systems, privacy-preserving computation techniques, and privacy attack methodologies. We synthesize findings from over 60 publications to contextualize our contributions and identify gaps addressed by this thesis.

## 2.1 Rare Disease Phenotyping and Standards

### 2.1.1 The Human Phenotype Ontology

The Human Phenotype Ontology (HPO) has emerged as the de facto standard for computational representation of clinical phenotypes in rare disease research. Originally developed by Robinson et al. (2008) at the Charité-Universitätsmedizin Berlin, HPO provides a structured vocabulary of phenotypic abnormalities organized as a directed acyclic graph (DAG).

The ontology has grown substantially since its inception. The 2021 release contains over 16,000 terms describing phenotypic abnormalities, with an additional 2,000+ terms for inheritance patterns, clinical modifiers, and frequency annotations (Köhler et al., 2021). Terms are organized hierarchically, with "Phenotypic abnormality" (HP:0000118) as the root, branching into major categories including "Abnormality of the nervous system," "Abnormality of the cardiovascular system," and 20 other organ-system domains.

HPO's design follows established ontology engineering principles. Each term has a unique identifier, label, definition, and synonyms. The subsumption hierarchy (IS-A relationships) enables both specific annotation and semantic inference—a patient annotated with "Focal clonic seizure" (HP:0002266) implicitly also has the more general "Seizure" (HP:0001250). This property is fundamental to semantic similarity computation.

The HPO consortium maintains active curation, with approximately 2,000 new terms added annually (Köhler et al., 2019). Integration with other ontologies—including the Gene Ontology (GO), Disease Ontology (DO), and Unified Medical Language System (UMLS)—enhances interoperability. Cross-references to OMIM, Orphanet, and DECIPHER link phenotypes to disease entities, enabling phenotype-driven differential diagnosis.

Several computational tools leverage HPO for clinical applications. Phenomizer uses semantic similarity to rank candidate diagnoses given a set of patient phenotypes (Köhler et al., 2009). Exomiser integrates HPO-based phenotype matching with variant prioritization for exome sequencing interpretation (Smedley et al., 2015). LIRICAL employs likelihood ratios to quantify phenotypic evidence supporting each candidate diagnosis (Robinson et al., 2020). These tools demonstrate the practical utility of standardized phenotype representation.

### 2.1.2 GA4GH Phenopackets

The Global Alliance for Genomics and Health (GA4GH) Phenopacket schema provides a comprehensive standard for representing clinical and genomic data (Jacobsen et al., 2022). Building on earlier work by the Monarch Initiative (Mungall et al., 2017), Phenopackets version 2.0 defines a hierarchical schema encompassing:

- **Individual**: demographic information (age, sex, vital status)
- **PhenotypicFeature**: observed or excluded HPO terms with modifiers
- **Disease**: diagnosed conditions with ontology identifiers
- **Biosample**: specimen information for molecular analysis
- **Interpretation**: clinical conclusions from genomic analysis
- **GenomicInterpretation**: variant-level findings

Phenopackets employ a "building block" architecture where complex clinical cases are constructed from reusable components. This modularity facilitates both human readability and machine processing. The schema is specified in Protocol Buffers with JSON serialization, enabling efficient transmission and storage.

Adoption has accelerated since the 2.0 release. The 100,000 Genomes Project uses Phenopackets for phenotype submission (Turnbull et al., 2018). ClinVar accepts Phenopacket-formatted submissions for variant interpretation. The Matchmaker Exchange network requires Phenopacket-compatible data formats. This growing ecosystem validates Phenopackets as the interoperability layer for rare disease genomics.

### 2.1.3 Disease Databases

Several authoritative databases provide curated disease-phenotype associations essential for phenotype matching.

**OMIM (Online Mendelian Inheritance in Man)** catalogs human genes and genetic disorders with detailed clinical descriptions (Amberger et al., 2019). Founded by Victor McKusick in 1966 as a print catalog, OMIM now contains entries for over 17,000 genes and 9,000 phenotypes with known molecular basis. Each entry includes clinical features, inheritance patterns, molecular genetics, and literature references. OMIM's gene-phenotype associations, curated from primary literature by expert geneticists, provide gold-standard annotations for computational analysis.

**Orphanet** serves as the European reference database for rare diseases (Rath et al., 2012). Established in 1997 by the French National Institute of Health and Medical Research (INSERM), Orphanet now operates as a consortium of 40 countries. The database contains information on over 6,000 rare diseases, including clinical descriptions, epidemiology, and management guidelines. Orphanet's systematic HPO annotations, derived from literature review and expert consensus, complement OMIM with broader European disease coverage.

**DECIPHER (Database of Chromosomal Imbalance and Phenotype in Humans using Ensembl Resources)** focuses on copy number variants and their phenotypic consequences (Firth et al., 2009). Developed at the Wellcome Sanger Institute, DECIPHER contains data from over 44,000 patients with chromosomal imbalances. Each record includes genomic coordinates, phenotype annotations, and—for consented cases—anonymized clinical details. DECIPHER's population-scale data enables genotype-phenotype correlation studies and CNV interpretation.

**ClinVar** aggregates variant interpretations from clinical laboratories worldwide (Landrum et al., 2018). Submissions include associated phenotypes, enabling phenotype-aware variant prioritization. With over 2 million variant records, ClinVar provides the largest public repository of clinically relevant genetic variation.

The HPO phenotype annotation file (phenotype.hpoa) integrates disease-phenotype associations from these sources. The current release contains 263,000+ annotations linking 13,000 diseases to 11,500 HPO terms. This comprehensive resource enables corpus-based information content computation and supports synthetic patient generation for algorithm evaluation.

## 2.2 Semantic Similarity Measures

Quantifying phenotypic similarity between patients requires principled approaches to comparing sets of ontology terms. Research in this area draws from information retrieval, computational linguistics, and biomedical informatics.

### 2.2.1 Set-Theoretic Measures

The simplest similarity measures treat phenotype profiles as sets and compute overlap statistics.

**Jaccard Index.** Introduced by Paul Jaccard in 1901 for botanical taxonomy, the Jaccard index measures the ratio of intersection to union cardinality:

$$J(A, B) = \frac{|A \cap B|}{|A \cup B|}$$

Jaccard similarity ranges from 0 (disjoint sets) to 1 (identical sets). Its simplicity and interpretability make it widely used, but it treats all terms as equally informative—a significant limitation when comparing rare versus common phenotypes.

**Dice Coefficient.** The Dice coefficient (Dice, 1945) provides an alternative formulation emphasizing intersection size:

$$D(A, B) = \frac{2|A \cap B|}{|A| + |B|}$$

Dice and Jaccard are monotonically related: $D = 2J/(1+J)$. Both require exact term matches, ignoring semantic relationships between related terms.

**Overlap Coefficient.** The overlap coefficient normalizes by the smaller set:

$$O(A, B) = \frac{|A \cap B|}{\min(|A|, |B|)}$$

This measure is useful when comparing sets of different sizes, as it achieves 1.0 when one set is a subset of the other.

### 2.2.2 Information-Theoretic Measures

Information content (IC) quantifies term specificity based on corpus frequency, enabling differential weighting of rare versus common phenotypes.

**Resnik Similarity.** Resnik (1995) proposed measuring similarity based on the information content of the most informative common ancestor (MICA):

$$sim_{Resnik}(t_1, t_2) = IC(MICA(t_1, t_2))$$

where $IC(t) = -\log P(t)$ and $P(t)$ is the probability of encountering term $t$ in a reference corpus. Terms sharing a specific common ancestor (high IC) are more similar than those sharing only a general ancestor (low IC).

Resnik similarity has several desirable properties: it is symmetric, reflexive, and bounded below by 0. However, it is not normalized—maximum similarity depends on the terms being compared—complicating cross-comparison interpretation.

**Lin Similarity.** Lin (1998) proposed a normalized variant based on information-theoretic principles:

$$sim_{Lin}(t_1, t_2) = \frac{2 \cdot IC(MICA(t_1, t_2))}{IC(t_1) + IC(t_2)}$$

Lin similarity ranges from 0 to 1, with 1 indicating identical terms. The normalization accounts for term specificity, making similarity scores comparable across term pairs.

**Jiang-Conrath Distance.** Jiang and Conrath (1997) proposed a distance measure (lower = more similar):

$$dist_{JC}(t_1, t_2) = IC(t_1) + IC(t_2) - 2 \cdot IC(MICA(t_1, t_2))$$

This can be converted to similarity via $sim_{JC} = 1 / (1 + dist_{JC})$. Jiang-Conrath often performs well in empirical evaluations (Pesquita et al., 2009).

### 2.2.3 Aggregation Strategies

Extending pairwise term similarity to sets of phenotypes requires aggregation. Several strategies have been proposed.

**Best-Match Average (BMA).** For each term in set A, find its maximum similarity to any term in B, then average across both directions:

$$BMA(A, B) = \frac{1}{2}\left(\frac{1}{|A|}\sum_{a \in A} \max_{b \in B} sim(a,b) + \frac{1}{|B|}\sum_{b \in B} \max_{a \in A} sim(a,b)\right)$$

BMA is symmetric and handles sets of different sizes. It performs well empirically and is used by tools including Phenomizer (Köhler et al., 2009).

**Maximum.** Take the maximum pairwise similarity across all term pairs. This strategy emphasizes the strongest match but ignores overall profile similarity.

**Average.** Average all pairwise similarities. This can be dominated by low-similarity pairs when sets are large.

**Graph-Based Information Content.** Schlicker et al. (2006) proposed GraSM, which computes IC based on graph structure rather than corpus frequency. This approach is useful when corpus statistics are unavailable but requires careful ontology design.

### 2.2.4 Intrinsic Information Content

When corpus-based IC is unavailable, intrinsic approaches estimate term specificity from ontology structure alone.

Seco et al. (2004) proposed IC based on hyponym (descendant) counts:

$$IC_{Seco}(t) = 1 - \frac{\log(|hypo(t)| + 1)}{\log(|C|)}$$

where $|hypo(t)|$ is the number of hyponyms and $|C|$ is the total concept count. Leaf nodes (no hyponyms) receive maximum IC.

Sánchez et al. (2011) extended this with additional structural features. Zhou et al. (2008) incorporated depth information. These measures enable similarity computation without requiring annotated corpora, though they may not reflect actual term usage patterns in clinical practice.

### 2.2.5 Evaluation of Similarity Measures

Pesquita et al. (2009) provide a comprehensive survey of semantic similarity measures for biomedical ontologies, comparing 13 measures across multiple evaluation tasks. Key findings include:

- IC-based measures consistently outperform edge-counting approaches
- BMA aggregation performs well across diverse datasets
- Corpus-based IC generally outperforms intrinsic IC when annotations are available
- Resnik and Lin similarity show strong performance for gene function prediction

For HPO specifically, Köhler et al. (2009) demonstrated that semantic similarity effectively identifies correct diagnoses among thousands of candidates, with Resnik-based measures achieving top performance.

## 2.3 Federated Patient Matching Systems

### 2.3.1 Matchmaker Exchange

The Matchmaker Exchange (MME) represents the most successful implementation of federated rare disease patient matching (Philippakis et al., 2015). Established under GA4GH auspices, MME connects seven databases:

- **GeneMatcher** (Johns Hopkins University): gene-centric matching (Sobreira et al., 2015)
- **PhenomeCentral** (University of Toronto): phenotype-first matching (Buske et al., 2015)
- **DECIPHER** (Wellcome Sanger Institute): CNV-phenotype matching
- **MyGene2** (University of Washington): family-driven matching
- **PatientMatcher** (Sweden): national rare disease network
- **seqr** (Broad Institute): exome analysis platform
- **Matchbox** (NIH UDP): Undiagnosed Diseases Program

MME operates via a standardized API where each node can query others for matching patients. The query format includes phenotypes (HPO terms), genomic features (genes, variants), and patient metadata. Matching algorithms are implemented locally by each node, enabling algorithmic diversity while maintaining interoperability.

Since 2015, MME has facilitated thousands of matches contributing to over 100 disease gene discoveries. Notably, MME enabled identification of TBCK-related intellectual disability syndrome through matching patients across three continents (Bhoj et al., 2016).

### 2.3.2 GeneMatcher

GeneMatcher pioneered the concept of gene-based patient matching (Sobreira et al., 2015). Researchers submit genes of interest, and the system identifies other submitters interested in the same genes. This "dating service for genes" has connected researchers on over 6,000 genes, directly contributing to novel disease gene publications.

The gene-centric approach complements phenotype matching: genes identify candidate collaborators, while phenotype comparison validates clinical similarity. GeneMatcher's success demonstrates demand for federated matching, though its gene-first model requires prior variant identification.

### 2.3.3 PhenomeCentral

PhenomeCentral implements phenotype-driven matching using semantic similarity (Buske et al., 2015). Patients are represented by HPO term profiles, and similarity is computed using the Resnik-based HRSS algorithm (Hybrid Relative Specificity Similarity). The system prioritizes matches where phenotype overlap is concentrated in specific, diagnostically informative features.

PhenomeCentral introduced several innovations: tiered access control (local, matchmaker, public), configurable similarity thresholds, and integration with variant data. Its architecture influenced the broader MME design.

### 2.3.4 Beacon Network

The GA4GH Beacon Network enables queries about variant presence across federated databases (Fiume et al., 2019). In its simplest form, a Beacon responds to the question "Do you have any genomes with variant X at position Y?" with a Boolean yes/no answer.

Beacon v2, released in 2022, extends the model substantially. Queries can include phenotype filters (HPO terms), and responses can provide counts or full records with appropriate authorization. The granularity model—boolean, count, record—enables institutions to participate at their comfort level.

Beacons demonstrated the feasibility of federated genomic queries but also exposed privacy vulnerabilities, as discussed in Section 2.5.

### 2.3.5 Limitations of Current Systems

Despite their success, current matching systems face limitations:

**Privacy concerns.** Participation requires sharing patient data with other institutions, creating regulatory and ethical barriers. Many institutions decline to participate due to privacy policies.

**Data governance.** Matching outputs reveal sensitive information about patient phenotypes. Even negative results can be informative (the patient doesn't match disease X).

**Scalability.** Full pairwise similarity computation scales quadratically with database size. Current systems mitigate this through indexing and caching but face challenges as databases grow.

**Semantic heterogeneity.** Despite HPO standardization, phenotype annotation practices vary across institutions, affecting match quality.

Our work addresses the first two limitations through privacy-preserving computation techniques.

## 2.4 Privacy-Preserving Computation

### 2.4.1 Private Set Intersection

Private Set Intersection (PSI) enables two parties to compute the intersection of their sets without revealing elements outside the intersection. PSI directly applies to phenotype matching: determining shared phenotypes without exposing non-shared terms.

**Diffie-Hellman PSI.** The earliest practical PSI protocols use Diffie-Hellman key exchange (Meadows, 1986; Huberman et al., 1999). Given sets $X$ and $Y$, the protocol proceeds:

1. Party A hashes and exponentiates each element: $H(x)^a$ for $x \in X$
2. Party B exponentiates received values and their own elements: $(H(x)^a)^b$ and $H(y)^b$
3. Party A exponentiates B's values: $(H(y)^b)^a$
4. Intersection: elements where $H(x)^{ab} = H(y)^{ab}$

Under the Decisional Diffie-Hellman assumption, this protocol is secure against semi-honest adversaries. Communication complexity is O(|X| + |Y|) group elements.

**OT-Based PSI.** Oblivious Transfer (OT) extensions enable more efficient PSI for large sets (Pinkas et al., 2014, 2015). The core idea uses OT to evaluate a pseudorandom function on set elements, then compare outputs. The PSZ protocol (Pinkas et al., 2018) achieves O(n) computation and O(n) communication for sets of size n, representing a significant improvement over naive approaches.

**Circuit-Based PSI.** Secure computation of Boolean or arithmetic circuits can implement PSI with stronger security guarantees (Huang et al., 2012). Garbled circuits achieve security against malicious adversaries but incur higher computational overhead. The Pinkas et al. (2018) circuit-based protocol using cuckoo hashing achieves practical performance for million-element sets.

**PSI Variants.** Several PSI variants address specific requirements:
- **Cardinality PSI**: reveals only $|X \cap Y|$, not intersection elements (De Cristofaro & Tsudik, 2010)
- **Threshold PSI**: reveals intersection only if $|X \cap Y| \geq t$ (Ghosh & Simkin, 2019)
- **Fuzzy PSI**: tolerates approximate matches (Indyk & Woodruff, 2006)
- **Multi-party PSI**: extends to n > 2 parties (Kolesnikov et al., 2017)

For phenotype matching, cardinality PSI enables Jaccard similarity computation ($|X \cap Y|$ and $|X \cup Y| = |X| + |Y| - |X \cap Y|$) without revealing which terms match.

### 2.4.2 Differential Privacy

Differential privacy provides a mathematical framework for quantifying privacy loss (Dwork et al., 2006; Dwork, 2006). A randomized mechanism $\mathcal{M}$ satisfies $(\varepsilon, \delta)$-differential privacy if for all datasets $D_1, D_2$ differing in one record and all output sets $S$:

$$\Pr[\mathcal{M}(D_1) \in S] \leq e^\varepsilon \cdot \Pr[\mathcal{M}(D_2) \in S] + \delta$$

The privacy parameter $\varepsilon$ (epsilon) bounds the multiplicative difference in output probabilities caused by any single record. Smaller $\varepsilon$ provides stronger privacy. The parameter $\delta$ allows for a small probability of greater leakage; pure DP sets $\delta = 0$.

**Laplace Mechanism.** For a numeric query $f$ with sensitivity $\Delta f = \max_{D_1, D_2} |f(D_1) - f(D_2)|$, adding Laplace noise achieves pure $\varepsilon$-DP:

$$\mathcal{M}(D) = f(D) + \text{Lap}(\Delta f / \varepsilon)$$

Similarity scores bounded in [0,1] have sensitivity 1, enabling straightforward application.

**Gaussian Mechanism.** For $(\varepsilon, \delta)$-DP with $\delta > 0$, Gaussian noise often provides better utility:

$$\mathcal{M}(D) = f(D) + \mathcal{N}(0, \sigma^2), \quad \sigma = \frac{\Delta f \sqrt{2\ln(1.25/\delta)}}{\varepsilon}$$

The Gaussian mechanism is particularly useful for high-dimensional outputs where Laplace noise compounds poorly.

**Exponential Mechanism.** For selecting among discrete options with utility function $u(D, r)$, the exponential mechanism samples proportionally to exponentiated utility (McSherry & Talwar, 2007):

$$\Pr[\mathcal{M}(D) = r] \propto \exp\left(\frac{\varepsilon \cdot u(D, r)}{2\Delta u}\right)$$

This is useful for private top-k selection in phenotype matching results.

**Composition.** When multiple DP mechanisms are applied, privacy degrades. Basic composition gives $\varepsilon_{total} = \sum_i \varepsilon_i$. Advanced composition (Dwork et al., 2010) provides tighter bounds:

$$\varepsilon_{total} = \sqrt{2k \ln(1/\delta')} \cdot \varepsilon + k\varepsilon(e^\varepsilon - 1)$$

for $k$ mechanisms each satisfying $(\varepsilon, \delta)$-DP. Privacy accounting (Abadi et al., 2016) enables precise tracking across complex pipelines.

**DP in Healthcare.** Differential privacy has seen increasing healthcare adoption. The US Census Bureau uses DP for disclosure avoidance (Abowd, 2018). Apple and Google deploy DP for telemetry collection (Erlingsson et al., 2014). Medical applications include DP-protected GWAS (Johnson & Shmatikov, 2013), clinical trial analysis (Dankar & El Emam, 2013), and genomic Beacons (Raisaro et al., 2017).

### 2.4.3 K-Anonymity and Generalization

K-anonymity requires that each released record be indistinguishable from at least $k-1$ others with respect to quasi-identifying attributes (Sweeney, 2002). Originally developed for de-identification of tabular data, k-anonymity concepts apply to phenotype protection.

**Quasi-Identifiers.** Phenotype combinations can uniquely identify individuals. A patient with "Tetralogy of Fallot," "Polydactyly," and "Microcephaly" may be the only such person in a database. These rare combinations serve as quasi-identifiers.

**Generalization.** K-anonymity is typically achieved through generalization (replacing specific values with more general ones) and suppression (removing values). In the phenotype domain, the HPO hierarchy enables natural generalization: "Tetralogy of Fallot" can be generalized to "Conotruncal defect" or further to "Abnormal heart morphology."

**Limitations.** K-anonymity has known weaknesses. L-diversity (Machanavajjhala et al., 2007) addresses homogeneity attacks where all k records share a sensitive value. T-closeness (Li et al., 2007) requires that sensitive attribute distributions in equivalence classes match the global distribution. These extensions can guide phenotype generalization strategies.

### 2.4.4 Secure Multi-Party Computation

Secure Multi-Party Computation (MPC) enables parties to jointly compute functions on private inputs without revealing those inputs (Yao, 1986; Goldreich et al., 1987). MPC provides a general framework encompassing PSI as a special case.

**Garbled Circuits.** Yao's garbled circuits protocol encrypts a Boolean circuit such that one party can evaluate it on encrypted inputs without learning intermediate values. Modern optimizations (half-gates, free-XOR) achieve practical performance for complex functions (Zahur et al., 2015).

**Secret Sharing.** Protocols based on secret sharing (Shamir, 1979) distribute computation across multiple parties, with privacy guaranteed as long as fewer than a threshold collude. The SPDZ protocol (Damgård et al., 2012) provides malicious security with preprocessing.

**Homomorphic Encryption.** Fully homomorphic encryption (FHE) enables computation on encrypted data (Gentry, 2009). While theoretically powerful, FHE incurs substantial overhead. Partial HE schemes (Paillier, ElGamal) enable specific operations more efficiently. Recent advances in BGV and CKKS schemes have improved practicality (Cheon et al., 2017).

For phenotype matching, MPC could enable private similarity computation without PSI's intersection revelation. However, the computational overhead currently limits practical deployment.

## 2.5 Privacy Attacks and Auditing

Understanding attack vectors is essential for designing robust privacy protections.

### 2.5.1 Genomic Beacon Attacks

Shringarpure and Bustamante (2015) demonstrated re-identification attacks against the Beacon Network. Despite revealing only Boolean responses ("yes, this variant exists"), repeated queries enable statistical inference. With approximately 5,000 queries, an adversary can determine whether a target individual's genome is in the database.

The attack exploits the statistical signature of an individual's genome. Querying variants known to be present in the target yields more "yes" responses than expected for a random individual. This differential enables high-confidence membership inference.

Raisaro et al. (2017) proposed DP-protected Beacons that add noise to responses. However, utility degrades substantially under strong privacy parameters. The Beacon community subsequently developed access control mechanisms limiting query rates.

### 2.5.2 Membership Inference

Membership inference attacks determine whether a specific record was used in training a model or generating a statistic (Shokri et al., 2017). The adversary trains shadow models on datasets with known membership, learning to distinguish members from non-members based on output patterns.

For phenotype matching, membership inference could reveal whether a patient with specific phenotypes exists in a database. The attack is particularly concerning for rare phenotype combinations that create distinctive query responses.

Defenses include DP mechanisms, output perturbation, and limiting query access. Our evaluation includes membership inference experiments to empirically measure privacy protection.

### 2.5.3 Attribute Inference

Attribute inference attacks predict sensitive attributes not directly revealed in outputs (Yeom et al., 2018). In phenotype matching, an adversary might infer undisclosed phenotypes from similarity scores with known patient profiles.

The attack proceeds by: (1) training a model relating observable outputs to sensitive attributes on auxiliary data, (2) applying this model to infer attributes of target individuals. Success depends on correlation between revealed and hidden information.

Our k-anonymity and rare term filtering mechanisms specifically target attribute inference by removing or generalizing identifying phenotypes before any computation.

### 2.5.4 Reconstruction Attacks

Reconstruction attacks recover individual records from aggregate statistics (Dinur & Nissim, 2003). This foundational result showed that answering too many linear queries enables full database reconstruction, motivating the development of differential privacy.

For phenotype databases, reconstruction could recover individual phenotype profiles from sufficient similarity queries. DP provides provable protection against reconstruction within privacy budget bounds.

### 2.5.5 Privacy Auditing

Jagielski et al. (2020) introduced privacy auditing frameworks that empirically estimate privacy leakage. Rather than relying solely on theoretical guarantees, auditing measures actual information leakage through controlled experiments.

Key metrics include:
- **Empirical privacy**: the observed ε in membership inference experiments
- **Adversarial advantage**: accuracy improvement over random guessing
- **Distinguishing probability**: confidence in distinguishing members from non-members

We adopt this empirical approach in our evaluation, complementing theoretical DP guarantees with measured attack success rates.

## 2.6 Related Privacy-Preserving Systems

### 2.6.1 Secure Genome Analysis

Ayday et al. (2013) surveyed privacy-preserving techniques for genomic data, categorizing approaches by threat model and functionality. Key techniques include:

- **Encrypted storage**: patient data encrypted at rest, decrypted only for authorized computation
- **Access control**: fine-grained permissions on genomic data elements
- **Secure computation**: MPC/HE for privacy-preserving analysis

The iDASH (integrating Data for Analysis, Anonymization, and Sharing) competition has driven development of practical secure genomic tools, with challenges including secure GWAS, private sequence comparison, and encrypted machine learning (Tang et al., 2020).

### 2.6.2 Privacy-Preserving GWAS

Genome-wide association studies (GWAS) aggregate genetic data across thousands of individuals, creating both scientific value and privacy risk. Johnson and Shmatikov (2013) proposed DP-GWAS using the exponential mechanism for SNP selection.

Subsequent work addressed specific GWAS computations: chi-squared tests (Uhlerop et al., 2013), regression analysis (Chen et al., 2019), and meta-analysis across institutions (Cho et al., 2018). These methods inform our approach to privacy-utility tradeoffs in aggregate phenotype analysis.

### 2.6.3 The Montgomery Lab and Rare Variant Expression

Relevant to our work, Montgomery and colleagues have developed approaches for identifying rare variants affecting gene expression while preserving privacy (GTEx Consortium et al., 2017). Their analysis of the GTEx dataset demonstrated that rare variants contribute substantially to expression variation, with outlier analysis providing a powerful discovery tool.

Frésard et al. (2019) extended this approach to rare disease diagnosis, using RNA-seq outlier detection to identify causal variants missed by DNA sequencing alone. This work demonstrates the value of multi-modal data integration for rare disease, while also highlighting the sensitivity of expression data that motivates privacy protection.

The GREGoR Consortium, which Montgomery co-leads, is generating multi-omic data on thousands of rare disease families with plans for broad data sharing (GREGoR Consortium, 2025). Our privacy-preserving matching framework could enable queries against this resource while protecting participant privacy.

## 2.7 Gaps and Contributions

Our review identifies several gaps in the literature that this thesis addresses:

1. **Integrated privacy framework.** While PSI, DP, and k-anonymity have been studied independently, their composition for phenotype matching has not been systematically evaluated.

2. **HPO-aware privacy.** Existing PSI protocols do not leverage ontology structure for semantic protection or similarity computation.

3. **Practical evaluation on rare diseases.** Privacy-preserving genomic tools have been evaluated primarily on common variant data; rare disease phenotype matching presents distinct challenges.

4. **Privacy-utility quantification.** Empirical characterization of tradeoffs across mechanism configurations, including combined approaches, remains limited.

5. **GA4GH compatibility.** Integration of privacy mechanisms with Phenopackets, Beacon, and MME standards enables practical deployment within existing infrastructure.

This thesis addresses these gaps through a modular, standards-compliant privacy framework with comprehensive evaluation on realistic rare disease data.

---

## References

Abadi, M., Chu, A., Goodfellow, I., McMahan, H. B., Mironov, I., Talwar, K., & Zhang, L. (2016). Deep Learning with Differential Privacy. In *ACM CCS* (pp. 308-318).

Abowd, J. M. (2018). The U.S. Census Bureau Adopts Differential Privacy. In *ACM KDD* (pp. 2867-2867).

Amberger, J. S., Bocchini, C. A., Scott, A. F., & Hamosh, A. (2019). OMIM.org: Leveraging Knowledge Across Phenotype-Gene Relationships. *Nucleic Acids Research*, 47(D1), D1038-D1043.

Ayday, E., De Cristofaro, E., Hubaux, J. P., & Tsudik, G. (2013). Whole Genome Sequencing: Revolutionary Medicine or Privacy Nightmare? *IEEE Computer*, 46(2), 58-66.

Bhoj, E. J., Li, D., Harr, M., et al. (2016). Mutations in TBCK, Encoding TBC1-Domain-Containing Kinase, Lead to a Recognizable Syndrome of Intellectual Disability and Hypotonia. *American Journal of Human Genetics*, 98(4), 782-788.

Buske, O. J., Girdea, M., Dumitriu, S., et al. (2015). PhenomeCentral: A Portal for Phenotypic and Genotypic Matchmaking of Patients with Rare Genetic Diseases. *Human Mutation*, 36(10), 931-940.

Chen, F., Wang, S., Jiang, X., et al. (2019). PRINCESS: Privacy-Protecting Rare Disease International Network Collaboration via Encryption through Software Guard Extensions. *Bioinformatics*, 35(5), 871-878.

Cheon, J. H., Kim, A., Kim, M., & Song, Y. (2017). Homomorphic Encryption for Arithmetic of Approximate Numbers. In *ASIACRYPT* (pp. 409-437).

Cho, H., Wu, D. J., & Berger, B. (2018). Secure Genome-Wide Association Analysis Using Multiparty Computation. *Nature Biotechnology*, 36(6), 547-551.

Damgård, I., Pastro, V., Smart, N., & Zakarias, S. (2012). Multiparty Computation from Somewhat Homomorphic Encryption. In *CRYPTO* (pp. 643-662).

Dankar, F. K., & El Emam, K. (2013). Practicing Differential Privacy in Health Care: A Review. *Transactions on Data Privacy*, 6(1), 35-67.

De Cristofaro, E., & Tsudik, G. (2010). Practical Private Set Intersection Protocols with Linear Complexity. In *Financial Cryptography* (pp. 143-159).

Dice, L. R. (1945). Measures of the Amount of Ecologic Association Between Species. *Ecology*, 26(3), 297-302.

Dinur, I., & Nissim, K. (2003). Revealing Information While Preserving Privacy. In *ACM PODS* (pp. 202-210).

Dwork, C. (2006). Differential Privacy. In *ICALP* (pp. 1-12). Springer.

Dwork, C., McSherry, F., Nissim, K., & Smith, A. (2006). Calibrating Noise to Sensitivity in Private Data Analysis. In *TCC* (pp. 265-284). Springer.

Dwork, C., Rothblum, G. N., & Vadhan, S. (2010). Boosting and Differential Privacy. In *IEEE FOCS* (pp. 51-60).

Erlingsson, Ú., Pihur, V., & Korolova, A. (2014). RAPPOR: Randomized Aggregatable Privacy-Preserving Ordinal Response. In *ACM CCS* (pp. 1054-1067).

Firth, H. V., Richards, S. M., Bevan, A. P., et al. (2009). DECIPHER: Database of Chromosomal Imbalance and Phenotype in Humans Using Ensembl Resources. *American Journal of Human Genetics*, 84(4), 524-533.

Fiume, M., Cupak, M., Keenan, S., et al. (2019). Federated Discovery and Sharing of Genomic Data Using Beacons. *Nature Biotechnology*, 37(3), 220-224.

Frésard, L., Smail, C., Ferraro, N. M., et al. (2019). Identification of Rare-Disease Genes Using Blood Transcriptome Sequencing and Large Control Cohorts. *Nature Medicine*, 25(6), 911-919.

Gentry, C. (2009). A Fully Homomorphic Encryption Scheme. *PhD Thesis*, Stanford University.

Ghosh, S., & Simkin, M. (2019). The Communication Complexity of Threshold Private Set Intersection. In *CRYPTO* (pp. 3-29).

Goldreich, O., Micali, S., & Wigderson, A. (1987). How to Play Any Mental Game. In *ACM STOC* (pp. 218-229).

GREGoR Consortium. (2025). GREGoR: Accelerating Genomics for Rare Diseases. *Nature*, in press.

GTEx Consortium, et al. (2017). Genetic Effects on Gene Expression Across Human Tissues. *Nature*, 550(7675), 204-213.

GTEx Consortium, et al. (2017). The Impact of Rare Variation on Gene Expression Across Tissues. *Nature*, 550(7675), 239-243.

Huang, Y., Evans, D., & Katz, J. (2012). Private Set Intersection: Are Garbled Circuits Better than Custom Protocols? In *NDSS*.

Huberman, B. A., Franklin, M., & Hogg, T. (1999). Enhancing Privacy and Trust in Electronic Communities. In *ACM EC* (pp. 78-86).

Indyk, P., & Woodruff, D. (2006). Polylogarithmic Private Approximations and Efficient Matching. In *TCC* (pp. 245-264).

Jacobsen, J. O. B., Baudis, M., Baynam, G. S., et al. (2022). The GA4GH Phenopacket Schema Defines a Computable Representation of Clinical Data. *Nature Biotechnology*, 40(6), 817-820.

Jagielski, M., Ullman, J., & Oprea, A. (2020). Auditing Differentially Private Machine Learning. In *USENIX Security* (pp. 1871-1888).

Jiang, J. J., & Conrath, D. W. (1997). Semantic Similarity Based on Corpus Statistics and Lexical Taxonomy. In *ROCLING* (pp. 19-33).

Johnson, A., & Shmatikov, V. (2013). Privacy-Preserving Data Exploration in Genome-Wide Association Studies. In *ACM KDD* (pp. 1079-1087).

Köhler, S., Carmody, L., Vasilevsky, N., et al. (2019). Expansion of the Human Phenotype Ontology (HPO) Knowledge Base and Resources. *Nucleic Acids Research*, 47(D1), D1018-D1027.

Köhler, S., Gargano, M., Matentzoglu, N., et al. (2021). The Human Phenotype Ontology in 2021. *Nucleic Acids Research*, 49(D1), D1207-D1217.

Köhler, S., Schulz, M. H., Krawitz, P., et al. (2009). Clinical Diagnostics in Human Genetics with Semantic Similarity Searches in Ontologies. *American Journal of Human Genetics*, 85(4), 457-464.

Kolesnikov, V., Matania, N., Pinkas, B., Rosulek, M., & Trieu, N. (2017). Practical Multi-Party Private Set Intersection from Symmetric-Key Techniques. In *ACM CCS* (pp. 1257-1272).

Landrum, M. J., Lee, J. M., Benson, M., et al. (2018). ClinVar: Improving Access to Variant Interpretations and Supporting Evidence. *Nucleic Acids Research*, 46(D1), D1062-D1067.

Li, N., Li, T., & Venkatasubramanian, S. (2007). t-Closeness: Privacy Beyond k-Anonymity and l-Diversity. In *IEEE ICDE* (pp. 106-115).

Lin, D. (1998). An Information-Theoretic Definition of Similarity. In *ICML* (pp. 296-304).

Machanavajjhala, A., Kifer, D., Gehrke, J., & Venkitasubramaniam, M. (2007). l-Diversity: Privacy Beyond k-Anonymity. *ACM TKDD*, 1(1), 3.

McSherry, F., & Talwar, K. (2007). Mechanism Design via Differential Privacy. In *IEEE FOCS* (pp. 94-103).

Meadows, C. (1986). A More Efficient Cryptographic Matchmaking Protocol for Use in the Absence of a Continuously Available Third Party. In *IEEE S&P* (pp. 134-137).

Mungall, C. J., McMurry, J. A., Köhler, S., et al. (2017). The Monarch Initiative: An Integrative Data and Analytic Platform Connecting Phenotypes to Genotypes Across Species. *Nucleic Acids Research*, 45(D1), D712-D722.

Pesquita, C., Faria, D., Falcão, A. O., Lord, P., & Couto, F. M. (2009). Semantic Similarity in Biomedical Ontologies. *PLoS Computational Biology*, 5(7), e1000443.

Philippakis, A. A., Azzariti, D. R., Birney, E., et al. (2015). The Matchmaker Exchange: A Platform for Rare Disease Gene Discovery. *Human Mutation*, 36(10), 915-921.

Pinkas, B., Schneider, T., Weinert, C., & Zohner, M. (2018). Efficient Circuit-Based PSI via Cuckoo Hashing. In *EUROCRYPT* (pp. 125-157).

Pinkas, B., Schneider, T., & Zohner, M. (2014). Faster Private Set Intersection Based on OT Extension. In *USENIX Security* (pp. 797-812).

Pinkas, B., Schneider, T., & Zohner, M. (2018). Scalable Private Set Intersection Based on OT Extension. *ACM TOPS*, 21(2), 1-35.

Raisaro, J. L., Tramèr, F., Ji, Z., et al. (2017). Addressing Beacon Re-Identification Attacks: Quantification and Mitigation of Privacy Risks. *JAMIA*, 24(4), 799-805.

Rath, A., Olry, A., Dhombres, F., et al. (2012). Representation of Rare Diseases in Health Information Systems: The Orphanet Approach to Serve a Wide Range of End Users. *Human Mutation*, 33(5), 803-808.

Resnik, P. (1995). Using Information Content to Evaluate Semantic Similarity in a Taxonomy. In *IJCAI* (pp. 448-453).

Robinson, P. N., Köhler, S., Bauer, S., et al. (2008). The Human Phenotype Ontology: A Tool for Annotating and Analyzing Human Hereditary Disease. *American Journal of Human Genetics*, 83(5), 610-615.

Robinson, P. N., Ravanmehr, V., Jacobsen, J. O. B., et al. (2020). Interpretable Clinical Genomics with a Likelihood Ratio Paradigm. *American Journal of Human Genetics*, 107(3), 403-417.

Sánchez, D., Batet, M., Isern, D., & Valls, A. (2011). Ontology-Based Semantic Similarity: A New Feature-Based Approach. *Expert Systems with Applications*, 39(9), 7718-7728.

Schlicker, A., Domingues, F. S., Rahnenführer, J., & Lengauer, T. (2006). A New Measure for Functional Similarity of Gene Products Based on Gene Ontology. *BMC Bioinformatics*, 7(1), 302.

Seco, N., Veale, T., & Hayes, J. (2004). An Intrinsic Information Content Metric for Semantic Similarity in WordNet. In *ECAI* (pp. 1089-1090).

Shamir, A. (1979). How to Share a Secret. *Communications of the ACM*, 22(11), 612-613.

Shokri, R., Stronati, M., Song, C., & Shmatikov, V. (2017). Membership Inference Attacks Against Machine Learning Models. In *IEEE S&P* (pp. 3-18).

Shringarpure, S. S., & Bustamante, C. D. (2015). Privacy Risks from Genomic Data-Sharing Beacons. *American Journal of Human Genetics*, 97(5), 631-646.

Smedley, D., Jacobsen, J. O. B., Jäger, M., et al. (2015). Next-Generation Diagnostics and Disease-Gene Discovery with the Exomiser. *Nature Protocols*, 10(12), 2004-2015.

Sobreira, N., Schiettecatte, F., Valle, D., & Hamosh, A. (2015). GeneMatcher: A Matching Tool for Connecting Investigators with an Interest in the Same Gene. *Human Mutation*, 36(10), 928-930.

Sweeney, L. (2002). k-Anonymity: A Model for Protecting Privacy. *International Journal of Uncertainty, Fuzziness and Knowledge-Based Systems*, 10(5), 557-570.

Tang, H., Jiang, X., Wang, X., et al. (2020). Protecting Genomic Data Analytics in the Cloud: State of the Art and Opportunities. *BMC Medical Genomics*, 9(1), 63.

Turnbull, C., Scott, R. H., Thomas, E., et al. (2018). The 100 000 Genomes Project: Bringing Whole Genome Sequencing to the NHS. *BMJ*, 361, k1687.

Uhlerop, C., Slavković, A., & Fienberg, S. E. (2013). Privacy-Preserving Data Sharing for Genome-Wide Association Studies. *Journal of Privacy and Confidentiality*, 5(1), 137-166.

Yao, A. C. (1986). How to Generate and Exchange Secrets. In *IEEE FOCS* (pp. 162-167).

Yeom, S., Giacomelli, I., Fredrikson, M., & Jha, S. (2018). Privacy Risk in Machine Learning: Analyzing the Connection to Overfitting. In *IEEE CSF* (pp. 268-282).

Zahur, S., Rosulek, M., & Evans, D. (2015). Two Halves Make a Whole: Reducing Data Transfer in Garbled Circuits Using Half Gates. In *EUROCRYPT* (pp. 220-250).

Zhou, Z., Wang, Y., & Gu, J. (2008). New Model of Semantic Similarity Measuring in WordNet. In *ICIS* (pp. 256-261).
