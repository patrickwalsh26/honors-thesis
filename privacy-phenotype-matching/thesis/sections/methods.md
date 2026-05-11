# Methods

## 3.1 System Overview

We present a privacy-preserving phenotype matching framework that enables federated rare disease patient discovery while protecting sensitive clinical information. Our system integrates three complementary privacy mechanisms—Private Set Intersection (PSI), differential privacy (DP), and k-anonymity—within a modular architecture compatible with established genomic data sharing standards.

The framework operates on patient records encoded as GA4GH Phenopackets (Jacobsen et al., 2022), with phenotypes represented using Human Phenotype Ontology (HPO) terms (Köhler et al., 2021). Given a query patient *Q* and a database of *n* patients *D = {P₁, P₂, ..., Pₙ}*, the system computes privacy-preserving similarity scores *s(Q, Pᵢ)* and returns the top-*k* most similar patients without revealing the precise phenotype overlap or non-matching terms.

Figure 1 illustrates the system architecture. A query phenopacket first undergoes rare term filtering to remove quasi-identifying phenotypes. The filtered query then enters the PSI protocol, which computes set intersection cardinality without revealing individual terms. Differential privacy noise is added to similarity scores, and finally, k-anonymity checks suppress results from insufficiently sized cohorts.

## 3.2 Data Representation

### 3.2.1 GA4GH Phenopackets

We adopt the Global Alliance for Genomics and Health (GA4GH) Phenopacket schema version 2.0 as our primary data representation (Jacobsen et al., 2022). Phenopackets provide a standardized, computable representation of clinical phenotype data that facilitates interoperability across institutions and computational systems.

Each phenopacket *P* contains:

- **Subject**: Patient identifier and demographic information (sex, age at encounter)
- **Phenotypic Features**: A set of observed clinical phenotypes *Φ(P) = {φ₁, φ₂, ..., φₘ}*, where each *φᵢ* is an HPO term
- **Diseases**: Diagnosed conditions with ontology identifiers (OMIM, Orphanet)
- **Metadata**: Provenance information including creation timestamp and contributing resources

Formally, we represent a patient's phenotype profile as the set of HPO term identifiers:

$$\Phi(P) = \{t \in \mathcal{H} : t \text{ is observed in } P\}$$

where $\mathcal{H}$ denotes the set of all HPO terms.

### 3.2.2 Human Phenotype Ontology

The Human Phenotype Ontology (HPO) provides a standardized vocabulary of 18,000+ phenotypic abnormality terms organized as a directed acyclic graph (DAG) (Köhler et al., 2021). Each term *t ∈ H* has:

- A unique identifier (e.g., HP:0001250 for "Seizure")
- A human-readable label and definition
- Hierarchical relationships to parent (more general) and child (more specific) terms

We denote the ancestor set of term *t* as *Anc(t)* and the descendant set as *Desc(t)*. The root term HP:0000001 ("All") is an ancestor of every term in the ontology.

The ontology structure enables semantic similarity computation and supports privacy-preserving generalization, where specific terms can be replaced with more general parent terms to reduce identifiability.

### 3.2.3 Information Content

Information content (IC) quantifies the specificity of an HPO term based on its frequency in a reference corpus (Resnik, 1995). We compute IC using two approaches:

**Corpus-based IC.** Given a corpus of phenopackets *C*, the probability of observing term *t* is:

$$P(t) = \frac{|\{P \in C : t \in \Phi(P) \lor \exists t' \in \Phi(P) : t \in Anc(t')\}|}{|C|}$$

The information content is then:

$$IC(t) = -\log_2 P(t)$$

Terms with lower corpus frequency have higher IC, indicating greater diagnostic specificity. For example, "Seizure" (HP:0001250) appearing in 2,530 of 12,974 diseases has lower IC than "Myoclonic seizure" (HP:0002123) appearing in 412 diseases.

**Intrinsic IC.** When corpus statistics are unavailable, we estimate IC from ontology structure using descendant counts (Sánchez et al., 2011):

$$IC_{intrinsic}(t) = 1 - \frac{\log(|Desc(t)| + 1)}{\log(|\mathcal{H}|)}$$

Leaf terms with no descendants receive maximum IC, while the root receives IC = 0.

## 3.3 Phenotype Similarity Metrics

We implement four similarity metrics spanning set-theoretic, vector-space, and semantic approaches. All metrics produce scores in [0, 1], where 1 indicates identical phenotype profiles.

### 3.3.1 Jaccard Similarity

The Jaccard index measures phenotype overlap as the ratio of intersection to union cardinality:

$$J(\Phi(Q), \Phi(P)) = \frac{|\Phi(Q) \cap \Phi(P)|}{|\Phi(Q) \cup \Phi(P)|}$$

Jaccard similarity treats all HPO terms as equally informative and requires exact term matches. While computationally efficient, it does not account for semantic relationships between related terms.

### 3.3.2 Cosine Similarity with IC Weighting

Cosine similarity represents phenotype profiles as vectors in a high-dimensional space indexed by HPO terms. We weight each dimension by information content to emphasize rare, diagnostically specific phenotypes:

$$\cos_{IC}(\Phi(Q), \Phi(P)) = \frac{\sum_{t \in \Phi(Q) \cap \Phi(P)} IC(t)^2}{\sqrt{\sum_{t \in \Phi(Q)} IC(t)^2} \cdot \sqrt{\sum_{t \in \Phi(P)} IC(t)^2}}$$

When IC weights are unavailable, we use uniform weighting (IC(t) = 1 for all t), reducing to standard cosine similarity on binary term vectors.

### 3.3.3 Resnik Similarity

Resnik similarity leverages the ontology structure to compare terms based on their most informative common ancestor (MICA) (Resnik, 1995). For two terms *t₁* and *t₂*:

$$sim_{Resnik}(t_1, t_2) = IC(MICA(t_1, t_2))$$

where:

$$MICA(t_1, t_2) = \arg\max_{t \in Anc(t_1) \cap Anc(t_2)} IC(t)$$

To extend pairwise term similarity to phenotype sets, we employ the Best-Match Average (BMA) strategy (Pesquita et al., 2009):

$$BMA(\Phi(Q), \Phi(P)) = \frac{1}{2}\left(\frac{1}{|\Phi(Q)|}\sum_{q \in \Phi(Q)} \max_{p \in \Phi(P)} sim(q, p) + \frac{1}{|\Phi(P)|}\sum_{p \in \Phi(P)} \max_{q \in \Phi(Q)} sim(p, q)\right)$$

BMA computes the average best match in both directions, ensuring symmetric similarity scores.

### 3.3.4 Simplified Resnik Similarity

For computational efficiency in privacy-preserving protocols, we implement a simplified Resnik variant that approximates semantic similarity using only exact term matches weighted by IC:

$$sim_{simple}(\Phi(Q), \Phi(P)) = \frac{\sum_{t \in \Phi(Q) \cap \Phi(P)} IC(t)}{\max\left(\sum_{t \in \Phi(Q)} IC(t), \sum_{t \in \Phi(P)} IC(t)\right)}$$

This formulation avoids ontology traversal during similarity computation while preserving the emphasis on rare phenotypes through IC weighting.

## 3.4 Privacy-Preserving Mechanisms

Our framework implements three complementary privacy mechanisms that can be composed to provide layered protection. Each mechanism addresses distinct privacy concerns and offers configurable parameters to balance utility and protection.

### 3.4.1 Private Set Intersection (PSI)

Private Set Intersection enables two parties to compute the intersection of their sets without revealing elements outside the intersection (Meadows, 1986). We implement Diffie-Hellman-based PSI using elliptic curve cryptography.

**Protocol Description.** Let *Q* denote the query party with phenotype set *Φ(Q)* and *S* denote the server with database phenotypes. The protocol proceeds as follows:

1. **Setup**: Both parties agree on elliptic curve parameters (NIST P-256/secp256r1) and a hash function *H: {0,1}* → G* mapping strings to curve points.

2. **Query Blinding**: The query party samples a random scalar *α ← Z_q* and computes blinded query elements:
   $$\forall t \in \Phi(Q): B_t = H(t)^\alpha$$

3. **Server Processing**: For each blinded element *B_t*, the server samples *β ← Z_q* and computes:
   $$\forall t \in \Phi(Q): C_t = B_t^\beta = H(t)^{\alpha\beta}$$

   The server also computes its own blinded elements:
   $$\forall t' \in \Phi(S): D_{t'} = H(t')^\beta$$

4. **Query Completion**: The query party raises server elements to its secret:
   $$\forall t' \in \Phi(S): E_{t'} = D_{t'}^\alpha = H(t')^{\alpha\beta}$$

5. **Intersection**: Elements appearing in both *{C_t}* and *{E_{t'}}* correspond to the intersection *Φ(Q) ∩ Φ(S)*.

**Security Guarantees.** Under the Decisional Diffie-Hellman (DDH) assumption, this protocol is secure in the semi-honest (honest-but-curious) adversary model. Neither party learns elements outside the intersection.

**Cardinality-Only Mode.** For enhanced privacy, we support a cardinality-only variant that reveals only |*Φ(Q) ∩ Φ(S)*| without disclosing which specific terms match. This enables Jaccard similarity computation:

$$J = \frac{|Φ(Q) \cap Φ(S)|}{|Φ(Q)| + |Φ(S)| - |Φ(Q) \cap Φ(S)|}$$

### 3.4.2 Differential Privacy

Differential privacy provides a rigorous mathematical framework for quantifying privacy loss (Dwork et al., 2006). A randomized mechanism *M* satisfies (ε, δ)-differential privacy if for all datasets *D₁*, *D₂* differing in one record and all output sets *S*:

$$\Pr[M(D_1) \in S] \leq e^\varepsilon \cdot \Pr[M(D_2) \in S] + \delta$$

The privacy parameter ε (epsilon) controls the privacy-utility tradeoff, with smaller ε providing stronger privacy but requiring more noise.

**Laplace Mechanism.** For a numeric query *f* with sensitivity *Δf = max_{D₁,D₂} |f(D₁) - f(D₂)|*, the Laplace mechanism achieves (ε, 0)-DP by adding noise:

$$M(D) = f(D) + Lap\left(\frac{\Delta f}{\varepsilon}\right)$$

where *Lap(b)* denotes a random variable drawn from the Laplace distribution with scale *b*.

We apply the Laplace mechanism to privatize similarity scores. For similarity functions bounded in [0, 1], the sensitivity is Δf = 1, yielding:

$$\tilde{s}(Q, P) = s(Q, P) + Lap\left(\frac{1}{\varepsilon}\right)$$

**Gaussian Mechanism.** For (ε, δ)-DP with δ > 0, the Gaussian mechanism adds normally distributed noise:

$$M(D) = f(D) + \mathcal{N}\left(0, \sigma^2\right)$$

where:

$$\sigma = \frac{\Delta f \cdot \sqrt{2\ln(1.25/\delta)}}{\varepsilon}$$

The Gaussian mechanism often provides better utility for the same privacy guarantee when δ is acceptably small (e.g., δ = 10⁻⁵).

**Exponential Mechanism.** For selecting among discrete options (e.g., top-k patients), the exponential mechanism samples outputs proportionally to their utility:

$$\Pr[M(D) = r] \propto \exp\left(\frac{\varepsilon \cdot u(D, r)}{2\Delta u}\right)$$

where *u(D, r)* is a utility function and *Δu* is its sensitivity.

**Privacy Accounting.** When multiple queries are issued, privacy degrades according to composition theorems. We implement a privacy accountant that tracks cumulative ε expenditure and enforces budget limits.

### 3.4.3 K-Anonymity and Rare Term Filtering

K-anonymity requires that each record be indistinguishable from at least *k-1* other records with respect to quasi-identifying attributes (Sweeney, 2002). In phenotype matching, rare HPO term combinations can serve as quasi-identifiers that uniquely identify patients.

**Rare Term Filtering.** We mitigate quasi-identifier risk by filtering phenotypes that appear infrequently in the reference corpus. For a prevalence threshold *τ* and corpus *C*:

$$\Phi'(P) = \{t \in \Phi(P) : freq(t, C) \geq \tau\}$$

where *freq(t, C)* is the number of patients in *C* with term *t* or any of its descendants.

We implement two filtering strategies:

1. **Suppression**: Remove rare terms entirely from the query
2. **Generalization**: Replace rare terms with their most specific ancestor having sufficient prevalence

**Result Suppression.** The k-anonymity guard suppresses query results when fewer than *k* patients match, preventing inference about small cohorts:

$$R(Q) = \begin{cases} \text{top-}k \text{ results} & \text{if } |matches| \geq k \\ \varnothing & \text{otherwise} \end{cases}$$

### 3.4.4 Composed Privacy Pipeline

Our framework composes the three mechanisms in a configurable pipeline:

1. **Rare Term Filtering** (pre-processing): Remove or generalize quasi-identifying phenotypes
2. **PSI Computation**: Securely compute phenotype overlap
3. **Similarity Scoring**: Calculate similarity metrics from PSI output
4. **DP Noise Addition**: Perturb similarity scores with calibrated noise
5. **K-Anonymity Check**: Suppress results below cohort size threshold

Each mechanism can be independently enabled or disabled via configuration. The privacy accountant tracks cumulative privacy loss across all mechanisms.

## 3.5 Evaluation Framework

### 3.5.1 Retrieval Metrics

We evaluate phenotype matching as an information retrieval task where, for each query patient *Q*, the system ranks all database patients by similarity. Ground truth relevance is defined as patients sharing the same underlying disease diagnosis.

**Precision at k (P@k).** The fraction of top-*k* retrieved patients that are relevant:

$$P@k = \frac{|\{P \in top_k : relevant(P)\}|}{k}$$

**Recall at k (R@k).** The fraction of all relevant patients retrieved in the top-*k*:

$$R@k = \frac{|\{P \in top_k : relevant(P)\}|}{|relevant|}$$

**Normalized Discounted Cumulative Gain (nDCG@k).** A ranking-aware metric that rewards placing relevant results higher:

$$DCG@k = \sum_{i=1}^{k} \frac{rel_i}{\log_2(i+1)}$$

$$nDCG@k = \frac{DCG@k}{IDCG@k}$$

where *IDCG@k* is the DCG of the ideal ranking.

**Mean Reciprocal Rank (MRR).** The average reciprocal position of the first relevant result:

$$MRR = \frac{1}{|Q|}\sum_{q \in Q} \frac{1}{rank_q}$$

### 3.5.2 Privacy Metrics

We quantify privacy protection through adversarial attack simulations.

**Membership Inference.** An adversary attempts to determine whether a target patient is in the database based on query responses (Shokri et al., 2017). We measure attack success rate as:

$$Adv_{MI} = Accuracy_{attack} - 0.5$$

where 0.5 represents random guessing.

**Attribute Inference.** An adversary attempts to infer undisclosed phenotypes from partial query responses. We measure:

$$Adv_{AI} = \frac{TP + TN}{TP + TN + FP + FN} - Prior$$

where *Prior* is the baseline accuracy from phenotype prevalence.

### 3.5.3 Privacy-Utility Frontier

We characterize the privacy-utility tradeoff by sweeping privacy parameters and plotting retrieval performance against privacy loss:

- **DP Sweep**: ε ∈ {0.1, 0.5, 1.0, 2.0, 5.0, 10.0, ∞}
- **K-Anonymity Sweep**: k ∈ {2, 5, 10, 20, 50}
- **Rare Term Threshold**: τ ∈ {0.001, 0.005, 0.01, 0.05}

The Pareto frontier identifies optimal configurations that maximize utility for a given privacy budget.

## 3.6 Experimental Dataset

### 3.6.1 HPO Annotations Corpus

We constructed our evaluation dataset from the HPO phenotype annotation file (phenotype.hpoa), which contains curated disease-phenotype associations from multiple authoritative sources (Köhler et al., 2021).

**Table 1: HPO Annotations Dataset Summary**

| Source | Diseases | Annotations | Mean Phenotypes/Disease |
|--------|----------|-------------|------------------------|
| OMIM | 8,592 | 174,519 | 20.3 |
| Orphanet | 4,335 | 88,051 | 20.3 |
| DECIPHER | 47 | 942 | 20.0 |
| **Total** | **12,974** | **263,512** | **20.3** |

The corpus contains 11,514 unique HPO terms spanning the full breadth of phenotypic abnormalities.

### 3.6.2 Synthetic Patient Generation

We generated synthetic patient cohorts by sampling from real disease phenotype profiles. This approach produces realistic phenotype distributions while avoiding the need for real patient data, which would require institutional review board approval and data use agreements.

**Generation Procedure:**

1. **Disease Selection**: Sample *d* diseases from HPOA, weighted by number of annotations
2. **Patient Assignment**: Assign *n/d* patients to each disease (balanced design)
3. **Phenotype Sampling**: For each patient, sample phenotypes from the disease profile with recall rate *ρ* (default: 0.75)
4. **Noise Addition**: Add *η* fraction (default: 0.10) of random phenotypes from the corpus
5. **Ground Truth**: Patients with the same underlying disease are labeled as relevant

**Table 2: Evaluation Cohort Characteristics**

| Parameter | Value |
|-----------|-------|
| Total patients | 500 |
| Diseases sampled | 100 |
| Patients per disease | 5 |
| Phenotype recall (ρ) | 0.75 |
| Noise rate (η) | 0.10 |
| Unique phenotypes | 1,145 |
| Mean phenotypes/patient | 10.4 ± 5.9 |
| Disease source distribution | 69% OMIM, 31% Orphanet |

The resulting cohort exhibits phenotype distributions consistent with real rare disease patients, with each patient having 3-24 observed phenotypes (median: 9).

## 3.7 Implementation

We implemented the privacy-preserving phenotype matching framework in Python 3.10. The system comprises approximately 4,500 lines of code organized into modular components.

**Core Dependencies:**

- `cryptography` (v41.0+): Elliptic curve operations for PSI (ECDH, HKDF)
- `pronto` (v2.5+): HPO ontology parsing and traversal
- `numpy` (v1.24+): Numerical computation and similarity matrices
- `scipy` (v1.11+): Statistical distributions for DP mechanisms

**Standards Compliance:**

- GA4GH Phenopackets v2.0: Native JSON schema support
- GA4GH Beacon v2: Query/response adapters for federation
- Matchmaker Exchange v1.1: Export compatibility for existing networks

**Computational Complexity:**

- Similarity matrix computation: O(n²m) for *n* patients with *m* average phenotypes
- PSI protocol: O(m) exponentiations per patient pair
- DP noise addition: O(1) per score
- Overall query processing: O(n · m) for ranking against *n* database patients

The complete implementation, evaluation scripts, and synthetic datasets are available at [repository URL] under an open-source license.

---

## References

Dwork, C., McSherry, F., Nissim, K., & Smith, A. (2006). Calibrating Noise to Sensitivity in Private Data Analysis. In *Theory of Cryptography Conference* (pp. 265-284). Springer.

Jacobsen, J. O. B., Baudis, M., Baynam, G. S., Beckmann, J. S., Beltran, S., Buske, O. J., ... & Robinson, P. N. (2022). The GA4GH Phenopacket schema defines a computable representation of clinical data. *Nature Biotechnology*, 40(6), 817-820.

Köhler, S., Gargano, M., Matentzoglu, N., Carmody, L. C., Lewis-Smith, D., Vasilevsky, N. A., ... & Robinson, P. N. (2021). The Human Phenotype Ontology in 2021. *Nucleic Acids Research*, 49(D1), D1207-D1217.

Meadows, C. (1986). A More Efficient Cryptographic Matchmaking Protocol for Use in the Absence of a Continuously Available Third Party. In *IEEE Symposium on Security and Privacy* (pp. 134-137).

Pesquita, C., Faria, D., Falcão, A. O., Lord, P., & Couto, F. M. (2009). Semantic Similarity in Biomedical Ontologies. *PLoS Computational Biology*, 5(7), e1000443.

Resnik, P. (1995). Using Information Content to Evaluate Semantic Similarity in a Taxonomy. In *Proceedings of the 14th International Joint Conference on Artificial Intelligence* (pp. 448-453).

Sánchez, D., Batet, M., Isern, D., & Valls, A. (2011). Ontology-based semantic similarity: A new feature-based approach. *Expert Systems with Applications*, 39(9), 7718-7728.

Shokri, R., Stronati, M., Song, C., & Shmatikov, V. (2017). Membership Inference Attacks Against Machine Learning Models. In *IEEE Symposium on Security and Privacy* (pp. 3-18).

Sweeney, L. (2002). k-Anonymity: A Model for Protecting Privacy. *International Journal of Uncertainty, Fuzziness and Knowledge-Based Systems*, 10(5), 557-570.
