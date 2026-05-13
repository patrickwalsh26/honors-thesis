# Literature Review

This chapter positions our contributions across four bodies of work: phenotype standards and disease databases (§2.1), semantic-similarity measures (§2.2), federated patient-matching systems (§2.3), and privacy-preserving computation with its attack surface (§2.4–§2.5). Section 2.7 names the gaps our work closes.

## 2.1 Rare-Disease Phenotyping and Standards

The **Human Phenotype Ontology** (Robinson et al., 2008; Köhler et al., 2021) provides the de facto vocabulary for computational rare-disease phenotyping. The current release contains roughly 19,000 phenotypic-abnormality terms in a DAG rooted at "Phenotypic abnormality" (HP:0000118), with auxiliary branches for inheritance, modifier, frequency, and clinical course. Ancestor closure is load-bearing for semantic similarity: a patient annotated "Focal clonic seizure" (HP:0002266) implicitly also matches "Seizure" (HP:0001250). HPO underlies clinical tools including Phenomizer (Köhler et al., 2009), Exomiser (Smedley et al., 2015), and LIRICAL (Robinson et al., 2020).

The **GA4GH Phenopacket** v2 schema (Jacobsen et al., 2022) standardises patient representation (subject, phenotypic features with observed/excluded annotation, diagnoses, biosample, interpretation) as a Protobuf/JSON message designed for federated exchange. Adoption spans the 100,000 Genomes Project, ClinVar, and the Matchmaker Exchange network. We adopt Phenopackets natively (§3.2.1).

**Disease databases.** OMIM (Amberger et al., 2019) and Orphanet (Rath et al., 2012) curate disease–phenotype associations from primary literature; DECIPHER (Firth et al., 2009) focuses on copy-number variants; ClinVar (Landrum et al., 2018) aggregates variant interpretations. The HPO project's `phenotype.hpoa` integrates these into roughly 263,000 disease–phenotype annotations linking 13,000 diseases to 11,500 HPO terms; this is the corpus from which both synthetic patients (§3.6.2) and corpus-IC priors (§3.2.3) are derived. For real-patient evaluation we additionally use the Monarch Phenopacket Store (Danis et al., 2025), a curated 9,588-phenopacket corpus drawn from published case reports with confirmed OMIM diagnoses.

## 2.2 Semantic-Similarity Measures

Phenotype matching reduces to a similarity computation over HPO term sets. Three classes of measure span the relevant design space.

**Set-theoretic measures** (Jaccard 1901; Dice 1945; overlap coefficient) treat term sets as sets and weight terms uniformly. They are computationally trivial and interpretable but ignore the diagnostic value differential between rare and common phenotypes.

**Information-theoretic measures** weight matches by term specificity. **Resnik similarity** (Resnik, 1995) sets pairwise similarity to the IC of the most informative common ancestor, $\mathrm{IC}(\text{MICA}(t_1, t_2))$; **Lin** (1998) and **Jiang–Conrath** (1997) provide normalised variants. IC may be computed from a reference corpus or estimated intrinsically from ontology structure (Seco et al., 2004; Sánchez et al., 2011). Aggregation to term sets uses **Best-Match Average** (BMA), the symmetric average of forward and reverse best-match similarities, which performs strongly in empirical surveys (Pesquita et al., 2009) and is used by Phenomizer.

**Vector-space measures** with IC weighting (Methods §3.3.2) trade ontology traversal cost for vector arithmetic. Our §4.6 real-cohort comparison shows that IC-weighted cosine is essentially tied with full Resnik+BMA over the HPO DAG when IC priors are estimated from the cohort itself.

## 2.3 Federated Patient-Matching Systems

The **Matchmaker Exchange** (MME; Philippakis et al., 2015) connects seven major rare-disease databases via a standardised query API: GeneMatcher (Sobreira et al., 2015), PhenomeCentral (Buske et al., 2015), DECIPHER, MyGene2, PatientMatcher, seqr, and Matchbox. MME has facilitated thousands of matches and contributed to over 100 disease-gene discoveries, including TBCK-related intellectual disability (Bhoj et al., 2016). PhenomeCentral implements semantic-similarity matching with a Resnik-derived HRSS measure; GeneMatcher uses gene-centric matching. The **GA4GH Beacon Network** (Fiume et al., 2019) extends the federated model to variant-presence and aggregate-count queries.

All current MME implementations exchange patient phenotype data in cleartext between participating institutions. This privacy boundary is the principal obstacle to broader institutional participation, and it is the design point of this thesis.

## 2.4 Privacy-Preserving Computation

### 2.4.1 Private Set Intersection

PSI lets two parties compute $A \cap B$ without revealing elements outside the intersection. The earliest practical protocols use Diffie–Hellman key exchange (Meadows, 1986; Huberman et al., 1999) and are secure under the Decisional Diffie–Hellman assumption against semi-honest adversaries, with $O(|A|+|B|)$ communication. Later OT-extension protocols (Pinkas et al., 2014, 2018) and circuit-based variants (Huang et al., 2012) achieve practical performance for million-element sets and stronger security models. **Cardinality-only PSI** (De Cristofaro & Tsudik, 2010) reveals only $|A \cap B|$, sufficient for Jaccard similarity. **Threshold PSI** (Ghosh & Simkin, 2019) and **multi-party PSI** (Kolesnikov et al., 2017) extend the primitive to richer query semantics. We use DH-PSI on NIST P-256 (§3.4.1) for its simplicity and adequate performance at our scale.

### 2.4.2 Differential Privacy

Differential privacy (Dwork et al., 2006) bounds the multiplicative change in any output's probability between adjacent datasets by $e^\varepsilon$ (plus an additive $\delta$). The **Laplace mechanism** adds noise of scale $\Delta f/\varepsilon$ to a numeric query with sensitivity $\Delta f$; the **Gaussian mechanism** gives $(\varepsilon, \delta)$-DP at typically lower variance for high-dimensional outputs. The **exponential mechanism** (McSherry & Talwar, 2007) samples from a discrete output set with probability $\propto \exp(\varepsilon \cdot u / 2\Delta u)$ for a utility function $u$; the choice of $u$ is the central design decision and dictates the mechanism's behaviour on compressed-signal regimes (§4.7). Composition theorems (Kairouz et al., 2015; Abadi et al., 2016) bound cumulative privacy loss across multiple mechanisms. DP has seen production deployment in the U.S. Census Bureau (Abowd, 2018) and consumer telemetry pipelines (Erlingsson et al., 2014), and biomedical applications including DP-GWAS (Johnson & Shmatikov, 2013) and DP-Beacons (Raisaro et al., 2017).

### 2.4.3 k-Anonymity and Generalisation

k-anonymity (Sweeney, 2002) requires each released record to be indistinguishable from at least $k-1$ others over quasi-identifying attributes; achievement uses generalisation (replacing values with more general ones) and suppression (removing them). The HPO DAG provides a natural generalisation lattice: "Tetralogy of Fallot" → "Conotruncal defect" → "Abnormal heart morphology." Refinements address known weaknesses: $\ell$-diversity (Machanavajjhala et al., 2007) for homogeneity attacks and $t$-closeness (Li et al., 2007) for distribution-based attacks.

### 2.4.4 Secure Multi-Party Computation

Generic MPC (Yao, 1986; Goldreich et al., 1987) computes arbitrary functions over private inputs and subsumes PSI as a special case. Modern instantiations (garbled circuits with half-gates and free-XOR, Zahur et al. 2015; secret-sharing protocols like SPDZ, Damgård et al. 2012; homomorphic-encryption schemes BGV and CKKS, Cheon et al. 2017) offer richer functionality at higher cost. We use the simpler PSI primitive because phenotype set intersection is the only computation our threat model demands.

## 2.5 Privacy Attacks and Auditing

**Beacon attacks** (Shringarpure & Bustamante, 2015) demonstrated that Boolean variant-presence queries enable membership inference with ~5,000 queries; Raisaro et al. (2017) responded with DP-protected Beacons. **Membership-inference attacks** more generally (Shokri et al., 2017) train shadow models on known-membership data and apply the learned classifier to a target. **Attribute-inference attacks** (Yeom et al., 2018) predict undisclosed attributes from query outputs and auxiliary data. **Reconstruction attacks** (Dinur & Nissim, 2003) recover full records from sufficiently many linear queries: the foundational result that motivated differential privacy.

**Empirical privacy auditing** (Jagielski et al., 2020) measures actual leakage through controlled experiments rather than relying on theoretical guarantees alone. We adopt this approach in §4.5, instantiating Yeom-threshold and Shokri-shadow MI attacks against our DP score-release oracle and a rare-term singling-out attack against the k-anonymity gate.

## 2.6 Privacy-Preserving Genomic Systems

Several systems apply privacy-preserving computation to specific genomic problems. iDASH challenges have driven secure GWAS, private sequence comparison, and encrypted ML (Tang et al., 2020). DP-GWAS uses the exponential mechanism for SNP selection (Johnson & Shmatikov, 2013); follow-on work addresses chi-squared tests (Uhlerop et al., 2013), regression (Chen et al., 2019), and federated meta-analysis (Cho et al., 2018). Encrypted-storage and access-control approaches are surveyed by Ayday et al. (2013). The GREGoR Consortium (2025), which our advisor co-leads, is generating multi-omic data on thousands of rare-disease families and represents a natural deployment target for the framework presented here.

## 2.7 Gaps Addressed by This Thesis

The above bodies of work establish the primitives (PSI, DP, k-anonymity), the data standards (HPO, Phenopackets), and the attack vectors (Shringarpure–Bustamante, Shokri, Yeom). Four gaps remain:

1. **No integrated composition with formal threat-model analysis.** Prior work studies PSI, DP, and k-anonymity in isolation; their composition for phenotype matching has not been specified with the per-step disclosure analysis that distinguishes a deployable system from a collection of primitives. We close this in §3.1.2.
2. **Synthetic-cohort over-optimism.** Existing privacy-utility evaluations rely overwhelmingly on disease-profile-sampled cohorts; we are aware of no prior work that benchmarks rare-disease matching privacy on a real published-patient cohort with confirmed diagnoses. We close this in §4.6.
3. **Score-magnitude pathology of standard DP mechanisms.** The default Laplace-on-similarity composition fails on real cohorts because same-disease score gaps are compressed; the literature does not measure this gap directly. We diagnose it in §5.1.5 and demonstrate the rank-utility exponential-mechanism fix in §4.7.
4. **No clinician-facing demonstration.** Privacy-preserving genomic primitives have not, to our knowledge, been packaged as an interactive deployment that a clinician could exercise without further engineering. We close this in §3.8 with the pilot system.
