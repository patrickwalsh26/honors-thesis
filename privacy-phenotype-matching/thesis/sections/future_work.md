# Future Work

This chapter outlines extensions and future research directions that build upon the foundation established in this thesis. We organize these into technical advances, system improvements, validation studies, and standards development.

## 6.1 Technical Extensions

### 6.1.1 Malicious-Secure Protocols

Our current PSI implementation assumes semi-honest adversaries who follow the protocol but attempt to learn from observations. Extending to malicious security—where adversaries may arbitrarily deviate from the protocol—would strengthen guarantees for high-stakes deployments.

Malicious-secure PSI protocols exist based on committed oblivious transfer (Pinkas et al., 2018) and cut-and-choose techniques (Lindell, 2017). These protocols incur 2–5× computational overhead but provide security against active attacks. Integration with our framework would involve replacing the Diffie-Hellman PSI component while preserving the overall pipeline architecture.

Additionally, verifiable computation techniques could enable clients to verify that servers correctly executed matching algorithms, detecting cheating without revealing inputs. Zero-knowledge proofs of correct execution are an active research area with improving practicality (Bünz et al., 2018).

### 6.1.2 Homomorphic Encryption for Similarity

Fully homomorphic encryption (FHE) enables computation on encrypted data, potentially allowing similarity computation without any information leakage beyond the final result. Recent FHE schemes (CKKS for approximate arithmetic, BGV for exact computation) have achieved practical performance for certain applications (Cheon et al., 2017).

An FHE-based phenotype matching system would:
1. Encrypt phenotype profiles under a common key
2. Compute similarity using homomorphic operations (additions, multiplications)
3. Decrypt only the final similarity scores

Challenges include encoding phenotype sets efficiently, implementing IC-weighted similarity within FHE constraints, and managing computational overhead. Preliminary work suggests that approximate similarity (e.g., locality-sensitive hashing) may be more tractable than exact computation. This direction merits dedicated investigation.

### 6.1.3 Federated Learning for Phenotype Matching

Machine learning approaches to phenotype matching could learn optimal similarity functions from data, potentially outperforming hand-crafted metrics. Federated learning enables model training across distributed datasets without centralizing raw data (McMahan et al., 2017).

A federated phenotype matching model would:
1. Train local models at each institution on their patient cohorts
2. Aggregate model updates centrally without sharing patient data
3. Deploy the global model for privacy-preserving inference

Differential privacy can be incorporated during aggregation (DP-SGD; Abadi et al., 2016), providing formal privacy guarantees for the training process. The resulting model could then be used for efficient, privacy-preserving matching at scale.

### 6.1.4 Multi-Party Computation for Coalition Queries

Current PSI operates between two parties. Extending to multi-party settings would enable coalition queries where a patient is matched against multiple institutions simultaneously without any institution learning about others' data.

Multi-party PSI protocols (Kolesnikov et al., 2017) can compute the intersection across n parties. Secure aggregation protocols can combine similarity scores from multiple sources without revealing individual contributions. These techniques would enable true federated matching across the Matchmaker Exchange network.

### 6.1.5 Ontology-Aware Privacy

The HPO ontology provides natural generalization pathways that could be more deeply integrated into privacy mechanisms. Semantic generalization—replacing specific terms with ancestors—could provide privacy while preserving similarity computation using the generalized terms.

Research questions include:
- Optimal generalization levels that balance privacy and utility
- Automated selection of generalization based on term frequencies
- Integration of ontology structure into DP mechanisms (e.g., exponential mechanism over ontology nodes)

This direction leverages the unique structure of phenotype data that distinguishes it from generic set matching.

## 6.2 System Improvements

### 6.2.1 Real-Time Privacy Budget Management

Production deployment requires sophisticated privacy budget management beyond our current implementation. Enhancements would include:

- **Fine-grained accounting**: Track budget per patient, per institution, and per query type
- **Adaptive mechanisms**: Adjust ε based on remaining budget and query importance
- **Budget forecasting**: Project budget consumption to inform policy decisions
- **Renewal policies**: Implement budget replenishment for long-term participation

Integration with institutional identity management systems would enable per-user budget allocation and audit trails.

### 6.2.2 Scalable Indexing

Our current implementation computes pairwise similarities, scaling quadratically with database size. For large databases (>10,000 patients), this becomes impractical.

Locality-sensitive hashing (LSH) provides approximate nearest neighbor search with sub-linear query time. Our MinHash-LSH implementation (Section 3.3) offers a foundation, but deeper integration with privacy mechanisms is needed:

- Privacy-preserving LSH index construction
- DP-protected candidate retrieval
- Accuracy-privacy tradeoffs for approximate matching

Tree-based indices (KD-trees, ball trees) and learned indices are additional directions.

### 6.2.3 Clinical Decision Support Integration

Effective deployment requires integration with clinical workflows. Future work includes:

- **EHR integration**: Plugins for Epic, Cerner, and other EHR systems
- **Phenopacket extraction**: Automated conversion from clinical notes to structured phenotypes via NLP
- **Result presentation**: Clinician-friendly interfaces for match review
- **Outcome tracking**: Integration with diagnostic workflows to measure clinical impact

User studies with clinicians would inform interface design and identify workflow integration points.

### 6.2.4 Multi-Modal Matching

Rare disease diagnosis increasingly relies on multi-omic data. Extending our framework to incorporate genomic, transcriptomic, and other data types would enhance matching power:

- **Gene-phenotype joint matching**: Combine phenotype similarity with gene/variant overlap
- **Expression outlier integration**: Incorporate RNA-seq outlier status (Frésard et al., 2019)
- **Imaging phenotypes**: Extend HPO-based matching to radiological and pathological features

Each modality introduces distinct privacy considerations requiring tailored mechanisms.

## 6.3 Validation Studies

### 6.3.1 Real Patient Cohort Evaluation

Validation on actual patient cohorts is essential to confirm clinical utility. This requires:

- **IRB approval**: Protocol development for privacy-preserving research
- **Institutional partnerships**: Collaboration with rare disease centers
- **Retrospective evaluation**: Apply matching to historical cohorts with known diagnoses
- **Blinded assessment**: Clinician evaluation of match quality without diagnostic labels

The Undiagnosed Diseases Network (UDN) and GREGoR Consortium represent potential validation partners with large, characterized cohorts.

### 6.3.2 Prospective Clinical Trial

Ultimately, demonstrating that privacy-preserving matching improves diagnostic outcomes requires prospective evaluation:

- **Study design**: Randomized comparison of standard care vs. matching-augmented diagnosis
- **Endpoints**: Time to diagnosis, diagnostic accuracy, healthcare utilization
- **Patient populations**: Undiagnosed rare disease patients at participating centers
- **Privacy monitoring**: Ongoing assessment of privacy mechanism performance

Such a trial would require multi-year commitment and substantial funding but would provide definitive evidence of clinical value.

### 6.3.3 Cross-Institutional Federation Pilot

Testing federated deployment across real institutions would validate practical feasibility:

- **Partner recruitment**: 3–5 institutions with rare disease programs
- **Infrastructure deployment**: Privacy-preserving matching nodes at each site
- **Query protocols**: Standardized procedures for cross-institutional matching
- **Governance framework**: Data use agreements, IRB coordination, dispute resolution

Lessons learned would inform broader deployment through Matchmaker Exchange.

### 6.3.4 Adversarial Red Team Evaluation

Security evaluation by external red teams would stress-test privacy mechanisms:

- **Attack simulation**: Sophisticated attacks beyond membership/attribute inference
- **Implementation audit**: Code review for cryptographic vulnerabilities
- **Protocol analysis**: Formal verification of security properties
- **Penetration testing**: Attempt to extract protected information

This evaluation would identify weaknesses before production deployment.

## 6.4 Standards and Interoperability

### 6.4.1 GA4GH Privacy Extensions

Proposing privacy extensions to GA4GH standards would promote ecosystem-wide adoption:

- **Phenopacket privacy fields**: Metadata indicating privacy mechanism parameters
- **Beacon v3 privacy modes**: Standardized DP and PSI query types
- **MME privacy handshake**: Protocol for negotiating privacy parameters between nodes

Engagement with GA4GH working groups would advance these proposals through the standards process.

### 6.4.2 FHIR Integration

HL7 FHIR (Fast Healthcare Interoperability Resources) is increasingly the standard for clinical data exchange. FHIR integration would enable:

- **Phenotype extraction**: Map FHIR Observations/Conditions to HPO terms
- **Matching as FHIR service**: Expose matching via FHIR REST APIs
- **Result representation**: Return matches as FHIR resources

This integration would facilitate adoption within healthcare settings already using FHIR.

### 6.4.3 International Regulatory Harmonization

Privacy regulations vary across jurisdictions. Future work should address:

- **GDPR compliance**: Detailed analysis of DP as "appropriate technical measures"
- **Cross-border data flows**: Privacy-preserving computation as transfer mechanism
- **Regulatory engagement**: Dialogue with FDA, EMA, and other authorities on privacy-preserving diagnostics

Harmonized guidance would reduce deployment uncertainty.

## 6.5 Summary

Future work spans multiple dimensions: strengthening cryptographic foundations, scaling to production workloads, validating clinical utility, and advancing interoperability standards. The most impactful near-term directions are:

1. **Real patient validation** to confirm clinical utility
2. **Multi-party protocols** to enable true federation
3. **Clinical integration** for practical deployment
4. **Standards proposals** for ecosystem adoption

These efforts would transform privacy-preserving phenotype matching from a research prototype to a production tool that advances rare disease diagnosis while protecting patient confidentiality.

---

## References

Abadi, M., Chu, A., Goodfellow, I., et al. (2016). Deep Learning with Differential Privacy. In *ACM CCS* (pp. 308-318).

Bünz, B., Bootle, J., Boneh, D., et al. (2018). Bulletproofs: Short Proofs for Confidential Transactions and More. In *IEEE S&P* (pp. 315-334).

Cheon, J. H., Kim, A., Kim, M., & Song, Y. (2017). Homomorphic Encryption for Arithmetic of Approximate Numbers. In *ASIACRYPT* (pp. 409-437).

Frésard, L., Smail, C., Ferraro, N. M., et al. (2019). Identification of Rare-Disease Genes Using Blood Transcriptome Sequencing and Large Control Cohorts. *Nature Medicine*, 25(6), 911-919.

Kolesnikov, V., Matania, N., Pinkas, B., Rosulek, M., & Trieu, N. (2017). Practical Multi-Party Private Set Intersection from Symmetric-Key Techniques. In *ACM CCS* (pp. 1257-1272).

Lindell, Y. (2017). How to Simulate It – A Tutorial on the Simulation Proof Technique. In *Tutorials on the Foundations of Cryptography* (pp. 277-346). Springer.

McMahan, B., Moore, E., Ramage, D., Hampson, S., & y Arcas, B. A. (2017). Communication-Efficient Learning of Deep Networks from Decentralized Data. In *AISTATS* (pp. 1273-1282).

Pinkas, B., Schneider, T., Weinert, C., & Zohner, M. (2018). Efficient Circuit-Based PSI via Cuckoo Hashing. In *EUROCRYPT* (pp. 125-157).
