# Future Work

Four directions follow directly from the results of this thesis, ordered by expected impact-per-effort.

## 6.1 Variance-Reduced and Composed Rank-Based Mechanisms

Section §4.7 evaluated the iterative exponential mechanism with rank utility and demonstrated 90% recovery of non-private nDCG@10 at ε = 5. Three refinements are immediately tractable on the same cohort. **Permute-and-Flip** (McKenna & Sheldon, 2020) is a variance-reduced alternative with the same ε-DP guarantee and tighter low-ε accuracy bounds; **one-shot top-k via GAP-K** (Bassily et al., 2021) avoids the ε/k budget split that iterative composition incurs; and **joint composition with the k-anonymity gate** of §4.5.2 has only been characterised on the synthetic cohort in §4.4.4. Each is expected to shift the rank-utility curve further toward ε ≤ 1, the empirically MI-defended regime.

## 6.2 EHR-Derived and Cross-Institutional Validation

The Phenopacket Store cohort (§4.6) is drawn from published case reports — better-phenotyped than typical EHRs. §5.1.5 predicts that the synthetic-to-real DP gap will widen further on EHR data; measuring this requires a cross-institutional retrospective on consented cohorts from the Undiagnosed Diseases Network (UDN), the GREGoR Consortium, or comparable rare-disease programmes. A **prospective clinical trial** comparing standard-of-care diagnosis with matching-augmented diagnosis on undiagnosed-disease patients would establish clinical utility; this is a multi-year effort but is the path to demonstrable outcome improvement. A **federation pilot** across 3–5 institutions running the framework as MME nodes would surface deployment issues (data-use agreements, IRB coordination, network failure modes) that the single-cohort evaluation cannot.

## 6.3 Cryptographic Strengthening

Our PSI implementation is semi-honest secure (§3.4.1). **Malicious-secure PSI** based on committed OT (Pinkas et al., 2018) or cut-and-choose (Lindell, 2017) raises the security model at a 2–5× computational cost. **Zero-knowledge proofs** of correct execution (Bünz et al., 2018) detect server-side cheating without revealing inputs. **Multi-party PSI** (Kolesnikov et al., 2017) extends the framework to coalition queries that the current two-party protocol does not support. **Fully homomorphic encryption** (Cheon et al., 2017) would permit similarity computation entirely on ciphertexts; today's CKKS/BGV implementations are within an order of magnitude of practicality for IC-weighted cosine.

## 6.4 Multi-Modal Matching and Workflow Integration

Rare-disease diagnosis increasingly couples phenotype with genomic, transcriptomic, and imaging signals. Extending the framework to **gene–phenotype joint matching**, **RNA-seq outlier integration** (Frésard et al., 2019), and **imaging phenotypes** introduces modality-specific privacy considerations that the current PSI/DP/k-anonymity composition does not address by default. Practical deployment additionally requires **FHIR integration** for EHR-resident phenotype extraction, a **richer pilot UI** (the current Streamlit prototype is a single-user demo, §3.8), and **institutional privacy-budget accounting** beyond the per-session tracker in the pilot. **Adversarial red-team evaluation** by an external team — going beyond our MI / singling-out attacks — should precede any production rollout.

## 6.5 GA4GH Standards Extensions

Three concrete proposals would let the broader rare-disease ecosystem benefit from the mechanisms here. (1) **Phenopacket privacy metadata** indicating the mechanism parameters used to generate or transmit a record. (2) **Beacon v3 privacy modes** standardising DP and PSI query types alongside the existing boolean / count / record granularities. (3) **MME privacy handshake** for negotiating privacy parameters between query and database nodes during the existing match protocol. Engagement with the relevant GA4GH working groups would advance these proposals through the standards process.
