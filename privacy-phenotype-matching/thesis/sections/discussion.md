# Discussion

The retrieval, attack, and rank-DP experiments of Chapter 4 jointly support a single deployment-level conclusion: privacy-preserving phenotype matching is technically feasible at ε-DP budgets that also defeat membership-inference attacks, provided the DP mechanism is calibrated to rank rather than score magnitude. Sections 5.1 interpret the supporting evidence in order. Section 5.2 contrasts the framework with prior federated systems; §5.3 gives concrete deployment recommendations; §5.4–5.5 catalogue limitations and ethical considerations.

## 5.1 Interpretation of Results

### 5.1.1 Effectiveness of Phenotype-Based Retrieval

The synthetic-cohort baseline achieves nDCG@10 = 0.997 (Resnik, IC-weighted), validating the premise of phenotype-based matching: patients with the same disease consistently share enough HPO terms to be retrieved together. On the real published-cohort benchmark (§4.6), Cosine-IC achieves MRR = 0.87 and nDCG@10 = 0.69 — within the 0.7–0.9 MRR band reported by Phenomizer and LIRICAL on comparable tasks (Köhler et al., 2009; Robinson et al., 2020). The synthetic-cohort numbers reflect what is achievable when phenotype profiles match disease templates faithfully; the real-cohort numbers reflect the genuine noise of published clinical phenotyping (atypical presentations, incomplete documentation, and inter-curator variability).

The consistent advantage of IC-weighted metrics over unweighted alternatives confirms prior work (Köhler et al., 2009; Pesquita et al., 2009). More surprising is that Cosine-IC marginally outperforms full Resnik+BMA over the HPO DAG on the real cohort (MRR 0.87 vs. 0.83): when corpus-IC priors capture enough of the term-discrimination signal, ontology-traversal cost buys little additional accuracy. This argues for the simpler metric in deployment.

### 5.1.2 Sparsity and Its Implications

Similarity-score distributions are sharply bimodal: median = 0, mean ≈ 0.01, with a heavy tail of high-similarity same-disease pairs. The result helps retrieval (true matches separate cleanly from background) and complicates privacy (the distinctive shape of a query's response — predominantly zero with a few high values — leaks information about the query patient's cohort size). The three privacy mechanisms address this jointly: DP noise smooths the response distribution, k-anonymity suppresses responses with too few non-zero entries, and rare-term filtering removes the quasi-identifiers that drive the heaviest tails before computation begins.

### 5.1.3 Privacy-Utility Tradeoffs

The tradeoffs separate cleanly along three axes — the DP mechanism, the k-anonymity gate, and the rare-term filter — and each has a distinct cohort-dependence pattern.

**DP under Laplace** degrades gracefully on the synthetic cohort (2.1% nDCG loss at ε = 5, 10.5% at ε = 1) but catastrophically on the real cohort (98% nDCG loss at ε = 1, requiring ε ≥ 20 for 80% retention). The cohort sensitivity is not a property of DP itself — the (ε, 0)-DP guarantee is unchanged across cohorts; only the utility cost varies. §5.1.5 explains why and §4.7 measures the rank-utility fix.

**Membership-inference defence is independent of cohort and robust at deployable ε.** Table 13 shows shadow-model attack AUC collapsing from 0.98 (no DP) to 0.50 (random) at ε ≤ 1 — the regime where the DP guarantee is empirically meaningful against the strongest attacker we evaluate. The MI defence is real, quantifiable, and the same on both cohorts.

**k-anonymity is comparatively cheap.** Re-identification probability against the rare-term singling-out adversary falls three orders of magnitude (0.42 → 0.005) between k = 1 and k = 10, and every unique-patient query is blocked already at k = 2. Composing k-anonymity with moderate DP gives defence in depth without the utility collapse of low-ε per-score DP alone.

**Rare-term filtering exhibits the sharpest internal tradeoff.** A 1% prevalence threshold retains 98% of utility; a 10% threshold removes the diagnostic phenotypes that make matching work in the first place (nDCG@10 drops to 0.847). Rare phenotypes are simultaneously the most identifying and the most diagnostic; aggressive filtering wastes the signal we are trying to preserve.

### 5.1.4 Mechanism Composition

No single mechanism is sufficient. PSI protects against server-side phenotype enumeration; DP bounds output leakage; k-anonymity prevents singling-out; rare-term filtering removes quasi-identifiers before computation. The §4.4.4 composed configuration (ε = 5, k = 5, 1% filter) preserves 96.5% of synthetic-cohort utility, but the synthetic-to-real gap (§5.1.5) means deployment configurations should re-tune ε upward on real-population data.

### 5.1.5 The Synthetic-to-Real Generalization Gap

The most consequential finding of our real-cohort evaluation is that the safe DP budget is 20–50× larger on real published patients than synthetic-cohort experiments suggest. The mechanism is straightforward: synthetic patients are sampled from a single disease's phenotype profile with 75% recall and 10% noise, producing well-separated similarity-score distributions where same-disease pairs score in the 0.6–0.9 range and different-disease pairs near 0. Real patients exhibit substantial within-disease phenotypic heterogeneity, co-morbidities, and atypical presentations, compressing same-disease scores into the 0.2–0.5 range with substantial overlap against different-disease pairs. A Laplace noise scale of 1/ε that is dwarfed by a 0.6 same-disease score becomes comparable to a 0.3 one, collapsing the rank signal.

This has two implications for the field. First, **published privacy-utility evaluations of rare-disease matching are systematically optimistic** to the extent they rely on disease-profile-sampled cohorts. Reviewers should expect deployment-time ε ≥ 10 with per-score Laplace even where synthetic-cohort experiments support ε ≤ 1. Second, **the principled fix is to replace the score utility with a rank utility under the exponential mechanism**. Results §4.7 evaluates this directly and confirms the claim empirically: rank-based exponential mechanism recovers 90% of non-private nDCG@10 at ε = 5 (vs. 13% for Laplace) and a usable MRR = 0.544 at ε = 1 (vs. 0.033 for Laplace), with the same ε-DP guarantee. Score-based exponential mechanism, by contrast, suffers the identical compression pathology as Laplace and offers no improvement — confirming that the mechanism *family* matters less than the *utility-function sensitivity-to-signal ratio*. Section 5.3.1 revises our deployment recommendation accordingly.

## 5.2 Comparison to Related Work

**Beacons.** The Shringarpure–Bustamante (2015) Beacon attack exploits Boolean variant-presence responses; Raisaro et al. (2017) propose DP-Beacons but report substantial utility loss. Our framework operates at the phenotype level, which (i) admits ontology-driven generalisation that genomic positions do not, (ii) composes three mechanisms rather than relying on DP alone, and (iii) restricts the released transcript to similarity scores or top-$k$ identifiers rather than the per-variant Boolean signature that drives the Shringarpure–Bustamante attack.

**Matchmaker Exchange.** MME nodes currently exchange phenotype data in cleartext between participating institutions. Our framework is API-compatible with MME but interposes PSI on the wire, enabling matching by institutions whose data-sharing restrictions preclude cleartext exchange today.

**Privacy-preserving genomics more broadly.** Prior work has emphasised variant-level protection (Ayday et al., 2013; Chen et al., 2019); phenotype privacy has received less attention despite the known re-identification potential of rare phenotype combinations (El Emam et al., 2011). Our rank-utility mechanism extends naturally to other tasks with compressed similarity-score distributions — for instance, RNA-seq outlier-based diagnosis (Frésard et al., 2019), where the protection target is multi-omic rather than purely phenotypic.

## 5.3 Practical Deployment Considerations

### 5.3.1 Parameter Selection

Our deployment recommendations differ from earlier drafts of this work because the real-cohort evaluation revealed the synthetic-to-real gap of §5.1.5. The values below are calibrated against the Phenopacket Store benchmark (§4.6), not the synthetic cohort.

| Parameter | Recommended Value | Rationale |
|-----------|-------------------|-----------|
| DP mechanism | **Iterative exponential mechanism, rank utility** (§4.7) | Noise tracks rank gaps (always O(n)) rather than compressed score magnitudes. 10× more budget-efficient than Laplace at matched ε-DP on real patients. |
| ε (per query, rank-exp) | **2–5** | Retains 77–90% of non-private nDCG@10 on the Phenopacket Store cohort. At ε = 1, still delivers MRR = 0.544 — usable for triage if MI-defended privacy is the priority. |
| ε (per query, Laplace fallback) | ≥ 20 | If a rank-utility implementation is unavailable, Laplace per-score requires ε in this range for comparable retention. Not recommended for new deployments. |
| k (anonymity) | 5–10 | Re-identification probability ≤ 0.05 at k = 10 against the rare-term singling-out adversary |
| Rare-term threshold | 1% | Preserves ~98% utility while suppressing extreme outliers |
| Similarity metric | Cosine-IC | Marginally outperforms Resnik+BMA on the real cohort at lower computational cost |

The empirical MI defense (§4.5.1) is unchanged: the DP guarantee still bounds attack AUC, regardless of ε's utility cost. Institutions choosing Laplace at ε = 20 obtain the same theoretical (ε, 0)-DP guarantee they would obtain at ε = 1 — they simply accept a weaker formal bound in exchange for usable retrieval. Whether ε = 20 constitutes "meaningful" privacy is a deployment-level judgment that compositions across many queries (§3.4.4 accountant) make tighter than the per-query value suggests.

### 5.3.2 Privacy Budget Management

For ongoing query systems, privacy budget management is essential. Each query consumes a portion of the cumulative ε budget. We recommend:

1. **Per-patient budgets**: Allocate a lifetime privacy budget (e.g., ε = 10) per patient, tracking consumption across queries.
2. **Query rate limiting**: Restrict query frequency to slow budget consumption and detect abuse.
3. **Budget renewal**: Consider periodic budget renewal (e.g., annually) for patients with updated phenotypes.
4. **Opt-out mechanisms**: Allow patients to withdraw, invalidating future queries against their data.

Our implementation includes a privacy accountant that tracks cumulative ε expenditure and can enforce budget limits.

### 5.3.3 Integration Pathways

We envision three deployment models:

**Institutional deployment.** A single institution deploys the system for internal matching. Privacy mechanisms protect against insider threats and audit compliance. This model requires minimal infrastructure change.

**Federated deployment.** Multiple institutions participate via MME-compatible interfaces. Each institution runs a privacy-preserving node that responds to queries without revealing raw phenotypes. This model requires coordinated deployment but enables cross-institutional matching.

**Centralized trusted service.** A trusted third party operates the matching service, with institutions submitting encrypted phenotypes. This model simplifies deployment but requires trust in the central operator.

Each model involves different trust assumptions and architectural requirements. Our modular design supports all three.

### 5.3.4 Regulatory Alignment

Privacy-preserving phenotype matching aligns with key regulatory frameworks:

**HIPAA.** The Privacy Rule permits de-identified data sharing. Our k-anonymity and rare term filtering support de-identification requirements. DP provides additional protection exceeding minimum standards.

**GDPR.** The right to data protection is balanced against research exemptions. Privacy-preserving computation demonstrates "privacy by design" and "data minimization" principles. Consent frameworks should disclose matching participation and privacy mechanisms.

**Common Rule.** For federally funded research in the US, IRB approval is required for human subjects research. Privacy-preserving matching may qualify for expedited review when risks are minimal, though IRB practices vary.

We recommend institutional consultation with compliance officers, as specific requirements depend on data sources, patient populations, and use cases.

## 5.4 Limitations

### 5.4.1 Evaluation Limitations

Our evaluation has several limitations that contextualize the results:

**Real-cohort scale.** The Phenopacket Store benchmark (§4.6) uses 1,500 patients across 100 OMIM diseases — a balanced subsample of the 8,343 filtered patients to keep all-vs-all retrieval tractable. The DP curve and the Resnik+BMA baseline would benefit from being re-run on the full filtered corpus; scaling experiments are out of scope for this submission.

**Curated case-report selection.** Phenopacket Store patients are extracted from peer-reviewed case reports and are consequently better-phenotyped than typical clinical records. EHR-derived cohorts (which we do not access in this work) would likely exhibit even greater within-disease heterogeneity and further widen the synthetic-to-real gap of §5.1.5.

**Single-institution simulation.** Queries run within a single corpus rather than across federated institutions with heterogeneous annotation practices. Cross-institutional matching faces additional challenges from differing phenotyping conventions, HPO version mismatches, and variable annotation depth.

**Privacy-attack scope.** Our membership-inference experiment (§4.5.1) implements the Yeom-threshold and Shokri-shadow attacks; stronger attackers (label-only MI, gradient-leakage, query-based reconstruction) are not evaluated. The k-anonymity ablation considers a single-rare-term adversary; multi-term quasi-identifiers may permit residual leakage at the recommended k = 5–10.

### 5.4.2 Technical Limitations

**Semi-honest adversary model.** Our PSI implementation assumes semi-honest (honest-but-curious) adversaries who follow the protocol but attempt to learn from observations. Malicious adversaries who deviate from the protocol are not addressed. Malicious-secure PSI protocols exist but incur substantial overhead.

**Computational overhead.** PSI adds ~12ms per patient pair, making real-time queries against large databases challenging. For a database of 100,000 patients, naive PSI would require ~20 minutes per query. Approximate methods and indexing could reduce this but merit separate evaluation.

**Phenotype-only matching.** We focus on phenotype matching without genomic data. Clinical matching often integrates both modalities. Extending privacy protection to combined phenotype-genotype matching requires additional mechanisms.

### 5.4.3 Scope Limitations

**Clinical validation.** We have not validated clinical utility with actual rare disease diagnoses. Whether privacy-preserving matching leads to diagnoses that would not otherwise occur requires prospective clinical studies.

**User interface.** Our implementation provides programmatic APIs but not clinician-facing interfaces. Effective deployment requires UX design for clinical workflows.

**Longitudinal phenotyping.** Phenotypes evolve as patients age and diseases progress. Our static matching does not address temporal phenotype dynamics.

## 5.5 Ethical Considerations

Privacy-preserving computation reduces but does not eliminate identification risk; residual risk depends on adversary capabilities that may evolve. Three deployment principles follow from our measured results and apply to any production rollout.

**Transparency in consent.** Informed-consent documents must state which mechanisms are used, what each protects against (e.g., membership inference at ε ≤ 1 reduces shadow-model AUC to random, §4.5.1), and what they explicitly do not protect against (out-of-scope threats, §3.1.2). Privacy claims should be calibrated to the real-cohort budget regime (§4.6, §4.7), not the more optimistic synthetic-cohort numbers historically reported in the literature.

**Patient agency.** Patients should be able to opt in or out at any time, with matching against their record disabled on withdrawal; receive optional notification when their record contributes to a successful match; and benefit from discoveries enabled by their data through equitable governance arrangements (Ramoni et al., 2017). The per-session privacy accountant we expose in the pilot (§3.8) is a foundation for the institutional-budget machinery these guarantees require.

**Equitable adoption.** Complex privacy mechanisms deployed only at well-resourced institutions risk creating a two-tier system in which underrepresented populations contribute disproportionately less to the rare-disease cohort. Open-source release, low-friction reproducibility (`make reproduce`, §3.7), and turn-key cloud deployment (§3.8) lower the barrier; institutional matching networks should explicitly fund participation by lower-resource sites. Secondary use beyond matching — commercial profiling, surveillance — must be excluded by governance policy rather than technical controls alone.

## 5.6 Summary

Privacy-preserving phenotype matching is technically feasible, but the privacy-utility tradeoff is sharper on real patients than synthetic-cohort evaluations imply. We make four contributions:

1. **Empirically validated retrieval on real published patients.** Non-private Cosine-IC achieves MRR = 0.87 / nDCG@10 = 0.69 on 1,500 Phenopacket Store patients across 100 OMIM diseases, placing the system within the Phenomizer/LIRICAL band.

2. **Empirically measured privacy defense.** Shadow-model MI attack AUC collapses from 0.98 (no DP) to 0.50 (random) at ε ≤ 1; k-anonymity at k = 10 cuts re-identification probability from 0.42 to 0.005. These numbers validate threat-model invariants I2 and I3 (§3.1.2).

3. **The synthetic-to-real privacy budget gap, and the rank-utility fix.** Per-score Laplace DP needs ε that is 20–50× larger on real cohorts than synthetic experiments suggest, because real similarity-score distributions are compressed. Replacing the score utility with a rank utility under the iterative exponential mechanism recovers 90% of non-private nDCG@10 at ε = 5 (vs. 13% for Laplace), with the same ε-DP guarantee — closing the gap empirically.

4. **A revised deployment configuration** that operationalizes the fix (§5.3.1): Cosine-IC similarity, rank-utility exponential mechanism with ε ∈ [2, 5], k ∈ [5, 10], 1% rare-term filtering.

---

