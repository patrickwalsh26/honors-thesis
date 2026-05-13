# Results

This chapter reports six experiments. Sections 4.1–4.4 evaluate retrieval and privacy-utility tradeoffs on a synthetic cohort sampled from the HPOA corpus. Section 4.5 measures empirical privacy via membership-inference and singling-out attacks against the mechanisms of §3.4. Section 4.6 evaluates retrieval on the Phenopacket Store real-patient cohort and exposes a 20–50× gap between synthetic and real DP budgets. Section 4.7 closes that gap with the rank-utility exponential mechanism.

## 4.1 Dataset Characteristics

### 4.1.1 HPO Annotations Corpus

Our evaluation leverages the HPO phenotype annotation file (phenotype.hpoa), which contains curated disease-phenotype associations from three authoritative sources. Table 1 summarizes the corpus composition.

**Table 1: HPO Annotations Corpus Summary**

| Source | Diseases | Phenotype Annotations | Mean Phenotypes/Disease |
|--------|----------|----------------------|------------------------|
| OMIM | 8,592 | 174,394 | 20.3 |
| Orphanet | 4,335 | 88,051 | 20.3 |
| DECIPHER | 47 | 954 | 20.3 |
| **Total** | **12,974** | **263,399** | **20.3** |

The corpus encompasses 11,514 unique HPO terms, representing comprehensive coverage of phenotypic abnormalities across rare disease categories. OMIM contributes the largest share of diseases (66.2%), reflecting its role as the primary reference for Mendelian disorders. Orphanet provides substantial coverage of European rare disease classifications (33.4%), while DECIPHER contributes a smaller set of chromosomal abnormality phenotypes (0.4%).

### 4.1.2 Evaluation Cohort

We generated a synthetic patient cohort by sampling from the HPO annotations corpus. Our generation procedure (described in Section 3.6.2) produces patients whose phenotype profiles realistically reflect the variability observed in clinical practice, including incomplete phenotyping and incidental findings.

**Table 2: Evaluation Cohort Characteristics**

| Parameter | Value |
|-----------|-------|
| Total patients | 500 |
| Diseases represented | 100 |
| Patients per disease | 5 |
| Disease source distribution | 69% OMIM, 31% Orphanet |
| Unique phenotypes in cohort | 1,145 |
| Phenotypes per patient (mean ± SD) | 10.4 ± 5.9 |
| Phenotypes per patient (median) | 9.0 |
| Phenotypes per patient (range) | 3–24 |
| Phenotype recall from disease profile | 75% |
| Noise phenotype rate | 10% |

The balanced design assigns exactly 5 patients to each of 100 sampled diseases, establishing clear ground truth for retrieval evaluation: for any query patient, the 4 other patients with the same underlying disease constitute the relevant set. This design enables precise measurement of retrieval metrics while controlling for disease prevalence effects.

The phenotype distribution (mean: 10.4, range: 3–24) aligns with clinical observations. Studies of rare disease cohorts report similar ranges, with patients typically presenting 5–15 observable phenotypic features at initial evaluation (Köhler et al., 2019). The 75% recall rate models incomplete phenotyping (not all disease-associated phenotypes are documented for every patient), while the 10% noise rate captures incidental or co-morbid findings unrelated to the primary diagnosis.

### 4.1.3 Information Content Distribution

We computed corpus-based information content (IC) for each HPO term using the evaluation cohort. Figure 1 (not shown) presents the IC distribution. The most informative terms (highest IC) are rare, highly specific phenotypes appearing in few patients:

**Table 3: Example High-IC Phenotypes**

| HPO Term | Label | IC (bits) | Patients |
|----------|-------|-----------|----------|
| HP:0000639 | Nystagmus | 7.2 | 3 |
| HP:0001274 | Agenesis of corpus callosum | 6.8 | 5 |
| HP:0002353 | EEG abnormality | 6.5 | 7 |

Conversely, common phenotypes have low IC:

**Table 4: Example Low-IC Phenotypes**

| HPO Term | Label | IC (bits) | Patients |
|----------|-------|-----------|----------|
| HP:0001263 | Global developmental delay | 2.1 | 89 |
| HP:0001252 | Hypotonia | 2.4 | 72 |
| HP:0001249 | Intellectual disability | 2.5 | 68 |

This IC distribution directly informs our similarity metrics: rare phenotypes contribute more heavily to similarity scores, appropriately reflecting their greater diagnostic specificity.

## 4.2 Baseline Retrieval Performance

We evaluated four similarity metrics on the retrieval task: given a query patient, rank all database patients by similarity and assess whether the 4 patients with the same underlying disease appear at the top of the ranking.

### 4.2.1 Overall Performance

Table 5 presents retrieval performance across metrics and evaluation cutoffs.

**Table 5: Baseline Retrieval Performance**

| Metric | P@1 | P@5 | R@5 | R@10 | nDCG@10 | Hit@5 | MRR |
|--------|-----|-----|-----|------|---------|-------|-----|
| Jaccard | 0.996 | 0.792 | 0.990 | 0.999 | 0.995 | 1.000 | 0.998 |
| Cosine (IC) | 0.998 | 0.792 | 0.991 | 0.999 | 0.996 | 1.000 | 0.999 |
| Cosine (unweighted) | 0.996 | 0.792 | 0.990 | 0.999 | 0.995 | 1.000 | 0.998 |
| **Resnik (simplified)** | **1.000** | **0.793** | **0.992** | **1.000** | **0.997** | **1.000** | **1.000** |

*P@k = Precision at k; R@k = Recall at k; nDCG = normalized Discounted Cumulative Gain; Hit@k = Hit Rate at k; MRR = Mean Reciprocal Rank. Best values in bold.*

All four metrics achieve near-saturated performance on the synthetic cohort: P@1 ≥ 0.996, R@10 ≥ 0.999, Hit@5 = 1.000 across all configurations. Simplified Resnik is best on every metric (P@1 = MRR = 1.000), confirming that IC weighting helps but the improvement is small because the synthetic cohort's clean disease-template sampling makes the task easy regardless of metric (this saturation is a feature of the cohort, not the system; §4.6 reports a markedly harder real-cohort regime where metric choice matters more). The 0.2–0.4% of queries that exhibit rank inversions involve patients with phenotypic overlap arising from phenotypically related diseases, a known limitation that the privacy mechanisms below do not introduce.

### 4.2.2 Precision-Recall Tradeoff

Figure 2 (not shown) plots precision versus recall at varying k values. The characteristic precision-recall tradeoff emerges: precision decreases as k increases (since fewer of the expanded result set are relevant), while recall increases (more relevant patients are retrieved).

**Table 6: Precision and Recall at Varying k (Resnik Similarity)**

| k | Precision@k | Recall@k |
|---|-------------|----------|
| 1 | 1.000 | 0.250 |
| 5 | 0.793 | 0.992 |
| 10 | 0.400 | 1.000 |
| 20 | 0.200 | 1.000 |
| 50 | 0.080 | 1.000 |

With 4 relevant patients per query, the theoretical maximum P@5 is 0.80 (achieved when all 4 relevant patients appear in the top 5, plus one irrelevant patient). Our observed P@5 = 0.793 approaches this bound, indicating that the top 5 results almost always consist of the 4 relevant patients plus at most one false positive.

## 4.3 Similarity Score Distributions

Understanding the distribution of similarity scores is essential for setting appropriate thresholds and interpreting results. Table 7 presents distributional statistics.

**Table 7: Similarity Score Distributions (All Patient Pairs)**

| Metric | Mean | SD | Median | Q25 | Q75 | Max |
|--------|------|----|--------|-----|-----|-----|
| Jaccard | 0.010 | 0.053 | 0.000 | 0.000 | 0.000 | 1.000 |
| Cosine (IC) | 0.011 | 0.063 | 0.000 | 0.000 | 0.000 | 1.000 |
| Resnik | 0.008 | 0.050 | 0.000 | 0.000 | 0.000 | 1.000 |

Several observations emerge:

**High sparsity.** The median similarity is 0.000 for all metrics, indicating that most patient pairs share no phenotypes. This reflects the diverse disease composition of the cohort: patients with different underlying diseases typically have non-overlapping phenotype sets.

**Heavy-tailed distribution.** While the mean similarity is low (0.008–0.011), the maximum reaches 1.000, representing identical phenotype profiles (patients with the same disease and similar phenotype sampling). The large gap between mean and maximum indicates a heavy-tailed distribution.

**Privacy implications.** The sparsity of similarity scores has important privacy implications. Most queries return predominantly zero-similarity results, with only a small cluster of similar patients. This pattern could itself be informative to an adversary, motivating our privacy mechanisms that protect result distributions.

### 4.3.1 Intra-Disease vs. Inter-Disease Similarity

To further characterize similarity structure, we partition patient pairs into intra-disease (same underlying disease) and inter-disease (different diseases) groups.

**Table 8: Similarity by Disease Relationship (Cosine IC)**

| Relationship | Mean | SD | Min | Max |
|--------------|------|----|-----|-----|
| Intra-disease | 0.547 | 0.218 | 0.112 | 1.000 |
| Inter-disease | 0.008 | 0.042 | 0.000 | 0.683 |

Intra-disease pairs exhibit substantially higher similarity (mean: 0.547) than inter-disease pairs (mean: 0.008), confirming that phenotype similarity effectively discriminates disease-sharing patients. The minimal overlap in distributions (inter-disease max: 0.683 vs. intra-disease min: 0.112) explains the excellent retrieval performance.

## 4.4 Privacy-Utility Tradeoff Analysis

We evaluated the impact of privacy mechanisms on retrieval performance. Each mechanism introduces a different form of protection with corresponding utility implications.

### 4.4.1 Differential Privacy Impact

We applied the Laplace mechanism to similarity scores with varying privacy parameter ε and measured retrieval performance degradation.

**Table 9: Retrieval Performance Under Differential Privacy**

| ε | nDCG@10 | MRR | P@5 | Δ nDCG@10 |
|---|---------|-----|-----|-----------|
| ∞ (no DP) | 0.997 | 1.000 | 0.793 | baseline |
| 10.0 | 0.989 | 0.994 | 0.784 | −0.8% |
| 5.0 | 0.976 | 0.981 | 0.769 | −2.1% |
| 2.0 | 0.941 | 0.952 | 0.738 | −5.6% |
| 1.0 | 0.892 | 0.908 | 0.691 | −10.5% |
| 0.5 | 0.821 | 0.843 | 0.628 | −17.7% |
| 0.1 | 0.634 | 0.671 | 0.482 | −36.4% |

*Δ nDCG@10 = relative change from non-private baseline.*

The results reveal a smooth privacy-utility tradeoff:

**Moderate privacy (ε = 2–5) preserves utility.** At ε = 5.0, nDCG@10 decreases by only 2.1% while providing meaningful privacy guarantees. At ε = 2.0, the 5.6% degradation remains acceptable for many applications.

**Strong privacy (ε ≤ 1) significantly impacts utility.** At ε = 1.0, widely considered a reasonable privacy threshold, nDCG@10 drops by 10.5%. At ε = 0.5, utility degradation reaches 17.7%.

**Very strong privacy (ε = 0.1) substantially degrades utility.** With ε = 0.1, nDCG@10 falls to 0.634, a 36.4% reduction. While this provides strong theoretical guarantees, the practical utility may be insufficient for clinical applications.

### 4.4.2 K-Anonymity Impact

K-anonymity enforcement suppresses results when the matching cohort contains fewer than k patients, protecting rare phenotype combinations from identification.

**Table 10: K-Anonymity Impact on Result Availability**

| k | Suppression Rate | Available Queries | P@5 (available) |
|---|------------------|-------------------|-----------------|
| 2 | 0.0% | 500/500 | 0.793 |
| 5 | 0.0% | 500/500 | 0.793 |
| 10 | 2.4% | 488/500 | 0.791 |
| 20 | 8.6% | 457/500 | 0.787 |
| 50 | 24.2% | 379/500 | 0.772 |

With k ≤ 5, no queries are suppressed because all diseases have at least 5 patients in our balanced cohort. At k = 10, 2.4% of queries return empty results due to insufficient cohort size. This rate increases with k, reaching 24.2% at k = 50.

For available (non-suppressed) queries, retrieval precision remains stable (P@5: 0.772–0.793), indicating that k-anonymity enforcement does not degrade result quality; it merely restricts result availability for rare cases.

### 4.4.3 Rare Term Filtering Impact

Rare term filtering removes phenotypes appearing in fewer than a threshold number of patients, reducing quasi-identifier risk.

**Table 11: Rare Term Filtering Impact**

| Prevalence Threshold | Terms Filtered | Mean Phenotypes Remaining | nDCG@10 |
|---------------------|----------------|---------------------------|---------|
| 0 (none) | 0 | 10.4 | 0.997 |
| 1% (5 patients) | 412 | 8.7 | 0.984 |
| 2% (10 patients) | 687 | 7.2 | 0.961 |
| 5% (25 patients) | 891 | 5.1 | 0.912 |
| 10% (50 patients) | 982 | 3.4 | 0.847 |

Aggressive filtering (10% threshold) removes the majority of unique phenotypes, substantially reducing discriminative power (nDCG@10: 0.847). Conservative filtering (1% threshold) removes only extremely rare terms while preserving 98.4% of baseline utility.

### 4.4.4 Combined Privacy Mechanisms

We evaluated the composition of multiple privacy mechanisms: the configuration most likely to be deployed in practice.

**Table 12: Combined Privacy Mechanism Performance**

| Configuration | nDCG@10 | MRR | Suppression |
|---------------|---------|-----|-------------|
| Baseline (no privacy) | 0.997 | 1.000 | 0% |
| ε = 5.0 only | 0.976 | 0.981 | 0% |
| k = 5 only | 0.997 | 1.000 | 0% |
| Rare filter (1%) only | 0.984 | 0.991 | 0% |
| ε = 5.0 + k = 5 | 0.976 | 0.981 | 0% |
| ε = 5.0 + k = 5 + filter (1%) | 0.962 | 0.971 | 0% |
| ε = 2.0 + k = 5 + filter (1%) | 0.928 | 0.942 | 0% |
| ε = 1.0 + k = 5 + filter (1%) | 0.873 | 0.894 | 0% |

The composed configuration with ε = 5.0, k = 5, and a 1% rare-term filter retains 96.5% of baseline nDCG@10. The synthetic-cohort utility cost of layered protection is therefore modest under Laplace DP, a conclusion that §4.6 will overturn for real patients.

### 4.4.5 Privacy-Utility Pareto Frontier

Figure 7 plots the privacy-utility Pareto frontier for the synthetic cohort. Three operating points span the practical range: ε = 5.0 / k = 5 / 1% filter retains 96% of nDCG@10; ε = 2.0 / k = 5 / 1% filter retains 93%; ε = 1.0 / k = 5 / 2% filter retains 85%. These numbers, generalised from synthetic data, motivated earlier published privacy-utility claims in this space; §4.6 reports the substantially harsher real-cohort behaviour that revises the recommendation.

## 4.5 Empirical Privacy: Adversarial Attack Evaluation

The §4.4 Pareto analysis reports privacy *parameters*; it does not show that those parameters defeat actual adversaries. We close this gap by instantiating the attacks from the threat model (§3.1.2) and measuring their success on our system: a membership-inference attack against the DP score-release oracle (testing invariant **I2**) and a singling-out attack against the k-anonymity gate (testing **I3**).

### 4.5.1 Membership-Inference Attack vs. Differential Privacy

Following the membership-inference (MI) literature, we evaluate two attackers against the IC-weighted cosine retrieval oracle wrapped in the Laplace mechanism:

- A **threshold attacker** (Yeom et al., 2018) submits the target patient as a query, observes the max returned similarity, and decides "member" if the score exceeds a threshold. The AUC of the resulting ROC measures the attacker's discriminative power without committing to a threshold.
- A **shadow-model attacker** (Shokri et al., 2017) trains a Random Forest classifier on features extracted from query-response distributions (max, mean, std, gap, percentile quantiles) over three independently sampled shadow databases (250 patients each, 80 in/80 out per shadow). The trained classifier is evaluated on a held-out shadow split.

We sweep $\varepsilon \in \{0.1, 0.5, 1.0, 2.0, 5.0, 10.0, +\infty\}$, where $\varepsilon = +\infty$ denotes the non-private baseline. The 500-patient HPOA cohort is split 50/50 into members and non-members.

**Table 13: Membership-inference attack ROC AUC vs. DP budget.** Higher AUC = stronger attack. AUC = 0.5 is random guessing.

| ε | Threshold attack AUC | Shadow-model attack AUC | Defender advantage* |
|---|---------------------|--------------------------|---------------------|
| ∞ (no DP) | 0.967 | 0.976 | 0.000 |
| 10 | 0.882 | 0.909 | 0.075 |
| 5 | 0.578 | 0.617 | 0.389 |
| 2 | 0.500 | 0.519 | 0.476 |
| 1 | 0.500 | 0.500 | 0.481 |
| 0.5 | 0.500 | 0.500 | 0.481 |
| 0.1 | 0.500 | 0.500 | 0.481 |

\*Defender advantage = (AUC$_{\varepsilon=\infty}$ − AUC$_{\varepsilon}$), averaged over both attackers.

![Membership-inference attack ROC AUC vs.\ Laplace privacy budget on the synthetic HPOA cohort. Both attackers collapse to random guessing at $\varepsilon \leq 2.0$.](figures/fig10_mia_pareto.pdf){#fig:mia width=85%}

Figure 10 plots the privacy-utility frontier. Three regimes are apparent: (i) **practical privacy** at $\varepsilon \leq 2.0$, where both attackers are reduced to random guessing despite full transcript access; (ii) a **transition band** at $2.0 < \varepsilon < 10.0$, where the shadow-model attacker recovers signal faster than the simple threshold attacker; and (iii) a **vulnerable** regime at $\varepsilon \geq 10.0$, where AUC exceeds 0.88, so DP at this budget should be regarded as non-protective in this deployment. The shadow-model attacker uniformly dominates the threshold attacker in the transition regime, illustrating that simple-score thresholds underestimate the adversary capability available to a real attacker.

Crucially, the regime where DP successfully defends against MI ($\varepsilon \leq 2.0$) overlaps with the regime where retrieval utility remains usable (nDCG@10 = 0.928 at $\varepsilon = 2.0$ per §4.4.1). This is the central operating-point claim of the thesis: $\varepsilon \approx 1$–$2$ simultaneously delivers strong empirical MI resistance and acceptable retrieval performance.

### 4.5.2 Singling-Out Attack and k-Anonymity Ablation

The DP analysis above protects against membership inference but does not by itself prevent singling-out via rare phenotype quasi-identifiers. We evaluate the k-anonymity gate of §3.4.3 against the rare-term singling-out adversary defined in §3.1.2.

**Threat realization.** The adversary holds a single rare HPO term $t$ known to apply to a target patient (e.g., inferred from a publication or registry). They query the database for all patients exhibiting $t$. Without protection, the system returns the cohort $S_t = \{P : t \in \Phi(P)\}$, and uniform-random selection inside $S_t$ identifies the target with probability $1/|S_t|$. The k-anonymity gate suppresses the response (returns $\bot$) when $|S_t| < k$.

We enumerated every HPO term in the cohort with cohort prevalence $1 \leq |S_t| \leq 50$ (1,142 terms; median $|S_t| = 4$, max $|S_t| = 46$), and measured the expected re-identification probability $\Pr_{\text{re-id}} = \mathbb{E}_t [\mathbb{1}\{|S_t| \geq k\} / |S_t|]$ across that distribution.

**Table 14: k-anonymity ablation against rare-term singling-out (n = 1,142 rare terms).**

| k | Suppression rate | Re-id probability (post-gate) | Singling-out rate (post-gate) | Mean released cohort size |
|---|-----------------|-------------------------------|-------------------------------|---------------------------|
| 1 (no gate) | 0.000 | 0.419 | 0.224 | 4.4 |
| 2 | 0.224 | 0.195 | 0.000 | 5.4 |
| 5 | 0.678 | 0.048 | 0.000 | 8.3 |
| 10 | 0.926 | 0.005 | 0.000 | 15.5 |
| 20 | 0.988 | 0.0005 | 0.000 | 26.9 |
| 50 | 1.000 | 0.000 | 0.000 | n/a |

![k-anonymity ablation on rare HPO-term singling-out. (a) Suppression rate vs.\ k. (b) Expected re-identification probability post-gate vs.\ k; horizontal dashed line is the unprotected baseline.](figures/fig11_kanon_ablation.pdf){#fig:kanon width=92%}

Figure 11 visualizes both axes of the tradeoff: suppression rate climbs steeply between $k = 2$ and $k = 10$ (22% → 93%), while re-identification probability collapses by nearly three orders of magnitude (0.419 → 0.005). Crucially, even at $k = 2$ the singling-out rate drops to 0: every rare term that originally identified a unique patient is now blocked. The persistent re-id probability at $k = 2$ (19.5%) comes entirely from queries that *are* released ($|S_t| \geq 2$), where the attacker still benefits from small cohort sizes; this is the gap that **rare-term filtering** (§4.4.3) closes by suppressing the query before it ever reaches the gate.

**Composition.** Combining $\varepsilon = 1.0$ DP (§4.5.1) with $k = 5$ k-anonymity reduces both attack surfaces simultaneously: MI attack AUC = 0.500 and rare-term re-identification probability ≤ 0.048, while retaining nDCG@10 = 0.870 (§4.4.4 composed configuration). This is the recommended deployment configuration for production use.

## 4.6 External Validation: Real Published Patient Cohort

The retrieval results in §4.2–4.4 use a synthetic cohort sampled from HPOA disease profiles. To check whether those numbers reflect real clinical phenotyping (known to be noisier, more variable, and less complete than disease-template sampling), we re-ran the evaluation on the **Monarch Phenopacket Store** (Danis et al., 2025), a curated collection of GA4GH phenopackets derived from published case reports with confirmed OMIM diagnoses and PMID provenance.

After filtering for patients with ≥3 observed phenotypes from diseases with ≥2 documented patients, the cohort contains 8,343 real patients across 586 OMIM diseases. For tractable all-vs-all retrieval we constructed a balanced subsample of the 100 most-populated diseases capped at 15 patients each, yielding 1,500 patients for the experiments below.

**Table 15: Retrieval performance on real published cohort (Phenopacket Store).** Best published rare-disease semantic-similarity systems (Phenomizer, LIRICAL) typically report MRR in the 0.7–0.9 range on similar tasks (Köhler et al., 2009; Robinson et al., 2020), placing our system within the established literature.

| System | MRR | nDCG@10 | Recall@5 | P@1 |
|--------|-----|---------|----------|-----|
| Cosine-IC (this work, non-private) | **0.869** | **0.689** | **0.738** | **0.821** |
| Resnik+BMA (Phenomizer-style, non-private) | 0.828 | 0.642 | 0.689 | 0.783 |

Two findings are notable: (i) absolute retrieval performance drops by roughly 30% versus the synthetic cohort (MRR 1.00 → 0.87, nDCG 0.99 → 0.69), confirming that synthetic-cohort numbers materially over-estimate real-world utility; (ii) Cosine-IC outperforms full Resnik+BMA over the HPO DAG, suggesting that the additional ontology-traversal cost of Resnik buys little above IC-weighted vector matching when the IC priors are estimated from a sufficiently large corpus.

### 4.6.1 Privacy-Utility Tradeoff on the Real Cohort

Re-running the Laplace ε sweep on the same 1,500-patient cohort reveals a markedly steeper privacy-utility curve than the synthetic-cohort analysis of §4.4.1 suggested.

**Table 16: Cosine-IC retrieval under the Laplace mechanism (Phenopacket Store cohort).**

| ε | MRR | nDCG@10 | Recall@5 | P@1 |
|---|-----|---------|----------|-----|
| ∞ (no DP) | 0.869 | 0.689 | 0.738 | 0.821 |
| 50 | 0.861 | 0.677 | 0.726 | 0.814 |
| 20 | 0.786 | 0.582 | 0.629 | 0.723 |
| 10 | 0.548 | 0.340 | 0.374 | 0.440 |
| 5 | 0.194 | 0.090 | 0.092 | 0.087 |
| 2 | 0.055 | 0.020 | 0.021 | 0.022 |
| 1 | 0.034 | 0.012 | 0.012 | 0.011 |
| 0.5 | 0.028 | 0.010 | 0.010 | 0.010 |

![Retrieval performance on the Monarch Phenopacket Store cohort (1,500 real published-case-report patients across 100 OMIM diseases). Cosine-IC nDCG@10 and MRR fall sharply for $\varepsilon < 20$; the dotted line marks the Resnik+BMA Phenomizer-style non-private baseline.](figures/fig12_phenopacket_benchmark.pdf){#fig:ppstore width=85%}

Figure 12 plots the resulting curve. On the synthetic cohort (§4.4.1), ε = 1.0 retained ≥85% of baseline nDCG; on the real cohort, only ε ≥ 50 retains the equivalent fraction. The gap arises because real-patient similarity-score distributions are compressed: between-disease gaps are smaller and more overlapping than in disease-template-sampled patients, so the Laplace noise scale of 1/ε needed to dominate the signal is substantially larger.

This finding has direct deployment implications. Privacy claims based on synthetic-cohort experiments are systematically optimistic: a deployed system targeting the empirical MI-defended regime of ε ≤ 1 (Table 13) would deliver retrieval performance close to random on real patients with this metric. The next subsection (§4.7) shows that the rank-based exponential mechanism resolves this gap empirically.

## 4.7 Rank-Based Differential Privacy on the Real Cohort

The Laplace results above (§4.6.1) reveal a pathology specific to *score-based* DP: when same-disease similarity-score gaps are compressed into a 0.2–0.3 range, Lap(1/ε) noise dominates the rank-determining signal. Two principled responses exist in the DP literature, both ε-DP for top-k release: the exponential mechanism (McSherry & Talwar, 2007) operating on a utility function, and Report-Noisy-Max (Dwork & Roth, 2014). For exponential-mechanism utility we evaluate two natural choices (the raw similarity score, and the candidate's *rank* in the true ordering) and compare against per-score Laplace at matched ε.

The rank utility is principled here because the rank function has sensitivity 1 under record addition/removal (adding or removing a single patient shifts any other patient's rank by at most 1), while the rank *gaps* are O(n) by construction (rank 1 vs. rank 100 differ by 99) regardless of how compressed the underlying similarities are. Rank-based selection decouples the noise scale from the score magnitude that score-based mechanisms inherit.

For top-k selection we use the standard iterative composition: ε is split into k rounds of ε/k each, sampling without replacement, yielding ε-DP overall by the composition theorem. The same composition applies to the score-based and rank-based variants.

**Table 17: Three DP top-10 release mechanisms on the Phenopacket Store cohort (1,500 patients, 100 OMIM diseases).** All mechanisms are ε-DP under the same privacy accounting.

| ε | Laplace (per-score) | Exp. mech. (score) | **Exp. mech. (rank)** |
|---|--------------------|--------------------|-----------------------|
|     | MRR / nDCG@10     | MRR / nDCG@10      | MRR / nDCG@10           |
| ∞ (no DP) | 0.869 / 0.689 | 0.869 / 0.689 | 0.869 / 0.689 |
| 50 | 0.862 / 0.678 | 0.067 / 0.024 | **0.865 / 0.688** |
| 20 | 0.786 / 0.582 | 0.039 / 0.013 | **0.853 / 0.682** |
| 10 | 0.551 / 0.342 | 0.037 / 0.012 | **0.834 / 0.667** |
| 5 | 0.186 / 0.087 | 0.029 / 0.010 | **0.782 / 0.620** |
| 2 | 0.054 / 0.020 | 0.033 / 0.011 | **0.665 / 0.466** |
| 1 | 0.033 / 0.012 | 0.029 / 0.010 | **0.544 / 0.333** |
| 0.5 | 0.032 / 0.011 | 0.034 / 0.011 | **0.396 / 0.199** |

![Three ε-DP top-k release mechanisms on the Phenopacket Store cohort. Rank-based exponential mechanism (green) recovers near-baseline retrieval at ε ≥ 5; per-score Laplace (blue) requires ε ≥ 50 for the same retention; score-based exponential mechanism (orange) collapses because compressed similarity scores yield indistinguishable utilities.](figures/fig14_rank_based_dp.pdf){#fig:rankdp width=95%}

Figure 14 visualizes the gap. Three findings:

- **Rank-based exponential mechanism recovers 90% of baseline nDCG@10 at ε = 5** (0.620 vs. 0.689), and 96% at ε = 10. Per-score Laplace requires ε ≥ 50 for comparable retention, a 10× budget efficiency advantage at matched ε-DP.
- **Score-based exponential mechanism collapses across the entire sweep.** This confirms the §5.1.5 hypothesis directly: the pathology is the score-utility-magnitude pairing, not the mechanism class. Replacing Laplace with the exponential mechanism without changing the utility function provides no relief.
- **The MI-defended regime (ε ≤ 1) is now usable.** At ε = 1, where Table 13 reports shadow-model MI attack AUC = 0.50 (random), rank-based selection delivers MRR = 0.544, a 16× improvement over Laplace MRR = 0.033 at the same ε. The deployment recommendation revised in §5.3.1 reflects this.

## 4.8 Computational Performance

We measured computational overhead for the privacy-preserving pipeline components.

**Table 18: Computational Overhead**

| Operation | Time (ms) | Notes |
|-----------|-----------|-------|
| Similarity matrix (500×500) | 1,847 | All pairwise similarities |
| Per-patient similarity | 3.7 | Single query vs. database |
| PSI per patient pair | 12.4 | EC operations, NIST P-256 |
| DP noise addition | 0.002 | Per similarity score |
| K-anonymity check | 0.1 | Per result set |
| IC computation (corpus) | 1,205 | One-time preprocessing |

*Measured on Apple M1 Pro, 16GB RAM, Python 3.10.*

The dominant cost is similarity computation (3.7 ms per query patient), which scales linearly with database size. PSI adds substantial overhead (12.4 ms per pair) due to elliptic curve operations, but this can be parallelized and is acceptable for federated matching where queries are infrequent. DP and k-anonymity checks impose negligible overhead.

For a database of 10,000 patients, per-query latency would be approximately 40 seconds for full similarity computation, or 2 minutes with PSI. Approximate methods (locality-sensitive hashing) could reduce this to sub-second latency at the cost of retrieval accuracy.

## 4.9 Summary of Key Findings

Our experimental evaluation yields the following principal findings:

1. **Phenotype-based retrieval is highly effective.** Baseline metrics achieve nDCG@10 > 0.99 and MRR = 1.0, demonstrating that phenotype similarity reliably identifies patients with the same underlying disease.

2. **IC weighting improves performance.** Information content weighting, which emphasizes rare phenotypes, provides consistent improvement over unweighted approaches (Resnik: MRR = 1.000 vs. Jaccard: MRR = 0.998).

3. **On the synthetic cohort, privacy is achievable with moderate utility cost.** Combined privacy mechanisms (ε = 5.0, k = 5, 1% filter) preserve 96.5% of baseline utility while providing layered protection; this is a synthetic-cohort result that should be read alongside findings 7–8 below.

4. **Differential privacy empirically thwarts membership inference at deployable budgets.** Shadow-model attack ROC AUC drops from 0.976 (no DP) to 0.500 (random guessing) at ε = 1.0, while nDCG@10 remains above 0.85. The threshold attack and shadow attack agree on this transition.

5. **k-anonymity collapses singling-out by orders of magnitude.** Expected re-identification probability against the rare-term adversary falls from 0.419 (no gate) to 0.005 at k = 10, with all unique-patient queries suppressed even at k = 2.

6. **Real-cohort validation places our system in the published literature's range.** On 1,500 real published-case-report patients from 100 OMIM diseases (Phenopacket Store), Cosine-IC retrieval achieves MRR = 0.87 and nDCG@10 = 0.69, within the 0.7–0.9 MRR band typical of Phenomizer/LIRICAL-class systems. Resnik+BMA over the full HPO DAG is essentially tied with the simpler Cosine-IC.

7. **Synthetic cohorts overestimate the safe Laplace-DP budget by 1–2 orders of magnitude.** With per-score Laplace noise, ε = 1 preserves >85% of synthetic-cohort nDCG but only ~2% of real-cohort nDCG; ε ≥ 20 is needed for comparable retention. Practical deployment cannot rely on synthetic-cohort privacy claims for score-based mechanisms.

8. **Rank-based exponential mechanism closes the gap.** Replacing per-score Laplace with iterative exponential mechanism on rank-utility recovers 90% of non-private nDCG@10 at ε = 5 on the real cohort (vs. 13% for Laplace at matched ε) under the same ε-DP guarantee, a 10× budget efficiency advantage. At ε = 1, where empirical MI attack AUC falls to 0.50 (random), rank-based retrieval achieves MRR = 0.544 (16× over Laplace).

9. **Tradeoffs are smooth and configurable.** The privacy-utility frontier allows institutions to select operating points matching their risk tolerance and regulatory requirements.

10. **Computational overhead is practical.** Per-query latency remains under 1 second for databases of 500 patients, with clear paths to scaling via approximation or parallelization.

---

