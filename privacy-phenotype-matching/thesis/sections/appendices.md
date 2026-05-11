# Appendices

## Appendix A: HPO Term Structure

The Human Phenotype Ontology (HPO) is distributed as an OBO-format file in which each term carries an identifier, label, definition, and parent links. We reproduce one canonical term for reference; the methods chapter (§3.2.2) describes how ancestor and descendant sets are derived from these links.

```
[Term]
id: HP:0001250
name: Seizure
def: "Seizures are an intermittent abnormality of the central nervous
     system due to a sudden, excessive, disorderly discharge of cerebral
     neurons and characterized clinically by some combination of
     disturbance of sensation, loss of consciousness, impairment of
     psychic function, or convulsive movements." [HPO:probinson]
is_a: HP:0012638 ! Abnormal nervous system physiology
```

The HPO release used in this thesis contains 19,906 terms (data/hpo_ontology/hp.obo) organized under five top-level branches: phenotypic abnormality (~16k terms), clinical modifier, mode of inheritance, clinical course, and frequency. Information content distributions for the evaluation cohorts are reported inline with the experiments in §4.1.3 and §4.6.

---

## Appendix B: Algorithm Pseudocode

We give pseudocode for the four mechanisms that are load-bearing for the privacy claims in §3.4 and the rank-based result in §4.7. The remaining mechanisms (k-anonymity gate, rare-term filter) are simple enough that a single sentence in Methods covers them.

### B.1 Diffie–Hellman Private Set Intersection (§3.4.1)

```
Algorithm: DH-PSI Protocol
Input:  Client set A = {a_1, ..., a_n}, Server set B = {b_1, ..., b_m}
Output: A ∩ B, revealed only to Client

// Phase 1: Client encoding
1. Client samples α ← Z_q
2. For each a_i ∈ A: H_A[i] ← H(a_i)^α
3. Client → Server : {H_A[i]}

// Phase 2: Server processing
4. Server samples β ← Z_q
5. For each received H_A[i]: H_A'[i] ← H_A[i]^β
6. For each b_j ∈ B:         H_B[j]  ← H(b_j)^β
7. Server → Client : {H_A'[i]}, {H_B[j]}  (server-side shuffled)

// Phase 3: Client intersection
8. For each H_B[j]: H_B'[j] ← H_B[j]^α
9. For each a_i ∈ A: if H_A'[i] ∈ {H_B'[j]} then include a_i in intersection

// Correctness: H_A'[i] = H(a_i)^(αβ) = H(b_j)^(βα) = H_B'[j] iff a_i = b_j
```

Security: semi-honest under the Decisional Diffie–Hellman assumption on the elliptic curve (NIST P-256 in our implementation).

### B.2 Laplace Mechanism (§3.4.2)

```
Algorithm: Laplace privatisation of a similarity score
Input:  Similarity s ∈ [0, 1], budget ε > 0, sensitivity Δ (= 1 here)
Output: Noisy score ŝ, (ε,0)-DP w.r.t. one record

1. b ← Δ / ε
2. η ← Laplace(0, b)
3. Return clamp(s + η, [0, 1])
```

### B.3 Resnik Similarity with Best-Match Average (§3.3.3)

```
Algorithm: Resnik+BMA over two phenotype sets
Input:  P1, P2: term sets; ic[·]; ancestors(t) for all t
Output: BMA similarity score

Function MICA_IC(p, q):
    return max(ic[a] for a in ancestors(p) ∩ ancestors(q))

1. forward  ← mean(max(MICA_IC(p, q) for q in P2)  for p in P1)
2. backward ← mean(max(MICA_IC(p, q) for p in P1)  for q in P2)
3. Return (forward + backward) / 2
```

### B.4 Rank-Utility Exponential Mechanism for Top-k (§4.7)

```
Algorithm: Iterative exponential mechanism, rank utility
Input:  True similarities sim[1..n], budget ε, top-k size k
Output: List of k selected indices, ε-DP overall

1. order ← argsort(sim, descending)         // true ranking
2. rank[i] ← position of i in order         // sensitivity = 1
3. remaining ← {1, ..., n}; selected ← []
4. ε_round ← ε / k
5. For round = 1 to k:
     utilities ← {-rank[i] : i ∈ remaining}     // sensitivity Δu = 1
     logits    ← (ε_round / (2·Δu)) · utilities
     probs     ← softmax(logits)
     i*        ← sample(remaining, probs)
     selected.append(i*); remaining.remove(i*)
6. Return selected
```

The same skeleton with `utilities ← sim[remaining]` reproduces the score-utility exponential mechanism that exhibits the compression pathology of §4.6.1.

---

## Appendix C: Reproducing the Reported Numbers

All retrieval and privacy-utility tables in Chapter 4 are produced by two experiment scripts and the data-loading utilities listed below. The scripts use a fixed random seed and pin the cohort-construction parameters used in the main text.

| Source of finding | Script / dataset | Wall time |
|---|---|---|
| §4.1–4.4 synthetic cohort | `experiments/evaluate_hpoa.py` + `data/hpoa_evaluation/cohort_phenopackets.json` | ~30 s |
| §4.5 MI + k-anonymity attacks | `experiments/evaluate_privacy_attacks.py` | ~5 s |
| §4.6 real-cohort retrieval | `experiments/evaluate_phenopacket_store.py` (incl. Resnik+BMA baseline) | ~4 min |
| §4.7 rank-based DP | `experiments/evaluate_rank_based_dp.py` | ~30 s |

The complete pipeline runs via `make reproduce` from a clean checkout (Dockerfile bundled). Raw CSV outputs are written to `results/` and figures to `figures/`; see the project README for paths.

---

## Appendix D: Software Documentation

The reference implementation, pilot system, evaluation scripts, and Docker build instructions are documented in the project README (`README.md`) and the per-component READMEs under `app/README.md` and `scripts/`. The pilot system is hosted at <https://honors-thesis-54tubqjkgwqjm4zyegglxw.streamlit.app/>. Source code is released under an open-source license at the GitHub repository linked from the title page.

---

## Appendix E: Glossary

| Term | Definition |
|------|------------|
| **BMA** | Best Match Average; aggregation strategy for semantic similarity over term sets. |
| **DDH** | Decisional Diffie–Hellman; the cryptographic hardness assumption underlying §3.4.1. |
| **DP** | Differential Privacy; the (ε, δ) framework introduced by Dwork et al. (2006). |
| **HPO / HPOA** | Human Phenotype Ontology and its annotation corpus (`phenotype.hpoa`). |
| **IC** | Information content of a term, $-\log_2 P(t)$ over a reference corpus. |
| **k-anonymity** | Privacy model in which each released record is indistinguishable from $k-1$ others over the quasi-identifiers. |
| **MICA** | Most Informative Common Ancestor; the highest-IC term in the ancestor-set intersection of two HPO terms. |
| **MME** | Matchmaker Exchange; the federated rare-disease patient-matching standard. |
| **MRR / nDCG** | Mean Reciprocal Rank / Normalised Discounted Cumulative Gain; ranked-retrieval metrics. |
| **Phenopacket** | GA4GH v2 schema for patient phenotype, diagnosis, and provenance. |
| **PSI** | Private Set Intersection; the two-party cryptographic primitive of §3.4.1. |
| **Sensitivity** | Maximum change in a query function's output when one record is added or removed; the calibration parameter for additive-noise DP mechanisms. |
