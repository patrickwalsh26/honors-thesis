# Appendices

This section provides supplementary materials supporting the main thesis, including HPO term examples, algorithm implementations, complete experimental results, and software documentation.

---

## Appendix A: Human Phenotype Ontology Examples

### A.1 HPO Term Structure

The Human Phenotype Ontology (HPO) organizes phenotypes in a hierarchical structure using directed acyclic graphs. Each term has a unique identifier, name, definition, and relationships to other terms.

**Example HPO Term Entry:**

```
[Term]
id: HP:0001250
name: Seizure
def: "Seizures are an intermittent abnormality of the central nervous
     system due to a sudden, excessive, disorderly discharge of cerebral
     neurons and characterized clinically by some combination of
     disturbance of sensation, loss of consciousness, impairment of
     psychic function, or convulsive movements." [HPO:probinson]
synonym: "Epileptic seizure" EXACT []
synonym: "Seizures" EXACT []
is_a: HP:0012638 ! Abnormal nervous system physiology
```

### A.2 HPO Hierarchical Organization

The HPO contains approximately 18,000 terms organized under five main branches:

| Root Category | Example Terms | Term Count |
|--------------|---------------|------------|
| **Phenotypic abnormality** (HP:0000118) | Seizure, Intellectual disability | ~16,000 |
| **Clinical modifier** (HP:0012823) | Severe, Progressive | ~300 |
| **Mode of inheritance** (HP:0000005) | Autosomal dominant | ~40 |
| **Clinical course** (HP:0031797) | Childhood onset | ~60 |
| **Frequency** (HP:0040279) | Frequent, Occasional | ~15 |

### A.3 Example Phenotype Profiles

**Example 1: Marfan Syndrome (OMIM:154700)**

| HPO ID | Phenotype | Information Content |
|--------|-----------|---------------------|
| HP:0001166 | Arachnodactyly | 7.23 |
| HP:0001382 | Joint hypermobility | 4.85 |
| HP:0000501 | Glaucoma | 5.92 |
| HP:0001083 | Ectopia lentis | 9.41 |
| HP:0001065 | Striae distensae | 8.17 |
| HP:0002705 | High, narrow palate | 6.34 |
| HP:0001634 | Mitral valve prolapse | 6.89 |
| HP:0004927 | Pulmonary artery dilatation | 8.93 |
| HP:0002616 | Aortic root aneurysm | 8.56 |
| HP:0001519 | Disproportionate tall stature | 7.78 |

**Example 2: Rett Syndrome (OMIM:312750)**

| HPO ID | Phenotype | Information Content |
|--------|-----------|---------------------|
| HP:0002376 | Developmental regression | 6.89 |
| HP:0001263 | Global developmental delay | 3.21 |
| HP:0001344 | Absent speech | 7.45 |
| HP:0012171 | Stereotypical hand wringing | 11.82 |
| HP:0002540 | Inability to walk | 5.67 |
| HP:0002871 | Central apnea | 9.23 |
| HP:0002069 | Generalized tonic-clonic seizures | 6.15 |
| HP:0002650 | Scoliosis | 4.92 |
| HP:0001257 | Spasticity | 5.34 |

### A.4 Information Content Distribution

Information content (IC) values in our evaluation corpus ranged from 0.5 to 14.2 bits:

| IC Range | Description | Example Terms | Frequency |
|----------|-------------|---------------|-----------|
| 0–2 bits | Very common | Growth abnormality | 8.2% |
| 2–4 bits | Common | Seizure, Intellectual disability | 24.5% |
| 4–6 bits | Moderate | Joint hypermobility | 31.2% |
| 6–8 bits | Uncommon | Ectopia lentis | 22.8% |
| 8–10 bits | Rare | Central apnea | 10.1% |
| >10 bits | Very rare | Stereotypical hand wringing | 3.2% |

---

## Appendix B: Algorithm Pseudocode

### B.1 Diffie-Hellman Private Set Intersection

```
Algorithm: DH-PSI Protocol
Input: Client set A = {a_1, ..., a_n}, Server set B = {b_1, ..., b_m}
Output: A ∩ B (revealed only to client)

// Phase 1: Client encoding
1. Client generates random secret key α ← Z_p
2. For each a_i ∈ A:
     H_A[i] ← H(a_i)^α    // Hash to curve, raise to α
3. Client sends {H_A[1], ..., H_A[n]} to Server

// Phase 2: Server processing
4. Server generates random secret key β ← Z_p
5. For each H_A[i] received:
     H_A'[i] ← H_A[i]^β   // Raise client's values to β
6. For each b_j ∈ B:
     H_B[j] ← H(b_j)^β    // Hash and raise server's values
7. Server sends {H_A'[1], ..., H_A'[n]} and {H_B[1], ..., H_B[m]} to Client
   (shuffled randomly)

// Phase 3: Client intersection
8. For each H_B[j] received:
     H_B'[j] ← H_B[j]^α   // Raise to client's secret
9. For each a_i ∈ A:
     If H_A'[i] ∈ {H_B'[1], ..., H_B'[m]}:
         Add a_i to intersection

// Correctness: H_A'[i] = H(a_i)^(αβ) = H(b_j)^(βα) = H_B'[j] iff a_i = b_j
Return intersection
```

**Security Properties:**
- Server learns only |A|, not individual elements
- Client learns only A ∩ B, not B \ A
- Security relies on Decisional Diffie-Hellman (DDH) assumption

### B.2 Laplace Mechanism for Differential Privacy

```
Algorithm: Laplace Mechanism
Input: Query function f, Database D, Privacy parameter ε, Sensitivity Δf
Output: Differentially private answer

1. Compute true answer: y ← f(D)
2. Compute noise scale: b ← Δf / ε
3. Sample noise: η ← Laplace(0, b)
4. Return noisy answer: ỹ ← y + η

// For similarity scores (sensitivity = 1.0):
Function PrivateSimilarity(sim, ε):
    b ← 1.0 / ε
    η ← Laplace(0, b)
    Return max(0, min(1, sim + η))  // Clamp to [0, 1]
```

**Privacy Guarantee:**
For any two adjacent databases D, D' differing in one record:
$$\Pr[M(D) = y] \leq e^\epsilon \cdot \Pr[M(D') = y]$$

### B.3 Gaussian Mechanism for (ε, δ)-DP

```
Algorithm: Gaussian Mechanism
Input: Query function f, Database D, Privacy parameters (ε, δ), Sensitivity Δf
Output: (ε, δ)-differentially private answer

1. Compute true answer: y ← f(D)
2. Compute noise scale: σ ← (Δf / ε) · √(2 · ln(1.25/δ))
3. Sample noise: η ← Normal(0, σ²)
4. Return noisy answer: ỹ ← y + η
```

### B.4 Resnik Similarity (Best Match Average)

```
Algorithm: Resnik Similarity (BMA)
Input: Phenotype sets P1 = {p_1, ..., p_n}, P2 = {q_1, ..., q_m}
       IC values: ic[term] for all terms
       Ontology ancestors: anc(term) for all terms
Output: Semantic similarity score in [0, 1]

// Compute information content of Most Informative Common Ancestor
Function MICA_IC(p, q):
    common_ancestors ← anc(p) ∩ anc(q)
    return max({ic[a] : a ∈ common_ancestors})

// Compute Best Match Average
1. For each p_i ∈ P1:
     best_match[i] ← max({MICA_IC(p_i, q_j) : q_j ∈ P2})
2. For each q_j ∈ P2:
     best_match_rev[j] ← max({MICA_IC(q_j, p_i) : p_i ∈ P1})

3. forward_avg ← (1/n) · Σ best_match[i]
4. reverse_avg ← (1/m) · Σ best_match_rev[j]

5. bma ← (forward_avg + reverse_avg) / 2

// Normalize to [0, 1]
6. max_ic ← max({ic[t] : t in ontology})
7. Return bma / max_ic
```

### B.5 K-Anonymity Filtering

```
Algorithm: K-Anonymity Query Suppression
Input: Query phenotypes Q, Database patients P, Threshold k
Output: Results or SUPPRESSED

// Count patients sharing quasi-identifying phenotypes
1. candidate_patients ← {}
2. For each patient p ∈ P:
     shared_terms ← Q ∩ phenotypes(p)
     If |shared_terms| > 0:
         Add p to candidate_patients

// Suppress if cohort too small
3. If |candidate_patients| < k:
     Return SUPPRESSED  // Query reveals too much

// Otherwise, proceed with matching
4. results ← []
5. For each p ∈ candidate_patients:
     sim ← compute_similarity(Q, phenotypes(p))
     Add (p, sim) to results

6. Sort results by similarity (descending)
7. Return results
```

### B.6 Rare Term Filtering

```
Algorithm: Rare Term Filter
Input: Phenotype set P, Term frequencies freq[t], Threshold τ
Output: Filtered phenotype set

1. filtered ← {}
2. For each term t ∈ P:
     If freq[t] ≥ τ:  // Term is common enough
         Add t to filtered
     Else:
         // Optionally generalize to parent
         parent ← get_parent(t)
         If parent exists and freq[parent] ≥ τ:
             Add parent to filtered

3. Return filtered
```

---

## Appendix C: Complete Experimental Results

### C.1 Dataset Characteristics

**Evaluation Cohort Statistics:**

| Metric | Value |
|--------|-------|
| Total patients | 500 |
| Unique diseases | 100 |
| Patients per disease | 5 |
| Unique phenotypes | 1,145 |
| Mean phenotypes/patient | 10.36 |
| Std phenotypes/patient | 5.92 |
| Median phenotypes/patient | 9.0 |
| Min phenotypes/patient | 3 |
| Max phenotypes/patient | 24 |

**Disease Source Distribution:**

| Source | Disease Count | Percentage |
|--------|---------------|------------|
| OMIM | 345 | 69% |
| Orphanet | 155 | 31% |

### C.2 Baseline Similarity Metric Performance

**Complete Retrieval Metrics by Similarity Function:**

| Metric | k | Jaccard | Cosine (IC) | Cosine (unweighted) | Simplified Resnik |
|--------|---|---------|-------------|---------------------|-------------------|
| **Precision** | 1 | 0.996 | 0.998 | 0.996 | 1.000 |
| | 5 | 0.792 | 0.792 | 0.792 | 0.793 |
| | 10 | 0.400 | 0.399 | 0.400 | 0.400 |
| | 20 | 0.200 | 0.200 | 0.200 | 0.200 |
| | 50 | 0.080 | 0.080 | 0.080 | 0.080 |
| **Recall** | 1 | 0.249 | 0.250 | 0.249 | 0.250 |
| | 5 | 0.990 | 0.991 | 0.990 | 0.992 |
| | 10 | 0.999 | 0.999 | 0.999 | 1.000 |
| | 20 | 1.000 | 1.000 | 1.000 | 1.000 |
| | 50 | 1.000 | 1.000 | 1.000 | 1.000 |
| **nDCG** | 1 | 0.996 | 0.998 | 0.996 | 1.000 |
| | 5 | 0.991 | 0.992 | 0.990 | 0.993 |
| | 10 | 0.995 | 0.996 | 0.995 | 0.997 |
| | 20 | 0.996 | 0.997 | 0.996 | 0.998 |
| | 50 | 0.996 | 0.997 | 0.996 | 0.998 |
| **Hit Rate** | 1 | 0.996 | 0.998 | 0.996 | 1.000 |
| | 5 | 1.000 | 1.000 | 1.000 | 1.000 |
| | 10 | 1.000 | 1.000 | 1.000 | 1.000 |
| **MRR** | - | 0.998 | 0.999 | 0.998 | 1.000 |

### C.3 Similarity Score Distributions

**Distribution Statistics by Metric:**

| Statistic | Jaccard | Cosine (IC) | Resnik |
|-----------|---------|-------------|--------|
| Mean | 0.0101 | 0.0111 | 0.0082 |
| Std Dev | 0.0525 | 0.0627 | 0.0495 |
| Median | 0.0 | 0.0 | 0.0 |
| Min | 0.0 | 0.0 | 0.0 |
| Max | 1.0 | 1.0 | 1.0 |
| 25th %ile | 0.0 | 0.0 | 0.0 |
| 75th %ile | 0.0 | 0.0 | 0.0 |

### C.4 Differential Privacy Results

**nDCG@10 Under Varying Privacy Budgets:**

| ε Value | Mean nDCG@10 | Std Dev | Relative to Baseline |
|---------|--------------|---------|---------------------|
| ∞ (none) | 0.9975 | 0.0152 | 100.0% |
| 10.0 | 0.9891 | 0.0234 | 99.2% |
| 5.0 | 0.9766 | 0.0312 | 97.9% |
| 2.0 | 0.9423 | 0.0489 | 94.5% |
| 1.0 | 0.8927 | 0.0634 | 89.5% |
| 0.5 | 0.8234 | 0.0812 | 82.5% |
| 0.1 | 0.6123 | 0.1234 | 61.4% |

**MRR Under Varying Privacy Budgets:**

| ε Value | Mean MRR | Std Dev |
|---------|----------|---------|
| ∞ (none) | 1.000 | 0.000 |
| 10.0 | 0.982 | 0.045 |
| 5.0 | 0.956 | 0.078 |
| 2.0 | 0.912 | 0.112 |
| 1.0 | 0.834 | 0.156 |
| 0.5 | 0.723 | 0.198 |

### C.5 K-Anonymity Results

**Query Availability by K Threshold:**

| k Value | Queries Answered | Suppression Rate | Mean nDCG@10 (answered) |
|---------|------------------|------------------|-------------------------|
| 1 | 500 | 0.0% | 0.9975 |
| 2 | 500 | 0.0% | 0.9975 |
| 5 | 500 | 0.0% | 0.9975 |
| 10 | 487 | 2.6% | 0.9972 |
| 15 | 462 | 7.6% | 0.9968 |
| 20 | 445 | 11.0% | 0.9964 |

### C.6 Rare Term Filtering Results

**Performance Under Rare Term Filtering:**

| Prevalence Threshold | Terms Removed (%) | Mean nDCG@10 | MRR |
|---------------------|-------------------|--------------|-----|
| 0% (none) | 0.0% | 0.9975 | 1.000 |
| 0.1% | 2.3% | 0.9971 | 0.998 |
| 0.5% | 8.7% | 0.9952 | 0.995 |
| 1.0% | 15.2% | 0.9823 | 0.989 |
| 2.0% | 24.6% | 0.9534 | 0.961 |
| 5.0% | 42.1% | 0.8912 | 0.912 |
| 10.0% | 61.3% | 0.7845 | 0.823 |

### C.7 Combined Mechanism Performance

**Performance with Layered Privacy (ε=5.0, k=5, 1% filtering):**

| Metric | Baseline | Combined | Relative |
|--------|----------|----------|----------|
| nDCG@10 | 0.9975 | 0.9623 | 96.5% |
| MRR | 1.000 | 0.9534 | 95.3% |
| Precision@1 | 1.000 | 0.956 | 95.6% |
| Recall@10 | 1.000 | 0.998 | 99.8% |
| Hit Rate@5 | 1.000 | 0.998 | 99.8% |

### C.8 Computational Performance

**Execution Time by Operation (500 patients):**

| Operation | Mean Time (ms) | Std Dev |
|-----------|---------------|---------|
| Jaccard similarity (per pair) | 0.12 | 0.03 |
| Cosine similarity (per pair) | 0.18 | 0.05 |
| Resnik similarity (per pair) | 2.34 | 0.45 |
| PSI protocol (per pair) | 12.45 | 1.23 |
| DP noise addition (per score) | 0.02 | 0.01 |
| K-anonymity check | 1.56 | 0.34 |
| Full query (100 candidates) | 287.3 | 45.2 |

**Memory Usage:**

| Component | Memory (MB) |
|-----------|-------------|
| HPO ontology (loaded) | 45.2 |
| IC values (11,514 terms) | 0.4 |
| 500 patient profiles | 2.1 |
| Similarity matrix (500×500) | 2.0 |

---

## Appendix D: Software Documentation

### D.1 Repository Structure

```
privacy-phenotype-matching/
├── src/
│   ├── similarity/
│   │   ├── __init__.py
│   │   ├── hpo_similarity.py      # Similarity metrics
│   │   ├── graph_similarity.py    # Graph-based measures
│   │   └── lsh_ann.py             # MinHash LSH indexing
│   ├── privacy/
│   │   ├── __init__.py
│   │   ├── psi.py                 # Private Set Intersection
│   │   ├── differential_privacy.py # DP mechanisms
│   │   ├── k_anonymity.py         # K-anonymity filtering
│   │   └── privacy_calculator.py  # Privacy accounting
│   ├── evaluation/
│   │   ├── __init__.py
│   │   ├── metrics.py             # Retrieval metrics
│   │   ├── privacy_utility.py     # Tradeoff analysis
│   │   └── leakage_audit.py       # Privacy auditing
│   ├── data_integration/
│   │   ├── __init__.py
│   │   ├── hpoa_patient_generator.py  # Synthetic patients
│   │   ├── beacon_v2.py           # Beacon adapter
│   │   ├── matchmaker_exchange.py # MME adapter
│   │   └── external_sources.py    # Data loaders
│   └── data_generation/
│       ├── __init__.py
│       └── synthetic_phenopackets.py  # Phenopacket generation
├── data/
│   ├── hpo/                       # HPO ontology files
│   └── hpo_annotations/           # HPOA annotations
├── results/
│   └── hpoa/                      # Evaluation results
├── thesis/
│   ├── sections/                  # Thesis markdown files
│   └── compile_thesis.py          # Compilation script
└── requirements.txt
```

### D.2 Installation

```bash
# Clone repository
git clone https://github.com/username/privacy-phenotype-matching.git
cd privacy-phenotype-matching

# Create virtual environment
python3 -m venv venv
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Download HPO data
python -c "from src.data_integration.external_sources import download_hpo_data; download_hpo_data()"
```

### D.3 Core API Reference

#### Similarity Module

```python
from src.similarity.hpo_similarity import (
    JaccardSimilarity,
    CosineSimilarity,
    ResnikSimilarity
)

# Initialize with IC values
ic_values = load_ic_values("data/hpo/ic_values.json")
resnik = ResnikSimilarity(ic_values=ic_values)

# Compute similarity between phenotype sets
patient1_terms = ["HP:0001250", "HP:0001263", "HP:0002376"]
patient2_terms = ["HP:0001250", "HP:0001344", "HP:0002540"]

similarity = resnik.compute_pairwise_similarity(patient1_terms, patient2_terms)
print(f"Resnik similarity: {similarity:.4f}")
```

#### Privacy Module

```python
from src.privacy.psi import DiffieHellmanPSI
from src.privacy.differential_privacy import LaplaceMechanism
from src.privacy.k_anonymity import KAnonymityFilter

# Private Set Intersection
psi = DiffieHellmanPSI()
client_set = {"HP:0001250", "HP:0001263"}
server_set = {"HP:0001250", "HP:0002376", "HP:0001344"}

intersection = psi.compute_intersection(client_set, server_set)
print(f"Intersection: {intersection}")

# Differential Privacy
dp = LaplaceMechanism(epsilon=5.0)
private_score = dp.privatize_similarity(0.85)
print(f"Privatized similarity: {private_score:.4f}")

# K-Anonymity
k_anon = KAnonymityFilter(k=5)
should_suppress = k_anon.check_suppression(query_phenotypes, database)
```

#### Evaluation Module

```python
from src.evaluation.metrics import RetrievalMetrics

metrics = RetrievalMetrics()

# Compute retrieval performance
results = metrics.evaluate(
    similarities=similarity_matrix,
    ground_truth=disease_labels,
    k_values=[1, 5, 10]
)

print(f"nDCG@10: {results['ndcg'][10]:.4f}")
print(f"MRR: {results['mrr']:.4f}")
```

### D.4 Phenopacket Format

Our system uses GA4GH Phenopackets v2.0 for patient representation:

```json
{
  "id": "patient-001",
  "subject": {
    "id": "patient-001",
    "timeAtLastEncounter": {
      "age": {
        "iso8601duration": "P25Y"
      }
    },
    "sex": "FEMALE"
  },
  "phenotypicFeatures": [
    {
      "type": {
        "id": "HP:0001250",
        "label": "Seizure"
      },
      "excluded": false,
      "onset": {
        "age": {
          "iso8601duration": "P5Y"
        }
      }
    },
    {
      "type": {
        "id": "HP:0001263",
        "label": "Global developmental delay"
      },
      "excluded": false
    }
  ],
  "diseases": [
    {
      "term": {
        "id": "OMIM:312750",
        "label": "Rett Syndrome"
      }
    }
  ],
  "metaData": {
    "created": "2026-04-20T12:00:00Z",
    "createdBy": "privacy-phenotype-matching",
    "phenopacketSchemaVersion": "2.0"
  }
}
```

### D.5 Configuration Options

```python
# Privacy configuration
PRIVACY_CONFIG = {
    "psi": {
        "enabled": True,
        "curve": "SECP256R1",
        "security_parameter": 128
    },
    "differential_privacy": {
        "enabled": True,
        "epsilon": 5.0,
        "delta": 1e-5,
        "mechanism": "laplace"  # or "gaussian"
    },
    "k_anonymity": {
        "enabled": True,
        "k": 5,
        "suppress_on_failure": True
    },
    "rare_term_filtering": {
        "enabled": True,
        "prevalence_threshold": 0.01,
        "generalize_to_parent": True
    }
}

# Similarity configuration
SIMILARITY_CONFIG = {
    "metric": "resnik",  # "jaccard", "cosine", "resnik"
    "ic_source": "corpus",  # "corpus" or "intrinsic"
    "aggregation": "bma"  # "bma" or "max"
}
```

### D.6 Running Evaluations

```bash
# Generate synthetic patients from HPOA
python -m src.data_integration.hpoa_patient_generator \
    --n_patients 500 \
    --n_diseases 100 \
    --output data/synthetic_patients.json

# Run baseline evaluation
python -m src.evaluation.privacy_utility \
    --input data/synthetic_patients.json \
    --output results/evaluation_results.json

# Run privacy-utility tradeoff analysis
python -m src.evaluation.privacy_utility \
    --epsilon_values 0.1,0.5,1.0,2.0,5.0,10.0 \
    --k_values 1,2,5,10,15,20 \
    --output results/privacy_tradeoff.json
```

### D.7 License and Citation

This software is released under the MIT License.

**Citation:**
```bibtex
@mastersthesis{walsh2026privacy,
  title={Privacy-Preserving Phenotype Matching for
         Rare Disease Cohort Discovery},
  author={Walsh, Patrick},
  year={2026},
  school={Stanford University},
  advisor={Montgomery, Stephen B.}
}
```

---

## Appendix E: Glossary

| Term | Definition |
|------|------------|
| **BMA** | Best Match Average; aggregation strategy for semantic similarity |
| **DDH** | Decisional Diffie-Hellman; cryptographic hardness assumption |
| **DP** | Differential Privacy; formal privacy framework with quantifiable guarantees |
| **HPO** | Human Phenotype Ontology; standardized vocabulary for clinical phenotypes |
| **HPOA** | HPO Annotations; database linking diseases to phenotype terms |
| **IC** | Information Content; measure of term specificity (-log probability) |
| **k-anonymity** | Privacy model ensuring each record is indistinguishable from k-1 others |
| **LSH** | Locality-Sensitive Hashing; approximate nearest neighbor technique |
| **MICA** | Most Informative Common Ancestor; lowest common ancestor with highest IC |
| **MME** | Matchmaker Exchange; federated network for rare disease patient matching |
| **MRR** | Mean Reciprocal Rank; retrieval metric emphasizing first relevant result |
| **nDCG** | Normalized Discounted Cumulative Gain; ranked retrieval metric |
| **Phenopacket** | GA4GH standard format for patient phenotype data exchange |
| **PSI** | Private Set Intersection; cryptographic protocol for set comparison |
| **Resnik** | IC-based semantic similarity using most informative common ancestor |
| **Sensitivity** | Maximum change in query output from adding/removing one record |

---

## Appendix F: Supplementary Figures

### F.1 System Architecture Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         Privacy-Preserving Pipeline                      │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                          │
│  ┌──────────────┐     ┌──────────────┐     ┌──────────────────────────┐ │
│  │   Patient    │     │  Rare Term   │     │   K-Anonymity           │ │
│  │  Phenopacket │────▶│   Filtering  │────▶│   Validation            │ │
│  │   (Input)    │     │              │     │                         │ │
│  └──────────────┘     └──────────────┘     └───────────┬──────────────┘ │
│                                                         │                │
│                                           ┌─────────────▼──────────────┐ │
│                                           │   Private Set             │ │
│                                           │   Intersection (PSI)      │ │
│                                           │                           │ │
│                                           │   Phenotype Intersection  │ │
│                                           │   Without Full Disclosure │ │
│                                           └─────────────┬──────────────┘ │
│                                                         │                │
│  ┌──────────────┐     ┌──────────────┐     ┌───────────▼──────────────┐ │
│  │   Private    │◀────│  Differential│◀────│   Similarity            │ │
│  │   Results    │     │   Privacy    │     │   Computation           │ │
│  │   (Output)   │     │   (DP Noise) │     │   (Resnik/Cosine)       │ │
│  └──────────────┘     └──────────────┘     └──────────────────────────┘ │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

### F.2 HPO Ontology Excerpt

```
                        HP:0000118
                   Phenotypic abnormality
                           │
           ┌───────────────┼───────────────┐
           │               │               │
    HP:0000707      HP:0000152      HP:0000924
Abnormality of    Head and neck  Abnormality of
nervous system    abnormality    skeletal system
           │               │               │
    HP:0012638      HP:0000234      HP:0000929
  Abnormal NS      Abnormality    Abnormal skull
  physiology       of head        morphology
           │               │
    HP:0001250      HP:0000252
     Seizure       Microcephaly
           │
    HP:0002069
  Generalized
tonic-clonic seizure
```

### F.3 Privacy-Utility Tradeoff Curve

```
nDCG@10
  │
1.0├───●────●────●                    ● Baseline (no privacy)
   │        ╲                         ● ε = 10
0.9├──────────●──────●                ● ε = 5
   │              ╲                   ● ε = 2
0.8├───────────────────●              ● ε = 1
   │                     ╲
0.7├────────────────────────●         ● ε = 0.5
   │
0.6├─────────────────────────────●    ● ε = 0.1
   │
   └──┬──────┬──────┬──────┬──────┬─────▶ Privacy (1/ε)
     0.1    0.2    0.5    1.0    2.0

Recommended Operating Range: ε ∈ [2.0, 5.0]
   - High utility (>93% nDCG@10)
   - Meaningful privacy guarantees
```

---

*End of Appendices*
