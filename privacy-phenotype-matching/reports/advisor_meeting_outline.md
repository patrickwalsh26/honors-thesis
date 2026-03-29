# Advisor Meeting Outline - Fall Quarter Progress

**Patrick Walsh | Privacy-Preserving Phenotype Matching**
**Meeting Date**: [To be scheduled]
**Duration**: 30 minutes

---

## Agenda

1. **Project Recap** (3 min)
2. **Fall Quarter Deliverables** (10 min)
3. **Live Demo** (7 min)
4. **Results & Evaluation** (5 min)
5. **Winter Quarter Plan** (3 min)
6. **Discussion & Feedback** (2 min)

---

## 1. Project Recap (3 min)

### Goal
Build privacy-preserving phenotype matching system for rare disease cohort discovery

### Key Innovation
Two-step reveal ladder: Private similarity test → Minimal data disclosure

### Privacy Mechanisms (Planned)
- Private Set Intersection (PSI)
- Differential Privacy (DP)
- k-anonymity guardrails
- Rare term suppression

---

## 2. Fall Quarter Deliverables (10 min)

### ✅ Infrastructure
- Professional Python codebase (~1,500 LOC)
- GA4GH Phenopackets v2.0 compliant
- HPO ontology integration (15K+ terms)
- Full documentation (README, quickstart, 15-page report)

### ✅ Synthetic Data Generation
**4 Rare Disease Profiles:**
- Marfan syndrome (1 in 5,000)
- Ehlers-Danlos syndrome (1 in 10,000)
- Achondroplasia (1 in 25,000)
- Progeria (1 in 2.5M)

**Features:**
- Realistic phenotypic variation (3-15 HPO terms/patient)
- Core, common, rare feature tiers
- Noise and negation modeling
- 200-patient test cohort generated

### ✅ Baseline Similarity Metrics
**Three approaches implemented:**
1. **Jaccard**: Simple set overlap
2. **Cosine (IC-weighted)**: Vector similarity with information content
3. **Simplified Resnik**: Ontology-aware (exact matches)

**Information Content:**
- Computed empirically from corpus
- 49 unique HPO terms
- IC range: 0.92 (common) to 4.61 (rare)

### ✅ Evaluation Framework
- Top-k patient retrieval
- Recall@k metrics (k=5, 10, 20)
- Similarity matrix computation
- Disease-stratified evaluation

---

## 3. Live Demo (7 min)

### Demo Script: `examples/demo_baseline.py`

**Show:**
1. **Data loading** - 200 patients, 4 diseases
2. **IC computation** - Term informativeness
3. **Similarity comparison** - Same vs different disease
4. **Retrieval** - Query Marfan patient, get top-10 matches
5. **Evaluation** - Recall@20 = 0.39 across diseases

**Expected output:**
```
Query: Marfan syndrome patient
Top-10 matches: All Marfan patients (perfect precision)
Similarity scores: 0.77-0.97
```

---

## 4. Results & Evaluation (5 min)

### Key Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Same-disease similarity | 0.66 | Strong within-disease coherence |
| Cross-disease similarity | 0.00 | Clear separation (synthetic) |
| Recall@20 | 0.39 | 39% relevant in top-20 |
| Top-10 precision | 1.00 | Perfect for demo query |

### Analysis

**Strengths:**
- IC weighting improves over basic Jaccard
- Clear disease separation with synthetic data
- Baseline provides targets for privacy-utility tradeoff

**Limitations:**
- Synthetic data has perfect labels (real data will be harder)
- Simplified Resnik (no ancestor graph yet)
- Small term vocabulary (49 terms)

**Next Steps:**
- Full Resnik with HPO hierarchy
- Larger, noisier synthetic cohorts
- Real data validation (pending STARR/MIMIC access)

---

## 5. Winter Quarter Plan (3 min)

### Phase 1: Privacy Mechanisms (Weeks 1-4)
**Priority: PSI Implementation**
- Private Set Intersection for secure overlap
- Differential privacy noise on counts
- k-anonymity filtering
- Rare term suppression

**Deliverable:** Privacy-utility curves (retrieval quality vs ε)

### Phase 2: Advanced Retrieval (Weeks 5-7)
- LSH-ANN approximate nearest neighbors
- Full Resnik with ancestor graph
- Compare PSI vs LSH-ANN performance

**Deliverable:** Protocol comparison report

### Phase 3: Integration (Weeks 8-10)
- MME adapter (Matchmaker Exchange)
- Beacon v2 integration
- Docker containerization
- Preliminary thesis chapter

**Deliverable:** Working containerized service + preliminary report

---

## 6. Discussion & Feedback (2 min)

### Questions for Advisor

1. **Privacy parameters**: What ε/k values are acceptable for clinical use?

2. **Real data access**:
   - Status of STARR-OMOP-deid application?
   - Should we pursue MIMIC-IV access?
   - Alternative datasets?

3. **Winter priorities**:
   - PSI vs LSH-ANN: Which first?
   - Should we implement both?
   - Focus on depth (one protocol well) vs breadth (both protocols)?

4. **Publication targets**:
   - Aim for conference (ACM BCB, AMIA, PSB)?
   - Or journal (npj Digital Medicine, JBI)?
   - Timeline for submission?

5. **Collaboration**:
   - Interest from other Montgomery Lab members?
   - Potential connections to clinical groups?
   - Real use case for testing?

6. **Scope management**:
   - Current scope appropriate for honors thesis?
   - Any features to add/remove?

---

## Materials Available

📄 **Full Report**: `reports/fall_2024_progress_report.md` (15 pages)
📄 **Executive Summary**: `reports/executive_summary.md` (3 pages)
💻 **Demo Script**: `examples/demo_baseline.py`
📊 **Code Repository**: 8 commits, clean history
📦 **Datasets**: 200 synthetic phenopackets

**Git Repository:**
```bash
git clone [repository URL]
cd privacy-phenotype-matching
python3 examples/demo_baseline.py
```

---

## Backup Slides (If Time Permits)

### Technical Deep Dive: Resnik Similarity

Resnik similarity uses information content of most informative common ancestor (MICA):

```
IC(term) = -log(P(term))
Resnik(t1, t2) = IC(MICA(t1, t2))
BMA(A, B) = avg of best matches for each term
```

**Example:**
- HP:0001519 (Dolichostenomelia): IC = 2.3
- HP:0001166 (Arachnodactyly): IC = 2.5
- Both → Marfan core features

### Privacy-Utility Tradeoff Theory

```
Privacy ↑ (larger ε, k) → Utility ↓
Privacy ↓ (smaller ε, k) → Utility ↑
```

**Goal:** Find optimal balance where:
- Privacy guarantees are clinically acceptable
- Retrieval quality enables useful cohort discovery

**Planned experiments:**
- Vary ε: [0.1, 0.5, 1.0, 2.0, 5.0]
- Vary k: [3, 5, 10, 20]
- Measure: Recall@20 vs privacy parameters

### Related Work Comparison

| System | Privacy | Matching | Standards |
|--------|---------|----------|-----------|
| **SHEPHERD** | None | Phenotype-guided | No |
| **MatchMaker Exchange** | Manual | Curated | Yes (MME) |
| **Beacon v2** | DP counts | Boolean | Yes (Beacon) |
| **Our system** | PSI+DP+k-anon | Automatic | Yes (GA4GH) |

**Key distinction:** Only system with multi-layered privacy for automated phenotype matching

---

## Action Items (Post-Meeting)

- [ ] Incorporate advisor feedback
- [ ] Prioritize Winter Quarter tasks
- [ ] Schedule check-in meetings (bi-weekly?)
- [ ] Update timeline based on data access status
- [ ] Begin PSI implementation (Week 1)

---

**Contact**: walshp26@stanford.edu
**Repository**: [Add after creating GitHub repo]
**Next Meeting**: [To be scheduled - 6 weeks?]
