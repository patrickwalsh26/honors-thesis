#!/usr/bin/env python3
"""
Demonstration of DECIPHER data integration for phenotype matching.

This script shows how to:
1. Load/generate DECIPHER patient data
2. Convert to GA4GH Phenopacket format
3. Compute phenotype similarities
4. Find similar patients
5. Evaluate matching performance

This serves as a template for using real DECIPHER data once access is obtained.

Usage:
    python examples/demo_decipher.py
"""

import sys
from pathlib import Path
from collections import defaultdict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_integration.decipher_loader import (
    DECIPHERLoader,
    DECIPHERSimulator,
    DECIPHERPatient,
)
from src.similarity.hpo_similarity import (
    JaccardSimilarity,
    CosineSimilarity,
    SimplifiedResnikSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
)


def print_header(text):
    """Print a formatted header."""
    print("\n" + "=" * 70)
    print(f"  {text}")
    print("=" * 70)


def print_patient_summary(patient: DECIPHERPatient):
    """Print a summary of a DECIPHER patient."""
    print(f"  ID: {patient.patient_id}")
    print(f"  Sex: {patient.sex}")
    print(f"  Phenotypes: {len(patient.phenotypes)} total")
    print(f"  Sample HPO terms: {', '.join(patient.phenotypes[:5])}...")
    if patient.genes:
        print(f"  Genes: {', '.join(patient.genes[:3])}...")


def print_phenopacket_summary(pp: dict):
    """Print a summary of a Phenopacket."""
    print(f"  ID: {pp['id']}")
    print(f"  Sex: {pp['subject']['sex']}")
    print(f"  Features: {len(pp['phenotypicFeatures'])} phenotypes")
    features = [f["type"]["id"] for f in pp["phenotypicFeatures"][:5]]
    print(f"  Sample terms: {', '.join(features)}...")


def demo_data_loading():
    """Demonstrate loading DECIPHER data."""
    print_header("1. Loading DECIPHER Data")

    # Since we likely don't have real DECIPHER data, use simulator
    print("\nGenerating simulated DECIPHER-like patient data...")
    simulator = DECIPHERSimulator(seed=42)
    patients = simulator.generate_patients(n_patients=200)

    print(f"\nGenerated {len(patients)} patients")

    # Show statistics
    phenotypes_per_patient = [len(p.phenotypes) for p in patients]
    print(f"\nPhenotype distribution:")
    print(f"  Mean phenotypes/patient: {sum(phenotypes_per_patient)/len(phenotypes_per_patient):.1f}")
    print(f"  Min: {min(phenotypes_per_patient)}, Max: {max(phenotypes_per_patient)}")

    # Show example patient
    print("\nExample patient:")
    print_patient_summary(patients[0])

    return patients


def demo_phenopacket_conversion(patients):
    """Demonstrate conversion to Phenopacket format."""
    print_header("2. Converting to Phenopacket Format")

    loader = DECIPHERLoader()
    phenopackets = loader.to_phenopackets(patients)

    print(f"\nConverted {len(phenopackets)} patients to Phenopacket format")

    print("\nExample Phenopacket:")
    print_phenopacket_summary(phenopackets[0])

    # Show metadata
    print("\nMetadata:")
    meta = phenopackets[0]["metaData"]
    print(f"  Schema version: {meta.get('phenopacket_schema_version')}")
    print(f"  Created by: {meta.get('created_by')}")
    if meta.get("external_references"):
        print(f"  External ref: {meta['external_references'][0]['id']}")

    return phenopackets


def demo_statistics(patients):
    """Demonstrate computing dataset statistics."""
    print_header("3. Dataset Statistics")

    loader = DECIPHERLoader()
    stats = loader.get_phenotype_statistics(patients)

    print(f"\nPatient Statistics:")
    print(f"  Total patients: {stats['total_patients']}")
    print(f"  Total phenotype annotations: {stats['total_phenotype_annotations']}")
    print(f"  Unique phenotypes: {stats['unique_phenotypes']}")

    print(f"\nPhenotypes per Patient:")
    ppp = stats['phenotypes_per_patient']
    print(f"  Mean: {ppp['mean']:.1f}")
    print(f"  Min: {ppp['min']}, Max: {ppp['max']}")

    print(f"\nSex Distribution:")
    for sex, count in stats['sex_distribution'].items():
        print(f"  {sex}: {count} ({100*count/stats['total_patients']:.1f}%)")

    print(f"\nRare Phenotypes (appearing once): {stats.get('rare_phenotypes_count', stats.get('rare_phenotypes', 'N/A'))}")

    print(f"\nMost Common Phenotypes:")
    for hpo_id, count in stats['most_common_phenotypes'][:5]:
        print(f"  {hpo_id}: {count} patients")


def demo_similarity_metrics(phenopackets):
    """Demonstrate similarity computation."""
    print_header("4. Computing Phenotype Similarities")

    # Compute IC values
    print("\nComputing information content from corpus...")
    ic_values = compute_empirical_ic(phenopackets)
    print(f"Computed IC for {len(ic_values)} unique HPO terms")

    # Show IC distribution
    ic_sorted = sorted(ic_values.items(), key=lambda x: x[1], reverse=True)
    print("\nMost informative terms (high IC = rare):")
    for term, ic in ic_sorted[:5]:
        print(f"  {term}: {ic:.4f}")

    print("\nLeast informative terms (low IC = common):")
    for term, ic in ic_sorted[-5:]:
        print(f"  {term}: {ic:.4f}")

    # Compare metrics
    metrics = {
        "Jaccard": JaccardSimilarity(),
        "Cosine (IC-weighted)": CosineSimilarity(ic_values),
        "Simplified Resnik": SimplifiedResnikSimilarity(ic_values)
    }

    # Compare patients 0 and 1
    print("\nComparing first two patients:")
    for name, metric in metrics.items():
        calc = PhenopacketSimilarityCalculator(metric)
        sim = calc.compute_similarity(phenopackets[0], phenopackets[1])
        print(f"  {name:25s}: {sim:.4f}")

    return ic_values, metrics


def demo_patient_retrieval(phenopackets, ic_values):
    """Demonstrate finding similar patients."""
    print_header("5. Finding Similar Patients")

    # Use Cosine similarity with IC weighting
    metric = CosineSimilarity(ic_values)
    calc = PhenopacketSimilarityCalculator(metric)

    # Use first patient as query
    query_pp = phenopackets[0]
    print("\nQuery patient:")
    print_phenopacket_summary(query_pp)

    # Find most similar
    print("\nTop 10 most similar patients (Cosine IC-weighted):")
    matches = calc.find_most_similar(query_pp, phenopackets, top_k=10)

    for rank, (idx, score) in enumerate(matches, 1):
        matched_pp = phenopackets[idx]
        print(f"{rank:2d}. Score: {score:.4f} | {matched_pp['id']}")

    # Analyze what makes top match similar
    print("\nAnalyzing top match:")
    top_match = phenopackets[matches[1][0]]  # Skip self-match

    query_terms = set(f["type"]["id"] for f in query_pp["phenotypicFeatures"])
    match_terms = set(f["type"]["id"] for f in top_match["phenotypicFeatures"])

    common = query_terms & match_terms
    print(f"  Query phenotypes: {len(query_terms)}")
    print(f"  Match phenotypes: {len(match_terms)}")
    print(f"  Common phenotypes: {len(common)}")
    print(f"  Overlap: {', '.join(list(common)[:5])}...")


def demo_batch_similarity(phenopackets, ic_values):
    """Demonstrate batch similarity computation."""
    print_header("6. Batch Similarity Analysis")

    # Use subset for speed
    subset = phenopackets[:50]
    print(f"\nComputing 50x50 similarity matrix...")

    metric = CosineSimilarity(ic_values)
    calc = PhenopacketSimilarityCalculator(metric)

    import numpy as np
    matrix = calc.compute_similarity_matrix(subset)

    print(f"\nSimilarity Matrix Statistics:")
    print(f"  Shape: {matrix.shape}")

    # Off-diagonal statistics
    off_diag = matrix[np.triu_indices_from(matrix, k=1)]
    print(f"  Mean similarity: {np.mean(off_diag):.4f}")
    print(f"  Std similarity: {np.std(off_diag):.4f}")
    print(f"  Min similarity: {np.min(off_diag):.4f}")
    print(f"  Max similarity: {np.max(off_diag):.4f}")

    # Find most similar pair
    i, j = np.unravel_index(
        np.argmax(matrix - np.eye(len(matrix))),
        matrix.shape
    )
    print(f"\nMost similar pair:")
    print(f"  {subset[i]['id']} <-> {subset[j]['id']}")
    print(f"  Similarity: {matrix[i, j]:.4f}")


def demo_filtering(patients):
    """Demonstrate patient filtering."""
    print_header("7. Filtering Patients")

    loader = DECIPHERLoader()

    print(f"\nOriginal patient count: {len(patients)}")

    # Filter by phenotype count
    filtered = loader.filter_patients(
        patients,
        min_phenotypes=5,
        max_phenotypes=15
    )
    print(f"After filtering (5-15 phenotypes): {len(filtered)}")

    # Create evaluation cohort
    cohort = loader.create_evaluation_cohort(
        patients,
        n_patients=100,
        balance_by_phenotype_count=True,
        random_seed=42
    )
    print(f"Balanced evaluation cohort: {len(cohort)}")

    # Show cohort distribution
    bins = defaultdict(int)
    for p in cohort:
        n = len(p.phenotypes)
        if n <= 3:
            bins["1-3"] += 1
        elif n <= 6:
            bins["4-6"] += 1
        elif n <= 10:
            bins["7-10"] += 1
        else:
            bins["11+"] += 1

    print("\nCohort phenotype distribution:")
    for bin_name, count in sorted(bins.items()):
        print(f"  {bin_name} phenotypes: {count}")


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 70)
    print("  DECIPHER Data Integration Demo")
    print("  Privacy-Preserving Phenotype Matching")
    print("=" * 70)

    # Load/generate data
    patients = demo_data_loading()

    # Convert to Phenopackets
    phenopackets = demo_phenopacket_conversion(patients)

    # Statistics
    demo_statistics(patients)

    # Similarity metrics
    ic_values, metrics = demo_similarity_metrics(phenopackets)

    # Patient retrieval
    demo_patient_retrieval(phenopackets, ic_values)

    # Batch analysis
    demo_batch_similarity(phenopackets, ic_values)

    # Filtering
    demo_filtering(patients)

    print("\n" + "=" * 70)
    print("  Demo complete!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Obtain DECIPHER data access: https://www.deciphergenomics.org/")
    print("  2. Load real patient data with DECIPHERLoader")
    print("  3. Run privacy experiments with evaluate_decipher.py")
    print("  4. Generate thesis figures and tables")
    print()


if __name__ == "__main__":
    main()
