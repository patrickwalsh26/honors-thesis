#!/usr/bin/env python3
"""
Demonstration of baseline phenotype matching capabilities.

This script shows:
1. Loading synthetic phenopackets
2. Computing various similarity metrics
3. Finding similar patients
4. Evaluating retrieval performance
"""

import sys
from pathlib import Path
import json
import numpy as np
from collections import defaultdict

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.similarity.hpo_similarity import (
    JaccardSimilarity,
    CosineSimilarity,
    SimplifiedResnikSimilarity,
    PhenopacketSimilarityCalculator,
    compute_empirical_ic,
    load_phenopackets
)


def print_header(text):
    """Print a formatted header."""
    print("\n" + "=" * 70)
    print(f"  {text}")
    print("=" * 70)


def print_phenopacket_summary(pp):
    """Print a summary of a phenopacket."""
    disease = pp["diseases"][0]["term"]["label"] if pp["diseases"] else "Unknown"
    n_features = len(pp["phenotypicFeatures"])
    features = [f["type"]["id"] for f in pp["phenotypicFeatures"][:5]]

    print(f"  ID: {pp['id']}")
    print(f"  Disease: {disease}")
    print(f"  Features: {n_features} total")
    print(f"  Sample terms: {', '.join(features)}...")


def demo_data_loading():
    """Demonstrate loading and examining phenopackets."""
    print_header("1. Loading Synthetic Phenopackets")

    # Load mixed cohort
    phenopackets = load_phenopackets("data/synthetic/mixed_cohort_200.json")
    print(f"\nLoaded {len(phenopackets)} phenopackets")

    # Group by disease
    by_disease = defaultdict(list)
    for pp in phenopackets:
        disease = pp["diseases"][0]["term"]["label"] if pp["diseases"] else "Unknown"
        by_disease[disease].append(pp)

    print("\nCohort composition:")
    for disease, pps in by_disease.items():
        print(f"  {disease}: {len(pps)} patients")

    # Show example phenopacket
    print("\nExample phenopacket:")
    print_phenopacket_summary(phenopackets[0])

    return phenopackets


def demo_similarity_metrics(phenopackets):
    """Demonstrate different similarity metrics."""
    print_header("2. Computing Similarity Metrics")

    # Compute IC values from the corpus
    print("\nComputing information content from corpus...")
    ic_values = compute_empirical_ic(phenopackets)
    print(f"Computed IC for {len(ic_values)} unique HPO terms")

    # Show top/bottom IC terms
    sorted_ic = sorted(ic_values.items(), key=lambda x: x[1], reverse=True)
    print("\nMost informative terms (high IC = rare):")
    for term, ic in sorted_ic[:5]:
        print(f"  {term}: {ic:.4f}")

    print("\nLeast informative terms (low IC = common):")
    for term, ic in sorted_ic[-5:]:
        print(f"  {term}: {ic:.4f}")

    # Initialize different similarity metrics
    metrics = {
        "Jaccard": JaccardSimilarity(),
        "Cosine (IC-weighted)": CosineSimilarity(ic_values),
        "Simplified Resnik": SimplifiedResnikSimilarity(ic_values)
    }

    # Compare metrics on sample pairs
    print("\nComparing metrics on sample patient pairs:")
    print(f"\nSame disease (Marfan #0 vs #1):")
    marfan_patients = [pp for pp in phenopackets if "Marfan" in pp["diseases"][0]["term"]["label"]]

    for name, metric in metrics.items():
        calc = PhenopacketSimilarityCalculator(metric)
        sim = calc.compute_similarity(marfan_patients[0], marfan_patients[1])
        print(f"  {name:25s}: {sim:.4f}")

    print(f"\nDifferent diseases (Marfan vs Achondroplasia):")
    achon_patients = [pp for pp in phenopackets if "Achondroplasia" in pp["diseases"][0]["term"]["label"]]

    for name, metric in metrics.items():
        calc = PhenopacketSimilarityCalculator(metric)
        sim = calc.compute_similarity(marfan_patients[0], achon_patients[0])
        print(f"  {name:25s}: {sim:.4f}")

    return ic_values, metrics


def demo_retrieval(phenopackets, metrics):
    """Demonstrate patient retrieval."""
    print_header("3. Finding Similar Patients")

    # Use Marfan patient as query
    query_pp = None
    for pp in phenopackets:
        if "Marfan" in pp["diseases"][0]["term"]["label"]:
            query_pp = pp
            break

    print("\nQuery patient:")
    print_phenopacket_summary(query_pp)

    # Search for similar patients using Cosine similarity
    ic_values = compute_empirical_ic(phenopackets)
    metric = CosineSimilarity(ic_values)
    calc = PhenopacketSimilarityCalculator(metric)

    print("\nTop 10 most similar patients (Cosine IC-weighted):")
    matches = calc.find_most_similar(query_pp, phenopackets, top_k=10)

    for rank, (idx, score) in enumerate(matches, 1):
        matched_pp = phenopackets[idx]
        disease = matched_pp["diseases"][0]["term"]["label"]
        print(f"{rank:2d}. Score: {score:.4f} | {matched_pp['id']:20s} | {disease}")


def demo_retrieval_evaluation(phenopackets):
    """Demonstrate retrieval evaluation metrics."""
    print_header("4. Evaluating Retrieval Performance")

    ic_values = compute_empirical_ic(phenopackets)
    metric = SimplifiedResnikSimilarity(ic_values)
    calc = PhenopacketSimilarityCalculator(metric)

    # For each disease, use one patient as query and measure recall
    diseases = defaultdict(list)
    for pp in phenopackets:
        disease = pp["diseases"][0]["term"]["label"]
        diseases[disease].append(pp)

    print("\nRecall@k for retrieving same-disease patients:")
    print(f"{'Disease':<40s} R@5    R@10   R@20")
    print("-" * 70)

    for disease, disease_patients in diseases.items():
        if len(disease_patients) < 2:
            continue

        # Use first patient as query
        query = disease_patients[0]
        relevant_ids = {pp["id"] for pp in disease_patients[1:]}

        # Find top matches in full cohort
        matches = calc.find_most_similar(query, phenopackets, top_k=20)

        # Compute recall@k
        recalls = {}
        for k in [5, 10, 20]:
            top_k_ids = {phenopackets[idx]["id"] for idx, _ in matches[:k]}
            retrieved_relevant = len(top_k_ids & relevant_ids)
            recall = retrieved_relevant / len(relevant_ids) if relevant_ids else 0
            recalls[k] = recall

        print(f"{disease:<40s} {recalls[5]:.3f}  {recalls[10]:.3f}  {recalls[20]:.3f}")

    # Overall statistics
    print("\nOverall baseline performance:")
    print("  - Similarity metrics: Jaccard, Cosine, Simplified Resnik")
    print("  - IC computation: Empirical from synthetic corpus")
    print("  - Next steps: Add privacy-preserving protocols (PSI, DP)")


def demo_similarity_matrix(phenopackets):
    """Demonstrate computing full similarity matrix."""
    print_header("5. Similarity Matrix Visualization")

    # Use small subset for speed
    subset = phenopackets[:20]

    ic_values = compute_empirical_ic(phenopackets)
    metric = CosineSimilarity(ic_values)
    calc = PhenopacketSimilarityCalculator(metric)

    print(f"\nComputing 20x20 similarity matrix...")
    matrix = calc.compute_similarity_matrix(subset)

    print(f"Matrix shape: {matrix.shape}")
    print(f"Mean similarity: {np.mean(matrix[np.triu_indices_from(matrix, k=1)]):.4f}")
    print(f"Max similarity (off-diagonal): {np.max(matrix - np.eye(len(matrix))):.4f}")
    print(f"Min similarity: {np.min(matrix):.4f}")

    # Show which patients are most similar on average
    mean_sims = matrix.mean(axis=1)
    most_similar_idx = np.argmax(mean_sims)
    print(f"\nMost 'typical' patient (highest avg similarity): {subset[most_similar_idx]['id']}")
    print(f"Average similarity: {mean_sims[most_similar_idx]:.4f}")


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 70)
    print("  Privacy-Preserving Phenotype Matching - Baseline Demo")
    print("  Stanford CS Honors Thesis 2024-2025")
    print("=" * 70)

    # Run demonstrations
    phenopackets = demo_data_loading()
    ic_values, metrics = demo_similarity_metrics(phenopackets)
    demo_retrieval(phenopackets, metrics)
    demo_retrieval_evaluation(phenopackets)
    demo_similarity_matrix(phenopackets)

    print("\n" + "=" * 70)
    print("  Demo complete!")
    print("=" * 70)
    print("\nNext steps:")
    print("  1. Implement privacy-preserving protocols (PSI)")
    print("  2. Add differential privacy mechanisms")
    print("  3. Evaluate privacy-utility tradeoffs")
    print("  4. Integrate with GA4GH standards (MME, Beacon)")
    print()


if __name__ == "__main__":
    main()
