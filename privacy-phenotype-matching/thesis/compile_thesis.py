#!/usr/bin/env python3
"""
Thesis Compilation Script

Compiles individual thesis sections into a single document.
Supports multiple output formats: Markdown, HTML, LaTeX, and PDF.

Usage:
    python compile_thesis.py [--format FORMAT] [--output OUTPUT]

Formats:
    md      - Single Markdown file (default)
    html    - HTML with styling
    latex   - LaTeX document
    pdf     - PDF via pandoc (requires pandoc and LaTeX)

Requirements:
    - pandoc (for HTML, LaTeX, PDF conversion)
    - pdflatex or xelatex (for PDF output)

Install pandoc:
    brew install pandoc        # macOS
    apt install pandoc         # Ubuntu/Debian
    choco install pandoc       # Windows
"""

import argparse
import subprocess
import sys
from pathlib import Path
from datetime import datetime


THESIS_DIR = Path(__file__).parent
SECTIONS_DIR = THESIS_DIR / "sections"
OUTPUT_DIR = THESIS_DIR / "output"

# Section order for compilation
SECTION_ORDER = [
    "introduction.md",
    "literature_review.md",
    "methods.md",
    "results.md",
    "discussion.md",
    "future_work.md",
    "conclusion.md",
    "bibliography.md",
    "appendices.md",
]

# Front matter template
FRONT_MATTER = """---
title: "Privacy-Preserving Phenotype Matching for Rare Disease Cohort Discovery"
author: "Patrick Walsh"
date: "April 2026"
institution: "Stanford University"
department: "Department of Computer Science"
advisor: "Professor Stephen B. Montgomery, Ph.D."
abstract: |
    Rare diseases affect approximately 300 million people worldwide, yet patients
    endure an average diagnostic odyssey of 4.8-7 years. We present a privacy-preserving
    phenotype matching framework that enables federated rare disease cohort discovery
    while protecting patient confidentiality through Private Set Intersection,
    differential privacy, and k-anonymity. Evaluated on 12,974 real disease profiles,
    our system achieves nDCG@10 = 99.7% while providing quantifiable privacy guarantees.
keywords:
    - rare diseases
    - phenotype matching
    - privacy-preserving computation
    - private set intersection
    - differential privacy
    - Human Phenotype Ontology
---

"""

# LaTeX preamble for academic formatting
LATEX_PREAMBLE = r"""
\documentclass[12pt, letterpaper]{article}
\usepackage[margin=1in]{geometry}
\usepackage{setspace}
\doublespacing
\usepackage{amsmath, amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{hyperref}
\usepackage{natbib}
\bibliographystyle{apalike}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhead[L]{\small Privacy-Preserving Phenotype Matching}
\fancyhead[R]{\small Walsh}
\usepackage{caption}
\captionsetup{font=small}
"""


def compile_markdown() -> str:
    """Compile all sections into a single Markdown document."""
    output = []

    # Add title page
    output.append("# Privacy-Preserving Phenotype Matching for Rare Disease Cohort Discovery\n")
    output.append("**Patrick Walsh**\n")
    output.append("*Honors Thesis*\n")
    output.append("Department of Computer Science, Stanford University\n")
    output.append("April 2026\n")
    output.append("\n**Thesis Advisor:** Professor Stephen B. Montgomery, Ph.D.\n")
    output.append("\n---\n")

    # Add each section
    for section_file in SECTION_ORDER:
        section_path = SECTIONS_DIR / section_file
        if section_path.exists():
            print(f"  Adding: {section_file}")
            content = section_path.read_text()
            output.append(content)
            output.append("\n\n---\n\n")
        else:
            print(f"  [MISSING] {section_file}")
            output.append(f"\n\n<!-- Section {section_file} not yet written -->\n\n")

    return "\n".join(output)


def compile_to_html(md_content: str, output_path: Path) -> bool:
    """Convert Markdown to HTML using pandoc."""
    try:
        # Write temp markdown file
        temp_md = output_path.with_suffix(".temp.md")
        temp_md.write_text(md_content)

        cmd = [
            "pandoc",
            str(temp_md),
            "-o", str(output_path),
            "--standalone",
            "--toc",
            "--toc-depth=3",
            "-c", "https://cdn.jsdelivr.net/npm/water.css@2/out/water.css",
            "--metadata", "title=Privacy-Preserving Phenotype Matching",
            "--mathjax",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        temp_md.unlink()  # Clean up

        if result.returncode != 0:
            print(f"Pandoc error: {result.stderr}")
            return False
        return True
    except FileNotFoundError:
        print("Error: pandoc not found. Install with: brew install pandoc")
        return False


def compile_to_latex(md_content: str, output_path: Path) -> bool:
    """Convert Markdown to LaTeX using pandoc."""
    try:
        temp_md = output_path.with_suffix(".temp.md")
        temp_md.write_text(FRONT_MATTER + md_content)

        cmd = [
            "pandoc",
            str(temp_md),
            "-o", str(output_path),
            "--standalone",
            "--toc",
            "-V", "geometry:margin=1in",
            "-V", "fontsize=12pt",
            "-V", "documentclass=article",
            "--natbib",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        temp_md.unlink()

        if result.returncode != 0:
            print(f"Pandoc error: {result.stderr}")
            return False
        return True
    except FileNotFoundError:
        print("Error: pandoc not found. Install with: brew install pandoc")
        return False


def compile_to_pdf(md_content: str, output_path: Path) -> bool:
    """Convert Markdown to PDF using pandoc."""
    try:
        temp_md = output_path.with_suffix(".temp.md")
        temp_md.write_text(FRONT_MATTER + md_content)

        cmd = [
            "pandoc",
            str(temp_md),
            "-o", str(output_path),
            "--toc",
            "-V", "geometry:margin=1in",
            "-V", "fontsize=12pt",
            "-V", "documentclass=article",
            "--pdf-engine=xelatex",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        temp_md.unlink()

        if result.returncode != 0:
            print(f"Pandoc error: {result.stderr}")
            return False
        return True
    except FileNotFoundError:
        print("Error: pandoc not found. Install with: brew install pandoc")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Compile thesis sections into a single document"
    )
    parser.add_argument(
        "--format", "-f",
        choices=["md", "html", "latex", "pdf", "all"],
        default="md",
        help="Output format (default: md)"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Output filename (default: thesis_YYYYMMDD.FORMAT)"
    )
    args = parser.parse_args()

    # Create output directory
    OUTPUT_DIR.mkdir(exist_ok=True)

    # Generate default filename
    date_str = datetime.now().strftime("%Y%m%d")

    print("=" * 60)
    print("  Thesis Compilation")
    print("=" * 60)
    print(f"\nSections directory: {SECTIONS_DIR}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    # Check which sections exist
    print("Checking sections...")
    existing = [s for s in SECTION_ORDER if (SECTIONS_DIR / s).exists()]
    missing = [s for s in SECTION_ORDER if not (SECTIONS_DIR / s).exists()]

    print(f"  Found: {len(existing)} sections")
    print(f"  Missing: {len(missing)} sections")
    if missing:
        print(f"  Missing: {', '.join(missing)}")

    # Compile markdown
    print("\nCompiling sections...")
    md_content = compile_markdown()

    # Output based on format
    formats = ["md", "html", "latex", "pdf"] if args.format == "all" else [args.format]

    for fmt in formats:
        if args.output:
            output_path = OUTPUT_DIR / f"{args.output}.{fmt}"
        else:
            output_path = OUTPUT_DIR / f"thesis_{date_str}.{fmt}"

        print(f"\nGenerating {fmt.upper()}...")

        if fmt == "md":
            output_path.write_text(md_content)
            print(f"  Saved: {output_path}")
        elif fmt == "html":
            if compile_to_html(md_content, output_path):
                print(f"  Saved: {output_path}")
        elif fmt == "latex":
            if compile_to_latex(md_content, output_path):
                print(f"  Saved: {output_path}")
        elif fmt == "pdf":
            if compile_to_pdf(md_content, output_path):
                print(f"  Saved: {output_path}")

    print("\n" + "=" * 60)
    print("  Compilation complete!")
    print("=" * 60)

    # Print status summary
    print(f"\nThesis Status:")
    print(f"  Sections complete: {len(existing)}/{len(SECTION_ORDER)}")
    if missing:
        print(f"\n  Still needed:")
        for m in missing:
            print(f"    - {m}")


if __name__ == "__main__":
    main()
