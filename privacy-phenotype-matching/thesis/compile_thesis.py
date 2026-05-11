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
date: "May 2026"
institution: "Stanford University"
department: "Department of Computer Science"
advisor: "Professor Stephen B. Montgomery, Ph.D."
abstract: |
    Rare diseases affect roughly 300 million people worldwide, and the diagnostic
    odyssey averages 4.8-7 years. Federated patient matching across institutions
    can accelerate diagnosis, but sharing phenotype data is risky: rare phenotype
    combinations act as quasi-identifiers. We present a privacy-preserving
    phenotype matching framework that composes Private Set Intersection (PSI),
    differential privacy (DP), and k-anonymity with rare-term filtering, evaluated
    against a formal two-party semi-honest threat model. On 1,500 real published
    case-report patients from the Monarch Phenopacket Store (Danis et al., 2025),
    non-private IC-weighted cosine retrieval achieves MRR = 0.87 and nDCG@10 = 0.69,
    placing the system within the Phenomizer/LIRICAL band. Shadow-model membership-
    inference attack ROC AUC drops from 0.98 (no DP) to 0.50 (random) at epsilon
    less than or equal to 1; k-anonymity at k = 10 reduces re-identification probability against
    the rare-term adversary from 0.42 to 0.005. The safe DP budget on real patients
    is 20-50x larger than synthetic-cohort experiments suggest, identifying rank-based
    mechanisms as the principled response.
keywords:
    - rare diseases
    - phenotype matching
    - differential privacy
    - private set intersection
    - k-anonymity
    - Human Phenotype Ontology
    - GA4GH Phenopackets
    - membership inference
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


def build_overleaf_bundle(md_content: str, bundle_dir: Path) -> bool:
    """Build an Overleaf-ready directory + zip with main.tex and figures/."""
    import re
    import shutil
    import zipfile

    if bundle_dir.exists():
        shutil.rmtree(bundle_dir)
    bundle_dir.mkdir(parents=True)
    fig_dst = bundle_dir / "figures"
    fig_dst.mkdir()

    # The compile_markdown() helper prepends a plain-markdown title block for
    # the md/html paths. For LaTeX we use \maketitle, so strip that prefix down
    # to the first real `# ` chapter heading to avoid a duplicate top section.
    lines = md_content.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if i > 0 and line.startswith("# ") and "Privacy-Preserving" not in line:
            md_content = "".join(lines[i:])
            break

    # Replace local figures/ references with figures/ (kept the same name) and
    # collect the set of referenced files so we copy only what's used.
    fig_src = THESIS_DIR.parent / "figures"
    referenced = set(re.findall(r"figures/([\w\-.]+\.(?:png|pdf))", md_content))
    copied = 0
    for name in sorted(referenced):
        src = fig_src / name
        if src.exists():
            shutil.copy2(src, fig_dst / name)
            copied += 1
        else:
            print(f"  [warn] figure not found: {name}")
    print(f"  Copied {copied}/{len(referenced)} referenced figures")

    # Always copy the canonical figure set even if not referenced in prose, so
    # the bundle is self-contained for the List of Figures.
    if fig_src.exists():
        for src in sorted(fig_src.glob("fig*.png")) + sorted(fig_src.glob("fig*.pdf")):
            dst = fig_dst / src.name
            if not dst.exists():
                shutil.copy2(src, dst)

    # Pandoc -> main.tex (no standalone wrapper; we provide our own preamble).
    temp_md = bundle_dir / "thesis.temp.md"
    temp_md.write_text(FRONT_MATTER + md_content)
    body_tex = bundle_dir / "body.tex"
    cmd = [
        "pandoc",
        str(temp_md),
        "-o", str(body_tex),
        "-V", "geometry:margin=1in",
        "-V", "fontsize=12pt",
        "--natbib",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    temp_md.unlink()
    if result.returncode != 0:
        print(f"  Pandoc error: {result.stderr}")
        return False

    # Wrap pandoc output in a documentclass + preamble so it builds standalone
    # on Overleaf without further tweaking.
    body = body_tex.read_text()
    body_tex.unlink()
    main_tex = bundle_dir / "main.tex"
    main_tex.write_text(
        r"""\documentclass[12pt, letterpaper]{article}
% Compile with pdfLaTeX. inputenc + fontenc handle the Unicode characters
% pandoc emits (epsilon, ≤, →, etc.). XeLaTeX/LuaLaTeX also work without
% these but may need a different fontenc.
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{textcomp}
\usepackage{lmodern}
\usepackage[margin=1in]{geometry}
\usepackage{setspace}
\onehalfspacing
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\graphicspath{{figures/}}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{calc}
\usepackage[hidelinks]{hyperref}
\usepackage{natbib}
\bibliographystyle{apalike}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhead[L]{\small Privacy-Preserving Phenotype Matching}
\fancyhead[R]{\small Walsh}
\usepackage{caption}
\captionsetup{font=small}
% Pandoc's longtable output uses these helpers.
\providecommand{\tightlist}{%
  \setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\newcommand{\real}[1]{#1}
\newlength{\cslhangindent}
\setlength{\cslhangindent}{1.5em}
\newenvironment{CSLReferences}%
  {}%
  {\par}
\newcommand{\passthrough}[1]{#1}
\title{Privacy-Preserving Phenotype Matching for Rare Disease Cohort Discovery}
\author{Patrick Walsh \\ Department of Computer Science \\ Stanford University}
\date{May 2026}
\begin{document}
\maketitle
\tableofcontents
\newpage
"""
        + body
        + "\n\\end{document}\n"
    )

    # Drop a tiny README so a fresh Overleaf project knows what to compile.
    (bundle_dir / "README.md").write_text(
        "Generated by privacy-phenotype-matching/thesis/compile_thesis.py.\n\n"
        "Overleaf setup:\n"
        "  1. Create a new project -> Upload Project -> select the zip.\n"
        "  2. Menu -> Settings -> Main document = main.tex.\n"
        "  3. Menu -> Settings -> Compiler = pdfLaTeX (default) works.\n"
        "     XeLaTeX/LuaLaTeX also work; if you switch, comment out the\n"
        "     `inputenc` and `fontenc` lines in main.tex.\n\n"
        "Figures live under figures/ (referenced via \\graphicspath).\n"
    )

    # Bundle as a zip alongside the directory for one-shot Overleaf upload.
    zip_path = bundle_dir.with_suffix(".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(bundle_dir.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(bundle_dir.parent))
    print(f"  Bundle directory: {bundle_dir}")
    print(f"  Overleaf zip:     {zip_path}")
    return True


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
        choices=["md", "html", "latex", "pdf", "overleaf", "all"],
        default="md",
        help="Output format (default: md). 'overleaf' produces a zip ready to "
             "upload to Overleaf (main.tex + figures/)."
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
    formats = ["md", "html", "latex", "pdf", "overleaf"] if args.format == "all" else [args.format]

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
        elif fmt == "overleaf":
            bundle_dir = OUTPUT_DIR / f"thesis_{date_str}_overleaf"
            if args.output:
                bundle_dir = OUTPUT_DIR / f"{args.output}_overleaf"
            build_overleaf_bundle(md_content, bundle_dir)

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
