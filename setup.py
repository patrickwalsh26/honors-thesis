"""Setup script for privacy-phenotype-matching package."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="privacy-phenotype-matching",
    version="0.1.0",
    author="Patrick Walsh",
    author_email="walshp26@stanford.edu",
    description="Privacy-Preserving Phenotype Matching for Rare Disease Cohorts",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/patrickwalsh/privacy-phenotype-matching",  # Update when repo created
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Security :: Cryptography",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.4.0",
            "pytest-cov>=4.1.0",
            "black>=23.7.0",
            "flake8>=6.0.0",
            "mypy>=1.4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "generate-phenopackets=src.data_generation.synthetic_phenopackets:main",
        ],
    },
)
