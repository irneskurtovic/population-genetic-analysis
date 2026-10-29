# Population Genetic Analysis

## Overview
This project explores genetic variation across populations using data from the 1000 Genomes Project. The analysis focuses on chromosome 22 and includes population-level allele frequencies, Principal Component Analysis (PCA), and pairwise population differentiation using FST.

## Objectives
- Calculate alternate allele frequencies across populations.
- Explore genetic structure using PCA.
- Compare genetic differentiation between populations using pairwise FST.
- Visualize and export the results for further analysis.

## Dataset
- **Source:** 1000 Genomes Project, Phase 3
- **Samples:** 2,504 individuals
- **Reference data:** A subset of chromosome 22 variants
- **Population metadata:** 1000 Genomes population panel

This is an exploratory analysis of a limited chromosome 22 dataset, not a genome-wide population genetics study.

## Methods
1. **Allele frequency analysis:** Calculate alternate allele frequencies for each variant and population.
2. **PCA:** Examine genetic variation and sample clustering using principal components.
3. **Pairwise FST:** Estimate genetic differentiation between population pairs.

## Repository Structure
- `data/` — input dataset, population metadata, and Python scripts.
- `final/` — exported tables, figures, and analysis reports.
- `requirements.txt` — Python dependencies.
- `.gitignore` — files excluded from version control.

## Reproducibility
See `data/README.md` for descriptions of the input files and analysis scripts. Install the dependencies listed in `requirements.txt` before running the scripts.

## Limitations
The analysis uses a limited subset of chromosome 22 variants. PCA results may be affected by variant filtering and linkage disequilibrium. Findings should therefore be interpreted as exploratory and not as a comprehensive representation of genome-wide population structure.

## Tools
Python, scikit-allel, pandas, NumPy, Matplotlib, and scikit-learn.
