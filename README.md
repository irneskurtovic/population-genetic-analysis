# Population Variant Analysis

An exploratory bioinformatics project using a subset of 1000 Genomes Project
Phase 3 chromosome 22 genotype data. The project compares alternate-allele
frequencies among populations, visualizes sample relationships with principal
component analysis (PCA), and estimates pairwise population differentiation
with Weir–Cockerham FST.

## Objectives

- Calculate REF/ALT allele counts and ALT frequencies per population.
- Explore sample structure using PCA on filtered biallelic SNPs.
- Estimate pairwise FST between populations and visualize the estimates.
- Keep the analysis reproducible with scripts that resolve paths relative to
  the project location.

## Dataset

The project uses files associated with the [1000 Genomes Project Phase 3 collection](https://www.internationalgenome.org/data-portal/data-collection/phase3),
maintained by the International Genome Sample Resource (IGSR). Phase 3 uses
the GRCh37 reference build. This project analyzes a limited subset of
chromosome 22, not the full chromosome or genome. Inspection of the local VCF
confirms 1,000 variant records and the header values
`##source=1000GenomesPhase3Pipeline` and `b37`. Its exact source VCF filename
and the method used to select these 1,000 records are unknown.

### Input files and provenance

- `data/chr22_first1000.vcf` is a local 1,000-record chromosome 22 subset. Its
  header identifies `1000GenomesPhase3Pipeline` and `b37`. The relevant
  collection and release information are available from the [IGSR Phase 3 collection](https://www.internationalgenome.org/data-portal/data-collection/phase3)
  and the [20130502 release directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
  The exact source VCF filename and selection method are not recorded in the
  local file or project documentation.
- `data/integrated_call_samples_v3.20130502.ALL.panel` is the Phase 3 sample
  panel format, with sample, population, super-population, and gender columns.
  Its filename matches the Phase 3 panel naming. The [IGSR Phase 3 collection](https://www.internationalgenome.org/data-portal/data-collection/phase3)
  and [20130502 release directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
  provide source context. The local panel contains 2,504 sample IDs, all
  matching the VCF sample IDs. A download record or checksum for this local
  copy is not available in the project.

### Reuse and distribution

IGSR states that 1000 Genomes Project data are available without embargo after
the final publication and asks users to cite the project. IGSR also cautions
that rights and restrictions can vary among data in its collections, and that
users must check applicable data reuse statements and third-party rights.
Public availability does not automatically establish permission to
redistribute these local input copies or derived outputs in a GitHub
repository. No redistribution licence is asserted here. Before publishing,
users must verify the applicable reuse terms for the source files and outputs.
See the [IGSR disclaimer](https://www.internationalgenome.org/IGSR_disclaimer/), the
[Phase 3 collection](https://www.internationalgenome.org/data-portal/data-collection/phase3),
the [20130502 release directory](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/),
and [IGSR citation guidance](https://www.internationalgenome.org/faq/how-do-I-cite-IGSR/) for source,
reuse, and citation information.

## Project structure

```text
population_variant_analysis/
├── data/
│   ├── analyze_population.py
│   ├── population_fst.py
│   ├── population_pca.py
│   ├── chr22_first1000.vcf
│   └── integrated_call_samples_v3.20130502.ALL.panel
├── final/
│   ├── Population_Analysis_Report.xlsx
│   ├── fst/
│   │   ├── pairwise_fst_heatmap.png
│   │   ├── pairwise_fst_matrix.csv
│   │   └── pairwise_fst_results.csv
│   └── pca/
│       ├── pca_sample_coordinates.csv
│       ├── pca_variants_used.csv
│       └── population_pca.png
├── .gitignore
├── README.md
└── requirements.txt
```

The three scripts read the VCF and panel from `data/`. They create or write
results beneath `final/`. The existing Excel report is a project artifact; it
is not written by these three scripts.

## Methods

### Population variant frequencies

`analyze_population.py` counts called REF and ALT alleles for records with a
single non-missing ALT allele in each population. It calculates per-population
ALT frequencies and aggregates allele counts for an overall frequency
summary. It also creates plots for a configured target variant (when present)
and the ten records with the lowest overall ALT frequency.

### Exploratory PCA

`population_pca.py` selects biallelic SNPs, applies a maximum 5% missing-call
rate and a 5% minimum minor allele frequency, mean-imputes remaining missing
dosages, standardizes variants, and computes two principal components. It does
not perform linkage-disequilibrium (LD) pruning.

### Pairwise FST

`population_fst.py` selects biallelic SNPs, filters on a maximum 5% missing-call
rate, then calculates pairwise Weir–Cockerham FST estimates between populations
and saves a matrix, pair table, and heatmap.

## Installation

From the project root, create and activate a virtual environment, then install
the listed dependencies:

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

On macOS or Linux, activate the environment with
`source .venv/bin/activate` after creating it with `python3 -m venv .venv`.

## Run the analyses

Run these commands from the project root after installing the dependencies:

```powershell
python data/analyze_population.py
python data/population_pca.py
python data/population_fst.py
```

The scripts use paths based on their own file locations, so the working
directory is not used to locate the input data. Running them from the project
root as shown keeps invocation consistent.

## Output files

| Script | Output path(s) |
|---|---|
| `analyze_population.py` | `final/population_variant_frequencies.csv`, `final/variant_frequency_summary.csv`, `final/lowest_10_variant_frequencies.png`, and, if the configured target is found, `final/rs529371290_population_frequency.png` |
| `population_pca.py` | `final/pca/pca_sample_coordinates.csv`, `final/pca/pca_variants_used.csv`, `final/pca/population_pca.png` |
| `population_fst.py` | `final/fst/pairwise_fst_matrix.csv`, `final/fst/pairwise_fst_results.csv`, `final/fst/pairwise_fst_heatmap.png` |

The scripts write these named outputs when run and may replace existing files
with the same names. Existing files in `final/` are retained in the repository.

## Limitations

- The supplied VCF contains only a subset of chromosome 22 variant records, so
  the results are not genome-wide estimates and may not represent the complete
  Phase 3 data.
- PCA is exploratory and is not genome-wide. Only two components are reported,
  and the script does not perform LD pruning.
- Frequency and FST estimates depend on the supplied subset and the scripts'
  filters; they should not be interpreted as comprehensive population
  estimates.
- The analysis assumes the sample identifiers in the VCF match those in the
  population panel.

## Reproducibility and future improvements

The scripts resolve input and output paths relative to their own locations,
validate key input fields and sample matching, and document important filters
in code. Install dependencies from `requirements.txt` and run each command
above to regenerate the corresponding outputs. For stronger reproducibility,
future work could record dependency versions, document exact source-file
accessions/checksums and subset-generation steps, add automated validation, use
the complete chromosome data with appropriate compute resources, and apply LD
pruning and sensitivity checks for PCA and FST.

## Data citation

Please cite the 1000 Genomes Project and follow the IGSR data reuse guidance
when reusing these data. IGSR recommends the final Phase 3 paper: The 1000
Genomes Project Consortium, “A global reference for human genetic variation,”
*Nature* 526, 68–74 (2015), doi:10.1038/nature15393. See the [Phase 3 collection page](https://www.internationalgenome.org/data-portal/data-collection/phase3)
and [IGSR citation guidance](https://www.internationalgenome.org/faq/how-do-I-cite-IGSR/) for release
information, data reuse policies, and associated publications.
