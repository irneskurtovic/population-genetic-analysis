"""Pairwise population FST analysis for a chromosome 22 VCF subset.

Inputs expected in population_variant_analysis/data/:
- chr22_first1000.vcf
- integrated_call_samples_v3.20130502.ALL.panel

Outputs are written to population_variant_analysis/final/fst/.
This is an exploratory analysis of a limited subset, not a genome-wide estimate.
"""
from itertools import combinations
from pathlib import Path

import allel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
VCF_PATH = DATA_DIR / "chr22_first1000.vcf"
PANEL_PATH = DATA_DIR / "integrated_call_samples_v3.20130502.ALL.panel"
OUTPUT_DIR = PROJECT_DIR / "final" / "fst"
MAX_MISSING_RATE = 0.05


def decode_text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def load_panel(panel_path, sample_ids):
    panel = pd.read_csv(panel_path, sep=r"\s+")
    required = {"sample", "pop"}
    missing_columns = required - set(panel.columns)
    if missing_columns:
        raise ValueError(f"Panel is missing required columns: {sorted(missing_columns)}")
    if panel["sample"].duplicated().any():
        raise ValueError("Panel contains duplicate sample IDs.")

    panel = panel.set_index("sample")
    missing_samples = sorted(set(sample_ids) - set(panel.index))
    if missing_samples:
        raise ValueError(
            f"{len(missing_samples)} VCF samples were not found in the panel. "
            f"Examples: {missing_samples[:10]}"
        )
    return panel.loc[sample_ids].reset_index()


def main():
    if not VCF_PATH.exists():
        raise FileNotFoundError(f"VCF file not found: {VCF_PATH}")
    if not PANEL_PATH.exists():
        raise FileNotFoundError(f"Population panel not found: {PANEL_PATH}")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    callset = allel.read_vcf(
        str(VCF_PATH),
        fields=["samples", "variants/REF", "variants/ALT", "calldata/GT"],
    )
    if callset is None or "calldata/GT" not in callset:
        raise ValueError("Could not read genotype data from the VCF.")

    sample_ids = [decode_text(x) for x in callset["samples"]]
    metadata = load_panel(PANEL_PATH, sample_ids)
    population_labels = metadata["pop"].astype(str).to_numpy()
    populations = sorted(metadata["pop"].astype(str).unique())

    refs = [decode_text(x) for x in callset["variants/REF"]]
    alts = callset["variants/ALT"]
    biallelic_snp_indices = []
    for i, (ref, alt_row) in enumerate(zip(refs, alts)):
        alt_values = [decode_text(a).strip() for a in alt_row]
        alt_values = [a for a in alt_values if a and a != "."]
        if len(ref) == 1 and len(alt_values) == 1 and len(alt_values[0]) == 1:
            biallelic_snp_indices.append(i)

    if not biallelic_snp_indices:
        raise ValueError("No biallelic SNPs found in the VCF subset.")

    genotypes = allel.GenotypeArray(callset["calldata/GT"])[biallelic_snp_indices]
    # scikit-allel returns a 2D (variants x samples) array from is_called().
    # Each value indicates whether that sample has a complete genotype call.
    called = genotypes.is_called()
    call_rate = called.mean(axis=1)
    keep = call_rate >= (1 - MAX_MISSING_RATE)
    genotypes = genotypes[keep]
    if genotypes.shape[0] < 2:
        raise ValueError("Fewer than two biallelic SNPs passed the call-rate filter.")

    # Compute pairwise Weir-Cockerham FST, aggregating variance components across sites.
    fst_matrix = pd.DataFrame(np.nan, index=populations, columns=populations, dtype=float)
    # Assign diagonal through pandas to avoid writing to a read-only NumPy view.
    for pop in populations:
        fst_matrix.loc[pop, pop] = 0.0
    pair_results = []
    sample_indices = {
        pop: np.flatnonzero(population_labels == pop) for pop in populations
    }

    for pop1, pop2 in combinations(populations, 2):
        subpops = [sample_indices[pop1], sample_indices[pop2]]
        a, b, c = allel.weir_cockerham_fst(genotypes, subpops, max_allele=1)
        numerator = np.nansum(a)
        denominator = np.nansum(a + b + c)
        fst = float(numerator / denominator) if denominator != 0 else np.nan
        fst_matrix.loc[pop1, pop2] = fst
        fst_matrix.loc[pop2, pop1] = fst
        pair_results.append({"population_1": pop1, "population_2": pop2, "fst": fst})

    matrix_path = OUTPUT_DIR / "pairwise_fst_matrix.csv"
    pairs_path = OUTPUT_DIR / "pairwise_fst_results.csv"
    plot_path = OUTPUT_DIR / "pairwise_fst_heatmap.png"
    fst_matrix.to_csv(matrix_path, index_label="population")
    pd.DataFrame(pair_results).to_csv(pairs_path, index=False)

    values = fst_matrix.to_numpy(dtype=float)
    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        raise ValueError("No finite pairwise FST estimates were produced.")
    limit = max(abs(float(finite_values.min())), abs(float(finite_values.max())), 0.01)

    fig, ax = plt.subplots(figsize=(12, 10))
    image = ax.imshow(values, cmap="coolwarm", vmin=-limit, vmax=limit, aspect="auto")
    ax.set_xticks(range(len(populations)))
    ax.set_yticks(range(len(populations)))
    ax.set_xticklabels(populations, rotation=90)
    ax.set_yticklabels(populations)
    ax.set_title("Pairwise Weir–Cockerham FST (chromosome 22 VCF subset)")
    fig.colorbar(image, ax=ax, label="FST estimate")
    fig.tight_layout()
    fig.savefig(plot_path, dpi=300)
    plt.close(fig)

    print("Pairwise FST analysis completed successfully.")
    print(f"Samples analysed: {len(sample_ids)}")
    print(f"Populations analysed: {len(populations)}")
    print(f"Biallelic SNPs initially selected: {len(biallelic_snp_indices)}")
    print(f"SNPs used after call-rate filtering: {genotypes.shape[0]}")
    print(f"Population pairs analysed: {len(pair_results)}")
    print(f"Results saved to: {OUTPUT_DIR}")
    print("Note: estimates use a limited chromosome 22 subset; interpret as exploratory, not genome-wide.")


if __name__ == "__main__":
    main()