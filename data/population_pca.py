"""Exploratory PCA of 1000 Genomes chromosome 22 genotype data.

Input files expected in the project's data/ directory:
- chr22_first1000.vcf
- integrated_call_samples_v3.20130502.ALL.panel

Outputs are written to final/pca/.
This is an exploratory analysis of a limited VCF subset, not a genome-wide
population-structure analysis. LD pruning is not performed by this script.
"""
from pathlib import Path

import allel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "final" / "pca"
VCF_PATH = DATA_DIR / "chr22_first1000.vcf"
PANEL_PATH = DATA_DIR / "integrated_call_samples_v3.20130502.ALL.panel"

MAX_MISSING_RATE = 0.05
MIN_MAF = 0.05


def decode_text(value):
    """Convert byte/string VCF fields to ordinary strings."""
    return value.decode() if isinstance(value, bytes) else str(value)


def load_panel(panel_path, sample_ids):
    panel = pd.read_csv(panel_path, sep=r"\s+")
    required = {"sample", "pop", "super_pop"}
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
        raise FileNotFoundError(f"Panel file not found: {PANEL_PATH}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    callset = allel.read_vcf(
        str(VCF_PATH),
        fields=[
            "samples",
            "variants/CHROM",
            "variants/POS",
            "variants/ID",
            "variants/REF",
            "variants/ALT",
            "calldata/GT",
        ],
    )
    if callset is None or "calldata/GT" not in callset:
        raise ValueError("Could not read genotype data from the VCF.")

    sample_ids = [decode_text(x) for x in callset["samples"]]
    sample_metadata = load_panel(PANEL_PATH, sample_ids)

    refs = [decode_text(x) for x in callset["variants/REF"]]
    alts = callset["variants/ALT"]
    keep_indices = []
    for i, (ref, alt_row) in enumerate(zip(refs, alts)):
        alt_values = [decode_text(a).strip() for a in alt_row]
        alt_values = [a for a in alt_values if a and a != "."]
        # PCA uses biallelic SNPs only; exclude indels and multiallelic sites.
        if len(ref) == 1 and len(alt_values) == 1 and len(alt_values[0]) == 1:
            keep_indices.append(i)

    if not keep_indices:
        raise ValueError("No biallelic SNPs found in the VCF subset.")

    gt = allel.GenotypeArray(callset["calldata/GT"])[keep_indices]
    # Dosage: 0, 1, or 2 ALT copies; missing genotypes become -1.
    dosage = gt.to_n_alt(fill=-1).astype(float)
    positions = np.asarray(callset["variants/POS"])[keep_indices]
    variant_ids = [decode_text(x) for x in np.asarray(callset["variants/ID"])[keep_indices]]
    chroms = [decode_text(x) for x in np.asarray(callset["variants/CHROM"])[keep_indices]]

    # Filter SNPs by call rate, then by minor allele frequency.
    called = dosage >= 0
    call_rate = called.mean(axis=1)
    keep_call_rate = call_rate >= (1 - MAX_MISSING_RATE)
    dosage = dosage[keep_call_rate]
    positions = positions[keep_call_rate]
    variant_ids = np.asarray(variant_ids)[keep_call_rate]
    chroms = np.asarray(chroms)[keep_call_rate]
    called = dosage >= 0

    if dosage.shape[0] == 0:
        raise ValueError("No SNPs passed the call-rate filter.")

    # Impute missing dosage values with each SNP's mean called dosage.
    means = np.array([row[row >= 0].mean() if np.any(row >= 0) else np.nan for row in dosage])
    valid_mean = np.isfinite(means)
    dosage = dosage[valid_mean]
    positions = positions[valid_mean]
    variant_ids = variant_ids[valid_mean]
    chroms = chroms[valid_mean]
    means = means[valid_mean]
    for i in range(dosage.shape[0]):
        dosage[i, dosage[i] < 0] = means[i]

    allele_freq = dosage.mean(axis=1) / 2
    maf = np.minimum(allele_freq, 1 - allele_freq)
    keep_maf = (maf >= MIN_MAF) & (maf < 0.5 + 1e-12)
    dosage = dosage[keep_maf]
    positions = positions[keep_maf]
    variant_ids = variant_ids[keep_maf]
    chroms = chroms[keep_maf]
    maf = maf[keep_maf]

    if dosage.shape[0] < 2:
        raise ValueError(
            f"Only {dosage.shape[0]} SNPs passed filtering; at least 2 are needed for PCA."
        )

    # Standardize SNPs by expected genotype standard deviation.
    p = dosage.mean(axis=1) / 2
    scale = np.sqrt(2 * p * (1 - p))
    usable = np.isfinite(scale) & (scale > 0)
    dosage = dosage[usable]
    positions = positions[usable]
    variant_ids = variant_ids[usable]
    chroms = chroms[usable]
    scale = scale[usable]
    p = p[usable]
    maf = maf[usable]
    standardized = (dosage - (2 * p)[:, None]) / scale[:, None]

    if standardized.shape[0] < 2:
        raise ValueError("Fewer than two informative SNPs remain after filtering.")

    pca = PCA(n_components=2)
    coordinates = pca.fit_transform(standardized.T)

    sample_metadata["PC1"] = coordinates[:, 0]
    sample_metadata["PC2"] = coordinates[:, 1]
    sample_metadata.to_csv(OUTPUT_DIR / "pca_sample_coordinates.csv", index=False)

    variant_table = pd.DataFrame({
        "chromosome": chroms,
        "position": positions,
        "variant_id": variant_ids,
        "minor_allele_frequency": maf,
    })
    variant_table.to_csv(OUTPUT_DIR / "pca_variants_used.csv", index=False)

    # Colour points by super-population to keep the legend readable.
    groups = sorted(sample_metadata["super_pop"].unique())
    cmap = plt.get_cmap("tab10", max(len(groups), 1))
    fig, ax = plt.subplots(figsize=(10, 7))
    for i, group in enumerate(groups):
        subset = sample_metadata[sample_metadata["super_pop"] == group]
        ax.scatter(subset["PC1"], subset["PC2"], s=16, alpha=0.75,
                   label=group, color=cmap(i), edgecolors="none")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)")
    ax.set_title("Exploratory PCA — chromosome 22 VCF subset")
    ax.legend(title="Super-population", frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "population_pca.png", dpi=300)
    plt.close(fig)

    print("PCA completed successfully.")
    print(f"Samples analysed: {len(sample_metadata)}")
    print(f"Biallelic SNPs initially selected: {len(keep_indices)}")
    print(f"SNPs used in PCA after filtering: {standardized.shape[0]}")
    print(f"PC1 variance explained: {pca.explained_variance_ratio_[0] * 100:.2f}%")
    print(f"PC2 variance explained: {pca.explained_variance_ratio_[1] * 100:.2f}%")
    print(f"Outputs saved to: {OUTPUT_DIR}")
    print("Note: this uses a limited subset of chromosome 22 and does not perform LD pruning.")


if __name__ == "__main__":
    main()
