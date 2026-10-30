"""Calculate chromosome 22 ALT-allele frequencies across 1000 Genomes populations.

Inputs (in ../data/):
    chr22_first1000.vcf
    integrated_call_samples_v3.20130502.ALL.panel

Outputs are written to ../final/. This is an analysis of the supplied VCF subset,
not a genome-wide analysis.
"""

from pathlib import Path

import allel
import matplotlib.pyplot as plt
import pandas as pd


# ---------------------------------------------------------------------------
# Paths and configuration
# ---------------------------------------------------------------------------

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"
FINAL_DIR = PROJECT_DIR / "final"

VCF_PATH = DATA_DIR / "chr22_first1000.vcf"
PANEL_PATH = DATA_DIR / "integrated_call_samples_v3.20130502.ALL.panel"

TARGET_VARIANT_ID = "rs529371290"
TARGET_VARIANT_POSITION = 16134585

POPULATION_RESULTS_PATH = FINAL_DIR / "population_variant_frequencies.csv"
SUMMARY_RESULTS_PATH = FINAL_DIR / "variant_frequency_summary.csv"
TARGET_PLOT_PATH = FINAL_DIR / "rs529371290_population_frequency.png"
LOWEST_VARIANTS_PLOT_PATH = FINAL_DIR / "lowest_10_variant_frequencies.png"


# ---------------------------------------------------------------------------
# Input loading and validation
# ---------------------------------------------------------------------------

def decode_text(value) -> str:
    """Convert VCF byte strings or other values to text."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def read_vcf_sample_ids(vcf_path: Path) -> list[str]:
    """Read sample IDs from the VCF #CHROM header."""
    with vcf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return line.rstrip("\r\n").split("\t")[9:]
    raise ValueError(f"No #CHROM header found in VCF: {vcf_path}")


def load_inputs():
    """Load panel and VCF data and validate sample alignment."""
    if not VCF_PATH.is_file():
        raise FileNotFoundError(f"VCF file not found: {VCF_PATH}")
    if not PANEL_PATH.is_file():
        raise FileNotFoundError(f"Population panel not found: {PANEL_PATH}")

    panel = pd.read_csv(PANEL_PATH, sep=r"\s+")
    required_columns = {"sample", "pop"}
    missing_columns = required_columns - set(panel.columns)
    if missing_columns:
        raise ValueError(
            f"Panel is missing required columns: {sorted(missing_columns)}"
        )
    if panel["sample"].duplicated().any():
        duplicates = panel.loc[panel["sample"].duplicated(), "sample"].tolist()
        raise ValueError(f"Duplicate sample IDs in panel: {duplicates[:10]}")

    vcf_samples = read_vcf_sample_ids(VCF_PATH)
    if len(vcf_samples) != len(set(vcf_samples)):
        raise ValueError("The VCF header contains duplicate sample IDs.")

    panel_sample_ids = set(panel["sample"].astype(str))
    missing_samples = sorted(set(vcf_samples) - panel_sample_ids)
    if missing_samples:
        raise ValueError(
            f"{len(missing_samples)} VCF samples were not found in the panel. "
            f"Examples: {missing_samples[:10]}"
        )

    callset = allel.read_vcf(
        str(VCF_PATH),
        fields=[
            "variants/CHROM",
            "variants/POS",
            "variants/ID",
            "variants/REF",
            "variants/ALT",
            "calldata/GT",
        ],
    )
    if callset is None or "calldata/GT" not in callset:
        raise ValueError("Could not load genotype data from the VCF.")

    genotypes = allel.GenotypeArray(callset["calldata/GT"])
    if genotypes.shape[1] != len(vcf_samples):
        raise ValueError(
            f"VCF header has {len(vcf_samples)} samples, but genotype data "
            f"has {genotypes.shape[1]} samples."
        )

    sample_to_population = panel.set_index("sample")["pop"]
    sample_populations = sample_to_population.loc[vcf_samples].to_numpy()
    populations = sorted(panel["pop"].astype(str).unique().tolist())

    print(f"VCF: {VCF_PATH}")
    print(f"Panel: {PANEL_PATH}")
    print(f"Samples in VCF: {len(vcf_samples)}")
    print(f"Samples matched to panel: {len(vcf_samples)}")
    print(f"Populations represented in panel: {len(populations)}")
    print(f"Variant records loaded: {genotypes.shape[0]}")

    return panel, callset, genotypes, populations, sample_populations


# ---------------------------------------------------------------------------
# Variant processing
# ---------------------------------------------------------------------------

def extract_single_alt(alt_values) -> str | None:
    """Return the sole non-missing ALT allele, or None for other records."""
    alleles = []
    for value in alt_values:
        allele = decode_text(value).strip()
        if allele and allele != ".":
            alleles.append(allele)
    return alleles[0] if len(alleles) == 1 else None


def calculate_population_frequencies(
    callset: dict,
    genotypes: allel.GenotypeArray,
    populations: list[str],
    sample_populations,
) -> pd.DataFrame:
    """Calculate REF/ALT counts and ALT frequency per variant and population."""
    indices_by_population = {
        population: [
            i for i, sample_population in enumerate(sample_populations)
            if sample_population == population
        ]
        for population in populations
    }

    rows = []
    for variant_index in range(genotypes.shape[0]):
        alt = extract_single_alt(callset["variants/ALT"][variant_index])
        if alt is None:
            continue

        ref = decode_text(callset["variants/REF"][variant_index])
        chromosome = decode_text(callset["variants/CHROM"][variant_index])
        position = int(callset["variants/POS"][variant_index])
        variant_id = decode_text(callset["variants/ID"][variant_index])

        for population, sample_indices in indices_by_population.items():
            population_genotypes = genotypes[
                variant_index:variant_index + 1, sample_indices, :
            ]
            counts = population_genotypes.count_alleles(max_allele=1)[0]
            ref_count = int(counts[0])
            alt_count = int(counts[1])
            total_count = ref_count + alt_count
            frequency = alt_count / total_count if total_count else float("nan")

            rows.append({
                "chromosome": chromosome,
                "position": position,
                "variant_id": variant_id,
                "ref": ref,
                "alt": alt,
                "population": population,
                "ref_allele_count": ref_count,
                "alt_allele_count": alt_count,
                "total_alleles": total_count,
                "alt_frequency": frequency,
            })

    if not rows:
        raise ValueError("No biallelic variant records were found in the VCF.")

    results = pd.DataFrame(rows)
    if not results["alt_frequency"].dropna().between(0, 1).all():
        raise ValueError("Some calculated ALT frequencies are outside [0, 1].")
    if not (
        results["ref_allele_count"] + results["alt_allele_count"]
        == results["total_alleles"]
    ).all():
        raise ValueError("Allele-count consistency check failed.")

    return results


def build_variant_summary(population_results: pd.DataFrame) -> pd.DataFrame:
    """Aggregate allele counts across populations for each distinct variant."""
    keys = ["chromosome", "position", "variant_id", "ref", "alt"]
    summary = (
        population_results.groupby(keys, dropna=False)[
            ["alt_allele_count", "total_alleles"]
        ]
        .sum()
        .reset_index()
    )
    summary["overall_alt_frequency"] = (
        summary["alt_allele_count"] / summary["total_alleles"]
    )
    return summary


# ---------------------------------------------------------------------------
# Output tables and plots
# ---------------------------------------------------------------------------

def save_target_variant_plot(population_results: pd.DataFrame) -> None:
    """Save population ALT frequencies for the selected variant."""
    data = population_results[
        (population_results["position"] == TARGET_VARIANT_POSITION)
        & (population_results["variant_id"] == TARGET_VARIANT_ID)
    ].copy()

    if data.empty:
        print(
            f"Target variant {TARGET_VARIANT_ID} at position "
            f"{TARGET_VARIANT_POSITION} not found; skipping its plot."
        )
        return

    data = data.sort_values("alt_frequency")
    plt.figure(figsize=(12, 6))
    plt.bar(data["population"], data["alt_frequency"])
    plt.title(
        "ALT allele frequency by population\n"
        f"Chromosome 22, position {TARGET_VARIANT_POSITION} "
        f"({TARGET_VARIANT_ID})"
    )
    plt.xlabel("Population")
    plt.ylabel("ALT allele frequency")
    plt.ylim(0, 1.05)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(TARGET_PLOT_PATH, dpi=300)
    plt.close()
    print(f"Saved plot: {TARGET_PLOT_PATH}")


def save_lowest_variants_plot(summary: pd.DataFrame) -> None:
    """Save a plot of the ten records with the lowest overall ALT frequency."""
    plot_data = summary.nsmallest(10, "overall_alt_frequency").copy()
    plot_data["label"] = (
        plot_data["variant_id"].astype(str)
        + " ("
        + plot_data["position"].astype(str)
        + ")"
    )

    plt.figure(figsize=(12, 6))
    plt.bar(plot_data["label"], plot_data["overall_alt_frequency"])
    plt.title("10 variant records with the lowest overall ALT frequency")
    plt.xlabel("Variant (ID and chromosome 22 position)")
    plt.ylabel("Overall ALT allele frequency")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(LOWEST_VARIANTS_PLOT_PATH, dpi=300)
    plt.close()
    print(f"Saved plot: {LOWEST_VARIANTS_PLOT_PATH}")


def print_summary(
    population_results: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    """Print concise result and validation summaries."""
    variant_keys = ["chromosome", "position", "variant_id", "ref", "alt"]
    print("\n--- Analysis summary ---")
    print(
        "Distinct biallelic variant records:",
        summary[variant_keys].drop_duplicates().shape[0],
    )
    print("Rows in population table:", len(population_results))
    print("Populations represented:", population_results["population"].nunique())
    print("Rows with missing frequency:", population_results["alt_frequency"].isna().sum())

    print("\nFive lowest overall ALT frequencies:")
    print(summary.nsmallest(5, "overall_alt_frequency").to_string(index=False))
    print("\nFive highest overall ALT frequencies:")
    print(summary.nlargest(5, "overall_alt_frequency").to_string(index=False))


def main() -> None:
    """Run the population-level variant-frequency analysis."""
    FINAL_DIR.mkdir(parents=True, exist_ok=True)

    _, callset, genotypes, populations, sample_populations = load_inputs()
    population_results = calculate_population_frequencies(
        callset, genotypes, populations, sample_populations
    )
    population_results.to_csv(POPULATION_RESULTS_PATH, index=False)
    print(f"Saved population results: {POPULATION_RESULTS_PATH}")

    summary = build_variant_summary(population_results)
    summary.to_csv(SUMMARY_RESULTS_PATH, index=False)
    print(f"Saved variant summary: {SUMMARY_RESULTS_PATH}")

    print_summary(population_results, summary)
    save_target_variant_plot(population_results)
    save_lowest_variants_plot(summary)


if __name__ == "__main__":
    main()
