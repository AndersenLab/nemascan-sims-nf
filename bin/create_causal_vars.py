#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys


def load_strain_set_variants(strain_set_variant_file):
    """
    Load the strain set variants from a plink .bim file
    """

    # just read in the marker column from the .bim file
    strain_set_variants = pd.read_csv(
        strain_set_variant_file,
        sep="\t",
        header=None,
        usecols=[1],
        names=["Marker"],
    )

    return strain_set_variants


def select_variants(markers, n_var, rng):
    """
    Select n_var variants from the list of markers
    """
    # Select n_var variants from the list of markers
    # Ensure markers is a 1D array by extracting the "Marker" column if it's a DataFrame
    if isinstance(markers, pd.DataFrame):
        markers = markers["Marker"].values

    selected_variants = rng.choice(markers, n_var, replace=False)
    return selected_variants


def simulate_effect_gamma(
    selected_variants, n_var, rng, effect_shape=0.4, effect_scale=1.66
):
    # pull the effect size from a gamma distribution
    effects = rng.gamma(effect_shape, effect_scale, n_var)

    directions = rng.choice([-1, 1], n_var)
    effects = effects * directions

    # Create a DataFrame from the selected variants array and effects
    result_df = pd.DataFrame({"Marker": selected_variants, "EFFECT": effects})

    return result_df


def simulate_effect_uniform(selected_variants, n_var, rng, low_end=0.1, high_end=0.5):
    # pull the effect size from a uniform distribution
    effects = rng.uniform(low_end, high_end, n_var)

    directions = rng.choice([-1, 1], n_var)
    effects = effects * directions

    # Create a DataFrame from the selected variants array and effects
    result_df = pd.DataFrame({"Marker": selected_variants, "EFFECT": effects})

    return result_df


def read_bed_genotypes(bed_file, bim_file, fam_file, selected_variant_ids):
    """
    Read per-strain genotypes at selected variant positions from PLINK .bed.

    Returns a long-format DataFrame: QTL (CHROM:POS), strain, allele
    Allele encoding: -1.0 (hom ref / A1), 1.0 (hom alt / A2), NaN (het or missing)
    Matches the encoding convention of the marker set genotype matrix.
    """
    fam = pd.read_csv(fam_file, sep=r'\s+', header=None,
                      names=["fid", "iid", "pid", "mid", "sex", "pheno"])
    bim = pd.read_csv(bim_file, sep='\t', header=None,
                      names=["chrom", "marker", "cm", "pos", "a1", "a2"])

    n_samples = len(fam)
    bytes_per_variant = (n_samples + 3) // 4

    with open(bed_file, 'rb') as f:
        header = f.read(3)
        assert header == b'\x6c\x1b\x01', "Not a valid variant-major PLINK .bed"
        raw = np.frombuffer(f.read(), dtype=np.uint8)

    selected_mask = bim['marker'].isin(selected_variant_ids)
    selected_bim = bim[selected_mask].reset_index(drop=True)
    selected_idx = np.where(selected_mask.values)[0]

    rows = []
    for bim_idx, row in zip(selected_idx, selected_bim.itertuples()):
        qtl = f"{row.chrom}:{row.pos}"
        start = bim_idx * bytes_per_variant
        variant_bytes = raw[start: start + bytes_per_variant]
        # Unpack 2 bits per sample, LSB-first within each byte (PLINK variant-major)
        genos = np.unpackbits(variant_bytes, bitorder='little').reshape(-1, 2)[:n_samples]
        low, high = genos[:, 0], genos[:, 1]
        # 00=hom A1=-1, 11=hom A2=+1, else (01=missing, 10=het)=NaN
        allele = np.where((low == 0) & (high == 0), -1.0,
                 np.where((low == 1) & (high == 1),  1.0, np.nan))
        for strain, a in zip(fam['iid'].values, allele):
            rows.append({"QTL": qtl, "strain": strain, "allele": a})

    return pd.DataFrame(rows)


if __name__ == "__main__":
    # Get command line arguments
    strain_set_variant_file = sys.argv[1]
    n_var = int(sys.argv[2])

    if sys.argv[3] == "gamma":
        effect_type = "gamma"
    elif "-" in sys.argv[3]:
        effect_type = "uniform"
        effect_range = sys.argv[3].split("-")
        low_end = float(effect_range[0])
        high_end = float(effect_range[1])
    else:
        print(
            f"Error: Invalid effect specification '{sys.argv[3]}'. Use 'gamma' or 'low_end-high_end'"
        )
        sys.exit(1)
    rep = int(sys.argv[4])
    rng = np.random.default_rng()
    # Read in annotated strain_set variants
    strain_var = load_strain_set_variants(strain_set_variant_file)

    # Pool size guard: fail fast with an informative message rather than a cryptic
    # numpy traceback when np.random.choice receives n_var > pool size.
    n_pool = len(strain_var)
    if n_var > n_pool:
        print(
            f"Error: requested nqtl={n_var} but CV pool contains only {n_pool} variants. "
            f"Reduce --nqtl or relax --cv_maf / --cv_ld thresholds."
        )
        sys.exit(1)

    # Select Causal variants for orthogroups
    causal_vars = select_variants(strain_var, n_var, rng)
    print("Selected Causal Variants")

    if effect_type == "gamma":
        print("Simulating effects using gamma distribution")
        causal_vars_effects = simulate_effect_gamma(causal_vars, n_var, rng)
    else:
        print(f"Simulating effects using uniform distribution: {low_end}-{high_end}")
        causal_vars_effects = simulate_effect_uniform(
            causal_vars, n_var, rng, low_end=low_end, high_end=high_end
        )
    print("Simulated effects")

    # Write output for trait simulations - just id and effect
    causal_vars_effects.to_csv(
        "causal_vars.txt",
        sep=" ",
        index=False,
        header=False,
    )

    # Extract per-strain genotypes at the selected causal positions from CV_TO_SIMS.bed
    bed_prefix = strain_set_variant_file.replace(".bim", "")
    causal_geno_df = read_bed_genotypes(
        bed_prefix + ".bed",
        strain_set_variant_file,
        bed_prefix + ".fam",
        causal_vars_effects["Marker"].values,
    )
    causal_geno_df.to_csv(
        f"causal_genotypes.sim.{n_var}.{rep}.tsv",
        sep="\t",
        index=False,
        header=True,
        na_rep="NA",
    )
