#!/usr/bin/env python3
import hashlib
import json
import math
import sys

import numpy as np
import pandas as pd


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


def unrank_combination(pos, N, k):
    """
    Lexicographic combination unranking.

    Return the ``pos``-th k-subset of ``range(N)`` (0-based rank), as a list of
    k indices sorted ascending. ``pos`` must be in ``[0, C(N, k))``. O(N * k).
    """
    result = []
    x = 0
    remaining = pos
    for i in range(k):
        # Greedily place the i-th (0-based) chosen index. C(N-1-x, k-1-i) counts
        # the combinations whose i-th element is exactly x; skip past blocks that
        # rank below `remaining`.
        while True:
            c = math.comb(N - 1 - x, k - 1 - i)
            if remaining < c:
                break
            remaining -= c
            x += 1
        result.append(x)
        x += 1
    return result


def prp(index, C, seed):
    """
    Stable seeded pseudo-random permutation (bijection) on ``[0, C)``.

    Cycle-walking balanced Feistel network sized to ceil(log2 C) bits (rounded
    up to an even width). A balanced Feistel is a permutation for any round
    function, and cycle-walking (re-encrypting until the image lands in
    ``[0, C)``) preserves bijectivity on the restricted domain. ``O(1)`` per
    call in expectation. The seed makes a capped subset of reps an unbiased
    sample of the combination space rather than the lexicographically-first
    cluster.
    """
    if C <= 1:
        return index

    total_bits = (C - 1).bit_length()
    if total_bits < 2:
        total_bits = 2
    if total_bits % 2:
        total_bits += 1
    half = total_bits // 2
    mask = (1 << half) - 1
    rounds = 4

    def feistel(x):
        left = (x >> half) & mask
        right = x & mask
        for r in range(rounds):
            h = hashlib.sha256(f"{seed}|{r}|{right}".encode()).digest()
            f = int.from_bytes(h[:8], "big") & mask
            left, right = right, left ^ f
        return (left << half) | right

    x = index
    while True:
        x = feistel(x)
        if x < C:
            return x


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
    # Positional args:
    #   bim  nqtl  effect(="gamma")  rep  filter_id  pool_hash
    strain_set_variant_file = sys.argv[1]
    n_var = int(sys.argv[2])

    # Effect distribution is fixed to gamma in production; the uniform branch is
    # retained so the script stays usable standalone / under tests.
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

    # filter_id (region label) and pool_hash (the deterministic-selection basis)
    # are supplied by the module in production. Defaulted here so the script is
    # self-contained for standalone / test invocation.
    filter_id = sys.argv[5] if len(sys.argv) > 5 else "genome"
    pool_hash = sys.argv[6] if len(sys.argv) > 6 else None

    # Read in annotated strain_set variants (in .bim row order — the unranking
    # is against this order, which PLINK --extract preserves).
    strain_var = load_strain_set_variants(strain_set_variant_file)
    markers = strain_var["Marker"].tolist()

    N = len(markers)
    k = n_var
    C = math.comb(N, k) if N >= k else 0
    index = rep - 1                          # 0-based position in the cell ordering

    # Safety net behind the upstream fanout cap: the cap should keep impossible
    # cells from ever being scheduled, so reaching this is a guard, not the
    # normal path. Structured JSON to stderr feeds the failure record.
    if C == 0 or index >= C:
        err = {
            "error": "pool_starvation",
            "filter_id": filter_id,
            "pool_size": int(N),
            "nqtl_requested": int(k),
            "rep": int(rep),
        }
        print(json.dumps(err), file=sys.stderr)
        sys.exit(1)

    # pool_hash is the single source for the selection seed. When not supplied,
    # fall back to the sha of the resolved marker list — the same definition the
    # upstream resolver uses, so the value never diverges from the
    # channel-supplied one in production.
    if not pool_hash:
        pool_hash = hashlib.sha256("\n".join(markers).encode()).hexdigest()

    # Cell seed EXCLUDES rep, effect, and h2: all reps of a (region, nqtl) cell
    # share one ordering (made distinct per rep by the PRP), identical across
    # executions (so a later --rep_start run continues the same sequence without
    # collisions) AND across h2 values (the causal set is shared; h2 only scales
    # the downstream simulation). causal_set_id is intentionally not computed
    # here — it is derived in R from this same pool_hash at trait-write time.
    cell_seed = hashlib.sha256(f"{pool_hash}|nqtl={k}".encode()).hexdigest()

    pos = prp(index, C, cell_seed)           # seeded bijection on [0, C)
    combo = unrank_combination(pos, N, k)    # k 0-based marker indices, lexicographic
    selected = [markers[i] for i in combo]
    print("Selected Causal Variants")

    # Effects: gamma-drawn, reproducible per (cell, rep) — seeded by (cell_seed,
    # rep), NOT per h2.
    eff_rng = np.random.default_rng(int(cell_seed[:16], 16) ^ rep)

    if effect_type == "gamma":
        print("Simulating effects using gamma distribution")
        causal_vars_effects = simulate_effect_gamma(selected, k, eff_rng)
    else:
        print(f"Simulating effects using uniform distribution: {low_end}-{high_end}")
        causal_vars_effects = simulate_effect_uniform(
            selected, k, eff_rng, low_end=low_end, high_end=high_end
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
