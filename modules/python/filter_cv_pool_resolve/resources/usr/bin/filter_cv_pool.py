#!/usr/bin/env python3
"""
Resolve a causal-variant region filter against a PLINK .bim and emit the
extract list, the resolved-interval sidecar, the resolved-pool hash, and a
replication-aware per-nqtl metrics table.

Look a region label up in the species-appropriate recombination-domain
reference table, resolve it against a PLINK .bim, and emit:
  - extract.list             (file)  marker IDs in any resolved interval
  - resolved_intervals.json  (file)  canonical sorted intervals applied
  - pool_hash.txt            (file)  sha256 of the resolved marker list
  - filter_metrics.tsv       (file)  ONE ROW PER nqtl (replication-aware)

The trait hash is label-based and derived in R, not here -- this script
computes no trait digest. `pool_hash` is the single source for both the
per-cell selection seed (create_causal_vars.py) and the R-side
generate_causal_set_id().
"""
import argparse
import hashlib
import json
import math
import sys
from functools import reduce

import pandas as pd

# Canonical chromosome key: Roman <-> numeric, so a Roman reference table and a
# numerically-recoded .bim (I->1 ... V->5, X->6) match. Caenorhabditis species
# all have five autosomes plus X, so X maps to 6.
_ROMAN_TO_NUM = {"I": "1", "II": "2", "III": "3", "IV": "4", "V": "5", "X": "6"}


def normalize_chrom(chrom):
    """Return a canonical chromosome key shared by Roman and numeric forms."""
    s = str(chrom).strip().upper()
    if s.startswith("CHR"):
        s = s[3:]
    return _ROMAN_TO_NUM.get(s, s)


def lookup_region(domain_file, species, filter_id):
    """Resolve a region label to its interval set for one species.

    "genome" -> empty interval set (select-all). A label absent from the
    reference (for this species) is a hard error -- a mis-typed region must
    never silently select the whole pool.
    Returns a list of (chrom, start, end) tuples in the reference's own
    chromosome naming (normalized later, at marker matching).
    """
    if filter_id == "genome":
        return []

    ref = pd.read_csv(domain_file, sep="\t", comment="#")
    rows = ref[(ref["species"] == species) & (ref["filter_id"] == filter_id)]
    if rows.empty:
        sys.exit(
            f"Error: region label '{filter_id}' not found for species "
            f"'{species}' in {domain_file}. Check the label spelling and that "
            f"the reference table covers this species."
        )
    return [
        (str(r.chrom), int(r.start), int(r.end))
        for r in rows.itertuples(index=False)
    ]


def apply_filter(bim, intervals):
    """Select marker IDs whose POS falls in any interval (union).

    An empty interval set selects every marker (the "genome" / select-all
    case). Chromosomes are matched on the normalized key so Roman intervals
    match a numeric .bim.
    """
    if not intervals:
        return bim["marker"].tolist()

    norm_chrom = bim["CHROM"].map(normalize_chrom)
    masks = []
    for chrom, start, end in intervals:
        masks.append(
            (norm_chrom == normalize_chrom(chrom))
            & (bim["POS"] >= start)
            & (bim["POS"] <= end)
        )
    union = reduce(lambda a, b: a | b, masks)
    return bim.loc[union, "marker"].tolist()


def compute_rep_plan(N, k, reps, rep_start):
    """Replication-aware combinatorics for one (pool, nqtl).

    C(N, k) bounds the number of distinct causal sets drawable without
    replacement. effective_reps caps the requested replicate range to what the
    pool can actually supply.
    """
    C = math.comb(N, k) if N >= k else 0
    requested_end = rep_start + reps - 1
    effective_end = min(requested_end, C)
    effective_reps = max(0, effective_end - rep_start + 1)

    if C == 0:
        status = "starved"
    elif rep_start > C:
        status = "exhausted"
    elif effective_reps < reps:
        status = "reps_capped"
    else:
        status = "ok"

    n_combinations = str(C) if C <= 1e9 else ">1e9"
    return C, n_combinations, effective_reps, status


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bim", required=True)
    p.add_argument("--filter_id", required=True,
                   help='region label; "genome" => select-all')
    p.add_argument("--species", required=True)
    p.add_argument("--domain_file", required=True,
                   help="shipped per-species reference (species,filter_id,chrom,start,end)")
    p.add_argument("--nqtl_file", required=True, help="one nqtl value per line")
    p.add_argument("--reps", type=int, required=True)
    p.add_argument("--rep_start", type=int, required=True)
    p.add_argument("--group", required=True)
    p.add_argument("--ms_maf", required=True)
    p.add_argument("--extract_out", required=True)
    p.add_argument("--resolved_out", required=True)
    p.add_argument("--pool_hash_out", required=True)
    p.add_argument("--metrics_out", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    # 1-2. Resolve the label to intervals (in reference chrom naming).
    intervals = lookup_region(args.domain_file, args.species, args.filter_id)

    # 3. Read the .bim.
    bim = pd.read_csv(
        args.bim, sep="\t", header=None,
        names=["CHROM", "marker", "cm", "POS", "A1", "A2"],
        dtype={"CHROM": str},
    )

    # 4-5. Mask + union -> selected marker IDs. N = pool size.
    selected = apply_filter(bim, intervals)
    N = len(selected)
    with open(args.extract_out, "w") as fh:
        fh.write("\n".join(selected))
        if selected:
            fh.write("\n")

    # 6. Resolved-interval sidecar: canonical sorted intervals ([] for genome).
    sidecar = sorted(
        ({"chrom": c, "start": s, "end": e} for c, s, e in intervals),
        key=lambda r: (normalize_chrom(r["chrom"]), r["start"], r["end"]),
    )
    with open(args.resolved_out, "w") as fh:
        json.dump(sidecar, fh)

    # 6b. pool_hash = sha256 of the resolved marker list -- single source for
    #     cell_seed (create_causal_vars.py) and R's generate_causal_set_id().
    pool_hash = hashlib.sha256("\n".join(selected).encode()).hexdigest()
    with open(args.pool_hash_out, "w") as fh:
        fh.write(pool_hash + "\n")

    # 7-8. One metrics row per nqtl (11-column schema).
    with open(args.nqtl_file) as fh:
        nqtls = [int(line.strip()) for line in fh if line.strip()]

    header = ["group", "ms_maf", "filter_id", "species", "nqtl", "pool_size",
              "n_combinations", "rep_start", "requested_reps", "effective_reps",
              "status"]
    with open(args.metrics_out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for k in nqtls:
            _, n_comb, eff_reps, status = compute_rep_plan(
                N, k, args.reps, args.rep_start)
            row = [args.group, args.ms_maf, args.filter_id, args.species,
                   str(k), str(N), n_comb, str(args.rep_start),
                   str(args.reps), str(eff_reps), status]
            fh.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
