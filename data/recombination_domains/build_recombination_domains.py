#!/usr/bin/env python3
"""
Build the shipped recombination-domain reference table from the per-species
chromosome-domain pickles in sources/.

Source: Rockman, M. V. & Kruglyak, L. (2009). "Recombinational Landscape and
Population Genomics of Caenorhabditis elegans." PLoS Genetics 5(3): e1000419.
The per-species pickles in sources/ encode each chromosome's five recombination
domains (left_tip, left_arm, center, right_arm, right_tip) as Start/Stop
coordinates on each species' own assembly.

This script collapses the five domains into the three canonical filter labels
used by FILTER_CV_POOL_RESOLVE (per example_filters.md):

    tip_<chrom>     = left_tip  U right_tip   (two intervals, unioned at resolve time)
    arm_<chrom>     = left_arm  U right_arm   (two intervals, unioned at resolve time)
    center_<chrom>  = center                  (one interval)

Output columns: species, filter_id, chrom, start, end (one row per interval;
arm_/tip_ labels span two rows). The resolver unions all rows sharing a
(species, filter_id).

Regenerate with:  python3 data/recombination_domains/build_recombination_domains.py
"""
import pickle
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = HERE / "sources"
OUT = HERE / "recombination_domains.tsv"

# pickle file -> pipeline species name
SPECIES = {
    "ce_chrom_dict.pkl": "c_elegans",
    "cb_chrom_dict.pkl": "c_briggsae",
    "ct_chrom_dict.pkl": "c_tropicalis",
}

CHROM_ORDER = ["I", "II", "III", "IV", "V", "X"]

# filter_id label -> the source domains it unions, in emit order
LABELS = [
    ("tip_{c}",    ["left_tip", "right_tip"]),
    ("arm_{c}",    ["left_arm", "right_arm"]),
    ("center_{c}", ["center"]),
]


def main():
    rows = [("species", "filter_id", "chrom", "start", "end")]
    for fname in ("ce_chrom_dict.pkl", "cb_chrom_dict.pkl", "ct_chrom_dict.pkl"):
        species = SPECIES[fname]
        with open(SOURCES / fname, "rb") as fh:
            d = pickle.load(fh)
        for chrom in CHROM_ORDER:
            for label_tmpl, domains in LABELS:
                filter_id = label_tmpl.format(c=chrom)
                for dom in domains:
                    seg = d[chrom][dom]
                    rows.append((species, filter_id, chrom,
                                 str(seg["Start"]), str(seg["Stop"])))

    with open(OUT, "w") as out:
        for r in rows:
            out.write("\t".join(r) + "\n")
    print(f"wrote {OUT} ({len(rows) - 1} interval rows across "
          f"{len(SPECIES)} species)")


if __name__ == "__main__":
    main()
