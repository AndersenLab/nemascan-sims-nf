# Recombination-domain reference table

`recombination_domains.tsv` is the shipped per-species reference table that
`FILTER_CV_POOL_RESOLVE` looks region labels up in. It is the default value of
`--cv_region_domainfile`.

## Schema

One row per resolved interval (`arm_*` / `tip_*` labels span two rows that the
resolver unions; `center_*` is a single interval):

| column | description |
|---|---|
| `species` | pipeline species name (`c_elegans`, `c_briggsae`, `c_tropicalis`) |
| `filter_id` | region label listed in `--cv_region_filterfile` |
| `chrom` | chromosome in Roman form (`I`–`V`, `X`); normalized against the `.bim` at resolve time |
| `start` | interval start (bp, 1-based, inclusive), on that species' own assembly |
| `end` | interval end (bp, inclusive) |

`filter_id` labels follow the canonical recombination-domain convention:

- `tip_<chrom>` = left tip ∪ right tip
- `arm_<chrom>` = left arm ∪ right arm
- `center_<chrom>` = the single central low-recombination domain

The literal label `genome` is **not** in the table — it is a sentinel handled
inside the resolver (empty interval set ⇒ select the whole pool).

## Provenance

Domain boundaries derive from Rockman, M. V. & Kruglyak, L. (2009).
"Recombinational Landscape and Population Genomics of *Caenorhabditis
elegans*." *PLoS Genetics* 5(3): e1000419.

The per-species boundary pickles in `sources/` were supplied alongside this
work (`{ce,cb,ct}_chrom_dict.pkl`); each encodes the five domains
(`left_tip`, `left_arm`, `center`, `right_arm`, `right_tip`) per chromosome on
that species' own assembly. `build_recombination_domains.py` collapses them
into the three canonical labels above.

> **Coordinate caveat.** Each species' coordinates are on its own assembly. If
> the VCF a run maps against uses a different assembly, an interval can resolve
> to an empty marker pool (silently). `filter_pool_metrics.tsv` surfaces this
> as `pool_size = 0` / `status = starved` — check it before launching a sweep.

## Regenerate

```bash
python3 data/recombination_domains/build_recombination_domains.py
```
