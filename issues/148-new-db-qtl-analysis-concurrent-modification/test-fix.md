# 2026-04-08: Test fix for concurrent modification in QTL analysis

The fix replaces independent wrap-combine-unwrap channel construction with a single
  merge-before-expansion pattern using .merge() + .multiMap(). This affects both the DB
  analysis path (ch_db_analysis) and the legacy assessment path (ch_intervals). Each
  emission now owns its own tuple — no shared ArrayList references.

## 1. Local stub-run

Verify that the channel wiring is working

```bash
NXF_VER=24.10.5 nextflow run main.nf -profile test -stub-run
```

**Suceeded**

Run stub check with the legacy path which was also impacted by the change

```bash
NXF_VER=24.10.5 nextflow run main.nf -profile test -stub-run --legacy_assess
```

**Suceeded**

## 2. Local full run

### Fixed architecture test (1 pop, 1 rep)

```bash
NXF_VER=24.10.5 nextflow run main.nf -profile test,docker
```

**Suceeded**

### Fixed architecture test (1 pop, 1 rep)

```bash
NXF_VER=24.10.5 nextflow run main.nf -profile test,docker --legacy_assess
```

**Suceeded**

## 3. HPC resume

```bash
nextflow run main.nf \
  -profile rockfish \
  --strainfile data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/strains_three_species.tsv \
  --nqtl data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/nqtl.csv \
  --h2 data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/h2.csv \
  --effect data/sims_ce-cb-ct_nqtl1-50_h2grid_50reps/effect_sizes.csv \
  --reps 25 \
  --cv_maf 0.05 \
  --cv_ld 0.8 \
  --output_dir Sims_ce-cb-ct_nqtl1-50_h2grid_50reps \
  -work-dir /scratch4/eande106/Ryan/nf-work-sims-ce-cb-ct5 \
  -resume
```

Started at 2026-04-08 13:09 

Suceeded!