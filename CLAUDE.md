# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NemaScan-sims-nf is a Nextflow pipeline for GWAS simulation in *C. elegans*. It generates synthetic phenotypes with known causal variants, runs GCTA-based genome-wide association mapping, detects QTL, and assesses detection power/precision against ground truth.

## Common Commands

### Running the pipeline
```bash
# Local test (fixed architecture, 13 strains, 3 chromosomes)
nextflow run main.nf -profile test

# Local test with Docker containers
nextflow run main.nf -profile test,docker

# Variable architecture (multi-parameter grid)
nextflow run main.nf -profile test_variable

# With legacy R-based assessment for cross-validation
nextflow run main.nf -profile test --legacy_assess

# Stub-run (wiring check, no computation)
nextflow run main.nf -profile test -stub-run

# HPC (Rockfish SLURM cluster with Singularity)
nextflow run main.nf -profile rockfish --strainfile strains.tsv --vcf 20231213
```

### Testing
```bash
# Run R unit + integration tests (from project root)
Rscript tests/run_tests.R

# With integration test data (after running collect_test_data.sh)
TEST_DB_DIR=tests/integration_data/db \
TEST_WORK_DIR=tests/.nf-work \
TEST_EXISTING_ASSESSMENT=tests/integration_data/simulation_assessment_results.tsv \
TEST_DB_ASSESSMENT=tests/integration_data/db_simulation_assessment_results.tsv \
Rscript tests/run_tests.R

# Collect integration test data from a pipeline run
bash tests/collect_test_data.sh              # run pipeline + collect
bash tests/collect_test_data.sh --collect-only  # collect from previous run
```

### Documentation
```bash
quarto render docs/   # renders to docs/_site/

# Publish to GitHub Pages (works from any branch)
quarto publish gh-pages docs/
```

## Architecture

### Pipeline Phases (main.nf)

1. **Prepare Marker Set**: VCF → strain extraction (bcftools) → LD-pruned PLINK binaries → genotype matrix eigenvalues
2. **Simulate Phenotypes**: Select random causal variants (Python) → GCTA phenotype simulation → variance check/upscaling
3. **GWAS Mapping**: GRM construction → association mapping in 4 mode×type combinations:
   - Mode: `inbred` (fastGWA-mlm-exact) vs `loco` (mlma-loco)
   - Type: `pca` (PC1 covariate) vs `nopca` (no covariate)
4. **Database Integration** (default): Write marker sets, long-format genotype matrices, GWA results, and per-trait data (metadata, causal variants, pre-upscaled phenotype) to Parquet DB → QTL detection → assessment
5. **Legacy Assessment** (optional, `--legacy_assess`): R-based QTL intervals → assessment (for cross-validation with DB path)

### Parallelization

The pipeline creates a cartesian product over: strain sets × MAF thresholds × nQTL × effect sizes × h2 values × replicates × mapping mode × mapping type. Each combination runs independently.

### Key Directories

- `modules/` — Nextflow process definitions, organized by tool (bcftools, plink, gcta, python, r, db_migration). Each module has `main.nf` and optionally `resources/` with scripts.
- `bin/` — Pipeline scripts (Python/R) available to all processes via `moduleBinaries = true`
- `R/` — Shared R library (9 files). `setup.R` loads 7 in dependency order:
  ```
  utils.R → io.R → database.R → queries.R → analysis.R → qtl_database.R → sim_performance.R
  ```
  `assessment.R` is **not** loaded by `setup.R` — it's sourced directly by `db_migration/assess_sims` and the legacy `Assess_Sims.R` scripts.
- `conf/` — Execution profiles: `rockfish.config` (SLURM/Singularity), `docker.config`
- `data/test/` — Test VCF (~825MB, 13 strains, chr I/II/V), parameter files, strain lists
- `tests/` — R/testthat test suite with fixtures and integration data

### Module Resources Pattern

Database migration modules (`modules/db_migration/`) use a two-layer pattern:
1. Each module's `main.nf` sets `R_SOURCE_DIR="${projectDir}/R"` as an environment variable
2. The executable R script lives in `resources/usr/bin/` and is auto-added to `$PATH` via `moduleBinaries = true`
3. The resource script uses `R_SOURCE_DIR` to source only the R library files it needs (not all of them)
4. All 7 `db_migration` processes emit `path "versions.yml", emit: versions` (R version captured via `Rscript --version`), consistent with the `modules/r/` module pattern. Stub blocks use `R: stub`. Do NOT add `when:` or `task.ext.args` to db_migration modules — they are intentionally absent; these modules use direct Nextflow parameter passing, not the `task.ext` config pattern.

Legacy modules in `modules/r/` use `bin/` scripts instead, which are available to all processes.

### Mapping Filename Convention

Mapping output files follow a strict naming pattern parsed by `parse_mapping_filename()` in `R/utils.R`:
```
{nqtl}_{rep}_{h2}_{maf}_{effect}_{population}_processed_LMM-EXACT-{INBRED|LOCO}[_PCA]_mapping.tsv
```
Example: `5_1_0.2_0.05_gamma_ct.fullpop.20210901_processed_LMM-EXACT-INBRED_PCA_mapping.tsv`

### Configuration

Requires Nextflow `>=24.10.0,<25.0.0`. Uses `nextflow.preview.output = true` for publish blocks and `nextflow.enable.moduleBinaries = true` so module `resources/` directories are added to `$PATH`.

### Output Structure

Primary output is `db_simulation_assessment_results.tsv`. The Parquet database lives under `{outputDir}/db/` with six directories:

| Directory | Content |
|-----------|---------|
| `marker_sets/` | `{pop}_{maf}_markers.parquet` + `{pop}_{maf}_genotypes.parquet` (long format) + `marker_set_metadata.parquet` |
| `mappings/` | Hive-partitioned GWA statistics: `population={pop}/mapping_id={id}/data.parquet` |
| `traits/` | Per-trait metadata: `{trait_id}.parquet` |
| `causal_variants/` | Causal variant parameters: `{trait_id}_causal.parquet` |
| `phenotypes/` | Pre-upscaled phenotype values: `{trait_id}_phenotype.parquet` |
| *(root)* | `mappings_metadata.parquet` |

`var.exp` was removed from the mappings schema in Phase 5 — it was always `NA` in the DB path. Raw data for offline variance-explained estimation is now in `traits/`, `causal_variants/`, and `phenotypes/`.

## Key Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--strainfile` | TSV with strain set names and lists (required) | — |
| `--vcf` | VCF path or CaeNDR release date (required) | — |
| `--reps` | Simulation replicates per parameter set | 2 |
| `--sthresh` | Significance threshold: "BF", "EIGEN", or numeric | "BF" |
| `--legacy_assess` | Also run legacy R assessment chain | false |
| `--group_qtl` | Distance (bp) for grouping nearby QTL | 1000 |
| `--ci_size` | SNVs flanking peak for CI definition | 150 |
| `--alpha` | Significance level for threshold calculation | 0.05 |
| `--db_output` | Parquet DB directory (absolute path) | `{outputDir}/db` |

## Documentation Site

The `docs/` directory is a Quarto website published to GitHub Pages. One `.qmd` page per pipeline phase plus `development.qmd` (testing/publishing guide) and `cross-validation.qmd` (legacy vs DB concordance). Publishing: `quarto publish gh-pages docs/` works from any branch.

## Development Notes

- The notes vault at `nemascan-sims-nf-notes/` (Obsidian, may be added as working directory) contains development plans, testing records, and technical notes. Its own `CLAUDE.md` describes navigation.
- Test VCF generation: `bash data/test/generate_test_vcf.sh /path/to/source.vcf.gz`
- GWA output formats differ by mode: `.fastGWA` (inbred) vs `.mlma` (loco) — `read_raw_gwa_file()` in `R/io.R` auto-detects format and normalizes column names.
- GCTA phenotype simulation can produce very small variance; `check_vp.py` upscales by 1000× if Vp < 1e-6.
- Log10P values are computed at query time via `safe_log10p()` in `R/analysis.R`, never stored in the Parquet database.
- `PYTHON_CHECK_VP` is defined in `modules/python/check_vp.nf` (flat file, not in a subdirectory like other modules).
- The test suite has two tiers: **unit tests** use static fixtures in `tests/fixtures/` and need no env vars; **integration tests** skip automatically when `TEST_DB_DIR` etc. are not set.
- **Phase 5 — Trait IDs and Mapping IDs:**
  - `generate_trait_id(group, maf, nqtl, effect, rep, h2)` → `MD5(paste(...))[1:12]` — 12-char hex string, stable across `-resume`, identifies all trait data files
  - `generate_mapping_id(trait_id, mode, type)` → `MD5(paste(...))[1:12]` — identifies a specific mapping run; used as Hive partition key in `mappings/population={pop}/mapping_id={id}/`
  - Both computed only in R — no Groovy counterpart.
- **Phase 5 — `var.exp` removed:** The `var.exp` column is gone from `mappings_schema()` and `prepare_mapping_data()`. `R/queries.R` already handled its absence (`NULL AS "var.exp"` fallback) and `R/assessment.R` uses `any_of("var.exp")` — no query-side changes needed.
- **Phase 5 — Parquet write safety:** All DB write functions default to `overwrite = TRUE` because `arrow::write_parquet()` is non-atomic. A crashed task leaves a partial file; retry with `overwrite = TRUE` replaces it cleanly.
- **Phase 5 — Narrow glob:** `R/queries.R` uses `_markers\\.parquet$` (not `\\.parquet$`) to exclude genotype files from the `markers` DuckDB view. The genotype files in `marker_sets/` have a different schema and must not be union-read with marker files.
- **Phase 5 — `.merge()` ordering:** `DB_MIGRATION_WRITE_TRAIT_DATA` uses `.merge()` to align `GCTA_SIMULATE_PHENOTYPES` output channels. `.merge()` aligns by emission order, not key — do not insert any reordering operator (`.filter()`, `.map()`, `.branch()`) between `GCTA_SIMULATE_PHENOTYPES.out.*` and the `.merge()` chain.
- **DB write barrier pattern:** `DB_MIGRATION_WRITE_MARKER_SET` emits `val true` on `done`; downstream `WRITE_GWA_TO_DB` gates on `.collect()` of all marker set completions combined with the GWA channel — ensures marker schemas exist before mappings are written.
- **Phase 5 — Pre-upscaled phenotype:** The phenotype stored in `phenotypes/` is the GCTA `.phen` output *before* `check_vp.py`. It may differ from the GWAS input by 1000× if upscaling was applied. This is intentional — the ANOVA SS ratio for offline variance-explained is scale-invariant.
- **GCTA thread pinning:** All GCTA steps except `--mlma-loco` use `--thread-num 1` to ensure deterministic floating-point results across runs. Multi-threaded BLAS introduces non-deterministic reduction order. `--mlma-loco` remains configurable via `GWA_THREADS=${task.cpus}` in `perform_gwa/main.nf`. The `gcta_make_grm` Rockfish label uses `cpus = 1`; `gcta_perform_gwa` uses `cpus = 4` to support loco parallelism.
- cSpell in the IDE flags many domain-specific terms as unknown (e.g. `nextflow`, `bcftools`, `gcta`, `Rscript`, `testthat`). These are false positives, not real errors.