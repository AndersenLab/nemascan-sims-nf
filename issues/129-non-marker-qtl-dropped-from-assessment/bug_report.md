# Bug Report: `score_causal_markers()` — non-marker causal variants silently dropped from DB assessment output

**Date:** 2026-03-22
**Nextflow version:** 24.10.4
**Executor:** local / SLURM (Rockfish)

---

## Description

When `cv_maf < ms_maf`, causal variants can be sampled from positions that are not in
the GWA marker set (`TO_SIMS.bim`). These non-marker causal variants have no entry in
the GWA output and therefore no `log10p` value. `score_causal_markers()` in
`R/assessment.R` ends with `dplyr::filter(!is.na(log10p))`, which silently drops all
such variants from `effects_scores`.

Because they never enter `effects_scores`, they are absent from the union built by
`build_assessment_union()` and do not appear in the assessment TSV at all — not even as
false-negative (`Simulated=TRUE, Detected=FALSE`) rows. As a result, the HPC
verification assertion that checks for non-NA `Simulated.QTL.VarExp` on non-marker QTL
positions finds zero matching rows and fails:

```
Non-marker QTL with non-NA var.exp: 0
Error: nrow(non_marker_with_varexp) > 0 is not TRUE
Execution halted
```

The local Docker run (Assertion 3 in the testing plan) passed because it only checks
`sum(!is.na(Simulated.QTL.VarExp)) > 0` — satisfied by the single marker-position
causal variant (`3:6026848`) that survives the filter. The HPC verification is stricter:
it specifically checks that non-marker QTL positions carry a non-NA var.exp entry.

## Command Used

```bash
# Pipeline run (HPC)
nextflow run main.nf \
  -profile test_cv_pool,rockfish \
  -work-dir /scratch4/eande106/Ryan/nf-work-cv-pool \
  --output_dir results_hpc_cv_pool

# Local verification (after rsync — run from repo root)
Rscript - << 'EOF'
library(arrow)
library(dplyr)

output_dir <- "results_hpc_cv_pool"

# Load causal variant positions for all traits
cv_files <- list.files(file.path(output_dir, "db/traits/causal_variants"),
                        pattern = "_causal\\.parquet$", full.names = TRUE)
causal <- bind_rows(lapply(cv_files, function(f) as.data.frame(read_parquet(f))))

# Load marker positions to identify non-marker causal variants
marker_files <- list.files(file.path(output_dir, "db/markers/marker_sets"),
                             pattern = "_markers\\.parquet$", full.names = TRUE)
markers <- bind_rows(lapply(marker_files, function(f) as.data.frame(read_parquet(f)))) %>%
  distinct(CHROM, POS)

non_marker_qtl <- causal %>%
  mutate(CHROM = as.character(CHROM), POS = as.integer(POS)) %>%
  anti_join(markers, by = c("CHROM", "POS")) %>%
  pull(QTL) %>% unique()

cat("Non-marker causal QTL positions:", length(non_marker_qtl), "\n")
stopifnot(length(non_marker_qtl) > 0)

# Check at least one has non-NA var.exp in the assessment output
source("R/setup.R")
source("R/assessment.R")
stub <- data.frame(
  nQTL=1L, simREP=1L, h2=0.5, maf=0.05,
  effect_distribution="gamma", strain_set_id="g",
  mode="inbred", type="pca", threshold="BF",
  algorithm_id="inbred_pca", alpha=0.05,
  ci_size=150L, snp_grouping=1000L,
  QTL="1:100", CHROM="1", POS=100L,
  Simulated="TRUE", Detected="TRUE",
  startPOS=1L, peakPOS=100L, endPOS=200L,
  interval_size=200L, min_log10p=5.0, peak_log10p=6.0,
  Simulated.QTL.VarExp=NA_real_, interval.var.exp=NA_real_,
  designation="TP", nDetected_QTL=1L,
  trait_id="abc", marker_set_id="def", mapping_id="ghi",
  stringsAsFactors=FALSE
)
col_order <- colnames(format_assessment_tsv(stub))

res <- read.table(
  file.path(output_dir, "db_simulation_assessment_results.tsv"),
  sep = "\t", header = FALSE,
  col.names = col_order,
  stringsAsFactors = FALSE
)

non_marker_with_varexp <- res %>%
  filter(QTL %in% non_marker_qtl, !is.na(Simulated.QTL.VarExp))
cat("Non-marker QTL with non-NA var.exp:", nrow(non_marker_with_varexp), "\n")
stopifnot(nrow(non_marker_with_varexp) > 0)
cat("PASS\n")
EOF
```

## Error Output

```
Non-marker QTL with non-NA var.exp: 0
Error: nrow(non_marker_with_varexp) > 0 is not TRUE
Execution halted
```

## Config / Profile

```groovy
// test_cv_pool profile — cv_maf=0.01 is lower than ms_maf=0.05,
// so the CV pool is broader than the marker set and non-marker
// causal variant positions are possible.
profiles {
    test_cv_pool {
        params.cv_maf = 0.01
        params.cv_ld  = 0.99
        // inherits all other params from the test profile
    }
}
```

## Additional Context

- Root function: `score_causal_markers()` in `R/assessment.R:153–173`
- Offending line: `dplyr::filter(!is.na(log10p))` at `R/assessment.R:170`
- Non-marker causal variants confirmed present in the run:
  `1:2217428`, `2:9853034`, `2:11624318`, `3:16593852` (4 of 5 simulated QTLs)
- Per-trait causal genotypes **are** written correctly to
  `db/traits/causal_genotypes/147b00054e04d5e7402e_causal_geno.parquet` (1000 rows).
  `assess_sims.R` reads them and passes them through `compute_var_exp_anova()` — but
  the computed `Simulated.QTL.VarExp` values are lost when the rows are filtered out.
- The same `!is.na(log10p)` filter exists in the legacy `bin/Assess_Sims.R:121`. The
  legacy path was not designed for CV pool runs, so this is a design difference, not a
  bug in the legacy path. See the RCA comment for concordance impact analysis.

---

*Generated with the Quarto GitHub Issue Framework. Submit via `gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
