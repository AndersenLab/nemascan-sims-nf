# test-cv_pool.R - Integration tests for the CV pool broader than marker set scenario
#
# Tests the case where cv_maf < ms_maf, producing causal variants that are absent from
# GWA output (non-marker causal variants). These appear as FN rows with NA log10p in
# the assessment output and require per-trait causal genotype storage for var.exp.
#
# Requires:
#   TEST_DB_DIR        - Path to pipeline database output (test_cv_pool profile)
#   TEST_DB_ASSESSMENT - Path to db_simulation_assessment_results.tsv
#   TEST_CV_POOL=true  - Enable cv_pool-specific assertions
#
# Do NOT set TEST_LEGACY_ASSESSMENT for this scenario. The legacy path is incompatible:
# non-marker causal variants are absent from GWA files and invisible to Assess_Sims.R,
# so a legacy TSV would silently misrepresent False Negative rates.
#
# Usage:
#   TEST_DB_DIR=tests/integration_data/test_cv_pool/db \
#   TEST_DB_ASSESSMENT=tests/integration_data/test_cv_pool/db_simulation_assessment_results.tsv \
#   TEST_CV_POOL=true \
#   Rscript tests/run_tests.R

db_dir             <- Sys.getenv("TEST_DB_DIR", unset = "")
db_assessment_path <- Sys.getenv("TEST_DB_ASSESSMENT", unset = "")
cv_pool_active     <- identical(Sys.getenv("TEST_CV_POOL", unset = "false"), "true")

skip_if_no_cv_pool <- function() {
  if (!cv_pool_active) {
    skip("TEST_CV_POOL=true not set — skipping cv_pool integration tests")
  }
  if (db_dir == "" || !dir.exists(db_dir)) {
    skip("TEST_DB_DIR not set or directory does not exist")
  }
}

skip_if_no_cv_pool_assessment <- function() {
  skip_if_no_cv_pool()
  if (db_assessment_path == "" || !file.exists(db_assessment_path)) {
    skip("TEST_DB_ASSESSMENT not set or file does not exist")
  }
}

# Parse the DB assessment TSV (no header, 35-column format from format_assessment_tsv())
read_cv_pool_assessment <- function(path) {
  col_names <- c(
    "QTL", "Simulated", "Detected", "CHROM", "POS", "RefAllele",
    "Frequency", "Effect", "Simulated.QTL.VarExp", "log10p", "significant",
    "BETA", "startPOS", "peakPOS", "endPOS", "detected.peak",
    "interval.log10p", "interval.var.exp", "interval.Frequency",
    "peak_id", "interval_size", "top.hit",
    "nQTL", "simREP", "h2", "maf", "effect_distribution",
    "strain_set_id", "mode", "type", "threshold", "algorithm_id",
    "alpha", "ci_size", "snp_grouping"
  )
  df <- tryCatch(
    data.table::fread(path,
      header = FALSE, col.names = col_names,
      sep = "\t", na.strings = c("NA", "")
    ) %>% as.data.frame(),
    error = function(e) {
      data.table::fread(path, header = TRUE, sep = "\t", na.strings = c("NA", "")) %>%
        as.data.frame()
    }
  )
  df$QTL       <- as.character(df$QTL)
  df$Simulated <- as.character(df$Simulated)
  df$Detected  <- as.character(df$Detected)
  df
}

# ── Causal Genotypes Directory ───────────────────────────────────────────────

test_that("traits/causal_genotypes/ directory exists and contains parquet files", {
  skip_if_no_cv_pool()

  cg_dir <- file.path(db_dir, "traits", "causal_genotypes")
  expect_true(dir.exists(cg_dir),
    label = "traits/causal_genotypes/ subdirectory exists for cv_pool scenario"
  )

  cg_files <- list.files(cg_dir, pattern = "\\.parquet$", full.names = TRUE)
  expect_gt(length(cg_files), 0,
    label = "at least one causal genotype parquet file exists"
  )
})

test_that("causal genotype parquets have correct schema", {
  skip_if_no_cv_pool()

  cg_dir   <- file.path(db_dir, "traits", "causal_genotypes")
  cg_files <- list.files(cg_dir, pattern = "\\.parquet$", full.names = TRUE)
  if (length(cg_files) == 0) skip("no causal genotype parquet files found")

  cg     <- arrow::read_parquet(cg_files[1], as_data_frame = FALSE)
  schema <- cg$schema

  expected_cols <- c("trait_id", "QTL", "CHROM", "POS", "strain", "allele")
  expect_true(all(expected_cols %in% names(schema)),
    label = paste(
      "causal genotype parquet missing columns:",
      paste(setdiff(expected_cols, names(schema)), collapse = ", ")
    )
  )
  expect_gt(nrow(cg), 0, label = "causal genotype parquet is non-empty")
})

# ── Non-Marker Causal Variant Representation ─────────────────────────────────

test_that("non-marker causal variants appear as FN rows with NA log10p in assessment", {
  skip_if_no_cv_pool_assessment()

  db_assess <- read_cv_pool_assessment(db_assessment_path)

  # Non-marker causal variants: Simulated=TRUE but absent from GWA output → log10p=NA.
  # Present when cv_maf < ms_maf so some causal variants have MAF below the GWA threshold.
  fn_non_marker <- db_assess[db_assess$Simulated == "TRUE" & is.na(db_assess$log10p), ]

  expect_gt(nrow(fn_non_marker), 0,
    label = paste(
      "non-marker causal variants (Simulated=TRUE, log10p=NA) appear in assessment;",
      "confirm pipeline was run with test_cv_pool profile (cv_maf=0.01, ms_maf=0.05)"
    )
  )
})

test_that("non-marker FN rows have significant=NA (not FALSE) in assessment", {
  skip_if_no_cv_pool_assessment()

  db_assess     <- read_cv_pool_assessment(db_assessment_path)
  fn_non_marker <- db_assess[db_assess$Simulated == "TRUE" & is.na(db_assess$log10p), ]
  if (nrow(fn_non_marker) == 0) skip("no non-marker FN rows found")

  # significant=NA means "not in GWA output" — distinct from significant=FALSE which
  # means "tested by GWA but below threshold". This distinction is important for
  # designate_qtl() to correctly classify these rows as Missed.CV.
  expect_true(all(is.na(fn_non_marker$significant)),
    label = "non-marker FN rows have significant=NA (not FALSE — they were never tested by GWA)"
  )
})

test_that("Simulated.QTL.VarExp is populated for non-marker FN rows", {
  skip_if_no_cv_pool_assessment()

  db_assess     <- read_cv_pool_assessment(db_assessment_path)
  fn_non_marker <- db_assess[db_assess$Simulated == "TRUE" & is.na(db_assess$log10p), ]
  if (nrow(fn_non_marker) == 0) skip("no non-marker FN rows found")

  # var.exp is computed via ANOVA from per-trait causal genotypes (traits/causal_genotypes/).
  # Non-marker causal variants still have var.exp because genotype data is stored
  # per-trait independently of the GWA marker set.
  expect_gt(
    sum(!is.na(fn_non_marker$Simulated.QTL.VarExp)), 0,
    label = "at least some non-marker FN rows have non-NA Simulated.QTL.VarExp (from causal_genotypes)"
  )
})

test_that("non-marker causal variants appear in causal_variants/ but not in markers DB view", {
  skip_if_no_cv_pool()

  cv_files <- list.files(
    file.path(db_dir, "traits", "causal_variants"),
    pattern = "_causal\\.parquet$", full.names = TRUE
  )
  if (length(cv_files) == 0) skip("no causal variant parquet files found")

  # Sample a subset of causal variant parquets
  all_causal <- do.call(rbind, lapply(head(cv_files, 5), arrow::read_parquet))
  causal_qtls <- unique(all_causal$QTL)

  # Query the markers DB view to get all marker IDs
  con <- open_mapping_db(db_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  marker_ids <- DBI::dbGetQuery(con, "SELECT DISTINCT marker FROM markers")$marker

  # At least one causal QTL should be absent from the markers view
  non_marker_causal <- causal_qtls[!causal_qtls %in% marker_ids]
  expect_gt(
    length(non_marker_causal), 0,
    label = paste(
      "at least one causal QTL is absent from the GWA marker set;",
      "expected for test_cv_pool profile (cv_maf=0.01 < ms_maf=0.05)"
    )
  )
})

# ── designate_qtl() Correctness with significant=NA ──────────────────────────

test_that("designate_qtl() classifies non-marker FN rows as Missed.CV", {
  skip_if_no_cv_pool_assessment()

  db_assess     <- read_cv_pool_assessment(db_assessment_path)
  fn_non_marker <- db_assess[db_assess$Simulated == "TRUE" & is.na(db_assess$log10p), ]
  if (nrow(fn_non_marker) == 0) skip("no non-marker FN rows found")

  designated <- designate_qtl(fn_non_marker)
  # Non-marker FN rows: Simulated=TRUE, Detected=FALSE, significant=NA.
  # After the designate_qtl() fix, these must be classified as Missed.CV so they
  # contribute to the Power denominator in calculate_simrep_performance().
  expect_true(all(designated$designation == "Missed.CV"),
    label = "non-marker FN rows (significant=NA) are classified as Missed.CV"
  )
})

# ── Standalone Power/FDR Check (no legacy comparison needed) ─────────────────

test_that("calculate_simrep_performance() runs correctly on cv_pool assessment data", {
  skip_if_no_cv_pool_assessment()

  db_assess <- read_cv_pool_assessment(db_assessment_path)

  perf <- db_assess %>%
    designate_qtl() %>%
    dplyr::group_by(nQTL, simREP, h2, maf, effect_distribution, strain_set_id, mode, type) %>%
    count_outcomes() %>%
    calculate_simrep_performance()

  expect_gt(nrow(perf), 0, label = "performance dataframe has rows")
  expect_true(all(c("Power", "FDR") %in% names(perf)))
  expect_true(all(perf$Power >= 0 & perf$Power <= 1),
    label = "Power values in [0, 1]"
  )
  expect_true(all(perf$FDR >= 0 & perf$FDR <= 1),
    label = "FDR values in [0, 1]"
  )

  # With non-marker causal variants, some simulated QTLs are structurally undetectable.
  # Power must be < 1 for at least one replicate (validates the cv_pool scenario was exercised).
  expect_true(any(perf$Power < 1),
    label = "Power < 1 for at least one replicate — non-marker causal variants reduce detectable power"
  )
})
