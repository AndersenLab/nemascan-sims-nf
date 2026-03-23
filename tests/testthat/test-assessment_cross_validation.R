# test-assessment_cross_validation.R - Cross-validation: DB-path vs legacy-path assessment
#
# Compares db_simulation_assessment_results.tsv (from DB_MIGRATION_ANALYZE_QTL + DB_MIGRATION_ASSESS_SIMS)
# against simulation_assessment_results.tsv (from R_GET_GCTA_INTERVALS + R_ASSESS_SIMS) to verify
# that the DB-path assessment reproduces the legacy pipeline's analytical results.
#
# Computes three per-replicate concordance scores (matching the analysis in
# tests/notebooks/assessment_cross_validation.qmd):
#   1. Detection Concordance  - Detected flag agreement per QTL row
#   2. Designation Concordance - designate_qtl() category agreement per QTL row
#   3. Interval Overlap Concordance - Jaccard overlap of [startPOS, endPOS] intervals
#
# A mapping entry is identified by (nQTL, simREP, h2, maf, strain_set_id, norm_alg).
# Each score is 1.0 for perfect concordance, 0.0 for no concordance.
#
# Requires pipeline output from BOTH paths:
#   TEST_LEGACY_ASSESSMENT - Path to simulation_assessment_results.tsv (legacy path)
#   TEST_DB_ASSESSMENT     - Path to db_simulation_assessment_results.tsv (DB path)
#
# Known differences (by design):
#   - interval.var.exp is NA in DB-path (raw GWA r² not stored in Phase 5 DB; permanent design difference)
#   - Simulated.QTL.VarExp is populated in DB-path via ANOVA (see test-assessment_var_exp.R)
#   - algorithm_id format differs (both encode the same info, normalized for joining)
#
# Usage:
#   TEST_LEGACY_ASSESSMENT=/path/to/simulation_assessment_results.tsv \
#   TEST_DB_ASSESSMENT=/path/to/db_simulation_assessment_results.tsv \
#   Rscript tests/run_tests.R

legacy_assessment_path <- Sys.getenv("TEST_LEGACY_ASSESSMENT", unset = "")
db_assessment_path <- Sys.getenv("TEST_DB_ASSESSMENT", unset = "")

skip_if_no_assessment_cross_validation <- function() {
  if (legacy_assessment_path == "" || !file.exists(legacy_assessment_path)) {
    skip("TEST_LEGACY_ASSESSMENT not set or file does not exist")
  }
  if (db_assessment_path == "" || !file.exists(db_assessment_path)) {
    skip("TEST_DB_ASSESSMENT not set or file does not exist")
  }
}

# Helper: parse the legacy assessment TSV (no header, column structure from Assess_Sims.R)
read_legacy_assessment <- function(path) {
  col_names <- c(
    "QTL", "Simulated", "Detected", "CHROM", "POS", "RefAllele",
    "Frequency", "Effect", "Simulated.QTL.VarExp", "log10p", "significant",
    "startPOS", "peakPOS", "endPOS",
    "detected.peak", "BETA", "interval.log10p",
    "peak_id", "interval_size",
    "interval.var.exp", "interval.Frequency",
    "top.hit", "nQTL", "simREP", "h2", "maf",
    "effect_distribution", "strain_set_id",
    "mode", "type", "threshold", "algorithm_id",
    "alpha", "ci_size", "snp_grouping"
  )
  # Deal with presence/absence of header in output files
  # needed to accommodate different versions of legacy assessment output (with or without header) and avoid test failures due to that
  df <- tryCatch(
    data.table::fread(path,
      header = FALSE, col.names = col_names,
      sep = "\t", na.strings = c("NA", "")
    ) %>%
      as.data.frame(),
    error = function(e) {
      data.table::fread(path,
        header = TRUE, sep = "\t",
        na.strings = c("NA", "")
      ) %>%
        as.data.frame()
    }
  )

  df$QTL <- as.character(df$QTL)
  df$Simulated <- as.character(df$Simulated)
  df$Detected <- as.character(df$Detected)
  df$startPOS <- as.integer(df$startPOS)
  df$peakPOS <- as.integer(df$peakPOS)
  df$endPOS <- as.integer(df$endPOS)
  df$nQTL <- as.character(df$nQTL)
  df$simREP <- as.character(df$simREP)
  df$h2 <- as.character(df$h2)
  df$maf <- as.character(df$maf)
  df$strain_set_id <- as.character(df$strain_set_id)

  df
}

# Helper: parse the DB assessment TSV (no header, from assess_sims.R via format_assessment_tsv())
read_db_assessment <- function(path) {
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
    ) %>%
      as.data.frame(),
    error = function(e) {
      data.table::fread(path,
        header = TRUE, sep = "\t",
        na.strings = c("NA", "")
      ) %>%
        as.data.frame()
    }
  )

  df$QTL <- as.character(df$QTL)
  df$Simulated <- as.character(df$Simulated)
  df$Detected <- as.character(df$Detected)
  df$startPOS <- as.integer(df$startPOS)
  df$peakPOS <- as.integer(df$peakPOS)
  df$endPOS <- as.integer(df$endPOS)
  df$nQTL <- as.character(df$nQTL)
  df$simREP <- as.character(df$simREP)
  df$h2 <- as.character(df$h2)
  df$maf <- as.character(df$maf)
  df$strain_set_id <- as.character(df$strain_set_id)

  df
}

# Helper: normalize algorithm_id across both paths for joining
# Both paths now use canonical algorithm values ("inbred"/"loco").
# Normalization strips no prefix but uppercases for comparison:
#   "inbred_nopca_EIGEN" -> "INBRED_NOPCA_EIGEN"
normalize_algorithm_id <- function(alg_id) {
  toupper(sub("^LMM-EXACT-", "", toupper(alg_id)))
}

# Helper: normalize and add join columns to a data frame for cross-path joining.
# All columns that discriminate rows (effect_distribution, norm_alg) must be present
# so that the join is 1-to-1 per (mapping, QTL). Dropping any discriminating column
# creates a many-to-many relationship because the same QTL position appears once per
# (effect_distribution × algorithm) combination.
add_join_cols <- function(df) {
  df %>%
    dplyr::mutate(norm_alg = normalize_algorithm_id(algorithm_id))
}

# Columns that uniquely identify one (mapping, QTL) row in both assessment outputs.
qtl_join_cols <- c("nQTL", "simREP", "h2", "maf", "effect_distribution",
                   "strain_set_id", "norm_alg", "QTL")

# Mapping group columns (replicate identity + method/threshold via norm_alg)
mapping_group <- c("nQTL", "simREP", "h2", "maf", "strain_set_id", "norm_alg")

# Helper: Jaccard overlap of two intervals [start_a, end_a] and [start_b, end_b].
# Returns 1 when start and end match perfectly, 0 when no overlap.
interval_jaccard <- function(start_a, end_a, start_b, end_b) {
  intersection <- pmax(0, pmin(end_a, end_b) - pmax(start_a, start_b))
  union <- pmax(end_a, end_b) - pmin(start_a, start_b)
  ifelse(union == 0, 1, intersection / union)
}

# Helper: join QTL rows from both paths on the typed multi-column key.
# Uses explicit column-by-column join (not a concatenated string) so that
# mismatches on any single dimension are surfaced clearly in test failures.
# relationship = "one-to-one" turns silent many-to-many into an error.
build_joined_assessment <- function() {
  legacy <- read_legacy_assessment(legacy_assessment_path) %>%
    designate_qtl() %>%
    add_join_cols()
  db_assess <- read_db_assessment(db_assessment_path) %>%
    designate_qtl() %>%
    add_join_cols()

  dplyr::inner_join(
    legacy %>%
      dplyr::select(dplyr::all_of(qtl_join_cols), norm_alg,
        nQTL, simREP, h2, maf, strain_set_id,
        designation_legacy = designation,
        Detected_legacy = Detected, Simulated_legacy = Simulated,
        startPOS_legacy = startPOS, peakPOS_legacy = peakPOS,
        endPOS_legacy = endPOS
      ),
    db_assess %>%
      dplyr::select(dplyr::all_of(qtl_join_cols),
        designation_db = designation,
        Detected_db = Detected, Simulated_db = Simulated,
        startPOS_db = startPOS, peakPOS_db = peakPOS,
        endPOS_db = endPOS
      ),
    by = qtl_join_cols,
    relationship = "one-to-one"
  )
}


test_that("generate_mapping_id() is deterministic across invocations", {
  skip_if_no_assessment_cross_validation()

  # Read known params from DB assessment to get a real mapping's parameters
  db_assess <- read_db_assessment(db_assessment_path)
  params_row <- db_assess[1, ]

  algorithm <- params_row$mode   # "inbred" or "loco" per canonical form
  pca       <- params_row$type == "pca"

  db_dir <- Sys.getenv("TEST_DB_DIR", unset = NA_character_)
  if (is.na(db_dir)) skip("TEST_DB_DIR not set")

  # First independent computation
  ms_meta1 <- read_marker_set_metadata(
    params_row$strain_set_id, as.numeric(params_row$maf), db_dir
  )
  if (is.null(ms_meta1)) skip("marker set metadata not found in TEST_DB_DIR")

  ms1 <- generate_marker_set_id(
    params_row$strain_set_id, as.numeric(params_row$maf),
    ms_meta1$species, ms_meta1$vcf_release_id, as.numeric(ms_meta1$ms_ld)
  )
  trait1 <- generate_trait_id(ms1$hash, as.integer(params_row$nQTL),
                              params_row$effect_distribution,
                              as.integer(params_row$simREP),
                              as.numeric(params_row$h2),
                              cv_maf_effective = as.numeric(ms_meta1$ms_maf),
                              cv_ld            = as.numeric(ms_meta1$ms_ld))
  map1   <- generate_mapping_id(trait1$hash, algorithm, pca)

  # Second independent computation with identical params — tests determinism
  ms_meta2 <- read_marker_set_metadata(
    params_row$strain_set_id, as.numeric(params_row$maf), db_dir
  )
  if (is.null(ms_meta2)) skip("marker set metadata not found in TEST_DB_DIR (second read)")

  ms2 <- generate_marker_set_id(
    params_row$strain_set_id, as.numeric(params_row$maf),
    ms_meta2$species, ms_meta2$vcf_release_id, as.numeric(ms_meta2$ms_ld)
  )
  trait2 <- generate_trait_id(ms2$hash, as.integer(params_row$nQTL),
                              params_row$effect_distribution,
                              as.integer(params_row$simREP),
                              as.numeric(params_row$h2),
                              cv_maf_effective = as.numeric(ms_meta2$ms_maf),
                              cv_ld            = as.numeric(ms_meta2$ms_ld))
  map2   <- generate_mapping_id(trait2$hash, algorithm, pca)

  expect_equal(map1$hash, map2$hash,
               label = "mapping_id is deterministic: same params → same hash (two independent calls)")
  expect_true(grepl("^[0-9a-f]{20}$", map1$hash),
              label = "mapping_id is 20-char lowercase hex")
})


# ── Output File Existence ────────────────────────────────────────────────────

test_that("both assessment output files exist and are non-empty", {
  skip_if_no_assessment_cross_validation()

  expect_true(file.exists(legacy_assessment_path),
    label = "legacy assessment file exists"
  )
  expect_true(file.exists(db_assessment_path),
    label = "DB assessment file exists"
  )

  legacy_size <- file.info(legacy_assessment_path)$size
  db_size <- file.info(db_assessment_path)$size

  expect_gt(legacy_size, 0, label = "legacy assessment is non-empty")
  expect_gt(db_size, 0, label = "DB assessment is non-empty")
})


# ── Detection Concordance ────────────────────────────────────────────────────

test_that("per-mapping detection concordance is 1.0 for all mappings", {
  skip_if_no_assessment_cross_validation()

  joined <- build_joined_assessment()
  expect_gt(nrow(joined), 0, label = "joined assessment has rows")

  detection_conc <- joined %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(mapping_group))) %>%
    dplyr::summarise(
      n_qtl = dplyr::n(),
      n_detect_agree = sum(Detected_legacy == Detected_db),
      detection_concordance = n_detect_agree / n_qtl,
      .groups = "drop"
    )

  imperfect <- detection_conc %>%
    dplyr::filter(detection_concordance < 1)

  expect_equal(
    nrow(imperfect), 0,
    label = sprintf(
      "all %d mappings have perfect detection concordance (found %d < 1.0)",
      nrow(detection_conc), nrow(imperfect)
    )
  )
})


# ── Designation Concordance ──────────────────────────────────────────────────

test_that("per-mapping designation concordance is 1.0 for all mappings", {
  skip_if_no_assessment_cross_validation()

  joined <- build_joined_assessment()
  expect_gt(nrow(joined), 0, label = "joined assessment has rows")

  designation_conc <- joined %>%
    dplyr::filter(!is.na(designation_legacy), !is.na(designation_db)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(mapping_group))) %>%
    dplyr::summarise(
      n_designated = dplyr::n(),
      n_desig_agree = sum(designation_legacy == designation_db),
      designation_concordance = n_desig_agree / n_designated,
      .groups = "drop"
    )

  expect_gt(nrow(designation_conc), 0,
    label = "at least one mapping has designated QTLs"
  )

  imperfect <- designation_conc %>%
    dplyr::filter(designation_concordance < 1)

  expect_equal(
    nrow(imperfect), 0,
    label = sprintf(
      "all %d mappings have perfect designation concordance (found %d < 1.0)",
      nrow(designation_conc), nrow(imperfect)
    )
  )
})


# ── Interval Overlap Concordance ─────────────────────────────────────────────

test_that("per-mapping interval Jaccard concordance is 1.0 for all mappings", {
  skip_if_no_assessment_cross_validation()

  joined <- build_joined_assessment()
  expect_gt(nrow(joined), 0, label = "joined assessment has rows")

  interval_conc_per_qtl <- joined %>%
    dplyr::filter(
      Detected_legacy == "TRUE", Detected_db == "TRUE",
      !is.na(startPOS_legacy), !is.na(startPOS_db),
      !is.na(endPOS_legacy), !is.na(endPOS_db)
    ) %>%
    dplyr::mutate(
      qtl_interval_jaccard = interval_jaccard(
        startPOS_legacy, endPOS_legacy, startPOS_db, endPOS_db
      )
    )

  if (nrow(interval_conc_per_qtl) == 0) {
    skip("no QTL intervals detected by both paths to compare")
  }

  interval_conc <- interval_conc_per_qtl %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(mapping_group))) %>%
    dplyr::summarise(
      n_intervals = dplyr::n(),
      n_perfect = sum(qtl_interval_jaccard == 1),
      interval_concordance = mean(qtl_interval_jaccard),
      .groups = "drop"
    )

  imperfect <- interval_conc %>%
    dplyr::filter(interval_concordance < 1)

  expect_equal(
    nrow(imperfect), 0,
    label = sprintf(
      "all %d mappings have perfect interval concordance (found %d < 1.0)",
      nrow(interval_conc), nrow(imperfect)
    )
  )

  # Also verify every individual QTL interval is a perfect match
  n_imperfect_qtls <- sum(interval_conc_per_qtl$qtl_interval_jaccard < 1)
  expect_equal(
    n_imperfect_qtls, 0,
    label = sprintf(
      "all %d QTL intervals have Jaccard = 1.0 (found %d < 1.0)",
      nrow(interval_conc_per_qtl), n_imperfect_qtls
    )
  )
})


# ── Variance-Explained Concordance ────────────────────────────────────────────
#
# interval.var.exp is NOT compared between paths:
#   Legacy: GCTA's LMM-adjusted r² at the peak SNP (accounts for polygenic background via GRM).
#   DB path: NA_real_ by design — raw per-marker GWA r² was not stored in Phase 5. Even if
#   attempted via ANOVA, peak GWAS markers are rarely the causal variant, so the join would
#   produce NA anyway. This is a permanent design difference, not a missing feature.
#
# Non-marker causal variants (DB path only, when cv_maf < ms_maf):
#   DB path: appear as Simulated=TRUE, Detected=FALSE FN rows with non-NA Simulated.QTL.VarExp.
#   Legacy path: entirely absent (bin/Assess_Sims.R applies the same !is.na(log10p) filter
#   and has no cv_maf < ms_maf support). The inner_join on QTL used for concordance comparison
#   naturally excludes these DB-only rows — this is intentional, not an oversight.

test_that("Simulated.QTL.VarExp is populated in DB assessment output", {
  skip_if_no_assessment_cross_validation()

  db_assess <- read_db_assessment(db_assessment_path)
  expect_gt(
    sum(!is.na(db_assess$Simulated.QTL.VarExp)), 0,
    label = "at least some rows have non-NA Simulated.QTL.VarExp in DB output"
  )
})


test_that("Simulated.QTL.VarExp is concordant between DB and legacy paths", {
  skip_if_no_assessment_cross_validation()

  legacy    <- read_legacy_assessment(legacy_assessment_path) %>%
    dplyr::mutate(norm_alg = normalize_algorithm_id(algorithm_id))
  db_assess <- read_db_assessment(db_assessment_path) %>%
    dplyr::mutate(norm_alg = normalize_algorithm_id(algorithm_id))

  # Join key must include effect_distribution and norm_alg to avoid many-to-many:
  # the same QTL position appears once per (effect_distribution × algorithm) combination,
  # so omitting either column causes one legacy row to match multiple DB rows.
  join_cols <- c("QTL", "nQTL", "simREP", "h2", "maf", "strain_set_id",
                 "effect_distribution", "norm_alg")

  joined <- dplyr::inner_join(
    legacy %>%
      dplyr::select(dplyr::all_of(join_cols),
                    varexp_legacy = Simulated.QTL.VarExp),
    db_assess %>%
      dplyr::select(dplyr::all_of(join_cols),
                    varexp_db     = Simulated.QTL.VarExp),
    by = join_cols,
    relationship = "one-to-one"
  ) %>%
    dplyr::filter(!is.na(varexp_legacy), !is.na(varexp_db))

  if (nrow(joined) == 0) skip("No non-NA Simulated.QTL.VarExp rows to compare")

  expect_equal(joined$varexp_legacy, joined$varexp_db, tolerance = 1e-6)
})
