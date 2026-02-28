# test-assessment_cross_validation.R - Cross-validation: DB-path vs existing-path assessment
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
#   TEST_EXISTING_ASSESSMENT - Path to simulation_assessment_results.tsv (existing path)
#   TEST_DB_ASSESSMENT       - Path to db_simulation_assessment_results.tsv (DB path)
#
# Known differences (by design):
#   - var.exp and Simulated.QTL.VarExp are NA in DB-path (no genotype matrix)
#   - algorithm_id format differs (both encode the same info, normalized for joining)
#
# Usage:
#   TEST_EXISTING_ASSESSMENT=/path/to/simulation_assessment_results.tsv \
#   TEST_DB_ASSESSMENT=/path/to/db_simulation_assessment_results.tsv \
#   Rscript tests/run_tests.R

existing_assessment_path <- Sys.getenv("TEST_EXISTING_ASSESSMENT", unset = "")
db_assessment_path <- Sys.getenv("TEST_DB_ASSESSMENT", unset = "")

skip_if_no_assessment_cross_validation <- function() {
  if (existing_assessment_path == "" || !file.exists(existing_assessment_path)) {
    skip("TEST_EXISTING_ASSESSMENT not set or file does not exist")
  }
  if (db_assessment_path == "" || !file.exists(db_assessment_path)) {
    skip("TEST_DB_ASSESSMENT not set or file does not exist")
  }
}

# Helper: parse the existing assessment TSV (no header, column structure from Assess_Sims.R)
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

# Helper: create a mapping key for joining rows across the two assessments
make_join_key <- function(df) {
  df %>%
    dplyr::mutate(
      norm_alg = normalize_algorithm_id(algorithm_id),
      join_key = paste(nQTL, simREP, h2, maf, effect_distribution,
        strain_set_id, norm_alg, QTL,
        sep = "|"
      )
    )
}

# Mapping group columns (replicate identity + method/threshold via norm_alg)
mapping_group <- c("nQTL", "simREP", "h2", "maf", "strain_set_id", "norm_alg")

# Helper: Jaccard overlap of two intervals [start_a, end_a] and [start_b, end_b].
# Returns 1 when start and end match perfectly, 0 when no overlap.
interval_jaccard <- function(start_a, end_a, start_b, end_b) {
  intersection <- pmax(0, pmin(end_a, end_b) - pmax(start_a, start_b))
  union <- pmax(end_a, end_b) - pmin(start_a, start_b)
  ifelse(union == 0, 1, intersection / union)
}

# Helper: join QTL rows from both paths on the composite key
build_joined_assessment <- function() {
  existing <- read_legacy_assessment(existing_assessment_path) %>%
    designate_qtl() %>%
    make_join_key()
  db_assess <- read_db_assessment(db_assessment_path) %>%
    designate_qtl() %>%
    make_join_key()

  dplyr::inner_join(
    existing %>%
      dplyr::select(join_key, norm_alg,
        nQTL, simREP, h2, maf, strain_set_id,
        designation_legacy = designation,
        Detected_legacy = Detected, Simulated_legacy = Simulated,
        startPOS_legacy = startPOS, peakPOS_legacy = peakPOS,
        endPOS_legacy = endPOS
      ),
    db_assess %>%
      dplyr::select(join_key,
        designation_db = designation,
        Detected_db = Detected, Simulated_db = Simulated,
        startPOS_db = startPOS, peakPOS_db = peakPOS,
        endPOS_db = endPOS
      ),
    by = "join_key"
  )
}


test_that("generate_mapping_id() is deterministic across invocations", {
  skip_if_no_assessment_cross_validation()

  # Read known params from DB assessment to get a real mapping's parameters
  db_assess <- read_db_assessment(db_assessment_path)
  params_row <- db_assess[1, ]

  algorithm <- params_row$mode   # "inbred" or "loco" per canonical form
  pca       <- params_row$type == "pca"

  ms1    <- generate_marker_set_id(params_row$strain_set_id, as.numeric(params_row$maf))
  trait1 <- generate_trait_id(ms1$hash, as.integer(params_row$nQTL),
                              params_row$effect_distribution,
                              as.integer(params_row$simREP),
                              as.numeric(params_row$h2))
  map1   <- generate_mapping_id(trait1$hash, algorithm, pca)

  ms2    <- generate_marker_set_id(params_row$strain_set_id, as.numeric(params_row$maf))
  trait2 <- generate_trait_id(ms2$hash, as.integer(params_row$nQTL),
                              params_row$effect_distribution,
                              as.integer(params_row$simREP),
                              as.numeric(params_row$h2))
  map2   <- generate_mapping_id(trait2$hash, algorithm, pca)

  expect_equal(map1$hash, map2$hash,
               label = "mapping_id is deterministic: same params → same hash")
  expect_true(grepl("^[0-9a-f]{20}$", map1$hash),
              label = "mapping_id is 20-char lowercase hex")
})


# ── Output File Existence ────────────────────────────────────────────────────

test_that("both assessment output files exist and are non-empty", {
  skip_if_no_assessment_cross_validation()

  expect_true(file.exists(existing_assessment_path),
    label = "existing assessment file exists"
  )
  expect_true(file.exists(db_assessment_path),
    label = "DB assessment file exists"
  )

  existing_size <- file.info(existing_assessment_path)$size
  db_size <- file.info(db_assessment_path)$size

  expect_gt(existing_size, 0, label = "existing assessment is non-empty")
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
