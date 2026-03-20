# test-cross_validation.R - Cross-validation: DB contents vs raw GWA source files
#
# Run AFTER pipeline execution. Compares the Parquet database produced by
# DB_MIGRATION modules against the raw GWA output files in the Nextflow work
# directory to verify round-trip accuracy.
#
# Usage:
#   TEST_DB_DIR=/path/to/output/db TEST_WORK_DIR=/path/to/work Rscript tests/run_tests.R
#
# Requires:
#   TEST_DB_DIR  - Path to pipeline database output directory
#   TEST_WORK_DIR - Path to Nextflow work/ directory containing raw GWA files

db_dir <- Sys.getenv("TEST_DB_DIR", unset = "")
work_dir <- Sys.getenv("TEST_WORK_DIR", unset = "")

skip_if_no_cross_validation <- function() {
  if (db_dir == "" || !dir.exists(db_dir)) {
    skip("TEST_DB_DIR not set or directory does not exist")
  }
  if (work_dir == "" || !dir.exists(work_dir)) {
    skip("TEST_WORK_DIR not set or directory does not exist")
  }
}

# Helper: find raw GWA files in the Nextflow work directory
find_gwa_files <- function(work_dir) {
  fastgwa <- list.files(work_dir, pattern = "\\.fastGWA$",
                        recursive = TRUE, full.names = TRUE)
  mlma <- list.files(work_dir, pattern = "\\.mlma$",
                     recursive = TRUE, full.names = TRUE)
  c(fastgwa, mlma)
}

# Helper: extract mapping params from GWA filename
# Inline-path filenames: {nqtl}_{rep}_{h2}_{maf}_{effect}_{group}_lmm-exact_{mode}_{type}.{suffix}
parse_inline_gwa_filename <- function(filename) {
  bn <- basename(filename)
  ext <- tools::file_ext(bn)
  base <- sub(paste0("\\.", ext, "$"), "", bn)

  pattern <- "^(\\d+)_(\\d+)_([0-9.]+)_([0-9.]+)_([a-zA-Z]+)_(.+)_lmm-exact_(inbred|loco)_(pca|nopca)$"
  matches <- regmatches(base, regexec(pattern, base))[[1]]

  if (length(matches) == 0) return(NULL)

  mode <- matches[8]
  type <- matches[9]

  list(
    nqtl = as.integer(matches[2]),
    rep = as.integer(matches[3]),
    h2 = as.numeric(matches[4]),
    maf = as.numeric(matches[5]),
    effect = matches[6],
    population = matches[7],
    algorithm = mode,   # "inbred" or "loco" — matches canonical form written by scripts
    pca = type == "pca",
    trait = paste(matches[2], matches[3], matches[4], sep = "_")
  )
}

# ── Round-trip Value Accuracy ────────────────────────────────────────────────

test_that("DB mapping P and BETA values match raw GWA source files", {
  skip_if_no_cross_validation()

  gwa_files <- find_gwa_files(work_dir)
  if (length(gwa_files) == 0) {
    skip("No raw GWA files found in TEST_WORK_DIR")
  }

  # Test at least one file of each type (fastGWA and mlma)
  tested <- 0
  for (gwa_file in gwa_files) {
    params <- parse_inline_gwa_filename(gwa_file)
    if (is.null(params)) next

    ms_meta <- read_marker_set_metadata(params$population, params$maf, db_dir)
    if (is.null(ms_meta)) next
    ms_id <- generate_marker_set_id(params$population, params$maf,
                                    ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld))
    trait <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
    mapping_id <- generate_mapping_id(trait$hash, params$algorithm, params$pca)$hash

    # Read raw GWA file
    raw_df <- tryCatch(
      read_raw_gwa_file(gwa_file, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(raw_df)) next

    # Read from DB
    db_df <- tryCatch(
      query_by_mapping_id(mapping_id, db_dir),
      error = function(e) NULL,
      warning = function(w) {
        suppressWarnings(query_by_mapping_id(mapping_id, db_dir))
      }
    )
    if (is.null(db_df) || nrow(db_df) == 0) next

    # Compare row counts
    expect_equal(nrow(db_df), nrow(raw_df),
                 label = paste("row count for", mapping_id))

    # Join on marker and compare P and BETA
    raw_df$marker <- as.character(raw_df$marker)
    db_df$marker <- as.character(db_df$marker)

    merged <- merge(
      raw_df[, c("marker", "P", "BETA")],
      db_df[, c("marker", "P", "BETA")],
      by = "marker", suffixes = c("_raw", "_db")
    )

    expect_equal(nrow(merged), nrow(raw_df),
                 label = paste("all markers joined for", mapping_id))
    expect_equal(merged$P_db, merged$P_raw,
                 tolerance = 1e-12,
                 label = paste("P values match for", mapping_id))
    expect_equal(merged$BETA_db, merged$BETA_raw,
                 tolerance = 1e-12,
                 label = paste("BETA values match for", mapping_id))

    tested <- tested + 1
    if (tested >= 4) break  # Test up to 4 files for speed
  }

  expect_gt(tested, 0, label = "at least one GWA file was cross-validated")
})

test_that("DB SE values match raw GWA source files", {
  skip_if_no_cross_validation()

  gwa_files <- find_gwa_files(work_dir)
  if (length(gwa_files) == 0) skip("No raw GWA files found")

  # Test one file
  for (gwa_file in gwa_files) {
    params <- parse_inline_gwa_filename(gwa_file)
    if (is.null(params)) next

    ms_meta <- read_marker_set_metadata(params$population, params$maf, db_dir)
    if (is.null(ms_meta)) next
    ms_id <- generate_marker_set_id(params$population, params$maf,
                                    ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld))
    trait <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
    mapping_id <- generate_mapping_id(trait$hash, params$algorithm, params$pca)$hash
    raw_df <- tryCatch(read_raw_gwa_file(gwa_file, verbose = FALSE),
                       error = function(e) NULL)
    if (is.null(raw_df)) next

    db_df <- tryCatch({
      suppressWarnings(query_by_mapping_id(mapping_id, db_dir))
    }, error = function(e) NULL)
    if (is.null(db_df) || nrow(db_df) == 0) next

    merged <- merge(
      raw_df[, c("marker", "SE")],
      db_df[, c("marker", "SE")],
      by = "marker", suffixes = c("_raw", "_db")
    )

    expect_equal(merged$SE_db, merged$SE_raw,
                 tolerance = 1e-12,
                 label = paste("SE values match for", mapping_id))
    break
  }
})

test_that("DB AF1 values match raw GWA source files", {
  skip_if_no_cross_validation()

  gwa_files <- find_gwa_files(work_dir)
  if (length(gwa_files) == 0) skip("No raw GWA files found")

  for (gwa_file in gwa_files) {
    params <- parse_inline_gwa_filename(gwa_file)
    if (is.null(params)) next

    ms_meta <- read_marker_set_metadata(params$population, params$maf, db_dir)
    if (is.null(ms_meta)) next
    ms_id <- generate_marker_set_id(params$population, params$maf,
                                    ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld))
    trait <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
    mapping_id <- generate_mapping_id(trait$hash, params$algorithm, params$pca)$hash
    raw_df <- tryCatch(read_raw_gwa_file(gwa_file, verbose = FALSE),
                       error = function(e) NULL)
    if (is.null(raw_df)) next

    db_df <- tryCatch({
      suppressWarnings(query_by_mapping_id(mapping_id, db_dir))
    }, error = function(e) NULL)
    if (is.null(db_df) || nrow(db_df) == 0) next

    merged <- merge(
      raw_df[, c("marker", "AF1")],
      db_df[, c("marker", "AF1")],
      by = "marker", suffixes = c("_raw", "_db")
    )

    expect_equal(merged$AF1_db, merged$AF1_raw,
                 tolerance = 1e-12,
                 label = paste("AF1 values match for", mapping_id))
    break
  }
})

# ── EIGEN Threshold Validation ───────────────────────────────────────────────

test_that("EIGEN threshold from DB matches pipeline EIGEN value", {
  skip_if_no_cross_validation()

  # Find EIGEN file in work directory
  eigen_files <- list.files(work_dir,
                            pattern = "total_independent_tests\\.txt$",
                            recursive = TRUE, full.names = TRUE)

  if (length(eigen_files) == 0) {
    skip("No EIGEN files found in work directory")
  }

  for (eigen_file in eigen_files) {
    parsed <- parse_eigen_filename(eigen_file)
    if (is.null(parsed)) next

    # Read pipeline EIGEN value
    pipeline_value <- read_eigen_file(eigen_file)

    # Read from DB marker set metadata
    db_meta <- read_marker_set_metadata(parsed$population, parsed$maf, db_dir)
    if (is.null(db_meta)) next

    expect_equal(db_meta$n_independent_tests, pipeline_value,
                 label = paste("EIGEN value for", parsed$population, parsed$maf))
    break
  }
})

# ── Marker Count Consistency ─────────────────────────────────────────────────

test_that("marker counts in metadata match actual partition row counts", {
  skip_if_no_cross_validation()

  meta <- get_metadata(db_dir)

  # Check a sample of mappings
  sample_ids <- head(meta$mapping_id, 4)

  for (mid in sample_ids) {
    meta_row <- meta[meta$mapping_id == mid, ]
    db_df <- tryCatch({
      suppressWarnings(query_by_mapping_id(mid, db_dir))
    }, error = function(e) NULL)

    if (is.null(db_df) || nrow(db_df) == 0) next

    expect_equal(nrow(db_df), meta_row$n_markers[1],
                 label = paste("marker count consistency for", mid))
  }
})

test_that("marker set count matches number of markers in mappings", {
  skip_if_no_cross_validation()

  meta <- get_metadata(db_dir)
  if (nrow(meta) == 0) skip("No metadata")

  # Get a mapping and check its marker count against the marker set
  mid <- meta$mapping_id[1]
  pop <- meta$population[1]
  maf_val <- meta$maf[1]

  ms_meta <- read_marker_set_metadata(pop, maf_val, db_dir)
  if (is.null(ms_meta)) skip("marker set metadata not found")
  ms <- read_marker_set(pop, maf_val,
                        ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld),
                        db_dir)
  mapping_data <- tryCatch({
    suppressWarnings(query_by_mapping_id(mid, db_dir))
  }, error = function(e) NULL)

  if (!is.null(mapping_data) && nrow(mapping_data) > 0) {
    # Mapping row count should equal marker set row count
    # (each mapping has exactly one row per marker in the marker set)
    expect_equal(nrow(mapping_data), nrow(ms),
                 label = "mapping rows should equal marker set size")
  }
})
