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

# Helper: look up mapping_id from metadata by filename-parsed parameters.
# Avoids reconstructing the full hash chain (which requires cv_maf_effective
# and cv_ld values that are not encoded in the GWA filename).
lookup_mapping_id <- function(params, meta) {
  matches <- meta[
    meta$population == params$population &
    meta$maf       == params$maf &
    meta$nqtl      == params$nqtl &
    meta$rep       == params$rep &
    meta$h2        == params$h2 &
    meta$effect    == params$effect &
    meta$algorithm == params$algorithm &
    meta$pca       == params$pca,
  ]
  if (nrow(matches) != 1) return(NULL)
  matches$mapping_id[1]
}

# ── Round-trip Value Accuracy ────────────────────────────────────────────────

test_that("DB mapping P and BETA values match raw GWA source files", {
  skip_if_no_cross_validation()

  gwa_files <- find_gwa_files(work_dir)
  if (length(gwa_files) == 0) {
    skip("No raw GWA files found in TEST_WORK_DIR")
  }

  meta <- get_metadata(db_dir)

  # Test at least one file of each type (fastGWA and mlma)
  tested <- 0
  for (gwa_file in gwa_files) {
    params <- parse_inline_gwa_filename(gwa_file)
    if (is.null(params)) next

    mapping_id <- lookup_mapping_id(params, meta)
    if (is.null(mapping_id)) next

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

  meta <- get_metadata(db_dir)

  # Test one file
  for (gwa_file in gwa_files) {
    params <- parse_inline_gwa_filename(gwa_file)
    if (is.null(params)) next

    mapping_id <- lookup_mapping_id(params, meta)
    if (is.null(mapping_id)) next

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

  meta <- get_metadata(db_dir)

  for (gwa_file in gwa_files) {
    params <- parse_inline_gwa_filename(gwa_file)
    if (is.null(params)) next

    mapping_id <- lookup_mapping_id(params, meta)
    if (is.null(mapping_id)) next

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

    # Skip broken symlinks (e.g. after rsync from HPC without --copy-unsafe-links)
    if (!file.exists(eigen_file)) next

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

# ── GWA File Column Assertions ───────────────────────────────────────────────

test_that("read_raw_gwa_file() returns normalized columns for inbred and loco formats", {
  skip_if_no_cross_validation()

  gwa_files <- find_gwa_files(work_dir)
  if (length(gwa_files) == 0) skip("no raw GWA files found")

  fastgwa_files <- gwa_files[grepl("\\.fastGWA$", gwa_files)]
  mlma_files    <- gwa_files[grepl("\\.mlma$", gwa_files)]

  # Normalized column set that read_raw_gwa_file() must produce for both formats
  expected_cols <- c("marker", "CHROM", "POS", "A1", "A2", "AF1", "BETA", "SE", "P")

  tested <- 0L

  if (length(fastgwa_files) > 0) {
    raw <- tryCatch(
      read_raw_gwa_file(fastgwa_files[1], verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(raw)) {
      missing <- setdiff(expected_cols, names(raw))
      expect_equal(length(missing), 0,
        label = paste("fastGWA missing normalized cols:", paste(missing, collapse = ", "))
      )
      expect_type(raw$P, "double")
      expect_gt(nrow(raw), 0, label = "fastGWA file is non-empty")
      tested <- tested + 1L
    }
  }

  if (length(mlma_files) > 0) {
    raw <- tryCatch(
      read_raw_gwa_file(mlma_files[1], verbose = FALSE),
      error = function(e) NULL
    )
    if (!is.null(raw)) {
      missing <- setdiff(expected_cols, names(raw))
      expect_equal(length(missing), 0,
        label = paste("mlma missing normalized cols:", paste(missing, collapse = ", "))
      )
      expect_type(raw$P, "double")
      expect_gt(nrow(raw), 0, label = "mlma file is non-empty")
      tested <- tested + 1L
    }
  }

  expect_gt(tested, 0L, label = "at least one GWA format was tested")
})

# ── Full Mapping Coverage ─────────────────────────────────────────────────────

test_that("query_by_mapping_id() returns non-empty result for every mapping_id in metadata", {
  skip_if_no_cross_validation()

  meta <- get_metadata(db_dir)
  expect_gt(nrow(meta), 0, label = "metadata has at least one mapping")

  for (mid in meta$mapping_id) {
    db_df <- tryCatch(
      suppressWarnings(query_by_mapping_id(mid, db_dir)),
      error = function(e) NULL
    )
    expect_false(is.null(db_df),
      label = paste("mapping queryable without error:", mid)
    )
    if (!is.null(db_df)) {
      expect_gt(nrow(db_df), 0,
        label = paste("non-empty mapping result:", mid)
      )
    }
  }
})

# ── Loco/Inbred Marker Count Invariant ───────────────────────────────────────

test_that("loco and inbred GWA have equal marker counts for same marker set", {
  skip_if_no_cross_validation()

  meta <- get_metadata(db_dir)
  if (!all(c("inbred", "loco") %in% unique(meta$algorithm))) {
    skip("both inbred and loco mappings required for this test")
  }

  # Check all unique populations — not just [1,].
  # After the LOCO marker count fix, both inbred and loco store n_markers equal
  # to the full LD-pruned marker set size in mappings metadata. The BF threshold
  # denominator in analyze_qtl.R uses nrow(mapping_data) at runtime regardless of
  # this stored value.
  populations <- unique(meta$population)
  any_checked <- FALSE

  for (pop in populations) {
    inbred_pop <- meta[meta$algorithm == "inbred" & meta$population == pop, ]
    loco_pop   <- meta[meta$algorithm == "loco"   & meta$population == pop, ]
    if (nrow(inbred_pop) == 0 || nrow(loco_pop) == 0) next

    # Use the first MAF available for this population
    maf_val <- inbred_pop$maf[1]
    inbred_row <- inbred_pop[inbred_pop$maf == maf_val, ][1, ]
    loco_row   <- loco_pop[loco_pop$maf == maf_val, ][1, ]
    if (is.na(loco_row$n_markers)) next

    expect_equal(
      loco_row$n_markers, inbred_row$n_markers,
      label = sprintf(
        "loco n_markers (%d) == inbred n_markers (%d) for population=%s maf=%s",
        loco_row$n_markers, inbred_row$n_markers, pop, maf_val
      )
    )
    any_checked <- TRUE
  }

  if (!any_checked) skip("no matching inbred+loco population pairs found")
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

# ── Direct vs Legacy Mapping Read Parity ─────────────────────────────────────
#
# query_mapping_direct() replaces query_for_threshold_analysis() in the
# analyze_qtl / assess_sims hot paths to avoid the union_by_name FD blow-up
# (see issues/new-db-qtl-analysis-too-many-files/rca.md). This test asserts
# the two code paths return content-equivalent dataframes.

test_that("query_mapping_direct() matches query_for_threshold_analysis() schema and content", {
  skip_if_no_cross_validation()

  meta <- get_metadata(db_dir)
  skip_if(nrow(meta) == 0, "mappings metadata is empty")

  # Find a mapping_id whose legacy query returns rows (defensively handle
  # empty partitions — unlikely in a healthy DB but keeps the test robust).
  checked <- 0L
  for (i in seq_len(min(nrow(meta), 4L))) {
    row        <- meta[i, ]
    mapping_id <- as.character(row$mapping_id)
    population <- as.character(row$population)
    maf        <- as.numeric(row$maf)

    ms_meta <- tryCatch(
      read_marker_set_metadata(population, maf, db_dir),
      error = function(e) NULL
    )
    if (is.null(ms_meta)) next

    legacy_df <- tryCatch(
      suppressWarnings(query_for_threshold_analysis(mapping_id, db_dir)),
      error = function(e) NULL
    )
    if (is.null(legacy_df) || nrow(legacy_df) == 0) next

    direct_df <- query_mapping_direct(mapping_id, population, ms_meta, db_dir)

    # Schema parity — column names and order must match exactly
    expect_equal(names(direct_df), names(legacy_df),
                 label = paste("column names for", mapping_id))
    expect_equal(nrow(direct_df), nrow(legacy_df),
                 label = paste("row count for", mapping_id))

    # Content parity (sort both to eliminate any residual ordering drift)
    legacy_sorted <- legacy_df %>% dplyr::arrange(CHROM, POS, marker)
    direct_sorted <- direct_df %>% dplyr::arrange(CHROM, POS, marker)

    expect_equal(direct_sorted$marker, legacy_sorted$marker,
                 label = paste("marker values for", mapping_id))
    expect_equal(direct_sorted$mapping_id, legacy_sorted$mapping_id,
                 label = paste("mapping_id values for", mapping_id))
    expect_equal(direct_sorted$population, legacy_sorted$population,
                 label = paste("population values for", mapping_id))
    expect_equal(as.character(direct_sorted$CHROM), as.character(legacy_sorted$CHROM),
                 label = paste("CHROM values for", mapping_id))
    expect_equal(as.integer(direct_sorted$POS), as.integer(legacy_sorted$POS),
                 label = paste("POS values for", mapping_id))
    expect_equal(direct_sorted$P, legacy_sorted$P,
                 tolerance = 1e-12,
                 label = paste("P values for", mapping_id))
    expect_equal(direct_sorted$BETA, legacy_sorted$BETA,
                 tolerance = 1e-12,
                 label = paste("BETA values for", mapping_id))
    expect_equal(direct_sorted$SE, legacy_sorted$SE,
                 tolerance = 1e-12,
                 label = paste("SE values for", mapping_id))
    expect_equal(direct_sorted$AF1, legacy_sorted$AF1,
                 tolerance = 1e-12,
                 label = paste("AF1 values for", mapping_id))
    expect_equal(as.numeric(direct_sorted$maf), as.numeric(legacy_sorted$maf),
                 label = paste("maf values for", mapping_id))

    # var.exp is NA_real_ in both paths under the Phase 5 DB (column removed
    # from mappings schema; legacy path takes the 'NULL AS "var.exp"' fallback
    # at R/queries.R:318, direct path sets NA_real_ explicitly).
    expect_true(all(is.na(direct_sorted$`var.exp`)),
                label = paste("direct var.exp all NA for", mapping_id))
    expect_true(all(is.na(legacy_sorted$`var.exp`)),
                label = paste("legacy var.exp all NA for", mapping_id))

    checked <- checked + 1L
    break  # one successful parity check is enough
  }

  expect_gt(checked, 0L,
            label = "at least one mapping was parity-checked across direct and legacy paths")
})
