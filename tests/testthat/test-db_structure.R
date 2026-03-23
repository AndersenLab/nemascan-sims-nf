# test-db_structure.R - Integration tests for database output structure
#
# Run AFTER pipeline execution. Validates the Parquet database produced by
# DB_MIGRATION modules has correct directory structure, metadata, and is queryable.
#
# Usage: TEST_DB_DIR=/path/to/output/db Rscript tests/run_tests.R
#
# Requires TEST_DB_DIR environment variable pointing to a completed pipeline
# database output directory.

# Skip entire file if TEST_DB_DIR is not set (unit-test-only runs)
db_dir <- Sys.getenv("TEST_DB_DIR", unset = "")
skip_if_no_db <- function() {
  if (db_dir == "" || !dir.exists(db_dir)) {
    skip("TEST_DB_DIR not set or directory does not exist — skipping integration tests")
  }
}

# Expected population for test profile (override with TEST_POPULATION env var)
expected_population <- Sys.getenv("TEST_POPULATION", unset = "ce.test.200strains")

# ── Directory Structure ──────────────────────────────────────────────────────

test_that("database root directory exists and has expected subdirectories", {
  skip_if_no_db()
  expect_true(dir.exists(db_dir))
  expect_true(dir.exists(file.path(db_dir, "markers")))
  expect_true(dir.exists(file.path(db_dir, "markers", "marker_sets")))
  expect_true(dir.exists(file.path(db_dir, "markers", "genotypes")))
  expect_true(dir.exists(file.path(db_dir, "mappings")))
})

test_that("marker set parquet files exist", {
  skip_if_no_db()
  markers_dir <- file.path(db_dir, "markers", "marker_sets")
  marker_files <- list.files(markers_dir, pattern = "_markers\\.parquet$")
  expect_gt(length(marker_files), 0, label = "at least one marker set file")

  # Derive expected filename from hash via metadata lookup (5-param v=2 signature)
  ms_meta <- read_marker_set_metadata(expected_population, 0.05, db_dir)
  if (is.null(ms_meta)) skip("marker set metadata not found in TEST_DB_DIR")
  expected_hash <- generate_marker_set_id(
    expected_population, 0.05,
    ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld)
  )$hash
  expected_file <- paste0(expected_hash, "_markers.parquet")
  expect_true(
    expected_file %in% marker_files,
    label = paste("expected marker set file:", expected_file)
  )
})

test_that("mapping partitions exist in Hive-style structure", {
  skip_if_no_db()
  mappings_dir <- file.path(db_dir, "mappings")

  # Check population partition directories
  pop_dirs <- list.dirs(mappings_dir, recursive = FALSE)
  pop_names <- grep("^population=", basename(pop_dirs), value = TRUE)
  expect_gt(length(pop_names), 0, label = "at least one population partition")

  expected_pop_dir <- paste0("population=", expected_population)
  expect_true(
    expected_pop_dir %in% pop_names,
    label = paste("expected population partition:", expected_pop_dir)
  )

  # Check mapping_id sub-partitions contain data.parquet
  data_files <- list.files(mappings_dir,
    pattern = "data\\.parquet$",
    recursive = TRUE, full.names = TRUE
  )
  expect_gt(length(data_files), 0, label = "at least one mapping data.parquet")
})

test_that("metadata files exist", {
  skip_if_no_db()
  expect_true(file.exists(file.path(db_dir, "mappings_metadata.parquet")))
  expect_true(file.exists(file.path(db_dir, "marker_set_metadata.parquet")))
})

# ── Marker Set Metadata ─────────────────────────────────────────────────────

test_that("marker_set_metadata.parquet has required columns and valid data", {
  skip_if_no_db()
  ms_meta <- arrow::read_parquet(file.path(db_dir, "marker_set_metadata.parquet"))

  # Required columns
  expected_cols <- c(
    "marker_set_id", "marker_set_hash_string", "hash_schema_version",
    "population", "maf", "n_markers", "n_independent_tests",
    "eigen_source_file", "strainfile_hash", "strain_list", "created_at"
  )
  expect_true(all(expected_cols %in% names(ms_meta)),
    label = paste(
      "marker_set_metadata columns:",
      paste(setdiff(expected_cols, names(ms_meta)), collapse = ", ")
    )
  )

  # At least one record

  expect_gt(nrow(ms_meta), 0)

  # n_markers should be positive integers
  expect_true(all(ms_meta$n_markers > 0))

  # n_independent_tests should be positive (not NA for inline path with EIGEN)
  expect_true(all(!is.na(ms_meta$n_independent_tests)))
  expect_true(all(ms_meta$n_independent_tests > 0))
})

test_that("marker_set_metadata contains strainfile_hash and strain_list", {
  skip_if_no_db()
  ms_meta <- arrow::read_parquet(file.path(db_dir, "marker_set_metadata.parquet"))

  expect_true("strainfile_hash" %in% names(ms_meta),
    label = "strainfile_hash column present"
  )
  expect_true("strain_list" %in% names(ms_meta),
    label = "strain_list column present"
  )

  # strainfile_hash must be populated (required, not NA) and 64-char SHA-256 hex
  expect_true(all(!is.na(ms_meta$strainfile_hash)),
    label = "strainfile_hash must not be NA"
  )
  expect_true(all(grepl("^[0-9a-f]{64}$", ms_meta$strainfile_hash)),
    label = "strainfile_hash should be 64-char lowercase hex SHA-256"
  )

  # strain_list must be populated and non-empty
  expect_true(all(!is.na(ms_meta$strain_list)),
    label = "strain_list must not be NA"
  )
  expect_true(all(nchar(ms_meta$strain_list) > 0),
    label = "strain_list should be non-empty"
  )
})

test_that("marker set parquet has correct schema", {
  skip_if_no_db()
  ms_meta <- read_marker_set_metadata(expected_population, 0.05, db_dir)
  if (is.null(ms_meta)) skip("marker set metadata not found in TEST_DB_DIR")
  ms <- read_marker_set(expected_population, 0.05,
                        ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld),
                        db_dir)

  expected_cols <- c("marker_set_id", "marker", "CHROM", "POS", "A1", "A2")
  expect_true(all(expected_cols %in% names(ms)))
  expect_gt(nrow(ms), 0)
  # marker_set_hash_string must NOT be in the data file (metadata only)
  expect_false("marker_set_hash_string" %in% names(ms),
    label = "marker_set_hash_string absent from data file"
  )
  expect_false("AF1" %in% names(ms), label = "AF1 absent from marker data schema")
  expect_false("population" %in% names(ms), label = "population absent from marker data schema")
  expect_false("maf" %in% names(ms), label = "maf absent from marker data schema")

  # CHROM should be character
  expect_type(ms$CHROM, "character")
  # POS should be integer
  expect_type(ms$POS, "integer")
})

# ── Mappings Metadata ────────────────────────────────────────────────────────

test_that("mappings_metadata.parquet has required columns", {
  skip_if_no_db()
  meta <- get_metadata(db_dir)

  expected_cols <- c(
    "mapping_id", "population", "maf", "nqtl", "rep",
    "h2", "effect", "algorithm", "pca", "trait",
    "n_markers", "source_file", "processed_at",
    "processing_version"
  )
  present <- expected_cols[expected_cols %in% names(meta)]
  # Allow some optional columns (source_file, processed_at, processing_version
  # may come from aggregate_metadata.R which uses slightly different column set)
  core_cols <- c(
    "mapping_id", "mapping_hash_string", "trait_id",
    "marker_set_id", "hash_schema_version",
    "population", "maf", "nqtl", "rep",
    "h2", "effect", "algorithm", "pca", "n_markers"
  )
  expect_true(all(core_cols %in% names(meta)),
    label = paste(
      "missing core columns:",
      paste(setdiff(core_cols, names(meta)), collapse = ", ")
    )
  )
})

test_that("hash_string columns are present in metadata files and absent from data files", {
  skip_if_no_db()

  # mapping_hash_string must be in mappings_metadata, absent from data partitions
  meta <- get_metadata(db_dir)
  expect_true("mapping_hash_string" %in% names(meta),
    label = "mapping_hash_string present in mappings_metadata.parquet"
  )

  data_files <- list.files(file.path(db_dir, "mappings"),
    pattern = "data\\.parquet$",
    recursive = TRUE, full.names = TRUE
  )
  if (length(data_files) > 0) {
    first_data <- arrow::read_parquet(data_files[1])
    expect_false("mapping_hash_string" %in% names(first_data),
      label = "mapping_hash_string absent from mapping data.parquet"
    )
    expect_false("mapping_id" %in% names(first_data),
      label = "mapping_id absent from mapping data.parquet (partition key only)"
    )
  }
})

test_that("mappings_metadata has correct row count for test profile", {
  skip_if_no_db()
  meta <- get_metadata(db_dir)

  # For -profile test (1 group, 1 maf, 1 nqtl=5, 1 h2=0.8, 1 rep, 1 effect):
  # 2 modes (inbred, loco) × 2 types (pca, nopca) = 4 mappings
  expected_count <- as.integer(Sys.getenv("TEST_EXPECTED_MAPPINGS", unset = "4"))
  expect_equal(nrow(meta), expected_count,
    label = paste("expected", expected_count, "mappings, got", nrow(meta))
  )
})

test_that("algorithm values are valid", {
  skip_if_no_db()
  meta <- get_metadata(db_dir)

  valid_algorithms <- c("inbred", "loco")
  expect_true(all(meta$algorithm %in% valid_algorithms),
    label = paste(
      "unexpected algorithms:",
      paste(setdiff(unique(meta$algorithm), valid_algorithms),
        collapse = ", "
      )
    )
  )

  # Both algorithms should be present for a full run
  expect_setequal(unique(meta$algorithm), valid_algorithms)
})

test_that("PCA flags are correct", {
  skip_if_no_db()
  meta <- get_metadata(db_dir)

  expect_type(meta$pca, "logical")
  # Both PCA and noPCA mappings should exist
  expect_true(any(meta$pca == TRUE), label = "at least one PCA mapping")
  expect_true(any(meta$pca == FALSE), label = "at least one noPCA mapping")
})

test_that("mapping IDs are 20-character lowercase hex hashes", {
  skip_if_no_db()
  meta <- get_metadata(db_dir)

  expect_true(all(grepl("^[0-9a-f]{20}$", meta$mapping_id)),
    label = "all mapping IDs should be 20-char lowercase hex"
  )
})

# ── Database Queryability ────────────────────────────────────────────────────

test_that("database can be opened and queried via open_mapping_db()", {
  skip_if_no_db()
  con <- open_mapping_db(db_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Should have markers view
  tables <- DBI::dbListTables(con)
  expect_true("markers" %in% tables)
  expect_true("mappings" %in% tables)

  # Basic query should work
  result <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM mappings")
  expect_gt(result$n[1], 0)
})

test_that("db_stats() returns valid statistics", {
  skip_if_no_db()
  stats <- db_stats(db_dir)

  expect_true(stats$exists)
  expect_gt(stats$n_marker_sets, 0)
  expect_gt(stats$n_populations, 0)
  expect_gt(stats$n_mappings, 0)
  expect_gt(stats$n_partitions, 0)
  expect_gt(stats$total_size_mb, 0)
})

test_that("query_mapping_data() returns joined data with marker info", {
  skip_if_no_db()
  data <- query_mapping_data(population = expected_population, base_dir = db_dir)

  expect_gt(nrow(data), 0)
  # Should have mapping columns
  expect_true("P" %in% names(data))
  expect_true("BETA" %in% names(data))
  # Should have joined marker columns
  expect_true("CHROM" %in% names(data))
  expect_true("POS" %in% names(data))
})

test_that("get_threshold_params() returns valid thresholds", {
  skip_if_no_db()
  tp <- get_threshold_params(expected_population, 0.05, base_dir = db_dir)

  expect_gt(tp$n_markers, 0)
  expect_gt(tp$bf_threshold, 0)
  expect_gt(tp$eigen_threshold, 0)
  # EIGEN threshold should be less strict than BF (fewer independent tests than markers)
  expect_lt(tp$eigen_threshold, tp$bf_threshold)
  expect_gt(tp$n_independent_tests, 0)
  expect_lt(tp$n_independent_tests, tp$n_markers)
})

# ── Mapping Data Integrity ───────────────────────────────────────────────────

test_that("all mapping partitions have non-empty data", {
  skip_if_no_db()
  data_files <- list.files(file.path(db_dir, "mappings"),
    pattern = "data\\.parquet$",
    recursive = TRUE, full.names = TRUE
  )

  for (f in data_files) {
    df <- arrow::read_parquet(f)
    expect_gt(nrow(df), 0,
      label = paste("non-empty parquet:", basename(dirname(f)))
    )
  }
})

test_that("P values are in valid range (0, 1]", {
  skip_if_no_db()
  con <- open_mapping_db(db_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  result <- DBI::dbGetQuery(con, "
    SELECT
      COUNT(*) as total,
      SUM(CASE WHEN P <= 0 THEN 1 ELSE 0 END) as invalid_low,
      SUM(CASE WHEN P > 1 THEN 1 ELSE 0 END) as invalid_high,
      SUM(CASE WHEN P IS NULL THEN 1 ELSE 0 END) as null_count
    FROM mappings
  ")

  expect_equal(result$invalid_low[1], 0, label = "no P <= 0")
  expect_equal(result$invalid_high[1], 0, label = "no P > 1")
  expect_equal(result$null_count[1], 0, label = "no NULL P values")
})

test_that("var.exp is absent or NA for inline-path mappings", {
  skip_if_no_db()
  con <- open_mapping_db(db_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # For inline path (nemascan-sims-nf integration), var.exp is not available
  # (no genotype matrix during raw GWA import). The column may be absent entirely
  # or present with all NA values depending on the write path.
  cols <- DBI::dbGetQuery(con, "SELECT column_name FROM information_schema.columns WHERE table_name = 'mappings'")
  if ("var.exp" %in% cols$column_name) {
    result <- DBI::dbGetQuery(con, "
      SELECT COUNT(*) as non_null_varexp
      FROM mappings
      WHERE \"var.exp\" IS NOT NULL
    ")
    expect_equal(result$non_null_varexp[1], 0,
      label = "var.exp should be NA for inline-path databases"
    )
  } else {
    succeed("var.exp column absent — expected for inline-path databases")
  }
})
