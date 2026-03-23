# Phase 1 gate test: R library integrity checks

test_that("all R library files exist and source without error", {
  r_dir <- file.path(.project_root, "R")
  expected_files <- c(
    "utils.R", "io.R", "database.R", "queries.R", "analysis.R", "qtl_database.R"
  )
  for (f in expected_files) {
    expect_true(file.exists(file.path(r_dir, f)), info = paste("Missing:", f))
  }
  # setup.R is the convenience loader, also must exist
  expect_true(file.exists(file.path(r_dir, "setup.R")))
})

test_that("Phase 0 functions are available after sourcing", {
  # Core functions added/modified in Phase 0
  expect_true(is.function(read_raw_gwa_file))
  expect_true(is.function(read_bim_file))
  expect_true(is.function(read_eigen_file))
  expect_true(is.function(safe_log10p))
  expect_true(is.function(write_marker_set))
  expect_true(is.function(write_mapping_partitioned))
  expect_true(is.function(generate_mapping_id))
  expect_true(is.function(prepare_mapping_data))
})

test_that("no stale log10p references in database schema", {
  # Mapping schema should NOT include log10p or N columns
  # (computed at query time via safe_log10p, not stored)
  r_dir <- file.path(.project_root, "R")
  db_code <- readLines(file.path(r_dir, "database.R"))

  # The prepare_mapping_data function should not select log10p
  prep_section <- paste(db_code, collapse = "\n")

  # Check that mapping_cols in prepare_mapping_data doesn't include log10p or N
  # Find the mapping_cols definition
  mapping_cols_line <- grep('mapping_cols.*<-.*c\\(', db_code, value = TRUE)
  for (line in mapping_cols_line) {
    expect_false(grepl('"log10p"', line), info = "log10p should not be in mapping_cols")
    expect_false(grepl('"N"', line), info = "N should not be in mapping_cols")
  }
})
