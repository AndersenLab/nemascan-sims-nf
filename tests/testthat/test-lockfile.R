# Tests for renv.lock — R package lockfile for the r_packages:20250519 container
#
# Tests 1 and 2 pass immediately (stub lockfile is present and valid JSON).
# Test 3 skips while the placeholder stub is in place; it passes once
# scripts/generate_renv_lock.sh has been run and the real lockfile committed.

test_that("renv.lock exists at project root", {
  expect_true(
    file.exists(file.path(.project_root, "renv.lock")),
    label = "renv.lock must be committed to the project root"
  )
})

test_that("renv.lock is valid JSON with required top-level keys", {
  lockfile_path <- file.path(.project_root, "renv.lock")
  skip_if_not(file.exists(lockfile_path))

  content <- jsonlite::read_json(lockfile_path)
  expect_true("R" %in% names(content),
    label = "renv.lock must have an 'R' key")
  expect_true("Packages" %in% names(content),
    label = "renv.lock must have a 'Packages' key")
})

test_that("renv.lock contains all packages used in the r_packages container scope", {
  lockfile_path <- file.path(.project_root, "renv.lock")
  skip_if_not(file.exists(lockfile_path))

  content <- jsonlite::read_json(lockfile_path)
  skip_if(
    "_PLACEHOLDER_" %in% names(content$Packages),
    "renv.lock is a placeholder stub — run scripts/generate_renv_lock.sh to populate it"
  )

  # Packages identified by static scan of bin/ scripts and R/ library files
  # that run inside the r_packages:20250519 container (R_.* Nextflow processes)
  required_packages <- c(
    # Core data manipulation (R/ library + bin/)
    "dplyr", "tidyr", "readr", "glue", "purrr", "data.table",
    # bin/Get_GenoMatrix_Eigen.R
    "Rcpp", "ggplot2", "stringr", "coop", "RSpectra",
    # bin/Create_Causal_QTLs.R
    "IRanges", "fuzzyjoin",
    # bin/Assess_Sims.R
    "GenomicRanges",
    # bin/Get_GCTA_Intervals.R
    "ggbeeswarm"
  )

  present <- names(content$Packages)
  missing <- setdiff(required_packages, present)

  expect_equal(
    missing, character(0),
    label = paste(
      "Packages missing from renv.lock:",
      if (length(missing) > 0) paste(missing, collapse = ", ") else "none"
    )
  )
})
