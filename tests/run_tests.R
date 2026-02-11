#!/usr/bin/env Rscript
# Run all tests. Usage: Rscript tests/run_tests.R
# IMPORTANT: Run from nemascan-sims-nf/ project root:
#   cd nemascan-sims-nf && Rscript tests/run_tests.R
#
# Requires: testthat package
# Exit code 0 = all pass, 1 = failures

library(testthat)

# Source all R library files
r_dir <- file.path(dirname(getwd()), "R")  # adjust if run from tests/
if (!dir.exists(r_dir)) r_dir <- file.path(getwd(), "R")
for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) {
  source(f)
}

# Run tests
results <- test_dir(
  file.path("tests", "testthat"),
  reporter = c("summary", "fail")
)

# Exit with appropriate code
if (!as.data.frame(results)$failed |> sum() == 0) quit(status = 1)
