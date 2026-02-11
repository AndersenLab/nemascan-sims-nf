test_that("safe_log10p matches stored log10p from processed mappings", {
  df <- data.table::fread(fixture_path("test_mapping.tsv"))
  computed <- safe_log10p(df$P)
  # Values should match within floating-point tolerance
  expect_equal(computed, df$log10p, tolerance = 1e-10)
})

test_that("significance calls are identical via computed vs stored log10p", {
  df <- data.table::fread(fixture_path("test_mapping.tsv"))
  n_eigen <- as.numeric(readLines(fixture_path("test_eigen.txt"), n = 1))
  n_markers <- nrow(df)

  # BF threshold
  bf_thresh <- -log10(0.05 / n_markers)
  sig_stored  <- df$log10p >= bf_thresh
  sig_computed <- safe_log10p(df$P) >= bf_thresh
  expect_identical(sig_stored, sig_computed)

  # EIGEN threshold
  eigen_thresh <- -log10(0.05 / n_eigen)
  sig_stored_e  <- df$log10p >= eigen_thresh
  sig_computed_e <- safe_log10p(df$P) >= eigen_thresh
  expect_identical(sig_stored_e, sig_computed_e)
})
