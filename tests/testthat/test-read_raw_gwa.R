test_that("read_raw_gwa_file reads .fastGWA format correctly", {
  df <- read_raw_gwa_file(fixture_path("test.fastGWA"), verbose = FALSE)
  expect_equal(names(df), c("marker", "CHROM", "POS", "A1", "A2", "AF1", "BETA", "SE", "P"))
  expect_type(df$CHROM, "character")
  expect_type(df$POS, "integer")
  expect_type(df$P, "double")
  expect_equal(nrow(df), 20)
  expect_false("N" %in% names(df))
  expect_false("log10p" %in% names(df))
})

test_that("read_raw_gwa_file reads .mlma format correctly", {
  df <- read_raw_gwa_file(fixture_path("test.mlma"), verbose = FALSE)
  expect_equal(names(df), c("marker", "CHROM", "POS", "A1", "A2", "AF1", "BETA", "SE", "P"))
  expect_equal(nrow(df), 20)
})

test_that("read_raw_gwa_file errors on unknown format", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines("x\ty", tmp)
  expect_error(read_raw_gwa_file(tmp), "Unknown GWA file format")
})
