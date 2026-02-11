test_that("read_bim_file parses .bim format correctly", {
  df <- read_bim_file(fixture_path("test.bim"))
  expect_equal(names(df), c("marker", "CHROM", "POS", "A1", "A2", "AF1"))
  expect_true(all(is.na(df$AF1)))
  expect_equal(nrow(df), 50)
  # marker should be "CHROM:POS"
  expect_equal(df$marker[1], paste(df$CHROM[1], df$POS[1], sep = ":"))
})
