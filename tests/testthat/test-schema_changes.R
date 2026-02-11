test_that("mappings_schema excludes N and log10p, includes AF1", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_true("AF1" %in% field_names)
  expect_false("N" %in% field_names)
  expect_false("log10p" %in% field_names)
})

test_that("prepare_mapping_data handles missing var.exp gracefully", {
  df <- read_raw_gwa_file(fixture_path("test.fastGWA"), verbose = FALSE)
  params <- list(nqtl = 5L, rep = 1L, h2 = 0.8, maf = 0.05, effect = "gamma",
                 population = "test_pop", algorithm = "INBRED", pca = TRUE)
  result <- prepare_mapping_data(df, params)
  # var.exp should not be in the result (not in raw GWA, not in available_cols)
  # OR if it is added as NA, that's also acceptable
  if ("var.exp" %in% names(result)) {
    expect_true(all(is.na(result$var.exp)))
  }
  # AF1 should be present
  expect_true("AF1" %in% names(result))
  # N and log10p should NOT be present
  expect_false("N" %in% names(result))
  expect_false("log10p" %in% names(result))
})

test_that("generate_mapping_id includes explicit PCA/noPCA suffix", {
  params_pca  <- list(nqtl = 5, rep = 1, h2 = 0.8, maf = 0.05, effect = "gamma",
                      population = "pop", algorithm = "INBRED", pca = TRUE)
  params_nopca <- modifyList(params_pca, list(pca = FALSE))
  id_pca  <- generate_mapping_id(params_pca)
  id_nopca <- generate_mapping_id(params_nopca)
  expect_match(id_pca, "_PCA$")
  expect_match(id_nopca, "_noPCA$")
  expect_false(id_pca == id_nopca)
})
