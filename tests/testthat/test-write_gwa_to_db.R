test_that("fastGWA data writes to database with correct schema and values", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  df <- read_raw_gwa_file(fixture_path("test.fastGWA"), verbose = FALSE)
  params <- list(nqtl = 5L, rep = 1L, h2 = 0.8, maf = 0.05, effect = "gamma",
                 population = "test_pop", algorithm = "inbred", pca = TRUE)
  ms_id   <- generate_marker_set_id(params$population, params$maf)
  trait   <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
  mapping <- generate_mapping_id(trait$hash, params$algorithm, params$pca)

  write_mapping_partitioned(df, params, ms_id, trait, db_dir)

  # Verify parquet was created
  partition_path <- get_partition_path("test_pop", mapping$hash, db_dir)
  parquet_path <- file.path(partition_path, "data.parquet")
  expect_true(file.exists(parquet_path))

  # Read back and verify
  result <- arrow::read_parquet(parquet_path)
  expect_equal(nrow(result), nrow(df))
  expect_true("AF1" %in% names(result))
  expect_false("N" %in% names(result))
  expect_false("log10p" %in% names(result))

  # Verify measurement values are preserved
  expect_equal(result$P, df$P)
  expect_equal(result$BETA, df$BETA)
  expect_equal(result$AF1, df$AF1)

  # Verify FK columns are present
  expect_true("marker_set_id" %in% names(result))
  expect_true("trait_id" %in% names(result))
  expect_match(unique(result$mapping_id), "^[0-9a-f]{20}$")
  expect_false("mapping_hash_string" %in% names(result))

  # Verify redundant params not stored per-row
  expect_false("nqtl" %in% names(result))
  expect_false("algorithm" %in% names(result))
  expect_false("trait" %in% names(result))
})

test_that("mlma data writes with correct LOCO metadata", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  df <- read_raw_gwa_file(fixture_path("test.mlma"), verbose = FALSE)
  params <- list(nqtl = 5L, rep = 1L, h2 = 0.8, maf = 0.05, effect = "gamma",
                 population = "test_pop", algorithm = "loco", pca = FALSE)
  ms_id   <- generate_marker_set_id(params$population, params$maf)
  trait   <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
  mapping <- generate_mapping_id(trait$hash, params$algorithm, params$pca)

  write_mapping_partitioned(df, params, ms_id, trait, db_dir)

  partition_path <- get_partition_path("test_pop", mapping$hash, db_dir)
  result <- arrow::read_parquet(file.path(partition_path, "data.parquet"))

  expect_match(unique(result$mapping_id), "^[0-9a-f]{20}$")
})

test_that("PCA=TRUE and PCA=FALSE produce different mapping_ids", {
  params_base <- list(nqtl = 5L, rep = 1L, h2 = 0.8, maf = 0.05, effect = "gamma",
                      population = "test_pop", algorithm = "inbred")
  ms_id  <- generate_marker_set_id(params_base$population, params_base$maf)
  trait  <- generate_trait_id(ms_id$hash, params_base$nqtl, params_base$effect,
                               params_base$rep, params_base$h2)

  mapping_pca    <- generate_mapping_id(trait$hash, params_base$algorithm, TRUE)
  mapping_nopca  <- generate_mapping_id(trait$hash, params_base$algorithm, FALSE)

  expect_false(mapping_pca$hash == mapping_nopca$hash)
})

test_that("marker set creation from bim file works end-to-end", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  marker_df <- read_bim_file(fixture_path("test.bim"))
  write_marker_set(marker_df, "test_pop", 0.05, db_dir,
                   overwrite = TRUE,
                   n_independent_tests = 1234)

  # Verify marker set
  ms <- read_marker_set("test_pop", 0.05, db_dir)
  expect_equal(nrow(ms), 50)
  expect_true("marker_set_id" %in% names(ms))
  expect_false("AF1" %in% names(ms))
  expect_false("population" %in% names(ms))
  expect_false("maf" %in% names(ms))

  # Verify metadata
  meta <- read_marker_set_metadata("test_pop", 0.05, db_dir)
  expect_equal(meta$n_markers, 50)
  expect_equal(meta$n_independent_tests, 1234)
  expect_equal(meta$marker_set_id, generate_marker_set_id("test_pop", 0.05)$hash)

  # Verify threshold calculation
  tp <- get_threshold_params("test_pop", 0.05, base_dir = db_dir)
  expect_equal(tp$n_markers, 50)
  expect_equal(tp$bf_threshold, -log10(0.05 / 50))
  expect_equal(tp$eigen_threshold, -log10(0.05 / 1234))
})

test_that("overwrite = TRUE replaces existing marker set on retry", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  marker_df <- read_bim_file(fixture_path("test.bim"))

  write_marker_set(marker_df, "test_pop", 0.05, db_dir,
                   overwrite = TRUE, n_independent_tests = 1234)
  meta1 <- read_marker_set_metadata("test_pop", 0.05, db_dir)
  expect_equal(meta1$n_independent_tests, 1234)

  write_marker_set(marker_df, "test_pop", 0.05, db_dir,
                   overwrite = TRUE, n_independent_tests = 5678)
  meta2 <- read_marker_set_metadata("test_pop", 0.05, db_dir)
  expect_equal(meta2$n_independent_tests, 5678)
})
