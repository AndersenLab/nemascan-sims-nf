# test-genotype_storage.R - Unit tests for genotype matrix storage functions
#
# Tests write_genotype_matrix(), genotype_matrix_exists(), and read_genotype_matrix()
# using a small synthetic fixture derived from real C. elegans genotype data.
#
# Fixture: tests/fixtures/test_genotype_matrix.tsv
#   - 5 markers (CHROM 1-5), 3 strains (BRC20263, JU397, XZ1734)
#   - One NA value (CHROM 3 / JU397) to verify NA preservation
#   - Subsetted from ce.96.allout15_irrepressible.grosbeak_0.05_Genotype_Matrix.tsv


# ── Round-trip write / read ───────────────────────────────────────────────────

test_that("write_genotype_matrix produces correct long-format row count", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)

  result <- read_genotype_matrix("test_pop", 0.05, db_dir)
  # 5 markers × 3 strains = 15 rows
  expect_equal(nrow(result), 15L)
})

test_that("write_genotype_matrix output has correct schema columns", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  result <- read_genotype_matrix("test_pop", 0.05, db_dir)

  expect_setequal(names(result), c("CHROM", "POS", "strain", "allele"))
})

test_that("write_genotype_matrix drops REF and ALT columns", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  result <- read_genotype_matrix("test_pop", 0.05, db_dir)

  expect_false("REF" %in% names(result))
  expect_false("ALT" %in% names(result))
})

test_that("write_genotype_matrix preserves allele values including NA", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  result <- read_genotype_matrix("test_pop", 0.05, db_dir)

  # All allele values should be -1, 1, or NA
  non_na_vals <- result$allele[!is.na(result$allele)]
  expect_true(all(non_na_vals %in% c(-1, 1)),
              label = "non-NA alleles are -1 or 1")

  # Exactly one NA (CHROM 3 / JU397 in the fixture)
  expect_equal(sum(is.na(result$allele)), 1L,
               label = "exactly one NA allele preserved from fixture")
})

test_that("write_genotype_matrix preserves strain names", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  result <- read_genotype_matrix("test_pop", 0.05, db_dir)

  expect_setequal(unique(result$strain), c("BRC20263", "JU397", "XZ1734"))
})


# ── Existence checks ──────────────────────────────────────────────────────────

test_that("genotype_matrix_exists returns FALSE before write, TRUE after", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  expect_false(genotype_matrix_exists("test_pop", 0.05, db_dir))

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)

  expect_true(genotype_matrix_exists("test_pop", 0.05, db_dir))
})


# ── Overwrite behavior ────────────────────────────────────────────────────────

test_that("write_genotype_matrix with overwrite = TRUE succeeds on second write", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  # Second write should not error
  expect_no_error(
    write_genotype_matrix(tsv, "test_pop", 0.05, db_dir, overwrite = TRUE)
  )

  # Data should still be correct after overwrite
  result <- read_genotype_matrix("test_pop", 0.05, db_dir)
  expect_equal(nrow(result), 15L)
})

test_that("write_genotype_matrix with overwrite = FALSE skips existing file", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  path_before <- get_genotype_matrix_path("test_pop", 0.05, db_dir)
  mtime_before <- file.mtime(path_before)

  Sys.sleep(0.01)
  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir, overwrite = FALSE)
  mtime_after <- file.mtime(path_before)

  expect_equal(mtime_before, mtime_after,
               label = "file not modified when overwrite = FALSE")
})


# ── Arrow schema enforcement ──────────────────────────────────────────────────

test_that("output Parquet has correct Arrow types", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  tsv <- fixture_path("test_genotype_matrix.tsv")

  write_genotype_matrix(tsv, "test_pop", 0.05, db_dir)
  path <- get_genotype_matrix_path("test_pop", 0.05, db_dir)

  tbl <- arrow::read_parquet(path, as_data_frame = FALSE)
  schema <- tbl$schema

  expect_equal(schema$GetFieldByName("CHROM")$type, arrow::utf8())
  expect_equal(schema$GetFieldByName("POS")$type, arrow::int32())
  expect_equal(schema$GetFieldByName("strain")$type, arrow::utf8())
  expect_equal(schema$GetFieldByName("allele")$type, arrow::float64())
})
