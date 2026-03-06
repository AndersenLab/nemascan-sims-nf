test_that("mappings_schema excludes N and log10p, includes AF1", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_true("AF1" %in% field_names)
  expect_false("N" %in% field_names)
  expect_false("log10p" %in% field_names)
})

test_that("mappings_schema has exactly 8 fields", {
  s <- mappings_schema()
  expect_equal(length(names(s)), 8L)
})

test_that("mappings_schema contains FK columns marker_set_id, trait_id, mapping_id", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_true("marker_set_id" %in% field_names)
  expect_true("trait_id" %in% field_names)
  expect_true("mapping_id" %in% field_names)
})

test_that("mappings_schema does not contain redundant per-row constant columns", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_false("nqtl" %in% field_names)
  expect_false("algorithm" %in% field_names)
  expect_false("population" %in% field_names)
  expect_false("trait" %in% field_names)
  expect_false("maf" %in% field_names)
  expect_false("h2" %in% field_names)
  expect_false("effect" %in% field_names)
  expect_false("pca" %in% field_names)
  expect_false("rep" %in% field_names)
})

test_that("metadata_schema does not contain trait column", {
  s <- metadata_schema()
  expect_false("trait" %in% names(s))
})

test_that("metadata_schema contains hash FK columns", {
  s <- metadata_schema()
  field_names <- names(s)
  expect_true("mapping_hash_string" %in% field_names)
  expect_true("trait_id" %in% field_names)
  expect_true("marker_set_id" %in% field_names)
  expect_true("hash_schema_version" %in% field_names)
})

test_that("trait_metadata_schema includes sim_seed field", {
  s <- trait_metadata_schema()
  expect_true("sim_seed" %in% names(s),
              label = "sim_seed present in trait_metadata_schema")
})

test_that("trait_metadata_schema sim_seed is int32 type", {
  s <- trait_metadata_schema()
  expect_equal(s$GetFieldByName("sim_seed")$type, arrow::int32())
})

