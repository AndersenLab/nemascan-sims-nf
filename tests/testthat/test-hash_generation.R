# Golden values computed 2026-02-27 after implementing SHA-256 hash functions.
# ms_hash   <- generate_marker_set_id("ce100", 0.05)$hash   → "879ca1a75fc7fb5183e7"
# trait_hash <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash → "2aa558d92d805fb85ee8"
# map_hash   <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash → "6c4379c7d9a216d77e49"

ms_hash    <- "879ca1a75fc7fb5183e7"
trait_hash <- "2aa558d92d805fb85ee8"
map_hash   <- "6c4379c7d9a216d77e49"

# ==============================================================================
# generate_marker_set_id()
# ==============================================================================

test_that("generate_marker_set_id returns list with hash and hash_string", {
  result <- generate_marker_set_id("ce100", 0.05)
  expect_type(result, "list")
  expect_true("hash" %in% names(result))
  expect_true("hash_string" %in% names(result))
})

test_that("generate_marker_set_id hash is 20-char lowercase hex", {
  result <- generate_marker_set_id("ce100", 0.05)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_marker_set_id hash_string starts with v=1| and contains required fields", {
  result <- generate_marker_set_id("ce100", 0.05)
  expect_true(startsWith(result$hash_string, "v=1|"))
  expect_true(grepl("population=", result$hash_string))
  expect_true(grepl("maf=", result$hash_string))
})

test_that("generate_marker_set_id float serialization uses 10 decimal places", {
  result <- generate_marker_set_id("ce100", 0.05)
  expect_true(grepl("maf=0.0500000000", result$hash_string))
})

test_that("generate_marker_set_id is deterministic", {
  h1 <- generate_marker_set_id("ce100", 0.05)$hash
  h2 <- generate_marker_set_id("ce100", 0.05)$hash
  expect_equal(h1, h2)
})

test_that("generate_marker_set_id is sensitive to population", {
  h_ce100 <- generate_marker_set_id("ce100", 0.05)$hash
  h_ce96  <- generate_marker_set_id("ce96",  0.05)$hash
  expect_false(h_ce100 == h_ce96)
})

test_that("generate_marker_set_id is sensitive to maf", {
  h_05 <- generate_marker_set_id("ce100", 0.05)$hash
  h_01 <- generate_marker_set_id("ce100", 0.01)$hash
  expect_false(h_05 == h_01)
})

test_that("generate_marker_set_id normalizes population to lowercase", {
  h_upper <- generate_marker_set_id("CE100", 0.05)$hash
  h_lower <- generate_marker_set_id("ce100", 0.05)$hash
  expect_equal(h_upper, h_lower)
})

test_that("generate_marker_set_id trims whitespace from population", {
  h_spaces <- generate_marker_set_id("  ce100  ", 0.05)$hash
  h_clean  <- generate_marker_set_id("ce100", 0.05)$hash
  expect_equal(h_spaces, h_clean)
})

test_that("generate_marker_set_id canonical form stores lowercase population", {
  result <- generate_marker_set_id("CE100", 0.05)
  expect_true(grepl("population=ce100", result$hash_string))
})

test_that("generate_marker_set_id golden value", {
  expect_equal(generate_marker_set_id("ce100", 0.05)$hash, ms_hash)
})

# ==============================================================================
# generate_trait_id()
# ==============================================================================

test_that("generate_trait_id returns list with hash and hash_string", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)
  expect_type(result, "list")
  expect_true("hash" %in% names(result))
  expect_true("hash_string" %in% names(result))
})

test_that("generate_trait_id hash is 20-char lowercase hex", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_trait_id hash_string contains parent hash verbatim", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)
  expect_true(grepl(paste0("parent=", ms_hash), result$hash_string))
})

test_that("generate_trait_id float serialization uses 10 decimal places for h2", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)
  expect_true(grepl("h2=0.8000000000", result$hash_string))
})

test_that("generate_trait_id is deterministic", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash
  expect_equal(h1, h2)
})

test_that("generate_trait_id is sensitive to nqtl", {
  h5 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash
  h3 <- generate_trait_id(ms_hash, 3, "gamma", 1, 0.8)$hash
  expect_false(h5 == h3)
})

test_that("generate_trait_id is sensitive to effect", {
  h_gamma  <- generate_trait_id(ms_hash, 5, "gamma",  1, 0.8)$hash
  h_normal <- generate_trait_id(ms_hash, 5, "normal", 1, 0.8)$hash
  expect_false(h_gamma == h_normal)
})

test_that("generate_trait_id is sensitive to rep", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 2, 0.8)$hash
  expect_false(h1 == h2)
})

test_that("generate_trait_id is sensitive to h2", {
  h_08 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash
  h_04 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.4)$hash
  expect_false(h_08 == h_04)
})

test_that("generate_trait_id is sensitive to parent marker_set_hash", {
  other_ms <- generate_marker_set_id("ce96", 0.05)$hash
  h_ce100  <- generate_trait_id(ms_hash,    5, "gamma", 1, 0.8)$hash
  h_ce96   <- generate_trait_id(other_ms,   5, "gamma", 1, 0.8)$hash
  expect_false(h_ce100 == h_ce96)
})

test_that("generate_trait_id accepts arbitrary 20-char hex verbatim", {
  arbitrary <- "abcdef1234abcdef1234"
  result <- generate_trait_id(arbitrary, 5, "gamma", 1, 0.8)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_trait_id normalizes effect to lowercase", {
  h_upper <- generate_trait_id(ms_hash, 5, "GAMMA", 1, 0.8)$hash
  h_lower <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash
  expect_equal(h_upper, h_lower)
})

test_that("generate_trait_id rejects non-hex or wrong-length marker_set_hash", {
  expect_error(generate_trait_id("not-a-hash", 5, "gamma", 1, 0.8))
  expect_error(generate_trait_id("879ca1a75fc7fb5183e7EXTRA", 5, "gamma", 1, 0.8))
})

test_that("generate_trait_id golden value", {
  expect_equal(generate_trait_id(ms_hash, 5, "gamma", 1, 0.8)$hash, trait_hash)
})

# ==============================================================================
# generate_mapping_id()
# ==============================================================================

test_that("generate_mapping_id returns list with hash and hash_string", {
  result <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)
  expect_type(result, "list")
  expect_true("hash" %in% names(result))
  expect_true("hash_string" %in% names(result))
})

test_that("generate_mapping_id hash is 20-char lowercase hex", {
  result <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_mapping_id hash_string contains parent trait_hash verbatim", {
  result <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)
  expect_true(grepl(paste0("parent=", trait_hash), result$hash_string))
})

test_that("generate_mapping_id hash_string does not contain simulation params", {
  result <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)
  expect_false(grepl("nqtl=",   result$hash_string))
  expect_false(grepl("rep=",    result$hash_string))
  expect_false(grepl("h2=",     result$hash_string))
  expect_false(grepl("effect=", result$hash_string))
})

test_that("generate_mapping_id is deterministic", {
  h1 <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash
  h2 <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash
  expect_equal(h1, h2)
})

test_that("generate_mapping_id is sensitive to algorithm", {
  h_inbred <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash
  h_loco   <- generate_mapping_id(trait_hash, "LMM-EXACT-LOCO",   TRUE)$hash
  expect_false(h_inbred == h_loco)
})

test_that("generate_mapping_id is sensitive to pca flag", {
  h_pca    <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash
  h_nopca  <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", FALSE)$hash
  expect_false(h_pca == h_nopca)
})

test_that("generate_mapping_id is sensitive to parent trait_hash", {
  other_trait <- generate_trait_id(ms_hash, 3, "gamma", 1, 0.8)$hash
  h1 <- generate_mapping_id(trait_hash,   "LMM-EXACT-INBRED", TRUE)$hash
  h2 <- generate_mapping_id(other_trait,  "LMM-EXACT-INBRED", TRUE)$hash
  expect_false(h1 == h2)
})

test_that("generate_mapping_id accepts arbitrary 20-char hex verbatim", {
  arbitrary <- "abcdef1234abcdef1234"
  result <- generate_mapping_id(arbitrary, "LMM-EXACT-INBRED", TRUE)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_mapping_id rejects string pca", {
  expect_error(generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", pca = "TRUE"))
})

test_that("generate_mapping_id rejects integer pca", {
  expect_error(generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", pca = 1L))
})

test_that("generate_mapping_id pca=TRUE serializes as 'pca=TRUE'", {
  result <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)
  expect_true(grepl("pca=TRUE", result$hash_string))
})

test_that("generate_mapping_id pca=FALSE serializes as 'pca=FALSE'", {
  result <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", FALSE)
  expect_true(grepl("pca=FALSE", result$hash_string))
})

test_that("generate_mapping_id normalizes algorithm to uppercase", {
  h_lower <- generate_mapping_id(trait_hash, "lmm-exact-inbred", TRUE)$hash
  h_upper <- generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash
  expect_equal(h_lower, h_upper)
})

test_that("generate_mapping_id rejects non-hex or wrong-length trait_hash", {
  expect_error(generate_mapping_id("not-a-hash", "LMM-EXACT-INBRED", TRUE))
})

test_that("generate_mapping_id golden value", {
  expect_equal(generate_mapping_id(trait_hash, "LMM-EXACT-INBRED", TRUE)$hash, map_hash)
})

# ==============================================================================
# Path stub (Step 2 scope — skip if not yet implemented)
# ==============================================================================

test_that("get_markers_path returns path ending in {20-char-hex}_markers.parquet", {
  skip("get_markers_path not yet updated for hash-based paths (Step 2)")
  result <- get_markers_path("ce100", 0.05, "data/db")
  expect_match(result, "[0-9a-f]{20}_markers\\.parquet$")
})
