# Golden values recomputed 2026-03-20 after generate_trait_id() → v=2 (adds cv_maf_effective, cv_ld).
# ms_hash   <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash → "db39e634ef0bfa76e1bd"
# trait_hash <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash            → "5a068bb1a67b5e6ba41d"
# map_hash   <- generate_mapping_id(trait_hash, "inbred", TRUE)$hash                       → "0b8f4ba93b6ac58d2fce"

ms_hash    <- "db39e634ef0bfa76e1bd"
trait_hash <- "5a068bb1a67b5e6ba41d"
map_hash   <- "0b8f4ba93b6ac58d2fce"

# ==============================================================================
# generate_marker_set_id()
# ==============================================================================

test_that("generate_marker_set_id returns list with hash and hash_string", {
  result <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)
  expect_type(result, "list")
  expect_true("hash" %in% names(result))
  expect_true("hash_string" %in% names(result))
})

test_that("generate_marker_set_id hash is 20-char lowercase hex", {
  result <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_marker_set_id hash_string starts with v=2| and contains required fields", {
  result <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)
  expect_true(startsWith(result$hash_string, "v=2|"))
  expect_true(grepl("population=", result$hash_string))
  expect_true(grepl("maf=", result$hash_string))
  expect_true(grepl("species=", result$hash_string))
  expect_true(grepl("vcf_release_id=", result$hash_string))
  expect_true(grepl("ms_ld=", result$hash_string))
})

test_that("generate_marker_set_id float serialization uses 10 decimal places", {
  result <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)
  expect_true(grepl("maf=0.0500000000", result$hash_string))
  expect_true(grepl("ms_ld=0.8000000000", result$hash_string))
})

test_that("generate_marker_set_id is deterministic", {
  h1 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  h2 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  expect_equal(h1, h2)
})

test_that("generate_marker_set_id is sensitive to population", {
  h_ce100 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  h_ce96  <- generate_marker_set_id("ce.96",   0.05, "c_elegans", "20220216", 0.8)$hash
  expect_false(h_ce100 == h_ce96)
})

test_that("generate_marker_set_id is sensitive to maf", {
  h_05 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  h_01 <- generate_marker_set_id("ce.test", 0.01, "c_elegans", "20220216", 0.8)$hash
  expect_false(h_05 == h_01)
})

test_that("generate_marker_set_id is sensitive to species", {
  h_ce <- generate_marker_set_id("ce.test", 0.05, "c_elegans",  "20220216", 0.8)$hash
  h_cb <- generate_marker_set_id("ce.test", 0.05, "c_briggsae", "20220216", 0.8)$hash
  expect_false(h_ce == h_cb)
})

test_that("generate_marker_set_id is sensitive to vcf_release_id", {
  h1 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  h2 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20210121", 0.8)$hash
  expect_false(h1 == h2)
})

test_that("generate_marker_set_id is sensitive to ms_ld", {
  h_08 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  h_05 <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.5)$hash
  expect_false(h_08 == h_05)
})

test_that("generate_marker_set_id normalizes population to lowercase", {
  h_upper <- generate_marker_set_id("CE.TEST", 0.05, "c_elegans", "20220216", 0.8)$hash
  h_lower <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  expect_equal(h_upper, h_lower)
})

test_that("generate_marker_set_id normalizes species to lowercase", {
  h_upper <- generate_marker_set_id("ce.test", 0.05, "C_ELEGANS", "20220216", 0.8)$hash
  h_lower <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)$hash
  expect_equal(h_upper, h_lower)
})

test_that("generate_marker_set_id trims whitespace from population", {
  h_spaces <- generate_marker_set_id("  ce.test  ", 0.05, "c_elegans", "20220216", 0.8)$hash
  h_clean  <- generate_marker_set_id("ce.test",     0.05, "c_elegans", "20220216", 0.8)$hash
  expect_equal(h_spaces, h_clean)
})

test_that("generate_marker_set_id canonical form stores lowercase population", {
  result <- generate_marker_set_id("CE.TEST", 0.05, "c_elegans", "20220216", 0.8)
  expect_true(grepl("population=ce.test", result$hash_string))
})

test_that("generate_marker_set_id canonical form stores lowercase species", {
  result <- generate_marker_set_id("ce.test", 0.05, "C_ELEGANS", "20220216", 0.8)
  expect_true(grepl("species=c_elegans", result$hash_string))
})

test_that("generate_marker_set_id v=2 golden value", {
  result <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)
  expect_equal(result$hash, ms_hash)
  expect_equal(result$hash_string,
    "v=2|population=ce.test|maf=0.0500000000|species=c_elegans|vcf_release_id=20220216|ms_ld=0.8000000000")
})

# ==============================================================================
# generate_trait_id()
# ==============================================================================

test_that("generate_trait_id returns list with hash and hash_string", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)
  expect_type(result, "list")
  expect_true("hash" %in% names(result))
  expect_true("hash_string" %in% names(result))
})

test_that("generate_trait_id hash is 20-char lowercase hex", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_trait_id hash_string contains parent hash verbatim", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)
  expect_true(grepl(paste0("parent=", ms_hash), result$hash_string))
})

test_that("generate_trait_id float serialization uses 10 decimal places for h2", {
  result <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)
  expect_true(grepl("h2=0.8000000000", result$hash_string))
})

test_that("generate_trait_id is deterministic", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  expect_equal(h1, h2)
})

test_that("generate_trait_id is sensitive to nqtl", {
  h5 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h3 <- generate_trait_id(ms_hash, 3, "gamma", 1, 0.8, 0.05, 0.8)$hash
  expect_false(h5 == h3)
})

test_that("generate_trait_id is sensitive to effect", {
  h_gamma  <- generate_trait_id(ms_hash, 5, "gamma",  1, 0.8, 0.05, 0.8)$hash
  h_normal <- generate_trait_id(ms_hash, 5, "normal", 1, 0.8, 0.05, 0.8)$hash
  expect_false(h_gamma == h_normal)
})

test_that("generate_trait_id accepts numeric range effect string (e.g. '0.2-0.3')", {
  result <- generate_trait_id(ms_hash, 5, "0.2-0.3", 1, 0.8, 0.05, 0.8)
  expect_match(result$hash, "^[0-9a-f]{20}$")
  expect_true(grepl("effect=0.2-0.3", result$hash_string))
})

test_that("generate_trait_id range effect differs from named distribution", {
  h_range <- generate_trait_id(ms_hash, 5, "0.2-0.3", 1, 0.8, 0.05, 0.8)$hash
  h_gamma <- generate_trait_id(ms_hash, 5, "gamma",   1, 0.8, 0.05, 0.8)$hash
  expect_false(h_range == h_gamma)
})

test_that("generate_trait_id is sensitive to rep", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 2, 0.8, 0.05, 0.8)$hash
  expect_false(h1 == h2)
})

test_that("generate_trait_id is sensitive to h2", {
  h_08 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h_04 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.4, 0.05, 0.8)$hash
  expect_false(h_08 == h_04)
})

test_that("generate_trait_id is sensitive to parent marker_set_hash", {
  other_ms <- generate_marker_set_id("ce.96", 0.05, "c_elegans", "20220216", 0.8)$hash
  h_ce100  <- generate_trait_id(ms_hash,  5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h_ce96   <- generate_trait_id(other_ms, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  expect_false(h_ce100 == h_ce96)
})

test_that("generate_trait_id accepts arbitrary 20-char hex verbatim", {
  arbitrary <- "abcdef1234abcdef1234"
  result <- generate_trait_id(arbitrary, 5, "gamma", 1, 0.8, 0.05, 0.8)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_trait_id normalizes effect to lowercase", {
  h_upper <- generate_trait_id(ms_hash, 5, "GAMMA", 1, 0.8, 0.05, 0.8)$hash
  h_lower <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  expect_equal(h_upper, h_lower)
})

test_that("generate_trait_id rejects non-hex or wrong-length marker_set_hash", {
  expect_error(generate_trait_id("not-a-hash", 5, "gamma", 1, 0.8, 0.05, 0.8))
  expect_error(generate_trait_id("879ca1a75fc7fb5183e7EXTRA", 5, "gamma", 1, 0.8, 0.05, 0.8))
})

test_that("generate_trait_id golden value", {
  expect_equal(generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash, trait_hash)
})

# ==============================================================================
# generate_mapping_id()
# ==============================================================================

test_that("generate_mapping_id returns list with hash and hash_string", {
  result <- generate_mapping_id(trait_hash, "inbred", TRUE)
  expect_type(result, "list")
  expect_true("hash" %in% names(result))
  expect_true("hash_string" %in% names(result))
})

test_that("generate_mapping_id hash is 20-char lowercase hex", {
  result <- generate_mapping_id(trait_hash, "inbred", TRUE)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_mapping_id hash_string contains parent trait_hash verbatim", {
  result <- generate_mapping_id(trait_hash, "inbred", TRUE)
  expect_true(grepl(paste0("parent=", trait_hash), result$hash_string))
})

test_that("generate_mapping_id hash_string does not contain simulation params", {
  result <- generate_mapping_id(trait_hash, "inbred", TRUE)
  expect_false(grepl("nqtl=",   result$hash_string))
  expect_false(grepl("rep=",    result$hash_string))
  expect_false(grepl("h2=",     result$hash_string))
  expect_false(grepl("effect=", result$hash_string))
})

test_that("generate_mapping_id is deterministic", {
  h1 <- generate_mapping_id(trait_hash, "inbred", TRUE)$hash
  h2 <- generate_mapping_id(trait_hash, "inbred", TRUE)$hash
  expect_equal(h1, h2)
})

test_that("generate_mapping_id is sensitive to algorithm", {
  h_inbred <- generate_mapping_id(trait_hash, "inbred", TRUE)$hash
  h_loco   <- generate_mapping_id(trait_hash, "loco",   TRUE)$hash
  expect_false(h_inbred == h_loco)
})

test_that("generate_mapping_id is sensitive to pca flag", {
  h_pca    <- generate_mapping_id(trait_hash, "inbred", TRUE)$hash
  h_nopca  <- generate_mapping_id(trait_hash, "inbred", FALSE)$hash
  expect_false(h_pca == h_nopca)
})

test_that("generate_mapping_id is sensitive to parent trait_hash", {
  other_trait <- generate_trait_id(ms_hash, 3, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h1 <- generate_mapping_id(trait_hash,   "inbred", TRUE)$hash
  h2 <- generate_mapping_id(other_trait,  "inbred", TRUE)$hash
  expect_false(h1 == h2)
})

test_that("generate_mapping_id accepts arbitrary 20-char hex verbatim", {
  arbitrary <- "abcdef1234abcdef1234"
  result <- generate_mapping_id(arbitrary, "inbred", TRUE)
  expect_match(result$hash, "^[0-9a-f]{20}$")
})

test_that("generate_mapping_id rejects string pca", {
  expect_error(generate_mapping_id(trait_hash, "inbred", pca = "TRUE"))
})

test_that("generate_mapping_id rejects integer pca", {
  expect_error(generate_mapping_id(trait_hash, "inbred", pca = 1L))
})

test_that("generate_mapping_id pca=TRUE serializes as 'pca=TRUE'", {
  result <- generate_mapping_id(trait_hash, "inbred", TRUE)
  expect_true(grepl("pca=TRUE", result$hash_string))
})

test_that("generate_mapping_id pca=FALSE serializes as 'pca=FALSE'", {
  result <- generate_mapping_id(trait_hash, "inbred", FALSE)
  expect_true(grepl("pca=FALSE", result$hash_string))
})

test_that("generate_mapping_id normalizes algorithm to lowercase", {
  h_upper <- generate_mapping_id(trait_hash, "INBRED", TRUE)$hash
  h_lower <- generate_mapping_id(trait_hash, "inbred", TRUE)$hash
  expect_equal(h_upper, h_lower)
})

test_that("generate_mapping_id rejects non-hex or wrong-length trait_hash", {
  expect_error(generate_mapping_id("not-a-hash", "inbred", TRUE))
})

test_that("generate_mapping_id golden value", {
  expect_equal(generate_mapping_id(trait_hash, "inbred", TRUE)$hash, map_hash)
})

# ==============================================================================
# Path functions (Step 2 scope)
# ==============================================================================

test_that("get_markers_path returns path ending in {20-char-hex}_markers.parquet", {
  result <- get_markers_path("ce.test", 0.05, "c_elegans", "20220216", 0.8, "data/db")
  expect_match(result, "[0-9a-f]{20}_markers\\.parquet$")
})

test_that("get_markers_path path includes markers/marker_sets subdir", {
  result <- get_markers_path("ce.test", 0.05, "c_elegans", "20220216", 0.8, "data/db")
  expect_true(grepl("markers/marker_sets", result, fixed = TRUE))
})

test_that("get_markers_path is deterministic", {
  p1 <- get_markers_path("ce.test", 0.05, "c_elegans", "20220216", 0.8, "data/db")
  p2 <- get_markers_path("ce.test", 0.05, "c_elegans", "20220216", 0.8, "data/db")
  expect_equal(p1, p2)
})

test_that("get_markers_path produces different paths for different params", {
  p_ce100 <- get_markers_path("ce.test", 0.05, "c_elegans", "20220216", 0.8, "data/db")
  p_ce96  <- get_markers_path("ce.96",   0.05, "c_elegans", "20220216", 0.8, "data/db")
  expect_false(p_ce100 == p_ce96)
})

test_that("init_database creates markers/marker_sets, markers/genotypes, traits/causal_variants, traits/phenotypes subdirs", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  expect_true(dir.exists(file.path(db_dir, "markers", "marker_sets")))
  expect_true(dir.exists(file.path(db_dir, "markers", "genotypes")))
  expect_true(dir.exists(file.path(db_dir, "traits", "causal_variants")))
  expect_true(dir.exists(file.path(db_dir, "traits", "phenotypes")))
  expect_true(dir.exists(file.path(db_dir, "mappings")))
})

# ==============================================================================
# T2 — CV pool hash tests
# ==============================================================================

# CV exclusion from generate_marker_set_id()
test_that("generate_marker_set_id rejects cv_maf argument", {
  expect_error(
    generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8, cv_maf = 0.01),
    regexp = "unused argument"
  )
})

# v=2 trait hash structure
test_that("generate_trait_id v=2 hash_string starts with v=2| and encodes cv params", {
  r <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)
  expect_true(startsWith(r$hash_string, "v=2|"))
  expect_true(grepl("cv_maf_effective=0.0500000000", r$hash_string))
  expect_true(grepl("cv_ld=0.8000000000", r$hash_string))
})

# v=2 sensitivity
test_that("generate_trait_id v=2 is sensitive to cv_maf_effective", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.01, 0.8)$hash
  expect_false(h1 == h2)
})

test_that("generate_trait_id v=2 is sensitive to cv_ld", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.99)$hash
  expect_false(h1 == h2)
})

test_that("generate_trait_id v=2 is deterministic", {
  h1 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  h2 <- generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash
  expect_equal(h1, h2)
})

# Golden value (uses trait_hash variable set at top of file)
test_that("generate_trait_id v=2 golden value", {
  expect_equal(
    generate_trait_id(ms_hash, 5, "gamma", 1, 0.8, 0.05, 0.8)$hash,
    trait_hash
  )
})
