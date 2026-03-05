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

test_that("marker_set_metadata_schema contains strainfile_hash and strain_list", {
  s <- marker_set_metadata_schema()
  field_names <- names(s)
  expect_true("strainfile_hash" %in% field_names)
  expect_true("strain_list" %in% field_names)
})

test_that("marker_set_metadata_schema strainfile_hash is utf8 type", {
  s <- marker_set_metadata_schema()
  expect_equal(s[["strainfile_hash"]]$type, arrow::utf8())
})

test_that("marker_set_metadata_schema strain_list is utf8 type", {
  s <- marker_set_metadata_schema()
  expect_equal(s[["strain_list"]]$type, arrow::utf8())
})

test_that("write_marker_set_metadata stores strainfile_hash and strain_list", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  ms_id <- generate_marker_set_id("test_pop", 0.05, "c_elegans", "20220216", 0.8)
  write_marker_set_metadata(
    population             = "test_pop",
    maf                    = 0.05,
    species                = "c_elegans",
    vcf_release_id         = "20220216",
    ms_ld                  = 0.8,
    n_markers              = 100L,
    marker_set_id          = ms_id$hash,
    marker_set_hash_string = ms_id$hash_string,
    strainfile_hash        = "a3b4c5d6e7f8a3b4c5d6e7f8a3b4c5d6e7f8a3b4c5d6e7f8a3b4c5d6e7f8a3b4",
    strain_list            = "N2,CB4856,MY23",
    base_dir               = db_dir
  )

  meta <- read_marker_set_metadata("test_pop", 0.05, db_dir)
  expect_equal(meta$strainfile_hash,
               "a3b4c5d6e7f8a3b4c5d6e7f8a3b4c5d6e7f8a3b4c5d6e7f8a3b4c5d6e7f8a3b4")
  expect_equal(meta$strain_list, "N2,CB4856,MY23")
})

test_that("digest::digest with file=TRUE produces 64-char SHA-256 hex", {
  # Validates the hash mechanism used by write_marker_set.R — previously computed in Groovy
  tmp <- tempfile()
  writeLines(c("group\tspecies\tvcf\tms_maf\tms_ld\tstrains",
               "pop1\tc_elegans\t20220216\t0.05\t0.8\tN2,CB4856"), tmp)
  on.exit(unlink(tmp))

  hash <- digest::digest(tmp, algo = "sha256", file = TRUE)
  expect_true(grepl("^[0-9a-f]{64}$", hash),
              label = "SHA-256 via digest::digest(file=TRUE) should be 64-char hex")
  # Deterministic: same content always produces the same hash
  expect_equal(hash, digest::digest(tmp, algo = "sha256", file = TRUE))
})

test_that("write_marker_set_metadata requires strainfile_hash and strain_list", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  ms_id <- generate_marker_set_id("test_pop2", 0.05, "c_elegans", "20220216", 0.8)
  # strainfile_hash and strain_list are required — omitting either must error
  expect_error(
    write_marker_set_metadata(
      population             = "test_pop2",
      maf                    = 0.05,
      species                = "c_elegans",
      vcf_release_id         = "20220216",
      ms_ld                  = 0.8,
      n_markers              = 50L,
      marker_set_id          = ms_id$hash,
      marker_set_hash_string = ms_id$hash_string,
      base_dir               = db_dir
    )
  )
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

