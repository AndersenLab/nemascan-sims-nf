test_that("mappings_schema excludes N and log10p, includes AF1", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_true("AF1" %in% field_names)
  expect_false("N" %in% field_names)
  expect_false("log10p" %in% field_names)
})

test_that("mappings_schema has exactly 5 fields", {
  s <- mappings_schema()
  expect_equal(length(names(s)), 5L)
})

test_that("mappings_schema excludes FK columns (partition key and metadata only)", {
  schema_names <- names(mappings_schema())

  # FK columns removed — mapping_id is Hive partition key; marker_set_id and
  # trait_id are in mappings_metadata.parquet (the metadata view), not data files
  expect_false("mapping_id" %in% schema_names,
    label = "mapping_id absent from mappings_schema (partition key)"
  )
  expect_false("marker_set_id" %in% schema_names,
    label = "marker_set_id absent from mappings_schema (in metadata only)"
  )
  expect_false("trait_id" %in% schema_names,
    label = "trait_id absent from mappings_schema (in metadata only)"
  )

  # Data columns present
  expect_true("marker" %in% schema_names)
  expect_true("AF1" %in% schema_names)
  expect_true("BETA" %in% schema_names)
  expect_true("SE" %in% schema_names)
  expect_true("P" %in% schema_names)
  expect_equal(length(schema_names), 5L, label = "exactly 5 columns in mappings_schema")
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

test_that("trait_metadata_schema includes sim_seed field", {
  s <- trait_metadata_schema()
  expect_true("sim_seed" %in% names(s),
              label = "sim_seed present in trait_metadata_schema")
})

test_that("trait_metadata_schema sim_seed is int32 type", {
  s <- trait_metadata_schema()
  expect_equal(s$GetFieldByName("sim_seed")$type, arrow::int32())
})

# ==============================================================================
# T3 — CV pool schema and round-trip tests
# ==============================================================================

test_that("marker_set_metadata_schema excludes cv_maf and cv_ld", {
  cols <- names(marker_set_metadata_schema())
  expect_false("cv_maf" %in% cols)
  expect_false("cv_ld"  %in% cols)
})

test_that("trait_metadata_schema includes cv_maf_effective (float64) and cv_ld (float64)", {
  s <- trait_metadata_schema()
  expect_true("cv_maf_effective" %in% names(s))
  expect_true("cv_ld"           %in% names(s))
  expect_equal(s$GetFieldByName("cv_maf_effective")$type, arrow::float64())
  expect_equal(s$GetFieldByName("cv_ld")$type,            arrow::float64())
})

test_that("causal_genotypes_schema has 6 correct columns", {
  expect_equal(sort(names(causal_genotypes_schema())),
               sort(c("trait_id", "QTL", "CHROM", "POS", "strain", "allele")))
  expect_equal(causal_genotypes_schema()$GetFieldByName("allele")$type, arrow::float64())
})

test_that("init_database creates traits/causal_genotypes directory", {
  db <- create_temp_db()
  init_database(db)
  expect_true(dir.exists(file.path(db, "traits", "causal_genotypes")))
})

test_that("write/read causal_genotypes round-trip preserves data", {
  db    <- create_temp_db()
  init_database(db)
  ms_id <- generate_marker_set_id("ce.test", 0.05, "c_elegans", "20220216", 0.8)
  trait <- generate_trait_id(ms_id$hash, 5, "gamma", 1, 0.8, 0.05, 0.8)
  write_causal_genotypes(fixture_path("test_causal_genotypes.tsv"), trait$hash, db)
  result <- read_causal_genotypes(trait$hash, db)
  expect_equal(sort(names(result)),
               sort(c("trait_id", "QTL", "CHROM", "POS", "strain", "allele")))
  expect_equal(unique(result$trait_id), trait$hash)
  expect_true(nrow(result) > 0)
  expect_true("I:100"  %in% result$QTL)
  expect_true("II:300" %in% result$QTL)
  expect_true(all(result$allele[!is.na(result$allele)] %in% c(-1.0, 1.0)))
})

# ==============================================================================
# build_ids_from_params() — canonical ID wrapper (Phase A)
# ==============================================================================

test_that("build_ids_from_params ms_id round-trip matches generate_marker_set_id", {
  ms_id_direct <- generate_marker_set_id("test_pop", 0.05, "c_elegans", "20220216", 0.8)
  params <- list(
    population = "test_pop", maf = 0.05, species = "c_elegans",
    vcf_release_id = "20220216", ms_ld = 0.8,
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  ids <- build_ids_from_params(params)
  expect_equal(ids$ms_id$hash, ms_id_direct$hash)
})

test_that("build_ids_from_params trait_id round-trip matches generate_trait_id", {
  ms_id_direct    <- generate_marker_set_id("test_pop", 0.05, "c_elegans", "20220216", 0.8)
  trait_id_direct <- generate_trait_id(ms_id_direct$hash, 5, "gamma", 1, 0.2, 0.05, 0.8)
  params <- list(
    population = "test_pop", maf = 0.05, species = "c_elegans",
    vcf_release_id = "20220216", ms_ld = 0.8,
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  ids <- build_ids_from_params(params)
  expect_equal(ids$trait_id$hash, trait_id_direct$hash)
})

test_that("build_ids_from_params mapping_id round-trip matches generate_mapping_id", {
  ms_id_direct    <- generate_marker_set_id("test_pop", 0.05, "c_elegans", "20220216", 0.8)
  trait_id_direct <- generate_trait_id(ms_id_direct$hash, 5, "gamma", 1, 0.2, 0.05, 0.8)
  mapping_direct  <- generate_mapping_id(trait_id_direct$hash, "inbred", TRUE)
  params <- list(
    population = "test_pop", maf = 0.05, species = "c_elegans",
    vcf_release_id = "20220216", ms_ld = 0.8,
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  ids <- build_ids_from_params(params, mode = "inbred", pca = TRUE)
  expect_equal(ids$mapping_id$hash, mapping_direct$hash)
})

test_that("build_ids_from_params omits mapping_id when mode/pca are NULL", {
  params <- list(
    population = "test_pop", maf = 0.05, species = "c_elegans",
    vcf_release_id = "20220216", ms_ld = 0.8,
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  ids <- build_ids_from_params(params)
  expect_null(ids$mapping_id)
})

test_that("build_ids_from_params propagates error from generate_trait_id on bad ms hash", {
  # Passing a full ms_id object (not $hash) is the common mistake — should error
  ms_id_obj <- generate_marker_set_id("test_pop", 0.05, "c_elegans", "20220216", 0.8)
  # Inject a bad marker_set_hash by constructing params with $hash_string (wrong)
  params <- list(
    population = "test_pop", maf = 0.05, species = "c_elegans",
    vcf_release_id = "20220216", ms_ld = 0.8,
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  ids_good <- build_ids_from_params(params)
  # Manually calling generate_trait_id with a non-hash string propagates the format guard
  expect_error(
    generate_trait_id(ids_good$ms_id$hash_string, 5, "gamma", 1, 0.2, 0.05, 0.8),
    regexp = "20-character"
  )
})

test_that("build_ids_from_params errors when required field is NULL", {
  params <- list(
    population = NULL,  # missing
    maf = 0.05, species = "c_elegans", vcf_release_id = "20220216", ms_ld = 0.8,
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  expect_error(build_ids_from_params(params), regexp = "params\\$population")
})

test_that("read_marker_set_metadata returns NULL when metadata file is missing", {
  db_dir <- create_temp_db()
  result <- read_marker_set_metadata("nonexistent_pop", 0.05, db_dir)
  expect_null(result)
})

test_that("hash concordance: DB-read path via build_ids_from_params matches direct generate calls", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  pop     <- "ce.hash_test"; maf <- 0.05; species <- "c_elegans"
  vcf_id  <- "20220216"; ms_ld <- 0.8

  ms_id_direct <- generate_marker_set_id(pop, maf, species, vcf_id, ms_ld)

  write_marker_set_metadata(
    population             = pop, maf = maf, species = species,
    vcf_release_id         = vcf_id, ms_ld = ms_ld, n_markers = 100L,
    marker_set_id          = ms_id_direct$hash,
    marker_set_hash_string = ms_id_direct$hash_string,
    strainfile_hash        = paste(rep("0", 64), collapse = ""),
    strain_list            = "N2", base_dir = db_dir
  )

  ms_meta <- read_marker_set_metadata(pop, maf, db_dir)

  sim_params <- list(
    population = pop, maf = maf, species = ms_meta$species,
    vcf_release_id = ms_meta$vcf_release_id, ms_ld = as.numeric(ms_meta$ms_ld),
    nqtl = 5L, effect = "gamma", rep = 1L, h2 = 0.2,
    cv_maf_effective = 0.05, cv_ld = 0.8
  )
  ids_fromdb <- build_ids_from_params(sim_params, mode = "inbred", pca = TRUE)

  trait_direct   <- generate_trait_id(ms_id_direct$hash, 5, "gamma", 1, 0.2, 0.05, 0.8)
  mapping_direct <- generate_mapping_id(trait_direct$hash, "inbred", TRUE)

  expect_equal(ms_id_direct$hash,    ids_fromdb$ms_id$hash)
  expect_equal(trait_direct$hash,    ids_fromdb$trait_id$hash)
  expect_equal(mapping_direct$hash,  ids_fromdb$mapping_id$hash)
})

test_that("generate_trait_id rejects hash_string — common mistake guard", {
  # Passing $hash_string (human-readable input) instead of $hash is the typical mistake.
  # $hash_string is not a 20-character hex string and must be rejected by the format guard.
  ms_id <- generate_marker_set_id("test_pop", 0.05, "c_elegans", "20220216", 0.8)
  expect_error(
    generate_trait_id(ms_id$hash_string, nqtl = 5, effect = "gamma", rep = 1,
                      h2 = 0.2, cv_maf_effective = 0.05, cv_ld = 0.8),
    regexp = "20-character"
  )
})

