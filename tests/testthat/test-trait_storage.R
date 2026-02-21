# test-trait_storage.R - Unit tests for trait data storage functions
#
# Tests generate_trait_id(), write/read_trait_metadata(), write/read_causal_variants(),
# write/read_phenotype_data(), phenotype_exists().
#
# Also validates Step 1 (glob narrowing), Step 4 (var.exp removal and overwrite default).
#
# Fixtures:
#   tests/fixtures/test_sims.par  - real causal variants (.par format)
#   tests/fixtures/test_sims.phen - real phenotype data (space-delimited, 3 strains)


# ── generate_trait_id() ───────────────────────────────────────────────────────

test_that("generate_trait_id is deterministic across calls", {
  id1 <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  id2 <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  expect_equal(id1, id2)
})

test_that("generate_trait_id produces different results for different params", {
  id_a <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  id_b <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 2, 0.8)
  id_c <- generate_trait_id("ce.test.200strains", 0.05, 10, "gamma", 1, 0.8)
  expect_false(id_a == id_b, label = "different rep -> different id")
  expect_false(id_a == id_c, label = "different nqtl -> different id")
  expect_false(id_b == id_c)
})

test_that("generate_trait_id output is 12-character lowercase hex", {
  id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  expect_equal(nchar(id), 12L)
  expect_match(id, "^[0-9a-f]{12}$")
})

test_that("generate_trait_id produces expected golden value", {
  # Golden-value test: pins generate_trait_id() output against a known reference.
  # If this test fails, the hash algorithm or input formatting has changed —
  # all stored trait data would become unreachable.
  #
  # Computed from: digest::digest("ce.test.200strains_0.05_5_gamma_1_0.8",
  #                               algo = "md5", serialize = FALSE)
  result <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  expect_equal(result, "41743948dbc0")
})


# ── Trait metadata round-trip ─────────────────────────────────────────────────

test_that("write_trait_metadata / read_trait_metadata round-trip preserves all fields", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  trait_id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)

  write_trait_metadata(
    trait_id  = trait_id,
    nqtl      = 5L,
    rep       = 1L,
    h2        = 0.8,
    maf       = 0.05,
    effect    = "gamma",
    population = "ce.test.200strains",
    base_dir  = db_dir
  )

  result <- read_trait_metadata(trait_id, db_dir)

  expect_equal(result$trait_id, trait_id)
  expect_equal(result$nqtl, 5L)
  expect_equal(result$rep, 1L)
  expect_equal(result$h2, 0.8)
  expect_equal(result$maf, 0.05)
  expect_equal(result$effect, "gamma")
  expect_equal(result$population, "ce.test.200strains")
})

test_that("write_trait_metadata created_at is ISO 8601 format string", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  trait_id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)

  write_trait_metadata(
    trait_id   = trait_id,
    nqtl       = 5L,
    rep        = 1L,
    h2         = 0.8,
    maf        = 0.05,
    effect     = "gamma",
    population = "ce.test.200strains",
    base_dir   = db_dir
  )

  result <- read_trait_metadata(trait_id, db_dir)
  expect_match(result$created_at, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}$")
})


# ── Causal variants round-trip ────────────────────────────────────────────────

test_that("write_causal_variants / read_causal_variants_data round-trip", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  trait_id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  par_file <- fixture_path("test_sims.par")

  write_causal_variants(par_file, trait_id, db_dir)

  result <- read_causal_variants_data(trait_id, db_dir)

  # Should have columns from load_causal_variants(): QTL, CHROM, POS, RefAllele, Frequency, Effect
  expect_true("QTL" %in% names(result))
  expect_true("CHROM" %in% names(result))
  expect_true("POS" %in% names(result))
  expect_gt(nrow(result), 0L)

  # Round-trip: values preserved
  original <- load_causal_variants(par_file)
  expect_equal(nrow(result), nrow(original))
  expect_equal(result$QTL, original$QTL)
})


# ── Phenotype data round-trip ─────────────────────────────────────────────────

test_that("write_phenotype_data / read_phenotype_data round-trip preserves values", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  trait_id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  phen_file <- fixture_path("test_sims.phen")

  write_phenotype_data(phen_file, trait_id, db_dir)

  result <- read_phenotype_data(trait_id, db_dir)

  expect_equal(nrow(result), 3L)
  expect_setequal(names(result), c("strain", "phenotype"))
  expect_setequal(result$strain, c("AB1", "BRC20067", "BRC20263"))

  # Values from the real fixture
  ab1_val <- result$phenotype[result$strain == "AB1"]
  expect_equal(ab1_val, 1.43652, tolerance = 1e-5)
})

test_that("read_phenotype_data output has correct Arrow types via Parquet schema", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  trait_id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  phen_file <- fixture_path("test_sims.phen")

  write_phenotype_data(phen_file, trait_id, db_dir)

  path <- get_phenotype_path(trait_id, db_dir)
  tbl <- arrow::read_parquet(path, as_data_frame = FALSE)
  schema <- tbl$schema

  expect_equal(schema$GetFieldByName("strain")$type, arrow::utf8())
  expect_equal(schema$GetFieldByName("phenotype")$type, arrow::float64())
})

test_that("phenotype_exists returns FALSE before write, TRUE after", {
  db_dir <- create_temp_db()
  init_database(db_dir)
  trait_id <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)

  expect_false(phenotype_exists(trait_id, db_dir))

  write_phenotype_data(fixture_path("test_sims.phen"), trait_id, db_dir)

  expect_true(phenotype_exists(trait_id, db_dir))
})


# ── Overwrite and re-run behavior ─────────────────────────────────────────────

test_that("deterministic re-run: same params produce same trait_id and clean overwrite", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  trait_id_1 <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  write_phenotype_data(fixture_path("test_sims.phen"), trait_id_1, db_dir)

  # Re-run with same params
  trait_id_2 <- generate_trait_id("ce.test.200strains", 0.05, 5, "gamma", 1, 0.8)
  expect_equal(trait_id_1, trait_id_2)

  expect_no_error(
    write_phenotype_data(fixture_path("test_sims.phen"), trait_id_2, db_dir,
                         overwrite = TRUE)
  )

  result <- read_phenotype_data(trait_id_1, db_dir)
  expect_equal(nrow(result), 3L)
})


# ── Step 1 validation: glob narrowing ────────────────────────────────────────

test_that("open_mapping_db markers view excludes genotype files", {
  db_dir <- create_temp_db()
  init_database(db_dir)

  # Write a genotype matrix (stored in markers/ dir as *_genotypes.parquet)
  write_genotype_matrix(
    fixture_path("test_genotype_matrix.tsv"),
    "test_pop", 0.05, db_dir
  )

  # Confirm the genotype file is in the markers/ directory
  markers_dir <- file.path(db_dir, "markers")
  all_files <- list.files(markers_dir)
  genotype_files <- grep("_genotypes\\.parquet$", all_files, value = TRUE)
  expect_equal(length(genotype_files), 1L,
               label = "genotype parquet file is in markers/ dir")

  # open_mapping_db uses pattern "_markers\.parquet$" — genotype files must not appear
  # in the markers view
  con <- open_mapping_db(db_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  tables <- DBI::dbListTables(con)
  if ("markers" %in% tables) {
    # markers view should only include files matching _markers.parquet
    # (no marker set was written, so the view may not exist — that's fine)
    result <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM markers")
    # The count should be 0 since we only wrote a genotype file (no marker set)
    expect_equal(result$n, 0L,
                 label = "markers view has 0 rows when only genotype file exists")
  } else {
    succeed("markers view not created when no *_markers.parquet files exist — correct")
  }
})


# ── Step 4 validation: var.exp removed from mappings schema ──────────────────

test_that("mappings_schema does not include var.exp column", {
  s <- mappings_schema()
  field_names <- names(s)
  expect_false("var.exp" %in% field_names,
               label = "var.exp absent from mappings_schema()")
})


# ── Step 4 validation: write_marker_set overwrite default is TRUE ─────────────

test_that("write_marker_set default overwrite argument is TRUE", {
  # Inspect formal arguments of write_marker_set
  args <- formals(write_marker_set)
  expect_true("overwrite" %in% names(args),
               label = "write_marker_set has overwrite parameter")
  expect_true(isTRUE(args$overwrite),
               label = "write_marker_set default overwrite is TRUE")
})
