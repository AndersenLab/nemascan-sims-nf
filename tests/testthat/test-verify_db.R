# test-verify_db.R - Unit tests for bin/verify_db.R
#
# Builds a tiny fixture database in a tempdir using the real database.R hash
# helpers, then exercises verify_db() in both expectation and discovery modes.
# These tests run as plain unit tests — they do NOT require TEST_DB_DIR or any
# integration pipeline output.

# Tiny null-coalescing helper (base R only gained %||% in 4.4)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# ── Setup ────────────────────────────────────────────────────────────────────

verify_db_script <- file.path(.project_root, "bin", "verify_db.R")
skip_if_no_verify_db <- function() {
  if (!file.exists(verify_db_script)) {
    skip("bin/verify_db.R not found")
  }
  if (!exists("verify_db", mode = "function")) {
    # source the script once per test run
    Sys.setenv(NEMASCAN_R_LIB = file.path(.project_root, "R"))
    source(verify_db_script, local = FALSE)
  }
}

# Helper: build a minimal but valid DB layout in `base_dir`.
# Writes a real mappings_metadata.parquet (queryable by DuckDB), per-group
# marker_set_metadata, and empty-but-valid parquet placeholder files for the
# marker_sets, genotypes, and traits artifacts. The returned list is the full
# expected architecture so tests can call verify_db() symmetrically.
build_fixture_db <- function(base_dir,
                             populations = c("ce.test", "cb.test"),
                             species     = c("c_elegans", "c_briggsae"),
                             vcf_release = c("20220216", "20210803"),
                             maf         = 0.05,
                             ms_ld       = 0.8,
                             nqtl        = c(1L, 5L),
                             h2          = c(0.3, 0.8),
                             effect      = "gamma",
                             reps        = 2L,
                             cv_maf      = 0.05,
                             cv_ld       = 0.8,
                             algorithm   = c("inbred", "loco"),
                             pca         = c(TRUE, FALSE)) {
  stopifnot(length(populations) == length(species),
            length(populations) == length(vcf_release))

  dir.create(file.path(base_dir, "markers", "marker_sets"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_dir, "markers", "genotypes"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_dir, "marker_set_metadata"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_dir, "traits", "causal_variants"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_dir, "traits", "phenotypes"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_dir, "traits", "causal_genotypes"),
             recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(base_dir, "mappings"),
             recursive = TRUE, showWarnings = FALSE)

  stub_parquet <- function(path) {
    arrow::write_parquet(data.frame(x = integer(0)), path)
  }

  ms_rows <- list()
  metadata_rows <- list()
  groups_info <- list()

  for (i in seq_along(populations)) {
    pop <- populations[i]
    sp  <- species[i]
    rel <- vcf_release[i]

    ms_id <- generate_marker_set_id(pop, maf, sp, rel, ms_ld)
    groups_info[[i]] <- list(
      population = pop, maf = maf, species = sp,
      vcf_release_id = rel, ms_ld = ms_ld, marker_set_id = ms_id$hash
    )

    # marker set metadata per-population parquet
    ms_rows[[length(ms_rows) + 1L]] <- data.frame(
      marker_set_id          = ms_id$hash,
      marker_set_hash_string = ms_id$hash_string,
      hash_schema_version    = "v=2",
      population             = pop,
      maf                    = as.numeric(maf),
      species                = sp,
      vcf_release_id         = rel,
      ms_ld                  = as.numeric(ms_ld),
      n_markers              = 100L,
      n_independent_tests    = 80,
      eigen_source_file      = NA_character_,
      strainfile_hash        = "deadbeef",
      strain_list            = "A,B,C",
      created_at             = Sys.time(),
      stringsAsFactors       = FALSE
    )
    arrow::write_parquet(
      ms_rows[[length(ms_rows)]],
      file.path(base_dir, "marker_set_metadata",
                paste0(pop, "_", maf, "_metadata.parquet"))
    )

    # marker_set + genotype placeholder files
    stub_parquet(file.path(base_dir, "markers", "marker_sets",
                           paste0(ms_id$hash, "_markers.parquet")))
    stub_parquet(file.path(base_dir, "markers", "genotypes",
                           paste0(ms_id$hash, "_genotypes.parquet")))

    for (nq in nqtl) {
      for (hval in h2) {
        for (eff in effect) {
          for (r in seq_len(reps)) {
            params <- list(
              population       = pop,
              maf              = as.numeric(maf),
              species          = sp,
              vcf_release_id   = rel,
              ms_ld            = as.numeric(ms_ld),
              nqtl             = as.integer(nq),
              effect           = eff,
              rep              = as.integer(r),
              h2               = as.numeric(hval),
              cv_maf_effective = as.numeric(cv_maf),
              cv_ld            = as.numeric(cv_ld)
            )
            trait <- generate_trait_id(ms_id$hash, nq, eff, r, hval, cv_maf, cv_ld)

            # trait artifact placeholders
            stub_parquet(file.path(base_dir, "traits",
                                   paste0(trait$hash, ".parquet")))
            stub_parquet(file.path(base_dir, "traits", "causal_variants",
                                   paste0(trait$hash, "_causal.parquet")))
            stub_parquet(file.path(base_dir, "traits", "phenotypes",
                                   paste0(trait$hash, "_phenotype.parquet")))
            stub_parquet(file.path(base_dir, "traits", "causal_genotypes",
                                   paste0(trait$hash, "_causal_geno.parquet")))

            for (alg in algorithm) {
              for (usePCA in pca) {
                mapping <- generate_mapping_id(trait$hash, alg, usePCA)

                part_dir <- file.path(
                  base_dir, "mappings",
                  paste0("population=", pop),
                  paste0("mapping_id=", mapping$hash)
                )
                dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
                arrow::write_parquet(
                  data.frame(
                    marker = character(0), AF1 = numeric(0),
                    BETA = numeric(0), SE = numeric(0), P = numeric(0)
                  ),
                  file.path(part_dir, "data.parquet")
                )
                arrow::write_parquet(
                  data.frame(mapping_id = mapping$hash),
                  file.path(part_dir, "meta.parquet")
                )

                metadata_rows[[length(metadata_rows) + 1L]] <- data.frame(
                  mapping_id          = mapping$hash,
                  mapping_hash_string = mapping$hash_string,
                  trait_id            = trait$hash,
                  marker_set_id       = ms_id$hash,
                  hash_schema_version = "v=2",
                  population          = pop,
                  maf                 = as.numeric(maf),
                  nqtl                = as.integer(nq),
                  rep                 = as.integer(r),
                  h2                  = as.numeric(hval),
                  effect              = eff,
                  algorithm           = alg,
                  pca                 = as.logical(usePCA),
                  n_markers           = 100L,
                  source_file         = NA_character_,
                  processed_at        = Sys.time(),
                  processing_version  = "test",
                  stringsAsFactors    = FALSE
                )
              }
            }
          }
        }
      }
    }
  }

  metadata_df <- dplyr::bind_rows(metadata_rows)
  arrow::write_parquet(
    metadata_df,
    file.path(base_dir, "mappings_metadata.parquet")
  )

  # Minimal strainfile TSV matching the production format:
  #   group\tspecies\tvcf\tms_maf\tms_ld\tstrains
  strainfile_path <- file.path(base_dir, "strains.tsv")
  strainfile_df <- data.frame(
    group   = populations,
    species = species,
    vcf     = paste0("/fake/", species, "/vcf.gz"),
    ms_maf  = as.character(maf),
    ms_ld   = as.character(ms_ld),
    strains = "A,B,C",
    stringsAsFactors = FALSE
  )
  readr::write_tsv(strainfile_df, strainfile_path)

  list(
    base_dir    = base_dir,
    strainfile  = strainfile_path,
    populations = populations,
    nqtl        = nqtl,
    h2          = h2,
    effect      = effect,
    reps        = reps,
    cv_maf      = cv_maf,
    cv_ld       = cv_ld,
    algorithm   = algorithm,
    pca         = pca,
    n_marker_sets = length(populations),
    n_traits    = length(populations) * length(nqtl) * length(h2) *
                  length(effect) * reps,
    n_mappings  = length(populations) * length(nqtl) * length(h2) *
                  length(effect) * reps * length(algorithm) * length(pca),
    metadata_df = metadata_df,
    groups_info = groups_info
  )
}

# ── Tests ────────────────────────────────────────────────────────────────────

test_that("verify_db returns a well-formed list for an empty DB (discovery)", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_empty_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  result <- verify_db(tmp)

  expected_names <- c("status", "mode", "base_dir", "run_at", "expected_params",
                      "observed_mappings", "observed_traits",
                      "mappings_delta", "traits_delta",
                      "missing_traits", "missing_mappings",
                      "summary", "marker_set_coverage", "errors")
  expect_true(all(expected_names %in% names(result)))
  expect_equal(result$mode, "discovery")
  expect_equal(result$status, "discovery")
  expect_null(result$mappings_delta)
  expect_null(result$missing_mappings)
  expect_s3_class(result$errors, "data.frame")
})

test_that("verify_db reports 'complete' on a fully-populated fixture", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_full_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fx <- build_fixture_db(tmp)

  result <- verify_db(
    base_dir   = fx$base_dir,
    strainfile = fx$strainfile,
    nqtl       = fx$nqtl,
    h2         = fx$h2,
    effect     = fx$effect,
    reps       = fx$reps,
    cv_maf     = fx$cv_maf,
    cv_ld      = fx$cv_ld
  )

  expect_equal(result$mode, "expectation")
  expect_equal(result$status, "complete")
  expect_false(is.null(result$mappings_delta))
  expect_equal(sum(result$mappings_delta$missing), 0L)
  expect_equal(nrow(result$missing_mappings %||% data.frame()), 0L)
  expect_equal(nrow(result$missing_traits %||% data.frame()), 0L)

  mappings_row <- result$summary[result$summary$component == "mappings_metadata_rows", ]
  expect_equal(mappings_row$expected, fx$n_mappings)
  expect_equal(mappings_row$observed, fx$n_mappings)
})

test_that("verify_db reports 'incomplete' and lists missing mappings when rows are dropped", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_partial_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fx <- build_fixture_db(tmp)

  # Drop the last 4 rows from mappings_metadata.parquet to simulate partial run
  meta_df <- fx$metadata_df
  dropped <- meta_df[(nrow(meta_df) - 3L):nrow(meta_df), ]
  truncated <- meta_df[seq_len(nrow(meta_df) - 4L), ]
  arrow::write_parquet(
    truncated,
    file.path(fx$base_dir, "mappings_metadata.parquet")
  )

  result <- verify_db(
    base_dir   = fx$base_dir,
    strainfile = fx$strainfile,
    nqtl       = fx$nqtl,
    h2         = fx$h2,
    effect     = fx$effect,
    reps       = fx$reps,
    cv_maf     = fx$cv_maf,
    cv_ld      = fx$cv_ld
  )

  expect_equal(result$status, "incomplete")
  expect_false(is.null(result$missing_mappings))
  expect_equal(nrow(result$missing_mappings), 4L)
  expect_true(all(dropped$mapping_id %in% result$missing_mappings$mapping_id))

  mappings_row <- result$summary[result$summary$component == "mappings_metadata_rows", ]
  expect_equal(mappings_row$missing, 4L)
})

test_that("verify_db discovery mode derives the observed parameter grid", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_disc_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fx <- build_fixture_db(tmp)

  result <- verify_db(base_dir = fx$base_dir)

  expect_equal(result$mode, "discovery")
  expect_equal(result$status, "discovery")
  expect_false(is.null(result$observed_mappings))
  expect_gt(nrow(result$observed_mappings), 0L)
  expect_null(result$mappings_delta)
  expect_null(result$missing_mappings)

  # observed_grid should reflect the fixture
  expect_setequal(result$observed_grid$populations, fx$populations)
  expect_setequal(result$observed_grid$nqtl,        fx$nqtl)
  expect_setequal(result$observed_grid$h2,          fx$h2)
  expect_setequal(result$observed_grid$algorithm,   fx$algorithm)
})

test_that("verify_db never throws on a corrupted mappings_metadata.parquet", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_corrupt_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Valid directory, but the metadata file is zero bytes → parquet parse fail
  writeLines(character(), file.path(tmp, "mappings_metadata.parquet"))

  expect_no_error(result <- verify_db(tmp))
  expect_true(result$status %in% c("error", "discovery", "incomplete"))
  expect_gte(nrow(result$errors), 1L)
})

test_that("verify_db honors max_enumeration_rows guard", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_guard_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fx <- build_fixture_db(tmp)

  result <- verify_db(
    base_dir             = fx$base_dir,
    strainfile           = fx$strainfile,
    nqtl                 = fx$nqtl,
    h2                   = fx$h2,
    effect               = fx$effect,
    reps                 = fx$reps,
    cv_maf               = fx$cv_maf,
    cv_ld                = fx$cv_ld,
    max_enumeration_rows = 1L
  )

  # Aggregate delta is still produced, but missing_mappings is NULL because
  # per-ID enumeration was skipped
  expect_false(is.null(result$mappings_delta))
  expect_null(result$missing_mappings)
  expect_true(any(grepl("max_enumeration_rows", result$errors$message)))
})

test_that("verify_db fails closed on partial architecture args", {
  skip_if_no_verify_db()
  tmp <- tempfile(pattern = "vdb_partial_args_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fx <- build_fixture_db(tmp)

  # Omit `reps` → expectation mode is requested but incomplete
  result <- verify_db(
    base_dir   = fx$base_dir,
    strainfile = fx$strainfile,
    nqtl       = fx$nqtl,
    h2         = fx$h2,
    effect     = fx$effect
  )

  expect_equal(result$mode, "expectation")
  expect_equal(result$status, "incomplete")
  expect_true(any(grepl("missing arguments", result$errors$message)))
})

test_that("plot helpers return ggplot objects and never throw", {
  skip_if_no_verify_db()
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    skip("ggplot2 not available")
  }
  tmp <- tempfile(pattern = "vdb_plot_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  fx <- build_fixture_db(tmp)
  result <- verify_db(
    base_dir   = fx$base_dir,
    strainfile = fx$strainfile,
    nqtl       = fx$nqtl,
    h2         = fx$h2,
    effect     = fx$effect,
    reps       = fx$reps,
    cv_maf     = fx$cv_maf,
    cv_ld      = fx$cv_ld
  )

  expect_s3_class(plot_n_mappings_by_param(result), "ggplot")
  expect_s3_class(plot_n_traits_by_param(result),   "ggplot")
  expect_s3_class(plot_completeness_heatmap(result), "ggplot")
  expect_s3_class(plot_marker_set_coverage(result), "ggplot")

  # Even with a degenerate (empty) result, plot helpers must not throw
  empty <- list(
    mode = "discovery",
    observed_mappings = NULL, observed_traits = NULL,
    mappings_delta = NULL,    marker_set_coverage = NULL,
    expected_params = list(reps = NULL, effect = NULL)
  )
  expect_s3_class(plot_n_mappings_by_param(empty),  "ggplot")
  expect_s3_class(plot_n_traits_by_param(empty),    "ggplot")
  expect_s3_class(plot_completeness_heatmap(empty), "ggplot")
  expect_s3_class(plot_marker_set_coverage(empty),  "ggplot")
})

# Utility shim: %||% isn't base R — use rlang if available, else a tiny local
`%||%` <- function(a, b) if (is.null(a)) b else a
