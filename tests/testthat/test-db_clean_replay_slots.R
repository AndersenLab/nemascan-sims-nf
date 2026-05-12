# test-db_clean_replay_slots.R - End-to-end test for the replay-slot cleanup script
#
# Builds a minimal fixture DB with two replicates of a single parameter cell,
# writes a synthetic replay manifest naming only one of them as failed, then
# invokes modules/local/db_clean_replay_slots/.../db_clean_replay_slots.R via
# Rscript. Verifies that exactly the failed slot's trait + mapping files are
# removed and the surviving slot is untouched.
#
# Runs as a plain unit test — no TEST_DB_DIR required.

# ── Constants matching the cleanup script's defaults ───────────────────────────
.cleanup_script <- file.path(
  .project_root, "modules", "local", "db_clean_replay_slots",
  "resources", "usr", "bin", "db_clean_replay_slots.R"
)
.replay_columns <- c(
  "session", "task_hash", "attempt", "max_retries",
  "species", "group", "maf", "nqtl", "effect", "h2", "rep", "mode", "type",
  "exit"
)
.fixture_species     <- "c_elegans"
.fixture_vcf_release <- "20220216"
.fixture_ms_ld       <- 0.8
.fixture_nqtl        <- "5"
.fixture_effect      <- "gamma"
.fixture_h2          <- 0.3
.fixture_cv_maf      <- 0.05
.fixture_cv_ld       <- 0.8

# ── Helpers ────────────────────────────────────────────────────────────────────

#' Compute the trait_id hash for a fixture parameter cell at a given rep.
.fixture_trait_hash <- function(group, maf, rep) {
  ms <- generate_marker_set_id(
    group, maf,
    .fixture_species, .fixture_vcf_release, .fixture_ms_ld
  )
  generate_trait_id(
    ms$hash, .fixture_nqtl, .fixture_effect, rep, .fixture_h2,
    .fixture_cv_maf, .fixture_cv_ld
  )$hash
}

#' Build the four trait-side parquet paths for a trait_id (mirrors the
#' cleanup script's trait_files_for_id()).
.trait_files <- function(db_root, trait_id) {
  c(
    fs::path(db_root, "traits",                  glue::glue("{trait_id}.parquet")),
    fs::path(db_root, "traits", "causal_variants",  glue::glue("{trait_id}_causal.parquet")),
    fs::path(db_root, "traits", "causal_genotypes", glue::glue("{trait_id}_causal_geno.parquet")),
    fs::path(db_root, "traits", "phenotypes",       glue::glue("{trait_id}_phenotype.parquet"))
  )
}

#' Build the eight mapping-side parquet paths for a trait_id (mirrors the
#' cleanup script's mapping_files_for_id()).
.mapping_files <- function(db_root, trait_id, group) {
  combos <- list(
    c("inbred", "TRUE"),  c("inbred", "FALSE"),
    c("loco",   "TRUE"),  c("loco",   "FALSE")
  )
  unlist(lapply(combos, function(mp) {
    m <- generate_mapping_id(trait_id, mp[[1]], as.logical(mp[[2]]))
    slot <- fs::path(
      db_root, "mappings",
      glue::glue("population={group}"),
      glue::glue("mapping_id={m$hash}")
    )
    c(fs::path(slot, "data.parquet"), fs::path(slot, "meta.parquet"))
  }))
}

#' Enumerate the trait + mapping parquets for a single (group, maf, rep) cell.
list_db_files_for_rep <- function(db_root, group, maf, rep) {
  trait_hash <- .fixture_trait_hash(group, maf, rep)
  list(
    trait   = .trait_files(db_root, trait_hash),
    mapping = .mapping_files(db_root, trait_hash, group)
  )
}

#' Materialize a minimal but schema-valid DB layout with marker_set_metadata
#' and one trait+mapping fan-out per rep.
setup_minimal_db_for_cleanup_test <- function(db_root, group, maf, reps) {
  fs::dir_create(fs::path(db_root, "marker_set_metadata"))
  fs::dir_create(fs::path(db_root, "traits", "causal_variants"))
  fs::dir_create(fs::path(db_root, "traits", "causal_genotypes"))
  fs::dir_create(fs::path(db_root, "traits", "phenotypes"))
  fs::dir_create(fs::path(db_root, "mappings"))

  ms <- generate_marker_set_id(
    group, maf,
    .fixture_species, .fixture_vcf_release, .fixture_ms_ld
  )

  write_marker_set_metadata(
    population             = group,
    maf                    = maf,
    species                = .fixture_species,
    vcf_release_id         = .fixture_vcf_release,
    ms_ld                  = .fixture_ms_ld,
    n_markers              = 100L,
    marker_set_id          = ms$hash,
    marker_set_hash_string = ms$hash_string,
    n_independent_tests    = 80,
    eigen_source_file      = NA_character_,
    strainfile_hash        = "deadbeef",
    strain_list            = "A,B,C",
    base_dir               = db_root
  )

  empty_df <- data.frame(x = integer(0))
  for (r in reps) {
    files <- list_db_files_for_rep(db_root, group, maf, r)
    purrr::walk(files$trait, function(p) {
      arrow::write_parquet(empty_df, p)
    })
    purrr::walk(files$mapping, function(p) {
      fs::dir_create(fs::path_dir(p))
      arrow::write_parquet(empty_df, p)
    })
  }

  invisible(db_root)
}

#' Build a one-row manifest tibble for a single failed slot (14-column schema).
make_single_slot_manifest <- function(group, rep, maf, species = .fixture_species) {
  tibble::tibble(
    session     = "00000000-0000-0000-0000-000000000001",
    task_hash   = "abcd1234",
    attempt     = 4L,
    max_retries = 3L,
    species     = species,
    group       = group,
    maf         = maf,
    nqtl        = .fixture_nqtl,
    effect      = .fixture_effect,
    h2          = .fixture_h2,
    rep         = as.integer(rep),
    mode        = "inbred",
    type        = "pca",
    exit        = 1L
  )
}

#' Write a manifest tibble to a fresh tempfile and return the path.
write_temp_manifest <- function(manifest) {
  path <- tempfile(pattern = "replay_", fileext = ".tsv")
  readr::write_tsv(manifest, path)
  path
}

#' Invoke an Rscript with R_SOURCE_DIR pre-set; raises on non-zero exit.
Rscript_call <- function(script_path, args = character()) {
  # suppressWarnings: system2 emits an R warning whenever the child returns
  # non-zero, even when the caller is intentionally probing failure modes
  # (we re-raise via stop() below using the captured status).
  output <- suppressWarnings(system2(
    "Rscript",
    args = c(script_path, args),
    env  = paste0("R_SOURCE_DIR=", file.path(.project_root, "R")),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("Rscript failed (status %d):\n%s",
                 status, paste(output, collapse = "\n")))
  }
  invisible(output)
}

skip_if_no_cleanup_script <- function() {
  if (!file.exists(.cleanup_script)) {
    skip("db_clean_replay_slots.R script not found")
  }
}

# ── Tests ──────────────────────────────────────────────────────────────────────

test_that("DB_CLEAN_REPLAY_SLOTS removes only the failed slot's files", {
  skip_if_no_cleanup_script()

  tmp_db <- tempfile(pattern = "replay_clean_")
  fs::dir_create(tmp_db)
  withr::defer(fs::dir_delete(tmp_db))

  group <- "ce.test"
  maf   <- 0.05

  setup_minimal_db_for_cleanup_test(tmp_db, group = group, maf = maf,
                                    reps = c(1L, 2L))

  manifest_path <- make_single_slot_manifest(group, rep = 2L, maf = maf) |>
    write_temp_manifest()
  withr::defer(unlink(manifest_path))

  Rscript_call(.cleanup_script,
               args = c("--replay_tsv", manifest_path, "--db_root", tmp_db))

  rep1 <- list_db_files_for_rep(tmp_db, group, maf, 1L)
  rep2 <- list_db_files_for_rep(tmp_db, group, maf, 2L)

  expect_true(all(fs::file_exists(rep1$trait)),
              label = "rep=1 trait files survive cleanup")
  expect_true(all(fs::file_exists(rep1$mapping)),
              label = "rep=1 mapping files survive cleanup")
  expect_false(any(fs::file_exists(rep2$trait)),
               label = "rep=2 trait files removed")
  expect_false(any(fs::file_exists(rep2$mapping)),
               label = "rep=2 mapping files removed")

  expect_true(fs::dir_exists(fs::path(tmp_db, "mappings",
                                      glue::glue("population={group}"))),
              label = "population= directory survives (no directory-level deletes)")
})

test_that("DB_CLEAN_REPLAY_SLOTS is idempotent", {
  skip_if_no_cleanup_script()

  tmp_db <- tempfile(pattern = "replay_clean_idem_")
  fs::dir_create(tmp_db)
  withr::defer(fs::dir_delete(tmp_db))

  group <- "ce.test"
  maf   <- 0.05

  setup_minimal_db_for_cleanup_test(tmp_db, group = group, maf = maf,
                                    reps = c(1L, 2L))

  manifest_path <- make_single_slot_manifest(group, rep = 2L, maf = maf) |>
    write_temp_manifest()
  withr::defer(unlink(manifest_path))

  expect_no_error(
    Rscript_call(.cleanup_script,
                 args = c("--replay_tsv", manifest_path, "--db_root", tmp_db))
  )
  expect_no_error(
    Rscript_call(.cleanup_script,
                 args = c("--replay_tsv", manifest_path, "--db_root", tmp_db))
  )

  # rep=1 must still be intact after the second invocation.
  rep1 <- list_db_files_for_rep(tmp_db, group, maf, 1L)
  expect_true(all(fs::file_exists(rep1$trait)))
  expect_true(all(fs::file_exists(rep1$mapping)))
})

test_that("DB_CLEAN_REPLAY_SLOTS errors loudly when marker_set_metadata is missing", {
  skip_if_no_cleanup_script()

  tmp_db <- tempfile(pattern = "replay_clean_missing_")
  fs::dir_create(tmp_db)
  withr::defer(fs::dir_delete(tmp_db))

  fs::dir_create(fs::path(tmp_db, "marker_set_metadata"))
  fs::dir_create(fs::path(tmp_db, "traits"))
  fs::dir_create(fs::path(tmp_db, "mappings"))

  manifest_path <- make_single_slot_manifest("ce.test", rep = 2L, maf = 0.05) |>
    write_temp_manifest()
  withr::defer(unlink(manifest_path))

  expect_error(
    Rscript_call(.cleanup_script,
                 args = c("--replay_tsv", manifest_path, "--db_root", tmp_db)),
    regexp = "marker_set_metadata|status"
  )
})
