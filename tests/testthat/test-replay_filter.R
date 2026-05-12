# test-replay_filter.R - Document type expectations for the replay tuple key
#
# main.nf builds a typed Set of 6-tuples (species, group, nqtl, effect, h2, rep)
# from replay.tsv via parseManifestTuple (main.nf:90-99). The Groovy filter
# uses Set.contains() with strict type equality — any drift between manifest
# coercion and the channel-tuple types silently empties the replay filter and
# produces a no-op rerun.
#
# These R-side expectations mirror the Groovy assertions at main.nf:243-250.
# When parseManifestTuple changes, update this test in the same commit.

first_manifest_row <- function() {
  readr::read_tsv(
    fixture_path("replay-tsv-good.tsv"),
    col_types = readr::cols(
      session     = readr::col_character(),
      task_hash   = readr::col_character(),
      attempt     = readr::col_integer(),
      max_retries = readr::col_integer(),
      species     = readr::col_character(),
      group       = readr::col_character(),
      maf         = readr::col_double(),
      nqtl        = readr::col_character(),
      effect      = readr::col_character(),
      h2          = readr::col_double(),
      rep         = readr::col_integer(),
      mode        = readr::col_character(),
      type        = readr::col_character(),
      exit        = readr::col_integer()
    )
  ) |>
    dplyr::slice(1L)
}

test_that("manifest row coerces to canonical R types matching parseManifestTuple", {
  row <- first_manifest_row()

  # Index map mirrors main.nf:88 — 0=species 1=group 2=nqtl 3=effect 4=h2 5=rep
  expect_type(row$species, "character")  # Groovy String
  expect_type(row$group,   "character")  # Groovy String
  expect_type(row$nqtl,    "character")  # Groovy String — kept as string in main.nf
  expect_type(row$effect,  "character")  # Groovy String
  expect_type(row$h2,      "double")     # Groovy Float (R has no Float; double is the analogue)
  expect_type(row$rep,     "integer")    # Groovy Integer
})

test_that("the 6-tuple filter key is well-formed for a single manifest row", {
  row <- first_manifest_row()

  # Build the same 6-tuple the channel filter uses (main.nf:462).
  tuple_key <- list(
    species = row$species,
    group   = row$group,
    nqtl    = row$nqtl,
    effect  = row$effect,
    h2      = row$h2,
    rep     = row$rep
  )

  expect_true(all(!vapply(tuple_key, is.null, logical(1L))))
  expect_true(all(!vapply(tuple_key, function(x) length(x) != 1L, logical(1L))))
  expect_false(any(is.na(unlist(tuple_key, use.names = FALSE))))
})

test_that("nqtl is preserved as a string (not coerced to integer)", {
  # main.nf:247 asserts the channel-tuple nqtl is String. If readr ever
  # silently coerced "5" → 5L, the Set.contains() check would miss every row.
  row <- first_manifest_row()
  expect_identical(row$nqtl, "5")
  expect_false(is.numeric(row$nqtl))
})

test_that("h2 round-trips as a double across the manifest", {
  row <- first_manifest_row()
  expect_equal(row$h2, 0.3, tolerance = 1e-12)
})
