# test-failure_manifest.R - Unit tests for replay.tsv schema and round-trip parsing
#
# Validates the column set, ordering, and column types emitted by main.nf's
# `workflow.onComplete` block (the `columnOrder` array around main.nf:1096).
# These tests pin the schema so that adding/removing fields in NF_TRAP_PAYLOAD
# without a corresponding test update fails loudly.
#
# Runs as a plain unit test — no TEST_DB_DIR required.

# Canonical column order, mirroring main.nf:1096-1100. Update both in lockstep.
.replay_columns <- c(
  "session", "task_hash", "attempt", "max_retries",
  "species", "group", "maf", "nqtl", "effect", "h2", "rep", "mode", "type",
  "exit"
)

# Column types expected by db_clean_replay_slots.R (read_replay_manifest).
.replay_col_types <- function() {
  readr::cols(
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
}

# Mirrors main.nf:227-232: a manifest is well-formed only if the first line
# starts with "session\t" and contains every required filter key. Defined
# inline so the unit test does not depend on the runtime Rscript executing.
validate_replay_header <- function(path) {
  con <- file(path, "r")
  on.exit(close(con), add = TRUE)
  first <- readLines(con, n = 1L, warn = FALSE)
  if (length(first) == 0L) {
    stop(sprintf("replay.tsv header missing — file is empty: %s", path))
  }
  if (!startsWith(first, "session\t")) {
    stop(sprintf(
      "replay.tsv header does not start with 'session\\t' — file may be malformed: %s",
      path
    ))
  }
  cols <- strsplit(first, "\t", fixed = TRUE)[[1]]
  required <- c("species", "group", "nqtl", "effect", "h2", "rep")
  missing <- setdiff(required, cols)
  if (length(missing) > 0L) {
    stop(sprintf(
      "replay.tsv header missing required columns: %s",
      paste(missing, collapse = ", ")
    ))
  }
  invisible(cols)
}

read_replay_fixture <- function(name) {
  readr::read_tsv(fixture_path(name), col_types = .replay_col_types())
}

test_that("replay.tsv good fixture parses with expected columns and types", {
  df <- read_replay_fixture("replay-tsv-good.tsv")

  expect_equal(colnames(df), .replay_columns)

  expect_type(df$session,     "character")
  expect_type(df$task_hash,   "character")
  expect_type(df$attempt,     "integer")
  expect_type(df$max_retries, "integer")
  expect_type(df$species,     "character")
  expect_type(df$group,       "character")
  expect_type(df$maf,         "double")
  expect_type(df$nqtl,        "character") # nqtl is intentionally a string
  expect_type(df$effect,      "character")
  expect_type(df$h2,          "double")
  expect_type(df$rep,         "integer")
  expect_type(df$mode,        "character")
  expect_type(df$type,        "character")
  expect_type(df$exit,        "integer")
})

test_that("validate_replay_header accepts a well-formed manifest", {
  cols <- validate_replay_header(fixture_path("replay-tsv-good.tsv"))
  expect_equal(cols, .replay_columns)
})

test_that("validate_replay_header rejects a manifest with a malformed header", {
  expect_error(
    validate_replay_header(fixture_path("replay-tsv-bad-header.tsv")),
    regexp = "header"
  )
})

test_that("replay.tsv good fixture round-trips through readr", {
  path <- fixture_path("replay-tsv-good.tsv")
  df <- read_replay_fixture("replay-tsv-good.tsv")

  tmp <- tempfile(pattern = "replay_roundtrip_", fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  readr::write_tsv(df, tmp)

  df2 <- readr::read_tsv(tmp, col_types = .replay_col_types())
  expect_equal(df, df2)
})

test_that("replay.tsv good fixture has at least one failed slot", {
  df <- read_replay_fixture("replay-tsv-good.tsv")
  expect_gt(nrow(df), 0L)
  expect_true(all(df$exit != 0L), label = "all rows are recorded as failures")
})
