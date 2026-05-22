#!/usr/bin/env Rscript

library(optparse)

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--db_root",             type = "character"),
    optparse::make_option("--expected_reps",       type = "integer"),
    optparse::make_option("--filter_pool_metrics", type = "character", default = NA_character_,
                          help = "filter_pool_metrics.tsv with per-cell effective_reps (cap-aware)"),
    optparse::make_option("--output_dir",          type = "character")
)))

r_source_dir <- Sys.getenv("R_SOURCE_DIR", unset = "")
if (r_source_dir == "") stop("R_SOURCE_DIR environment variable must be set")
source(file.path(r_source_dir, "utils.R"))
source(file.path(r_source_dir, "database.R"))
source(file.path(r_source_dir, "queries.R"))

#' Query parameter cells with fewer distinct reps than a single global expected count.
#' Legacy fallback used only when no filter_pool_metrics is available (e.g. databases written
#' before per-cell rep metrics existed).
find_short_cells <- function(con, expected_reps) {
    DBI::dbGetQuery(con, "
        SELECT population, maf, nqtl, effect, h2, algorithm, pca,
               COUNT(DISTINCT rep) AS rep_count
        FROM   metadata
        GROUP  BY population, maf, nqtl, effect, h2, algorithm, pca
        HAVING rep_count <> ?
    ", params = list(expected_reps)) |>
        tibble::as_tibble()
}

#' Is a usable filter_pool_metrics path available?
have_metrics <- function(metrics_path) {
    !is.null(metrics_path) && !is.na(metrics_path) && nzchar(metrics_path) &&
        fs::file_exists(metrics_path)
}

#' Cap-aware short-cell detection.
#'
#' Compares each (group, ms_maf, filter_id, nqtl, h2) trait cell's distinct-rep count
#' against its per-cell effective_reps from the metrics TSV (the single source of truth —
#' combinatorics are never recomputed here) instead of one global integer. Both the trait
#' metadata (traits/*.parquet — the only table carrying cv_region_filter) and the metrics
#' TSV are read by DuckDB so the join keys coerce consistently. Returns cells with
#' rep_count < effective_reps, annotated with the resolve status.
find_short_cells_capaware <- function(db_root, metrics_path) {
    traits_dir <- fs::path(db_root, "traits")
    if (!fs::dir_exists(traits_dir) ||
        length(fs::dir_ls(traits_dir, glob = "*.parquet", fail = FALSE)) == 0L) {
        return(tibble::tibble())
    }
    traits_glob <- fs::path(traits_dir, "*.parquet")
    con <- duckdb::dbConnect(duckdb::duckdb(), read_only = FALSE)
    on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)
    sql <- sprintf("
        WITH actual AS (
            SELECT CAST(population AS VARCHAR)        AS \"group\",
                   CAST(maf AS DOUBLE)                AS ms_maf,
                   CAST(cv_region_filter AS VARCHAR)  AS filter_id,
                   CAST(nqtl AS BIGINT)               AS nqtl,
                   CAST(h2 AS DOUBLE)                 AS h2,
                   COUNT(DISTINCT rep)                AS rep_count
            FROM read_parquet('%s', union_by_name = true)
            GROUP BY 1, 2, 3, 4, 5
        ),
        plan AS (
            SELECT CAST(\"group\" AS VARCHAR)         AS \"group\",
                   CAST(ms_maf AS DOUBLE)             AS ms_maf,
                   CAST(filter_id AS VARCHAR)         AS filter_id,
                   CAST(nqtl AS BIGINT)               AS nqtl,
                   CAST(effective_reps AS BIGINT)     AS effective_reps,
                   CAST(status AS VARCHAR)            AS status
            FROM read_csv_auto('%s', delim = '\t', header = true)
        )
        SELECT a.\"group\", a.ms_maf, a.filter_id, a.nqtl, a.h2,
               a.rep_count, p.effective_reps, p.status
        FROM actual a
        JOIN plan p USING (\"group\", ms_maf, filter_id, nqtl)
        WHERE a.rep_count < p.effective_reps
        ORDER BY 1, 2, 3, 4, 5
    ", traits_glob, metrics_path)
    DBI::dbGetQuery(con, sql) |> tibble::as_tibble()
}

#' Read intentionally-capped cells (status != ok) from the metrics TSV for an
#' informational note. Returns an empty tibble on any read error.
read_capped_cells <- function(metrics_path) {
    tryCatch({
        m <- utils::read.delim(metrics_path, sep = "\t", stringsAsFactors = FALSE)
        tibble::as_tibble(m[!is.na(m$status) & m$status != "ok", , drop = FALSE])
    }, error = function(e) tibble::tibble())
}

#' Identify mapping partitions whose data.parquet is missing or empty.
#'
#' Glob the expected data.parquet files and read each Parquet footer's row
#' count via arrow::ParquetFileReader. The footer is the last few KB of the
#' file; this avoids materializing any row data, schema unification, or
#' Dataset object construction. Files with num_rows == 0L (or that fail to
#' open) are flagged — matches the legacy "data.parquet present but empty
#' → bad" semantics and returns real file paths. We deliberately bypass
#' arrow::open_dataset()'s Hive partition parsing because, when called with
#' a vector of file paths rather than a directory root, R-arrow does not
#' surface partition columns from the path segments — group_by(population,
#' mapping_id) then fails with "Column not found" (Patch v2 regression).
find_corrupt_parquet <- function(db_root) {
    mappings_root <- fs::path(db_root, "mappings")
    if (!fs::dir_exists(mappings_root)) return(character(0))

    data_files <- as.character(fs::dir_ls(
        mappings_root, recurse = TRUE, glob = "*/data.parquet"
    ))
    if (length(data_files) == 0L) return(character(0))

    n_rows <- vapply(data_files, function(f) {
        tryCatch(
            as.integer(arrow::ParquetFileReader$create(f)$num_rows),
            error = function(e) 0L
        )
    }, integer(1))

    data_files[n_rows == 0L]
}

#' Format a human-readable failure report
build_failure_report <- function(short_cells, corrupt, expected_reps, output_dir) {
    msgs <- character(0)

    if (nrow(short_cells) > 0L) {
        msgs <- c(
            msgs,
            glue::glue("{nrow(short_cells)} parameter cell(s) short of {expected_reps} reps:"),
            utils::capture.output(print(short_cells, n = Inf))
        )
    }

    if (length(corrupt) > 0L) {
        msgs <- c(
            msgs,
            glue::glue("{length(corrupt)} corrupt or empty mapping partition(s):"),
            paste0("  ", corrupt)
        )
    }

    c(msgs, glue::glue(
        "Run with --replay {output_dir}/replay.tsv -resume to recover, then re-validate."
    ))
}

#' Write the REPLAY_REQUIRED marker file and exit non-zero, or clear it and exit 0
write_outcome <- function(ok, marker, report, expected_reps) {
    if (ok) {
        if (fs::file_exists(marker)) fs::file_delete(marker)
        message(glue::glue(
            "Replication check OK: all data.parquet files readable.",
            " ({expected_reps} reps expected per parameter cell.)"
        ))
        quit(status = 0L)
    } else {
        writeLines(report, marker)
        message(paste(report, collapse = "\n"))
        quit(status = 1L)
    }
}

#' Build the recovery hint shown in REPLAY_REQUIRED when the DB is empty/absent.
#' Branches on whether replay.tsv was actually produced — pointing `--replay` at
#' a non-existent manifest is unactionable.
build_recovery_hint <- function(output_dir) {
    replay_tsv <- fs::path(output_dir, "replay.tsv")
    if (fs::file_exists(replay_tsv)) {
        glue::glue("Run with --replay {replay_tsv} -resume to recover, then re-validate.")
    } else {
        paste(
            "No replay manifest was produced — `.failures/` is empty.",
            glue::glue("Inspect upstream task logs in {output_dir}/.command.* (or the workdir under workflow.workDir) to determine why the failure trap did not fire.")
        )
    }
}

# --- main ---

# Guard: open_mapping_db() raises "Database directory not found" if db_root is
# missing — opaque to the user. The directory is genuinely absent when every
# upstream task fails before any DB write completes (e.g. all
# GCTA_SIMULATE_PHENOTYPES tasks failing with errorStrategy='ignore', leaving
# no path to WRITE_GWA_TO_DB / WRITE_TRAIT_DATA / AGGREGATE_METADATA), or when
# a resumed run with cleared output finds no published artifacts. Treat this
# the same as the metadata-table-missing case below.
if (!fs::dir_exists(opts$db_root)) {
    hint <- build_recovery_hint(opts$output_dir)
    message(glue::glue("ERROR: database directory not found: {opts$db_root}"))
    message("No DB writes occurred — likely all upstream tasks failed before any data could be written.")
    message(hint)
    marker <- fs::path(opts$output_dir, "REPLAY_REQUIRED")
    writeLines(c(
        glue::glue("Database directory not found: {opts$db_root}"),
        "No DB writes occurred — all upstream tasks may have failed before any marker/trait/mapping data could be written.",
        hint
    ), marker)
    quit(status = 1L)
}

con <- open_mapping_db(opts$db_root)

if (!DBI::dbExistsTable(con, "metadata")) {
    hint <- build_recovery_hint(opts$output_dir)
    message("ERROR: metadata view missing — mappings_metadata.parquet was never written.")
    message("This indicates all upstream mapping tasks failed. Treating all replication cells as incomplete.")
    message(hint)
    marker <- fs::path(opts$output_dir, "REPLAY_REQUIRED")
    writeLines(c(
        "metadata view missing — no mapping data written to database.",
        "All upstream tasks may have failed with errorStrategy='ignore'.",
        hint
    ), marker)
    quit(status = 1L)
}

# Cap-aware replication check when filter_pool_metrics is available: each cell is compared
# to its own effective_reps (the number of reps that could be drawn given the available
# causal-variant combinations) rather than one global target, so cells that were
# intentionally limited are not mistaken for genuinely-dropped reps. Falls back to the
# legacy single-target check when no metrics file is present (e.g. databases written before
# the metrics file existed).
corrupt <- find_corrupt_parquet(opts$db_root)

if (have_metrics(opts$filter_pool_metrics)) {
    # Intentionally-limited cells (fewer reps than requested because the region/nqtl pool is
    # too small): report them so the operator knows replication was capped on purpose, not lost.
    capped <- read_capped_cells(opts$filter_pool_metrics)
    if (nrow(capped) > 0L) {
        message(glue::glue(
            "NOTE: {nrow(capped)} cell(s) were intentionally rep-capped (status != ok). ",
            "Reduced replication here is expected, not a dropped task — see {opts$filter_pool_metrics}."
        ))
        message(paste(utils::capture.output(print(capped, n = Inf)), collapse = "\n"))
    }

    short_cells <- tryCatch(
        find_short_cells_capaware(opts$db_root, opts$filter_pool_metrics),
        error = function(e) {
            message("WARNING: cap-aware replication check failed (", conditionMessage(e),
                    "); skipping per-cell rep comparison.")
            tibble::tibble()
        }
    )
    if (nrow(short_cells) > 0L) {
        message(glue::glue(
            "WARNING: {nrow(short_cells)} trait cell(s) have fewer reps than their per-cell ",
            "effective_reps target. Treating as a warning — pipeline will continue."
        ))
        message(paste(utils::capture.output(print(short_cells, n = Inf)), collapse = "\n"))
    }
} else {
    short_cells <- find_short_cells(con, opts$expected_reps)
    if (nrow(short_cells) > 0L) {
        message(glue::glue(
            "WARNING: {nrow(short_cells)} parameter cell(s) have fewer than {opts$expected_reps} reps. ",
            "Treating as a warning — pipeline will continue."
        ))
        message(paste(utils::capture.output(print(short_cells, n = Inf)), collapse = "\n"))
    }
}

ok     <- length(corrupt) == 0L
marker <- fs::path(opts$output_dir, "REPLAY_REQUIRED")
report <- build_failure_report(tibble::tibble(), corrupt, opts$expected_reps, opts$output_dir)

write_outcome(ok, marker, report, opts$expected_reps)
