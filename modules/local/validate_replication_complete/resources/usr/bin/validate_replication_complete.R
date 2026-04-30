#!/usr/bin/env Rscript

library(optparse)

opts <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--db_root",       type = "character"),
    optparse::make_option("--expected_reps", type = "integer"),
    optparse::make_option("--output_dir",    type = "character")
)))

r_source_dir <- Sys.getenv("R_SOURCE_DIR", unset = "")
if (r_source_dir == "") stop("R_SOURCE_DIR environment variable must be set")
source(file.path(r_source_dir, "utils.R"))
source(file.path(r_source_dir, "database.R"))
source(file.path(r_source_dir, "queries.R"))

#' Query parameter cells with fewer distinct reps than expected
find_short_cells <- function(con, expected_reps) {
    DBI::dbGetQuery(con, "
        SELECT group_, maf, nqtl, effect, h2, mode, type,
               COUNT(DISTINCT rep) AS rep_count
        FROM   mappings_metadata
        GROUP  BY group_, maf, nqtl, effect, h2, mode, type
        HAVING rep_count <> ?
    ", params = list(expected_reps)) |>
        tibble::as_tibble()
}

#' Identify unreadable or empty parquet files under a directory
find_corrupt_parquet <- function(db_root) {
    data_files <- fs::dir_ls(
        fs::path(db_root, "mappings"),
        recurse = TRUE,
        glob    = "*/data.parquet"
    )

    is_bad <- purrr::map_lgl(data_files, function(f) {
        tryCatch(
            nrow(arrow::open_dataset(f)) == 0L,
            error = function(e) TRUE
        )
    })

    data_files[is_bad]
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
            glue::glue("{length(corrupt)} corrupt or empty data.parquet file(s):"),
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
            "Replication check OK: every parameter cell has {expected_reps} reps; ",
            "all parquet files readable."
        ))
        quit(status = 0L)
    } else {
        writeLines(report, marker)
        message(paste(report, collapse = "\n"))
        quit(status = 1L)
    }
}

# --- main ---

con         <- open_mapping_db(opts$db_root)
short_cells <- find_short_cells(con, opts$expected_reps)
corrupt     <- find_corrupt_parquet(opts$db_root)

ok     <- nrow(short_cells) == 0L && length(corrupt) == 0L
marker <- fs::path(opts$output_dir, "REPLAY_REQUIRED")
report <- build_failure_report(short_cells, corrupt, opts$expected_reps, opts$output_dir)

write_outcome(ok, marker, report, opts$expected_reps)
