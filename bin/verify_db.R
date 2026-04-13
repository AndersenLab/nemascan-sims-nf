#!/usr/bin/env Rscript
# verify_db.R — Verify completeness of the nemascan-sims-nf Parquet database.
#
# Two modes:
#   expectation: user supplies pipeline architecture args (strainfile, nqtl, h2,
#     effect, reps, cv_maf, cv_ld). The script compares expected file/row counts
#     to what is present in {outputDir}/db/ and enumerates missing combinations.
#   discovery:   no architecture args supplied. Grid is derived from
#     mappings_metadata.parquet and only observed distributions are reported.
#
# verify_db() never errors out. Failures are collected in result$errors and
# result$status is one of "complete" / "incomplete" / "discovery" / "error".
# Plot helpers also never throw — a failed plot returns a ggplot whose title
# describes the error, so qmd chunks that unconditionally print never break.
#
# CLI:
#   Rscript bin/verify_db.R --base_dir path/to/db \
#     [--strainfile s.tsv --nqtl 1,2,3 --h2 0.2,0.5 --effect gamma \
#      --reps 25 --cv_maf 0.05 --cv_ld 0.8]
#
# Interactive use (from R or a qmd):
#   source("bin/verify_db.R")
#   result <- verify_db(base_dir, strainfile = "...", nqtl = c(1, 2, 5), ...)
#   print(plot_n_mappings_by_param(result))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(glue)
  library(DBI)
  library(duckdb)
  library(arrow)
  library(ggplot2)
  library(readr)
})

# ---------------------------------------------------------------------------
# R library sourcing — skipped when verify_db is sourced from a qmd where
# R/setup.R was already loaded.
# ---------------------------------------------------------------------------

.vdb_locate_r_library <- function() {
  r_lib <- Sys.getenv("NEMASCAN_R_LIB", unset = "")
  if (nzchar(r_lib) && dir.exists(r_lib)) return(r_lib)

  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
                        error = function(e) NULL)
  if (is.null(this_file) || !nzchar(this_file) || !file.exists(this_file)) {
    cli_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", cli_args, value = TRUE)
    if (length(file_arg) == 1L) {
      this_file <- tryCatch(normalizePath(sub("^--file=", "", file_arg), mustWork = FALSE),
                            error = function(e) NULL)
    }
  }
  if (!is.null(this_file) && file.exists(this_file)) {
    candidate <- normalizePath(file.path(dirname(this_file), "..", "R"), mustWork = FALSE)
    if (dir.exists(candidate)) return(candidate)
  }
  ""
}

.vdb_source_r_library <- function() {
  if (exists("build_ids_from_params", mode = "function") &&
      exists("open_mapping_db",        mode = "function") &&
      exists("get_all_marker_set_metadata", mode = "function")) {
    return(invisible(NULL))
  }
  r_lib <- .vdb_locate_r_library()
  if (!nzchar(r_lib)) {
    stop("Cannot locate R library; set NEMASCAN_R_LIB=<repo>/R before sourcing verify_db.R")
  }
  setup_file <- file.path(r_lib, "setup.R")
  if (file.exists(setup_file)) {
    source(setup_file)
  } else {
    for (f in c("utils.R", "io.R", "database.R", "queries.R")) {
      p <- file.path(r_lib, f)
      if (file.exists(p)) source(p)
    }
  }
  invisible(NULL)
}
.vdb_source_r_library()


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

.VDB_DEFAULT_ENUM_LIMIT <- 500000L
.VDB_STATUS_LEVELS <- c("complete", "incomplete", "discovery", "error")


# ---------------------------------------------------------------------------
# Error collection helpers
# ---------------------------------------------------------------------------

.vdb_new_errors <- function() {
  e <- new.env(parent = emptyenv())
  e$rows <- list()
  e
}

.vdb_record <- function(errors_env, step, message, level = "error") {
  errors_env$rows[[length(errors_env$rows) + 1L]] <- tibble::tibble(
    step      = step,
    level     = level,
    message   = as.character(message),
    timestamp = Sys.time()
  )
  invisible(NULL)
}

.vdb_errors_to_tibble <- function(errors_env) {
  if (length(errors_env$rows) == 0L) {
    return(tibble::tibble(
      step      = character(),
      level     = character(),
      message   = character(),
      timestamp = as.POSIXct(character())
    ))
  }
  dplyr::bind_rows(errors_env$rows)
}

.vdb_safely <- function(expr, step, errors_env, default = NULL) {
  tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        .vdb_record(errors_env, step, conditionMessage(w), "warning")
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      .vdb_record(errors_env, step, conditionMessage(e), "error")
      default
    }
  )
}


# ---------------------------------------------------------------------------
# Strainfile parsing
# ---------------------------------------------------------------------------

.vdb_read_strainfile <- function(path, errors_env) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) {
    .vdb_record(errors_env, "read_strainfile",
                glue("Strainfile not found: {path}"), "error")
    return(NULL)
  }
  df <- .vdb_safely(
    readr::read_tsv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      progress  = FALSE
    ),
    "read_strainfile",
    errors_env
  )
  if (is.null(df) || nrow(df) == 0L) {
    .vdb_record(errors_env, "read_strainfile",
                "Strainfile is empty or unreadable", "error")
    return(NULL)
  }

  # Accept either {group, ms_maf} or {population, maf} naming
  if ("group" %in% names(df) && !"population" %in% names(df)) {
    df$population <- df$group
  }
  if ("ms_maf" %in% names(df) && !"maf" %in% names(df)) {
    df$maf <- df$ms_maf
  }
  if (!all(c("population", "maf") %in% names(df))) {
    .vdb_record(errors_env, "read_strainfile",
                "Strainfile missing required columns 'group'/'population' and 'ms_maf'/'maf'",
                "error")
    return(NULL)
  }

  df$maf <- suppressWarnings(as.numeric(df$maf))
  if ("ms_ld" %in% names(df)) {
    df$ms_ld <- suppressWarnings(as.numeric(df$ms_ld))
  }

  keep_cols <- intersect(
    c("population", "maf", "species", "vcf", "ms_ld", "strains"),
    names(df)
  )
  df[, keep_cols, drop = FALSE]
}


# ---------------------------------------------------------------------------
# Merge strainfile groups with marker_set_metadata to hydrate
# species / vcf_release_id / ms_ld needed for hash generation
# ---------------------------------------------------------------------------

.vdb_hydrate_groups <- function(strainfile_df, base_dir, errors_env) {
  if (is.null(strainfile_df) || nrow(strainfile_df) == 0L) return(NULL)

  all_meta <- .vdb_safely(
    get_all_marker_set_metadata(base_dir),
    "get_all_marker_set_metadata",
    errors_env,
    default = NULL
  )

  needed <- c("population", "maf", "species", "vcf_release_id", "ms_ld",
              "marker_set_id", "strainfile_hash", "strain_list",
              "n_markers", "n_independent_tests")
  if (is.null(all_meta) || !is.data.frame(all_meta) || nrow(all_meta) == 0L) {
    meta_sel <- data.frame(
      population          = character(),
      maf                 = numeric(),
      species             = character(),
      vcf_release_id      = character(),
      ms_ld               = numeric(),
      marker_set_id       = character(),
      strainfile_hash     = character(),
      strain_list         = character(),
      n_markers           = integer(),
      n_independent_tests = numeric(),
      stringsAsFactors    = FALSE
    )
  } else {
    for (col in needed) {
      if (!col %in% names(all_meta)) all_meta[[col]] <- NA
    }
    meta_sel <- all_meta[, needed, drop = FALSE]
    meta_sel$maf   <- suppressWarnings(as.numeric(meta_sel$maf))
    meta_sel$ms_ld <- suppressWarnings(as.numeric(meta_sel$ms_ld))
  }

  # Strainfile wins for join keys; metadata fills in species/vcf_release_id/ms_ld.
  overlap <- intersect(names(strainfile_df),
                       setdiff(needed, c("population", "maf", "ms_ld")))
  strainfile_df <- dplyr::select(strainfile_df, -dplyr::any_of(overlap))

  has_strainfile_ms_ld <- "ms_ld" %in% names(strainfile_df)
  joined <- dplyr::left_join(strainfile_df, meta_sel,
                             by = c("population", "maf"),
                             suffix = c("", ".db"))

  if (has_strainfile_ms_ld && "ms_ld.db" %in% names(joined)) {
    joined$ms_ld <- ifelse(is.na(joined$ms_ld), joined$ms_ld.db, joined$ms_ld)
    joined$ms_ld.db <- NULL
  } else if (!has_strainfile_ms_ld && "ms_ld.db" %in% names(joined)) {
    joined$ms_ld <- joined$ms_ld.db
    joined$ms_ld.db <- NULL
  }

  joined$marker_set_metadata_present <- !is.na(joined$marker_set_id)
  tibble::as_tibble(joined)
}


# ---------------------------------------------------------------------------
# DB open + observed queries
# ---------------------------------------------------------------------------

.vdb_open_db <- function(base_dir, errors_env) {
  if (!dir.exists(base_dir)) {
    .vdb_record(errors_env, "open_db",
                glue("Database directory not found: {base_dir}"), "error")
    return(NULL)
  }
  .vdb_safely(
    open_mapping_db(base_dir, read_only = TRUE),
    "open_mapping_db",
    errors_env
  )
}

.vdb_has_view <- function(con, view_name, errors_env) {
  if (is.null(con)) return(FALSE)
  views <- .vdb_safely(
    DBI::dbGetQuery(
      con,
      "SELECT table_name FROM information_schema.tables WHERE table_schema = 'main'"
    ),
    "list_views",
    errors_env,
    default = data.frame(table_name = character())
  )
  view_name %in% views$table_name
}

.vdb_query_observed <- function(con, errors_env) {
  if (is.null(con)) return(NULL)
  if (!.vdb_has_view(con, "metadata", errors_env)) {
    .vdb_record(errors_env, "query_observed",
                "metadata view not present (mappings_metadata.parquet missing)",
                "warning")
    return(NULL)
  }
  .vdb_safely(
    DBI::dbGetQuery(con, "
      SELECT
        population,
        nqtl,
        h2,
        effect,
        algorithm,
        pca,
        COUNT(*)                 AS n_mappings,
        COUNT(DISTINCT trait_id) AS n_distinct_traits,
        COUNT(DISTINCT rep)      AS n_distinct_rep
      FROM metadata
      GROUP BY population, nqtl, h2, effect, algorithm, pca
    "),
    "query_observed_mappings",
    errors_env
  )
}

.vdb_query_observed_ids <- function(con, errors_env) {
  if (is.null(con)) return(NULL)
  if (!.vdb_has_view(con, "metadata", errors_env)) return(NULL)
  .vdb_safely(
    DBI::dbGetQuery(con, "SELECT DISTINCT trait_id, mapping_id FROM metadata"),
    "query_observed_ids",
    errors_env
  )
}


# ---------------------------------------------------------------------------
# Derive an observed parameter grid from the aggregated observed table
# ---------------------------------------------------------------------------

.vdb_derive_grid_from_db <- function(observed) {
  if (is.null(observed) || nrow(observed) == 0L) {
    return(list(
      populations = character(0),
      nqtl        = integer(0),
      h2          = numeric(0),
      effect      = character(0),
      algorithm   = character(0),
      pca         = logical(0),
      reps_max    = NA_integer_
    ))
  }
  list(
    populations = sort(unique(observed$population)),
    nqtl        = sort(unique(as.integer(observed$nqtl))),
    h2          = sort(unique(as.numeric(observed$h2))),
    effect      = sort(unique(as.character(observed$effect))),
    algorithm   = sort(unique(as.character(observed$algorithm))),
    pca         = sort(unique(as.logical(observed$pca))),
    reps_max    = max(as.integer(observed$n_distinct_rep), na.rm = TRUE)
  )
}


# ---------------------------------------------------------------------------
# Expected enumeration (expectation mode only)
# ---------------------------------------------------------------------------

.vdb_enumerate_expected <- function(groups_df, nqtl, h2, effect, reps,
                                    cv_maf, cv_ld, algorithm, pca,
                                    max_rows, errors_env) {
  n_groups <- nrow(groups_df)
  if (n_groups == 0L) return(NULL)

  n_rows_estimate <- n_groups * length(nqtl) * length(h2) *
                     length(effect) * as.integer(reps) *
                     length(algorithm) * length(pca)
  if (n_rows_estimate > max_rows) {
    .vdb_record(errors_env, "enumerate_expected",
      glue("Cartesian product too large ({n_rows_estimate} rows > ",
           "max_enumeration_rows={max_rows}); skipping per-ID enumeration. ",
           "Aggregate-level checks still produced."),
      "warning")
    return(NULL)
  }

  groups_df$.vdb_group_row <- seq_len(n_groups)
  base_grid <- tidyr::crossing(
    .vdb_group_row = groups_df$.vdb_group_row,
    nqtl           = as.integer(nqtl),
    h2             = as.numeric(h2),
    effect         = as.character(effect),
    rep            = seq_len(as.integer(reps)),
    algorithm      = as.character(algorithm),
    pca            = as.logical(pca)
  )

  join_cols <- intersect(
    c("population", "maf", "species", "vcf_release_id", "ms_ld",
      "marker_set_id", "marker_set_metadata_present"),
    names(groups_df)
  )
  grid <- dplyr::left_join(
    base_grid,
    groups_df[, c(".vdb_group_row", join_cols), drop = FALSE],
    by = ".vdb_group_row"
  )
  grid$.vdb_group_row <- NULL

  trait_ids   <- character(nrow(grid))
  mapping_ids <- character(nrow(grid))
  n_failed    <- 0L
  failed_msgs <- character()

  for (i in seq_len(nrow(grid))) {
    row <- grid[i, , drop = FALSE]
    if (isFALSE(row$marker_set_metadata_present)) {
      trait_ids[i]   <- NA_character_
      mapping_ids[i] <- NA_character_
      next
    }
    params <- list(
      population       = row$population,
      maf              = row$maf,
      species          = row$species,
      vcf_release_id   = row$vcf_release_id,
      ms_ld            = row$ms_ld,
      nqtl             = row$nqtl,
      effect           = row$effect,
      rep              = row$rep,
      h2               = row$h2,
      cv_maf_effective = cv_maf,
      cv_ld            = cv_ld
    )
    ids <- tryCatch(
      build_ids_from_params(params, mode = row$algorithm, pca = row$pca),
      error = function(e) {
        n_failed   <<- n_failed + 1L
        if (length(failed_msgs) < 3L) {
          failed_msgs <<- c(failed_msgs, conditionMessage(e))
        }
        NULL
      }
    )
    if (is.null(ids)) {
      trait_ids[i]   <- NA_character_
      mapping_ids[i] <- NA_character_
    } else {
      trait_ids[i]   <- ids$trait_id$hash
      mapping_ids[i] <- ids$mapping_id$hash
    }
  }

  if (n_failed > 0L) {
    .vdb_record(errors_env, "enumerate_expected",
      glue("{n_failed} grid rows failed ID computation. ",
           "First error(s): {paste(failed_msgs, collapse=' | ')}"),
      "warning")
  }

  grid$trait_id   <- trait_ids
  grid$mapping_id <- mapping_ids
  grid <- dplyr::select(grid,
    dplyr::any_of(c("population", "maf", "species", "vcf_release_id", "ms_ld",
                    "nqtl", "h2", "effect", "rep", "algorithm", "pca",
                    "trait_id", "mapping_id", "marker_set_metadata_present"))
  )
  tibble::as_tibble(grid)
}


# ---------------------------------------------------------------------------
# Expected/observed join producing delta tables
# ---------------------------------------------------------------------------

.vdb_aggregate_expected <- function(groups_df, nqtl, h2, effect, reps,
                                    algorithm, pca) {
  if (is.null(groups_df) || nrow(groups_df) == 0L) return(NULL)
  expected <- tidyr::crossing(
    population = groups_df$population,
    nqtl       = as.integer(nqtl),
    h2         = as.numeric(h2),
    effect     = as.character(effect),
    algorithm  = as.character(algorithm),
    pca        = as.logical(pca)
  )
  expected$expected <- as.integer(reps)
  tibble::as_tibble(expected)
}

.vdb_build_mappings_delta <- function(expected_agg, observed) {
  if (is.null(expected_agg)) return(NULL)
  obs <- if (is.null(observed)) {
    data.frame(population = character(), nqtl = integer(), h2 = numeric(),
               effect = character(), algorithm = character(), pca = logical(),
               n_mappings = integer())
  } else {
    observed[, c("population", "nqtl", "h2", "effect", "algorithm", "pca",
                 "n_mappings"), drop = FALSE]
  }
  # Coerce types to match
  obs$population <- as.character(obs$population)
  obs$nqtl       <- as.integer(obs$nqtl)
  obs$h2         <- as.numeric(obs$h2)
  obs$effect     <- as.character(obs$effect)
  obs$algorithm  <- as.character(obs$algorithm)
  obs$pca        <- as.logical(obs$pca)

  delta <- dplyr::full_join(
    expected_agg, obs,
    by = c("population", "nqtl", "h2", "effect", "algorithm", "pca")
  )
  delta$expected     <- tidyr::replace_na(delta$expected,   0L)
  delta$observed     <- tidyr::replace_na(delta$n_mappings, 0L)
  delta$missing      <- pmax(delta$expected - delta$observed, 0L)
  delta$pct_complete <- ifelse(delta$expected > 0,
                               delta$observed / delta$expected,
                               NA_real_)
  delta <- dplyr::select(delta,
    population, nqtl, h2, effect, algorithm, pca,
    expected, observed, missing, pct_complete
  )
  dplyr::arrange(delta, population, nqtl, h2, effect, algorithm, pca)
}

.vdb_build_traits_delta <- function(groups_df, nqtl, h2, effect, reps,
                                    observed) {
  if (is.null(groups_df) || nrow(groups_df) == 0L) return(NULL)
  expected <- tidyr::crossing(
    population = groups_df$population,
    nqtl       = as.integer(nqtl),
    h2         = as.numeric(h2),
    effect     = as.character(effect)
  )
  expected$expected <- as.integer(reps)

  if (is.null(observed) || nrow(observed) == 0L) {
    obs <- data.frame(population = character(), nqtl = integer(), h2 = numeric(),
                      effect = character(), observed = integer())
  } else {
    obs <- observed %>%
      dplyr::group_by(population, nqtl, h2, effect) %>%
      dplyr::summarise(observed = dplyr::first(n_distinct_traits),
                       .groups = "drop")
    obs$observed <- as.integer(obs$observed)
  }
  obs$population <- as.character(obs$population)
  obs$nqtl       <- as.integer(obs$nqtl)
  obs$h2         <- as.numeric(obs$h2)
  obs$effect     <- as.character(obs$effect)

  delta <- dplyr::full_join(
    expected, obs,
    by = c("population", "nqtl", "h2", "effect")
  )
  delta$expected     <- tidyr::replace_na(delta$expected, 0L)
  delta$observed     <- tidyr::replace_na(delta$observed, 0L)
  delta$missing      <- pmax(delta$expected - delta$observed, 0L)
  delta$pct_complete <- ifelse(delta$expected > 0,
                               delta$observed / delta$expected,
                               NA_real_)
  tibble::as_tibble(delta)
}


# ---------------------------------------------------------------------------
# File-system coverage (marker sets, trait artifacts, mapping partitions)
# ---------------------------------------------------------------------------

.vdb_count_files <- function(dir_path, pattern, recursive = FALSE) {
  if (!dir.exists(dir_path)) return(0L)
  length(list.files(dir_path, pattern = pattern,
                    recursive = recursive, full.names = FALSE))
}

.vdb_marker_set_coverage <- function(base_dir, groups_df, errors_env) {
  if (is.null(groups_df) || nrow(groups_df) == 0L) return(NULL)
  config <- .make_db_config(base_dir)
  markers_dir   <- file.path(config$base_dir, config$markers_dir,
                             config$marker_sets_subdir)
  genotypes_dir <- file.path(config$base_dir, config$markers_dir,
                             config$genotypes_subdir)

  rows <- lapply(seq_len(nrow(groups_df)), function(i) {
    row <- groups_df[i, , drop = FALSE]
    ms_id <- row$marker_set_id
    mf <- if (!is.na(ms_id)) file.path(markers_dir,   paste0(ms_id, "_markers.parquet"))   else NA_character_
    gf <- if (!is.na(ms_id)) file.path(genotypes_dir, paste0(ms_id, "_genotypes.parquet")) else NA_character_
    tibble::tibble(
      population             = row$population,
      maf                    = row$maf,
      marker_set_id          = as.character(ms_id),
      metadata_row_exists    = !is.na(ms_id),
      markers_file_exists    = !is.na(mf) && file.exists(mf),
      genotypes_file_exists  = !is.na(gf) && file.exists(gf),
      n_markers              = if ("n_markers" %in% names(row)) row$n_markers else NA_integer_,
      n_independent_tests    = if ("n_independent_tests" %in% names(row)) row$n_independent_tests else NA_real_,
      strainfile_hash        = if ("strainfile_hash" %in% names(row)) row$strainfile_hash else NA_character_
    )
  })
  dplyr::bind_rows(rows)
}

.vdb_component_counts <- function(base_dir, errors_env) {
  config <- .make_db_config(base_dir)
  bd <- config$base_dir

  mappings_dir <- file.path(bd, config$mappings_dir)
  data_files <- .vdb_safely(
    list.files(mappings_dir, pattern = "^data\\.parquet$",
               recursive = TRUE, full.names = FALSE),
    "list_mapping_data_files",
    errors_env,
    default = character()
  )
  meta_files <- .vdb_safely(
    list.files(mappings_dir, pattern = "^meta\\.parquet$",
               recursive = TRUE, full.names = FALSE),
    "list_mapping_meta_files",
    errors_env,
    default = character()
  )

  list(
    markers_files          = .vdb_count_files(
      file.path(bd, config$markers_dir, config$marker_sets_subdir),
      "_markers\\.parquet$"
    ),
    genotypes_files        = .vdb_count_files(
      file.path(bd, config$markers_dir, config$genotypes_subdir),
      "_genotypes\\.parquet$"
    ),
    marker_set_metadata    = .vdb_count_files(
      file.path(bd, config$marker_set_metadata_dir),
      "_metadata\\.parquet$"
    ),
    traits_files           = .vdb_count_files(
      file.path(bd, config$traits_dir),
      "^[0-9a-f]{20}\\.parquet$"
    ),
    causal_variants_files  = .vdb_count_files(
      file.path(bd, config$traits_dir, config$causal_variants_subdir),
      "_causal\\.parquet$"
    ),
    phenotypes_files       = .vdb_count_files(
      file.path(bd, config$traits_dir, config$phenotypes_subdir),
      "_phenotype\\.parquet$"
    ),
    causal_genotypes_files = .vdb_count_files(
      file.path(bd, config$traits_dir, config$causal_genotypes_subdir),
      "_causal_geno\\.parquet$"
    ),
    mapping_data_files     = length(data_files),
    mapping_meta_files     = length(meta_files),
    mappings_metadata_exists = file.exists(file.path(bd, config$metadata_file))
  )
}


# ---------------------------------------------------------------------------
# Component summary table
# ---------------------------------------------------------------------------

.vdb_build_summary <- function(mode, counts, expected_counts, observed_mapping_rows) {
  components <- c(
    "marker_sets_markers", "marker_sets_genotypes", "marker_set_metadata",
    "traits_metadata_files", "causal_variants_files",
    "phenotypes_files", "causal_genotypes_files",
    "mapping_data_files", "mapping_meta_files", "mappings_metadata_rows"
  )

  observed <- c(
    counts$markers_files,
    counts$genotypes_files,
    counts$marker_set_metadata,
    counts$traits_files,
    counts$causal_variants_files,
    counts$phenotypes_files,
    counts$causal_genotypes_files,
    counts$mapping_data_files,
    counts$mapping_meta_files,
    if (is.null(observed_mapping_rows)) NA_integer_ else as.integer(observed_mapping_rows)
  )

  if (mode == "discovery" || is.null(expected_counts)) {
    expected <- rep(NA_integer_, length(components))
  } else {
    expected <- c(
      expected_counts$n_marker_sets,
      expected_counts$n_marker_sets,
      expected_counts$n_marker_sets,
      expected_counts$n_traits,
      expected_counts$n_traits,
      expected_counts$n_traits,
      expected_counts$n_traits,
      expected_counts$n_mappings,
      expected_counts$n_mappings,
      expected_counts$n_mappings
    )
  }

  missing <- ifelse(is.na(expected), NA_integer_,
                    pmax(as.integer(expected) - as.integer(observed), 0L))
  pct <- ifelse(is.na(expected) | expected == 0, NA_real_,
                pmin(as.integer(observed) / as.integer(expected), 1))

  tibble::tibble(
    component    = components,
    expected     = as.integer(expected),
    observed     = as.integer(observed),
    missing      = as.integer(missing),
    pct_complete = pct
  )
}


# ---------------------------------------------------------------------------
# Top-level verify_db()
# ---------------------------------------------------------------------------

#' Verify completeness of a nemascan-sims-nf Parquet database.
#'
#' @param base_dir  Path to the {outputDir}/db directory.
#' @param strainfile Optional path to the pipeline strainfile TSV.
#' @param nqtl,h2,effect Vectors of expected architecture values. May be NULL.
#' @param reps Integer number of replicates.
#' @param cv_maf,cv_ld Causal-variant pool MAF/LD scalars used at pipeline time.
#' @param algorithm Character vector of mapping algorithms expected.
#' @param pca Logical vector of PCA covariate states expected.
#' @param max_enumeration_rows Guard that disables per-ID enumeration above
#'   this cartesian-product size. Aggregate delta tables still computed.
#' @return A list (see inline docs). `verify_db()` never throws.
verify_db <- function(base_dir,
                      strainfile = NULL,
                      nqtl       = NULL,
                      h2         = NULL,
                      effect     = NULL,
                      reps       = NULL,
                      cv_maf     = NULL,
                      cv_ld      = NULL,
                      algorithm  = c("inbred", "loco"),
                      pca        = c(TRUE, FALSE),
                      max_enumeration_rows = .VDB_DEFAULT_ENUM_LIMIT) {

  errors_env <- .vdb_new_errors()
  run_at     <- Sys.time()

  architecture_args <- list(
    strainfile = strainfile, nqtl = nqtl, h2 = h2, effect = effect,
    reps = reps, cv_maf = cv_maf, cv_ld = cv_ld
  )
  any_supplied <- any(!vapply(architecture_args, is.null, logical(1)))
  all_supplied <- all(!vapply(architecture_args, is.null, logical(1)))
  mode <- if (any_supplied) "expectation" else "discovery"

  if (mode == "expectation" && !all_supplied) {
    missing_args <- names(architecture_args)[
      vapply(architecture_args, is.null, logical(1))]
    .vdb_record(errors_env, "verify_db",
                glue("Expectation mode missing arguments: ",
                     "{paste(missing_args, collapse=', ')}. ",
                     "Falling back to incomplete expectation — per-ID ",
                     "enumeration skipped for missing dimensions."),
                "warning")
  }

  empty_result <- list(
    status              = "error",
    mode                = mode,
    base_dir            = base_dir,
    run_at              = run_at,
    expected_params     = list(groups = NULL, nqtl = nqtl, h2 = h2,
                               effect = effect, reps = reps,
                               cv_maf = cv_maf, cv_ld = cv_ld,
                               algorithm = algorithm, pca = pca),
    observed_mappings   = NULL,
    observed_traits     = NULL,
    mappings_delta      = NULL,
    traits_delta        = NULL,
    missing_traits      = NULL,
    missing_mappings    = NULL,
    summary             = NULL,
    marker_set_coverage = NULL,
    errors              = .vdb_errors_to_tibble(errors_env)
  )

  if (!dir.exists(base_dir)) {
    .vdb_record(errors_env, "verify_db",
                glue("base_dir does not exist: {base_dir}"), "error")
    empty_result$errors <- .vdb_errors_to_tibble(errors_env)
    return(empty_result)
  }

  # ---- Strainfile → hydrated groups (expectation mode only) ---------------
  groups_df <- NULL
  if (!is.null(strainfile)) {
    sf_df   <- .vdb_read_strainfile(strainfile, errors_env)
    groups_df <- .vdb_hydrate_groups(sf_df, base_dir, errors_env)
  }

  # ---- Counts + observed query --------------------------------------------
  counts <- .vdb_safely(
    .vdb_component_counts(base_dir, errors_env),
    "component_counts", errors_env,
    default = list(
      markers_files = 0L, genotypes_files = 0L, marker_set_metadata = 0L,
      traits_files = 0L, causal_variants_files = 0L, phenotypes_files = 0L,
      causal_genotypes_files = 0L, mapping_data_files = 0L,
      mapping_meta_files = 0L, mappings_metadata_exists = FALSE
    )
  )

  con <- .vdb_open_db(base_dir, errors_env)
  on.exit({
    if (!is.null(con)) {
      try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    }
  }, add = TRUE)

  observed         <- .vdb_query_observed(con, errors_env)
  observed_ids     <- .vdb_query_observed_ids(con, errors_env)
  observed_n_rows  <- if (is.null(observed_ids)) NULL else nrow(observed_ids)

  # ---- Discovery mode: short-circuit and return observed only -------------
  if (mode == "discovery") {
    grid <- .vdb_derive_grid_from_db(observed)
    observed_traits <- if (is.null(observed)) NULL else .vdb_observed_traits(observed)
    summary_tbl <- .vdb_build_summary("discovery", counts, NULL, observed_n_rows)
    result <- list(
      status              = "discovery",
      mode                = "discovery",
      base_dir            = base_dir,
      run_at              = run_at,
      expected_params     = list(groups = NULL, nqtl = NULL, h2 = NULL,
                                 effect = NULL, reps = NULL,
                                 cv_maf = NULL, cv_ld = NULL,
                                 algorithm = algorithm, pca = pca),
      observed_mappings   = observed,
      observed_traits     = observed_traits,
      observed_grid       = grid,
      mappings_delta      = NULL,
      traits_delta        = NULL,
      missing_traits      = NULL,
      missing_mappings    = NULL,
      summary             = summary_tbl,
      marker_set_coverage = NULL,
      errors              = .vdb_errors_to_tibble(errors_env)
    )
    return(result)
  }

  # ---- Expectation mode (possibly partial) -------------------------------
  observed_traits <- if (is.null(observed)) NULL else .vdb_observed_traits(observed)

  can_aggregate <- !is.null(groups_df) && nrow(groups_df) > 0L &&
                   !is.null(nqtl) && !is.null(h2) && !is.null(effect) &&
                   !is.null(reps)

  expected_counts <- NULL
  mappings_delta  <- NULL
  traits_delta    <- NULL
  missing_traits  <- NULL
  missing_mappings <- NULL

  if (can_aggregate) {
    g <- nrow(groups_df)
    n_traits   <- as.integer(g * length(nqtl) * length(h2) * length(effect) *
                             as.integer(reps))
    n_mappings <- as.integer(n_traits * length(algorithm) * length(pca))
    expected_counts <- list(
      n_marker_sets = as.integer(g),
      n_traits      = n_traits,
      n_mappings    = n_mappings
    )

    expected_agg <- .vdb_aggregate_expected(groups_df, nqtl, h2, effect, reps,
                                            algorithm, pca)
    mappings_delta <- .vdb_build_mappings_delta(expected_agg, observed)
    traits_delta   <- .vdb_build_traits_delta(groups_df, nqtl, h2, effect, reps,
                                              observed)
  }

  # ---- Per-ID enumeration (expensive; only when all CV args present) -----
  can_enumerate <- can_aggregate && !is.null(cv_maf) && !is.null(cv_ld)
  if (can_enumerate) {
    enum <- .vdb_enumerate_expected(
      groups_df, nqtl, h2, effect, reps, cv_maf, cv_ld,
      algorithm, pca, max_enumeration_rows, errors_env
    )
    if (!is.null(enum)) {
      valid <- !is.na(enum$mapping_id)
      enum_valid <- enum[valid, , drop = FALSE]
      obs_ids <- if (is.null(observed_ids)) {
        data.frame(trait_id = character(), mapping_id = character())
      } else {
        observed_ids
      }

      missing_mappings <- dplyr::anti_join(
        enum_valid,
        dplyr::distinct(obs_ids[, "mapping_id", drop = FALSE]),
        by = "mapping_id"
      ) %>%
        dplyr::select(dplyr::any_of(c("population", "nqtl", "h2", "effect",
                                      "rep", "algorithm", "pca", "mapping_id")))

      missing_traits <- dplyr::anti_join(
        dplyr::distinct(enum_valid,
                        population, nqtl, h2, effect, rep, trait_id),
        dplyr::distinct(obs_ids[, "trait_id", drop = FALSE]),
        by = "trait_id"
      )
    }
  }

  marker_set_coverage <- .vdb_marker_set_coverage(base_dir, groups_df, errors_env)

  summary_tbl <- .vdb_build_summary("expectation", counts, expected_counts,
                                    observed_n_rows)

  status <- .vdb_resolve_status(
    mode      = "expectation",
    errors_df = .vdb_errors_to_tibble(errors_env),
    observed  = observed,
    summary   = summary_tbl,
    can_aggregate = can_aggregate
  )

  list(
    status              = status,
    mode                = "expectation",
    base_dir            = base_dir,
    run_at              = run_at,
    expected_params     = list(
      groups = groups_df, nqtl = nqtl, h2 = h2, effect = effect,
      reps = reps, cv_maf = cv_maf, cv_ld = cv_ld,
      algorithm = algorithm, pca = pca
    ),
    observed_mappings   = observed,
    observed_traits     = observed_traits,
    mappings_delta      = mappings_delta,
    traits_delta        = traits_delta,
    missing_traits      = missing_traits,
    missing_mappings    = missing_mappings,
    summary             = summary_tbl,
    marker_set_coverage = marker_set_coverage,
    errors              = .vdb_errors_to_tibble(errors_env)
  )
}

.vdb_observed_traits <- function(observed) {
  observed %>%
    dplyr::group_by(population, nqtl, h2, effect) %>%
    dplyr::summarise(
      n_traits       = dplyr::first(n_distinct_traits),
      n_distinct_rep = dplyr::first(n_distinct_rep),
      .groups        = "drop"
    )
}

.vdb_resolve_status <- function(mode, errors_df, observed, summary, can_aggregate) {
  n_errors <- sum(errors_df$level == "error", na.rm = TRUE)
  if (is.null(observed) || nrow(observed) == 0L) {
    if (n_errors > 0L) return("error")
  }
  if (mode == "discovery") return("discovery")
  if (!can_aggregate) return("incomplete")
  any_missing <- !is.null(summary) &&
                 any(!is.na(summary$missing) & summary$missing > 0L)
  if (any_missing) return("incomplete")
  if (n_errors > 0L) return("incomplete")
  "complete"
}


# ---------------------------------------------------------------------------
# Plotting helpers — pure functions of `result`, never throw
# ---------------------------------------------------------------------------

.vdb_plot_error <- function(msg) {
  ggplot2::ggplot() +
    ggplot2::labs(title = paste("verify_db plot error:", msg)) +
    ggplot2::theme_void()
}

.vdb_plot_safely <- function(expr, label) {
  tryCatch(force(expr),
           error = function(e) .vdb_plot_error(paste(label, "-", conditionMessage(e))))
}

plot_n_mappings_by_param <- function(result, facet_by = "population") {
  .vdb_plot_safely({
    if (is.null(result$observed_mappings) || nrow(result$observed_mappings) == 0L) {
      return(.vdb_plot_error("no observed mappings"))
    }
    df <- result$observed_mappings %>%
      dplyr::mutate(
        nQTL = factor(nqtl, levels = sort(unique(nqtl))),
        h2   = factor(h2,   levels = sort(unique(h2))),
        ap   = paste(algorithm,
                     ifelse(as.logical(pca), "pca", "nopca"),
                     sep = "/")
      )
    ref_line <- NA_real_
    if (result$mode == "expectation" &&
        !is.null(result$expected_params$reps) &&
        !is.null(result$expected_params$effect)) {
      ref_line <- length(result$expected_params$effect) *
                  as.integer(result$expected_params$reps)
    }
    p <- ggplot2::ggplot(df,
      ggplot2::aes(x = nQTL, y = n_mappings, fill = h2)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9)) +
      ggplot2::facet_grid(
        stats::as.formula(paste("ap ~", facet_by))
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(angle = 90, hjust = 1),
        axis.title.x    = ggplot2::element_text(face = "bold"),
        axis.title.y    = ggplot2::element_text(face = "bold"),
        legend.position = "right"
      ) +
      ggplot2::labs(
        x = "Number of QTL",
        y = "Observed mappings",
        fill = expression(bold(italic(h)^{2})),
        title = "Mappings per (nQTL, h2) cell",
        subtitle = if (is.na(ref_line)) NULL else
                   glue("Dashed line = expected per bar ({ref_line})")
      )
    if (!is.na(ref_line)) {
      p <- p + ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed",
                                   colour = "firebrick")
    }
    p
  }, "plot_n_mappings_by_param")
}

plot_n_traits_by_param <- function(result, facet_by = "population") {
  .vdb_plot_safely({
    obs_t <- result$observed_traits
    if (is.null(obs_t) || nrow(obs_t) == 0L) {
      return(.vdb_plot_error("no observed traits"))
    }
    df <- obs_t %>%
      dplyr::mutate(
        nQTL = factor(nqtl, levels = sort(unique(nqtl))),
        h2   = factor(h2,   levels = sort(unique(h2)))
      )
    ref_line <- NA_real_
    if (result$mode == "expectation" && !is.null(result$expected_params$reps)) {
      ref_line <- as.integer(result$expected_params$reps)
    }
    p <- ggplot2::ggplot(df,
      ggplot2::aes(x = nQTL, y = n_traits, fill = h2)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.9)) +
      ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by))) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
        axis.title.x = ggplot2::element_text(face = "bold"),
        axis.title.y = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(
        x = "Number of QTL",
        y = "Distinct traits",
        fill = expression(bold(italic(h)^{2})),
        title = "Traits per (nQTL, h2) cell"
      )
    if (!is.na(ref_line)) {
      p <- p + ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed",
                                   colour = "firebrick")
    }
    p
  }, "plot_n_traits_by_param")
}

plot_completeness_heatmap <- function(result) {
  .vdb_plot_safely({
    if (result$mode == "expectation") {
      df <- result$mappings_delta
      fill_col <- "pct_complete"
      title <- "Completeness (observed / expected)"
      scale <- ggplot2::scale_fill_viridis_c(limits = c(0, 1),
                                             na.value = "grey80")
    } else {
      df <- result$observed_mappings
      fill_col <- "n_mappings"
      title <- "Observed mapping count (discovery mode)"
      scale <- ggplot2::scale_fill_viridis_c(na.value = "grey80")
    }
    if (is.null(df) || nrow(df) == 0L) {
      return(.vdb_plot_error("no delta/observed data"))
    }
    df <- df %>%
      dplyr::mutate(
        nqtl = factor(nqtl, levels = sort(unique(nqtl))),
        h2   = factor(h2,   levels = sort(unique(h2))),
        ap   = paste(algorithm,
                     ifelse(as.logical(pca), "pca", "nopca"),
                     sep = "/")
      )
    ggplot2::ggplot(df,
      ggplot2::aes(x = nqtl, y = h2, fill = .data[[fill_col]])) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::facet_grid(ap ~ population) +
      scale +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
      ) +
      ggplot2::labs(
        x = "Number of QTL",
        y = expression(bold(italic(h)^{2})),
        fill = fill_col,
        title = title
      )
  }, "plot_completeness_heatmap")
}

plot_marker_set_coverage <- function(result) {
  .vdb_plot_safely({
    ms <- result$marker_set_coverage
    if (is.null(ms) || nrow(ms) == 0L) {
      return(.vdb_plot_error("no marker set coverage"))
    }
    long <- ms %>%
      dplyr::mutate(label = paste(population, maf, sep = "_")) %>%
      dplyr::select(label, metadata_row_exists, markers_file_exists,
                    genotypes_file_exists) %>%
      tidyr::pivot_longer(-label, names_to = "artifact", values_to = "present")
    ggplot2::ggplot(long,
      ggplot2::aes(x = label, y = artifact, fill = present)) +
      ggplot2::geom_tile(colour = "white") +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "#1b7837",
                                            `FALSE` = "#c51b7d"),
                                 na.value = "grey70") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      ) +
      ggplot2::labs(
        x = "Population / MAF",
        y = "Artifact",
        fill = "Present",
        title = "Marker set coverage"
      )
  }, "plot_marker_set_coverage")
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

.vdb_parse_values <- function(raw) {
  if (is.null(raw)) return(NULL)
  if (is.character(raw) && length(raw) == 1L && file.exists(raw)) {
    lines <- readLines(raw, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    return(lines)
  }
  if (is.character(raw) && length(raw) == 1L && grepl(",", raw, fixed = TRUE)) {
    parts <- strsplit(raw, ",", fixed = TRUE)[[1]]
    return(trimws(parts))
  }
  raw
}

.vdb_parse_cli <- function(argv) {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("Package 'optparse' is required for the verify_db CLI")
  }
  option_list <- list(
    optparse::make_option("--base_dir",   type = "character",
      help = "Path to {outputDir}/db (required)"),
    optparse::make_option("--strainfile", type = "character", default = NULL),
    optparse::make_option("--nqtl",       type = "character", default = NULL,
      help = "CSV file or inline comma-separated nqtl values"),
    optparse::make_option("--h2",         type = "character", default = NULL,
      help = "CSV file or inline comma-separated h2 values"),
    optparse::make_option("--effect",     type = "character", default = NULL,
      help = "CSV file or inline comma-separated effect strings"),
    optparse::make_option("--reps",       type = "integer",   default = NULL),
    optparse::make_option("--cv_maf",     type = "double",    default = NULL),
    optparse::make_option("--cv_ld",      type = "double",    default = NULL),
    optparse::make_option("--max_enumeration_rows", type = "integer",
      default = .VDB_DEFAULT_ENUM_LIMIT),
    optparse::make_option("--quiet",      action = "store_true", default = FALSE)
  )
  optparse::parse_args(
    optparse::OptionParser(option_list = option_list),
    args = argv, positional_arguments = FALSE
  )
}

.vdb_print_summary <- function(result) {
  cat("verify_db: ", result$base_dir, "\n", sep = "")
  cat("mode: ", result$mode, "   status: ", result$status, "\n\n", sep = "")

  if (!is.null(result$summary)) {
    fmt <- function(x) ifelse(is.na(x), "—", as.character(x))
    fmt_pct <- function(x) ifelse(is.na(x), "—",
                                  sprintf("%5.1f%%", x * 100))
    lines <- sprintf("  %-28s %10s %10s %10s %10s",
                     result$summary$component,
                     fmt(result$summary$expected),
                     fmt(result$summary$observed),
                     fmt(result$summary$missing),
                     fmt_pct(result$summary$pct_complete))
    cat(sprintf("  %-28s %10s %10s %10s %10s\n",
                "component", "expected", "observed", "missing", "pct"))
    cat(paste(lines, collapse = "\n"), "\n\n", sep = "")
  }

  if (!is.null(result$mappings_delta)) {
    gaps <- result$mappings_delta %>%
      dplyr::filter(.data$missing > 0) %>%
      dplyr::arrange(dplyr::desc(.data$missing))
    if (nrow(gaps) > 0) {
      cat("missing mapping cells (top 10 by count):\n")
      print(utils::head(gaps, 10L))
      cat("\n")
    }
  }

  n_err <- sum(result$errors$level == "error",   na.rm = TRUE)
  n_warn <- sum(result$errors$level == "warning", na.rm = TRUE)
  cat("errors: ", n_err, "   warnings: ", n_warn, "\n", sep = "")
  if (nrow(result$errors) > 0) {
    for (i in seq_len(nrow(result$errors))) {
      cat(sprintf("  [%s] %s: %s\n",
                  result$errors$level[i],
                  result$errors$step[i],
                  result$errors$message[i]))
    }
  }
  invisible(NULL)
}

.vdb_exit_code <- function(status) {
  switch(status,
    complete   = 0L,
    discovery  = 0L,
    incomplete = 1L,
    error      = 2L,
    3L)
}

.vdb_main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opt <- tryCatch(.vdb_parse_cli(argv),
                  error = function(e) {
                    message("verify_db: ", conditionMessage(e))
                    return(NULL)
                  })
  if (is.null(opt) || is.null(opt$base_dir)) {
    message("verify_db: --base_dir is required")
    return(2L)
  }

  result <- verify_db(
    base_dir             = opt$base_dir,
    strainfile           = opt$strainfile,
    nqtl                 = as.integer(.vdb_parse_values(opt$nqtl)),
    h2                   = as.numeric(.vdb_parse_values(opt$h2)),
    effect               = as.character(.vdb_parse_values(opt$effect)),
    reps                 = opt$reps,
    cv_maf               = opt$cv_maf,
    cv_ld                = opt$cv_ld,
    max_enumeration_rows = opt$max_enumeration_rows
  )

  if (!isTRUE(opt$quiet)) .vdb_print_summary(result)
  .vdb_exit_code(result$status)
}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

.vdb_is_cli <- function() {
  cli_args <- commandArgs(trailingOnly = FALSE)
  any(grepl("--file=", cli_args, fixed = TRUE)) &&
    any(grepl("verify_db\\.R$", cli_args))
}

if (.vdb_is_cli() && !interactive()) {
  quit(status = .vdb_main(), save = "no")
}
