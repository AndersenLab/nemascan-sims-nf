#!/usr/bin/env Rscript
# assess_panel_performance.R — Assess GWA mapping panel performance across
# populations, mapping methods, and significance thresholds from a Parquet
# simulation database.
#
# Reads all mappings from the database, applies BF and EIGEN significance
# thresholds, scores each against its known causal variants, and writes two
# output TSVs: per-replicate Power/FDR and per-panel mean Power/FDR.
#
# Usage:
#   Rscript scripts/assess_panel_performance.R \
#     --db_path /path/to/db \
#     --output_dir 2026-07-08_mapping-panel-size-tests/data/proc \
#     [--ci_size 150] [--snp_grouping 1000] [--alpha 0.05] [--workers 7]


# ==============================================================================
# 1. CLI arguments
# ==============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--db_path",
              type = "character", default = NULL,
              help = "Path to Parquet simulation database [required]"),
  make_option("--output_dir",
              type = "character", default = ".",
              help = "Directory for output TSVs [default: cwd]"),
  make_option("--ci_size",
              type = "integer", default = 150L,
              help = "Markers flanking peak for CI definition [default: 150]"),
  make_option("--snp_grouping",
              type = "integer", default = 1000L,
              help = "Marker index distance for QTL grouping [default: 1000]"),
  make_option("--alpha",
              type = "double", default = 0.05,
              help = "Significance level [default: 0.05]"),
  make_option("--workers",
              type = "integer", default = NULL,
              help = "Parallel workers [default: detectCores() - 1]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$db_path))        stop("--db_path is required")
if (!dir.exists(opt$db_path))    stop("--db_path directory not found: ", opt$db_path)

db_path      <- opt$db_path
output_dir   <- opt$output_dir
ci_size      <- opt$ci_size
snp_grouping <- opt$snp_grouping
alpha        <- opt$alpha

detected   <- parallel::detectCores(logical = TRUE)
n_workers  <- if (!is.null(opt$workers)) opt$workers else max(1L, detected - 1L, na.rm = TRUE)


# ==============================================================================
# 2. Library loading
# ==============================================================================

# Resolve R library dir from the script's own location.
# scripts/ and R/ are siblings under the project root, so R/ is always ../R
# relative to this file. commandArgs() --file= is set by Rscript at launch.
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[startsWith(script_args, "--file=")])
r_lib_dir   <- normalizePath(file.path(dirname(script_path), "..", "R"),
                              mustWork = TRUE)

suppressPackageStartupMessages({
  Sys.setenv(NEMASCAN_R_LIB = r_lib_dir)
  source(file.path(r_lib_dir, "setup.R"))      # utils, io, database, queries, analysis,
                                               # qtl_database, sim_performance
  source(file.path(r_lib_dir, "assessment.R")) # compile_full_assessment()
  library(memoise)
  library(furrr)
  library(readr)
})


# ==============================================================================
# 3. Memoised readers
#
# Each worker receives a serialized copy of these closures (empty cache).
# Within a worker's population batch, repeated calls with the same arguments
# are served from the in-process cache rather than re-reading Parquet files.
# ==============================================================================

marker_meta_c <- memoise(read_marker_set_metadata)  # keyed on (population, maf, base_dir)
causal_c      <- memoise(read_causal_variants_data) # keyed on (trait_id, base_dir)


# ==============================================================================
# 4. Per-mapping helper functions
# ==============================================================================

#' Resolve metadata row and marker set metadata for one mapping
#'
#' @param mapping_id 20-char hex mapping ID
#' @param meta Metadata dataframe from get_metadata()
#' @param base_dir Database root directory
#' @return Named list: row (single metadata row), population, maf, ms_meta
get_mapping_context <- function(mapping_id, meta, base_dir) {
  m <- dplyr::filter(meta, .data$mapping_id == .env$mapping_id)
  if (nrow(m) != 1) {
    stop(sprintf("Expected one metadata row for %s, got %d", mapping_id, nrow(m)))
  }
  pop     <- m$population
  maf     <- as.numeric(m$maf)
  ms_meta <- marker_meta_c(pop, maf, base_dir)
  if (is.null(ms_meta)) {
    stop(sprintf("No marker set metadata found for population=%s maf=%s", pop, maf))
  }
  list(row = m, population = pop, maf = maf, ms_meta = ms_meta)
}

#' Compute the significance threshold for one method
#'
#' Returns NULL when EIGEN data is absent so the caller can skip that threshold
#' rather than error.
#'
#' @param threshold_method Character: "BF" or "EIGEN"
#' @param mapping_data Dataframe of mapping rows; row count drives the BF denominator
#' @param ms_meta Named list from read_marker_set_metadata(); must contain n_independent_tests
#' @param alpha Significance level
#' @return list(threshold_value, threshold_method) or NULL
resolve_threshold <- function(threshold_method, mapping_data, ms_meta, alpha) {
  if (threshold_method == "BF") {
    calculate_threshold("BF", n_markers = nrow(mapping_data), alpha = alpha)
  } else if (threshold_method == "EIGEN") {
    n_ind <- as.numeric(ms_meta$n_independent_tests)
    if (is.na(n_ind)) return(NULL)
    calculate_threshold("EIGEN", n_independent = n_ind, alpha = alpha)
  } else {
    stop(sprintf("Unsupported threshold_method: %s", threshold_method))
  }
}

#' Assemble the mapping_params list consumed by compile_full_assessment()
#'
#' @param ctx Context list from get_mapping_context()
#' @param threshold_method Character: "BF" or "EIGEN"
#' @param ci_size Marker CI size
#' @param snp_grouping Marker grouping distance
#' @param alpha Significance level
#' @return Named list of simulation and assessment parameters
build_mapping_params <- function(ctx, threshold_method, ci_size, snp_grouping, alpha) {
  m <- ctx$row
  list(
    population       = ctx$population,
    maf              = ctx$maf,
    nqtl             = m$nqtl,
    effect           = m$effect,
    rep              = m$rep,
    h2               = m$h2,
    algorithm        = m$algorithm,
    pca              = m$pca,
    threshold_method = threshold_method,
    mode             = m$algorithm,
    type             = if (isTRUE(m$pca)) "pca" else "nopca",
    alpha            = alpha,
    ci_size          = ci_size,
    snp_grouping     = snp_grouping
  )
}

#' Assess one mapping across all significance thresholds
#'
#' Reads the mapping partition once, then loops over thresholds. For each:
#' applies the threshold, extracts QTL regions, and scores against causal
#' variants via compile_full_assessment().
#'
#' @param mapping_id 20-char hex mapping ID
#' @param meta Metadata dataframe (rows for this population)
#' @param base_dir Database root directory
#' @param thresholds Character vector of threshold methods
#' @param ci_size Marker CI size
#' @param snp_grouping Marker grouping distance
#' @param alpha Significance level
#' @return Dataframe of per-QTL assessment rows for all thresholds, or NULL
assess_mapping <- function(mapping_id, meta, base_dir,
                           thresholds   = c("BF", "EIGEN"),
                           ci_size      = 150L,
                           snp_grouping = 1000L,
                           alpha        = 0.05) {
  ctx          <- get_mapping_context(mapping_id, meta, base_dir)
  mapping_data <- query_mapping_direct(mapping_id, ctx$population, ctx$ms_meta, base_dir)
  if (nrow(mapping_data) == 0) return(NULL)
  causal_vars  <- causal_c(ctx$row$trait_id, base_dir)

  purrr::map(thresholds, function(tm) {
    thr <- resolve_threshold(tm, mapping_data, ctx$ms_meta, alpha)
    if (is.null(thr)) return(NULL)
    processed   <- analyze_mapping(mapping_data,
                                   thr$threshold_value, thr$threshold_method,
                                   ci_size      = ci_size,
                                   snp_grouping = snp_grouping,
                                   verbose      = FALSE)
    qtl_regions <- extract_qtl_regions(processed)
    params      <- build_mapping_params(ctx, tm, ci_size, snp_grouping, alpha)
    compile_full_assessment(processed, qtl_regions, causal_vars, params)
  }) |>
    dplyr::bind_rows()
}


# ==============================================================================
# 5. Metadata and panel sizes
# ==============================================================================

cat("Loading metadata from:", db_path, "\n")
meta <- get_metadata(db_path)
cat(sprintf("  %d mappings across %d populations\n",
            nrow(meta), length(unique(meta$population))))

panel_size <- get_all_marker_set_metadata(base_dir = db_path) |>
  dplyr::transmute(
    strain_set_id = population,
    n_strains     = lengths(strsplit(strain_list, ","))
  ) |>
  dplyr::distinct()


# ==============================================================================
# 6. Output filename prefix
# ==============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

date_tag    <- format(Sys.Date(), "%Y%m%d")
db_tag      <- basename(db_path)
params_tag  <- sprintf("ci%d_snp%d_a%s", ci_size, snp_grouping, alpha)
file_prefix <- paste(date_tag, db_tag, params_tag, sep = "_")


# ==============================================================================
# 7. Temp directory for incremental worker output
# ==============================================================================

tmp_dir <- tempfile("assess_panels_")
dir.create(tmp_dir)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)


# ==============================================================================
# 8. Parallel assessment
# ==============================================================================

plan(multisession, workers = n_workers)
cat(sprintf("Assessing %d populations with %d workers...\n",
            length(unique(meta$population)), n_workers))

pop_groups <- dplyr::group_split(meta, population)

# Explicit globals: helper functions, memoised readers, and scalar params.
# Workers source the R library in their init block to satisfy the underlying
# function dependencies (read_marker_set_metadata, compile_full_assessment, etc.).
worker_globals <- list(
  r_lib_dir            = r_lib_dir,
  db_path              = db_path,
  ci_size              = ci_size,
  snp_grouping         = snp_grouping,
  alpha                = alpha,
  tmp_dir              = tmp_dir,
  get_mapping_context  = get_mapping_context,
  resolve_threshold    = resolve_threshold,
  build_mapping_params = build_mapping_params,
  assess_mapping       = assess_mapping,
  marker_meta_c        = marker_meta_c,
  causal_c             = causal_c
)

result_paths <- furrr::future_map(
  pop_groups,
  function(pop_meta) {
    # Source R library to make underlying functions available in this process.
    suppressPackageStartupMessages({
      Sys.setenv(NEMASCAN_R_LIB = r_lib_dir)
      source(file.path(r_lib_dir, "setup.R"))
      source(file.path(r_lib_dir, "assessment.R"))
    })

    rows <- purrr::map(
      pop_meta$mapping_id,
      ~ assess_mapping(.x,
                       meta         = pop_meta,
                       base_dir     = db_path,
                       thresholds   = c("BF", "EIGEN"),
                       ci_size      = ci_size,
                       snp_grouping = snp_grouping,
                       alpha        = alpha)
    ) |>
      dplyr::bind_rows()

    out_file <- file.path(tmp_dir, paste0(pop_meta$population[[1]], ".tsv"))
    readr::write_tsv(rows, out_file)
    out_file
  },
  .options = furrr::furrr_options(seed = TRUE, globals = worker_globals)
)

plan(sequential)


# ==============================================================================
# 9. Bind and clean up
# ==============================================================================

cat("Combining population results...\n")
assess <- purrr::map(result_paths, readr::read_tsv, show_col_types = FALSE) |>
  dplyr::bind_rows()
cat(sprintf("  %d assessment rows\n", nrow(assess)))


# ==============================================================================
# 10. per_rep aggregation
# ==============================================================================

cat("Computing per-replicate Power/FDR...\n")
per_rep <- assess |>
  dplyr::group_by(strain_set_id, mode, type, threshold,
                  nQTL, h2, effect_distribution, maf, simREP) |>
  designate_qtl() |>
  dplyr::group_modify(~ count_outcomes(.x)) |>
  calculate_simrep_performance(stringent = FALSE)


# ==============================================================================
# 11. per_panel aggregation
# ==============================================================================

cat("Computing per-panel mean Power/FDR...\n")
per_panel <- per_rep |>
  dplyr::group_by(strain_set_id, mode, type, threshold) |>
  dplyr::summarise(
    mean_power = mean(Power, na.rm = TRUE),
    mean_fdr   = mean(FDR,   na.rm = TRUE),
    n_reps     = dplyr::n(),
    .groups    = "drop"
  ) |>
  dplyr::left_join(panel_size, by = "strain_set_id")


# ==============================================================================
# 12. Write outputs
# ==============================================================================

per_rep_file   <- file.path(output_dir, paste0(file_prefix, "_per_rep_assessment.tsv"))
per_panel_file <- file.path(output_dir, paste0(file_prefix, "_per_panel_assessment.tsv"))

readr::write_tsv(per_rep,   per_rep_file)
readr::write_tsv(per_panel, per_panel_file)

cat("Done. Outputs written:\n")
cat(" ", per_rep_file,   "\n")
cat(" ", per_panel_file, "\n")
