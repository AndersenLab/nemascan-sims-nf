#!/usr/bin/env Rscript
# compute_causal_detection.R — Per-causal-variant detection status across all
# mappings and significance thresholds for every simulated trait in a Parquet DB.
#
# Detection is defined as the causal marker's own log10p exceeding the
# significance threshold (direct p-value test, not interval-based). Covers all
# 4 mode×type combinations (inbred/loco × pca/nopca) at BF and EIGEN thresholds.
#
# Output: one row per (causal × mapping × threshold) with detected = TRUE/FALSE.
# Expected size: nrow(var_exp) × 4 mappings × 2 thresholds.
#
# Usage:
#   Rscript scripts/compute_causal_detection.R \
#     --db_path /path/to/db \
#     --output_dir 2026-07-08_mapping-panel-size-tests/data/proc \
#     [--workers 7]


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
              help = "Directory for output TSV [default: cwd]"),
  make_option("--workers",
              type = "integer", default = NULL,
              help = "Parallel workers [default: detectCores() - 1]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$db_path))     stop("--db_path is required")
if (!dir.exists(opt$db_path)) stop("--db_path directory not found: ", opt$db_path)

db_path    <- opt$db_path
output_dir <- opt$output_dir

detected  <- parallel::detectCores(logical = TRUE)
n_workers <- if (!is.null(opt$workers)) opt$workers else max(1L, detected - 1L, na.rm = TRUE)


# ==============================================================================
# 2. Library loading
# ==============================================================================

script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[startsWith(script_args, "--file=")])
r_lib_dir   <- normalizePath(file.path(dirname(script_path), "..", "R"),
                             mustWork = TRUE)

suppressPackageStartupMessages({
  Sys.setenv(NEMASCAN_R_LIB = r_lib_dir)
  source(file.path(r_lib_dir, "setup.R"))  # utils, io, database, queries, analysis,
                                           # qtl_database, sim_performance
  library(furrr)
  library(readr)
})


# ==============================================================================
# 3. Per-trait worker
# ==============================================================================

#' Query causal detection for one trait across all mappings and thresholds
#'
#' Calls `query_causal_pvalues()` inside a tryCatch so a single bad trait
#' does not abort the parallel batch. Returns NULL on any error.
#'
#' @param trait_id 20-char hex trait ID
#' @param population Strain-set / panel identifier the trait belongs to
#' @param base_dir Database root directory
#' @param meta Pre-loaded mapping metadata (from `get_metadata()`), passed
#'   through so the full metadata table is not re-read from disk per trait.
#' @return Long data frame (causal × mapping × threshold rows), or NULL on error
detection_for_trait <- function(trait_id, population, base_dir, meta = NULL) {
  tryCatch({
    query_causal_pvalues(trait_id, base_dir = base_dir, meta = meta)
  }, error = function(e) {
    message(sprintf("  [skip] trait %s (pop %s): %s", trait_id, population, conditionMessage(e)))
    NULL
  })
}


# ==============================================================================
# 4. Enumerate traits
# ==============================================================================

cat("Loading metadata from:", db_path, "\n")
meta <- get_metadata(db_path)

traits <- meta %>%
  dplyr::distinct(trait_id, population)

cat(sprintf("  %d traits across %d populations\n",
            nrow(traits), length(unique(traits$population))))


# ==============================================================================
# 5. Output filename
# ==============================================================================

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

run_tag  <- format(Sys.time(), "%Y%m%d_%H%M%S")
db_tag   <- basename(db_path)
out_file <- file.path(output_dir,
                      paste(run_tag, db_tag, "causal_detection.tsv", sep = "_"))


# ==============================================================================
# 6. Parallel computation
# ==============================================================================

plan(multisession, workers = n_workers)
cat(sprintf("Computing causal detection for %d traits with %d workers...\n",
            nrow(traits), n_workers))

# `meta` (loaded once above) is shipped to each worker as a global so the full
# mappings_metadata.parquet is read a single time per worker process instead of
# once per trait — see query_causal_pvalues()'s `meta` argument.
worker_globals <- list(
  r_lib_dir          = r_lib_dir,
  db_path            = db_path,
  meta               = meta,
  detection_for_trait = detection_for_trait
)

results <- furrr::future_map2(
  traits$trait_id,
  traits$population,
  function(trait_id, population) {
    if (!exists("query_causal_pvalues", mode = "function")) {
      suppressPackageStartupMessages({
        Sys.setenv(NEMASCAN_R_LIB = r_lib_dir)
        source(file.path(r_lib_dir, "setup.R"))
      })
    }
    detection_for_trait(trait_id, population, base_dir = db_path, meta = meta)
  },
  .options = furrr::furrr_options(seed = TRUE, globals = worker_globals)
)

plan(sequential)


# ==============================================================================
# 7. Bind and write
# ==============================================================================

cat("Combining results...\n")
detection <- dplyr::bind_rows(results)

cat(sprintf("  %d rows across %d traits\n",
            nrow(detection), dplyr::n_distinct(detection$trait_id)))

readr::write_tsv(detection, out_file)

cat("Done. Output written:\n")
cat(" ", out_file, "\n")
