# qtl_database.R - QTL database storage functions
#
# Provides functions for storing QTL analysis results in Parquet format.
# QTL results are stored separately from the mapping database to allow
# re-running analysis with different parameters.
#
# Database structure:
#   {qtl_dir}/
#   ├── qtl_regions/
#   │   └── {population}_{algorithm}_qtl_regions.parquet
#   ├── analysis_summary.parquet
#   └── analysis_metadata.parquet

library(arrow)
library(dplyr)
library(glue)

# ==============================================================================
# Database Configuration
# ==============================================================================

.qtl_constants <- list(
  qtl_regions_dir = "qtl_regions",
  summary_file = "analysis_summary.parquet",
  metadata_file = "analysis_metadata.parquet",
  compression = "snappy"
)


# ==============================================================================
# Schema Definitions
# ==============================================================================

#' Define Arrow schema for QTL regions
#'
#' @return Arrow schema object
qtl_regions_schema <- function() {
  arrow::schema(
    mapping_id = arrow::utf8(),
    threshold_method = arrow::utf8(),
    peak_id = arrow::int32(),
    CHROM = arrow::utf8(),
    startPOS = arrow::int32(),
    peakPOS = arrow::int32(),
    endPOS = arrow::int32(),
    interval_size = arrow::int32(),
    n_sig_markers = arrow::int32(),
    max_log10p = arrow::float64(),
    peak_marker = arrow::utf8(),
    sig_threshold_value = arrow::float64(),
    ci_size = arrow::int32(),
    snp_grouping = arrow::int32(),
    algorithm = arrow::utf8(),
    population = arrow::utf8(),
    maf = arrow::float64()
  )
}


#' Define Arrow schema for analysis summary
#'
#' @return Arrow schema object
analysis_summary_schema <- function() {
  arrow::schema(
    mapping_id = arrow::utf8(),
    threshold_method = arrow::utf8(),
    population = arrow::utf8(),
    maf = arrow::float64(),
    algorithm = arrow::utf8(),
    n_markers = arrow::int32(),
    n_significant = arrow::int32(),
    pct_significant = arrow::float64(),
    n_qtl = arrow::int32(),
    max_log10p = arrow::float64(),
    threshold_value = arrow::float64(),
    ci_size = arrow::int32(),
    snp_grouping = arrow::int32()
  )
}


#' Define Arrow schema for analysis metadata
#'
#' @return Arrow schema object
analysis_metadata_schema <- function() {
  arrow::schema(
    run_id = arrow::utf8(),
    alpha = arrow::float64(),
    ci_size = arrow::int32(),
    snp_grouping = arrow::int32(),
    n_populations = arrow::int32(),
    n_mappings = arrow::int32(),
    n_qtl_bf = arrow::int32(),
    n_qtl_eigen = arrow::int32(),
    created_at = arrow::timestamp("us")
  )
}


# ==============================================================================
# Database Initialization
# ==============================================================================

#' Initialize QTL database directory structure
#'
#' @param qtl_dir QTL database root directory
#' @return Invisibly returns the qtl_dir path
init_qtl_database <- function(qtl_dir) {
  qtl_regions_dir <- file.path(qtl_dir, .qtl_constants$qtl_regions_dir)

  for (dir in c(qtl_dir, qtl_regions_dir)) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      log_msg(glue::glue("Created directory: {dir}"))
    }
  }

  invisible(qtl_dir)
}


# ==============================================================================
# Write Functions
# ==============================================================================

#' Get path to QTL regions file for a population/algorithm
#'
#' @param population Population identifier
#' @param algorithm Algorithm (INBRED or LOCO)
#' @param qtl_dir QTL database root directory
#' @return Path to QTL regions Parquet file
get_qtl_regions_path <- function(population, algorithm, qtl_dir) {
  file.path(
    qtl_dir,
    .qtl_constants$qtl_regions_dir,
    glue::glue("{population}_{algorithm}_qtl_regions.parquet")
  )
}


#' Write QTL regions to database
#'
#' Writes QTL regions for a (population, algorithm) batch.
#' Overwrites existing file if present.
#'
#' @param df Dataframe with QTL regions
#' @param population Population identifier
#' @param algorithm Algorithm (INBRED or LOCO)
#' @param qtl_dir QTL database root directory
#' @return Invisibly returns the output path
write_qtl_regions <- function(df, population, algorithm, qtl_dir) {
  init_qtl_database(qtl_dir)
  output_path <- get_qtl_regions_path(population, algorithm, qtl_dir)

  if (nrow(df) == 0) {
    log_msg(glue::glue("No QTL regions to write for {population}_{algorithm}"))
    return(invisible(output_path))
  }

  # Ensure proper column types
  df <- df %>%
    dplyr::mutate(
      mapping_id = as.character(mapping_id),
      threshold_method = as.character(threshold_method),
      peak_id = as.integer(peak_id),
      CHROM = as.character(CHROM),
      startPOS = as.integer(startPOS),
      peakPOS = as.integer(peakPOS),
      endPOS = as.integer(endPOS),
      interval_size = as.integer(interval_size),
      n_sig_markers = as.integer(n_sig_markers),
      max_log10p = as.numeric(max_log10p),
      peak_marker = as.character(peak_marker),
      sig_threshold_value = as.numeric(sig_threshold_value),
      ci_size = as.integer(ci_size),
      snp_grouping = as.integer(snp_grouping),
      algorithm = as.character(algorithm),
      population = as.character(population),
      maf = as.numeric(maf)
    )

  arrow::write_parquet(df, output_path, compression = .qtl_constants$compression)
  log_msg(glue::glue("Wrote {nrow(df)} QTL regions to: {basename(output_path)}"))

  invisible(output_path)
}


#' Write analysis summary for a batch
#'
#' Writes per-mapping summary statistics for a (population, algorithm) batch.
#'
#' @param df Dataframe with summary statistics
#' @param population Population identifier
#' @param algorithm Algorithm (INBRED or LOCO)
#' @param qtl_dir QTL database root directory
#' @return Invisibly returns the output path
write_analysis_summary_batch <- function(df, population, algorithm, qtl_dir) {
  init_qtl_database(qtl_dir)
  output_path <- file.path(
    qtl_dir,
    glue::glue("analysis_summary_{population}_{algorithm}.parquet")
  )

  if (nrow(df) == 0) {
    log_msg(glue::glue("No summary data to write for {population}_{algorithm}"))
    return(invisible(output_path))
  }

  # Ensure proper column types
  df <- df %>%
    dplyr::mutate(
      mapping_id = as.character(mapping_id),
      threshold_method = as.character(threshold_method),
      population = as.character(population),
      maf = as.numeric(maf),
      algorithm = as.character(algorithm),
      n_markers = as.integer(n_markers),
      n_significant = as.integer(n_significant),
      pct_significant = as.numeric(pct_significant),
      n_qtl = as.integer(n_qtl),
      max_log10p = as.numeric(max_log10p),
      threshold_value = as.numeric(threshold_value),
      ci_size = as.integer(ci_size),
      snp_grouping = as.integer(snp_grouping)
    )

  arrow::write_parquet(df, output_path, compression = .qtl_constants$compression)
  log_msg(glue::glue("Wrote {nrow(df)} summary rows to: {basename(output_path)}"))

  invisible(output_path)
}


#' Merge batch summary files into single analysis_summary.parquet
#'
#' @param qtl_dir QTL database root directory
#' @return Invisibly returns the output path
merge_analysis_summaries <- function(qtl_dir) {
  summary_pattern <- "^analysis_summary_.*\\.parquet$"
  summary_files <- list.files(qtl_dir, pattern = summary_pattern, full.names = TRUE)

  if (length(summary_files) == 0) {
    log_msg("No batch summary files found to merge", level = "WARN")
    return(invisible(NULL))
  }

  log_msg(glue::glue("Merging {length(summary_files)} batch summary files..."))

  # Read and combine all batch files
  all_summaries <- lapply(summary_files, arrow::read_parquet)
  combined <- dplyr::bind_rows(all_summaries)

  # Write combined summary
  output_path <- file.path(qtl_dir, .qtl_constants$summary_file)
  arrow::write_parquet(combined, output_path, compression = .qtl_constants$compression)
  log_msg(glue::glue("Wrote combined summary: {nrow(combined)} rows"))

  # Clean up batch files
  file.remove(summary_files)
  log_msg(glue::glue("Removed {length(summary_files)} batch files"))

  invisible(output_path)
}


#' Write analysis metadata
#'
#' Records the parameters used for QTL analysis.
#'
#' @param alpha Significance level
#' @param ci_size CI size parameter
#' @param snp_grouping SNP grouping parameter
#' @param n_populations Number of populations processed
#' @param n_mappings Total number of mappings processed
#' @param n_qtl_bf Total QTL found with BF threshold
#' @param n_qtl_eigen Total QTL found with EIGEN threshold
#' @param qtl_dir QTL database root directory
#' @return Invisibly returns the output path
write_analysis_metadata <- function(alpha, ci_size, snp_grouping,
                                     n_populations, n_mappings,
                                     n_qtl_bf, n_qtl_eigen,
                                     qtl_dir) {
  init_qtl_database(qtl_dir)
  output_path <- file.path(qtl_dir, .qtl_constants$metadata_file)

  # Generate run ID from timestamp
  run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

  metadata <- data.frame(
    run_id = run_id,
    alpha = as.numeric(alpha),
    ci_size = as.integer(ci_size),
    snp_grouping = as.integer(snp_grouping),
    n_populations = as.integer(n_populations),
    n_mappings = as.integer(n_mappings),
    n_qtl_bf = as.integer(n_qtl_bf),
    n_qtl_eigen = as.integer(n_qtl_eigen),
    created_at = Sys.time(),
    stringsAsFactors = FALSE
  )

  arrow::write_parquet(metadata, output_path, compression = .qtl_constants$compression)
  log_msg(glue::glue("Wrote analysis metadata: run_id={run_id}"))

  invisible(output_path)
}


# ==============================================================================
# Read Functions
# ==============================================================================

#' Read QTL regions for a population/algorithm
#'
#' @param population Population identifier
#' @param algorithm Algorithm (INBRED or LOCO)
#' @param qtl_dir QTL database root directory
#' @return Dataframe with QTL regions
read_qtl_regions <- function(population, algorithm, qtl_dir) {
  input_path <- get_qtl_regions_path(population, algorithm, qtl_dir)

  if (!file.exists(input_path)) {
    warning(glue::glue("QTL regions file not found: {input_path}"))
    return(data.frame())
  }

  arrow::read_parquet(input_path) %>%
    as.data.frame()
}


#' Read all QTL regions from database
#'
#' @param qtl_dir QTL database root directory
#' @return Dataframe with all QTL regions
read_all_qtl_regions <- function(qtl_dir) {
  qtl_regions_dir <- file.path(qtl_dir, .qtl_constants$qtl_regions_dir)

  if (!dir.exists(qtl_regions_dir)) {
    warning(glue::glue("QTL regions directory not found: {qtl_regions_dir}"))
    return(data.frame())
  }

  files <- list.files(qtl_regions_dir, pattern = "\\.parquet$", full.names = TRUE)

  if (length(files) == 0) {
    return(data.frame())
  }

  all_regions <- lapply(files, arrow::read_parquet)
  dplyr::bind_rows(all_regions) %>%
    as.data.frame()
}


#' Read analysis summary
#'
#' @param qtl_dir QTL database root directory
#' @return Dataframe with analysis summary
read_analysis_summary <- function(qtl_dir) {
  input_path <- file.path(qtl_dir, .qtl_constants$summary_file)

  if (!file.exists(input_path)) {
    warning(glue::glue("Analysis summary file not found: {input_path}"))
    return(data.frame())
  }

  arrow::read_parquet(input_path) %>%
    as.data.frame()
}


#' Read analysis metadata
#'
#' @param qtl_dir QTL database root directory
#' @return Dataframe with analysis metadata
read_analysis_metadata <- function(qtl_dir) {
  input_path <- file.path(qtl_dir, .qtl_constants$metadata_file)

  if (!file.exists(input_path)) {
    warning(glue::glue("Analysis metadata file not found: {input_path}"))
    return(data.frame())
  }

  arrow::read_parquet(input_path) %>%
    as.data.frame()
}


# ==============================================================================
# Utility Functions
# ==============================================================================
#' Get QTL database statistics
#'
#' @param qtl_dir QTL database root directory
#' @return List with database statistics
qtl_db_stats <- function(qtl_dir) {
  if (!dir.exists(qtl_dir)) {
    return(list(
      exists = FALSE,
      n_qtl_files = 0,
      n_qtl_regions = 0,
      n_summary_rows = 0
    ))
  }

  qtl_regions_dir <- file.path(qtl_dir, .qtl_constants$qtl_regions_dir)
  qtl_files <- if (dir.exists(qtl_regions_dir)) {
    list.files(qtl_regions_dir, pattern = "\\.parquet$", full.names = TRUE)
  } else {
    character()
  }

  n_qtl_regions <- 0
  if (length(qtl_files) > 0) {
    n_qtl_regions <- sum(sapply(qtl_files, function(f) {
      nrow(arrow::read_parquet(f))
    }))
  }

  summary_path <- file.path(qtl_dir, .qtl_constants$summary_file)
  n_summary_rows <- if (file.exists(summary_path)) {
    nrow(arrow::read_parquet(summary_path))
  } else {
    0
  }

  list(
    exists = TRUE,
    n_qtl_files = length(qtl_files),
    n_qtl_regions = n_qtl_regions,
    n_summary_rows = n_summary_rows,
    qtl_files = basename(qtl_files)
  )
}
