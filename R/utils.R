# utils.R - Utility functions for mapping reprocessing
# Filename parsing, validation helpers, and logging

#' Parse mapping filename to extract simulation parameters
#'
#' Extracts nqtl, rep, h2, maf, effect, population, algorithm, and pca status
#' from standard NemaScan mapping filenames.
#'
#' @param filename Character string of the mapping filename (basename, not full path)
#' @return Named list with parsed parameters, or NULL if parsing fails
#'
#' @examples
#' parse_mapping_filename("5_1_0.2_0.05_gamma_ct.fullpop.20210901_processed_LMM-EXACT-INBRED_PCA_mapping.tsv")
parse_mapping_filename <- function(filename) {
  # Remove path if present
  filename <- basename(filename)

  # Normalize filename: remove _reprocessed suffix and _qtl_regions variant
  # This allows parsing both original and reprocessed files
  normalized <- filename
  normalized <- sub("_qtl_regions_reprocessed\\.tsv$", "_mapping.tsv", normalized)
  normalized <- sub("_mapping_reprocessed\\.tsv$", "_mapping.tsv", normalized)
  normalized <- sub("_reprocessed\\.tsv$", ".tsv", normalized)

  # Expected pattern:
  # {nqtl}_{rep}_{h2}_{maf}_{effect}_{population}_processed_LMM-EXACT-{INBRED/LOCO}_{PCA}_mapping.tsv
  # or without PCA:
  # {nqtl}_{rep}_{h2}_{maf}_{effect}_{population}_processed_LMM-EXACT-{INBRED/LOCO}_mapping.tsv

  # Check for PCA in filename
  has_pca <- grepl("_PCA_mapping\\.tsv$", normalized)

  # Build regex pattern
  if (has_pca) {
    pattern <- "^(\\d+)_(\\d+)_([0-9.]+)_([0-9.]+)_([a-zA-Z]+)_(.+)_processed_LMM-EXACT-(INBRED|LOCO)_PCA_mapping\\.tsv$"
  } else {
    pattern <- "^(\\d+)_(\\d+)_([0-9.]+)_([0-9.]+)_([a-zA-Z]+)_(.+)_processed_LMM-EXACT-(INBRED|LOCO)_mapping\\.tsv$"
  }

  matches <- regmatches(normalized, regexec(pattern, normalized))[[1]]

  if (length(matches) == 0) {
    warning(glue::glue("Could not parse filename: {filename}"))
    return(NULL)
  }

  list(
    nqtl = as.integer(matches[2]),
    rep = as.integer(matches[3]),
    h2 = as.numeric(matches[4]),
    maf = as.numeric(matches[5]),
    effect = matches[6],
    population = matches[7],
    algorithm = matches[8],
    pca = has_pca,
    filename = filename
  )
}


#' Extract population and MAF from filename for EIGEN file matching
#'
#' @param filename Mapping filename
#' @return Named list with population and maf
extract_population_maf <- function(filename) {
  parsed <- parse_mapping_filename(filename)
  if (is.null(parsed)) {
    return(NULL)
  }
  list(
    population = parsed$population,
    maf = parsed$maf
  )
}


#' Build expected EIGEN filename from population and MAF
#'
#' @param population Population identifier
#' @param maf Minor allele frequency
#' @return Expected EIGEN filename pattern
build_eigen_filename <- function(population, maf) {
  glue::glue("{population}_{maf}_total_independent_tests.txt")
}


#' Validate mapping dataframe has required columns
#'
#' @param df Dataframe to validate
#' @param required_cols Character vector of required column names
#' @return TRUE if valid, stops with error otherwise
validate_mapping_columns <- function(df, required_cols = c("marker", "CHROM", "POS")) {
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(glue::glue(
      "Missing required columns in mapping file: {paste(missing_cols, collapse = ', ')}"
    ))
  }

  TRUE
}


#' Log message with timestamp
#'
#' @param msg Message to log
#' @param level Log level: "INFO", "WARN", "ERROR", "DEBUG"
#' @param verbose If FALSE, only show WARN and ERROR
log_msg <- function(msg, level = "INFO", verbose = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  if (!verbose && level %in% c("INFO", "DEBUG")) {
    return(invisible(NULL))
  }

  prefix <- switch(level,
    "INFO" = "[INFO]",
    "WARN" = "[WARN]",
    "ERROR" = "[ERROR]",
    "DEBUG" = "[DEBUG]",
    "[INFO]"
  )

  message(glue::glue("{timestamp} {prefix} {msg}"))
}


#' Create output filename from input filename
#'
#' @param input_filename Input mapping filename
#' @param suffix Suffix to add before extension (default: "_reprocessed")
#' @return Output filename
create_output_filename <- function(input_filename, suffix = "_reprocessed") {
  # Remove .tsv extension, add suffix, re-add extension
  base <- sub("\\.tsv$", "", basename(input_filename))
  glue::glue("{base}{suffix}.tsv")
}


#' Get unique trait name from mapping dataframe
#'
#' @param df Mapping dataframe
#' @return Trait name string, or NA if not found
get_trait_name <- function(df) {
  if ("trait" %in% names(df)) {
    traits <- unique(df$trait)
    if (length(traits) == 1) {
      return(traits[1])
    } else if (length(traits) > 1) {
      warning("Multiple traits found in mapping file, using first")
      return(traits[1])
    }
  }
  return(NA_character_)
}
