# io.R - File I/O functions for database migration
# Reading mapping files, raw GWA output, BIM files, and EIGEN file handling

library(readr)
library(yaml)
library(glue)
library(data.table)


#' Read raw GWA output file with auto-detection
#'
#' Auto-detects format from file extension (.fastGWA vs .mlma) and normalizes
#' columns to the database schema.
#'
#' @param file_path Path to GWA output file
#' @param verbose Whether to log messages (default TRUE)
#' @return Dataframe with columns: marker, CHROM, POS, A1, A2, AF1, BETA, SE, P
read_raw_gwa_file <- function(file_path, verbose = TRUE) {
  ext <- tools::file_ext(file_path)
  df <- data.table::fread(file_path)

  if (ext == "fastGWA") {
    df <- df %>%
      dplyr::rename(CHROM = CHR, marker = SNP) %>%
      dplyr::select(marker, CHROM, POS, A1, A2, AF1, BETA, SE, P)
  } else if (ext == "mlma") {
    df <- df %>%
      dplyr::rename(CHROM = Chr, marker = SNP, POS = bp,
                     AF1 = Freq, BETA = b, SE = se, P = p) %>%
      dplyr::select(marker, CHROM, POS, A1, A2, AF1, BETA, SE, P)
  } else {
    stop(paste("Unknown GWA file format:", ext))
  }

  # Ensure CHROM is character type
  df$CHROM <- as.character(df$CHROM)

  if (verbose) {
    log_msg(glue::glue("Read {nrow(df)} markers from {basename(file_path)} ({ext} format)"))
  }

  # No log10p computation — derived at analysis time via safe_log10p()
  # No N column — removed from schema (not in .mlma files, unused downstream)
  df
}


#' Read PLINK .bim file for marker set creation
#'
#' Reads 6-column no-header .bim format and creates marker identifiers.
#' AF1 is not included — allele frequency is not available in .bim files
#' and is not part of the markers schema (stored per-mapping in mappings only).
#'
#' @param bim_path Path to .bim file
#' @return Dataframe with columns: marker, CHROM, POS, A1, A2
read_bim_file <- function(bim_path) {
  df <- data.table::fread(bim_path, header = FALSE,
    col.names = c("CHROM", "SNP", "cM", "POS", "A1", "A2"))
  df %>%
    dplyr::mutate(
      CHROM = as.character(CHROM),
      marker = paste(CHROM, POS, sep = ":")
    ) %>%
    dplyr::select(marker, CHROM, POS, A1, A2)
}


#' Load YAML configuration file
#'
#' @param config_path Path to YAML config file
#' @return List containing configuration values
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(glue::glue("Config file not found: {config_path}"))
  }

  config <- yaml::read_yaml(config_path)

  # Validate required sections exist
  required_sections <- c("significance", "intervals")
  missing <- setdiff(required_sections, names(config))
  if (length(missing) > 0) {
    stop(glue::glue("Config missing required sections: {paste(missing, collapse = ', ')}"))
  }

  config
}


#' Read mapping file with validation
#'
#' Reads a *_mapping.tsv file and validates required columns are present.
#' Uses data.table::fread for fast reading.
#'
#' @param file_path Path to mapping TSV file
#' @param required_cols Columns that must be present
#' @param verbose Whether to log messages (default TRUE, set FALSE for batch operations)
#' @return Dataframe with mapping data
read_mapping_file <- function(file_path,
                              required_cols = c("marker", "CHROM", "POS"),
                              verbose = TRUE) {
  if (!file.exists(file_path)) {
    stop(glue::glue("Mapping file not found: {file_path}"))
  }

  if (verbose) {
    log_msg(glue::glue("Reading mapping file: {basename(file_path)}"))
  }

  df <- data.table::fread(
    file_path,
    sep = "\t",
    header = TRUE,
    colClasses = c(CHROM = "character", marker = "character", POS = "integer"),
    showProgress = FALSE
  )

  # Convert to data.frame for compatibility with existing code
  df <- as.data.frame(df)

  if (verbose) {
    log_msg(glue::glue("Loaded {nrow(df)} markers from mapping file"))
  }

  # Validate required columns
  validate_mapping_columns(df, required_cols)

  df
}


#' Read EIGEN independent tests file
#'
#' @param file_path Path to *_total_independent_tests.txt file
#' @return Numeric value of independent tests
read_eigen_file <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(glue::glue("EIGEN file not found: {file_path}"))
  }

  log_msg(glue::glue("Reading EIGEN file: {basename(file_path)}"))

  # File contains a single numeric value
  value <- readr::read_lines(file_path, n_max = 1)
  independent_tests <- as.numeric(value)

  if (is.na(independent_tests)) {
    stop(glue::glue("Could not parse EIGEN value from file: {file_path}"))
  }

  log_msg(glue::glue("EIGEN independent tests: {independent_tests}"))

  independent_tests
}


#' Find EIGEN file matching population and MAF
#'
#' Searches for the appropriate EIGEN file based on population and MAF
#' extracted from the mapping filename.
#'
#' @param mapping_filename Mapping filename to match
#' @param search_dir Directory to search for EIGEN files
#' @return Path to matching EIGEN file, or NULL if not found
find_eigen_file <- function(mapping_filename, search_dir) {
  pop_maf <- extract_population_maf(mapping_filename)

  if (is.null(pop_maf)) {
    warning("Could not extract population/MAF from filename")
    return(NULL)
  }

  expected_name <- build_eigen_filename(pop_maf$population, pop_maf$maf)
  expected_path <- file.path(search_dir, expected_name)

  if (file.exists(expected_path)) {
    log_msg(glue::glue("Found EIGEN file: {expected_name}"))
    return(expected_path)
  }

  # Search recursively if not found directly
  found <- list.files(
    search_dir,
    pattern = glob2rx(expected_name),
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(found) > 0) {
    log_msg(glue::glue("Found EIGEN file: {found[1]}"))
    return(found[1])
  }

  warning(glue::glue("EIGEN file not found: {expected_name}"))
  NULL
}


#' Find all EIGEN independent tests files in a directory
#'
#' @param search_dir Directory to search
#' @param recursive Search subdirectories (default TRUE)
#' @return Character vector of EIGEN file paths
find_eigen_files <- function(search_dir, recursive = TRUE) {
  if (!dir.exists(search_dir)) {
    return(character())
  }

  files <- list.files(
    search_dir,
    pattern = "_total_independent_tests\\.txt$",
    recursive = recursive,
    full.names = TRUE
  )

  log_msg(glue::glue("Found {length(files)} EIGEN files in {search_dir}"))

  files
}


#' Parse EIGEN filename to extract population and MAF
#'
#' @param filename EIGEN filename (basename or full path)
#' @return Named list with population and maf, or NULL if parsing fails
parse_eigen_filename <- function(filename) {
  # Remove path if present
  filename <- basename(filename)

  # Expected pattern: {population}_{maf}_total_independent_tests.txt
  pattern <- "^(.+)_([0-9.]+)_total_independent_tests\\.txt$"

  matches <- regmatches(filename, regexec(pattern, filename))[[1]]

  if (length(matches) == 0) {
    warning(glue::glue("Could not parse EIGEN filename: {filename}"))
    return(NULL)
  }

  list(
    population = matches[2],
    maf = as.numeric(matches[3]),
    filename = filename
  )
}


#' Build EIGEN lookup table from directory
#'
#' Scans directory for EIGEN files and builds a lookup table keyed by
#' population and MAF. This should be called once at the start of migration
#' to enable marker sets to be created with their EIGEN values.
#'
#' @param search_dir Directory to search for EIGEN files
#' @return Named list keyed by "{population}_{maf}" with n_independent_tests values
build_eigen_lookup <- function(search_dir) {
  eigen_files <- find_eigen_files(search_dir)

  if (length(eigen_files) == 0) {
    log_msg("No EIGEN files found in search directory", level = "WARN")
    return(list())
  }

  lookup <- list()

  for (eigen_path in eigen_files) {
    parsed <- parse_eigen_filename(eigen_path)

    if (is.null(parsed)) {
      log_msg(glue::glue("Skipping unparseable EIGEN file: {basename(eigen_path)}"),
        level = "WARN"
      )
      next
    }

    tryCatch(
      {
        n_tests <- read_eigen_file(eigen_path)
        key <- paste0(parsed$population, "_", parsed$maf)
        lookup[[key]] <- list(
          n_independent_tests = n_tests,
          source_file = basename(eigen_path)
        )
      },
      error = function(e) {
        log_msg(glue::glue("Error reading EIGEN file {basename(eigen_path)}: {e$message}"),
          level = "WARN"
        )
      }
    )
  }

  log_msg(glue::glue("Built EIGEN lookup with {length(lookup)} marker sets"))
  lookup
}


#' Get EIGEN value from lookup table
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param eigen_lookup Lookup table from build_eigen_lookup()
#' @return Named list with n_independent_tests and source_file, or NULL if not found
get_eigen_from_lookup <- function(population, maf, eigen_lookup) {
  key <- paste0(population, "_", maf)
  eigen_lookup[[key]]
}
