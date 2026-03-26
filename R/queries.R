# queries.R - Database query functions for simulation mapping data
#
# This file contains all read-only query operations for the mapping database.
# Write operations remain in database.R.
#
# Dependencies: Requires database.R to be sourced first for:
#   - .make_db_config(), .db_constants
#   - get_db_connection(), get_threshold_params()
#   - read_marker_set(), get_marker_set_metadata_dir()
#   - get_all_marker_set_metadata(), list_marker_sets()
#   - log_msg() from utils.R

library(arrow)
library(duckdb)
library(dplyr)
library(glue)


# ==============================================================================
# Database Connection
# ==============================================================================

#' Open DuckDB connection to mapping database
#'
#' Creates a DuckDB connection with markers and mappings views registered.
#' Uses partitioned storage: mappings/population={pop}/mapping_id={id}/data.parquet
#'
#' @param base_dir Database root directory
#' @param read_only Open in read-only mode (default TRUE)
#' @return DuckDB connection object
open_mapping_db <- function(base_dir = "data/db", read_only = TRUE) {
  config <- .make_db_config(base_dir)
  if (!dir.exists(config$base_dir)) {
    stop(glue::glue("Database directory not found: {config$base_dir}"))
  }

  con <- DBI::dbConnect(duckdb::duckdb(), ":memory:", read_only = read_only)

  # Register markers view
  markers_dir  <- file.path(config$base_dir, config$markers_dir, config$marker_sets_subdir)
  marker_files <- list.files(markers_dir, pattern = "_markers\\.parquet$", full.names = TRUE)

  if (length(marker_files) > 0) {
    file_list <- paste0("'", marker_files, "'", collapse = ", ")
    DBI::dbExecute(con, glue::glue("
      CREATE VIEW markers AS
      SELECT * FROM read_parquet([{file_list}])
    "))
  }

  # Register mappings view from partitioned structure
  mappings_dir <- file.path(config$base_dir, config$mappings_dir)
  n_sources <- 0

  if (dir.exists(mappings_dir)) {
    partition_files <- list.files(mappings_dir, pattern = "data\\.parquet$",
                                   recursive = TRUE, full.names = TRUE)
    n_sources <- length(partition_files)

    if (n_sources > 0) {
      glob_pattern <- file.path(mappings_dir, "**", "data.parquet")
      DBI::dbExecute(con, glue::glue("
        CREATE VIEW mappings AS
        SELECT * FROM read_parquet('{glob_pattern}', hive_partitioning = true, union_by_name = true)
      "))
    }
  }

  # Register metadata view
  metadata_file <- file.path(config$base_dir, config$metadata_file)
  if (file.exists(metadata_file)) {
    DBI::dbExecute(con, glue::glue("
      CREATE VIEW metadata AS
      SELECT * FROM read_parquet('{metadata_file}')
    "))
  }

  # Register marker set metadata view (per-population files)
  ms_metadata_dir <- get_marker_set_metadata_dir(base_dir)
  ms_metadata_glob <- file.path(ms_metadata_dir, "*_metadata.parquet")
  ms_metadata_files <- Sys.glob(ms_metadata_glob)
  if (length(ms_metadata_files) > 0) {
    DBI::dbExecute(con, glue::glue("
      CREATE VIEW marker_set_metadata AS
      SELECT * FROM read_parquet('{ms_metadata_glob}', union_by_name = true)
    "))
  }

  log_msg(glue::glue("Opened database with {length(marker_files)} marker sets, {n_sources} mapping partitions"))

  con
}


# ==============================================================================
# Basic Query Functions
# ==============================================================================

#' Query marker set
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @return Dataframe with marker data
query_markers <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  read_marker_set(population, maf, species, vcf_release_id, ms_ld, base_dir)
}


#' Query mappings by population
#'
#' Works with both partitioned and non-partitioned storage structures.
#'
#' @param population Population identifier
#' @param base_dir Database root directory
#' @return Dataframe with mapping data
query_by_population <- function(population, base_dir = "data/db") {
  con <- open_mapping_db(base_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  result <- DBI::dbGetQuery(con, glue::glue("
    SELECT * FROM mappings
    WHERE population = '{population}'
  "))

  if (nrow(result) == 0) {
    stop(glue::glue("No mapping data found for population: {population}"))
  }

  result
}


#' Query specific mapping by ID
#'
#' @param mapping_id Unique mapping identifier
#' @param base_dir Database root directory
#' @return Dataframe with mapping data
query_by_mapping_id <- function(mapping_id, base_dir = "data/db") {
  con <- open_mapping_db(base_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  result <- DBI::dbGetQuery(con, glue::glue("
    SELECT * FROM mappings
    WHERE mapping_id = '{mapping_id}'
  "))

  if (nrow(result) == 0) {
    warning(glue::glue("No data found for mapping_id: {mapping_id}"))
  }

  result
}


#' Get mapping metadata summary
#'
#' @param base_dir Database root directory
#' @return Dataframe with all mapping metadata
get_metadata <- function(base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  metadata_file <- file.path(config$base_dir, config$metadata_file)

  if (!file.exists(metadata_file)) {
    stop("Metadata file not found. Have you added any mappings to the database?")
  }

  arrow::read_parquet(metadata_file) %>%
    as.data.frame()
}


#' List available populations in database
#'
#' Returns populations from partitioned storage structure.
#'
#' @param base_dir Database root directory
#' @return Character vector of population names
list_populations <- function(base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  mappings_dir <- file.path(config$base_dir, config$mappings_dir)

  if (!dir.exists(mappings_dir)) {
    return(character())
  }

  # Get populations from partition directories (population=X)
  partition_dirs <- list.dirs(mappings_dir, recursive = FALSE)
  partition_pops <- grep("^population=", basename(partition_dirs), value = TRUE)

  sub("^population=", "", partition_pops)
}


# ==============================================================================
# Standardized Query Functions for Downstream Processing
# ==============================================================================

#' Query full mapping data with markers joined
#'
#' Returns complete mapping data with marker information joined.
#' Primary function for downstream analysis in Qmd scripts.
#'
#' @param population Optional population filter
#' @param maf Optional MAF filter
#' @param h2 Optional heritability filter
#' @param nqtl Optional nqtl filter
#' @param algorithm Optional algorithm filter (INBRED or LOCO)
#' @param base_dir Database root directory
#' @return Dataframe with complete mapping data including marker info
query_mapping_data <- function(population = NULL, maf = NULL, h2 = NULL,
                               nqtl = NULL, algorithm = NULL,
                               base_dir = "data/db") {
  con <- open_mapping_db(base_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  conditions <- character()

  if (!is.null(population)) {
    conditions <- c(conditions, glue::glue("m.population = '{population}'"))
  }
  if (!is.null(maf)) {
    conditions <- c(conditions, glue::glue("mm.maf = {maf}"))
  }
  if (!is.null(h2)) {
    conditions <- c(conditions, glue::glue("mm.h2 = {h2}"))
  }
  if (!is.null(nqtl)) {
    conditions <- c(conditions, glue::glue("mm.nqtl = {nqtl}"))
  }
  if (!is.null(algorithm)) {
    conditions <- c(conditions, glue::glue("mm.algorithm = '{algorithm}'"))
  }

  where_clause <- if (length(conditions) > 0) {
    paste("WHERE", paste(conditions, collapse = " AND "))
  } else {
    ""
  }

  query <- glue::glue("
    SELECT
      m.*,
      mk.CHROM,
      mk.POS,
      mk.A1,
      mk.A2
    FROM mappings m
    LEFT JOIN metadata mm ON m.mapping_id = mm.mapping_id
    LEFT JOIN markers mk
      ON mm.marker_set_id = mk.marker_set_id
      AND m.marker = mk.marker
    {where_clause}
    ORDER BY m.mapping_id, mk.CHROM, mk.POS
  ")

  DBI::dbGetQuery(con, query)
}


#' Query simulation run summary
#'
#' Returns aggregated summary statistics for a simulation run.
#'
#' @param population Optional population filter
#' @param base_dir Database root directory
#' @return Dataframe with simulation run summary
query_simulation_summary <- function(population = NULL, base_dir = "data/db") {
  con <- open_mapping_db(base_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  where_clause <- if (!is.null(population)) {
    glue::glue("WHERE population = '{population}'")
  } else {
    ""
  }

  DBI::dbGetQuery(con, glue::glue("
    SELECT
      population,
      maf,
      nqtl,
      h2,
      algorithm,
      pca,
      COUNT(DISTINCT rep) as n_replicates,
      COUNT(DISTINCT mapping_id) as n_mappings
    FROM metadata
    {where_clause}
    GROUP BY population, maf, nqtl, h2, algorithm, pca
    ORDER BY population, maf, nqtl, h2, algorithm, pca
  "))
}


#' Query markers for significance threshold analysis
#'
#' Returns mapping statistics needed to apply significance thresholds.
#' Note: log10p is NOT stored; compute via safe_log10p(P) after query.
#'
#' @param mapping_id Mapping ID to query
#' @param base_dir Database root directory
#' @param con Optional existing DuckDB connection (for batch operations)
#' @return Dataframe with marker statistics for threshold analysis
query_for_threshold_analysis <- function(mapping_id, base_dir = "data/db", con = NULL) {
  own_con <- is.null(con)
  if (own_con) {
    con <- get_db_connection(base_dir, use_cache = TRUE)
  }

  # Check if var.exp column exists in mappings (absent in inline-path databases)
  mapping_cols <- DBI::dbGetQuery(con, "SELECT column_name FROM information_schema.columns WHERE table_name = 'mappings'")
  has_var_exp <- "var.exp" %in% mapping_cols$column_name

  var_exp_col <- if (has_var_exp) 'm."var.exp"' else 'NULL AS "var.exp"'

  result <- DBI::dbGetQuery(con, glue::glue("
    SELECT
      m.marker,
      m.mapping_id,
      m.P,
      m.BETA,
      m.SE,
      {var_exp_col},
      mk.CHROM,
      mk.POS,
      m.AF1,
      m.population,
      mm.maf
    FROM mappings m
    LEFT JOIN metadata mm ON m.mapping_id = mm.mapping_id
    LEFT JOIN markers mk
      ON mm.marker_set_id = mk.marker_set_id
      AND m.marker = mk.marker
    WHERE m.mapping_id = '{mapping_id}'
    ORDER BY mk.CHROM, mk.POS, m.marker
  "))

  # Verify CHROM:POS uniqueness (LOCO files may have duplicates if not deduped at write time)
  if (nrow(result) > 0) {
    dup_check <- result %>% dplyr::group_by(CHROM, POS) %>% dplyr::filter(dplyr::n() > 1)
    if (nrow(dup_check) > 0) {
      n_dups <- nrow(dup_check)
      example <- paste0(dup_check$CHROM[1], ":", dup_check$POS[1])
      stop(glue::glue(
        "CHROM:POS uniqueness violation in mapping {mapping_id}: ",
        "{n_dups} duplicate rows (e.g. {example}). ",
        "Check for LOCO deduplication in write_gwa_to_db.R"
      ))
    }
  }

  result
}


#' Bulk query multiple mappings for threshold analysis
#'
#' Efficiently queries multiple mappings in a single database call.
#' Much faster than calling query_for_threshold_analysis() in a loop.
#'
#' @param mapping_ids Character vector of mapping IDs to query
#' @param base_dir Database root directory
#' @param con Optional existing DuckDB connection
#' @return Dataframe with marker statistics for all mappings
query_bulk_for_threshold_analysis <- function(mapping_ids, base_dir = "data/db", con = NULL) {
  if (length(mapping_ids) == 0) {
    return(data.frame())
  }

  own_con <- is.null(con)
  if (own_con) {
    con <- get_db_connection(base_dir, use_cache = TRUE)
  }

  # Check if var.exp column exists in mappings (absent in inline-path databases)
  mapping_cols <- DBI::dbGetQuery(con, "SELECT column_name FROM information_schema.columns WHERE table_name = 'mappings'")
  has_var_exp <- "var.exp" %in% mapping_cols$column_name

  var_exp_col <- if (has_var_exp) 'm."var.exp"' else 'NULL AS "var.exp"'

  # Build IN clause with proper quoting
  ids_quoted <- paste0("'", mapping_ids, "'", collapse = ", ")

  result <- DBI::dbGetQuery(con, glue::glue("
    SELECT
      m.marker,
      m.mapping_id,
      m.P,
      m.BETA,
      m.SE,
      {var_exp_col},
      mk.CHROM,
      mk.POS,
      m.AF1,
      m.population,
      mm.maf
    FROM mappings m
    LEFT JOIN metadata mm ON m.mapping_id = mm.mapping_id
    LEFT JOIN markers mk
      ON mm.marker_set_id = mk.marker_set_id
      AND m.marker = mk.marker
    WHERE m.mapping_id IN ({ids_quoted})
    ORDER BY m.mapping_id, mk.CHROM, mk.POS, m.marker
  "))

  # Verify CHROM:POS uniqueness per mapping_id
  if (nrow(result) > 0) {
    dup_check <- result %>%
      dplyr::group_by(mapping_id, CHROM, POS) %>%
      dplyr::filter(dplyr::n() > 1)
    if (nrow(dup_check) > 0) {
      n_dups <- nrow(dup_check)
      example_id <- dup_check$mapping_id[1]
      example_pos <- paste0(dup_check$CHROM[1], ":", dup_check$POS[1])
      stop(glue::glue(
        "CHROM:POS uniqueness violation in bulk query: ",
        "{n_dups} duplicate rows (e.g. {example_id} at {example_pos}). ",
        "Check for LOCO deduplication in write_gwa_to_db.R"
      ))
    }
  }

  result
}


#' Get unique marker count for a mapping
#'
#' Returns the number of unique markers, used for Bonferroni correction.
#'
#' @param mapping_id Mapping ID to query
#' @param base_dir Database root directory
#' @return Integer count of unique markers
get_marker_count <- function(mapping_id, base_dir = "data/db") {
  con <- open_mapping_db(base_dir)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  result <- DBI::dbGetQuery(con, glue::glue("
    SELECT COUNT(DISTINCT marker) as n
    FROM mappings
    WHERE mapping_id = '{mapping_id}'
  "))

  result$n[1]
}


# ==============================================================================
# Enhanced Query Functions
# ==============================================================================

#' Query mappings by method (algorithm and PCA)
#'
#' Returns all mappings for a population using a specific mapping method.
#'
#' @param population Population identifier
#' @param algorithm Mapping algorithm ("INBRED" or "LOCO")
#' @param pca Optional PCA filter (TRUE/FALSE/NULL for any)
#' @param base_dir Database root directory
#' @return Dataframe with mapping metadata matching the method
query_by_method <- function(population, algorithm, pca = NULL, base_dir = "data/db") {
  metadata <- get_metadata(base_dir)

  result <- metadata %>%
    dplyr::filter(
      population == !!population,
      algorithm == !!algorithm
    )

  if (!is.null(pca)) {
    result <- result %>%
      dplyr::filter(pca == !!pca)
  }

  if (nrow(result) == 0) {
    pca_str <- if (is.null(pca)) "" else paste0(", pca=", pca)
    log_msg(glue::glue(
      "No mappings found for {population} with {algorithm}{pca_str}"
    ), level = "WARN")
  } else {
    log_msg(glue::glue(
      "Found {nrow(result)} mappings for {population} with {algorithm}"
    ))
  }

  result
}


#' Query all mappings across all populations
#'
#' Returns metadata for all mappings in the database.
#'
#' @param base_dir Database root directory
#' @return Dataframe with all mapping metadata
query_all_populations <- function(base_dir = "data/db") {
  metadata <- get_metadata(base_dir)

  n_pops <- length(unique(metadata$population))
  log_msg(glue::glue(
    "Retrieved {nrow(metadata)} mappings across {n_pops} populations"
  ))

  metadata
}


#' Query mappings by genetic architecture
#'
#' Returns mappings matching specified heritability and/or number of QTL.
#'
#' @param h2 Optional heritability filter
#' @param nqtl Optional number of QTL filter
#' @param population Optional population filter
#' @param maf Optional MAF filter
#' @param base_dir Database root directory
#' @return Dataframe with matching mapping metadata
query_by_architecture <- function(h2 = NULL, nqtl = NULL, population = NULL,
                                  maf = NULL, base_dir = "data/db") {
  metadata <- get_metadata(base_dir)

  result <- metadata

  if (!is.null(population)) {
    result <- result %>% dplyr::filter(population == !!population)
  }
  if (!is.null(maf)) {
    result <- result %>% dplyr::filter(maf == !!maf)
  }
  if (!is.null(h2)) {
    result <- result %>% dplyr::filter(h2 == !!h2)
  }
  if (!is.null(nqtl)) {
    result <- result %>% dplyr::filter(nqtl == !!nqtl)
  }

  # Build filter description for logging
  filters <- character()
  if (!is.null(h2)) filters <- c(filters, paste0("h2=", h2))
  if (!is.null(nqtl)) filters <- c(filters, paste0("nqtl=", nqtl))
  if (!is.null(population)) filters <- c(filters, paste0("pop=", population))
  if (!is.null(maf)) filters <- c(filters, paste0("maf=", maf))

  filter_str <- if (length(filters) > 0) paste(filters, collapse = ", ") else "none"

  log_msg(glue::glue(
    "Query by architecture ({filter_str}): {nrow(result)} mappings"
  ))

  result
}


#' Get architecture summary across database
#'
#' Returns an overview of available h2/nQTL combinations.
#'
#' @param base_dir Database root directory
#' @return Dataframe summarizing genetic architectures
query_architecture_summary <- function(base_dir = "data/db") {
  metadata <- get_metadata(base_dir)

  summary <- metadata %>%
    dplyr::group_by(population, maf, h2, nqtl) %>%
    dplyr::summarise(
      n_algorithms = dplyr::n_distinct(algorithm),
      n_replicates = dplyr::n_distinct(rep),
      n_mappings = dplyr::n(),
      algorithms = paste(unique(algorithm), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::arrange(population, maf, h2, nqtl)

  log_msg(glue::glue(
    "Architecture summary: {nrow(summary)} unique combinations"
  ))

  summary
}


#' Query for analysis with threshold parameters
#'
#' Enhanced query that returns mapping data with threshold parameters included.
#' Combines query_for_threshold_analysis with get_threshold_params.
#'
#' @param mapping_id Mapping ID to query
#' @param alpha Significance level for threshold calculation (default: 0.05)
#' @param base_dir Database root directory
#' @return List with mapping_data dataframe and threshold_params
query_for_analysis <- function(mapping_id, alpha = 0.05, base_dir = "data/db") {
  # Get mapping data
  mapping_data <- query_for_threshold_analysis(mapping_id, base_dir)

  if (nrow(mapping_data) == 0) {
    stop(glue::glue("No data found for mapping_id: {mapping_id}"))
  }

  # Extract population and maf from data
  population <- unique(mapping_data$population)[1]
  maf <- unique(mapping_data$maf)[1]

  # Get threshold parameters
  threshold_params <- get_threshold_params(population, maf, alpha, base_dir)

  list(
    mapping_data = mapping_data,
    threshold_params = threshold_params
  )
}


# ==============================================================================
# Database Statistics
# ==============================================================================

#' Get database statistics
#'
#' Returns statistics for partitioned storage structure.
#'
#' @param base_dir Database root directory
#' @return List with database statistics
db_stats <- function(base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  if (!dir.exists(config$base_dir)) {
    return(list(
      exists = FALSE,
      n_marker_sets = 0,
      n_marker_sets_with_eigen = 0,
      n_populations = 0,
      n_mappings = 0,
      n_partitions = 0,
      total_size_mb = 0
    ))
  }

  markers_dir <- file.path(config$base_dir, config$markers_dir, config$marker_sets_subdir)
  marker_files <- if (dir.exists(markers_dir)) {
    list.files(markers_dir, pattern = "_markers\\.parquet$", full.names = TRUE)
  } else {
    character()
  }

  mappings_dir <- file.path(config$base_dir, config$mappings_dir)
  partition_files <- if (dir.exists(mappings_dir)) {
    list.files(mappings_dir, pattern = "data\\.parquet$", recursive = TRUE, full.names = TRUE)
  } else {
    character()
  }

  all_files <- c(marker_files, partition_files)
  metadata_file <- file.path(config$base_dir, config$metadata_file)
  if (file.exists(metadata_file)) {
    all_files <- c(all_files, metadata_file)
  }

  # Include marker set metadata files in size calculation
  ms_metadata_dir <- get_marker_set_metadata_dir(base_dir)
  if (dir.exists(ms_metadata_dir)) {
    ms_meta_files <- list.files(ms_metadata_dir, pattern = "_metadata\\.parquet$", full.names = TRUE)
    all_files <- c(all_files, ms_meta_files)
  }

  total_bytes <- sum(file.info(all_files)$size, na.rm = TRUE)

  n_mappings <- 0
  if (file.exists(metadata_file)) {
    n_mappings <- nrow(arrow::read_parquet(metadata_file))
  }

  # Count marker sets with EIGEN data
  n_with_eigen <- 0
  marker_set_metadata <- get_all_marker_set_metadata(base_dir)
  if (nrow(marker_set_metadata) > 0) {
    n_with_eigen <- sum(!is.na(marker_set_metadata$n_independent_tests))
  }

  populations <- list_populations(base_dir)

  list(
    exists = TRUE,
    n_marker_sets = length(marker_files),
    n_marker_sets_with_eigen = n_with_eigen,
    n_populations = length(populations),
    n_mappings = n_mappings,
    n_partitions = length(partition_files),
    total_size_mb = round(total_bytes / 1024^2, 2),
    marker_sets = list_marker_sets(base_dir),
    populations = populations
  )
}
