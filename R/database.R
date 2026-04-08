# database.R - Parquet/DuckDB storage for simulation mapping data
#
# Database structure reflects simulation data hierarchy:
#   1. Mapping population - strain sets used for GWAS
#   2. Marker set - variants selected by MAF threshold (population + MAF)
#   3. Simulated trait - trait architecture (nQTL, h2, effect, rep)
#   4. Mapping - association statistics (BETA, SE, P) per marker
#
# Storage organization:
#   - markers/{population}_{maf}_markers.parquet - static marker data per marker set
#   - mappings/population={pop}/mapping_id={id}/data.parquet - Hive-partitioned mappings
#   - mappings_metadata.parquet - summary table for quick queries
#   - marker_set_metadata.parquet - EIGEN values, marker counts
#
# Significance thresholds and QTL intervals are computed at query time
# via safe_log10p(P) in analysis.R, not stored in the database.

library(arrow)
library(duckdb)
library(dplyr)
library(glue)

# ==============================================================================
# Database Configuration
# ==============================================================================

# Module-level constants for database structure
# These are fixed and never change - only base_dir is configurable
.db_constants <- list(
  markers_dir            = "markers",
  marker_sets_subdir     = "marker_sets",
  genotypes_subdir       = "genotypes",
  traits_dir             = "traits",
  causal_variants_subdir   = "causal_variants",
  causal_genotypes_subdir  = "causal_genotypes",
  phenotypes_subdir        = "phenotypes",
  mappings_dir           = "mappings",
  markers_pattern        = "{marker_set_id}_markers.parquet",
  genotypes_pattern      = "{marker_set_id}_genotypes.parquet",
  mappings_pattern       = "{population}_mappings.parquet",
  traits_pattern         = "{trait_id}.parquet",
  causal_variants_pattern  = "{trait_id}_causal.parquet",
  causal_genotypes_pattern = "{trait_id}_causal_geno.parquet",
  phenotypes_pattern       = "{trait_id}_phenotype.parquet",
  metadata_file          = "mappings_metadata.parquet",
  marker_set_metadata_dir  = "marker_set_metadata",
  compression            = "snappy"
)

# Default database location
.default_base_dir <- "data/db"

#' Construct full database config from base_dir
#'
#' @param base_dir Database root directory
#' @return List with all config fields
.make_db_config <- function(base_dir) {
  c(list(base_dir = base_dir), .db_constants)
}


# ==============================================================================
# Connection Caching
# ==============================================================================

# Module-level cache for database connections and metadata
.db_cache <- new.env(hash = TRUE, parent = emptyenv())


#' Get a cached database connection
#'
#' Returns a cached DuckDB connection, creating one if needed.
#' Use this for batch operations to avoid repeated connection overhead.
#'
#' @param base_dir Database root directory
#' @param use_cache If TRUE (default), use cached connection; if FALSE, create new
#' @return DuckDB connection object
get_db_connection <- function(base_dir = "data/db", use_cache = TRUE) {
  config <- .make_db_config(base_dir)
  cache_key <- config$base_dir

  # Check if cached connection exists and is valid
  if (use_cache && exists(cache_key, envir = .db_cache)) {
    cached <- get(cache_key, envir = .db_cache)
    # Test if connection is still valid
    tryCatch({
      DBI::dbGetQuery(cached$con, "SELECT 1")
      return(cached$con)
    }, error = function(e) {
      # Connection is stale, remove it
      rm(list = cache_key, envir = .db_cache)
    })
  }

  # Create new connection
  con <- open_mapping_db(base_dir, read_only = TRUE)

  if (use_cache) {
    assign(cache_key, list(con = con, config = config), envir = .db_cache)
  }

  con
}


#' Close cached database connection
#'
#' Closes and removes the cached connection for a database.
#'
#' @param base_dir Database root directory (optional, closes all if NULL)
#' @return Invisibly returns TRUE
close_cached_connection <- function(base_dir = NULL) {
  if (is.null(base_dir)) {
    # Close all cached connections
    for (key in ls(envir = .db_cache)) {
      cached <- get(key, envir = .db_cache)
      tryCatch(
        DBI::dbDisconnect(cached$con, shutdown = TRUE),
        error = function(e) NULL
      )
    }
    rm(list = ls(envir = .db_cache), envir = .db_cache)
  } else {
    cache_key <- base_dir
    if (exists(cache_key, envir = .db_cache)) {
      cached <- get(cache_key, envir = .db_cache)
      tryCatch(
        DBI::dbDisconnect(cached$con, shutdown = TRUE),
        error = function(e) NULL
      )
      rm(list = cache_key, envir = .db_cache)
    }
  }

  invisible(TRUE)
}


#' Execute expression with database connection
#'
#' Provides a connection and ensures cleanup on exit.
#' Useful for one-off queries where caching isn't needed.
#'
#' @param base_dir Database root directory
#' @param expr Expression to evaluate with connection available as `con`
#' @return Result of expression evaluation
with_db_connection <- function(base_dir = "data/db", expr) {
  config <- .make_db_config(base_dir)
  con <- open_mapping_db(base_dir, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
  eval(substitute(expr), envir = list(con = con), enclos = parent.frame())
}


# ==============================================================================
# Cached Metadata Lookups
# ==============================================================================

#' Get threshold parameters for a marker set
#'
#' Returns all threshold-related parameters in a single lookup.
#' Results are cached for batch operations.
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param alpha Significance level (default: 0.05)
#' @param base_dir Database root directory
#' @return Named list with n_markers, n_independent_tests, bf_threshold, eigen_threshold
get_threshold_params <- function(population, maf, alpha = 0.05, base_dir = "data/db") {
  # Get marker set metadata (includes n_independent_tests from EIGEN)
  ms_meta <- read_marker_set_metadata(population, maf, base_dir)

  if (is.null(ms_meta)) {
    stop(glue::glue("Marker set metadata not found for {population}_{maf}"))
  }

  n_markers <- ms_meta$n_markers
  n_independent <- ms_meta$n_independent_tests

  # NOTE: bf_threshold here uses the marker-set n_markers (from the PLINK .bim file).
  # This may differ from the actual GWA output row count when GCTA internally excludes
  # near-singular markers. Callers that need concordance with the legacy BF computation
  # (which uses sum(log10p > 0) from the GWA output) must override n_markers with
  # nrow(mapping_data) — see analyze_qtl.R and assess_sims.R Step 3.
  # Calculate thresholds
  bf_threshold <- -log10(alpha / n_markers)
  eigen_threshold <- if (!is.na(n_independent) && n_independent > 0) {
    -log10(alpha / n_independent)
  } else {
    NA_real_
  }

  list(
    population = population,
    maf = maf,
    n_markers = n_markers,
    n_independent_tests = n_independent,
    bf_threshold = bf_threshold,
    eigen_threshold = eigen_threshold,
    alpha = alpha
  )
}


#' Get thresholds for all MAF values in a population
#'
#' Returns threshold parameters for all marker sets in a population.
#'
#' @param population Population identifier
#' @param alpha Significance level (default: 0.05)
#' @param base_dir Database root directory
#' @return Dataframe with threshold info for each MAF value
get_population_thresholds <- function(population, alpha = 0.05, base_dir = "data/db") {
  all_meta <- get_all_marker_set_metadata(base_dir)

  if (nrow(all_meta) == 0) {
    return(data.frame())
  }

  pop_meta <- all_meta %>%
    dplyr::filter(population == !!population) %>%
    dplyr::mutate(
      bf_threshold = -log10(alpha / n_markers),
      eigen_threshold = dplyr::if_else(
        !is.na(n_independent_tests) & n_independent_tests > 0,
        -log10(alpha / n_independent_tests),
        NA_real_
      ),
      alpha = alpha
    )

  pop_meta
}


# ==============================================================================
# Schema Definitions
# ==============================================================================

#' Define Arrow schema for marker data
#'
#' Marker-specific data that is static across all mappings for a marker set.
#' Marker sets are defined by population + MAF threshold.
#'
#' @return Arrow schema object
markers_schema <- function() {
  arrow::schema(
    marker_set_id = arrow::utf8(),
    marker        = arrow::utf8(),
    CHROM         = arrow::utf8(),
    POS           = arrow::int32(),
    A1            = arrow::utf8(),
    A2            = arrow::utf8()
  )
}


#' Define Arrow schema for mapping data
#'
#' Mapping-specific data that varies for each GWAS mapping.
#' Note: N and log10p removed from schema — log10p is computed at analysis
#' time via safe_log10p(P); N is not in .mlma files and unused downstream.
#' AF1 stored per-mapping from raw GWA output (see AF1 Dual Storage design note).
#'
#' @return Arrow schema object
mappings_schema <- function() {
  arrow::schema(
    marker        = arrow::utf8(),
    AF1           = arrow::float64(),
    BETA          = arrow::float64(),
    SE            = arrow::float64(),
    P             = arrow::float64()
  )
}


#' Define schema for mapping metadata table
#'
#' @return Arrow schema object
metadata_schema <- function() {
  arrow::schema(
    mapping_id          = arrow::utf8(),
    mapping_hash_string = arrow::utf8(),
    trait_id            = arrow::utf8(),
    marker_set_id       = arrow::utf8(),
    hash_schema_version = arrow::utf8(),
    population          = arrow::utf8(),
    maf                 = arrow::float64(),
    nqtl                = arrow::int32(),
    rep                 = arrow::int32(),
    h2                  = arrow::float64(),
    effect              = arrow::utf8(),
    algorithm           = arrow::utf8(),
    pca                 = arrow::boolean(),
    n_markers           = arrow::int32(),
    source_file         = arrow::utf8(),
    processed_at        = arrow::timestamp("us"),
    processing_version  = arrow::utf8()
  )
}


#' Define schema for marker set metadata table
#'
#' Stores per-marker-set metadata including EIGEN independent tests values.
#' Marker sets are uniquely identified by population + MAF.
#'
#' @return Arrow schema object
marker_set_metadata_schema <- function() {
  arrow::schema(
    marker_set_id          = arrow::utf8(),
    marker_set_hash_string = arrow::utf8(),
    hash_schema_version    = arrow::utf8(),  # "v=2": v=1 hashes incompatible — regenerate DB
    population             = arrow::utf8(),
    maf                    = arrow::float64(),
    species                = arrow::utf8(),
    vcf_release_id         = arrow::utf8(),
    ms_ld                  = arrow::float64(),
    n_markers              = arrow::int32(),
    n_independent_tests    = arrow::float64(),
    eigen_source_file      = arrow::utf8(),
    strainfile_hash        = arrow::utf8(),
    strain_list            = arrow::utf8(),
    created_at             = arrow::timestamp("us")
  )
}


#' Define Arrow schema for genotype matrix data (long format)
#'
#' Genotype matrix stored in long format with fixed 4-column schema.
#' Allele encoding: -1 (hom ref), 1 (hom alt), NA (het/missing).
#'
#' @return Arrow schema object
genotype_matrix_schema <- function() {
  arrow::schema(
    marker_set_id = arrow::utf8(),
    CHROM         = arrow::utf8(),
    POS           = arrow::int32(),
    strain        = arrow::utf8(),
    allele        = arrow::float64()
  )
}


#' Define Arrow schema for trait metadata
#'
#' @return Arrow schema object
trait_metadata_schema <- function() {
  arrow::schema(
    trait_id          = arrow::utf8(),
    trait_hash_string = arrow::utf8(),
    marker_set_id     = arrow::utf8(),
    nqtl              = arrow::int32(),
    rep               = arrow::int32(),
    sim_seed          = arrow::int32(),
    h2                = arrow::float64(),
    maf               = arrow::float64(),
    effect            = arrow::utf8(),
    population        = arrow::utf8(),
    cv_maf_effective  = arrow::float64(),
    cv_ld             = arrow::float64(),
    created_at        = arrow::utf8()
  )
}


#' Define Arrow schema for phenotype data
#'
#' @return Arrow schema object
phenotype_schema <- function() {
  arrow::schema(
    strain    = arrow::utf8(),
    phenotype = arrow::float64()
  )
}


#' Define Arrow schema for per-trait causal genotype data
#'
#' Stores per-strain genotypes at causal variant positions, keyed by trait_id.
#' Enables var.exp computation for non-marker causal variants (positions absent
#' from the marker set genotype matrix) via inner_join on CHROM/POS.
#'
#' Allele encoding: -1.0 (hom ref / A1), 1.0 (hom alt / A2), NA (het or missing).
#' Matches the encoding convention of the marker set genotype matrix.
#'
#' @return Arrow schema object
causal_genotypes_schema <- function() {
  arrow::schema(
    trait_id = arrow::utf8(),
    QTL      = arrow::utf8(),   # "CHROM:POS" string — matches .par file format
    CHROM    = arrow::utf8(),
    POS      = arrow::int32(),
    strain   = arrow::utf8(),
    allele   = arrow::float64() # -1.0 (hom ref), 1.0 (hom alt), NA (het/missing)
  )
}


# ==============================================================================
# Database Initialization
# ==============================================================================

#' Initialize database directory structure
#'
#' @param base_dir Database root directory
#' @return Invisibly returns the base directory path
init_database <- function(base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  markers_dir               <- file.path(base_dir, config$markers_dir)
  marker_sets_subdir        <- file.path(markers_dir, config$marker_sets_subdir)
  genotypes_subdir          <- file.path(markers_dir, config$genotypes_subdir)
  traits_dir                <- file.path(base_dir, config$traits_dir)
  causal_variants_subdir    <- file.path(traits_dir, config$causal_variants_subdir)
  causal_genotypes_subdir   <- file.path(traits_dir, config$causal_genotypes_subdir)
  phenotypes_subdir         <- file.path(traits_dir, config$phenotypes_subdir)
  mappings_dir              <- file.path(base_dir, config$mappings_dir)
  for (dir in c(base_dir, markers_dir, marker_sets_subdir, genotypes_subdir,
                traits_dir, causal_variants_subdir, causal_genotypes_subdir,
                phenotypes_subdir, mappings_dir)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(base_dir)
}


#' Generate marker set ID from population, MAF, species, VCF release ID, and LD threshold
#'
#' @param population Population identifier (normalized to lowercase-trimmed). Canonical: "ce100", "ce96"
#' @param maf MAF threshold as numeric. Do not pre-format — pass raw numeric.
#' @param species Species identifier (normalized to lowercase-trimmed). E.g. "c_elegans", "c_briggsae"
#' @param vcf_release_id 8-digit CaeNDR release date string, e.g. "20220216"
#' @param ms_ld LD R² pruning threshold as numeric. Do not pre-format.
#' @return list with $hash (20-char lowercase hex) and $hash_string (human-readable v=2 input string)
#'
#' Pass $hash (not $hash_string) to generate_trait_id().
#' hash_schema_version prefix "v=2": v=1 hashes are incompatible — regenerate DB.
#'
#' ⚠ CV parameters (cv_maf, cv_ld) are NOT included in this hash. Two runs that share the same
#' (population, maf, species, vcf_release_id, ms_ld) but differ in cv_maf/cv_ld will produce
#' identical marker_set_id values. See README Marker Set ID section for advisory.
generate_marker_set_id <- function(population, maf, species, vcf_release_id, ms_ld) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for generate_marker_set_id()")
  }
  population     <- tolower(trimws(population))
  species        <- tolower(trimws(species))
  hash_string    <- paste0(
    "v=2|population=",      population,
    "|maf=",                sprintf("%.10f", as.numeric(maf)),
    "|species=",            species,
    "|vcf_release_id=",     vcf_release_id,
    "|ms_ld=",              sprintf("%.10f", as.numeric(ms_ld))
  )
  list(
    hash_string = hash_string,
    hash        = substr(digest::digest(hash_string, algo = "sha256", serialize = FALSE), 1, 20)
  )
}


#' Generate deterministic trait ID from parent marker set hash and trait parameters
#'
#' @param marker_set_hash 20-char lowercase hex string from generate_marker_set_id()$hash.
#'   Passing $hash_string by mistake is caught by the format guard.
#' @param nqtl Number of simulated QTLs (integer)
#' @param effect Effect size string — named distribution ("gamma") or numeric range ("0.2-0.3"),
#'   matching a row from the --effect_size CSV. Normalized to lowercase; passed verbatim to hash.
#' @param rep Simulation replicate number (integer)
#' @param h2 Heritability (numeric)
#' @param cv_maf_effective Effective MAF threshold for the CV pool (numeric). Two runs with
#'   identical MS params but different CV pool configs select from different variant universes
#'   and deserve different trait IDs.
#' @param cv_ld LD R² pruning threshold for the CV pool (numeric).
#' @return list with $hash (20-char lowercase hex) and $hash_string (human-readable input)
#'
#' hash_schema_version prefix "v=2": encodes cv_maf_effective and cv_ld; v=1 hashes incompatible
generate_trait_id <- function(marker_set_hash, nqtl, effect, rep, h2,
                               cv_maf_effective, cv_ld) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for generate_trait_id()")
  }
  if (!grepl("^[0-9a-f]{20}$", marker_set_hash)) {
    stop("marker_set_hash must be a 20-character lowercase hex string (pass $hash, not $hash_string)")
  }
  effect      <- tolower(trimws(effect))
  hash_string <- paste0(
    "v=2|parent=", marker_set_hash,
    "|nqtl=",              as.integer(nqtl),
    "|effect=",            effect,
    "|rep=",               as.integer(rep),
    "|h2=",                sprintf("%.10f", as.numeric(h2)),
    "|cv_maf_effective=",  sprintf("%.10f", as.numeric(cv_maf_effective)),
    "|cv_ld=",             sprintf("%.10f", as.numeric(cv_ld))
  )
  list(
    hash_string = hash_string,
    hash        = substr(digest::digest(hash_string, algo = "sha256", serialize = FALSE), 1, 20)
  )
}


#' Generate unique mapping ID from parent trait hash and mapping parameters
#'
#' @param trait_hash 20-char lowercase hex string from generate_trait_id()$hash
#' @param algorithm GWA mode string (normalized to lowercase). Canonical values: "inbred", "loco"
#'   — matches val(mode) emitted by GCTA_PERFORM_GWA in main.nf.
#' @param pca Logical scalar (TRUE/FALSE). Derived from opt$type == "pca" in write_gwa_to_db.R
#' @return list with $hash (20-char lowercase hex) and $hash_string (human-readable input)
#'
#' @note Callers must pass a scalar trait_hash string (from generate_trait_id()$hash), not a
#'   params list. To resolve a mapping ID from a filename-parsed params object, first call
#'   read_marker_set_metadata() -> generate_marker_set_id() -> generate_trait_id() ->
#'   generate_mapping_id(trait$hash, params$algorithm, params$pca)$hash.
#'
#' hash_schema_version prefix "v=1": increment if hash construction rules change; existing DBs regenerated
generate_mapping_id <- function(trait_hash, algorithm, pca) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package 'digest' is required for generate_mapping_id()")
  }
  if (!grepl("^[0-9a-f]{20}$", trait_hash)) {
    stop("trait_hash must be a 20-character lowercase hex string (pass $hash, not $hash_string)")
  }
  stopifnot(is.logical(pca), length(pca) == 1, !is.na(pca))
  algorithm <- tolower(trimws(algorithm))
  pca_str   <- if (pca) "TRUE" else "FALSE"
  hash_string <- paste0(
    "v=1|parent=", trait_hash,
    "|algorithm=", algorithm,
    "|pca=",       pca_str
  )
  list(
    hash_string = hash_string,
    hash        = substr(digest::digest(hash_string, algo = "sha256", serialize = FALSE), 1, 20)
  )
}


#' Build all simulation IDs from a named parameter list
#'
#' Canonical entry point for ID computation in all four DB migration scripts.
#' Adding a new simulation dimension requires updating only this function's
#' generate_*() calls and the sim_params list construction at each call site.
#'
#' @param params Named list with fields: population, maf, species, vcf_release_id,
#'   ms_ld, nqtl, effect, rep, h2, cv_maf_effective, cv_ld.
#'   Use \code{population} (not \code{group}) — the CLI \code{--group} arg maps to
#'   \code{params$population} at all four call sites to match the DB schema column name.
#'   species/vcf_release_id/ms_ld come from \code{read_marker_set_metadata()} at runtime.
#' @param mode GWA mode string ("inbred" or "loco"). If NULL, mapping_id is omitted.
#' @param pca Logical scalar (TRUE/FALSE). Derived from \code{opt$type == "pca"}.
#'   If NULL (or mode is NULL), mapping_id is omitted from the result.
#' @return Named list with \code{ms_id}, \code{trait_id}, and (if both mode and pca
#'   are non-NULL) \code{mapping_id}. Each element is the full ID object from the
#'   corresponding \code{generate_*()} call (use \code{$hash} for the 20-char string).
build_ids_from_params <- function(params, mode = NULL, pca = NULL) {
  if (is.null(params$population) || is.na(params$population))   stop("params$population is NA/NULL")
  if (is.null(params$maf)        || is.na(params$maf))          stop("params$maf is NA/NULL")
  if (is.null(params$species)    || is.na(params$species))      stop("params$species is NA/NULL — check marker set metadata read")
  if (is.null(params$vcf_release_id) || is.na(params$vcf_release_id)) stop("params$vcf_release_id is NA/NULL — check marker set metadata read")
  if (is.null(params$ms_ld) || is.na(as.numeric(params$ms_ld))) stop("params$ms_ld is NA/NULL — check marker set metadata read")
  if (is.null(params$nqtl)       || is.na(params$nqtl))         stop("params$nqtl is NA/NULL")
  if (is.null(params$effect)     || is.na(params$effect))       stop("params$effect is NA/NULL")
  if (is.null(params$rep)        || is.na(params$rep))          stop("params$rep is NA/NULL")
  if (is.null(params$h2)         || is.na(params$h2))           stop("params$h2 is NA/NULL")
  if (is.null(params$cv_maf_effective) || is.na(params$cv_maf_effective)) stop("params$cv_maf_effective is NA/NULL")
  if (is.null(params$cv_ld)      || is.na(params$cv_ld))        stop("params$cv_ld is NA/NULL")

  ms_id    <- generate_marker_set_id(params$population, params$maf, params$species,
                                     params$vcf_release_id, params$ms_ld)
  trait_id <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep,
                                params$h2, params$cv_maf_effective, params$cv_ld)
  result   <- list(ms_id = ms_id, trait_id = trait_id)
  if (!is.null(mode) && !is.null(pca)) {
    result$mapping_id <- generate_mapping_id(trait_id$hash, mode, pca)
  }
  result
}


# ==============================================================================
# Marker Set Operations
# ==============================================================================

.markers_path_from_hash <- function(ms_id_hash, base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  file.path(config$base_dir, config$markers_dir, config$marker_sets_subdir,
            paste0(ms_id_hash, "_markers.parquet"))
}

.genotypes_path_from_hash <- function(ms_id_hash, base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  file.path(config$base_dir, config$markers_dir, config$genotypes_subdir,
            paste0(ms_id_hash, "_genotypes.parquet"))
}


#' Get path to marker set file
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @return Path to marker set Parquet file
get_markers_path <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  result <- generate_marker_set_id(population, maf, species, vcf_release_id, ms_ld)
  # ⚠ CV params (cv_maf, cv_ld) are NOT included in this hash — runs differing only in CV params
  # produce identical marker_set_id values. See README Marker Set ID section for advisory.
  .markers_path_from_hash(result$hash, base_dir)
}


#' Check if marker set exists in database
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @return TRUE if marker set exists
marker_set_exists <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  file.exists(get_markers_path(population, maf, species, vcf_release_id, ms_ld, base_dir))
}


#' Write marker set to database
#'
#' Extracts and stores marker-specific data (marker, CHROM, POS, A1, A2, AF1).
#' Also creates/updates marker set metadata.
#'
#' @param df Mapping dataframe with marker data
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @param overwrite If TRUE, replace existing marker set
#' @param n_independent_tests Number of EIGEN independent tests (optional)
#' @param eigen_source_file Source filename for EIGEN value (optional)
#' @return Invisibly returns the output path
write_marker_set <- function(df, population, maf, species, vcf_release_id, ms_ld,
                             base_dir = "data/db",
                             overwrite = TRUE, n_independent_tests = NA_real_,
                             eigen_source_file = NA_character_,
                             strainfile_hash,
                             strain_list) {
  init_database(base_dir)
  config <- .make_db_config(base_dir)

  ms_id        <- generate_marker_set_id(population, maf, species, vcf_release_id, ms_ld) # compute once
  # ⚠ CV params (cv_maf, cv_ld) are NOT included in this hash — runs differing only in CV params
  # produce identical marker_set_id values. See README Marker Set ID section for advisory.
  markers_path <- .markers_path_from_hash(ms_id$hash, base_dir)  # no second hash call

  if (file.exists(markers_path) && !overwrite) {
    log_msg(glue::glue("Marker set {population}_{maf} already exists - skipping"))
    return(invisible(markers_path))
  }

  marker_cols <- c("marker", "CHROM", "POS", "A1", "A2")  # AF1/population/maf removed

  markers_df <- df %>%
    dplyr::select(dplyr::any_of(marker_cols)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(marker_set_id = ms_id$hash)  # short hash only — no hash_string in data

  markers_df <- prepare_markers_for_parquet(markers_df)

  log_msg(glue::glue("Writing {nrow(markers_df)} markers to: {basename(markers_path)}"))

  arrow::write_parquet(
    markers_df,
    markers_path,
    compression = config$compression
  )

  # Also write marker set metadata
  write_marker_set_metadata(
    population             = population,
    maf                    = maf,
    species                = species,
    vcf_release_id         = vcf_release_id,
    ms_ld                  = ms_ld,
    n_markers              = nrow(markers_df),
    marker_set_id          = ms_id$hash,         # thread through
    marker_set_hash_string = ms_id$hash_string,  # stored once in metadata only
    n_independent_tests    = n_independent_tests,
    eigen_source_file      = eigen_source_file,
    strainfile_hash        = strainfile_hash,
    strain_list            = strain_list,
    base_dir               = base_dir
  )

  invisible(markers_path)
}


#' Read marker set from database
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @return Dataframe with marker data
read_marker_set <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  markers_path <- get_markers_path(population, maf, species, vcf_release_id, ms_ld, base_dir)

  if (!file.exists(markers_path)) {
    stop(glue::glue("Marker set not found: {population}_{maf}"))
  }

  arrow::read_parquet(markers_path) %>%
    as.data.frame()
}


#' List all marker sets in database
#'
#' @param base_dir Database root directory
#' @return Dataframe with population, maf, and marker_set_id columns
list_marker_sets <- function(base_dir = "data/db") {
  meta <- get_all_marker_set_metadata(base_dir)
  if (nrow(meta) == 0) {
    return(data.frame(population = character(), maf = numeric(),
                      marker_set_id = character(), stringsAsFactors = FALSE))
  }
  # intersect() returns only columns present in the metadata file.
  # All expected columns are present after Step 3.
  meta[, intersect(c("population", "maf", "marker_set_id"), names(meta))]
}


# ==============================================================================
# Genotype Matrix Operations
# ==============================================================================

#' Get path to genotype matrix file
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @return Path to genotype matrix Parquet file
get_genotype_matrix_path <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  result <- generate_marker_set_id(population, maf, species, vcf_release_id, ms_ld)
  # ⚠ CV params (cv_maf, cv_ld) are NOT included in this hash — runs differing only in CV params
  # produce identical marker_set_id values. See README Marker Set ID section for advisory.
  .genotypes_path_from_hash(result$hash, base_dir)
}


#' Check if genotype matrix exists in database
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database root directory
#' @return TRUE if genotype matrix exists
genotype_matrix_exists <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  file.exists(get_genotype_matrix_path(population, maf, species, vcf_release_id, ms_ld, base_dir))
}


#' Write genotype matrix to database in long format
#'
#' Reads a wide-format TSV genotype matrix, pivots to long format, and writes
#' as Parquet. The wide-to-long pivot uses data.table::melt() for performance.
#'
#' @param genotype_tsv Path to wide-format TSV from
#'   BCFTOOLS_CREATE_GENOTYPE_MATRIX
#' @param population Population/strain group identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database base directory
#' @param overwrite Overwrite existing file (default TRUE for retry safety)
#' @return Invisible path to written Parquet file
write_genotype_matrix <- function(genotype_tsv, population, maf, species, vcf_release_id, ms_ld,
                                  base_dir = "data/db", overwrite = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for write_genotype_matrix()")
  }

  ms_id    <- generate_marker_set_id(population, maf, species, vcf_release_id, ms_ld) # compute once
  # ⚠ CV params (cv_maf, cv_ld) are NOT included in this hash — runs differing only in CV params
  # produce identical marker_set_id values. See README Marker Set ID section for advisory.
  out_path <- .genotypes_path_from_hash(ms_id$hash, base_dir) # no second hash call

  if (file.exists(out_path) && !overwrite) {
    message("Genotype matrix already exists: ", out_path)
    return(invisible(out_path))
  }

  # Read wide-format TSV
  geno_wide <- data.table::fread(genotype_tsv)

  # Drop REF/ALT (redundant with markers file)
  drop_cols <- intersect(c("REF", "ALT"), names(geno_wide))
  if (length(drop_cols) > 0) {
    geno_wide[, (drop_cols) := NULL]
  }

  # Pivot to long format: CHROM, POS are id vars; strain columns become rows
  id_vars <- c("CHROM", "POS")
  geno_long <- data.table::melt(
    geno_wide,
    id.vars = id_vars,
    variable.name = "strain",
    value.name = "allele"
  )

  # Cast schema
  geno_long[, CHROM := as.character(CHROM)]
  geno_long[, POS := as.integer(POS)]
  geno_long[, strain := as.character(strain)]
  geno_long[, allele := as.numeric(allele)]
  geno_long[, marker_set_id := ms_id$hash]  # join key (hash only)
  # Reorder columns to match schema: marker_set_id must be first
  data.table::setcolorder(geno_long, c("marker_set_id", "CHROM", "POS", "strain", "allele"))

  # Write Parquet with schema enforcement
  arrow::write_parquet(
    arrow::as_arrow_table(geno_long, schema = genotype_matrix_schema()),
    sink = out_path
  )

  invisible(out_path)
}


#' Read genotype matrix from database
#'
#' @param population Population/strain group identifier
#' @param maf MAF threshold
#' @param species Species identifier (e.g. "c_elegans")
#' @param vcf_release_id VCF release date string (e.g. "20220216")
#' @param ms_ld LD R² pruning threshold
#' @param base_dir Database base directory
#' @return data.frame with columns: CHROM, POS, strain, allele (long format)
read_genotype_matrix <- function(population, maf, species, vcf_release_id, ms_ld, base_dir = "data/db") {
  path <- get_genotype_matrix_path(population, maf, species, vcf_release_id, ms_ld, base_dir)
  if (!file.exists(path)) {
    stop("Genotype matrix not found: ", path)
  }
  as.data.frame(arrow::read_parquet(path))
}


# ==============================================================================
# Trait Metadata Operations
# ==============================================================================

#' Get path to trait metadata file
#'
#' @param trait_id Deterministic 20-character hex hash from generate_trait_id()
#' @param base_dir Database base directory
#' @return Path to trait metadata Parquet file
get_trait_metadata_path <- function(trait_id, base_dir) {
  traits_dir <- file.path(base_dir, .db_constants$traits_dir)
  filename <- gsub("\\{trait_id\\}", trait_id, .db_constants$traits_pattern)
  file.path(traits_dir, filename)
}


#' Write trait metadata to database
#'
#' @param trait_id Deterministic 20-character hex hash from generate_trait_id()
#' @param trait_hash_string Human-readable hash input string from generate_trait_id()$hash_string
#' @param marker_set_id Parent marker set hash from generate_marker_set_id()$hash
#' @param nqtl Number of QTLs
#' @param rep Replicate number
#' @param sim_seed Integer seed used for the stochastic simulation
#' @param h2 Heritability
#' @param maf MAF threshold
#' @param effect Effect distribution
#' @param population Strain group ID
#' @param base_dir Database base directory
#' @param overwrite Overwrite existing file (default TRUE for retry safety)
write_trait_metadata <- function(trait_id, trait_hash_string, marker_set_id,
                                 nqtl, rep, sim_seed, h2, maf, effect,
                                 population, cv_maf_effective, cv_ld,
                                 base_dir, overwrite = TRUE) {
  out_path <- get_trait_metadata_path(trait_id, base_dir)
  if (file.exists(out_path) && !overwrite) {
    return(invisible(out_path))
  }
  df <- data.frame(
    trait_id          = trait_id,
    trait_hash_string = as.character(trait_hash_string),
    marker_set_id     = as.character(marker_set_id),
    nqtl             = as.integer(nqtl),
    rep              = as.integer(rep),
    sim_seed         = as.integer(sim_seed),
    h2               = as.numeric(h2),
    maf              = as.numeric(maf),
    effect           = as.character(effect),
    population       = as.character(population),
    cv_maf_effective = as.numeric(cv_maf_effective),
    cv_ld            = as.numeric(cv_ld),
    created_at       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    stringsAsFactors = FALSE
  )
  arrow::write_parquet(
    arrow::as_arrow_table(df, schema = trait_metadata_schema()),
    sink = out_path
  )
  invisible(out_path)
}


#' Read trait metadata from database
#'
#' Provides API symmetry with read_genotype_matrix(), read_causal_variants_data(),
#' and read_phenotype_data(). Required for future offline workflows that enumerate
#' or inspect stored traits. (Addresses DS review L2)
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return data.frame with trait metadata columns
read_trait_metadata <- function(trait_id, base_dir) {
  path <- get_trait_metadata_path(trait_id, base_dir)
  if (!file.exists(path)) {
    stop("Trait metadata not found: ", path)
  }
  as.data.frame(arrow::read_parquet(path))
}


# ==============================================================================
# Causal Variant Operations
# ==============================================================================

#' Get path to causal variants file
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return Path to causal variants Parquet file
get_causal_variants_path <- function(trait_id, base_dir) {
  cv_dir <- file.path(base_dir, .db_constants$traits_dir,
                      .db_constants$causal_variants_subdir)
  filename <- gsub("\\{trait_id\\}", trait_id,
                   .db_constants$causal_variants_pattern)
  file.path(cv_dir, filename)
}


#' Check if causal variants exist in database
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return TRUE if causal variants file exists
causal_variants_exist <- function(trait_id, base_dir) {
  file.exists(get_causal_variants_path(trait_id, base_dir))
}


#' Write causal variants to database
#'
#' @param par_file Path to .par file from GCTA simulation
#' @param trait_id Deterministic trait ID
#' @param base_dir Database base directory
#' @param overwrite Overwrite existing file (default TRUE)
write_causal_variants <- function(par_file, trait_id, base_dir,
                                  overwrite = TRUE) {
  out_path <- get_causal_variants_path(trait_id, base_dir)
  if (file.exists(out_path) && !overwrite) {
    return(invisible(out_path))
  }
  # load_causal_variants() is defined in assessment.R
  cv <- load_causal_variants(par_file)
  arrow::write_parquet(cv, sink = out_path)
  invisible(out_path)
}


#' Read causal variants from database
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return data.frame with causal variant columns
read_causal_variants_data <- function(trait_id, base_dir) {
  path <- get_causal_variants_path(trait_id, base_dir)
  if (!file.exists(path)) {
    stop("Causal variants not found: ", path)
  }
  as.data.frame(arrow::read_parquet(path))
}


# ==============================================================================
# Phenotype Operations
# ==============================================================================

#' Get path to phenotype data file
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return Path to phenotype Parquet file
get_phenotype_path <- function(trait_id, base_dir) {
  pheno_dir <- file.path(base_dir, .db_constants$traits_dir,
                         .db_constants$phenotypes_subdir)
  filename <- gsub("\\{trait_id\\}", trait_id,
                   .db_constants$phenotypes_pattern)
  file.path(pheno_dir, filename)
}


#' Check if phenotype data exists in database
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return TRUE if phenotype file exists
phenotype_exists <- function(trait_id, base_dir) {
  file.exists(get_phenotype_path(trait_id, base_dir))
}


#' Write phenotype data to database
#'
#' IMPORTANT: Stores the pre-upscaled phenotype values from GCTA_SIMULATE_PHENOTYPES,
#' captured BEFORE the pipeline flow continues through PLINK_UPDATE_BY_H2 ->
#' GCTA_MAKE_GRM (which may iteratively scale values by 1000x per round when Vp < 1e-4) ->
#' GCTA_PERFORM_GWA. The stored values may differ from the values used in
#' the actual GWAS mapping by a factor of 1000^Nx.
#'
#' Rationale: For future variance explained estimation, the ANOVA SS ratio
#' (SS_between / SS_total) is scale-invariant — the ratio is identical whether
#' computed on pre-upscaled or post-upscaled phenotypes. Storing pre-upscaled
#' values is preferred because they represent the direct output of GCTA
#' simulation before any pipeline-specific transformations.
#'
#' @param pheno_file Path to .phen file (space-delimited: FID IID value, no header)
#' @param trait_id Deterministic trait ID
#' @param base_dir Database base directory
#' @param overwrite Overwrite existing file (default TRUE)
write_phenotype_data <- function(pheno_file, trait_id, base_dir,
                                 overwrite = TRUE) {
  out_path <- get_phenotype_path(trait_id, base_dir)
  if (file.exists(out_path) && !overwrite) {
    return(invisible(out_path))
  }
  phen <- read.table(pheno_file, header = FALSE,
                     col.names = c("FID", "IID", "value"))
  df <- data.frame(
    strain    = phen$FID,
    phenotype = phen$value,
    stringsAsFactors = FALSE
  )
  arrow::write_parquet(
    arrow::as_arrow_table(df, schema = phenotype_schema()),
    sink = out_path
  )
  invisible(out_path)
}


#' Read phenotype data from database
#'
#' Returns pre-upscaled phenotype values from GCTA simulation. These may differ
#' from the values used in GWAS mapping if variance upscaling was applied by
#' GCTA_MAKE_GRM iterative REML scaling. The ANOVA SS ratio (SS_between / SS_total) is scale-invariant,
#' so these values are appropriate for variance explained estimation.
#'
#' @param trait_id Deterministic 12-character hex hash
#' @param base_dir Database base directory
#' @return data.frame with columns: strain, phenotype
read_phenotype_data <- function(trait_id, base_dir) {
  path <- get_phenotype_path(trait_id, base_dir)
  if (!file.exists(path)) {
    stop("Phenotype data not found: ", path)
  }
  as.data.frame(arrow::read_parquet(path))
}


# ==============================================================================
# Causal Genotype Operations
# ==============================================================================

#' Get path to per-trait causal genotype file
#'
#' @param trait_id Deterministic 20-character hex hash from generate_trait_id()
#' @param base_dir Database base directory
#' @return Path to causal genotype Parquet file
get_causal_genotypes_path <- function(trait_id, base_dir) {
  cg_dir <- file.path(base_dir, .db_constants$traits_dir,
                      .db_constants$causal_genotypes_subdir)
  filename <- gsub("\\{trait_id\\}", trait_id,
                   .db_constants$causal_genotypes_pattern)
  file.path(cg_dir, filename)
}


#' Check if per-trait causal genotype file exists in database
#'
#' @param trait_id Deterministic 20-character hex hash
#' @param base_dir Database base directory
#' @return TRUE if causal genotype file exists
causal_genotypes_exist <- function(trait_id, base_dir) {
  file.exists(get_causal_genotypes_path(trait_id, base_dir))
}


#' Write per-trait causal genotype data to database
#'
#' Reads a TSV file with QTL, strain, allele columns (written by
#' create_causal_vars.py), parses CHROM/POS from the QTL "CHROM:POS" string,
#' adds trait_id, casts types to schema, and writes Parquet.
#'
#' @param causal_geno_file Path to TSV file with QTL, strain, allele columns
#' @param trait_id Deterministic 20-character hex hash from generate_trait_id()
#' @param base_dir Database base directory
#' @param overwrite Overwrite existing file (default TRUE for retry safety)
#' @return Invisibly returns the output path
write_causal_genotypes <- function(causal_geno_file, trait_id, base_dir,
                                   overwrite = TRUE) {
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package 'tidyr' is required for write_causal_genotypes()")
  }
  out_path <- get_causal_genotypes_path(trait_id, base_dir)
  if (file.exists(out_path) && !overwrite) {
    return(invisible(out_path))
  }
  raw <- read.table(causal_geno_file, header = TRUE, sep = "\t",
                    na.strings = "NA", stringsAsFactors = FALSE)
  df <- raw %>%
    tidyr::separate(QTL, into = c("CHROM", "POS"), sep = ":", remove = FALSE,
                    convert = FALSE) %>%
    dplyr::mutate(
      trait_id = trait_id,
      CHROM    = as.character(CHROM),
      POS      = as.integer(POS),
      strain   = as.character(strain),
      allele   = as.numeric(allele)
    ) %>%
    dplyr::select(trait_id, QTL, CHROM, POS, strain, allele)
  arrow::write_parquet(
    arrow::as_arrow_table(df, schema = causal_genotypes_schema()),
    sink = out_path
  )
  invisible(out_path)
}


#' Read per-trait causal genotype data from database
#'
#' Returns genotypes at causal variant positions for a specific trait.
#' The returned data frame has CHROM/POS columns enabling direct
#' inner_join(by = c("CHROM", "POS")) with the marker set genotype matrix.
#'
#' @param trait_id Deterministic 20-character hex hash
#' @param base_dir Database base directory
#' @return data.frame with columns: trait_id, QTL, CHROM, POS, strain, allele
read_causal_genotypes <- function(trait_id, base_dir) {
  path <- get_causal_genotypes_path(trait_id, base_dir)
  if (!file.exists(path)) {
    stop("Causal genotypes not found: ", path)
  }
  as.data.frame(arrow::read_parquet(path))
}


# ==============================================================================
# Marker Set Metadata Operations
# ==============================================================================

#' Get path to marker set metadata file
#'
#' Get per-population marker set metadata file path
#'
#' Returns the path to an individual population's metadata Parquet file.
#' Each DB_MIGRATION_WRITE_MARKER_SET task writes its own file, eliminating
#' the concurrent read-modify-write race on a shared file.
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param base_dir Database root directory
#' @return Path to per-population metadata Parquet file
get_per_population_metadata_path <- function(population, maf, base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  metadata_dir <- file.path(config$base_dir, config$marker_set_metadata_dir)
  file.path(metadata_dir, paste0(population, "_", maf, "_metadata.parquet"))
}


#' Get the marker set metadata directory path
#'
#' @param base_dir Database root directory
#' @return Path to the per-population metadata directory
get_marker_set_metadata_dir <- function(base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  file.path(config$base_dir, config$marker_set_metadata_dir)
}


#' Write or update marker set metadata
#'
#' Stores metadata for a marker set including the number of independent tests
#' from EIGEN decomposition. Updates existing record or creates new one.
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param n_markers Number of markers in the set
#' @param n_independent_tests Number of independent tests from EIGEN (optional)
#' @param eigen_source_file Source filename for EIGEN value (optional)
#' @param base_dir Database root directory
#' @return Invisibly returns the metadata file path
write_marker_set_metadata <- function(population, maf, species, vcf_release_id, ms_ld,
                                       n_markers,
                                       marker_set_id, marker_set_hash_string,
                                       hash_schema_version = "v=2",
                                       n_independent_tests = NA_real_,
                                       eigen_source_file = NA_character_,
                                       strainfile_hash,
                                       strain_list,
                                       base_dir = "data/db") {
  config <- .make_db_config(base_dir)

  # Per-population metadata file — each task writes its own file,
  # eliminating the concurrent read-modify-write race condition that
  # caused lost population records on HPC (SLURM) with multiple
  # DB_MIGRATION_WRITE_MARKER_SET tasks running in parallel.
  metadata_path <- get_per_population_metadata_path(population, maf, base_dir)
  metadata_dir <- dirname(metadata_path)
  if (!dir.exists(metadata_dir)) {
    dir.create(metadata_dir, recursive = TRUE)
  }

  new_record <- data.frame(
    marker_set_id          = as.character(marker_set_id),
    marker_set_hash_string = as.character(marker_set_hash_string),
    hash_schema_version    = as.character(hash_schema_version),
    population             = as.character(population),
    maf                    = as.numeric(maf),
    species                = as.character(species),
    vcf_release_id         = as.character(vcf_release_id),
    ms_ld                  = as.numeric(ms_ld),
    n_markers              = as.integer(n_markers),
    n_independent_tests    = as.numeric(n_independent_tests),
    eigen_source_file      = as.character(eigen_source_file),
    strainfile_hash        = as.character(strainfile_hash),
    strain_list            = as.character(strain_list),
    created_at             = Sys.time(),
    stringsAsFactors       = FALSE
  )

  arrow::write_parquet(new_record, metadata_path, compression = config$compression)
  log_msg(glue::glue("Wrote marker set metadata: {population}_{maf}"))

  invisible(metadata_path)
}


#' Read marker set metadata for a specific marker set
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param base_dir Database root directory
#' @return Named list with marker set metadata, or NULL if not found
read_marker_set_metadata <- function(population, maf, base_dir = "data/db") {
  per_pop_path <- get_per_population_metadata_path(population, maf, base_dir)
  if (!file.exists(per_pop_path)) {
    return(NULL)
  }

  record <- as.data.frame(arrow::read_parquet(per_pop_path))
  if (nrow(record) == 0) {
    return(NULL)
  }

  as.list(record[1, ])
}


#' Get all marker set metadata
#'
#' @param base_dir Database root directory
#' @return Dataframe with all marker set metadata
get_all_marker_set_metadata <- function(base_dir = "data/db") {
  metadata_dir <- get_marker_set_metadata_dir(base_dir)
  if (!dir.exists(metadata_dir)) {
    return(data.frame(
      population = character(),
      maf = numeric(),
      n_markers = integer(),
      n_independent_tests = numeric(),
      eigen_source_file = character(),
      created_at = as.POSIXct(character())
    ))
  }

  files <- list.files(metadata_dir, pattern = "_metadata\\.parquet$", full.names = TRUE)
  if (length(files) == 0) {
    return(data.frame(
      population = character(),
      maf = numeric(),
      n_markers = integer(),
      n_independent_tests = numeric(),
      eigen_source_file = character(),
      created_at = as.POSIXct(character())
    ))
  }

  records <- lapply(files, function(f) as.data.frame(arrow::read_parquet(f)))
  dplyr::bind_rows(records)
}


#' Get number of independent tests for a marker set
#'
#' Returns the EIGEN-derived number of independent tests for computing
#' the EIGEN significance threshold.
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param base_dir Database root directory
#' @return Numeric value of independent tests, or NA if not available
get_n_independent_tests <- function(population, maf, base_dir = "data/db") {
  metadata <- read_marker_set_metadata(population, maf, base_dir)

  if (is.null(metadata)) {
    warning(glue::glue("No metadata found for marker set: {population}_{maf}"))
    return(NA_real_)
  }

  n_tests <- metadata$n_independent_tests

  if (is.na(n_tests)) {
    warning(glue::glue("No EIGEN data available for marker set: {population}_{maf}"))
  }

  n_tests
}


#' Update EIGEN value for existing marker set
#'
#' @param population Population identifier
#' @param maf MAF threshold
#' @param n_independent_tests Number of independent tests from EIGEN
#' @param eigen_source_file Source filename for EIGEN value
#' @param base_dir Database root directory
#' @return Invisibly returns TRUE if updated, FALSE if marker set not found
update_marker_set_eigen <- function(population, maf, n_independent_tests,
                                     eigen_source_file = NA_character_,
                                     base_dir = "data/db") {
  stop("update_marker_set_eigen() is deprecated — use write_marker_set() directly (write_marker_set.R module).")

  existing <- read_marker_set_metadata(population, maf, base_dir)

  if (is.null(existing)) {
    warning(glue::glue("Marker set not found: {population}_{maf}"))
    return(invisible(FALSE))
  }

  write_marker_set_metadata(
    population = population,
    maf = maf,
    n_markers = existing$n_markers,
    n_independent_tests = n_independent_tests,
    eigen_source_file = eigen_source_file,
    base_dir = base_dir
  )

  invisible(TRUE)
}


#' Prepare markers dataframe for Parquet storage
#'
#' @param df Markers dataframe
#' @return Dataframe with proper types
prepare_markers_for_parquet <- function(df) {
  char_cols <- c("marker", "CHROM", "A1", "A2", "population")
  for (col in char_cols) {
    if (col %in% names(df)) {
      df[[col]] <- as.character(df[[col]])
    }
  }

  if ("POS" %in% names(df)) {
    df$POS <- as.integer(df$POS)
  }

  num_cols <- c("AF1", "maf")
  for (col in num_cols) {
    if (col %in% names(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }

  df
}


# ==============================================================================
# Mapping Operations
# ==============================================================================

#' Get path to mappings file for a population
#'
#' @param population Population identifier
#' @param base_dir Database root directory
#' @return Path to mappings Parquet file
get_mappings_path <- function(population, base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  file.path(
    config$base_dir,
    config$mappings_dir,
    glue::glue(config$mappings_pattern, population = population)
  )
}


#' Write mapping to database
#'
#' Stores mapping-specific data (BETA, SE, P) for a single GWAS mapping.
#' Also creates/updates the marker set if needed.
#' Deduplicates by CHROM:POS before storage (LOCO files have duplicate rows).
#'
#' @param df Mapping dataframe (raw, without significance processing)
#' @param source_file Original source filename
#' @param base_dir Database root directory
#' @param overwrite If TRUE, replace existing data for this mapping_id
#' @return Invisibly returns the output file path
write_mapping_to_db <- function(df, source_file, base_dir = "data/db",
                                 overwrite = TRUE) {
  stop("write_mapping_to_db() is deprecated — use the write_gwa_to_db.R module (modules/db_migration/write_gwa_to_db).")

  config <- .make_db_config(base_dir)
  init_database(base_dir)

  params <- parse_mapping_filename(source_file)
  if (is.null(params)) {
    stop("Could not parse mapping filename - cannot determine parameters for storage")
  }

  ms_id       <- generate_marker_set_id(params$population, params$maf)
  trait_obj   <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
  mapping_obj <- generate_mapping_id(trait_obj$hash, params$algorithm, params$pca)
  mapping_id  <- mapping_obj$hash
  population  <- params$population
  maf         <- params$maf

  # Deduplicate by CHROM:POS (LOCO files have duplicate rows per marker)
  n_before <- nrow(df)
  df <- df %>%
    dplyr::distinct(CHROM, POS, .keep_all = TRUE)
  n_after <- nrow(df)
  if (n_before != n_after) {
    log_msg(glue::glue("Deduplicated: {n_before} -> {n_after} unique markers (removed {n_before - n_after} duplicates)"))
  }

  log_msg(glue::glue("Writing mapping {mapping_id} to database"))

  # 1. Write marker set (only if not already present)
  write_marker_set(df, population, maf, base_dir, overwrite = FALSE)

  # 2. Prepare mapping data (statistics only)
  mapping_df <- prepare_mapping_data(df, params)

  # 3. Write/merge mapping data
  mappings_path <- get_mappings_path(population, base_dir)

  if (file.exists(mappings_path)) {
    merge_mapping_data(mappings_path, mapping_df, mapping_id, overwrite, base_dir)
  } else {
    arrow::write_parquet(
      mapping_df,
      mappings_path,
      compression = config$compression
    )
    log_msg(glue::glue("Created new mappings file with {nrow(mapping_df)} rows"))
  }

  # 4. Update metadata table
  update_metadata_table(df, params, source_file, base_dir, ms_id, trait_obj, mapping_obj)

  invisible(mappings_path)
}


#' Prepare mapping data for storage
#'
#' Selects mapping statistics columns and adds metadata columns.
#' Note: AF1 included, N and log10p removed from schema.
#'
#' @param df Raw mapping dataframe
#' @param params Parsed filename parameters
#' @return Dataframe ready for storage
prepare_mapping_data <- function(df, params) {
  mapping_cols   <- c("marker", "AF1", "BETA", "SE", "P")
  available_cols <- intersect(mapping_cols, names(df))
  df %>%
    dplyr::select(dplyr::all_of(available_cols)) %>%
    prepare_mappings_for_parquet()
}


#' Prepare mappings dataframe for Parquet storage
#'
#' @param df Mappings dataframe
#' @return Dataframe with proper types
prepare_mappings_for_parquet <- function(df) {
  char_cols <- c("marker", "mapping_id", "marker_set_id", "trait_id")
  for (col in char_cols) {
    if (col %in% names(df)) {
      df[[col]] <- as.character(df[[col]])
    }
  }

  num_cols <- c("AF1", "BETA", "SE", "P")
  for (col in num_cols) {
    if (col %in% names(df)) {
      df[[col]] <- as.numeric(df[[col]])
    }
  }

  df
}


#' Merge new mapping data with existing Parquet file
#'
#' @param parquet_file Path to existing Parquet file
#' @param new_data Dataframe with new mapping data
#' @param mapping_id Mapping ID to merge
#' @param overwrite If TRUE, replace existing data for this mapping_id
#' @param base_dir Database root directory
merge_mapping_data <- function(parquet_file, new_data, mapping_id, overwrite, base_dir) {
  config <- .make_db_config(base_dir)
  existing <- as.data.frame(arrow::read_parquet(parquet_file))

  if (overwrite) {
    existing_filtered <- existing[existing$mapping_id != mapping_id, , drop = FALSE]
    result <- dplyr::bind_rows(existing_filtered, new_data)
    log_msg(glue::glue("Replaced mapping {mapping_id} - now {nrow(result)} total rows"))
  } else {
    existing_count <- sum(existing$mapping_id == mapping_id)

    if (existing_count > 0) {
      log_msg(glue::glue("Mapping {mapping_id} already exists ({existing_count} rows) - skipping"),
              level = "WARN")
      return(invisible(NULL))
    }

    result <- dplyr::bind_rows(existing, new_data)
    log_msg(glue::glue("Appended mapping {mapping_id} - now {nrow(result)} total rows"))
  }

  arrow::write_parquet(
    result,
    parquet_file,
    compression = config$compression
  )
}


#' Update the mappings metadata table
#'
#' @param df Mapping dataframe
#' @param params Parsed filename parameters
#' @param source_file Original source filename
#' @param base_dir Database root directory
update_metadata_table <- function(df, params, source_file, base_dir, ms_id, trait_obj, mapping_obj) {
  config        <- .make_db_config(base_dir)
  metadata_file <- file.path(config$base_dir, config$metadata_file)

  meta_record <- data.frame(
    mapping_id          = mapping_obj$hash,
    mapping_hash_string = mapping_obj$hash_string,
    trait_id            = trait_obj$hash,
    marker_set_id       = ms_id$hash,
    hash_schema_version = "v=1",
    population          = params$population,
    maf                 = as.numeric(params$maf),
    nqtl                = as.integer(params$nqtl),
    rep                 = as.integer(params$rep),
    h2                  = as.numeric(params$h2),
    effect              = params$effect,
    algorithm           = params$algorithm,
    pca                 = as.logical(params$pca),
    n_markers           = nrow(df),
    source_file         = basename(source_file),
    processed_at        = Sys.time(),
    processing_version  = "2.0.0",
    stringsAsFactors    = FALSE
  )

  if (file.exists(metadata_file)) {
    existing_meta <- arrow::read_parquet(metadata_file) %>%
      dplyr::filter(mapping_id != !!mapping_obj$hash)
    updated_meta <- dplyr::bind_rows(existing_meta, meta_record)
  } else {
    updated_meta <- meta_record
  }

  arrow::write_parquet(updated_meta, metadata_file)
  log_msg(glue::glue("Updated metadata table ({nrow(updated_meta)} mappings)"))
}


# ==============================================================================
# Batch Operations
# ==============================================================================

#' Create a metadata record for a mapping
#'
#' @param df Mapping dataframe
#' @param params Parsed filename parameters
#' @param source_file Original source filename
#' @param ms_id Marker set ID object from generate_marker_set_id()
#' @param trait_obj Trait ID object from generate_trait_id()
#' @param mapping_obj Mapping ID object from generate_mapping_id()
#' @return Dataframe with single metadata row
create_metadata_record <- function(df, params, source_file, ms_id, trait_obj, mapping_obj) {
  data.frame(
    mapping_id          = mapping_obj$hash,
    mapping_hash_string = mapping_obj$hash_string,
    trait_id            = trait_obj$hash,
    marker_set_id       = ms_id$hash,
    hash_schema_version = "v=1",
    population          = params$population,
    maf                 = as.numeric(params$maf),
    nqtl                = as.integer(params$nqtl),
    rep                 = as.integer(params$rep),
    h2                  = as.numeric(params$h2),
    effect              = params$effect,
    algorithm           = params$algorithm,
    pca                 = as.logical(params$pca),
    n_markers           = nrow(df),
    source_file         = basename(source_file),
    processed_at        = Sys.time(),
    processing_version  = "2.1.0",
    stringsAsFactors    = FALSE
  )
}


#' Batch write all metadata records at once
#'
#' @param metadata_list List of metadata dataframes
#' @param base_dir Database root directory
batch_write_metadata <- function(metadata_list, base_dir) {
  if (length(metadata_list) == 0) return(invisible(NULL))

  config <- .make_db_config(base_dir)
  metadata_file <- file.path(config$base_dir, config$metadata_file)
  new_metadata <- dplyr::bind_rows(metadata_list)
  new_mapping_ids <- unique(new_metadata$mapping_id)

  log_msg(glue::glue("Writing metadata for {length(new_mapping_ids)} mappings"))

  if (file.exists(metadata_file)) {
    existing <- as.data.frame(arrow::read_parquet(metadata_file))
    # Remove any existing entries for mapping_ids we're adding
    existing_filtered <- existing[!existing$mapping_id %in% new_mapping_ids, ]
    result <- dplyr::bind_rows(existing_filtered, new_metadata)
  } else {
    result <- new_metadata
  }

  arrow::write_parquet(result, metadata_file)
  log_msg(glue::glue("Metadata table updated: {nrow(result)} total mappings"))

  invisible(metadata_file)
}


# ==============================================================================
# Partitioned Dataset Operations (True Parallel Writes)
# ==============================================================================

#' Get path to partitioned mapping directory
#'
#' Returns the path for a specific mapping in the partitioned structure:
#' mappings/population={pop}/mapping_id={id}/
#'
#' @param population Population identifier
#' @param mapping_id Mapping identifier
#' @param base_dir Database root directory
#' @return Path to partition directory
get_partition_path <- function(population, mapping_id, base_dir = "data/db") {
  config <- .make_db_config(base_dir)
  file.path(
    config$base_dir,
    config$mappings_dir,
    glue::glue("population={population}"),
    glue::glue("mapping_id={mapping_id}")
  )
}


#' Write a single mapping to partitioned storage
#'
#' Writes mapping data to its own partition directory. This enables true parallel
#' writes since each mapping goes to a unique location with no coordination needed.
#'
#' Directory structure: mappings/population={pop}/mapping_id={id}/data.parquet
#'
#' @param df Mapping dataframe (raw, deduplicated)
#' @param params Parsed filename parameters
#' @param base_dir Database root directory
#' @return Invisibly returns the partition path
write_mapping_partitioned <- function(df, params, ms_id, trait_id, base_dir = "data/db") {
  config     <- .make_db_config(base_dir)
  mapping    <- generate_mapping_id(trait_id$hash, params$algorithm, params$pca)
  mapping_id <- mapping$hash

  mapping_df <- prepare_mapping_data(df, params)

  partition_path <- get_partition_path(params$population, mapping_id, base_dir)
  dir.create(partition_path, recursive = TRUE, showWarnings = FALSE)

  arrow::write_parquet(mapping_df, file.path(partition_path, "data.parquet"),
                       compression = config$compression)

  invisible(partition_path)
}


#' Write per-mapping metadata sidecar
#'
#' Writes a single-row metadata parquet to {partition_path}/meta.parquet.
#' Contains all metadata_schema() columns. Called by write_gwa_to_db.R
#' in parallel; safe because each mapping writes to its own directory.
#'
#' @param params List with population, maf, nqtl, effect, rep, h2, algorithm, pca
#' @param ms_id  Marker set ID object from generate_marker_set_id()
#' @param trait_id Trait ID object from generate_trait_id()
#' @param n_markers Integer, number of markers in the mapping
#' @param base_dir Database root directory
write_mapping_metadata <- function(params, ms_id, trait_id, n_markers, base_dir = "data/db") {
  config   <- .make_db_config(base_dir)
  mapping  <- generate_mapping_id(trait_id$hash, params$algorithm, params$pca)
  part_dir <- get_partition_path(params$population, mapping$hash, base_dir)
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)

  meta <- data.frame(
    mapping_id          = mapping$hash,
    mapping_hash_string = mapping$hash_string,
    trait_id            = trait_id$hash,
    marker_set_id       = ms_id$hash,
    hash_schema_version = "v=1",
    population          = params$population,
    maf                 = as.numeric(params$maf),
    nqtl                = as.integer(params$nqtl),
    rep                 = as.integer(params$rep),
    h2                  = as.numeric(params$h2),
    effect              = params$effect,
    algorithm           = params$algorithm,
    pca                 = as.logical(params$pca),
    n_markers           = as.integer(n_markers),
    source_file         = "",
    processed_at        = Sys.time(),
    processing_version  = "2.1.0-nf",
    stringsAsFactors    = FALSE
  )

  arrow::write_parquet(
    arrow::as_arrow_table(meta, schema = metadata_schema()),
    sink = file.path(part_dir, "meta.parquet"),
    compression = config$compression
  )
  invisible(file.path(part_dir, "meta.parquet"))
}
