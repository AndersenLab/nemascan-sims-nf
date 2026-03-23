#!/usr/bin/env Rscript
# aggregate_metadata.R - Build metadata table by scanning parquet files
#
# Usage: Rscript aggregate_metadata.R <base_dir>
# Output: mappings_metadata.parquet and aggregation_summary.txt in base_dir

library(arrow)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript aggregate_metadata.R <base_dir>")
}

base_dir <- args[1]

# Source R modules - R_SOURCE_DIR env var is required (set by Nextflow process script block)
# For standalone usage: R_SOURCE_DIR=./R Rscript aggregate_metadata.R <args>
r_source_dir <- Sys.getenv("R_SOURCE_DIR", unset = "")
if (r_source_dir == "") {
  stop("R_SOURCE_DIR environment variable must be set. Example: R_SOURCE_DIR=./R Rscript aggregate_metadata.R <args>")
}
if (!dir.exists(r_source_dir)) {
  stop(paste("R_SOURCE_DIR does not exist:", r_source_dir))
}

source(file.path(r_source_dir, "utils.R"))
source(file.path(r_source_dir, "database.R"))

# Find all meta sidecar files (one per mapping, written by write_gwa_to_db.R)
mappings_dir <- file.path(base_dir, "mappings")
meta_files <- list.files(
  mappings_dir,
  pattern = "meta\\.parquet$",
  full.names = TRUE,
  recursive = TRUE
)

cat("Found", length(meta_files), "mapping metadata sidecars\n")

if (length(meta_files) == 0) {
  # Fallback: -resume runs skip WRITE_GWA_TO_DB (cached), so no meta.parquet
  # sidecars are written. Fall back to the legacy data.parquet scan.
  # Simulation parameters (maf, nqtl, h2, etc.) cannot be recovered from
  # data.parquet — they will be NA in the output. DB_MIGRATION_ANALYZE_QTL
  # will still fail for those mappings. Re-run without -resume to fix.
  warning(
    "No meta.parquet sidecars found — falling back to data.parquet scan. ",
    "Simulation parameters (maf, nqtl, h2, etc.) and FK references (marker_set_id, trait_id) will be NA. ",
    "Re-run without -resume to generate complete metadata."
  )
  data_files <- list.files(
    mappings_dir,
    pattern = "data\\.parquet$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(data_files) == 0) {
    cat("No data.parquet or meta.parquet files found. Nothing to aggregate.\n")
    writeLines("No mappings found", "aggregation_summary.txt")
    quit(status = 0)
  }
  metadata_list <- lapply(data_files, function(f) {
    tryCatch({
      partition_dir <- dirname(f)
      mapping_id    <- sub("^mapping_id=", "", basename(partition_dir))
      n_markers     <- nrow(as.data.frame(arrow::read_parquet(f)))
      data.frame(
        mapping_id          = mapping_id,
        marker_set_id       = NA_character_,
        trait_id            = NA_character_,
        n_markers           = as.integer(n_markers),
        processed_at        = Sys.time(),
        processing_version  = "2.1.0-nf",
        maf                 = NA_real_,
        nqtl                = NA_integer_,
        rep                 = NA_integer_,
        h2                  = NA_real_,
        effect              = NA_character_,
        algorithm           = NA_character_,
        pca                 = NA,
        population          = NA_character_,
        hash_schema_version = NA_character_,
        mapping_hash_string = NA_character_,
        source_file         = NA_character_,
        stringsAsFactors    = FALSE
      )
    }, error = function(e) { warning(paste("Failed to read:", f, "-", e$message)); NULL })
  })
} else {
  # Primary path: combine all sidecars — each is already a complete metadata row
  metadata_list <- lapply(meta_files, function(f) {
    tryCatch(
      as.data.frame(arrow::read_parquet(f)),
      error = function(e) { warning(paste("Failed to read:", f, "-", e$message)); NULL }
    )
  })
}

metadata_list <- Filter(Negate(is.null), metadata_list)
cat("Successfully read", length(metadata_list), "metadata entries\n")

if (length(metadata_list) == 0) {
  cat("No valid metadata files - nothing to aggregate\n")
  writeLines("No valid mappings found", "aggregation_summary.txt")
  quit(status = 0)
}

metadata_df <- bind_rows(metadata_list)

# Write metadata parquet
metadata_path <- file.path(base_dir, "mappings_metadata.parquet")
arrow::write_parquet(metadata_df, metadata_path, compression = "snappy")
cat("Wrote metadata to:", metadata_path, "\n")

# Write aggregation summary to the declared process output file.
# Warnings written only to stderr are captured in .command.err in the work
# directory but are NOT visible in the Nextflow run log for processes that
# exit 0. Writing to aggregation_summary.txt ensures operators can see them
# after the run without digging into work directories.
summary_lines <- c(
  paste("Aggregation complete"),
  paste("Total mappings:", nrow(metadata_df)),
  paste("Unique mapping IDs:", length(unique(metadata_df$mapping_id)))
)

# Orphan detection: partitions with data.parquet but no meta.parquet.
# Indicates write_gwa_to_db crashed after writing data but before writing
# the sidecar. Only meaningful in the primary path — in the fallback, no
# meta.parquet exists at all and every partition would appear orphaned.
if (length(meta_files) > 0) {
  data_dirs <- unique(dirname(list.files(
    mappings_dir, pattern = "data\\.parquet$", full.names = TRUE, recursive = TRUE
  )))
  meta_dirs <- unique(dirname(meta_files))
  orphan_dirs <- setdiff(data_dirs, meta_dirs)

  if (length(orphan_dirs) > 0) {
    orphan_msg <- paste0(
      "WARN: ", length(orphan_dirs), " orphaned partition(s) ",
      "(data.parquet without meta.parquet):\n",
      paste(orphan_dirs, collapse = "\n")
    )
    warning(orphan_msg)
    summary_lines <- c(summary_lines, orphan_msg)
  }
}

writeLines(summary_lines, "aggregation_summary.txt")

cat("\nAggregation complete:\n")
cat("  Total mappings:", nrow(metadata_df), "\n")
