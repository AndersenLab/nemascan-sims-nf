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

# Find all mapping parquet files
mappings_dir <- file.path(base_dir, "mappings")
parquet_files <- list.files(
  mappings_dir,
  pattern = "data\\.parquet$",
  full.names = TRUE,
  recursive = TRUE
)

cat("Found", length(parquet_files), "mapping parquet files\n")

if (length(parquet_files) == 0) {
  cat("No parquet files found - nothing to aggregate\n")
  writeLines("No mappings found", "aggregation_summary.txt")
  quit(status = 0)
}

# Extract metadata from each parquet file (read just one row for efficiency)
metadata_list <- lapply(parquet_files, function(pq_file) {
  tryCatch({
    # Read first row to get metadata columns
    df <- arrow::read_parquet(pq_file) %>%
      head(1) %>%
      select(any_of(c(
        "mapping_id", "population", "maf", "nqtl", "rep",
        "h2", "effect", "algorithm", "pca", "trait"
      )))

    # Get marker count from full file
    n_markers <- arrow::read_parquet(pq_file) %>% nrow()

    df %>%
      mutate(
        n_markers = n_markers,
        source_file = basename(dirname(dirname(pq_file))),
        processed_at = Sys.time(),
        processing_version = "2.1.0-nf"
      )
  }, error = function(e) {
    warning(paste("Failed to read:", pq_file, "-", e$message))
    NULL
  })
})

# Remove failed reads and combine
metadata_list <- Filter(Negate(is.null), metadata_list)
cat("Successfully read", length(metadata_list), "parquet files\n")

if (length(metadata_list) == 0) {
  cat("No valid parquet files - nothing to aggregate\n")
  writeLines("No valid mappings found", "aggregation_summary.txt")
  quit(status = 0)
}

# Combine into single dataframe
metadata_df <- bind_rows(metadata_list)

# Write metadata parquet
metadata_path <- file.path(base_dir, "mappings_metadata.parquet")
arrow::write_parquet(metadata_df, metadata_path, compression = "snappy")
cat("Wrote metadata to:", metadata_path, "\n")

# Write summary
summary_text <- paste0(
  "Aggregation complete\n",
  "Total mappings: ", nrow(metadata_df), "\n",
  "Populations: ", paste(unique(metadata_df$population), collapse = ", "), "\n",
  "MAF values: ", paste(unique(metadata_df$maf), collapse = ", "), "\n"
)
writeLines(summary_text, "aggregation_summary.txt")

cat("\nAggregation complete:\n")
cat("  Total mappings:", nrow(metadata_df), "\n")
cat("  Populations:", length(unique(metadata_df$population)), "\n")
