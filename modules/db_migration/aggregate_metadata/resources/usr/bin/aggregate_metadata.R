#!/usr/bin/env Rscript
# aggregate_metadata.R - Build metadata table by scanning parquet files
#
# Usage: Rscript aggregate_metadata.R <base_dir>
# Output: mappings_metadata.parquet and aggregation_summary.txt in base_dir

library(arrow)
library(dplyr)
library(duckdb)

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

# Enumerate partition directories
mappings_dir   <- file.path(base_dir, "mappings")
partition_dirs <- Sys.glob(file.path(mappings_dir, "population=*", "mapping_id=*"))

# Check for meta.parquet and data.parquet files in each partition directory
has_meta <- file.exists(file.path(partition_dirs, "meta.parquet"))
has_data <- file.exists(file.path(partition_dirs, "data.parquet"))

# Identify mappings with valid metadata sidecars and orphaned partitions (data.parquet without meta.parquet)
good_meta_files <- file.path(partition_dirs[has_meta], "meta.parquet")
orphan_dirs     <- partition_dirs[has_data & !has_meta]

cat("Found", length(good_meta_files), "mapping metadata sidecars\n")

# Fail-fast when no sidecars exist
if (length(good_meta_files) == 0) {
  if (any(has_data)) {
    msg <- paste0(
      "Zero meta.parquet sidecars found but ", sum(has_data), " data.parquet partition(s) exist.\n",
      "This indicates WRITE_GWA_TO_DB crashed after writing data.parquet but before writing meta.parquet.\n",
      "Affected partitions:\n",
      paste(partition_dirs[has_data], collapse = "\n")
    )
    writeLines(c("FATAL: no meta sidecars", msg), "aggregation_summary.txt")
    stop(msg)
  } else {
    cat("No data.parquet or meta.parquet files found. Nothing to aggregate.\n")
    writeLines("No mappings found", "aggregation_summary.txt")
    quit(status = 0)
  }
}

# Read all sidecars in one multithreaded Arrow scan. All sidecars are written
# with as_arrow_table(schema = metadata_schema()) so they are schema-identical;
# The explicit vetted vector used to differentiate between 
# meta.parquet and data.parquet in the same directory tree.
arrow::set_cpu_count(as.integer(Sys.getenv("NF_TASK_CPUS", "4")))

ds  <- arrow::open_dataset(good_meta_files, format = "parquet")
tbl <- dplyr::compute(ds)

metadata_path <- file.path(base_dir, "mappings_metadata.parquet")
arrow::write_parquet(tbl, metadata_path, compression = "snappy")
cat("Wrote metadata to:", metadata_path, "\n")


summary_lines <- c("Aggregation complete")

if (length(orphan_dirs) > 0) {
  orphan_msg <- paste0(
    "WARN: ", length(orphan_dirs), " orphaned partition(s) ",
    "(data.parquet without meta.parquet):\n",
    paste(orphan_dirs, collapse = "\n")
  )
  warning(orphan_msg)
  summary_lines <- c(summary_lines, orphan_msg)
}

# Compute summary stats in one DuckDB pass over the already-written output;
# avoids re-reading 50k sidecars for counts.
con <- duckdb::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
on.exit(duckdb::dbDisconnect(con, shutdown = TRUE), add = TRUE)

stats <- duckdb::dbGetQuery(con, sprintf(
  "SELECT COUNT(*) AS n_rows, COUNT(DISTINCT mapping_id) AS n_unique
   FROM read_parquet('%s')",
  metadata_path
))

summary_lines <- c(
  summary_lines,
  paste("Total mappings:", stats$n_rows),
  paste("Unique mapping IDs:", stats$n_unique)
)

writeLines(summary_lines, "aggregation_summary.txt")

cat("\nAggregation complete:\n")
cat("  Total mappings:", stats$n_rows, "\n")
cat("  Unique mapping IDs:", stats$n_unique, "\n")
