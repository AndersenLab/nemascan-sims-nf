#!/usr/bin/env Rscript
# write_gwa_to_db.R - Write raw GWA output to Parquet database
#
# Reads a GCTA GWA output file (.fastGWA or .mlma), constructs mapping params
# from Nextflow channel metadata, and writes to Hive-partitioned Parquet.
#
# Usage: write_gwa_to_db.R --group <group> --maf <maf> --nqtl <nqtl>
#            --effect <effect> --rep <rep> --h2 <h2> --mode <mode>
#            --type <type> --gwa_file <file> --base_dir <db_dir>

library(optparse)

option_list <- list(
  make_option("--group", type = "character", help = "Population/strain group identifier"),
  make_option("--maf", type = "character", help = "MAF threshold"),
  make_option("--nqtl", type = "integer", help = "Number of simulated QTLs"),
  make_option("--effect", type = "character", help = "Effect size distribution (gamma/uniform)"),
  make_option("--rep", type = "integer", help = "Simulation replicate number"),
  make_option("--h2", type = "double", help = "Heritability"),
  make_option("--mode", type = "character", help = "GWA mode (inbred/loco)"),
  make_option("--type", type = "character", help = "PCA type (pca/nopca)"),
  make_option("--gwa_file", type = "character", help = "Path to GWA output file"),
  make_option("--base_dir", type = "character", help = "Database output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c("group", "maf", "nqtl", "effect", "rep", "h2", "mode", "type", "gwa_file", "base_dir")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}

# Source R library via R_SOURCE_DIR (set by Nextflow process script block)
r_source_dir <- Sys.getenv("R_SOURCE_DIR", unset = "")
if (r_source_dir == "") {
  stop("R_SOURCE_DIR environment variable must be set")
}
source(file.path(r_source_dir, "utils.R"))
source(file.path(r_source_dir, "io.R"))
source(file.path(r_source_dir, "database.R"))

# Read raw GWA output (auto-detects .fastGWA vs .mlma format)
log_msg(paste("Reading GWA file:", opt$gwa_file))
gwa_df <- read_raw_gwa_file(opt$gwa_file)

# Deduplicate by marker (CHROM:POS) — LOCO .mlma files may contain duplicate rows
# Matches dedup pattern in write_mapping_to_db() (database.R)
n_before <- nrow(gwa_df)
gwa_df <- gwa_df %>% dplyr::distinct(marker, .keep_all = TRUE)
n_after <- nrow(gwa_df)
if (n_before != n_after) {
  log_msg(paste("Deduplicated:", n_before, "->", n_after, "unique markers (removed",
                n_before - n_after, "duplicates)"))
}

# Construct mapping params from Nextflow channel metadata
# algorithm = opt$mode ("inbred" or "loco") — canonical form per Step 1
#   type "pca"    → pca = TRUE
#   type "nopca"  → pca = FALSE
pca <- opt$type == "pca"

params <- list(
  population = opt$group,
  maf = as.numeric(opt$maf),
  nqtl = as.integer(opt$nqtl),
  effect = opt$effect,
  rep = as.integer(opt$rep),
  h2 = as.numeric(opt$h2),
  algorithm = opt$mode,   # "inbred" or "loco" — canonical form per Step 1
  pca = pca
)

ms_id   <- generate_marker_set_id(params$population, params$maf)
trait   <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
mapping <- generate_mapping_id(trait$hash, params$algorithm, params$pca)

log_msg(paste("Marker set ID:", ms_id$hash))
log_msg(paste("Trait ID:",      trait$hash))
log_msg(paste("Writing mapping:", mapping$hash))

# Write to Hive-partitioned Parquet
write_mapping_partitioned(
  df       = gwa_df,
  params   = params,
  ms_id    = ms_id,
  trait_id = trait,
  base_dir = opt$base_dir
)

write_mapping_metadata(
  params    = params,
  ms_id     = ms_id,
  trait_id  = trait,
  n_markers = nrow(gwa_df),
  base_dir  = opt$base_dir
)
log_msg(paste("Metadata sidecar written for mapping:", mapping$hash))

log_msg(paste("Mapping written successfully:", mapping$hash))
