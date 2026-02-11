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

# Construct mapping params from Nextflow channel metadata
# Maps GWA mode/type to database algorithm/pca fields:
#   mode "inbred" → algorithm "LMM-EXACT-INBRED"
#   mode "loco"   → algorithm "LMM-EXACT-LOCO"
#   type "pca"    → pca = TRUE
#   type "nopca"  → pca = FALSE
algorithm <- if (opt$mode == "inbred") "LMM-EXACT-INBRED" else "LMM-EXACT-LOCO"
pca <- opt$type == "pca"

params <- list(
  population = opt$group,
  maf = as.numeric(opt$maf),
  nqtl = as.integer(opt$nqtl),
  effect = opt$effect,
  rep = as.integer(opt$rep),
  h2 = as.numeric(opt$h2),
  algorithm = algorithm,
  pca = pca,
  trait = paste(opt$nqtl, opt$rep, opt$h2, sep = "_")
)

mapping_id <- generate_mapping_id(params)
log_msg(paste("Writing mapping:", mapping_id))

# Write to Hive-partitioned Parquet
# var.exp is NA for inline-path databases (no genotype matrix available)
write_mapping_partitioned(
  df = gwa_df,
  params = params,
  base_dir = opt$base_dir
)

log_msg(paste("Mapping written successfully:", mapping_id))
