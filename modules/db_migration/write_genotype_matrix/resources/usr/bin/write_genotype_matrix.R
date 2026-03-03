#!/usr/bin/env Rscript

# Write genotype matrix to database in long-format Parquet
#
# Input: wide-format TSV from BCFTOOLS_CREATE_GENOTYPE_MATRIX
# Output: long-format Parquet in {db_dir}/markers/genotypes/{marker_set_id}_genotypes.parquet

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--group",            type = "character"),
  make_option("--maf",              type = "character"),
  make_option("--genotype_matrix",  type = "character"),
  make_option("--species",          type = "character"),
  make_option("--vcf_release_id",   type = "character"),
  make_option("--ms_ld",            type = "double"),
  make_option("--base_dir",         type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Source R functions
r_source_dir <- Sys.getenv("R_SOURCE_DIR")
source(file.path(r_source_dir, "utils.R"))
source(file.path(r_source_dir, "io.R"))
source(file.path(r_source_dir, "database.R"))

# Initialize database directories
init_database(opt$base_dir)

# Write genotype matrix
# ⚠ cv_maf/cv_ld are NOT included in marker_set_id — see README Marker Set ID section
message("Writing genotype matrix: ", opt$group, "_", opt$maf)
write_genotype_matrix(
  genotype_tsv   = opt$genotype_matrix,
  population     = opt$group,
  maf            = opt$maf,
  species        = opt$species,
  vcf_release_id = opt$vcf_release_id,
  ms_ld          = as.numeric(opt$ms_ld),
  base_dir       = opt$base_dir
)
message("Done: genotype matrix written")
