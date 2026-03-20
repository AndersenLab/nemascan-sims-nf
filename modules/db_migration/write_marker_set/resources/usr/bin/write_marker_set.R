#!/usr/bin/env Rscript
# write_marker_set.R - Write marker set to Parquet database from .bim file
#
# Reads a PLINK .bim file and EIGEN independent tests file, then writes the
# marker set to the database. This is the inline-path version that reads .bim
# directly (not processed mapping TSVs).
#
# Usage: write_marker_set.R --group <group> --maf <maf> --bim <bim_file>
#                            --eigen <eigen_file> --species <species>
#                            --vcf_release_id <vcf_release_id> --ms_ld <ms_ld>
#                            --base_dir <db_dir> --strainfile_path <path>
#                            --strains <comma_separated_strains>

library(optparse)

option_list <- list(
  make_option("--group",          type = "character", help = "Population/strain group identifier"),
  make_option("--maf",            type = "character", help = "MAF threshold"),
  make_option("--bim",            type = "character", help = "Path to PLINK .bim file"),
  make_option("--eigen",          type = "character", help = "Path to EIGEN independent tests file"),
  make_option("--species",        type = "character", help = "Species (c_elegans, c_briggsae, c_tropicalis)"),
  make_option("--vcf_release_id", type = "character", help = "VCF release date (e.g. 20220216)"),
  make_option("--ms_ld",          type = "double",    help = "LD R² threshold for marker SNP selection"),
  make_option("--base_dir",        type = "character", help = "Database output directory"),
  make_option("--strainfile_path", type = "character", help = "Path to the strainfile (for SHA-256 provenance hash)"),
  make_option("--strains",         type = "character", help = "Comma-separated strain names for this group")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c("group", "maf", "bim", "eigen", "species", "vcf_release_id", "ms_ld", "base_dir",
              "strainfile_path", "strains")
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

# Read .bim file for marker data (AF1 will be NA — only available from GWA output)
log_msg(paste("Reading .bim file:", opt$bim))
markers_df <- read_bim_file(opt$bim)

# Read EIGEN independent tests value from single-line text file
log_msg(paste("Reading EIGEN file:", opt$eigen))
n_indep_tests <- read_eigen_file(opt$eigen)

# Compute SHA-256 of the strainfile — consistent with other hash generation in R/database.R
# digest::digest() with file=TRUE reads raw bytes, same behaviour as the pipeline's Java MessageDigest
log_msg(paste("Computing strainfile SHA-256:", opt$strainfile_path))
strainfile_hash <- digest::digest(opt$strainfile_path, algo = "sha256", file = TRUE)

# Write marker set to database
# ⚠ cv_maf/cv_ld are NOT included in marker_set_id — see README Marker Set ID section
log_msg(paste("Writing marker set:", opt$group, opt$maf))
write_marker_set(
  df = markers_df,
  population = opt$group,
  maf = as.numeric(opt$maf),
  species = opt$species,
  vcf_release_id = opt$vcf_release_id,
  ms_ld = as.numeric(opt$ms_ld),
  base_dir = opt$base_dir,
  overwrite = TRUE,
  n_independent_tests = n_indep_tests,
  eigen_source_file = basename(opt$eigen),
  strainfile_hash = strainfile_hash,
  strain_list = opt$strains
)

log_msg(paste("Marker set written successfully:", opt$group, opt$maf))
