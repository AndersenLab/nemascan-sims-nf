#!/usr/bin/env Rscript
# write_gwa_to_db.R - Write raw GWA output to Parquet database
#
# Reads a GCTA GWA output file (.fastGWA or .mlma), reads species/vcf_release_id/ms_ld
# from marker_set_metadata.parquet (written by DB_MIGRATION_WRITE_MARKER_SET), computes
# IDs via build_ids_from_params(), and writes to Hive-partitioned Parquet.
#
# Usage: write_gwa_to_db.R --group <group> --maf <maf> --nqtl <nqtl>
#            --effect <effect> --rep <rep> --h2 <h2> --mode <mode>
#            --type <type> --cv_maf_effective <num> --cv_ld <num>
#            --gwa_file <file> --base_dir <db_dir>

library(optparse)

option_list <- list(
  make_option("--group",          type = "character", help = "Population/strain group identifier"),
  make_option("--maf",            type = "character", help = "MAF threshold"),
  make_option("--nqtl",           type = "integer",   help = "Number of simulated QTLs"),
  make_option("--effect",         type = "character", help = "Effect size distribution (gamma/uniform)"),
  make_option("--rep",            type = "integer",   help = "Simulation replicate number"),
  make_option("--h2",             type = "double",    help = "Heritability"),
  make_option("--mode",           type = "character", help = "GWA mode (inbred/loco)"),
  make_option("--type",           type = "character", help = "PCA type (pca/nopca)"),
  make_option("--cv_maf_effective", type = "double",  help = "Effective MAF threshold for the CV pool"),
  make_option("--cv_ld",            type = "double",  help = "CV LD pruning threshold"),
  make_option("--gwa_file",         type = "character", help = "Path to GWA output file"),
  make_option("--base_dir",         type = "character", help = "Database output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c("group", "maf", "nqtl", "effect", "rep", "h2", "mode", "type",
              "cv_maf_effective", "cv_ld", "gwa_file", "base_dir")
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
n_before <- nrow(gwa_df)
gwa_df <- gwa_df %>% dplyr::distinct(marker, .keep_all = TRUE)
n_after <- nrow(gwa_df)
if (n_before != n_after) {
  log_msg(paste("Deduplicated:", n_before, "->", n_after, "unique markers (removed",
                n_before - n_after, "duplicates)"))
}

pca <- opt$type == "pca"

# Read species/vcf_release_id/ms_ld from marker_set_metadata.parquet.
# DB_MIGRATION_WRITE_MARKER_SET is guaranteed to have completed before this task
# starts (via ch_marker_barrier in main.nf), so the metadata file always exists here.
ms_meta <- read_marker_set_metadata(opt$group, as.numeric(opt$maf), opt$base_dir)
if (is.null(ms_meta)) {
  stop(paste0(
    "Marker set metadata not found for population='", opt$group,
    "', maf=", opt$maf, " in ", opt$base_dir,
    ". Ensure DB_MIGRATION_WRITE_MARKER_SET completed before resuming. ",
    "If VCF/species parameters changed, a full re-run (not -resume) is required."
  ))
}
if (is.na(as.numeric(ms_meta$ms_ld))) {
  stop("ms_ld field is NA in marker set metadata — DB may be corrupt")
}

# Compute all IDs via canonical wrapper
sim_params <- list(
  population       = opt$group,
  maf              = as.numeric(opt$maf),
  species          = ms_meta$species,
  vcf_release_id   = ms_meta$vcf_release_id,
  ms_ld            = as.numeric(ms_meta$ms_ld),
  nqtl             = opt$nqtl,
  effect           = opt$effect,
  rep              = opt$rep,
  h2               = opt$h2,
  cv_maf_effective = as.numeric(opt$cv_maf_effective),
  cv_ld            = as.numeric(opt$cv_ld)
)
ids     <- build_ids_from_params(sim_params, mode = opt$mode, pca = pca)
ms_id   <- ids$ms_id
trait   <- ids$trait_id
mapping <- ids$mapping_id

log_msg(paste("Marker set ID:", ms_id$hash))
log_msg(paste("Trait ID:",      trait$hash))
log_msg(paste("Writing mapping:", mapping$hash))

# params list for write functions (algorithm/pca needed; species/vcf/ms_ld are not)
params <- list(
  population = opt$group,
  maf        = as.numeric(opt$maf),
  nqtl       = as.integer(opt$nqtl),
  effect     = opt$effect,
  rep        = as.integer(opt$rep),
  h2         = as.numeric(opt$h2),
  algorithm  = opt$mode,
  pca        = pca
)

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
