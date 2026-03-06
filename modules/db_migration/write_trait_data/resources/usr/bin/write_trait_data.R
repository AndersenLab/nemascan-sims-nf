#!/usr/bin/env Rscript

# Write trait data (metadata, causal variants, phenotype) to database
#
# Trait ID is computed internally via generate_trait_id() — no cross-language
# hash computation. The Nextflow .multiMap{} block passes raw simulation
# parameters; this script computes the hash in R only.
#
# The phenotype stored is the pre-upscaled value from GCTA_SIMULATE_PHENOTYPES
# (before check_vp.py). This is intentional — ANOVA SS ratio is scale-invariant.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--group",      type = "character"),
  make_option("--maf",        type = "character"),
  make_option("--nqtl",       type = "integer"),
  make_option("--effect",     type = "character"),
  make_option("--rep",        type = "integer"),
  make_option("--h2",         type = "double"),
  make_option("--pheno_file", type = "character"),
  make_option("--par_file",   type = "character"),
  make_option("--base_dir",   type = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Source R functions
r_source_dir <- Sys.getenv("R_SOURCE_DIR")
source(file.path(r_source_dir, "utils.R"))
source(file.path(r_source_dir, "io.R"))
source(file.path(r_source_dir, "database.R"))
source(file.path(r_source_dir, "assessment.R"))  # for load_causal_variants()

# Initialize database directories
init_database(opt$base_dir)

# --- Runtime validation of .merge() synchronization (M4) ---
# Verify .phen filename matches expected parameters (catches .merge() misalignment)
expected_prefix <- paste(opt$nqtl, opt$rep, opt$h2, opt$maf,
                         opt$effect, opt$group, sep = "_")
actual_basename <- basename(opt$pheno_file)
if (!grepl(expected_prefix, actual_basename)) {
  stop(paste("Phenotype file mismatch! Expected prefix:", expected_prefix,
             "but got:", actual_basename))
}

# Compute marker set ID (parent of trait)
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
ms_id <- generate_marker_set_id(
  opt$group, as.numeric(opt$maf),
  ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld)
)

# Compute trait ID (child of marker set)
trait    <- generate_trait_id(ms_id$hash, opt$nqtl, opt$effect, opt$rep, opt$h2)
trait_id <- trait$hash

message("Writing trait data: trait_id=", trait_id,
        " group=", opt$group, " nqtl=", opt$nqtl,
        " rep=", opt$rep, " h2=", opt$h2)

# Write trait metadata
write_trait_metadata(
  trait_id          = trait_id,
  trait_hash_string = trait$hash_string,
  marker_set_id     = ms_id$hash,
  nqtl       = opt$nqtl,
  rep        = opt$rep,
  sim_seed   = as.integer(opt$rep),
  h2         = opt$h2,
  maf        = opt$maf,
  effect     = opt$effect,
  population = opt$group,
  base_dir   = opt$base_dir
)

# Write causal variants
write_causal_variants(
  par_file = opt$par_file,
  trait_id = trait_id,
  base_dir = opt$base_dir
)

# Write phenotype data (pre-upscaled — scale-invariant for ANOVA SS)
write_phenotype_data(
  pheno_file = opt$pheno_file,
  trait_id   = trait_id,
  base_dir   = opt$base_dir
)

message("Done: trait data written for trait_id=", trait_id)
