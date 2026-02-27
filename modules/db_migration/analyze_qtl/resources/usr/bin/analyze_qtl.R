#!/usr/bin/env Rscript
# analyze_qtl.R - Query DB and detect QTL intervals
#
# Step 1 of the split DB-path analysis (formerly part of analyze_and_assess.R):
#   1. Query mapping data from Parquet database
#   2. Calculate significance threshold (BF or EIGEN)
#   3. Detect QTL intervals via analyze_mapping()
#   4. Write QTL regions TSV
#
# Does NOT use the .par file — it is staged as a pass-through for ASSESS_SIMS.
#
# Usage: analyze_qtl.R --group <group> --maf <maf> --nqtl <nqtl>
#            --effect <effect> --rep <rep> --h2 <h2> --mode <mode>
#            --type <type> --threshold <BF|EIGEN>
#            --base_dir <db_dir> --ci_size <int> --snp_grouping <int>
#            --alpha <num>

library(optparse)

option_list <- list(
  make_option("--group", type = "character", help = "Population/strain group identifier"),
  make_option("--maf", type = "character", help = "MAF threshold"),
  make_option("--nqtl", type = "integer", help = "Number of simulated QTLs"),
  make_option("--effect", type = "character", help = "Effect size distribution"),
  make_option("--rep", type = "integer", help = "Simulation replicate number"),
  make_option("--h2", type = "double", help = "Heritability"),
  make_option("--mode", type = "character", help = "GWA mode (inbred/loco)"),
  make_option("--type", type = "character", help = "PCA type (pca/nopca)"),
  make_option("--threshold", type = "character", help = "Threshold method (BF/EIGEN)"),
  make_option("--base_dir", type = "character", help = "Database output directory"),
  make_option("--ci_size", type = "integer", default = 150, help = "CI size in markers"),
  make_option("--snp_grouping", type = "integer", default = 1000, help = "SNP grouping distance"),
  make_option("--alpha", type = "double", default = 0.05, help = "Significance level")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c(
  "group", "maf", "nqtl", "effect", "rep", "h2", "mode", "type",
  "threshold", "base_dir"
)
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
source(file.path(r_source_dir, "queries.R"))
source(file.path(r_source_dir, "analysis.R"))

# Construct mapping params
algorithm <- if (opt$mode == "inbred") "LMM-EXACT-INBRED" else "LMM-EXACT-LOCO"

# set log pca flag based on type (pca/nopca)
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
  trait = paste(opt$nqtl, opt$rep, opt$h2, sep = "_"),
  threshold_method = toupper(opt$threshold)
)

ms_id      <- generate_marker_set_id(params$population, params$maf)
trait      <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
mapping    <- generate_mapping_id(trait$hash, params$algorithm, params$pca)
mapping_id <- mapping$hash
log_msg(paste("Analyzing mapping:", mapping_id, "with threshold:", opt$threshold))

# Step 1: Query mapping data from database
mapping_data <- query_for_threshold_analysis(mapping_id, opt$base_dir)

if (nrow(mapping_data) == 0) {
  stop(paste("No data found in database for mapping_id:", mapping_id))
}

log_msg(paste("Queried", nrow(mapping_data), "markers from database"))

# Step 2: Get threshold parameters
threshold_params <- get_threshold_params(opt$group, as.numeric(opt$maf), opt$alpha, opt$base_dir)

# Step 3: Calculate threshold
threshold_method <- toupper(opt$threshold)
if (threshold_method == "BF") {
  threshold <- calculate_threshold("BF", n_markers = threshold_params$n_markers, alpha = opt$alpha)
} else if (threshold_method == "EIGEN") {
  threshold <- calculate_threshold("EIGEN", n_independent = threshold_params$n_independent_tests, alpha = opt$alpha)
} else {
  threshold <- calculate_threshold(as.numeric(opt$threshold), alpha = opt$alpha)
}

# Step 4: Run analysis pipeline (flag significant, group intervals, define CIs)
processed <- analyze_mapping(
  df = mapping_data,
  threshold_value = threshold$threshold_value,
  threshold_method = threshold$threshold_method,
  ci_size = opt$ci_size,
  snp_grouping = opt$snp_grouping
)

# Step 5: Extract QTL regions
qtl_regions <- extract_qtl_regions(processed)

# Step 6: Write QTL regions TSV
output_file <- paste0(
  opt$nqtl, "_", opt$rep, "_", opt$h2, "_", opt$maf, "_",
  opt$effect, "_", opt$group, "_", opt$mode, "_", opt$type, "_",
  opt$threshold, "_qtl_regions.tsv"
)

write.table(
  qtl_regions,
  file = output_file,
  sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE
)

log_msg(paste("QTL regions written:", output_file, "(", nrow(qtl_regions), "regions)"))
