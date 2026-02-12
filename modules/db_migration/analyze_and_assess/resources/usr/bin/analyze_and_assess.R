#!/usr/bin/env Rscript
# analyze_and_assess.R - Query DB, detect QTL intervals, assess against simulated truth
#
# Implements the DB-path equivalent of R_GET_GCTA_INTERVALS + R_ASSESS_SIMS:
#   1. Query mapping data from Parquet database
#   2. Calculate significance threshold (BF or EIGEN)
#   3. Detect QTL intervals via analyze_mapping()
#   4. Load causal variants from .par file
#   5. Assess detection (simulated vs detected)
#   6. Write assessment TSV (column-compatible with existing output)
#
# Usage: analyze_and_assess.R --group <group> --maf <maf> --nqtl <nqtl>
#            --effect <effect> --rep <rep> --h2 <h2> --mode <mode>
#            --type <type> --threshold <BF|EIGEN> --par_file <file>
#            --base_dir <db_dir> --ci_size <int> --snp_grouping <int>

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
  make_option("--par_file", type = "character", help = "Path to .par file with causal variants"),
  make_option("--base_dir", type = "character", help = "Database output directory"),
  make_option("--ci_size", type = "integer", default = 150, help = "CI size in markers"),
  make_option("--snp_grouping", type = "integer", default = 1000, help = "SNP grouping distance")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c("group", "maf", "nqtl", "effect", "rep", "h2", "mode", "type",
               "threshold", "par_file", "base_dir")
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
source(file.path(r_source_dir, "assessment.R"))

# Construct mapping params (same as write_gwa_to_db.R)
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
  trait = paste(opt$nqtl, opt$rep, opt$h2, sep = "_"),
  threshold_method = toupper(opt$threshold)
)

mapping_id <- generate_mapping_id(params)
log_msg(paste("Analyzing mapping:", mapping_id, "with threshold:", opt$threshold))

# Step 1: Query mapping data from database
mapping_data <- query_for_threshold_analysis(mapping_id, opt$base_dir)

if (nrow(mapping_data) == 0) {
  stop(paste("No data found in database for mapping_id:", mapping_id))
}

log_msg(paste("Queried", nrow(mapping_data), "markers from database"))

# Step 2: Get threshold parameters
threshold_params <- get_threshold_params(opt$group, as.numeric(opt$maf), 0.05, opt$base_dir)

# Step 3: Calculate threshold
threshold_method <- toupper(opt$threshold)
if (threshold_method == "BF") {
  threshold <- calculate_threshold("BF", n_markers = threshold_params$n_markers)
} else if (threshold_method == "EIGEN") {
  threshold <- calculate_threshold("EIGEN", n_independent = threshold_params$n_independent_tests)
} else {
  threshold <- calculate_threshold(as.numeric(opt$threshold))
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

# Step 6: Load causal variants from .par file
causal_variants <- load_causal_variants(opt$par_file)

# Step 7: Compile full assessment
assessment <- compile_full_assessment(
  mapping_data = processed,
  qtl_regions = qtl_regions,
  causal_variants = causal_variants,
  mapping_params = params
)

# Step 8: Format and write output
if (nrow(assessment) > 0) {
  formatted <- format_assessment_tsv(assessment)
} else {
  # Even with no QTLs, write an empty file with correct structure
  formatted <- assessment
}

output_file <- paste0(
  opt$nqtl, "_", opt$rep, "_", opt$h2, "_", opt$maf, "_",
  opt$effect, "_", opt$group, "_", opt$mode, "_", opt$type, "_",
  opt$threshold, "_db_assessment.tsv"
)

write.table(
  formatted,
  file = output_file,
  sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE
)

log_msg(paste("Assessment written:", output_file, "(", nrow(formatted), "rows)"))
