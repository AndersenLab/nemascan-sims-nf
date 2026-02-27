#!/usr/bin/env Rscript
# assess_sims.R - Assess QTL detection against simulated truth
#
# Step 2 of the split DB-path analysis (formerly part of analyze_and_assess.R):
#   1. Read QTL regions TSV (from analyze_qtl.R)
#   2. Load causal variants from .par file
#   3. Query mapping data for causal marker scores
#   4. Compile full assessment (simulated vs detected)
#   5. Write assessment TSV
#
# Usage: assess_sims.R --group <group> --maf <maf> --nqtl <nqtl>
#            --effect <effect> --rep <rep> --h2 <h2> --mode <mode>
#            --type <type> --threshold <BF|EIGEN> --par_file <file>
#            --qtl_regions <file> --base_dir <db_dir>
#            --ci_size <int> --snp_grouping <int> --alpha <num>

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
  make_option("--qtl_regions", type = "character", help = "Path to QTL regions TSV from analyze_qtl.R"),
  make_option("--base_dir", type = "character", help = "Database output directory"),
  make_option("--ci_size", type = "integer", default = 150, help = "CI size in markers"),
  make_option("--snp_grouping", type = "integer", default = 1000, help = "SNP grouping distance"),
  make_option("--alpha", type = "double", default = 0.05, help = "Significance level")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c("group", "maf", "nqtl", "effect", "rep", "h2", "mode", "type",
               "threshold", "par_file", "qtl_regions", "base_dir")
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

# Construct mapping params
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
  threshold_method = toupper(opt$threshold),
  mode = opt$mode,
  type = opt$type,
  alpha = opt$alpha,
  ci_size = opt$ci_size,
  snp_grouping = opt$snp_grouping
)

ms_id      <- generate_marker_set_id(params$population, params$maf)
trait      <- generate_trait_id(ms_id$hash, params$nqtl, params$effect, params$rep, params$h2)
mapping    <- generate_mapping_id(trait$hash, params$algorithm, params$pca)
mapping_id <- mapping$hash
log_msg(paste("Assessing mapping:", mapping_id))

# Step 1: Read QTL regions from analyze_qtl.R output
qtl_regions <- data.table::fread(opt$qtl_regions, header = TRUE) %>% as.data.frame()
log_msg(paste("Read", nrow(qtl_regions), "QTL regions"))

# Step 2: Load causal variants from .par file
causal_variants <- load_causal_variants(opt$par_file)

# Step 3: Query mapping data for causal marker scores + peak info
mapping_data <- query_for_threshold_analysis(mapping_id, opt$base_dir)
if (nrow(mapping_data) == 0) {
  stop(paste("No data found in database for mapping_id:", mapping_id))
}

# Compute threshold (needed for analyze_mapping to add peak annotations)
threshold_params <- get_threshold_params(opt$group, as.numeric(opt$maf), opt$alpha, opt$base_dir)
threshold_method <- toupper(opt$threshold)
if (threshold_method == "BF") {
  threshold <- calculate_threshold("BF", n_markers = threshold_params$n_markers, alpha = opt$alpha)
} else if (threshold_method == "EIGEN") {
  threshold <- calculate_threshold("EIGEN", n_independent = threshold_params$n_independent_tests, alpha = opt$alpha)
} else {
  threshold <- calculate_threshold(as.numeric(opt$threshold), alpha = opt$alpha)
}

# Re-run analysis to get the processed mapping data with peak annotations
processed <- analyze_mapping(
  df = mapping_data,
  threshold_value = threshold$threshold_value,
  threshold_method = threshold$threshold_method,
  ci_size = opt$ci_size,
  snp_grouping = opt$snp_grouping
)

# Step 4: Compile full assessment
assessment <- compile_full_assessment(
  mapping_data = processed,
  qtl_regions = qtl_regions,
  causal_variants = causal_variants,
  mapping_params = params
)

# Step 5: Format and write output
if (nrow(assessment) > 0) {
  formatted <- format_assessment_tsv(assessment)
} else {
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
