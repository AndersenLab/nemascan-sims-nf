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
  make_option("--cv_maf_effective", type = "double", help = "Effective CV MAF threshold"),
  make_option("--cv_ld", type = "double", help = "CV LD pruning threshold"),
  make_option("--ci_size", type = "integer", default = 150, help = "CI size in markers"),
  make_option("--snp_grouping", type = "integer", default = 1000, help = "SNP grouping distance"),
  make_option("--alpha", type = "double", default = 0.05, help = "Significance level")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Validate required args
required <- c("group", "maf", "nqtl", "effect", "rep", "h2", "mode", "type",
               "threshold", "par_file", "qtl_regions", "base_dir",
               "cv_maf_effective", "cv_ld")
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
pca <- opt$type == "pca"

params <- list(
  population = opt$group,
  maf = as.numeric(opt$maf),
  nqtl = as.integer(opt$nqtl),
  effect = opt$effect,
  rep = as.integer(opt$rep),
  h2 = as.numeric(opt$h2),
  algorithm = opt$mode,   # "inbred" or "loco" — canonical form per Step 1
  pca = pca,
  trait = paste(opt$nqtl, opt$rep, opt$h2, sep = "_"),
  threshold_method = toupper(opt$threshold),
  mode = opt$mode,
  type = opt$type,
  alpha = opt$alpha,
  ci_size = opt$ci_size,
  snp_grouping = opt$snp_grouping
)

ms_meta <- read_marker_set_metadata(params$population, as.numeric(params$maf), opt$base_dir)
if (is.null(ms_meta)) {
  stop(paste0(
    "Marker set metadata not found for population='", params$population,
    "', maf=", params$maf, " in ", opt$base_dir,
    ". Ensure DB_MIGRATION_WRITE_MARKER_SET completed before resuming. ",
    "If VCF/species parameters changed, a full re-run (not -resume) is required."
  ))
}
if (is.na(as.numeric(ms_meta$ms_ld))) {
  stop("ms_ld field is NA in marker set metadata — DB may be corrupt")
}
sim_params <- list(
  population       = params$population,
  maf              = as.numeric(params$maf),
  species          = ms_meta$species,
  vcf_release_id   = ms_meta$vcf_release_id,
  ms_ld            = as.numeric(ms_meta$ms_ld),
  nqtl             = params$nqtl,
  effect           = params$effect,
  rep              = params$rep,
  h2               = params$h2,
  cv_maf_effective = as.numeric(opt$cv_maf_effective),
  cv_ld            = as.numeric(opt$cv_ld)
)
ids        <- build_ids_from_params(sim_params, mode = opt$mode, pca = opt$type == "pca")
mapping_id <- ids$mapping_id$hash
trait_hash <- ids$trait_id$hash
log_msg(paste("Assessing mapping:", mapping_id))

# Step 1: Read QTL regions from analyze_qtl.R output
qtl_regions <- data.table::fread(opt$qtl_regions, header = TRUE) %>% as.data.frame()
log_msg(paste("Read", nrow(qtl_regions), "QTL regions"))

# Step 2: Load causal variants from .par file
causal_variants <- load_causal_variants(opt$par_file)

# Step 3a: Compute variance explained from DB-stored genotype and phenotype data.
# var.exp is a per-trait property; ASSESS_SIMS runs once per (trait × mode × type × threshold).
# All invocations for a trait compute identical values — redundancy is accepted.
# var.exp is scale-invariant (SS ratio); read is the full genotype Parquet (~400–600 MB at
# production scale, <8 GB task allocation). Accepted redundancy across mode/type invocations.
genotype_matrix <- tryCatch(
  read_genotype_matrix(
    params$population, as.numeric(params$maf),
    ms_meta$species, ms_meta$vcf_release_id, as.numeric(ms_meta$ms_ld),
    opt$base_dir
  ),
  error = function(e) stop(paste0(
    "Failed to read genotype matrix for population='", params$population,
    "', maf=", params$maf, ": ", conditionMessage(e)))
)
phenotype_data <- tryCatch(
  read_phenotype_data(trait_hash, opt$base_dir),
  error = function(e) stop(paste0(
    "Failed to read phenotype data for trait='", trait_hash, "': ",
    conditionMessage(e)))
)

# Read per-trait causal genotypes (covers non-marker positions when cv_maf < ms_maf).
causal_geno <- tryCatch(
  read_causal_genotypes(trait_hash, opt$base_dir),
  error = function(e) {
    warning("Could not read causal genotypes for trait='", trait_hash,
            "': ", conditionMessage(e))
    NULL
  }
)

# Detect non-marker causal positions and warn so users can distinguish sources.
non_marker_count <- causal_variants %>%
  dplyr::anti_join(
    genotype_matrix %>%
      dplyr::mutate(CHROM = as.character(CHROM), POS = as.integer(POS)) %>%
      dplyr::select(CHROM, POS) %>% dplyr::distinct(),
    by = c("CHROM", "POS")
  ) %>% nrow()

if (non_marker_count > 0) {
  message(non_marker_count,
          " non-marker causal variant(s): var.exp computed from per-trait causal genotypes")
}

# Merge genotype sources: causal genotypes are authoritative for non-marker positions;
# marker genotype matrix supplies the rest.
merged_geno <- if (!is.null(causal_geno) && nrow(causal_geno) > 0) {
  cg_cols <- causal_geno %>%
    dplyr::mutate(CHROM = as.character(CHROM), POS = as.integer(POS)) %>%
    dplyr::select(CHROM, POS, strain, allele)
  marker_only <- genotype_matrix %>%
    dplyr::mutate(CHROM = as.character(CHROM), POS = as.integer(POS)) %>%
    dplyr::select(CHROM, POS, strain, allele) %>%
    dplyr::anti_join(
      cg_cols %>% dplyr::select(CHROM, POS) %>% dplyr::distinct(),
      by = c("CHROM", "POS")
    )
  dplyr::bind_rows(marker_only, cg_cols)
} else {
  genotype_matrix %>%
    dplyr::mutate(CHROM = as.character(CHROM), POS = as.integer(POS)) %>%
    dplyr::select(CHROM, POS, strain, allele)
}

var_exp_df <- compute_var_exp_anova(merged_geno, phenotype_data, causal_variants)

# Augment causal_variants with Simulated.QTL.VarExp before compile_full_assessment().
# build_assessment_union() picks it up via score_causal_markers() join (any_of select).
causal_variants <- dplyr::left_join(causal_variants, var_exp_df, by = "QTL")

# interval.var.exp stays NA_real_ by design — GCTA's LMM-adjusted r² at the peak marker
# cannot be recovered from the DB (raw per-marker GWA r² was not stored in Phase 5).

# Step 3: Query mapping data for causal marker scores + peak info
mapping_data <- query_for_threshold_analysis(mapping_id, opt$base_dir)
if (nrow(mapping_data) == 0) {
  stop(paste("No data found in database for mapping_id:", mapping_id))
}

# Compute threshold (needed for analyze_mapping to add peak annotations)
# BF uses nrow(mapping_data): actual GCTA output count, matching legacy Get_GCTA_Intervals.R.
# See analyze_qtl.R Step 3 comment for the full explanation.
threshold_params <- get_threshold_params(opt$group, as.numeric(opt$maf), opt$alpha, opt$base_dir)
threshold_method <- toupper(opt$threshold)
if (threshold_method == "BF") {
  threshold <- calculate_threshold("BF", n_markers = nrow(mapping_data), alpha = opt$alpha)
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

# Ensure Simulated.QTL.VarExp is present even if all causal variants were filtered
# (e.g. no mapping peak matched). format_assessment_tsv() handles the NA default,
# but guard here prevents downstream confusion if assessment has no matching rows.
if (!"Simulated.QTL.VarExp" %in% names(assessment)) {
  assessment$Simulated.QTL.VarExp <- NA_real_
}

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
