# assessment.R - Simulation assessment functions for DB-path QTL analysis
#
# Implements the assessment logic from bin/Assess_Sims.R without
# GenomicRanges or genotype matrix dependencies:
#   - load_causal_variants(): read .par file
#   - score_causal_markers(): match causal variants to mapping data
#   - find_peak_causal_overlaps(): match peaks to causal variants
#   - build_assessment_union(): combine simulated + detected into unified df
#   - compile_full_assessment(): orchestrate the above helpers
#   - format_assessment_tsv(): format output with standardized column schema
#
# Dependencies: analysis.R (analyze_mapping, extract_qtl_regions, safe_log10p)
#               utils.R (log_msg)

library(dplyr)
library(tidyr)
library(glue)


# =============================================================================
# Load Causal Variants from .par File
# =============================================================================

#' Read simulation parameter (.par) file and extract causal variant positions
#'
#' The .par file has columns: QTL, RefAllele, Frequency, Effect
#' where QTL is in CHROM:POS format.
#'
#' @param par_file Path to .par file
#' @return Dataframe with columns: QTL, CHROM, POS, RefAllele, Frequency, Effect
load_causal_variants <- function(par_file) {
  if (!file.exists(par_file)) {
    stop(glue::glue("Par file not found: {par_file}"))
  }

  effects <- data.table::fread(par_file, header = TRUE) %>%
    as.data.frame() %>%
    tidyr::separate(QTL, c("CHROM", "POS"), sep = ":", remove = FALSE) %>%
    dplyr::mutate(
      CHROM = as.character(CHROM),
      POS = as.integer(POS)
    )

  log_msg(glue::glue("Loaded {nrow(effects)} causal variants from {basename(par_file)}"))

  effects
}


# =============================================================================
# Assessment Helper Functions
# =============================================================================

#' Score causal variant markers in the mapping data
#'
#' Joins causal variant markers against mapping data by QTL == marker to obtain
#' log10p and significance status for each causal variant.
#'
#' Causal variants are sampled post-MAF-filter in PYTHON_SIMULATE_EFFECTS_GLOBAL,
#' so they should always be present in GWA output. The filter(!is.na(log10p))
#' is defensive; both paths have identical behavior.
#'
#' @param mapping_data Processed mapping dataframe with marker, P, significant columns
#' @param causal_variants Dataframe from load_causal_variants()
#' @return Dataframe with causal variant info plus log10p and significant columns
score_causal_markers <- function(mapping_data, causal_variants) {
  # Compute log10p on mapping data if not present
  if (!"log10p" %in% names(mapping_data)) {
    mapping_data <- mapping_data %>%
      dplyr::mutate(log10p = safe_log10p(P))
  }

  # Get scores for causal variant markers in the mapping data
  causal_marker_scores <- mapping_data %>%
    dplyr::mutate(QTL = marker) %>%
    dplyr::filter(QTL %in% causal_variants$QTL) %>%
    dplyr::select(QTL, log10p, significant) %>%
    dplyr::filter(!duplicated(QTL))

  # Add log10p and significance to causal variants
  effects_scores <- causal_variants %>%
    dplyr::left_join(causal_marker_scores, by = "QTL") %>%
    dplyr::filter(!is.na(log10p))

  effects_scores
}


#' Find overlaps between detected peaks and causal variants
#'
#' For each detected peak, finds causal variants with POS in [startPOS, endPOS]
#' on the same chromosome. Peaks with no overlapping causal variant are
#' included as false positives.
#'
#' @param peak_info Dataframe with one row per peak: CHROM, marker, POS, AF1, BETA,
#'   log10p, startPOS, peakPOS, endPOS, peak_id, interval_size, detected.peak
#' @param effects_scores Dataframe from score_causal_markers(): QTL, CHROM, POS,
#'   RefAllele, Frequency, Effect, log10p, significant
#' @return Dataframe with overlap info for detected peaks + false positive rows
find_peak_causal_overlaps <- function(peak_info, effects_scores) {
  if (nrow(peak_info) == 0 || nrow(effects_scores) == 0) {
    return(data.frame(QTL = character(), stringsAsFactors = FALSE))
  }

  has_var_exp <- "var.exp" %in% names(peak_info)

  # Prepare peak-level interval info with renamed columns to avoid collisions
  peaks <- peak_info %>%
    dplyr::transmute(
      peak_id, peak_CHROM = CHROM, peak_POS = POS,
      startPOS, peakPOS, endPOS,
      detected.peak, peak_BETA = BETA,
      interval.log10p = log10p,
      interval.var.exp = if (has_var_exp) var.exp else NA_real_,
      interval.Frequency = AF1,
      interval_size
    )

  # Cross-join peaks with causal variants on same chromosome, filter by position
  matched <- peaks %>%
    dplyr::inner_join(
      effects_scores,
      by = c("peak_CHROM" = "CHROM"),
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(POS >= startPOS, POS <= endPOS) %>%
    dplyr::mutate(CHROM = peak_CHROM, BETA = peak_BETA) %>%
    dplyr::select(-peak_CHROM, -peak_POS, -peak_BETA)

  # False-positive peaks: no causal variant within interval
  false_positives <- peaks %>%
    dplyr::filter(!peak_id %in% matched$peak_id) %>%
    dplyr::transmute(
      QTL = detected.peak,
      CHROM = peak_CHROM,
      POS = peak_POS,
      RefAllele = NA_character_,
      Frequency = NA_real_,
      Effect = NA_real_,
      log10p = NA_real_,
      significant = NA_integer_,
      startPOS, peakPOS, endPOS,
      detected.peak, BETA = peak_BETA,
      interval.log10p, interval.var.exp, interval.Frequency,
      peak_id, interval_size
    )

  dplyr::bind_rows(matched, false_positives)
}


#' Build the Simulated/Detected union dataframe
#'
#' Combines causal variant scores and peak overlap info into a unified assessment
#' dataframe with simulation metadata and classification flags.
#'
#' @param effects_scores Dataframe from score_causal_markers()
#' @param overlap_df Dataframe from find_peak_causal_overlaps()
#' @param mapping_params Named list with: nqtl, rep, h2, maf, effect, population,
#'   algorithm, pca, threshold_method. Optional: mode, type, alpha, ci_size, snp_grouping
#' @return Assessment dataframe with Simulated/Detected classification
build_assessment_union <- function(effects_scores, overlap_df, mapping_params) {
  all_qtl_ids <- unique(c(effects_scores$QTL, overlap_df$QTL))

  if (length(all_qtl_ids) == 0) {
    log_msg("No QTLs to assess (no causal variants matched and no intervals detected)")
    return(data.frame())
  }

  all_qtl <- data.frame(QTL = all_qtl_ids, stringsAsFactors = FALSE) %>%
    dplyr::filter(!duplicated(QTL)) %>%
    dplyr::mutate(
      QTL = as.character(QTL),
      Simulated = QTL %in% effects_scores$QTL,
      Detected = QTL %in% overlap_df$QTL
    )

  # Join causal variant info (CHROM, POS, RefAllele, Frequency, Effect, log10p, significant)
  if (nrow(effects_scores) > 0) {
    all_qtl <- all_qtl %>%
      dplyr::left_join(
        effects_scores %>%
          dplyr::select(QTL, CHROM, POS, RefAllele, Frequency, Effect,
                        log10p, significant) %>%
          dplyr::filter(!duplicated(QTL)),
        by = "QTL"
      )
  }

  # Join overlap/interval info
  if (nrow(overlap_df) > 0) {
    overlap_join_cols <- intersect(
      c("QTL", "startPOS", "peakPOS", "endPOS", "detected.peak", "BETA",
        "interval.log10p", "interval.var.exp", "interval.Frequency",
        "peak_id", "interval_size"),
      names(overlap_df)
    )
    overlap_join <- overlap_df %>%
      dplyr::select(dplyr::all_of(overlap_join_cols)) %>%
      dplyr::filter(!duplicated(QTL))

    # For false-positive peaks, also get CHROM/POS from overlap_df
    fp_coords <- overlap_df %>%
      dplyr::filter(!QTL %in% effects_scores$QTL) %>%
      dplyr::select(QTL, CHROM_ov = CHROM, POS_ov = POS) %>%
      dplyr::filter(!duplicated(QTL))

    all_qtl <- all_qtl %>%
      dplyr::left_join(overlap_join, by = "QTL") %>%
      dplyr::left_join(fp_coords, by = "QTL") %>%
      dplyr::mutate(
        CHROM = dplyr::coalesce(CHROM, CHROM_ov),
        POS = dplyr::coalesce(as.integer(POS), as.integer(POS_ov))
      ) %>%
      dplyr::select(-dplyr::any_of(c("CHROM_ov", "POS_ov")))
  }

  # Build algorithm_id and add metadata
  algorithm_id <- paste0(
    mapping_params$algorithm,
    if (isTRUE(mapping_params$pca)) "_PCA" else "_noPCA",
    "_", mapping_params$threshold_method
  )

  all_qtl <- all_qtl %>%
    dplyr::mutate(
      Simulated.QTL.VarExp = NA_real_,
      top.hit = if ("detected.peak" %in% names(.)) QTL == detected.peak else NA,
      nQTL = as.character(mapping_params$nqtl),
      simREP = as.character(mapping_params$rep),
      h2 = as.character(mapping_params$h2),
      maf = as.character(mapping_params$maf),
      effect_distribution = as.character(mapping_params$effect),
      strain_set_id = as.character(mapping_params$population),
      mode = if (!is.null(mapping_params$mode)) as.character(mapping_params$mode) else NA_character_,
      type = if (!is.null(mapping_params$type)) as.character(mapping_params$type) else NA_character_,
      threshold = as.character(mapping_params$threshold_method),
      algorithm_id = algorithm_id,
      alpha = if (!is.null(mapping_params$alpha)) as.numeric(mapping_params$alpha) else NA_real_,
      ci_size = if (!is.null(mapping_params$ci_size)) as.integer(mapping_params$ci_size) else NA_integer_,
      snp_grouping = if (!is.null(mapping_params$snp_grouping)) as.integer(mapping_params$snp_grouping) else NA_integer_,
      Simulated = factor(Simulated, levels = c("TRUE", "FALSE")),
      Detected = factor(Detected, levels = c("TRUE", "FALSE"))
    )

  n_simulated <- sum(all_qtl$Simulated == "TRUE")
  n_detected <- sum(all_qtl$Detected == "TRUE")
  n_total <- nrow(all_qtl)
  log_msg(glue::glue(
    "Assessment complete: {n_total} QTLs ({n_simulated} simulated, {n_detected} detected)"
  ))

  all_qtl
}


# =============================================================================
# Full Assessment Compilation
# =============================================================================

#' Compile full assessment combining simulated and detected QTLs
#'
#' Orchestrates the assessment pipeline: score causal markers -> build peak info ->
#' find overlaps -> build union dataframe.
#'
#' Produces a dataframe containing:
#'   - True positives: causal variants that fall within detected intervals
#'   - False negatives: causal variants not within any detected interval
#'   - False positives: detected peaks that don't contain any causal variant
#'
#' @param mapping_data Processed mapping dataframe (from analyze_mapping())
#'   with columns: marker, CHROM, POS, P, log10p, significant, peak_id,
#'   startPOS, peakPOS, endPOS, AF1, BETA, var.exp
#' @param qtl_regions Dataframe from extract_qtl_regions()
#' @param causal_variants Dataframe from load_causal_variants()
#' @param mapping_params Named list with simulation parameters:
#'   nqtl, rep, h2, maf, effect, population, algorithm, pca, threshold_method.
#'   Optional: mode, type, alpha, ci_size, snp_grouping
#' @return Dataframe with one row per QTL (simulated + detected, deduplicated)
compile_full_assessment <- function(mapping_data, qtl_regions, causal_variants,
                                    mapping_params) {

  # Step 1: Score causal variant markers
  effects_scores <- score_causal_markers(mapping_data, causal_variants)

  # Step 2: Build peak info from QTL regions (one row per peak: highest log10p marker)
  if (nrow(qtl_regions) == 0) {
    peak_info <- data.frame(
      CHROM = character(), marker = character(), POS = integer(),
      AF1 = numeric(), BETA = numeric(), log10p = numeric(),
      startPOS = integer(), peakPOS = integer(), endPOS = integer(),
      peak_id = integer(), interval_size = integer(),
      detected.peak = character(),
      stringsAsFactors = FALSE
    )
  } else {
    # Compute log10p if not present
    if (!"log10p" %in% names(mapping_data)) {
      mapping_data <- mapping_data %>%
        dplyr::mutate(log10p = safe_log10p(P))
    }
    peak_info <- mapping_data %>%
      dplyr::filter(!is.na(peak_id)) %>%
      dplyr::group_by(peak_id) %>%
      dplyr::slice_max(log10p, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(
        CHROM, marker, POS, AF1, BETA, log10p,
        startPOS, peakPOS, endPOS, peak_id, interval_size,
        dplyr::any_of("var.exp")
      ) %>%
      dplyr::mutate(detected.peak = marker)
  }

  # Step 3: Find peak-causal overlaps
  overlap_df <- find_peak_causal_overlaps(peak_info, effects_scores)

  # Step 4: Build assessment union
  build_assessment_union(effects_scores, overlap_df, mapping_params)
}


# =============================================================================
# Output Formatting
# =============================================================================

#' Format assessment dataframe for TSV output
#'
#' Selects and orders columns to match the standardized output schema shared
#' between DB-path and legacy-path assessments.
#'
#' @param assessment_df Dataframe from compile_full_assessment()
#' @return Dataframe with standardized column ordering
format_assessment_tsv <- function(assessment_df) {
  if (nrow(assessment_df) == 0) return(assessment_df)

  # Standardized schema: columns shared between DB-path and legacy-path
  expected_cols <- c(
    "QTL", "Simulated", "Detected", "CHROM", "POS", "RefAllele", "Frequency", "Effect",
    "Simulated.QTL.VarExp", "log10p", "significant", "BETA",
    "startPOS", "peakPOS", "endPOS", "detected.peak",
    "interval.log10p", "interval.var.exp", "interval.Frequency",
    "peak_id", "interval_size", "top.hit",
    "nQTL", "simREP", "h2", "maf", "effect_distribution", "strain_set_id",
    "mode", "type", "threshold", "algorithm_id", "alpha", "ci_size", "snp_grouping"
  )

  for (col in expected_cols) {
    if (!col %in% names(assessment_df)) {
      assessment_df[[col]] <- NA
    }
  }

  assessment_df %>%
    dplyr::select(dplyr::all_of(expected_cols))
}
