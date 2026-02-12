# assessment.R - Simulation assessment functions for DB-path QTL analysis
#
# Implements the assessment logic from bin/Assess_Sims.R without
# GenomicRanges or genotype matrix dependencies:
#   - load_causal_variants(): read .par file
#   - assess_qtl_detection(): match causal variants to detected QTL intervals
#   - compile_full_assessment(): build the full Simulated/Detected dataframe
#   - format_assessment_tsv(): format output to match existing pipeline output
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
# QTL Detection Assessment
# =============================================================================

#' Check which causal variants fall within detected QTL intervals
#'
#' For each causal variant, checks if its position falls within any detected
#' QTL interval (startPOS <= POS <= endPOS on the same chromosome).
#' Uses simple position-based overlap (equivalent to GenomicRanges for
#' point-in-interval queries).
#'
#' @param qtl_regions Dataframe from extract_qtl_regions() with columns:
#'   peak_id, CHROM, startPOS, peakPOS, endPOS, peak_marker, max_log10p, etc.
#' @param causal_variants Dataframe from load_causal_variants() with columns:
#'   QTL, CHROM, POS
#' @return Dataframe with one row per (causal_variant, matched_interval) pair,
#'   including unmatched causal variants (left join behavior)
assess_qtl_detection <- function(qtl_regions, causal_variants) {

  if (nrow(qtl_regions) == 0) {
    # No detected intervals — all causal variants are undetected
    log_msg("No QTL regions detected — all causal variants are undetected")
    result <- causal_variants %>%
      dplyr::mutate(
        detected_peak_id = NA_integer_,
        detected_peak_marker = NA_character_,
        startPOS = NA_integer_,
        peakPOS = NA_integer_,
        endPOS = NA_integer_
      )
    return(result)
  }

  if (nrow(causal_variants) == 0) {
    log_msg("No causal variants — returning empty assessment", level = "WARN")
    return(data.frame())
  }

  # For each causal variant, find overlapping QTL intervals
  # Simple position-based overlap: causal POS within [startPOS, endPOS] on same CHROM
  overlap <- causal_variants %>%
    dplyr::left_join(
      qtl_regions %>%
        dplyr::select(peak_id, CHROM, startPOS, peakPOS, endPOS, peak_marker),
      by = "CHROM",
      relationship = "many-to-many"
    ) %>%
    dplyr::filter(
      !is.na(startPOS) & POS >= startPOS & POS <= endPOS
    ) %>%
    dplyr::rename(
      detected_peak_id = peak_id,
      detected_peak_marker = peak_marker
    )

  # Left join back to get unmatched causal variants
  matched_qtls <- unique(overlap$QTL)
  unmatched <- causal_variants %>%
    dplyr::filter(!QTL %in% matched_qtls) %>%
    dplyr::mutate(
      detected_peak_id = NA_integer_,
      detected_peak_marker = NA_character_,
      startPOS = NA_integer_,
      peakPOS = NA_integer_,
      endPOS = NA_integer_
    )

  result <- dplyr::bind_rows(overlap, unmatched)

  n_detected <- length(matched_qtls)
  n_total <- nrow(causal_variants)
  log_msg(glue::glue(
    "Assessment: {n_detected}/{n_total} causal variants detected in QTL intervals"
  ))

  result
}


# =============================================================================
# Full Assessment Compilation
# =============================================================================

#' Compile full assessment combining simulated and detected QTLs
#'
#' Produces a dataframe matching the structure of Assess_Sims.R output,
#' containing:
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
#'   nqtl, rep, h2, maf, effect, population, algorithm, pca, threshold_method
#' @return Dataframe with one row per QTL (simulated + detected, deduplicated)
compile_full_assessment <- function(mapping_data, qtl_regions, causal_variants,
                                    mapping_params) {

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

  # Build peak info from QTL regions (detected intervals)
  if (nrow(qtl_regions) == 0) {
    peak_info <- data.frame(
      CHROM = character(), marker = character(), POS = integer(),
      AF1 = numeric(), BETA = numeric(), log10p = numeric(),
      startPOS = integer(), peakPOS = integer(), endPOS = integer(),
      peak_id = integer(), interval_size = integer(), var.exp = numeric(),
      detected.peak = character(),
      stringsAsFactors = FALSE
    )
  } else {
    # One row per peak: the marker with highest log10p in each peak group
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

  # Match causal variants to detected intervals using position overlap
  if (nrow(peak_info) > 0 && nrow(effects_scores) > 0) {
    # For each detected peak, check which causal variants overlap
    overlap_results <- lapply(seq_len(nrow(peak_info)), function(i) {
      peak <- peak_info[i, ]
      overlapping <- effects_scores %>%
        dplyr::filter(CHROM == peak$CHROM,
                      POS >= peak$startPOS,
                      POS <= peak$endPOS)
      if (nrow(overlapping) > 0) {
        overlapping %>%
          dplyr::mutate(
            startPOS = peak$startPOS,
            peakPOS = peak$peakPOS,
            endPOS = peak$endPOS,
            detected.peak = peak$detected.peak,
            interval.log10p = peak$log10p,
            interval.var.exp = if ("var.exp" %in% names(peak)) peak$var.exp else NA_real_,
            interval.Frequency = peak$AF1,
            peak_id = peak$peak_id,
            interval_size = peak$interval_size
          )
      } else {
        # No causal variant overlaps — this is a false positive peak
        data.frame(
          QTL = peak$detected.peak,
          CHROM = peak$CHROM,
          POS = peak$POS,
          RefAllele = NA_character_,
          Frequency = NA_real_,
          Effect = NA_real_,
          log10p = NA_real_,
          significant = NA_integer_,
          startPOS = peak$startPOS,
          peakPOS = peak$peakPOS,
          endPOS = peak$endPOS,
          detected.peak = peak$detected.peak,
          interval.log10p = peak$log10p,
          interval.var.exp = if ("var.exp" %in% names(peak)) peak$var.exp else NA_real_,
          interval.Frequency = peak$AF1,
          peak_id = peak$peak_id,
          interval_size = peak$interval_size,
          stringsAsFactors = FALSE
        )
      }
    })
    overlap_df <- dplyr::bind_rows(overlap_results)
  } else {
    overlap_df <- data.frame(
      QTL = character(),
      stringsAsFactors = FALSE
    )
  }

  # Build the all.QTL dataframe
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

  # Join effect scores (causal variant info)
  if (nrow(effects_scores) > 0) {
    join_cols <- intersect(
      c("QTL", "log10p", "significant"),
      names(effects_scores)
    )
    all_qtl <- all_qtl %>%
      dplyr::left_join(
        effects_scores %>% dplyr::select(dplyr::all_of(join_cols)),
        by = "QTL"
      )
  }

  # Join overlap info (detected interval info)
  if (nrow(overlap_df) > 0) {
    overlap_join_cols <- intersect(
      c("QTL", "startPOS", "peakPOS", "endPOS", "detected.peak",
        "interval.log10p", "interval.var.exp", "interval.Frequency",
        "peak_id", "interval_size"),
      names(overlap_df)
    )
    all_qtl <- all_qtl %>%
      dplyr::left_join(
        overlap_df %>%
          dplyr::select(dplyr::all_of(overlap_join_cols)) %>%
          dplyr::filter(!duplicated(QTL)),
        by = "QTL"
      )
  }

  # Add top.hit flag and simulation metadata
  algorithm_id <- paste0(
    mapping_params$algorithm,
    if (isTRUE(mapping_params$pca)) "_PCA" else "_noPCA",
    "_", mapping_params$threshold_method
  )

  all_qtl <- all_qtl %>%
    dplyr::mutate(
      top.hit = if ("detected.peak" %in% names(.)) QTL == detected.peak else NA,
      nQTL = as.character(mapping_params$nqtl),
      simREP = as.character(mapping_params$rep),
      h2 = as.character(mapping_params$h2),
      maf = as.character(mapping_params$maf),
      effect_distribution = as.character(mapping_params$effect),
      strain_set_id = as.character(mapping_params$population),
      algorithm_id = algorithm_id,
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
# Output Formatting
# =============================================================================

#' Format assessment dataframe for TSV output
#'
#' Selects and orders columns to match the output format of Assess_Sims.R.
#' Note: Simulated.QTL.VarExp and var.exp are NA in DB-path assessments
#' since the genotype matrix is not available.
#'
#' @param assessment_df Dataframe from compile_full_assessment()
#' @return Dataframe with columns ordered to match existing output
format_assessment_tsv <- function(assessment_df) {
  if (nrow(assessment_df) == 0) return(assessment_df)

  # Ensure all expected columns exist (add NAs for missing)
  expected_cols <- c(
    "QTL", "Simulated", "Detected", "log10p", "significant",
    "startPOS", "peakPOS", "endPOS", "detected.peak",
    "interval.log10p", "interval.var.exp", "interval.Frequency",
    "peak_id", "interval_size", "top.hit",
    "nQTL", "simREP", "h2", "maf", "effect_distribution",
    "strain_set_id", "algorithm_id"
  )

  for (col in expected_cols) {
    if (!col %in% names(assessment_df)) {
      assessment_df[[col]] <- NA
    }
  }

  # Select and order columns
  assessment_df %>%
    dplyr::select(dplyr::all_of(expected_cols))
}
