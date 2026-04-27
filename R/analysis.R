# analysis.R - Query-time analysis functions for threshold and interval calculation
# These functions work with dataframes from database queries

library(dplyr)
library(glue)


# =============================================================================
# Safe log10p Computation
# =============================================================================

#' Maximum value for -log10(p) to prevent infinite values
LOG10P_MAX <- 300

#' Compute -log10(p) with safe handling of zero and extreme values
#'
#' Computes -log10(p), replacing infinite values (from p = 0) with LOG10P_MAX.
#' This replaces the previously-stored log10p column in the database schema.
#'
#' @param p Numeric vector of p-values
#' @return Numeric vector of -log10(p) values, capped at LOG10P_MAX
safe_log10p <- function(p) {
  log10p <- -log10(p)
  log10p[is.infinite(log10p)] <- LOG10P_MAX
  log10p
}


# =============================================================================
# Threshold Calculation Functions
# =============================================================================

#' Calculate significance threshold value
#'
#' Calculates the -log10(p) threshold for calling significant markers.
#'
#' @param method Threshold method: "BF", "EIGEN", or numeric value
#' @param n_markers Number of unique markers (required for BF method)
#' @param n_independent Number of independent tests from EIGEN (required for EIGEN method)
#' @param alpha Significance level (default: 0.05)
#' @return Named list with threshold_value and threshold_method
calculate_threshold <- function(method,
                                n_markers = NULL,
                                n_independent = NULL,
                                alpha = 0.05) {

  if (is.numeric(method)) {
    # User-specified threshold
    threshold_value <- method
    threshold_method <- "USER"
    log_msg(glue::glue("Using user-specified threshold: {threshold_value}"))

  } else if (toupper(method) == "BF") {
    # Bonferroni correction
    if (is.null(n_markers) || n_markers <= 0) {
      stop("BF method requires positive n_markers parameter")
    }
    threshold_value <- -log10(alpha / n_markers)
    threshold_method <- "BF"
    log_msg(glue::glue(
      "Bonferroni threshold: -log10({alpha}/{n_markers}) = {round(threshold_value, 4)}"
    ))

  } else if (toupper(method) == "EIGEN") {
    # EIGEN-based threshold
    if (is.null(n_independent) || n_independent <= 0) {
      stop("EIGEN method requires positive n_independent parameter")
    }
    threshold_value <- -log10(alpha / n_independent)
    threshold_method <- "EIGEN"
    log_msg(glue::glue(
      "EIGEN threshold: -log10({alpha}/{n_independent}) = {round(threshold_value, 4)}"
    ))

  } else {
    stop(glue::glue("Unknown significance method: {method}"))
  }

  list(
    threshold_value = threshold_value,
    threshold_method = threshold_method
  )
}


#' Flag markers above significance threshold
#'
#' Adds significance columns to mapping dataframe.
#' Computes log10p from P values using safe_log10p() at query time.
#'
#' @param df Mapping dataframe with P column
#' @param threshold_value Numeric threshold on -log10 scale
#' @param threshold_method Character string describing method used (optional)
#' @param verbose Whether to log messages (default: TRUE)
#' @return Dataframe with added columns: log10p, sig_threshold_value, sig_threshold_method, significant
flag_significant_markers <- function(df,
                                     threshold_value,
                                     threshold_method = "USER",
                                     verbose = TRUE) {

  result <- df %>%
    dplyr::mutate(
      log10p = safe_log10p(P),
      sig_threshold_value = threshold_value,
      sig_threshold_method = threshold_method,
      significant = as.integer(log10p >= threshold_value)
    )

  if (verbose) {
    n_sig <- sum(result$significant, na.rm = TRUE)
    n_total <- nrow(result)
    log_msg(glue::glue(
      "Found {n_sig} significant markers out of {n_total} ({round(100*n_sig/n_total, 2)}%)"
    ))
  }

  result
}


# =============================================================================
# QTL Interval Functions
# =============================================================================

#' Group significant markers into QTL intervals
#'
#' Assigns peak IDs to significant markers based on proximity.
#' Markers within snp_grouping distance on the same chromosome are grouped together.
#'
#' @param df Mapping dataframe with significant column
#' @param snp_grouping Number of markers - if two significant markers are within this
#'   index distance, they are assigned to the same peak
#' @param verbose Whether to log messages (default: TRUE)
#' @return Dataframe with peak_id and marker_index columns added
group_markers_to_intervals <- function(df, snp_grouping = 1000, verbose = TRUE) {
  n_sig <- sum(df$significant, na.rm = TRUE)

  if (n_sig == 0) {
    if (verbose) log_msg("No significant markers found, skipping interval grouping")
    return(df %>%
             dplyr::mutate(
               marker_index = NA_integer_,
               peak_id = NA_integer_
             ))
  }

  if (verbose) log_msg(glue::glue("Grouping {n_sig} significant markers into QTL intervals..."))

  # Add row ID for tracking and create marker index within each chromosome
  indexed_df <- df %>%
    dplyr::mutate(.row_id = dplyr::row_number()) %>%
    dplyr::group_by(CHROM) %>%
    dplyr::arrange(POS, .by_group = TRUE) %>%
    dplyr::mutate(marker_index = dplyr::row_number()) %>%
    dplyr::ungroup()

  # Extract significant markers for peak assignment
  sig_markers <- indexed_df %>%
    dplyr::filter(significant == 1) %>%
    dplyr::arrange(CHROM, POS)

  # Assign peak IDs based on proximity
  if (nrow(sig_markers) == 1) {
    sig_markers$peak_id <- 1L
  } else {
    sig_markers$peak_id <- 1L

    for (i in 2:nrow(sig_markers)) {
      prev_chrom <- sig_markers$CHROM[i - 1]
      curr_chrom <- sig_markers$CHROM[i]
      prev_index <- sig_markers$marker_index[i - 1]
      curr_index <- sig_markers$marker_index[i]

      # Same peak if same chromosome and within grouping distance
      if (curr_chrom == prev_chrom &&
          abs(curr_index - prev_index) < snp_grouping) {
        sig_markers$peak_id[i] <- sig_markers$peak_id[i - 1]
      } else {
        sig_markers$peak_id[i] <- sig_markers$peak_id[i - 1] + 1L
      }
    }
  }

  n_peaks <- length(unique(sig_markers$peak_id))
  if (verbose) log_msg(glue::glue("Identified {n_peaks} distinct QTL intervals"))

  # Join peak IDs back using row ID (avoids duplicate marker issues)
  result <- indexed_df %>%
    dplyr::left_join(
      sig_markers %>% dplyr::select(.row_id, peak_id),
      by = ".row_id"
    ) %>%
    dplyr::select(-.row_id)

  result
}


#' Define confidence intervals for QTL peaks
#'
#' For each QTL peak, defines start and end positions based on CI size.
#' Computes log10p from P values using safe_log10p().
#'
#' @param df Mapping dataframe with peak_id, marker_index, and P columns
#' @param ci_size Number of markers to left and right of peak for CI boundary
#' @param verbose Whether to log messages (default: TRUE)
#' @return Dataframe with startPOS, peakPOS, endPOS, interval_size columns added
define_confidence_intervals <- function(df, ci_size = 150, verbose = TRUE) {
  # Check if any peaks exist
  if (all(is.na(df$peak_id))) {
    if (verbose) log_msg("No QTL peaks to define intervals for")
    return(df %>%
             dplyr::mutate(
               startPOS = NA_integer_,
               peakPOS = NA_integer_,
               endPOS = NA_integer_,
               interval_size = NA_integer_
             ) %>%
             dplyr::select(-marker_index))
  }

  if (verbose) log_msg(glue::glue("Defining confidence intervals (CI size: {ci_size} markers)..."))

  # Compute log10p from P
  df <- df %>% dplyr::mutate(log10p = safe_log10p(P))

  # Get chromosome-wise index boundaries
  chrom_bounds <- df %>%
    dplyr::group_by(CHROM) %>%
    dplyr::summarise(
      min_index = min(marker_index),
      max_index = max(marker_index),
      .groups = "drop"
    )

  # Build position lookup (index -> POS for each chromosome)
  pos_lookup <- df %>%
    dplyr::select(CHROM, marker_index, POS) %>%
    dplyr::distinct()

  # For each peak, find the peak marker and define CI boundaries
  peak_info <- df %>%
    dplyr::filter(!is.na(peak_id)) %>%
    dplyr::group_by(CHROM, peak_id) %>%
    dplyr::summarise(
      # Peak is the marker with highest log10p in the group
      peak_index = marker_index[which.max(log10p)],
      peakPOS = POS[which.max(log10p)],
      # CI boundaries based on min/max index of significant markers in group
      min_sig_index = min(marker_index),
      max_sig_index = max(marker_index),
      .groups = "drop"
    ) %>%
    # Calculate CI start/end indices
    dplyr::mutate(
      start_index = min_sig_index - ci_size,
      end_index = max_sig_index + ci_size
    ) %>%
    # Join chromosome bounds to constrain indices
    dplyr::left_join(chrom_bounds, by = "CHROM") %>%
    # Constrain to valid index range
    dplyr::mutate(
      start_index = pmax(start_index, min_index),
      end_index = pmin(end_index, max_index)
    )

  # Look up positions for start/end indices
  peak_info <- peak_info %>%
    dplyr::left_join(
      pos_lookup %>% dplyr::rename(start_index = marker_index, startPOS = POS),
      by = c("CHROM", "start_index")
    ) %>%
    dplyr::left_join(
      pos_lookup %>% dplyr::rename(end_index = marker_index, endPOS = POS),
      by = c("CHROM", "end_index")
    ) %>%
    dplyr::mutate(interval_size = endPOS - startPOS) %>%
    dplyr::select(CHROM, peak_id, startPOS, peakPOS, endPOS, interval_size)

  # Join interval info back to full dataframe
  result <- df %>%
    dplyr::left_join(peak_info, by = c("CHROM", "peak_id")) %>%
    dplyr::select(-marker_index)

  result
}


# =============================================================================
# Combined Analysis Pipeline
# =============================================================================

#' Analyze mapping with significance threshold and QTL intervals
#'
#' Main analysis function that applies threshold and defines QTL intervals.
#'
#' @param df Mapping dataframe from database query (must have CHROM, POS, P)
#' @param threshold_value Numeric threshold on -log10 scale
#' @param threshold_method Character describing method (default: "USER")
#' @param ci_size Number of markers for CI boundary (default: 150)
#' @param snp_grouping Distance for grouping markers into same peak (default: 1000)
#' @param verbose Whether to log messages (default: TRUE)
#' @return Fully processed mapping dataframe with QTL annotations
analyze_mapping <- function(df,
                            threshold_value,
                            threshold_method = "USER",
                            ci_size = 150,
                            snp_grouping = 1000,
                            verbose = TRUE) {

  if (verbose) {
    log_msg("Starting mapping analysis...")
    log_msg(glue::glue("  Threshold: {round(threshold_value, 4)} ({threshold_method})"))
    log_msg(glue::glue("  CI size: {ci_size}"))
    log_msg(glue::glue("  SNP grouping: {snp_grouping}"))
  }

  # Step 1: Flag significant markers (computes log10p from P)
  processed <- flag_significant_markers(
    df = df,
    threshold_value = threshold_value,
    threshold_method = threshold_method,
    verbose = verbose
  )

  # Step 2: Group significant markers into intervals
  processed <- group_markers_to_intervals(
    df = processed,
    snp_grouping = snp_grouping,
    verbose = verbose
  )

  # Step 3: Define confidence intervals
  processed <- define_confidence_intervals(
    df = processed,
    ci_size = ci_size,
    verbose = verbose
  )

  # Summary
  if (verbose) {
    n_sig <- sum(processed$significant, na.rm = TRUE)
    n_peaks <- length(unique(na.omit(processed$peak_id)))
    log_msg(glue::glue(
      "Analysis complete: {n_sig} significant markers in {n_peaks} QTL intervals"
    ))
  }

  processed
}


# =============================================================================
# Summary and Extraction Functions
# =============================================================================

#' Extract QTL regions from processed mapping
#'
#' Returns a summary dataframe with one row per QTL interval.
#' Computes log10p from P values using safe_log10p().
#'
#' @param df Processed mapping dataframe with QTL annotations and P column
#' @return Dataframe with QTL region summaries
extract_qtl_regions <- function(df) {
  if (!"peak_id" %in% names(df) || all(is.na(df$peak_id))) {
    log_msg("No QTL regions to extract")
    return(data.frame())
  }

  # Compute log10p from P
  df <- df %>% dplyr::mutate(log10p = safe_log10p(P))

  regions <- df %>%
    dplyr::filter(!is.na(peak_id)) %>%
    dplyr::group_by(peak_id) %>%
    dplyr::summarise(
      CHROM = unique(CHROM)[1],
      startPOS = unique(startPOS)[1],
      peakPOS = unique(peakPOS)[1],
      endPOS = unique(endPOS)[1],
      interval_size = unique(interval_size)[1],
      n_sig_markers = sum(significant),
      max_log10p = max(log10p),
      peak_marker = marker[which.max(log10p)],
      sig_threshold_value = unique(sig_threshold_value)[1],
      sig_threshold_method = unique(sig_threshold_method)[1],
      .groups = "drop"
    ) %>%
    dplyr::arrange(CHROM, peakPOS)

  log_msg(glue::glue("Extracted {nrow(regions)} QTL regions"))

  regions
}


#' Summarize detection results
#'
#' Returns summary statistics for a processed mapping.
#' Computes log10p from P values using safe_log10p().
#'
#' @param df Processed mapping dataframe with P column
#' @return Named list with detection summary statistics
summarize_detection <- function(df) {
  n_markers <- nrow(df)

  # Compute log10p from P
  df <- df %>% dplyr::mutate(log10p = safe_log10p(P))

  n_sig <- sum(df$significant, na.rm = TRUE)
  n_qtl <- length(unique(na.omit(df$peak_id)))
  max_log10p <- max(df$log10p, na.rm = TRUE)

  # Get threshold info if available
  threshold_value <- if ("sig_threshold_value" %in% names(df)) {
    unique(df$sig_threshold_value)[1]
  } else {
    NA_real_
  }

  threshold_method <- if ("sig_threshold_method" %in% names(df)) {
    unique(df$sig_threshold_method)[1]
  } else {
    NA_character_
  }

  list(
    n_markers = n_markers,
    n_significant = n_sig,
    pct_significant = round(100 * n_sig / n_markers, 2),
    n_qtl = n_qtl,
    max_log10p = round(max_log10p, 4),
    threshold_value = threshold_value,
    threshold_method = threshold_method
  )
}


# =============================================================================
# Batch Processing Functions
# =============================================================================

#' Process a single mapping from database with full pipeline
#'
#' Convenience function that queries a mapping from the database, calculates
#' the appropriate threshold, and runs the full analysis pipeline.
#'
#' @param mapping_id Mapping ID to process
#' @param method Threshold method: "BF", "EIGEN", or numeric value
#' @param base_dir Database root directory
#' @param alpha Significance level (default: 0.05)
#' @param ci_size Number of markers for CI boundary (default: 150)
#' @param snp_grouping Distance for grouping markers into same peak (default: 1000)
#' @param verbose Whether to log progress messages (default: TRUE)
#' @return List with processed data and summary statistics
process_single_mapping <- function(mapping_id,
                                   method = "EIGEN",
                                   base_dir = "data/db",
                                   alpha = 0.05,
                                   ci_size = 150,
                                   snp_grouping = 1000,
                                   verbose = TRUE) {

  if (verbose) log_msg(glue::glue("Processing mapping: {mapping_id}"))

  # Step 1: Query mapping data from database
  mapping_data <- query_for_threshold_analysis(mapping_id, base_dir)

  if (nrow(mapping_data) == 0) {
    stop(glue::glue("No data found for mapping_id: {mapping_id}"))
  }

  # Extract population and MAF
  population <- unique(mapping_data$population)[1]
  maf <- unique(mapping_data$maf)[1]

  # Step 2: Get threshold parameters
  params <- get_threshold_params(population, maf, alpha, base_dir)

  # Step 3: Calculate threshold based on method
  if (is.numeric(method)) {
    threshold <- calculate_threshold(method)
  } else if (toupper(method) == "BF") {
    threshold <- calculate_threshold("BF", n_markers = params$n_markers, alpha = alpha)
  } else if (toupper(method) == "EIGEN") {
    if (is.na(params$n_independent_tests)) {
      stop(glue::glue("EIGEN data not available for {population}_{maf}"))
    }
    threshold <- calculate_threshold("EIGEN", n_independent = params$n_independent_tests, alpha = alpha)
  } else {
    stop(glue::glue("Unknown threshold method: {method}"))
  }

  # Step 4-6: Run analysis pipeline
  processed <- analyze_mapping(
    df = mapping_data,
    threshold_value = threshold$threshold_value,
    threshold_method = threshold$threshold_method,
    ci_size = ci_size,
    snp_grouping = snp_grouping
  )

  # Step 7: Extract results
  qtl_regions <- extract_qtl_regions(processed)
  summary_stats <- summarize_detection(processed)

  list(
    mapping_id = mapping_id,
    population = population,
    maf = maf,
    processed_data = processed,
    qtl_regions = qtl_regions,
    summary = summary_stats,
    threshold_params = params
  )
}
