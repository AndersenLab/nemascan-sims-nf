# sim_performance.R - Simulation performance analysis functions
#
# Adapted from Caenorhabditis GWAS Manuscript analysis pipeline.
# Provides reusable functions for classifying QTL detection outcomes
# and calculating power/FDR metrics.
#
# Functions:
#   designate_qtl()                - Classify each QTL row into 4 outcome categories
#   count_outcomes()               - Pivot designation counts per group
#   calculate_simrep_performance() - Compute Power and FDR per simulation replicate
#   calculate_varexp_performance() - Compute Power without FDR (for variance-explained analysis)
#
# Dependencies: dplyr, tidyr (both part of tidyverse, already available)

library(dplyr)
library(tidyr)


#' Classify each QTL row into detection outcome categories
#'
#' Adds a `designation` column with one of four values:
#'   - Detected.CV: Simulated=TRUE, Detected=TRUE, significant=TRUE or NA (true positive)
#'   - Missed.CV: Simulated=TRUE, Detected=FALSE, significant=FALSE or NA (false negative)
#'   - CV.Not.Significant.In.Interval: Simulated=TRUE, Detected=TRUE, significant=FALSE
#'   - False.Discovery: Simulated=FALSE, Detected=TRUE, significant=TRUE (false positive)
#'
#' significant=NA rows arise when a causal variant is absent from GWA output entirely
#' (e.g., when cv_maf < ms_maf produces non-marker causal variants). When Detected=TRUE,
#' the variant fell inside a detected QTL interval — it was correctly localized despite
#' having no direct GWA test, so it counts as Detected.CV (true positive). When
#' Detected=FALSE, it counts as Missed.CV (false negative for Power calculation).
#'
#' Auto-detects `significant` (Phase 4 schema) vs `aboveBF` (legacy schema) column.
#' The significance column is coerced to logical for consistent case_when evaluation.
#'
#' @param df A dataframe with columns: Simulated, Detected, and either
#'   `significant` or `aboveBF`
#' @return The input dataframe with an additional `designation` column
designate_qtl <- function(df) {
  # Auto-detect significance column
  if ("significant" %in% names(df)) {
    sig_col <- "significant"
  } else if ("aboveBF" %in% names(df)) {
    sig_col <- "aboveBF"
  } else {
    stop("designate_qtl(): dataframe must have a 'significant' or 'aboveBF' column")
  }

  # Coerce Simulated/Detected to logical
  df <- df %>%
    dplyr::mutate(
      .sim = as.logical(Simulated),
      .det = as.logical(Detected),
      .sig = as.logical(.data[[sig_col]])
    )

  df %>%
    dplyr::mutate(
      designation = dplyr::case_when(
        # significant=NA with Detected=TRUE: non-marker variant inside a detected interval
        # (correctly localized despite no direct GWA test → true positive)
        .sim == TRUE  & .det == TRUE  & (.sig == TRUE | is.na(.sig))  ~ "Detected.CV",
        # significant=NA with Detected=FALSE: non-marker variant not in any detected interval
        # (structurally undetectable → false negative for Power denominator)
        .sim == TRUE  & .det == FALSE & (.sig == FALSE | is.na(.sig)) ~ "Missed.CV",
        .sim == TRUE  & .det == TRUE  & .sig == FALSE ~ "CV.Not.Significant.In.Interval",
        .sim == FALSE & .det == TRUE  & .sig == TRUE  ~ "False.Discovery"
      )
    ) %>%
    dplyr::select(-dplyr::all_of(c(".sim", ".det", ".sig")))
}


#' Count detection outcomes for a grouped dataframe
#'
#' Pivots designation counts per group. The input should already be grouped
#' (e.g., by method, replicate, etc.) before calling via group_modify().
#'
#' Ensures all four designation columns exist in the output, filling with 0
#' for any categories absent in the data.
#'
#' @param grouped_designation_df A (possibly grouped) dataframe with a
#'   `designation` column produced by designate_qtl()
#' @return An ungrouped dataframe with columns: Detected.CV, Missed.CV,
#'   CV.Not.Significant.In.Interval, False.Discovery (plus any grouping columns)
count_outcomes <- function(grouped_designation_df) {
  summary_for_group <- grouped_designation_df %>%
    dplyr::count(designation, name = "n_counts") %>%
    tidyr::pivot_wider(
      names_from = designation,
      values_from = n_counts,
      values_fill = 0
    )

  expected_cols <- c(
    "Detected.CV",
    "Missed.CV",
    "CV.Not.Significant.In.Interval",
    "False.Discovery"
  )

  for (col_name in expected_cols) {
    if (!col_name %in% names(summary_for_group)) {
      summary_for_group <- summary_for_group %>%
        dplyr::mutate(!!rlang::sym(col_name) := 0)
    }
  }

  dplyr::ungroup(summary_for_group)
}


#' Calculate Power and FDR for simulation replicates
#'
#' Computes Power and FDR from designation count columns.
#'
#' Power definitions:
#'   - Non-stringent (default): (Detected.CV + CV.Not.Significant.In.Interval) / Simulated
#'   - Stringent: Detected.CV / Simulated
#'
#' FDR: False.Discovery / (Detected.CV + CV.Not.Significant.In.Interval + False.Discovery)
#'
#' @param count_df A dataframe with columns: Detected.CV, Missed.CV,
#'   CV.Not.Significant.In.Interval, False.Discovery
#' @param stringent Logical. If TRUE, only count significant detections as power.
#'   Defaults to FALSE.
#' @return The input dataframe with added columns: Detected, Simulated, Power, FDR
calculate_simrep_performance <- function(count_df, stringent = FALSE) {
  if (stringent) {
    result <- count_df %>%
      dplyr::mutate(
        Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery,
        Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
        Power = dplyr::if_else(Simulated == 0, 0, Detected.CV / Simulated),
        FDR = dplyr::if_else(Detected == 0, 0, False.Discovery / Detected)
      )
  } else {
    result <- count_df %>%
      dplyr::mutate(
        Detected = Detected.CV + CV.Not.Significant.In.Interval + False.Discovery,
        Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
        Power = dplyr::if_else(
          Simulated == 0, 0,
          (Detected.CV + CV.Not.Significant.In.Interval) / Simulated
        ),
        FDR = dplyr::if_else(Detected == 0, 0, False.Discovery / Detected)
      )
  }
  result
}


#' Calculate Power for variance-explained analysis (no FDR)
#'
#' Similar to calculate_simrep_performance() but excludes False.Discovery
#' from the calculation, as variance-explained analysis focuses on
#' detection sensitivity without false positive rate.
#'
#' @param count_df A dataframe with columns: Detected.CV, Missed.CV,
#'   CV.Not.Significant.In.Interval, False.Discovery
#' @return The input dataframe with added columns: Detected, Simulated, Power
#'   (False.Discovery column removed)
calculate_varexp_performance <- function(count_df) {
  result <- count_df %>%
    dplyr::select(-False.Discovery) %>%
    dplyr::mutate(
      Detected = Detected.CV + CV.Not.Significant.In.Interval,
      Simulated = Detected.CV + CV.Not.Significant.In.Interval + Missed.CV,
      Power = dplyr::if_else(
        Simulated == 0, 0,
        (Detected.CV + CV.Not.Significant.In.Interval) / Simulated
      )
    )

  dplyr::ungroup(result)
}
