# test-analyze_mapping.R - Unit tests for analyze_mapping() and sub-functions in R/analysis.R
#
# Covers the behavioral change from removing check_significant_proportion():
# analyze_mapping() now always attempts QTL interval definition regardless of
# what fraction of markers are significant.
#
# Organized into four sections:
#   1. flag_significant_markers()
#   2. group_markers_to_intervals()
#   3. define_confidence_intervals()
#   4. analyze_mapping() — regression tests for removed guard + integration

# ---- Helpers ----------------------------------------------------------------

# Minimal mapping dataframe with required columns: CHROM, POS, P, marker
make_mapping_df <- function(n = 20, chrom = "I", p_null = 0.5) {
  pos <- as.integer(seq(1000, by = 1000, length.out = n))
  data.frame(
    CHROM  = chrom,
    POS    = pos,
    P      = rep(p_null, n),
    marker = paste0(chrom, ":", pos),
    stringsAsFactors = FALSE
  )
}

# Combine two single-chrom frames into a multi-chrom frame
two_chrom_df <- function(n = 10) {
  rbind(make_mapping_df(n, chrom = "I"), make_mapping_df(n, chrom = "II"))
}


# =============================================================================
# 1. flag_significant_markers()
# =============================================================================

test_that("flag_significant_markers adds required columns", {
  df     <- make_mapping_df(10)
  result <- flag_significant_markers(df, threshold_value = 5, verbose = FALSE)
  required <- c("log10p", "sig_threshold_value", "sig_threshold_method", "significant")
  expect_true(all(required %in% names(result)), label = "all flag columns present")
  expect_equal(nrow(result), 10L, label = "row count preserved")
})

test_that("flag_significant_markers marks markers >= threshold as significant", {
  df       <- make_mapping_df(5)
  df$P[3]  <- 1e-7   # log10p = 7, above threshold 5
  result   <- flag_significant_markers(df, threshold_value = 5, verbose = FALSE)
  expect_equal(result$significant[3], 1L, label = "marker at p=1e-7 is significant")
  expect_equal(sum(result$significant), 1L, label = "exactly one significant marker")
})

test_that("flag_significant_markers treats exactly-at-threshold as significant (>= not >)", {
  df          <- make_mapping_df(3)
  threshold   <- 5
  df$P[2]     <- 10^(-threshold)   # log10p exactly equals threshold
  result      <- flag_significant_markers(df, threshold_value = threshold, verbose = FALSE)
  expect_equal(result$significant[2], 1L, label = "boundary marker is significant")
})

test_that("flag_significant_markers caps p=0 at LOG10P_MAX=300", {
  df       <- make_mapping_df(3)
  df$P[1]  <- 0
  result   <- flag_significant_markers(df, threshold_value = 5, verbose = FALSE)
  expect_equal(result$log10p[1], 300, label = "p=0 maps to LOG10P_MAX")
  expect_equal(result$significant[1], 1L, label = "p=0 marker is significant")
})

test_that("flag_significant_markers returns all-zero significant when nothing exceeds threshold", {
  df     <- make_mapping_df(10)    # all p=0.5 → log10p ≈ 0.3, well below 5
  result <- flag_significant_markers(df, threshold_value = 5, verbose = FALSE)
  expect_equal(sum(result$significant), 0L, label = "no markers significant")
})

test_that("flag_significant_markers records threshold value and method in output", {
  df     <- make_mapping_df(5)
  result <- flag_significant_markers(df, threshold_value = 7.3, threshold_method = "BF",
                                     verbose = FALSE)
  expect_equal(unique(result$sig_threshold_value), 7.3, label = "threshold value stored")
  expect_equal(unique(result$sig_threshold_method), "BF", label = "threshold method stored")
})


# =============================================================================
# 2. group_markers_to_intervals()
# =============================================================================

test_that("group_markers_to_intervals returns all-NA peak_id when no significant markers", {
  df            <- make_mapping_df(10)
  df$significant <- 0L
  result        <- group_markers_to_intervals(df, verbose = FALSE)
  expect_true(all(is.na(result$peak_id)), label = "all peak_ids NA when n_sig=0")
})

test_that("group_markers_to_intervals assigns peak_id=1 for a single significant marker", {
  df            <- make_mapping_df(10)
  df$significant <- 0L
  df$significant[5] <- 1L
  result        <- group_markers_to_intervals(df, verbose = FALSE)
  sig_peak      <- result$peak_id[result$significant == 1]
  expect_equal(sig_peak, 1L, label = "single sig marker gets peak_id=1")
})

test_that("group_markers_to_intervals groups adjacent markers into the same peak", {
  df            <- make_mapping_df(20)
  df$significant <- 0L
  df$significant[c(5, 6)] <- 1L   # adjacent (index distance = 1 < snp_grouping=1000)
  result        <- group_markers_to_intervals(df, snp_grouping = 1000, verbose = FALSE)
  sig_peaks     <- result$peak_id[result$significant == 1]
  expect_equal(length(unique(sig_peaks)), 1L,
               label = "adjacent sig markers share one peak_id")
})

test_that("group_markers_to_intervals splits distant markers into separate peaks", {
  df            <- make_mapping_df(20)
  df$significant <- 0L
  df$significant[c(1, 20)] <- 1L   # index distance = 19 > snp_grouping=10
  result        <- group_markers_to_intervals(df, snp_grouping = 10, verbose = FALSE)
  sig_peaks     <- result$peak_id[result$significant == 1]
  expect_equal(length(unique(sig_peaks)), 2L,
               label = "distant sig markers get different peak_ids")
})

test_that("group_markers_to_intervals always separates markers on different chromosomes", {
  df            <- two_chrom_df(10)
  df$significant <- 0L
  # First marker on each chrom — index distance within each chrom is 1,
  # so if same-chrom logic applied they'd merge; different chroms must separate them.
  df$significant[c(1, 11)] <- 1L
  result <- group_markers_to_intervals(df, snp_grouping = 1000, verbose = FALSE)
  sig_peaks <- result$peak_id[result$significant == 1]
  expect_equal(length(unique(sig_peaks)), 2L,
               label = "markers on different chroms always get distinct peak_ids")
})

test_that("group_markers_to_intervals handles all-significant markers on one chromosome", {
  df            <- make_mapping_df(20)
  df$significant <- 1L   # 100% significant
  result        <- group_markers_to_intervals(df, snp_grouping = 1000, verbose = FALSE)
  n_peaks       <- length(unique(na.omit(result$peak_id)))
  expect_equal(n_peaks, 1L,
               label = "all sig markers on one chrom collapse to one peak with default grouping")
  expect_equal(sum(!is.na(result$peak_id)), 20L,
               label = "all markers receive a peak_id")
})


# =============================================================================
# 3. define_confidence_intervals()
# =============================================================================

test_that("define_confidence_intervals returns NA interval columns when no peaks exist", {
  df             <- make_mapping_df(10)
  df$significant <- 0L
  df$peak_id     <- NA_integer_
  df$marker_index <- seq_len(nrow(df))
  result         <- define_confidence_intervals(df, verbose = FALSE)
  expect_true(all(is.na(result$startPOS)),     label = "startPOS all NA when no peaks")
  expect_true(all(is.na(result$peakPOS)),      label = "peakPOS all NA when no peaks")
  expect_true(all(is.na(result$endPOS)),       label = "endPOS all NA when no peaks")
  expect_true(all(is.na(result$interval_size)), label = "interval_size all NA when no peaks")
})

test_that("define_confidence_intervals clamps CI at left chromosome boundary", {
  # sig marker at index 1 — ci_size=5 would try to extend to index -4
  df            <- make_mapping_df(10)
  df$significant <- 0L
  df$significant[1] <- 1L
  df            <- group_markers_to_intervals(df, verbose = FALSE)
  result        <- define_confidence_intervals(df, ci_size = 5, verbose = FALSE)
  # startPOS must be >= first position on the chromosome (1000)
  sig_start <- result$startPOS[result$significant == 1][1]
  expect_equal(sig_start, 1000L,
               label = "startPOS clamped at chromosome start position")
})

test_that("define_confidence_intervals clamps CI at right chromosome boundary", {
  # sig marker at last index — ci_size=5 would extend past the chromosome end
  df            <- make_mapping_df(10)
  df$significant <- 0L
  df$significant[10] <- 1L
  df            <- group_markers_to_intervals(df, verbose = FALSE)
  result        <- define_confidence_intervals(df, ci_size = 5, verbose = FALSE)
  sig_end <- result$endPOS[result$significant == 1][1]
  expect_equal(sig_end, 10000L,
               label = "endPOS clamped at chromosome end position")
})

test_that("define_confidence_intervals sets peakPOS to position of max log10p marker in peak", {
  df            <- make_mapping_df(10)
  df$significant <- 0L
  # Two adjacent sig markers; marker at index 5 has much lower p-value (higher log10p)
  df$P[4]  <- 1e-6
  df$P[5]  <- 1e-10   # this should be the peak
  df$significant[c(4, 5)] <- 1L
  df        <- group_markers_to_intervals(df, verbose = FALSE)
  result    <- define_confidence_intervals(df, verbose = FALSE)
  peak_pos  <- unique(na.omit(result$peakPOS))
  expect_equal(peak_pos, 5000L, label = "peakPOS is position of max log10p marker")
})

test_that("define_confidence_intervals produces separate intervals for multiple peaks", {
  df            <- make_mapping_df(30)
  df$significant <- 0L
  df$significant[c(5, 25)] <- 1L   # well-separated: index distance = 20 > snp_grouping=10
  df            <- group_markers_to_intervals(df, snp_grouping = 10, verbose = FALSE)
  result        <- define_confidence_intervals(df, ci_size = 2, verbose = FALSE)
  n_peaks       <- length(unique(na.omit(result$peak_id)))
  expect_equal(n_peaks, 2L, label = "two distinct peaks produce two intervals")
  start_positions <- result$startPOS[result$significant == 1]
  expect_equal(length(unique(start_positions)), 2L,
               label = "each peak has a different startPOS")
})


# =============================================================================
# 4. analyze_mapping() — regression tests for removed guard + integration
# =============================================================================

# --- Core regression: removed check_significant_proportion() guard ---

test_that("analyze_mapping returns intervals when >15% of markers are significant [regression]", {
  # Before fix: check_significant_proportion() returned FALSE → all-NA intervals
  # After fix:  always proceeds to group_markers_to_intervals() → real intervals
  df       <- make_mapping_df(20)
  df$P[1:10] <- 1e-10   # 50% significant, well above old 15% threshold
  result   <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  expect_false(all(is.na(result$peak_id)),
               label = "peak_id is not all-NA even with 50% significant markers")
  expect_false(all(is.na(result$startPOS)),
               label = "startPOS is not all-NA even with 50% significant markers")
  expect_false(all(is.na(result$endPOS)),
               label = "endPOS is not all-NA even with 50% significant markers")
})

test_that("analyze_mapping returns intervals when 100% of markers are significant [regression]", {
  df       <- make_mapping_df(20)
  df$P     <- 1e-10    # all significant
  result   <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  expect_false(all(is.na(result$peak_id)),
               label = "peak_id assigned when all markers are significant")
  # All on one chrom with default snp_grouping → one peak
  expect_equal(length(unique(na.omit(result$peak_id))), 1L,
               label = "all sig markers on one chrom collapse to one peak")
})

# --- Unchanged behavior: 0 significant markers still returns all-NA ---

test_that("analyze_mapping returns all-NA intervals when no markers are significant [unchanged]", {
  df     <- make_mapping_df(20)   # all p=0.5, log10p ≈ 0.3, below threshold=5
  result <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  expect_true(all(is.na(result$peak_id)),
              label = "no sig markers → all-NA peak_id")
  expect_true(all(is.na(result$startPOS)),
              label = "no sig markers → all-NA startPOS")
})

# --- Unchanged behavior: <15% significant still works ---

test_that("analyze_mapping defines intervals when <15% of markers are significant [unchanged]", {
  df       <- make_mapping_df(20)
  df$P[2]  <- 1e-10    # 1/20 = 5%, below old 15% threshold
  result   <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  n_peaks  <- length(unique(na.omit(result$peak_id)))
  expect_equal(n_peaks, 1L, label = "one peak for one sig marker")
  sig_start <- result$startPOS[result$significant == 1][1]
  expect_false(is.na(sig_start), label = "startPOS defined for the sig marker")
})

# --- Output column completeness ---

test_that("analyze_mapping always produces the full set of output columns", {
  df       <- make_mapping_df(10)
  df$P[3]  <- 1e-10
  result   <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  required <- c(
    "log10p", "sig_threshold_value", "sig_threshold_method", "significant",
    "peak_id", "startPOS", "peakPOS", "endPOS", "interval_size"
  )
  expect_true(all(required %in% names(result)), label = "all required output columns present")
})

test_that("analyze_mapping records threshold value and method in every output row", {
  df       <- make_mapping_df(10)
  df$P[3]  <- 1e-10
  result   <- analyze_mapping(df, threshold_value = 5, threshold_method = "BF", verbose = FALSE)
  expect_equal(unique(result$sig_threshold_value), 5, label = "threshold value in every row")
  expect_equal(unique(result$sig_threshold_method), "BF", label = "threshold method in every row")
})

# --- Multi-chromosome behavior ---

test_that("analyze_mapping assigns peaks per-chromosome, not across chromosomes", {
  df       <- two_chrom_df(10)
  df$P[3]  <- 1e-10   # chrom I, position 3000
  df$P[13] <- 1e-10   # chrom II, position 3000
  result   <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  sig_rows <- result[result$significant == 1, ]
  expect_equal(nrow(sig_rows), 2L, label = "two significant markers found")
  expect_equal(length(unique(sig_rows$peak_id)), 2L,
               label = "one peak per chromosome")
  expect_true(setequal(sig_rows$CHROM, c("I", "II")),
              label = "peaks on correct chromosomes")
})

test_that("analyze_mapping with >15% sig on one chrom and 0% on another works correctly", {
  df        <- two_chrom_df(10)
  # Make 80% of chrom I significant, chrom II untouched
  chrom_i   <- df$CHROM == "I"
  df$P[chrom_i][1:8] <- 1e-10
  result    <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  chrom_i_peaks <- unique(na.omit(result$peak_id[result$CHROM == "I"]))
  chrom_ii_peaks <- result$peak_id[result$CHROM == "II"]
  expect_gte(length(chrom_i_peaks), 1L,
             label = "at least one peak on heavily-significant chromosome")
  expect_true(all(is.na(chrom_ii_peaks)),
              label = "no peaks on chromosome with no significant markers")
})

# --- Edge cases ---

test_that("analyze_mapping treats p=0 markers as significant (log10p capped at 300)", {
  df       <- make_mapping_df(10)
  df$P[5]  <- 0
  result   <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  expect_equal(result$significant[5], 1L, label = "p=0 marker is significant")
  expect_false(is.na(result$peak_id[5]), label = "p=0 marker receives a peak_id")
})

test_that("analyze_mapping with dense significant region collapses to one peak per default snp_grouping", {
  df        <- make_mapping_df(100)
  df$P[25:75] <- 1e-10   # 51 consecutive significant markers (51% of 100)
  result    <- analyze_mapping(df, threshold_value = 5, snp_grouping = 1000, verbose = FALSE)
  n_peaks   <- length(unique(na.omit(result$peak_id)))
  expect_equal(n_peaks, 1L,
               label = "dense consecutive sig region collapses to one peak")
})

test_that("analyze_mapping with two separated significant clusters produces two peaks", {
  df        <- make_mapping_df(100)
  df$P[c(5, 95)] <- 1e-10   # two markers far apart (index distance 90 > snp_grouping=10)
  result    <- analyze_mapping(df, threshold_value = 5, snp_grouping = 10, verbose = FALSE)
  n_peaks   <- length(unique(na.omit(result$peak_id)))
  expect_equal(n_peaks, 2L,
               label = "two well-separated sig markers produce two peaks")
})

test_that("analyze_mapping CI extends ci_size markers from peak boundary", {
  df        <- make_mapping_df(30)
  df$P[15]  <- 1e-10    # single sig marker at index 15, position 15000
  result    <- analyze_mapping(df, threshold_value = 5, ci_size = 3, snp_grouping = 1000,
                               verbose = FALSE)
  sig_row   <- result[result$significant == 1, ]
  # startPOS = position at index 15-3=12 = 12000; endPOS = position at index 15+3=18 = 18000
  expect_equal(sig_row$startPOS[1], 12000L,
               label = "startPOS is ci_size markers left of sig marker")
  expect_equal(sig_row$endPOS[1], 18000L,
               label = "endPOS is ci_size markers right of sig marker")
})

test_that("analyze_mapping does not drop original input columns", {
  df          <- make_mapping_df(10)
  df$P[3]     <- 1e-10
  df$extra_col <- "keep_me"
  result      <- analyze_mapping(df, threshold_value = 5, verbose = FALSE)
  expect_true("extra_col" %in% names(result), label = "extra input columns are preserved")
  expect_true("marker" %in% names(result),    label = "marker column is preserved")
})
