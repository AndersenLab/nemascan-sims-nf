# test-assessment.R - Unit tests for R/assessment.R functions

# =============================================================================
# load_causal_variants()
# =============================================================================

test_that("load_causal_variants reads .par file correctly", {
  cv <- load_causal_variants(fixture_path("test_sims.par"))

  expect_s3_class(cv, "data.frame")
  expect_true(nrow(cv) == 4)

  # Check required columns
  expect_true(all(c("QTL", "CHROM", "POS", "RefAllele", "Frequency", "Effect") %in% names(cv)))

  # QTL should be in CHROM:POS format
 expect_equal(cv$QTL[1], "2:7459092")

  # CHROM and POS parsed correctly
  expect_equal(cv$CHROM[1], "2")
  expect_equal(cv$POS[1], 7459092L)

  # Types
  expect_type(cv$CHROM, "character")
  expect_type(cv$POS, "integer")
  expect_type(cv$Effect, "double")
})

test_that("load_causal_variants stops on missing file", {
  expect_error(
    load_causal_variants("/nonexistent/path/test.par"),
    "not found"
  )
})


# =============================================================================
# assess_qtl_detection()
# =============================================================================

test_that("assess_qtl_detection finds overlapping causal variants", {
  causal <- data.frame(
    QTL = c("2:100", "3:500"),
    CHROM = c("2", "3"),
    POS = c(100L, 500L),
    RefAllele = c("A", "G"),
    Frequency = c(0.1, 0.2),
    Effect = c(1.0, -0.5),
    stringsAsFactors = FALSE
  )

  regions <- data.frame(
    peak_id = c(1L, 2L),
    CHROM = c("2", "4"),
    startPOS = c(50L, 400L),
    peakPOS = c(100L, 450L),
    endPOS = c(150L, 600L),
    peak_marker = c("2:100", "4:450"),
    max_log10p = c(5.0, 3.0),
    n_sig_markers = c(10L, 5L),
    interval_size = c(100L, 200L),
    sig_threshold_value = c(4.0, 4.0),
    sig_threshold_method = c("EIGEN", "EIGEN"),
    stringsAsFactors = FALSE
  )

  result <- assess_qtl_detection(regions, causal)

  expect_s3_class(result, "data.frame")
  # 2:100 is in region on CHROM 2 (50-150), 3:500 is NOT in region on CHROM 4
  expect_equal(nrow(result), 2)

  detected <- result %>% dplyr::filter(!is.na(detected_peak_id))
  undetected <- result %>% dplyr::filter(is.na(detected_peak_id))

  expect_equal(nrow(detected), 1)
  expect_equal(detected$QTL, "2:100")
  expect_equal(nrow(undetected), 1)
  expect_equal(undetected$QTL, "3:500")
})

test_that("assess_qtl_detection handles no QTL regions", {
  causal <- data.frame(
    QTL = c("2:100"),
    CHROM = c("2"),
    POS = c(100L),
    RefAllele = c("A"),
    Frequency = c(0.1),
    Effect = c(1.0),
    stringsAsFactors = FALSE
  )

  empty_regions <- data.frame(
    peak_id = integer(),
    CHROM = character(),
    startPOS = integer(),
    peakPOS = integer(),
    endPOS = integer(),
    peak_marker = character(),
    stringsAsFactors = FALSE
  )

  result <- assess_qtl_detection(empty_regions, causal)

  expect_equal(nrow(result), 1)
  expect_true(is.na(result$detected_peak_id[1]))
})

test_that("assess_qtl_detection handles empty causal variants", {
  regions <- data.frame(
    peak_id = 1L,
    CHROM = "2",
    startPOS = 50L,
    peakPOS = 100L,
    endPOS = 150L,
    peak_marker = "2:100",
    stringsAsFactors = FALSE
  )

  empty_causal <- data.frame(
    QTL = character(),
    CHROM = character(),
    POS = integer(),
    stringsAsFactors = FALSE
  )

  result <- assess_qtl_detection(regions, empty_causal)
  expect_equal(nrow(result), 0)
})


# =============================================================================
# compile_full_assessment()
# =============================================================================

test_that("compile_full_assessment produces correct structure", {
  # Create minimal mapping data with QTL intervals
  mapping_data <- data.frame(
    marker = c("2:50", "2:100", "2:150", "3:500", "4:200"),
    CHROM = c("2", "2", "2", "3", "4"),
    POS = c(50L, 100L, 150L, 500L, 200L),
    P = c(0.001, 0.0001, 0.001, 0.5, 0.5),
    AF1 = c(0.1, 0.1, 0.1, 0.2, 0.3),
    BETA = c(0.5, 1.0, 0.5, 0.01, 0.01),
    significant = c(1L, 1L, 1L, 0L, 0L),
    peak_id = c(1L, 1L, 1L, NA_integer_, NA_integer_),
    startPOS = c(50L, 50L, 50L, NA_integer_, NA_integer_),
    peakPOS = c(100L, 100L, 100L, NA_integer_, NA_integer_),
    endPOS = c(150L, 150L, 150L, NA_integer_, NA_integer_),
    interval_size = c(100L, 100L, 100L, NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )

  qtl_regions <- data.frame(
    peak_id = 1L,
    CHROM = "2",
    startPOS = 50L,
    peakPOS = 100L,
    endPOS = 150L,
    interval_size = 100L,
    n_sig_markers = 3L,
    max_log10p = 4.0,
    peak_marker = "2:100",
    sig_threshold_value = 3.0,
    sig_threshold_method = "EIGEN",
    stringsAsFactors = FALSE
  )

  causal <- data.frame(
    QTL = c("2:100", "3:500"),
    CHROM = c("2", "3"),
    POS = c(100L, 500L),
    RefAllele = c("A", "G"),
    Frequency = c(0.1, 0.2),
    Effect = c(1.0, -0.5),
    stringsAsFactors = FALSE
  )

  params <- list(
    nqtl = 2, rep = 1, h2 = 0.5, maf = 0.05,
    effect = "gamma", population = "testpop",
    algorithm = "LMM-EXACT-INBRED", pca = TRUE,
    threshold_method = "EIGEN"
  )

  result <- compile_full_assessment(mapping_data, qtl_regions, causal, params)

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)

  # Check required columns exist
  expected_cols <- c("QTL", "Simulated", "Detected", "nQTL", "simREP",
                     "h2", "maf", "effect_distribution", "strain_set_id",
                     "algorithm_id")
  for (col in expected_cols) {
    expect_true(col %in% names(result), label = paste("Column", col, "exists"))
  }

  # 2:100 should be Simulated=TRUE, Detected=TRUE (in QTL interval)
  qtl_100 <- result %>% dplyr::filter(QTL == "2:100")
  expect_equal(as.character(qtl_100$Simulated), "TRUE")
  expect_equal(as.character(qtl_100$Detected), "TRUE")

  # 3:500 should be Simulated=TRUE, Detected=FALSE (no interval on CHROM 3)
  qtl_500 <- result %>% dplyr::filter(QTL == "3:500")
  expect_equal(as.character(qtl_500$Simulated), "TRUE")
  expect_equal(as.character(qtl_500$Detected), "FALSE")

  # Simulation metadata columns
  expect_equal(unique(as.character(result$nQTL)), "2")
  expect_equal(unique(as.character(result$strain_set_id)), "testpop")
})

test_that("compile_full_assessment handles no detected QTLs", {
  mapping_data <- data.frame(
    marker = c("2:100", "3:500"),
    CHROM = c("2", "3"),
    POS = c(100L, 500L),
    P = c(0.5, 0.5),
    AF1 = c(0.1, 0.2),
    BETA = c(0.01, 0.01),
    significant = c(0L, 0L),
    peak_id = c(NA_integer_, NA_integer_),
    startPOS = c(NA_integer_, NA_integer_),
    peakPOS = c(NA_integer_, NA_integer_),
    endPOS = c(NA_integer_, NA_integer_),
    interval_size = c(NA_integer_, NA_integer_),
    stringsAsFactors = FALSE
  )

  qtl_regions <- data.frame(
    peak_id = integer(), CHROM = character(),
    startPOS = integer(), peakPOS = integer(), endPOS = integer(),
    peak_marker = character(), stringsAsFactors = FALSE
  )

  causal <- data.frame(
    QTL = c("2:100"),
    CHROM = c("2"),
    POS = c(100L),
    RefAllele = c("A"),
    Frequency = c(0.1),
    Effect = c(1.0),
    stringsAsFactors = FALSE
  )

  params <- list(
    nqtl = 1, rep = 1, h2 = 0.5, maf = 0.05,
    effect = "gamma", population = "testpop",
    algorithm = "LMM-EXACT-INBRED", pca = FALSE,
    threshold_method = "BF"
  )

  result <- compile_full_assessment(mapping_data, qtl_regions, causal, params)

  # All causal variants should be Simulated=TRUE, Detected=FALSE
  expect_true(all(result$Simulated == "TRUE"))
  expect_true(all(result$Detected == "FALSE"))
})

test_that("compile_full_assessment identifies false positives", {
  # Mapping with detected peak on CHROM 4 that does NOT contain the causal variant
  # Causal variant 2:100 IS in the mapping data (as it would be in a real marker set)
  # but the detected peak is on CHROM 4
  mapping_data <- data.frame(
    marker = c("2:100", "4:200", "4:250", "4:300"),
    CHROM = c("2", "4", "4", "4"),
    POS = c(100L, 200L, 250L, 300L),
    P = c(0.5, 0.0001, 0.00001, 0.0001),
    AF1 = c(0.1, 0.3, 0.3, 0.3),
    BETA = c(0.01, 0.5, 1.0, 0.5),
    significant = c(0L, 1L, 1L, 1L),
    peak_id = c(NA_integer_, 1L, 1L, 1L),
    startPOS = c(NA_integer_, 200L, 200L, 200L),
    peakPOS = c(NA_integer_, 250L, 250L, 250L),
    endPOS = c(NA_integer_, 300L, 300L, 300L),
    interval_size = c(NA_integer_, 100L, 100L, 100L),
    stringsAsFactors = FALSE
  )

  qtl_regions <- data.frame(
    peak_id = 1L,
    CHROM = "4",
    startPOS = 200L,
    peakPOS = 250L,
    endPOS = 300L,
    interval_size = 100L,
    n_sig_markers = 3L,
    max_log10p = 5.0,
    peak_marker = "4:250",
    sig_threshold_value = 3.0,
    sig_threshold_method = "EIGEN",
    stringsAsFactors = FALSE
  )

  # Causal variant on CHROM 2 — not within the CHROM 4 detected interval
  causal <- data.frame(
    QTL = c("2:100"),
    CHROM = c("2"),
    POS = c(100L),
    RefAllele = c("A"),
    Frequency = c(0.1),
    Effect = c(1.0),
    stringsAsFactors = FALSE
  )

  params <- list(
    nqtl = 1, rep = 1, h2 = 0.5, maf = 0.05,
    effect = "gamma", population = "testpop",
    algorithm = "LMM-EXACT-INBRED", pca = TRUE,
    threshold_method = "EIGEN"
  )

  result <- compile_full_assessment(mapping_data, qtl_regions, causal, params)

  # Should have the causal variant (simulated, not detected) and
  # the false positive peak (detected, not simulated)
  expect_true(nrow(result) >= 2)

  fp <- result %>% dplyr::filter(Detected == "TRUE", Simulated == "FALSE")
  fn <- result %>% dplyr::filter(Simulated == "TRUE", Detected == "FALSE")

  expect_equal(nrow(fp), 1, label = "one false positive peak")
  expect_equal(fp$QTL, "4:250")
  expect_equal(nrow(fn), 1, label = "one false negative (undetected causal)")
  expect_equal(fn$QTL, "2:100")
})


# =============================================================================
# format_assessment_tsv()
# =============================================================================

test_that("format_assessment_tsv produces correct column set", {
  input <- data.frame(
    QTL = "2:100",
    Simulated = factor("TRUE", levels = c("TRUE", "FALSE")),
    Detected = factor("TRUE", levels = c("TRUE", "FALSE")),
    log10p = 5.0,
    significant = 1L,
    startPOS = 50L,
    peakPOS = 100L,
    endPOS = 150L,
    detected.peak = "2:100",
    interval.log10p = 5.0,
    interval.var.exp = NA_real_,
    interval.Frequency = 0.1,
    peak_id = 1L,
    interval_size = 100L,
    top.hit = TRUE,
    nQTL = "1",
    simREP = "1",
    h2 = "0.5",
    maf = "0.05",
    effect_distribution = "gamma",
    strain_set_id = "testpop",
    algorithm_id = "LMM-EXACT-INBRED_PCA_EIGEN",
    stringsAsFactors = FALSE
  )

  result <- format_assessment_tsv(input)

  expected_cols <- c(
    "QTL", "Simulated", "Detected", "log10p", "significant",
    "startPOS", "peakPOS", "endPOS", "detected.peak",
    "interval.log10p", "interval.var.exp", "interval.Frequency",
    "peak_id", "interval_size", "top.hit",
    "nQTL", "simREP", "h2", "maf", "effect_distribution",
    "strain_set_id", "algorithm_id"
  )

  expect_equal(names(result), expected_cols)
})

test_that("format_assessment_tsv handles empty dataframe", {
  result <- format_assessment_tsv(data.frame())
  expect_equal(nrow(result), 0)
})

test_that("format_assessment_tsv adds missing columns as NA", {
  input <- data.frame(
    QTL = "2:100",
    Simulated = factor("TRUE", levels = c("TRUE", "FALSE")),
    Detected = factor("FALSE", levels = c("TRUE", "FALSE")),
    stringsAsFactors = FALSE
  )

  result <- format_assessment_tsv(input)
  expect_true("startPOS" %in% names(result))
  expect_true(is.na(result$startPOS[1]))
})
