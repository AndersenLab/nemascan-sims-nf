# test-assessment_var_exp.R - Unit tests for compute_var_exp_anova()
#
# Tests ANOVA-based variance explained computation for the DB analysis path.
# Uses long-format DB schema fixtures (test_genotype_matrix_anova.tsv,
# test_phenotypes_anova.tsv, test_causal_anova.par).
#
# Golden values pinned from aov(trait.value ~ allele, data=data_clean) with allele
# as factor, run on the fixture data (4 strains, phenotypes A=2.5, B=3.8, C=2.1, D=4.0):
#   QTL 1:100 (balanced 2/2 split): SSallele=2.56, SSresid=0.10, SST=2.66
#   QTL 2:300 (balanced 2/2 split): SSallele=0.01, SSresid=2.65, SST=2.66
#
# ANOVA SS_allele = Σ_g n_g * (ȳ_g - ȳ)² (between-group sum of squares).
# For a balanced two-group case this coincidentally equals (Σx·y)²/Σ(x²),
# but the correct general formula is the group mean deviation formula.
#
# R library is sourced by tests/testthat/setup.R — do NOT re-source here.

# ── Fixture helpers ────────────────────────────────────────────────────────────

load_geno_fixture <- function() {
  data.table::fread(fixture_path("test_genotype_matrix_anova.tsv"),
    header = TRUE, sep = "\t"
  ) %>%
    as.data.frame() %>%
    dplyr::mutate(CHROM = as.character(CHROM), POS = as.integer(POS))
}

load_pheno_fixture <- function() {
  data.table::fread(fixture_path("test_phenotypes_anova.tsv"),
    header = TRUE, sep = "\t"
  ) %>%
    as.data.frame()
}

load_causal_fixture <- function() {
  load_causal_variants(fixture_path("test_causal_anova.par"))
}

# Golden values pinned from actual aov() output on fixture data:
#   4 strains, phenotypes: A=2.5, B=3.8, C=2.1, D=4.0
#   QTL 1:100: alleles A=-1, B=1, C=-1, D=1 → SSallele=2.56, SSresid=0.10, SST=2.66
#   QTL 2:300: alleles A=-1, B=-1, C=1, D=1 → SSallele=0.01, SSresid=2.65, SST=2.66
GOLDEN_1_100 <- 0.962406015037594   # aov(): 2.56 / 2.660000000000001
GOLDEN_2_300 <- 0.003759398496241   # aov(): 0.01 / 2.66


# ── Test cases ─────────────────────────────────────────────────────────────────

test_that("compute_var_exp_anova: balanced QTL 1:100 matches golden aov() value", {
  geno   <- load_geno_fixture()
  pheno  <- load_pheno_fixture()
  causal <- load_causal_fixture()

  result <- compute_var_exp_anova(geno, pheno, causal)

  row_1 <- result[result$QTL == "1:100", ]
  expect_equal(nrow(row_1), 1L)
  expect_equal(row_1$Simulated.QTL.VarExp, GOLDEN_1_100, tolerance = 1e-6)
})


test_that("compute_var_exp_anova: balanced QTL 2:300 matches golden aov() value", {
  geno   <- load_geno_fixture()
  pheno  <- load_pheno_fixture()
  causal <- load_causal_fixture()

  result <- compute_var_exp_anova(geno, pheno, causal)

  row_2 <- result[result$QTL == "2:300", ]
  expect_equal(nrow(row_2), 1L)
  expect_equal(row_2$Simulated.QTL.VarExp, GOLDEN_2_300, tolerance = 1e-6)
})


test_that("compute_var_exp_anova: no position overlap yields all NA", {
  geno  <- load_geno_fixture()
  pheno <- load_pheno_fixture()

  # Causal variants at positions not in genotype matrix
  causal_miss <- data.frame(
    QTL = c("99:999", "88:888"),
    CHROM = c("99", "88"),
    POS = c(999L, 888L),
    RefAllele = c("A", "G"),
    Frequency = c(0.5, 0.5),
    Effect = c(1.0, -1.0),
    stringsAsFactors = FALSE
  )

  result <- compute_var_exp_anova(geno, pheno, causal_miss)
  expect_true(all(is.na(result$Simulated.QTL.VarExp)),
    label = "all Simulated.QTL.VarExp are NA when no positions overlap")
})


test_that("compute_var_exp_anova: near-monomorphic QTL 3:500 (3:1 allele split) returns NA", {
  # QTL 3:500 in the fixture has alleles A=-1, B=-1, C=-1, D=1 (minority class has 1 strain).
  # The guard `any(class_counts < 2)` fires because the +1 group has only 1 member.
  # This exercises the near-monomorphic path via compute_var_exp_anova() → var_exp_one().
  geno   <- load_geno_fixture()
  pheno  <- load_pheno_fixture()
  causal <- load_causal_fixture()  # now includes 3:500

  result <- compute_var_exp_anova(geno, pheno, causal)

  row_3 <- result[result$QTL == "3:500", ]
  expect_equal(nrow(row_3), 1L)
  expect_true(is.na(row_3$Simulated.QTL.VarExp),
    label = "near-monomorphic QTL 3:500 (3:1 split) returns NA (min-strains-per-class guard)")
})


test_that("compute_var_exp_anova: monomorphic QTL (all same allele) returns NA", {
  pheno <- load_pheno_fixture()

  # Monomorphic: all strains have allele = -1 (only one class — guard triggers)
  geno_mono <- data.frame(
    CHROM  = "1",
    POS    = 100L,
    strain = c("strain_A", "strain_B", "strain_C", "strain_D"),
    allele = c(-1, -1, -1, -1),
    stringsAsFactors = FALSE
  )
  causal_mono <- data.frame(
    QTL = "1:100",
    CHROM = "1",
    POS = 100L,
    RefAllele = "A",
    Frequency = 0.0,
    Effect = 1.5,
    stringsAsFactors = FALSE
  )

  result <- compute_var_exp_anova(geno_mono, pheno, causal_mono)
  expect_equal(nrow(result), 1L)
  expect_true(is.na(result$Simulated.QTL.VarExp),
    label = "monomorphic QTL returns NA (minimum strains per group guard)")
})


test_that("compute_var_exp_anova: all-NA alleles (heterozygote-only) returns NA", {
  pheno <- load_pheno_fixture()

  # All alleles are NA (heterozygotes / missing)
  geno_na <- data.frame(
    CHROM  = "1",
    POS    = 100L,
    strain = c("strain_A", "strain_B", "strain_C", "strain_D"),
    allele = c(NA_real_, NA_real_, NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )
  causal_na <- data.frame(
    QTL = "1:100",
    CHROM = "1",
    POS = 100L,
    RefAllele = "A",
    Frequency = 0.5,
    Effect = 1.5,
    stringsAsFactors = FALSE
  )

  result <- compute_var_exp_anova(geno_na, pheno, causal_na)
  expect_equal(nrow(result), 1L)
  expect_true(is.na(result$Simulated.QTL.VarExp),
    label = "all-NA alleles returns NA (nrow < 2 guard after filtering)")
})


test_that("format_assessment_tsv: absent Simulated.QTL.VarExp column is added as NA", {
  # Verifies that format_assessment_tsv() adds Simulated.QTL.VarExp as NA when
  # causal_variants has no Simulated.QTL.VarExp column (canonical NA fallback in the loop).

  mapping_data <- data.frame(
    marker        = "1:100",
    CHROM         = "1",
    POS           = 100L,
    P             = 0.001,
    log10p        = 3.0,
    significant   = 1L,
    AF1           = 0.5,
    BETA          = 1.5,
    SE            = 0.2,
    peak_id       = NA_integer_,
    startPOS      = NA_integer_,
    peakPOS       = NA_integer_,
    endPOS        = NA_integer_,
    interval_size = NA_integer_,
    stringsAsFactors = FALSE
  )
  qtl_regions <- data.frame(stringsAsFactors = FALSE)

  # causal_variants deliberately has NO Simulated.QTL.VarExp column
  causal_variants_no_varexp <- data.frame(
    QTL        = "1:100",
    CHROM      = "1",
    POS        = 100L,
    RefAllele  = "A",
    Frequency  = 0.5,
    Effect     = 1.5,
    stringsAsFactors = FALSE
  )

  mapping_params <- list(
    population = "test_pop", maf = 0.05, nqtl = 1L, effect = "gamma",
    rep = 1L, h2 = 0.5, algorithm = "inbred", pca = FALSE,
    threshold_method = "BF", mode = "inbred", type = "nopca",
    alpha = 0.05, ci_size = 150L, snp_grouping = 1000L
  )

  result    <- compile_full_assessment(mapping_data, qtl_regions,
                                       causal_variants_no_varexp, mapping_params)
  formatted <- format_assessment_tsv(result)

  expect_true("Simulated.QTL.VarExp" %in% names(formatted),
    label = "Simulated.QTL.VarExp column present even when absent from causal_variants input")
  sim_rows <- formatted[as.character(formatted$Simulated) == "TRUE", ]
  expect_gt(nrow(sim_rows), 0L, label = "at least one simulated QTL row")
  expect_true(all(is.na(sim_rows$Simulated.QTL.VarExp)),
    label = "Simulated.QTL.VarExp is NA when not provided in causal_variants")
})


test_that("compile_full_assessment: Simulated.QTL.VarExp flows from causal_variants to output", {
  # Verify that Simulated.QTL.VarExp in causal_variants is preserved through
  # compile_full_assessment() via the any_of() select in build_assessment_union().

  # Minimal mapping data with one marker at the causal position.
  # significant and log10p are required by score_causal_markers().
  mapping_data <- data.frame(
    marker        = "1:100",
    CHROM         = "1",
    POS           = 100L,
    P             = 0.001,
    log10p        = 3.0,
    significant   = 1L,
    AF1           = 0.5,
    BETA          = 1.5,
    SE            = 0.2,
    peak_id       = NA_integer_,
    startPOS      = NA_integer_,
    peakPOS       = NA_integer_,
    endPOS        = NA_integer_,
    interval_size = NA_integer_,
    stringsAsFactors = FALSE
  )

  qtl_regions <- data.frame(stringsAsFactors = FALSE)  # no detected peaks

  causal_variants <- data.frame(
    QTL        = "1:100",
    CHROM      = "1",
    POS        = 100L,
    RefAllele  = "A",
    Frequency  = 0.5,
    Effect     = 1.5,
    Simulated.QTL.VarExp = 0.962,
    stringsAsFactors = FALSE
  )

  mapping_params <- list(
    population = "test_pop",
    maf = 0.05,
    nqtl = 1L,
    effect = "gamma",
    rep = 1L,
    h2 = 0.5,
    algorithm = "inbred",
    pca = FALSE,
    threshold_method = "BF",
    mode = "inbred",
    type = "nopca",
    alpha = 0.05,
    ci_size = 150L,
    snp_grouping = 1000L
  )

  result    <- compile_full_assessment(
    mapping_data    = mapping_data,
    qtl_regions     = qtl_regions,
    causal_variants = causal_variants,
    mapping_params  = mapping_params
  )
  formatted <- format_assessment_tsv(result)

  expect_true("Simulated.QTL.VarExp" %in% names(formatted),
    label = "Simulated.QTL.VarExp column present in formatted output")
  sim_rows <- formatted[as.character(formatted$Simulated) == "TRUE", ]
  expect_gt(nrow(sim_rows), 0L, label = "at least one simulated QTL row")
  expect_equal(sim_rows$Simulated.QTL.VarExp[1], 0.962, tolerance = 1e-6,
    label = "Simulated.QTL.VarExp value preserved through compile_full_assessment()")
})


# ==============================================================================
# T4 — Non-marker causal variant var.exp tests
# ==============================================================================

test_that("compute_var_exp_anova: non-marker position computes var.exp when genotype provided", {
  pheno <- load_pheno_fixture()
  # Synthetic genotype at non-marker position 9:999 — same allele pattern as 1:100
  # (A=-1, B=1, C=-1, D=1) so it should yield the same golden var.exp as 1:100
  geno_nonmarker <- data.frame(
    CHROM  = "9",
    POS    = 999L,
    strain = c("strain_A", "strain_B", "strain_C", "strain_D"),
    allele = c(-1, 1, -1, 1),
    stringsAsFactors = FALSE
  )
  causal_9 <- data.frame(
    QTL       = "9:999",
    CHROM     = "9",
    POS       = 999L,
    RefAllele = "A",
    Frequency = 0.5,
    Effect    = 1.5,
    stringsAsFactors = FALSE
  )
  result <- compute_var_exp_anova(geno_nonmarker, pheno, causal_9)
  expect_equal(result$Simulated.QTL.VarExp, GOLDEN_1_100, tolerance = 1e-6,
               label = "non-marker position with supplied genotype matches expected var.exp")
})

test_that("compute_var_exp_anova: non-marker without genotype still returns NA", {
  # When only the marker genotype matrix is passed (no per-trait merge),
  # a non-marker position is absent from the geno data → inner_join yields NA
  result <- compute_var_exp_anova(
    load_geno_fixture(),
    load_pheno_fixture(),
    data.frame(
      QTL       = "9:999",
      CHROM     = "9",
      POS       = 999L,
      RefAllele = "A",
      Frequency = 0.5,
      Effect    = 1.0,
      stringsAsFactors = FALSE
    )
  )
  expect_true(is.na(result$Simulated.QTL.VarExp),
    label = "non-marker position absent from geno matrix returns NA")
})
