# Bug Report

2026-03-20

# Bug Report: LOCO and INBRED mapping modes produce different marker counts from identical PLINK input

**Date:** 2026-03-20 **Nextflow version:** version 24.10.4, build 5934
(20-01-2025 16:47 UTC) **Executor:** local

------------------------------------------------------------------------

## Description

`GCTA_PERFORM_GWA` produces a different number of tested markers
depending on mapping mode, even though both modes receive the same
PLINK-filtered input files. With `--mlma-loco` (LOCO mode), GCTA
re-computes allele frequencies internally and excludes 414 markers that
PLINK had already passed at the same MAF threshold (0.05), resulting in
6602 markers tested. With `--fastGWA-mlm-exact` (INBRED mode), GCTA
trusts the PLINK pre-filtering and tests all 7016 markers.

This discrepancy was the root cause of the interval concordance failure
documented in `issues/114-integration-cross-validation-failures/` and
resolved in PR \#121 (branch `mckeowr1/issue111`). The fix applied there
(`nrow(mapping_data)` instead of `marker_set_metadata$n_markers` as the
BF threshold denominator) correctly accounts for this behaviour at the
analysis level. This report documents the underlying GCTA behaviour for
reference.

## Command Used

``` bash
nextflow run main.nf -profile test,docker
```

The discrepancy is observed in the `GCTA_PERFORM_GWA` process. Work
directories for the integration test run:

| Mode   | Type  | Workdir                                            |
|--------|-------|----------------------------------------------------|
| LOCO   | noPCA | `tests/.nf-work/19/65278cd3e96ef555dfcfb8752c1c5d` |
| LOCO   | PCA   | `tests/.nf-work/65/56ba6301af38c6c02bf263261bec38` |
| INBRED | noPCA | `tests/.nf-work/da/644766a1e8219a4d3cdab5090c57f9` |
| INBRED | PCA   | `tests/.nf-work/b0/65c6763c4c11287b7a5eff796d74ff` |

All four runs used: population `ce.test.200strains`, nqtl=5, rep=1,
h2=0.8, maf=0.05, effect=gamma.

## Error Output

Both modes receive identical PLINK BED/BIM/FAM files produced by
`PLINK_RECODE_VCF` (`modules/plink/recode_vcf/main.nf`) with
`--maf 0.05`. The BIM file contains 7016 markers.

**LOCO log**
(`tests/.nf-work/19/65278cd3e96ef555dfcfb8752c1c5d/.command.log`, lines
52–58):

    Reading PLINK BIM file from [TO_SIMS_5_1_0.8_0.05_gamma_ce.test.200strains.bim].
    7016 SNPs to be included from [TO_SIMS_5_1_0.8_0.05_gamma_ce.test.200strains.bim].
    Reading PLINK BED file from [TO_SIMS_5_1_0.8_0.05_gamma_ce.test.200strains.bed] in SNP-major format ...
    Genotype data for 200 individuals and 7016 SNPs to be included from [TO_SIMS_5_1_0.8_0.05_gamma_ce.test.200strains.bed].
    Calculating allele frequencies ...
    Filtering SNPs with MAF > 0.05 ...
    After filtering SNPs with MAF > 0.05, there are 6602 SNPs (414 SNPs with MAF < 0.05).

Per-chromosome breakdown after LOCO filtering (lines 65–76):

    Chr 1:  1357 SNPs on chromosome 1 are included in the analysis.
    Chr 2:  3162 SNPs on chromosome 2 are included in the analysis.
    Chr 3:  2083 SNPs on chromosome 3 are included in the analysis.

**INBRED log**
(`tests/.nf-work/da/644766a1e8219a4d3cdab5090c57f9/.command.log`, lines
55–96):

    Reading PLINK BIM file from [TO_SIMS_5_1_0.8_0.05_gamma_ce.test.200strains.bim]...
    7016 SNPs to be included from BIM file(s).
    Threshold to filter variants: MAF > 0.050000.
    Reading the sparse GRM file from [5_1_0.8_0.05_gamma_ce.test.200strains_sparse_grm_inbred]...
    After matching all the files, 200 individuals to be included in the analysis.
    ...
    7016 SNPs have been processed.
    Saved 7016 SNPs.

**Marker count summary:**

| Source | Marker count | Notes |
|----|:--:|----|
| PLINK `.bim` file | 7016 | After `--maf 0.05` in `PLINK_RECODE_VCF` |
| INBRED (fastGWA) output | 7016 | No additional filtering |
| LOCO (mlma-loco) output | 6602 | 414 markers re-filtered by GCTA |

## Config / Profile

Both mapping modes ran with the `test,docker` profile. The `--maf`
argument passed to `gcta64` in `GCTA_PERFORM_GWA` is the same value used
by `PLINK_RECODE_VCF`:

``` groovy
// modules/gcta/perform_gwa/main.nf (line 61)
--maf ${maf}
```

``` groovy
// modules/plink/recode_vcf/main.nf (lines 30, 40)
--maf ${maf}
```

## Additional Context

### Why the filtering differs between GCTA commands

`--mlma-loco` performs its own allele frequency calculation from the raw
genotype matrix before association testing. The GCTA log explicitly
shows the steps: “Calculating allele frequencies…” followed by
“Filtering SNPs with MAF \> 0.05”. This recalculation is internal to
`--mlma-loco` and produces slightly different per-marker MAF estimates
than PLINK’s calculation — enough to reclassify 414 markers as falling
below the 0.05 boundary.

`--fastGWA-mlm-exact` reports the MAF threshold (“Threshold to filter
variants: MAF \> 0.050000”) but does not perform a separate frequency
calculation step on the genotype matrix. It processes all 7016 markers
that were passed in.

This appears to be an inherent difference in how the two GCTA algorithms
handle the MAF parameter: LOCO re-derives frequencies to guard against
cross-chromosome GRM mismatches, while fastGWA trusts the upstream
filtering.

### Downstream impact on BF threshold

The Bonferroni threshold is `−log10(α / n_tests)`. When `n_tests` is
drawn from `marker_set_metadata$n_markers` (7016, from the `.bim` file),
the DB threshold is higher than the legacy threshold (which uses
`sum(log10p > 0)` = actual GWA output row count = 6602):

``` r
bf_legacy <- -log10(0.05 / 6602)   # 5.1207
bf_db     <- -log10(0.05 / 7016)   # 5.1471
```

The 0.026 difference excluded the borderline marker `3:16538229` (log10p
= 5.123) in the DB path, shifting the `startPOS` of the `loco_nopca_bf`
interval by 79,774 bp.

### Fix applied (PR \#121, branch `mckeowr1/issue111`)

`analyze_qtl.R` and `assess_sims.R` now use `nrow(mapping_data)` as the
BF denominator instead of `threshold_params$n_markers`. This correctly
matches the legacy path’s `sum(log10p > 0)` because both equal the
actual number of rows in the GWA output file (6602 for LOCO, 7016 for
INBRED — which happen to agree with `.bim` for INBRED, making INBRED
unaffected by the fix).

------------------------------------------------------------------------

*Generated with the Quarto GitHub Issue Framework. Submit via
`gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
