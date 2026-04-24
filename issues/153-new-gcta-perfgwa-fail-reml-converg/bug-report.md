

# Bug Report: GCTA-REML fails to converge in GWAS analysis

## Description

Several process failures for two NF processes `GCTA_MAKE GRM` and `GCTA_PERFORM_GWA` 

Pipeline failed after reaching the max number of errors for the `GCTA_PERFORM_GWA` process for a particular input parameter combination. The error message indicates that the fastGWA-REML method used in GCTA is unable to converge for the given data and parameters.

## Command Used

*(rockfish production run; full command not captured — add from .nextflow.log)*

## Error output

### Execution Summary

```
[e9/5c2619] Submitted process > GCTA_PERFORM_GWA (5 48 0.8 gamma inbred nopca ce.n300.r24_0.05)                                                                                                                                
executor >  local (251), slurm (58236)                                                                                                                                                                                         
[be/d0c7c4] LOCAL_GET_CONTIG_INFO (ce.n100.r01)                                      [100%] 1 of 1 ✔                                                                                                                           
[c8/2bc47a] BCFTOOLS_EXTRACT_STRAINS (ce.n500.r20)                                   [100%] 250 of 250 ✔                                                                                                                       
[1b/7342d4] BCFTOOLS_RENAME_CHROMS (ce.n500.r13)                                     [100%] 250 of 250 ✔                                                                                                                       
[62/581265] PLINK_RECODE_MS_VCF (ce.n500.r20_0.05)                                   [100%] 250 of 250 ✔                                                                                                                       
[e4/579828] PLINK_RECODE_CV_VCF (ce.n500.r20_0.05)                                   [100%] 250 of 250 ✔                                                                                                                       
[95/91729b] BCFTOOLS_CREATE_GENOTYPE_MATRIX (ce.n500.r37_0.05)                       [100%] 250 of 250 ✔                                                                                                                       
[89/250b3a] R_FIND_GENOTYPE_MATRIX_EIGEN (2 ce.n500.r29_0.05)                        [100%] 1500 of 1500 ✔                                                                                                                     
[62/bb3d7e] LOCAL_COMPILE_EIGENS (ce.n500.r20 0.05)                                  [100%] 250 of 250 ✔                                                                                                                       
[11/73ace9] PYTHON_SIMULATE_EFFECTS_GLOBAL (5 28 gammace.n100.r17_0.05)              [100%] 12500 of 12500 ✔                                                                                                                   
[01/fbef31] GCTA_SIMULATE_PHENOTYPES (5 32 0.8 gamma ce.n300.r29_0.05)               [ 85%] 10700 of 12500                                                                                                                     
[c9/96f4ed] PLINK_UPDATE_BY_H2 (5 31 0.8 gamma ce.n200.r50_0.05)                     [ 62%] 6700 of 10700                                                                                                                      
[09/1fbb04] GCTA_MAKE_GRM (5 19 0.8 gamma inbred ce.n100.r37_0.05)                   [ 57%] 7702 of 13404, failed: 4, retries: 4                                                                                               
[fc/cb0551] GCTA_PERFORM_GWA (5 48 0.8 gamma inbred pca ce.n300.r24_0.05)            [ 48%] 7202 of 14804, failed: 5, retries: 4                                                                                               
[b9/7727d7] DB_MIGRATION_WRITE_MARKER_SET (ce.n100.r19_0.05)                         [100%] 250 of 250 ✔                                                                                                                       
[b5/c22e9c] DB_MIGRATION_WRITE_GENOTYPE_MATRIX (ce.n500.r34_0.05)                    [100%] 250 of 250 ✔                                                                                                                       
[3a/3e1b3b] DB_MIGRATION_WRITE_TRAIT_DATA (ce.n100.r41_5_36_0.8)                     [ 71%] 7600 of 10700                                                                                                                      
[33/f97ca8] DB_MIGRATION_WRITE_GWA_TO_DB (5 7 0.8 gamma loco nopca ce.n300.r18_0.05) [ 32%] 2300 of 7000                                                                                                                       
[-        ] DB_MIGRATION_AGGREGATE_METADATA                                          -                                                                                                                                         
[-        ] DB_MIGRATION_ANALYZE_QTL                                                 -                                                                                                                                         
[-        ] DB_MIGRATION_ASSESS_SIMS                                                 -                                                                                                                                        
```

### Caused by

```
  Process `GCTA_PERFORM_GWA (5 32 0.8 gamma inbred nopca ce.n100.r06_0.05)` terminated with an error exit status (1)
```

### Command output:

```
  Overall computational time: 0.05 sec.
  *******************************************************************
  * Genome-wide Complex Trait Analysis (GCTA)
  * version v1.94.1 Linux
  * Built at Nov 15 2022 21:14:25, by GCC 8.5
  * (C) 2010-present, Yang Lab, Westlake University
  * Please report bugs to Jian Yang <jian.yang@westlake.edu.cn>
  *******************************************************************
  Analysis started at 17:25:12 EDT on Wed Apr 22 2026.
  Hostname: c618

  Options:

  --fastGWA-mlm-exact
  --bfile TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06
  --grm-sparse 5_32_0.8_0.05_gamma_ce.n100.r06_sparse_grm_inbred
  --out 5_32_0.8_0.05_gamma_ce.n100.r06_lmm-exact_inbred_nopca
  --pheno 5_32_0.8_0.05_gamma_ce.n100.r06_sims.pheno
  --extract plink_snplist.txt
  --thread-num 1

  The program will be running with up to 1 threads.
  Reading PLINK FAM file from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.fam]...
  100 individuals to be included from FAM file.
  Reading phenotype data from [5_32_0.8_0.05_gamma_ce.n100.r06_sims.pheno]...
  100 overlapping individuals with non-missing data to be included from the phenotype file.
  100 individuals to be included. 0 males, 0 females, 100 unknown.
  Reading PLINK BIM file from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.bim]...
  38597 SNPs to be included from BIM file(s).
  Get 38597 SNPs from list [plink_snplist.txt].
  After extracting SNP, 38597 SNPs remain.
  Reading the sparse GRM file from [5_32_0.8_0.05_gamma_ce.n100.r06_sparse_grm_inbred]...
  After matching all the files, 100 individuals to be included in the analysis.
  Estimating the genetic variance (Vg) by fastGWA-REML (grid search)...
  Iteration 1, step size: 1.8458, logL: -245.408. Vg: 114.44, searching range: 112.594 to 116.286
  Iteration 2, step size: 0.246107, logL: -244.137. Vg: 115.547, searching range: 115.301 to 115.793
  Iteration 3, step size: 0.0328143, logL: -243.816. Vg: 115.629, searching range: 115.596 to 115.662
  Iteration 4, step size: 0.00437524, logL: -243.764. Vg: 115.649, searching range: 115.645 to 115.653
  Iteration 5, step size: 0.000583365, logL: -243.764. Vg: 115.648, searching range: 115.647 to 115.648
  Iteration 6, step size: 7.7782e-05, logL: -243.764. Vg: 115.648, searching range: 115.647 to 115.648
  Iteration 7, step size: 1.03709e-05, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Iteration 8, step size: 1.38279e-06, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Iteration 9, step size: 1.84372e-07, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Iteration 10, step size: 2.45829e-08, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Iteration 11, step size: 3.27773e-09, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Iteration 12, step size: 4.3703e-10, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Iteration 13, step size: 5.82702e-11, logL: -243.764. Vg: 115.648, searching range: 115.648 to 115.648
  Best guess Vg range: 115.647573108802 to 115.647573108919, Vp: 115.362665270647
  Error: fastGWA-REML can't converge.
```

## Analysis

### Error classification

This is a **different failure mode from issue #143**. Issue #143 was `"the Vp is below 1e-5"` (phenotype variance too small for GCTA's internal threshold, fixed by iterative Vp scaling in `GCTA_MAKE_GRM`). The current error is `"fastGWA-REML can't converge"`.

### Root cause: Vg > Vp (h² > 1)

The GCTA log reveals the underlying condition:

```
Best guess Vg range: 115.647573..., Vp: 115.362665...
```

**Vg (115.648) > Vp (115.363)** — the REML grid search converges numerically (step size reaches ~6e-11 by iteration 13) but the MLE solution has estimated genetic variance exceeding total phenotypic variance, implying h² ≈ 1.0025. Since h² must lie in [0, 1], the solution is outside the valid parameter space. GCTA detects this and hard-fails rather than clamping to the boundary (h² = 1).

### Retry behavior

The execution summary shows `failed: 4, retries: 4` for `GCTA_MAKE_GRM` and `failed: 5, retries: 4` for `GCTA_PERFORM_GWA`. The `rockfish.config` `maxRetries = 3` setting was designed for transient SLURM failures. This GCTA convergence failure is **deterministic** — the same phenotype with the same sparse GRM will always fail. All retries are wasted compute.

---

## Debugging steps

The following commands should be run from the pipeline launch directory on the cluster (where `.nextflow.log` and the `nextflow` binary are available). Use `nextflow log` to retrieve specific task work directories — do not search scratch space broadly.

### Step 0: Identify the run name

Activate the conda environment used for the run to ensure `nextflow` is available

```bash
conda activate /data/eande106/software/conda_envs/nf24_env
```

```bash
nextflow log
```

```
2026-04-22 12:51:50     4h 38m 23s      spontaneous_swirles     ERR     04477aefd0      39a52330-c538-40ee-9c9d-b13ea7d406be    nextflow run main.nf -profile rockfish --strainfile data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/strains_panelsize_ce.tsv --nqtl data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/nqtl.csv --h2 data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/h2.csv --effect data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/effect_sizes.csv --reps 50 --cv_maf 0.05 --cv_ld 0.8 --output_dir Sims_panelsize_ce_nqtl5_h2-0.8_50reps -work-dir /scratch4/eande106/Ryan_tmp/nf-work-panelsize
```
The run name is `spontaneous_swirles`.

### Step 1: Confirm the error is deterministic across retries

Retrieve all failed `GCTA_PERFORM_GWA` task work directories and check that every attempt shows the same convergence failure (not a SLURM OOM or timeout):

```bash
nextflow log spontaneous_swirles \
  -filter 'process == "GCTA_PERFORM_GWA" && status == "FAILED"' \
  -f name,workdir,exit,attempt \
  | head -40
```

```
GCTA_PERFORM_GWA (5 32 0.8 gamma inbred nopca ce.n100.r06_0.05) /scratch4/eande106/Ryan_tmp/nf-work-panelsize/f3/747893c194d638e6a37318ff11f5d4 1       1
GCTA_PERFORM_GWA (5 32 0.8 gamma inbred nopca ce.n100.r06_0.05) /scratch4/eande106/Ryan_tmp/nf-work-panelsize/9b/29062ca47de46d427a0d1b5a72be8c 1       2
GCTA_PERFORM_GWA (5 32 0.8 gamma inbred nopca ce.n100.r06_0.05) /scratch4/eande106/Ryan_tmp/nf-work-panelsize/9b/a1a1c1301ffef5b06d739e5d0c35c3 1       3
GCTA_PERFORM_GWA (5 32 0.8 gamma inbred pca ce.n100.r06_0.05)   /scratch4/eande106/Ryan_tmp/nf-work-panelsize/70/cebf607bfd422b62fe128015d5591e 1       1
GCTA_PERFORM_GWA (5 32 0.8 gamma inbred nopca ce.n100.r06_0.05) /scratch4/eande106/Ryan_tmp/nf-work-panelsize/91/5765b40e98558fb17f709f97d61486 1       4
```

Check each failed task's `.command.out`  by iterating over the `nextflow log` output above:

The loop:

1. Re-runs the same nextflow log query.
2. Reads each output line.
3. Extracts the workdir field using awk from the end of the line.
4. Greps the exact convergence error string from each `.command.out` 


Lets check the files as well, since the GCTA output may be there instead 

```bash
nextflow log spontaneous_swirles \
  -filter 'process == "GCTA_PERFORM_GWA" && status == "FAILED"' \
  -f name,workdir,exit,attempt \
  | head -40 \
  | while IFS= read -r line; do
      workdir="$(awk '{print $(NF-2)}' <<< "$line")"
      echo "=== $workdir ==="
      grep -F "Error: fastGWA-REML can't converge." "$workdir/.command.out" || true
      echo
    done
```

```
=== /scratch4/eande106/Ryan_tmp/nf-work-panelsize/f3/747893c194d638e6a37318ff11f5d4 ===
Error: fastGWA-REML can't converge.

=== /scratch4/eande106/Ryan_tmp/nf-work-panelsize/9b/29062ca47de46d427a0d1b5a72be8c ===
Error: fastGWA-REML can't converge.

=== /scratch4/eande106/Ryan_tmp/nf-work-panelsize/9b/a1a1c1301ffef5b06d739e5d0c35c3 ===
Error: fastGWA-REML can't converge.

=== /scratch4/eande106/Ryan_tmp/nf-work-panelsize/70/cebf607bfd422b62fe128015d5591e ===
Error: fastGWA-REML can't converge.

=== /scratch4/eande106/Ryan_tmp/nf-work-panelsize/91/5765b40e98558fb17f709f97d61486 ===
Error: fastGWA-REML can't converge.
```

Every retry attempt shows `"fastGWA-REML can't converge"` and the failure is deterministic.

### Step 2: Confirm the corresponding GCTA_MAKE_GRM completed with Vp well above threshold

For the same parameter combination (`5 32 0.8 gamma inbred ce.n100.r06_0.05`), find the matching `GCTA_MAKE_GRM` work directory and inspect its REML output:

```bash
nextflow log spontaneous_swirles \
  -filter 'process == "GCTA_MAKE_GRM" && tag == "5 32 0.8 gamma inbred ce.n100.r06_0.05"' \
  -f tag,workdir,status
```

Output:

```
5 32 0.8 gamma inbred ce.n100.r06_0.05  /scratch4/eande106/Ryan_tmp/nf-work-panelsize/d2/0fddf7b2eec655177227293c460c5e COMPLETED
```

Now check the `check_vp.hsq` file in that work directory, which contains the REML estimates from `GCTA_MAKE_GRM`:

```bash
cat /scratch4/eande106/Ryan_tmp/nf-work-panelsize/d2/0fddf7b2eec655177227293c460c5e/check_vp.hsq
```

Output:

```
Source  Variance        SE
V(G)    107.446470      16.102063
V(e)    0.000115        1.189889
Vp      107.446585      15.850732
V(G)/Vp 0.999999        0.011074
logL    -246.057
logL0   -286.833
LRT     81.551
df      1
Pval    0.0000e+00
n       100
```

E`GCTA_MAKE_GRM` completed successfully and `check_vp.hsq` shows Vp well above 1e-4 (e.g. Vp = 107.446585). This rules out the issue #143 Vp-too-small failure mode and confirms the two errors are distinct. The `GCTA_MAKE_GRM` REML converged to a valid solution with h² ≈ 1, but the `GCTA_PERFORM_GWA` REML with the sparse GRM fails to converge because it finds an invalid solution with h² > 1.

#### Step 2a: Check if the guard against issue #143 is working as intended for this parameter combination

Check the `GCTA_MAKE_GRM` logs for this parameter combination to confirm that the Vp scaling loop is working as intended (it should not have triggered since Vp was already large):

Check the `.command.out` for the `GCTA_MAKE_GRM` task for this parameter combination:

```bash
cat /scratch4/eande106/Ryan_tmp/nf-work-panelsize/d2/0fddf7b2eec655177227293c460c5e/.command.out
```

Output:

```
*******************************************************************
* Genome-wide Complex Trait Analysis (GCTA)
* version v1.94.1 Linux
* Built at Nov 15 2022 21:14:25, by GCC 8.5
* (C) 2010-present, Yang Lab, Westlake University
* Please report bugs to Jian Yang <jian.yang@westlake.edu.cn>
*******************************************************************
Analysis started at 15:28:05 EDT on Wed Apr 22 2026.
Hostname: c640

Accepted options:
--bfile TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06
--autosome
--extract plink_snplist.txt
--make-grm-inbred
--out TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred
--thread-num 1

Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to speed up the computation if there are multiple processors in your machine.

Reading PLINK FAM file from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.fam].
100 individuals to be included from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.fam].
Reading PLINK BIM file from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.bim].
38597 SNPs to be included from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.bim].
Reading a list of SNPs from [plink_snplist.txt].
38597 SNPs are extracted from [plink_snplist.txt].
38597 SNPs from chromosome 1 to chromosome 22 are included in the analysis.
Reading PLINK BED file from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.bed] in SNP-major format ...
Genotype data for 100 individuals and 38597 SNPs to be included from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06.bed].
Calculating allele frequencies ...
Recoding genotypes (individual major mode) ...

Calculating the genetic relationship matrix (GRM) ... (Note: default speed-optimized mode, may use huge RAM)

Summary of the GRM:
Mean of diagonals = 2
Variance of diagonals = 1.11747
Mean of off-diagonals = -0.0202019
Variance of off-diagonals = 0.220803
GRM of 100 individuals has been saved in the file [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.bin] (in binary format).
Number of SNPs to calculate the genetic relationship between each pair of individuals has been saved in the file [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.N.bin] (in binary format).
IDs for the GRM file [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.bin] have been saved in the file [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.id].

Analysis finished at 15:28:06 EDT on Wed Apr 22 2026
Overall computational time: 0.26 sec.
*******************************************************************
* Genome-wide Complex Trait Analysis (GCTA)
* version v1.94.1 Linux
* Built at Nov 15 2022 21:14:25, by GCC 8.5
* (C) 2010-present, Yang Lab, Westlake University
* Please report bugs to Jian Yang <jian.yang@westlake.edu.cn>
*******************************************************************
Analysis started at 15:28:06 EDT on Wed Apr 22 2026.
Hostname: c640

Accepted options:
--grm TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred
--pheno working_pheno.txt
--reml
--out check_vp
--thread-num 1

Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to speed up the computation if there are multiple processors in your machine.

Reading IDs of the GRM from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.id].
100 IDs are read from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.id].
Reading the GRM from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.bin].
GRM for 100 individuals are included from [TO_SIMS_5_32_0.8_0.05_gamma_ce.n100.r06_gcta_grm_inbred.grm.bin].
Reading phenotypes from [working_pheno.txt].
Non-missing phenotypes of 100 individuals are included from [working_pheno.txt].

100 individuals are in common in these files.

Performing  REML analysis ... (Note: may take hours depending on sample size).
100 observations, 1 fixed effect(s), and 2 variance component(s)(including residual variance).
Calculating prior values of variance components by EM-REML ...
Updated prior values: 53.7589 39.5189
logL: -258.533
Running AI-REML algorithm ...
Iter.   logL    V(G)    V(e)
1       -253.48 67.24130        22.01094
2       -249.66 76.49509        14.18447
3       -248.39 83.21566        9.63918
4       -247.81 88.59034        6.55820
5       -247.47 102.95575       0.00012 (1 component(s) constrained)
6       -246.01 107.25889       0.00012 (1 component(s) constrained)
7       -246.05 107.31806       0.00012 (1 component(s) constrained)
8       -246.05 107.44632       0.00012 (1 component(s) constrained)
9       -246.06 107.44647       0.00012 (1 component(s) constrained)
Log-likelihood ratio converged.

Calculating the logLikelihood for the reduced model ...
(variance component 1 is dropped from the model)
Calculating prior values of variance components by EM-REML ...
Updated prior values: 115.36267
logL: -286.83258
Running AI-REML algorithm ...
Iter.   logL    V(e)
1       -286.83 115.36267
Log-likelihood ratio converged.

Summary result of REML analysis:
Source  Variance        SE
V(G)    107.446470      16.102063
V(e)    0.000115        1.189889
Vp      107.446585      15.850732
V(G)/Vp 0.999999        0.011074

Sampling variance/covariance of the estimates of variance components:
2.592764e+02    -4.723275e+00
-4.723275e+00   1.415835e+00

Summary result of REML analysis has been saved in the file [check_vp.hsq].

Analysis finished at 15:28:06 EDT on Wed Apr 22 2026
Overall computational time: 0.02 sec.
Vp=107.446585 >= 0.0001 at round 1; done
```

It looks like no Vp scaling iterations were triggered, which is expected since Vp was already large. This confirms that the issue #143 guard is working as intended and that this convergence failure is a distinct issue.

**However:** The REML analysis identies the same h² ≈ 1 solution (V(G)/Vp = 0.999999) as the `GCTA_PERFORM_GWA` REML

### Step 3: Characterize which parameter combinations fail

Check whether failures are concentrated in small-n, high-h² inbred tasks (and absent from loco tasks with the same parameters):

```bash
# All GCTA_PERFORM_GWA failures with their task names
nextflow log spontaneous_swirles \
  -filter 'process == "GCTA_PERFORM_GWA" && status == "FAILED"' \
  -f name \
  | sort | uniq -c | sort -rn | head -40
```

Output:
```
4 GCTA_PERFORM_GWA (5 32 0.8 gamma inbred nopca ce.n100.r06_0.05)
1 GCTA_PERFORM_GWA (5 32 0.8 gamma inbred pca ce.n100.r06_0.05)
```

---

## Proposed fix direction

The fix must handle the case where `--fastGWA-mlm-exact` with `--grm-sparse` converges to h² > 1. Two options:

**Option A — Graceful skip with empty output (preferred):** When GCTA fails with `"fastGWA-REML can't converge"`, 

Could write a header-only `.fastGWA` file and exit 0. Downstream `DB_MIGRATION_WRITE_GWA_TO_DB` would ingest an empty mapping — zero association statistics for this trait. The pipeline continues; downstream QTL assessment simply finds no QTL for this replicate.
