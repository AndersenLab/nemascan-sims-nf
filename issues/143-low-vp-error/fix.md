# Fix: Iterative Vp Scaling — Add loop in GCTA_MAKE_GRM & Remove PYTHON_CHECK_VP
 
## Context
 
GCTA's `--fastGWA-mlm-exact` (inbred mode) hard-fails when internal Vp < 1e-5
(error: `"the Vp is below 1e-5"`). The current single-pass guard
(`PYTHON_CHECK_VP` / `bin/check_vp.py`) scales phenotypes x1000 once if
full-GRM REML Vp < 1e-4. This is insufficient for two reasons:
 
1. **Vp can be extremely small** — REML reports `Vp = 0.000003` (3e-6). One
   x1000 pass would suffice here, but for even smaller values (REML reports
   `0.000000`), a single pass leaves Vp still below threshold.
 
2. **Full-GRM vs sparse-GRM discrepancy** — inbred mode sparsifies the GRM
   (`--make-bK-sparse 0.05`), which can yield a lower internal Vp than the
   full-GRM REML estimate. A phenotype that passes the REML check can still
   fail inside `--fastGWA-mlm-exact`.
 
The fix replaces the single-pass Python guard with an iterative REML loop in
`GCTA_MAKE_GRM` (primary) plus a catch-and-retry in `GCTA_PERFORM_GWA`
(safety net for the sparse-GRM edge case).
  
---
 
## Changes
 
### 1. `modules/gcta/make_grm/main.nf` — iterative REML loop
 
**Output declaration** (line 16) — replace `pheno_hsq_and_par` tuple:
 
```
# OLD:
tuple path("..._sims.phen"), path("check_vp.hsq"), path("..._sims.par"), emit: pheno_hsq_and_par
 
# NEW:
tuple path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"), path("${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par"), emit: pheno
```
 
- Extension `.phen` → `.pheno` (matches what `GCTA_PERFORM_GWA` already expects)
- `check_vp.hsq` removed from output (still generated internally)
- Emit name `pheno_hsq_and_par` → `pheno`
 
**Script block** (lines 38-45) — replace single REML call + fallback with iterative loop:
 
```bash
# Iterative REML-scale loop: run REML, check Vp, scale x1000 per round
# until full-GRM Vp >= 1e-4 (max 4 rounds).
VP_THRESHOLD="0.0001"
VP_MAX_ROUNDS=4
PHEN_INPUT="${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen"
GRM_PREFIX="TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_gcta_grm_${mode}"
 
cp -L "\${PHEN_INPUT}" working_pheno.txt   # copy from symlink
 
for vp_round in \$(seq 1 \${VP_MAX_ROUNDS}); do
    if ! gcta64 --grm \${GRM_PREFIX} \
                --pheno working_pheno.txt \
                --reml --out check_vp \
                --thread-num 1; then
        echo "REML failed (round \${vp_round}); falling back to Vp=1.0"
        printf 'Source\\tVariance\\tSE\\nVp\\t1.0\\tNA\\n' > check_vp.hsq
        break
    fi
 
    VP_VALUE=\$(awk -F'\\t' '\$1 == "Vp" {print \$2}' check_vp.hsq)
    VP_OK=\$(awk -v vp="\${VP_VALUE}" -v thresh="\${VP_THRESHOLD}" \
             'BEGIN { print (vp+0 >= thresh+0) ? "1" : "0" }')
 
    if [ "\${VP_OK}" -eq 1 ]; then
        echo "Vp=\${VP_VALUE} >= \${VP_THRESHOLD} at round \${vp_round}; done"
        break
    fi
 
    if [ "\${vp_round}" -eq "\${VP_MAX_ROUNDS}" ]; then
        echo "WARNING: Vp=\${VP_VALUE} still below \${VP_THRESHOLD} after \${VP_MAX_ROUNDS} rounds"
        break
    fi
 
    echo "Vp=\${VP_VALUE} < \${VP_THRESHOLD} at round \${vp_round}; scaling x1000"
    awk '{printf "%s %s %.10g\\n", \$1, \$2, \$3 * 1000}' working_pheno.txt > scaled_pheno.tmp
    mv scaled_pheno.tmp working_pheno.txt
done
 
cp working_pheno.txt "${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"
```
 
**Stub block** (line 58) — change touch from `.phen` to `.pheno`, remove `check_vp.hsq` touch:
 
```bash
# OLD (lines 58-60):
touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.phen
touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par
touch "check_vp.hsq"
 
# NEW:
touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno
touch ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.par
```
 
### 2. `main.nf` — remove PYTHON_CHECK_VP, rewire
 
**Remove include** (line 21):
```groovy
// DELETE this line:
include { PYTHON_CHECK_VP               } from './modules/python/check_vp.nf'
```
 
**Remove PYTHON_CHECK_VP call and rewire** (lines 517-545) — replace the entire block from the comments through the GWA channel definitions:
 
```groovy
    // Simulate GWA using output from GCTA_MAKE_GRM
    ch_type = Channel.of(
        "pca",
        "nopca"
        )
 
    ch_gwa_params = GCTA_MAKE_GRM.out.params.combine(ch_type)
    ch_gwa_grm =    GCTA_MAKE_GRM.out.grm.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    ch_gwa_plink =  GCTA_MAKE_GRM.out.plink.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
    ch_gwa_pheno =  GCTA_MAKE_GRM.out.pheno.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
```
 
This replaces all `PYTHON_CHECK_VP.out.*` references with `GCTA_MAKE_GRM.out.*`.
 
Also **remove** the versions line (529):
```groovy
// DELETE:
ch_versions = ch_versions.mix(PYTHON_CHECK_VP.out.versions)
```
 
### 3. Delete files
 
- `modules/python/check_vp.nf`
- `bin/check_vp.py`
 
### 5. `conf/rockfish.config` — remove label
 
Delete lines 143-150 (the `withLabel: python_check_vp` block).
 
`conf/docker.config` `withName: 'PYTHON_.*'` stays — still needed for `PYTHON_SIMULATE_EFFECTS_GLOBAL`.
 

---
 
## Fix rationale
 
- **DB storage unaffected**: Pre-upscaled phenotype is captured from `GCTA_SIMULATE_PHENOTYPES` output, upstream of `GCTA_MAKE_GRM`
- **Scale invariance**: P-values, ANOVA SS ratios, QTL detection are all scale-invariant
- **Channel structure preserved**: `GCTA_MAKE_GRM.out.{params, grm, plink, pheno}` emit tuples with identical structure to what `PYTHON_CHECK_VP.out.*` 
