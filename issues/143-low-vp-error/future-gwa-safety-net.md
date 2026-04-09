# Future: Catch-and-Retry Safety Net in GCTA_PERFORM_GWA

## Status: Deferred

The iterative REML loop in `GCTA_MAKE_GRM` (implemented in this branch)
handles the primary Vp scaling problem. However, there is a secondary edge
case where the fix may be insufficient:

## The sparse-GRM edge case

`GCTA_MAKE_GRM` estimates Vp using the **full dense GRM**. In inbred mode,
`GCTA_PERFORM_GWA` sparsifies the GRM (`--make-bK-sparse 0.05`) before
running `--fastGWA-mlm-exact`. The sparse GRM can yield a lower internal Vp
than the full-GRM REML estimate, causing GCTA to hard-fail with:

```
the Vp is below 1e-5
```

A phenotype that passes the REML check in `GCTA_MAKE_GRM` could still fail
inside `GCTA_PERFORM_GWA` if sparsification reduces the effective Vp below
GCTA's internal threshold.

## Implementation plan

### File: `modules/gcta/perform_gwa/main.nf`

Replace the GWA command block (the single `gcta64` call after `plink_snplist.txt`
extraction, currently around line 63) with mode-branched logic:

```bash
awk '{print \$2}' TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}.bim > plink_snplist.txt

if [[ ${mode} == "inbred" ]]; then
    # Safety net: catch "Vp is below 1e-5" from sparse-GRM edge case
    VP_RETRY_MAX=3
    VP_RETRY=0
    PHENO_FILE="${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno"
    GWA_OUT="${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type}"

    cp -L "\${PHENO_FILE}" gwa_pheno_working.txt

    while true; do
        set +e
        gcta64 \${COMMAND} \\
            --bfile TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group} \\
            \${GRM_OPTION} \${GRM_PREFIX} \\
            \${COVAR} \\
            --out \${GWA_OUT} \\
            --pheno gwa_pheno_working.txt \\
            --extract plink_snplist.txt \\
            --thread-num \${GWA_THREADS}
        GWA_EXIT=\$?
        set -e

        if [ \${GWA_EXIT} -eq 0 ]; then
            break
        fi

        if grep -q "the Vp is below 1e-5" "\${GWA_OUT}.log" 2>/dev/null; then
            VP_RETRY=\$((VP_RETRY + 1))
            if [ \${VP_RETRY} -ge \${VP_RETRY_MAX} ]; then
                echo "ERROR: Vp still below 1e-5 after \${VP_RETRY_MAX} retries"
                exit 1
            fi
            echo "Vp below 1e-5 (attempt \${VP_RETRY}); scaling pheno x1000"
            awk '{printf "%s %s %.10g\\n", \$1, \$2, \$3 * 1000}' gwa_pheno_working.txt > gwa_scaled.tmp
            mv gwa_scaled.tmp gwa_pheno_working.txt
        else
            echo "GWA failed with non-Vp error (exit \${GWA_EXIT})"
            exit \${GWA_EXIT}
        fi
    done
else
    gcta64 \${COMMAND} \\
        --bfile TO_SIMS_${nqtl}_${rep}_${h2}_${maf}_${effect}_${group} \\
        \${GRM_OPTION} \${GRM_PREFIX} \\
        \${COVAR} \\
        --out ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_lmm-exact_${mode}_${type} \\
        --pheno ${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno \\
        --extract plink_snplist.txt \\
        --thread-num \${GWA_THREADS}
fi
```

The existing loco `mv` rename block (lines 72-75) stays unchanged after this.

### Key details

- **Inbred only:** Loco mode uses the full (non-sparse) GRM, so the
  full-GRM REML estimate in `GCTA_MAKE_GRM` matches what `--mlma-loco` sees.
  Only inbred mode (`--fastGWA-mlm-exact` with `--grm-sparse`) is affected.

- **`set +e` / `set -e`:** The retry loop must temporarily disable `errexit`
  to capture the GCTA exit code without killing the script.

- **Log file grep:** GCTA writes the "Vp is below 1e-5" message to
  `${GWA_OUT}.log`. The grep distinguishes Vp failures (retryable) from other
  GCTA errors (non-retryable, propagated immediately).

- **Max 3 retries:** Each retry scales phenotypes x1000 (variance x10^6).
  Three retries handle Vp down to ~1e-23, which is far below any realistic
  simulation output.

- **Symlink safety:** `cp -L` copies the input phenotype file before
  modification, avoiding mutation of upstream task outputs.

### Documentation updates needed

If implemented, update:

1. `CLAUDE.md` line 150 — change "future enhancement" to present tense
2. `docs/gwas-mapping.qmd` Section 2 (GCTA_PERFORM_GWA) — add a note about
   the inbred-mode safety net after the association mapping step description
3. `docs/gwas-mapping.qmd` inputs table — note that the pheno input may be
   re-scaled internally for inbred mode

### Verification

1. Stub-run: `nextflow run main.nf -profile test -stub-run`
2. Local test: `nextflow run main.nf -profile test` (exercises both modes)
3. HPC: reduced-reps production run; check `.nextflow.log` for retry messages
   and zero Vp failures in inbred mode
