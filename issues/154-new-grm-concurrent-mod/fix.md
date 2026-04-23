## Root Cause

Same class of bug as issue #148 — the wrap-combine-unwrap anti-pattern in channel construction:

```groovy
// Location 1 (main.nf:511-512 before fix):
ch_grm_plink = PLINK_UPDATE_BY_H2.out.plink.map{ it: [it] }.combine(ch_mode).map{ it: it[0] }
ch_grm_pheno = PLINK_UPDATE_BY_H2.out.pheno.map{ it: [it] }.combine(ch_mode).map{ it: it[0] }

// Location 2 (main.nf:529-531 before fix):
ch_gwa_grm   = GCTA_MAKE_GRM.out.grm.map{   it: [it] }.combine(ch_type).map{ it: it[0] }
ch_gwa_plink = GCTA_MAKE_GRM.out.plink.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
ch_gwa_pheno = GCTA_MAKE_GRM.out.pheno.map{ it: [it] }.combine(ch_type).map{ it: it[0] }
```

`.combine(ch_mode)` produces two emissions from each upstream tuple (one for "inbred", one for "loco"). `.map{ it: it[0] }` returns the **same underlying ArrayList reference** for both emissions. Under SLURM job array batching + Singularity, two submission threads (one per mode/type) iterated the shared ArrayList simultaneously in `BashWrapperBuilder.createContainerBuilder` → `ConcurrentModificationException` inside `TaskArrayCollector.collect`.

The exception named `GCTA_MAKE_GRM` as the failing process. The submission-level race is most likely in Location 2 (the output fan-out that feeds `GCTA_PERFORM_GWA`'s job array), but Location 1 is equally vulnerable.

## Fix Applied

Replaced both wrap-combine-unwrap blocks with merge-before-expansion + `.multiMap{}`, consistent with the #148 fix. Each multiMap emission owns its own independent tuple — no shared ArrayList references between the inbred/loco or pca/nopca submissions in the same SLURM job array batch.

**Location 1 (PLINK_UPDATE_BY_H2 → GCTA_MAKE_GRM):**
```groovy
PLINK_UPDATE_BY_H2.out.params
    .merge(PLINK_UPDATE_BY_H2.out.plink)
    .merge(PLINK_UPDATE_BY_H2.out.pheno)
    .combine(ch_mode)
    .multiMap { group, maf, nqtl, effect, rep, h2,
                bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                pheno, par, mode, suffix ->
        params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix)
        plink:  tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
        pheno:  tuple(pheno, par)
    }
    .set { ch_grm }
```

**Location 2 (GCTA_MAKE_GRM → GCTA_PERFORM_GWA):**
```groovy
GCTA_MAKE_GRM.out.params
    .merge(GCTA_MAKE_GRM.out.grm)
    .merge(GCTA_MAKE_GRM.out.plink)
    .merge(GCTA_MAKE_GRM.out.pheno)
    .combine(ch_type)
    .multiMap { group, maf, nqtl, effect, rep, h2, mode, suffix,
                grm_bin, grm_n, grm_id,
                bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests,
                pheno, par, type ->
        params: tuple(group, maf, nqtl, effect, rep, h2, mode, suffix, type)
        grm:    tuple(grm_bin, grm_n, grm_id)
        plink:  tuple(bed, bim, fam, plink_map, nosex, ped, plink_log, gm, n_indep_tests)
        pheno:  tuple(pheno, par)
    }
    .set { ch_gwa }
```

## Verification

`NXF_VER=24.10.5 nextflow run main.nf -profile test -stub-run` — passed (`Success: true`, 4.2s)
