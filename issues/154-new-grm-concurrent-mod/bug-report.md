
# Bug Report: Concurrent Modification Issue when executing `GCTA_MAKE_GRM`

## Description

Resumed failed run from [#153](https://github.com/AndersenLab/nemascan-sims-nf/issues/153)


## Command used

```bash
nextflow run main.nf \
  -profile rockfish \
  --strainfile data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/strains_panelsize_ce.tsv \
  --nqtl data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/nqtl.csv \
  --h2 data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/h2.csv \
  --effect data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/effect_sizes.csv \
  --reps 50 \
  --cv_maf 0.05 \
  --cv_ld 0.8 \
  --output_dir Sims_panelsize_ce_nqtl5_h2-0.8_50reps \
  -work-dir /scratch4/eande106/Ryan_tmp/nf-work-panelsize \
  -resume
```

## Error output

```
ERROR ~ Error executing process > 'GCTA_MAKE_GRM (5 42 0.8 gamma inbred ce.n100.r05_0.05)
Caused by:
  java.util.ConcurrentModificationException
 -- Check '.nextflow.log' file for details                                  
```

## `.nextflow.log`

```bash
grep -nE "GCTA_MAKE_GRM|ConcurrentModificationException|ERROR ~" .nextflow.log
```

Selected output from the command:


```
10306:  task: name=GCTA_MAKE_GRM (5 42 0.8 gamma inbred ce.n100.r05_0.05); work-dir=/scratch4/eande106/Ryan_tmp/nf-work-panelsize/29/f1428d9f8ab9d002dc1c97e424541c
10307:  error [java.util.ConcurrentModificationException]: java.util.ConcurrentModificationException
10309:Apr-23 18:32:48.538 [Actor Thread 26] ERROR nextflow.processor.TaskProcessor - Error executing process > 'GCTA_MAKE_GRM (5 42 0.8 gamma inbred ce.n100.r05_0.05)'
10312:  java.util.ConcurrentModificationException
10315:java.util.ConcurrentModificationException: null
10340:Apr-23 18:32:48.635 [Actor Thread 26] DEBUG nextflow.Session - Session aborted -- Cause: java.util.ConcurrentModificationException
10373:[process] GCTA_MAKE_GRM
11704:  task: name=GCTA_MAKE_GRM; work-dir=null
11733:    Error report: Error executing process > 'GCTA_MAKE_GRM (5 42 0.8 gamma inbred ce.n100.r05_0.05)'
11736:  java.util.ConcurrentModificationException
```

## Notes

We have encountered some concurrent modification errors before Issue #141 and Issue #148

In #148 the error was in the previous method used to create the `ch_db_analysis_pheno` input channel for database ingestion.

In #141 The same error was triggered when executing the `PLINK_UPDATE_BY_H2` process. 

Each output of `GCTA_SIMULATE_PHENOTYPES` emited an output channel that was consumed by two operators without explicit forking. 
