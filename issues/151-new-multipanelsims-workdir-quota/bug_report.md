

# Bug Report: Large Multipanel Simulations Exceed Workdir Quota

## Description

Pipeline execution failed with the above error message. This issue is related to the new multipanelsims workdir quota being filled up, which is causing the pipeline to fail when it tries to write output files.

Tried to include 750 mapping panels in a single Nextflow run, which exceeded the storage quota for the work directory. This is likely due to the large number of intermediate files generated during the pipeline execution.

With current workdir output this test run exceeds current capacities.

## Command Used

```bash
nextflow run main.nf \
>   -profile rockfish \
>   --strainfile data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/strains_panelsize.tsv \
>   --nqtl data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/nqtl.csv \
>   --h2 data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/h2.csv \
>   --effect data/sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps/effect_sizes.csv \
>   --reps 50 \
>   --cv_maf 0.05 \
>   --cv_ld 0.8 \
>   --output_dir Sims_panelsize_ce-cb-ct_nqtl5_h2-0.8_50reps \
>   -work-dir /scratch4/eande106/Ryan_tmp/nf-work-panelsize
```

Breakdown of the mapping panels included in the run:

```
        species n_strains  n
1    c_briggsae       100 50
2    c_briggsae       200 50
3    c_briggsae       300 50
4    c_briggsae       400 50
5    c_briggsae       500 50
6     c_elegans       100 50
7     c_elegans       200 50
8     c_elegans       300 50
9     c_elegans       400 50
10    c_elegans       500 50
11 c_tropicalis       100 50
12 c_tropicalis       200 50
13 c_tropicalis       300 50
14 c_tropicalis       400 50
15 c_tropicalis       500 50
```

## Error Output

```
Command exit status:
executor >  local (1), slurm (8173)
[ef/dabf1d] LOCAL_GET_CONTIG_INFO (ce.n100.r01)                   [100%] 1 of 1 ✔
[1b/b97f73] BCFTOOLS_EXTRACT_STRAINS (ct.n200.r29)                [100%] 750 of 750 ✔
[c2/db0bae] BCFTOOLS_RENAME_CHROMS (cb.n200.r04)                  [100%] 750 of 750 ✔
[da/725444] PLINK_RECODE_MS_VCF (cb.n500.r15_0.05)                [ 97%] 786 of 809, failed: 63, retries: 59
[82/7423fa] PLINK_RECODE_CV_VCF (cb.n500.r25_0.05)                [100%] 750 of 750 ✔
[68/02c61b] BCFTOOLS_CREATE_GENOTYPE_MATRIX (cb.n500.r39_0.05)    [100%] 700 of 700
[09/9a7c9d] R_FIND_GENOTYPE_MATRIX_EIGEN (6 ce.n100.r10_0.05)     [ 86%] 3659 of 4221, failed: 21, retries: 21
[-        ] LOCAL_COMPILE_EIGENS                                  -
[-        ] PYTHON_SIMULATE_EFFECTS_GLOBAL                        -
[-        ] GCTA_SIMULATE_PHENOTYPES                              -
[-        ] PLINK_UPDATE_BY_H2                                    -
[-        ] GCTA_MAKE_GRM                                         -
[-        ] GCTA_PERFORM_GWA                                      -
[-        ] DB_MIGRATION_WRITE_MARKER_SET                         -
[32/bdbf4f] DB_MIGRATION_WRITE_GENOTYPE_MATRIX (ct.n500.r38_0.05) [ 89%] 626 of 700
[-        ] DB_MIGRATION_WRITE_TRAIT_DATA                         -
[-        ] DB_MIGRATION_WRITE_GWA_TO_DB                          -
[-        ] DB_MIGRATION_AGGREGATE_METADATA                       -
[-        ] DB_MIGRATION_ANALYZE_QTL                              -
[-        ] DB_MIGRATION_ASSESS_SIMS                              -
ERROR ~ Error executing process > 'PLINK_RECODE_MS_VCF (cb.n400.r02_0.05)'

Caused by:
  Process `PLINK_RECODE_MS_VCF (cb.n400.r02_0.05)` terminated with an error exit status (1)


Command executed:

  cp cb.n400.r02_renamed.vcf.gz recoded.vcf.gz
  cp cb.n400.r02_renamed.vcf.gz.tbi recoded.vcf.gz.tbi
  
  plink --vcf recoded.vcf.gz \
      --snps-only \
      --biallelic-only \
      --maf 0.05 \
      --set-missing-var-ids @:# \
      --indep-pairwise 50 10 0.8 \
      --geno \
      --not-chr 7 \
      --allow-extra-chr
  
  plink --vcf recoded.vcf.gz \
      --make-bed \
      --snps-only \
      --biallelic-only \
      --maf 0.05 \
      --set-missing-var-ids @:# \
      --extract plink.prune.in \
      --geno \
      --recode \
      --out TO_SIMS \
      --allow-extra-chr
  
  awk -F":" '$1=$1' OFS="\t" plink.prune.in | \
      sort -k1,1d -k2,2n > markers.txt
  
  cat <<-END_VERSIONS > versions.yml
  "PLINK_RECODE_MS_VCF":
      plink: $( plink --version |& head -n 1 | cut -f 2 )
  END_VERSIONS

Command exit status:
  1

Command output:
  (empty)

Command error:
  cp: error writing 'recoded.vcf.gz': Quota exceeded

Work dir:
  /scratch4/eande106/Ryan_tmp/nf-work-panelsize/30/8943e77e3ff05ba11a0ecee7d56a5c

Container:
  /vast/eande106/singularity/andersenlab-plink-1.9.img

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

---

*Generated with the Quarto GitHub Issue Framework. Submit via `gh issue create --title "<title>" --label "bug" --body-file <this_file>.md`*
