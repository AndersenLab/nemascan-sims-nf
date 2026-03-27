During HPC testing with the three-species HPC test profile there were errors with the `DB_MIGRATION_WRITE_GWA_TO_DB` process and the `DB_MIGRATION_WRITE_TRAIT_DATA` process

Branch: `fix/merge-unit-with-integration`
Commit ID: `d01dbc9`

## Commands Run: 

### Generate Input Test Data

Generate test VCFs

```bash
SIMG=/vast/eande106/singularity/quay.io-biocontainers-bcftools-1.16--hfe4b78e_1.img

# C. elegans
singularity exec --bind /vast/eande106 --bind "$PWD" $SIMG \
  bash data/test/generate_test_ce_vcf.sh \
  /vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.isotype.vcf.gz

# C. briggsae
singularity exec --bind /vast/eande106 --bind "$PWD" $SIMG \
  bash data/test/generate_test_cb_vcf.sh \
  /vast/eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz

# C. tropicalis
singularity exec --bind /vast/eande106 --bind "$PWD" $SIMG \
  bash data/test/generate_test_ct_vcf.sh \
  /vast/eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz
```

Generate the HPC strain files

```bash
SIMG=/vast/eande106/singularity/quay.io-biocontainers-bcftools-1.16--hfe4b78e_1.img

singularity exec --bind /vast/eande106 --bind "$PWD" $SIMG \
  bash data/test/generate_hpc_strainfile.sh \
  /vast/eande106/data/c_elegans/WI/variation/20250625/vcf/WI.20250625.hard-filter.isotype.vcf.gz \
  /vast/eande106/data/c_briggsae/WI/variation/20250626/vcf/WI.20250626.hard-filter.isotype.vcf.gz \
  /vast/eande106/data/c_tropicalis/WI/variation/20250627/vcf/WI.20250627.hard-filter.isotype.vcf.gz
```

### Run the pipeline with the test profile w/o legacy assess profile

Set up the environment and singularity cache for the three-species profile on a login node

```bash
conda activate /data/eande106/software/conda_envs/nf24_env
export NXF_SINGULARITY_CACHEDIR=/vast/eande106/singularity
```

Run the pipeline

```bash
nextflow run main.nf \
  -profile test_hpc_three_species,rockfish \
  -work-dir /scratch4/eande106/Ryan/nf-work-test-hpc-three-species
```

## Nextflow Log Infomation

```
2026-03-26 14:43:48     46m 36s         crazy_lalande           ERR     8ed69d4cc9      53e0962f-5c18-4458-ac4f-acbd6a378ca3      nextflow run main.nf -profile test_hpc_three_species,rockfish -work-dir /scratch4/eande106/Ryan/nf-work-test-hpc-three-species          
```

Run-name: `crazy_lalande`

Check log for failures:

### `DB_MIGRATION_WRITE_GWA_TO_DB`
```bash
nextflow log crazy_lalande -f name,status,exit,workdir | grep -i "WRITE_GWA"
```

Returns:

```log
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca ce.hpc.popA_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/4e/3de34c95dc9830daa8598c37179443
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca ce.hpc.popA_0.05)      COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/b8/b73040673ae6c0c5f60015e718e346
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca ce.hpc_0.05)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/21/ff61768eb6faffb140234143f958ff
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca cb.hpc.popB_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/3d/38793ba7d69256ea4d69ffdf8513f8
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca ce.hpc_0.05)       COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/2b/eb20bdb8ee6635b898dd8c85f71c78
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca ce.hpc.popA_0.05)  COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/09/cc2ad2fd6413781bc90e7d577d6264
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca ce.hpc.popA_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/b4/6f10a28aa6be67c050638426b39b8d
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca cb.hpc.popB_0.05)      COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/ea/583115b4de66ae455b7f7ee50adeed
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca cb.hpc.popB_0.05)  COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/3a/08bbd364e4834ec59d0a14f08e1f0f
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca ct.hpc.popB_0.05)        FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/ff/7b1f70b7568507121be79e69881e7f
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca cb.hpc.popB_0.05)        COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/b5/7c9fb0b5da9f645ca62293f69f3066
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca ct.hpc.popB_0.05)      FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/9a/2eaf3f5ef920055dfb4f5efcdd6bdc
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca ct.hpc.popA_0.05)        ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/5d/16ed829fc8e0786aa0c80c7d05f41f
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca ct.hpc.popA_0.05)      ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/e2/d57e8772666f5a615264d46f5e34ff
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca ce.hpc.popB_0.05)        ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/3d/58738009ccec1795e14a322bd15abf
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca ce.hpc.popB_0.05)        ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/8e/664753d48e369a59aefc31348a8900
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca ce.hpc.popB_0.05)  ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/3e/251cc14b3db1c3d1e61833927fdd45
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca ce.hpc.popB_0.05)      ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/8d/af82145d3ca8d986cf605375c05f48
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca ct.hpc.popA_0.05)  ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/d5/c2165e40ef1240aea4e45c862543e9
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca ct.hpc.popA_0.05)        ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/d8/79c3d1ae0852fd730446cb352c04cc
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca cb.hpc.popA_0.05)        ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/28/e6c650350c7c6c0eff6c5d359d6f34
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca cb.hpc.popA_0.05)      ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/bf/c2d1e2bf5763a9546d19e35015c660
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca ct.hpc_0.05)     ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/39/c15e59d7e72f950a382c32807b109f
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca ct.hpc_0.05)       ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/21/9f4a2577b2278017b84edf06b477d1
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca cb.hpc.popA_0.05)        ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/e3/51ed7c6537409ef92f1e5865a8b407
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca cb.hpc.popA_0.05)  ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/49/bd9578038c693cf492d72599fc7c22
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca ce.hpc_0.05)   ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/1c/558e370f37dd4d9e6f17be323f4ea3
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco nopca cb.hpc_0.05)     ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/1b/d37c38dd4af92f2a7e14b9a86753c8
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma loco pca cb.hpc_0.05)       ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/40/9ac5eb3e26351203445168f5884c54
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca ce.hpc_0.05)     ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/39/5592cf5dd544fc3baca9729c12d956
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca cb.hpc_0.05)   ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/23/1ade5fd390d311d9da57bd41ceaa68
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca cb.hpc_0.05)     ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/5b/9c5616947ff020e159e3eb6304fafa
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred pca ct.hpc_0.05)     ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/a8/a7553608cde7e7f3ec31456da8f45b
DB_MIGRATION_WRITE_GWA_TO_DB (5 1 0.8 gamma inbred nopca ct.hpc_0.05)   ABORTED -       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/8f/2a80d8b2be70a5ca15bc36ac16f4f3
```

### `DB_MIGRATION_WRITE_TRAIT_DATA`

```bash
nextflow log crazy_lalande -f name,status,exit,workdir | grep -i "WRITE_TRAIT"
```

```log
DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/7e/03277f518f2fafcfcbf8e768df1a76
DB_MIGRATION_WRITE_TRAIT_DATA (ce.hpc_5_1_0.8)  COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/c9/6e7a0063c064d05a5f3e3f2a079f61
DB_MIGRATION_WRITE_TRAIT_DATA (ce.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/d3/ff879edba1d1e50f600da81937142f
DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc_5_1_0.8)  COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/07/28eeff4458708b5e0e6c06313b94b8
DB_MIGRATION_WRITE_TRAIT_DATA (cb.hpc.popA_5_1_0.8)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/19/1f09d37c5eb9a97de6420480a900da
DB_MIGRATION_WRITE_TRAIT_DATA (ce.hpc.popA_5_1_0.8)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/26/f14e743352d5470c88f684b2363608
DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popA_5_1_0.8)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/1b/3dcdbc6f03c24873a3814650af460a
DB_MIGRATION_WRITE_TRAIT_DATA (cb.hpc.popB_5_1_0.8)     COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/c0/c9624f2a8f5c809920411120a1515f
DB_MIGRATION_WRITE_TRAIT_DATA (cb.hpc_5_1_0.8)  COMPLETED       0       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/07/42130d4bd69a72ff053ce497138b97
DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/ac/2310270dc5a6a156325028cb762f0e
DB_MIGRATION_WRITE_TRAIT_DATA (ce.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/6d/fa66cc0d3b8624a0f796ceee1988c1
DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/2e/e9a57d49bd412dc7c2590f7dc2f160
DB_MIGRATION_WRITE_TRAIT_DATA (ce.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/cc/bfba655a2006894110a972bad29aeb
DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/e7/587444dbd1580dfe51c0679ca7d9a0
DB_MIGRATION_WRITE_TRAIT_DATA (ce.hpc.popB_5_1_0.8)     FAILED  1       /scratch4/eande106/Ryan/nf-work-test-hpc-three-species/f7/11b434ace5e205c896cceb23869f5f
```

### Failed process workdirs

## Error messages: 

### `DB_MIGRATION_WRITE_TRAIT_DATA`

```
ERROR ~ Error executing process > 'DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popB_5_1_0.8)'

Caused by:
  Process `DB_MIGRATION_WRITE_TRAIT_DATA (ct.hpc.popB_5_1_0.8)` terminated with an error exit status (1)


Command executed:

  export R_SOURCE_DIR="/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/R"
  write_trait_data.R --group ct.hpc.popB --maf 0.05 --nqtl 5         --effect gamma --rep 1 --h2 0.8         --pheno_file 5_1_0.8_0.05_gamma_ct.hpc.popB_sims.phen --par_file 5_1_0.8_0.05_gamma_ct.hpc.popB_sims.par --base_dir /vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Analysis_Results-20260326/db         --causal_geno_file causal_genotypes.sim.5.1.tsv         --cv_maf_effective 0.05         --cv_ld 0.8
  
  cat <<-END_VERSIONS > versions.yml
  "DB_MIGRATION_WRITE_TRAIT_DATA":
      R: $( Rscript --version |& cut -f 4 )
  END_VERSIONS

Command exit status:
  1

Command output:
  (empty)

Command error:
  WARNING: While bind mounting '/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf:/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf': destination is already in the mount point list
  
  Attaching package: ‘arrow’
  
  The following object is masked from ‘package:utils’:
  
      timestamp
  
  Loading required package: DBI
  
  Attaching package: ‘dplyr’
  
  The following objects are masked from ‘package:data.table’:
  
      between, first, last
  
  The following objects are masked from ‘package:stats’:
  
      filter, lag
  
  The following objects are masked from ‘package:base’:
  
      intersect, setdiff, setequal, union
  
  
  Attaching package: ‘purrr’
  
  The following object is masked from ‘package:data.table’:
  
      transpose
  
  Error: Marker set metadata not found for population='ct.hpc.popB', maf=0.05 in /vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Analysis_Results-20260326/db. Ensure DB_MIGRATION_WRITE_MARKER_SET completed before resuming. If VCF/species parameters changed, a full re-run (not -resume) is required.
  Execution halted
```

### `DB_MIGRATION_WRITE_GWA_TO_DB`

```
WARNING: While bind mounting '/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf:/vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf': destination is already in the mount point list

Attaching package: ‘arrow’

The following object is masked from ‘package:utils’:

    timestamp

Loading required package: DBI

Attaching package: ‘dplyr’

The following objects are masked from ‘package:data.table’:

    between, first, last

The following objects are masked from ‘package:stats’:

    filter, lag

The following objects are masked from ‘package:base’:

    intersect, setdiff, setequal, union

2026-03-26 15:30:16 [INFO] Reading GWA file: 5_1_0.8_0.05_gamma_ct.hpc.popB_lmm-exact_inbred_pca.fastGWA
2026-03-26 15:30:16 [INFO] Read 5008 markers from 5_1_0.8_0.05_gamma_ct.hpc.popB_lmm-exact_inbred_pca.fastGWA (fastGWA format)
Error: Marker set metadata not found for population='ct.hpc.popB', maf=0.05 in /vast/eande106/projects/Ryan/simulation_pipeline/nemascan-sims-nf/Analysis_Results-20260326/db. Ensure DB_MIGRATION_WRITE_MARKER_SET completed before resuming. If VCF/species parameters changed, a full re-run (not -resume) is required.
Execution halted
```