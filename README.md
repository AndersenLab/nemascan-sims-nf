## Usage

nextflow andersenlab/nemascan-sim-nf --strainfile /path/to/strainfile --vcf /path/to/vcf -output-dir my-results

Mandatory argument (General):

    --strainfile      File               A TSV file with two columns: the first is a name for the strain set and the second is a comma-separated strain list without spaces
    --vcf             File               Generally a CaeNDR release date (i.e. 20231213). Can also provide a user-specified VCF with index in same folder

Optional arguments (General):

    --nqtl            File               A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: data/simulate_nqtl.csv)
    --h2              File               A CSV file with phenotype heritability, one value per line (Default is located: data/simulate_h2.csv)
    --reps             Integer            The number of replicates to simulate per number of QTL and heritability (Default: 2)
    --maf             File               A CSV file where each line is a minor allele frequency threshold to test for simulations (Default: data/simulate_maf.csv)
    --effect          File               A CSV file where each line is an effect size range (e.g. 0.2-0.3) to test for simulations (Default: data/simulate_effect_sizes.csv)
    --qtlloc          File               A BED file with three columns: chromosome name (numeric 1-6), start postion, end postion. The genomic range specified is where markers will be pulled from to simulate QTL (Default: null [which defaults to using the whole genome to randomly simulate a QTL])
    --sthresh         String             Significance threshold for QTL - Options: BF - for bonferroni correction, EIGEN - for SNV eigen value correction, or another number e.g. 4
    --group_qtl       Integer            If two QTL are less than this distance from each other, combine the QTL into one, (DEFAULT = 1000)
    --ci_size         Integer            Number of SNVs to the left and right of the peak marker used to define the QTL confidence interval, (DEFAULT = 150)
    --sparse_cut      Decimal            Any off-diagonal value in the genetic relatedness matrix greater than this is set to 0 (Default: 0.05)
    --simulate_qtlloc Boolean            Whether to simulate QTLs in specific genomic regions (Default: false)
    -output-dir       String             Name of folder that will contain the results (Default: Simulations_{date})

# Simulations

## Preparing marker sets

The marker set generation phase transforms a multi-sample VCF into LD-pruned PLINK binary filesets and a numeric genotype matrix through four processes: `BCFTOOLS_EXTRACT_STRAINS`, `BCFTOOLS_RENAME_CHROMS`, `PLINK_RECODE_VCF`, and `BCFTOOLS_CREATE_GENOTYPE_MATRIX`. See the [Marker Set Generation](docs/marker-set-generation.qmd) documentation for detailed process descriptions, commands, and parameter references.

## Simulated traits

The trait simulation phase selects causal variants, simulates quantitative phenotypes with GCTA, and prepares phenotype-annotated PLINK filesets through three processes: `PYTHON_SIMULATE_EFFECTS_GLOBAL`, `GCTA_SIMULATE_PHENOTYPES`, and `PLINK_UPDATE_BY_H2`. See the [Trait Simulation](docs/trait-simulation.qmd) documentation for detailed process descriptions, commands, and parameter references.

## GWAS Mappings

The GWAS mapping phase constructs genetic relatedness matrices, verifies phenotypic variance, and performs association mapping under four mode/type conditions (`inbred`/`loco` x `pca`/`nopca`) through three processes: `GCTA_MAKE_GRM`, `PYTHON_CHECK_VP`, and `GCTA_PERFORM_GWA`. See the [GWAS Mapping](docs/gwas-mapping.qmd) documentation for detailed process descriptions, commands, and parameter references.

### QTL Detection & Assessment

The final phase detects QTL intervals from GWA results and assesses detection
performance against simulated truth. Two parallel paths produce equivalent
output: the DB path (default) writes results to a Parquet database and queries
via DuckDB, while the legacy path (opt-in via `--legacy_assess`) operates
directly on intermediate files. See the [QTL Detection & Assessment](docs/qtl-detection-assessment.qmd)
documentation for detailed process descriptions, commands, and output schema.

## Test Data

### Generating the Test VCF

The test VCF (`data/test/test.vcf.gz`) is not committed to the repository due to its size (~71 MB). The script `data/test/generate_test_vcf.sh` creates it by subsetting the full CaeNDR isotype reference VCF to the 13 strains in `data/test/test_strains.txt` across chromosomes I, II, and V.

**Requirements:** bcftools (>= 1.16), tabix, and the source VCF

```bash
# Download the source VCF from CaeNDR (~7.8 GB)
# https://elegansvariation.org/data/release/latest
# File: WI.20220216.hard-filter.isotype.vcf.gz

# Generate the test VCF (takes 5-10 minutes)
./data/test/generate_test_vcf.sh /path/to/WI.20220216.hard-filter.isotype.vcf.gz
```

This produces:
- `data/test/test.vcf.gz` — BGZF-compressed VCF (13 samples, chromosomes I/II/V, monomorphic sites removed)
- `data/test/test.vcf.gz.tbi` — tabix index

The script also strips `##contig` headers for absent chromosomes (III, IV, X). This is required because `LOCAL_GET_CONTIG_INFO` parses contig headers to build the chromosome mapping, and headers for chromosomes without data cause downstream failures in `R_FIND_GENOTYPE_MATRIX_EIGEN`.

### Generating Integration Test Data

The script `tests/collect_test_data.sh` runs the pipeline with the `test` profile and collects outputs needed by integration tests. It runs the actual pipeline (not a separate reimplementation), so the test data is always consistent with the current code.

**Requirements:** Docker, Nextflow (NXF_VER=24.10.4 or any 24.10.x)

```bash
# Run pipeline and collect outputs (~20-45 min on first run)
bash tests/collect_test_data.sh

# Collect outputs from a previous pipeline run without rerunning
bash tests/collect_test_data.sh --collect-only

# Remove previous results before running
bash tests/collect_test_data.sh --clean
```

The script runs `nextflow run main.nf -profile test,docker --legacy_assess`, which produces both DB-path and legacy assessment outputs. It then copies the results into `tests/integration_data/`.

After collection, the script prints the exact command to run integration tests:

```bash
TEST_DB_DIR=tests/integration_data/db \
TEST_WORK_DIR=tests/.nf-work \
TEST_EXISTING_ASSESSMENT=tests/integration_data/simulation_assessment_results.tsv \
TEST_DB_ASSESSMENT=tests/integration_data/db_simulation_assessment_results.tsv \
Rscript tests/run_tests.R
```

Subsequent runs benefit from Nextflow's `-resume` behavior — only processes affected by code changes are rerun.

### Unit Tests

Unit tests do not require pipeline output. They use static fixtures in `tests/fixtures/` and run with:

```bash
Rscript tests/run_tests.R
```

Integration test files skip automatically when the required environment variables are not set.

## Development Notes

### Rendering Documentation

The `docs/` directory contains a [Quarto](https://quarto.org/) website project with static documentation describing the cross-validation framework, concordance scores, and how to run the analysis. No pipeline output is required to render.

```bash
quarto render docs/
```

The rendered site is output to `docs/_site/`. To deploy as a GitHub Page, configure the repository to serve from that directory.

## HPC Testing Guide (Rockfihs)

### 1. Environment Setup

```bash
git clone https://github.com/AndersenLab/nemascan-sims-nf .
cd nemascan-sims-nf
```

### 2. Generate Test Data

Request an interactive session on RF - equivlent to SRUN command on other SLURM managed HPCs

```bash
interact -n 1 -c 1 -a eande106 -m 64G -p queue-name -t “30”
```
Load `bcftools` and `tabix` versions that are pre-installed for simple strain and marker filtering operations to get test data. The defaults on rockfish are listed in the sample command below

```bash
$ module load bcftools
$ bcftools --version
>bcftools 1.15.1
>Using htslib 1.15.1-15-ge51f72f
>Copyright (C) 2022 Genome Research Ltd.
```

```bash
$ module load tabix
$ tabix --version
> tabix (htslib) 1.13+ds
> Copyright (C) 2021 Genome Research Ltd.
```

Set the path to the source VCF. In this example we will point to the C. elegans WI-hard-filter.isotype.vcf from the 20220216 CaeNDR release.
```bash
SOURCE_VCF=/vast/eande106/data/c_elegans/WI/variation/20220216/vcf/WI.20220216.hard-filter.isotype.vcf.gz
```

Run the script to generate the test data from the source VCF

[ ] - Currently an error when sourcing bcftools with `module load` command on RF

```bash
./data/test/generate_test_vcf.sh $SOURCE_VCF
```

Verify output:

```bash
ls -lh data/test/test.vcf.gz       # expect ~825MB
ls -lh data/test/test.vcf.gz.tbi   # expect ~37KB
```

### 3. Stub-Run Validation (Quick Wiring Check)

Stub-runs verify process wiring without executing real computations. These run in seconds on the login node.

Prior to running any form of the pipeline prepare the env

```bash
# soruce NF settings from bash profile
source ~/.bash_profile

# load the nextflow conda 
conda activate /data/eande106/software/conda_envs/nf24_env
```

#### 3.1 Fixed architecture (test profile)
Simulations with one mapping panel and one trait architecture 

```bash
nextflow run main.nf -profile test -stub-run
```

#### 3.2 Variable architecture (test_variable profile)
Simulations with one mapping panel and multipe trait architectures 

```bash
nextflow run main.nf -profile test_variable -stub-run
```

#### 3.3 With `--legacy_assess` flag 
Test parallel analysis to enable legacy post-processing modules for comparisons

```bash
nextflow run main.nf -profile test --legacy_assess -stub-run
```

### 4. End-to-End Test Run

Uses SLURM + Singularity on real test data.

```bash
nextflow run main.nf -profile test,rockfish
```

#### 4.1 End-to-End Test Run with `--legacy_assess` flag

```bash
nextflow run main.nf -profile test,rockfish --legacy_assess
```