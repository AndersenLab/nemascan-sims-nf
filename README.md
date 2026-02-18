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

### Define QTL Regions of Interest
* QTL regions are defined with the NF process `R_GET_GCTA_INTERVALS`
* This process runs an Rscript `bin/Get_GCTA_Intervals.R` which defines QTL intervals from the mapping outputs of `GCTA_PERFORM_GWA` and several pipeline parameters
* Essentially the script takes the raw mapping results, applies significance criteria, groups significant markers into QTLs, defines confidence intervals for those QTL regions and estimates their effect size.
* The script has several processing steps
    1. Load the libraries and command line arguments given by the pipeline
    2. Load the input data
        * Phenotype Data
        * GCTA Mapping Data
        * Genotype matrix
    3. Set the significance threshold
        * The script accepts a commandline argument (argument #10) which specifies the significance threshold to be applied to markers.
        * This argument can one of the following:
            * `BF` - Bonferroni threshold
            * `EIGEN` - Defined by the number of independent tests from Eigen decomposition of the genotype matrix
            * A user defined numeric value that is used as the threshold.
        * The `--sthresh` argument supplied to the pipeline sets the input for this argument to the process.
            * The default setting for the pipeline `--sthresh` argument is the `BF` threshold in the `nextflow.config` file.
    4. Process Mapping Data w/ `process_mapping_df()` function
#### `process_mapping_df()`
This is the core function of the script and performs many operations to define QTL intervals. The function returns the variable `Processed` which contains the original mapping data and these additional columns
- `strain`
- `value`
- `allele`
- `var.exp`
- `startPOS`: the starting position of the QTL interval
- `peakPOS`: The position of the peak marker of the QTL interval
- `endPOS`: the end position of the QTL interval
- `peak_id`: the id of the QTL interval. Is `1` if there is just one QTL identified for the trait or `2`..`Inf` if there are multiple QTL identified for the trait.
- `interval_size`: The number of bases spanned by the QTL interval
1. Threshold application
    * Step calculates the significance threshold and identifies marker SNPs exceeding that threshold.
        1. First the mapping df is grouped by trait (in the case that multiple mappings of different traits occurred)
        2. Depending on the threshold set, SNPs are flagged as being above `1` or below `0` the significance threshold in a newly created column `aboveBF`
    * Note: The function uses an externally defined variable `QTL_cutoff` which is not passed as an argument to the script.
    * The column to denote if a SNP is above the significance threshold is named `aboveBF` regardless of the significance threshold that is applied. This is likely required so that the outputs have standard formatting for later processing steps.
2. Filtering
    * After applying the significance threshold to flag SNPs as either above (`1`) or below (`0`) the significance threshold in the column `aboveBF` there are three possible next steps
    1. If more than 15% of the total SNPs are above the significance threshold all columns added by the mapping function `process_mapping_df()` are set to `NA`.
    2. If there are no significant SNPs the columns added by the `process_mapping_df()` are also set to `NA`
3. Variant effect calculation for significant SNPs
    * This step adds the `var.exp` column to the processed mapping result by correlating phenotype values with genotype values at the significant SNPs.
    * Uses pearsons correlation R2 between the phenotype values and allelic state (REF/ALT).
4. QTL interval definition
    * Identifies the most significant SNP in a QTL region of interest
The output is a processed dataframe containing the original mapping data augmented with QTL interval information (start, peak, end positions, peak ID, interval size, and Variance explained)

## Assessing Mapping Performance
The process `R_ASSESS_SIMS` runs the Rscript `Assess_Sims.R` to evaluate the performance of GWAS simulations.

It loads the simulated trait outputs, mapping outputs, and a number of simulation pipeline parameters.

This final process outputs a `simulation_assessment_results.tsv` to the analysis directory. Each row represents a QTL simulated or Detected with the following columns:

- `QTL`: Peak marker ID for the QTL interval
- `Simulated`: TRUE/FALSE if the QTL was simulated
- `Detected`: TRUE/FALSE if the QTL was detected in mapping
- `CHROM`: Chromosome of the QTL (numeric ID e.g., 1 = I, 2 = II, etc.)
- `POS`: Position of the marker
- `RefAllele`: Reference allele for the marker
- `Frequency`: Allele frequency of the marker
- `Effect`: Effect size of the marker
- `Simulated.QTL.Var.Exp`: Variance explained by the simulated QTL
- `log10p`: -log10(p-value) of the marker from mapping
- `aboveBF`: TRUE/FALSE if the peak marker is above the significance threshold (see `algorithm_id` column)
- `startPOS`: Start position of the QTL interval
- `peakPOS`: Peak position of the QTL interval
- `endPOS`: End position of the QTL interval
- `detected.peak`: TRUE/FALSE if the marker is the detected peak in mapping
- `interval.Frequency`: Allele frequency of the peak marker in the QTL interval
- `BETA`: Effect size estimate from mapping for the peak marker
- `interval.log10p`: -log10(p-value) of the peak marker in the QTL interval
- `peak_id`: Numeric ID of the QTL interval (e.g 1, 2, ... n, where N is the total number of QTL detected)
- `interval_size`: Size of the QTL interval in base pairs
- `interval.var.exp`: Variance explained by the peak marker in the QTL interval
- `top.hit`: TRUE/FALSE if the marker is the top hit in QTL interval
- `nQTL`: Number of QTL simulated for the trait
- `simREP`: Replicate number of the simulation
- `h2`: Heritability of the simulated trait
- `maf`: Minor allele frequency threshold used in simulation
- `effect_distribution`: Effect size range used in simulation
- `strain_set_id`: Name of the strain set used in simulation
- `algorithm_id`: Mapping method (Inbred, Loco, Inbred + PCA, LOCO + PCA) and significance threshold (e.g. `inbred_pca_EIGEN`, or `inbred_pca_BF`)

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