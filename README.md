## Usage

nextflow andersenlab/nemascan-sim-nf --strainfile /path/to/strainfile --vcf /path/to/vcf -output-dir my-results

Mandatory argument (General):

    --strainfile      File               A TSV file with two columns: the first is a name for the strain set and the second is a comma-separated strain list without spaces
    --vcf             File               Generally a CaeNDR release date (i.e. 20231213). Can also provide a user-specified VCF with index in same folder

Optional arguments (General):

    --nqtl            File               A CSV file with the number of QTL to simulate per phenotype, one value per line (Default is located: data/simulate_nqtl.csv)
    --h2              File               A CSV file with phenotype heritability, one value per line (Default is located: data/simulate_h2.csv)
    --rep             Integer            The number of replicates to simulate per number of QTL and heritability (Default: 2)
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
[ ] - Add
## Simulated traits
Trait simulation involves the following sequential steps:
1.  Select causal variants.
2.  Simulate phenotypes based on these variants.
3.  Update PLINK fileset with simulated phenotypes.
4.  Create a Genetic Relatedness Matrix (GRM) and estimate phenotypic variance.
5.  Verify and adjust phenotypic variance (Vp).
### 1. Selecting Causal Variants
Causal variants are chosen from the available marker set by the `PYTHON_SIMULATE_EFFECTS_GLOBAL` process, which executes the `bin/create_causal_vars.py` script.


This script takes the following inputs:
*   A `.bim` file (PLINK binary marker information).
*   The desired number of causal variants (`nQTL`).
*   The effect range, specified either as a numeric range (e.g., `0.4-0.9`) or as `gamma`.

The selection process involves two main steps:
1.  **Variant Selection**: The `select_variants()` function randomly chooses `nQTL` variants without replacement from the markers listed in the `.bim` file.
2.  **Effect Size Assignment**:
    *   If the effect range is `gamma`, the `simulate_effect_gamma()` function assigns effect sizes. These sizes are drawn from a gamma distribution (`gamma(effect_shape=0.4, effect_scale=1.66)`), and each variant is randomly assigned a direction (positive or negative effect, i.e., `1` or `-1`).
    *   If a numeric range (e.g., `0.4-0.9`) is provided, the `simulate_effect_uniform()` function assigns effect sizes. These are drawn from a uniform distribution spanning the specified `low_end` to `high_end`, and a direction (`1` or `-1`) is also randomly assigned.

The script outputs a file named `causal_variants.txt` in the designated output directory. This file lists the selected causal variant IDs and their assigned effect sizes.

```causal_variants.txt
14266 -0.46737319855194537
741541 -0.41074062864697813
4210220 -0.41437796351367484
14314458 -0.47016528365409455
627759 0.47718905922917665
14211523 -0.41997813721985267
2278072 -0.4364813139721949
3481620 0.43486733972778335
8059998 -0.45881269191927904
16687130 -0.43335413403568435
```
### 2. Simulating Phenotypes with `GCTA_SIMULATE_PHENOTYPES`
The `causal_variants.txt` file (generated in the previous step) is used by the `GCTA_SIMULATE_PHENOTYPES` process. This process employs the `gcta64 --simu-qt` command to simulate quantitative traits based on the selected causal variants. (Refer to the [GCTA GWAS Simulation documentation](https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation) for more details).

Key parameters for `gcta64 --simu-qt`:
*   `--simu-causal-loci`: This parameter takes the `causal_variants.txt` file, which provides the SNP IDs and their effect sizes for the simulation.
*   `--simu-hsq`: This specifies the target heritability (h²) of the trait. The value for this parameter is taken from the file provided to the `--h2` pipeline parameter.


This process generates two primary output files:
*   `{prefix}.par`: A parameter file with a header, detailing:
    *   `QTL`: SNP ID of the causal variant.
    *   `RefAllele`: Reference allele.
    *   `Frequency`: Allele frequency.
    *   `Effect size`: The effect size used in the simulation for that QTL.
*   `{prefix}.phen`: A phenotype file without a header, containing:
    *   Column 1: Family ID.
    *   Column 2: Individual ID.
    *   Column 3: Simulated phenotype value.

### 3. Updating PLINK Fileset with Simulated Phenotypes (`PLINK_UPDATE_BY_H2`)
This step, handled by the `PLINK_UPDATE_BY_H2` process, integrates the simulated phenotypes (from `{prefix}.phen`) into the primary PLINK fileset used for simulation (referred to as the "TO_SIMS" fileset). This associates the newly generated phenotype data with the existing genetic data for each individual, preparing it for downstream analyses such as association mapping.

The filtering parameters that are applied should remain consistent across all plink commands

[ ] - Evaluate if this process is needed. Would be useful if we wanted to select causal variants from a different pool of markers.


### 4. Creating GRM and Estimating Phenotypic Variance (`GCTA_MAKE_GRM`)
The `GCTA_MAKE_GRM` process performs two sequential GCTA operations:

1.  **Genetic Relatedness Matrix (GRM) Construction**:
    *   Command: `gcta64 --make-grm-inbred` or `gcta64 --make-grm`.
    *   Purpose: This command builds a GRM using the genetic data from the PLINK fileset. The GRM quantifies the genetic similarity between pairs of individuals.

2.  **Restricted Maximum Likelihood (REML) Analysis**:
    *   Command: `gcta64 --reml`.
    *   Purpose: This analysis uses the GRM (from step 1) and the simulated phenotype data (`{prefix}.phen`) to estimate the proportion of phenotypic variance (Vp) explained by the SNPs included in the GRM. (See [GCTA GREML analysis documentation](https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis)).

The output of the `gcta64 --reml` analysis is a plain text file with the `*.hsq` extension 

[ ] - example *.hsq output

The `PYTHON_CHECK_VP` process, which runs the `bin/check_vp.py` script, inspects the estimated phenotypic variance (Vp) from the `*.hsq` file (produced by GCTA REML). This step ensures that the simulated phenotypes exhibit sufficient variance.

Inputs to `bin/check_vp.py`:
*   The `*.hsq` file (containing Vp estimates).
*   The current phenotype file (e.g., `{prefix}.phen` from the GCTA simulation step).

Script Logic:
1.  The script parses the `*.hsq` file to extract the `Vp` value.
2.  It then checks the `Vp` against a threshold:
    *   If `Vp` is less than `0.000001`: The script scales up the phenotype values for all individuals by multiplying them by `1000`. This adjustment aims to increase the phenotypic variance.
    *   If `Vp` is greater than or equal to `0.000001`: The original phenotype values are retained without modification.
3.  The script writes the (potentially modified) phenotype data to a temporary file named `new_phenos.temp`.
4.  Finally, this `new_phenos.temp` file is renamed to reflect the specific simulation parameters, following a pattern like `${nqtl}_${rep}_${h2}_${maf}_${effect}_${group}_sims.pheno`. This becomes the final phenotype file for this simulation iteration.