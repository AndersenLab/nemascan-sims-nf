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
Traits are simulated through a series of steps.
1) Select causal variants 
2) Simulate phenotypes
3) Update by H2
4) Make genetic relatedness matrix (GRM) & Check phenotypic variance
5) Check & Update Vp
### Selecting causal variants
Causal variants are selected from the marker set in the process `PYTHON_SIMULATE_EFFECTS_GLOBAL` which runs the script `bin/create_causal_vars.py`. 

The `.bim` file containing the markers is input into the script, along with the number of causal variants (`nQTL`) and the effect range. The effect range is either a numeric range (e.g. `0.4-0.9`) or `gamma`. 

First the script randomly selects N variants without replacement from the marker set with the `select_variants()` function. 

Next, the script assigns the selected variants an effect size.

If the effect range is `gamma` the script pulls the effect sizes from a gamma distribution `gamma(effect_shape=0.4, effect_scale=1.66)` and randomly assigns the variants a direction (either `-1` or `1`). These steps are performed by the function `simulate_effect_gamma()` which is defined in the script.

If the effect range is numeric (e.g. `0.4-0.9`) a uniform distribution is created ranging from the `low_end` of the distribution (`0.4`) to the high_end of the distribution (`0.9`) and effects are randomly drawn from the distribution and a direction (either `-1` or `1`) is assigned. These steps are perfromed by the `simulate_effect_uniform()` function defined in the script.

Finally, the script saves the causal_variants to a file `causal_variants.txt` in the output directory.

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
### Simulating phenotypes
The outputs are then passed to the `GCTA_SIMULATE_PHENOTYPES` process which runs the command `gcta64 --simu-qt` to simulate a quantitative trait [GWAS Simulation documentation](https://yanglab.westlake.edu.cn/software/gcta/#GWASSimulation)

The parameter `--simu-causal-loci` supplies the causal variants selected in the prior step. GCTA expects the input to have two columns (SNP ID and effect size)

The parameter `--simu-hsq` specifies the heritability of the trait. This is passed to the command by the values in the file suppled to the pipeline parameter `--h2`

There are two outputs from this process:
`{prefix}.par`
  - The parfile has a header and the following columns:
    - QTL: SNPid of the causal variant
    - RefAllele: Reference allele 
    - Frequency
    - Effect size
`{prefix}.phen`
  - This is the simulated phenotype data. The file has no header and multiple columns:
    - family ID
    - individual ID
    - simulated phenotypes
