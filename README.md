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
    -output-dir       String             Name of folder that will contain the results (Default: Simulations_{date})
