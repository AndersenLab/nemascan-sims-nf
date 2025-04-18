process LOCAL_GET_CONTIG_INFO {

    tag "${meta.id}"
    label "local_get_contig_info"

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    path "contigs.txt",          emit: contigs
    path "contig_mapping.tsv",   emit: mapping

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    zcat ${vcf} | head -n 200 | grep "##contig" | \\
        awk 'BEGIN{COUNT=0}{split(\$1,A,"="); split(A[3],B,","); CHROM=B[1]; printf "%s\\t%i\\n", CHROM, COUNT; COUNT=COUNT+1}' > contig_mapping.tsv
    cut -f 1 contig_mapping.tsv > contigs.txt
    """

    stub:
    """
    touch contigs.txt
    touch contig_lengths.tsv
    touch genome_partition.txt
    """
}
