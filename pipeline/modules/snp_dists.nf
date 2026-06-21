process SNP_DISTS {
    publishDir "${params.outdir}/snippy", mode: 'copy'

    input:
    path core_aln

    output:
    path "pairwise_snp_matrix.tsv", emit: matrix

    script:
    """
    snp-dists ${core_aln} > pairwise_snp_matrix.tsv
    """
}
