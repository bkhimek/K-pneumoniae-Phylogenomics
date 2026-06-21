process GUBBINS {
    publishDir "${params.outdir}/gubbins", mode: 'copy'

    input:
    path full_aln

    output:
    path "gubbins.filtered_polymorphic_sites.fasta", emit: filtered_aln
    path "gubbins.final_tree.tre",                   emit: tree
    path "gubbins.*",                                emit: all

    script:
    """
    run_gubbins.py \\
        --threads 4 \\
        --prefix gubbins \\
        --tree-builder fasttree \\
        --mar \\
        ${full_aln}
    """
}
