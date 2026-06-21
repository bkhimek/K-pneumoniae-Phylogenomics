process IQTREE {
    publishDir "${params.outdir}/iqtree", mode: 'copy'

    input:
    path filtered_aln

    output:
    path "iqtree.treefile", emit: treefile
    path "iqtree.iqtree",   emit: log
    path "iqtree.*",        emit: all

    script:
    """
    iqtree2 \\
        -s ${filtered_aln} \\
        -m GTR+G \\
        -bb 1000 \\
        -alrt 1000 \\
        -T AUTO \\
        --threads-max ${task.cpus} \\
        -pre iqtree \\
        -redo
    """
}
