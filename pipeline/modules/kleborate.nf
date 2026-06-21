process KLEBORATE {
    tag "${meta.id}"
    publishDir "${params.outdir}/kleborate", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_kleborate.tsv"), emit: results
    path "${meta.id}_kleborate.tsv",                   emit: tsv

    script:
    """
    kleborate \\
        -a ${assembly} \\
        -o ${meta.id}_kleborate.tsv \\
        --preset kpsc
    """
}
