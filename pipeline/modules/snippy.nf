process SNIPPY {
    tag "${meta.id}"
    publishDir "${params.outdir}/snippy", mode: 'copy'

    input:
    tuple val(meta), path(assembly)
    path reference

    output:
    tuple val(meta), path("${meta.id}"), emit: results
    path "${meta.id}",                   emit: outdir

    script:
    """
    snippy \\
        --cpus ${task.cpus} \\
        --outdir ${meta.id} \\
        --ref ${reference} \\
        --ctgs ${assembly} \\
        --force
    """
}

process SNIPPY_CORE {
    publishDir "${params.outdir}/snippy", mode: 'copy'

    input:
    path snippy_dirs
    path reference

    output:
    path "core.full.aln", emit: full_aln
    path "core.aln",      emit: core_aln
    path "core.*",        emit: all

    script:
    """
    snippy-core \\
        --ref ${reference} \\
        --prefix core \\
        ${snippy_dirs}
    """
}
