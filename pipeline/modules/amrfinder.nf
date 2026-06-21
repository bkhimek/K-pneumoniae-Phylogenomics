process AMRFINDER {
    tag "${meta.id}"
    publishDir "${params.outdir}/amrfinder", mode: 'copy'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("${meta.id}_amrfinder.tsv"), emit: results
    path "${meta.id}_amrfinder.tsv",                   emit: tsv

    script:
    """
    amrfinder \\
        --nucleotide ${assembly} \\
        --organism Klebsiella_pneumoniae \\
        --plus \\
        --name ${meta.id} \\
        --output ${meta.id}_amrfinder.tsv \\
        --threads ${task.cpus}
    """
}
