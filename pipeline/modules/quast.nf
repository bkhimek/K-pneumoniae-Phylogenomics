process QUAST {
    tag "${meta.id}"
    container 'quay.io/biocontainers/quast:5.2.0--py39pl5321h2add14b_1'
    publishDir "${params.outdir}/quast", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_transposed_report.tsv"), emit: report

    script:
    """
    quast.py ${fasta} \
        --output-dir quast_out \
        --threads ${task.cpus} \
        --fast

    cp quast_out/transposed_report.tsv ${meta.id}_transposed_report.tsv
    """
}
