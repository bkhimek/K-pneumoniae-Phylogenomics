#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { QUAST                  } from './modules/quast'
include { AMRFINDER              } from './modules/amrfinder'
include { KLEBORATE              } from './modules/kleborate'
include { SNIPPY; SNIPPY_CORE   } from './modules/snippy'
include { GUBBINS                } from './modules/gubbins'
include { IQTREE                 } from './modules/iqtree'
include { SNP_DISTS } from './modules/snp_dists'

workflow {

    ch_genomes = Channel
        .fromPath(params.input, type: 'file')
        .map { fasta ->
            def meta = [id: fasta.simpleName
                            .replaceAll(/_genomic$/, '')
                            .replaceAll(/_ASM.*$/, '')]
            tuple(meta, fasta)
        }

    reference = file(params.reference)

    // QC
    QUAST(ch_genomes)

    // AMR & typing
    AMRFINDER(ch_genomes)
    KLEBORATE(ch_genomes)

    // Per-sample WGS alignment
    SNIPPY(ch_genomes, reference)

    // Core SNP alignment
    ch_snippy_dirs = SNIPPY.out.outdir.collect()
    SNIPPY_CORE(ch_snippy_dirs, reference)

    // Recombination removal → ML tree
    GUBBINS(SNIPPY_CORE.out.full_aln)
    IQTREE(GUBBINS.out.filtered_aln)
    SNP_DISTS(SNIPPY_CORE.out.core_aln)
}
