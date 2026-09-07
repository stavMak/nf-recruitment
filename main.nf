#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input  = 'samplesheet.csv'
params.outdir = 'results'

process FASTQC {
    tag "$sample"
    conda 'bioconda::fastqc=0.11.9'
    publishDir "${params.outdir}/1_fastqc/${sample}", mode: 'copy'
    cpus 4

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_fastqc.{zip,html}"), emit: qc

    script:
    """
    fastqc ${reads} --threads ${task.cpus}
    """
}

workflow {
    // Read the samplesheet, one row per sample, and build (sample, [R1, R2]) tuples
    reads_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            tuple(row.sample, [file(row.fastq_1), file(row.fastq_2)])
        }

    FASTQC(reads_ch)
}
