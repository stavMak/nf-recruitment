#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input  = 'samplesheet.csv'
params.outdir = 'results'

process FASTQC {
    tag "${sample}"
    conda 'bioconda::fastqc=0.11.9'
    publishDir { "${params.outdir}/1_fastqc/${sample}" }, mode: 'copy'
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

process COUNT_READS {
    tag "${sample}"

    input:
    tuple val(sample), path(reads)

    output:
    path("${sample}_count.tsv"), emit: count

    script:
    """
    n_lines=\$(zcat ${reads[0]} | wc -l)
    n_reads=\$((n_lines / 4))
    printf "%s\\t%s\\n" "${sample}" "\${n_reads}" > ${sample}_count.tsv
    """
}

process FASTP {
    tag "${sample}"
    conda 'bioconda::fastp=0.23.4'
    publishDir { "${params.outdir}/2_fastp/${sample}" }, mode: 'copy'
    cpus 4

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_trimmed_1.fq.gz"), path("${sample}_trimmed_2.fq.gz"), emit: trimmed
    path("${sample}_fastp.json"), emit: json
    path("${sample}_fastp.html"), emit: html

    script:
    """
    fastp \\
        -i ${reads[0]} -I ${reads[1]} \\
        -o ${sample}_trimmed_1.fq.gz -O ${sample}_trimmed_2.fq.gz \\
        --json ${sample}_fastp.json \\
        --html ${sample}_fastp.html \\
        --thread ${task.cpus}
    """
}

process MULTIQC_FASTQC {
    conda 'bioconda::multiqc=1.19'
    publishDir "${params.outdir}/3_multiqc_fastqc", mode: 'copy'

    input:
    path(fastqc_files)

    output:
    path("multiqc_report.html")
    path("multiqc_data")

    script:
    """
    multiqc .
    """
}

process MULTIQC_FASTP {
    conda 'bioconda::multiqc=1.19'
    publishDir "${params.outdir}/4_multiqc_fastp", mode: 'copy'

    input:
    path(fastp_files)

    output:
    path("multiqc_report.html")
    path("multiqc_data")

    script:
    """
    multiqc .
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
    COUNT_READS(reads_ch)
    FASTP(reads_ch)

    // Merge every sample's small count file into one norm.tsv
    COUNT_READS.out.count
        .collectFile(name: 'norm.tsv', storeDir: params.outdir)

    // MultiQC on the raw FastQC reports
    MULTIQC_FASTQC(FASTQC.out.qc.map { sample, files -> files }.collect())

    // MultiQC on the fastp reports (json + html)
    MULTIQC_FASTP(FASTP.out.json.mix(FASTP.out.html).collect())
}