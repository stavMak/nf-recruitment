#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input            = 'samplesheet.csv'
params.outdir           = 'results'
params.reference        = '/media/penbio24/sata2/stavroula/20251027_spiking/20260909_FZB24_reference_genome/Gunter/FZB42.fasta'
params.min_identity     = 90
params.min_coverage     = 90
params.min_cov_metapop  = 10
params.run_mapping      = true   // set to true (--run_mapping true) to run mapping/filtering + MetaPop after checking QC

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

process BUILD_INDEX {
    conda 'bioconda::strobealign=0.13.0'
    cpus 8

    input:
    path(reference)

    output:
    tuple path(reference), path("${reference}.r*"), emit: indexed

    script:
    """
    strobealign -t ${task.cpus} -r 150 --create-index ${reference}
    """
}

process MAP_FILTER {
    tag "${sample}"
    conda 'bioconda::strobealign=0.13.0 bioconda::samtools=1.19 conda-forge::gawk'
    publishDir { "${params.outdir}/5_mapping/${sample}" }, mode: 'copy'
    cpus 8

    input:
    tuple val(sample), path(r1), path(r2)
    tuple path(reference), path(index_files)

    output:
    tuple val(sample), path("${sample}_filtered.bam"), path("${sample}_filtered.bam.bai"), emit: filtered
    path("${sample}_strobealign_log.txt"), emit: log

    script:
    """
    strobealign -t ${task.cpus} --use-index ${reference} ${r1} ${r2} 2> ${sample}_strobealign_log.txt | \\
        samtools sort -@ ${task.cpus} -o ${sample}_sorted.bam
    samtools index ${sample}_sorted.bam

    samtools view -h ${sample}_sorted.bam | \\
    gawk -v min_id=${params.min_identity} -v min_cov=${params.min_coverage} '
    BEGIN {OFS="\\t"}
    /^@/ {print; next}
    {
        seq      = \$10
        read_len = length(seq)
        cigar    = \$6
        aligned  = 0
        while (match(cigar, /([0-9]+)([M=X])/, a)) {
            aligned += a[1]
            cigar    = substr(cigar, RSTART+RLENGTH)
        }
        nm = -1
        for (j=12; j<=NF; j++) if (\$j ~ /^NM:i:/) { nm = substr(\$j,6); break }
        if (read_len > 0 && nm >= 0) {
            id  = (read_len - nm) / read_len * 100
            cov = aligned / read_len * 100
            if (id >= min_id && cov >= min_cov) print
        }
    }' | samtools view -b -o ${sample}_filtered.bam

    samtools index ${sample}_filtered.bam

    rm ${sample}_sorted.bam ${sample}_sorted.bam.bai
    """
}

process METAPOP {
    conda '/home/penbio24/miniconda3/envs/metapop'
    publishDir "${params.outdir}/6_metapop", mode: 'copy'
    cpus 32

    input:
    path(bam_files)
    path(norm_tsv)
    path(reference)

    output:
    path("MetaPop")

    script:
    """
    mkdir bam_dir
    for f in *.bam *.bam.bai; do
        ln -s "\$(realpath "\$f")" bam_dir/
    done

    mkdir ref_dir
    ln -s "\$(realpath ${reference})" ref_dir/

    metapop \\
        --input_samples bam_dir/ \\
        --threads ${task.cpus} \\
        --reference ref_dir/ \\
        --norm ${norm_tsv} \\
        --min_cov ${params.min_cov_metapop}
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

    // Mapping/filtering + MetaPop only run once you've checked QC and set --run_mapping true
    if (params.run_mapping) {
        BUILD_INDEX(file(params.reference))
        MAP_FILTER(FASTP.out.trimmed, BUILD_INDEX.out.indexed)

        // Collect all filtered BAMs (+ their index files) from every sample
        all_bams_ch = MAP_FILTER.out.filtered
            .flatMap { sample, bam, bai -> [bam, bai] }
            .collect()

        METAPOP(
            all_bams_ch,
            file("${params.outdir}/norm.tsv"),
            file(params.reference)
        )
    }
}