#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input            = 'samplesheet.csv'
params.outdir           = 'results'
params.reference        = null
params.metapop_env      = null
params.min_identity     = 90
params.min_coverage     = 90
params.min_cov_metapop  = 20
params.id_min_metapop   = 90

// ---- Per-step switches: turn any stage on/off independently ----
// Downstream stages read their inputs from disk (results/) when the
// upstream stage is switched off, so you can run each step on its own
// AND delete work/ between runs.
params.run_fastqc       = true   // per-lane FastQC on raw reads
params.run_fastp        = true   // merge lanes + fastp trimming + read counts (norm.tsv)
params.run_multiqc      = true   // aggregate FastQC + fastp reports (fresh + from disk)
params.run_mapping      = true   // strobealign + filtering -> filtered BAMs
params.run_metapop      = true   // MetaPop on all BAMs (fresh + from disk)

process MERGE_LANES {
    tag "${sample}"
    publishDir { "${params.outdir}/0_merged/${sample}" }, mode: 'copy'

    input:
    tuple val(sample), path(reads1), path(reads2)

    output:
    tuple val(sample), path("${sample}_R1.fq.gz"), path("${sample}_R2.fq.gz"), emit: merged

    script:
    """
    cat \$(printf '%s\\n' ${reads1} | sort) > ${sample}_R1.fq.gz
    cat \$(printf '%s\\n' ${reads2} | sort) > ${sample}_R2.fq.gz
    """
}

process FASTQC {
    tag "${lane_id}"
    conda 'bioconda::fastqc=0.11.9'
    publishDir { "${params.outdir}/1_fastqc/${sample}" }, mode: 'copy'
    cpus 4

    input:
    tuple val(sample), val(lane_id), path(r1), path(r2)

    output:
    tuple val(sample), path("*_fastqc.{zip,html}"), emit: qc

    script:
    """
    fastqc ${r1} ${r2} --threads ${task.cpus}
    """
}

process COUNT_READS {
    tag "${sample}"
    publishDir "${params.outdir}/0_counts", mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    path("${sample}_count.tsv"), emit: count

    script:
    """
    n_lines=\$(zcat ${reads[0]} | wc -l)
    n_reads=\$((n_lines / 4))
    # Name must match the BAM prefix MetaPop keys on: <sample>_filtered.bam -> <sample>_filtered
    printf "%s\\t%s\\n" "${sample}_filtered" "\${n_reads}" > ${sample}_count.tsv
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
    conda params.metapop_env
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
        --output . \\
        --threads ${task.cpus} \\
        --reference ref_dir/ \\
        --norm ${norm_tsv} \\
        --min_cov ${params.min_cov_metapop} \\
        --id_min ${params.id_min_metapop}
    """
}

workflow {

    // ---------------- Inputs ----------------
    // Samplesheet has ONE ROW PER LANE. A sample sequenced over several
    // lanes appears on several rows sharing the same sample name.
    rows_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true)

    // Per-lane tuples: used for raw FastQC and for lane merging.
    // lane_id = the R1 filename without the _1.fq.gz suffix (unique per lane).
    per_lane_ch = rows_ch.map { row ->
        def lane_id = file(row.fastq_1).name.replaceAll(/_1\.fq\.gz$/, '')
        tuple(row.sample, lane_id, file(row.fastq_1), file(row.fastq_2))
    }

    // ---------------- Per-lane FastQC (raw reads) ----------------
    if (params.run_fastqc) {
        FASTQC(per_lane_ch)
    }

    // ---------------- Merge lanes -> trim (fastp) + read counts ----------------
    if (params.run_fastp) {
        // Collect every lane of each sample into one group
        grouped_ch = per_lane_ch
            .map { sample, lane_id, r1, r2 -> tuple(sample, r1, r2) }
            .groupTuple()   // -> (sample, [r1_laneA, r1_laneB, ...], [r2_laneA, r2_laneB, ...])

        MERGE_LANES(grouped_ch)

        // Standard (sample, [R1, R2]) shape for the downstream QC processes
        reads_ch = MERGE_LANES.out.merged.map { sample, r1, r2 -> tuple(sample, [r1, r2]) }

        COUNT_READS(reads_ch)   // per-sample counts are published to 0_counts/

        FASTP(reads_ch)
    }

    // ---------------- MultiQC (this run's reports + whatever is on disk) ----------------
    if (params.run_multiqc) {
        // FastQC reports
        fresh_fastqc = params.run_fastqc ? FASTQC.out.qc.map { s, files -> files } : Channel.empty()
        prev_fastqc  = Channel.fromPath("${params.outdir}/1_fastqc/**/*_fastqc.{zip,html}")
        MULTIQC_FASTQC(
            fresh_fastqc.mix(prev_fastqc).flatten().unique { it.name }.collect()
        )

        // fastp reports
        fresh_fastp = params.run_fastp ? FASTP.out.json.mix(FASTP.out.html) : Channel.empty()
        prev_fastp  = Channel.fromPath("${params.outdir}/2_fastp/**/*_fastp.{json,html}")
        MULTIQC_FASTP(
            fresh_fastp.mix(prev_fastp).unique { it.name }.collect()
        )
    }

    // ---------------- Mapping / filtering ----------------
    fresh_bams_ch = Channel.empty()
    if (params.run_mapping) {
        // Trimmed reads: fresh from FASTP this run, else read from disk (2_fastp)
        if (params.run_fastp) {
            trimmed_ch = FASTP.out.trimmed
        } else {
            trimmed_ch = Channel
                .fromPath("${params.outdir}/2_fastp/*/*_trimmed_{1,2}.fq.gz")
                .map { f -> tuple(f.name.replaceAll(/_trimmed_[12]\.fq\.gz$/, ''), f) }
                .groupTuple()
                .map { sample, files -> def s = files.sort(); tuple(sample, s[0], s[1]) }
        }

        BUILD_INDEX(file(params.reference))
        // .first() makes the index a value channel so it is reused for every sample
        MAP_FILTER(trimmed_ch, BUILD_INDEX.out.indexed.first())

        fresh_bams_ch = MAP_FILTER.out.filtered.flatMap { sample, bam, bai -> [bam, bai] }
    }

    // ---------------- MetaPop (all BAMs + all counts: fresh + from disk) ----------------
    if (params.run_metapop) {
        // norm.tsv assembled from ALL per-sample counts:
        // fresh this run (if fastp ran) + every count already on disk (0_counts/).
        // collectFile returns the merged file as a channel, so METAPOP now waits
        // for it (proper dependency) instead of reading a path that may not exist yet.
        fresh_counts_ch = params.run_fastp ? COUNT_READS.out.count : Channel.empty()
        prev_counts_ch  = Channel.fromPath("${params.outdir}/0_counts/*_count.tsv")

        norm_ch = fresh_counts_ch
            .mix(prev_counts_ch)
            .unique { it.name }
            .collectFile(name: 'norm.tsv', storeDir: params.outdir, sort: true)

        // filtered BAMs: fresh this run + any already on disk from previous runs
        prev_bams_ch = Channel.fromPath("${params.outdir}/5_mapping/**/*_filtered.bam*")

        all_bams_ch = fresh_bams_ch
            .mix(prev_bams_ch)
            .unique { it.name }
            .collect()

        METAPOP(
            all_bams_ch,
            norm_ch,
            file(params.reference)
        )
    }
}