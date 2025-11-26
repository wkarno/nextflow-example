#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
========================================================================================
    FASTQ Processing Pipeline
========================================================================================
    FastQC -> fastp -> BWA-MEM -> Sort -> Stats -> MultiQC
========================================================================================
*/

// Print pipeline header
log.info """\
    =============================================
    F A S T Q   P R O C E S S I N G   P I P E L I N E
    =============================================
    samplesheet  : ${params.samplesheet}
    outdir       : ${params.outdir}
    reference    : ${params.reference}
    =============================================
    """
    .stripIndent()

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/fastqc/raw", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_results

    script:
    """
    fastqc -t ${task.cpus} -q ${reads}
    """
}

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/fastp", mode: 'copy'
    container 'biocontainers/fastp:v0.23.4_cv1'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R{1,2}.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.html", emit: html
    path "${sample_id}_fastp.json", emit: json

    script:
    def read1 = reads[0]
    def read2 = reads[1]
    """
    fastp \\
        -i ${read1} \\
        -I ${read2} \\
        -o ${sample_id}_trimmed_R1.fastq.gz \\
        -O ${sample_id}_trimmed_R2.fastq.gz \\
        --thread ${task.cpus} \\
        --qualified_quality_phred ${params.fastp_qualified_quality} \\
        --length_required ${params.fastp_min_length} \\
        --cut_mean_quality ${params.fastp_cut_mean_quality} \\
        --html ${sample_id}_fastp.html \\
        --json ${sample_id}_fastp.json \\
        ${params.fastp_extra_args}
    """
}

process FASTQC_TRIMMED {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/fastqc/trimmed", mode: 'copy'
    container 'biocontainers/fastqc:v0.11.9_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}", emit: fastqc_results

    script:
    """
    fastqc -t ${task.cpus} -q ${reads}
    """
}

process BWA_INDEX {
    tag "$fasta"
    publishDir "${params.outdir}/reference", mode: 'copy'
    container 'biocontainers/bwa:v0.7.17_cv1'

    input:
    path fasta

    output:
    path "${fasta}.*", emit: index

    script:
    """
    bwa index ${fasta}
    """
}

process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'
    container 'biocontainers/bwa:v0.7.17_cv1'

    input:
    tuple val(sample_id), path(reads)
    path reference
    path index

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam

    script:
    """
    bwa mem \\
        -t ${task.cpus} \\
        -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \\
        ${reference} \\
        ${reads[0]} ${reads[1]} \\
        | samtools view -Sb - > ${sample_id}.bam
    """
}

process SORT_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'
    container 'biocontainers/samtools:v1.19.2_cv1'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), emit: sorted_bam
    tuple val(sample_id), path("${sample_id}.sorted.bam.bai"), emit: bai

    script:
    """
    samtools sort \\
        -@ ${task.cpus} \\
        -o ${sample_id}.sorted.bam \\
        ${bam}

    samtools index ${sample_id}.sorted.bam
    """
}

process ALIGNMENT_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/stats", mode: 'copy'
    container 'biocontainers/samtools:v1.19.2_cv1'

    input:
    tuple val(sample_id), path(bam)

    output:
    path "${sample_id}.flagstat.txt", emit: flagstat
    path "${sample_id}.stats.txt", emit: stats
    path "${sample_id}.idxstats.txt", emit: idxstats

    script:
    """
    samtools flagstat ${bam} > ${sample_id}.flagstat.txt
    samtools stats ${bam} > ${sample_id}.stats.txt
    samtools idxstats ${bam} > ${sample_id}.idxstats.txt
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'ewels/multiqc:v1.19'

    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . -o . --filename multiqc_report.html
    """
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

def parseSamplesheet(samplesheet) {
    def samples = []

    // Read and parse CSV
    file(samplesheet).eachLine { line ->
        // Skip header
        if (line.startsWith('sampleID') || line.startsWith('#')) {
            return
        }

        def fields = line.split(',')
        if (fields.size() >= 3) {
            def sample_id = fields[0].trim()
            def fastq1 = file(fields[1].trim())
            def fastq2 = file(fields[2].trim())

            // Validate files exist
            if (!fastq1.exists()) {
                error "FASTQ file does not exist: ${fastq1}"
            }
            if (!fastq2.exists()) {
                error "FASTQ file does not exist: ${fastq2}"
            }

            samples.add([sample_id, [fastq1, fastq2]])
        }
    }

    if (samples.isEmpty()) {
        error "No valid samples found in samplesheet: ${samplesheet}"
    }

    return samples
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    // Parse samplesheet
    samples_ch = Channel.fromList(parseSamplesheet(params.samplesheet))

    // Reference genome
    reference = file(params.reference)

    // Index reference if index files don't exist
    if (!file("${params.reference}.bwt").exists()) {
        BWA_INDEX(reference)
        index_ch = BWA_INDEX.out.index
    } else {
        index_ch = Channel.fromPath("${params.reference}.*").collect()
    }

    // FastQC on raw reads
    FASTQC_RAW(samples_ch)

    // Trim reads with fastp
    FASTP(samples_ch)

    // FastQC on trimmed reads
    FASTQC_TRIMMED(FASTP.out.trimmed_reads)

    // Align with BWA-MEM
    BWA_MEM(FASTP.out.trimmed_reads, reference, index_ch)

    // Sort BAM files
    SORT_BAM(BWA_MEM.out.bam)

    // Calculate alignment statistics
    ALIGNMENT_STATS(SORT_BAM.out.sorted_bam)

    // Collect all QC outputs for MultiQC
    multiqc_input = Channel.empty()
        .mix(FASTQC_RAW.out.fastqc_results.collect().ifEmpty([]))
        .mix(FASTP.out.json.collect().ifEmpty([]))
        .mix(FASTQC_TRIMMED.out.fastqc_results.collect().ifEmpty([]))
        .mix(ALIGNMENT_STATS.out.flagstat.collect().ifEmpty([]))
        .mix(ALIGNMENT_STATS.out.stats.collect().ifEmpty([]))
        .collect()

    // Run MultiQC
    MULTIQC(multiqc_input)
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """\
        =============================================
        Pipeline completed at: ${workflow.complete}
        Duration            : ${workflow.duration}
        Success             : ${workflow.success}
        Exit status         : ${workflow.exitStatus}
        Output directory    : ${params.outdir}
        =============================================
        """
        .stripIndent()
}
