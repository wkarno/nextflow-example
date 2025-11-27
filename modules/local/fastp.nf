process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/fastp", mode: 'copy'
    container 'biocontainers/fastp:v0.20.1_cv1'

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
