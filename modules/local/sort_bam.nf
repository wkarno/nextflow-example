process SORT_BAM {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

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
