process BWA_MEM {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/alignment", mode: 'copy'
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0'

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
