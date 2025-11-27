process ALIGNMENT_STATS {
    tag "$sample_id"
    publishDir "${params.outdir}/${sample_id}/stats", mode: 'copy'
    container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

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
