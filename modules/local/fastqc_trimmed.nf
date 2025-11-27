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
