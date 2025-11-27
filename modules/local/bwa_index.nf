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
