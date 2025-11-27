process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    container 'ewels/multiqc:v1.19'

    input:
    path '*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data

    script:
    """
    multiqc . -o . --filename multiqc_report.html
    """
}
