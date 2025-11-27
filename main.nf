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
    IMPORT MODULES
========================================================================================
*/

include { FASTQC_RAW       } from './modules/local/fastqc_raw'
include { FASTP             } from './modules/local/fastp'
include { FASTQC_TRIMMED    } from './modules/local/fastqc_trimmed'
include { BWA_INDEX         } from './modules/local/bwa_index'
include { BWA_MEM           } from './modules/local/bwa_mem'
include { SORT_BAM          } from './modules/local/sort_bam'
include { ALIGNMENT_STATS   } from './modules/local/alignment_stats'
include { MULTIQC           } from './modules/local/multiqc'

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
