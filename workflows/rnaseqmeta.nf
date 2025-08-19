/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { STAR_GENOMEGENERATE    } from '../modules/nf-core/star/genomegenerate/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { SALMON_INDEX           } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT           } from '../modules/nf-core/salmon/quant/main'
include { SUBREAD_FEATURECOUNTS  } from '../modules/nf-core/subread/featurecounts/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BATCH_CORRECTION       } from '../modules/local/batch_correction'
include { DIFFERENTIAL_EXPRESSION} from '../modules/local/differential_expression'
include { META_ANALYSIS          } from '../modules/local/meta_analysis'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseqmeta_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQMETA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Prepare reference files
    //
    ch_fasta = params.fasta ? Channel.fromPath(params.fasta, checkIfExists: true).collect() : Channel.empty()
    ch_gtf = params.gtf ? Channel.fromPath(params.gtf, checkIfExists: true).collect() : Channel.empty()
    ch_transcript_fasta = params.transcript_fasta ? Channel.fromPath(params.transcript_fasta, checkIfExists: true).collect() : Channel.empty()

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run FASTP for read trimming and filtering
    //
    FASTP (
        ch_samplesheet,
        [],  // adapter_fasta
        false, // discard_trimmed_pass
        false, // save_trimmed_fail
        false  // save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // Conditional workflow: STAR alignment or Salmon quantification
    //
    if (params.aligner == 'star') {
        //
        // MODULE: Generate STAR genome index
        //
        STAR_GENOMEGENERATE (
            ch_fasta.map { fasta -> [[:], fasta] },
            ch_gtf.map { gtf -> [[:], gtf] }
        )
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)

        //
        // MODULE: Align reads with STAR
        //
        STAR_ALIGN (
            FASTP.out.reads,
            STAR_GENOMEGENERATE.out.index,
            ch_gtf.map { gtf -> [[:], gtf] },
            false, // star_ignore_sjdbgtf
            '', // seq_platform
            ''  // seq_center
        )
        ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})
        ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)

        //
        // MODULE: Count features with featureCounts
        //
        SUBREAD_FEATURECOUNTS (
            STAR_ALIGN.out.bam.map { meta, bam -> [meta, bam, params.gtf] }
        )
        ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})
        ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)

    } else {
        //
        // MODULE: Generate Salmon index
        //
        SALMON_INDEX (
            ch_fasta,
            ch_transcript_fasta
        )
        ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)

        //
        // MODULE: Quantify with Salmon
        //
        SALMON_QUANT (
            FASTP.out.reads,
            SALMON_INDEX.out.index,
            ch_gtf,
            ch_transcript_fasta,
            false, // alignment_mode
            params.lib_type ?: 'A'
        )
        ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.results.collect{it[1]})
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)
    }

    //
    // Prepare count matrices for meta-analysis
    //
    ch_count_matrices = params.aligner == 'star' ? 
        SUBREAD_FEATURECOUNTS.out.counts : 
        SALMON_QUANT.out.results

    //
    // MODULE: Batch effect correction
    //
    if (params.perform_batch_correction) {
        BATCH_CORRECTION (
            ch_count_matrices,
            params.sample_info
        )
        ch_corrected_counts = BATCH_CORRECTION.out.corrected_counts
        ch_versions = ch_versions.mix(BATCH_CORRECTION.out.versions)
    } else {
        ch_corrected_counts = ch_count_matrices
    }

    //
    // MODULE: Differential expression analysis per study
    //
    DIFFERENTIAL_EXPRESSION (
        ch_corrected_counts,
        params.sample_info
    )
    ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_EXPRESSION.out.plots.collect{it[1]})
    ch_versions = ch_versions.mix(DIFFERENTIAL_EXPRESSION.out.versions)

    //
    // MODULE: Meta-analysis across studies
    //
    if (params.perform_meta_analysis) {
        // Group DE results by comparison
        ch_grouped_results = DIFFERENTIAL_EXPRESSION.out.results
            .map { meta, results -> [meta.comparison ?: 'default', results] }
            .groupTuple()
            .map { comparison, results_list -> 
                [[id: comparison], results_list] 
            }

        META_ANALYSIS (
            ch_grouped_results,
            params.sample_info
        )
        ch_multiqc_files = ch_multiqc_files.mix(META_ANALYSIS.out.plots.collect{it[1]})
        ch_versions = ch_versions.mix(META_ANALYSIS.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'rnaseqmeta_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
