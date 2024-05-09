#!/usr/bin/env nextflow
/* 
 * Enables Nextflow DSL 2
 */

nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    =======================================================
    SampoAbfFinemap v${workflow.manifest.version}
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
        nextflow run main.nf \
        --InputDir \
        --SnpRefFile \
        -profile singularity,slurm\
        -resume

    Mandatory arguments:
    --InputDir              Directory with parquet files.
    --SnpRefFile            SNP reference file with alleles information.

    Optional arguments:
    --OutputDir             Output directory.
    """.stripIndent()
}

// Define parameters for input and output directories

// Include module files
include { FinemapAbf } from './modules/AbfFinemapping.nf'

// Create Channel for initial input data files
Channel
    .fromPath("${params.inputDir}/*.parquet.snappy", checkIfExists: true)
    .ifEmpty { exit 1, "Sumstats directory is empty!" }
    .set { input_files_ch }
Channel
    .fromPath("${params.SnpRefFile}/*", checkIfExists: true)
    .ifEmpty { exit 1, "SNP reference file is missing!" }
    .set { snp_ref_ch }


log.info """=======================================================
SampoAbfFinemap v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']            = 'SampoAbfFinemap'
summary['Pipeline Version']         = workflow.manifest.version
summary['Working dir']              = workflow.workDir
summary['Container Engine']         = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']             = "$HOME"
summary['Current user']             = "$USER"
summary['Current path']             = "$PWD"
summary['Working dir']              = workflow.workDir
summary['Script dir']               = workflow.projectDir
summary['Config Profile']           = workflow.profile
summary['Sumstats folder']          = params.inputDir
summary['Output folder']            = params.outputDir
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "========================================="


// Define the workflow
workflow {
    FinemapAbf(input_files_ch.combine(snp_ref_ch))
    FinemapAbf.out.map { it[1] }.collectFile(name: 'FinemappingResultsSummary.txt', keepHeader: true, sort: true, storeDir: "${params.outputDir}")
    FinemapAbf.out.map { it[0] }.collectFile(storeDir: "${params.outputDir}")
} 


