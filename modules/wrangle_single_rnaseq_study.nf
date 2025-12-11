#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleRnaSeqStudy {
    container "veupathdb/study-wrangler:1.0.27"

//    publishDir params.outputDir + "/${study}", mode: 'copy'

    input:
    tuple val(study), path(dat)

    output:
    tuple val(study), path("install.json"), path("*.cache")

    script:
    """
    wrangleRNASeq.R 
    """
}
