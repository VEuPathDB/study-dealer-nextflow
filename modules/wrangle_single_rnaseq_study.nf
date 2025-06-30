#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleRnaSeqStudy {
    container "veupathdb/study-wrangler:latest"

    publishDir params.outputDir + "/${study}", mode: 'copy'

    input:
    tuple val(study), path(dat)

    output:
    tuple path("install.json"), path("*.cache")

    script:
    """
    wrangleRNASeq.R 
    """
}
