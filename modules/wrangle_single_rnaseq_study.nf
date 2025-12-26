#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleRnaSeqStudy {
    container "veupathdb/study-wrangler:${params.studyWranglerTag}"

    publishDir params.outputDir + "/${study}", mode: 'copy'

    input:
    tuple val(study), path(dat), val(extDbNames)

    output:
    tuple val(study), path("install.json"), path("*.cache")

    script:
    """
    wrangleRNASeq.R 
    """
}
