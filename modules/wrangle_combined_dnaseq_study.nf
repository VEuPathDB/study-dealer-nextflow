#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleCombinedDnaseqStudy {
    container "veupathdb/study-wrangler:${params.studyWranglerTag}"

    publishDir { "${params.outputDir}/${study}" }, mode: 'copy'

    maxForks 4

    input:
    tuple val(study), val(orgAbbrev), path(dat), val(extDbNames)

    output:
    tuple val(study), path("install.json"), path("*.cache"), val(extDbNames)

    script:
    """
    combinedDnaseqStudyWrangle.R $study $orgAbbrev
    """
}
