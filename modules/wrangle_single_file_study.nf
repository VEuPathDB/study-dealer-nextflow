#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleFileStudy {
    container "veupathdb/study-wrangler:1.0.16"

    publishDir params.outputDir + "/" + params.datasetName, mode: 'copy'
    

    input:
    path(dat)
    path(customWrangleScript)

    output:
    tuple val(params.datasetName), path("install.json"), path("*.cache")

    script:
    """
    singleFileCustomWrangle.R $customWrangleScript
    """
}
