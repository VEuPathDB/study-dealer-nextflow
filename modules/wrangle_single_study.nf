#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleStudy {
    container "veupathdb/study-wrangler:${params.studyWranglerTag}"

    publishDir params.outputDir + "/" + params.datasetName, mode: 'copy'
    

    input:
    path(dat)

    output:
    tuple val(params.datasetName), path("install.json"), path("*.cache")

    script:
    """
    singleStudyWrangle.R $params.datasetName $params.mode
    """
}
