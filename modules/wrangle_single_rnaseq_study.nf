#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleRnaSeqStudy {
    container "veupathdb/study-wrangler:latest"

    input:
    tuple val(study), path(dat)

//    output:
//    tuple path("install.json"), path("*.cache")

    script:
    """
    echo $study
    """
}
