#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleFileStudy {
    container "veupathdb/study-wrangler:latest"

    input:
    path(dat)
    path(customWrangleScript)

    output:
    tuple path("install.json"), path("*.cache")

    script:
    """
    singleFileCustomWrangle.R $customWrangleScript
    """
}
