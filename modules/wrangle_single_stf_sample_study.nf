#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleStfSampleStudy {
    container "veupathdb/study-wrangler:${params.studyWranglerTag}"
    errorStrategy 'ignore'
    maxForks 10

    input:
    tuple val(study), path(dat), val(extDbNames)

    output:
    tuple val(study), path("install.json"), path("*.cache"), val(extDbNames)

    // NOTE:  for now we can assume study=datasetName
    script:
    """
    singleStfSampleStudyWrangle.R $study
    """
}

