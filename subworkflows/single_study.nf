#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { wrangleSingleFileStudy } from '../modules/wrangle_single_file_study'

include { wrangleSingleRnaSeqStudy } from '../modules/wrangle_single_rnaseq_study'


workflow single_study {
    take:
    dataFile
    customWrangleScript

    main:
    artifacts = wrangleSingleFileStudy(dataFile, customWrangleScript)

    // TODO:  maybe use the vdi container with the gus environment.  We need ApiCommonData built
    //loadVdiArtifacts(artifacts)
}


workflow single_rnaseq_study {
    take:
    obj

    main:
    artifacts = wrangleSingleRnaSeqStudy(obj)

    // TODO:  maybe use the vdi container with the gus environment.  We need ApiCommonData built
    //loadVdiArtifacts(artifacts)
}
