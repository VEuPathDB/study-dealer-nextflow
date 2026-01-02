#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { wrangleSingleStudy } from '../modules/wrangle_single_study'

include { wrangleSingleRnaSeqStudy } from '../modules/wrangle_single_rnaseq_study'


include { loadVdiArtifacts  } from '../modules/load_vdi_artifacts'


workflow single_study {
    take:
    obj

    main:
    artifacts = wrangleSingleStudy(obj)
}



workflow single_rnaseq_study {
    take:
    obj

    main:
    artifacts = wrangleSingleRnaSeqStudy(obj)

    // TODO:  maybe use the vdi container with the gus environment.  We need ApiCommonData built
    //loadVdiArtifacts(artifacts)

    emit:
    strandedness_summary = wrangleSingleRnaSeqStudy.out.strandedness_summary
}
