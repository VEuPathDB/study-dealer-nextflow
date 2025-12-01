#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { multiple_rnaseq_studies } from './workflows/multiple_rnaseq_studies'

include { single_study } from './subworkflows/single_study'


workflow {

    main:

    if(params.mode == "phenotype") {
        phenotypeFiles = Channel.fromPath(params.filePatterns['phenotype'])
        single_study(phenotypeFiles.collect())
    }

    if(params.mode == "antibodyArray") {
        antibodyArrayFiles = Channel.fromPath(params.filePatterns['antibodyArray'])
        single_study(antibodyArrayFiles.collect())
    }

    if(params.mode == "cellularLocalization") {
        cellularLocalizationFiles = Channel.fromPath(params.filePatterns['cellularLocalization'])
        single_study(cellularLocalizationFiles.collect())
    }


    if(params.mode == "rnaseq") {
        multiple_rnaseq_studies()
    }

    if(params.mode == "dnaseq_chipChip") {
        // TODO 
    }
    if(params.mode == "dnaseq_chipSeq") {
        // TODO
    }
    if(params.mode == "dnaseq_SNP_CNV") {
        // TODO
    }


}
