#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { multiple_rnaseq_studies } from './workflows/multiple_rnaseq_studies'

include { single_study } from './subworkflows/single_study'


workflow {

    main:



    if(params.mode == "phenotype" ||
       params.mode == "antibodyArray" ||
       params.mode == "cellularLocalization" ||
       params.mode == "rflp") {
        singleStudyFiles = Channel.fromPath(params.filePatterns[params.mode])
        single_study(singleStudyFiles.collect())
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
