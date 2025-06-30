#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { single_phenotype_study } from './workflows/single_phenotype_study'

include { multiple_rnaseq_studies } from './workflows/multiple_rnaseq_studies'

workflow {

    main:

    if(params.mode == "phenotype") {
        single_phenotype_study()
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
