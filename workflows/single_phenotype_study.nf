#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { single_study } from '../subworkflows/single_study'

/*
 *  We can grab phenotype datasets from the workflow datadir.  Each dataset will have one txt or tab file and
 *  an R code.  This r code must include a "wrangle" function (ie a custom function to preprocess this phenotype dataset)
 */
workflow single_phenotype_study {

    main:
    phenotypeFile = Channel.fromPath(params.filePatterns['phenotype'])

//    customWrangleScript = Channel.fromPath(params.filePatterns['phenotypeScript'])

    single_study(phenotypeFile)

}
