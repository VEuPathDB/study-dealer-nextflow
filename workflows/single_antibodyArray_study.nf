#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { single_study } from '../subworkflows/single_study'

/*
 *  We can grab  datasets from the workflow datadir.  Each dataset will have one txt or tab file. */
/*  The R is from this repo code must include a "wrangle" function (ie a custom function to preprocess this phenotype dataset)
 */
workflow single_antibodyArray_study {

    main:
    antibodyArrayFiles = Channel.fromPath(params.filePatterns['antibodyArray'])
    
    single_study(antibodyArrayFiles.collect())

}
