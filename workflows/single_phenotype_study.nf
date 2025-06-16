#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { single_study } from '../subworkflows/single_study'

/*
 *  We can grab phenotype datasets from the workflow datadir.  Each dataset will have one txt or tab file and
 *  an R code.  This r code must include a "wrangle" function (ie a custom function to preprocess this phenotype dataset)
 *  we can grab the orgAbbrev,projectName, and datasetName from the path to be stored as meta data. (the datasetName is needed 
 *  when calling the GUS plugin 
 */
workflow single_phenotype_study {

    main:
    phenotypeFile = Channel.fromPath(params.inputDir + "/*/*/" + params.mode + "/" + params.datasetName + "/*.{txt,tab}")
        .map { file -> def matcher = (file.toString() =~ /${params.inputDir}\/(.+?)\/(.+?)\/${params.mode}\/(.+?)\/(.+)/)
              [ [projectName:matcher[0][1], organismAbbrev:matcher[0][2], datasetName:matcher[0][3] ], file ]
        }

    customWrangleScript = Channel.fromPath(params.inputDir + "/*/*/" + params.mode + "/" + params.datasetName + "/*.{R,r}")

    single_study(phenotypeFile, customWrangleScript)

}
