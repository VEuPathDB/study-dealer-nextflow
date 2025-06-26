#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { single_rnaseq_study } from '../subworkflows/single_study'

include { addOrganismPrefix } from '../modules/utils'

def slurpJson(jsonFilePath) {
    def jsonSlurper = new JsonSlurper()
    def jsonArray = jsonSlurper.parse(new File(jsonFilePath))

    def rv = [:]

    for(obj in jsonArray) {
        def study = obj.get("study");
        def datasets = obj.get("datasets");
        for(d in datasets) {
            rv.put(d, study)
        }
    }
    return(rv)
} 

def addFileMetadataToCounts(file, datasetToStudyMap) {

    // need to grab the datasetName and the organismAbbrev from the file path
    def matcher = (file.toString() =~ /${params.workflowDataDir}\/(.+?)\/(.+?)\/${params.mode}\/(.+?)\//)

    def projectName = matcher[0][1];
    def organismAbbrev = matcher[0][2];
    def datasetName = matcher[0][3];

    def studyMetadata = datasetName; // if we've only got one dataset for a study, use that name for the study
    if(datasetToStudyMap.containsKey(datasetName)) {
        studyMetadata = datasetToStudyMap.get(datasetName)
    }

    return [ studyMetadata, organismAbbrev, file ]
}



workflow multiple_rnaseq_studies {

    main:
    datasetToStudyMap = slurpJson("${params.multiDatasetStudies}")

    // mix counts files and add to ai sample meta data;  group result by "study"
    inputs = Channel.fromPath(params.filePatterns['ebiRnaSeqCounts'])
        .mix(Channel.fromPath(params.filePatterns['rnaSeqCounts']))
        .map { file -> addFileMetadataToCounts(file, datasetToStudyMap)}
        .mix(Channel.fromPath(params.filePatterns['rnaseqAiMetadata'])
             .map { file -> [ file.baseName, "SAMPLE_DETAILS", file ] })


    renamedAndGrouped = addOrganismPrefix(inputs)
        .groupTuple(by:0)

    single_rnaseq_study(renamedAndGrouped)
}
