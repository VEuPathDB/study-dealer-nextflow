#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

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

def addFileMetadata(file, datasetToStudyMap) {

    def matcher = (file.toString() =~ /${params.inputDir}\/(.+?)\/(.+?)\/${params.mode}\/(.+?)\//)

    def projectName = matcher[0][1];
    def organismAbbrev = matcher[0][2];
    def datasetName = matcher[0][3];

    def studyMetadata = datasetName; // if we've only got one dataset for a study, use that name for the study
    if(datasetToStudyMap.containsKey(datasetName)) {
        studyMetadata = datasetToStudyMap.get(datasetName)
    }

    // return a Tuple of Meta object + the file
    return [ [study:studyMetadata, projectName:projectName, organismAbbrev:organismAbbrev, datasetName:datasetName ], file ]
}


workflow multiple_rnaseq_studies {

    main:
    datasetToStudyMap = slurpJson("${params.multiDatasetStudies}")

    // Get all raw counts (mix ebi and  non ebi datasets)
    countsFiles = Channel.fromPath(params.inputDir + "/*/*/" + params.mode + "/*/bulkrnaseq/analysisDir/analysis_output/countsForEda*")
        .mix(Channel.fromPath(params.inputDir + "/*/*/" + params.mode + "/*/analysis_output/countsForEda*"))
        .map { file -> addFileMetadata(file, datasetToStudyMap)}


    // TODO:  make a channel for the AI generated metadata files

    // TODO:  join channels.  metadata for the counts should match the file names of the AI files
    
    countsFiles.view()
    


}
