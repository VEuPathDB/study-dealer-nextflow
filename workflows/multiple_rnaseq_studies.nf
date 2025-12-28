#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

import groovy.json.JsonSlurper

include { single_rnaseq_study } from '../subworkflows/single_study'

include { addOrganismPrefixAndFilterRows } from '../modules/utils'
include { dumpEdaExternalDatabaseNames } from '../modules/utils'
include { dumpAllExternalDatabaseNames } from '../modules/utils'
include { fileMatcher as fileMatcher_one } from '../modules/utils'
include { fileMatcher as fileMatcher_two } from '../modules/utils'
include { reportStudiesWithoutSampleDetails } from '../modules/utils'

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

    return [ studyMetadata, organismAbbrev, file, datasetName ]
}


def addFileMetadataSampleDetails(file) {

    // need to grab the study/dataset from the directory name
    def matcher = (file.toString() =~ /\/([^\/]+?)\/entity-sample.+/)

    def studyOrDatasetDirName = matcher[0][1];

    return [ studyOrDatasetDirName, "SAMPLE_DETAILS", file, "" ]
}


workflow multiple_rnaseq_studies {

    main:
    datasetToStudyMap = slurpJson("${params.multiDatasetStudies}")

    // mix counts files and add to ai sample meta data;  group result by "study"
    inputs = fileMatcher_one(params.filePatterns['ebiRnaSeqCounts'])
        .splitCsv()
        .mix(fileMatcher_two(params.filePatterns['rnaSeqCounts']).splitCsv())
        .map { file -> addFileMetadataToCounts(file, datasetToStudyMap)}
        .mix(Channel.fromPath(params.filePatterns['rnaseqAiMetadata'])
             .map { file -> addFileMetadataSampleDetails(file)  })
//      .mix(Channel.fromPath("/home/jbrestel/tempAiRnaSeqMetaData/*/*.{tsv,yaml}")  


    renamedAndGrouped = addOrganismPrefixAndFilterRows(inputs)
        .groupTuple(by:0)

    // Split studies based on whether they have SAMPLE_DETAILS
    splitStudies = renamedAndGrouped
        .branch {
            studyName, files, datasetNames ->
            def hasSampleDetails = files.any { it.name.endsWith("entity-sample.yaml") }
            withSampleDetails: hasSampleDetails
            withoutSampleDetails: !hasSampleDetails
        }

    // Report studies without SAMPLE_DETAILS
    studiesMissingSampleDetails = reportStudiesWithoutSampleDetails(splitStudies.withoutSampleDetails)

    studiesMissingSampleDetails.collectFile(name: "studies_missing_sample_details.txt", storeDir: params.outputDir)

    // Continue processing only studies with SAMPLE_DETAILS
    studiesWithSampleDetails = splitStudies.withSampleDetails

    // Dump EDA external database names from the database (already in EDA)
    edaDatabaseNamesSet = dumpEdaExternalDatabaseNames()
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
        .collect()
        .map { it.toSet() }

    // Dump all external database names (all databases with releases)
    allDatabaseNamesSet = dumpAllExternalDatabaseNames()
        .splitText()
        .map { it.trim() }
        .filter { it != "" }
        .collect()
        .map { it.toSet() }


    // Filter the channel:
    // 1. Exclude studies that are already in the EDA database
    // 2. Include only studies that have database names in all_external_database_names
    filtered = studiesWithSampleDetails
        .combine(edaDatabaseNamesSet)
        .combine(allDatabaseNamesSet)
        .map { studyName, files, databaseNames, edaDbNamesSet, allDbNamesSet ->
            // Clean up databaseNames list: remove nulls, empty strings, and get unique values
            def cleanedNames = databaseNames.findAll { it != null && it != "" }.unique()
            tuple(studyName, files, cleanedNames, edaDbNamesSet, allDbNamesSet)
        }
        .filter { studyName, files, cleanedNames, edaDbNamesSet, allDbNamesSet ->
            // Check if any of the database names are already in the EDA database
            def matchingEdaNames = cleanedNames.findAll { edaDbNamesSet.contains(it) }

            if (matchingEdaNames) {
                log.info "Skipping study '${studyName}' - database name(s) already in EDA: ${matchingEdaNames.join(', ')}"
                return false
            }

            // Check if any of the database names are in the all external database names
            def matchingAllNames = cleanedNames.findAll { allDbNamesSet.contains(it) }

            if (!matchingAllNames) {
                log.info "Skipping study '${studyName}' - database name(s) not found in external database releases: ${cleanedNames.join(', ')}"
                return false
            }

            // Log that we're processing this study
            //log.info "Processing study '${studyName}' - database name(s) found in external database releases: ${matchingAllNames.join(', ')}"
            return true
        }
        .map { studyName, files, cleanedNames, edaDbNamesSet, allDbNamesSet ->
            tuple(studyName, files, cleanedNames)
        }

    single_rnaseq_study(filtered)
 }
