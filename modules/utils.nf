#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def addFileMetadataSampleDetails(file) {

    // need to grab the study/dataset from the directory name
    def matcher = (file.toString() =~ /\/([^\/]+?)\/entity-sample.+/)

    def studyOrDatasetDirName = matcher[0][1];

    return [ studyOrDatasetDirName, "SAMPLE_DETAILS", file, "" ]
}



process addOrganismPrefixAndFilterRows {
    container "veupathdb/alpine_bash:latest"

    maxForks 10

    input:
    tuple val(study), val(orgAbbrev), path(file), val(datasetName)


    output:
    tuple val(study), path("${orgAbbrev}_$file"), val(datasetName)

    script:
    """
    grep -v "^__" $file >${orgAbbrev}_$file
    """

}


process fileMatcher {
    container "veupathdb/alpine_bash:latest"

    input:
    val(globPattern)

    output:
    path("matchedFiles.txt")    

    script:
    """
    ls $globPattern >matchedFiles.txt
    """

}



process dumpEdaExternalDatabaseNames {
    container "veupathdb/vdi-plugin-wrangler:latest"
    
    output:
    path "external_database_names.txt"

    script:
    """
    dumpEdaExternalDatabaseNames.pl \
        --gusConfigFile ${params.gusConfigFile} \
        --outputFile external_database_names.txt \
        --edaOnly \
        --verbose
    """
}


process dumpAllExternalDatabaseNames {
    container "veupathdb/vdi-plugin-wrangler:latest"

    output:
    path "all_external_database_names.txt"

    script:
    """
    dumpEdaExternalDatabaseNames.pl \
        --gusConfigFile ${params.gusConfigFile} \
        --outputFile all_external_database_names.txt \
        --verbose
    """
}


process reportStudiesWithoutSampleDetails {
    container "veupathdb/alpine_bash:latest"

    maxForks 10

    input:
    tuple val(studyName), path(files), val(datasetNames)

    output:
    path "no_sample_details.txt"

    script:
    def datasetList = datasetNames.findAll { it != "" }.unique().join(", ")
    """
    echo "Study: ${studyName}" >> no_sample_details.txt
    echo "Datasets: ${datasetList}" >> no_sample_details.txt
    echo "Files: ${files.join(', ')}" >> no_sample_details.txt
    echo "---" >> no_sample_details.txt
    """
}
