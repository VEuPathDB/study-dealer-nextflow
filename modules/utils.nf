#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


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
    container "veupathdb/vdi-plugin-wrangler"

    output:
    path "external_database_names.txt"

    script:
    """
    dumpEdaExternalDatabaseNames.pl \
        --gusConfigFile ${params.gusHomeDir}/config/gus.config \
        --outputFile external_database_names.txt \
        --edaOnly \
        --verbose
    """
}


process dumpAllExternalDatabaseNames {
    container "veupathdb/vdi-plugin-wrangler"

    output:
    path "all_external_database_names.txt"

    script:
    """
    dumpEdaExternalDatabaseNames.pl \
        --gusConfigFile ${params.gusHomeDir}/config/gus.config \
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
