#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process addOrganismPrefixAndFilterRows {
    container "veupathdb/alpine_bash:latest"

    input:
    tuple val(study), val(orgAbbrev), path(file), val(datasetName)


    output:
    tuple val(study), path("${orgAbbrev}_$file"), val(datasetName)

    script:
    """
    grep -v "^__" $file >${orgAbbrev}_$file
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
