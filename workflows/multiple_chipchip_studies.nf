#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { single_stf_sample_study } from '../subworkflows/single_study'

include { addFileMetadataSampleDetails } from '../modules/utils'
include { dumpEdaExternalDatabaseNames } from '../modules/utils'
include { dumpAllExternalDatabaseNames } from '../modules/utils'

workflow multiple_chipchip_studies {

    main:
    
    inputs = Channel.fromPath(params.filePatterns['chipChipMetadata'])
        .map { file -> addFileMetadataSampleDetails(file)  }
        .groupTuple(by:0)


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
    filtered = inputs
        .combine(edaDatabaseNamesSet)
        .combine(allDatabaseNamesSet)
        .map { studyName, prefixes, files, databaseNames, edaDbNamesSet, allDbNamesSet ->
             // Clean up databaseNames list: remove nulls, empty strings, and get unique values
            tuple(studyName, files, [studyName], edaDbNamesSet, allDbNamesSet)
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


    single_stf_sample_study(filtered)
 }
