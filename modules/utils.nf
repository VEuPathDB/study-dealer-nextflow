#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process addOrganismPrefixAndFilterRows {
    container "veupathdb/alpine_bash:latest"

    input:
    tuple val(study), val(orgAbbrev), path(file)


    output:
    tuple val(study), path("${orgAbbrev}_$file")
    
    script:        
    """
    grep -v "^__" $file >${orgAbbrev}_$file 
    """
    
}
