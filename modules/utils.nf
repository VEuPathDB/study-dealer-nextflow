#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process addOrganismPrefix {
    container "veupathdb/alpine_bash:latest"

    input:
    tuple val(study), val(orgAbbrev), path(file)


    output:
    tuple val(study), path("${orgAbbrev}_$file")
    
    script:        
    """
    ln -s $file ${orgAbbrev}_$file 
    """
    
}
