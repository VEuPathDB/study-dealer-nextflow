#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process loadVdiArtifacts {
    container "veupathdb/vdi-plugin-wrangler:latest"

    maxForks = 5
    
    input:
    tuple val(study), path(installJson), path(cache), val(extDbNames)
    path(strandSummary), optional true
    
    script:
    """
    ga ApiCommonData::Load::Plugin::InsertEdaEntityTypeAndStudy \
        --installJson install.json  \
        --studyFile study.cache \
        --entityTypeGraphFile entitytypegraph.cache  \
        --extDbRlsSpec "${extDbNames.collect{ "${it}|%" }.join(',')}" \
        --gusConfigFile ${params.gusConfigFile} \
        --commit
    """

}
