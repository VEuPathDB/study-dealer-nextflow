#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process loadVdiArtifacts {
    container "veupathdb/vdi-plugin-wrangler:latest"

    maxForks = 5
    
    input:
    tuple val(study), path(installJson), path(cache), val(extDbNames)
    
    script:
    """
    ga ApiCommonData::Load::Plugin::InsertEdaStudyFromArtifacts \
        --inputDirectory \$PWD  \
        --outputDirectory \$PWD/loadArtifactsOut  \
        --extDbRlsSpec "${extDbNames.collect{ "${it}|%" }.join(',')}" \
        --gusConfigFile ${params.gusConfigFile} \
        --commit
    """

}
