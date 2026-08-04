#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process loadVdiArtifacts {
    container "veupathdb/vdi-plugin-wrangler:latest"

    tag "plugin"
    
    maxForks = 5
    
    input:
    tuple val(study), path(installJson), path(cache), val(extDbNames)

    // params.dryRunLoad echoes the command instead of running it, so a new mode
    // can be validated end to end -- artifacts built, ext db spec assembled --
    // without touching the database.
    script:
    def dryRun = params.dryRunLoad ? "echo " : ""
    """
    ${dryRun}ga ApiCommonData::Load::Plugin::InsertEdaStudyFromArtifacts \
        --inputDirectory \$PWD  \
        --outputDirectory \$PWD/loadArtifactsOut  \
        --extDbRlsSpec "${extDbNames.collect{ "${it}|%" }.join(',')}" \
        --gusConfigFile ${params.gusConfigFile} \
        --commit
    """

}
