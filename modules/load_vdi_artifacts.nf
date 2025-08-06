#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process loadVdiArtifacts {
    container "jbrestel/vdi-isasimple"

    input:
    tuple val(datasetName), path(installJson), path(cache)


    script:
    """
    ga ApiCommonData::Load::Plugin::InsertEdaEntityTypeAndStudy \
        --installJson install.json  \
        --studyFile study.cache \
        --entityTypeGraphFile entitytypegraph.cache  \
        --extDbRlsSpec "$datasetName|%" \
        --gusConfigFile $params.gusConfigFile \
        --commit
    """

}
