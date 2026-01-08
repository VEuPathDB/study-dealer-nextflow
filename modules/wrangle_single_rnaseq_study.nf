#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process wrangleSingleRnaSeqStudy {
    container "veupathdb/study-wrangler:${params.studyWranglerTag}"
    errorStrategy 'ignore'
    maxForks 10

    publishDir params.outputDir + "/${study}", mode: 'copy'

    input:
    tuple val(study), path(dat), val(extDbNames)

    output:
    tuple val(study), path("install.json"), path("*.cache"), val(extDbNames)
    path "strandedness_summary.tsv", optional: true, emit: strandedness_summary

    script:
    """
    wrangleRNASeq.R

    # Extract strandedness summary if report exists (skip the SUMMARY prefix field)
    if [ -f strandedness_report.txt ]; then
        grep "^SUMMARY" strandedness_report.txt | awk -F'\t' -v study="${study}" 'BEGIN{OFS="\t"} {print study, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11}' > strandedness_summary.tsv
    fi
    """
}
