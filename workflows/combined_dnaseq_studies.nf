#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { combined_dnaseq_study } from '../subworkflows/single_study'


// One study per REFERENCE ORGANISM.
//
// Unlike rnaseq, there is no multiDatasetStudy.json: the variables of every
// dnaseq experiment aligned to a given reference organism have already been
// harmonized into a single stf, and the directory layout
// <sampleDetailsDir>/<project>/<organismAbbrev>/ IS the study grouping. Nothing
// to keep in sync by hand, and nothing that can silently disagree with the data.
//
// This workflow's only job is to pair each organism's stf with that same
// organism's per-sample alignment stats and hand both to the wrangler.

def keyByProjectAndOrganism(file, root) {
    def matcher = (file.toString() =~ /${root}\/([^\/]+)\/([^\/]+)\//)

    if(!matcher) {
        throw new IllegalArgumentException("cannot parse project/organismAbbrev out of ${file} (root ${root})")
    }

    def projectName = matcher[0][1];
    def organismAbbrev = matcher[0][2];

    return [ [projectName, organismAbbrev], file ]
}


// The alignment stats path carries the datasetName, and for dnaseq the dataset
// name IS the external database name (see the dnaSeqExperimentFrom* dataset
// classes in EbrcModelCommon: the loader registers %DATASET_NAME% as its own
// external database). The stf files deliberately do not carry it -- dataset
// membership is recorded there as the DS_ digest in the dataset_id variable --
// so the paths are where the load spec has to come from.
def keyAlignmentStatsFile(file) {
    def matcher = (file.toString() =~ /${params.workflowDataDir}\/([^\/]+)\/([^\/]+)\/${params.mode}\/([^\/]+)\//)

    if(!matcher) {
        throw new IllegalArgumentException("cannot parse project/organismAbbrev/datasetName out of ${file}")
    }

    def projectName = matcher[0][1];
    def organismAbbrev = matcher[0][2];
    def datasetName = matcher[0][3];

    return [ [projectName, organismAbbrev], file, datasetName ]
}


workflow combined_dnaseq_studies {

    main:

    // <sampleDetailsDir>/<project>/<organismAbbrev>/entity-sample.{tsv,yaml}
    stfFiles = Channel.fromPath(params.filePatterns['dnaseqStf'])
        .map { file -> keyByProjectAndOrganism(file, params.sampleDetailsDir) }
        .groupTuple()

    // <workflowDataDir>/<project>/<organismAbbrev>/dnaseq/<datasetName>/dnaseqNextflow/analysisDir/results/<sample>/<sample>_alignment_stats.tsv
    alignmentStats = Channel.fromPath(params.filePatterns[params.mode])
        .map { file -> keyAlignmentStatsFile(file) }
        .groupTuple(by: 0)

    // The alignment stats are the authority on what exists: they are what the
    // workflow actually aligned. The two half-matches are therefore not
    // symmetric.
    //
    //   aligned but no stf -> fatal. We would be dropping data the workflow
    //       produced because its metadata is missing, which is a gap to fix, not
    //       to skip past.
    //   stf but nothing aligned -> report and skip. A test database legitimately
    //       holds a subset of the aligned data, and only the metadata that
    //       corresponds to it should be processed.
    // Destructured as a single row rather than as named arguments because
    // join(remainder: true) does not pad: an unmatched left entry arrives as
    // [key, stfFiles, null] with no slot for the datasetNames the right hand
    // side would have carried. A four argument closure fails on those rows.
    paired = stfFiles.join(alignmentStats, remainder: true)
        .branch { row ->
            complete: row[1] != null && row[2] != null
            orphanedStats: row[1] == null
            unaligned: true
        }

    paired.orphanedStats.subscribe { row ->
        def key = row[0]
        throw new IllegalStateException("${key[0]}/${key[1]} has ${row[2].size()} aligned samples but no stf in ${params.sampleDetailsDir} - the alignment stats are the authority on which samples exist, so this is missing metadata, not an optional input")
    }

    paired.unaligned.subscribe { row ->
        def key = row[0]
        log.warn "Skipping ${key[0]}/${key[1]} - stf present but none of its samples were aligned"
    }

    studies = paired.complete
        .map { row ->
            def (key, stf, stats, datasetNames) = row
            def (projectName, organismAbbrev) = key

            // One study per organism, but it is derived from every dnaseq dataset
            // aligned to that organism, so every one of their external databases
            // has to be on the load spec -- EDA.StudyExternalDatabaseRelease is
            // many-to-one with the study. Taken from the alignment stats paths
            // rather than from a list somewhere: the datasets that contributed
            // samples are exactly the datasets whose stats we just read.
            def extDbNames = datasetNames.unique().sort()

            def studyName = "${organismAbbrev}_dnaSeqVariations"

            log.info "${projectName}/${organismAbbrev}: ${stats.size()} samples from ${extDbNames.size()} dataset(s)"

            tuple(studyName, organismAbbrev, stf + stats, extDbNames)
        }

    combined_dnaseq_study(studies)
}
