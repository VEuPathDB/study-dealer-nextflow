#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { combined_dnaseq_study } from '../subworkflows/single_study'


// One study per REFERENCE ORGANISM.
//
// Unlike rnaseq, there is no multiDatasetStudy.json: the variables of every
// dnaseq experiment aligned to a given reference organism have already been
// harmonized into a single stf, and the directory layout
// <dnaseqStfDir>/<project>/<organismAbbrev>/ IS the study grouping. Nothing to
// keep in sync by hand, and nothing that can silently disagree with the data.
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


workflow combined_dnaseq_studies {

    main:

    // <dnaseqStfDir>/<project>/<organismAbbrev>/entity-sample.{tsv,yaml}
    stfFiles = Channel.fromPath(params.dnaseqStfPattern)
        .map { file -> keyByProjectAndOrganism(file, params.dnaseqStfDir) }
        .groupTuple()

    // <workflowDataDir>/<project>/<organismAbbrev>/dnaseq/<datasetName>/dnaseqNextflow/analysisDir/results/<sample>/<sample>_alignment_stats.tsv
    alignmentStats = Channel.fromPath(params.filePatterns[params.mode])
        .map { file -> keyByProjectAndOrganism(file, params.workflowDataDir) }
        .groupTuple()

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
    paired = stfFiles.join(alignmentStats, remainder: true)
        .branch { key, stf, stats ->
            complete: stf != null && stats != null
            orphanedStats: stf == null
            unaligned: true
        }

    paired.orphanedStats.subscribe { key, stf, stats ->
        throw new IllegalStateException("${key[0]}/${key[1]} has ${stats.size()} aligned samples but no stf in ${params.dnaseqStfDir} - the alignment stats are the authority on which samples exist, so this is missing metadata, not an optional input")
    }

    paired.unaligned.subscribe { key, stf, stats ->
        log.warn "Skipping ${key[0]}/${key[1]} - stf present but none of its samples were aligned"
    }

    studies = paired.complete
        .map { key, stf, stats ->
            def (projectName, organismAbbrev) = key

            // The study name and the external database name are deliberately the
            // same string: there is exactly one dnaseq study per organism, so
            // giving it a second identity would only create a mapping to
            // maintain. e.g. pvinvinckeiCY_dnaSeqVariations
            def studyName = "${organismAbbrev}_dnaSeqVariations"

            tuple(studyName, organismAbbrev, stf + stats, [studyName])
        }

    combined_dnaseq_study(studies)
}
