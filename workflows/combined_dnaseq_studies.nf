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

    // An organism is only wrangled when we have both halves. An stf with no
    // alignment stats has nothing to load, and alignment stats with no stf would
    // produce a study with no sample characteristics. Either way it is a missing
    // input, not something to quietly proceed with, so the half-matches are
    // warned about instead of being dropped invisibly by an inner join.
    paired = stfFiles.join(alignmentStats, remainder: true)
        .branch { key, stf, stats ->
            complete: stf != null && stats != null
            incomplete: true
        }

    paired.incomplete.subscribe { key, stf, stats ->
        def what = stf == null ? "alignment stats but no stf" : "stf but no alignment stats"
        log.warn "Skipping ${key[0]}/${key[1]} - found ${what}"
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
