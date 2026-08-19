# study-dealer-nextflow

A Nextflow pipeline that discovers raw VEuPathDB study data in a workflow's data directory and "wrangles" it into VDI-ready study artifacts.

## Overview

This pipeline locates study data files produced elsewhere in a VEuPathDB genomic data workflow run and hands them to the `study-wrangler` R package (`veupathdb/study-wrangler`) to be validated and converted into the `install.json` + cache-file artifacts consumed by VDI (VEuPathDB's dataset loading pipeline), which are then loaded into the EDA database via the `InsertEdaStudyFromArtifacts` GUS plugin. In effect, it is the bridge between a workflow's raw analysis output and a loadable EDA study/dataset.

Processing branches on `params.mode`, which selects one of three processing paths:

- **Single-study modes** — `phenotype`, `antibodyArray`, `cellularLocalization`, and `rflp` each wrangle one dataset directory at a time (selected via `params.datasetName`) using a mode-specific custom R wrangling script.
- **`rnaseq`** — discovers RNA-seq count matrices (EBI-style or standard `countsForEda` output) and AI-generated sample metadata across the whole workflow data directory, optionally groups related per-organism datasets into multi-dataset studies via a JSON mapping file, and wrangles each resulting study.
- **`chipChip`** — discovers ChIP-chip sample metadata across the workflow data directory, filters out studies that are already loaded or whose external database isn't registered, and wrangles the remaining studies.

Two additional modes, `dnaseq_chipSeq` and `dnaseq_SNP_CNV`, are recognized by the mode dispatcher but not yet implemented.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Singularity (the default `nextflow.config` includes `conf/singularity.config`; processes run in the `veupathdb/study-wrangler` and `veupathdb/vdi-plugin-wrangler` container images)
- A GUS config file (for `--gusConfigFile`), needed by the artifact-loading step

## Usage

```
nextflow run VEuPathDB/study-dealer-nextflow -r main \
  --mode rnaseq \
  --workflowDataDir /path/to/workflow/data \
  --gusConfigFile /path/to/gus.config \
  --outputDir /path/to/results \
  -resume
```

```
nextflow run VEuPathDB/study-dealer-nextflow -r main \
  --mode phenotype \
  --workflowDataDir /path/to/workflow/data \
  --datasetName tgonGT1_crisprPhenotype_CrisprScreen_RSRC \
  --gusConfigFile /path/to/gus.config \
  --outputDir /path/to/results \
  -resume
```

The pipeline has a single (default) entry point that branches on `params.mode`:

- **`phenotype` / `antibodyArray` / `cellularLocalization` / `rflp`** locate the dataset's data file(s) under `<workflowPath>/<datasetName>` and run them through `wrangleSingleStudy` (`singleStudyWrangle.R <datasetName> <mode>`), which sources the mode/dataset-specific `wrangle()` script, validates the resulting study, and exports it to VDI format.
- **`rnaseq`** globs RNA-seq count files and AI-generated sample metadata across the workflow data directory, tags each file with its study/organism/dataset (using `multiDatasetStudies` to group related per-organism datasets into a single study), filters out studies already present in the EDA database, and wrangles each study group before loading its artifacts via `loadVdiArtifacts`.
- **`chipChip`** globs ChIP-chip sample metadata across the workflow data directory, groups files by study, filters out studies already loaded into EDA or whose external database isn't registered, and wrangles each remaining study group.

## Key Parameters

| Parameter | Default | Description |
|---|---|---|
| `mode` | `rnaseq` | Which processing path to run: `phenotype`, `antibodyArray`, `cellularLocalization`, `rflp`, `rnaseq`, or `chipChip` |
| `workflowDataDir` | `$baseDir/data` | Root data directory being scanned, mirroring a ReFlow workflow's data directory layout |
| `workflowPath` | `${workflowDataDir}/*/*/${mode}` | Glob pattern used to locate per-project/per-organism directories for the selected mode |
| `datasetName` | `""` | Name of the dataset directory to process; required for the single-study modes, optional filter for others |
| `gusConfigFile` | `/path/to/gusConfigFile` | Path to the GUS config file used by the artifact-loading step |
| `multiDatasetStudies` | `data/rnaseq_sample_reannotation/multiDatasetStudy.json` | JSON file mapping individual RNA-seq dataset names to the multi-dataset study they belong to |
| `filePatterns` | (see `nextflow.config`) | Map of glob patterns used to locate each mode's required input files |
| `studyWranglerTag` | `1.0.27` | Docker/Singularity image tag for the `study-wrangler` container (override for local development against an unreleased wrangler version) |
| `outputDir` | `$launchDir/results` | Directory the wrangled study artifacts are published to (one subdirectory per dataset/study) |

## Output

For each dataset or study processed: an `install.json` file and one or more `.cache` files, published to `outputDir` under a per-dataset/study subdirectory. These are the VDI-ready artifacts produced by `study-wrangler`'s `export_to_vdi()`; the `rnaseq` and `chipChip` modes go on to load these artifacts directly into the EDA database via the `InsertEdaStudyFromArtifacts` GUS plugin.
