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

### Basic Usage

```bash
# Run RNA-seq mode (default)
nextflow run main.nf

# Run phenotype mode
nextflow run main.nf --mode phenotype

# Process a specific dataset
nextflow run main.nf --mode phenotype --datasetName "PlasmoDB_Rod_Mal_Phenotype_RSRC"

# Custom output directory
nextflow run main.nf --outputDir /path/to/results
```

### Configuration

Ensure your GUS home directory is set and contains `config/gus.config`:

```bash
nextflow run main.nf --gusHomeDir /path/to/gusHome
```

The pipeline expects the GUS config file at `{gusHomeDir}/config/gus.config`.

## Processing Modes

### RNA-seq Mode (`--mode rnaseq`)

Processes RNA sequencing studies with support for:
- EBI RNA-seq count data
- Standard RNA-seq count matrices
- WGCNA eigengene data
- AI-generated sample metadata
- Multi-organism studies

**Features:**
- Automatic strandedness detection (sense/antisense/unstranded)
- Multi-dataset study grouping via JSON configuration
- Database validation to prevent reprocessing

**Example:**
```bash
nextflow run main.nf --mode rnaseq
```

### Phenotype Mode (`--mode phenotype`)

Processes phenotype datasets using custom R wrangling scripts.

**Example:**
```bash
nextflow run main.nf --mode phenotype --datasetName "PlasmoDB_Rod_Mal_Phenotype_RSRC"
```

### Other Modes

- **Antibody Array** (`--mode antibodyArray`)
- **Cellular Localization** (`--mode cellularLocalization`)
- **RFLP** (`--mode rflp`)

All use the same pattern: custom R scripts in `lib/R/{mode}/{DatasetName}.R`

## Configuration Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `mode` | Processing mode (rnaseq, phenotype, etc.) | `rnaseq` |
| `gusHomeDir` | Path to GUS environment root (must contain config/gus.config) | `/path/to/gusHome` |
| `workflowDataDir` | Base directory for input data | `$baseDir/data` |
| `outputDir` | Results output directory | `$launchDir/results` |
| `datasetName` | Filter to specific dataset | `""` |
| `multiDatasetStudies` | JSON file for study mapping | `${params.sampleDetailsDir}/multiDatasetStudy.json` |
| `studyWranglerTag` | Docker image version | `1.0.27` |

## Input Data Organization

The pipeline expects data organized in a specific directory structure:

```
data/
└── {projectName}/          # e.g., PlasmoDB, HostDB
    └── {organismAbbrev}/   # e.g., pfal3D7, hsapREF
        └── {mode}/         # e.g., rnaseq, phenotype
            └── {datasetName}/
                └── {data files}
```

**Example:**
```
data/PlasmoDB/pfal3D7/rnaseq/pfal3D7_Lee_Gambian_ebi_rnaSeq_RSRC/analysis_output/countsForEda_firststrand.txt
```

## Multi-Dataset Studies (RNA-seq)

For RNA-seq studies spanning multiple organisms, define study groupings in `data/sample_details/multiDatasetStudy.json`:

```json
[
  {
    "study": "HPI_Lee_Gambian",
    "datasets": [
      "hsapREF_Lee_Gambian_ebi_rnaSeq_RSRC",
      "pfal3D7_Lee_Gambian_ebi_rnaSeq_RSRC"
    ]
  }
]
```

Datasets not in the JSON are processed as individual studies.

## Output Structure

Results are published to `{outputDir}/{study_name}/`:

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
