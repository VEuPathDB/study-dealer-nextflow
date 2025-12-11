# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline called "study-dealer-nextflow" that processes biological study data for the VEuPathDB project. It "deals studies to the wrangler" by taking various genomics datasets and processing them through standardized wrangling workflows.

## Architecture

The pipeline follows a mode-based architecture where different data types are processed through specialized workflows:

### Entry Point
- `main.nf`: Main workflow dispatcher that routes execution based on `params.mode`
  - `mode = "rnaseq"`: Processes RNA sequencing studies via `multiple_rnaseq_studies` workflow
  - `mode = "phenotype"`: Processes phenotype datasets via `single_phenotype_study` workflow
  - Additional modes planned: `dnaseq_chipChip`, `dnaseq_chipSeq`, `dnaseq_SNP_CNV`

### Workflow Structure
- `workflows/`: Contains high-level workflow definitions
  - `multiple_rnaseq_studies.nf`: Handles multi-dataset RNA-seq studies, groups datasets by study using JSON mapping
  - `single_phenotype_study.nf`: Processes individual phenotype datasets with custom R wrangling scripts
- `subworkflows/`: Contains reusable workflow components
  - `single_study.nf`: Defines both `single_study` and `single_rnaseq_study` sub-workflows
- `modules/`: Contains individual process definitions
  - `utils.nf`: Contains `addOrganismPrefixAndFilterRows` process
  - `wrangle_single_rnaseq_study.nf`: Main RNA-seq processing using study-wrangler container
  - `wrangle_single_file_study.nf`: Generic file processing
  - `load_vdi_artifacts.nf`: VDI artifact loading (currently unused)

### Data Processing Architecture

The pipeline uses a sophisticated file pattern matching system defined in `nextflow.config`:
- Input data organized by: `{projectName}/{organismAbbrev}/{mode}/{datasetName}/`
- Supports multiple file types: counts files, metadata, WGCNA eigengenes, phenotype data
- Uses JSON mapping (`multiDatasetStudy.json`) to group datasets into studies for RNA-seq mode

### Key Processing Steps for RNA-seq
1. Collect files via glob patterns from `params.filePatterns`
2. Extract metadata (study, organism, dataset) from file paths using regex
3. Mix counts files with AI-generated sample metadata
4. Group by study using `addOrganismPrefixAndFilterRows`
5. Process each study group through `single_rnaseq_study` workflow

## Development Commands

### Running the Pipeline
```bash
# Run RNA-seq mode (default)
nextflow run main.nf

# Run specific mode
nextflow run main.nf --mode phenotype

# Filter to specific dataset
nextflow run main.nf --datasetName "specific_dataset_name"

# Custom output directory
nextflow run main.nf --outputDir /path/to/output
```

### Configuration
- Main config: `nextflow.config`
- Docker config: `conf/docker.config` (enabled by default)
- Singularity config: `conf/singularity.config` (available but disabled)
- Pipeline runs with Docker by default with user mapping and host network access

### Key Parameters
- `gusConfigFile`: Path to GUS configuration file (default: `$launchDir/input/gus.config`)
- `workflowDataDir`: Base directory for input data (default: `$baseDir/data`)
- `mode`: Processing mode (rnaseq, phenotype, etc.)
- `outputDir`: Results output directory (default: `$launchDir/results`)
- `datasetName`: Filter to specific dataset (optional)
- `multiDatasetStudies`: JSON file mapping datasets to studies

### Docker Containers Used
- `veupathdb/study-wrangler:1.0.27`: Main RNA-seq processing container
- `veupathdb/alpine_bash:latest`: Utility container for file filtering

### R Scripts
- `bin/wrangleRNASeq.R`: Main RNA-seq data processing script using study.wrangler library
- `bin/singleFileCustomWrangle.R`: Custom wrangling for phenotype datasets
- Scripts expect specific function signatures (e.g., `wrangle` function for phenotype processing)

## Data Directory Structure

```
data/
├── HostDB/
│   └── hsapREF/
│       └── rnaseq/
├── PlasmoDB/
│   ├── pfal3D7/
│   │   └── rnaseq/
│   └── pberANKA/
│       └── phenotype/
└── rnaseq_sample_reannotation/
    └── multiDatasetStudy.json
```

## Important Implementation Notes

- Pipeline uses Nextflow DSL2
- File path parsing relies on specific regex patterns that match the data directory structure
- RNA-seq mode supports both EBI RNA-seq counts and regular RNA-seq counts patterns
- Multi-dataset studies are defined via JSON configuration for flexible study grouping
- Process isolation uses maxForks = 2 to limit concurrent execution
- All file processing includes organism prefix addition and filtering of rows starting with "__"