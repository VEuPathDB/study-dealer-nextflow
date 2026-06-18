# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow workflow system for processing and "dealing" studies to data wranglers. The system handles multiple data types (RNA-seq, phenotype data, and planned support for DNA-seq variants) and processes them through containerized workflows.

## Essential Commands

### Running Workflows
```bash
# Run RNA-seq workflow (default mode)
nextflow run main.nf

# Run phenotype workflow for specific dataset
nextflow run main.nf --mode phenotype --datasetName tgonGT1_crisprPhenotype_CrisprScreen_RSRC

# Run with custom output directory
nextflow run main.nf --outputDir /path/to/output

# Resume failed workflow
nextflow run main.nf -resume
```

### Development Commands
```bash
# Clean work directory
nextflow clean -f

# View workflow execution report
nextflow log

# Check Nextflow version
nextflow -version
```

## Architecture

### Workflow Structure
- **main.nf**: Entry point that routes to different workflows based on `params.mode`
- **workflows/**: Contains main workflow definitions
  - `multiple_rnaseq_studies.nf`: Processes multiple RNA-seq studies with metadata mapping
  - `single_phenotype_study.nf`: Handles individual phenotype datasets
- **subworkflows/**: Reusable workflow components
- **modules/**: Individual process definitions and utilities

### Data Organization
The system expects data in `data/` directory organized as:
```
data/
├── [ProjectDB]/[organism]/[datatype]/[dataset]/
└── rnaseq_sample_reannotation/
    ├── multiDatasetStudy.json  # Maps datasets to studies
    └── [study]/entity-sample.{tsv,yaml}  # Sample metadata
```

### Key Configuration Parameters
- `params.mode`: Workflow type ("rnaseq", "phenotype", "dnaseq_*")
- `params.workflowDataDir`: Root data directory (default: "$baseDir/data")
- `params.multiDatasetStudies`: JSON file mapping datasets to studies
- `params.filePatterns`: Glob patterns for different file types

### Container Support
Uses Docker by default with user ID mapping for file permissions. Singularity/Apptainer disabled in config.

### Study-Dataset Relationship
The system handles the complexity where multiple datasets can belong to the same study through the `multiDatasetStudy.json` mapping file, allowing proper grouping of related data.