# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline called "study-dealer-nextflow" that processes biological study data for the VEuPathDB project. It "deals studies to the wrangler" by taking various genomics datasets and processing them through standardized wrangling workflows.

## Architecture

The pipeline follows a mode-based architecture where different data types are processed through specialized workflows:

### Entry Point
- `main.nf`: Main workflow dispatcher that routes execution based on `params.mode`
  - `mode = "rnaseq"`: Processes RNA sequencing studies via `multiple_rnaseq_studies` workflow
  - `mode = "phenotype"`, `antibodyArray`, `cellularLocalization`, `rflp`: Process via `single_study` subworkflow
  - Additional modes planned: `dnaseq_chipChip`, `dnaseq_chipSeq`, `dnaseq_SNP_CNV`

### Workflow Structure
- `workflows/`: Contains high-level workflow definitions
  - `multiple_rnaseq_studies.nf`: Handles multi-dataset RNA-seq studies, groups datasets by study using JSON mapping
- `subworkflows/`: Contains reusable workflow components
  - `single_study.nf`: Defines both `single_study` and `single_rnaseq_study` sub-workflows
- `modules/`: Contains individual process definitions
  - `utils.nf`: Contains `addOrganismPrefixAndFilterRows` process
  - `wrangle_single_rnaseq_study.nf`: Main RNA-seq processing using study-wrangler container
  - `wrangle_single_study.nf`: Generic file processing with mode-specific R script resolution
  - `load_vdi_artifacts.nf`: VDI artifact loading (currently unused)
- `lib/R/`: Custom R wrangling scripts organized by mode
  - `{mode}/{datasetName}.R`: Each must implement a `wrangle()` function
  - Scripts are dynamically loaded via MY_R_LIB environment variable
  - Examples: `phenotype/`, `antibodyArray/`, `cellularLocalization/`, `rflp/`

### Data Processing Architecture

The pipeline uses a sophisticated file pattern matching system defined in `nextflow.config`:
- Input data organized by: `{projectName}/{organismAbbrev}/{mode}/{datasetName}/`
- Supports multiple file types: counts files, metadata, WGCNA eigengenes, phenotype data
- Uses JSON mapping (`multiDatasetStudy.json`) to group datasets into studies for RNA-seq mode

### Key Processing Steps for RNA-seq
1. Collect files via glob patterns from `params.filePatterns`
2. Extract metadata (study, organism, dataset) from file paths using regex patterns
3. Mix counts files with AI-generated sample metadata
4. Add organism prefix and filter rows starting with "__"
5. Group by study name using `groupTuple(by:0)`
6. Validate against database names (only process if NOT in EDA but IS in external database releases)
7. Process each study group through `single_rnaseq_study` workflow

### Important Architectural Patterns

#### Regex-Based Path Parsing
File paths are parsed to extract metadata using regex that matches directory structure:
```groovy
${params.workflowDataDir}/(.+?)/(.+?)/${params.mode}/(.+?)/
// Captures: projectName, organismAbbrev, datasetName
```
This decouples metadata from filenames and enables flexible file layouts.

#### JSON-Based Study Mapping
Multi-dataset studies are defined in `multiDatasetStudy.json`:
```json
[{"study": "HPI_Lee_Gambian", "datasets": ["hsapREF_Lee_Gambian_ebi_rnaSeq_RSRC", "pfal3D7_Lee_Gambian_ebi_rnaSeq_RSRC"]}]
```
Default behavior: dataset name = study name if not in JSON.

#### Dynamic R Script Resolution
Custom wrangling scripts are loaded via environment variable pattern:
- `MY_R_LIB` environment variable points to `lib/R/`
- Scripts located at: `${MY_R_LIB}/{mode}/{datasetName}.R`
- Each script runs in isolated environment and must implement `wrangle()` function
- Enables adding new datasets without code changes

#### Database Name Validation
Recent addition (PR #5) to prevent reprocessing:
- Dumps EDA database names and all external database names
- Filters studies: only process if database name NOT in EDA AND IS in external databases
- Prevents redundant wrangling of already-processed studies

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

# Use custom GUS config
nextflow run main.nf --gusHomeDir /path/to/gus
```

### Configuration
- Main config: `nextflow.config`
- Docker config: `conf/docker.config` (enabled by default)
- Singularity config: `conf/singularity.config` (available but disabled)
- Pipeline runs with Docker by default with user mapping (`-u $(id -u):$(id -g)`) and host network access

### Key Parameters
- `gusHomeDir`: Path to GUS environment root (must contain config/gus.config) (default: `/path/to/gusHome`)
- `workflowDataDir`: Base directory for input data (default: `$baseDir/data`)
- `mode`: Processing mode (rnaseq, phenotype, etc.) (default: `rnaseq`)
- `outputDir`: Results output directory (default: `$launchDir/results`)
- `datasetName`: Filter to specific dataset (default: `""`)
- `multiDatasetStudies`: JSON file mapping datasets to studies
- `studyWranglerTag`: Docker image version (default: `1.0.27`)

### Docker Containers Used
- `veupathdb/study-wrangler:{params.studyWranglerTag}`: Main processing container with R and study.wrangler library
- `veupathdb/alpine_bash:latest`: Utility container for file filtering
- `veupathdb/vdi-plugin-wrangler`: For dumping external database names

### R Scripts
- `bin/wrangleRNASeq.R`: Main RNA-seq data processing script using study.wrangler library
  - Handles strandedness detection (sense/antisense/unstranded)
  - Supports WGCNA eigengene processing
  - Includes benchmarking utilities
- `bin/singleStudyWrangle.R`: Wrapper for phenotype processing
  - Dynamically loads mode-specific scripts from `lib/R/{mode}/{datasetName}.R`
  - Expects `wrangle()` function that returns a study object
- `lib/R/{mode}/{datasetName}.R`: Custom wrangling scripts
  - Must implement `wrangle()` function
  - Uses study.wrangler library for entity/study creation
  - Validates with "baseline" and "eda" profiles
  - Exports to VDI format (install.json + *.cache files)

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

Expected file structure: `{workflowDataDir}/{projectName}/{organismAbbrev}/{mode}/{datasetName}/{files}`

## Important Implementation Notes

### Core Concepts
- Pipeline uses Nextflow DSL2
- File path parsing relies on specific regex patterns that match the data directory structure
- RNA-seq mode supports both EBI RNA-seq counts and regular RNA-seq counts patterns
- Multi-dataset studies are defined via JSON configuration for flexible study grouping
- Process isolation uses `maxForks = 2` to limit concurrent execution (prevents memory overload)
- All file processing includes organism prefix addition and filtering of rows starting with "__"

### Channel Operations
RNA-seq mode uses complex channel mixing:
1. Mix multiple file patterns (EBI counts, standard counts, WGCNA, AI metadata)
2. Map files to tuples: `[study_name, organism_prefix, file_path, dataset_name]`
3. Add organism prefix and filter rows
4. Group by study using `groupTuple(by:0)`
5. Filter based on database name validation

### Study Wrangler Integration
All modes use the study.wrangler R library pattern:
1. Read data into entity (tibble-like object)
2. Set entity metadata (name, display_name)
3. Configure columns as IDs/variables
4. Set variable metadata (display_name, data_shape, etc.)
5. Create study object from entities
6. Validate with "baseline" and "eda" profiles
7. Export to VDI format (install.json + *.cache files)

### Adding New Datasets

**For phenotype-like modes:**
1. Add file pattern to `nextflow.config` if not already present
2. Create `lib/R/{mode}/{DatasetName}.R` with `wrangle()` function
3. Run: `nextflow run main.nf --mode {mode} --datasetName {DatasetName}`

**For RNA-seq mode:**
- Data files should match existing patterns (countsForEda, eigengenes, sample metadata)
- Multi-dataset studies: add entry to `multiDatasetStudy.json`

### Output Structure
```
results/{study_name}/
├── install.json              # Study definition
├── study.cache              # Serialized study object
├── entitytypegraph.cache    # Entity relationships
└── *.cache                  # Additional cache files
```

### Performance Considerations
- Recent work (PR #5) added benchmarking to identify bottlenecks
- `maxForks = 2` prevents resource exhaustion but slows throughput
- Database validation prevents redundant processing
- Strandedness detection optimized for large count matrices
