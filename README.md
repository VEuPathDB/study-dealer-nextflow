# study-dealer-nextflow

A Nextflow pipeline that processes biological study data for the VEuPathDB project. The pipeline "deals studies to the wrangler" by taking various genomics datasets and processing them through standardized wrangling workflows.

## Overview

This pipeline automates the processing of diverse genomics datasets into standardized formats for the VEuPathDB ecosystem. It supports multiple data types including RNA sequencing, phenotype data, antibody arrays, cellular localization, and RFLP datasets.

### Key Features

- **Mode-based architecture**: Process different data types (RNA-seq, phenotype, etc.) through specialized workflows
- **Multi-dataset study support**: Automatically groups related datasets into unified studies
- **Flexible configuration**: JSON-based study mapping and glob pattern file discovery
- **Containerized execution**: Docker-based processing with reproducible environments
- **Validation & quality control**: Built-in validation against baseline and EDA profiles
- **Smart processing**: Avoids reprocessing already-wrangled studies

## Prerequisites

- [Nextflow](https://www.nextflow.io/) (version 20.01.0 or later)
- Docker (recommended) or Singularity/Apptainer
- Access to VEuPathDB data directories
- GUS home directory with configuration file at `{gusHomeDir}/config/gus.config`

## Quick Start

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
| `multiDatasetStudies` | JSON file for study mapping | `$baseDir/data/rnaseq_sample_reannotation/multiDatasetStudy.json` |
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

For RNA-seq studies spanning multiple organisms, define study groupings in `data/rnaseq_sample_reannotation/multiDatasetStudy.json`:

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
results/{study_name}/
├── install.json              # Study definition (VDI format)
├── study.cache              # Serialized study object
├── entitytypegraph.cache    # Entity relationships
└── *.cache                  # Additional cache files
```

## Adding New Datasets

### For Phenotype-like Modes

1. **Create wrangling script**: `lib/R/{mode}/{DatasetName}.R`

```R
library(study.wrangler)

wrangle <- function() {
  # Read data
  entity <- entity_from_file("data.txt")

  # Configure entity
  entity <- entity %>%
    set_entity_metadata(name = "samples", display_name = "Samples") %>%
    redetect_column_as_id("sample_id") %>%
    set_variable_metadata("phenotype", display_name = "Phenotype")

  # Create and validate study
  study <- study("MyStudy", entity)
  if(!validate(study, profiles=c("baseline", "eda"))) {
    stop("Study invalid")
  }

  # Export
  export_to_vdi(study, getwd())
  return(study)
}
```

2. **Run pipeline**:
```bash
nextflow run main.nf --mode {mode} --datasetName {DatasetName}
```

### For RNA-seq Mode

Ensure data files match expected patterns:
- Count files: `countsForEda_*.txt`
- Eigengenes: `*eigengenes*.txt`
- Sample metadata: `entity-sample.{tsv,yaml}`

For multi-organism studies, add entry to `multiDatasetStudy.json`.

## Docker Configuration

The pipeline uses Docker by default with the following settings:

- **User mapping**: Runs containers as current user for proper file permissions
- **Network**: Host network access enabled
- **Environment variables**:
  - `MY_R_LIB`: Points to custom R scripts (`lib/R/`)
  - `GUS_HOME`: Database configuration location

### Containers Used

- `veupathdb/study-wrangler:1.0.27`: Main processing with R and study.wrangler library
- `veupathdb/alpine_bash:latest`: Utility operations
- `veupathdb/vdi-plugin-wrangler`: Database name extraction

### Using Singularity

To use Singularity instead of Docker, edit `nextflow.config`:

```groovy
includeConfig 'conf/singularity.config'  // Instead of docker.config
```

## Troubleshooting

### Common Issues

**Pipeline skips my study:**
- Check that database name exists in external databases but NOT in EDA
- Review logs for filtering messages: `log.info "Skipping study..."`

**Validation errors:**
- Ensure `wrangle()` function returns valid study object
- Check that all required metadata is set
- Verify variable data shapes match actual data

**File not found errors:**
- Verify data directory structure matches expected pattern
- Check file patterns in `nextflow.config`
- Ensure `workflowDataDir` parameter points to correct location

**GUS config errors:**
- Verify GUS config file exists at `{gusHomeDir}/config/gus.config`
- Check `--gusHomeDir` parameter points to correct location

### Debug Mode

Run with Nextflow's debug options:

```bash
nextflow run main.nf -with-trace -with-report -with-dag flowchart.html
```

## Performance

- **Concurrent processes**: Limited to 2 (`maxForks = 2`) to prevent resource exhaustion
- **Benchmarking**: RNA-seq mode includes timing utilities for performance analysis
- **Database validation**: Recent addition prevents redundant reprocessing

## Architecture

### Components

- **workflows/**: High-level orchestration (RNA-seq multi-dataset handling)
- **subworkflows/**: Reusable components (single study processing)
- **modules/**: Atomic processes (wrangling, filtering, validation)
- **bin/**: Executable R scripts (wrangleRNASeq.R, singleStudyWrangle.R)
- **lib/R/**: Custom dataset-specific wrangling scripts

### Processing Flow (RNA-seq)

1. Discover files via glob patterns
2. Extract metadata from directory paths
3. Mix counts, WGCNA, and sample metadata
4. Add organism prefixes and filter system rows
5. Group datasets by study
6. Validate against database names
7. Wrangle each study group
8. Validate output
9. Export to VDI format

## Contributing

When adding new functionality:

1. Follow the mode-based dispatcher pattern in `main.nf`
2. Use study.wrangler library for consistency
3. Implement validation with "baseline" and "eda" profiles
4. Add file patterns to `nextflow.config`
5. Document custom wrangling scripts

## License

VEuPathDB Project

## Contact

For issues or questions, please contact the VEuPathDB development team.
