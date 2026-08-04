#!/usr/local/bin/Rscript

## Build the final per-REFERENCE-ORGANISM dnaseq study.
##
## Two inputs, both staged into the working directory by nextflow:
##
##   entity-sample.{tsv,yaml}  the harmonized sample characteristics for every
##                             dnaseq experiment aligned to this reference
##                             organism (union of their variables)
##   *_alignment_stats.tsv     one file per sample, emitted by the dnaseq
##                             workflow that produced the alignments
##
## The alignment stats land on the sample entity rather than on a child assay
## entity: this study is scoped to a single reference organism, so a sample has
## exactly one alignment against it. One row per sample, no fan-out.
##
## usage:  combinedDnaseqStudyWrangle.R <studyName> <organismAbbrev>

suppressPackageStartupMessages({
  library(tidyverse)
  library(study.wrangler)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: combinedDnaseqStudyWrangle.R <studyName> <organismAbbrev>")

study_name <- args[1]
organism_abbrev <- args[2]

## The alignment stats columns, in the order they should appear in the UI.
## `sample` is the join key and is dropped after joining.
##
## NOTE on the two *_mapped_pct columns: despite the names, the workflow writes
## fractions, not percentages, and raw_mapped_pct exceeds 1.0 for paired data
## (mapped mates over raw read pairs). The display names say ratio so nobody
## reads 0.09 as 9%.
ALIGNMENT_VARS <- list(
  raw_fastq_reads = list(
    display_name = "Raw reads",
    definition = "Number of reads in the submitted fastq file(s), before trimming."
  ),
  trimmed_total_sequences = list(
    display_name = "Trimmed reads",
    definition = "Number of sequences surviving adapter and quality trimming, as counted by samtools stats."
  ),
  reads_mapped = list(
    display_name = "Reads mapped",
    definition = "Number of trimmed reads mapped to the reference genome."
  ),
  trimmed_mapped_pct = list(
    display_name = "Mapped / trimmed reads (ratio)",
    definition = "Reads mapped divided by trimmed reads. A fraction between 0 and 1, not a percentage."
  ),
  raw_mapped_pct = list(
    display_name = "Mapped / raw reads (ratio)",
    definition = "Reads mapped divided by raw fastq reads. A fraction, which can exceed 1 for paired data because both mates of a pair are counted as mapped reads while the raw count is of read pairs."
  ),
  average_read_length = list(
    display_name = "Average read length",
    definition = "Mean length in base pairs of the mapped reads."
  ),
  mean_coverage = list(
    display_name = "Mean coverage",
    definition = "Mean depth of coverage across the reference genome."
  )
)

ALIGNMENT_CATEGORY <- "alignment_statistics"
ALIGNMENT_CATEGORY_DISPLAY <- "Alignment statistics"

## ---------------------------------------------------------------- read inputs

sample_entity <- entity_from_stf("entity-sample.tsv", "entity-sample.yaml")

stats_files <- list.files(".", pattern = "_alignment_stats\\.tsv$")
if (length(stats_files) == 0) stop("no *_alignment_stats.tsv files were staged")

message("reading ", length(stats_files), " alignment stats files")

stats <- stats_files %>%
  map(~ read_tsv(.x, col_types = cols(sample = col_character(), .default = col_double()))) %>%
  list_rbind()

## --------------------------------------------------------------- sanity checks
## sample_id in the stf is byte-identical to the `sample` column of these files
## by construction -- that is the whole reason sample_id was not "improved"
## during harmonization. If that ever stops being true we want a hard failure
## here, not a study full of empty alignment columns.

dupes <- stats$sample[duplicated(stats$sample)]
if (length(dupes) > 0)
  stop("duplicate samples in alignment stats: ", paste(unique(dupes), collapse = ", "))

stf_ids <- sample_entity %>% get_data() %>% pull(sample_id)

missing_stats <- setdiff(stf_ids, stats$sample)
missing_meta <- setdiff(stats$sample, stf_ids)

if (length(missing_meta) > 0)
  stop(length(missing_meta), " samples have alignment stats but no row in the stf: ",
       paste(head(missing_meta, 10), collapse = ", "))

if (length(missing_stats) > 0)
  stop(length(missing_stats), " samples are in the stf but have no alignment stats: ",
       paste(head(missing_stats, 10), collapse = ", "))

message("joined ", nrow(stats), " samples of alignment stats onto ", length(stf_ids),
        " stf samples for ", organism_abbrev)

## ----------------------------------------------------------------- merge them

alignment_cols <- names(ALIGNMENT_VARS)

stats <- stats %>%
  rename(sample_id = sample) %>%
  select(all_of(c("sample_id", alignment_cols)))

## modify_data's argument is the right hand side of a pipeline whose left hand
## side is the entity's data, so left_join receives it as `x` implicitly.
sample_entity <- sample_entity %>%
  modify_data(left_join(stats, by = "sample_id"))

## New columns arrive with no metadata at all: sync_variable_metadata() creates
## the rows, redetect_columns_as_variables() then guarantees they are treated as
## variables rather than as extra id columns (a numeric column whose values
## happen to be unique per sample would otherwise be detected as an id and make
## the entity invalid). The entity's id column is re-asserted last because the
## data changed underneath it.
sample_entity <- sample_entity %>%
  sync_variable_metadata() %>%
  redetect_columns_as_variables(alignment_cols) %>%
  redetect_column_as_id("sample_id") %>%
  set_variable_metadata("sample_id", entity_name = "sample")

for (v in alignment_cols) {
  sample_entity <- sample_entity %>%
    set_variable_metadata(v,
                          display_name = ALIGNMENT_VARS[[v]]$display_name,
                          definition = ALIGNMENT_VARS[[v]]$definition)
}

sample_entity <- sample_entity %>%
  create_variable_category(ALIGNMENT_CATEGORY,
                           display_name = ALIGNMENT_CATEGORY_DISPLAY,
                           children = alignment_cols)

## ------------------------------------------------------------ validate + emit

study <- study(name = study_name, sample_entity)

if (!validate(study, profiles = c("baseline", "eda"))) {
  stop("Stopping....Study is not valid")
}

export_to_vdi(study, getwd())
