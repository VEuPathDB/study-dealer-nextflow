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
##
## The alignment stats are the authority on which samples exist. They are what
## the workflow actually aligned, and sample_id in the stf is byte-identical to
## their `sample` column by construction -- that is the whole reason sample_id
## was not "improved" during harmonization.
##
## So the two directions are NOT symmetric:
##
##   sample aligned but not in the stf  -> hard failure. The workflow produced
##       data we have no metadata for, which is a real gap in the metadata.
##   sample in the stf but not aligned  -> report and drop it. A test database
##       legitimately holds a subset of the aligned data, and we only want to
##       process the metadata that corresponds to it.

dupes <- stats$sample[duplicated(stats$sample)]
if (length(dupes) > 0)
  stop("duplicate samples in alignment stats: ", paste(unique(dupes), collapse = ", "))

stf_ids <- sample_entity %>% get_data() %>% pull(sample_id)

missing_meta <- setdiff(stats$sample, stf_ids)
if (length(missing_meta) > 0)
  stop(length(missing_meta), " samples have alignment stats but no row in the stf: ",
       paste(head(missing_meta, 10), collapse = ", "))

not_aligned <- setdiff(stf_ids, stats$sample)
if (length(not_aligned) > 0) {
  message("NOTE: dropping ", length(not_aligned), " of ", length(stf_ids),
          " stf samples with no alignment stats (subset of aligned data): ",
          paste(head(not_aligned, 10), collapse = ", "),
          if (length(not_aligned) > 10) ", ..." else "")

  sample_entity <- sample_entity %>%
    modify_data(filter(!sample_id %in% not_aligned))

  ## Dropping samples can empty a variable out completely -- one recorded by
  ## only the experiments that are absent from this subset. An all-empty
  ## variable is noise in the UI, and a category left with no children is not
  ## valid, so both go.
  empty_vars <- sample_entity %>%
    get_data() %>%
    select(-sample_id) %>%
    select(where(~ all(is.na(.)))) %>%
    names()

  if (length(empty_vars) > 0) {
    message("NOTE: dropping ", length(empty_vars),
            " variables left with no values by the sample subset: ",
            paste(empty_vars, collapse = ", "))

    sample_entity <- sample_entity %>%
      modify_data(select(-all_of(empty_vars))) %>%
      sync_variable_metadata()
  }

  ## Looped because pruning a category can orphan its parent.
  repeat {
    meta <- sample_entity %>% get_variable_and_category_metadata()

    childless <- meta %>%
      filter(data_type == "category", !variable %in% meta$parent_variable) %>%
      pull(variable)

    if (length(childless) == 0) break

    message("NOTE: dropping empty variable categories: ", paste(childless, collapse = ", "))
    for (category in childless) {
      sample_entity <- sample_entity %>% delete_variable_category(category)
    }
  }
}

message("joined ", nrow(stats), " samples of alignment stats onto ",
        length(stf_ids) - length(not_aligned), " stf samples for ", organism_abbrev)

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
