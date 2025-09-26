#!/usr/local/bin/Rscript

library(tidyverse)

# use study.wrangler v1.0.14 or above for speed
library(study.wrangler)

# Benchmarking utilities
benchmark_results <- list()

benchmark <- function(name, expr) {
  cat("Benchmarking:", name, "...\n")
  result <- system.time({
    value <- force(expr)
  })
  benchmark_results[[name]] <<- result
  cat("  Elapsed time:", result["elapsed"], "seconds\n")
  cat("  User time:", result["user"], "seconds\n")
  cat("  System time:", result["sys"], "seconds\n\n")
  return(value)
}

print_benchmark_summary <- function() {
  cat("\n=== BENCHMARK SUMMARY ===\n")
  total_elapsed <- 0
  for (name in names(benchmark_results)) {
    elapsed <- benchmark_results[[name]]["elapsed"]
    total_elapsed <- total_elapsed + elapsed
    cat(sprintf("%-30s: %8.3f seconds\n", name, elapsed))
  }
  cat(sprintf("%-30s: %8.3f seconds\n", "TOTAL", total_elapsed))
  cat("========================\n\n")
}


read_rnaseq_data <- function(filename) {
  cat("Reading RNA-seq data from:", basename(filename), "\n")

  # read in as all-character
  data <- benchmark(paste("read_tsv", basename(filename)), {
    read_tsv(
      filename,
      col_names = FALSE,
      col_types = cols(.default = "c")
    )
  })

  # and transpose
  data <- benchmark(paste("transpose", basename(filename)), {
    data %>%
      as.matrix() %>%
      t()
  })

  # make colnames manually to avoid pairwise compute
  colnames(data) <- sprintf("V%i", seq_len(ncol(data)))

  # wrap it in a tibble with NO name repair overhead
  data <- benchmark(paste("as_tibble", basename(filename)), {
    as_tibble(data, .name_repair = "minimal")
  })

  # Extract the header row
  headers <- data %>% slice_head(n = 1) %>% unlist(use.names = FALSE)
  headers[1] <- 'sample.ID'
  # Drop header row
  data <- data %>% slice_tail(n = nrow(data) - 1)
  # Name the columns and we're back to a proper tibble
  colnames(data) <- headers

  # convert the counts columns to integer
  # and make a trivial assay ID column
  data <- benchmark(paste("type_conversion", basename(filename)), {
    data %>%
      mutate(
        across(
          -sample.ID,
          as.integer
        ),
        assay.ID = sample.ID
      ) %>%
      relocate(assay.ID, .after = sample.ID) %>%
      select(-starts_with('__'))
  })

  return(data)
}

#
# figures out which countsForEda file is sense/antisense
#
# returns a named ("sense", "antisense" or "unstranded") list of tibbles
#
rename_counts_by_strandedness <- function(counts_list, orgAbbrev) {
  if (length(counts_list) == 1) {
    # Only one file → label it unstranded
    names(counts_list) <- paste0("unstranded_", orgAbbrev)
    return(counts_list)
  }
  
  stopifnot(length(counts_list) == 2)
  
  # Extract tibbles
  x <- counts_list[[1]]
  y <- counts_list[[2]]
  
  # Remove metadata columns (assumes first 2 cols are not genes)
  x_counts <- x %>% select(-1, -2)
  y_counts <- y %>% select(-1, -2)
  
  # Compute per-sample sums (across each row)
  x_sums <- rowSums(x_counts)
  y_sums <- rowSums(y_counts)

  if (all(x_sums > y_sums)) {
    names(counts_list) <- c("Sense", "Antisense")
  } else if (all(x_sums < y_sums)) {
    names(counts_list) <- c("Antisense", "Sense")
  } else {
    stop("Strandedness could not be consistently determined for ", names(counts_list))
  }

  prefix = paste0(orgAbbrev, "_")

  names(counts_list) = paste0(prefix, names(counts_list)) 
  
  return(counts_list);
}

counts_to_entity <- function(tbl, name, orgAbbrev) {

  # strip off the prefix
  name = sub(glue("{orgAbbrev}_"), "", name); 
  
  assays <- entity_from_tibble(
    tbl,
    name = glue("{orgAbbrev}_{name}_counts"),
    display_name = glue("{orgAbbrev} {name} htseq counts"),
    display_name_plural = glue("{orgAbbrev} {name} htseq counts"),
    skip_type_convert = TRUE
  ) %>%
    set_parents('sample', 'sample.ID') %>% 
    set_variable_display_names_from_provider_labels() %>%
    set_variable_metadata('sample.ID', display_name = "Sample ID", hidden=list('variableTree')) %>%
    set_variable_metadata('assay.ID', display_name = "HTSeq Count", hidden=list('variableTree'))
    
  
  assays <- assays %>%
    create_variable_category(
      category_name = 'gene_counts',
      children = assays %>% get_variable_metadata() %>% pull(variable),
      display_name = 'Gene counts',
      definition = 'Counts per gene from RNA-Seq'
    ) %>%
    create_variable_collection(
      'gene_counts',
      member = 'gene',
      member_plural = 'genes',
      is_proportion = FALSE,
      is_compositional = TRUE
    )

  assays
}


wgcna_to_entity <- function(tbl, name, orgAbbrev) {

  assays <- entity_from_tibble(
    tbl,
    name = glue("{orgAbbrev}_eigengene"),
    display_name = glue("{orgAbbrev} Eigengene (wgcna)"),
    display_name_plural = glue("{orgAbbrev} Eigengenes (wgcna)"),
    skip_type_convert = TRUE
    #TO DO, description = ???
  ) %>%
    set_parents(glue("{orgAbbrev}_{name}_counts"), 'assay.ID') %>%
    set_variable_display_names_from_provider_labels() %>%
    set_variable_metadata('sample.ID', display_name = "Sample ID", hidden=list('variableTree')) %>%
    set_variable_metadata('assay.ID', display_name = "Assay ID", hidden=list('variableTree'))

  assays <- assays %>%
    create_variable_category(
      category_name = 'wgcna',
      children = assays %>% get_variable_metadata() %>% pull(variable),
      display_name = 'Eigengenes',
      definition = 'Eigengene from RNA-Seq'
    ) %>%
    create_variable_collection(
      'wgcna',
      member = 'eigengene',
      member_plural = 'eigengenes',
      stable_id =  "EUPATH_0005051",
      is_proportion =  FALSE,
      is_compositional = FALSE
    )

  assays
}


# returns a named list
# of filename character vectors
group_files_by_prefix <- function(filenames, pattern) {
  prefixes <- sub(glue("{pattern}.*"), "", filenames) # extract text up to first “_”
  split(filenames, prefixes)
}

#
# the main event... wrangle()
#
wrangle <- function() {

  countsPattern = "_countsForEda_"
  wgcnaPattern = '_merged-0.25-eigengenes_'
  
  # find the sample STF file
  sample_tsv_path <- file.path(
    '*entity-sample.tsv'
  )

  sample_tsv_file = Sys.glob(sample_tsv_path)
  
  
  # find the countsForEda_*.txt files
  counts_file_glob <- file.path(
    glue("*{countsPattern}*")
  )

  counts_filenames <- Sys.glob(counts_file_glob)
  count_files_by_prefix <- group_files_by_prefix(counts_filenames, countsPattern)

  # find the countsForEda_*.txt files
  wgcna_file_glob <- file.path(
    glue("*{wgcnaPattern}*")
  
  )

  wgcna_filenames <- Sys.glob(wgcna_file_glob)
  wgcna_files_by_prefix <- group_files_by_prefix(wgcna_filenames, wgcnaPattern)

  # call countsData() on each group, passing (files, prefix)  
  all_counts_entities <- count_files_by_prefix %>%
    imap(~ countsData(.x, .y)) %>%
    flatten()


  # call countsData() on each group, passing (files, prefix)  
  all_wgcna_entities <- wgcna_files_by_prefix %>%
    imap(~ wgcnaData(.x, .y)) %>%
    flatten()

  
  # the entity will have the name 'sample'
  samples <- entity_from_stf(sample_tsv_file) %>%
    set_variable_metadata('sample.ID', display_name = "Sample ID", hidden=list('variableTree'))
  
  study <- study_from_entities(c(samples, all_counts_entities, all_wgcna_entities), name = "RNA-Seq study")
  
  study
}


countsData <- function(counts_filenames, orgAbbrev) {

  counts_data <- counts_filenames %>%
    set_names() %>%
    map(read_rnaseq_data)  %>%
    rename_counts_by_strandedness(., orgAbbrev) %>%
    imap(counts_to_entity, orgAbbrev)

  return(counts_data)
}

wgcnaData <- function(counts_filenames, orgAbbrev) {

  counts_data <- counts_filenames %>%
    set_names() %>%
    map(read_rnaseq_data)  %>%
    imap(wgcna_to_entity, orgAbbrev)
    #imap(purrr::partial(wgcna_to_entity, orgAbbrev=orgAbbrev))

  return(counts_data)
}



## study = wrangle()
## if(!validate(study, profiles=c("baseline", "eda")) {
##   stop("Stopping....Study is not valid");
## }

## export_to_vdi(study, getwd());
