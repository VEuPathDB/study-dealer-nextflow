#!/usr/local/bin/Rscript

library(tidyverse)

# use study.wrangler v1.0.14 or above for speed
library(study.wrangler)


read_counts_data <- function(filename) {

  # read in as all-character
  data <- read_tsv(
    filename,
    col_names = FALSE,
    col_types = cols(.default = "c")
  ) %>%
    # and transpose
    as.matrix() %>%
    t()

  # make colnames manually to avoid pairwise compute
  colnames(data) <- sprintf("V%i", seq_len(ncol(data)))

  # wrap it in a tibble with NO name repair overhead
  data <- as_tibble(data, .name_repair = "minimal")
  
  # Extract the header row
  headers <- data %>% slice_head(n = 1) %>% unlist(use.names = FALSE)
  headers[1] <- 'sample.ID'
  # Drop header row
  data <- data %>% slice_tail(n = nrow(data) - 1)
  # Name the columns and we're back to a proper tibble
  colnames(data) <- headers

  # convert the counts columns to integer
  # and make a trivial assay ID column
  data <- data %>%
    mutate(
      across(
        -sample.ID,
        as.integer
      ),
      assay.ID = sample.ID
    ) %>%
    relocate(assay.ID, .after = sample.ID)
  
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
    names(counts_list) <- c("sense", "antisense")
  } else if (all(x_sums < y_sums)) {
    names(counts_list) <- c("antisense", "sense")
  } else {
    stop("Strandedness could not be consistently determined for ", names(counts_list))
  }

  suffix = paste0("_", orgAbbrev)

  names(counts_list) = paste0(names(counts_list), suffix) 
  
  return(counts_list);
}

counts_to_entity <- function(tbl, name) {
  assays <- entity_from_tibble(
    tbl,
    name = paste0(name, "_assay"),
    display_name = paste(str_to_title(name), 'assay'),
    display_name_plural = paste(str_to_title(name), 'assays'),
    skip_type_convert = TRUE
    #TO DO, description = ???
  ) %>%
    set_parents('sample', 'sample.ID')
  
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
      member_plural = 'genes'
      # TO DO: normalization_method? is_compositional?
    )
  
  assays
}

#
# the main event... wrangle()
#
wrangle <- function() {
  
  # find the sample STF file
  sample_tsv_path <- file.path(
    '*entity-sample.tsv'
  )

  sample_tsv_file = Sys.glob(sample_tsv_path)
  
  
  # find the countsForEda_*.txt files
  counts_file_glob <- file.path(
    '*_countsForEda*.txt'
  )

  counts_filenames <- Sys.glob(counts_file_glob)

  prefixes = strsplit(counts_filenames, "_")
  uniquePrefixes = unique(sapply(prefixes, function(x) x[1]))


  all_counts_data = c()

  for(i in 1:length(uniquePrefixes)) {
    prefix = uniquePrefixes[i];

    filtered_counts_filenames = counts_filenames[startsWith(counts_filenames, prefix)]

    filteredCountsData = countsData(filtered_counts_filenames, prefix)

    all_counts_data[i] = filteredCountsData 
  }

  
  # the entity will have the name 'sample'
  samples <- entity_from_stf(sample_tsv_file)
  
  study <- study_from_entities(c(samples, all_counts_data), name = "RNA-Seq study")
  
  study
}


countsData <- function(counts_filenames, orgAbbrev) {

  counts_data <- counts_filenames %>% 
    set_names() %>%
    map(read_counts_data)  %>%
    rename_counts_by_strandedness(., orgAbbrev) %>%
    imap(counts_to_entity)

  return(counts_data)
}


study = wrangle()

# TODO:  will this stop if validation fails?
validate(study)

export_to_vdi(study, getwd());
