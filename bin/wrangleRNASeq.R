#!/usr/local/bin/Rscript

library(tidyverse)

library(glue)

# use study.wrangler v1.0.14 or above for speed
library(study.wrangler)

# Benchmarking utilities
benchmark_results <- list()

benchmark <- function(name, expr) {
  cat("Running:", name, "...\n")
  result <- system.time({
    value <- force(expr)
  })
  benchmark_results[[name]] <<- result
  ## cat("  Elapsed time:", result["elapsed"], "seconds\n")
  ## cat("  User time:", result["user"], "seconds\n")
  ## cat("  System time:", result["system"], "seconds\n\n")
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


read_wgcna_data <- function(filename) {
  col_specs <- cols(
    .default = 'd', 
    "Module" = col_character()
  )

  data <- read_rnaseq_data_default(filename, col_specs) %>%
    mutate(wgcna.ID = sample.ID) %>%
    relocate(wgcna.ID, .after = sample.ID) %>%
    rename(assay.ID = sample.ID)

  
  return(data)
}

read_rnaseq_data <- function(filename) {

  col_specs <- cols(
    .default = col_integer(),
    "...1" = col_character()
  )

  data <- read_rnaseq_data_default(filename, col_specs) %>%
    mutate(assay.ID = sample.ID) %>%
    relocate(assay.ID, .after = sample.ID)


  return(data)
}


read_rnaseq_data_default <- function(filename, columnSpec) {
  cat("Reading RNA-seq data from:", basename(filename), "\n")

  data <- benchmark(paste("read_tsv", basename(filename)), {
    read_tsv(
      filename,
      col_names = TRUE,
      col_types= columnSpec
    )
  })

  headers <- colnames(data)
  headers[1] <- 'sample.ID'
  colnames(data) <- headers

  data <- data %>%
    pivot_longer(-sample.ID, names_to = "column_name", values_to = "value") %>%
    pivot_wider(names_from = sample.ID, values_from = value) %>%
    rename(sample.ID = column_name)

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

  # Compute per-sample sums (across each row - samples are rows after transpose)
  x_sums <- rowSums(x_counts)
  y_sums <- rowSums(y_counts)

  # Count how many samples support each hypothesis
  x_higher_count <- sum(x_sums > y_sums)
  y_higher_count <- sum(x_sums < y_sums)
  total_samples <- length(x_sums)

  # Also check gene-level consistency (genes are columns after transpose)
  # For each gene, sum across all samples
  x_gene_sums <- colSums(x_counts)
  y_gene_sums <- colSums(y_counts)
  x_gene_higher_count <- sum(x_gene_sums > y_gene_sums)
  y_gene_higher_count <- sum(x_gene_sums < y_gene_sums)
  total_genes <- length(x_gene_sums)

  # Write diagnostic report
  report_file <- "strandedness_report.txt"
  cat("=== STRANDEDNESS DETECTION REPORT ===\n", file = report_file)
  cat(sprintf("Files: %s\n\n", paste(names(counts_list), collapse = ", ")), file = report_file, append = TRUE)

  cat("SAMPLE-LEVEL ANALYSIS:\n", file = report_file, append = TRUE)
  cat(sprintf("  Total samples: %d\n", total_samples), file = report_file, append = TRUE)
  cat(sprintf("  Samples where file1 > file2: %d (%.1f%%)\n",
              x_higher_count, 100 * x_higher_count / total_samples), file = report_file, append = TRUE)
  cat(sprintf("  Samples where file2 > file1: %d (%.1f%%)\n",
              y_higher_count, 100 * y_higher_count / total_samples), file = report_file, append = TRUE)

  # Identify sample outliers
  majority_direction <- if (x_higher_count > y_higher_count) "first" else "second"
  sample_outliers <- if (majority_direction == "first") {
    which(x_sums < y_sums)
  } else {
    which(x_sums > y_sums)
  }

  if (length(sample_outliers) > 0) {
    cat("\n  SAMPLE OUTLIERS (opposite from majority):\n", file = report_file, append = TRUE)
    for (idx in sample_outliers) {
      ratio <- if (majority_direction == "first") {
        y_sums[idx] / x_sums[idx]
      } else {
        x_sums[idx] / y_sums[idx]
      }
      cat(sprintf("    Sample: %s\n", x$sample.ID[idx]), file = report_file, append = TRUE)
      cat(sprintf("      File1 total: %d\n", x_sums[idx]), file = report_file, append = TRUE)
      cat(sprintf("      File2 total: %d\n", y_sums[idx]), file = report_file, append = TRUE)
      cat(sprintf("      Ratio: %.3f\n", ratio), file = report_file, append = TRUE)
    }
  }

  cat("\nGENE-LEVEL ANALYSIS:\n", file = report_file, append = TRUE)
  cat(sprintf("  Total genes: %d\n", total_genes), file = report_file, append = TRUE)
  cat(sprintf("  Genes where file1 > file2: %d (%.1f%%)\n",
              x_gene_higher_count, 100 * x_gene_higher_count / total_genes), file = report_file, append = TRUE)
  cat(sprintf("  Genes where file2 > file1: %d (%.1f%%)\n",
              y_gene_higher_count, 100 * y_gene_higher_count / total_genes), file = report_file, append = TRUE)

  # Identify gene outliers (top 20 most discrepant)
  gene_direction <- if (x_gene_higher_count > y_gene_higher_count) "first" else "second"
  gene_outlier_indices <- if (gene_direction == "first") {
    which(x_gene_sums < y_gene_sums)
  } else {
    which(x_gene_sums > y_gene_sums)
  }

  if (length(gene_outlier_indices) > 0) {
    # Calculate fold-change for outliers
    gene_ratios <- if (gene_direction == "first") {
      y_gene_sums[gene_outlier_indices] / pmax(x_gene_sums[gene_outlier_indices], 1)
    } else {
      x_gene_sums[gene_outlier_indices] / pmax(y_gene_sums[gene_outlier_indices], 1)
    }

    # Sort by ratio and take top 20
    top_indices <- head(order(gene_ratios, decreasing = TRUE), 20)
    gene_names <- colnames(x_counts)[gene_outlier_indices[top_indices]]

    cat("\n  TOP 20 GENE OUTLIERS (highest fold-change opposite from majority):\n", file = report_file, append = TRUE)
    for (i in seq_along(top_indices)) {
      orig_idx <- gene_outlier_indices[top_indices[i]]
      cat(sprintf("    %s: file1=%d, file2=%d, ratio=%.2f\n",
                  gene_names[i],
                  x_gene_sums[orig_idx],
                  y_gene_sums[orig_idx],
                  gene_ratios[top_indices[i]]), file = report_file, append = TRUE)
    }
  }

  cat("\n", file = report_file, append = TRUE)

  # Use majority vote with at least 80% agreement threshold
  threshold <- 0.8

  # Collect sample outlier names for summary
  sample_outlier_names <- if (length(sample_outliers) > 0) {
    paste(x$sample.ID[sample_outliers], collapse = ",")
  } else {
    "none"
  }

  if (x_higher_count >= threshold * total_samples) {
    result <- sprintf("Strandedness determined by majority vote: %d/%d samples (%.1f%%) show first > second",
                x_higher_count, total_samples, 100 * x_higher_count / total_samples)
    cat(result, "\n")
    cat(result, "\n", file = report_file, append = TRUE)
    cat(sprintf("Decision: File 1 = Sense, File 2 = Antisense\n"), file = report_file, append = TRUE)

    # Add parseable summary line for Nextflow
    cat(sprintf("\nSUMMARY\t%s\t%d\t%d\t%d\t%.1f\t%s\t%d\t%d\t%.1f\t%s\n",
                orgAbbrev,
                total_samples, x_higher_count, y_higher_count,
                100 * x_higher_count / total_samples,
                "Sense/Antisense",
                total_genes, length(gene_outlier_indices),
                100 * length(gene_outlier_indices) / total_genes,
                sample_outlier_names), file = report_file, append = TRUE)

    cat(sprintf("\nReport written to: %s\n", report_file), file = report_file, append = TRUE)
    names(counts_list) <- c("Sense", "Antisense")
  } else if (y_higher_count >= threshold * total_samples) {
    result <- sprintf("Strandedness determined by majority vote: %d/%d samples (%.1f%%) show second > first",
                y_higher_count, total_samples, 100 * y_higher_count / total_samples)
    cat(result, "\n")
    cat(result, "\n", file = report_file, append = TRUE)
    cat(sprintf("Decision: File 1 = Antisense, File 2 = Sense\n"), file = report_file, append = TRUE)

    # Add parseable summary line for Nextflow
    cat(sprintf("\nSUMMARY\t%s\t%d\t%d\t%d\t%.1f\t%s\t%d\t%d\t%.1f\t%s\n",
                orgAbbrev,
                total_samples, x_higher_count, y_higher_count,
                100 * max(x_higher_count, y_higher_count) / total_samples,
                "Antisense/Sense",
                total_genes, length(gene_outlier_indices),
                100 * length(gene_outlier_indices) / total_genes,
                sample_outlier_names), file = report_file, append = TRUE)

    cat(sprintf("\nReport written to: %s\n", report_file), file = report_file, append = TRUE)
    names(counts_list) <- c("Antisense", "Sense")
  } else {
    # Add parseable summary line even for failures
    cat(sprintf("\nSUMMARY\t%s\t%d\t%d\t%d\t%.1f\t%s\t%d\t%d\t%.1f\t%s\n",
                orgAbbrev,
                total_samples, x_higher_count, y_higher_count,
                100 * max(x_higher_count, y_higher_count) / total_samples,
                "UNDETERMINED",
                total_genes, length(gene_outlier_indices),
                100 * length(gene_outlier_indices) / total_genes,
                sample_outlier_names), file = report_file, append = TRUE)

    cat(sprintf("\nReport written to: %s\n", report_file), file = report_file, append = TRUE)
    stop(sprintf("Strandedness could not be determined for %s: only %d/%d samples (%.1f%%) agree",
                 paste(basename(names(counts_list)), collapse=", "),
                 max(x_higher_count, y_higher_count),
                 total_samples,
                 100 * max(x_higher_count, y_higher_count) / total_samples))
  }

  prefix = paste0(orgAbbrev, "_")

  names(counts_list) = paste0(prefix, names(counts_list))

  return(counts_list);
}

counts_to_entity <- function(tbl, name, orgAbbrev) {

  # strip off the prefix
  name = sub(glue("{orgAbbrev}_"), "", name);

  
  assays <- benchmark(paste("entity_from_tibble", name), {
    entity_from_tibble(
      tbl,
      name = glue("{orgAbbrev} {name} counts"),
      display_name = glue("{orgAbbrev} {name} htseq counts"),
      display_name_plural = glue("{orgAbbrev} {name} htseq counts"),
      skip_type_convert = TRUE
    ) %>%
      set_parents('sample', 'sample.ID') %>%
  ##    set_variable_display_names_from_provider_labels() %>%
      set_variable_metadata('sample.ID', display_name = "Sample ID", hidden=list('variableTree')) %>%
      set_variable_metadata('assay.ID', display_name = "HTSeq Count", hidden=list('variableTree'))
  })

  assays <- benchmark(paste("set_stable_ids", name), {
    # Directly mutate the variables tibble instead of using reduce()
    gene_idx <- assays@variables$data_type != 'id'
    assays@variables$stable_id[gene_idx] <- assays@variables$variable[gene_idx]
    assays@variables$display_name[gene_idx] <- assays@variables$variable[gene_idx]
    assays
  })


  
  assays <- benchmark(paste("create_variable_metadata", name), {
    assays %>%
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
  })

  assays
}


wgcna_to_entity <- function(tbl, name, orgAbbrev) {
  
  assays <- benchmark(paste("wgcna_entity_from_tibble", name), {
    entity_from_tibble(
      tbl,
      name = glue("{orgAbbrev} eigengene"),
      display_name = glue("{orgAbbrev} Eigengene (wgcna)"),
      display_name_plural = glue("{orgAbbrev} Eigengenes (wgcna)"),
      skip_type_convert = TRUE
      #TO DO, description = ???
    ) %>%
      set_parents(glue("{orgAbbrev} Sense counts"), 'assay.ID') %>%
      set_variable_display_names_from_provider_labels() %>%
      set_variable_metadata('assay.ID', display_name = "Assay ID", hidden=list('variableTree')) %>%
      set_variable_metadata('wgcna.ID', display_name = "WGCNA ID", hidden=list('variableTree'))
  })

  assays <- benchmark(paste("wgcna_variable_metadata", name), {
    assays %>%
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
  })

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

  file_discovery <- benchmark("file_globbing", {
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

    list(
      sample_tsv_file = sample_tsv_file,
      count_files_by_prefix = count_files_by_prefix,
      wgcna_files_by_prefix = wgcna_files_by_prefix
    )
  })

  # call countsData() on each group, passing (files, prefix)
  all_counts_entities <- benchmark("process_counts_data", {
    file_discovery$count_files_by_prefix %>%
      imap(~ countsData(.x, .y)) %>%
      flatten()
  })

  # call wgcnaData() on each group, passing (files, prefix)
  all_wgcna_entities <- benchmark("process_wgcna_data", {
    file_discovery$wgcna_files_by_prefix %>%
      imap(~ wgcnaData(.x, .y)) %>%
      flatten()
  })

  # the entity will have the name 'sample'
  samples <- benchmark("process_sample_metadata", {
    entity_from_stf(file_discovery$sample_tsv_file) %>%
      set_variable_metadata('sample.ID', display_name = "Sample ID", hidden=list('variableTree'))
  })


  #return(list(samples = samples, all_counts_entities = all_counts_entities, all_wgcna_entities = all_wgcna_entities))
  
  study <- benchmark("create_study", {
    study_from_entities(c(samples, all_counts_entities, all_wgcna_entities), name = "RNA-Seq study")
  })

  print_benchmark_summary()

  return(study);
}


countsData <- function(counts_filenames, orgAbbrev) {

  counts_data = benchmark(paste("read_counts_files", orgAbbrev), {
    counts_filenames %>%
      set_names() %>%
      map(read_rnaseq_data)
  })
  
  counts_data <- benchmark(paste("rename_strandedness", orgAbbrev), {
    rename_counts_by_strandedness(counts_data, orgAbbrev)
  })

  counts_data <- benchmark(paste("convert_to_entities", orgAbbrev), {
    counts_data %>%
      imap(counts_to_entity, orgAbbrev)
  })

  return(counts_data)
}

wgcnaData <- function(counts_filenames, orgAbbrev) {

  counts_data <- benchmark(paste("read_wgcna_files", orgAbbrev), {
    counts_filenames %>%
      set_names() %>%
      map(read_wgcna_data)
  })

  counts_data <- benchmark(paste("convert_wgcna_to_entities", orgAbbrev), {
    counts_data %>%
      imap(wgcna_to_entity, orgAbbrev)
  })

  return(counts_data)
}

study = wrangle()

if(!validate(study, profiles=c("baseline", "eda"))) {
  stop("Stopping....Study is not valid");
}

 export_to_vdi(study, getwd());
