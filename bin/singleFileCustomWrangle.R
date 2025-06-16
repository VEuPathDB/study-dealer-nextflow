#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(study.wrangler)


# Check if at least two arguments are provided
if (length(args) < 1) {
  stop("usage:  singleFileCustomWrangle.R custom.R")
}

wrangleScript <- args[1]

source(wrangleScript)

study = wrangle()


# TODO:  will this stop if validation fails?
validate(study)

export_to_vdi(study, getwd());
