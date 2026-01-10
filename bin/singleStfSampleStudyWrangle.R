#!/usr/local/bin/Rscript

library(study.wrangler)

args <- commandArgs(trailingOnly = TRUE)

my_r_lib <- Sys.getenv("MY_R_LIB")
# Check if at least one argument provided
if (length(args) < 1) {
  stop("usage:  wrangleAntibodyArray.R dataset")
}

datasetName = args[1];

sampleEntity = entity_from_stf("entity-sample.tsv", "entity-sample.yaml")

study = study(name=datasetName, sampleEntity)

if(!validate(study, profiles=c("baseline", "eda"))) {
  stop("Stopping....Study is not valid");
}

export_to_vdi(study, getwd());


