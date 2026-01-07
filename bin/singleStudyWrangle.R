#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


my_r_lib <- Sys.getenv("MY_R_LIB")
# Check if at least two arguments are provided
if (length(args) < 2) {
  stop("usage:  wrangleAntibodyArray.R dataset type")
}

datasetName = args[1];
datasetType = args[2];

datasetScript <- paste0(my_r_lib, "/", datasetType, "/", datasetName, ".R");
dataset_env <- new.env()

source(datasetScript, local=dataset_env)

# Always assign the studyName to be the dataset name for uniqueness
study <- dataset_env$wrangle() %>%
  dataset_env$set_study_name(datasetName)



if(!validate(study, profiles=c("baseline", "eda"))) {
  stop("Stopping....Study is not valid");
}

export_to_vdi(study, getwd());


