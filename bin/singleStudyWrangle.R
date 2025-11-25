#!/usr/local/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)


my_r_lib <- Sys.getenv("MY_R_LIB")
# Check if at least two arguments are provided
if (length(args) < 2) {
  stop("usage:  wrangleAntibodyArray.R dataset type")
}

datasetScript <- paste0(my_r_lib, "/", args[2], "/", args[1], ".R");
dataset_env <- new.env()

source(datasetScript, local=dataset_env)

study = dataset_env$wrangle()

if(!validate(study, profiles=c("baseline", "eda"))) {
  stop("Stopping....Study is not valid");
}

export_to_vdi(study, getwd());


