library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  my_r_lib <- Sys.getenv("MY_R_LIB")
  shared_env <- new.env()
  source(paste0(my_r_lib, "/phenotype/PHI-base_curated_phenotype_shared.R"), local = shared_env)

  shared_env$wrangle_phi_base("phenotypes_TriTrypDB.txt")
}
