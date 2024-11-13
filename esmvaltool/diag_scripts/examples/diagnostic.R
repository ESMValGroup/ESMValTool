# #############################################################################
# diagnostic.R
# Authors:       E. Arnone (ISAC-CNR, Italy)
# #############################################################################
# Description
# Example of ESMValTool diagnostic written in R
#
# Required
#
# Optional
#
# Caveats
#
# Modification history
#    20180620-A_arnone_e: written for v2.0
#
# ############################################################################
library(tools)
library(yaml)
# get path to script and source subroutines (if needed)
diag_scripts_dir <- Sys.getenv("diag_scripts")
# source paste0(diag_scripts_dir,"/subroutine.r")
print(file.path("source ", diag_scripts_dir, "subroutine.r"))
# read settings and metadata files (assuming one variable only)
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}
metadata <- yaml::read_yaml(settings$input_files)
# get name of climofileis for first variable and list
# associated to first climofile
climofiles <- names(metadata)
climolist <- get(climofiles[1], metadata)
# get diagnostic name from metadata file
diag_base <- climolist$diagnostic
print(paste0(diag_base, ": starting routine"))
# create work and plot directories if they do not exist
print(paste0(diag_base, ": creating work and plot directories"))
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(plot_dir, recursive = T, showWarnings = F)
# extract metadata
models_name <- unname(sapply(metadata, "[[", "dataset"))
reference_model <-
  unname(sapply(metadata, "[[", "reference_dataset"))[1]
models_start_year <- unname(sapply(metadata, "[[", "start_year"))
models_end_year <- unname(sapply(metadata, "[[", "end_year"))
models_experiment <- unname(sapply(metadata, "[[", "exp"))
models_ensemble <- unname(sapply(metadata, "[[", "ensemble"))
## Loop through input models
for (model_idx in c(1:(length(models_name)))) {
  # Setup parameters and path
  model <- models_name[model_idx]
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]
  infile <- climofiles[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  print(paste0(diag_base, ": working on file ", infile))
  print(paste0(diag_base, ": calling diagnostic with following parameters"))
  print(paste(model, model_exp, model_ens, year1, year2))
  ## Call actual diagnostic
  print(paste0(diag_base, ": I am your R diagnostic"))
}
