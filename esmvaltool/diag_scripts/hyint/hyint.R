# #############################################################################
# hyint.R
# Authors:       E. Arnone (ISAC-CNR, Italy)
#                J. von Hardenberg (ISAC-CNR, Italy)
# #############################################################################
# Description
# HyInt is a tool for calculation of the HY-INT index (Giorgi et al. 2011)
# and additional hydroclimatic indices (Giorgi et al. 2014)
# which allow an estimate of the overall behaviour of the hydroclimatic cycle.
# The tool calculates also timeseries and trends over selected regions and
# produces a variety of types of plots including maps and timeseries. The timeseries/trend
# and plotting modules handle also ETCCDI indices data calculated with the climdex library through
# an ad hoc pre-processing.
#
# Details
# The following indices are calculated based on input daily precipitation data:
# PRY = mean annual precipitation
# INT = mean annual precipitation intensity (intensity during wet days, or simple precipitation intensity index SDII)
# WSL = mean annual wet spell length (number of consecutive days during each wet spell)
# DSL = mean annual dry spell lenght (number of consecutive days during each dry spell)
# PA  = precipitation area (area over which of any given day precipitation occurs)
# R95 = heavy precipitation index (percent of total precipitation above the 95% percentile of the reference distribution)
# HY-INT = hydroclimatic intensity. HY-INT = normalized(INT) x normalized(DSL).
#
# For EC-Earth data and then extended to any model and observational data, producing plots
# of data vs. a reference dataset (e.g. ERA-INTERIM). Indices are normalized over a reference
# period. Both absolute and normalized values are made available: users can select the indices
# to be stored and plotted. The tool makes extensives use of the cfg_hyint configuration file
# for user selectable options and ease feeding needed inputs (e.g. region boundaries for timeseries
# or value ranges and labels for figures).
#
# Required
# It reads daily precipitation data through ESMValTool. If requested, input precipitation data are pre-processed
# interpolating on a common grid set by the user in the hyint_parameters file.
# R libraries:"tools","PCICt","ncdf4","maps"
#
# Optional
# Several options can be selected via the configuration file, e.g. the provision of an
# external normalization functions for the indices; a reference climatology for the R95 index; type of plots; etc.
#
# Caveats
# Spatial data selection based on elevation works only with regridding at 320x160 (or by producing by hand grid files at needed resolution)
#
# Modification history
#    20181001-A_arno_en: converted to latest v2.0
#    20180302-A_arno_en: converted to ESMValTool2
#    20171206-A_arno_en: modularized version accepting climdex indices
#    20171010-A_arno_en: modularized version
#    20170901-A_arno_en: 1st github version
#
# ############################################################################

library(tools)
library(yaml)
library(ncdf4)

# get path to script and source subroutines (if needed)
args <- commandArgs(trailingOnly = FALSE)
spath <- paste0(dirname(unlist(strsplit(grep("--file", args,
                                             value = TRUE), "="))[2]), "/")

source(paste0(spath, "hyint_functions.R"))
source(paste0(spath, "hyint_metadata.R"))
source(paste0(spath, "hyint_preproc.R"))
source(paste0(spath, "hyint_diagnostic.R"))
source(paste0(spath, "hyint_etccdi_preproc.R"))
source(paste0(spath, "hyint_trends.R"))
source(paste0(spath, "hyint_plot_maps.R"))
source(paste0(spath, "hyint_plot_trends.R"))
source(paste0(spath, "hyint_parameters.R"))

diag_script_cfg <- paste0(spath, "hyint_parameters.R")

# Read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings_file <- args[1]
settings <- yaml::read_yaml(settings_file)

# load data from settings
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}

metadata <- yaml::read_yaml(settings$input_files)

# get name of climofile for selected variable and 
# list associated to first climofile
climofiles <- names(metadata)
climolist <- get(climofiles[1], metadata)
list0 <- climolist
climolist0 <- list0

# get variable name
varname <- paste0("'", climolist$short_name, "'")
var0 <- varname
var0 <- "pr"

diag_base <- climolist0$diagnostic
print(paste0(diag_base, ": starting routine"))

if (length(etccdi_dir) != 1) {
  etcddi_dir <- work_dir
}
regridding_dir <- run_dir
dir.create(plot_dir, recursive = T, showWarnings = F)
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)

# extract metadata
models_name <- unname(sapply(metadata, "[[", "dataset"))
reference_model <- unname(sapply(metadata, "[[", "reference_dataset"))[1]
models_start_year <- unname(sapply(metadata, "[[", "start_year"))
models_end_year <- unname(sapply(metadata, "[[", "end_year"))
models_experiment <- unname(sapply(metadata, "[[", "exp"))
models_ensemble <- unname(sapply(metadata, "[[", "ensemble"))

# select reference dataset; if not available, use last of list
ref_idx <- which(models_name == reference_model)
if (length(ref_idx) == 0) {
  ref_idx <- length(models_name)
}

## Run regridding and diagnostic
if (write_netcdf) {

  # loop through models
  for (model_idx in c(1:(length(models_name)))) {

    # Setup filenames
    climofile <- climofiles[model_idx]
    sgrid <- "noregrid"
    if (rgrid != F) {
      sgrid <- rgrid
    }
    regfile <- getfilename_regridded(regridding_dir, sgrid, var0, model_idx)

    # If needed, pre-process file and adding absolute time axis
    if (run_regridding) {
      if (!file.exists(regfile) | force_regridding) {
        dummy <- hyint_preproc(work_dir, model_idx, climofile, regfile)
      } else {
        gridfile <- getfilename_indices(work_dir, diag_base, model_idx,
                                        grid = T)
        grid_command <- paste("cdo griddes ", regfile, " > ", gridfile)
        system(grid_command)
        print(paste0(diag_base, ": data file exists: ", regfile))
        print(paste0(diag_base, ": corresponding grid: ", gridfile))
      }
    }

    if (run_diagnostic) {
      # Loop through seasons and call diagnostic
      for (seas in seasons) {
        hyint_diagnostic(work_dir, regfile, model_idx, seas,
                         rewrite = force_diagnostic)
      }
    }
  }
}

## Preprocess ETCCDI input files and merge them with HyInt indices
if (write_netcdf & etccdi_preproc) {
  for (model_idx in c(1:(length(models_name)))) {
    gridfile <- getfilename_indices(work_dir, diag_base, model_idx, grid = T)
    dummy <- hyint_etccdi_preproc(work_dir, etccdi_dir, etccdi_list_import,
                                  gridfile, model_idx, "ALL", yrmon = "yr")
  }
}

## Calculate timeseries and trends
if (write_netcdf & run_timeseries) {
  for (model_idx in c(1:(length(models_name)))) {
    for (seas in seasons) {
      hyint_trends(work_dir, model_idx, seas)
    }
  }
}

## Create figures
if (write_plots) {
  plot_type_list <- plot_type
   # select reference dataset; if not available, use last of list
  ref_idx <- which(models_name == reference_model)
  if (length(ref_idx) == 0) {
    ref_idx <- length(models_name)
  }
  for (plot_type in plot_type_list) {
    print(paste0("***** PLOT TYPE: ", plot_type, " *********"))
    for (seas in seasons) {
      if (plot_type <= 10) {
        # Plot maps
        hyint_plot_maps(work_dir, plot_dir, work_dir, ref_idx, seas)
      } else {
        # Plot timeseries and trends
        hyint_plot_trends(work_dir, plot_dir, work_dir, ref_idx, seas)
      }
    }
  }
}

print(paste0(diag_base, ": done."))
