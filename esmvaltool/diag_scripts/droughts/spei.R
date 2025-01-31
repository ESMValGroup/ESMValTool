# DESCRIPTION:
#
# This file is based on code from diag_save_spei_all.R and should replace
# diag_spei.R at some point. It is more flexible and provides the same
# functionality as the previous diag_spei.R as a special case.
# Calculation of PET and plotting of results is covered in seperate diagnostics.
#
# Authors: [Peter Berg, Katja Weigel, Lukas Ruhe]
#
# CHANGELOG:
#
# NOTE: modified PET to work with ancestors. Maybe this breaks internal PET,
# when ancestors are not set.
#
# NOTE: with given ancestors the metadata files are not longer part of
# settings.yml. Variables from ancestors can be added in recipe (i.e. pet/pr)
#
# NOTE: In the output folder are now (ESMValTool update) additinal xml and yml
# files which need to be skipped when manually loaded. 
# UPDATE: reduced this diag to work only with ancestors producing metadata.yml 
# or preprocessed variables.
# 
# NOTE: Is the reference dataset used to apply the mask (NaNs) to all other
# datasets? What should be chosen for reference? and why?
# UPDATE: refrence_dataset will be read from script settings instead of dataset
# UPDATE2: individual masks and refperiods should be fine. Reference dataset
# completly removed.
# 
# NOTE: All metadata is kept in a named list with model keys. Removed some loops 
# and variables. In most functions meta just replaces yml[m][1].
#
# NOTE: A similar correction for time unit was hardcoded for each
# variable. Changed it to a function that is called multiple times. #DRY
# 
# UPDATE: cleaned up getnc/getpetnc -> get_var_from_nc() getpetnc 
# UPDATE: gettimenc is another special case of get_var_from_nc() #DRY
#
# NOTE: latitude is taken from reference dataset, but also for each ds during
# calculation. 
# 
# NOTE: Common functions moved to utils.R this includes general read and write 
# functions for nc files that replace ncwrite, ncwritespei, ncwritepet, getpetnc
# gettimenc... and general utility functions like default values for lists
# 
# NOTE: move all if conditions for each refperiod param into a seperate function
# fill_refperiod. #DRY
#
# NOTE: added optional parameter `short_name_pet` to use variables other than 
# evspsblpot from recipe.
# 
# NOTE: set log-Logistic as default distribution if nothing is given in the 
# recipe. Missing distribution raised an unclear error before.
#
# OPTIONS:
#
# write_coeffs: boolean, default FALSE
#   write xi, alpha and kappa to netcdf files
# write_wb: boolean, default FALSE
#   write water balance to netcdf file
# short_name_pet: string, default "evspsblpot"
#   short name of the variable to use as PET (i.e. ET)
# distributionn: string, default "log-Logistic"
#   type of distribution used for SPEI calibration. 
#   Options: "Gamma", "log-Logistic", "Pearson III"
# refstart_year: integer, default first year of time series
# refstart_month: integer, default 1
# refend_year: integer, default last year of time series
# refend_month: integer, default 12

library(yaml)
library(ncdf4)
library(SPEI)
library(R.utils)

setwd(dirname(commandArgs(asValues=TRUE)$file))
source("utils.R")

fill_refperiod <- function(params, tsvec) {
  params <- list_default(params, "refstart_year", tsvec[1])
  params <- list_default(params, "refstart_month", tsvec[2])
  params <- list_default(params, "refend_year", tsvec[3])
  params <- list_default(params, "refend_month", tsvec[4])
}

# ---------------------------------------------------------------------------- #
# Script starts here --------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
params <- read_yaml(commandArgs(trailingOnly = TRUE)[1])
params <- list_default(params, "write_coeffs", FALSE)
params <- list_default(params, "write_wb", FALSE)
params <- list_default(params, "short_name_pet", "evspsblpot")
params <- list_default(params, "distribution", "log-Logistic")
dir.create(params$work_dir, recursive = TRUE)
dir.create(params$plot_dir, recursive = TRUE)
fillfloat <- 1.e+20
as.single(fillfloat)
provenance_file <- paste0(params$run_dir, "/", "diagnostic_provenance.yml")
meta_file <- paste0(params$work_dir, "/metadata.yml")
provenance <- list()
meta <- list()  # collect output files metadata

print("--- Load input meta")
ancestor_meta <- list()
for (anc in params$input_files){
  if (!endsWith(anc, 'metadata.yml')){anc <- paste0(anc, "/metadata.yml")}
  ancestor_meta <- append(ancestor_meta, read_yaml(anc))
}
grouped_meta <- group_meta(ancestor_meta)

print("--- Process each dataset")
for (dataset in names(grouped_meta)){
  print(paste("-- processing dataset:", dataset))
  metas <- grouped_meta[[dataset]]  # list of files for this dataset
  mask <- get_merged_mask(metas) # use mask per dataset
  # START SPEI CALC
  pr_meta <- select_var(metas, "pr")
  pr <- get_var_from_nc(pr_meta)
  lat <- get_var_from_nc(pr_meta, custom_var="lat")
  tsvec <- get_var_from_nc(pr_meta, custom_var="time")
  pet <- get_var_from_nc(select_var(metas, params$short_name_pet))
  pme <- pr - pet
  if (params$write_wb) {
    filename_wb <- write_nc_file_like(params, pr_meta, pme, fillfloat, short_name="wb")
  }

  fill_refperiod(params, tsvec)
  pme_spei <- pme * NA
  coeffs <- array(numeric(), c(3, dim(pme)[1], dim(pme)[2], 12)) # xi, alpha, kappa
  for (i in 1:dim(pme)[1]){
    wh <- which(!is.na(mask[i,]))
    if (length(wh) > 1){
      tmp <- pme[i, wh,]
      ts_data <- ts(t(tmp), freq=12, start=c(tsvec[1], tsvec[2]))
      spei_results <- spei(
        ts_data,
        params$smooth_month,
        na.rm = TRUE,
        distribution = params$distribution,
        ref.start = c(params$refstart_year, params$refstart_month),
        ref.end = c(params$refend_year, params$refend_month)
      )
      coeffs[, i, wh, ] <- spei_results$coefficients
      pme_spei[i, wh, ] <- t(spei_results$fitted)
    }
  }
  pme_spei[pme_spei > 10000] <- NA  # replaced with fillfloat in write function
  filename <- write_nc_file_like(params, pr_meta, pme_spei, fillfloat)
  if (params$write_coeffs) {
    filename_xi <- write_nc_file_like(
      params, pr_meta, coeffs[1, , ,], fillfloat, 
      short_name = "xi", 
      moty=TRUE)
    meta[[filename_xi]] <- list(
      filename=filename_xi, 
      short_name="xi", 
      dataset=dataset)
    filename_alpha <- write_nc_file_like(
      params, pr_meta, coeffs[2, , ,], fillfloat, 
      short_name = "alpha", 
      moty = TRUE)
    meta[[filename_alpha]] <- list(
      filename=filename_xi, 
      short_name="alpha", 
      dataset=dataset)
    filename_kappa <- write_nc_file_like(
      params, pr_meta, coeffs[3, , ,], fillfloat, 
      short_name = "kappa", 
      moty=TRUE)
    meta[[filename_kappa]] <- list(
      filename=filename_kappa, 
      short_name="kappa", 
      dataset=dataset)
  }

  print("-- prepare metadata for output")
  meta[[filename]] <- list(
    filename = filename,
    short_name = "spei",
    dataset = dataset)
  xprov$caption <- "SPEI index per grid point."
  # generate metadata.yml
  input_meta = select_var(metas, "pr")
  input_meta$filename = filename
  input_meta$short_name = "spei"
  input_meta$long_name = "SPEI"
  input_meta$units = "1"
  input_meta$index = "SPEI"
  meta[[filename]] <- input_meta
  # meta[[filename]][["index"]] <- "SPEI"
  for (t in 1:dim(pme)[3]) {
    tmp <- pme_spei[, , t]
    tmp[is.na(mask)] <- NA
    pme_spei[, , t] <- tmp
  }
}  # end of big dataset loop

provenance[[filename]] <- xprov
write_yaml(provenance, provenance_file)
write_yaml(meta, meta_file)