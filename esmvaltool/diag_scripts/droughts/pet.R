#' @importFrom checkmate makeAssertCollection
#' @importFrom stats cycle end frequency start ts optim
#' @importFrom zoo as.yearmon rollapply

# This file is used to calculate potential evapotranspiration (PET) based on
# the SPEI.R package. It can be used as ancestor for the SPEI calculation.
# Documentation can be found at:
# https://docs.esmvaltool.org/en/latest/recipes/recipe_droughts.html
# https://cran.r-project.org/web/packages/SPEI/SPEI.pdf
#
# Authors: [Peter Berg, Katja Weigel, Lukas Lindenlaub]

library(yaml)
library(ncdf4)
library(SPEI)
library(R.utils)

setwd(dirname(commandArgs(asValues=TRUE)$file))
source("utils.R")


calculate_hargreaves <- function(metas, xprov, use_pr=FALSE) {
  data <- list(tasmin=NULL, tasmax=NULL, rsdt=NULL, pr=NULL)
  for (meta in metas) {
    if (meta$short_name %in% names(data) && !(meta$short_name == "pr" && !use_pr)) {
      print(paste("read", meta$short_name))
      data[[meta$short_name]] <- get_var_from_nc(meta)
      if (meta$short_name == "tasmin") {
        data$lat <- get_var_from_nc(meta, custom_var="lat")
      }
      print(paste("Shape of", meta$short_name))
      print(dim(data[[meta$short_name]]))
    } else {
      print(paste("Variable", meta$short_name, "not used for hargreaves"))
    }
  }
  dpet <- data$tasmin * NA
  for (i in 1:dim(dpet)[2]) {
    pet_tmp <- hargreaves(
      t(data$tasmin[, i, ]),
      t(data$tasmax[, i, ]),
      lat = rep(data$lat[i], dim(dpet)[1]),
      Pre = t_or_null(data$pr[, i, ]),
      Ra = t_or_null(data$rsdt[, i, ]),
      na.rm = TRUE
    )
    d2 <- dim(pet_tmp)
    pet_tmp <- as.numeric(pet_tmp)
    dim(pet_tmp) <- d2
    dpet[, i, ] <- t(pet_tmp)
  }
  return(dpet)
}

calculate_thornthwaite <- function(metas, xprov) {
  for (meta in metas) {
    if (meta$short_name == "tas") {
      tas = get_var_from_nc(meta)
      lat = get_var_from_nc(meta, custom_var="lat")
      xprov$ancestors <- append(xprov$ancestors, meta$filename)
    } else {
      print(paste("Variable", meta$short_name, "not used for Thornthwaite."))
    }
  }
  dpet <- tas * NA
  for (i in 1:dim(dpet)[2]){
    tmp <- tas[, i, ]
    pet_tmp <- thornthwaite(
      t(tmp), rep(lat[i], dim(dpet)[1]), na.rm = TRUE
    )
    d2 <- dim(pet_tmp)
    pet_tmp <- as.numeric(pet_tmp)
    dim(pet_tmp) <- d2
    dpet[, i, ] <- t(pet_tmp)
  }
  return(dpet)
}


calculate_penman <- function(metas, xprov, method="ICID", crop="tall") {
  data <- list(
    tasmin = NULL,
    tasmax = NULL,
    clt = NULL,
    sfcWind = NULL,
    ps = NULL,
    psl = NULL,
    hurs = NULL,
    rsds = NULL,
    rsdt=NULL)
  # load relevant variables
  for(meta in metas){
    if (meta$short_name %in% names(data)){
      data[[meta$short_name]] <- get_var_from_nc(meta)
      xprov$ancestors <- append(xprov$ancestors, meta$filename)
      if(meta$short_name == "tasmin"){
        data$lat <- get_var_from_nc(meta, custom_var="lat")
      }
    } else {
      print("variable not used for penman:")
      print(meta$short_name)
    }
  }
  # calculate pet for each latitude
  dpet <- data$tasmin * NA
  for (i in 1:dim(dpet)[2]){
    pet_tmp <- penman(
      t(data$tasmin[, i, ]),
      t(data$tasmax[, i, ]),
      t(data$sfcWind[, i, ]),
      lat = rep(data$lat[i], dim(dpet)[1]),
      CC = t_or_null(data$clt[, i, ]),
      P = t_or_null(data$ps[, i, ]),
      P0 = t_or_null(data$psl[, i, ]),
      RH = t_or_null(data$hurs[, i, ]),
      Ra = t_or_null(data$rsdt[, i, ]),
      Rs = t_or_null(data$rsds[, i, ]),
      na.rm = TRUE,
      method = method,
      crop = crop
    )
    d2 <- dim(pet_tmp)
    pet_tmp <- as.numeric(pet_tmp)
    dim(pet_tmp) <- d2
    dpet[, i, ] <- t(pet_tmp)
  }
  return(dpet)
}


# ---------------------------------------------------------------------------- #
# Script starts here --------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

params <- read_yaml(commandArgs(trailingOnly = TRUE)[1])
params$method = ifelse(!is.null(params$method), params$method, "ICID")
params$crop = ifelse(!is.null(params$crop), params$crop, "tall")
dir.create(params$work_dir, recursive = TRUE)
dir.create(params$plot_dir, recursive = TRUE)
fillfloat <- 1.e+20
as.single(fillfloat)
provenance_file <- paste0(params$run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()
meta <- list()  # collect output files metadata
meta_file <- paste0(params$work_dir, "/metadata.yml")

print("--- Load metadata")
ancestor_meta <- list()
for (anc in params$input_files){
  if (!endsWith(anc, 'metadata.yml')){anc <- paste0(anc, "/metadata.yml")}
  ancestor_meta <- append(ancestor_meta, read_yaml(anc))
}
grouped_meta <- group_meta(ancestor_meta)

print("--- Process each dataset")
for (dataset in names(grouped_meta)){
  metas <- grouped_meta[[dataset]]  # list of files for this dataset
  xprov$ancestors <- list()
  switch(params$pet_type,
    Penman = {pet <- calculate_penman(metas, xprov,
      method=params$method, crop=params$crop)},
    Thornthwaite = {pet <- calculate_thornthwaite(metas, xprov)},
    Hargreaves = {pet <- calculate_hargreaves(metas, xprov)},
    stop("pet_type must be one of: Penman, Hargreaves, Thornthwaite")
  )
  mask <- get_merged_mask(metas)
  for (t in 1:dim(pet)[3]){
    tmp <- pet[, , t]
    tmp[is.na(mask)] <- NA
    pet[, , t] <- tmp
  }
  pet <- monthly_to_daily(pet, dim=3)

  # write PET to file
  first_meta = metas[[names(metas)[1]]]
  filename <- write_nc_file_like(
    params, first_meta, pet, fillfloat,
    short_name="evspsblpot",
    long_name="Potential Evapotranspiration",
    units="mm day-1")
  input_meta <- select_var(metas, "tasmin")
  input_meta$filename <- filename
  input_meta$short_name <- "evspsblpot"
  input_meta$long_name <- "Potential Evapotranspiration"
  input_meta$units <- "mm day-1"
  meta[[filename]] <- input_meta

  xprov$caption <- "PET per grid point."
  provenance[[filename]] <- xprov
}

write_yaml(provenance, provenance_file)
write_yaml(meta, meta_file)
