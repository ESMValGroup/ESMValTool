# OPTIONS:
#
# write_coeffs: boolean, default FALSE
#   write xi, alpha and kappa to netcdf files
#   TODO: hardcoded for this 3 parameters. Only works for log-Logistic yet.
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
#
# Authors: [Peter Berg, Katja Weigel, Lukas Lindenlaub]

library(yaml)
library(ncdf4)
library(SPEI)
library(R.utils)

setwd(dirname(commandArgs(asValues=TRUE)$file))
source("utils.R")

fill_refperiod <- function(cfg, tsvec) {
  cfg <- list_default(cfg, "refstart_year", tsvec[1])
  cfg <- list_default(cfg, "refstart_month", tsvec[2])
  cfg <- list_default(cfg, "refend_year", tsvec[3])
  cfg <- list_default(cfg, "refend_month", tsvec[4])
  return(cfg)
}

# ---------------------------------------------------------------------------- #
# Script starts here --------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
cfg <- read_yaml(commandArgs(trailingOnly = TRUE)[1])
cfg <- list_default(cfg, "write_coeffs", FALSE)
cfg <- list_default(cfg, "write_wb", FALSE)
cfg <- list_default(cfg, "short_name_pet", "evspsblpot")
cfg <- list_default(cfg, "distribution", "log-Logistic")
dir.create(cfg$work_dir, recursive = TRUE)
dir.create(cfg$plot_dir, recursive = TRUE)
fillfloat <- 1.e+20
as.single(fillfloat)
provenance_file <- paste0(cfg$run_dir, "/", "diagnostic_provenance.yml")
meta_file <- paste0(cfg$work_dir, "/metadata.yml")
provenance <- list()
meta <- list()  # collect output files metadata

print("--- Load input meta")
ancestor_meta <- list()
for (anc in cfg$input_files){
  if (!endsWith(anc, 'metadata.yml')){anc <- paste0(anc, "/metadata.yml")}
  ancestor_meta <- append(ancestor_meta, read_yaml(anc))
}
grouped_meta <- group_meta(ancestor_meta)

print("--- Process each dataset")
for (dataset in names(grouped_meta)){
  print(paste("-- processing dataset:", dataset))
  metas <- grouped_meta[[dataset]]  # list of files for this dataset
  mask <- get_merged_mask(metas) # use mask per dataset
  pr_meta <- select_var(metas, "pr")
  pr <- get_var_from_nc(pr_meta)
  lat <- get_var_from_nc(pr_meta, custom_var="lat")
  tsvec <- get_var_from_nc(pr_meta, custom_var="time")
  pet_meta <- select_var(metas, cfg$short_name_pet, strict=FALSE)
  if (is.null(pet_meta)) {
    cfg$indexname <- "SPI"
    cfg$long_name <- "Standardized Precipitation Index"
    pme <- pr
  } else {
    cfg$indexname <- "SPEI"
    cfg$long_name <- "Standardized Precipitation Evapotranspiration Index"
    pet <- get_var_from_nc(pet_meta)
    pme <- pr - pet
  }
  if (cfg$write_wb) {
    filename_wb <- write_nc_file_like(cfg, pr_meta, pme, fillfloat, short_name="wb")
    meta[[filename_wb]] <- list(filename=filename_wb, short_name="wb", dataset = dataset)
    provenance[[filename_wb]] <- list(caption="Water balance per grid point.")
  }

  cfg <- fill_refperiod(cfg, tsvec)
  pme_spei <- pme * NA
  coeffs <- list()
  for (i in 1:dim(pme)[1]){
    wh <- which(!is.na(mask[i,]))
    if (length(wh) > 1){
      tmp <- pme[i, wh,]
      ts_data <- ts(t(tmp), freq=12, start=c(tsvec[1], tsvec[2]))
      spei_results <- spei(
        ts_data,
        cfg$smooth_month,
        na.rm = TRUE,
        distribution = cfg$distribution,
        ref.start = c(cfg$refstart_year, cfg$refstart_month),
        ref.end = c(cfg$refend_year, cfg$refend_month)
      )
      pme_spei[i, wh, ] <- t(spei_results$fitted)
      if (cfg$write_coeffs) {
        for (c_name in rownames(spei_results$coefficients)){
          if (is.null(coeffs[[c_name]])){
            coeffs[[c_name]] <- array(NA, dim=c(dim(pme)[1], dim(pme)[2], 12))
          }
          coeffs[[c_name]][i, wh, ] <- spei_results$coefficients[c_name, ,]
        }
      }
    }
  }
  pme_spei[pme_spei > 10000] <- NA  # replaced with fillfloat in write function
  # apply mask
  for (t in 1:dim(pme)[3]) {
    tmp <- pme_spei[, , t]
    tmp[is.na(mask)] <- NA
    pme_spei[, , t] <- tmp
  }
  filename <- write_nc_file_like(cfg, pr_meta, pme_spei, fillfloat, short_name=cfg$indexname)
  new_meta = list(filename=filename, short_name=cfg$indexname,
      long_name=cfg$indexname, units="1", dataset=dataset)
  meta[[filename]] <- modifyList(pr_meta, new_meta)
  provenance[[filename]] = list(caption=paste(cfg$indexname, " index per grid point."))
  if (cfg$write_coeffs) {
    for (c_name in names(coeffs)){
      filename_c <- write_nc_file_like(cfg, pr_meta, coeffs[[c_name]], fillfloat, short_name=c_name, moty=TRUE)
      meta[[filename_c]] <- list(filename=filename_c, short_name=c_name, dataset=dataset)
      provenance[[filename_c]] <- list(caption=paste(c_name, " (fitting parameter) per grid point."))
    }
  }
}  # end of dataset loop

write_yaml(provenance, provenance_file)
write_yaml(meta, meta_file)
