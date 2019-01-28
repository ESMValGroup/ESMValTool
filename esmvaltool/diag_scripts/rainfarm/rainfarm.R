# #############################################################################
# rainfarm.R
# Authors:       E. Arnone (ISAC-CNR, Italy)
# 	         J. von Hardenberg (ISAC-CNR, Italy)
# #############################################################################
# Description
# ESMValTool diagnostic calling the RainFARM library written in Julia (by von Hardenberg, ISAC-CNR, Italy).
# RainFARM is a stochastic precipitation downscaling method, further adapted for climate downscaling.
#
# Required
# CDO
# Julia language: https://julialang.org
# RainFARM Julia library: https://github.com/jhardenberg/RainFARM.jl
#
# Optional
#
# Caveats
#
# Modification history
#    20181210 hard_jo: cleanup and using juliacall
#    20180508-A_arnone_e: Conversion to v2.0
#    20170908-A_arnone_e: 1st github version
#
# ############################################################################

library(tools)
library(yaml)
library(JuliaCall)  #nolint

julia_setup()
julia_library("RainFARM")

diag_scripts_dir <- Sys.getenv("diag_scripts")

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}
metadata <- yaml::read_yaml(settings$input_files)

# get name of climofile for selected variable and list associated to first climofile
climofiles <- names(metadata)
climolist <- get(climofiles[1], metadata)

varname <- climolist$short_name

diag_base <- climolist$diagnostic
print(paste(diag_base, ": starting routine"))

# create working dirs if they do not exist
work_dir <- settings$work_dir
regridding_dir <- settings$run_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)

# switch to working directory
setwd(work_dir)

# setup provenance file and list
provenance_file <- paste0(regridding_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

# extract metadata
models_name <- unname(sapply(metadata, "[[", "dataset"))
reference_model <- unname(sapply(metadata, "[[", "reference_dataset"))[1]
models_start_year <- unname(sapply(metadata, "[[", "start_year"))
models_end_year <- unname(sapply(metadata, "[[", "end_year"))
models_experiment <- unname(sapply(metadata, "[[", "exp"))
models_ensemble <- unname(sapply(metadata, "[[", "ensemble"))

# Loop through input models, apply pre-processing and call RainFARM
for (model_idx in c(1:(length(models_name)))) {

  exp <- models_name[model_idx]
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]
  infile <- climofiles[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  infilename <- file_path_sans_ext(basename(infile))

  print(paste0(diag_base, ": calling rainfarm"))
  outfilename <- paste0(work_dir, "/", infilename, "_downscaled")

  # call rfarm
  ans <- julia_call("read_netcdf2d", infile, varname, need_return = "R" )
  pr <- ans[[1]]
  lon_mat <- ans[[2]]
  lat_mat <- ans[[3]]
  ans <- julia_call("lon_lat_fine", lon_mat, lat_mat, nf, need_return = "R" )
  lon_f <- ans[[1]]
  lat_f <- ans[[2]]
  if (slope == 0){
      ans <- julia_call("fft3d", pr, need_return = "R"  )
      fxp <- ans[[1]]
      ftp <- ans[[2]]
      sx <- julia_call("fitslopex", fxp, kmin = 1, need_return = "R" )
      print(paste0("Computed spatial spectral slope: ", sx))
  } else {
      sx <- slope
      print(paste0("Fixed spatial spectral slope: ", sx))
  }
  if (weights_climo != F) {
      print(paste0("Using external climatology for weights: ", weights_climo))
      fileweights <- paste0(work_dir, "/", infilename, "_w.nc")
      ww <- julia_call("rfweights", weights_climo, infile, nf,
                      weightsfn = fileweights, varname = varname,
                      fsmooth = conserv_smooth, need_return = "R");
  } else {
    print("Not using weights")
    ww <- 1.
  }
  if (conserv_glob) {
     print("Conserving global field")
  } else if (conserv_smooth) {
     print("Smooth conservation")
  } else {
     print("Box conservation")
  }
  for (iens in 1:nens) {
      print(paste0("Realization ", iens))
      rd <- julia_call("rainfarm", pr, sx, nf, ww,
                              fglob = conserv_glob,
                              fsmooth = conserv_smooth,
                              verbose = T, need_return = "R")
      fname <- sprintf("%s_%04d.nc", outfilename, iens)
      julia_call("write_netcdf2d", fname, rd, lon_f, lat_f, varname, infile )

      # Set provenance for this output file
      caption <- paste0("RainFARM precipitation downscaling")
      xbase <- list(ancestors = list(infile),
                    authors = list("arno_en", "hard_jo"),
                    references = list("donofrio14jh", "rebora06jhm",
                                      "terzago18nhess"),
                    projects = list("c3s-magic"),
                    caption = caption,
                    statistics = list("other"),
                    realms = list("atmos"),
                    themes = list("phys"),
                    domains = list("reg"))

      # Store provenance in main provenance list
      provenance[[fname]] <- xbase
  }
  # Write provenance to file
  write_yaml(provenance, provenance_file)
}
