# #############################################################################
# diagnostic.R
# Authors:       Irene Cionni (ENEA, Italy)
# #############################################################################
# Description
# This script modifies the wind capacity factor diagnostic wrote for the
# MAGIC project from BSC, see also esmvaltool/diag_scripts/magic_bsc/.
#
# Required
# season: String to include shortcut for season in plot title
#         and name (e.g. "djf"). It will be converted to upper case.
#         This season should be the one set in the preprocessor, since it is
#         only used as a string and does not affect the data in the diagnostic.
#         In the default recipe this is solved through a node anchor.
#
# Optional
# maxval_colorbar: Optional upper limit for the colorbar.
#
# Caveats
#
# Modification history
#    20210401-cionni_irene: written for v2.0
#    20210401-weigel_katja: changed to allow multiple models
#    20210621-weigel_katja: formattiong updates
#
# ############################################################################
library(abind)
library(climdex.pcic)
library(ggplot2)
library(multiApply) # nolint
library(ncdf4)
library(RColorBrewer) # nolint
library(s2dverification)
library(yaml)

# Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(
  file_arg_name, "",
  initial_options[grep(file_arg_name, initial_options)]
)
script_dirname <- dirname(script_name)

source(file.path(script_dirname, "PV_CF.R"))
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)
input_files_per_var <- yaml::read_yaml(params$input_files[1])
var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x) x$dataset)
model_names <- unname(model_names)
var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
input_files_per_var1 <- yaml::read_yaml(params$input_files[2])
var_names1 <- names(input_files_per_var1)
var1 <- lapply(input_files_per_var1, function(x) x$short_name)
fullpath_filenames1 <- names(var1)
var0 <- unname(var0)[1]
var1 <- unname(var1)[1]
start_years <- lapply(input_files_per_var, function(x) x$start_year)
start_years <- unname(start_years)
end_years <- lapply(input_files_per_var, function(x) x$end_year)
end_years <- unname(end_years)
seasons <- toupper(params$season)

var0 <- unlist(var0)
for (i in seq(1, length(model_names))) {
  start_year <- c(unlist(start_years[i]))
  end_year <- c(unlist(end_years[i]))
  no_of_years <- length(seq(start_year, end_year, 1))
  data_nc <- nc_open(fullpath_filenames[i])
  data <- ncvar_get(data_nc, var0)
  names(dim(data)) <- c("lon", "lat", "time")
  lat <- ncvar_get(data_nc, "lat")
  lon <- ncvar_get(data_nc, "lon")
  units <- ncatt_get(data_nc, var0, "units")$value
  calendar <- ncatt_get(data_nc, "time", "calendar")$value
  long_names <- ncatt_get(data_nc, var0, "long_name")$value
  time <- ncvar_get(data_nc, "time")
  start_date <- as.POSIXct(substr(ncatt_get(
    data_nc, "time",
    "units"
  )$value, 11, 29))
  nc_close(data_nc)
  time <- as.Date(time,
    origin = substr(start_date, 1, 10),
    calendar = calendar
  )
  time <- as.POSIXct(time, format = "%Y-%m-%d")
  time_dim <- which(names(dim(data)) == "time")
  time <- as.PCICt(time, cal = calendar)
  if (calendar != "360_day" & calendar != "365_day") {
    time <- as.character(time)
    jdays <- as.numeric(strftime(time, format = "%j"))
    pos <- which(substr(time, 6, 10) == "02-29")
    if (length(pos) > 0) {
      time <- time[-pos]
      data <- apply(
        data, c(seq(1, length(dim(data)), 1))[-time_dim],
        function(x) {
          x[-pos]
        }
      )
      data <- aperm(data, c(2, 3, 1))
      names(dim(data)) <- c("lon", "lat", "time")
    }
  }
  dims <- dim(data)
  dims <- append(dims[-time_dim], c(no_of_years, dims[time_dim] /
    no_of_years), after = 2)
  dim(data) <- dims
  data <- aperm(data, c(3, 4, 2, 1))
  names(dim(data)) <- c("year", "day", "lat", "lon")
  ######## var1#########################################
  var1 <- unlist(var1)
  data_nc1 <- nc_open(fullpath_filenames1[i])
  data1 <- ncvar_get(data_nc1, var1)

  names(dim(data1)) <- c("lon", "lat", "time")
  lat <- ncvar_get(data_nc1, "lat")
  lon <- ncvar_get(data_nc1, "lon")
  units <- ncatt_get(data_nc1, var1, "units")$value
  calendar <- ncatt_get(data_nc1, "time", "calendar")$value
  long_names <- ncatt_get(data_nc1, var1, "long_name")$value
  time <- ncvar_get(data_nc1, "time")
  start_date <- as.POSIXct(substr(ncatt_get(
    data_nc1, "time",
    "units"
  )$value, 11, 29))
  nc_close(data_nc1)
  time <- as.Date(time,
    origin = substr(start_date, 1, 10),
    calendar = calendar
  )
  time <- as.POSIXct(time, format = "%Y-%m-%d")
  time_dim <- which(names(dim(data1)) == "time")
  time <- as.PCICt(time, cal = calendar)
  if (calendar != "360_day" & calendar != "365_day") {
    time <- as.character(time)
    jdays <- as.numeric(strftime(time, format = "%j"))
    pos <- which(substr(time, 6, 10) == "02-29")
    if (length(pos) > 0) {
      time <- time[-pos]
      data1 <- apply(
        data1, c(seq(1, length(dim(data1)), 1))[-time_dim],
        function(x) {
          x[-pos]
        }
      )
      data1 <- aperm(data1, c(2, 3, 1))
      names(dim(data1)) <- c("lon", "lat", "time")
    }
  }
  dims1 <- dim(data1)
  dims1 <- append(dims1[-time_dim], c(no_of_years, dims1[time_dim] /
    no_of_years), after = 2)
  dim(data1) <- dims1
  data1 <- aperm(data1, c(3, 4, 2, 1))
  names(dim(data1)) <- c("year", "day", "lat", "lon")

  #####################################
  # CF modelC
  ####################################

  seas_data <- Mean1Dim(data, 2)
  data_cf1 <- rsds2cf(data1, data)
  dim(data_cf1) <- dim(data)
  #---------------------------
  # Aggregate daily data to seasonal means
  #---------------------------

  seas_data_cf1 <- Mean1Dim(data_cf1, 2)
  ##############################
  # Make some plots
  ##############################
  #---------------------------
  # Prepare data, labels and colorscales
  #---------------------------
  p <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  q <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  years <- seq(start_year, end_year)
  turb_types <- c("PVCF")

  seas_data_cf_all <- seas_data_cf1

  mean_data_cf_all <- Mean1Dim(seas_data_cf_all, 1)

  anom_data_cf_all <- seas_data_cf_all - InsertDim( # nolint
    Mean1Dim(seas_data_cf_all, 1), 1, dim(data)[1]
  ) # nolint
  pct_anom_data_cf_all <- (seas_data_cf_all / InsertDim( # nolint
    Mean1Dim(seas_data_cf_all, 1), 1, dim(data)[1]
  )) - 1 # nolint
  #---------------------------
  # Plot seasonal CF maps
  #---------------------------
  filepng <- paste0(
    plot_dir, "/", "capacity_factor_", model_names[i], "_",
    start_year, "-", end_year, "_", seasons, ".png"
  )
  title <- paste0(
    seasons, " CF from ", model_names[i],
    " (", start_year, "-", end_year, ")"
  )

  # Optional upper limit for the color bar set in recipe
  if (length(params$maxval_colorbar) == 1) {
    maxval_colorbar <- params$maxval_colorbar
  } else {
    maxval_colorbar <- max(seas_data_cf_all, na.rm = TRUE)
  }

  PlotEquiMap(Mean1Dim(seas_data_cf_all, 1), lon, lat,
    filled.continents = F,
    brks = seq(
      from = 0, to = maxval_colorbar,
      length.out = 11
    ),
    color_fun = clim.palette("yellowred"),
    toptitle = title,
    height = 6,
    fileout = filepng
  )
  filencdf <- paste0(
    work_dir, "/", "capacity_factor_", model_names[i], "_",
    start_year, "-", end_year, ".nc"
  )
  dimlon <- ncdim_def(
    name = "lon", units = "degrees_east",
    vals = as.vector(lon), longname = "longitude"
  )
  dimlat <- ncdim_def(
    name = "lat", units = "degrees_north",
    vals = as.vector(lat), longname = "latitude"
  )
  dimtime <- ncdim_def(
    name = "season", units = "season",
    vals = start_year:end_year,
    longname = "season of the year: DJF, MAM, JJA, SON"
  )
  dimcurve <- ncdim_def(
    name = "curve", units = "name", vals = seq(1, 5, 1),
    longname = "Power curves of considered turbines"
  )
  names(dim(seas_data_cf_all)) <- c("time", "lat", "lon")
  defdata <- ncvar_def(
    name = "CapacityFactor", units = "%",
    dim = list(
      dimtime,
      lat = dimlat,
      lon = dimlon
    ),
    longname = paste(
      "Capacity Factor of PV",
      "based on rsds and tas"
    )
  )
  file <- nc_create(filencdf, list(defdata))
  ncvar_put(file, defdata, seas_data_cf_all)
  nc_close(file)

  # Set provenance for output files
  xprov <- list(
    ancestors = list(
      fullpath_filenames[i]
    ),
    authors = list(
      "cionni_irene"
    ),
    projects = list("crescendo"),
    caption = title,
    statistics = list("other"),
    realms = list("atmos"),
    themes = list("phys")
  )
  provenance[[filepng]] <- xprov
  provenance[[filencdf]] <- xprov
}

# Write provenance to file
write_yaml(provenance, provenance_file)
