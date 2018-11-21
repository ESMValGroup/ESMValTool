####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4

Sys.setenv(TAR = "/bin/tar")
library(s2dverification)
#library(startR)
#library(startR)

library(multiApply)
#library(ggplot2)
library(yaml)
library(ncdf4)

library(parallel)
##Until integrated into current version of s2dverification
#library(magic.bsc, lib.loc = "/home/Earth/nperez/git/magic.bsc.Rcheck/")
library(ClimProjDiags)
#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)
input_files_per_var <- yaml::read_yaml(params$input_files)
var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x) x$model)
model_names <- unname(model_names)
var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]
experiment <- lapply(input_files_per_var, function(x) x$exp)
experiment <- unlist(unname(experiment))


reference_files <- which(unname(experiment) == "historical")
projection_files <- which(unname(experiment) != "historical")

rcp_scenario <- unique(experiment[projection_files])
model_names <-  lapply(input_files_per_var, function(x) x$dataset)
model_names <- unlist(unname(model_names))[projection_files]

start_reference <- lapply(input_files_per_var, function(x) x$start_year)
start_reference <- c(unlist(unname(start_reference))[reference_files])[1]
end_reference <- lapply(input_files_per_var, function(x) x$end_year)
end_reference <- c(unlist(unname(end_reference))[reference_files])[1]

start_projection <- lapply(input_files_per_var, function(x) x$start_year)
start_projection <- c(unlist(unname(start_projection))[projection_files])[1]
end_projection <- lapply(input_files_per_var, function(x) x$end_year)
end_projection <- c(unlist(unname(end_projection))[projection_files])[1]


## Do not print warnings
#options(warn=-1)


# Which metric to be computed
metric <- params$metric
var0 <- unlist(var0)
projection <- "NULL"
reference_filenames <-  fullpath_filenames[reference_files]
hist_nc <- nc_open(reference_filenames)
historical_data <- ncvar_get(hist_nc, var0)

names(dim(historical_data)) <- rev(names(hist_nc$dim))[-1]
lat <- ncvar_get(hist_nc,"lat")
lon <- ncvar_get(hist_nc,"lon")
units <- ncatt_get(hist_nc, var0, "units")$value
calendar <- ncatt_get(hist_nc, "time", "calendar")$value
long_names <-  ncatt_get(hist_nc,var0,"long_name")$value
time <-  ncvar_get(hist_nc,"time")
start_date <- as.POSIXct(substr(ncatt_get(hist_nc, "time", "units")$value,11, 29 ))

time <- as.Date(time, origin = start_date, calendar = calendar)
#projection <- attr(historical_data, "Variables")$common$tas$coordinates

# nolint start
#hist_names <- names(dim(historical_data))
#jpeg(paste0(plot_dir, "/plot1.jpg"))
#PlotEquiMap(historical_data[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
# ------------------------------
# Provisional solution to error in dimension order:
# nolint end
if ((end_reference-start_reference + 1) * 12 == length(time)) {
  time <- seq(
    as.Date(
      paste(start_reference, "01", "01", sep = "-"),
      format = "%Y-%m-%d"
    ),
    as.Date(
      paste(end_reference, "12", "01", sep = "-"),
      format = "%Y-%m-%d"
    ),
    "month"
  )
}
historical_data <- as.vector(historical_data)
dim(historical_data) <- c(
  model = 1,
  var = 1,
  lon = length(lon),
  lat = length(lat),
  time = length(time)
)
historical_data <- aperm(historical_data, c(1, 2, 5, 3, 4))
# nolint start
# ------------------------------
#jpeg(paste0(plot_dir, "/plot2.jpg"))
#PlotEquiMap(historical_data[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
# nolint end

names(dim(historical_data)) <- c("model", "var", "time", "lon", "lat")
time_dimension <- which(names(dim(historical_data)) == "time")

attributes(lon) <- NULL
attributes(lat) <- NULL

dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))
model_dim <- which(names(dim(historical_data)) == "model")
###Compute the quantiles and standard deviation for the historical period.

if (var0 == "tasmin") {
  #metric <- "t10p"
  quantile <- 0.1
} else if (var0 == "tasmax") {
  #metric <- "t90p"
  quantile <- 0.9
} else if (var0 == "sfcWind") {
  #metric <- "Wx"
  historical_data <- 0.5 * 1.23 * (historical_data ** 3)
  quantile <- 0.9
} else if (var0 == "pr") {
  #metric <- c("cdd", "rx5day")
  historical_data <- historical_data * 60 * 60 * 24

}
attr(historical_data, "Variables")$dat1$time <- time
#names(dim(historical_data)) <- hist_names
base_sd <- base_sd_historical <- base_mean <- list()
for (m in 1 : length(metric)) {
  if (var0 != "pr") {
    thresholds <- Threshold(
      historical_data,
      calendar = calendar,
      qtiles = quantile,
      ncores = detectCores() -1
    )
  str(thresholds)
    base_index <- Climdex(
      data = historical_data,
      calendar = calendar,
      metric = metric[m],
      threshold = thresholds,
      ncores = detectCores() - 1
    )
  } else {
    base_index <- Climdex(
      data = historical_data,
      calendar = calendar,
      metric = metric[m],
      ncores = detectCores() - 1
    )
  }
  base_sd[[m]] <- Apply(
    list(base_index$result),
    target_dims = list(c(1)),
     "sd"
  )$output1
  base_sd_historical[[m]] <- InsertDim(
    base_sd[[m]], 1, dim(base_index$result)[1]
  )

  if (var0 != "pr") {
    base_mean[[m]] <- 10
    base_mean_historical <- 10
  } else {
    base_mean[[m]] <- Apply(
      list(base_index$result),
      target_dims = list(c(1)),
       "mean"
    )$output1
    base_mean_historical <- InsertDim(
      base_mean[[m]], 1, dim(base_index$result)[1]
    )
  }
}
# Compute the time series of the relevant index, using the quantiles
# and standard deviation from the index
projection_filenames <-  fullpath_filenames[projection_files]

for (i in 1 : length(projection_filenames)) {
  proj_nc <- nc_open(projection_filenames[i])
  projection_data <- ncvar_get(proj_nc, var0)
  time <-  ncvar_get(proj_nc,"time")
  start_date <- as.POSIXct(substr(ncatt_get(proj_nc, "time", "units")$value,11, 29 ))
  calendar <- ncatt_get(hist_nc, "time", "calendar")$value
  time <- as.Date(time, origin = start_date, calendar = calendar)

  #proj_names <- names(dim(projection_data))
  # ------------------------------
  #jpeg(paste0(plot_dir, "/plot3.jpg"))
  #PlotEquiMap(projection_data[1,1,1,,], lon = lon, lat = lat, filled = F)
  #dev.off()
      # ------------------------------
  # Provisional solution to error in dimension order:
  if ((end_projection-start_projection + 1) * 12 == length(time)) {
    time <- seq(
      as.Date(
        paste(start_projection, "01", "01", sep = "-"),
        format = "%Y-%m-%d"
      ),
      as.Date(
        paste(end_projection, "12", "01", sep = "-"),
        format = "%Y-%m-%d"
      ),
      "month"
    )
  }
  projection_data <- as.vector(projection_data)
  dim(projection_data) <- c(
    model = 1,
    var = 1,
    lon = length(lon),
    lat = length(lat),
    time = length(time)
  )
  projection_data <- aperm(projection_data, c(1,2,5,3,4))
  attr(projection_data, "Variables")$dat1$time <- time
  names(dim(projection_data)) <- c("model", "var", "time", "lon", "lat")
  # nolint start
  # ------------------------------
  #jpeg(paste0(plot_dir, "/plot4.jpg"))
  #PlotEquiMap(projection_data[1,1,1,,], lon = lon, lat = lat, filled = F)
  #dev.off()
  # nolint end

  if (var0 == "pr") {
    projection_data <- projection_data * 60 * 60 * 24
  } else if (var0 == "sfcWind") {
    projection_data <- 0.5 * 1.23 * (projection_data ** 3)
  }
    #names(dim(projection_data)) <- proj_names
  for (m in 1 : length(metric)) {

    if (var0 != "pr") {
      projection_index <- Climdex(data = projection_data, metric = metric[m], calendar = calendar,
                                  threshold = thresholds, ncores = detectCores() - 1)
      projection_mean <- 10
    } else {
      projection_index <- Climdex(data = projection_data, metric = metric[m], calendar = calendar,
                                  ncores = detectCores() - 1)
      projection_mean <- InsertDim(base_mean[[m]], 1, dim(projection_index$result)[1])
    }

    base_sd_proj <- InsertDim(base_sd[[m]], 1, dim(projection_index$result)[1])
    projection_index_standardized <- (projection_index$result - projection_mean) / base_sd_proj
    #model_dim <- which(proj_names == "model")
    for (mod in 1 : dim(projection_data)[model_dim]) {
      data <- drop(Mean1Dim(projection_index_standardized, 1))
      print(paste(
        "Attribute projection from climatological data is saved and,",
        "if it's correct, it can be added to the final output:",
        projection
      ))
      dimlon <- ncdim_def(
        name = "lon",
        units = "degrees_east",
        vals = as.vector(lon),
        longname = "longitude"
      )
      dimlat <- ncdim_def(
        name = "lat",
        units = "degrees_north",
        vals = as.vector(lat),
        longname = "latitude"
      )
      defdata <- ncvar_def(
        name = "data",
        units = units,
        dim = list(lat = dimlat, lon = dimlon),
        longname = paste("Mean", metric[m], long_names)
      )
      file <- nc_create(
        paste0(
          plot_dir, "/", var0, "_", metric[m], "_risk_insurance_index_",
          model_names, "_", start_projection, "_", end_projection, "_",
          start_reference, "_", end_reference, ".nc"
        ),
        list(defdata)
      )
      ncvar_put(file, defdata, data)
      nc_close(file)

      title <- paste0(
        "Index for  ", metric[m], " ", substr(start_projection, 1, 4), "-",
        substr(end_projection, 1, 4), " ", " (", rcp_scenario[i],
        " ", model_names ,")")

      breaks <- seq(-1 * ceiling(max(abs(data))), ceiling(max(abs(data))), 2 * ceiling(max(abs(data))) / 16)
      PlotEquiMap(
        data,
        lon = lon,
        lat = lat,
        filled.continents = FALSE,
        toptitle = title,
        brks = breaks,
        fileout = paste0(
          plot_dir, "/", metric[m], "_",model_names[mod],"_", rcp_scenario[i],
          "_", start_projection, "_", end_projection, ".png"
        )
      )
    }
  }

}
