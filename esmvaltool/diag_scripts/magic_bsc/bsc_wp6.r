####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4
#install.packages('yaml')
#install.packages('devtools')
#library(devtools)
#Sys.setenv(TAR = '/bin/tar')
#install_git('https://earth.bsc.es/gitlab/es/startR', branch = 'develop-hotfixes-0.0.2')
#install_git('https://earth.bsc.es/gitlab/es/easyNCDF', branch = 'master')
Sys.setenv(TAR = '/bin/tar')
library(s2dverification)
library(startR, lib.loc='/home/Earth/ahunter/R/x86_64-unknown-linux-gnu-library/3.2/')
library(multiApply)
library(ggplot2)
library(yaml)

##Until integrated into current version of s2dverification
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
source("https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/CombineIndices.R")
source("https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/SelBox.R")



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
rcp_scenario <- experiment

model_names <-  lapply(input_files_per_var, function(x) x$model)
model_names <- unlist(unname(model_names))

start_year <- lapply(input_files_per_var, function(x) x$start_year)
start_year <- c(unlist(unname(start_year)))[1]
end_year <- lapply(input_files_per_var, function(x) x$end_year)
end_year <- c(unlist(unname(end_year)))[1]


## Do not print warnings
#options(warn=-1)

#Parameters for Season() function
monini <- 1
moninf <- params$moninf
monsup <- params$monsup
region <- params$region
mask <- NULL ### How can we pass masks from ESMValTool?
if (!is.null(params$lat_max)) {
    region <- c(params$lon_min, params$lon_max, params$lat_min, params$lat_max)
} else {
  region <- "NULL"
}
running_mean <- params$running_mean
multi_year_average <- params$multi_year_average
weights <- params$weights
### Load data


data <- Start(model = fullpath_filenames,
                        var = var0,
                        var_var = 'var_names',
                        time = 'all',
                        lat = 'all',
                        lon = 'all',
                        lon_var = 'lon',
                        return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                        retrieve = TRUE)

lat <- attr(data, "Variables")$dat1$lat
lon <- attr(data, "Variables")$dat1$lon
attributes(lon) <- NULL
attributes(lat) <- NULL
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))
time_dim <- which(names(dim(data)) == "time")
timestamp <- ""

if (is.null(moninf)) {
  time <- attributes(data)$Variables$dat1$time

} else {
    day <- "01"
    month <- moninf
    if (length(moninf) == 1) {
      month <- paste0(0, month)
    }
    month <- paste0("-", month, "-")
    time <- start_year : end_year
    time <- as.POSIXct(paste0(time, month, day), tz = "CET")
    time <- julian(time, origin = as.POSIXct("1970-01-01"))
}

attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
attr(time, 'variables') <- metadata


if (!is.null(region)) {
  dim_names <- names(dim(data))
  londim <- which(names(dim(data)) == "lon")
  latdim <- which(names(dim(data)) == "lat")
  data <- WeightedMean(data, lon = lon, lat = lat, region = region, mask = NULL)
  names(dim(data)) <- dim_names[-c(londim,latdim)]
  time_dim <- which(names(dim(data)) == "time")
}

if (!is.null(running_mean)) {
  data <- Smoothing(data, runmeanlen = running_mean, numdimt = time_dim)
  timestamp <- paste0(running_mean, "-month-running-mean-")
}

if (!is.null(moninf)) {
  months <- paste0(month.abb[moninf],"-", month.abb[monsup])
  dims <- dim(data)
  dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
  dims <- dims[-time_dim]
  dim(data) <- dims
  names(dim(data))[c(time_dim, time_dim + 1)] <- c("month", "year")
  margins <- list(c(1 : length(dim(data)))[-c(time_dim + 1)])
  data <- Season(data, posdim = time_dim, monini = monini, moninf = moninf,
                                    monsup = monsup)
  margins <- list(c(1 : length(dim(data)))[-c(time_dim + 1)])
  years_dim <- which(names(dim(data)) == "year")
  timestamp <- "Mar-May"
  month_dim <- which(names(dim(data)) == "month")
  data <- adrop(data, month_dim)
  if (multi_year_average == "TRUE") {
    time <- time[1]
    data <- Mean1Dim(data, c(years_dim -1))
  }
}

if (!is.null(weights)) {
  indices_dim <- which(names(dim(data)) == "model")
  indices <- list()
  for (i in 1 : dim(data)[indices_dim]) {
    indices[[i]] <- Subset(data, along = indices_dim, indices = i)
  }
  print(weights)
  if (!is.numeric(weights)) {
    weights <- "NULL"
    data <- CombineIndices(indices, weights = NULL)
  } else {
    data <- CombineIndices(indices, weights = NULL)
  }
}



if (!is.null(region)) {
  data <- data[1,1, ]
  attributes(data) <- NULL
  dim(data) <-  c(time = length(data))
  metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec = 'double'))))
  names(metadata)[1] <- var0
  attr(data, 'variables') <- metadata
  variable_list <- list(variable = data, time = time)
  names(variable_list)[1] <- var0
  ArrayToNetCDF(variable_list,
                paste0(plot_dir, "/", var0, "_",paste(model_names, sep="", collapse="_"), "_", timestamp, "_", rcp_scenario, "_", start_year, "_", end_year, "_", ".nc"))
} else {
  data <- data[1,1, , ,]
  data <- aperm(data, c(3,2,1))
  names(dim(data)) <- c("lon", "lat", "time")
  metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec = 'double'))))
  names(metadata)[1] <- var0
  attr(data, 'variables') <- metadata
  variable_list <- list(variable = data, lat = lat, lon = lon, time = time)
  names(variable_list)[1] <- var0
  ArrayToNetCDF(variable_list,
                paste0(plot_dir, "/", var0, "_",paste(model_names, sep="", collapse="_"), "_", timestamp, "_", rcp_scenario, "_", start_year, "_", end_year, "_", ".nc"))
}





