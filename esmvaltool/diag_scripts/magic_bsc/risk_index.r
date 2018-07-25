####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4

Sys.setenv(TAR = '/bin/tar')
library(s2dverification)
library(startR)#, lib.loc='/home/Earth/ahunter/R/x86_64-unknown-linux-gnu-library/3.2/')
#library(startR, lib.loc = '/home/Earth/nmanuben/tmp/startR_mod/startR.Rcheck')
library(multiApply)
library(ggplot2)
library(yaml)
library(climdex.pcic)
library(devtools)
library(parallel)

##Until integrated into current version of s2dverification
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Climdex.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Threshold.R')


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

print(start_reference)
print(end_reference)

## Do not print warnings
#options(warn=-1)


# Which metric to be computed
metric <- params$metric
detrend <- params$detrend

reference_filenames <-  fullpath_filenames[reference_files]

historical_data <- Start(model = reference_filenames,
                         var = var0,
                         var_var = 'var_names',
                         time = 'all',
                         lat = 'all',
                         lon = 'all',
                         lon_var = 'lon',
    lon_reorder = CircularSort(0, 360),
                         return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                         retrieve = TRUE)
# ------------------------------
# Provisional solution to error in dimension order:
 lon <- attr(historical_data, "Variables")$dat1$lon
 lat <- attr(historical_data, "Variables")$dat1$lat
 time <- attr(historical_data, "Variables")$dat1$time
    historical_data <- as.vector(historical_data)
    dim(historical_data) <- c(model = 1, var = 1, lat = length(lat), lon = length(lon), time = length(time))
historical_data <- aperm(historical_data, c(1,2,5,3,4))
    attr(historical_data, "Variables")$dat1$time <- time
# ------------------------------

#print("real")
#print(historical_data[1,1,1:13,1:2,1:2])
#x <- as.vector(historical_data[1,1,,,])
#print("vec")
#print(x[1:9])
#dim(x) <- c(lat = 12, lon = 12, time = 7300)
#print("tras")
#print(x[1:4, 1:5 ,1])
#PlotEquiMap(x[,,2], lon=lon, lat=lat, filled =FALSE, fileout = paste0(plot_dir, "/Plots.png"))

time_dimension <- which(names(dim(historical_data)) == "time")

#lon[lon > 180] <- lon[lon > 180] - 360
#lon_order <- sort(lon, index.return = TRUE)
#historical_data <- Subset(historical_data, "lon", lon_order$ix)
#lon <- lon_order$x

PlotEquiMap(historical_data[1,1,1,,], lon=lon, lat=lat, filled = FALSE, fileout = paste0(plot_dir, "/Plots1.png"))
attributes(lon) <- NULL
attributes(lat) <- NULL
# attributes(years) <- NULL
# dim(years) <- c(length(years))
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))
model_dim <- which(names(dim(historical_data)) == "model")
###Compute the quantiles and standard deviation for the historical period.
if (detrend == TRUE) {
    historical_data <- Trend(historical_data,)
}
if (var0 == "tasmin") {
  metric <- "tx10p"
  quantile <- 0.1
} else if (var0 == "tasmax") {
  metric <- "tx90p"
  quantile <- 0.9
} else if (var0 == "sfcWind") {
  metric <- "Wx"
  historical_data <- 0.5 * 1.23 * (historical_data ** 3)  # Convert to wind power
  quantile <- 0.9
} else if (var0 == "pr") {
  metric <- c("cdd", "rx5day")
  historical_data <- historical_data * 60 * 60 * 24

}

base_sd <- base_sd_historical <- base_mean <- list()
for (m in 1 : length(metric)) {
  base_range <- as.numeric(c(substr(start_reference, 1, 4), substr(end_reference, 1, 4)))

  #Compute the 90th percentile for the historical period
  if (var0 != "pr") {
    thresholds <- Threshold(historical_data, base_range=as.numeric(base_range),
                            qtiles = quantile, ncores = detectCores() -1)
    base_index <- Climdex(data = historical_data, metric = metric[m],
                          quantiles = thresholds, ncores = detectCores() - 1)
  } else {
    base_index <- Climdex(data = historical_data, metric = metric[m], ncores = detectCores() - 1)
  }

  base_sd[[m]] <- Apply(list(base_index$result), target_dims = list(c(1)), AtomicFun = "sd")$output1
  base_sd_historical[[m]] <- InsertDim(base_sd[[m]], 1, dim(base_index$result)[1])

  if (var0 != "pr") {
    base_mean[[m]] <- 10
    base_mean_historical <- 10
  } else {
    base_mean[[m]] <-  Apply(list(base_index$result), target_dims = list(c(1)), AtomicFun = "mean")$output1
    base_mean_historical <- InsertDim(base_mean[[m]], 1, dim(base_index$result)[1])
  }
  historical_index_standardized <- (base_index$result - base_mean_historical) / base_sd_historical[[m]]
 #   print(str(historical_index_standardized))
 #   print(dim(historical_index_standardized))
 #   print(attributes(historical_index_standardized))
  for (mod in 1 : dim(historical_data)[model_dim]) {
    historical_index_standardized <- aperm(historical_index_standardized[,mod,1, ,], c(3,2,1))
    names(dim(historical_index_standardized)) <- c("lon", "lat", "time")
    metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec='double'))))
    names(metadata)[1] <- var0
    attr(historical_index_standardized, 'variables') <- metadata
    day <- "01"
    month <- "-01-"
    time <- as.numeric(base_index$years)
    time <- as.POSIXct(paste0(time, month, day), tz = "CET")
    time <- julian(time, origin = as.POSIXct("1970-01-01"))

    attributes(time) <- NULL
    dim(time) <- c(time = length(time))
    metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
    attr(time, "variables") <- metadata
         print(str(historical_index_standardized))
    print(dim(historical_index_standardized))
    print(attributes(historical_index_standardized))
    variable_list <- list(index = historical_index_standardized, lat = lat, lon = lon, time = time)
    names(variable_list)[1] <- var0

print("AS")

    ArrayToNetCDF(variable_list,
                  paste0(plot_dir, "/", metric[m], "_",model_names[mod],"_", "historical", "_", start_reference, "_", end_reference, ".nc"))
  }
}


#Compute the time series of the relevant index, using the quantiles and standard deviation from the index
projection_filenames <-  fullpath_filenames[projection_files]

for (i in 1 : length(projection_filenames)) {
    projection_data <- Start(model = projection_filenames[i],
                             var = var0,
                             var_var = 'var_names',
                             time = 'all',
                             lat = 'all',
                             lon = 'all',
                             lon_var = 'lon',
                        lon_reorder = CircularSort(0, 360),
                             return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                             retrieve = TRUE)
    # ------------------------------
# Provisional solution to error in dimension order:
 lon <- attr(projection_data, "Variables")$dat1$lon
 lat <- attr(projection_data, "Variables")$dat1$lat
 time <- attr(projection_data, "Variables")$dat1$time
    projection_data <- as.vector(projection_data)
    dim(projection_data) <- c(model = 1, var = 1, lat = length(lat), lon = length(lon), time = length(time))
    projection_data <- aperm(projection_data, c(1,2,5,3,4))
     attr(projection_data, "Variables")$dat1$time <- time
    print(dim(projection_data))
# ------------------------------


    #projection_data <- Subset(projection_data, "lon", lon_order$ix)
    #lon <- lon_order$x
    if (var0 == "pr") {
      projection_data <- projection_data * 60 * 60 * 24
    } else if (var0 == "sfcWind") {
      projection_data <- 0.5 * 1.23 * (projection_data ** 3)
    }

  for (m in 1 : length(metric)) {

    if (var0 != "pr") {
      projection_index <- Climdex(data = projection_data, metric = metric[m],
                                  quantiles = thresholds, ncores = detectCores() - 1)
      projection_mean <- 10
    } else {
      projection_index <- Climdex(data = projection_data, metric = metric[m],
                                  ncores = detectCores() - 1)
      projection_mean <- InsertDim(base_mean[[m]], 1, dim(projection_index$result)[1])

    }

    base_sd_proj <- InsertDim(base_sd[[m]], 1, dim(projection_index$result)[1])
    projection_index_standardized <- (projection_index$result - projection_mean) / base_sd_proj
#print("PROJ")
#print(str(projection_index_standardized))
 #   print(dim(projection_index_standardized))#
      # print(attributes(projection_index_standardized))
    for (mod in 1 : dim(projection_data)[model_dim]) {
        projection_index_standardized <- aperm(projection_index_standardized[,mod,1, ,], c(3,2,1))
      names(dim(projection_index_standardized)) <- c("lon", "lat", "time")

      metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec = 'double'))))
      names(metadata)[1] <- var0
        attr(projection_index_standardized, 'variables') <- metadata

      day <- "01"
      month <- "-01-"
        time <- as.numeric(projection_index$years)
      time <- as.POSIXct(paste0(time, month, day), tz = "CET")
      time <- as.numeric(julian(time, origin = as.POSIXct("1970-01-01")))

      attributes(time) <- NULL
        #attributes(lat) <- NULL
        #attributes(lon) <- NULL

      dim(time) <- c(time = length(time))
      metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
      attr(time, "variables") <- metadata
        names(variable_list)[1] <- var0


      variable_list <- list(index = projection_index_standardized, lat = lat, lon = lon, time = time)

print("dr")

   #   ArrayToNetCDF(variable_list,
    #                paste0(plot_dir, "/", metric[m], "_",model_names[mod],"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))

      title <- paste0("Index for  ", metric[m], " ", substr(start_projection, 1, 4), "-",
                      substr(end_projection, 1, 4), " ",
                      " (",rcp_scenario[i], ")")
      breaks <- -8 : 8
      PlotEquiMap(Mean1Dim(projection_index_standardized, 3), lon = lon, lat = lat, filled.continents = FALSE,
                  toptitle = title, brks = breaks,
                  fileout = paste0(plot_dir, "/", metric[m], "_",model_names[mod],"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".pdf"))
      }
  }

}











