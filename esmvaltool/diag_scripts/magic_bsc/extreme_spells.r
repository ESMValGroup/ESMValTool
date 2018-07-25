library(yaml)
library(s2dverification)
library(startR)
library(multiApply)
library(devtools)
library(climdex.pcic)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Threshold.R')
library(parallel)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Heatwaves.R')



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
model_names <- lapply(input_files_per_var, function(x) x$dataset)
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
op <- as.character(params$operator)
print(paste("OP:", op))
qtile = params$quantile
print(paste("QT:", qtile))
spell_length <- params$min_duration

reference_filenames <-  fullpath_filenames[reference_files]
print(reference_filenames)
print(var0)
historical_data  <- Start(model = reference_filenames,
                         var = var0,
                         var_var = 'var_names',
                         time = 'all',
                         lat = 'all',
                         lon = 'all',
                         lon_var = 'lon',
  # lon_reorder = CircularSort(0, 360),
                         return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                         retrieve = TRUE)

    # ------------------------------
# Provisional solution to error in dimension order:
 lon <- attr(historical_data, "Variables")$dat1$lon
 lat <- attr(historical_data, "Variables")$dat1$lat
 time <- attr(historical_data, "Variables")$dat1$time
    historical_data <- as.vector(historical_data)
    dim(historical_data) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    historical_data <- aperm(historical_data, c(1,2,5,4,3))
     attr(historical_data, "Variables")$dat1$time <- time
    print(dim(historical_data))
# ------------------------------





base_range <- c(as.numeric(substr(start_reference, 1, 4)), as.numeric(substr(end_reference, 1, 4)))
threshold <- Threshold(historical_data, base_range = base_range, qtiles = qtile, ncores = NULL)



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
    dim(projection_data) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    projection_data <- aperm(projection_data, c(1,2,5,4,3))
     attr(projection_data, "Variables")$dat1$time <- time
    print(dim(projection_data))
# ------------------------------


 heatwave <- Heatwave(projection_data, threshold, op = op, spell_length = spell_length, by.seasons = TRUE, ncores = NULL)


  if (var0 == "tasmax") {
    ## Select summer season

    heatwave_season <- heatwave$result[seq(2, dim(heatwave$result)[1] - 2, by = 4), 1, 1, , ]
    years <-  as.numeric(substr(start_projection, 1, 4)) : as.numeric(substr(end_projection, 1, 4))
  } else {
    ##Select winter season
    heatwave_season <- heatwave$result[seq(1, dim(heatwave$result)[1] - 2, by = 4), 1, 1, , ]
    years <-  as.numeric(substr(start_projection, 1, 4)) : as.numeric(substr(end_projection, 1, 4))
  }


    #print(lon)
  #  print(lon)
  #lon[lon > 180] <- lon[lon > 180] - 360
  #lon_order <- sort(lon, index.return = TRUE)

  data <- heatwave_season
  #  print(dim(data))
  data <- aperm(data, c(3,2,1))
  names(dim(data)) <- c("lon", "lat", "time")

  #data <- Subset(data, "lon", lon_order$ix)
    #print(length(lon))


  #lon <- lon_order$x
  attributes(lon) <- NULL
  attributes(lat) <- NULL
  dim(lon) <-  c(lon = length(lon))
  dim(lat) <- c(lat = length(lat))
  metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec='double'))))
  attr(data, 'variables') <- metadata
  time <- years
  attributes(time) <- NULL
  dim(time) <- c(time = length(time))
  metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'years since 0-0-0 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
  attr(time, "variables") <- metadata
  ArrayToNetCDF(list(metric= data, lat = lat, lon = lon, time = time),
                paste0(plot_dir, "/" ,var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))
  brks <- seq(0, 40, 4)
  title <- paste0("Days JJA ", var0, " ", substr(start_projection, 1, 4), "-",
                  substr(end_projection, 1, 4), " ",op, " the ", substr(as.character(qtile), 3, 4),
                  "th quantile for ",substr(start_reference, 1, 4), "-", substr(end_reference, 1, 4),
                  " (",rcp_scenario[i], ")")
#print(lon)
    #print(lat)
   dat <- Mean1Dim(data,3)

  #  PlotEquiMap(dat, lat = 12:1, lon = 23:1, fill = FALSE, fileout = paste0(plot_dir,"/New_", var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i],
    #                            "_", start_projection, "_", end_projection, ".pdf"))

  PlotEquiMap(Mean1Dim(data, 3), lat = lat, lon = lon, fill = FALSE,
              brks = brks, color_fun = clim.palette("yellowred"),
              units = paste0("Days" ), toptitle = title,
              fileout = paste0(plot_dir,"/", var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i],
                               "_", start_projection, "_", end_projection, ".pdf"),
              title_scale = 0.5)

}


## Plots
#source("https://earth.bsc.es/gitlab/rserrano/colorbar/raw/master/int_breaks.R")
#new_breaks <- int_breaks(data[1, , ], method = "equal", zero_centered = FALSE, n=10)







