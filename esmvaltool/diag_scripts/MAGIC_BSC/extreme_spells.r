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
model_names <-  lapply(input_files_per_var, function(x) x$model)
model_names <- unlist(unname(model_names))[projection_files]

start_reference <- lapply(input_files_per_var, function(x) x$start_year)
start_reference <- c(unlist(unname(start_reference))[reference_files])[1]
end_reference <- lapply(input_files_per_var, function(x) x$end_year)
end_reference <- c(unlist(unname(end_reference))[reference_files])[1]

start_projection <- lapply(input_files_per_var, function(x) x$start_year)
start_projection <- c(unlist(unname(start_projection))[projection_files])[1]
end_projection <- lapply(input_files_per_var, function(x) x$end_year)
end_projection <- c(unlist(unname(end_projection))[projection_files])[1]

print(start_climatology)
print(end_projection)

## Do not print warnings
#options(warn=-1)


# Which metric to be computed
op <- as.character(params$operator)
qtile = params$quantile
spell_length <- params$min_duration

reference_filenames <-  fullpath_filenames[reference_files]

historical_data <-  <- Start(model = reference_filenames,
                         var = var0,
                         var_var = 'var_names',
                         time = 'all',
                         lat = 'all',
                         lon = 'all',
                         lon_var = 'lon',
                         return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                         retrieve = TRUE)


base_range <- c(as.numeric(substr(start_reference, 1, 4)), as.numeric(substr(end_reference, 1, 4)))
threshold <- Threshold(historical_data, base_range = base_range, qtiles = qtile, ncores = NULL)

projection_filenames <-  projection_filenames[projection_files]

for (i in 1 : length(projection_filenames)) {
    projection_data <- Start(model = projection_filenames[i],
                             var = var0,
                             var_var = 'var_names',
                             time = 'all',
                             lat = 'all',
                             lon = 'all',
                             lon_var = 'lon',
                             return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                             retrieve = TRUE)

 heatwave <- Heatwave(rcp_data, threshold, op = op, spell_length = spell_length, by.seasons = TRUE, ncores = NULL)
  if (var0 == "tasmax") {
    ## Select summer season
    heatwave_season <- heatwave$result[seq(2, dim(heatwave$result)[1] - 2, by = 4), 1, 1, , ]
    years <-  as.numeric(substr(start_projection, 1, 4)) : as.numeric(substr(end_projection, 1, 4))
  } else {
    ##Select winter season
    heatwave_season <- heatwave$result[seq(1, dim(heatwave$result)[1] - 2, by = 4), 1, 1, , ]
    years <-  as.numeric(substr(start_projection, 1, 4)) : as.numeric(substr(end_projection, 1, 4))
  }

  lat <- attr(rcp_data, "Variables")$dat1$lat
  lon <- attr(rcp_data, "Variables")$dat1$lon
  lon[lon > 180] <- lon[lon > 180] - 360
  lon_order <- sort(lon, index.return = TRUE)
  data <- heatwave_season
  data <- aperm(data, c(3,2,1))
  names(dim(data)) <- c("lon", "lat", "time")
  data <- Subset(data, "lon", lon_order$ix)
  lon <- lon_order$x
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
                paste0(var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))
  brks <- seq(0, 40, 4)
  title <- paste0("Days JJA ", var0, " ", substr(start_projection, 1, 4), "-",
                  substr(end_projection, 1, 4), " ",op, " the ", substr(as.character(qtile), 3, 4),
                  "th quantile for ",substr(start_reference, 1, 4), "-", substr(end_reference, 1, 4),
                  " (",rcp_scenario[i], ")")
  PlotEquiMap(Mean1Dim(data, 3), lon = lon, lat = lat, fill = FALSE,
              brks = brks, color_fun = clim.palette("yellowred"),
              units = paste0("Days" ), toptitle = title,
              fileout = paste0(var0, "_extreme_spell_duration", "_",model_names,"_", rcp_scenario[i],
                               "_", start_projection, "_", end_projection, ".pdf"),
              title_scale = 0.5)
}


## Plots
#source("https://earth.bsc.es/gitlab/rserrano/colorbar/raw/master/int_breaks.R")
#new_breaks <- int_breaks(data[1, , ], method = "equal", zero_centered = FALSE, n=10)







