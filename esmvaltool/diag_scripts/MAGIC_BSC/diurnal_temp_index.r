## Insurance products
args <- commandArgs(trailingOnly = TRUE)
params <- yaml::read_yaml(args[1])

plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)

input_files_per_var <- yaml::read_yaml(params$input_files)
var_names <- names(input_files_per_var)
input_files <- lapply(var_names, function(x) names(input_files_per_var[[x]]))
names(input_files) <- var_names
model_names <- lapply(input_files_per_var, function(x) unname(sapply(x, '[[', 'model')))

## Do not print warnings
#options(warn=-1)


#Var considered
var0 <- var_names[1]

#Region considered
lat.max <- params$lat_max
lat.min <- params$lat_min
lon.max <- params$lon_max
lon.min <- params$lon_min


#Start and end periods for the historical and projection periods
start_historical <- as.POSIXct(params$start_historical)
end_historical <- as.POSIXct(params$end_historical)
start_projection <- as.POSIXct(params$start_projection)
end_projection <- as.POSIXct(params$end_projection)

#Regime parameters
metric <- params$metric
rcp8.5 <- params$rcp8.5
rcp2.6 <- params$rcp2.6
rcp_scenario <- c(rcp8.5, rcp2.6)

library(s2dverification)
library(startR)
library(multiApply)
library(devtools)
library(climdex.pcic)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Climdex.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Threshold.R')
library(parallel)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/DTR.R')


historical_tasmax <- Start(model = fullpath_filenames_historical_tasmax,
                                var = "tasmax",
                                var_var = 'var_names',
                                #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
                                time = values(list(as.POSIXct(start_historical), 
                                                   as.POSIXct(end_historical))),
                                time_tolerance = as.difftime(0, units = 'days'), 
                                lat = values(list(lat.min, lat.max)),
                                lon = values(list(lon.min, lon.max)),
                                lon_var = 'lon',
                                lon_reorder = CircularSort(0, 360),
                                return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                                retrieve = TRUE)

historical_tasmin <- Start(model = fullpath_filenames_historical_tasmin,
                                var = "tasmin",
                                var_var = 'var_names',
                                #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
                                time = values(list(as.POSIXct(start_historical), 
                                                   as.POSIXct(end_historical))),
                                time_tolerance = as.difftime(0, units = 'days'), 
                                lat = values(list(lat.min, lat.max)),
                                lon = values(list(lon.min, lon.max)),
                                lon_var = 'lon',
                                lon_reorder = CircularSort(0, 360),
                                return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                                retrieve = TRUE)

dtr_base <- DTR_ref(tmax = historical_tasmax, tmin = historical_tasmin, by.seasons = TRUE, ncores = NULL)

for (i in 1 : length(fullpath_filenames_projection_tasmax)){
  rcp_tasmax <- Start(model = fullpath_filenames_projection_tasmax[[i]],
                      var = 'tasmax',
                      var_var = 'var_names',
                      time = values(list(as.POSIXct(start_projection), 
                                         as.POSIXct(end_projection))),
                      time_tolerance = as.difftime(0, units = 'days'), 
                      lat = values(list(lat.min, lat.max)),
                      lon = values(list(lon.min, lon.max)),
                      lon_var = 'lon',
                      lon_reorder = CircularSort(0, 360),
                      return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                      retrieve = TRUE)
  
  
  rcp_tasmin <- Start(model = fullpath_filenames_projection_tasmin[[i]],
                      var = 'tasmin',
                      var_var = 'var_names',
                      time = values(list(as.POSIXct(start_projection), 
                                         as.POSIXct(end_projection))),
                      time_tolerance = as.difftime(0, units = 'days'), 
                      lat = values(list(lat.min, lat.max)),
                      lon = values(list(lon.min, lon.max)),
                      lon_var = 'lon',
                      lon_reorder = CircularSort(0, 360),
                      return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                      retrieve = TRUE)
  dtr_indicator <- DTR_indicator(rcp_tasmax, rcp_tasmin, ref = dtr_base, by.seasons = TRUE, ncores = NULL)
  
  lat <- attr(rcp_tasmax, "Variables")$dat1$lat
  lon <- attr(rcp_tasmax, "Variables")$dat1$lon
  lon[lon > 180] <- lon[lon > 180] - 360
  lon_order <- sort(lon, index.return = TRUE)
  dtr_indicator$indicator <- Subset(dtr_indicator$indicator, "lon", lon_order$ix)
  lon <- lon_order$x
  season <-   as.factor(dtr_indicator$season)# dtr_indicator$season
  year <- dtr_indicator$year
  attributes(lon) <- NULL
  attributes(lat) <- NULL
  dim(lon) <-  c(lon = length(lon))
  dim(lat) <- c(lat = length(lat))
  data <- dtr_indicator$indicator
  
  for (j in 1 : length(levels(season))) {
    data_subset <- data[,j,,,,]
    data_subset <- aperm(data_subset, c(3,2,1))
    names(dim(data_subset)) <- c("lon", "lat", "time")
    metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec='double'))))
    attr(data_subset, 'variables') <- metadata
    day <- "01"
    if (levels(season)[j] == "DJF") {
      month <- "-01-"
    } else if (levels(season)[j] == "MAM") {
      month <- "-03-"
    } else if (levels(season)[j] == "JJA") {
      month <- "-06-"
    } else {
      month <- "-09-"
    }
    time <- as.POSIXct(paste0(year, month, day))
    time <- julian(time, origin = as.POSIXct("1970-01-01"))
    
    attributes(time) <- NULL
    dim(time) <- c(time = length(time))
    metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
    attr(time, "variables") <- metadata
    ArrayToNetCDF(list(metric= data_subset, lat = lat, lon = lon, time = time), 
                  paste0("dtr_indicator","_",levels(season)[j] ,"_" ,model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))
    title <- paste0("Number of days in ",levels(season)[j]  , " exceeding the mean diurnal temperature range by 5 degrees ", " ", substr(start_projection, 1, 4), "-", 
                    substr(end_projection, 1, 4)) 
    print(time)
    breaks <- seq(0, max(data_subset),5)
    PlotEquiMap(Mean1Dim(data_subset, 3), lon = lon, lat = lat, filled.continents = FALSE,
                units = "Days", title_scale = 0.5, 
                toptitle = title, brks = breaks, color_fun = clim.palette("yellowred"),
                fileout =  paste0("dtr_indicator_",levels(season)[j], "_",model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".pdf"))
    
  }
  
  
  
 
}


### Plots 

### SON
#dtr_diff <- Mean1Dim(dtr_indicator$indicator[,3,1,1,,], 1) - Mean1Dim(dtr_indicator_reference$indicator[,3,1,1,,], 1)

#PlotLayout(PlotEquiMap,c(1,2),lon=lon,lat=lat,var=aperm(dtr_all_seasons,c(3,2,1)),
#           titles=c('DJF', 'MAM', 'JJA', 'SON'), toptitle = "Change in DTR indicator (2020-2040) - (1960-1990)",
#           filled.continents=F, units = "Days",
#           axelab=F,draw_separators = T,subsampleg = 1,brks=seq(-16,16,by=2),
#           bar_extra_labels = c(2,0,0,0),fileout='rcp85.png')



