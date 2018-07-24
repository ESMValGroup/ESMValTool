library(yaml)

## Insurance products
args <- commandArgs(trailingOnly = TRUE)
#params <- yaml::read_yaml(args[1])
params <- read_yaml(args[1])

plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)
print(plot_dir)
#FOR THE FIRST METADATA.yml
input_files_tasmax <- yaml::read_yaml(params$input_files[1])
model_names <- input_files_tasmax[[1]]$dataset
#print("AJH")
#print(str(input_files_tasmax))
#print(paste("MODEL:", model_names))
var_names_tmax <- input_files_tasmax[[1]]$short_name
#print(paste("var_names:", var_names_tmax))
experiment <- lapply(input_files_tasmax, function(x){x$exp})
#c(input_files_tasmax[[1]]$exp, input_files_tasmax[[2]]$exp)
filename_tasmax <- lapply(input_files_tasmax, function(x){x$filename})
#c(input_files_tasmax[[1]]$filename, input_files_tasmax[[2]]$filename)
#print(paste("EXP:", experiment))
#    print(paste("FIL:", filename_tasmax))


#lat.max <- input_files_tasmax[[1]]$lat_max
#lat.min <- input_files_tasmax[[1]]$lat_min
#lon.max <- input_files_tasmax[[1]]$lon_max
#lon.min <- input_files_tasmax[[1]]$lon_min
#print('lon')
#print(names(input_files_tasmax[[1]]))
#print(c(lat.max,lat.min, lon.max,lon.min))
#FOR THE SECOND METADA.yml
#print('tmin')
input_files_tasmin <- yaml::read_yaml(params$input_files[2])
var_names_tmin <- input_files_tasmin[[1]]$short_name
#print(paste("var_names:", var_names_tmin))
filename_tasmin <- lapply(input_files_tasmin, function(x){x$filename})
 #   print(paste("FIL:", filename_tasmin))

#print(params$input_files)
#input_files_per_var <- yaml::read_yaml(params$input_files[1])
#print(str(input_files_per_var))
#var_names <- names(input_files_per_var)
#print(paste("var_names:", var_names))
#input_files <- lapply(var_names, function(x) names(input_files_per_var[[x]]))
#print(paste("input_files:", input_files))
#names(input_files) <- var_names

#model_names <- lapply(input_files_per_var, function(x) x$dataset)
#model_names <- unlist(unname(model_names))
#print(paste("model_names:", model_names))

## Do not print warnings
#options(warn=-1)


#Var considered
#var0 <- var_names[1]

#Region considered
#lat.max <- params$lat_max
#lat.min <- params$lat_min
#lon.max <- params$lon_max
#lon.min <- params$lon_min

#print(names(params))
#print(str(lat.max))
#Start and end periods for the historical and projection periods
#start_historical <- as.POSIXct(params$start_historical)
#end_historical <- as.POSIXct(params$end_historical)
#start_projection <- as.POSIXct(params$start_projection)
#end_projection <- as.POSIXct(params$end_projection)

#experiment <- lapply(input_files_per_var, function(x) x$exp)
#experiment <- unlist(unname(experiment))

reference_files <- which(experiment == "historical")
#print(paste("refe",reference_files))
projection_files <- which(experiment != "historical")
#print(paste("rproj",projection_files))

start_historical <- input_files_tasmax[[reference_files]]$start_year
end_historical <- input_files_tasmax[[reference_files]]$end_year
start_projection <- input_files_tasmax[[projection_files[1]]]$start_year
end_projection <- input_files_tasmax[[projection_files[1]]]$end_year
#print(c(start_historical,end_historical,start_projection,end_projection))

#start_reference <- lapply(input_files_per_var, function(x) x$start_year)
#start_reference <- c(unlist(unname(start_reference))[reference_files])[1]
#end_reference <- lapply(input_files_per_var, function(x) x$end_year)
#end_reference <- c(unlist(unname(end_reference))[reference_files])[1]

#start_projection <- lapply(input_files_per_var, function(x) x$start_year)
#start_projection <- c(unlist(unname(start_projection))[projection_files])[1]
#end_projection <- lapply(input_files_per_var, function(x) x$end_year)
#end_projection <- c(unlist(unname(end_projection))[projection_files])[1]


#Regime parameters
metric <- params$metric
#print(paste("met: ", metric))
rcp8.5 <- params$rcp8.5
rcp2.6 <- params$rcp2.6
rcp_scenario <- c(rcp8.5, rcp2.6)

library(s2dverification)
library(startR) #, lib.loc = '/home/Earth/nmanuben/tmp/startR_mod/startR.Rcheck')
library(multiApply)
library(devtools)
library(climdex.pcic)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Climdex.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/Threshold.R')
library(parallel)
#source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Climdex/R/DTR.R')
source('https://earth.bsc.es/gitlab/nperez/magic.bsc/raw/master/R/DTR_Indicator.R')
source('https://earth.bsc.es/gitlab/nperez/magic.bsc/raw/master/R/DTR_Ref.R')
#print(str(input_files_per_var))
#var0 <- lapply(input_files_per_var, function(x) x$short_name)
#print(var0)
#fullpath_filenames <- names(var0)
#print(fullpath_filenames)

fullpath_filenames_historical_tasmax <- filename_tasmax[[reference_files]]
historical_tasmax <- Start(model = fullpath_filenames_historical_tasmax,
                                var = "tasmax",
                                var_var = 'var_names',
                                #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
   time = 'all',
   #                            time = values(list(as.POSIXct(paste(start_historical, "01", "01", sep = "-")),
    #                                              as.POSIXct(paste(end_historical, "12", "31", sep = "-")))),
#                                time_tolerance = as.difftime(0, units = 'days'),
                              #  lat = values(list(lat.min, lat.max)),
                               # lon = values(list(lon.min, lon.max)),
    lon='all', lat='all',
                                lon_var = 'lon',
                                lon_reorder = CircularSort(0, 360),
                                return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                                retrieve = TRUE)

fullpath_filenames_historical_tasmin <- filename_tasmin[[reference_files]]
#print(fullpath_filenames_historical_tasmin)
historical_tasmin <- Start(model = fullpath_filenames_historical_tasmin,
                                var = "tasmin",
                                var_var = 'var_names',
                                #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
time = 'all',
#  time = values(list(as.POSIXct(paste(start_historical, "01", "01", sep = "-")),
#                                                   as.POSIXct(paste(end_historical, "12", "31", sep = "-")))),
#                                time_tolerance = as.difftime(0, units = 'days'),
                          #      lat = values(list(lat.min, lat.max)),
                         #       lon = values(list(lon.min, lon.max)),
    lon='all', lat='all',
                                lon_var = 'lon',
                                lon_reorder = CircularSort(0, 360),
                                return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                                retrieve = TRUE)
print(dim(historical_tasmax))
print(dim(historical_tasmin))


dtr_base <- DTR_Ref(tmax = historical_tasmax, tmin = historical_tasmin, by.seasons = TRUE, ncores = NULL)
#print("EO")

#print(class(projection_files))
for (i in 1 : length(projection_files)){
   fullpath_filenames_projection_tasmax <- filename_tasmax[[projection_files[i]]]
  rcp_tasmax <- Start(model = fullpath_filenames_projection_tasmax,
                      var = 'tasmax',
                      var_var = 'var_names',
    time = 'all',
                      #time = values(list(as.POSIXct(start_projection),
                      #                   as.POSIXct(end_projection))),
                      #time_tolerance = as.difftime(0, units = 'days'),
                      #lat = values(list(lat.min, lat.max)),
                      #lon = values(list(lon.min, lon.max)),
    lon='all', lat='all',
                                lon_var = 'lon',
                      lon_reorder = CircularSort(0, 360),
                      return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                      retrieve = TRUE)

fullpath_filenames_projection_tasmin <- filename_tasmin[[projection_files[i]]]
  rcp_tasmin <- Start(model = fullpath_filenames_projection_tasmin,
                      var = 'tasmin',
                      var_var = 'var_names',
    time = 'all',
                      #time = values(list(as.POSIXct(start_projection),
                      #                   as.POSIXct(end_projection))),
                      #time_tolerance = as.difftime(0, units = 'days'),
                      #lat = values(list(lat.min, lat.max)),
                      #lon = values(list(lon.min, lon.max)),
     lon='all', lat='all',
                                lon_var = 'lon',
                      lon_reorder = CircularSort(0, 360),
                      return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                      retrieve = TRUE)


  dtr_indicator <- DTR_Indicator(rcp_tasmax, rcp_tasmin, ref = dtr_base, by.seasons = TRUE, ncores = NULL)

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

    print(levels(season))
print(dim(data))
  for (j in 1 : length(levels(season))) {
    data_subset <- data[,j,,,,]
      print(dim(data_subset))
    data_subset <- aperm(data_subset, c(2,3,1))
    names(dim(data_subset)) <- c("lon", "lat", "time")
      print(dim(data_subset))
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
                  paste0(plot_dir, "/", "dtr_indicator","_",levels(season)[j] ,"_" ,model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".nc"))
    title <- paste0("Number of days in ",levels(season)[j]  , " exceeding the mean diurnal temperature range by 5 degrees ", " ", substr(start_projection, 1, 4), "-",
                    substr(end_projection, 1, 4))
    print(time)
    breaks <- seq(0, max(data_subset),5)
    PlotEquiMap(Mean1Dim(data_subset, 3), lon = lon, lat = lat, filled.continents = FALSE,
                units = "Days", title_scale = 0.5,
                toptitle = title, brks = breaks, color_fun = clim.palette("yellowred"),
                fileout =  paste0(plot_dir, "/", "dtr_indicator_",levels(season)[j], "_",model_names,"_", rcp_scenario[i], "_", start_projection, "_", end_projection, ".pdf"))

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



