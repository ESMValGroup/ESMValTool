library(yaml)
library(s2dverification)
library(startR)
library(multiApply)
library(climdex.pcic)
library(magic.bsc, lib.loc = '/home/Earth/nperez/git/magic.bsc.Rcheck')
library(parallel)
library(ncdf4)

## Insurance products
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])

plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)


#FOR THE FIRST METADATA.yml
input_files_tasmax <- yaml::read_yaml(params$input_files[1])
model_names <- input_files_tasmax[[1]]$dataset
var_names_tmax <- input_files_tasmax[[1]]$short_name
experiment <- lapply(input_files_tasmax, function(x){x$exp})
filename_tasmax <- lapply(input_files_tasmax, function(x){x$filename})


input_files_tasmin <- yaml::read_yaml(params$input_files[2])
var_names_tmin <- input_files_tasmin[[1]]$short_name
filename_tasmin <- lapply(input_files_tasmin, function(x){x$filename})


reference_files <- which(experiment == "historical")
projection_files <- which(experiment != "historical")


start_historical <- input_files_tasmax[[reference_files]]$start_year
end_historical <- input_files_tasmax[[reference_files]]$end_year
start_projection <- input_files_tasmax[[projection_files[1]]]$start_year
end_projection <- input_files_tasmax[[projection_files[1]]]$end_year



#Regime parameters
metric <- params$metric
rcp8.5 <- params$rcp8.5
rcp2.6 <- params$rcp2.6
rcp_scenario <- c(rcp8.5, rcp2.6)




fullpath_filenames_historical_tasmax <- filename_tasmax[[reference_files]]
historical_tasmax <- Start(model = fullpath_filenames_historical_tasmax,
                                var = "tasmax",
                                var_var = 'var_names',
                                time = 'all',
                                lon='all', lat='all',
                                lon_var = 'lon',
                                lon_reorder = CircularSort(0, 360),
                                return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                                retrieve = TRUE)

fullpath_filenames_historical_tasmin <- filename_tasmin[[reference_files]]
historical_tasmin <- Start(model = fullpath_filenames_historical_tasmin,
                                var = "tasmin",
                                var_var = 'var_names',
                                time = 'all',
                                lon='all', lat='all',
                                lon_var = 'lon',
                                lon_reorder = CircularSort(0, 360),
                                return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                                retrieve = TRUE)
  lat <- attr(historical_tasmax, "Variables")$dat1$lat
  lon <- attr(historical_tasmax, "Variables")$dat1$lon

#jpeg(paste0(plot_dir, "/plot1tasmax.jpg"))
#PlotEquiMap(historical_tasmax[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
# ------------------------------------------------------------
# Provisional solution to error in dimension order and time values:
 time <- attr(historical_tasmin, "Variables")$dat1$time
 calendar <- attributes(time)$variables$time$calendar
time_his = time
 if ((end_historical-start_historical + 1) * 12 == length(time)) {
     time <-  seq(as.Date(paste(start_historical, '01', '01', sep = "-"), format = "%Y-%m-%d"), as.Date(paste(end_historical, '12', '01', sep = "-"), format = "%Y-%m-%d"), "day")
 }

    historical_tasmin <- as.vector(historical_tasmin)
    historical_tasmax <- as.vector(historical_tasmax)
    dim(historical_tasmin) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    dim(historical_tasmax) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    historical_tasmin <- aperm(historical_tasmin, c(1,2,5,4,3))
    historical_tasmax <- aperm(historical_tasmax, c(1,2,5,4,3))
     attr(historical_tasmin, "Variables")$dat1$time <- time
     attr(historical_tasmax, "Variables")$dat1$time <- time
# ------------------------------------------------------------
#jpeg(paste0(plot_dir, "/plot2tasmax.jpg"))
#PlotEquiMap(historical_tasmax[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()

long_names <- attr(historical_tasmin, "Variables")$common$tas$long_name
projection <- attr(historical_tasmin, "Variables")$common$tas$coordinates
units <- (attr(historical_tasmin,"Variables")$common)[[2]]$units

dtr_base <- DTRRef(tmax = historical_tasmax, tmin = historical_tasmin, by.seasons = TRUE, ncores = NULL)



#print(str(dtr_base))
#jpeg(paste0(plot_dir, "/plotBASE.jpg"))
#PlotEquiMap(dtr_base$dtr.ref[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()


for (i in 1 : length(projection_files)){
    fullpath_filenames_projection_tasmax <- filename_tasmax[[projection_files[i]]]
    rcp_tasmax <- Start(model = fullpath_filenames_projection_tasmax,
                      var = 'tasmax',
                      var_var = 'var_names',
                      time = 'all',
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
                      lon='all', lat='all',
                      lon_var = 'lon',
                      lon_reorder = CircularSort(0, 360),
                      return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                      retrieve = TRUE)
    lat <- attr(rcp_tasmax, "Variables")$dat1$lat
    lon <- attr(rcp_tasmax, "Variables")$dat1$lon

#jpeg(paste0(plot_dir, "/plot3tasmin.jpg"))
#PlotEquiMap(rcp_tasmin[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
# ------------------------------------------------------------
# Provisional solution to error in dimension order and time values:
 time <- attr(rcp_tasmin, "Variables")$dat1$time
 calendar <- attributes(time)$variables$time$calendar
 if ((end_projection-start_projection + 1) * 12 == length(time)) {
     time <-  seq(as.Date(paste(start_projection, '01', '01', sep = "-"), format = "%Y-%m-%d"), as.Date(paste(end_projection, '12', '01', sep = "-"), format = "%Y-%m-%d"), "day")
 }

    rcp_tasmin <- as.vector(rcp_tasmin)
    rcp_tasmax <- as.vector(rcp_tasmax)
    dim(rcp_tasmin) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    dim(rcp_tasmax) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    rcp_tasmin <- aperm(rcp_tasmin, c(1,2,5,4,3))
    rcp_tasmax <- aperm(rcp_tasmax, c(1,2,5,4,3))
     attr(rcp_tasmin, "Variables")$dat1$time <- time
     attr(rcp_tasmax, "Variables")$dat1$time <- time
# ------------------------------------------------------------
#jpeg(paste0(plot_dir, "/plot4tasmin.jpg"))
#PlotEquiMap(rcp_tasmin[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()

    dtr_indicator <- DTRIndicator(rcp_tasmax, rcp_tasmin, ref = dtr_base, by.seasons = TRUE, ncores = NULL)


#jpeg(paste0(plot_dir, "/plotIndicator.jpg"))
# PlotEquiMap(dtr_indicator$indicator[1,1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
}


### Plots

### SON

dtr_rcp <- array(dim = c(4, length(lon), length(lat)))
for (j in 1 : 4){
    dtr_rcp[j, , ] <- Mean1Dim(dtr_indicator$indicator[, j, , , , ], 1)
}
names(dim(dtr_rcp)) <- c("season", "lon", "lat")
PlotLayout(PlotEquiMap, plot_dims = c('lon', 'lat'), var = dtr_rcp,
       lon = lon, lat = lat ,
	   titles = c('DJF', 'MAM', 'JJA', 'SON'),
	   toptitle = paste("Number of days exceeding the DTR in 5 degrees during the period", start_projection, "-", end_projection),
       filled.continents = FALSE, units = "Days",
       axelab = FALSE, draw_separators = TRUE, subsampleg = 1,
	    brks = seq(0, max(dtr_rcp), 2), color_fun = clim.palette("yellowred"),
       bar_extra_labels = c(2, 0, 0, 0), title_scale = 0.7,
	   fileout = paste0(plot_dir, '/rcp85.png'))

print(paste("Attribute projection from climatological data is saved and, if it's correct, it can be added to the final output:", projection))

dimlon <- ncdim_def(name = "lon", units = "degrees_east", vals = as.vector(lon), longname = "longitude" )
dimlat <- ncdim_def(name = "lat", units = "degrees_north", vals = as.vector(lat), longname = "latitude")
dimseason <- ncdim_def(name = "season", units = "season", vals = 1 : 4, longname = "season of the year: DJF, MAM, JJA, SON")
defdata <- ncvar_def(name = "VulnerabilityIndex", units = "number_of_days", dim = list(season = dimseason, lat = dimlat, lon = dimlon), longname = "Number of days exceeding in 5 degrees the Diurnal Temeprature Range for the reference period")

file <- nc_create(paste0(plot_dir, "/Seasonal_DTRindicator_", model_names, "_", start_projection, "_", end_projection,"_",
                        start_historical, "_", end_historical, ".nc"), list(defdata))
ncvar_put(file, defdata, dtr_rcp)
#ncatt_put(file, 0, "Conventions", "CF-1.5")
nc_close(file)



