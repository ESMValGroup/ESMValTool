## Regimes namelist
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
model_names <- lapply(input_files_per_var, function(x) x$dataset)
model_names <- unlist(unname(model_names))


## Do not print warnings
#options(warn=-1)


#Var considered
var0 <- lapply(input_files_per_var, function(x) x$short_name)


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
ncenters <- params$ncenters
cluster_method <- params$cluster_method
EOFS <- params$EOFS


library(s2dverification)
library(ggplot2)
library(multiApply)
library(devtools)
library(startR)
library(magic.bsc, lib.loc = '/home/Earth/nperez/git/magic.bsc.Rcheck/')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Regimes/R/WeatherRegime.R')
#source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
#source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-debug-plot-ts/R/PlotTimeSeries.R')


fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]
data <- Start(model = fullpath_filenames,
              var = var0,
              var_var = 'var_names',
              #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
             # time = values(list(as.POSIXct(start_historical),
              #                   as.POSIXct(end_historical))),
    time ='all',
              time_tolerance = as.difftime(15, units = 'days'),
              lat = values(list(lat.min, lat.max)),
              lon = values(list(lon.min, lon.max)),
    #lat = 'all', lon = 'all',
              lon_var = 'lon',
              lon_reorder = CircularSort(0, 360),
              return_vars = list(time = 'model', lon = 'model', lat = 'model'),
              retrieve = TRUE)
      # ------------------------------
# Provisional solution to error in dimension order:
 lon <- attr(data, "Variables")$dat1$lon
 lat <- attr(data, "Variables")$dat1$lat
 time <- attr(data, "Variables")$dat1$time
    data <- as.vector(data)
    dim(data) <- c(model = 1, var = 1, lon = length(lon), lat = length(lat), time = length(time))
    data <- aperm(data, c(1,2,5,4,3))
     attr(data, "Variables")$dat1$time <- time
    print(dim(data))
# ------------------------------



time_dim <- which(names(dim(data)) == "time")
dims <- dim(data)
dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
dims <- dims[-time_dim]
dim(data) <- dims
names(dim(data))[c(time_dim, time_dim + 1)] <- c("ftime", "sdate")
Loess<-function(clim,loess_span){
  data<-data.frame(ensmean=clim,day=1:length(clim))
  loess_filt<-loess(ensmean~day,data,span=loess_span)
  output<-predict(loess_filt)
  return(output)
}

clim <- Clim(var_exp = data,var_obs = data,memb=T)

anom_obs <- Ano(data, clim$clim_obs)

#anom_exp<-Ano(data$mod,clim_smoothed_exp)

#WR_obs <- WeatherRegime(data = anom_obs, EOFS = FALSE, lat_weights = TRUE, lat = lat, lon = lon,
#                        ncenters = ncenters, method = cluster_method)


WR_obs <- WeatherRegime(data = anom_obs, EOFS = FALSE, lat = lat, lon = lon,# lat_weights = FALSE,
                        ncenters = ncenters, method = cluster_method)

#
titles<-paste0('freq = ', round(WR_obs$frequency, 1), '%')

# -----------------------------
# WeatherRegime.R returns an altered order of dimensions:
names(dim(WR_obs$composite)) <- c("lat", "lon", "Cluster", "Mod", "exp")
# -----------------------------



clim_frequencies<-paste0('freq = ',round(Mean1Dim(WR_obs$frequency,1),1),'%')


cosa <- aperm(drop(WR_obs$composite), c(3,1,2))

PlotLayout(PlotEquiMap, c(2, 3), lon = lon, lat = lat, var = cosa/100,
           titles = paste0(paste0('Cluster ', 1:4),' (',clim_frequencies,' )'), filled.continents = FALSE,
           axelab = FALSE, draw_separators = TRUE, subsampleg = 1, brks = seq(-16, 16, by = 2),
           bar_extra_labels = c(2, 0, 0, 0), fileout = paste0(plot_dir, '/observed_regimes.png'))



#reference<-drop(WR_obs$composite)
#names(dim(reference))<-c('lat','lon','nclust')
#WR_exp<-RegimesAssign(var_ano=anom_exp,ref_maps=reference,lats=lats)

