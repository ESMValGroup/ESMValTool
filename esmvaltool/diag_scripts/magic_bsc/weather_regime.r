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
ncenters <- params$ncenters
cluster_method <- params$cluster_method
EOFS <- params$EOFS


library(s2dverification)
library(startR)
library(ggplot2)
library(multiApply)
library(devtools)
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Regimes/R/WeatherRegime.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
#source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-debug-plot-ts/R/PlotTimeSeries.R')

fullpath_filenames <- input_files[[var0]]
data <- Start(model = fullpath_filenames,
              var = var0,
              var_var = 'var_names',
              #  sdate = paste0(seq(1970, 2000, by = 1), "0101"),
              time = values(list(as.POSIXct(start_historical), 
                                 as.POSIXct(end_historical))),
              time_tolerance = as.difftime(15, units = 'days'), 
              lat = values(list(lat.min, lat.max)),
              lon = values(list(lon.min, lon.max)),
              lon_var = 'lon',
              lon_reorder = CircularSort(0, 360),
              return_vars = list(time = 'model', lon = 'model', lat = 'model'),
              retrieve = TRUE)

print(dim(data))

lat <- attr(data, "Variables")$dat1$lat
lon <- attr(data, "Variables")$dat1$lon

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

WR_obs <- WeatherRegime(data = anom_obs, EOFS = FALSE, lat_weights = TRUE, lat = lat, lon = lon, 
                        ncenters = ncenters, method = cluster_method)

titles<-paste0('freq = ', round(WR_obs$frequency, 1), '%')
PlotLayout(PlotEquiMap,c(1,2),lon=lon,lat=lat,var=WR_obs$composite/100,
           titles=paste0(paste0('Cluster ',1:4),' (',paste0('freq = ',round(WR_obs$frequency,1),'%'),' )'),filled.continents=F,
           axelab=F,draw_separators = T,subsampleg = 1,brks=seq(-16,16,by=2),
           bar_extra_labels = c(2,0,0,0),fileout= paste0(plot_dir, '/observed_regimes.png'))


#reference<-drop(WR_obs$composite)
#names(dim(reference))<-c('lat','lon','nclust')
#WR_exp<-RegimesAssign(var_ano=anom_exp,ref_maps=reference,lats=lats)

