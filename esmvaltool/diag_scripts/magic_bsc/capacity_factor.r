### To do: extrapolate surface wind to 100m wind (different for land and sea)

####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4

#R package dependencies installation script
#install.packages('yaml')
#install.packages('devtools')
#library(devtools)
#Sys.setenv(TAR = '/bin/tar')
#install_git('https://earth.bsc.es/gitlab/es/startR', branch = 'develop-hotfixes-0.0.2')
#install_git('https://earth.bsc.es/gitlab/es/easyNCDF', branch = 'master')


Sys.setenv(TAR = '/bin/tar')

library(startR)
library(multiApply)
library(ggplot2)
library(yaml)
library(s2dverification)
#library(lubridate)
#source("~/PycharmProjects/ESMValTool/esmvaltool/diag_scripts/magic_bsc/PC.r")
#source("/home/Earth/nperez/git/ESMValTool/esmvaltool/diag_scripts/magic_bsc/PC.r")
library(climdex.pcic)


#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
rdir  <- sub("\\/recipes/.*", "", params$work_dir)
source(paste0(rdir,"/diag_scripts/magic_bsc/PC.r"))
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
start_year <- lapply(input_files_per_var, function(x) x$start_year)
start_year <- c(unlist(unname(start_year)))[1]
end_year <- lapply(input_files_per_var, function(x) x$end_year)
end_year <- c(unlist(unname(end_year)))[1]
seasons <- params$seasons
power_curves <- params$power_curves

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
no_of_years <- length(start_year : end_year)

time_dim <- which(names(dim(data)) == "time")

## TODO extrapolate from surface wind to 100m height
#---------------------------
# We assume power law and s sheaering exponents:
# land: 0.143
# sea: 0.11
# spd100 = spd10*(100/10)^0.11 = spd10*1.29
# spd100 = spd10*(100/10)^0.143 = spd10*1.39

# ratio <- ifelse(landmask > 50,1.39,1.29)

days <- as.Date(attr(data, "Variables")$dat1$time)
print(dim(data))
print(length(days))
print(no_of_years)
dims <- dim(data)
dims <- append(dims[-time_dim], c(no_of_years, dims[time_dim] / no_of_years), after = 1)
#print(data[1,1,1:5, 1,1])
print("CC")
print(dims)
dims <- dims[-c(1, 4)]

dim(data) <- dims
data <- aperm(data, c(2,1,3,4))
names(dim(data)) <- c("year", "day", "lat", "lon")
#print(data[1,1:5,1,1])

#####################################
# Cross with PC
#####################################

#---------------------------
# Load PC to use and compute CF for 6h values
#---------------------------
seas_data <- Mean1Dim(data,2)

#for (power_curve in power_curves){
#    pc = read_xml_pc(power_curve)
#    data_cf <- wind2CF(data,pc1)
#    dim(data_cf) <- dim(data)
#    seas_data_cf <- Mean1Dim(data_cf,2)

#}

pc1 <- read_xml_pc("/home/Earth/llledo/Documents/Power_Curves/Windographer_library/Enercon_E70_2.3MW.wtp")
pc2 <- read_xml_pc("/home/Earth/llledo/Documents/Power_Curves/Windographer_library/Gamesa_G80_2.0MW.wtp")
pc3 <- read_xml_pc("/home/Earth/llledo/Documents/Power_Curves/Windographer_library/Gamesa_G87_2.0MW.wtp")
pc4 <- read_xml_pc("/home/Earth/llledo/Documents/Power_Curves/Windographer_library/Vestas_V100_2.0MW.wtp")
pc5 <- read_xml_pc("/home/Earth/llledo/Documents/Power_Curves/Windographer_library/Vestas_V110_2.0MW.wtp")
#pc4 <- read_pc("/home/Earth/llledo/Repos/WindPower/archive/V100-2MW.txt")

data_cf1 <- wind2CF(data,pc1)
dim(data_cf1) <- dim(data)
data_cf2 <- wind2CF(data,pc2)
dim(data_cf2) <- dim(data)
data_cf3 <- wind2CF(data,pc3)
dim(data_cf3) <- dim(data)
data_cf4 <- wind2CF(data,pc4)
dim(data_cf4) <- dim(data)
data_cf5 <- wind2CF(data,pc5)
dim(data_cf5) <- dim(data)

#---------------------------
# Aggregate daily data to seasonal means
#---------------------------

seas_data_cf1 <- Mean1Dim(data_cf1,2)
seas_data_cf2 <- Mean1Dim(data_cf2,2)
seas_data_cf3 <- Mean1Dim(data_cf3,2)
seas_data_cf4 <- Mean1Dim(data_cf4,2)
seas_data_cf5 <- Mean1Dim(data_cf5,2)

#save(seas_erai_gwa,seas_erai_cf1,seas_erai_cf2,seas_erai_cf3,seas_erai_cf4,seas_erai_cf5,lats,lons,variable,seasons,first_year,last_year,bbox,nsdates,nleadtime,nlat,nlon,ratio,file="/esnas/scratch/llledo/PC_sensitivity/PC_sensitivity.Rdata")

##############################
# Make some plots
##############################
#---------------------------
# Prepare data, labels and colorscales
#---------------------------
library(RColorBrewer)
library(abind)
p <- colorRampPalette(brewer.pal(9,"YlOrRd"))
q <- colorRampPalette(rev(brewer.pal(11,"RdBu")))
years <- seq(start_year,end_year)
turb_types <- c("IEC I","IEC I/II","IEC II","IEC II/III","IEC III")

seas_data_cf_all <- abind(seas_data_cf1,seas_data_cf2,seas_data_cf3,seas_data_cf4,seas_data_cf5,along=0)
mean_data_cf_all <- Mean1Dim(seas_data_cf_all,2)
anom_data_cf_all <- seas_data_cf_all-InsertDim(Mean1Dim(seas_data_cf_all,2),2,dim(data)[1])
pct_anom_data_cf_all <- (seas_data_cf_all/InsertDim(Mean1Dim(seas_data_cf_all,2),2,dim(data)[1]))-1

#---------------------------
# Plot seasonal CF maps
#---------------------------

PlotLayout(PlotEquiMap,c(3,2),Mean1Dim(seas_data_cf_all, 2), lon, lat, filled.continents=F,toptitle=paste0(seasons, " CF from ", model_names, " (", start_year, "-", end_year, ")"),
           fileout = paste0(plot_dir, "/", "capacity_factor_", model_names,  "_", start_year, "-", end_year, ".png"))

#---------------------------
# Plot seasonal CF anomalies maps
#---------------------------

PlotLayout(PlotEquiMap,c(3,2),Mean1Dim(anom_data_cf_all, 2),lon, lat, filled.continents=F,toptitle=paste0(seasons, " CF Anomaly from ", model_names, " (", start_year, "-", end_year, ")")
           ,col_titles=turb_types,color_fun=q,brks=seq(-0.25,0.25,0.05),bar_scale=0.5,title_scale=0.7,axelab = F,
           fileout = paste0(plot_dir, "/", "capacity_factor_anomaly_", model_names,  "_", start_year, "-", end_year, ".png"))


#---------------------------
# Scatterplot wind vs CF
#---------------------------
#---------------------------
# Correlation maps
#---------------------------
cor12 <- apply(seas_data_cf_all,c(3,4),function(x) {cor(x[1,],x[2,])})
cor13 <- apply(seas_data_cf_all,c(3,4),function(x) {cor(x[1,],x[3,])})
cor14 <- apply(seas_data_cf_all,c(3,4),function(x) {cor(x[1,],x[4,])})
cor15 <- apply(seas_data_cf_all,c(3,4),function(x) {cor(x[1,],x[5,])})
cor24 <- apply(seas_data_cf_all,c(3,4),function(x) {cor(x[2,],x[4,])})
cor35 <- apply(seas_data_cf_all,c(3,4),function(x) {cor(x[3,],x[5,])})


PlotLayout(PlotEquiMap,c(1,2),list(cor13^2,cor35^2,cor24^2,cor15^2),lon,lat,nrow=2,ncol=2,filled.continents=F,toptitle="Seasonal CF determination coef.",
           titles=c("between cf1 and cf3","between cf3 and cf5","between cf2 and cf4","between cf1 and cf5"),brks=c(0.,0.3,0.5,0.6,0.7,0.8,0.9,0.93,0.96,0.98,0.99,1),bar_scale=0.5,
           title_scale=0.7,axelab=F,color_fun=p, fileout = paste0(plot_dir, "/", "capacity_factor_correlation_maps_", model_names,  "_", start_year, "-", end_year, ".png"))


#---------------------------
# RMSE maps
#---------------------------
rmse <- function(x,y) { sqrt(mean((x-y)^2,na.rm=T)) }
rmse12 <- apply(anom_data_cf_all,c(3,4),function(x) {rmse(x[1,],x[2,])})
rmse13 <- apply(anom_data_cf_all,c(3,4),function(x) {rmse(x[1,],x[3,])})
rmse35 <- apply(anom_data_cf_all,c(3,4),function(x) {rmse(x[3,],x[5,])})
rmse15 <- apply(anom_data_cf_all,c(3,4),function(x) {rmse(x[1,],x[5,])})
rmse24 <- apply(anom_data_cf_all,c(3,4),function(x) {rmse(x[2,],x[4,])})

PlotLayout(PlotEquiMap,c(1,2),list(rmse13,rmse35,rmse24,rmse15),lon,lat,nrow=2,ncol=2,filled.continents=F,toptitle="Seasonal CF RMSE",
           titles=c("between cf1 and cf3","between cf3 and cf5","between cf2 and cf4","between cf1 and cf5"),brks=seq(0,0.08,0.01),bar_scale=0.5,
           title_scale=0.7,axelab=F,color_fun=p, fileout = paste0(plot_dir, "/", "capacity_factor_rmse_maps_", model_names,  "_", start_year, "-", end_year, ".png"))

