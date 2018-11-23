### To do: extrapolate surface wind to 100m wind (different for land and sea)

# nolint start
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

# nolint end

Sys.setenv(TAR = "/bin/tar") # nolint

library(multiApply) # nolint
library(ggplot2)
library(yaml)
library(s2dverification)
library(climdex.pcic)
library(ncdf4)

#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
initial.options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(
    file_arg_name, "", initial.options[grep(file_arg_name, initial.options)]
)
script_dirname <- dirname(script_name)
source(file.path(script_dirname, "PC.r"))

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
start_year <- lapply(input_files_per_var, function(x) x$start_year)
start_year <- c(unlist(unname(start_year)))[1]
end_year <- lapply(input_files_per_var, function(x) x$end_year)
end_year <- c(unlist(unname(end_year)))[1]
seasons <- params$seasons
power_curves <- params$power_curves
power_curves_folder <- params$power_curves_folder

no_of_years <- length(start_year : end_year)
var0 <- unlist(var0)
data_nc <- nc_open(fullpath_filenames)
data <- ncvar_get(data_nc, var0)

names(dim(data)) <- c("lon", "lat", "time")
lat <- ncvar_get(data_nc,"lat")
lon <- ncvar_get(data_nc,"lon")
units <- ncatt_get(data_nc, var0, "units")$value
calendar <- ncatt_get(data_nc, "time", "calendar")$value
long_names <-  ncatt_get(data_nc,var0,"long_name")$value
time <-  ncvar_get(data_nc,"time")
start_date <- as.POSIXct(substr(ncatt_get(data_nc, "time", "units")$value,11, 29 ))
nc_close(data_nc)
time <- as.Date(time, origin = start_date, calendar = calendar)



time_dim <- which(names(dim(data)) == "time")

# nolint start
## TODO extrapolate from surface wind to 100m height
#---------------------------
# We assume power law and s sheaering exponents:
# land: 0.143
# sea: 0.11
# spd100 = spd10*(100/10)^0.11 = spd10*1.29
# spd100 = spd10*(100/10)^0.143 = spd10*1.39

# ratio <- ifelse(landmask > 50,1.39,1.29)
# nolint end

days <- time
print(dim(data))
print(length(days))
print(no_of_years)
dims <- dim(data)
dims <- append(
    dims[-time_dim], c(no_of_years, dims[time_dim] / no_of_years), after = 1
)
print("CC")
print(dims)
#dims <- dims[-c(1, 4)]

dim(data) <- dims
data <- aperm(data, c(2, 1, 3, 4))
names(dim(data)) <- c("year", "day", "lat", "lon")

#####################################
# Cross with PC
#####################################

#---------------------------
# Load PC to use and compute CF for 6h values
#---------------------------
seas_data <- Mean1Dim(data, 2)

pc1 <- read_xml_pc(file.path(power_curves_folder, "Enercon_E70_2.3MW.wtp"))
pc2 <- read_xml_pc(file.path(power_curves_folder, "Gamesa_G80_2.0MW.wtp"))
pc3 <- read_xml_pc(file.path(power_curves_folder, "Gamesa_G87_2.0MW.wtp"))
pc4 <- read_xml_pc(file.path(power_curves_folder, "Vestas_V100_2.0MW.wtp"))
pc5 <- read_xml_pc(file.path(power_curves_folder, "Vestas_V110_2.0MW.wtp"))

data_cf1 <- wind2CF(data, pc1)
dim(data_cf1) <- dim(data)
data_cf2 <- wind2CF(data, pc2)
dim(data_cf2) <- dim(data)
data_cf3 <- wind2CF(data, pc3)
dim(data_cf3) <- dim(data)
data_cf4 <- wind2CF(data, pc4)
dim(data_cf4) <- dim(data)
data_cf5 <- wind2CF(data, pc5)
dim(data_cf5) <- dim(data)

#---------------------------
# Aggregate daily data to seasonal means
#---------------------------

seas_data_cf1 <- Mean1Dim(data_cf1, 2)
seas_data_cf2 <- Mean1Dim(data_cf2, 2)
seas_data_cf3 <- Mean1Dim(data_cf3, 2)
seas_data_cf4 <- Mean1Dim(data_cf4, 2)
seas_data_cf5 <- Mean1Dim(data_cf5, 2)

# nolint start
# save(
#   seas_erai_gwa,
#   seas_erai_cf1, seas_erai_cf2, seas_erai_cf3,
#   seas_erai_cf4, seas_erai_cf5, lats, lons, variable, seasons, first_year,
#   last_year, bbox, nsdates, nleadtime, nlat, nlon, ratio,
#   file="/esnas/scratch/llledo/PC_sensitivity/PC_sensitivity.Rdata"
# )
# nolint end

##############################
# Make some plots
##############################
#---------------------------
# Prepare data, labels and colorscales
#---------------------------
library(RColorBrewer) # nolint
library(abind)
p <- colorRampPalette(brewer.pal(9, "YlOrRd"))
q <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
years <- seq(start_year, end_year)
turb_types <- c("IEC I", "IEC I/II", "IEC II", "IEC II/III", "IEC III")

seas_data_cf_all <- abind(
    seas_data_cf1, seas_data_cf2, seas_data_cf3, seas_data_cf4, seas_data_cf5,
    along = 0
)
mean_data_cf_all <- Mean1Dim(seas_data_cf_all, 2)
anom_data_cf_all <- seas_data_cf_all - InsertDim( # nolint
    Mean1Dim(seas_data_cf_all, 2), 2, dim(data)[1] # nolint
)
pct_anom_data_cf_all <- (seas_data_cf_all / InsertDim( # nolint
    Mean1Dim(seas_data_cf_all, 2), 2, dim(data)[1] # nolint
)) - 1

#---------------------------
# Plot seasonal CF maps
#---------------------------

PlotLayout( # nolint
    PlotEquiMap, # nolint
    c(3, 2),
    Mean1Dim(seas_data_cf_all, 2),
    lon,
    lat,
    filled.continents = F,
    toptitle = paste0(
        seasons, " CF from ",
        model_names, " (", start_year, "-", end_year, ")"
    ),
    fileout = paste0(
        plot_dir, "/", "capacity_factor_",
        model_names,  "_", start_year, "-", end_year, ".png"
    )
)

#---------------------------
# Plot seasonal CF anomalies maps
#---------------------------

PlotLayout( #nolint
    PlotEquiMap, # nolint
    c(3, 2),
    Mean1Dim(anom_data_cf_all, 2),
    lon,
    lat,
    filled.continents = F,
    toptitle = paste0(
        seasons, " CF Anomaly from ", model_names,
        " (", start_year, "-", end_year, ")"
    ),
    col_titles = turb_types,
    color_fun = q,
    brks = seq(-0.25, 0.25, 0.05),
    bar_scale = 0.5,
    title_scale = 0.7,
    axelab = F,
    fileout = paste0(
        plot_dir, "/", "capacity_factor_anomaly_", model_names,
        "_", start_year, "-", end_year, ".png"
    )
)


#---------------------------
# Scatterplot wind vs CF
#---------------------------
#---------------------------
# Correlation maps
#---------------------------
corr <- function(data, i, j){
  apply(
    data,
    c(3, 4),
    function(x) {
      cor(x[i, ], x[j, ])
    }
  )
}
cor12 <- corr(data, 1, 2)
cor12 <- corr(data, 1, 3)
cor14 <- corr(data, 1, 4)
cor15 <- corr(data, 1, 5)
cor24 <- corr(data, 2, 4)
cor35 <- corr(data, 3, 5)

PlotLayout( # nolint
    PlotEquiMap, # nolint
    c(1, 2),
    list(cor13 ^ 2, cor35 ^ 2, cor24 ^ 2, cor15 ^ 2),
    lon, lat, nrow = 2, ncol = 2,
    filled.continents = F,
    toptitle = "Seasonal CF determination coef.",
    titles = c(
        "between cf1 and cf3",
        "between cf3 and cf5",
        "between cf2 and cf4",
        "between cf1 and cf5"
    ),
    brks = c(0., 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 0.93, 0.96, 0.98, 0.99, 1),
    bar_scale = 0.5,
    title_scale = 0.7,
    axelab = F,
    color_fun = p,
    fileout = paste0(
        plot_dir, "/", "capacity_factor_correlation_maps_", model_names,
        "_", start_year, "-", end_year, ".png"
    )
)

#---------------------------
# RMSE maps
#---------------------------
rmse <- function(data, i, j){
  apply(
    data,
    c(3, 4),
    function(x, y){
      sqrt(mean( (x[i, ] - x[j, ]) ^ 2, na.rm = T))
    }
  )
}
rmse12 <- rmse(anom_data_cf_all, 1, 2)
rmse13 <- rmse(anom_data_cf_all, 1, 3)
rmse35 <- rmse(anom_data_cf_all, 3, 5)
rmse15 <- rmse(anom_data_cf_all, 1, 5)
rmse24 <- rmse(anom_data_cf_all, 2, 4)

PlotLayout( # nolint
    PlotEquiMap, # nolint
    c(1, 2),
    list(rmse13, rmse35, rmse24, rmse15),
    lon,
    lat,
    nrow = 2,
    ncol = 2,
    filled.continents = F,
    toptitle = "Seasonal CF RMSE",
    titles = c(
        "between cf1 and cf3",
        "between cf3 and cf5",
        "between cf2 and cf4",
        "between cf1 and cf5"
    ),
    brks = seq(0, 0.08, 0.01),
    bar_scale = 0.5,
    title_scale = 0.7,
    axelab = F,
    color_fun = p,
    fileout = paste0(
        plot_dir, "/", "capacity_factor_rmse_maps_", model_names,
        "_", start_year, "-", end_year, ".png"
    )
)
