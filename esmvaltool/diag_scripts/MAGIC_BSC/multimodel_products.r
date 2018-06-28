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
library(s2dverification)
library(startR, lib.loc='/home/Earth/ahunter/R/x86_64-unknown-linux-gnu-library/3.2/')
library(multiApply)
library(ggplot2)
library(yaml)

##Until integrated into current version of s2dverification
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-MagicWP5/R/AnoAgree.R')
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
##Until integrated into current version of s2dverification
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')
source("https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/SelBox.R")



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

climatology_class <- params$climatology_class
anomaly_class <- params$anomaly_class
climatology_files <- which(unname(experiment) == as.character(climatology_class))
anomaly_files <- which(unname(experiment) == as.character(anomaly_class))

model_names <-  lapply(input_files_per_var, function(x) x$model)
model_names <- unlist(unname(model_names))[anomaly_files]

start_climatology <- lapply(input_files_per_var, function(x) x$start_year)
start_climatology <- c(unlist(unname(start_climatology))[climatology_files])[1]
end_climatology <- lapply(input_files_per_var, function(x) x$end_year)
end_climatology <- c(unlist(unname(end_climatology))[climatology_files])[1]

start_anomaly <- lapply(input_files_per_var, function(x) x$start_year)
start_anomaly <- c(unlist(unname(start_anomaly))[anomaly_files])[1]
end_anomaly <- lapply(input_files_per_var, function(x) x$end_year)
end_anomaly <- c(unlist(unname(end_anomaly))[anomaly_files])[1]

agreement_threshold <- params$agreement_threshold
print(start_climatology)
print(end_anomaly)

font_size <- 12

## Do not print warnings
#options(warn=-1)

#Parameters for Season() function
monini <- params$monini
moninf <- params$moninf
monsup <- params$monsup
time_series_plot <- params$time_series_plot

### Load data and compute climatologies and anomalies

climatology_filenames <- fullpath_filenames[climatology_files]
reference_data <- Start(model = climatology_filenames,
                        var = var0,
                        var_var = 'var_names',
                        time = 'all',
                        lat = 'all',
                        lon = 'all',
                        lon_var = 'lon',
                        return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                        retrieve = TRUE)

lat <- attr(reference_data, "Variables")$dat1$lat
lon <- attr(reference_data, "Variables")$dat1$lon

if (!is.null(moninf)) {
  months <- paste0(month.abb[moninf],"-", month.abb[monsup])
} else {
  months <- ""
}


attributes(lon) <- NULL
attributes(lat) <- NULL
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))

time_dim <- which(names(dim(reference_data)) == "time")
dims <- dim(reference_data)

dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
dims <- dims[-time_dim]
dim(reference_data) <- dims
names(dim(reference_data))[c(time_dim, time_dim + 1)] <- c("month", "year")

## If moninf is not specified, compute the climatolgy and mean anomaly for each month
if (!is.null(moninf)){
reference_seasonal_mean <- Season(reference_data, posdim = time_dim, monini = monini, moninf = moninf,
                                   monsup = monsup)
} else {
 reference_seasonal_mean <- reference_data
}

margins <- list(c(1 : length(dim(reference_seasonal_mean)))[-c(time_dim + 1)])
years_dim <- which(names(dim(reference_seasonal_mean)) == "year")
climatology <- Mean1Dim(reference_seasonal_mean, years_dim)

anomaly_filenames <- fullpath_filenames[anomaly_files]
rcp_data <- Start(model = anomaly_filenames,
                  var = var0,
                  var_var = 'var_names',
                  time = 'all',
                  lat = 'all',
                  lon = 'all',
                  lon_var = 'lon',
                  lon_reorder = CircularSort(0, 360),
                  return_vars = list(time = 'model', lon = 'model',
                                     lat = 'model'),
                  retrieve = TRUE)



time_dim <- which(names(dim(rcp_data)) == "time")
dims <- dim(rcp_data)
dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
dims <- dims[-time_dim]
dim(rcp_data) <- dims
names(dim(rcp_data))[c(time_dim, time_dim + 1)] <- c("month", "year")
#
attr <- ((attr(rcp_data,"Variables")$common))[2]
units <- attr[[1]]$units

if (!is.null(moninf)) {
  proj_seasonal_mean <- Season(rcp_data, posdim = time_dim, monini = monini, moninf = moninf,
                             monsup = monsup)
   years_dim <- which(names(dim(proj_seasonal_mean)) == "year")
  multi_year_anomaly <- Mean1Dim(proj_seasonal_mean, years_dim) - climatology
  climatology <- InsertDim(climatology, years_dim, lendim = dim(proj_seasonal_mean)[years_dim])
} else {
  proj_seasonal_mean <- rcp_data
}
anomaly <- proj_seasonal_mean - climatology
#


if (!is.null(moninf)){
time <- seq(start_anomaly, end_anomaly, by = 1)
month <- moninf
  if (month <= 9) {
    month <- paste0(as.character(0), as.character(month))
  }
  month <- paste0("-", month, "-")
  day <- "01"
  time <- as.POSIXct(paste0(time, month, day), tz = "CET")
  time <- julian(time, origin = as.POSIXct("1970-01-01"))
} else {
 day <- "01"
 time <- sapply(start_anomaly : end_anomaly,rep, 12)
 month <- c("-01-", "-02-", "-03-", "-04-", "-05-", "-06-", "-07-", "-08-", "-09-", "-10-", "-11-", "-12-")
 time <- as.POSIXct(paste0(time, month, day), tz = "CET")
}



attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
attr(time, "variables") <- metadata
#Save the single model anomalies
for (mod in 1 : length(model_names)) {
  if (!is.null(moninf)) {
     data <- anomaly[mod,1,1, , ,]
  } else {
    data <- anomaly[mod,1,,1 , ,]
  }
  data <- anomaly[mod,1,1, , ,]
  data <- aperm(data, c(2,3,1))
  names(dim(data)) <- c("lat", "lon", "time")
  metadata <- list(variable = list(dim = list(list(name='time', unlim = FALSE)), units = units ))
  names(metadata)[1] <- var0
  attr(data, 'variables') <- metadata
  variable_list <- list(variable = data, lat = lat, lon = lon, time = time)
  names(variable_list)[1] <- var0
  ArrayToNetCDF(variable_list,
                paste0(plot_dir,  "/", var0, "_", months, "_anomaly_",model_names[mod],"_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".nc"))
}


#Compute and save the multi-model anomalies
#data <- Mean1Dim(anomaly, 1)
#if (!is.null(moninf)) {
#     data <- anomaly[mod,1,1, , ,]
#} else {
#  data <- anomaly[mod,1,,1 , ,]
#  months <- "annual"
#}
#data <- aperm(data, c(2,3,1))


### Plot timeseries

model_anomalies <- WeightedMean(anomaly, lon = lon, lat = lat, mask = NULL)

 if (!is.null(params$running_mean)) {
  model_anomalies <- Smoothing(model_anomalies, runmeanlen = params$running_mean, numdimt = 4)
}

data_frame <- as.data.frame.table(t(model_anomalies[,1,1,]))
years <- rep(start_anomaly : end_anomaly, dim(model_anomalies)[1])
data_frame$Year <- c(years)
names(data_frame)[2] <- "Model"
for (i in 1 : length(levels(data_frame$Model))) {
levels(data_frame$Model)[i] <- model_names[i]
}


if (!is.null(time_series_plot)) {
   g <- ggplot(data_frame, aes(x = Year, y = Freq, color = Model)) + theme_bw() +
        geom_line() + ylab(paste0("Anomaly (", units, ")")) +  xlab("Year") + theme(text=element_text(size = font_size),legend.text=element_text(size = font_size),
                                                                                    axis.title=element_text(size = font_size)) +
         stat_summary(data =  data_frame, fun.y= "mean", mapping = aes(x = data_frame$Year, y = data_frame$Freq, group = interaction(data_frame[2,3]),
                                                                       color = data_frame$Model), geom = "line", size = 1) +
         ggtitle(paste0(months, " ", var0, " anomaly (", start_anomaly, "-", end_anomaly,") - ", "(", start_climatology, "-", end_climatology,")"))
} else {
  g <- ggplot(data_frame, aes(x = Year, y = Freq)) + theme_bw() +
        ylab(paste0("Anomaly (", units, ")")) +  xlab("Year") + theme(text=element_text(size = font_size),legend.text=element_text(size = font_size),
                                                                      axis.title=element_text(size = font_size)) +
        stat_summary(data =  data_frame, fun.y= "mean", mapping = aes(x = data_frame$Year, y = data_frame$Freq, group = interaction(data_frame[2,3]),
                                                                       color = data_frame$Model), geom = "line", size = 0.8) +
        stat_summary(data =  data_frame, geom = "ribbon", fun.ymin = "min", fun.ymax = "max", mapping = aes(x = data_frame$Year, y = data_frame$Freq, group = interaction(data_frame[2,3])), alpha = 0.3, color = "red", fill = "red") +
         ggtitle(paste0(months, " ", var0, " anomaly (", start_anomaly, "-", end_anomaly,") - ", "(", start_climatology, "-", end_climatology,")"))
}


     ggsave(filename = paste0(plot_dir, "/", "Area-averaged ", var0, "_",months, "_multimodel-anomaly_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".png"), g, device = NULL)


##Plot maps

if (!is.null(agreement_threshold)) {

  model_dim <- which(names(dim(multi_year_anomaly)) == "model")
  agreement <- AnoAgree(multi_year_anomaly + rnorm(3*17*19), members_dim = model_dim)
} else {
  agreement_threshold <- 1000
  agreement <- NULL
}

colorbar_lim <- ceiling(max(abs(max(multi_year_anomaly)),abs(min(data))))
brks <- seq(-colorbar_lim, colorbar_lim, length.out = 11)
title <- paste0(months, " ", var0, " anomaly (", start_anomaly, "-", end_anomaly, ") - (", start_climatology, "-", end_climatology, ")")
data <- drop(Mean1Dim(multi_year_anomaly, model_dim))
PlotEquiMap(data, lat = lat, lon = lon, brks = brks, units =units, toptitle = title, filled.continents = FALSE,
            dots = drop(agreement) >= agreement_threshold,
            fileout = paste0(plot_dir, "/", var0, "_",months, "_multimodel-anomaly_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".png"))
data <- InsertDim(data, 3, 1)
names(dim(data)) <- c("lat", "lon", "time")
metadata <- list(variable = list(dim = list(list(name='time', unlim = FALSE)), units = units))
names(metadata)[1] <- var0
attr(data, 'variables') <- metadata
model_names_filename <- paste(model_names, collapse = '_')
agreement <- adrop(agreement,1)
agreement <- aperm(agreement, c(2,3,1))
names(dim(agreement)) <- c("lat", "lon", "time")
metadata <- list(variable = list(dim = list(list(name='time', unlim = FALSE)), units = "%"))
names(metadata)[1] <- "agreement"
attr(agreement, 'variables') <- metadata
time <- time[1]
attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
attr(time, "variables") <- metadata
variable_list <- list(variable = data, agreement = agreement, lat = lat, lon = lon, time = time)
names(variable_list)[1] <- var0

ArrayToNetCDF(variable_list,  paste0(plot_dir, "/", var0, "_",months, "_multimodel-anomaly_",
              model_names_filename,"_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".nc"))


