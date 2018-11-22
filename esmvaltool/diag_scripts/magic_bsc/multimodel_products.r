# nolint start
####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4

#R package dependencies installation script
#install.packages("yaml")
#install.packages("devtools")
#library(devtools)
#Sys.setenv(TAR = "/bin/tar")
#install_git("https://earth.bsc.es/gitlab/es/startR", branch = "develop-hotfixes-0.0.2")
#install_git("https://earth.bsc.es/gitlab/es/easyNCDF", branch = "master")
# nolint end


Sys.setenv(TAR = "/bin/tar") # nolint
library(s2dverification)
library(startR) # nolint
library(ClimProjDiags) # nolint
library(abind)
library(ggplot2)
library(yaml)
library(ncdf4)
library(multiApply)

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
model_names <- unique(unlist(unname(model_names)))

var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]
experiment <- lapply(input_files_per_var, function(x) x$exp)
experiment <- unlist(unname(experiment))

climatology_class <- params$climatology_class
anomaly_class <- params$anomaly_class
climatology_files <- which(
  unname(experiment) == as.character(climatology_class)
)
anomaly_files <- which(unname(experiment) == as.character(anomaly_class))

agreement_threshold <- params$agreement_threshold
start_climatology <- params$climatology_start_year
end_climatology <- params$climatology_end_year
start_anomaly <- params$anomaly_start_year
end_anomaly <- params$anomaly_end_year
font_size <- 12

#Parameters for Season() function
monini <- 1
moninf <- params$moninf
monsup <- params$monsup
if (is.null(moninf) & !is.null(monsup)){
    moninf <- monsup
} else if (!is.null(moninf) & is.null(monsup)) {
    monsup <- moninf
}
month_names <- c(
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AGO", "SEP",
  "OCT", "NOV", "DEC"
)
if (moninf == monsup) {
  months <- month_names[moninf]
} else {
  months <- paste0(month_names[moninf], month_names[monsup], sep = "-")
}

time_series_plot <- params$time_series_plot
### Load data and compute climatologies and anomalies

climatology_filenames <- fullpath_filenames[climatology_files]
reference_data <- Start(
  model = climatology_filenames,
  var = var0,
  var_var = "var_names",
  time = "all",
  lat = "all",
  lon = "all",
  lon_var = "lon",
  return_vars = list(time = "model", lon = "model", lat = "model"),
  retrieve = TRUE
)

lat <- attr(reference_data, "Variables")$dat1$lat
lon <- attr(reference_data, "Variables")$dat1$lon
long_names <- attr(reference_data, "Variables")$common$tas$long_name
projection <- attr(reference_data, "Variables")$common$tas$coordinates

# nolint start
#jpeg(paste0(plot_dir, "/plot.jpg"))
#PlotEquiMap(reference_data[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
# ------------------------------
# Provisional solution to error in dimension order and time values:
# nolint end
time <- attr(reference_data, "Variables")$dat1$time
calendar <- attributes(time)$variables$time$calendar
if ((end_climatology-start_climatology + 1) * 12 == length(time)) {
  time <- seq(
    as.Date(
      paste(start_climatology, "01", "01", sep = "-"),
      format = "%Y-%m-%d"
    ),
    as.Date(
      paste(end_climatology, "12", "01", sep = "-"),
      format = "%Y-%m-%d"
    ),
    "month"
  )
}

num_models <- dim(reference_data)[which(names(dim(reference_data))=="model")]
reference_data <- as.vector(reference_data)
dim(reference_data) <- c(
  num_models,
  var = 1,
  lon = length(lon),
  lat = length(lat),
  time = length(time)
)
reference_data <- aperm(reference_data, c(1,2,5,4,3))
attr(reference_data, "Variables")$dat1$time <- time

# nolint start
# ------------------------------
#jpeg(paste0(plot_dir, "/plot1.jpg"))
#PlotEquiMap(reference_data[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
#---------------------------------------------
# MONTHLY - SEASONAL - ANNUAL
# MONTH:  moninf = monsup
# SEASONAL: specify the moninf and monsup;
#           if winter: moninf = 12 monsup = 2;
#           any other moninf > monsup allowed
#---------------------------------------------
# nolint end

time_dim <- which(names(dim(reference_data)) == "time")
dims <- dim(reference_data)

if (moninf <= monsup) {
  dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
  dims <- dims[-time_dim]
  dim(reference_data) <- dims
  names(dim(reference_data))[c(time_dim, time_dim + 1)] <- c("month", "year")
  reference_seasonal_mean <- Season(
    reference_data,
    posdim = time_dim,
    monini = monini,
    moninf = moninf,
    monsup = monsup
  )
  reference_seasonal_mean <- adrop(adrop(reference_seasonal_mean, 2), 2)
} else {
  if (monsup == 2 & moninf == 12) {
    reference_seasonal_mean <- SeasonSelect(
      reference_data, season = "DJF", dates = time, calendar = calendar
    )$data
            # Adding one NA december at the begining
    time_dim <- which(names(dim(reference_seasonal_mean)) == "time")
    dims <- dim(reference_seasonal_mean)
    empty_array <- rep(NA, prod(dims[-time_dim]))
    dims[time_dim] = 1
    dim(empty_array) <- dims[-time_dim]
    nom <- names(dim(reference_seasonal_mean))
    reference_seasonal_mean <- abind(
      reference_seasonal_mean, empty_array, along = time_dim
    )
    # and removing the last december
    names(dim(reference_seasonal_mean)) <- nom
    dimensiones <- 1 : length(dim(reference_seasonal_mean))
    reference_seasonal_mean <- Apply(
      reference_seasonal_mean,
      target_dims = time_dim,
      fun = function(x){x[1 : (length(x) -1)]}
    )$output1
    dims <- dim(reference_seasonal_mean)
    time_dim <- which(names(dim(reference_seasonal_mean)) == "time")
    dims <- append(dims, c(3, dims[time_dim] / 3), after = time_dim)
    dims <- dims[-time_dim]
    dim(reference_seasonal_mean) <- dims
    names(dim(reference_seasonal_mean))[c(time_dim, time_dim + 1)] <- c("season", "year")
    reference_seasonal_mean <- Mean1Dim(reference_seasonal_mean, posdim = time_dim)
  }
}

margins <- list(c(1 : length(dim(reference_seasonal_mean)))[-c(time_dim + 1)])
years_dim <- which(names(dim(reference_seasonal_mean)) == "year")
climatology <- Mean1Dim(reference_seasonal_mean, years_dim)

anomaly_filenames <- fullpath_filenames[anomaly_files]
rcp_data <- Start(
  model = anomaly_filenames,
  var = var0,
  var_var = "var_names",
  time = "all",
  lat = "all",
  lon = "all",
  lon_var = "lon",
  #lon_reorder = CircularSort(0, 360),
  return_vars = list(time = "model", lon = "model", lat = "model"),
  retrieve = TRUE
)
# nolint start
#jpeg(paste0(plot_dir, "/plot2.jpg"))
#PlotEquiMap(rcp_data[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()
# ------------------------------
# Provisional solution to error in dimension order:
# nolint end
lon <- attr(rcp_data, "Variables")$dat1$lon
lat <- attr(rcp_data, "Variables")$dat1$lat
time <- attr(rcp_data, "Variables")$dat1$time
if (attributes(time)$variables$time$calendar != calendar) {
  print("Different calendars between climatology and anomaly.")
}
if ( (end_anomaly-start_anomaly + 1) * 12 == length(time)) {
  time <- seq(
    as.Date(
      paste(start_anomaly, "01", "01", sep = "-"),
      format = "%Y-%m-%d"
    ),
    as.Date(
      paste(end_anomaly, "12", "01", sep = "-"),
      format = "%Y-%m-%d"
    ),
    "month"
  )
}
num_models <- dim(rcp_data)[which(names(dim(rcp_data))=="model")]
units <- (attr(rcp_data,"Variables")$common)[[2]]$units
rcp_data <- as.vector(rcp_data)
dim(rcp_data) <- c(
  num_models,
  var = 1,
  lon = length(lon),
  lat = length(lat),
  time = length(time)
)
rcp_data <- aperm(rcp_data, c(1,2,5,4,3))
attr(rcp_data, "Variables")$dat1$time <- time

# nolint start
# ------------------------------
#jpeg(paste0(plot_dir, "/plot3.jpg"))
#PlotEquiMap(rcp_data[1,1,1,,], lon = lon, lat = lat, filled = F)
#dev.off()


#---------------------------------------------
# MONTHLY - SEASONAL - ANNUAL
# MONTH:  moninf = monsup
# SEASONAL: specify the moninf and monsup;
#           if winter: moninf = 12 monsup = 2;
#           any other moninf > monsup allowed
#---------------------------------------------
# nolint end

time_dim <- which(names(dim(rcp_data)) == "time")
dims <- dim(rcp_data)
mes <- as.numeric(substr(time, 6, 7))

if (moninf <= monsup) {
  dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
  dims <- dims[-time_dim]
  dim(rcp_data) <- dims
  names(dim(rcp_data))[c(time_dim, time_dim + 1)] <- c("month", "year")
  rcp_seasonal_mean <- Season(
    rcp_data,
    posdim = time_dim,
    monini = monini,
    moninf = moninf,
    monsup = monsup
  )
  rcp_seasonal_mean <- adrop(adrop(rcp_seasonal_mean, 2), 2)
} else {
  if (monsup == 2 & moninf == 12) {
    rcp_seasonal_mean <- SeasonSelect(
      rcp_data, season = "DJF",
      dates = time,
      calendar = calendar
    )$data
    time_dim <- which(names(dim(rcp_seasonal_mean)) == "time")
    dims <- dim(rcp_seasonal_mean)
    empty_array <- rep(NA, prod(dims[-time_dim]))
    dims[time_dim] = 1
    dim(empty_array) <- dims[-time_dim]
    nom <- names(dim(rcp_seasonal_mean))
    rcp_seasonal_mean <- abind(
      rcp_seasonal_mean, empty_array, along = time_dim
    )
    borrar <- dim(rcp_seasonal_mean)[time_dim]
    names(dim(rcp_seasonal_mean)) <- nom
    dimensiones <- 1 : length(dim(rcp_seasonal_mean))
    rcp_seasonal_mean <- Apply(
      rcp_seasonal_mean,
      target_dims = time_dim,
      fun = function(x){x[1 : (length(x) -1)]}
    )$output1
    dims <- dim(rcp_seasonal_mean)
    time_dim <- which(names(dim(rcp_seasonal_mean)) == "time")
    dims <- append(dims, c(3, dims[time_dim] / 3), after = time_dim)
    dims <- dims[-time_dim]
    dim(rcp_seasonal_mean) <- dims
    names(dim(rcp_seasonal_mean))[c(time_dim, time_dim + 1)] <- c("season", "year")
    rcp_seasonal_mean <- Mean1Dim(rcp_seasonal_mean, posdim = time_dim)
    rcp_seasonal_mean <- aperm(rcp_seasonal_mean, c(2, 1, 3, 4))
  }
}
years_dim <- which(names(dim(rcp_seasonal_mean)) == "year")
climatology <- InsertDim(
  climatology,
  years_dim,
  lendim = dim(rcp_seasonal_mean)[years_dim]
)
anomaly <- rcp_seasonal_mean - climatology
multi_year_anomaly <- Mean1Dim(anomaly, years_dim)

time <- seq(start_anomaly, end_anomaly, by = 1)
month <- moninf
  if (month <= 9) {
    month <- paste0(as.character(0), as.character(month))
  }
  month <- paste0("-", month, "-")
  day <- "01"
  time <- as.POSIXct(paste0(time, month, day), tz = "CET")
  time <- julian(time, origin = as.POSIXct("1970-01-01"))

attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(
  standard_name = "time",
  long_name = "time",
  units = "days since 1970-01-01 00:00:00",
  prec = "double",
  dim = list(list(name="time", unlim = FALSE))
))
attr(time, "variables") <- metadata

#Save the single model anomalies
for (mod in 1 : length(model_names)) {
    data <- anomaly[mod, , , ]
  #data <- anomaly[mod,1,1, , ,]
  data <- aperm(data, c(2, 3, 1))
  names(dim(data)) <- c("lat", "lon", "time")
  metadata <- list(variable = list(
    dim = list(list(name="time", unlim = FALSE)),
    units = units
  ))
  names(metadata)[1] <- var0
  attr(data, "variables") <- metadata
  attributes(lat) <- NULL
  attributes(lon) <- NULL
  dim(lat) <- c(lat = length(lat))
  dim(lon) <- c(lon = length(lon))
  variable_list <- list(variable = data, lat = lat, lon = lon, time = time)
  names(variable_list)[1] <- var0

  ArrayToNetCDF(
    variable_list,
    paste0(
      plot_dir,  "/", var0, "_", months, "_anomaly_", model_names[mod],
      "_", start_anomaly, "_", end_anomaly,"_", start_climatology, "_",
      end_climatology, ".nc"
    )
  )
}

model_anomalies <- WeightedMean(
  anomaly,
  lon = as.vector(lon),
  lat = as.vector(lat),
  mask = NULL
)
if (!is.null(params$running_mean)) {
  model_anomalies <- Smoothing(
    model_anomalies,
    runmeanlen = params$running_mean,
    numdimt = 2
  )
}
data_frame <- as.data.frame.table(t(model_anomalies[,]))
years <- rep(start_anomaly : end_anomaly, dim(model_anomalies)[1])
data_frame$Year <- c(years)
names(data_frame)[2] <- "Model"

for (i in 1 : length(levels(data_frame$Model))) {
    levels(data_frame$Model)[i] <- model_names[i]
}

if (!is.null(time_series_plot)) {
  g <- ggplot(
    data_frame,
    aes(x = Year, y = Freq, color = Model)) + theme_bw() +
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
ggsave(
  filename = paste0(
    plot_dir, "/", "Area-averaged ", var0, "_",months, "_multimodel-anomaly_",
    start_anomaly, "_", end_anomaly,"_", start_climatology, "_",
    end_climatology, ".png"),
  g,
  device = NULL
)

if (!is.null(agreement_threshold)) {
  model_dim <- which(names(dim(multi_year_anomaly)) == "model")
  agreement <- AnoAgree(
    multi_year_anomaly + rnorm(length(unique(model_names)) * length(lat) * length(lon)),
    membersdim = model_dim
  )
} else {
  agreement_threshold <- 1000
  agreement <- NULL
}

colorbar_lim <- ceiling(max(abs(max(multi_year_anomaly)),abs(min(data))))
brks <- seq(-colorbar_lim, colorbar_lim, length.out = 21)
title <- paste0(
  months, " ", var0, " anomaly (", start_anomaly, "-", end_anomaly,
  ") - (", start_climatology, "-", end_climatology, ")"
)
data <- drop(Mean1Dim(multi_year_anomaly, model_dim))

PlotEquiMap(
  data,
  lat = lat,
  lon = lon,
  brks = brks,
  units = units,
  toptitle = title,
  filled.continents = FALSE,
  dots = drop(agreement) >= agreement_threshold,
  fileout = paste0(
    plot_dir, "/", var0, "_",months, "_multimodel-anomaly_", start_anomaly,
    "_", end_anomaly,"_", start_climatology, "_", end_climatology, ".png"
  )
)
model_names_filename <- paste(model_names, collapse = "_")
print(paste(
  "Attribute projection from climatological data is saved and,",
  "if it's correct, it can be added to the final output:",
  projection
))

dimlon <- ncdim_def(
  name = "lon",
  units = "degrees_east",
  vals = as.vector(lon),
  longname = "longitude"
)
dimlat <- ncdim_def(
  name = "lat",
  units = "degrees_north",
  vals = as.vector(lat),
  longname = "latitude"
)
defdata <- ncvar_def(
  name = "data",
  units = units,
  dim = list(lat = dimlat, lon = dimlon),
  longname = paste("Mean",long_names)
)
defagreement <- ncvar_def(
  name = "agreement",
  units = "%",
  dim = list(lat = dimlat, lon = dimlon),
  longname = "Agremeent between models")

file <- nc_create(
  paste0(
    plot_dir, "/", var0, "_",months, "_multimodel-anomaly_",
    model_names_filename,"_", start_anomaly, "_", end_anomaly,"_",
    start_climatology, "_", end_climatology, ".nc"),
  list(defdata, defagreement)
)
ncvar_put(file, defdata, data)
ncvar_put(file, defagreement, agreement)
nc_close(file)
