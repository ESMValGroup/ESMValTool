library(s2dverification)
library(ggplot2)
library(multiApply) #nolint
library(ncdf4)
library(gridExtra) #nolint
library(ClimProjDiags) #nolint
library(yaml)

initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "", initial_options[grep(file_arg_name,
    initial_options)])
script_dirname <- dirname(script_name)
source(file.path(script_dirname, "WeatherRegime.r"))
source(file.path(script_dirname, "RegimesAssign.r"))

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

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

input_files_per_var <- yaml::read_yaml(params$input_files)


model_names <- lapply(input_files_per_var, function(x) x$dataset)
model_names <- unique(unlist(unname(model_names)))


var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]

experiment <- lapply(input_files_per_var, function(x) x$exp)
experiment <- unlist(unname(experiment))

reference_files <- which(unname(experiment) == "historical")
projection_files <- which(unname(experiment) != "historical")


#Region considered to select the plot
region <- params$region


#Start and end periods for the historical and projection periods
start_historical <- as.POSIXct(params$start_historical)
end_historical <- as.POSIXct(params$end_historical)
start_projection <- as.POSIXct(params$start_projection)
end_projection <- as.POSIXct(params$end_projection)

#Regime parameters
ncenters <- params$ncenters
cluster_method <- params$cluster_method
if (cluster_method != "kmeans") {
    cluster_method <- "complete"
}
EOFS <- params$EOFS
frequency <- params$frequency
detrend_order <- params$detrend_order

if (grepl("OImon", fullpath_filenames[1])) {
  var_lat <- "rlat"
  var_lon <- "rlon"
} else {
  var_lat <- "lat"
  var_lon <- "lon"
}
# ---------------------------
# Reading and formating
# ---------------------------
ref_nc <- nc_open(fullpath_filenames[reference_files])
var0 <- unlist(var0)
reference_data <- ncvar_get(ref_nc, var0)
print(str(reference_data))
names(dim(reference_data)) <- rev(names(ref_nc$dim))[-1]
lat <- ncvar_get(ref_nc, var_lat)
lon <- ncvar_get(ref_nc, var_lon)
units <- ncatt_get(ref_nc, var0, "units")$value
calendario <- ncatt_get(ref_nc, "time", "calendar")$value
long_names <-  ncatt_get(ref_nc, var0, "long_name")$value
time <-  ncvar_get(ref_nc, "time")
start_date <- as.POSIXct(substr(ncatt_get(ref_nc, "time",
                                          "units")$value, 11, 29))
nc_close(ref_nc)

data_type <- ifelse(grepl("day", fullpath_filenames[1]), "day", "month")
dates_historical <- seq(start_historical, end_historical, data_type)

print(calendario)
if (dim(reference_data)[3] != length(dates_historical)) {
  if (
    calendario == "365" | calendario == "365_days" |
    calendario == "365_day" | calendario == "noleap"
  ) {
dates_historical <-
dates_historical[-which(substr(dates_historical, 6, 10) == "02-29")]#nolint
  }
}

reference_data <- as.vector(reference_data)
dim(reference_data) <- c(
  model = 1,
  var = 1,
  lon = length(lon),
  lat = length(lat),
  time = length(dates_historical)
)

reference_data <- aperm(reference_data, c(1, 2, 5, 3, 4))
attr(reference_data, "Variables")$dat1$time <- dates_historical


names(dim(reference_data)) <- c("model", "var", "time", "lon", "lat")
time_dimension <- which(names(dim(reference_data)) == "time")


# -------------------------------
## Selecting the season or month
# -------------------------------
time_dim <- which(names(dim(reference_data)) == "time")

months <- c(
  "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP",
  "OCT", "NOV", "DEC"
)
seasons <- c("DJF", "MAM", "JJA", "SON")
mes <- match(frequency, months)
sea <- match(frequency, seasons)
print(mes)
print(frequency)
if (!is.na(mes)) {
  print("MONTHLY")
  dims <- dim(reference_data)
  ind <- which(as.numeric(substr(dates_historical, 6, 7)) == mes)
  years <- unique(as.numeric(substr(dates_historical, 1, 4)))
  reference_data <- reference_data[ , , ind , , ] #nolint
    print(dims)
  dims <- append(
    dims, c(length(ind) / length(years), length(years)), after = time_dim
  )
    print(dims)
} else if (!is.na(sea)) {
  print("Seasonal")
  reference_data <- SeasonSelect( #nolint
    reference_data,
    season = frequency,
    dates = dates_historical,
    calendar = calendario
  )
  time <- reference_data$dates
  years <- unique(as.numeric(substr(time, 1, 4)))
  reference_data <- reference_data$data
  reference_data <- InsertDim(reference_data, posdim = 1, lendim = 1) #nolint
  reference_data <- InsertDim(reference_data, posdim = 1, lendim = 1) #nolint
  names(dim(reference_data))[c(1, 2)] <- c("model", "var")
  dims <- dim(reference_data)
dims <- append(
    dims, c(length(time) / length(years), length(years)), after = time_dim
  )
 }


dims <- dims[-time_dim]
dim(reference_data) <- dims

names(dim(reference_data))[c(time_dim, time_dim + 1)] <- c("sdate", "ftime")

Loess <- function(clim, loess_span) {
  if (sum(is.na(clim)) != length(clim)) {
  data <- data.frame(ensmean = clim, day = 1 : length(clim))
  loess_filt <- loess(
    ensmean ~ day, data, span = loess_span, degree = detrend_order
  )
  output <- predict(loess_filt)
  } else {
   output <- rep(NA, length(clim))
  }
  return(output)
}
# -------------------------------
## Computing the WR_obs
# -------------------------------

clim_obs <- array(
  apply(reference_data, c(1, 2, 3, 5, 6), mean),
  dim = dim(reference_data)[-4]
)
print(sea)
print(data_type)
save(clim_obs, clim_obs, file = file.path(plot_dir, "ClimObs.RData"))
print(data_type == "day" | !is.na(sea))
if (data_type == "day" | !is.na(sea)) {
clim_obs <- aperm(
  apply(
      clim_obs,
      c(1 : length(dim(clim_obs)))[-which(names(dim(clim_obs)) == "sdate")],
      Loess,
      loess_span = 1),
  c(2, 3, 1, 4, 5)
)
}

anom_obs <- Ano(reference_data, clim_obs)
print(head(anom_obs[which(anom_obs > 0)]))
print(dim(anom_obs))
WR_obs <- WeatherRegime( #nolint
  data = anom_obs,
  EOFS = EOFS,
  lat = lat,
  lon = lon,
  ncenters = ncenters,
  method = cluster_method
)
print(str(WR_obs))
print(cluster_method)
if (cluster_method != "kmeans" & data_type == "day") {
    WR_obs$composite <- aperm(WR_obs$composite, c(2, 1, 3, 4, 5))
}
names(dim(WR_obs$composite)) <- c("lat", "lon", "Cluster", "Mod", "exp")
names(dim(WR_obs$cluster))[1] <- "Evolution"
 print(dim(WR_obs$composite))
# -------------------------------
## Plotting the WR_obs output
# -------------------------------

if (cluster_method == "kmeans") {
    clim_frequencies <- paste0(
    "freq = ",
    round(Mean1Dim(WR_obs$frequency, 1), 1), "%")
} else {
    freq_clus <- NULL
    for (i in 1 : ncenters) {
        freq_clus <- c(freq_clus, round(sum(WR_obs$cluster[, 1, 1] == i) /
            length(WR_obs$cluster[ , 1, 1]) * 100, 1))
    }
    clim_frequencies <- paste0(
    "freq = ", freq_clus, "%")
    WR_obs$frequency <- freq_clus
}

cosa <- aperm(drop(WR_obs$composite), c(3, 1, 2))
print(dim(cosa))
lim <- max(abs(cosa / 100), na.rm = TRUE)
print(lim)
if (lim < 1) {
    x <- floor(log10(lim)) + 1
    lim <- 10 ^ x
} else {
    lim <- ceiling(lim)
}

if (region == "Polar") {
    filepng1 <- paste0(
        plot_dir, "/", frequency, "-", var0, "_observed_regimes.png")
    title <- paste0(paste0("Cluster ", 1 : 4), " (", clim_frequencies, " )")
  PlotLayout( #nolint
    PlotStereoMap, #nolint
    c(2, 3),
    lon = lon,
    lat = lat,
    var = cosa / 100,
    titles = title,
    filled.continents = FALSE,
    axelab = FALSE,
    draw_separators = TRUE,
    subsampleg = 1,
    brks = seq(-1 * lim, lim, by = lim / 10),
    fileout = filepng1)
} else {
    filepng1 <- paste0(
      plot_dir, "/", frequency, "-", var0, "_observed_regimes.png")
    title <- paste0(paste0("Cluster ", 1 : 4), " (", clim_frequencies, " )")
  PlotLayout( #nolint
    PlotEquiMap, #nolint
    c(2, 3),
    lon = lon,
    lat = lat,
    var = cosa / 100,
    titles = title,
    filled.continents = FALSE,
    axelab = FALSE,
    draw_separators = TRUE,
    subsampleg = 1,
    brks = seq(-1 * lim, lim, by = lim / 10),
    fileout = filepng1)
}

# -------------------------------
## Save the WR_obs output to ncdf
# -------------------------------
time <- dates_historical
time <- julian(time, origin = as.POSIXct("1970-01-01"))
attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(
  standard_name = "time", long_name = "time",
  units = "days since 1970-01-01 00:00:00", prec = "double",
  dim = list(list(name = "time", unlim = FALSE))))
attr(time, "variables") <- metadata

attributes(lon) <- NULL
attributes(lat) <- NULL
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))
metadata <- list(variable = list(dim = list(list(name = "time",
                                                 unlim = FALSE))))
dim(WR_obs$frequency) <- c(frequency = length(WR_obs$frequency))
dim(WR_obs$pvalue) <- c(pvalue = length(WR_obs$pvalue))
dim(WR_obs$cluster) <- c(cluster = length(WR_obs$cluster))
print("dimensions")
print(dim(WR_obs$composite))
print(dim(WR_obs$pvalue))
print(dim(WR_obs$cluster))
print(dim(WR_obs$frequency))
print(dim(lat))
print(cluster_method)
if (cluster_method != "kmeans") {
print("IN")
 names(dim(WR_obs$composite))[c(1, 2)] <- c("lon", "lat")
print(dim(WR_obs$composite))
}
variable_list <- list(
  variable = WR_obs$composite,
  pvalue = WR_obs$pvalue,
  cluster = WR_obs$cluster,
  frequency = WR_obs$frequency,
  lat = lat,
  lon = lon,
  time = time)

names(variable_list)[1] <- var0
attributes(variable_list) <- NULL
filencdf1 <-   paste0(
    work_dir, "/", var0, "_", frequency, "_WR_obs_", model_names, "_",
    start_projection, "_", end_projection, "_", start_historical, "_",
    end_historical, ".nc")

ArrayToNetCDF( #nolint
  variable_list, filencdf1)
   # Set provenance for output files
    xprov <- list(ancestors = list(fullpath_filenames[reference_files]),
                  authors = list("torr_ve", "fuck_ne","cort_ni", "guem_vi",
                                 "hunt_al", "manu_ni"),
                  projects = list("c3s-magic"),
                  caption = "Observed modes of variability",
                  statistics = list("eof"),
                  ncenters = params$ncenters,
                  cluster_method = cluster_method,
                  EOFS = params$EOFS,
                  frequency = params$frequency,
                  detrend_order = params$detrend_order,
                  realms = list("atmos"),
                  themes = list("phys"),
                  plot_file = filepng1)

      provenance[[filencdf1]] <- xprov
# ---------------------------
# ---------------------------
# Reading and formating
# ---------------------------
proj_nc <- nc_open(fullpath_filenames[projection_files])
projection_data <- ncvar_get(proj_nc, var0)
names(dim(projection_data)) <- rev(names(proj_nc$dim))[-1]
time <-  ncvar_get(proj_nc, "time")
start_date <- as.POSIXct(substr(ncatt_get(proj_nc, "time",
                                          "units")$value, 11, 29))
nc_close(proj_nc)

dates_projection <- seq(start_projection, end_projection, data_type)

print(calendario)
if (dim(projection_data)[3] != length(dates_projection)) {
  if (
    calendario == "365" | calendario == "365_days" |
    calendario == "365_day" | calendario == "noleap"
  ) {
dates_projection <-
dates_projection[-which(substr(dates_projection, 6, 10) == "02-29")]#nolint
  }
}



data <- as.vector(projection_data)
dim(projection_data) <- c(
  model = 1,
  var = 1,
  lon = length(lon),
  lat = length(lat),
  time = length(time)
)
print(dim(projection_data))
projection_data <- aperm(projection_data, c(1, 2, 5, 3, 4))
attr(projection_data, "Variables")$dat1$time <- dates_projection
names(dim(projection_data)) <- c("model", "var", "time", "lon", "lat")
time_dimension <- which(names(dim(projection_data)) == "time")

# ---------------------------
# Selecting the period
# ---------------------------
time_dim <- which(names(dim(projection_data)) == "time")

if (!is.na(mes)) {
  print("MONTHLY")
  dims <- dim(projection_data)
  ind <- which(as.numeric(substr(dates_projection, 6, 7)) == mes)
  years <- unique(as.numeric(substr(dates_projection, 1, 4)))
  projection_data <- projection_data[ , , ind , , ] #nolint
  dims <- append(
    dims,
    c(length(ind) / length(years), length(years)),
    after = time_dim)
} else if (!is.na(sea)) {
  projection_data <- SeasonSelect( #nolint
  projection_data,
  season = frequency,
  dates = dates_projection,
  calendar = calendario
  )
  time <- projection_data$dates
  years <- unique(as.numeric(substr(time, 1, 4)))
  projection_data <- projection_data$data
  projection_data <- InsertDim(projection_data, posdim = 1, lendim = 1)#nolint
  projection_data <- InsertDim(projection_data, posdim = 1, lendim = 1)#nolint

  names(dim(projection_data))[c(1, 2)] <- c("model", "var")
  dims <- dim(projection_data)
  dims <- append(
    dims,
    c(length(time) / length(years), length(years)),
    after = time_dim
  )
}
dims <- dims[-time_dim]
dim(projection_data) <- dims
names(dim(projection_data))[c(time_dim, time_dim + 1)] <- c("sdate", "ftime")

clim_ref <- array(
  apply(projection_data, c(1, 2, 3, 5, 6), mean),
  dim = dim(projection_data)[-4]
)


if (data_type == "day" | !is.na(sea)) {
clim_ref <- aperm(
  apply(
    clim_ref,
    c(1 : length(dim(clim_ref)))[-which(names(dim(clim_ref)) == "sdate")],
    Loess,
    loess_span = 1),
    c(2, 3, 1, 4, 5))
}
print(dim(clim_ref))
print("NOHERE")
anom_exp <- Ano(projection_data, clim_ref)
print(dim(anom_exp))
reference <- drop(WR_obs$composite)
print(dim(reference))
if (cluster_method == "kmeans") {
  names(dim(reference)) <- c("lat", "lon", "nclust")

  reference <- aperm(reference, c(3, 2, 1))
  names(dim(reference)) <- c("nclust", "lon", "lat")
} else {
  names(dim(reference)) <- c("lon", "lat", "nclust")
}
print(head(reference[which(reference > 0)]))
print(head(anom_exp[which(anom_exp > 0)]))

WR_exp <- RegimesAssign( #nolint
  var_ano = anom_exp, ref_maps = reference, lats = lat, method = "distance"
)


# ---------------------------
# Plotting WR projection:
# ---------------------------

if (cluster_method == "kmeans") {
  cosa <- aperm(WR_exp$composite, c(3, 2, 1))
  names(dim(WR_exp$composite))[3] <- "nclust"
} else {
  names(dim(WR_exp$composite))[3] <- "nclust"
}
print("AAET")
print(dim(cosa))
print(dim(WR_exp$composite))
lim <- max(abs(cosa / 100), na.rm = TRUE)
if (lim < 1) {
    x <- floor(log10(lim)) + 1
    lim <- 10 ^ x
} else {
    lim <- ceiling(lim)
}
if (region == "Polar") {
    filepng2 <- paste0(
      plot_dir, "/", frequency, "-", var0, "_predicted_regimes.png")
    title <- paste0(paste0("Cluster ", 1 : dim(cosa)[1], " (",
        paste0("freq = ", round(WR_exp$frequency, 1), "%"), " )"))
  PlotLayout( #nolint
    PlotStereoMap, #nolint
    c(2, 3),
    lon = lon,
    lat = lat,
    var = cosa / 100,
    titles = title,
    filled.continents = FALSE,
    draw_separators = TRUE, subsampleg = 1,
    brks = seq(-1 * lim, lim, by = lim / 10),
    fileout = filepng2)
} else {
    filepng2 <- paste0(
      plot_dir, "/", frequency, "-", var0, "_predicted_regimes.png")
    title <- paste0(paste0("Cluster ", 1 : dim(cosa)[1], " (",
        paste0("freq = ", round(WR_exp$frequency, 1), "%"), " )"))
  PlotLayout( #nolint
    PlotEquiMap, #nolint
    c(2, 3),
    lon = lon,
    lat = lat,
    var = cosa / 100,
    titles = title,
    filled.continents = FALSE,
    axelab = FALSE, draw_separators = TRUE, subsampleg = 1,
    brks = seq(-1 * lim, lim, by = lim / 10),
    fileout = filepng2)
}


# -------------------------------
## Save the WR_exp output to ncdf
# -------------------------------
time <- dates_projection
time <- julian(time, origin = as.POSIXct("1970-01-01"))
attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = "time", long_name = "time",
                units = "days since 1970-01-01 00:00:00", prec = "double",
                dim = list(list(name = "time", unlim = FALSE))))
attr(time, "variables") <- metadata

attributes(lon) <- NULL
attributes(lat) <- NULL
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))
metadata <- list(variable = list(dim = list(list(name = "time",
                                                unlim = FALSE))))

dim(WR_exp$frequency) <- c(frequency = length(WR_exp$frequency))
dim(WR_exp$pvalue) <- c(pvalue = length(WR_exp$pvalue))
dim(WR_exp$cluster) <- c(cluster = length(WR_exp$cluster))

variable_list <- list(
  variable = WR_exp$composite,
  pvalue = WR_exp$pvalue,
  cluster = WR_exp$cluster,
  frequency = WR_exp$frequency,
  lat = lat,
  lon = lon,
  time = time)
names(variable_list)[1] <- var0

attributes(variable_list) <- NULL
filencdf2 <-   paste0(
    work_dir, "/", var0, "_", frequency, "_WR_exp_", model_names, "_",
    start_projection, "_", end_projection, "_", start_historical, "_",
    end_historical, ".nc")
ArrayToNetCDF(variable_list, filencdf2) #nolint

   # Set provenance for output files
    xprov <- list(ancestors = list(filencdf1,
                              fullpath_filenames[projection_files]),
                  authors = list("torr_ve", "fuck_ne","cort_ni", "guem_vi",
                                 "hunt_al", "manu_ni"),
                  projects = list("c3s-magic"),
                  caption = "Predicted modes of variability",
                  statistics = list("other"),
                  ncenters = params$ncenters,
                  cluster_method = cluster_method,
                  EOFS = params$EOFS,
                  frequency = params$frequency,
                  detrend_order = params$detrend_order,
                  realms = list("atmos"),
                  themes = list("phys"),
                  plot_file = filepng2)

      provenance[[filencdf2]] <- xprov
# ---------------------------
# Computing the RMSE:
# ---------------------------
 if (cluster_method != "kmeans") {
 print("DOME")
  cosa <- aperm(cosa, c(2, 3, 1))
  reference <- aperm(reference, c(3, 1, 2))
} else {
  cosa <- aperm(WR_exp$composite, c(2, 1, 3))
}
print("JEH")
print(dim(cosa))
print(dim(reference))
rmse <- NULL
for (i in 1 : ncenters) {
  for (j in 1 : ncenters) {
      rmse <- c(rmse, sqrt(mean( (reference[i, , ] - cosa[, , j]) ^ 2, #nolint
                na.rm = T)))
  }
}
dim(rmse) <- c(ncenters, ncenters)
print(rmse)

dimpattern <- ncdim_def(
  name = "pattern",
  units = "undim",
  vals = 1 : ncenters,
  longname = "Pattern"
)
title <- paste0(
    "Root Mean Squared Error between observed and ",
    "future projected patterns")
defrmse <- ncvar_def(
  name = "rmse",
  units = "undim",
  dim = list(observed = dimpattern, experiment = dimpattern),
  longname = title)
filencdf3 <- paste0(work_dir, "/", var0, "_", frequency, "_rmse_",
                         model_names, "_", start_projection, "_",
                         end_projection, "_", start_historical, "_",
                         end_historical, ".nc")
file <- nc_create(filencdf3, list(defrmse))
ncvar_put(file, defrmse, rmse)

nc_close(file)

colnames(rmse) <- paste("Obs", 1 : ncenters)
rownames(rmse) <- paste("Pre", 1 : ncenters)
filepng3 <- paste0(file.path(plot_dir, "Table_"), var0, "_", frequency,
            "_rmse_", model_names,
           "_", start_projection, "_", end_projection, "_", start_historical,
           "_", end_historical, ".png")
png(filepng3, height = 6, width = 18, units = "cm", res = 100)
grid.table(round(rmse, 2))
dev.off()

    # Set provenance for output files
    xprov <- list(ancestors = list(filencdf1, filencdf2),
                  authors = list("torr_ve", "fuck_ne","cort_ni", "guem_vi",
                                 "hunt_al", "manu_ni"),
                  projects = list("c3s-magic"),
                  caption = title,
                  statistics = list("rmsd"),
                  ncenters = params$ncenters,
                  cluster_method = cluster_method,
                  EOFS = params$EOFS,
                  frequency = params$frequency,
                  detrend_order = params$detrend_order,
                  realms = list("atmos"),
                  themes = list("phys"),
                  plot_file = filepng3)

      provenance[[filencdf3]] <- xprov

# Write provenance to file
write_yaml(provenance, provenance_file)
