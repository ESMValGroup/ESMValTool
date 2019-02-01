# nolint start
####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# nolint end
Sys.setenv(TAR = "/bin/tar") # nolint
library(s2dverification)
library(multiApply) # nolint
library(ggplot2)
library(yaml)
library(ncdf4)
library(ClimProjDiags) #nolint
library(abind)

#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
print(params)
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir

## Create working dirs if they do not exist
dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)


# setup provenance file and list
provenance_file <- paste0(plot_dir, "/", "combined_provenance.yml")
provenance <- list()

input_files_per_var <- yaml::read_yaml(params$input_files[1])
var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x) x$dataset)
model_names <- unique(unlist(unname(model_names)))

var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]

var0 <- unlist(var0)


start_year <- lapply(input_files_per_var, function(x) x$start_year)
start_year <- c(unlist(unname(start_year)))[1]
end_year <- lapply(input_files_per_var, function(x) x$end_year)
end_year <- c(unlist(unname(end_year)))[1]

#Parameters for Season() function
monini <- 1
moninf <- params$moninf
monsup <- params$monsup
region <- params$region
print("MESES")
print(moninf)
print(monsup)
mask <- NULL ### How can we pass masks from ESMValTool?
running_mean <- params$running_mean
multi_year_average <- params$multi_year_average
weights <- params$weights
### Load data
print(running_mean)
print(multi_year_average)
print(weights)
data_nc <- nc_open(fullpath_filenames)
lat <- ncvar_get(data_nc, "lat")
lon <- ncvar_get(data_nc, "lon")
units <- ncatt_get(data_nc, var0, "units")$value
calendar <- ncatt_get(data_nc, "time", "calendar")$value
long_names <-  ncatt_get(data_nc, var0, "long_name")$value
time <-  ncvar_get(data_nc, "time")
print(dim(data_nc))
data <- InsertDim(ncvar_get(data_nc, var0), 1, 1) # nolint
print(dim(data))
start_date <- as.POSIXct(substr(ncatt_get(data_nc, "time",
                                          "units")$value, 11, 29))
time <- as.Date(time, origin = start_date, calendar = calendar)
projection <- "NULL"
nc_close(data_nc)
if (length(params$input_files) >= 2) {
  for (i in 2 : length(params$input_files)) {
    input_files_per_var <- yaml::read_yaml(params$input_files[i])
    var_names <- names(input_files_per_var)
    model_names <- lapply(input_files_per_var, function(x) x$dataset)
    model_names <- unique(unlist(unname(model_names)))
    var0 <- lapply(input_files_per_var, function(x) x$short_name)
    fullpath_filenames <- names(var0)
    var0 <- unname(var0)[1]
    var0 <- unlist(var0)
    data_nc <- nc_open(fullpath_filenames)
    data <- abind(data,
    InsertDim(ncvar_get(data_nc, var0), 1, 1), along = 1) # nolint
    nc_close(data_nc)
  }
}
names(dim(data)) <- c("model", "lon", "lat", "time")
region <- c(min(lon), max(lon), min(lat), max(lat))
attributes(lon) <- NULL
attributes(lat) <- NULL
dim(lon) <-  c(lon = length(lon))
dim(lat) <- c(lat = length(lat))
time_dim <- which(names(dim(data)) == "time")
timestamp <- ""
attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(
  standard_name = "time",
  long_name = "time",
  units = "days since 1970-01-01 00:00:00",
  prec = "double",
  dim = list(list(name = "time", unlim = FALSE))
))
attr(time, "variables") <- metadata
if (!is.null(region)) {
  dim_names <- names(dim(data))
  londim <- which(names(dim(data)) == "lon")
  latdim <- which(names(dim(data)) == "lat")
  data <- WeightedMean( # nolint
    data,
    lon = as.vector(lon),
    lat = as.vector(lat),
    region = region,
    mask = NULL
  )
  names(dim(data)) <- dim_names[-c(londim, latdim)]
  time_dim <- which(names(dim(data)) == "time")
}

if (!is.null(running_mean)) {
    data <- Smoothing(data, runmeanlen = running_mean, numdimt = time_dim)
    timestamp <- paste0(running_mean, "-month-running-mean-")
}

print(paste("moninf", moninf))
if (!is.null(moninf)) {
  months <- paste0(month.abb[moninf], "-", month.abb[monsup])
  print(months)
  dims <- dim(data)
  dims <- append(dims, c(12, dims[time_dim] / 12), after = time_dim)
  dims <- dims[-time_dim]
  dim(data) <- dims
  names(dim(data))[c(time_dim, time_dim + 1)] <- c("month", "year")
  margins <- list(c(1 : length(dim(data)))[-c(time_dim + 1)])
  data <- Season(
    data,
    posdim = time_dim,
    monini = monini,
    moninf = moninf,
    monsup = monsup
  )
  margins <- list(c(1 : length(dim(data)))[-c(time_dim + 1)])
  years_dim <- which(names(dim(data)) == "year")
  timestamp <- "Mar-May"
  month_dim <- which(names(dim(data)) == "month")
  data <- adrop(data, month_dim)
  if (multi_year_average == "TRUE") {
    time <- time[1]
    data <- Mean1Dim(data, c(years_dim - 1))
  }
}

if (!is.null(weights)) {
  indices_dim <- which(names(dim(data)) == "model")
  indices <- list()
  print(indices_dim)
  for (i in 1 : dim(data)[indices_dim]) {
    indices[[i]] <- Subset(data, along = indices_dim, indices = i)
  }
    print(str(indices))
  if (!is.numeric(weights)) {
    weights <- "NULL"
    print("AQUI")
    data <- CombineIndices(indices, weights = NULL) # nolint
    print(dim(data))
  } else {
    data <- CombineIndices(indices, weights = weights) # nolint
  }
}
print(region)
if (!is.null(region)) {
  data <- data[1,  ]
  attributes(data) <- NULL
  dim(data) <-  c(time = length(data))
  metadata <- list(
    index = list(
      dim = list(list(name = "time", unlim = FALSE, prec = "double"))
    )
  )
  names(metadata)[1] <- var0
  attr(data, "variables") <- metadata
  variable_list <- list(variable = data, time = time)
  names(variable_list)[1] <- var0
  # nolint start
  # ArrayToNetCDF(
  #   variable_list,
  #   paste0(
  #     plot_dir, "/", var0, "_", paste(model_names, sep="", collapse="_"),
  #     "_", timestamp, "_", rcp_scenario, "_", start_year, "_", end_year,
  #     "_", ".nc"
  #   )
  # )
  # nolint end
 model_names_filename <- paste(model_names, collapse = "_")

  print(
    paste(
      "Attribute projection from climatological data is saved and,",
      "if it's correct, it can be added to the final output:",
      projection
    )
  )
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
  dimtime <- ncdim_def(
    name = "time",
    units = "days since 1970-01-01 00:00:00",
    vals = as.vector(time),
    longname = "time"
  )
  defdata <- ncvar_def(
    name = "data",
    units = units,
    dim = list(time = dimtime),
    longname = paste("Combination", long_names)
  )

  file <- nc_create(
    paste0(
      plot_dir, "/", var0, "_", paste0(model_names, collapse = "_"),
      "_", timestamp, "_", model_names_filename, "_", start_year,
      "_", end_year, "_", ".nc"),
    list(defdata)
  )
  ncvar_put(file, defdata, data)
  nc_close(file)
} else {
  data <- data[1, 1, , , ] # nolint
  data <- aperm(data, c(3, 2, 1))
  names(dim(data)) <- c("lon", "lat", "time")
  metadata <- list(
    index = list(
      dim = list(list(name = "time", unlim = FALSE, prec = "double"))))
  names(metadata)[1] <- var0
  attr(data, "variables") <- metadata
  variable_list <- list(variable = data, lat = lat, lon = lon, time = time)
  names(variable_list)[1] <- var0

  print(paste(
    "Attribute projection from climatological data is saved and,",
    "if it's correct, it can be added to the final output:",
    projection)
  )

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
  dimtime <- ncdim_def(
    name = "time",
    units = "days since 1970-01-01 00:00:00",
    vals = as.vector(time),
    longname = "time"
  )
  defdata <- ncvar_def(
    name = "data",
    units = units,
    dim = list(time = dimtime),
    longname = paste("Combination", long_names)
  )
  file <- nc_create(
    paste0(
      plot_dir, "/", var0, "_", paste0(model_names, collapse = "_"),
      "_", timestamp, "_", model_names_filename, "_", start_year, "_",
      end_year, "_", ".nc"),
    list(defdata)
  )
  ncvar_put(file, defdata, data)
  nc_close(file)
}
