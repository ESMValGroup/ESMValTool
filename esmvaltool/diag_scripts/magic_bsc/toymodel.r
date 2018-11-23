####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4
#install.packages("yaml")
#install.packages("devtools")
#library(devtools)
#Sys.setenv(TAR = "/bin/tar")
#install_git("https://earth.bsc.es/gitlab/es/startR", branch = "develop-hotfixes-0.0.2")
#install_git("https://earth.bsc.es/gitlab/es/easyNCDF", branch = "master")
Sys.setenv(TAR = "/bin/tar")
library(s2dverification)
library(ncdf4)
library(multiApply)
library(yaml)
library(abind)
library(ClimProjDiags)
library(RColorBrewer)

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
a <- 1
b <- params$beta
g <- 0.1
nm <- params$number_of_members
nstartd <- 1
nleadt = params$no_of_lead_times

#fullpath_filenames <- "/scratch/Earth/ahunter/esmvaltool_input/CMIP5/Tier2/NCEP/OBS_NCEP_reanaly_1_T3M_ta_200001-200212.nc"
var0 <- unlist(var0)
data_nc <- nc_open(fullpath_filenames)
data <- ncvar_get(data_nc, var0)
data <- InsertDim(InsertDim(data, 1, 1),1,1)
names(dim(data)) <- c("model", "var","lon", "lat", "time")
lat <- ncvar_get(data_nc,"lat")
lon <- ncvar_get(data_nc,"lon")
lon <- unlist(lon)
lat <- unlist(lat)
print(lon)
print(lat)
attributes(lon) <- NULL
attributes(lat) <- NULL
units <- ncatt_get(data_nc, var0, "units")$value
calendar <- ncatt_get(data_nc, "time", "calendar")$value
long_names <-  ncatt_get(data_nc,var0,"long_name")$value
time <-  ncvar_get(data_nc,"time")
start_date <- as.POSIXct(substr(ncatt_get(data_nc, "time", "units")$value,11, 29 ))
nc_close(data_nc)
time <- as.Date(time, origin = start_date, calendar = calendar)

dim_names <- names(dim(data))
lon_dim <- which(names(dim(data)) == "lon")
lat_dim <- which(names(dim(data)) == "lat")
data <- WeightedMean(data, lat = lat, lon = lon, londim = lon_dim, latdim = lat_dim)
names(dim(data)) <- dim_names[-c(lon_dim, lat_dim)]
time_dim <- which(names(dim(data)) == "time")

ToyModel <- function (
  alpha = 0.1, beta = 0.4, gamma = 1, sig = 1, trend = 0,
  nstartd = 30, nleadt = 4, nmemb = 10, obsini = NULL, fxerr = NULL
){
  if (any(!is.numeric(c(alpha, beta, gamma, sig, trend, nstartd,
      nleadt, nmemb)))) {
      stop(paste("Parameters alpha, beta, gamma, sig, trend, nstartd,",
          "nleadt and nmemb must be numeric."))
  }
  nstartd <- round(nstartd)
  nleadt <- round(nleadt)
  nmemb <- round(nmemb)
  if (!is.null(obsini)) {
    if (!is.numeric(obsini) || !is.array(obsini)) {
      stop("Parameter obsini must be a numeric array.")
    }
    if (length(dim(obsini)) != 4) {
      stop(paste(
        "Parameter obsini must be an array with dimensions",
        "c(1, 1, nleadt, nstartd)."
      ))
    }
    if (dim(obsini)[3] != nstartd || dim(obsini)[4] != nleadt) {
      stop(paste0(
        "The dimensions of parameter obsini and the parameters ",
        "nleadt and nstartd must match:\n  dim(obsini) = c(",
        dim(obsini)[3], ", ", dim(obsini)[4], ")\n  nstartd = ",
        nstartd, "  nleadt = ", nleadt
      ))
    }
  }
  if (!is.null(fxerr)) {
    if (!is.numeric(fxerr)) {
        stop("Parameter fxerr must be numeric.")
    }
  }
  time <- seq(1, nstartd)
  if (nstartd < 0) {
      stop("Number of start dates must be positive")
  }
  if (nleadt < 0) {
      stop("Number of lead-times must be positive")
  }
  if (nmemb < 0) {
      stop("Number of members must be positive")
  }
  if (!is.null(obsini)) {
    obs_ano <- obsini
  }
  else {
    obs_ano <- array(rnorm(nleadt * nstartd, mean = 0, sd = sig),
        dim = c(1, 1, nstartd, nleadt))
    obs_trend <- array(t(time) * rep(trend, times = nstartd),
        , dim = c(1, 1, nstartd, nleadt))
    obs <- obs_ano + obs_trend
    trend <- rep(c(trend), times = nleadt)
    sig <- rep(c(sig), times = nleadt)
  }
  forecast <- array(dim = c(length(gamma), nmemb, nstartd, nleadt))
  for (j in 1:nstartd) {
    for (f in 1:nleadt) {
      for (g in 1:length(gamma)) {
        auto_term <-  obs_ano[1, 1, j, f]
        if (is.numeric(fxerr)) {
          conf_term <- fxerr
        }
        else {
          conf_term <- rnorm(1, mean = 0, sd = beta)
        }
        trend_term <- gamma[g] * trend * j
        var_corr <- rnorm(
          nmemb,
          mean = 0,
          sd = sqrt(sig - alpha ^ 2 - beta ^ 2)
        )
        forecast[g, , j, f] <- matrix(auto_term, c(nmemb,1)) + matrix(conf_term, c(nmemb, 1)) + matrix(trend_term, c(nmemb, 1)) 
      }
    }
  }
  list(mod = forecast, obs = obs_ano)
}

forecast <- ToyModel(
  alpha = a,
  beta = b,
  gamma = g,
  nmemb = nm,
  obsini = InsertDim(data, 1, 1), # nolint
  nstartd = 1,
  nleadt = dim(data)[time_dim]
)
print(  min(c(forecast$obs, forecast$mod)))
print( max(c(forecast$obs, forecast$mod)))
#Quick plot of results
print(brewer.pal(n = nm, name = "Reds"))
jpeg(
  paste0(
    plot_dir, "/", "synthetic_", gsub(".nc", "",
    basename(fullpath_filenames)), ".jpg"
  ),
  height = 460,
  width = 600
)
plot(time, forecast$obs, type = "l",
 # ylim = c(
   # min(c(forecast$obs, forecast$mod), rm.na = TRUE),
  #  max(c(forecast$obs, forecast$mod), rm.na = TRUE)
 # ),
  ylab = paste(var0, "(", units, ")"),
  main = paste(nm, "synthetic members generated"),
  bty = "n"
)
matlines(
  time,
  t(forecast$mod[1, , 1, ]),
  col = brewer.pal(n = nm, name = "Blues")
)
lines(time, forecast$obs, lwd = 2)
dev.off()

#data <- Apply(list(data), target_dims = list(time_dim), AtomicFun = "Obs2ToyForecast")$mod
obs_data <- forecast$obs
data <- forecast$mod[1, , 1, ]
names(dim(data))[c(1, 2)] <- c("number", "time")

attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = "time", long_name = "time", units = "days since 1970-01-01 00:00:00", prec = "double", dim = list(list(name="time", unlim = FALSE))))
attr(time, "variables") <- metadata
metadata <- list(index = list(dim = list(list(name="time", unlim = FALSE, prec = "double"))))
names(metadata)[1] <- var0
attr(data, "variables") <- metadata
variable_list <- list(variable = data, time = time)
names(variable_list)[1] <- var0
print(str(data))
ArrayToNetCDF(variable_list,
              paste0(plot_dir, "/", "synthetic_", basename(fullpath_filenames)))
