####REQUIRED SYSTEM LIBS
####ŀibssl-dev
####libnecdf-dev
####cdo

# conda install -c conda-forge r-ncdf4
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
library(yaml)
library(abind)
library(magic.bsc, lib.loc = '/home/Earth/nperez/git/magic.bsc.Rcheck/')
#source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')

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
a <- params$alpha
b <- params$beta
g <- params$gamma
nm <- params$number_of_members
nstartd <- 1
nleadt = params$no_of_lead_times
print(args)
print(var0)
print(params)
print(fullpath_filenames)
#fullpath_filenames <-  "/scratch/Earth/ahunter/esmvaltool_input/CMIP5/Tier2/NCEP/OBS_NCEP_reanaly_1_T3M_ta_200001-200212.nc"
data <- Start(model = fullpath_filenames,
                        var = var0,
                        var_var = 'var_names',
                        time = 'all',
                        lat = 'all',
                        lon = 'all',
#                        plev = "all",
                        lon_var = 'lon',
                        return_vars = list(time = 'model', lon = 'model', lat = 'model'),
                        retrieve = TRUE)


dim_names <- names(dim(data))
time <- attributes(data)$Variables$dat1$time
lon_dim <- which(names(dim(data)) == "lon")
lat_dim <- which(names(dim(data)) == "lat")
lat <- attr(data, "Variables")$dat1$lat
lon <- attr(data, "Variables")$dat1$lon
attributes(lon) <- NULL
attributes(lat) <- NULL
data <- WeightedMean(data, lat = lat, lon = lon, londim = lon_dim, latdim = lat_dim)
names(dim(data)) <- dim_names[-c(lon_dim, lat_dim)]
time_dim <- which(names(dim(data)) == "time")

ToyModel <- function (alpha = 0.1, beta = 0.4, gamma = 1, sig = 1, trend = 0, 
    nstartd = 30, nleadt = 4, nmemb = 10, obsini = NULL, fxerr = NULL) 
{
    if (any(!is.numeric(c(alpha, beta, gamma, sig, trend, nstartd, 
        nleadt, nmemb)))) {
        stop(paste("Parameters 'alpha', 'beta', 'gamma', 'sig', 'trend', 'nstartd',", 
            "nleadt and nmemb must be numeric."))
    }
    nstartd <- round(nstartd)
    nleadt <- round(nleadt)
    nmemb <- round(nmemb)
    if (!is.null(obsini)) {
        if (!is.numeric(obsini) || !is.array(obsini)) {
            stop("Parameter 'obsini' must be a numeric array.")
        }
        if (length(dim(obsini)) != 4) {
            stop("Parameter 'obsini' must be an array with dimensions c(1, 1, nleadt, nstartd).")
        }
        if (dim(obsini)[3] != nstartd || dim(obsini)[4] != nleadt) {
            stop(paste0("The dimensions of parameter 'obsini' and the parameters 'nleadt' and 'nstartd' must match:\n  dim(obsini) = c(", 
                dim(obsini)[3], ", ", dim(obsini)[4], ")\n  nstartd = ", 
                nstartd, "  nleadt = ", nleadt))
        }
    }
    if (!is.null(fxerr)) {
        if (!is.numeric(fxerr)) {
            stop("Parameter 'fxerr' must be numeric.")
        }
    }
    time <- seq(1, nstartd)
    if (!(sig^2 - alpha^2 - beta^2 > 0)) {
        stop("Model variability not constrained: respect condition \"sig^2-alpha^2-beta^2 > 0\")")
    }
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
    forecast <- array(dim = c(length(gamma), nmemb, nstartd, 
        nleadt))
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
                var_corr <- rnorm(nmemb, mean = 0, sd = sqrt(sig - 
                  alpha^2 - beta^2))
                forecast[g, , j, f] <- matrix(auto_term, c(nmemb, 
                  1)) + matrix(conf_term, c(nmemb, 1)) + matrix(trend_term, 
                  c(nmemb, 1)) + var_corr
            }
        }
    }
    list(mod = forecast, obs = obs_ano)
}
#Obs2ToyForecast <- function(x) {
#nleadt <- length(x)
#x <- array(x, c(1, 1, 1, length(x)))
forecast <- ToyModel(alpha = a, beta = b, gamma = g, nmemb = nm,
              obsini = InsertDim(data,3,1), nstartd = 1, nleadt = dim(data)[time_dim])
#x$obs <- adrop(x$obs, c(1,3))
#x$mod <- adrop(x$mod, c(1,3))
#return(x)
#}
#Quick plot of results
jpeg(paste0(plot_dir, "/syntheticforecast.jpg"))
 plot(time,forecast$obs, type = "l")

for (i in 1:nm){
 lines(time,forecast$mod[1,i,,], col = 2)
}
dev.off()

#data <- Apply(list(data), target_dims = list(time_dim), AtomicFun = "Obs2ToyForecast")$mod
obs_data <- forecast$obs
data <- forecast$mod[1,,1,]
names(dim(data))[c(1, 2)] <- c("number", "time")

attributes(time) <- NULL
dim(time) <- c(time = length(time))
metadata <- list(time = list(standard_name = 'time', long_name = 'time', units = 'days since 1970-01-01 00:00:00', prec = 'double', dim = list(list(name='time', unlim = FALSE))))
attr(time, 'variables') <- metadata
metadata <- list(index = list(dim = list(list(name='time', unlim = FALSE, prec = 'double'))))
names(metadata)[1] <- var0
attr(data, 'variables') <- metadata
variable_list <- list(variable = data, time = time)
names(variable_list)[1] <- var0
print(str(data))
ArrayToNetCDF(variable_list,
              paste0(plot_dir, "/", "synthetic_", basename(fullpath_filenames)))





