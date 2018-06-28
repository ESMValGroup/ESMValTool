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
source('https://earth.bsc.es/gitlab/es/s2dverification/raw/develop-Magic_WP6/R/WeightedMean.R')

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

#fullpath_filenames <-  "/scratch/Earth/ahunter/esmvaltool_input/CMIP5/Tier2/NCEP/OBS_NCEP_reanaly_1_T3M_ta_200001-200212.nc"
data <- Start(model = fullpath_filenames,
                        var = var0,
                        var_var = 'var_names',
                        time = 'all',
                        lat = 'all',
                        lon = 'all',
                        plev = "all",
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


Obs2ToyForecast <- function(x) {
nleadt <- length(x)
x <- array(x, c(1, 1, 1, length(x)))
x <- ToyModel(alpha = a, beta = b, gamma = g, nmemb = nm,
              obsini = x, nstartd = 1, nleadt = nleadt)
x$obs <- adrop(x$obs, c(1,3))
x$mod <- adrop(x$mod, c(1,3))
return(x)
}



data <- Apply(list(data), target_dims = list(time_dim), AtomicFun = "Obs2ToyForecast")$mod
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
ArrayToNetCDF(variable_list,
              paste0(plot_dir, "/", "synthetic_", basename(fullpath_filenames)))





