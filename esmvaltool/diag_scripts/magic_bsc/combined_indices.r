Sys.setenv(TAR = "/bin/tar") # nolint
library(s2dverification)
library(multiApply) # nolint
library(ggplot2)
library(yaml)
library(ncdf4)
library(ClimProjDiags) #nolint
library(abind)
library(climdex.pcic)

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


# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
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
starting <- c(unlist(unname(start_year)))[1]
end_year <- lapply(input_files_per_var, function(x) x$end_year)
ending <- c(unlist(unname(end_year)))[1]
start_year <- as.POSIXct(as.Date(paste0(starting, "-01-01"), "%Y-%m-%d"))
end_year <- as.POSIXct(as.Date(paste0(ending, "-12-31"), "%Y-%m-%d"))

#Parameters for Season() function
monini <- 1
moninf <- params$moninf
monsup <- params$monsup
months <- ""
region <- params$region
running_mean <- params$running_mean
timestamp <- ""
multi_year_average <- params$multi_year_average
### Settings:
print("Settings")
print(running_mean)
print(multi_year_average)
print(moninf)
print(monsup)
print(start_year)
print(end_year)
print(region)

if (region == "Nino3") {
  lon_min <- 360 - 150
  lon_max <- 360 - 90
  lat_min <- -5
  lat_max <- 5
} else if (region == "Nino3.4") {
  lon_min <- 360 - 170
  lon_max <- 360 - 120
  lat_min <- -5
  lat_max <- 5
} else if (region == "Nino3.4") {
  lon_min <- 360 - 160
  lon_max <- 360 - 150
  lat_min <- -5
  lat_max <- 5
} else if (region == "NAO") {
  lon_min <- 360 + c(-90, -90)
  lon_max <- c(40, 40)
  lat_min <- c(25, 60)
  lat_max <- c(45, 80)
} else if (region == "SOI") {
  lon_min <- c(90, 360 - 130)
  lon_max <- c(140, 360 - 80)
  lat_min <- c(-5, -5)
  lat_max <- c(5, 5)
}
### Load data 
data_nc <- nc_open(fullpath_filenames)
lat <- as.vector(ncvar_get(data_nc, "lat"))
lon <- as.vector(ncvar_get(data_nc, "lon"))
units <- ncatt_get(data_nc, var0, "units")$value
long_names <-  ncatt_get(data_nc, var0, "long_name")$value

data <- InsertDim(ncvar_get(data_nc, var0), 1, 1) # nolint
names(dim(data)) <- c("model", "lon", "lat", "time")
time <- seq(start_year, end_year, "month")
nc_close(data_nc)

jpeg(paste0(plot_dir, "/", "plot0.jpg"))
PlotEquiMap(data[1, , , 1], lon = lon, lat = lat, filled = F ) #nolint
dev.off()


jpeg(paste0(plot_dir, "/", "time.jpg"))
plot(1 : length(time), data[1, 1, 1,], type = "l")

if (!is.null(running_mean)) {
    data <- Smoothing(data, runmeanlen = running_mean, numdimt = 4) #nolint
    timestamp <- paste0(running_mean, "-month-running-mean-")
}
lines(1 : length(time), data[1, 1, 1,], col = "blue")

if (!is.null(moninf)) {
  data <- Season(data, posdim = 4, monini = monini, #nolint
                 moninf = moninf, monsup = monsup)
  months <- paste0(month.abb[moninf], "-", month.abb[monsup])
}
lines(seq(1, length(time), 12), data[1, 1, 1,], col = "red")
dev.off()

jpeg(paste0(plot_dir, "/", "plot1.jpg"))
PlotEquiMap(data[1, , , 2], lon = lon, lat = lat, filled = F) #nolint
dev.off()

if (length(lon_min) == 1) {
#nolint start
# print(str(data))
# jpeg(paste0(plot_dir, "/", "plot2.jpg"))
# PlotEquiMap(data$data[1, , , 2], lon = data$lon, lat = data$lat, filled = F)
# dev.off()
#nolint end
  data <- WeightedMean(data, lon = lon, lat = lat, #nolint
                       region = c(lon_min, lon_max, lat_min, lat_max),
                       londim = 2, latdim = 3, mask = NULL)
print("HERE3")
  data <- drop(data)
print(str(data))
print(dim(data))
} else {
print("NOHERE")
#nolint start
#    data1_aux <- SelBox(data, lon = lon, lat = lat,
#                   region = c(lon_min[1], lon_max[1], lat_min[1], lat_max[1]),
#                   londim = 2, latdim = 3)
#    data2_aux <- SelBox(data, lon = lon, lat = lat,
#                   region = c(lon_min[2], lon_max[2], lat_min[2], lat_max[2]),
#                   londim = 2, latdim = 3)
# print(str(data1))
# jpeg(paste0(plot_dir, "/plot3.jpg"))
# PlotEquiMap(data1$data[1, , , 1], lon = data1$lon, lat = data1$lat, filled = F)
# dev.off()
# print(str(data2))
# jpeg(paste0(plot_dir, "/plot4.jpg"))
# PlotEquiMap(data2$data[1, , , 1], lon = data2$lon, lat = data2$lat, filled = F)
# dev.off()
#nolint end
    data1 <- WeightedMean(data, lon = lon, lat = lat, #nolint
                    region = c(lon_min[1], lon_max[1], lat_min[1], lat_max[1]),
                    londim = 2, latdim = 3, mask = NULL)
    data2 <- WeightedMean(data, lon = lon, lat = lat, #nolint
                    region = c(lon_min[2], lon_max[2], lat_min[2], lat_max[2]),
                    londim = 2, latdim = 3, mask = NULL)
    data1 <- drop(data1)
    data2 <- drop(data2)
print(str(data1))
print(str(data2))
    data <- CombineIndices(list(data1, data2), weights = c(1, -1), #nolint
                           operation = "add")
}

dimtime <- ncdim_def(name = "Time", units = "years",
                     vals = starting : ending,
                     longname = "Time")
defdata <- ncvar_def(name = "data", units = units, dim = list(time = dimtime),
               longname = paste("Index for region", region, "Variable", var0))
filencdf <- paste0(work_dir, "/", var0, "_", timestamp, "_", months, "_",
                   starting, ending, "_", ".nc")
file <- nc_create(filencdf, list(defdata))
ncvar_put(file, defdata, data)
nc_close(file)

print(summary(data))
png(paste0(plot_dir, "/", "Index_", region, ".png"), width = 300, height = 300)
plot(starting : ending, data, type = "l", col = "purple", lwd = 2, bty = "n",
     xlab = "Time (years)", ylab = "Index", main = paste("Region", region))
dev.off()
# Set provenance for output files
xprov <- list(ancestors = list(fullpath_filenames),
              authors = list("hunt_al", "manu_ni"),
              projects = list("c3s-magic"),
              caption = "Combined selection",
              statistics = list("other"),
              moninf = params$moninf,
              monsup = params$monsup,
              region = list(params$region),
              running_mean = params$running_mean,
              realms = list("atmos"),
              themes = list("phys"))
provenance[[filencdf]] <- xprov

# Write provenance to file
write_yaml(provenance, provenance_file)
