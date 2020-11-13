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

input_files_per_var <- yaml::read_yaml(params$input_files)
var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x) x$dataset)
model_names <- unique(unlist(unname(model_names)))
var0 <- lapply(input_files_per_var, function(x) x$short_name)
scenario <- unlist(lapply(input_files_per_var, function(x) x$exp))
ensemble <- unlist(lapply(input_files_per_var, function(x) x$ensemble))
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
standardized <- params$standardized
data_frame <- NULL
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
} else if (region == "Nino4") {
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
for (i in seq(1, length(model_names), 1)) {
    data_nc <- nc_open(fullpath_filenames[i])
    lat <- as.vector(ncvar_get(data_nc, "lat"))
    lon <- as.vector(ncvar_get(data_nc, "lon"))
    units <- ncatt_get(data_nc, var0, "units")$value
    long_names <-  ncatt_get(data_nc, var0, "long_name")$value

    data <- InsertDim(ncvar_get(data_nc, var0), 1, 1) # nolint
    names(dim(data)) <- c("model", "lon", "lat", "time")
    time <- seq(start_year, end_year, "month")
    nc_close(data_nc)
    if (standardized) {
        data <- Apply(list(data), target_dims = c("time"),
                     fun = function(x) {(x - mean(x)) / sqrt(var(x))}) #nolint
        data <- aperm(data$output1, c(2, 3, 4, 1))
        names(dim(data)) <- c("model", "lon", "lat", "time")
    }
    if (!is.null(running_mean)) {
        data <- Smoothing(data, runmeanlen = running_mean, numdimt = 4) #nolint
        timestamp <- paste0(running_mean, "-month-running-mean-")
    }
    if (!is.null(moninf)) {
        data <- Season(data, posdim = 4, monini = monini, #nolint
                       moninf = moninf, monsup = monsup)
        months <- paste0(month.abb[moninf], "-", month.abb[monsup])
    }
    if (length(lon_min) == 1) {
        data <- WeightedMean(data, lon = lon, lat = lat, #nolint
                             region = c(lon_min, lon_max, lat_min, lat_max),
                             londim = 2, latdim = 3, mask = NULL)

        data <- drop(data)
    } else {
        data1 <- WeightedMean(data, lon = lon, lat = lat, #nolint
                    region = c(lon_min[1], lon_max[1], lat_min[1], lat_max[1]),
                    londim = 2, latdim = 3, mask = NULL)
        data2 <- WeightedMean(data, lon = lon, lat = lat, #nolint
                    region = c(lon_min[2], lon_max[2], lat_min[2], lat_max[2]),
                    londim = 2, latdim = 3, mask = NULL)
        data1 <- drop(data1)
        data2 <- drop(data2)
        data <- CombineIndices(list(data1, data2), weights = c(1, -1), #nolint
                               operation = "add")
    }
    if (moninf > monsup) {
        period <- (starting : ending)[-1]
    } else {
        period <- starting : ending
    }
       data_frame <- cbind(data_frame, data)
} # close for loop i
# convert to data frame for ggplot
data_frame_plot <- as.data.frame.table(data_frame)
data_frame_plot$year <- rep(period, length(model_names))
names(data_frame_plot)[2] <- "model"
data_frame_plot$model <- as.factor(sort(rep(seq(1, length(model_names), 1),
                                             length(period))))
for (i in seq(1, length(levels(data_frame_plot$model)), 1)) {
    levels(data_frame_plot$model)[i] <- paste(model_names[i], scenario[i],
                                              ensemble[i])
}
font_size <- 12
g <- ggplot(data_frame_plot, aes(x = year, y = Freq, color = model)) +
     theme_bw() +
     geom_line() + ylab(paste0("Anomaly (", units, ")")) + xlab("Year") +
     theme(text = element_text(size = font_size),
           legend.text = element_text(size = font_size),
           axis.title = element_text(size = font_size)) +
     stat_summary(data =  data_frame_plot, fun.y = "mean",
           mapping = aes(x = data_frame_plot$year, y = data_frame_plot$Freq,
                  group = interaction(data_frame_plot[2, 3]),
                  color = data_frame_plot$model), geom = "line", size = 1) +
      ggtitle(paste0(region, " index for ", var0, " on ", months,
                     " (", starting, "-", ending, ")"))

filepng <-  paste0(plot_dir, "/", region, "_", var0, "_", months,
                   "_running-mean_", running_mean, "_",
                   starting, "-", ending, ".png")
ggsave(filename = filepng, g, device = NULL)
# Save ncdf data_frame
dimtime <- ncdim_def(name = "Time", units = "years",
                     vals = period, longname = "Time")
dimmodel <- ncdim_def(name = "Models", units = "names",
                      vals = seq(1, length(model_names), 1),
                      longname = paste(model_names, scenario))
defdata <- ncvar_def(name = "data", units = units,
                     dim = list(time = dimtime, model = dimmodel),
                     longname = paste("Index for region", region,
                                      "Variable", var0))
filencdf <- paste0(work_dir, "/", var0, "_", timestamp, "_", months, "_",
                   paste(model_names, collapse = ""),
                   starting, ending, "_", ".nc")
file <- nc_create(filencdf, list(defdata))
ncvar_put(file, defdata, data_frame)
nc_close(file)

# Set provenance for output files
xprov <- list(ancestors = fullpath_filenames,
              authors = list("perez-zanon_nuria", "hunter_alasdair",
                             "manubens_nicolau"),
              projects = list("c3s-magic"),
              caption = "Combined selection",
              statistics = list("other"),
              realms = list("atmos"),
              themes = list("phys"))
provenance[[filencdf]] <- xprov

# Write provenance to file
write_yaml(provenance, provenance_file)
