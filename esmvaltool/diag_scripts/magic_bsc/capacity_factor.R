library(abind)
library(climdex.pcic)
library(ggplot2)
library(multiApply) # nolint
library(ncdf4)
library(RColorBrewer) # nolint
library(s2dverification)
library(yaml)

# Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
initial.options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(
  file_arg_name, "",
  initial.options[grep(file_arg_name, initial.options)]
)
script_dirname <- dirname(script_name)

source(file.path(script_dirname, "PC.R"))
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

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


no_of_years <- length(seq(start_year, end_year, 1))
var0 <- unlist(var0)
for (i in seq(1, length(model_names), 1)) {
  data_nc <- nc_open(fullpath_filenames[i])
  data <- ncvar_get(data_nc, var0)

  names(dim(data)) <- c("lon", "lat", "time")
  lat <- ncvar_get(data_nc, "lat")
  lon <- ncvar_get(data_nc, "lon")
  units <- ncatt_get(data_nc, var0, "units")$value
  calendar <- ncatt_get(data_nc, "time", "calendar")$value
  long_names <- ncatt_get(data_nc, var0, "long_name")$value
  time <- ncvar_get(data_nc, "time")
  start_date <- as.POSIXct(substr(ncatt_get(
    data_nc, "time",
    "units"
  )$value, 11, 29))
  nc_close(data_nc)
  time <- as.Date(time, origin = substr(start_date, 1, 10),
                  calendar = calendar)
  time <- as.POSIXct(time, format = "%Y-%m-%d")
  time_dim <- which(names(dim(data)) == "time")
  time <- as.PCICt(time, cal = calendar)
  if (calendar != "360_day" & calendar != "365_day") {
    time <- as.character(time)
    jdays <- as.numeric(strftime(time, format = "%j"))
    pos <- which(substr(time, 6, 10) == "02-29")
    if (length(pos) > 0) {
      time <- time[-pos]
      data <- apply(
        data, c(seq(1, length(dim(data)), 1))[-time_dim],
        function(x) {
          x[-pos]
        }
      )
      data <- aperm(data, c(2, 3, 1))
      names(dim(data)) <- c("lon", "lat", "time")
    }
  }
  dims <- dim(data)
  dims <- append(dims[-time_dim], c(no_of_years, dims[time_dim] /
    no_of_years), after = 2)
  dim(data) <- dims
  # Convert to 100 m wind:
  data <- data * 1.39
  data <- aperm(data, c(3, 4, 2, 1))
  names(dim(data)) <- c("year", "day", "lat", "lon")
  #####################################
  # Cross with PC
  ####################################

  #---------------------------
  # Load PC to use and compute CF for 6h values
  #---------------------------
  seas_data <- Mean1Dim(data, 2)
  pc1 <- read_pc(file.path(script_dirname, power_curves[1]))
  pc2 <- read_pc(file.path(script_dirname, power_curves[2]))
  pc3 <- read_pc(file.path(script_dirname, power_curves[3]))
  pc4 <- read_pc(file.path(script_dirname, power_curves[4]))
  pc5 <- read_pc(file.path(script_dirname, power_curves[5]))


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

  ##############################
  # Make some plots
  ##############################
  #---------------------------
  # Prepare data, labels and colorscales
  #---------------------------
  p <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  q <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  years <- seq(start_year, end_year)
  turb_types <- c("IEC I", "IEC I/II", "IEC II", "IEC II/III", "IEC III")

  seas_data_cf_all <- abind(seas_data_cf1, seas_data_cf2, seas_data_cf3,
    seas_data_cf4, seas_data_cf5,
    along = 0
  )
  mean_data_cf_all <- Mean1Dim(seas_data_cf_all, 2)
  anom_data_cf_all <- seas_data_cf_all - InsertDim( # nolint
    Mean1Dim(seas_data_cf_all, 2), 2, dim(data)[1]
  ) # nolint
  pct_anom_data_cf_all <- (seas_data_cf_all / InsertDim( # nolint
    Mean1Dim(seas_data_cf_all, 2), 2, dim(data)[1]
  )) - 1 # nolint
  #---------------------------
  # Plot seasonal CF maps
  #---------------------------
  filepng <- paste0(
    plot_dir, "/", "capacity_factor_", model_names[i], "_",
    start_year, "-", end_year, ".png"
  )
  title <- paste0(
    seasons, " CF from ", model_names[i],
    " (", start_year, "-", end_year, ")"
  )

  PW_names <- c(
    "Enercon E70", "Gamesa G80", "Gamesa G87",
    "Vestas V100", "Vestas V110"
  )
  PlotLayout(PlotEquiMap, # nolint
    c(3, 2), Mean1Dim(seas_data_cf_all, 2), lon, lat,
    colNA = "white",
    brks = seq(
      from = 0, to = max(seas_data_cf_all, na.rm = TRUE),
      length.out = 10
    ), color_fun = clim.palette("yellowred"),
    filled.continents = FALSE, toptitle = title,
    titles = PW_names, fileout = filepng
  )

  filencdf <- paste0(
    work_dir, "/", "capacity_factor_", model_names[i], "_",
    start_year, "-", end_year, ".nc"
  )
  dimlon <- ncdim_def(
    name = "lon", units = "degrees_east",
    vals = as.vector(lon), longname = "longitude"
  )
  dimlat <- ncdim_def(
    name = "lat", units = "degrees_north",
    vals = as.vector(lat), longname = "latitude"
  )
  dimtime <- ncdim_def(
    name = "season", units = "season",
    vals = start_year:end_year,
    longname = "season of the year: DJF, MAM, JJA, SON"
  )
  dimcurve <- ncdim_def(
    name = "curve", units = "name", vals = seq(1, 5, 1),
    longname = "Power curves of considered turbines"
  )
  names(dim(seas_data_cf_all)) <- c("curve", "time", "lat", "lon")
  defdata <- ncvar_def(
    name = "CapacityFactor", units = "%",
    dim = list(
      season = dimcurve, dimtime, lat = dimlat,
      lon = dimlon
    ),
    longname = paste(
      "Capacity Factor of wind on",
      "different turbines"
    )
  )
  file <- nc_create(filencdf, list(defdata))
  ncvar_put(file, defdata, seas_data_cf_all)
  nc_close(file)

  # Set provenance for output files
  xprov <- list(
    ancestors = list(
      fullpath_filenames[i],
      file.path(script_dirname, power_curves[1]),
      file.path(script_dirname, power_curves[2]),
      file.path(script_dirname, power_curves[3]),
      file.path(script_dirname, power_curves[4]),
      file.path(script_dirname, power_curves[5])
    ),
    authors = list(
      "hunter_alasdair", "perez-zanon_nuria",
      "manubens_nicolau", "lledo_llorenc",
      "caron_louis-philippe", "bojovic_dragana",
      "gonzalez-reviriego_nube"
    ),
    projects = list("c3s-magic"),
    caption = title,
    statistics = list("other"),
    realms = list("atmos"),
    themes = list("phys"),
    plot_file = filepng
  )

  provenance[[filencdf]] <- xprov
}

# Write provenance to file
write_yaml(provenance, provenance_file)
