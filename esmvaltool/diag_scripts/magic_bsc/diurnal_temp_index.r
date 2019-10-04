library(yaml)
library(s2dverification)
library(multiApply) # nolint
library(climdex.pcic)
library(ClimProjDiags) # nolint
library(parallel)
library(ncdf4)

## Insurance products
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

# FOR THE FIRST METADATA.yml
input_files_tasmax <- yaml::read_yaml(params$input_files[1])
model_names <- input_files_tasmax[[1]]$dataset
var_names_tmax <- input_files_tasmax[[1]]$short_name
experiment <- lapply(input_files_tasmax, function(x) {
  x$exp
}) # nolint
filename_tasmax <-
  lapply(input_files_tasmax, function(x) {
    x$filename
  }) # nolint

input_files_tasmin <- yaml::read_yaml(params$input_files[2])
var_names_tmin <- input_files_tasmin[[1]]$short_name
filename_tasmin <-
  lapply(input_files_tasmin, function(x) {
    x$filename
  }) # nolint

reference_files <- which(experiment == "historical")
projection_files <- which(experiment != "historical")

start_historical <- input_files_tasmax[[reference_files]]$start_year
end_historical <- input_files_tasmax[[reference_files]]$end_year
start_projection <-
  input_files_tasmax[[projection_files[1]]]$start_year
end_projection <- input_files_tasmax[[projection_files[1]]]$end_year


fullpath_hist_tasmax <- filename_tasmax[[reference_files]]
file <- nc_open(fullpath_hist_tasmax)
historical_tasmax <- ncvar_get(file, "tasmax")
names(dim(historical_tasmax)) <- rev(names(file$dim))[-1]
lat <- ncvar_get(file, "lat")
lon <- ncvar_get(file, "lon")
units <- ncatt_get(file, "tasmax", "units")$value
calendario <- ncatt_get(file, "time", "calendar")$value
long_names <- ncatt_get(file, "tasmax", "long_name")$value
time <- ncvar_get(file, "time")
start_date <- as.POSIXct(substr(ncatt_get(
  file, "time",
  "units"
)$value, 11, 29))
nc_close(file)

fullpath_hist_tasmin <- filename_tasmin[[reference_files]]
file <- nc_open(fullpath_hist_tasmin)
historical_tasmin <- ncvar_get(file, "tasmin")
names(dim(historical_tasmin)) <- rev(names(file$dim))[-1]
lat <- ncvar_get(file, "lat")
lon <- ncvar_get(file, "lon")
units <- ncatt_get(file, "tasmin", "units")$value
calendario <- ncatt_get(file, "time", "calendar")$value
long_names <- ncatt_get(file, "tasmin", "long_name")$value
tunits <- ncatt_get(file, "time", "units")$value
time <- ncvar_get(file, "time")
start_date <- as.POSIXct(substr(ncatt_get(
  file, "time",
  "units"
)$value, 11, 29))
nc_close(file)
dia <- as.Date(strsplit(tunits, " ")[[1]][3], format = "%Y-%m-%d")
time <- time + dia


dtr_base <- DTRRef(
  tmax = historical_tasmax,
  # nolint
  tmin = historical_tasmin,
  by.seasons = TRUE,
  ncores = NULL,
  dates = time,
  calendar = calendario
)

for (i in 1:length(projection_files)) {
  fullpath_projection_tasmax <- filename_tasmax[[projection_files[i]]]
  file <- nc_open(fullpath_projection_tasmax)
  rcp_tasmax <- ncvar_get(file, "tasmax")
  names(dim(rcp_tasmax)) <- rev(names(file$dim))[-1]
  lat <- ncvar_get(file, "lat")
  lon <- ncvar_get(file, "lon")
  units <- ncatt_get(file, "tasmax", "units")$value
  calendario <- ncatt_get(file, "time", "calendar")$value
  long_names <- ncatt_get(file, "tasmax", "long_name")$value
  time <- ncvar_get(file, "time")
  start_date <- as.POSIXct(substr(ncatt_get(
    file, "time",
    "units"
  )$value, 11, 29))
  nc_close(file)

  fullpath_projection_tasmin <-
    filename_tasmin[[projection_files[i]]]
  file <- nc_open(fullpath_projection_tasmin)
  rcp_tasmin <- ncvar_get(file, "tasmin")
  names(dim(rcp_tasmin)) <- rev(names(file$dim))[-1]
  lat <- ncvar_get(file, "lat")
  lon <- ncvar_get(file, "lon")
  units <- ncatt_get(file, "tasmin", "units")$value
  calendario <- ncatt_get(file, "time", "calendar")$value
  long_names <- ncatt_get(file, "tasmin", "long_name")$value
  tunits <- ncatt_get(file, "time", "units")$value
  time <- ncvar_get(file, "time")
  start_date <- as.POSIXct(substr(ncatt_get(
    file, "time",
    "units"
  )$value, 11, 29))
  nc_close(file)

  dia <-
    as.Date(strsplit(tunits, " ")[[1]][3], format = "%Y-%m-%d")
  time <- time + dia

  dtr_indicator <-
    DTRIndicator(
      rcp_tasmax,
      rcp_tasmin,
      ref = dtr_base,
      by.seasons = TRUE,
      ncores = NULL,
      dates = time,
      calendar = calendario
    )

  dtr_rcp <- array(dim = c(4, length(lon), length(lat)))
  for (j in 1:4) {
    dtr_rcp[j, , ] <-
      Mean1Dim(dtr_indicator$indicator[, j, , ], 1)
  }
  names(dim(dtr_rcp)) <- c("season", "lon", "lat")
  title <- paste0(
    "Number of days exceeding the DTR by 5 degrees\n",
    "during the period ",
    start_projection,
    "-",
    end_projection
  )
  PlotLayout(
    PlotEquiMap,
    plot_dims = c("lon", "lat"),
    # nolint
    var = dtr_rcp,
    colNA = "white",
    lon = lon,
    lat = lat,
    titles = c("DJF", "MAM", "JJA", "SON"),
    toptitle = title,
    filled.continents = FALSE,
    units = "Days",
    axelab = FALSE,
    draw_separators = TRUE,
    subsampleg = 1,
    brks = seq(0, max(dtr_rcp, na.rm = TRUE), 2),
    color_fun = clim.palette("yellowred"),
    extra_margin = c(0, 0, 1, 0),
    bar_extra_labels = c(2, 0, 0, 0),
    title_scale = 0.7,
    fileout = file.path(
      plot_dir,
      paste0(
        "Seasonal_DTRindicator_",
        model_names,
        "_",
        start_projection,
        "_",
        end_projection,
        "_",
        start_historical,
        "_",
        end_historical,
        ".png"
      )
    ),
    col_inf = "white",
    col_sup = "darkred"
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
  dimseason <-
    ncdim_def(
      name = "season",
      units = "season",
      vals = 1:4,
      longname = "season of the year: DJF, MAM, JJA, SON"
    )
  defdata <-
    ncvar_def(
      name = "VulnerabilityIndex",
      units = "number_of_days",
      dim = list(
        season = dimseason,
        lat = dimlat,
        lon = dimlon
      ),
      longname = paste0(
        "Number of days exceeding in 5 degrees ",
        "the Diurnal Temeprature Range for ",
        "the reference period"
      )
    )
  filencdf <-
    paste0(
      work_dir,
      "/",
      "Seasonal_DTRindicator_",
      model_names,
      "_",
      start_projection,
      "_",
      end_projection,
      "_",
      start_historical,
      "_",
      end_historical,
      ".nc"
    )
  file <- nc_create(filencdf, list(defdata))
  ncvar_put(file, defdata, dtr_rcp)
  nc_close(file)


  # Set provenance for output files
  xprov <-
    list(
      ancestors = list(
        filename_tasmin[[reference_files]],
        filename_tasmax[[reference_files]],
        filename_tasmin[[projection_files[i]]],
        filename_tasmax[[projection_files[i]]]
      ),
      authors = list(
        "hunter_alasdair",
        "manubens_nicolau",
        "caron_louis-philippe"
      ),
      projects = list("c3s-magic"),
      caption = title,
      statistics = list("other"),
      realms = list("atmos"),
      themes = list("phys"),
      plot_file = file.path(
        plot_dir,
        paste0(
          "Seasonal_DTRindicator_",
          model_names,
          "_",
          start_projection,
          "_",
          end_projection,
          "_",
          start_historical,
          "_",
          end_historical,
          ".png"
        )
      )
    )

  provenance[[filencdf]] <- xprov
}

# Write provenance to file
write_yaml(provenance, provenance_file)
