library(yaml)
library(s2dverification)
library(multiApply) # nolint
library(ncdf4)
library(climdex.pcic)
library(parallel)
library(ClimProjDiags) # nolint
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir

dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

input_files_per_var <- yaml::read_yaml(params$input_files)
var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x)
  x$dataset)
model_names <- unname(model_names)
var0 <- lapply(input_files_per_var, function(x)
  x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]

experiment <- lapply(input_files_per_var, function(x)
  x$exp)
experiment <- unlist(unname(experiment))

reference_files <- which(unname(experiment) == "historical")
projection_files <- which(unname(experiment) != "historical")

rcp_scenario <- unique(experiment[projection_files])
model_names <- lapply(input_files_per_var, function(x)
  x$dataset)
model_names <- unlist(unname(model_names))[projection_files]

start_reference <-
  lapply(input_files_per_var, function(x)
    x$start_year)
start_reference <-
  c(unlist(unname(start_reference))[reference_files])[1]
end_reference <- lapply(input_files_per_var, function(x)
  x$end_year)
end_reference <-
  c(unlist(unname(end_reference))[reference_files])[1]

start_projection <-
  lapply(input_files_per_var, function(x)
    x$start_year)
start_projection <-
  c(unlist(unname(start_projection))[projection_files])[1]
end_projection <-
  lapply(input_files_per_var, function(x)
    x$end_year)
end_projection <-
  c(unlist(unname(end_projection))[projection_files])[1]


op <- as.character(params$operator)
qtile <- params$quantile
spell_length <- params$min_duration
season <- params$season

reference_filenames <- fullpath_filenames[reference_files]
projection <- "NULL"
reference_filenames <- fullpath_filenames[reference_files]
hist_nc <- nc_open(reference_filenames)
var0 <- unlist(var0)
historical_data <- ncvar_get(hist_nc, var0)

names(dim(historical_data)) <- rev(names(hist_nc$dim))[-1]
lat <- ncvar_get(hist_nc, "lat")
lon <- ncvar_get(hist_nc, "lon")
units <- ncatt_get(hist_nc, var0, "units")$value
calendar <- ncatt_get(hist_nc, "time", "calendar")$value
long_names <- ncatt_get(hist_nc, var0, "long_name")$value
time <- ncvar_get(hist_nc, "time")
start_date <- as.POSIXct(substr(ncatt_get(
  hist_nc, "time",
  "units"
)$value, 11, 29))
nc_close(hist_nc)
time <- as.Date(time, origin = start_date, calendar = calendar)


historical_data <- as.vector(historical_data)
dim(historical_data) <- c(
  model = 1,
  var = 1,
  lon = length(lon),
  lat = length(lat),
  time = length(time)
)
historical_data <- aperm(historical_data, c(1, 2, 5, 4, 3))
attr(historical_data, "Variables")$dat1$time <- time
names(dim(historical_data)) <-
  c("model", "var", "time", "lon", "lat")
time_dimension <- which(names(dim(historical_data)) == "time")

base_range <- c(
  as.numeric(substr(start_reference, 1, 4)),
  as.numeric(substr(end_reference, 1, 4))
)
threshold <-
  Threshold(
    historical_data,
    base.range = base_range,
    # nolint
    calendar = calendar,
    qtiles = qtile,
    ncores = NULL,
    na.rm = TRUE
  )

projection_filenames <- fullpath_filenames[projection_files]
for (i in 1:length(projection_filenames)) {
  proj_nc <- nc_open(projection_filenames[i])
  projection_data <- ncvar_get(proj_nc, var0)
  time <- ncvar_get(proj_nc, "time")
  start_date <- as.POSIXct(substr(ncatt_get(
    proj_nc, "time",
    "units"
  )$value, 11, 29))
  calendar <- ncatt_get(hist_nc, "time", "calendar")$value
  time <- as.Date(time, origin = start_date, calendar = calendar)
  nc_close(proj_nc)
  projection_data <- as.vector(projection_data)
  dim(projection_data) <- c(
    model = 1,
    var = 1,
    lon = length(lon),
    lat = length(lat),
    time = length(time)
  )
  projection_data <- aperm(projection_data, c(1, 2, 5, 4, 3))
  attr(projection_data, "Variables")$dat1$time <- time
  names(dim(projection_data)) <-
    c("model", "var", "time", "lon", "lat")
  # ------------------------------
  heatwave <- WaveDuration(
    projection_data,
    threshold,
    # nolint
    calendar = calendar,
    op = op,
    spell.length = spell_length,
    by.seasons = TRUE,
    ncores = NULL
  )
  if (season == "summer") {
    heatwave_season <-
      heatwave$result[seq(2, dim(heatwave$result)[1] - 2,
        by = 4
      ), 1, 1, , ] # nolint
    years <-
      heatwave$years[seq(2, length(heatwave$years) - 2, by = 4)]
  } else if (season == "winter") {
    heatwave_season <-
      heatwave$result[seq(1, dim(heatwave$result)[1] - 2,
        by = 4
      ), 1, 1, , ] # nolint
    years <-
      heatwave$years[seq(1, length(heatwave$years) - 1, by = 4)]
  } else if (season == "spring") {
    heatwave_season <-
      heatwave$result[seq(3, dim(heatwave$result)[1] - 2,
        by = 4
      ), 1, 1, , ] # nolint
    years <-
      heatwave$years[seq(3, length(heatwave$years) - 2, by = 4)]
  } else {
    heatwave_season <-
      heatwave$result[seq(4, dim(heatwave$result)[1] - 2,
        by = 4
      ), 1, 1, , ] # nolint
    years <-
      heatwave$years[seq(4, length(heatwave$years) - 2, by = 4)]
  }

  data <- heatwave_season
  names(dim(data)) <- c("time", "lon", "lat")
  attributes(lon) <- NULL
  attributes(lat) <- NULL
  dim(lon) <- c(lon = length(lon))
  dim(lat) <- c(lat = length(lat))
  time <- as.numeric(substr(years, 1, 4))
  attributes(time) <- NULL
  dim(time) <- c(time = length(time))
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
  dimtime <-
    ncdim_def(
      name = "time",
      units = "years since 0-0-0 00:00:00",
      vals = time,
      longname = "time"
    )
  defdata <- ncvar_def(
    name = "duration",
    units = "days",
    dim = list(
      season = dimtime,
      lat = dimlat,
      lon = dimlon
    ),
    longname = paste(
      "Number of days during the period",
      start_projection,
      "-",
      end_projection,
      "for",
      season,
      "in which",
      var0,
      "is",
      op,
      "than the",
      qtile,
      "quantile obtained from",
      start_reference,
      "-",
      end_reference
    )
  )
  filencdf <-
    paste0(
      work_dir,
      "/",
      var0,
      "_extreme_spell_duration",
      season,
      "_",
      model_names,
      "_",
      rcp_scenario[i],
      "_",
      start_projection,
      "_",
      end_projection,
      ".nc"
    )
  file <- nc_create(filencdf, list(defdata))
  ncvar_put(file, defdata, data)
  nc_close(file)

  # Check dimension order:
  if (length(lat) != dim(data)["lat"] |
    length(lon) != dim(data)["lon"]) {
    if (length(lat) == dim(data)["lon"] &
      length(lon) == dim(data)["lat"]) {
      poslat <- which(names(dim(data)) == "lat")
      poslon <- which(names(dim(data)) == "lon")
      names(dim(data))[poslat] <- "lon"
      names(dim(data))[poslon] <- "lat"
    }
  }
  timeseries <- WeightedMean(data,
    lon = as.vector(lon),
    # nolint
    lat = as.vector(lat),
    mask = NULL
  )
  data_frame <- data.frame(Experiment = timeseries)
  years <- rep(start_projection:end_projection)
  data_frame$Year <- c(years)

  title <-
    paste(
      "Days ",
      season,
      var0,
      paste0(
        substr(start_projection, 1, 4),
        "-",
        substr(end_projection, 1, 4)
      ),
      op,
      "the",
      qtile * 100,
      "th quantile for",
      paste0(
        substr(start_reference, 1, 4),
        "-",
        substr(end_reference, 1, 4),
        " (",
        rcp_scenario[i],
        ")"
      )
    )
  filepng1 <- file.path(
    plot_dir,
    paste0(
      "Time_",
      var0,
      "_extreme_spell_duration",
      season,
      "_",
      model_names,
      "_",
      rcp_scenario[i],
      "_",
      start_projection,
      "_",
      end_projection,
      ".png"
    )
  )

  g <- ggplot(data_frame, aes(x = Year, y = Experiment)) +
    theme_bw() + geom_line() + ylab(paste0("Number of Days")) +
    xlab("Year") + theme(
      text = element_text(size = 12),
      legend.text = element_text(size = 12),
      axis.title = element_text(size = 12)
    ) +
    ggtitle(title)
  ggsave(
    filename = filepng1,
    g,
    device = NULL,
    dpi = "print",
    width = 10,
    height = 6
  )

  brks <- seq(0, max(Mean1Dim(data, 1), na.rm = TRUE), 4)
  filepng <-
    paste0(
      plot_dir,
      "/",
      var0,
      "_extreme_spell_duration",
      season,
      "_",
      model_names,
      "_",
      rcp_scenario[i],
      "_",
      start_projection,
      "_",
      end_projection,
      ".png"
    )
  PlotEquiMap(
    Mean1Dim(data, 1),
    lat = lat,
    lon = lon,
    # nolint
    filled.continents = FALSE,
    brks = brks,
    color_fun = clim.palette("yellowred"),
    units = "Days",
    toptitle = title,
    fileout = filepng,
    title_scale = 0.5
  )

  # Set provenance for output files

  xprov <-
    list(
      ancestors = list(projection_filenames, reference_filenames),
      authors = list(
        "hunter_alasdair",
        "manubens_nicolau",
        "caron_louis-philippe"
      ),
      projects = list("c3s-magic"),
      caption = title,
      statistics = list("other"),
      op = as.character(params$operator),
      qtile = params$quantile,
      spell_length = params$min_duration,
      season = params$season,
      realms = list("atmos"),
      themes = list("phys"),
      plot_file = filepng
    )

  provenance[[filencdf]] <- xprov
}

# Write provenance to file
write_yaml(provenance, provenance_file)
