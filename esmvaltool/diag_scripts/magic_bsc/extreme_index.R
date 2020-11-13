library(s2dverification)
library(multiApply) # nolint
library(yaml)
library(ncdf4)
library(abind)
library(parallel)
library(ClimProjDiags) # nolint

# function to flatten nested lists
flatten_lists <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten_lists)), recursive = FALSE))
}

args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir

dir.create(plot_dir, recursive = TRUE)
dir.create(run_dir, recursive = TRUE)
dir.create(work_dir, recursive = TRUE)

weights <- c(
  t90p = params$weight_t90p,
  t10p = params$weight_t10p,
  Wx = params$weight_Wx,
  rx5day = params$weight_rx5day,
  cdd = params$weight_cdd
)
running_mean <- params$running_mean

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

wdata <- NULL
for (j in 1:4) { # nolint
  input_files_per_var <- yaml::read_yaml(params$input_files[j])
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
  end_reference <-
    lapply(input_files_per_var, function(x)
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

  var0 <- unlist(var0)
  projection <- "NULL"
  reference_filenames <- fullpath_filenames[reference_files]
  hist_nc <- nc_open(reference_filenames)
  historical_data <- ncvar_get(hist_nc, var0)

  names(dim(historical_data)) <- rev(names(hist_nc$dim))[-1]
  lat <- ncvar_get(hist_nc, "lat")
  lon <- ncvar_get(hist_nc, "lon")
  units <- ncatt_get(hist_nc, var0, "units")$value
  calendar <- ncatt_get(hist_nc, "time", "calendar")$value
  long_names <- ncatt_get(hist_nc, var0, "long_name")$value
  time <- ncvar_get(hist_nc, "time")
  # Time correction:
  start_date <-
    as.POSIXct(paste0(start_reference, "-01-01"),
      tz = "UTC",
      format = "%Y-%m-%d"
    )
  end_date <-
    as.POSIXct(paste0(end_reference, "-12-31"),
      tz = "UTC",
      format = "%Y-%m-%d"
    )
  nc_close(hist_nc)
  time <- seq(start_date, end_date, "day")
  if (calendar == "noleap" |
    calendar == "365_day" | calendar == "365") {
    time <- time[format(time, "%m-%d") != "02-29"]
  } else if (calendar == "360_day" | calendar == "360") {
    time <- time[format(time, "%m-%d") != "02-29"]
    time <- time[format(time, "%m-%d") != "01-31"]
    time <- time[format(time, "%m-%d") != "05-31"]
    time <- time[format(time, "%m-%d") != "07-31"]
    time <- time[format(time, "%m-%d") != "10-31"]
    time <- time[format(time, "%m-%d") != "12-31"]
  }
  # nolint start
  # hist_names <- names(dim(historical_data))
  # jpeg(paste0(plot_dir, "/plot1.jpg"))
  # PlotEquiMap(historical_data[1,1,1,,], lon = lon, lat = lat, filled = F)
  # dev.off()
  # ------------------------------
  # Provisional solution to error in dimension order:
  # nolint end
  historical_data <- as.vector(historical_data)
  dim(historical_data) <- c(
    model = 1,
    var = 1,
    lon = length(lon),
    lat = length(lat),
    time = length(time)
  )
  historical_data <- aperm(historical_data, c(1, 2, 5, 3, 4))
  # nolint start
  # ------------------------------
  # jpeg(paste0(plot_dir, "/plot2.jpg"))
  # PlotEquiMap(historical_data[1,1,1,,], lon = lon, lat = lat, filled = F)
  # dev.off()
  # nolint end
  names(dim(historical_data)) <-
    c("model", "var", "time", "lon", "lat")
  time_dimension <- which(names(dim(historical_data)) == "time")
  attributes(lon) <- NULL
  attributes(lat) <- NULL

  dim(lon) <- c(lon = length(lon))
  dim(lat) <- c(lat = length(lat))
  model_dim <- which(names(dim(historical_data)) == "model")
  ### Compute the quantiles and standard deviation for the historical period.
  if (var0 == "tasmin") {
    quantile <- 0.1
    metric <- "t10p"
  } else if (var0 == "tasmax") {
    quantile <- 0.9
    metric <- "t90p"
  } else if (var0 == "sfcWind") {
    historical_data <- 0.5 * 1.23 * (historical_data**3)
    quantile <- 0.9
    metric <- "Wx"
  } else if (var0 == "pr") {
    historical_data <- historical_data * 60 * 60 * 24
    metric <- c("rx5day", "cdd")
  }
  attr(historical_data, "Variables")$dat1$time <- time

  base_sd <- base_sd_historical <- base_mean <- list()
  for (m in seq_along(metric)) {
    if (var0 != "pr") {
      thresholds <- Threshold(
        # nolint
        historical_data,
        calendar = calendar,
        qtiles = quantile,
        ncores = detectCores() - 1,
        na.rm = TRUE
      )
      str(thresholds)
      base_index <- Climdex(
        # nolint
        data = historical_data,
        calendar = calendar,
        metric = metric[m],
        threshold = thresholds,
        ncores = detectCores() - 1
      )
    } else {
      base_index <- Climdex(
        # nolint
        data = historical_data,
        calendar = calendar,
        metric = metric[m],
        ncores = detectCores() - 1
      )
    }
    base_sd[[m]] <- Apply( # nolint
      list(base_index$result),
      target_dims = list(c(1)),
      "sd"
    )$output1
    base_sd_historical[[m]] <- InsertDim( # nolint
      base_sd[[m]], 1, dim(base_index$result)[1]
    )

    if (var0 != "pr") {
      base_mean[[m]] <- 10
      base_mean_historical <- 10
    } else {
      base_mean[[m]] <- Apply( # nolint
        list(base_index$result),
        target_dims = list(c(1)),
        "mean"
      )$output1
      base_mean_historical <- InsertDim( # nolint
        base_mean[[m]], 1, dim(base_index$result)[1]
      )
    }
  }
  # Compute the time series of the relevant index, using the quantiles
  # and standard deviation from the index
  projection_filenames <- fullpath_filenames[projection_files]

  for (i in seq_along(projection_filenames)) {
    proj_nc <- nc_open(projection_filenames[i])
    projection_data <- ncvar_get(proj_nc, var0)
    time <- ncvar_get(proj_nc, "time")
    # Time correction:
    start_date <-
      as.POSIXct(paste0(start_projection, "-01-01"),
        tz = "UTC",
        format = "%Y-%m-%d"
      )
    end_date <-
      as.POSIXct(paste0(end_projection, "-12-31"),
        tz = "UTC",
        format = "%Y-%m-%d"
      )
    nc_close(proj_nc)
    time <- seq(start_date, end_date, "day")
    if (calendar == "noleap" |
      calendar == "365_day" | calendar == "365") {
      time <- time[format(time, "%m-%d") != "02-29"]
    } else if (calendar == "360_day" | calendar == "360") {
      time <- time[format(time, "%m-%d") != "02-29"]
      time <- time[format(time, "%m-%d") != "01-31"]
      time <- time[format(time, "%m-%d") != "05-31"]
      time <- time[format(time, "%m-%d") != "07-31"]
      time <- time[format(time, "%m-%d") != "10-31"]
      time <- time[format(time, "%m-%d") != "12-31"]
    }
    projection_data <- as.vector(projection_data)
    dim(projection_data) <- c(
      model = 1,
      var = 1,
      lon = length(lon),
      lat = length(lat),
      time = length(time)
    )
    projection_data <- aperm(projection_data, c(1, 2, 5, 3, 4))
    attr(projection_data, "Variables")$dat1$time <- time
    names(dim(projection_data)) <-
      c("model", "var", "time", "lon", "lat")
    num_model <- dim(projection_data)["model"]
    print(num_model)
    # nolint start
    # ------------------------------
    # jpeg(paste0(plot_dir, "/plot4.jpg"))
    # PlotEquiMap(projection_data[1,1,1,,], lon = lon, lat = lat, filled = F)
    # dev.off()
    # nolint end
    if (var0 == "pr") {
      projection_data <- projection_data * 60 * 60 * 24
    } else if (var0 == "sfcWind") {
      projection_data <- 0.5 * 1.23 * (projection_data**3)
    }

    for (m in seq_along(metric)) {
      if (var0 != "pr") {
        projection_index <-
          Climdex(
            data = projection_data,
            metric = metric[m],
            calendar = calendar,
            threshold = thresholds,
            ncores = detectCores() - 1
          )
        projection_mean <- 10
      } else {
        # nolint
        projection_index <-
          Climdex(
            data = projection_data,
            metric = metric[m],
            calendar = calendar,
            ncores = detectCores() - 1
          )
        projection_mean <- InsertDim(
          base_mean[[m]], 1, # nolint
          dim(projection_index$result)[1]
        )
      }
      base_sd_proj <- InsertDim(
        base_sd[[m]], 1, # nolint
        dim(projection_index$result)[1]
      )
      projection_index_standardized <-
        (projection_index$result - projection_mean) / base_sd_proj
      for (mod in 1:num_model) {
        model_dim <-
          which(names(dim(projection_index_standardized)) == "model")
        if (length(model_dim) == 0) {
          data <- drop(projection_index_standardized)
        } else {
          print(dim(projection_index_standardized))
          data <- drop(projection_index_standardized[, mod, , , ])
        }
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
          units = "Years",
          vals = start_projection:end_projection,
          longname = "Time in years"
        )
        defdata <- ncvar_def(
          name = "data",
          units = units,
          dim = list(
            year = dimtime,
            lon = dimlon,
            lat = dimlat
          ),
          longname = paste("Annual", metric[m], long_names)
        )
        filencdf <- paste0(
          work_dir,
          "/",
          var0,
          "_",
          metric[m],
          "_risk_insurance_index_",
          model_names,
          "_",
          start_projection,
          "_",
          end_projection,
          "_",
          start_reference,
          "_",
          end_reference,
          ".nc"
        )
        file <- nc_create(filencdf, list(defdata))
        ncvar_put(file, defdata, projection_index_standardized)
        nc_close(file)

        # Plottings
        data <- drop(Mean1Dim(projection_index_standardized, 1))
        title <- paste0(
          "Index for  ",
          metric[m],
          " ",
          substr(start_projection, 1, 4),
          "-",
          substr(end_projection, 1, 4),
          " ",
          " (",
          rcp_scenario[i],
          " ",
          model_names,
          ")"
        )

        breaks <-
          seq(
            -1 * ceiling(max(abs(data), na.rm = TRUE)),
            ceiling(max(abs(data), na.rm = TRUE)),
            2 * ceiling(max(abs(data), na.rm = TRUE)) / 16
          )
        filepng <- paste0(
          plot_dir,
          "/",
          metric[m],
          "_",
          model_names[mod],
          "_",
          rcp_scenario[i],
          "_",
          start_projection,
          "_",
          end_projection,
          ".png"
        )
        PlotEquiMap(
          # nolint
          data,
          lon = lon,
          lat = lat,
          filled.continents = FALSE,
          toptitle = title,
          sizetit = 0.5,
          brks = breaks,
          fileout = filepng,
          colNA = "white"
        )
        # Set provenance for output files
        xprov <-
          list(
            ancestors = flatten_lists(list(projection_filenames, reference_filenames)),
            authors = list(
              "hunter_alasdair",
              "manubens_nicolau",
              "caron_louis-philippe"
            ),
            projects = list("c3s-magic"),
            caption = title,
            statistics = list("other"),
            weight = weights[j + (m - 1)],
            realms = list("atmos"),
            themes = list("phys"),
            plot_file = filepng
          )
        provenance[[filencdf]] <- xprov
        # compute weights in the data
        lon <- as.vector(lon)
        lat <- as.vector(lat)
        temporal <-
          WeightedMean(projection_index_standardized,
            # nolint
            lon = lon,
            lat = lat
          )
        time_dim <- which(names(dim(temporal)) == "year")
        if (!is.null(running_mean)) {
          temporal <- Smoothing(temporal,
            runmeanlen = running_mean,
            # nolint
            numdimt = time_dim
          )
          timestamp <-
            paste0(running_mean, "-month-running-mean-")
        }
        wdata[[j + (m - 1)]] <- temporal
      } # model index
    } # metric index
  } # number of projections
} # variable index
if (!is.numeric(weights)) {
  data <- CombineIndices(wdata, weights = NULL) # nolint
} else {
  data <- CombineIndices(wdata, weights = weights) # nolint
}

# Plotting time series:
data <- drop(data)
if (length(data) >= 5) {
  png(paste0(plot_dir, "/", "CombinedIndices.png"))
  plot(
    start_projection:end_projection,
    data,
    type = "l",
    lwd = 2,
    col = "darkblue",
    xlab = "Time (years)",
    ylab = "Combined indices"
  )
  dev.off()
}
dimtime <- ncdim_def(
  name = "time",
  units = "years",
  vals = start_projection:end_projection,
  longname = "time"
)
defdata <- ncvar_def(
  name = "data",
  units = "adimensional",
  dim = list(time = dimtime),
  longname = paste("Combination", long_names)
)
filencdf <-
  paste0(
    work_dir, "/", "_", paste0(model_names, collapse = "_"),
    ".nc"
  )
file <- nc_create(filencdf, list(defdata))
ncvar_put(file, defdata, data)
nc_close(file)


xprov <- list(
  ancestors = fullpath_filenames,
  authors = list(
    "hunter_alasdair",
    "manubens_nicolau",
    "perez-zanon_nuria"
  ),
  projects = list("c3s-magic"),
  plot_file = paste0(plot_dir, "/", "CombinedIndices.png"),
  caption = "Combined selection",
  statistics = list("other"),
  realms = list("atmos"),
  themes = list("phys")
)
provenance[[filencdf]] <- xprov

# Write provenance to file
write_yaml(provenance, provenance_file)
