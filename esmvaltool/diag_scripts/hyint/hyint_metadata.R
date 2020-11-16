# HyInt metadata

##########################################################
#----------------Metadata functions------------------------#
##########################################################

getmetadata_indices <- function(var, sfile) {
  ncfile <- nc_open(sfile)
  long_name <- (ncatt_get(ncfile, var, "long_name"))$value
  units <- (ncatt_get(ncfile, var, "units"))$value
  missval <- (ncatt_get(ncfile, var, "missing_value"))$value
  if (units == 0) {
    units <- ""
  }
  nc_close(ncfile)
  metadata <-
    list(
      long_name = long_name,
      units = units,
      missval = missval
    )

  return(metadata)
}

setmetadata_indices <- function(var) { # nolint
  longvar <- ""
  unit <- ""

  # name of the var
  if (var == "pry") {
    longvar <- "Annual mean precipitation"
    unit <- "days"
  }
  if (var == "pry_mean") {
    longvar <-
      "Normalization function: Annual mean precipitation averaged over
    available XX century data"
    unit <- "days"
  }
  if (var == "pry_mean_sd") {
    longvar <-
      "Normalization function: Standard deviation of the annual mean
    precipitation averaged over available XX century data"
    unit <- "days"
  }
  if (var == "pry_norm") {
    longvar <- "Normalized annual mean precipitation"
    unit <- ""
  }

  if (var == "dsl") {
    longvar <- "Annual mean dry spell length"
    unit <- "days"
  }
  if (var == "dsl_mean") {
    longvar <-
      "Normalization function: Annual mean dry spell length averaged
    over available XX century data"
    unit <- "days"
  }
  if (var == "dsl_mean_sd") {
    longvar <-
      "Normalization function: Standard deviation of the annual mean
    dry spell length averaged over available XX century data"
    unit <- "days"
  }
  if (var == "dsl_norm") {
    longvar <- "Normalized annual mean dry spell length"
    unit <- ""
  }


  if (var == "dsl_tseries") {
    longvar <- "dsl timeseries over selected regions"
    unit <- ""
  }
  if (var == "dsl_tseries_sd") {
    longvar <- "standard deviation about the mean dsl timeseries"
    unit <- ""
  }
  if (var == "dsl_trend") {
    longvar <- "trend coefficients over selected time period"
    unit <- ""
  }
  if (var == "dsl_trend_stat") {
    longvar <- "statistics of trend over selected time period"
    unit <- ""
  }

  if (var == "wsl") {
    longvar <- "Annual mean wet spell length"
    unit <- "days"
  }
  if (var == "wsl_mean") {
    longvar <-
      "Normalization function: Annual mean wet spell length averaged
    over available XX century data"
    unit <- "days"
  }
  if (var == "wsl_mean_sd") {
    longvar <-
      "Normalization function: Standard deviation of the annual mean
    wet spell length averaged over available XX century data"
    unit <- "days"
  }
  if (var == "wsl_norm") {
    longvar <- "Normalized annual mean wet spell length"
    unit <- ""
  }

  if (var == "int") {
    longvar <- "Annual mean precipitation intensity"
    unit <- "mm day-1"
  }
  if (var == "int_mean") {
    longvar <-
      "Normalization function: Annual mean precipitation intensity
    averaged over available XX century data"
    unit <- "mm day-1"
  }
  if (var == "int_mean_sd") {
    longvar <-
      "Normalization function: Standard deviation of the annual mean
    precipitation intensity averaged over available XX century data"
    unit <- "mm day-1"
  }
  if (var == "int_norm") {
    longvar <- "Normalized annual mean precipitation intensity"
    unit <- ""
  }

  if (var == "pa") {
    longvar <- "Precipitation area: area over which of any given day
    precipitation occurs."
    unit <- "mm day-1 km2"
  }
  if (var == "pa_mean") {
    longvar <- "Normalization function: Precipitation ara averaged over
    available XX century data"
    unit <- "mm day-1"
  }
  if (var == "pa_mean_sd") {
    longvar <- "Normalization function: Standard deviation of the
    precipitation area averaged over available XX century data"
    unit <- "mm day-1"
  }
  if (var == "pa_norm") {
    longvar <- "Normalized precipitation area"
    unit <- ""
  }

  if (var == "hyint") {
    longvar <- "Hydroclimatic intensity index"
    unit <- ""
  }

  metadata <- list(longvar = longvar, unit = unit)

  return(metadata)
}
