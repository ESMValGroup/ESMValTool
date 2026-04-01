# This file contains utility functions that are available in all
# R diagnostics related to drought indices.

# NOTE: Several codeblocks applying masks to combinations of variables are
# replaced by reusable get_merged_mask function.
#
# NOTE: Variables are converted to Units required by SPEI library when loaded
#   time -> start_year, start_month, end_year, end_month
#   lat -> fails > 90
#   tas/min/max -> kelvin to celsius
#   pr -> to mm month-1
#   sfcWind -> U10 to U2
#   psl -> hPa to kPa
#   rsdt/rsds -> Wm-2 to MJm-2d-1
#
# NOTE: `fillfloat` and `fillvalue` are variables that are used both for the
# same purpose but can have different values? fillvalue is read from reference
# data _FillValue but not used anywhere. Now fillfloat is used to write
# _FillValue attribute of the new created file.

# ---------------------------------------------------------------------------- #
# Variables ------------------------------------------------------------------ #
# ---------------------------------------------------------------------------- #

xprov <- list(
  ancestors = list(""),
  authors = list("berg_peter", "weigel_katja", "lindenlaub_lukas"),
  references = list("vicente10jclim"),
  projects = list("c3s-magic"),
  caption = "",
  statistics = list("other"),
  realms = list("atmos"),
  themes = list("phys"),
  domains = list("global")
)

# ---------------------------------------------------------------------------- #
# funcitons ------------------------------------------------------------------ #
# ---------------------------------------------------------------------------- #

convert_to_monthly <- function(id, v) {
  # converts kg/m2/s to mm month-1 depending on the calendar of the nc files
  # time coordinate. Assuming a density of 1000 kg m-3
  tcal <- ncatt_get(id, "time", attname = "calendar")
  if (tcal$value == "360_day") return(v* 30 * 24 * 3600.)
  time <- ncvar_get(id, "time")
  tunits <- ncatt_get(id, "time", attname = "units")
  tustr <- strsplit(tunits$value, " ")
  stdate <- as.Date(time[1], origin = unlist(tustr)[3])
  nddate <- as.Date(time[length(time)], origin = unlist(tustr)[3])
  if (tcal$value == "365_day") {  # Correct for missing leap years in nddate
    diff <- as.numeric(nddate - stdate, units = "days")
    dcorr <- floor( (diff / 365 - diff / 365.25) * 365.25)
    nddate <- nddate + dcorr
  }
  cnt <- 1
  monarr <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  date <- stdate
  while (date <= nddate){
    year <- as.numeric(substr(date, 1, 4))
    lpyear <- leap_year(year)
    month <- as.numeric(substr(date, 6, 7))
    mdays <- monarr[month]
    pdays <- mdays
    if (month == 2 & lpyear == TRUE){
      pdays <- 29
      mdays <- if (tcal$value != "365_day") 29 else 28
    }
    v[,,cnt] <- v[,,cnt] * mdays * 24 * 3600.
    date <- date + pdays
    cnt <- cnt + 1
  }
  return(v)
}

daily_to_monthly <- function(v, dim=1) {
  # multiplies by number of days (i.e. mm day-1 -> mm/mon)
  # ignores leap years to be compatible to SPEI library
  mlen <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  mlen_rep <- rep(mlen, length.out = dim(v)[dim])
  for (i in 1:dim(v)[dim]){
    if (dim == 1) {
      v[i,,] <- v[i,,] * mlen_rep[i]
    } else if (dim == 2) {
      v[,i,] <- v[,i,] * mlen_rep[i]
    } else if (dim == 3) {
      v[,,i] <- v[,,i] * mlen_rep[i]
    } else {
      stop("time dimension must be 1, 2 or 3")
    }
  }
  return(v)
}

monthly_to_daily <- function(v, dim=1){
  # divide by number of days (i.e. mm/mon -> mm day-1)
  # ignores leap years to be compatible to SPEI library
  mlen <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  mlen_rep <- rep(mlen, length.out = dim(v)[dim])
  for (i in 1:dim(v)[dim]){
    if (dim == 1) {
      v[i,,] <- v[i,,] / mlen_rep[i]
    } else if (dim == 2) {
      v[,i,] <- v[,i,] / mlen_rep[i]
    } else if (dim == 3) {
      v[,,i] <- v[,,i] / mlen_rep[i]
    } else {
      stop("time dimension must be 1, 2 or 3")
    }
  }
  return(v)
}


convert_to_cf <- function(id, v) {
  # NOTE: use mm day-1 converted by preprocessor if possible and convert
  # it to/from month using daily_to_monthly and monthly_to_daily instead.
  #
  # converts mm/mon back to kg/m2/s assuming the input ncfile to have the
  # correct calendar set. (inverse of `convert_to_monthly`)
  tcal <- ncatt_get(id, "time", attname = "calendar")
  if (tcal$value == "360_day") return(v / 30 / 24 / 3600.)
  time <- ncvar_get(id, "time")
  tunits <- ncatt_get(id, "time", attname = "units")
  tustr <- strsplit(tunits$value, " ")
  stdate <- as.Date(time[1], origin = unlist(tustr)[3])
  nddate <- as.Date(time[length(time)], origin = unlist(tustr)[3])
  if (tcal$value == "365_day") {  # Correct for missing leap years in nddate
    diff <- as.numeric(nddate - stdate, units = "days")
    dcorr <- floor( (diff / 365 - diff / 365.25) * 365.25)
    nddate <- nddate + dcorr
  }
  cnt <- 1
  monarr <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  date <- stdate
  while (date <= nddate){
    year <- as.numeric(substr(date, 1, 4))
    lpyear <- leap_year(year)
    month <- as.numeric(substr(date, 6, 7))
    mdays <- monarr[month]
    pdays <- mdays
    if (month == 2 & lpyear == TRUE){
      pdays <- 29
      mdays <- if (tcal$value != "365_day") 29 else 28
    }
    v[,,cnt] <- v[,,cnt] / mdays / 24 / 3600.
    date <- date + pdays
    cnt <- cnt + 1
  }
  return(v)
}

get_dims_from_nc <- function(meta) {
  id <- nc_open(meta$filename)
  xdim <- ncvar_dim(id, "lon")
  if (length(xdim)==0) xdim <- ncvar_dim(id, "longitude")
  ydim <- ncvar_dim(id, "lat")
  if (length(ydim)==0) ydim <- ncvar_dim(id, "latitude")
  tdim <- ncvar_dim(id, "time")
  nc_close(id)
  return(list(xdim, ydim, tdim))
}

get_attrs_from_nc <- function(meta) {
  id <- nc_open(meta$filename)
  allatt <- ncatt_get(id, meta$short_name)
  nc_close(id)
  return(allatt)
}

get_merged_mask <- function(metas, variables="all") {
  first <- TRUE
  for(meta in metas){
    var <- get_var_from_nc(meta)
    msk <- apply(var, c(1,2), FUN = mean, na.rm = TRUE)
    msk[msk > 1e29] <- NA
    if(first){
      mask <- msk
      first <- FALSE
    } else {
      mask[is.na(msk)] <- NA
    }
  }
  mask[mask > 1e29] <- NA  # why again?
  mask[!is.na(mask)] <- 1
  return(mask)
}


get_var_from_nc <- function(meta, custom_var=FALSE) {
  # read data for a specific variable from a nc file
  # NOTE: for some variables the data is converted to match units for SPEI
  if (custom_var != FALSE) var <- custom_var
  else var <- meta$short_name
  id <- nc_open(meta$filename, readunlim=FALSE)
  data <- ncvar_get(id, var)
  dim_names <- attributes(id$dim)$names
  time_dim <- which(sapply(dim_names, function(x) x == "time"))
  # convert to required units
  if (var == "time") {
    tcal <- ncatt_get(id, "time", attname = "calendar")
    tunits <- ncatt_get(id, "time", attname = "units")
    tustr <- strsplit(tunits$value, " ")
    stdate <- as.Date(data[1], origin = unlist(tustr)[3])
    nddate <- as.Date(data[length(data)], origin = unlist(tustr)[3])
    syear <- as.numeric(substr(stdate, 1, 4))
    smonth <- as.numeric(substr(stdate, 6, 7))
    eyear <- as.numeric(substr(nddate, 1, 4))
    emonth <- as.numeric(substr(nddate, 6, 7))
    return(c(syear,smonth,eyear,emonth))
  } else if (var == "lat") {
    if (max(data) > 90){
      print(paste0("Latitude must be [-90,90]: min=",
      min(data), " max=", max(data)))
      stop("Aborting!")
    }
  } else if (var %in% list("tas", "tasmin", "tasmax")) {
    data <- data - 273.15
  } else if (var %in% list("rsdt", "rsds")) {
    data <- data * (86400.0 / 1e6) # W/(m2) to MJ/(m2 d)
  } else if (var %in% list("pr", "evspsbl", "evspsblpot")) {
    if (meta$units == "mm day-1") {
      data <- daily_to_monthly(data, dim=time_dim)
    } else if (meta$units == "mm month-1") {  # do nothing
    } else {
      stop(paste(var,
        " is expected to be in mm day-1 or mm month-1 not in ", meta$units))
    }
  } else if (var == "sfcWind") {
    if (meta$units != "m s-1") {stop("sfcWind is expected to at 10m in m/s")}
    data <- data * (4.87/(log(67.9 * 10.0 - 5.42))) # U10m to U2m (*0.74778)
  } else if (var %in% list("psl", "ps")) {
    data <- data * 0.001 # Pa to kPa
  } else {
    print(paste("No conversion for", var))
  }
  nc_close(id)
  return(data)
}


group_meta <- function(metadata) {
  grouped <- list()
  for (meta in metadata) {
    if (! meta$dataset %in% names(grouped)) {
      grouped[[meta$dataset]] <- list()
    }
    entry <- list(meta)
    names(entry) <- list(meta$filename)
    entry[meta$filename] <- list(meta)
    grouped[[meta$dataset]] <- append(grouped[[meta$dataset]], entry)
  }
  return(grouped)
}


leap_year <- function(year) {
  return(
    ifelse(
      (year %% 4 == 0 & year %% 100 != 0) | year %% 400 == 0,
      TRUE, FALSE
    )
  )
}

list_default <- function(list, key, default) {
  if (!key %in% names(list)) {
    print("key not existing")
    print(key)
    list[[key]] <- default
  }
  return(list)
}

select_var <- function(metas, short_name, strict=TRUE) {
  # NOTE: metas should be a list of metadata entries, but R seems to handle it
  # differently if this list ist of length 1.
  for (meta in metas) {
    if (meta$short_name == short_name) {
      return(meta)
    }
  }
  if (strict) {
    stop(paste("No metadata found for", short_name, "in", meta$dataset))
  } else {
    print(paste("No metadata found for", short_name, "in", meta$dataset))
    return(NULL)
  }
}


t_or_null <- function(arr) {
  if (is.null(arr)) {
    return(NULL)
  } else {
    return(t(arr))
  }
}

update_short_name <- function(params, meta, short_name) {
  # NOTE: this expects meta from reference data and replaces short_name and
  # its occurance in the filename by string replacement.
  # A more general get_output_filename(meta) that reconstruct fname
  # would be better.
  fname <- strsplit(meta$filename, "/")[[1]]
  pcs <- strsplit(fname[length(fname)], "_")[[1]]
  pcs[which(pcs == meta$short_name)] <- short_name
  onam <- paste0(params$work_dir, "/", paste(pcs, collapse = "_"))
  new_meta <- meta
  new_meta$filename <- onam
  new_meta$short_name <- short_name
  return(new_meta)
}


whfcn <- function(x, ilow, ihigh){ # TODO: rename function
  return(length(which(x >= ilow & x < ihigh)))
}


write_nc_file_like <- function(
  params, meta, data, fillfloat,
  short_name="spei",
  units="1",
  long_name="Standardized Precipitation-Evapotranspiration Index",
  moty=FALSE
) {
  # loads a nc file from meta as reference to copy dimensions and attributes
  # from. Data, short_name and units from arguments will be used to update and
  # save the file. evspsblpot will be converted from monthly to cf units.
  # NOTE: cf convert skipped for debugging
  new_meta <- update_short_name(params, meta, short_name)
  print(paste("create new file:", new_meta$filename))
  ncid_in <- nc_open(meta$filename)
  xdim <- ncid_in$dim[["lon"]]
  if (length(xdim)==0) xdim <- ncid_in$dim[["longitude"]]
  ydim <- ncid_in$dim[["lat"]]
  if (length(ydim)==0) ydim <- ncid_in$dim[["latitude"]]
  if (moty) {
    tdim = ncdim_def("month", units, 1:12, longname="Month of the Year")
  } else {
    tdim <- ncid_in$dim[["time"]]
  }
  allatt <- ncatt_get(ncid_in, meta$short_name)
  var_data <- ncvar_def(short_name, units, list(xdim, ydim, tdim), fillfloat)
  idw <- nc_create(new_meta$filename, var_data, force_v4=TRUE)
  data[is.infinite(data)] <- fillfloat
  data[is.na(data)] <- fillfloat
  ncvar_put(idw, short_name, data)
  # NOTE: updated loop to use names instead of counter
  globat <- ncatt_get(ncid_in, 0)
  for (thisattname in names(globat)){
    ncatt_put(idw, 0, thisattname, globat[[thisattname]])
  }
  nc_close(idw)
  nc_close(ncid_in)
  return(new_meta$filename)
}
