# #############################################################################
# hyint_functions.R
#
# Author: Enrico Arnone (ISAC-CNR, Italy)
#
# #############################################################################
# Description
#     Functions used in HyInt routines
#
# Modification history
#    20170901-arnone_enrico: adapted to HyInt and extended
#    20170522-davini_paolo: creation for MiLES
# #############################################################################

# basis functions

##########################################################
#------------------------Packages------------------------#
##########################################################

# loading packages
library("maps")
library("ncdf4")
library("PCICt")

# check if fast linear fit is operative (after R 3.1):
# 3x faster than lm.fit, 36x faster than lm
if (exists(".lm.fit")) {
  lin.fit <- .lm.fit
} else {
  lin.fit <- lm.fit
}


##########################################################
#----------------Naming functions------------------------#
##########################################################

# function to flatten nested lists
flatten_lists <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten_lists)), recursive = FALSE))
}

getfilename_regridded <- function(spath, rgrid, var0, model_idx) {
  exp <- models_name[model_idx]
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  filename <- paste0(
    spath,
    "/",
    exp,
    "_",
    model_exp,
    "_",
    model_ens,
    "_",
    toString(year1),
    "-",
    toString(year2),
    "_",
    var0,
    "_",
    rgrid,
    ".nc"
  )
  return(filename)
}

getfilename_indices <-
  function(spath,
             label,
             model_idx,
             season,
             hist = F,
             hist_years = hist_years,
             grid = F,
             topo = F) {
    exp <- models_name[model_idx]
    model_exp <- models_experiment[model_idx]
    model_ens <- models_ensemble[model_idx]
    if (grid) {
      filename <- paste0(
        spath,
        "/",
        label,
        "_",
        exp,
        "_",
        model_exp,
        "_",
        model_ens,
        ".grid"
      )
    } else if (topo) {
      filename <- paste0(
        spath,
        "/",
        label,
        "_",
        exp,
        "_",
        model_exp,
        "_",
        model_ens,
        "_topo.nc"
      )
    } else {
      year1 <- models_start_year[model_idx]
      year2 <- models_end_year[model_idx]
      if (hist) {
        model_exp <- "historical"
        year1 <- hist_years[1]
        year2 <- hist_years[2]
      }
      filename <- paste0(
        spath,
        "/",
        label,
        "_",
        exp,
        "_",
        model_exp,
        "_",
        model_ens,
        "_",
        toString(year1),
        "_",
        toString(year2),
        "_",
        season,
        ".nc"
      )
    }
    return(filename)
  }

getfilename_etccdi <-
  function(spath, var, model_idx, yrmon = "yr") {
    # Function to get names of files of ETCCDI indices
    # If input 'var' is an array of names, 'filename' an array will be as well

    filename <- ""
    for (svar in var) {
      exp <- models_name[model_idx]
      model_exp <- models_experiment[model_idx]
      model_ens <- models_ensemble[model_idx]
      year1 <- toString(models_start_year[model_idx])
      year2 <- toString(models_end_year[model_idx])
      if (yrmon == "mon") {
        year1 <- paste0(year1, "01")
        year2 <- paste0(year2, "12")
      }
      filenametmp <- paste0(
        spath,
        "/",
        svar,
        "_",
        yrmon,
        "_",
        exp,
        "_",
        model_exp,
        "_",
        model_ens,
        "_",
        year1,
        "-",
        year2,
        ".nc"
      )
      filename <- c(filename, filenametmp)
    }
    filename <- filename[2:length(filename)]
    return(filename)
  }

getfilename_trends <- function(spath, label, model_idx, season) {
  exp <- models_name[model_idx]
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  filename <- paste0(
    spath,
    "/",
    diag_base,
    "_",
    exp,
    "_",
    model_exp,
    "_",
    model_ens,
    "_",
    toString(year1),
    "_",
    toString(year2),
    "_",
    season,
    "_tseries_",
    label,
    ".nc"
  )
  return(filename)
}

getfilename_figure <-
  function(spath,
             var,
             year1,
             year2,
             model_idx,
             season,
             syears,
             sregion,
             label,
             map,
             output_file_type,
             multimodel = F) {
    if (nchar(var) > 10) {
      var <- substr(var, 1, 10)
    }
    exp <- models_name[model_idx]
    year1 <- models_start_year[model_idx]
    year2 <- models_end_year[model_idx]
    model_exp <- models_experiment[model_idx]
    model_ens <- models_ensemble[model_idx]
    model_tag <- paste(exp, model_exp, model_ens, sep = "_")
    if (multimodel) {
      model_tag <- "multimodel"
    }
    figname <- paste0(
      spath,
      "/",
      paste(
        var,
        model_tag,
        paste(year1, year2, sep = "-"),
        season,
        syears,
        sregion,
        map,
        sep = "_"
      ),
      ".",
      output_file_type
    )
    if (!(label == "") & !(label == F)) {
      figname <- paste0(
        spath,
        "/",
        paste(
          var,
          model_tag,
          paste(year1, year2, sep = "-"),
          season,
          syears,
          sregion,
          label,
          map,
          sep = "_"
        ),
        ".",
        output_file_type
      )
    }
    return(figname)
  }


##########################################################
#-----------------Basic functions------------------------#
##########################################################

# normalize a time series
standardize <- function(timeseries) {
  out <-
    (timeseries - mean(timeseries, na.rm = T)) / sd(timeseries, na.rm = T)
  return(out)
}


# detect ics ipsilon lat-lon
whicher <- function(axis, number) {
  out <- which.min(abs(axis - number))
  return(out)
}


# area of longitude/latitude rectangle
area_lonlat <- function(lon1, lon2, lat1, lat2) {
  R <- 6378
  return(2 * pi * R^2 * abs(sin(lat1 / 180. * pi) - sin(lat2 / 180. * pi))
    * abs(lon1 - lon2) / 360)
}


# produce a 2d matrix of area size for given longitude/latitude grid points
area_size <- function(ics,
                      ipsilon,
                      resolution = NA,
                      norm = F) {
  if (is.na(resolution) &
    (length(ics) == 1) & (length(ipsilon) == 1)) {
    stop("Provide either resolution or two adjacent elements")
  }
  if (is.na(resolution) & (length(ics) != 1)) {
    resolution <- ics[2] - ics[1]
  }
  field <- array(NA, dim = c(length(ics), length(ipsilon)))
  for (j in seq_along(ipsilon)) {
    field[, j] <- area_lonlat(
      0,
      resolution,
      ipsilon[j] - 0.5 * resolution,
      ipsilon[j] + 0.5 * resolution
    )
  }
  if (norm) {
    field <- field / sum(field)
  }

  return(field)
}



# produce a 2d matrix of area weight
area_weight <- function(ics,
                        ipsilon,
                        root = T,
                        norm = F) {
  field <- array(NA, dim = c(length(ics), length(ipsilon)))
  if (root == T) {
    for (j in seq_along(ipsilon)) {
      field[, j] <- sqrt(cos(pi / 180 * ipsilon[j]))
    }
  }

  if (root == F) {
    for (j in 1:length(ipsilon)) {
      field[, j] <- cos(pi / 180 * ipsilon[j])
    }
  }
  if (norm) {
    field <- field / mean(field)
  }
  return(field)
}

# normalize a 2D or 3D field by a 2d matrix of area weight
area_weight_norm <-
  function(ics,
             ipsilon,
             field,
             root = T,
             norm = F) {
    timedim <- dim(field)[length(dim(field))]
    weights <- replicate(timedim, area_weight(ics, ipsilon,
      root = root,
      norm = norm
    ))
    field <- field * weights
    return(field)
  }

##########################################################
#--------------Time Based functions----------------------#
##########################################################

# check number of days for each month

number_days_month <- function(datas) {
  # evaluate the number of days in a defined month of a year
  datas <- as.Date(datas)
  m <- format(datas, format = "%m")
  while (format(datas, format = "%m") == m) {
    datas <- datas + 1
  }
  return(as.integer(format(datas - 1, format = "%d")))
}


# to convert season charname to months number
season2timeseason <- function(season) {
  if (season == "ALL") {
    timeseason <- 1:12
  }
  if (season == "JJA") {
    timeseason <- 6:8
  }
  if (season == "DJF") {
    timeseason <- c(1, 2, 12)
  }
  if (season == "MAM") {
    timeseason <- 3:5
  }
  if (season == "SON") {
    timeseason <- 9:11
  }
  if (!exists("timeseason")) {
    stop("wrong season selected!")
  }
  return(timeseason)
}

# leap year treu/false function
is_leapyear <- function(year) {
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

power_date_new <- function(datas) {
  whichdays <- as.numeric(format(datas, "%m"))
  # create a "season" for continuous time
  seas <- whichdays * 1
  ss <- 1
  for (i in 1:(length(whichdays) - 1)) {
    if (diff(whichdays)[i] > 1) {
      ss <- ss + 1
    }
    seas[i + 1] <- ss
  }
  etime <- list(
    day = as.numeric(format(datas, "%d")),
    month = as.numeric(format(datas, "%m")),
    year = as.numeric(format(datas, "%Y")),
    data = datas,
    season = seas
  )
  return(etime)
}

power_date <- function(season, ANNO1, ANNO2) {
  # evalute the number of days that will analyze in order
  # to create arrays of the needed dimensions

  # create continous calendar
  p1 <- as.Date(paste0(ANNO1, "-01-01"))
  p2 <- as.Date(paste0(ANNO2, "-12-31"))
  datas <- seq(p1, p2, by = "day")

  # select only days correspondeing to the needed season
  timeseason <- season2timeseason(season)
  month <- as.numeric(format(datas, "%m"))
  whichdays <- which(month %in% timeseason)

  # create a "season" for continuous time, used by persistance tracking
  seas <- whichdays * 1
  ss <- 1
  for (i in 1:(length(whichdays) - 1)) {
    if (diff(whichdays)[i] > 1) {
      ss <- ss + 1
    }
    seas[i + 1] <- ss
  }
  # produce a final timeseries of dates
  datas <- datas[whichdays]
  dataline <- list(
    day = as.numeric(format(datas, "%d")),
    month = as.numeric(format(datas, "%m")),
    year = as.numeric(format(datas, "%Y")),
    season = seas,
    data = datas
  )
  print("Time Array Built")
  print(paste("Length:", length(seas), "days for", season, "season"))
  print(paste("From", datas[1], "to", datas[length(seas)]))

  return(dataline)
}

power_date_no_leap <- function(season, ANNO1, ANNO2) {
  # apply to power_date object to clean out elements for leap years
  e <- power_date(season, ANNO1, ANNO2)
  leap.days <- which(e$month == 2 & e$day == 29)
  dataline.leap <- list(
    day = e$day[-leap.days],
    month = e$month[-leap.days],
    year = e$year[-leap.days],
    season = e$season[-leap.days],
    data = e$data[-leap.days]
  )
  print("FIXED FOR NO LEAP CALENDAR: Time Array Built")
  print(paste(
    "Length:",
    length(dataline.leap$season),
    "days for",
    season,
    "season"
  ))
  print(paste(
    "From", dataline.leap$data[1], "to",
    dataline.leap$data[length(dataline.leap$season)]
  ))
  return(dataline.leap)
}

power_date_30day <- function(season, ANNO1, ANNO2) {
  # apply to power.date object to clean out elements for leap years
  nmonths <- length(season2timeseason(season))
  nyears <- as.numeric(ANNO2) - as.numeric(ANNO1) + 1
  dd <- rep(seq(1, 30), nmonths * nyears)
  mm <- rep(rep(season2timeseason(season), each = 30), nyears)
  # create a "season" for continuous time, used by persistance tracking
  seas <- mm * 0 + 1
  ss <- 1
  for (i in 1:(length(mm) - 1)) {
    if (diff(mm)[i] > 1) {
      ss <- ss + 1
    }
    seas[i + 1] <- ss
  }
  dataline_30day <- list(
    day = dd,
    month = mm,
    season = seas
  )
  print("SIMPLIFIED CALENDAR FOR 30-day CALENDAR: Time Array Built")
  print(paste(
    "Length:",
    length(dataline.30day$season),
    "days for",
    season,
    "season"
  ))
  return(dataline_30day)
}

calc_region_timeseries <-
  function(x,
             y,
             indata,
             region,
             calc_sd = F,
             weighted_mean = T,
             root = F,
             norm = T,
             ..) {
    # This function subsets a lon/lat/time array based on an input
    # region(lon1,loni2,lat1,lat2) and returns its timeseries.
    # Area weights are applied if requested. The function returns also
    # the standard deviation of the averaging elements
    # (currently excluding weights)

    idimtimedata <- length(dim(indata))
    dimtimedata <- (dim(indata))[idimtimedata]
    retx <- which(region[1] <= x & x <= region[2])
    rety <- which(region[3] <= y & y <= region[4])
    if (!calc_sd) {
      print(paste(
        "Calc.region.timeseries: ",
        length(retx) * length(rety)
      ))
    }
    if (is.na(retx[1]) | is.na(rety[1])) {
      cat(
        "calc.region.timeseries: no data in selected region.",
        "Returning NA."
      )
      outdata <- array(dim = dimtimedata)
    } else {
      retdata <- indata[retx, rety, , drop = F]
      if (weighted_mean & !calc_sd) {
        retdata <- area_weight_norm(x[retx], y[rety], retdata,
          root = root, norm = norm
        )
      }
      outdata <- apply(retdata, idimtimedata, mean, na.rm = T)
      if (calc_sd) {
        outdata <- apply(retdata, idimtimedata, sd, na.rm = T)
      }
    }
    return(outdata)
  }

##################
#--------Data preprocessing
#################

#
# Method to create an asci grid file for
# to use to regrid on
# @param idx_dir path of directory containing
# files from which to create the grid
# Adapted from 20170920-sandstad_marit
#
create_grid <- function(ref_file = "./reffile",
                        path = idx_dir,
                        out_file = "./gridDef") {
  ## Picking the grid found in reference file to regrid over
  if (!file.exists(ref_file)) {
    ## Picking the grid found in the first file to regrid over
    ref_file <-
      list.files(path, pattern = "*.nc", full.names = TRUE)[1]
  }
  cdo("griddes", input = ref_file, stdout = out_file)
}

#
# Method to create a landSeaMask on a suitable grid
# @param regrid name w/path of gridfile to use
# to put the landdseamask on
# Adapted from 20170920-sandstad_marit
#
create_landseamask <-
  function(regrid = "./gridDef",
             ref_file = ref_file,
             loc = "./",
             regridded_topo = paste0("./", "regridded_topo.nc"),
             landmask = "./landSeaMask.nc",
             topo_only = F) {
    # Test if gridfile exists
    # otherwise call function to generate one
    if (!file.exists(regrid)) {
      if (length(ref_file) == 0) {
        print("Unable to access grid file")
        stop
      }
      create_grid(ref_file = ref_file, out_file = regrid)
    }

    ## Making topographic map
    ftopo <- cdo("topo", options = "-f nc")

    ## Regridding the topographic map to chosen grid
    cdo(
      "remapcon2",
      args = paste0("'", regrid, "'"),
      input = ftopo,
      output = regridded_topo
    )

    if (!topo_only) {
      # Set above sea-level gridpoints to missing
      ftopomiss1 <-
        cdo("setrtomiss", args = "0,9000", input = regridded_topo)

      # Set above sea-level gridpoints to 1
      ftopo1pos <- cdo("setmisstoc", args = "1", input = ftopomiss1)

      # Set below sea-level gridpoints to missing
      cdo("setrtomiss",
        args = "-9000,0",
        input = ftopo1pos,
        output = landmask
      )
      unlink(c(ftopomiss1, ftopo1pos))
    }
    unlink(ftopo)
  }

##
## Read seaLandElevationMask and mask data
##
apply_elevation_mask <- function(rfield,
                                 relevation,
                                 el_threshold,
                                 reverse = F) {
  if (!reverse) {
    if (el_threshold >= 0) {
      # mountains
      relevation[relevation < el_threshold] <- NA
      relevation <- relevation * 0 + 1
    } else {
      # oceans
      relevation[relevation > el_threshold] <- NA
      relevation <- relevation * 0 + 1
    }
  } else {
    if (el_threshold >= 0) {
      # mountains
      relevation[relevation > el_threshold] <- NA
      relevation <- relevation * 0 + 1
    } else {
      # oceans
      relevation[relevation < el_threshold] <- NA
      relevation <- relevation * 0 + 1
    }
  }
  itimedim <- dim(rfield)[length(dim(rfield))]
  myear_relevation <- replicate(itimedim, relevation)
  if (dim(myear_relevation) != dim(rfield)) {
    stop(
      "STOP - dimension of topography does not match
      dimension of field: remove old topography files if needed"
    )
  }
  rfield <- rfield * myear_relevation

  return(rfield)
}



##########################################################
#-------------------Data analysis------------------------#
##########################################################


###################################
# Function: Annual mean spell length
#
# About:  This function calculates the annual mean spell length of a given
#         field (lon x lat x time) reporting 1's for active parameter and
#         0's for non active parameter. In order to reduce memory usage only
#         the annual mean spell length is returned. E.g. calculation of dry
#         spell length needs input fields with 1 for dry days, 0 for wet ones.
#
# Author: E. Arnone ( ISAC-CNR, Torino)
# Last update: 14 June 2017

mean_spell_length <- function(m) {
  # Setup useful arrays and parameters
  nlon <- dim(m)[1]
  nlat <- dim(m)[2]
  ntime <- dim(m)[3]
  mean_spell_length_year <- m[, , 1] * NA

  # Loop through grid points
  for (ilon in 1:nlon) {
    for (ilat in 1:nlat) {
      spell_point <- (m[ilon, ilat, ])
      # Look for variations along time axis
      diff_spell_point <-
        spell_point[2:ntime] - spell_point[1:ntime - 1]
      # select when variation is positive (starting spell)
      spell_start <- which(diff_spell_point == 1) + 1
      if (!is.na(spell_point[1])) {
        if (spell_point[1] == 1) {
          spell_start <- c(1, spell_start)
        }
      } # if first day is active add it to list
      # select when variation is negative (ending spell)
      spell_stop <- which(diff_spell_point == -1)
      if (!is.na(spell_point[ntime])) {
        if (spell_point[ntime] == 1) {
          spell_stop <- c(spell_stop, ntime)
        }
      } # if last day is active add it to list
      # difference between stop and start gives spell length
      spell_length <- spell_stop - spell_start + 1
      # assign annual mean spell length to output array
      mean_spell_length_year[ilon, ilat] <-
        mean(spell_length, na.rm = T)
    }
  }
  return(mean_spell_length_year)
}

get_elevation <-
  function(filename = NULL,
             elev_range = c(-1000, 10000),
             mask = F,
             elev_plot = F) {
    # get elevation data from a high resolution topography file.

    funlink <- F
    if (is.null(filename)) {
      filename <- cdo("topo", options = "-f nc")
      funlink <- T
    }
    elevation <- ncdf_opener(
      filename,
      namevar = "elevation",
      namelon = "longitude",
      namelat = "latitude",
      rotate = "no"
    )
    lon_el <-
      ncdf_opener(filename, namevar = "longitude", rotate = "no")
    lat_el <-
      ncdf_opener(filename, namevar = "latitude", rotate = "no")
    elevation[which(elevation < elev_range[1] |
      elevation > elev_range[2])] <- NA
    if (mask) {
      elevation[which(elevation >= elev_range[1] &
        elevation <= elev_range[2])] <- 1
    }
    if (elev_plot) {
      filled_contour3(lon_el, lat_el, elevation, color.palette = rainbow)
      map(
        "world",
        regions = ".",
        interior = F,
        exact = F,
        boundary = T,
        add = T,
        col = "gray",
        lwd = 1.5
      )
    }
    el_list <-
      list(
        elevation = elevation,
        lon_el = lon_el,
        lat_el = lat_el
      )
    if (funlink) {
      unlink(filename)
    }
    return(el_list)
  }


##########################################################
#--------------NetCDF loading function-------------------#
##########################################################

# universal function to open a single var 3D (x,y,time) ncdf files: it includes
# rotation, y-axis filpping, time selection and CDO-based interpolation
# to replace both ncdf.opener.time and ncdf.opener (deprecated and removed)
# automatically rotate matrix to place greenwich at the center (flag "rotate")
# and flip the latitudes in order to have increasing
# if required (flag "interp2grid") additional interpolation with CDO is used.
# "grid" can be used to specify the target grid name
# time selection based on package PCICt must be specifed with both "tmonths"
#  and "tyears" flags. It returns a list including its own dimensions
ncdf_opener_universal <- # nolint
  function(namefile,
             namevar = NULL,
             namelon = NULL,
             namelat = NULL,
             tmonths = NULL,
             tyears = NULL,
             rotate = "full",
             interp2grid = F,
             grid = "r144x73",
             remap_method = "remapcon2",
             exportlonlat = TRUE,
             verbose = F) {
    # load package
    require(ncdf4)

    # verbose-only printing function
    printv <- function(value) {
      if (verbose) {
        print(value)
      }
    }

    # check if timeflag is activated or full file must be loaded
    if (is.null(tyears) | is.null(tmonths)) {
      timeflag <- FALSE
      printv("No time and months specified, loading all the data")
    } else {
      timeflag <- TRUE
      printv("tyears and tmonths are set!")
      require(PCICt)
    }

    if (rotate == "full") {
      rot <- T
      move1 <- move2 <- 1 / 2
    } # 180 degrees rotation of longitude
    if (rotate == "half") {
      rot <- T
      move1 <- 1 / 4
      move2 <- 3 / 4
    } # 90 degree rotation (useful for TM90)
    if (rotate == "no") {
      rot <- F
    } # keep as it is

    # interpolation made with CDO: second order conservative remapping
    if (interp2grid) {
      print(paste("Remapping with CDO on", grid, "grid"))
      if (is.null(namevar)) {
        namefile <- cdo(remap_method,
          args = paste0("'", grid, "'"),
          input = namefile
        )
      } else {
        selectf <- cdo("selvar", args = namevar, input = namefile)
        gridf <- tempfile()
        cdo("griddes", input = grid, stdout = gridf)
        namefile <- cdo(remap_method, args = gridf, input = selectf)
        unlink(c(selectf, gridf))
      }
    }

    # define rotate function (faster than with apply)
    rotation <- function(line) {
      vettore <- line
      dims <- length(dim(vettore))
      # for longitudes
      if (dims == 1) {
        ll <- length(line)
        line[(ll * move1):ll] <- vettore[1:(ll * move2 + 1)]
        line[1:(ll * move1 - 1)] <- vettore[(ll * move2 + 2):ll] - 360
      }
      # for x,y data
      if (dims == 2) {
        ll <- length(line[, 1])
        line[(ll * move1):ll, ] <- vettore[1:(ll * move2 + 1), ]
        line[1:(ll * move1 - 1), ] <- vettore[(ll * move2 + 2):ll, ]
      }
      # for x,y,t data
      if (dims == 3) {
        ll <- length(line[, 1, 1])
        line[(ll * move1):ll, , ] <- vettore[1:(ll * move2 + 1), , ]
        line[1:(ll * move1 - 1), , ] <-
          vettore[(ll * move2 + 2):ll, , ]
      }
      return(line)
    }

    # define flip function ('cos rev/apply is not working)
    flipper <- function(field) {
      dims <- length(dim(field))
      if (dims == 2) {
        ll <- length(field[1, ])
        field <- field[, ll:1]
      } # for x,y data
      if (dims == 3) {
        ll <- length(field[1, , 1])
        field <- field[, ll:1, ]
      } # for x,y,t data
      return(field)
    }

    # opening file: getting variable (if namevar is given, that variable
    # is extracted)
    printv(paste("opening file:", namefile))
    a <- nc_open(namefile)

    # if no name provided load the only variable available
    if (is.null(namevar)) {
      namevar <- names(a$var)
      if (length(namevar) > 1) {
        print(namevar)
        stop("More than one var in the files, please select it
             with namevar=yourvar")
      }
    }

    # load axis: updated version, looking for dimension directly stored
    # inside the variable
    naxis <-
      unlist(lapply(a$var[[namevar]]$dim, function(x)
        x["name"]))
    for (axis in naxis) {
      assign(axis, ncvar_get(a, axis))
      printv(paste(axis, ":", length(get(axis)), "records"))
    }

    if (timeflag) {
      printv("selecting years and months")

      # based on preprocessing of CDO time format: get calendar type and
      # use PCICt package for irregular data
      caldata <- ncatt_get(a, "time", "calendar")$value
      timeline <-
        as.PCICt(as.character(time), format = "%Y%m%d", cal = caldata)

      # break if the calendar has not been recognized
      if (any(is.na(timeline))) {
        stop("Calendar from NetCDF is unsupported or not present. Stopping!!!")
      }

      # break if the data requested is not there
      lastday_base <- paste0(max(tyears), "-", max(tmonths), "-28")
      maxdays <- number_days_month(lastday_base)
      if (caldata == "360_day") {
        maxdays <- 30
      }
      # uses number_days_month, which loops to get the month change
      lastday <- as.PCICt(paste0(
        max(tyears), "-", max(tmonths), "-",
        maxdays
      ),
      cal = caldata,
      format = "%Y-%m-%d"
      )
      firstday <-
        as.PCICt(paste0(min(tyears), "-", min(tmonths), "-01"),
          cal = caldata,
          format = "%Y-%m-%d"
        )
      if (max(timeline) < lastday | min(timeline) > firstday) {
        stop("You requested a time interval that is not present in the NetCDF")
      }
    }

    # time selection and variable loading
    printv("loading full field...")
    field <- ncvar_get(a, namevar)

    if (timeflag) {
      # select data we need
      select <- which(as.numeric(format(timeline, "%Y")) %in% tyears &
        as.numeric(format(timeline, "%m")) %in% tmonths)
      field <- field[, , select]
      time <- timeline[select]

      printv(paste("This is a", caldata, "calendar"))
      printv(paste(
        length(time), "days selected from", time[1],
        "to", time[length(time)]
      ))

      printv(paste("Months that have been loaded are.. "))
      printv(unique(format(time, "%Y-%m")))
    }

    # check for dimensions (presence or not of time dimension)
    dimensions <- length(dim(field))

    # if dimensions are multiple, get longitude, latitude
    # if needed, rotate and flip the array
    xlist <- c("lon", "Lon", "longitude", "Longitude")
    ylist <- c("lat", "Lat", "latitude", "Latitude")
    if (dimensions > 1) {
      # assign ics and ipsilon
      if (is.null(namelon)) {
        if (any(xlist %in% naxis)) {
          ics <- get(naxis[naxis %in% xlist], a$dim)$vals
        } else {
          print("WARNING: No lon found")
          ics <- NA
        }
      } else {
        ics <- ncvar_get(a, namelon)
      }
      if (is.null(namelat)) {
        if (any(ylist %in% naxis)) {
          ipsilon <- get(naxis[naxis %in% ylist], a$dim)$vals
        } else {
          print("WARNING: No lat found")
          ipsilon <- NA
        }
      } else {
        ipsilon <- ncvar_get(a, namelat)
      }

      # longitute rotation around Greenwich
      if (rot) {
        printv("rotating...")
        ics <- rotation(ics)
        field <- rotation(field)
      }
      if (ipsilon[2] < ipsilon[1] & length(ipsilon) > 1) {
        if (length(ics) > 1) {
          print("flipping...")
          ipsilon <- sort(ipsilon)
          field <- flipper(field)
        }
      }

      # exporting variables to the main program
      if (exportlonlat) {
        assign("ics", ics, envir = .GlobalEnv)
        assign("ipsilon", ipsilon, envir = .GlobalEnv)
      }
      assign(naxis[naxis %in% c(xlist, namelon)], ics)
      assign(naxis[naxis %in% c(ylist, namelat)], ipsilon)
    }

    if (dimensions > 3) {
      stop("This file is more than 3D file")
    }

    # close connection
    nc_close(a)

    # remove interpolated file
    if (interp2grid) {
      unlink(namefile)
    }

    # showing array properties
    printv(paste(dim(field)))
    if (timeflag) {
      printv(paste("From", time[1], "to", time[length(time)]))
    }

    # returning file list
    return(mget(c("field", naxis)))
  }

# ncdf.opener is a simplified wrapper for ncdf.opener.universal which returns
# only the field, ignoring the list
ncdf_opener <- function(namefile,
                        namevar = NULL,
                        namelon = NULL,
                        namelat = NULL,
                        tmonths = NULL,
                        tyears = NULL,
                        rotate = "full",
                        interp2grid = F,
                        grid = "r144x73",
                        remap_method = "remapcon2",
                        exportlonlat = T) {
  field <-
    ncdf_opener_universal(
      namefile,
      namevar,
      namelon,
      namelat,
      tmonths,
      tyears,
      rotate,
      interp2grid,
      grid,
      remap_method,
      exportlonlat = exportlonlat
    )
  return(field$field)
}


# function to open ncdf files (much more refined, with CDO-based interpolation)
ncdf_opener_time <- # nolint
  function(namefile,
             namevar = NULL,
             namelon = NULL,
             namelat = NULL,
             tmonths = NULL,
             tyears = NULL,
             ics = ics,
             ipsilon = ipsilon,
             rotate = "full",
             interp2grid = F,
             grid = "r144x73",
             remap_method = "remapcon2") {
    # function to open netcdf files. It uses ncdf4 library
    # time selection of month and years needed automatically rotate matrix
    # to place greenwich at the center (flag "rotate")
    # and flip the latitudes in order to have increasing
    # if require (flag "interp2grid") additional interpolation with CDO
    # can be used. "grid" can be used to specify the grid name
    require(ncdf4)
    require(PCICt)

    if (is.null(tyears) | is.null(tmonths)) {
      stop("Please specify both months and years to load")
    }

    if (rotate == "full") {
      rot <- T
      move1 <- move2 <- 1 / 2
    } # 180 degrees rotation of longitude
    if (rotate == "half") {
      rot <- T
      move1 <- 1 / 4
      move2 <- 3 / 4
    } # 90 degree rotation (useful for TM90)
    if (rotate == "no") {
      rot <- F
    } # keep as it is, breaking at Greemwich

    # interpolation made with CDO: second order conservative remapping
    if (interp2grid) {
      print(paste("Remapping with CDO on", grid, "grid"))
      namefile <- cdo(remap_method,
        args = paste0("'", grid, "'"),
        input = namefile
      )
    }

    # define rotate function (faster than with apply)
    rotation <- function(line) {
      vettore <- line
      dims <- length(dim(vettore))
      if (dims == 1) {
        # for longitudes
        ll <- length(line)
        line[(ll * move1):ll] <- vettore[1:(ll * move2 + 1)]
        line[1:(ll * move1 - 1)] <- vettore[(ll * move2 + 2):ll] - 360
      }
      if (dims == 2) {
        # for x,y data
        ll <- length(line[, 1])
        line[(ll * move1):ll, ] <- vettore[1:(ll * move2 + 1), ]
        line[1:(ll * move1 - 1), ] <- vettore[(ll * move2 + 2):ll, ]
      }
      if (dims == 3) {
        # for x,y,t data
        ll <- length(line[, 1, 1])
        line[(ll * move1):ll, , ] <- vettore[1:(ll * move2 + 1), , ]
        line[1:(ll * move1 - 1), , ] <-
          vettore[(ll * move2 + 2):ll, , ]
      }
      return(line)
    }

    # define flip function ('cos rev/apply is not working)
    flipper <- function(field) {
      dims <- length(dim(field))
      if (dims == 2) {
        ll <- length(field[1, ])
        field <- field[, ll:1]
      } # for x,y data
      if (dims == 3) {
        ll <- length(field[1, , 1])
        field <- field[, ll:1, ]
      } # for x,y,t data
      return(field)
    }


    # opening file: getting variable (if namevar is given,
    # that variable is extracted)
    print(paste("opening file:", namefile))
    a <- nc_open(namefile)

    # load axis: old version, loading the variable dimensions with a max of
    # 4 dimensions. It showed some issues with the time_bnds variable appearing
    # in some NetCDF file. naxis=names(a$dim)[1:min(c(4,length(a$dim)))]
    # load axis: updated version, looking for dimension directly stored inside
    #  the variable
    naxis <-
      unlist(lapply(a$var[[namevar]]$dim, function(x)
        x["name"]))
    for (axis in naxis) {
      print(axis)
      assign(axis, ncvar_get(a, axis))
    }
    # based on preprocessing of CDO time format: get calendar type and
    # use PCICt package for irregular data
    caldata <- ncatt_get(a, "time", "calendar")$value
    timeline <-
      as.PCICt(as.character(time), format = "%Y%m%d", cal = caldata)
    str(timeline)

    # break if the calendar has not been recognized
    if (any(is.na(timeline))) {
      stop("Calendar from NetCDF is unsupported or not present. Stopping!!!")
    }

    # break if the data requested is not there
    lastday_base <- paste0(max(tyears), "-", max(tmonths), "-28")
    # uses number_days_month, which loops to get the month change
    lastday <- as.PCICt(
      paste0(
        max(tyears),
        "-",
        max(tmonths),
        "-",
        number_days_month(lastday_base)
      ),
      cal = caldata,
      format = "%Y-%m-%d"
    )
    firstday <-
      as.PCICt(paste0(min(tyears), "-", min(tmonths), "-01"),
        cal = caldata,
        format = "%Y-%m-%d"
      )
    if (max(timeline) < lastday | min(timeline) > firstday) {
      stop("You requested a time interval that is not present in the NetCDF")
    }

    # time selection and variable loading
    # if no name provided load the only variable available
    if (is.null(namevar)) {
      namevar <- names(a$var)
    }
    field <- ncvar_get(a, namevar)

    # select data we need
    select <- which(as.numeric(format(timeline, "%Y")) %in% tyears &
      as.numeric(format(timeline, "%m")) %in% tmonths)

    field <- field[, , select]
    time <- timeline[select]

    # check for dimensions (presence or not of time dimension)
    dimensions <- length(dim(field))

    # if dimensions are multiple, get longitude, latitude
    # if needed, rotate and flip the array
    if (dimensions > 1) {
      # assign ics and ipsilon
      if (is.null(namelon)) {
        xlist <- c("lon", "Lon", "longitude", "Longitude")
        if (any(xlist %in% naxis)) {
          ics <- get(naxis[(naxis %in% xlist)], a$dim)$vals
        } else {
          stop("No lon found")
        }
      } else {
        ics <- ncvar_get(a, namelon)
      }
      if (is.null(namelat)) {
        ylist <- c("lat", "Lat", "latitude", "Latitude")
        if (any(ylist %in% naxis)) {
          ipsilon <- get(naxis[(naxis %in% ylist)], a$dim)$vals
        } else {
          stop("No lon found")
        }
      } else {
        ipsilon <- ncvar_get(a, namelat)
      }

      print("flipping and rotating")
      # longitute rotation around Greenwich
      if (rot) {
        ics <- rotation(ics)
        field <- rotation(field)
      }
      if (ipsilon[2] < ipsilon[1] & length(ipsilon) > 1) {
        if (length(ics) > 1) {
          ipsilon <- sort(ipsilon)
          field <- flipper(field)
        }
      }

      # exporting variables to the main program
      assign("ics", ics, envir = .GlobalEnv)
      assign("ipsilon", ipsilon, envir = .GlobalEnv)
      assign(naxis[naxis %in% xlist], ics)
      assign(naxis[naxis %in% ylist], ipsilon)
    }


    if (dimensions > 3) {
      stop("This file is more than 3D file")
    }

    # close connection
    nc_close(a)

    # remove interpolated file
    if (interp2grid) {
      unlink(namefile)
    }

    # showing array properties
    print(paste(dim(field)))
    print(paste("From", time[1], "to", time[length(time)]))

    return(mget(c("field", naxis)))
  }


##########################################################
#--------------Plotting functions------------------------#
##########################################################


# Figure functions
scale_figure <- function(plot_type,
                         diag_script_cfg,
                         nfields,
                         npancol,
                         npanrow) {
  source(diag_script_cfg)
  if (plot_type == 1 || plot_type == 11) {
    npancol <- 1
    npanrow <- 1
  }
  if (plot_type == 2) {
    npancol <- 1
    npanrow <- 3
  }
  if (plot_type == 3) {
    npancol <- 3
    napnrow <- nfields
  }
  npanels <- npancol * npanrow
  if (npancol > 1) {
    png_width <- png_width_multi * npancol
    pdf_width <- pdf_width_multi * npancol
    x11_width <- x11_width_multi * npancol
  }
  png_width <- png_width * figure_rel_width[plot_type]
  pdf_width <- pdf_width * figure_rel_width[plot_type]
  x11_width <- x11_width * figure_rel_width[plot_type]

  figure_aspect_ratio[plot_type] <- (figure_aspect_ratio[plot_type]
  * npancol / npanrow)

  plot_size <-
    c(png_width, png_width / figure_aspect_ratio[plot_type])
  if (tolower(output_file_type) == "pdf") {
    plot_size[1] <- pdf_width
    plot_size[2] <- pdf_width / figure_aspect_ratio[plot_type]
  } else if ((tolower(output_file_type) == "eps") |
    (tolower(output_file_type) == "epsi") |
    (tolower(output_file_type) == "ps")) {
    plot_size[1] <- pdf_width
    plot_size[2] <- pdf_width / figure_aspect_ratio[plot_type]
  } else if (tolower(output_file_type) == "x11") {
    plot_size[1] <- x11_width
    plot_size[2] <- x11_width / figure_aspect_ratio[plot_type]
  }
  print(plot_size)
  return(plot_size)
}

graphics_startup <- function(figname, output_file_type, plot_size) {
  source(diag_script_cfg)
  # choose output format for figure - by JvH
  if (tolower(output_file_type) == "png") {
    png(
      filename = figname,
      width = plot_size[1],
      height = plot_size[2]
    )
  } else if (tolower(output_file_type) == "pdf") {
    pdf(
      file = figname,
      width = plot_size[1],
      height = plot_size[2],
      onefile = T
    )
  } else if ((tolower(output_file_type) == "eps") |
    (tolower(output_file_type) == "epsi") |
    (tolower(output_file_type) == "ps")) {
    setEPS(
      width = plot_size[1],
      height = plot_size[2],
      onefile = T,
      paper = "special"
    )
    postscript(figname)
  } else if (tolower(output_file_type) == "x11") {
    x11(width = plot_size[1], height = plot_size[2])
  }
  return()
}

graphics_close <- function(figname) {
  print(figname)
  dev.off()
  return()
}

# extensive filled.contour function
filled_contour3 <- # nolint
  function(x = seq(0, 1, length.out = nrow(z)),
             y = seq(0, 1, length.out = ncol(z)),
             z,
             xlim = range(x, finite = TRUE),
             ylim = range(y, finite = TRUE),
             zlim = range(z, finite = TRUE),
             levels = pretty(zlim, nlevels),
             nlevels = 20,
             color.palette = cm.colors,
             col = color.palette(length(levels) - 1),
             extend = TRUE,
             plot.title,
             plot.axes,
             key.title,
             key.axes,
             asp = NA,
             xaxs = "i",
             yaxs = "i",
             las = 1,
             axes = TRUE,
             frame.plot = axes,
             mar,
             ...) {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    # modification to allow plot outside boundaries

    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else {
        stop("no 'z' matrix specified")
      }
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) {
      stop("increasing 'x' and 'y' values expected")
    }

    # trim extremes for nicer plots
    if (extend) {
      z[z < min(levels)] <- min(levels)
      z[z > max(levels)] <- max(levels)
    }

    plot.new()
    plot.window(xlim,
      ylim,
      "",
      xaxs = xaxs,
      yaxs = yaxs,
      asp = asp
    )
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) {
      stop("no proper 'z' matrix specified")
    }
    if (!is.double(z)) {
      storage.mode(z) <- "double"
    }
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
      col = col
    )
    if (missing(plot.axes)) {
      if (axes) {
        title(
          main = "",
          xlab = "",
          ylab = ""
        )
        Axis(x, side = 1, ...)
        Axis(y, side = 2, ...)
      }
    }
    else {
      plot.axes
    }
    if (frame.plot) {
      box()
    }
    if (missing(plot.title)) {
      title(...)
    } else {
      plot.title
    }
    invisible()
  }

image_scale3 <- function(z,
                         levels,
                         color.palette = heat.colors,
                         colorbar.label = "image.scale",
                         extend = T,
                         line.label = 2,
                         line.colorbar = 0,
                         cex.label = 1,
                         cex.colorbar = 1,
                         colorbar.width = 1,
                         new_fig_scale = c(-0.07, -0.03, 0.1, -0.1),
                         ...) {
  # save properties from main plotting region
  old.par <- par(no.readonly = TRUE)
  mfg.save <- par()$mfg
  old.fig <- par()$fig

  # defining plotting region with proper scaling
  xscal <- (old.fig[2] - old.fig[1])
  yscal <- (old.fig[4] - old.fig[3])
  lw <- colorbar.width
  lp <- line.colorbar / 100
  new.fig <- c(
    old.fig[2] + new_fig_scale[1] * xscal * lw - lp,
    old.fig[2] + new_fig_scale[2] * xscal - lp,
    old.fig[3] + new_fig_scale[3] * yscal,
    old.fig[4] + new_fig_scale[4] * yscal
  )

  if (missing(levels)) {
    levels <- seq(min(z), max(z), , 12)
  }
  # fixing color palette
  col <- color.palette(length(levels) - 1)

  # starting plot
  par(
    mar = c(1, 1, 1, 1),
    fig = new.fig,
    new = TRUE
  )

  # creating polygons for legend
  poly <- vector(mode = "list", length(col))
  for (i in seq(poly)) {
    poly[[i]] <- c(levels[i], levels[i + 1], levels[i + 1], levels[i])
  }

  xlim <- c(0, 1)
  if (extend) {
    longer <- 1.5
    dl <- diff(levels)[1] * longer
    ylim <- c(min(levels) - dl, max(levels) + dl)
  } else {
    ylim <- range(levels)
  }
  plot(
    1,
    1,
    t = "n",
    ylim = ylim,
    xlim = xlim,
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i",
    ...
  )
  for (i in seq(poly)) {
    polygon(c(0, 0, 1, 1), poly[[i]], col = col[i], border = NA)
  }
  if (extend) {
    polygon(c(0, 1, 1 / 2),
      c(levels[1], levels[1], levels[1] - dl),
      col = col[1],
      border = NA
    )
    polygon(c(0, 1, 1 / 2),
      c(
        levels[length(levels)], levels[length(levels)],
        levels[length(levels)] + dl
      ),
      col = col[length(col)],
      border = NA
    )
    polygon(
      c(0, 0, 1 / 2, 1, 1, 1 / 2),
      c(
        levels[1], levels[length(levels)], levels[length(levels)] + dl,
        levels[length(levels)], levels[1], levels[1] - dl
      ),
      border = "black",
      lwd = 2
    )
    ylim0 <- range(levels)
    prettyspecial <- pretty(ylim0)
    prettyspecial <- prettyspecial[prettyspecial <= max(ylim0) &
      prettyspecial >= min(ylim0)]
    axis(
      4,
      las = 1,
      cex.axis = cex.colorbar,
      at = prettyspecial,
      labels = prettyspecial,
      ...
    )
  } else {
    box()
    axis(4, las = 1, cex.axis = cex.colorbar, ...)
  }

  # box, axis and leged
  mtext(colorbar.label,
    line = line.label,
    side = 4,
    cex = cex.label,
    ...
  )

  # resetting properties for starting a new plot (mfrow style)
  par(old.par)
  par(mfg = mfg.save, new = FALSE)
  invisible()
}

cdo <-
  function(command,
             args = "",
             input = "",
             options = "",
             output = "",
             stdout = "",
             noout = F) {
    if (args != "") {
      args <- paste0(",", args)
    }
    if (stdout != "") {
      stdout <- paste0(" > '", stdout, "'")
      noout <- T
    }
    if (input[1] != "") {
      for (i in seq_along(input)) {
        input[i] <- paste0("'", input[i], "'")
      }
      input <- paste(input, collapse = " ")
    }
    output0 <- output
    if (output != "") {
      output <- paste0("'", output, "'")
    } else if (!noout) {
      output <- tempfile()
      output0 <- output
    }
    argstr <-
      paste0(
        options, " ", command, args, " ", input, " ", output,
        " ", stdout
      )
    print(paste("cdo", argstr))
    ret <- system2("cdo", args = argstr)
    if (ret != 0) {
      stop(paste("Failed (", ret, "): cdo", argstr))
    }
    return(output0)
  }
