# basis functions


##########################################################
#------------------------Packages------------------------#
##########################################################

# loadin packages
library("maps")
library("ncdf4")
library("PCICt")

# check if fast linear fit is operative (after R 3.1): 3x faster
# than lm.fit, 36x faster than lm
if (exists(".lm.fit")) {
  lin.fit <- .lm.fit
} else {
  lin.fit <- lm.fit
}

# check R version as numeric
R_version <-
  as.numeric(R.Version()$major) + as.numeric(R.Version()$minor) / 10

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

# produce a 2d matrix of area weight
area_weight <- function(ics, ipsilon, root = T) {
  field <- array(NA, dim = c(length(ics), length(ipsilon)))
  if (root == T) {
    for (j in seq_along(ipsilon)) {
      field[, j] <- sqrt(cos(pi / 180 * ipsilon[j]))
    }
  }

  if (root == F) {
    for (j in seq_along(ipsilon)) {
      field[, j] <- cos(pi / 180 * ipsilon[j])
    }
  }

  return(field)
}

# sector details for blocking extra diagnostics and EOFs sectors
sector_details <- function(SECTOR) {
  if (SECTOR == "Euro") {
    lons <- c(-15, 25)
    lats <- c(50, 65)
    namesec <- "Central Europe"
  }
  if (SECTOR == "Azores") {
    lons <- c(-70, -10)
    lats <- c(30, 40)
    namesec <- "Central Atlantic"
  }
  if (SECTOR == "Greenland") {
    lons <- c(-65, -15)
    lats <- c(62.5, 72.5)
    namesec <- "Greenland"
  }
  if (SECTOR == "FullPacific") {
    lons <- c(130, -150)
    lats <- c(60, 75)
    namesec <- "North Pacific"
  }
  if (SECTOR == "FullPacific2") {
    lons <- c(130, 210)
    lats <- c(60, 75)
    namesec <- "North Pacific"
  }

  left1 <- which.min(abs(ics - lons[1]))
  right1 <- which.min(abs(ics - lons[2]))
  low1 <- which.min(abs(ipsilon - lats[1]))
  high1 <- which.min(abs(ipsilon - lats[2]))

  latssel <- low1:high1
  if (SECTOR == "FullPacific") {
    lonssel <- c(left1:length(ics), 1:right1)
  } else {
    lonssel <- left1:right1
  }
  out <- list(
    lons = lons,
    lonssel = lonssel,
    lats = lats,
    latssel = latssel,
    name = namesec
  )
  return(out)
}

# weighted correlation
weighted_cor <- function(x, y, w) {
  w_mean_x <- sum(w * x) / sum(w)
  w_mean_y <- sum(w * y) / sum(w)

  w_cov_xy <- sum(w * (x - w_mean_x) * (y - w_mean_y)) / sum(w)
  w_var_y <- sum(w * (y - w_mean_y) * (y - w_mean_y)) / sum(w)
  w_var_x <- sum(w * (x - w_mean_x) * (x - w_mean_x)) / sum(w)

  corr <- w_cov_xy / sqrt(w_var_x * w_var_y)
  return(corr)
}

# weighted standard deviations
weighted_sd <- function(x, w) {
  w_mean <- sum(w * x) / sum(w)
  v1 <- sum(w)
  v2 <- sum(w^2)
  var <- v1 / (v1^2 - v2) * sum(w * (x - w_mean)^2)
  sdd <- sqrt(var)
  return(sdd)
}

# info string creator
info_builder <-
  function(dataset, expid, ens, year1, year2, season) {
    # loop on descriptors that are concatenated to create info string
    descriptors <-
      c(dataset, expid, ens, paste0(year1, "-", year2), season)
    info <- NULL
    for (dcode in descriptors) {
      if (length(dcode) > 0) {
        info <- paste(info, dcode)
      }
    }
    return(info)
  }

# basic switch to create NetCDF file names and folders
# (use recursive structure from v0.6)
file_builder <-
  function(DATADIR,
             dir_name,
             file_name,
             dataset,
             expid,
             ens,
             year1,
             year2,
             season) {
    # loop on descriptors that are concatenated to create dir and file name
    descriptors <-
      c(dataset, expid, ens, paste0(year1, "-", year2), season)
    for (dcode in descriptors) {
      if (length(dcode) > 0) {
        DATADIR <- file.path(DATADIR, dcode)
        file_name <- paste0(file_name, "_", dcode)
      }
    }

    # add directory name descriptor
    DATADIR <- file.path(DATADIR, dir_name)

    # actually dir.exists is in devtools only for R < 3.2,
    # then is included in base package
    if (exists("dir.exists")) {
      if (!dir.exists(DATADIR)) {
        dir.create(DATADIR, recursive = T)
      }
    } else {
      dir.create(DATADIR, recursive = T, showWarnings = F)
    }
    return(file.path(DATADIR, paste0(file_name, ".nc")))
  }

# basic switch to create figures names and folders
# (use recursive structure from v0.6)
fig_builder <- function(FIGDIR,
                        dir_name,
                        file_name,
                        dataset,
                        expid,
                        ens,
                        year1,
                        year2,
                        season,
                        output_file_type) {
  # loop on descriptors that are concatenated to create dir and file name
  descriptors <-
    c(dataset, expid, ens, paste0(year1, "-", year2), season)
  for (dcode in descriptors) {
    if (dcode != "NO") {
      FIGDIR <- file.path(FIGDIR, dcode)
      file_name <- paste0(file_name, "_", dcode)
    }
  }

  # add directory name descriptor
  FIGDIR <- file.path(FIGDIR, dir_name)

  # actually dir.exists is in devtools only for R < 3.2,
  # then is included in base package
  if (exists("dir.exists")) {
    if (!dir.exists(FIGDIR)) {
      dir.create(FIGDIR, recursive = T)
    }
  } else {
    dir.create(FIGDIR, recursive = T, showWarnings = F)
  }

  return(file.path(FIGDIR, paste0(file_name, ".", output_file_type)))
}

# progression bar
progression_bar <- function(index, total_length, each = 10) {
  if (any(index == round(seq(0, total_length, , each + 1)))) {
    progression <- paste("--->", round(index / total_length * 100), "%")
    print(progression)
  }
}


##########################################################
#--------------Time Based functions----------------------#
##########################################################

# to convert season charname to months number
season2timeseason <- function(season) {
  if (nchar(season) == 3 & toupper(season) == season) {
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
  } else {
    charseason <- strsplit(season, "_")[[1]]
    print(charseason)
    if (mean(nchar(charseason)) == 3) {
      timeseason <- which(charseason == month.abb)
    } else {
      timeseason <- which(charseason == month.name)
    }
  }
  print(timeseason)
  if (length(timeseason) == 0 |
    min(timeseason) < 0 | max(timeseason) > 13) {
    stop("wrong season selected!")
  }
  return(timeseason)
}

# leap year treu/false function
is_leapyear <- function(year) {
  return(((year %% 4 == 0) &
    (year %% 100 != 0)) | (year %% 400 == 0))
}

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

power_date_new <- function(datas) {
  whichdays <- as.numeric(format(datas, "%m"))
  # create a "season" for continuous time, used by persistance tracking
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
  print("Time Array Built")
  print(paste("Length:", length(seas)))
  print(paste("From", datas[1], "to", datas[length(seas)]))
  return(etime)
}

##########################################################
#--------------NetCDF loading function-------------------#
##########################################################

# universal function to open a single var 3D (x,y,time) ncdf files:
# it includes rotation, y-axis filpping, possible time selection and
# CDO-based interpolation to replace both ncdf.opener.time and ncdf.opener
# (deprecated and removed)
# automatically rotate matrix to place greenwich at the center (flag "rotate")
# and flip the latitudes in order to have increasing
# if required (flag "interp2grid") additional interpolation with CDO can be
# used. "grid" can be used to specify the target grid name
# time selection based on package PCICt must be specifed with
# both "tmonths" and "tyears" flags
# it returns a list including its own dimensions
ncdf_opener_universal <- # nolint
  function(namefile,
             namevar = NULL,
             namelon = NULL,
             namelat = NULL,
             tmonths = NULL,
             tyears = NULL,
             rotate = "full",
             interp2grid = FALSE,
             fillmiss = FALSE,
             grid = "r144x73",
             remap_method = "remapcon2",
             exportlonlat = TRUE,
             verbose = TRUE) {
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

    if (fillmiss) {
        namefile <- cdo("fillmiss", input = namefile)
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
      # for x,y,z,t data
      if (dims == 4) {
        ll <- length(line[, 1, 1, 1])
        line[(ll * move1):ll, , , ] <-
          vettore[1:(ll * move2 + 1), , , ]
        line[1:(ll * move1 - 1), , , ] <-
          vettore[(ll * move2 + 2):ll, , , ]
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
      if (dims == 4) {
        ll <- length(field[1, , 1, 1])
        field <- field[, ll:1, , ]
      } # for x,y,z,t data
      return(field)
    }

    # opening file: getting variable (if namevar is given, that variable
    # is extracted)
    printv(paste("opening file:", namefile))
    a <- nc_open(namefile)
    print(paste("Loading", namevar, "..."))

    # if no name provided load the only variable available
    if (is.null(namevar)) {
      namevar <- names(a$var)
      if (length(namevar) > 1) {
        print(namevar)
        stop(
          paste("More than one var in the files, please select it with"),
          " namevar=yourvar",
          sep = ""
        )
      }
    }

    # load axis: updated version, looking for dimension directly
    # stored inside the variable
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
      units <- ncatt_get(a, "time", "units")$value
      caldata <- ncatt_get(a, "time", "calendar")$value
      if (grepl("day as", units, fixed = TRUE) |
        grepl("days as", units, fixed = TRUE)) {
        timeline <- as.PCICt(as.character(time),
          format = "%Y%m%d",
          cal = caldata
        )
      } else if (grepl("day since", units, fixed = TRUE) |
        grepl("days since", units, fixed = TRUE)) {
        origin <- unlist(strsplit(units, "[a-zA-Z ]+"))[2]
        origin.pcict <-
          as.PCICt(origin, cal = caldata, format = "%Y-%m-%d")
        timeline <- origin.pcict + (floor(time) * 86400)
      } else {
        printv(units)
        stop("Time units from NetCDF unsupported. Stopping!!!")
      }

      # break if the calendar has not been recognized
      if (any(is.na(timeline))) {
        stop("Calendar from NetCDF is unsupported or not present. Stopping!!!")
      }

      # break if the data requested is not there
      lastday_base <- paste0(max(tyears), "-", max(tmonths), "-28")
      # uses number.days.month, which loops to get the month change
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
        print(firstday)
        print(lastday)
        print(min(timeline))
        print(max(timeline))
        stop("You requested a time interval that is not present in the NetCDF")
      }
    }

    # time selection and variable loading
    printv("loading full field...")
    field <- ncvar_get(a, namevar)

    if (timeflag) {
      # select data we need
      select <- which(as.numeric(format(timeline, "%Y")) %in%
        tyears & as.numeric(format(timeline, "%m")) %in% tmonths)
      field <- field[, , select]
      time <- timeline[select]

      printv(paste("This is a", caldata, "calendar"))
      printv(paste(
        length(time), "days selected from", time[1], "to",
        time[length(time)]
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
          printv("WARNING: No lon found")
          ics <- NA
        }
      } else {
        ics <- ncvar_get(a, namelon)
      }
      if (is.null(namelat)) {
        if (any(ylist %in% naxis)) {
          ipsilon <- get(naxis[naxis %in% ylist], a$dim)$vals
        } else {
          printv("WARNING: No lat found")
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
      # if ics and ipsilon exists, assign the rearranged values
      if (!is.na(ics[1])) {
        assign(naxis[naxis %in% c(xlist, namelon)], ics)
      }
      if (!is.na(ipsilon[1])) {
        assign(naxis[naxis %in% c(ylist, namelat)], ipsilon)
      }
    }

    if (dimensions > 4) {
      stop("This file is more than 4D file")
    }

    # close connection
    nc_close(a)

    # remove interpolated file
    if (interp2grid) {
      system2("rm", tempfile)
    }

    # showing array properties
    printv(paste(dim(field)))
    if (timeflag) {
      printv(paste("From", time[1], "to", time[length(time)]))
    }

    # returning file list
    return(mget(c("field", naxis)))
  }

# ncdf_opener is a simplified wrapper for ncdf_opener_universal which returns
# only the field, ignoring the list and no verbosity
ncdf_opener <- function(namefile,
                        namevar = NULL,
                        namelon = NULL,
                        namelat = NULL,
                        tmonths = NULL,
                        tyears = NULL,
                        rotate = "full",
                        interp2grid = FALSE,
                        fillmiss = FALSE,
                        grid = "r144x73",
                        remap_method = "remapcon2",
                        exportlonlat = TRUE,
                        verbose = FALSE) {
  field <- ncdf_opener_universal(
    namefile,
    namevar,
    namelon,
    namelat,
    tmonths,
    tyears,
    rotate,
    interp2grid,
    fillmiss,
    grid,
    remap_method,
    exportlonlat,
    verbose
  )
  return(field$field)
}

##########################################################
#--------------Plotting functions------------------------#
##########################################################

# function to open devices
open_plot_device <-
  function(figname, output_file_type, special = FALSE) {
    # Choose output format for figure
    output_file_type <- tolower(output_file_type)
    if (special == FALSE) {
      if (output_file_type == "png") {
        png(
          filename = figname,
          width = png_width,
          height = png_height
        )
      } else if (output_file_type == "pdf") {
        pdf(
          file = figname,
          width = pdf_width,
          height = pdf_height,
          onefile = T
        )
      } else if ((output_file_type == "eps") |
        (output_file_type == "epsi") |
        (output_file_type == "ps")) {
        setEPS(
          width = pdf_width,
          height = pdf_height,
          onefile = T,
          paper = "special"
        )
        postscript(figname)
      }
    } else {
      if (output_file_type == "png") {
        png(
          filename = figname,
          width = png_width / af,
          height = png_height * af / 2
        )
      } else if (output_file_type == "pdf") {
        pdf(
          file = figname,
          width = pdf_width / af,
          height = pdf_height * af / 2,
          onefile = T
        )
      } else if ((output_file_type == "eps") |
        (output_file_type == "epsi") |
        (output_file_type == "ps")) {
        setEPS(
          width = pdf_width / af,
          height = pdf_height * af / 2,
          onefile = T,
          paper = "special"
        )
        postscript(figname)
      }
    }
  }


# extensive filled_contour function
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
    # modification by Ian Taylor of the filled_contour function
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
    old.fig[2] - 0.07 * xscal * lw - lp,
    old.fig[2] - 0.03 * xscal - lp,
    old.fig[3] + 0.1 * yscal,
    old.fig[4] - 0.1 * yscal
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
        levels[1], levels[length(levels)],
        levels[length(levels)] + dl, levels[length(levels)], levels[1],
        levels[1] - dl
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

# function for interpolation and projection of a 2D field on a
# mapproj R projection
proj_plot <-
  function(lon,
             lat,
             field,
             lmin = NULL,
             proj = "azequalarea",
             param = NULL,
             orient = c(90, 0, 0),
             npoints = 201) {
    # default is azimuthal equal area map

    # required packages
    require(mapproj)
    require(akima)

    # it provides lower latitude limit for plots
    if (is.null(lmin)) {
      lmin <- min(lat)
    }

    # build grids
    lon.grid <- rep(lon, length(lat))
    lat.grid <- sort(rep(ipsilon, length(lon)))

    # project grid
    proj.grid <- mapproject(
      lon.grid,
      lat.grid,
      projection = proj,
      parameters = param,
      orientation = orient
    )

    # provide limits for future plots (for polar projection)
    limiter <- mapproject(c(0, 90, 180, 270),
      rep(lmin, 4),
      proj = "",
      orientation = orient
    )
    xlims <- sort(c(limiter$x[2], limiter$x[4]))
    ylims <- sort(c(limiter$y[1], limiter$y[3]))

    # plot grid
    lon.plot <-
      seq(min(proj.grid$x, na.rm = T),
        max(proj.grid$x, na.rm = T),
        length.out = npoints
      )
    lat.plot <-
      seq(min(proj.grid$y, na.rm = T),
        max(proj.grid$y, na.rm = T),
        length.out = npoints
      )

    # interpolation (akima needed)
    good <-
      is.finite(field) & is.finite(proj.grid$x) & is.finite(proj.grid$y)
    projected <-
      interp(proj.grid$x[good],
        proj.grid$y[good],
        field[good],
        lon.plot,
        lat.plot,
        duplicate = "strip"
      )
    return(projected = list(
      x = projected$x,
      y = projected$y,
      z = projected$z,
      xlim = xlims,
      ylim = ylims
    ))
  }

# addland function based on map which can handle projections
proj_addland <- function(proj = "no",
                         orient = c(90, 0, 0),
                         param = NULL,
                         color = "black") {
  # required packages
  require(maps)
  require(mapproj)

  if (proj == "no") {
    map(
      "world",
      regions = ".",
      interior = F,
      exact = F,
      boundary = T,
      add = T
    )
  } else {
    # get map, project and do the lines
    box()
    map(
      "world",
      add = T,
      projection = proj,
      orientation = orient,
      parameter = param,
      interior = F,
      exact = F,
      boundary = T
    )

    # default lines for northern hemisphere
    for (i in seq(-80, 80, 20)) {
      x0 <- ics
      y0 <- rep(i, length(ics))
      p <- mapproject(x0, y0, proj = "", orientation = orient)
      lines(p, lty = 3)
    }

    # default circles for northern hemisphere
    for (i in c(seq(-360, 360, 30))) {
      y0 <- seq(min(ipsilon), max(ipsilon), , 90)
      x0 <- rep(i, 90)
      p <- mapproject(x0, y0, proj = "", orientation = orient)
      lines(p, lty = 3)
    }
  }
}

# rearrange arrays for use both standard plotting and proj_plot
plot_prepare <- function(ics, ipsilon, field, proj, lat_lim) {
  if (proj == "no") {
    outfile <- list(
      x = ics,
      y = ipsilon,
      z = field,
      xlim = range(ics),
      ylim = lat_lim,
      xlab = "Longitude",
      ylab = "Latitude",
      axes = T
    )
  } else {
    field[is.na(field)] <- 0
    p <- proj_plot(
      ics,
      ipsilon,
      field,
      lmin = lat_lim[1],
      proj = proj,
      param = NULL,
      orient = c(90, 0, 0),
      npoints = 80
    )
    outfile <- list(
      x = p$x,
      y = p$y,
      z = p$z,
      xlim = p$xlim,
      ylim = p$ylim,
      xlab = "",
      ylab = "",
      axes = F
    )
  }
  return(outfile)
}

# function that provides labels and names for Blocking Plots
field_details <- function(field) {
  # default value
  legend_distance <- 3
  lev_hist <- NULL

  # case specific
  if (field == "TM90") {
    color_field <- c("dodgerblue", "darkred")
    color_diff <- NULL
    lev_field <- c(0, 30)
    lev_diff <- NULL
    legend_unit <- "Blocked Days (%)"
    title_name <- "TM90 Instantaneous Blocking"
  }

  if (field == "InstBlock") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 36, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    legend_unit <- "Blocked Days (%)"
    title_name <- "Instantaneous Blocking frequency:"
  }

  if (field == "ExtraBlock") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 36, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    legend_unit <- "Blocked Days (%)"
    title_name <-
      "Instantaneous Blocking frequency (GHGS2 condition):"
  }

  if (field == "BlockEvents") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 27, 3)
    lev_diff <- seq(-10.5, 10.5, 1)
    lev_hist <- c(0, 16)
    legend_unit <- "Blocked Days (%)"
    title_name <- "Blocking Events frequency:"
  }

  if (field == "LongBlockEvents") {
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 16, 2)
    lev_diff <- seq(-5.25, 5.25, .5)
    legend_unit <- "Blocked Days (%)"
    title_name <- "10-day Blocking Events frequency:"
  }

  if (field == "DurationEvents") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(5, 11.5, .5)
    lev_diff <- seq(-2.1, 2.1, .2)
    lev_hist <- c(6, 8)
    legend_unit <- "Duration (days)"
    title_name <- "Duration of Blocking Events:"
  }

  if (field == "NumberEvents") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(0, 100, 10)
    lev_diff <- seq(-42.5, 42.5, 5)
    lev_hist <- c(0, 60)
    legend_unit <- ""
    title_name <- "Number of Blocking Events:"
  }

  if (field == "Z500") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(4800, 6000, 50)
    lev_diff <- seq(-310, 310, 20)
    legend_unit <- "Geopotential Height (m)"
    title_name <- "Z500:"
    legend_distance <- 4
  }

  if (field == "BI") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(1, 6, 0.25)
    lev_diff <- seq(-2.1, 2.1, .2)
    legend_unit <- "BI index"
    title_name <- "Blocking Intensity (BI):"
  }

  if (field == "MGI") {
    color_field <- palette0
    color_diff <- palette2
    lev_field <- seq(0, 15, 1)
    lev_diff <- seq(-5.25, 5.25, .5)
    legend_unit <- "MGI Index"
    title_name <- "Meridional Gradient Inversion (MGI):"
  }

  if (field == "ACN" | field == "CN") {
    if (field == "ACN") {
      title_name <- "Anticyclonic Rossby wave breaking frequency:"
    }
    if (field == "CN") {
      title_name <- "Cyclonic Rossby wave breaking frequency:"
    }
    color_field <- palette1
    color_diff <- palette2
    lev_field <- seq(0, 20, 2)
    lev_diff <- seq(-5.25, 5.25, .5)
    legend_unit <- "RWB frequency (%)"
  }


  out <- list(
    color_field = color_field,
    color_diff = color_diff,
    lev_field = lev_field,
    lev_diff = lev_diff,
    lev_hist = lev_hist,
    legend_unit = legend_unit,
    legend_distance = legend_distance,
    title_name = title_name
  )
  return(out)
}

##########################################################
#------------Blocking Tracking Functions-----------------#
##########################################################

# time persistence (used for longitude filter too)
time_persistence <- function(timeseries, persistence = 5) {
  rr <- rle(timeseries)
  rr$values[which(rr$values == 1 & rr$length < persistence)] <- 0
  nn <- rep(rr$values, rr$length)
  return(nn)
}


# blocking 5 days tracking
blocking_persistence <-
  function(field, minduration = 5, time.array) {
    # function for persistence
    pers2 <- function(timeseries, persistence, time.array) {
      dd <- min(time.array$season):max(time.array$season)
      nn <- sapply(dd, function(x) {
        time_persistence(
          timeseries[which(time.array$season == x)],
          persistence
        )
      })
      xx <- c(unlist(nn))
      return(xx)
    }

    # check for etime
    if (length(time.array$month) != length(field[1, 1, ])) {
      stop("Wrong time array! Exiting...")
    }

    print("Time filtering...")
    newfield <- apply(field, c(1, 2), function(x)
      pers2(x,
        persistence = minduration, time.array
      ))
    newfield <- aperm(newfield, c(2, 3, 1))
    print("Mean field...")
    meanfield <- apply(newfield, c(1, 2), mean, na.rm = T) * 100


    print("Events detection...")
    maxdim <- max(apply(
      newfield, c(1, 2),
      function(x)
        length(rle(x)$length[which(rle(x)$values == 1)])
    ))
    events <- apply(
      newfield, c(1, 2),
      function(x)
        c(
          rle(x)$lengths[which(rle(x)$values == 1)],
          rep(NA, maxdim - length(rle(
            x
          )$length[which(rle(x)$values == 1)]))
        )
    )
    events <- aperm(events, c(2, 3, 1))
    print("Mean Duration...")
    duration <- apply(events, c(1, 2), mean, na.rm = T)
    print("Number of Events...")
    nevents <-
      apply(events, c(1, 2), function(x)
        length(x[!is.na(x)]))

    out <- list(
      track = newfield,
      percentage = meanfield,
      duration = duration,
      events = events,
      nevents = nevents
    )
    print(quantile(meanfield))
    print(min(duration, na.rm = T))
    return(out)
  }


# large scale extension with further implementation
largescale_extension_if <- function(ics, ipsilon, field) {
  print("Large Scale Extension based on fixed angle")
  fimin <- 30 # southern latitude to be analyzed
  fimax <- 75 # northern latitude to be analyzed
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]
  passo <- 5 / xreso # horizontal movemenent
  vertical <- 2.5 / yreso # vertical movement
  time <-
    which(apply(field, 3, max) != 0) # elements length of the dataset
  # (removing not blocked days)

  print(paste(
    "Box dimension:",
    passo * 2 * xreso,
    "° lon x ",
    vertical * 2 * yreso,
    "° lat"
  ))

  short <- function(ics, ipsilon, field, passo, vertical) {
    control <- field
    range <-
      which.min(abs(ipsilon - fimin)):which.min(abs(ipsilon - fimax))
    # check range for latitude excursion
    # reduce range considering border effect
    new <-
      rbind(field, field, field) # bind domain for cross-date line
    for (i in seq_along(ics)) {
      ii <- i + length(ics)
      # check to speed up
      if (!all(new[(ii - passo):(ii + passo), ] == 0)) {
        for (j in range) {
          control[i, j] <- mean(new[
            (ii - passo):(ii + passo),
            (j - vertical):(j + vertical)
          ], na.rm = T)
        }
      }
    }
    control[control > 0] <- 1
    return(control)
  }


  tt <- length(time)
  for (t in time) {
    progression_bar(t, tt)
    field[, , t] <-
      short(ics, ipsilon, field[, , t], passo, vertical)
  }
  return(field)
}


# Longitude filter for minimum extension
longitude_filter <- function(ics, ipsilon, field) {
  print("Longitude filter based on fixed angle")
  yreso <- ipsilon[2] - ipsilon[1]
  xreso <- ics[2] - ics[1]
  startipsilon <- which.min(abs(ipsilon - 30))
  estension <- (75 - 30) / yreso
  passo <- 15 / xreso

  print(paste("Continous longitude contrain", passo * xreso, "° lon"))

  tt <- length(field[1, 1, ])
  for (t in 1:tt) {
    progression_bar(t, tt)

    new <- rbind(field[, , t], field[, , t], field[, , t])
    for (j in startipsilon:((startipsilon + estension))) {
      new[, j] <- time_persistence(new[, j], persistence = passo)
    }
    field[, , t] <- new[length(ics) + (seq_along(ics)), ]
  }
  return(field)
}


##########################################################
#------------EOFs and regims functions-------------------#
##########################################################

eofs <-
  function(lon,
             lat,
             field,
             neof = 4,
             xlim,
             ylim,
             method = "SVD",
             do_standardize = F,
             do_regression = F) {
    # R tool for computing EOFs based on Singular Value Decomposition
    # ("SVD", default)
    # or with the eigenvectors of the covariance matrix ("covariance", slower)
    # If requested, computes linear regressions and standardizes the PCs
    # If you want to use the regressions, remember to standardize the PCs
    # Take as input a 3D anomaly field.
    # Requires "personal" functions area_weight, whicher and standardize

    # area weighting, based on the root of cosine
    print("Area Weighting...")
    ww <- area_weight(lon, lat, root = T)
    wwfield <- sweep(field, c(1, 2), ww, "*")

    # selection of the box
    box <- wwfield[
      whicher(lon, xlim[1]):whicher(lon, xlim[2]),
      whicher(lat, ylim[1]):whicher(lat, ylim[2]),
    ]
    slon <- lon[whicher(lon, xlim[1]):whicher(lon, xlim[2])]
    slat <- lat[whicher(lat, ylim[1]):whicher(lat, ylim[2])]

    # transform 3D field in a matrix
    new_box <-
      array(box, dim = c(dim(box)[1] * dim(box)[2], dim(box)[3]))

    # calling SVD
    if (method == "SVD") {
      print("Calling SVD...")
      SVD <- svd(new_box, nu = neof, nv = neof)

      # extracting EOFs (loading pattern),
      # expansions coefficient and variance explained
      pattern <- array(SVD$u, dim = c(dim(box)[1], dim(box)[2], neof))
      coefficient <- SVD$v
      variance <- (SVD$d[1:neof])^2 / sum((SVD$d)^2)
      if (do_standardize) {
        coefficient <- apply(coefficient, c(2), standardize)
      } else {
        coefficient <- sweep(coefficient, c(2), sqrt(variance), "*")
      }
    }

    # calling covariance matrix
    if (method == "covariance") {
      print("Calling eigenvectors of the covariance matrix...")
      covma <- cov(t(new_box))
      eig <- eigen(covma)
      coef <- (t(new_box) %*% eig$vector)[, 1:neof]
      pattern <- array(eig$vectors, dim = c(
        dim(box)[1], dim(box)[2],
        dim(box)[3]
      ))[, , 1:neof]
      variance <- eig$values[1:neof] / sum(eig$values)
      if (do_standardize) {
        coefficient <- apply(coef, c(2), standardize)
      } else {
        coefficient <- coef
      }
    }

    # linear regressions on anomalies
    regression <- NULL
    if (do_regression) {
      print("Linear Regressions (it can takes a while)... ")
      regression <- array(NA, dim = c(length(lon), length(lat), neof))
      for (i in 1:neof) {
        regression[, , i] <- apply(
          field, c(1, 2),
          function(x)
            lin.fit(as.matrix(coefficient[, i],
              ncol = 1
            ), x)$coefficients
        )
      }
    }

    # preparing output
    print("Finalize...")
    pattern <- list(x = slon, y = slat, z = pattern)
    out <- list(
      pattern = pattern,
      coeff = coefficient,
      variance = variance,
      regression = regression
    )
    return(out)
  }

eofs_coeff <-
  function(lon,
             lat,
             field,
             eof_object,
             do_standardize = F) {
    # Computes expansion coefficient (i.e. PCs) of a given dataset on the
    # loading pattern of EOF previously computed
    # Works only on eof_object obtained with "eofs" function

    # Area weighting, based on the root of cosine
    print("Area Weighting...")
    ww <- area_weight(lon, lat, root = T)
    wwfield <- sweep(field, c(1, 2), ww, "*")

    # selection of the box
    xlim <- c(min(eof_object$pattern$x), max(eof_object$pattern$x))
    ylim <- c(min(eof_object$pattern$y), max(eof_object$pattern$y))
    box <- wwfield[
      whicher(lon, xlim[1]):whicher(lon, xlim[2]),
      whicher(lat, ylim[1]):whicher(lat, ylim[2]),
    ]

    # transform 3D field in a matrix
    new_box <-
      array(box, dim = c(dim(box)[1] * dim(box)[2], dim(box)[3]))
    new_pattern <- array(eof_object$pattern$z,
      dim = c(
        dim(eof_object$pattern$z)[1] * dim(eof_object$pattern$z)[2],
        dim(eof_object$pattern$z)[3]
      )
    )

    # projects the coefficients
    coef <- (t(new_box) %*% new_pattern)

    # standardize
    if (do_standardize) {
      coefficient <- apply(coef, c(2), standardize)
    } else {
      coefficient <- coef
    }

    print("Finalize...")
    return(coefficient)
  }


regimes <- function(lon,
                    lat,
                    field,
                    ncluster = 4,
                    ntime = 1000,
                    neof = 10,
                    xlim,
                    ylim,
                    alg = "Hartigan-Wong") {
  # R tool to compute cluster analysis based on k-means.
  # Requires "personal" function eofs
  # Take as input a 3D anomaly field

  # Reduce the phase space with EOFs: use SVD and do not standardize PCs
  print("Launching EOFs...")
  reducedspace <- eofs(
    lon,
    lat,
    field,
    neof = neof,
    xlim = xlim,
    ylim = ylim,
    method = "SVD",
    do_regression = F,
    do_standardize = F
  )

  # extract the principal components
  PC <- reducedspace$coeff
  print(str(PC))

  # k-means computation repeat for ntime to find best solution.
  print("Computing k-means...")
  print(str(ncluster))
  regimes <- kmeans(
    PC,
    as.numeric(ncluster),
    nstart = ntime,
    iter.max = 1000,
    algorithm = alg
  )

  # Extract regimes frequencyr and timeseries of occupation
  cluster <- regimes$cluster
  frequencies <- regimes$size / dim(field)[3] * 100
  print(frequencies[order(frequencies, decreasing = T)])

  print("Creating Composites...")
  compose <-
    aperm(apply(field, c(1, 2), by, cluster, mean), c(2, 3, 1))

  # sorting from the more frequent to the less frequent
  kk <- order(frequencies, decreasing = T)
  cluster <- cluster + 10
  for (ss in 1:ncluster) {
    cluster[cluster == (ss + 10)] <- which(kk == ss)
  }

  # prepare output
  print("Finalize...")
  out <- list(
    cluster = cluster,
    frequencies = frequencies[kk],
    regimes = compose[, , kk],
    tot.withinss = regimes$tot.withinss
  )
  return(out)
}


regimes2 <-
  function(lon,
             lat,
             field,
             ncluster = 4,
             ntime = 1000,
             minvar = 0.8,
             xlim,
             ylim,
             alg = "Hartigan-Wong") {
    # R tool to compute cluster analysis based on k-means.
    # Requires "personal" function eofs (see above)
    # Take as input a 3D anomaly field

    # Reduce the phase space with EOFs: use SVD and do not standardize PCs
    print("Launching EOFs...")
    reducedspace <- eofs(
      lon,
      lat,
      field,
      neof = 20,
      xlim = xlim,
      ylim = ylim,
      method = "SVD",
      do_regression = F,
      do_standardize = F
    )
    reqpc <- which(cumsum(reducedspace$variance) > minvar)[1]
    print(
      paste(
        "Retaining",
        reqpc,
        "PCs to fullfil minimum explained variance required (",
        minvar * 100,
        "%)"
      )
    )

    # extract the principal components
    PC <- reducedspace$coeff[, 1:reqpc]
    print(str(PC))

    # k-means computation repeat for ntime to find best solution.
    print("Computing k-means...")
    print(str(ncluster))
    regimes <- kmeans(
      PC,
      as.numeric(ncluster),
      nstart = ntime,
      iter.max = 100,
      algorithm = alg
    )

    # Extract regimes frequencyr and timeseries of occupation
    cluster <- regimes$cluster
    frequencies <- regimes$size / dim(field)[3] * 100
    print(frequencies[order(frequencies, decreasing = T)])

    print("Creating Composites...")
    compose <-
      aperm(apply(field, c(1, 2), by, cluster, mean), c(2, 3, 1))

    # sorting from the more frequent to the less frequent
    kk <- order(frequencies, decreasing = T)
    cluster <- cluster + 10
    for (ss in 1:ncluster) {
      cluster[cluster == (ss + 10)] <- which(kk == ss)
    }

    # prepare output
    print("Finalize...")
    out <- list(
      cluster = cluster,
      frequencies = frequencies[kk],
      regimes = compose[, , kk],
      tot.withinss = regimes$tot.withinss
    )
    return(out)
  }

##########################################################
#-------------------Time Avg functions-------------------#
##########################################################

# fast function for monthly mean, using preallocation,
# vectorization and rowMeans
monthly_mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$month, etime$year)
  monthly <- array(NA, dim = c(
    length(ics), length(ipsilon),
    length(unique(condition))
  ))
  for (t in unique(condition)) {
    monthly[, , which(t == unique(condition))] <- rowMeans(field[
      , ,
      t == condition
    ], dims = 2)
  }
  return(monthly)
}

# introduce running mean
run_mean <- function(field, n = 5) {
  nn <- floor(n / 2)
  newfield <- field
  for (t in (1 + nn):(length(field) - nn)) {
    newfield[t] <- mean(field[(t - nn):(t + nn)])
  }
  return(newfield)
}

# improve running mean
# use vectorization for a 5 day running mean ad-hoc function
# (to be generalized!)
# about 10 times faster that a standard running mean function based
# on for loop
run_mean5 <- function(field) {
  newfield <- rowMeans(cbind(
    c(field[3:length(field)], NA, NA),
    c(field[2:length(field)], NA),
    field,
    c(NA, field[1:(length(field) - 1)]),
    c(NA, NA, field[1:(length(field) - 2)])
  ),
  na.rm = T
  )
  return(newfield)
}

# function for daily anomalies, use array predeclaration
# and rowMeans (40 times faster!)
daily_anom_mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(
    length(ics), length(ipsilon),
    length(unique(condition))
  ))
  anom <- field * NA
  for (t in unique(condition)) {
    daily[, , which(t == unique(condition))] <-
      rowMeans(field[, , t == condition], dims = 2)
    anom[, , which(t == condition)] <-
      sweep(
        field[, , which(t == condition)], c(1, 2),
        daily[, , which(t == unique(condition))], "-"
      )
  }
  return(anom)
}

# beta function for daily anomalies plus running mean
# (only 50% slower that standard daily avg)
daily_anom_run_mean <- function(ics, ipsilon, field, etime) {
  condition <- paste(etime$day, etime$month)
  daily <- array(NA, dim = c(
    length(ics), length(ipsilon),
    length(unique(condition))
  ))
  for (t in unique(condition)) {
    daily[, , which(t == unique(condition))] <-
      rowMeans(field[, , t == condition], dims = 2)
  }
  anom <- field * NA
  for (t in unique(condition)) {
    anom[, , which(t == condition)] <-
      sweep(
        field[, , which(t == condition)], c(1, 2),
        daily[, , which(t == unique(condition))], "-"
      )
  }
  return(anom)
}
