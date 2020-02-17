# #############################################################################
# quantilebias_functions.R
#
# Author: Enrico Arnone (ISAC-CNR, Italy)
#
# #############################################################################
# Description
#    Originally developed as functions to be used in HyInt routines
#
# Modification history
#    20170901-A_arno_en: adapted to HyInt and extended
#    20170522-A_davi_pa: Creation for MiLES
# #############################################################################

# basis functions

##########################################################
#------------------------Packages------------------------#
##########################################################

# loading packages
library("maps")
library("ncdf4")
library("PCICt")

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
      unlist(lapply(a$var[[namevar]]$dim, function(x) {
        x["name"]
      }))
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

##########################################################
#--------------Plotting functions------------------------#
##########################################################

graphics_startup <- function(figname, output_file_type, plot_size) {
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
                         col = col,
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
  if (missing(col)) {
    col <- color.palette(length(levels) - 1)
  }
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
