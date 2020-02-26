######################################################
#-----EOFs routines computation for MiLES--------#
#-------------P. Davini (Feb 2018)-------------------#
######################################################
miles_eofs_fast <- # nolint
  function(dataset,
             expid,
             ens,
             year1,
             year2,
             season,
             tele,
             z500filename,
             FILESDIR,
             PROGDIR,
             doforce) {
    # standard defined 4 EOFs
    neofs <- 4

    # t0
    t0 <- proc.time()

    # setting up time domain
    years <- year1:year2
    timeseason <- season2timeseason(season)

    # define folders using file.builder function (takes care of ensembles)
    print(".....")
    print(dataset)
    print(expid)
    print(ens)
    savefile1 <- file_builder(
      FILESDIR,
      paste0("EOFs/", tele),
      "EOFs",
      dataset,
      expid,
      ens,
      year1,
      year2,
      season
    )

    # select teleconnection region
    if (tele == "NAO") {
      xlim <- c(-90, 40)
      ylim <- c(20, 85)
      rotation <- "full"
    } else if (tele == "AO") {
      xlim <- c(-180, 180)
      ylim <- c(20, 85)
      rotation <- "full"
    } else if (tele == "PNA") {
      xlim <- c(140, 280)
      ylim <- c(20, 85)
      rotation <-
        "no" # 140E-80W: use trick of rotation for cross-dateline
    } else {
      # use non standard region, detect region with strsplit
      splitter <- as.numeric(strsplit(tele, "_")[[1]])
      if (length(splitter) == 4) {
        xlim <- c(splitter[1], splitter[2])
        ylim <- c(splitter[3], splitter[4])
        if (xlim[2] > 180) {
          rotation <- "no"
        } else {
          rotation <- "full"
        }
      } else {
        stop("Wrong teleconnection region!")
      }
    }

    # check if data is already there to avoid re-run
    if (file.exists(savefile1)) {
      print("Actually requested EOFs data is already there!")
      print(savefile1)
      if (doforce == TRUE) {
        print("Running with doforce=true... re-run!")
      } else {
        print("Skipping... activate doforce=true if you want to re-run it")
        q()
      }
    }

    # new file opening
    nomefile <- z500filename
    fieldlist <- ncdf_opener_universal(
      nomefile,
      "zg",
      tmonths = timeseason,
      tyears = years,
      rotate = rotation
    )
    print(str(fieldlist))

    # extract calendar and time unit from the original file
    tcal <- attributes(fieldlist$time)$cal
    tunit <- attributes(fieldlist$time)$units

    # time array
    etime <- power_date_new(fieldlist$time)

    # declare variable
    z500 <- fieldlist$field

    # monthly averaging
    print("monthly mean...")

    # new faster monthly mean function
    z500monthly <- monthly_mean(ics, ipsilon, z500, etime)

    # climatology
    print("climatological mean...")
    z500clim <-
      apply(z500monthly, c(1, 2), ave, rep(timeseason, length(years)))
    z500clim <- aperm(z500clim, c(2, 3, 1))

    # monthly anomalies
    print("anomalies...")
    z500anom <- z500monthly - z500clim

    # compute EOFs
    print("EOFs...")
    EOFS <- eofs(
      ics,
      ipsilon,
      z500anom,
      neof = neofs,
      xlim,
      ylim,
      method = "SVD",
      do_standardize = T,
      do_regression = T
    )
    # COEFF=eofs.coeff(ics,ipsilon,z500anom,EOFS,
    # do_standardize=T) #do we really need this?

    # flip signs of patterns and regressions for NAO and AO
    print("checking signs...")
    for (i in 1:neofs) {
      posreg <- NULL

      # define regions for sign control: boxes where values should be positive
      if (tele == "NAO") {
        if (i == 1) {
          posreg <- c(-30, 30, 40, 50)
        } # NAO
        if (i == 2) {
          posreg <- c(-60, 0, 40, 60)
        } # East Atlantic Pattern
        if (i == 3) {
          posreg <- c(-30, 30, 50, 70)
        } # Scandinavian Blocking
      }

      if (tele == "AO") {
        if (i == 1) {
          posreg <- c(-180, 180, 20, 50)
        } # Arctic Oscillation
        if (i == 2) {
          posreg <- c(-120, -60, 40, 60)
        } # PNA
      }

      # if definition of region exists
      if (!is.null(posreg)) {
        # convert into indices
        xbox <- whicher(EOFS$pattern$x, posreg[1]):whicher(
          EOFS$pattern$x,
          posreg[2]
        )
        ybox <- whicher(EOFS$pattern$y, posreg[3]):whicher(
          EOFS$pattern$y,
          posreg[4]
        )
        valuereg <- mean(EOFS$pattern$z[xbox, ybox, i])

        # if negative in the box, flip all signs!
        if (valuereg < 0) {
          EOFS$pattern$z[, , i] <- -EOFS$pattern$z[, , i]
          EOFS$regression <- -EOFS$regression
        }
      }
    }

    # expand EOF pattern to save it
    expanded_pattern <- EOFS$regression * NA
    expanded_pattern[
      whicher(ics, xlim[1]):whicher(ics, xlim[2]),
      whicher(ipsilon, ylim[1]):whicher(ipsilon, ylim[2]),
    ] <-
      EOFS$pattern$z

    t1 <- proc.time() - t0
    print(t1)


    ##########################################################
    #------------------------Save to NetCDF------------------#
    ##########################################################

    # saving output to netcdf files
    print("saving NetCDF climatologies...")
    print(savefile1)

    # monthly specific time
    monthtime <- as.numeric(etime$data[etime$day == 15])

    # dimensions definition
    TIME <- paste(tunit, " since ", year1, "-", timeseason[1],
      "-01 00:00:00",
      sep = ""
    )
    LEVEL <- 50000
    x <- ncdim_def("lon", "degrees_east", ics, longname = "longitude")
    y <-
      ncdim_def("lat", "degrees_north", ipsilon, longname = "latitude")
    z <- ncdim_def("plev", "Pa", LEVEL, longname = "pressure")
    ef <- ncdim_def("PC", "-", 1:neofs)
    t <- ncdim_def(
      "time",
      TIME,
      monthtime,
      calendar = tcal,
      longname = "time",
      unlim = T
    )

    # defining vars
    unit <- "m"
    longvar <- "EOFs Loading Pattern"
    pattern_ncdf <-
      ncvar_def(
        "Patterns",
        unit,
        list(x, y, z, ef),
        -999,
        longname = longvar,
        prec = "single",
        compression = 1
      )

    unit <- "m"
    longvar <- "EOFs Linear Regressions"
    regression_ncdf <-
      ncvar_def(
        "Regressions",
        unit,
        list(x, y, z, ef),
        -999,
        longname = longvar,
        prec = "single",
        compression = 1
      )

    unit <- paste0("0-", neofs)
    longvar <- "PCs timeseries"
    pc_ncdf <- ncvar_def(
      "PCs",
      unit,
      list(ef, t),
      -999,
      longname = longvar,
      prec = "single",
      compression = 1
    )

    unit <- "%"
    longvar <- "EOFs variance"
    variance_ncdf <- ncvar_def(
      "Variances",
      unit,
      list(ef),
      -999,
      longname = longvar,
      prec = "single",
      compression = 1
    )

    # saving files
    ncfile1 <- nc_create(
      savefile1,
      list(pattern_ncdf, pc_ncdf, variance_ncdf, regression_ncdf)
    )
    ncvar_put(
      ncfile1,
      "Patterns",
      expanded_pattern,
      start = c(1, 1, 1, 1),
      count = c(-1, -1, -1, -1)
    )
    ncvar_put(
      ncfile1,
      "Regressions",
      EOFS$regression,
      start = c(1, 1, 1, 1),
      count = c(-1, -1, -1, -1)
    )
    ncvar_put(ncfile1,
      "PCs",
      EOFS$coeff,
      start = c(1, 1),
      count = c(-1, -1)
    )
    ncvar_put(ncfile1,
      "Variances",
      EOFS$variance,
      start = c(1),
      count = c(-1)
    )
    nc_close(ncfile1)
    return(savefile1)
  }
