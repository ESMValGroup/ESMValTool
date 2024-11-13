######################################################
#-----Blocking routines computation for MiLES--------#
#-------------P. Davini (Oct 2014)-------------------#
######################################################
miles_block_fast <- # nolint
  function(dataset,
             expid,
             ens,
             year1,
             year2,
             season,
             z500filename,
             FILESDIR,
             doforce) {
    t0 <- proc.time()

    # setting up time domain
    years <- year1:year2
    timeseason <- season2timeseason(season)

    # define folders using file.builder function (takes care of ensembles)
    savefile1 <- file_builder(
      FILESDIR,
      "Block",
      "BlockClim",
      dataset,
      expid,
      ens,
      year1,
      year2,
      season
    )
    savefile2 <- file_builder(
      FILESDIR,
      "Block",
      "BlockFull",
      dataset,
      expid,
      ens,
      year1,
      year2,
      season
    )

    # check if data is already there to avoid re-run
    if (file.exists(savefile1) & file.exists(savefile2)) {
      print("Actually requested blocking data is already there!")
      print(savefile1)
      print(savefile2)
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
      namevar = "zg",
      tmonths = timeseason,
      tyears = years,
      rotate = "full",
      fillmiss = TRUE
    )
    print(str(fieldlist))

    # extract calendar and time unit from the original file
    tcal <- attributes(fieldlist$time)$cal
    tunit <- attributes(fieldlist$time)$units

    # time array to simplify time filtering
    etime <- power_date_new(fieldlist$time)
    totdays <- length(fieldlist$time)

    # declare variable
    z500 <- fieldlist$field

    # grid resolution
    yreso <- ipsilon[2] - ipsilon[1]
    xreso <- ics[2] - ics[1]

    # reso checks: this are not needed with default 2.5 grid,
    # but they may be relevant with
    # future envisaged power up to finer grids
    # xcritical factor is due to RWB longitudinal jump of 7.5
    # ycritical factor is due to Large Scale Extension of 2.5
    xcritical <- 2.5
    ycritical <- 2.5
    if (ycritical %% yreso != 0) {
      stop("Latitudinal resolution is not a factor of 5 deg")
    }

    if (xcritical %% xreso != 0) {
      stop("Longitudinal resolution is not a factor of 5 deg")
    }

    ##########################################################
    #--------------Tibaldi and Molteni 1990------------------#
    ##########################################################

    print("Tibaldi and Molteni (1990) index...")
    # TM90: parametres for blocking detection
    tm90_fi0 <- 60 # central_lat
    tm90_fin <- tm90_fi0 + 20
    tm90_fis <- tm90_fi0 - 20 # south and north lat, 80N and 40N
    tm90_central <- whicher(ipsilon, tm90_fi0)
    tm90_south <- whicher(ipsilon, tm90_fis)
    tm90_north <- whicher(ipsilon, tm90_fin)
    tm90_range <- seq(-5, 5, yreso) / yreso # 5 degrees to the north,
    #  5 to the south (larger than TM90 or D'Andrea et al 1998)

    # TM90: beta version, the amazing power of R vectorization!
    # 6 lines to get the climatology
    tm90_ghgn <- (z500[, tm90_north + tm90_range, ] -
      z500[, tm90_central + tm90_range, ]) / (tm90_fin - tm90_fi0)
    tm90_ghgs <- (z500[, tm90_central + tm90_range, ] -
      z500[, tm90_south + tm90_range, ]) / (tm90_fi0 - tm90_fis)
    tm90_check <-
      (tm90_ghgs > 0 & tm90_ghgn < (-10)) # TM90 conditions
    tm90_check[tm90_check == T] <- 1
    tm90_check[tm90_check == F] <- 0
    tottm90 <- apply(tm90_check, c(1, 3), max, na.rm = T)
    tm90 <- apply(tottm90, 1, mean) * 100
    print("Done!")

    ##########################################################
    #--------------Davini et al. 2012------------------------#
    ##########################################################

    # decleare main variables to be computed (considerable speed up!)
    totrwb <- totmeridional <- totbi <- z500 * NA
    totblocked <- totblocked2 <- z500 * 0

    # Davini et al. 2012: parameters to be set for blocking detection
    fi0 <- 30 # lowest latitude to be analyzed
    jump <- 15 # distance on which compute gradients
    step0 <- jump / yreso # number of grid points to be used
    central <-
      which.min(abs(ipsilon - fi0)) # lowest starting latitude
    north <- central + step0 # lowest north latitude
    south <- central - step0 # lowest sourth latitude
    maxsouth <- central - 2 * step0
    if (maxsouth <= 0) {
      stop("Latitude=0 has to be included in the domain.")
    }
    fin <- ipsilon[north]
    fis <- ipsilon[south]
    range <- (90 - fi0 - jump) / yreso # escursion to the north for
    # computing blocking (from 30 up to 75)

    print("--------------------------------------------------")
    print("Davini et al. (2012) index and diagnostics...")
    print(c("distance for gradients:", step0 * diff(ics)[1]))
    print(paste("range of latitudes ", fi0, "-", 90 - step0 * diff(ics)[1],
      " N",
      sep = ""
    ))

    ##########################################################
    #--------------Istantaneous Blocking---------------------#
    ##########################################################

    #----COMPUTING BLOCKING INDICES-----
    for (t in 1:totdays) {
      progression_bar(t, totdays)

      # multidim extension
      new_field <- rbind(z500[, , t], z500[, , t], z500[, , t])

      # computing blocking for different latitudes
      for (delta in 0:range) {
        ghgn <- (z500[, north + delta, t] -
          z500[, central + delta, t]) / (fin - fi0)
        ghgs <- (z500[, central + delta, t] -
          z500[, south + delta, t]) / (fi0 - fis)
        gh2gs <- (z500[, south + delta, t] -
          z500[, maxsouth + delta, t]) / (fi0 - fis)
        check1 <- which(ghgs > 0 & ghgn < (-10))
        check2 <- which(ghgs > 0 & ghgn < (-10) & gh2gs < (-5))
        # supplementary condition

        if (length(check2) > 0) {
          totblocked2[check2, central + delta, t] <- 1
        }

        if (length(check1) > 0) {
          # 1-MATRIX FOR INSTANTANEOUS BLOCKING
          totblocked[check1, central + delta, t] <- 1


          # 2-PART ON COMPUTATION OF ROSSBY WAVEBREAKING
          r <- check1 + length(ics)
          rwb_jump <- jump / 2
          steprwb <- rwb_jump / xreso
          rwb_west <-
            new_field[(r - steprwb), south + delta + steprwb]
          rwb_east <-
            new_field[(r + steprwb), south + delta + steprwb]
          fullgh <- (rwb_west - rwb_east)

          totrwb[check1[fullgh < 0], central + delta, t] <- (-10)
          # gradient decreasing: cyclonic RWB
          totrwb[check1[fullgh > 0], central + delta, t] <- 10
          # gradient increasing: anticyclonic RWB

          # 4-part about adapted version of blocking intensity
          # by Wiedenmann et al. (2002)
          step <- 60 / xreso
          ii <- check1 + length(ics)
          zu <- zd <- NULL
          for (ll in ii) {
            zu <- c(zu, min(new_field[(ll - step):ll, central + delta]))
            zd <-
              c(zd, min(new_field[ll:(ll + step), central + delta]))
          }
          mz <- z500[check1, central + delta, t]
          rc <- 0.5 * ((zu + mz) / 2 + (zd + mz) / 2)
          totbi[check1, central + delta, t] <- 100 * (mz / rc - 1)

          # 5 - part about meridional gradient index
          totmeridional[check1, central + delta, t] <- ghgs[check1]
        }
      }
    }

    print(paste("Total # of days:", t))
    print("-------------------------")

    ##########################################################
    #--------------------Mean Values-------------------------#
    ##########################################################

    # compute mean values (use rowMeans, faster when there are no NA values)
    frequency <- rowMeans(totblocked, dims = 2) * 100
    # frequency of Instantaneous Blocking days
    frequency2 <- rowMeans(totblocked2, dims = 2) * 100
    # frequency of Instantaneous Blocking days with GHGS2
    z500mean <- rowMeans(z500, dims = 2) # Z500 mean value
    bi <- apply(totbi, c(1, 2), mean, na.rm = T)
    # Blocking Intensity Index as Wiedenmann et al. (2002)
    mgi <- apply(totmeridional, c(1, 2), mean, na.rm = T)
    # Value of meridional gradient inversion

    # anticyclonic and cyclonic averages RWB
    cn <-
      apply(totrwb, c(1, 2), function(x)
        sum(x[x == (-10)], na.rm = T)) /
        (totdays) * (-10)
    acn <-
      apply(totrwb, c(1, 2), function(x)
        sum(x[x == (10)], na.rm = T)) /
        (totdays) * (10)

    t1 <- proc.time() - t0
    print(t1)

    print("Instantaneous blocking and diagnostics done!")

    ##########################################################
    #--------------------Time filtering----------------------#
    ##########################################################

    # spatial filtering on fixed longitude distance
    spatial <- longitude_filter(ics, ipsilon, totblocked)

    # large scale extension on 10x5 box
    large <- largescale_extension_if(ics, ipsilon, spatial)

    # 5-day persistence filter
    block <-
      blocking_persistence(large, minduration = 5, time.array = etime)

    # 10-day persistence for extreme long block
    longblock <- blocking_persistence(large,
      minduration = 10,
      time.array = etime
    )

    tf <- proc.time() - t1
    print(tf)


    ##########################################################
    #------------------------Save to NetCDF------------------#
    ##########################################################

    # saving output to netcdf files
    print("saving NetCDF climatologies...")

    # which fieds to plot/save
    fieldlist <- c(
      "TM90",
      "InstBlock",
      "ExtraBlock",
      "Z500",
      "MGI",
      "BI",
      "CN",
      "ACN",
      "BlockEvents",
      "LongBlockEvents",
      "DurationEvents",
      "NumberEvents"
    )
    full_fieldlist <- c(
      "TM90",
      "InstBlock",
      "ExtraBlock",
      "Z500",
      "MGI",
      "BI",
      "CN",
      "ACN",
      "BlockEvents",
      "LongBlockEvents"
    )

    # dimensions definition
    fulltime <- as.numeric(etime$data) - as.numeric(etime$data)[1]
    TIME <- paste(tunit, " since ", year1, "-", timeseason[1],
      "-01 00:00:00",
      sep = ""
    )
    LEVEL <- 50000
    x <- ncdim_def("lon", "degrees_east", ics, longname = "longitude")
    y <-
      ncdim_def("lat", "degrees_north", ipsilon, longname = "latitude")
    z <- ncdim_def("plev", "Pa", LEVEL, longname = "pressure")
    t1 <- ncdim_def(
      "time",
      TIME,
      0,
      unlim = T,
      calendar = tcal,
      longname = "time"
    )
    t2 <- ncdim_def(
      "time",
      TIME,
      fulltime,
      unlim = T,
      calendar = tcal,
      longname = "time"
    )

    for (var in fieldlist) {
      # name of the var
      if (var == "TM90") {
        longvar <- "Tibaldi-Molteni 1990 Instantaneous Blocking frequency"
        unit <- "%"
        field <- tm90
        full_field <- tottm90
      }
      if (var == "InstBlock") {
        longvar <- "Instantaneous Blocking frequency"
        unit <- "%"
        field <- frequency
        full_field <- totblocked
      }
      if (var == "ExtraBlock") {
        longvar <- "Instantaneous Blocking frequency (GHGS2)"
        unit <- "%"
        field <- frequency2
        full_field <- totblocked2
      }
      if (var == "Z500") {
        longvar <- "Geopotential Height"
        unit <- "m"
        field <- z500mean
        full_field <- z500
      }
      if (var == "BI") {
        longvar <- "BI index"
        unit <- ""
        field <- bi
        full_field <- totbi
      }
      if (var == "MGI") {
        longvar <- "MGI index"
        unit <- ""
        field <- mgi
        full_field <- totmeridional
      }
      if (var == "ACN") {
        longvar <- "Anticyclonic RWB frequency"
        unit <- "%"
        field <- acn
        full_field <- totrwb / 10
        full_field[full_field == (-1)] <- NA
      }
      if (var == "CN") {
        longvar <- "Cyclonic RWB frequency"
        unit <- "%"
        field <- cn
        full_field <- totrwb / 10
        full_field[full_field == (1)] <- NA
      }
      if (var == "BlockEvents") {
        longvar <- "Blocking Events frequency"
        unit <- "%"
        field <- block$percentage
        full_field <- block$track
      }
      if (var == "LongBlockEvents") {
        longvar <- "10-day Blocking Events frequency"
        unit <- "%"
        field <- longblock$percentage
        full_field <- longblock$track
      }
      if (var == "DurationEvents") {
        longvar <- "Blocking Events duration"
        unit <- "days"
        field <- block$duration
      }
      if (var == "NumberEvents") {
        longvar <- "Blocking Events number"
        unit <- ""
        field <- block$nevents
      }

      # fix eventual NaN
      field[is.nan(field)] <- NA

      # variable definitions
      if (var == "TM90") {
        var_ncdf <- ncvar_def(
          var,
          unit,
          list(x, t = t1),
          -999,
          longname = longvar,
          prec = "single",
          compression = 1
        )
        full_var_ncdf <- ncvar_def(
          var,
          unit,
          list(x, t = t2),
          -999,
          longname = longvar,
          prec = "single",
          compression = 1
        )
      } else {
        var_ncdf <- ncvar_def(
          var,
          unit,
          list(x, y, z, t = t1),
          -999,
          longname = longvar,
          prec = "single",
          compression = 1
        )
        full_var_ncdf <-
          ncvar_def(
            var,
            unit,
            list(x, y, z, t = t2),
            -999,
            longname = longvar,
            prec = "single",
            compression = 1
          )
      }

      assign(paste0("var", var), var_ncdf)
      assign(paste0("full_var", var), full_var_ncdf)
      assign(paste0("field", var), field)
      assign(paste0("full_field", var), full_field)
    }

    # Climatologies Netcdf file creation
    print(savefile1)
    namelist1 <- paste0("var", fieldlist)
    nclist1 <- mget(namelist1)
    ncfile1 <- nc_create(savefile1, nclist1)
    for (var in fieldlist) {
      # put variables into the ncdf file
      ndims <- get(paste0("var", var))$ndims
      ncvar_put(
        ncfile1,
        var,
        get(paste0("field", var)),
        start = rep(1, ndims),
        count = rep(-1, ndims)
      )
    }
    nc_close(ncfile1)

    # Fullfield Netcdf file creation
    print(savefile2)
    namelist2 <- paste0("full_var", full_fieldlist)
    nclist2 <- mget(namelist2)
    ncfile2 <- nc_create(savefile2, nclist2)
    for (var in full_fieldlist) {
      # put variables into the ncdf file
      ndims <- get(paste0("full_var", var))$ndims
      ncvar_put(
        ncfile2,
        var,
        get(paste0("full_field", var)),
        start = rep(1, ndims),
        count = rep(-1, ndims)
      )
    }
    nc_close(ncfile2)
    return(c(savefile1, savefile2))
  }
