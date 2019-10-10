######################################################
#-----Hydroclimatic Intensity (HyInt) diagnostic-----#
#-------------E. Arnone (June 2017)------------------#
######################################################
hyint_diagnostic <- function(work_dir, # nolint
                             infile,
                             model_idx,
                             season,
                             prov_info,
                             rewrite = FALSE) {
  # setting up path and parameters
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]

  outfile <-
    getfilename_indices(work_dir, diag_base, model_idx, season)

  # If diagnostic output file already exists skip calculation
  if (file.exists(outfile) & !rewrite) {
    print(paste0(diag_base, ": output file already exists:", outfile))
    print(paste0(diag_base, ": skipping calculation"))
    return()
  }

  # Test input file exists
  print(infile)
  if (!file.exists(infile)) {
    stop("HyInt: missing regridded input file. Run HyInt pre-processing.")
  }

  # setting up time domain
  years <- year1:year2
  timeseason <- season2timeseason(season)

  # file opening
  pr_list <- ncdf_opener_universal(
    infile,
    namevar = "pr",
    tmonths = timeseason,
    tyears = years,
    rotate = rotlongitude
  )

  # extract calendar and time unit from the original file
  tcal <- attributes(pr_list$time)$cal
  tunit <- attributes(pr_list$time)$units

  etime <- power_date_new(pr_list$time)

  # declare and convert variable
  pr <- pr_list$field * 86400. # convert (Kg m-2 s-1) to (mm day-1)

  #############################################################
  #--------HyInt calculation (Giorgi et al. 2011/14)----------#
  #############################################################

  # Setup useful arrays and parameters
  nyear <- length(years)
  pry <-
    pr[, , 1:nyear] * NA # annual mean precipitation (over all days)
  int <-
    pr[, , 1:nyear] * NA # mean prec. intensity (over wet days == SDII)
  dsl <- pr[, , 1:nyear] * NA # mean dry spell length (DSL)
  wsl <- pr[, , 1:nyear] * NA # mean wet spell length (WSL)
  pa <- pr[, , 1:nyear] * NA # precipitation area (PA)
  r95 <- pr[, , 1:nyear] * NA # heavy precipitation index (R95)
  pry_norm <- pry
  int_norm <- int
  dsl_norm <- dsl
  wsl_norm <- wsl
  pa_norm <- pa
  r95_norm <- r95

  # Evaluate r95_threshold over normalization period (or load it if requested)
  if (external_r95[1] == F) {
    r95_threshold <-
      apply(pr, c(1, 2), quantile, probs = 0.95, na.rm = T)
  } else {
    # if required, use HyInt file from historical period
    if (external_r95[1] == "HIST") {
      external_r95 <- getfilename_indices(
        work_dir,
        diag_base,
        model_idx,
        season,
        hist = T,
        hist_years = norm_years
      )
    }
    r95_idx <-
      model_idx # assume each model has its r95_threshold file
    if (length(external_r95) == 1) {
      # if list of files with r95_threshold has only 1 entry,
      # use that for all models
      r95_idx <- 1
    }
    print(paste(
      diag_base,
      ": loading external r95_threshold data from ",
      external_r95[r95_idx]
    ))
    r95_threshold <-
      ncdf_opener(external_r95[r95_idx], "r95_threshold",
        rotate = "no"
      )
  }
  r95_threshold360 <- replicate(360, r95_threshold)
  r95_threshold365 <- replicate(365, r95_threshold)
  r95_threshold366 <- replicate(366, r95_threshold)

  # Calculate indices
  print(paste0(diag_base, ": calculating indices"))
  for (iyear in 1:nyear) {
    ret_year <- which(etime$year == years[iyear])
    pr_year <- pr[, , ret_year]

    r95_thresh_year <- r95_threshold365
    if (length(pr_year[1, 1, ]) == 360) {
      r95_thresh_year <- r95_threshold360
    }
    if (length(pr_year[1, 1, ]) == 366) {
      r95_thresh_year <- r95_threshold366
    }


    # Identify dry and wet days (Salinger and Griffiths 2001)
    ret_dry <- (pr_year < 1) # Dry days when pr < 1 mm
    ret_wet <- (pr_year >= 1) # Rainy days when pr >= 1 mm
    ret_below_r95 <-
      (pr_year < r95_thresh_year) # Rainy days when pr <
    # reference 95% quantile
    pr_year_dry <- pr_year * 0.
    pr_year_dry[ret_dry] <- 1 # mask with 1 for dry day
    pr_year_wet <- pr_year * 0.
    pr_year_wet[ret_wet] <- 1 # mask with 1 for rainy day
    pr_year_int <- pr_year
    pr_year_int[ret_dry] <-
      NA # actual precipitation but with NA on dry days
    pr_year_r95 <- pr_year
    pr_year_r95[ret_below_r95] <-
      NA # actual precipitation but with NA on
    # days with pr < reference 95% quantile

    # Mean annual precipitation
    pry_year <- apply(pr_year, c(1, 2), mean, na.rm = T)

    # Mean annual precipitation intensity (INT/SDII; intensity during wet days)
    int_year <- apply(pr_year_int, c(1, 2), mean, na.rm = T)

    # Mean annual dry spell length (DSL:
    # number of consecutive dry days during each dry spell).
    dsl_year <- mean_spell_length(pr_year_dry)

    # Mean annual wet  spell length (WSL:
    # number of consecutive wet days during each wet spell).
    wsl_year <- mean_spell_length(pr_year_wet)

    # Precipitation area (PA: number of rainy days * area of grid box)
    area_size <- area_size(ics, ipsilon)
    pa_year <-
      (apply(pr_year_wet, c(1, 2), sum, na.rm = T)) * area_size

    # Heavy precipitation Index (R95: percent of total precipitation above the
    # 95% percentile of the reference distribution);
    r95_year <- apply(pr_year_r95, c(1, 2), sum, na.rm = T) /
      apply(pr_year, c(1, 2), sum, na.rm = T) * 100.

    # Assign in-loop variables to storage array
    pry[, , iyear] <- pry_year
    dsl[, , iyear] <- dsl_year
    wsl[, , iyear] <- wsl_year
    int[, , iyear] <- int_year
    pa[, , iyear] <- pa_year
    r95[, , iyear] <- r95_year
  }

  # remove desert areas if required
  # (mean annual precipitation <0.5 mm, Giorgi et al. 2014)
  if (removedesert) {
    retdes <- which(pry < 0.5)
    pry[retdes] <- NA
    # create mask with NAs for deserts and 1's for not-desert
    retdes2D <- apply(pry * 0, c(1, 2), sum) + 1
    retdes3D <- replicate(nyear, retdes2D)
    pry <- pry * retdes3D
    dsl <- dsl * retdes3D
    wsl <- wsl * retdes3D
    int <- int * retdes3D
    pa <- pa * retdes3D
    r95 <- r95 * retdes3D
  }

  # Normalize to available data in reference period
  # NOTE: take care of normalization by 0: when the normalizing function is 0
  #       (e.g. short dataset in reference period), the resulting normalized
  #       index will be NA.

  # calculate normalization function
  if (external_norm[1] == F) {
    ret_years <- which(years >= norm_years[1] & years <= norm_years[2])
    if (length(ret_years) == 0) {
      stop(
        paste0(
          diag_base,
          ": no data over selected normalization period,
          unable to normalize"
        )
      )
    }
    pry_mean <- apply(pry[, , ret_years], c(1, 2), mean, na.rm = T)
    dsl_mean <- apply(dsl[, , ret_years], c(1, 2), mean, na.rm = T)
    wsl_mean <- apply(wsl[, , ret_years], c(1, 2), mean, na.rm = T)
    int_mean <- apply(int[, , ret_years], c(1, 2), mean, na.rm = T)
    pa_mean <- apply(pa[, , ret_years], c(1, 2), mean, na.rm = T)
    r95_mean <- apply(r95[, , ret_years], c(1, 2), mean, na.rm = T)
    pry_mean_sd <- apply(pry[, , ret_years], c(1, 2), sd, na.rm = T)
    dsl_mean_sd <- apply(dsl[, , ret_years], c(1, 2), sd, na.rm = T)
    wsl_mean_sd <- apply(wsl[, , ret_years], c(1, 2), sd, na.rm = T)
    int_mean_sd <- apply(int[, , ret_years], c(1, 2), sd, na.rm = T)
    pa_mean_sd <- apply(pa[, , ret_years], c(1, 2), sd, na.rm = T)
    r95_mean_sd <- apply(r95[, , ret_years], c(1, 2), sd, na.rm = T)
  } else {
    # load normalization data from file
    mean_idx <-
      model_idx # assume each model has its normalization file
    if (external_norm[1] == "HIST") {
      # if required, use HyInt file from historical period
      external_norm <-
        getfilename_indices(
          work_dir,
          diag_base,
          model_idx,
          season,
          hist = T,
          hist_years = norm_years
        )
    }
    if (length(external_norm) == 1) {
      mean_idx <- 1
    }
    # if list of files with normalization functions has only 1 entry,
    # use that for all models
    print(paste(
      diag_base,
      ": loading external normalization data from ",
      external_norm[mean_idx]
    ))
    pry_mean <-
      ncdf_opener(external_norm[mean_idx], "pry_mean", rotate = "no")
    dsl_mean <-
      ncdf_opener(external_norm[mean_idx], "dsl_mean", rotate = "no")
    wsl_mean <-
      ncdf_opener(external_norm[mean_idx], "wsl_mean", rotate = "no")
    int_mean <-
      ncdf_opener(external_norm[mean_idx], "int_mean", rotate = "no")
    pa_mean <-
      ncdf_opener(external_norm[mean_idx], "pa_mean", rotate = "no")
    r95_mean <-
      ncdf_opener(external_norm[mean_idx], "r95_mean", rotate = "no")
    pry_mean_sd <-
      ncdf_opener(external_norm[mean_idx], "pry_mean_sd",
        rotate = "no"
      )
    dsl_mean_sd <-
      ncdf_opener(external_norm[mean_idx], "dsl_mean_sd",
        rotate = "no"
      )
    wsl_mean_sd <-
      ncdf_opener(external_norm[mean_idx], "wsl_mean_sd",
        rotate = "no"
      )
    int_mean_sd <-
      ncdf_opener(external_norm[mean_idx], "int_mean_sd",
        rotate = "no"
      )
    pa_mean_sd <- ncdf_opener(external_norm[mean_idx], "pa_mean_sd",
      rotate = "no"
    )
    r95_mean_sd <-
      ncdf_opener(external_norm[mean_idx], "r95_mean_sd",
        rotate = "no"
      )
  }

  # remove 0s from normalizing functions
  pry_mean[pry_mean == 0] <- NA
  dsl_mean[dsl_mean == 0] <- NA
  wsl_mean[wsl_mean == 0] <- NA
  int_mean[int_mean == 0] <- NA
  pa_mean[pa_mean == 0] <- NA
  r95_mean[r95_mean == 0] <- NA

  # perform normalization
  for (iyear in 1:nyear) {
    pry_norm[, , iyear] <- pry[, , iyear] / pry_mean
    dsl_norm[, , iyear] <- dsl[, , iyear] / dsl_mean
    wsl_norm[, , iyear] <- wsl[, , iyear] / wsl_mean
    int_norm[, , iyear] <- int[, , iyear] / int_mean
    pa_norm[, , iyear] <- pa[, , iyear] / pa_mean
    r95_norm[, , iyear] <- r95[, , iyear] / r95_mean
  }

  # Calculate HY-INT index
  hyint <- dsl_norm * int_norm


  # Calculate mean and mean_sd for hyint
  if (external_norm[1] == F) {
    # calculate or load hyint_mean from file for consistency with other indice
    ret_years <-
      which(years >= norm_years[1] & years <= norm_years[2])
    hyint_mean <-
      apply(hyint[, , ret_years], c(1, 2), mean, na.rm = T)
    hyint_mean_sd <-
      apply(hyint[, , ret_years], c(1, 2), sd, na.rm = T)
  } else {
    # load normalization data from file
    mean_idx <-
      model_idx # assume each model has its normalization file
    if (length(external_norm) == 1) {
      # if list of files with normalization functions has only 1 entry,
      # use that for all models
      mean_idx <- 1
    }
    hyint_mean <- ncdf_opener(external_norm[mean_idx], "hyint_mean",
      rotate = "no"
    )
    hyint_mean_sd <-
      ncdf_opener(external_norm[mean_idx], "hyint_mean_sd",
        rotate = "no"
      )
  }

  # HyInt list
  hyint_list <- list(
    pry = pry,
    dsl = dsl,
    wsl = wsl,
    int = int,
    pa = pa,
    r95 = r95,
    hyint = hyint,
    pry_mean = pry_mean,
    dsl_mean = dsl_mean,
    wsl_mean = wsl_mean,
    int_mean = int_mean,
    pa_mean = pa_mean,
    r95_mean = r95_mean,
    hyint_mean = hyint_mean,
    pry_mean_sd = pry_mean_sd,
    dsl_mean_sd = dsl_mean_sd,
    wsl_mean_sd = wsl_mean_sd,
    int_mean_sd = int_mean_sd,
    pa_mean_sd = pa_mean_sd,
    r95_mean_sd = r95_mean_sd,
    hyint_mean_sd = hyint_mean_sd,
    pry_norm = pry_norm,
    dsl_norm = dsl_norm,
    wsl_norm = wsl_norm,
    int_norm = int_norm,
    pa_norm = pa_norm,
    r95_norm = r95_norm,
    r95_threshold = r95_threshold
  )

  print(
    paste(
      diag_base,
      ": calculation done. Returning mean precipitation,
      sdii, dsl, wsl, pa, r95 (absolute and normalized values)
      and hyint indices"
    )
  )


  ##########################################################
  #------------------------Save to NetCDF------------------#
  ##########################################################

  # saving output to netcdf files
  print(paste0(diag_base, ": saving data to NetCDF file:"))

  # define fieds to be saved
  field_list <- c(
    "pry",
    "dsl",
    "wsl",
    "int",
    "pa",
    "r95",
    "hyint",
    "pry_mean",
    "dsl_mean",
    "wsl_mean",
    "int_mean",
    "pa_mean",
    "r95_mean",
    "hyint_mean",
    "pry_mean_sd",
    "dsl_mean_sd",
    "wsl_mean_sd",
    "int_mean_sd",
    "pa_mean_sd",
    "r95_mean_sd",
    "hyint_mean_sd",
    "pry_norm",
    "dsl_norm",
    "wsl_norm",
    "int_norm",
    "pa_norm",
    "r95_norm",
    "r95_threshold"
  )

  TIME <-
    paste(tunit, " since ", year1, "-", timeseason[1], "-01 00:00:00",
      sep = ""
    )

  # dimensions definition
  x <- ncdim_def("lon", "degrees_east", ics, longname = "longitude")
  y <-
    ncdim_def("lat", "degrees_north", ipsilon, longname = "latitude")
  t <- ncdim_def(
    timedimname,
    "years",
    years,
    unlim = T,
    calendar = tcal,
    longname = timedimname
  )
  # timedim <- ncdim_def( timedimname,
  #            "years since 1950-01-01 00:00:00", (years-1950),unlim=T)

  # t <- ncdim_def( timedimname, TIME, years,
  #      unlim=T, calendar=tcal, longname=timedimname)

  for (var in field_list) {
    field <- get(var, hyint_list)
    field[is.nan(field)] <- NA
    metadata <- setmetadata_indices(var)
    longvar <- metadata$longvar
    unit <- metadata$unit
    # variable definitions
    var_ncdf <- ncvar_def(
      var,
      unit,
      list(x, y, t),
      -999,
      longname = longvar,
      prec = "single",
      compression = 1
    )
    if ((var == "pry_mean") |
      (var == "int_mean") | (var == "dsl_mean") |
      (var == "wsl_mean") |
      (var == "pa_mean") | (var == "r95_mean") |
      (var == "hyint_mean") | (var == "pry_mean_sd") |
      (var == "int_mean_sd") | (var == "dsl_mean_sd") |
      (var == "wsl_mean_sd") | (var == "pa_mean_sd") |
      (var == "r95_mean_sd") | (var == "hyint_mean_sd") |
      (var == "r95_threshold")) {
      var_ncdf <- ncvar_def(
        var,
        unit,
        list(x, y),
        -999,
        longname = longvar,
        prec = "single",
        compression = 1
      )
    }
    assign(paste0("var", var), var_ncdf)
    assign(paste0("field", var), field)
  }

  # Netcdf file creation
  print(paste(diag_base, ": saving output to ", outfile))
  namelist <- paste0("var", field_list)
  nclist <- mget(namelist)
  ncfile <- nc_create(outfile, nclist)
  for (var in field_list) {
    # put variables into the ncdf file
    ndims <- get(paste0("var", var))$ndims
    ncvar_put(
      ncfile,
      var,
      get(paste0("field", var)),
      start = rep(1, ndims),
      count = rep(-1, ndims)
    )
  }
  nc_close(ncfile)

  # Set provenance for this output file
  caption <-
    paste0(
      "Hyint indices  for years ",
      year1,
      " to ",
      year2,
      " according to ",
      models_name[model_idx]
    )
  xprov <- list(
    ancestors = list(infile),
    model_idx = list(model_idx),
    caption = caption
  )

  # Store provenance in main provenance list
  prov_info[[outfile]] <- xprov

  print(paste(diag_base, ": diagnostic netCDF files saved"))
  return(prov_info)
}
