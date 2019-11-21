######################################################
#-------------Trends routine for HyInt---------------#
#-------------E. Arnone (June 2017)------------------#
######################################################


# MAIN TRENDS FUNCTION
hyint_trends <- function(work_dir, model_idx, season, prov_info) { # nolint
  # setup useful strings
  var_type <- c("tseries", "tseries-sd", "trend", "trend-stat")
  var_type_long <- c(
    "Timeseries",
    "St.dev of timeseries",
    "Trend coeff. for two intervals ",
    "Trend statistics for trend 1 (Estimate, Std. Error, t value, Pr(>|t|))"
  )

  # setup parameters
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]

  # set main paths
  outfile <- getfilename_trends(work_dir, label, model_idx, season)

  # Define regions to be used
  nregions <- length(selregions)

  # Define fields to be used (main list loaded from cfg_file)
  if (selfields[1] != F) {
    field_names <- field_names[selfields]
  }

  # Years to be considered based on namelist and cfg_file
  years <- year1:year2
  nyears <- length(years)

  print(paste0(diag_base, ": starting timeseries calculation"))

  #-----------------Loading data-----------------------#
  # open experiment field

  gridfile <-
    getfilename_indices(work_dir, diag_base, model_idx, grid = T)
  infile <-
    getfilename_indices(work_dir, diag_base, model_idx, season)
  # test if file contains all requested variables and
  # keep file open for reading attributes
  nc <- nc_open(infile)
  nc_att_glob <- ncatt_get(nc, 0)
  if (!all(is.element(field_names, names(nc$var)))) {
    missing <- (field_names)[!is.element(field_names, names(nc$var))]
    print(paste(
      "HyInt_trends: missing variable in input indices file: ",
      missing
    ))
    nc_close(nc)
    stop("HyInt: check field_names list in configuration file")
  }
  nc_close(nc)

  # Get seaLandElevation mask
  if (masksealand) {
    topofile <-
      getfilename_indices(work_dir, diag_base, ref_idx, topo = T)
    if (!file.exists(topofile)) {
      create_landseamask(
        regrid = gridfile,
        ref_file = infile,
        loc = run_dir,
        regridded_topo = topofile,
        topo_only = T
      )
    }
    relevation <- ncdf_opener(topofile, "topo", rotate = "no")
  }

  # remove desert areas if required
  # (mean annual precipitation <0.5 mm, Giorgi et al. 2014)
  if (removedesert) {
    pry <- ncdf_opener(infile, "pry", rotate = "no")
    retdes <- which(pry < 0.5)
    pry[retdes] <- NA
    # create mask with NAs for deserts and 1's for not-desert
    retdes2D <- apply(pry * 0, c(1, 2), sum) + 1
    retdes3D <- replicate(dim(pry)[length(dim(pry))], retdes2D)
  }

  for (var in field_names) {
    rfield <- ncdf_opener(infile, var, rotate = "no")
    print("===========================================")
    print(paste(infile, var))

    if (removedesert) {
      rfield <- rfield * retdes3D
    }
    if (masksealand) {
      rfield <- apply_elevation_mask(rfield, relevation, sealandelevation,
        reverse = reverse_masksealand
      )
    }
    # store size of time array
    ntime <- length(rfield[1, 1, ])

    #-----------------Calculating timeseries and trends-----------------------#

    # TIMESERIES:
    # - select required region and calculate timeseries
    # - timeseries are temporarily stored as a "region x time" matrix
    # - trends are temporarily stored as a "region x coefficient" matrix
    tfield <- matrix(nrow = nregions, ncol = ntime)
    tfield_sd <- matrix(nrow = nregions, ncol = ntime)
    rtrend <- matrix(nrow = nregions, ncol = 4)
    rtrend_stat <- matrix(nrow = nregions, ncol = 4)
    for (ireg in 1:nregions) {
      iselreg <- selregions[ireg]
      # extract data and perform averages
      print(paste("Working on ", region_names[iselreg]))

      tfield[ireg, ] <- calc_region_timeseries(ics, ipsilon, rfield,
        regions[iselreg, ],
        weighted_mean = weight_tseries
      )
      tfield_sd[ireg, ] <-
        calc_region_timeseries(ics, ipsilon, rfield,
          regions[iselreg, ],
          calc_sd = T
        )
    }

    # setup time array
    times <- as.numeric(year1) + 1:ntime - 1
    rettimes <- 1:length(times)
    if (trend_years[1] != F) {
      # apply trend to limited time interval if required
      rettimes <- which((times >= trend_years[1]) &
        times <= trend_years[2])
      if (length(trend_years) == 4) {
        # apply trend also to second time interval if required
        rettimes2 <- which((times >= trend_years[3]) &
          times <= trend_years[4])
      }
    }

    # LOOP through regions to calculate trends as required
    for (ireg in 1:nregions) {
      iselreg <- selregions[ireg]
      if (lm_trend) {
        # linear regression
        print("-----------------------------------------------------")
        print(paste(var, region_names[iselreg]))
        temp.tfield <- tfield[ireg, rettimes]
        if (length(which(!is.na(temp.tfield))) < 2) {
          print("less than 2 points in selected region - skipping")
        } else {
          lm_fit <- lm(temp.tfield ~ times[rettimes])
          lm_sum <- summary(lm_fit)
          # store trend coefficients (intercept and linear coef.)
          rtrend[ireg, 1:2] <- lm_fit$coefficients
          # store trend coef., standard error, t value, Pr(>|t|)
          rtrend_stat[ireg, ] <- lm_sum$coefficients[2, ]
          print(lm_sum$coefficients[2, ])
          if (length(trend_years) == 4) {
            # apply trend also to second time interval if required
            temp_tfield2 <- tfield[ireg, rettimes2]
            if (length(which(!is.na(temp.tfield))) < 2) {
              print("less than 2 points in second trend over selected region
                    - skipping")
            } else {
              lm_fit2 <- lm(temp_tfield2 ~ times[rettimes2])
              # store 2nd interval trend coefficients
              rtrend[ireg, 3:4] <- lm_fit2$coefficients
            }
          }
        }
      }
    }

    # assign timeseries and trends to named field variables
    assign(paste0(var, "_tseries"), tfield)
    assign(paste0(var, "_tseries-sd"), tfield_sd)
    assign(paste0(var, "_trend"), rtrend)
    assign(paste0(var, "_trend-stat"), rtrend_stat)
  } # close loop over fields

  # store field variables in named lists
  stseries_list <- c(
    paste0(field_names, "_tseries"),
    paste0(field_names, "_tseries-sd"),
    paste0(field_names, "_trend"),
    paste0(field_names, "_trend-stat")
  )
  rtseries_list <- mget(stseries_list)
  names(rtseries_list) <- stseries_list

  ##########################################################
  #------------------------Save to NetCDF------------------#
  ##########################################################

  # saving output to netcdf files
  print(paste0(diag_base, "_timeseries: saving data to NetCDF file:"))


  # dimensions definition
  var_region <- 1:nregions
  regiondim <- ncdim_def("region", "number", var_region)
  coeffdim <- ncdim_def("coefficients", "number", 1:4)
  boundarydim <- ncdim_def("boundaries", "degrees", 1:4)
  timedim <-
    ncdim_def(timedimname,
      "years since 1950-01-01 00:00:00",
      (years - 1950),
      unlim = T
    )

  # variables definition
  for (var in field_names) {
    for (itype in 1:length(var_type)) {
      svar <- paste0(var, "_", var_type[itype])
      rfield <- get(svar, rtseries_list)
      rfield[is.nan(rfield)] <- NA
      # copy and update attributes
      metadata <- getmetadata_indices(var, infile)
      long_name <- metadata$long_name
      description <<-
        paste0(var_type_long[itype], " of ", metadata$long_name)
      units <- metadata$units
      missval <- metadata$missing_value
      # variable definitions
      var_ncdf <- ncvar_def(
        svar,
        units,
        list(regiondim, timedim),
        missval,
        longname = long_name,
        prec = "single",
        compression = 1
      )
      if (itype > 2) {
        # trends
        var_ncdf <- ncvar_def(
          svar,
          units,
          list(regiondim, coeffdim),
          missval,
          longname = long_name,
          prec = "single",
          compression = 1
        )
      }
      assign(paste0("var", svar), var_ncdf)
      assign(paste0("field", svar), rfield)
      assign(paste0(svar, "_", "description"), description)
    }
  }

  varregions <- ncvar_def(
    "regions",
    "degrees",
    list(regiondim, boundarydim),
    -999,
    "region boundaries",
    prec = "single",
    compression = 1
  )
  regions_description <- "regions over which averages are performed"
  fieldregions <- regions[selregions, ]
  fieldregion_names <- region_names[selregions]
  fieldregion_codes <- region_codes[selregions]

  # Netcdf file creation
  print(paste0(diag_base, ": saving output to ", outfile))
  namelist <- c("regions", stseries_list)
  varnamelist <- paste0("var", c(namelist))
  nclist <- mget(varnamelist)
  ncfile <- nc_create(outfile, nclist)

  # put variables into the ncdf file
  for (var in namelist) {
    ndims <- get(paste0("var", var))$ndims
    tmp.field <- get(paste0("field", var))
    ncvar_put(ncfile,
      var,
      tmp.field,
      start = rep(1, ndims),
      count = rep(-1, ndims)
    )
    ncatt_put(ncfile, var, "description", get(paste0(var, "_description")))
  }

  # put additional attributes into dimension and data variables
  ncatt_put(
    ncfile,
    "regions",
    "regionnames",
    paste(fieldregion_names,
      collapse = " "
    )
  )
  ncatt_put(
    ncfile,
    "regions",
    "region_codes",
    paste(fieldregion_codes, collapse = " ")
  )

  nc_close(ncfile)

  # Set provenance for this output file
  caption <- paste0(
    "Hyint timeseries and trends for years ",
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

  print(paste(diag_base, ": timeseries netCDF file saved"))
  return(prov_info)
}
