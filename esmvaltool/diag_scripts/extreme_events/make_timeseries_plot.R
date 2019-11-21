# #############################################################################
# make_timeseries_plot.R
#
# Author: Marit Sandstad (CICERO, Norway)
#       : Christian W. Mohr (CICERO, Norway)
#
# #############################################################################
# Description
#    Code to plot a timeseries plot for a set of climdex indices
#
# Modification history
#    20190506-vonhardenberg_jost:    conversion to ESMValTool2
#    20180816-mohr_christianwilhelm: adding input procedure and plotting for/of
#                                    observation data
#    20180725-mohr_christianwilhelm: modification of time croppin
#    20180618-mohr_christianwilhelm: alpha levels for polygon plotting,
#                                    second y-axis,
#    20180131-lauer_axel:            clean-up of code, adaptation to ESMValTool
#                                    standards, added tagging, bugfixes: time
#                                    axis, cdo, filenames
#    2017 0920-sandstad_marit:       creation
#
# #############################################################################

##
##
## Method to call all preprocessing and loop through
## all models and indices then call plotting script
## to produce time series plots for a list of indices
## @param path is the path to where the original indices
## are stored
## @param idx_list lists the indices to be considered in
## this run. Defaults are the indices from the IPCC
## report.
##

###########################


timeseries_main <- function(path = "../work/extreme_events",
                            idx_list,
                            model_list,
                            obs_list,
                            plot_dir = "./plot",
                            normalize = FALSE,
                            start_yr = 2000,
                            end_yr = 2006) {
  ## For file structure and files
  tsgrid <- paste(path, "/tsGridDef", sep = "") # nolint
  time_cropped <- paste(path, "/timeCropped", sep = "") # nolint
  landmask <- paste(path, "/landSeaMask.nc", sep = "") # nolint
  regridded <- paste(path, "/regridded", sep = "") # nolint
  land <- paste(path, "/Land", sep = "") # nolint

  # Initial nc-file time crop, regrid, land and plot purge
  unlink(c(time_cropped, regridded, land, landmask, tsgrid),
    recursive = TRUE
  )

  # Initial grid and landmask creation reset
  gridandlandmask <- TRUE

  ## Loop over the indices to produce a plot for each index
  plotfiles <- list()
  idx <- idx_list[1]
  for (idx in idx_list) {
    ## Combine the list of models and observations
    modelobs_list <- unique(c(model_list, obs_list))

    ## Find the model files
    modelandobs <- basename(Sys.glob(file.path(
      path,
      paste(
        idx, "*.nc",
        sep = ""
      )
    )))

    if (ts_data) {
      ## Time crop
      returnvalue <- set_time_for_files_equal(
        path = path,
        idx = idx,
        model_list = modelobs_list,
        time_cropped = time_cropped
      ) # This is a temporary solution

      max_start <- returnvalue[1]
      min_end <- returnvalue[2]

      ## If there is no overlap in the files the index
      ## should be skipped
      if (max_start >= min_end) {
        print(paste("No time overlap in files for index", idx))
        break
      }

      ## Find the new model files after time cropping
      modelandobs <- basename(Sys.glob(file.path(
        time_cropped,
        paste0(idx, "*.nc")
      )))

      # !New Grid and landseamask for each idx
      # !(or just the first idx set) should be
      # !produced here
      if (gridandlandmask) {
        create_grid(path = path, loc = tsgrid)
        create_land_sea_mask(
          regrid = tsgrid,
          loc = path,
          landmask = landmask
        )
        gridandlandmask <- FALSE
      }

      ## Loop over each file so it can be regridded
      ## and landseaMasked
      for (m in modelandobs) {
        print(paste(time_cropped, "/", m, sep = ""))
        regrid_and_land_sea_mask(
          idx_raw = paste0(time_cropped, "/", m),
          regrid = tsgrid,
          landmask = landmask,
          regridded = regridded,
          land = land,
          loc = path
        )
      }

      ## Then do the preprocessing

      time_series_preprocessing(
        land = land,
        idx = idx,
        model_list = model_list,
        obs_list = obs_list,
        plot_dir = plot_dir,
        work_dir = path,
        normalize = normalize
      )
    }

    ## Produce plot for this index
    fname <- timeseries_plot(
      plot_dir = plot_dir,
      idx = idx,
      obs_list = obs_list,
      start_yr = start_yr,
      end_yr = end_yr,
      normalize = normalize
    )
    plotfiles <- c(plotfiles, fname)
  }

  # Final cleanup
  unlink(c(time_cropped, regridded, land, landmask, tsgrid),
    recursive = TRUE
  )

  return(plotfiles)
}

#
# Method that preprocesses idx-files for a single index
# in order to get data to plot time series plot
# for this index
# @param path is the path to index file location
# @param idx is the index to be processed.
#
time_series_preprocessing <- # nolint
  function(land = "./Land",
             idx = "tnnETCCDI_yr",
             model_list = model_list,
             obs_list = obs_list,
             plot_dir = "./plot",
             work_dir = "./work",
             normalize = FALSE) {
    tseriesdir <- paste0(work_dir, "/timeseries")
    if (!file.exists(tseriesdir)) {
      dir.create(tseriesdir)
    }

    ## List of indices which are never normalized:
    pidx <- c(
      "tn10pETCCDI_yr",
      "tx10pETCCDI_yr",
      "tn90pETCCDI_yr",
      "tx90pETCCDI_yr",
      "csdiETCCDI_yr",
      "wsdiETCCDI_yr",
      "tn10pETCCDI_mon",
      "tx10pETCCDI_mon",
      "tn90pETCCDI_mon",
      "tx90pETCCDI_mon",
      "csdiETCCDI_mon",
      "wsdiETCCDI_mon"
    )

    # Getting a list of all the files for the index
    modelsandobs <-
      basename(Sys.glob(file.path(land, paste0(idx, "*.nc"))))
    modelsandobssplitlist <- strsplit(modelsandobs, split = "_")
    modelsandobssplit <-
      unlist(lapply(modelsandobssplitlist, function(x) {
        x[3]
      }))

    # Extracting only the model files
    models <- modelsandobs[which(modelsandobssplit %in% model_list)]
    print("These are the models:")
    print(models)

    ## Extracting only the observation files
    obs_order <- which(modelsandobssplit %in% obs_list)
    obs <- modelsandobs[obs_order]
    print("These are the observations:")
    print(obs)

    #### NORMALIZE VALUES ####
    if (normalize) {
      # File string to be filled with file names that can be used
      # For the aggregated statistics
      # (ensmean, enspctl)
      file_string_models <- ""
      m <- models[1]
      for (m in models) {
        print(m)
        if (idx %in% pidx) {
          # Fieldmean results
          cdo(
            "fldmean",
            input = paste0(land, "/", m),
            output = paste0(land, "/", "fldm_", m),
            options = "-O"
          )

          ## add the preprocessed file to the filestring
          file_string_models <- paste(file_string_models, land,
            "/fldm_", m, " ",
            sep = ""
          ) # nolint
        } else {
          # Subtracting timemeans from land files:
          cdo(
            "sub",
            input = c(paste0(land, "/", m), paste0(land, "/tm_", m)),
            # nolint
            output = paste0(land, "/", "norm_", m),
            options = "-O"
          )

          # Detrended results:
          cdo(
            "detrend",
            input = paste0(land, "/norm_", m),
            # nolint
            output = paste0(land, "/", "detrend_", m),
            options = "-O"
          )

          # Timstd of detrend
          cdo(
            "timstd",
            input = paste0(land, "/detrend_", m),
            # nolint
            output = paste0(land, "/", "detrend_std_", m),
            options = "-O"
          )

          # Divide normalized by timstded detrend
          cdo(
            "div",
            input = c(
              paste0(land, "/norm", m),
              # nolint
              paste0(land, "/detrend_std_", m) # nolint
            ),
            output = paste0(land, "/detrend_standard_", m),
            options = "-O" # nolint
          )

          # Fieldmean results
          cdo(
            "fldmean",
            input = paste0(land, "/detrend_standard_", m),
            # nolint
            output = paste0(land, "/detrend_std_fldm_", m),
            options = "-O" # nolint
          )

          ## add the preprocessed file to the filestring
          file_string_models <- paste(file_string_models,
            land,
            "/detrend_std_fldm_",
            m,
            " ",
            sep = ""
          ) # nolint
        }
      }
      # Find model ensemble mean
      cdo(
        "ensmean",
        input = file_string_models,
        output = paste0(tseriesdir, "/", idx, "_ensmean_for_timeseries.nc"),
        options = "-O"
      )

      # Find ensemble 25th percentile
      cdo(
        "enspctl",
        args = "25",
        input = file_string_models,
        output = paste0(tseriesdir, "/", idx, "_25enspctl_for_timeseries.nc"),
        options = "-O"
      )

      # Find ensemble 75th percentile
      cdo(
        "enspctl",
        args = "75",
        input = file_string_models,
        output = paste0(tseriesdir, "/", idx, "_75enspctl_for_timeseries.nc"),
        options = "-O"
      )

      n <- 0
      for (o in obs) {
        print(o)

        if (idx %in% pidx) {
          # Fieldmean results
          cdo(
            "fldmean",
            input = paste0(land, "/", o),
            output = paste0(land, "/", "fldm_", o),
            options = "-O"
          )

          # Copy obs file to plot
          n <- n + 1
          file.copy(
            paste0(land, "/fldm_", o),
            # nolint
            paste0(
              tseriesdir,
              "/",
              idx,
              "_",
              modelsandobssplit[obs_order[n]],
              "_for_timeseries.nc"
            )
          )
        } else {
          # Subtracting timemeans from land files:
          cdo(
            "sub",
            input = c(paste0(land, "/", o), paste0(land, "/tm_", o)),
            # nolint
            output = paste0(land, "/norm_", o),
            options = "-O" # nolint
          )

          # Detrended results:
          cdo(
            "detrend",
            input = paste0(land, "/norm_", o),
            # nolint
            output = paste0(land, "/detrend_", o),
            options = "-O" # nolint
          )

          # Timstd of detrend
          cdo(
            "timstd",
            input = paste0(land, "/detrend_", o),
            # nolint
            output = paste0(land, "/detrend_std_", o),
            options = "-O" # nolint
          )

          # Divide normalized by timstded detrend
          cdo(
            "div",
            input = c(
              paste0(land, "/norm", o),
              # nolint
              paste0(land, "/detrend_std_", o) # nolint
            ),
            output = paste0(land, "/detrend_standard_", o),
            options = "-O" # nolint
          )

          # Fieldmean results
          cdo(
            "fldmean",
            input = paste0(land, "/detrend_standard_", o),
            # nolint
            output = paste0(land, "/detrend_std_fldm_", o),
            options = "-O" # nolint
          )

          # Copy obs file to plot
          n <- n + 1
          file.copy(
            paste0(land, "/detrend_std_fldm_", o),
            # nolint
            paste0(
              tseriesdir,
              "/",
              idx,
              "_",
              modelsandobssplit[obs_order[n]],
              "_for_timeseries.nc"
            )
          )
        }
      }
    }

    # ABSOLUTE VALUES ####
    # Non-normalized values fieldmeans
    if (!normalize) {
      file_string_models <- ""
      m <- models[1]
      for (m in models) {
        print(m)
        # Fieldmean results
        cdo(
          "fldmean",
          input = paste0(land, "/", m),
          output = paste0(land, "/fldm_", m),
          options = "-O" # nolint
        )

        ## add the preprocessed file to the filestring
        file_string_models <- paste(file_string_models, land,
          "/fldm_", m, " ",
          sep = ""
        ) # nolint
      }
      # Find model ensemble mean
      cdo(
        "ensmean",
        input = file_string_models,
        output = paste0(tseriesdir, "/", idx, "_ensmean_for_timeseries.nc"),
        options = "-O"
      )

      # Find ensemble 25th percentile
      cdo(
        "enspctl",
        args = "25",
        input = file_string_models,
        output = paste0(tseriesdir, "/", idx, "_25enspctl_for_timeseries.nc"),
        options = "-O"
      )

      # Find ensemble 75th percentile
      cdo(
        "enspctl",
        args = "75",
        input = file_string_models,
        output = paste0(tseriesdir, "/", idx, "_75enspctl_for_timeseries.nc"),
        options = "-O"
      )

      ## Extracting only the observation files
      obs_order <- which(modelsandobssplit %in% obs_list)
      obs <- modelsandobs[obs_order]

      print("These are the observations:")
      print(obs)
      n <- 0
      for (o in obs) {
        print(o)
        # Fieldmean results
        cdo(
          "fldmean",
          input = paste0(land, "/", o),
          output = paste0(land, "/fldm_", o),
          options = "-O" # nolint
        )
        # Copy obs file to plot
        n <- n + 1
        file.copy(
          paste0(land, "/fldm_", o),
          # nolint
          paste0(
            tseriesdir,
            "/",
            idx,
            "_",
            modelsandobssplit[obs_order[n]],
            "_for_timeseries.nc"
          )
        )
      }
    }
  }

#
#
# Method to plot the time series plot
# of single idx for already preprocessed data
# yearly data is assumed
# @param path - path to directory containing ensemble mean
# and percentile data.
# @param idx name of index to be processed
# @param start_yr start year for data to be used to convert
# values from days after start year format to
# year.
#

timeseries_plot <-
  function(plot_dir = "./plot",
             idx = "tn10pETCCDI_yr",
             obs_list,
             start_yr = 2006,
             end_yr = 2010,
             normalize = FALSE) {
    # Drawing parameters
    leg_names <- c(mip_name, obs_list)

    ## Reading the netcdf data files into R
    ## First ensemble mean file
    ensm <- nc_open(paste(
      work_dir,
      "/timeseries/",
      idx,
      "_ensmean_for_timeseries.nc",
      sep = ""
    ))

    ## Then 25th percentile file
    enspctl25 <- nc_open(paste(
      work_dir,
      "/timeseries/",
      idx,
      "_25enspctl_for_timeseries.nc",
      sep = ""
    ))
    ## Finally 75th percentile file
    enspctl75 <- nc_open(paste(
      work_dir,
      "/timeseries/",
      idx,
      "_75enspctl_for_timeseries.nc",
      sep = ""
    ))

    ## Reading in time variable and converting to years:
    ts <- nc.get.time.series(ensm) # nolint
    time_conv <- format(ts, "%Y") # extract years

    ## Stripping off the _yr tail to the index name
    idx_no <- which(idx_df$idx_etccdi_time == idx)
    idx_name <- paste(idx_df$idx_etccdi[idx_no], "ETCCDI", sep = "")

    # Reading in the y-variables to be plotted
    # First the ensemble mean
    idx_ensm <- ncvar_get(ensm, idx_name)
    # Then the 25th percentile
    idx_ens25 <- ncvar_get(enspctl25, idx_name)
    # Finally the 75th percentile
    idx_ens75 <- ncvar_get(enspctl75, idx_name)

    # Maximum and minimum x and y values
    max.x <- end_yr
    min.x <- start_yr
    irange <- ((time_conv >= min.x) & (time_conv <= max.x))
    max.y <-
      max(idx_ensm[irange], idx_ens25[irange], idx_ens75[irange])
    min.y <-
      min(idx_ensm[irange], idx_ens25[irange], idx_ens75[irange])

    # Reading in the observations and plotting via a loop
    obsdata_list <- list()
    n <- 0
    for (o in obs_list) {
      n <- n + 1
      nc_obs <- nc_open(paste(
        work_dir,
        "/timeseries/",
        idx,
        "_",
        o,
        "_for_timeseries.nc",
        sep = ""
      ))
      ts_obs <- nc.get.time.series(nc_obs) # nolint
      time_conv_obs <- format(ts_obs, "%Y") # extract years
      idx_obs <- ncvar_get(nc_obs, idx_name)
      nc_close(nc_obs)
      obsdata_list[[n]] <- list(o, as.numeric(time_conv_obs), idx_obs)
      irange <- ((time_conv_obs >= min.x) & (time_conv_obs <= max.x))
      max.y <- max(max.y, idx_obs[irange])
      min.y <- min(min.y, idx_obs[irange])

      if (n > length(ts_col_list)) {
        print(
          paste(
            "Error: There are more observations,",
            "than available color plotting parameters."
          )
        )
        print("Update cfg_ExtermeEvents.r file.")
        dev.off()
        break
      }
    }

    # Setting the x- and y-range limits for plotting
    xrng <- as.numeric(c(min.x, max.x))
    yrng <- c(min.y, max.y)
    print(xrng)
    print(yrng)

    ## Making name string for the plot
    plotname <- paste(plot_dir,
      "/",
      idx,
      "_",
      length(obs_list),
      "-obs_ensmean_timeseriesplot",
      sep = ""
    )

    ## Setting device to write the plot to
    figure_filename <- paste(plotname, output_file_type, sep = ".")

    ## Chose output format for figure
    if (tolower(output_file_type) == "png") {
      png(
        filename = figure_filename,
        width = ts_png_width,
        height = ts_png_height,
        units = ts_png_units,
        pointsize = ts_png_pointsize,
        bg = ts_png_bg
      )
    } else if (tolower(output_file_type) == "pdf") {
      pdf(file <- figure_filename)
    } else if (tolower(output_file_type) == "eps") {
      setEPS()
      postscript(figure_filename)
    }

    n <- 1
    # Parameters for plot
    par(mfrow = c(1, 1), mar = c(4.5, 4.5, 2, 3))
    # Plotting first the ensemblemean
    plot(
      time_conv,
      idx_ensm,
      type = "l",
      col = ts_col_list[n],
      lty = ts_lty_list[n],
      xlim = xrng,
      ylim = yrng,
      lwd = ts_lwd_list[n],
      ann = FALSE,
      xaxs = "i",
      yaxt = "n"
    )
    # Then making a transparent polygon between the 25th and 75 percentile
    polygon(
      c(time_conv, rev(time_conv)),
      c(idx_ens75, rev(idx_ens25)),
      col = alpha(ts_col_list[n], 0.1),
      border = NA
    )

    # Plotting observations and plotting via a loop
    n <- 0
    for (o in obs_list) {
      n <- n + 1
      lines(
        obsdata_list[[n]][[2]],
        obsdata_list[[n]][[3]],
        col = ts_col_list[n + 1],
        lty = ts_lty_list[n + 1],
        lwd = ts_lwd_list[n + 1]
      ) # plot observation
    }

    # Produce a legend
    legend(
      "top",
      legend = leg_names,
      col = ts_col_list,
      lty = ts_lty_list,
      lwd = ts_lwd_list,
      bty = "n",
      ncol = 3
    )

    # Produce a first y-axis
    axis(side = 2, at = pretty(yrng, 5))
    axis(side = 2, at = pretty(yrng, 5))
    pretty(yrng, 10)

    axis(side = 2, at = pretty(yrng, 5))

    # Produce a second y-axis
    axis(side = 4, at = pretty(yrng, 5))

    # Producing a title from info in netcdf file
    title(main = idx_df$name[idx_no], font.main = 2)

    # Choosing x-label
    title(xlab = "Year")

    # Chosing y-label from idx_ylab list
    title(ylab = idx_ylab[idx_no])
    # Resetting plotting device to default
    dev.off()

    # Close Ensemble files
    nc_close(ensm)
    nc_close(enspctl25)
    nc_close(enspctl75)

    return(figure_filename)
  }
