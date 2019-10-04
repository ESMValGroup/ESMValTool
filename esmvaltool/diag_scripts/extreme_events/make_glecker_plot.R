# ############################################################################
# make_glecker_plot.R
#
# Author: Christian W. Mohr (CICERO, Norway)
#         Marit Sandstad (CICERO, Norway)
#
#
# ############################################################################
# Description:
# Code to plot Glecker polygon diagram to compare climdex index
# 	performance between models and reanalysis.
#
# Modification history
#
#    20190506-vonhardenberg_jost: conversion to ESMValTool2
#    20180601-mohr_christianwilhelm: re-creation (complete new script
#                                    incorparating segments from
#                                    "make_timeseries_plot.r" and
#                                    "make_Glecker_plot.r")
#
# ############################################################################

gleckler_main <-
  function(path = "./",
             idx_list,
             model_list,
             obs_list,
             plot_dir = "../plot/extreme_events/",
             promptinput = promptinput,
             start_yr = 2000,
             end_yr = 2009) {
    #### CLIMDEX PREPROCESSING ####

    ## For file structure and files
    tsgrid <- paste(path, "/tsGridDef", sep = "") # nolint
    time_cropped <- paste(path, "/timeCropped", sep = "") # nolint
    landmask <- paste(path, "/landSeaMask.nc", sep = "") # nolint
    regridded <- paste(path, "/regridded", sep = "") # nolint
    land <- paste(path, "/Land", sep = "") # nolint

    nmodel <- length(model_list) # number of models
    nidx <- length(idx_list) # number of indices
    nobs <- length(obs_list) # number of observations

    if (file.exists(
      paste0(
        path,
        "/gleckler/Gleckler-Array_",
        nidx,
        "-idx_",
        nmodel,
        "-models_",
        nobs,
        "-obs",
        ".RDS"
      )
    )) {
      promptinput <- "n"
    }

    if (promptinput == "y") {
      # Initial nc-file time crop, regrid, land and plot purge
      unlink(c(
        time_cropped, regridded, land,
        landmask, tsgrid
      ),
      recursive = TRUE
      )

      ## Initial grid and landmask creation reset
      grid_and_landmask <- TRUE

      ## Combine model and observation list
      modelandobs_list <- unique(c(model_list, obs_list))

      ## Loop over the indices to produce a plot for each index
      for (idx in idx_list) {
        ## Time crop
        returnvalue <- set_time_for_files_equal(
          path = path,
          idx = idx,
          model_list = modelandobs_list,
          time_cropped = time_cropped,
          max_start = start_yr,
          min_end = end_yr
        )

        max_start <- returnvalue[1]
        min_end <- returnvalue[2]

        ## If there is no overlap in the files the index
        ## should be skipped
        if (max_start >= min_end) {
          print(paste("No time overlap in files for index", idx))
          break
        }

        ## Find the new model and observation names (after time cropping)
        modelsandobs <- basename(Sys.glob(file.path(
          time_cropped,
          paste0(idx, "*.nc")
        )))
        split_modelsandobs <- strsplit(modelsandobs, split = "_")
        modelsandobs_index <-
          unlist(lapply(split_modelsandobs, function(x) {
            x[3]
          }))

        ## new models
        models <-
          modelsandobs[which(modelsandobs_index %in% model_list)]

        ## new observations
        obs <- modelsandobs[which(modelsandobs_index %in% obs_list)]

        ## Find the start year (to be used in plotting)
        # start_yr <- strtoi(substr(models[1], nchar(models[1]) - 11,
        #                          nchar(models[1]) - 8))
        # !New Grid and landseamask for each idx
        # !(or just the first idx set) should be
        # !produced here
        if (grid_and_landmask) {
          create_grid(path = path, loc = tsgrid)
          create_land_sea_mask(
            regrid = tsgrid,
            loc = path,
            landmask = landmask
          )
          grid_and_landmask <- FALSE
        }

        ## Loop over each file so it can be regridded
        ## and landseaMasked
        for (mo in modelsandobs) {
          print(paste(time_cropped, "/", mo, sep = ""))
          regrid_and_land_sea_mask(
            idx_raw = paste(time_cropped, "/", mo, sep = ""),
            regrid = tsgrid,
            landmask = landmask,
            regridded = regridded,
            land = land,
            loc = path
          )
        }
      }

      #### Gleckler Array Processing ####
      rmserelarr <- gleckler_array(
        path = land,
        idx_list = idx_list,
        model_list = model_list,
        obs_list = obs_list
      )

      ## Save Array
      glecdir <- paste0(path, "/gleckler") # nolint
      if (!file.exists(glecdir)) {
        dir.create(glecdir)
      }
      saveRDS(
        object = rmserelarr,
        file = paste0(
          path,
          "/gleckler/Gleckler-Array_",
          nidx,
          "-idx_",
          nmodel,
          "-models_",
          nobs,
          "-obs",
          ".RDS"
        )
      )
      saveRDS(
        object = returnvalue,
        file = paste0(
          path,
          "/gleckler/Gleckler-years.RDS"
        )
      )

      # Final cleanup
      unlink(c(
        time_cropped, regridded, land,
        landmask, tsgrid
      ),
      recursive = TRUE
      )
    }

    #### Gleckler Plotting ####
    rmserelarr <- readRDS(
      file = paste0(
        path,
        "/gleckler/Gleckler-Array_",
        nidx,
        "-idx_",
        nmodel,
        "-models_",
        nobs,
        "-obs",
        ".RDS"
      )
    )
    year_range <- readRDS(file = paste0(
      path,
      "/gleckler/Gleckler-years.RDS"
    ))

    plotfile <-
      gleckler_plotting(
        arr = rmserelarr,
        idx_list = idx_list,
        model_list = model_list,
        obs_list = obs_list,
        plot_dir = plot_dir,
        syear = year_range[1],
        eyear = year_range[2]
      )
    return(plotfile)
  }

#### Computing the RMSEs ####

gleckler_array <- function(path = land,
                           idx_list = gleckler_idx,
                           model_list = model_list,
                           obs_list = obs_list) {
  ## Produce an array to hold all the model and reanalysis means

  ## Input data for testing the plotting routine
  nidx <- length(idx_list) # number of indices
  nmodel <- length(model_list) # number of models
  nobs <- length(obs_list) # number of reanalyses

  ## Check point for reanalysis data
  if (nobs == 0) {
    print("No reanalysis datasets provided")
    break
  }

  ## Function to calculate area mean
  area.mean <- function(x, lat) {
    nlon <- dim(x)[1]
    nlat <- dim(x)[2]

    meanlat <- apply(x, 2, function(x) {
      mean(x, na.rm = TRUE)
    })

    fi <- lat * 3.14159 / 180

    wgt.prod <- meanlat * cos(fi)

    # At some latitudes there is no land and therfore no data.
    nan.check <- is.nan(wgt.prod)
    # The mean of missing data is not a number, and hench results in NaNs.
    # These NaN must be removed in order to calculate the correct area mean.

    gl <- sum(wgt.prod[!nan.check])
    sumcos <- sum(cos(fi)[!nan.check])
    ar.m <- gl / sumcos
    return(ar.m)
  }

  ## Function to calculate the RMSE between the model and
  ## observed climatology (RMSExy)
  ## Equation 1, from Sillmann et. al 2013
  RMSE <-
    function(model = tm_model_idx,
                 obs = tm_obs_idx,
                 lat = model_lat) {
      RMSE <- sqrt(area.mean((model - obs)^2, lat))
      return(RMSE)
    }

  # Array for the RMSE spaces in the array are created so that the
  # RSMEall, ENSmean, ENSmedian and CMIP RMSE can be created
  rmsearr <- array(NA, dim = c(nidx + 1, nmodel + 3, nobs))
  rmserelarr <- rmsearr
  ensmodel_list <- list()

  i <- 2
  m <- 1
  o <- 1
  lat_collect <- TRUE

  for (i in seq_along(idx_list)) {
    for (m in seq_along(model_list)) {
      ## Read in model annual climatology

      tm_model <- nc_open(Sys.glob(file.path(
        path, paste0(
          "tm_", idx_list[i],
          "_", model_list[m],
          "*.nc"
        )
      )))
      idxs <- unlist(strsplit(idx_list[i], split = "_"))[1]
      tm_model_idx <- ncvar_get(tm_model, idxs)

      # extract latitudes for area mean calculations
      if (lat_collect) {
        model_lat <- ncvar_get(tm_model, "lat")
        lat_collect <- FALSE
      }

      nc_close(tm_model)
      ensmodel_list[[m]] <- tm_model_idx
    }
    ## Create a new array for adding the time mean model matices
    ensarr <- array(NA, dim = c(
      nrow(tm_model_idx),
      ncol(tm_model_idx),
      length(ensmodel_list) + 2
    ))

    # Copy each matrix from the multimodel list to the array "ensarr".
    # Notice the "+2" on the 3rd dimension. This is so later the model
    # ensemble mean and median matrices can be added to the array.
    for (n in seq_along(ensmodel_list)) {
      ensarr[, , n + 2] <- ensmodel_list[[n]]
    }

    ## Calculate the ensemble mean and median of
    ## all the model time mean matrices
    ensmean <- apply(ensarr, c(1, 2), function(x) {
      mean(na.omit(x))
    })
    ensmedian <- apply(ensarr, c(1, 2), function(x) {
      median(na.omit(x))
    })

    # Place the ensemble model mean and medians into the
    # first two matrices (3-dimention) of the array "ensarr"
    ensarr[, , 1] <- ensmean
    ensarr[, , 2] <- ensmedian

    j <- 1
    ## Calculate the RMSE for all the models and the ensemble mean and median
    for (j in 1:dim(ensarr)[3]) {
      ## Read in reannalysis annual climatology
      for (o in seq_along(obs_list)) {
        tm_obs <- nc_open(Sys.glob(file.path(
          path, paste0(
            "tm_", idx_list[i],
            "_", obs_list[o],
            "*.nc"
          )
        )))
        tm_obs_idx <- ncvar_get(tm_obs, idxs)
        nc_close(tm_obs)
        rmsearr[i + 1, j, o] <- RMSE(
          model = ensarr[, , j],
          obs = tm_obs_idx,
          lat = model_lat
        ) # Calculate each RMSE and place value in RMSE-array

        ## Calculate the model standard deviation.
        ## Later used for calculating the rmsemedian,std.
        ## Denominator in equation 3, from Sillmann et. al 2013
        rmsearr[i + 1, ncol(rmsearr), o] <-
          sqrt(area.mean((
            tm_obs_idx - area.mean(tm_obs_idx,
              lat = model_lat
            )
          )^2,
          lat = model_lat
          ))
      }
    }
  }

  ## Calculate the RMSE median for the models
  tmprmsearr <- rmsearr[, -c(1, 2, ncol(rmsearr)), ]
  if (length(dim(tmprmsearr)) == 3) {
    rmsemed <- apply(tmprmsearr, c(1, 3), function(x) {
      median(x, na.rm = TRUE)
    })
  } else {
    rmsemed <- apply(tmprmsearr, 1, function(x) {
      median(x, na.rm = TRUE)
    })
  }

  ## Function to calculate the relative RMSE (RMSE'xy)
  ## between the model and observed climatology
  ## Equation 2, from Sillmann et. al 2013
  rmserel <- function(rmse, rmsemed) {
    rmserel <- (rmse - rmsemed) / rmsemed
    return(rmserel)
  }

  ## Calculating the relative RMSE (RMSE'xy)
  m <- 1
  for (m in 1:(ncol(rmsearr) - 1)) {
    rmserelarr[, m, ] <-
      rmserel(rmse = rmsearr[, m, ], rmsemed = rmsemed)
  }

  ## Calculating the RMSE median,std. Equation 3, from Sillmann et. al 2013
  rmserelarr[, ncol(rmserelarr), ] <-
    rmsemed / rmsearr[, ncol(rmsearr), ]

  ## Calculating the RSME mean
  tmprmsearr <- rmserelarr[, -ncol(rmserelarr), ]
  if (length(dim(tmprmsearr)) == 3) {
    rmserelarr[1, -ncol(rmserelarr), ] <- apply(
      tmprmsearr, c(2, 3),
      function(x) {
        mean(x, na.rm = TRUE)
      }
    )
  } else {
    rmserelarr[1, -ncol(rmserelarr), ] <- apply(
      tmprmsearr, c(2),
      function(x) {
        mean(x, na.rm = TRUE)
      }
    )
  }
  print(rmserelarr)
  return(rmserelarr)
}

#### Plotting Routine ####
gleckler_plotting <-
  function(arr = rmserelarr,
             idx_list,
             model_list,
             obs_list,
             plot_dir = "../plots/extreme_events/",
             syear = max_start,
             eyear = min_end) {
    nidx <- length(idx_list) # number of indices
    nmodel <- length(model_list) # number of models
    nobs <- length(obs_list) # number of reanalyses

    ## Numbers for color scale
    sclseq <- seq(-0.55, 0.55, 0.1)

    ## Colour scale
    glc <- brewer.pal(length(sclseq) - 2, "RdYlBu") # nolint
    glc <- c("#662506", glc, "#3f007d")
    glc <- rev(glc)

    # Numbers for black & white scale
    sclseq_bw <- seq(0.05, 1.15, 0.1)
    sclseq_bw
    glbw <- gray(seq(0, 1, length.out = length(sclseq_bw)))
    glbw <- rev(glbw)

    ## Determining what shapes should be plotted, based on number of
    ## observations
    if (nobs == 1) {
      # One reanalysis references
      x1 <- c(0, 1, 1, 0)
      y1 <- c(0, 0, 1, 1)
      xs <- list(x1)
      ys <- list(y1)

      # text coordinates
      xtx <- 0.50
      ytx <- -0.25
      rotx <- 0 # text rotation in degrees
    }

    if (nobs == 2) {
      # Two reanalysis references
      x1 <- c(0, 1, 1) # lower triangle
      y1 <- c(0, 0, 1) # lower triangle
      x2 <- c(0, 1, 0) # upper triangle
      y2 <- c(0, 1, 1) # upper triangle

      xs <- list(x1, x2)
      ys <- list(y1, y2)

      # text coordinates
      xtx <- c(0.75, 0.25)
      ytx <- c(-0.25, 1.25)
      rotx <- c(0, 0) # text rotation in degrees
    }

    if (nobs == 3) {
      # Three reanalysis references
      x1 <- c(0, 0.5, 0.5, 0) # bottom left
      y1 <- c(0, 0, 0.5, 1) # bottom left
      x2 <- c(0.5, 1, 1, 0.5) # bottom right
      y2 <- c(0, 0, 1, 0.5) # bottom right
      x3 <- c(0, 0, 0.5, 1, 1) # top
      y3 <- c(1, 0.75, 0.5, 0.75, 1) # top

      xs <- list(x1, x2, x3)
      ys <- list(y1, y2, y3)

      # text coordinates
      xtx <- c(-0.25, 1.25, 0.5)
      ytx <- c(0.25, 0.25, 1.25)
      rotx <- c(90, 90, 0) # text rotation in degrees
    }

    if (nobs == 4) {
      # Four reanalysis references
      x1 <- c(0, 0.5, 1) # bottom triangle
      y1 <- c(0, 0.5, 0) # bottom triangle
      x2 <- c(0, 0.5, 0) # left triangle
      y2 <- c(0, 0.5, 1) # left triangle
      x3 <- c(0, 0.5, 1) # top triangle
      y3 <- c(1, 0.5, 1) # top triangle
      x4 <- c(1, 0.5, 1) # right triangle
      y4 <- c(1, 0.5, 0) # right triangle

      xs <- list(x1, x2, x3, x4)
      ys <- list(y1, y2, y3, y4)

      # text coordinates
      xtx <- c(0.5, -0.25, 0.5, 1.25)
      ytx <- c(-0.25, 0.5, 1.25, 0.5)
      rotx <- c(0, 90, 0, 90) # text rotation in degrees
    }

    if (!(nobs %in% c(1, 2, 3, 4))) {
      if (nobs == 0) {
        print("No reanalysis dataset provided")
        break
      } else {
        print(
          paste(
            "Too many reanalysis datasets provided.",
            "Please choose between 1 and 4 datasets"
          )
        )
        break
      }
    }

    print("--- Creating Gleckler plot ---")
    img.adj <- gl_mar_par * 0.05
    width.fct <- ((nmodel + 3) / (nidx + 1)) + sum(img.adj[c(2, 4)])
    height.fct <- 1 + sum(img.adj[c(1, 3)])

    figure_filename <-
      paste(
        plot_dir,
        "/Gleckler_",
        mip_name,
        "_",
        nmodel,
        "-models_",
        nidx,
        "-idx_",
        nobs,
        "-obs_",
        syear,
        "-",
        eyear,
        ".",
        output_file_type,
        sep = ""
      )

    ## Chose output format for figure
    if (tolower(output_file_type) == "png") {
      png(
        filename = figure_filename,
        width = gl_png_res * (width.fct / height.fct),
        height = gl_png_res,
        units = gl_png_units,
        pointsize = gl_png_pointsize,
        bg = gl_png_bg
      )
    } else if (tolower(output_file_type) == "pdf") {
      pdf(file <- figure_filename)
    } else if (tolower(output_file_type) == "eps") {
      setEPS()
      postscript(figure_filename)
    }

    par(
      mfrow = c(1, 1),
      mar = gl_mar_par,
      xpd = FALSE,
      oma = rep(0, 4)
    )
    plot(
      x = c(0, 1 + gl_rmsespacer),
      y = c(0, 1),
      type = "n",
      ann = FALSE,
      xaxs = "i",
      yaxs = "i",
      bty = "n",
      xaxt = "n",
      yaxt = "n"
    )

    ## Array dimentions
    xn <- ncol(arr)
    yn <- nrow(arr)

    ## Testing array plotting
    xi <- 1 # model
    yj <- 2 # index
    zk <- 1 # obs

    ## Plotting RMSE of models, ensemble mean and median and RSMEall
    for (xi in 1:(xn - 1)) {
      for (yj in 1:yn) {
        for (zk in 1:nobs) {
          polygon(
            x = (xs[[zk]] / xn) + ((xi - 1) / xn),
            y = (ys[[zk]] / yn) + ((yn - yj) / yn),
            col = glc[which.min(abs(sclseq - arr[yj, xi, zk]))]
          )
        }
      }
    }

    ## Plotting RMSE median standard diviation
    for (yj in 2:yn) {
      for (zk in 1:nobs) {
        polygon(
          x = (xs[[zk]] / xn) + ((xn - 1) / xn) + gl_rmsespacer,
          y = (ys[[zk]] / yn) + ((yn - yj) / yn),
          col = glbw[which.min(abs(sclseq_bw - arr[yj, xn, zk]))]
        )
      }
    }

    ## Produce the borders for the Glecker plot
    par(xpd = TRUE)
    rect(
      xleft = 0,
      ybottom = 0,
      xright = (1 - 1 / xn),
      ytop = (1 - 1 / yn),
      density = NULL,
      angle = 45,
      col = NA,
      border = 1,
      lty = par("lty"),
      lwd = 4
    )
    rect(
      xleft = 0,
      ybottom = (1 - 1 / yn),
      xright = (1 - 1 / xn),
      ytop = 1,
      density = NULL,
      angle = 45,
      col = NA,
      border = 1,
      lty = par("lty"),
      lwd = 4
    )

    ## Scale for Gleckler plot
    gleckler_scale <- function(sclseq,
                                   glc,
                                   xn,
                                   scaling_factor,
                                   text.scaling_factor,
                                   xscale_spacer) {
      par(xpd = TRUE)
      ## Square legend
      sqrxs <- c(0, 1, 1, 0)
      sqrys <- c(0, 0, 1, 1)

      # up-triangle legend
      utrixs <- c(0, 1, 0.5)
      utriys <- c(0, 0, 1)

      # down-triangle legend
      dtrixs <- c(0.5, 1, 0)
      dtriys <- c(0, 1, 1)

      # Legend number shifter
      seq_shift <-
        mean(diff(sclseq) / 2) # Shifts the legend numbers so that
      # they represent the border values

      # y-scale spacer
      yscale_spacer <- (1 - scaling_factor) / 2

      exlen <- length(glc)
      for (a in 1:exlen) {
        if (a == 1) {
          xtmp <- scaling_factor * (dtrixs / xn) + 1 + xscale_spacer / xn
          ytmp <-
            (scaling_factor * (dtriys / exlen + (a - 1) / exlen) +
              yscale_spacer)
          polygon(
            x = xtmp,
            y = ytmp,
            col = glc[a]
          )
          text(
            x = max(xtmp),
            y = max(ytmp),
            round(sclseq[a] + seq_shift, 1),
            cex = text.scaling_factor,
            pos = 4
          )
        } else if (a == exlen) {
          xtmp <- scaling_factor * (utrixs / xn) + 1 + xscale_spacer / xn
          ytmp <-
            (scaling_factor * (utriys / exlen + (a - 1) / exlen) +
              yscale_spacer)
          polygon(
            x = xtmp,
            y = ytmp,
            col = glc[a]
          )
        } else {
          xtmp <- scaling_factor * (sqrxs / xn) + 1 + xscale_spacer / xn
          ytmp <-
            (scaling_factor * (sqrys / exlen + (a - 1) / exlen) +
              yscale_spacer)
          polygon(
            x = xtmp,
            y = ytmp,
            col = glc[a]
          )
          text(
            x = max(xtmp),
            y = max(ytmp),
            round(sclseq[a] + seq_shift, 1),
            cex = text.scaling_factor,
            pos = 4
          )
        }
      }
    }

    ## Plot scales
    gleckler_scale(
      sclseq,
      glc,
      xn,
      scaling_factor = gl_scaling_factor,
      text.scaling_factor = gl_text_scaling_factor,
      xscale_spacer = gl_xscale_spacer_rmse
    )

    gleckler_scale(
      sclseq_bw,
      glbw,
      xn,
      scaling_factor = gl_scaling_factor,
      text.scaling_factor = gl_text_scaling_factor,
      xscale_spacer = gl_xscale_spacer_rmsestd
    )

    ## Plotting symbol legend
    exlen <- length(glc)
    xsym1 <-
      gl_scaling_factor * (0.5 / xn) + 1 + gl_xscale_spacer_rmse / xn
    exlen <- length(glbw)
    xsym2 <-
      gl_scaling_factor * (0.5 / xn) + 1 + gl_xscale_spacer_rmsestd / xn
    x.max_adj <- max(gl_symb_scaling_factor * (xs[[zk]] / xn))
    x.min_adj <- min(gl_symb_scaling_factor * (xs[[zk]] / xn))
    xmidadj <- (x.max_adj - x.min_adj) / 2

    gl_symb_xshift <- (xsym1 + xsym2) / 2 - xmidadj

    for (zk in 1:nobs) {
      xsym <- gl_symb_scaling_factor * (xs[[zk]] / xn) + gl_symb_xshift
      ysym <- (gl_symb_scaling_factor * (ys[[zk]] / xn)
        - gl_symb_yshift / xn) * width.fct / height.fct
      print(paste("xs:", xsym))
      print(paste("ys:", ysym))
      polygon(
        x = xsym,
        y = ysym,
        col = "white",
        border = 1
      )

      xtxsym <-
        gl_symb_scaling_factor * (xtx[[zk]] / xn) + gl_symb_xshift
      ytxsym <- (gl_symb_scaling_factor * (ytx[[zk]] / xn)
        - gl_symb_yshift / xn) * width.fct / height.fct

      text(
        x = xtxsym,
        y = ytxsym,
        labels = obs_list[zk],
        adj = 0.5,
        cex = gl_text_symb_scaling_factor,
        srt = rotx[zk]
      )
    }

    ## Label adjusting parameters
    axlabsize <- 0.8
    lineadj <- -0.5

    ## Add model labels
    col_names <- c("ENSMEAN", "ENSMEDIAN", model_list)
    xtcks1 <- seq((0.5 / xn), ((xn - 1) / xn), by = (1 / xn))
    axis(
      side = 1,
      at = xtcks1,
      labels = col_names,
      las = 2,
      cex.axis = axlabsize,
      tick = FALSE,
      line = lineadj
    )

    xtcks2 <- ((xn - 1) / xn) + gl_rmsespacer + (0.5 / xn)
    axis(
      side = 1,
      at = xtcks2,
      labels = expression("RMSE"["std"]),
      las = 2,
      cex.axis = axlabsize,
      tick = FALSE,
      line = lineadj
    )

    ## Add index labels
    row_names <-
      vector(mode = "character", length = length(idx_list))
    for (i in seq_along(idx_list)) {
      row_names[i] <- idx_df$idx_etccdi[which(idx_df$idx_etccdi_time
        %in% idx_list[i])]
    }
    row_names <- rev(c(expression("RSME"["all"]), row_names))
    ytcks1 <- seq((1 / yn) * 0.5, 1, by = (1 / yn))
    axis(
      side = 2,
      at = ytcks1,
      labels = row_names,
      las = 2,
      cex.axis = axlabsize,
      tick = FALSE,
      line = lineadj
    )

    mtext(
      text = paste(mip_name, " global land ", syear, "-", eyear, sep = ""),
      side = 3,
      line = 1,
      font = 2,
      cex = 1.1
    )

    dev.off()
    return(figure_filename)
  }
