######################################################
#------Regimes routines figures for MiLES------------#
#-------------P. Davini (May 2017)-------------------#
######################################################

miles_regimes_figures <- function(dataset,
                                  expid,
                                  ens,
                                  year1,
                                  year2,
                                  dataset_ref,
                                  expid_ref,
                                  ens_ref,
                                  year1_ref,
                                  year2_ref,
                                  season,
                                  FIGDIR,
                                  FILESDIR,
                                  REFDIR,
                                  nclusters) {
  if (nclusters != 4 | season != "DJF") {
    stop("Beta version: unsupported season and/or number of clusters")
  }

  ##########################################################
  #-----------------Loading datasets-----------------------#
  ##########################################################

  # loading anomalies and variances of experiment
  nomefile <-
    file_builder(
      FILESDIR,
      "Regimes",
      "RegimesPattern",
      dataset,
      expid,
      ens,
      year1,
      year2,
      season
    )
  frequencies_exp <- ncdf_opener(nomefile, "Frequencies")
  regimes_exp <-
    ncdf_opener(nomefile, namevar = "Regimes", rotate = "no")

  # loading names
  p <- nc_open(nomefile)
  names_exp <- ncvar_get(p, "Names")
  nc_close(p)
  print(names_exp)

  # loading reference field
  # check for REFDIR==FILESDIR, i.e. if we are using the climatology
  # provided by MiLES or another dataset MiLES-generated
  if (REFDIR != FILESDIR) {
    nomefile_ref <-
      paste0(
        file.path(REFDIR, "Regimes"),
        "/RegimesPattern_",
        dataset_ref,
        "_",
        year1_ref,
        "_",
        year2_ref,
        "_",
        season,
        ".nc"
      )
  } else {
    # use file.builder to create the path of the blocking files
    nomefile_ref <-
      file_builder(
        FILESDIR,
        "Regimes",
        "RegimesPattern",
        dataset_ref,
        expid_ref,
        ens_ref,
        year1_ref,
        year2_ref,
        season
      )
  }

  frequencies_ref <- ncdf_opener(nomefile_ref, "Frequencies")
  regimes_ref <-
    ncdf_opener(nomefile_ref, namevar = "Regimes", rotate = "no")

  # loading names
  p <- nc_open(nomefile_ref)
  names_ref <- ncvar_get(p, "Names")
  nc_close(p)
  print(names_ref)

  # plot properties
  lev_field <- seq(-250, 250, 20)
  lev_diff <- seq(-150, 150, 20)

  # standard properties
  info_exp <-
    info_builder(dataset, expid, ens, year1, year2, season)
  info_ref <- info_builder(
    dataset_ref, expid_ref, ens_ref,
    year1_ref, year2_ref, season
  )

  filenames <- c()
  kk0 <- 1
  # loop on regimes
  for (name in names_ref) {
    #-----plotting-------#
    # a bit complicated but it is used to compare similar regimes
    # even if they do not equal percentage of occurrence (using names)
    ii <- which(name == names_exp)
    jj <- which(name == names_ref)
    print(ii)
    print(jj)
    if (length(ii) == 0) {
      ii <- which(setdiff(names_exp, names_ref)[kk0] == names_exp)
      kk0 <- kk0 + 1
    }
    print(name)

    # final plot production
    figname <- fig_builder(
      FIGDIR,
      "Regimes",
      paste0("Regime", ii),
      dataset,
      expid,
      ens,
      year1,
      year2,
      season,
      output_file_type
    )
    print(figname)
    filenames <- c(filenames, figname)

    # Chose output format for figure - by JvH
    open_plot_device(figname, output_file_type)

    # where to plot frequencies values
    if (map_projection == "no") {
      varpoints <- c(120, 85)
    } else {
      varpoints <- c(0, 0.7)
    }

    # plot properties
    par(plotpar)

    im <- plot_prepare(ics,
      ipsilon,
      regimes_exp[, , ii],
      proj = map_projection,
      lat_lim = lat_lim
    )
    filled_contour3(
      im$x,
      im$y,
      im$z,
      xlab = im$xlab,
      ylab = im$ylab,
      main = paste(info_exp),
      levels = lev_field,
      color.palette = palette3,
      xlim = im$xlim,
      ylim = im$ylim,
      axes = im$axes
    )
    mtext(
      name,
      side = 3,
      line = .5,
      outer = TRUE,
      cex = 2,
      font = 2
    )
    proj_addland(proj = map_projection)
    text(varpoints[1],
      varpoints[2],
      paste("Frequencies: ", round(frequencies_exp[ii], 2), "%", sep = ""),
      cex = 2
    )

    im <- plot_prepare(ics,
      ipsilon,
      regimes_ref[, , jj],
      proj = map_projection,
      lat_lim = lat_lim
    )
    filled_contour3(
      im$x,
      im$y,
      im$z,
      xlab = im$xlab,
      ylab = im$ylab,
      main = paste(info_ref),
      levels = lev_field,
      color.palette = palette3,
      xlim = im$xlim,
      ylim = im$ylim,
      axes = im$axes
    )
    proj_addland(proj = map_projection)
    text(varpoints[1],
      varpoints[2],
      paste("Frequencies: ",
        round(frequencies_ref[ii], 2), "%",
        sep = ""
      ),
      cex = 2
    )
    image_scale3(
      volcano,
      levels = lev_field,
      color.palette = palette3,
      colorbar.label = "m",
      cex.colorbar = imgscl_colorbar,
      cex.label = imgscl_label,
      colorbar.width = 1 * af,
      line.label = imgscl_line
    )

    # delta field plot
    im <-
      plot_prepare(ics,
        ipsilon,
        regimes_exp[, , ii] - regimes_ref[, , jj],
        proj = map_projection,
        lat_lim = lat_lim
      )
    filled_contour3(
      im$x,
      im$y,
      im$z,
      xlab = im$xlab,
      ylab = im$ylab,
      main = paste("Difference"),
      levels = lev_diff,
      color.palette = palette2,
      xlim = im$xlim,
      ylim = im$ylim,
      axes = im$axes
    )
    proj_addland(proj = map_projection)
    image_scale3(
      volcano,
      levels = lev_diff,
      color.palette = palette2,
      colorbar.label = "m",
      cex.colorbar = imgscl_colorbar,
      cex.label = imgscl_label,
      colorbar.width = 1 * af,
      line.label = imgscl_line
    )

    dev.off()
  }
  return(list(figs = filenames, mod = nomefile, ref = nomefile_ref))
}
