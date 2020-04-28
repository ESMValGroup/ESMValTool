######################################################
#------Blocking routines plotting for MiLES----------#
#-------------P. Davini (May 2017)-------------------#
######################################################

miles_block_figures <- function(dataset,
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
                                REFDIR) {
  # which fieds to load/plot
  fieldlist <- c(
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
    "NumberEvents",
    "TM90"
  )

  ##########################################################
  #-----------------Loading datasets-----------------------#
  ##########################################################

  # open field
  for (field in fieldlist) {
    # use file.builder function
    nomefile <- file_builder(
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
    field_exp <-
      ncdf_opener(nomefile, namevar = field, rotate = "no")
    assign(paste(field, "_exp", sep = ""), field_exp)
  }

  # open reference field
  for (field in fieldlist) {
    # check for REFDIR==FILESDIR, i.e. if we are using the climatology
    # provided by MiLES or another dataset MiLES-generated
    if (REFDIR != FILESDIR) {
      nomefile_ref <- paste0(
        file.path(REFDIR, "Block"),
        "/BlockClim_",
        # nolint
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
      nomefile_ref <- file_builder(
        FILESDIR,
        "Block",
        "BlockClim",
        dataset_ref,
        expid_ref,
        ens_ref,
        year1_ref,
        year2_ref,
        season
      )
    }

    field_ref <-
      ncdf_opener(nomefile_ref, namevar = field, rotate = "no")
    assign(paste(field, "_ref", sep = ""), field_ref)
  }

  ##########################################################
  #-----------------Produce figures------------------------#
  ##########################################################

  # standard properties
  info_exp <-
    info_builder(dataset, expid, ens, year1, year2, season)
  info_ref <- info_builder(
    dataset_ref, expid_ref, ens_ref, year1_ref,
    year2_ref, season
  )

  filenames <- c()
  # loop on fields
  for (field in fieldlist) {
    # define field-dependent properties
    fp <- field_details(field)

    # get fields
    field_ref <- get(paste(field, "_ref", sep = ""))
    field_exp <- get(paste(field, "_exp", sep = ""))

    # create figure names with ad-hoc function
    figname <- fig_builder(
      FIGDIR,
      "Block",
      field,
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

    # special treatment for TM90: it is a 1D field!
    if (field == "TM90") {
      open_plot_device(figname, output_file_type, special = TRUE)

      # panels option
      par(
        cex.main = 2,
        cex.axis = 1.5,
        cex.lab = 1.5,
        mar = c(5, 5, 4, 3),
        oma = c(0, 0, 0, 0)
      )

      # rotation to simplify the view (90 deg to the west)
      n <- (-length(ics) / 4)
      ics2 <- c(tail(ics, n), head(ics, -n) + 360)
      field_exp2 <- c(tail(field_exp, n), head(field_exp, -n))
      field_ref2 <- c(tail(field_ref, n), head(field_ref, -n))

      # plot properties
      lwdline <- 4
      tm90cols <- fp$color_field
      plot(
        ics2,
        field_exp2,
        type = "l",
        lwd = lwdline,
        ylim = fp$lev_field,
        main = fp$title_name,
        xlab = "Longitude",
        ylab = fp$legend_unit,
        col = tm90cols[1]
      )
      points(
        ics2,
        field_ref2,
        type = "l",
        lwd = lwdline,
        lty = 1,
        col = tm90cols[2]
      )
      grid()
      legend(
        100,
        30,
        legend = c(info_exp, info_ref),
        lwd = lwdline,
        lty = c(1, 1),
        col = tm90cols,
        bg = "white",
        cex = 1.
      )

      dev.off()

      # skip other part of the script
      next
    }

    # Choose output format for figure - by JvH
    open_plot_device(figname, output_file_type)

    # plot options
    par(plotpar)

    # main experiment plot
    im <- plot_prepare(ics,
      ipsilon,
      field_exp,
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
      levels = fp$lev_field,
      color.palette = fp$color_field,
      xlim = im$xlim,
      ylim = im$ylim,
      axes = im$axes
    )
    mtext(
      fp$title_name,
      side = 3,
      line = .5,
      outer = TRUE,
      cex = 2,
      font = 2
    )
    proj_addland(proj = map_projection)

    # reference field plot
    im <- plot_prepare(ics,
      ipsilon,
      field_ref,
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
      levels = fp$lev_field,
      color.palette = fp$color_field,
      xlim = im$xlim,
      ylim = im$ylim,
      axes = im$axes
    )
    proj_addland(proj = map_projection)
    image_scale3(
      volcano,
      levels = fp$lev_field,
      color.palette = fp$color_field,
      colorbar.label = fp$legend_unit,
      cex.colorbar = imgscl_colorbar,
      cex.label = imgscl_label,
      colorbar.width = 1 * af,
      line.label = fp$legend_distance
    )

    # delta field plot
    im <- plot_prepare(ics,
      ipsilon,
      field_exp - field_ref,
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
      levels = fp$lev_diff,
      color.palette = fp$color_diff,
      xlim = im$xlim,
      ylim = im$ylim,
      axes = im$axes
    )
    proj_addland(proj = map_projection)
    image_scale3(
      volcano,
      levels = fp$lev_diff,
      color.palette = fp$color_diff,
      colorbar.label = fp$legend_unit,
      cex.colorbar = imgscl_colorbar,
      cex.label = imgscl_label,
      colorbar.width = 1 * af,
      line.label = fp$legend_distance
    )

    dev.off()
  }
  return(list(figs = filenames, mod = nomefile, ref = nomefile_ref))
}
