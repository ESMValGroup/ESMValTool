######################################################
#--------Trend plotting routine for HyInt------------#
#-------------E. Arnone (September 2017)-------------#
######################################################

hyint_plot_trends <- function(work_dir, # nolint
                              plot_dir,
                              ref_idx,
                              season,
                              prov_info) {
  #  Define subscripts for variable names
  var_type <- c("tseries", "tseries-sd", "trend", "trend-stat")

  # Set main paths
  dir.create(plot_dir, recursive = T)

  # Load palette
  palette(palette_ts)

  # Number of models
  nmodels <- length(models_name)

  # Define regions to be used
  nregions <- length(selregions)
  if ((plot_type == 13) | (plot_type == 15)) {
    # if plotting multiple models use only first region of list
    nregions <- 1
  }

  #  Check whether enough panels are allocated for plotting the
  # indices requested. In not, drop extra indices
  npanels <- npancol * npanrow
  if (npanels < length(selfields)) {
    selfields <- selfields[1:npanels]
  }

  # Update number of panels and columns if selfields has one element only
  if (length(selfields) == 1) {
    npancol <- 1
    npanrow <- 1
  }

  # Define fields to be used (note that the routine is
  # optimized for 6 fields in 3x2 panels per multi-panel figures)
  if (selfields[1] != F) {
    field_names <- field_names[selfields, drop = F]
    levels_m <- levels_m[selfields, , drop = F]
    tlevels_m <- tlevels_m[selfields, , drop = F]
    title_unit_m <- title_unit_m[selfields, , drop = F]
  }

  # Define field label for filenames
  field_label <- "multiindex"
  if (length(selfields) == 1) {
    field_label <- field_names
  }

  # Remove preset range of values for plotting if needed
  nyears <- models_end_year[ref_idx] - models_start_year[ref_idx]
  if (autolevels) {
    tlevels_m[] <- NA
    levels_m[] <- NA
  }

  # If on, switch shade option off and lines on for plot_type 13
  if (plot_type == 13 & add_trend_sd_shade) {
    add_trend_sd_shade <- F
    add_trend_sd_lines <- T
  }

  # Define array to store plotting limits for each panel of multi-panel figures
  plot_limits <- array(NaN, c(4, length(field_names)))

  # Load parameters for reference dataset
  year1_ref <- models_start_year[ref_idx]
  year2_ref <- models_end_year[ref_idx]

  # Handle label tag when overplotting data from tseries
  # files with different labels in plot_type 14 and 15
  label_figname <- label[1]
  if (length(label) > 1 & plot_type >= 10) {
    label_figname <- paste0(label[1], "-plus")
  }

  # Set figure dimensions
  plot_size <-
    scale_figure(
      plot_type,
      diag_script_cfg,
      length(selfields),
      npancol,
      npanrow
    )

  # Startup graphics for multi-model timeseries or trends
  plot_type_now <- (plot_type == 13) | (plot_type == 15)
  if (plot_type_now == T) {
    tseries_trend_tag <- "timeseries"
    if (plot_type == 15) {
      tseries_trend_tag <- "trend_summary"
    }
    figname <- getfilename_figure(
      plot_dir,
      field_label,
      year1_ref,
      year2_ref,
      ref_idx,
      season,
      "",
      region_codes[selregions[1]],
      label_figname,
      tseries_trend_tag,
      output_file_type,
      multimodel = T
    )
    graphics_startup(figname, output_file_type, plot_size)
    par(
      mfrow = c(npanrow, npancol),
      cex.main = 1.3,
      cex.axis = 1.2,
      cex.lab = 1.2,
      mar = c(5, 5, 5, 5),
      oma = c(1, 1, 1, 1)
    )
  }

  n_noplot <- 1
  if ((plot_type == 13 | plot_type == 15) & autolevels) {
    n_noplot <- 2
  }
  minmax_levels <- c(NA, NA)

  # if requested, loop twice over all models to get range of values for plots
  for (noplot in n_noplot:1) {
    # Loop over models
    for (model_idx in 1:nmodels) {
      # setting up path and parameters
      year1 <- models_start_year[model_idx]
      year2 <- models_end_year[model_idx]

      # Years to be considered based on namelist and cfg_file
      years <- year1:year2
      years_ref <- year1_ref:year2_ref
      if (ryearplot[1] == "ALL") {
        years <- year1:year2
      } else if (ryearplot[1] == "FIRST") {
        years <- year1
      } else {
        years <- years[match(ryearplot, years)]
        years <- years[!is.na(years)]
      }
      if (plot_type >= 14) {
        add_trend <- F
      } # do not plot trend line for plot 14 or 15

      # Startup graphics for multi-region timeseries
      if (plot_type == 12) {
        figname <- getfilename_figure(
          plot_dir,
          field_label,
          year1,
          year2,
          model_idx,
          season,
          "",
          "multiregion",
          label_figname,
          "timeseries",
          output_file_type
        )
        graphics_startup(figname, output_file_type, plot_size)
        par(
          mfrow = c(npanrow, npancol),
          cex.main = 1.3,
          cex.axis = 1.2,
          cex.lab = 1.2,
          mar = c(5, 5, 5, 5),
          oma = c(1, 1, 1, 1)
        )
      }
      #  Startup graphics for bar plot of trend coefficients
      if (plot_type == 14) {
        figname <- getfilename_figure(
          plot_dir,
          field_label,
          year1,
          year2,
          model_idx,
          season,
          "",
          "multiregion",
          label_figname,
          "trend_summary",
          output_file_type
        )
        graphics_startup(figname, output_file_type, plot_size)
        par(
          mfrow = c(npanrow, npancol),
          cex.main = 1.3,
          cex.axis = 1.2,
          cex.lab = 1.2,
          mar = c(8, 8, 2, 2),
          oma = c(1, 1, 1, 1)
        )
      }

      if (model_idx == 1) {
        store_label <- label
      }

      # ----- Loop over label when plotting more files in the same panel ----
      for (ilabel in seq_along(store_label)) {
        label <- store_label[ilabel]
        #-----------------Loading data-----------------------#

        # open timeseries and trends for exp and ref
        infile <-
          getfilename_trends(work_dir, label, model_idx, season)
        print(paste("HyInt_trends: reading file ", infile))
        field_long_names <- array(NaN, length(field_names))
        field_units <- array(NaN, length(field_names))

        if ((plot_type == 13) | (plot_type == 15)) {
          # Store data for provenance
          caption <-
            paste0(
              "Hyint timeseries for selected indices and regions ",
              "according to selected datasets"
            )
          if (plot_type == 15) {
            caption <- paste0(
              "Hyint trends for multiple indices and regions ",
              "according to selected datasets"
            )
          }
          if (length(prov_info[[figname]]) == 0) {
            anc_list <- flatten_lists(prov_info[[infile]]$ancestors)
            prov_fig_now <- list(
              figname = figname,
              caption = caption,
              model_idx = list(model_idx),
              ancestors = anc_list
            )
            prov_info[[figname]] <- prov_fig_now
          } else {
            if (is.na(match(model_idx, prov_info[[figname]]$model_idx))) {
              prov_info[[figname]]$model_idx <-
                c(prov_info[[figname]]$model_idx, model_idx)
              prov_info[[figname]]$ancestors <-
                c(prov_info[[figname]]$ancestors,
                  prov_info[[infile]]$ancestors)
            }
          }
        }

        for (var in field_names) {
          ivar <- which(field_names == var)
          for (stype in var_type[1:2]) {
            svar <- paste0(var, "_", stype)
            rfield <- ncdf_opener(infile, svar, "region", timedimname,
              rotate = "no"
            )
            assign(svar, rfield) # assign field data to field name
            nc <- nc_open(infile)
            dlname <- ncatt_get(nc, svar, "long_name")
            dunits <- ncatt_get(nc, svar, "units")
            field_long_names[ivar] <- dlname$value
            field_units[ivar] <- dunits$value
            nc_close(nc)
          }
          for (stype in var_type[3:4]) {
            svar <- paste0(var, "_", stype)
            rfield <-
              ncdf_opener(infile, svar, "region", "coefficients",
                rotate = "no"
              )
            assign(svar, rfield) # assign field data to field name
          }
        }

        # store size of time and region arrays
        time <- ncdf_opener(infile, timedimname,
          timedimname,
          rotate = "no"
        ) + 1950
        regions <-
          ncdf_opener(infile, "regions", "region", "boundaries",
            rotate = "no"
          )
        # setup time selection for trends
        rettimes <- which(!is.na(time))
        if (trend_years[1] != F) {
          # apply trend to limited time interval if required
          rettimes_tmp <-
            (time >= trend_years[1]) & time <= trend_years[2]
          rettimes <- which(rettimes_tmp)
          if (length(trend_years) == 4) {
            # apply trend also to second time interval if required
            rettime2_tmp <-
              (time >= trend_years[3]) & time <= trend_years[4]
            rettimes2 <- which(rettime2_tmp)
          }
        }
        xlim <- c(min(time), max(time))
        if (trend_years_only & (trend_years[1] != F)) {
          xlim <- trend_years[1:2]
        }


        #-----------------Producing figures------------------------#

        print(paste0(diag_base, ": starting figures"))

        # LOOP over fields
        for (field in field_names) {
          ifield <- which(field == field_names)

          if (noplot == 2 & model_idx == 1) {
            minmax_levels <- c(NA, NA)
            minmax_tlevels <- c(NA, NA)
            assign(paste0(field, "_levels"), minmax_levels)
            assign(paste0(field, "_tlevels"), minmax_tlevels)
          }

          if (anyNA(title_unit_m[ifield, 1:3])) {
            title_unit_m[ifield, 1] <- field_names[ifield]
            title_unit_m[ifield, 2:3] <- field_long_names[ifield]
            title_unit_m[ifield, 4] <- field_units[ifield]
          }

          # TIMESERIES: get timeseries and trends
          tfield_exp <- get(paste0(field, "_", var_type[1]))
          tfield_exp_sd <- get(paste0(field, "_", var_type[2]))
          trend_exp <- get(paste0(field, "_", var_type[3]))
          trend_exp_stat <- get(paste0(field, "_", var_type[4]))

          if (length(dim(tfield_exp)) < 2) {
            # reshape data to matrix if regions has only one element
            tfield_exp <- array(tfield_exp, c(1, length(tfield_exp)))
            tfield_exp_sd <-
              array(tfield_exp_sd, c(1, length(tfield_exp_sd)))
            trend_exp <- array(trend_exp, c(1, length(trend_exp)))
            trend_exp_stat <-
              array(trend_exp_stat, c(1, length(trend_exp_stat)))
          }
          if (plot_type == 13 | plot_type == 15) {
            # get only first region if working on multimodel
            tfield_exp <- tfield_exp[1, , drop = F]
            tfield_exp_sd <- tfield_exp_sd[1, , drop = F]
            trend_exp <- trend_exp[1, , drop = F]
            trend_exp_stat <- trend_exp_stat[1, , drop = F]
          }

          if (is.na(levels_m[ifield, 1]) |
            is.na(levels_m[ifield, 2])) {
            print("No value for range: assigning min and max")
            tmp.levels <- c(
              min(tfield_exp, na.rm = T),
              max(tfield_exp, na.rm = T)
            )
            if (add_trend_sd | add_trend_sd_shade) {
              tmp.levels <- c(
                min(tfield_exp - tfield_exp_sd, na.rm = T),
                max(tfield_exp + tfield_exp_sd, na.rm = T)
              )
            }
          } else {
            tmp.levels <- c(levels_m[ifield, 1], levels_m[ifield, 2])
          }

          if (nyears < 20 & (!autolevels)) {
            levrange <- max(tmp.levels, na.rm = T) - min(tmp.levels, na.rm = T)
            meanrange <- mean(tmp.levels, na.rm = T)
            tmp.levels <- c(
              meanrange - levrange * 1.5,
              meanrange + levrange * 1.5
            )
          }

          #  Startup graphics for one timeseries in one figure
          if (plot_type == 11) {
            figname <- getfilename_figure(
              plot_dir,
              field,
              year1,
              year2,
              model_idx,
              season,
              "",
              "multiregion",
              label_figname,
              "timeseries_single",
              output_file_type
            )
            graphics_startup(figname, output_file_type, plot_size)
            par(
              cex.main = 1.3,
              cex.axis = 1.2,
              cex.lab = 1.2,
              mar = c(4, 4, 2, 2),
              oma = c(1, 1, 1, 1)
            )
          }

          # Actual plotting
          if ((plot_type == 11) |
            (plot_type == 12) | (plot_type == 13)) {
            if (plot_type != 11) {
              # set active panel
              par_row <- (ifield - 1) %/% npancol + 1
              par_col <- (ifield - 1) %% npancol + 1
              par(mfg = c(par_row, par_col, npanrow, npancol))
            }

            # scale autolevels if required
            if (autolevels && (autolevels_scale != 1)) {
              autorange <- max(tmp.levels, nat.rm = T)
              -min(tmp.levels, na.rm = T)
              meanrange <- mean(tmp.levels, na.rm = T)
              tmp.levels <-
                c(
                  meanrange - autorange * autolevels_scale,
                  meanrange + autorange * autolevels_scale
                )
            }

            if (noplot == 2 & autolevels & plot_type == 13) {
              # Recursively store min and max values to be plotted
              # NOTE: this works as long as only one region at the time is used
              minmax_levels <- get(paste0(field, "_levels"))
              minmax_levels[1] <-
                min(c(minmax_levels[1], tmp.levels[1]),
                  na.rm = T
                )
              minmax_levels[2] <-
                max(c(minmax_levels[2], tmp.levels[2]),
                  na.rm = T
                )
              assign(paste0(field, "_levels"), minmax_levels)
              next
            }

            if (noplot == 1 & autolevels & plot_type == 13) {
              tmp.levels[1] <- (get(paste0(field, "_levels")))[1]
              tmp.levels[2] <- (get(paste0(field, "_levels")))[2]
            }

            # Base plot
            if (!(plot_type == 13 & model_idx > 1) & ilabel == 1) {
              ylab <- paste0(title_unit_m[ifield, 1])
              if (title_unit_m[ifield, 4] != "") {
                ylab <- paste0(ylab, title_unit_m[ifield, 4])
              }
              plot(
                time,
                type = "n",
                ylim = c(tmp.levels[1], tmp.levels[2]),
                xlim = xlim,
                xlab = "Year",
                ylab = ylab,
                main = title_unit_m[ifield, 3],
                xaxs = "i"
              )
              # store panel plot limits
              plot_limits[, ifield] <- par("usr")
            }

            # Update plot limits in case panel has changed
            par(usr = plot_limits[, ifield])

            # LOOP over regions to plot timeseries
            if (add_trend_sd_shade) {
              for (ireg in 1:nregions) {
                iselreg <- selregions[ireg]
                shade_area <- c(
                  tfield_exp[ireg, ] + tfield_exp_sd[ireg, ],
                  rev(tfield_exp[ireg, ] - tfield_exp_sd[ireg, ])
                )
                shade_area[shade_area < tmp.levels[1]] <-
                  tmp.levels[1]
                polygon(c(time, rev(time)),
                  shade_area,
                  col = "grey95",
                  border = NA
                )
              }
            }
            for (ireg in 1:nregions) {
              iselreg <- selregions[ireg]
              col_ts <- ireg
              if (length(label) > 1) {
                col_ts <- c(
                  "dodgerblue4",
                  "darkseagreen4",
                  "goldenrod4",
                  "coral4",
                  "grey",
                  "mediumorchid1",
                  "black"
                )[ilabel]
              }
              if (plot_type == 13) {
                col_ts <- model_idx
              }
              if (add_trend_sd) {
                lines(time,
                  tfield_exp[ireg, ] + tfield_exp_sd[ireg, ],
                  lty = 3,
                  col = col_ts
                )
                lines(time,
                  tfield_exp[ireg, ] - tfield_exp_sd[ireg, ],
                  lty = 3,
                  col = col_ts
                )
              }
              if (add_tseries_lines) {
                lines(time, tfield_exp[ireg, ], col = col_ts)
              }
              points(time, tfield_exp[ireg, ], col = col_ts)
              if (add_trend) {
                lines(
                  time[rettimes],
                  trend_exp[ireg, 1] + trend_exp[ireg, 2]
                  * time[rettimes],
                  col = col_ts,
                  lwd = 2
                )
                if (length(trend_years) == 4) {
                  # apply trend also to second time interval if required
                  lines(
                    time[rettimes2],
                    trend_exp[ireg, 3] + trend_exp[ireg, 4]
                    * time[rettimes2],
                    col = col_ts,
                    lwd = 2
                  )
                }
              }
            }
            if (abs(add_legend) & ((plot_type == 11) |
              (plot_type == 12)) & (ifield == 1)) {
              pos_legend <- c(
                plot_limits[1, ifield] + (plot_limits[2, ifield]
                - plot_limits[1, ifield]) * xy_legend[1],
                plot_limits[3, ifield] + (plot_limits[4, ifield]
                - plot_limits[3, ifield]) * xy_legend[2]
              )
              ncol <- 1
              if (add_legend < 0) {
                ncol <- nregions
              }
              if (add_legend > 1) {
                ncol <- add_legend
              }
              legend(
                pos_legend[1],
                pos_legend[2],
                region_codes[selregions],
                text.col = (1:nregions),
                ncol = ncol
              )
            }
            box(lwd = 2)
            if (plot_type == 11) {
              graphics_close(figname)
              # Store data for provenance
              caption <- paste0(
                "Hyint timeseries for index ",
                field,
                " over selected regions according to ",
                models_name[model_idx]
              )
              anc_list <- flatten_lists(prov_info[[infile]]$ancestors)
              prov_fig_now <- list(
                figname = figname,
                caption = caption,
                model_idx = list(model_idx),
                ancestors = anc_list
              )
              prov_info[[figname]] <- prov_fig_now
            }
          }
          if ((plot_type == 14) | (plot_type == 15)) {
            # plot trend coefficients for different regions,
            # one panel per field
            if (anyNA(tlevels_m[ifield, ])) {
              print("No value for range: assigning min and max")
              ylim <- c(
                min(trend_exp_stat[, 1] - trend_exp_stat[, 2], na.rm = T),
                max(trend_exp_stat[, 1] + trend_exp_stat[, 2], na.rm = T)
              )
              # scale autolevels if required
              if (autolevels && (autolevels_scale_t != 1)) {
                autorange <- max(ylim, na.rm = T) - min(ylim, na.rm = T)
                meanrange <- mean(ylim, na.rm = T)
                ylim <- c(
                  meanrange - autorange * autolevels_scale_t,
                  meanrange + autorange * autolevels_scale_t
                )
              }
            } else {
              ylim <- tlevels_m[ifield, ]
            }

            if (trend_years[1] != F) {
              xlim <- trend_years[1:2]
            }
            ylab <- paste0("Avg trend")
            # change y scale to % and 1/100 years
            if (scalepercent & (field != "hyint")) {
              trend_exp <- trend_exp * 100 # trend coefficients
              trend_exp_stat[, 2] <-
                trend_exp_stat[, 2] * 100 # standard error
              ylab <- paste0(ylab, " (%)")
              ylim <- ylim * 100
            }
            if (scale100years) {
              trend_exp <- trend_exp * 100 # trend coefficients
              trend_exp_stat[, 2] <-
                trend_exp_stat[, 2] * 100 # standard error
              ylab <- paste0(ylab, " (1/100 years)")
              ylim <- ylim * 100
            }
            nx <- nregions
            xlab <- "Regions"
            xlabels <- region_codes[selregions]
            if (plot_type == 15) {
              nx <- nmodels
              xlab <- "" # "Models"
              xlabels <- models_name
            }
            # hereafter xregions is the x which also holds models
            # for plottype 15
            xregions <- 1:nx

            # Actual plotting
            # set active panel
            par_row <- (ifield - 1) %/% npancol + 1
            par_col <- (ifield - 1) %% npancol + 1
            par(mfg = c(par_row, par_col, npanrow, npancol))

            if (noplot == 2 & autolevels & plot_type == 15) {
              # Recursively store min and max values to be plotted
              # NOTE: this works as long as only one region at the time is used
              minmax_tlevels <- get(paste0(field, "_tlevels"))
              minmax_tlevels[1] <- min(c(minmax_tlevels[1], ylim[1]),
                na.rm = T
              )
              minmax_tlevels[2] <- max(c(minmax_tlevels[2], ylim[2]),
                na.rm = T
              )
              assign(paste0(field, "_tlevels"), minmax_tlevels)
              next
            }
            if (noplot == 1 & autolevels & plot_type == 15) {
              ylim[1] <- (get(paste0(field, "_tlevels")))[1]
              ylim[2] <- (get(paste0(field, "_tlevels")))[2]
            }

            # Base plot
            if (!(plot_type == 15 & model_idx > 1) & ilabel == 1) {
              plot(
                xregions,
                xregions,
                type = "n",
                pch = 22,
                axes = F,
                xlab = xlab,
                ylab = ylab,
                ylim = ylim,
                main = (
                  paste0(
                    title_unit_m[ifield, 1],
                    " trend (", xlim[1], "-", xlim[2], ")"
                  )
                )
              )
              box()
              # store panel plot limits
              plot_limits[, ifield] <- par("usr")
            }

            # Update plot limits in case panel has changed
            par(usr = plot_limits[, ifield])
            for (ireg in 1:nregions) {
              iregion <- ireg
              ixregion <- ireg
              if (plot_type == 15) {
                ixregion <- model_idx
              }
              # add errorbar (standard error)
              if (!anyNA(trend_exp_stat[iregion, ])) {
                arrows(
                  xregions[ixregion],
                  trend_exp[iregion, 2] -
                    trend_exp_stat[iregion, 2],
                  xregions[ixregion],
                  trend_exp[iregion, 2] + trend_exp_stat[iregion, 2],
                  length = 0.05,
                  angle = 90,
                  code = 3
                )
                points(
                  xregions[ixregion],
                  trend_exp[iregion, 2],
                  pch = 22,
                  col = "grey40",
                  bg = "white",
                  cex = 2
                )
                # add filled points for significant (95% level)
                col90 <- "grey70"
                col95 <- "dodgerblue3"
                if (length(label) > 1) {
                  col90 <- c(
                    "dodgerblue3",
                    "darkseagreen3",
                    "goldenrod3",
                    "coral3",
                    "grey",
                    "mediumorchid1",
                    "black"
                  )
                  col95 <-
                    c(
                      "dodgerblue4",
                      "darkseagreen4",
                      "goldenrod4",
                      "coral4",
                      "grey",
                      "mediumorchid1",
                      "black"
                    )
                }
                if (trend_exp_stat[iregion, 4] <= 0.1) {
                  points(
                    xregions[ixregion],
                    trend_exp[iregion, 2],
                    pch = 22,
                    col = col90[ilabel],
                    bg = col90[ilabel],
                    cex = 2
                  )
                }
                if (trend_exp_stat[iregion, 4] <= 0.05) {
                  points(
                    xregions[ixregion],
                    trend_exp[iregion, 2],
                    pch = 22,
                    col = col95[ilabel],
                    bg = col95[ilabel],
                    cex = 2
                  )
                }
              } else {
                print(paste(
                  "MISSING VALUES in index ",
                  field,
                  ", region ",
                  region_codes[iregion]
                ))
                print(trend_exp_stat[iregion, ])
              }
            }
            if (length(label) > 1) {
              retsig90 <- which(trend_exp_stat[, 4] < 0.1)
              if (!is.na(retsig90[1])) {
                points(
                  xregions[retsig90],
                  trend_exp[retsig90, 2],
                  pch = 22,
                  col = "grey70",
                  bg = "grey70",
                  cex = 2
                )
              }
              retsig95 <- which(trend_exp_stat[, 4] < 0.05)
              if (!is.na(retsig95[1])) {
                points(
                  xregions[retsig95],
                  trend_exp[retsig95, 2],
                  pch = 22,
                  col = "dodgerblue3",
                  bg = "dodgerblue3",
                  cex = 2
                )
              }
            }
            box()
            if (!((plot_type == 15) & (model_idx > 1))) {
              if (add_zeroline & (ylim[1] != 0)) {
                lines(
                  c(-1, nx + 1),
                  c(0, 0),
                  lty = 2,
                  lwd = 1.5,
                  col = "grey40"
                )
              }
              las <- 1
              cex.axis <- 1
              if (plot_type == 15) {
                las <- 2
                cex.axis <- 0.8
              }
              axis(
                1,
                labels = xlabels,
                at = xregions,
                las = las,
                cex.axis = cex.axis
              )
              axis(2)
            }
          } # close if on plot_type 14 and 15
        } # close loop over field
      } # close loop over label
      if ((plot_type == 12) | (plot_type == 14)) {
        graphics_close(figname)
        # Store data for provenance
        caption <-
          paste0(
            "Hyint timeseries for selected indices and regions ",
            "according to ",
            models_name[model_idx]
          )
        if (plot_type == 14) {
          caption <- paste0(
            "Hyint trends for multiple indices and regions ",
            "according to ",
            models_name[model_idx]
          )
        }
        anc_list <- flatten_lists(prov_info[[infile]]$ancestors)
        prov_fig_now <- list(
          figname = figname,
          caption = caption,
          model_idx = list(model_idx),
          ancestors = anc_list
        )
        prov_info[[figname]] <- prov_fig_now
      }
    } # close loop over model
  } # close miniloop over noplot

  # Legend for plot_type 13
  if (abs(add_legend) & (plot_type == 13)) {
    ncol <- 1
    if (add_legend > 1) {
      ncol <- add_legend
    }
    if (add_legend < 0) {
      ncol <- nmodels
    }
    #  for (ifield in 1:length(field_names)) {
    ifield <- 1
    # set active panel
    par_row <- (ifield - 1) %/% npancol + 1
    par_col <- (ifield - 1) %% npancol + 1
    par(
      mfg = c(par_row, par_col, npanrow, npancol),
      usr = plot_limits[, ifield]
    )
    pos_legend <- c(
      plot_limits[1, ifield] + (plot_limits[2, ifield] -
        plot_limits[1, ifield]) * xy_legend[1],
      plot_limits[3, ifield] + (plot_limits[4, ifield] -
        plot_limits[3, ifield]) * xy_legend[2]
    )
    legend_label <- ""
    if (tag_legend[1]) {
      legend_label <- models_name
    }
    if (tag_legend[2]) {
      legend_label <- paste(legend_label,
        models_experiments,
        sep = " "
      )
    }
    if (tag_legend[3]) {
      legend_label <- paste(legend_label,
        models_ensemble,
        sep = " "
      )
    }
    legend(
      pos_legend[1],
      pos_legend[2],
      legend = legend_label,
      text.col = (1:nmodels),
      ncol = ncol,
      cex = 0.9
    )
    print(legend_label)
    print("legend_label")
  }
  if ((plot_type == 13) | (plot_type == 15)) {
    graphics_close(figname)
  }
  return(prov_info)
} # close function
