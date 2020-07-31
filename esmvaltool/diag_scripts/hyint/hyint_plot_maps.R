######################################################
#---------Maps plotting routine for HyInt------------#
#-------------E. Arnone (September 2017)-------------#
######################################################

# DECLARING THE FUNCTION: EXECUTION IS AT THE BOTTOM OF THE SCRIPT

hyint_plot_maps <- # nolint
  function(work_dir,
             plot_dir,
             ref_dir,
             ref_idx,
             season,
             prov_info) {
    # setting up path and parameters
    dataset_ref <- models_name[ref_idx]
    year1_ref <- models_start_year[ref_idx]
    year2_ref <- models_end_year[ref_idx]
    years_ref <- year1_ref:year2_ref

    # set main paths
    work_dir_exp <- work_dir
    plot_dir_exp <- plot_dir
    dir.create(plot_dir_exp, recursive = T)

    # Define fields to be used
    if (selfields[1] != F) {
      field_names <- field_names[selfields, drop = F]
      levels_m <- levels_m[selfields, , drop = F]
      title_unit_m <- title_unit_m[selfields, , drop = F]
    }
    nfields <- length(field_names)

    # Define quantity (exp, ref, exp-ref) to be plotted depending on plot_type
    # 1=exp_only, 2=ref_only, 3=exp/ref/exp-ref
    nquantity <- c(1, 3, 3, 1)

    # Define regions to be used
    nregions <- length(selregions)
    if (nregions > dim(regions)[1]) {
      stop(paste(diag_base, ": requesting regions outside list"))
    }

    if (autolevels) {
      levels_m[] <- NA
    }

    # ------- loading reference data ----------
    # load topography if needed
    if (masksealand) {
      topofile <-
        getfilename_indices(work_dir, diag_base, ref_idx, topo = T)
      gridfile <-
        getfilename_indices(work_dir, diag_base, ref_idx, grid = T)
      if (!file.exists(topofile)) {
        create_landseamask(
          regrid = gridfile,
          loc = run_dir,
          regridded_topo = topofile,
          topo_only = T
        )
      }
      relevation <- ncdf_opener(topofile, "topo", "lon", "lat",
        rotate = "no"
      )
    }
    if (highreselevation) {
      highresel <- get_elevation(elev_range = c(highreselevation, 9000))
    }

    # produce desert areas map if required from reference file
    # (mean annual precipitation <0.5 mm, Giorgi et al. 2014)
    if (removedesert) {
      # reference model
      ref_filename <-
        getfilename_indices(ref_dir, diag_base, ref_idx, season)
      pry <-
        ncdf_opener(ref_filename, "pry", "lon", "lat", rotate = "no")
      retdes <- which(pry < 0.5)
      pry[retdes] <- NA
      # create mask with NAs for deserts and 1's for non-desert
      ref_retdes2D <- apply(pry * 0, c(1, 2), sum) + 1
      ref_retdes3D <-
        replicate(dim(pry)[length(dim(pry))], ref_retdes2D)
    }

    # open reference field
    ref_filename <-
      getfilename_indices(ref_dir, diag_base, ref_idx, season)
    print(paste("Reading reference ", ref_filename))
    for (field in field_names) {
      field_ref <-
        ncdf_opener(ref_filename, field, "lon", "lat", rotate = "no")
      ics_ref <- ics
      ipsilon_ref <- ipsilon

      if (removedesert) {
        field_ref <- field_ref * ref_retdes3D
      }
      if (masksealand) {
        field_ref <- apply_elevation_mask(
          field_ref, relevation,
          sealandelevation
        )
      }
      # if requested calculate multiyear average and store at time=1
      # in this case skip multi-year plot_type 4
      if (rmultiyear_mean) {
        if (plot_type == 4) {
          print("skipping multi-year plot_type 4 with multiyear mean")
          return(0)
        }
        # exclude normalization years from multiyear mean
        retyears <- seq_along(years_ref)
        skipyears <- which(as.logical(match(
          years_ref,
          norm_years[1]:norm_years[2]
        )))
        retyears[skipyears] <- NA
        retyears <- retyears[which(is.finite(retyears))]
        field_ref[, , 1] <- apply(field_ref[, , retyears],
          c(1, 2), mean,
          na.rm = T
        )
      }
      assign(paste(field, "_ref", sep = ""), field_ref)
    }

    # Loop over models
    for (model_idx in c(1:(length(models_name)))) {
      # Do not compare reference with itself
      if ((model_idx == ref_idx) &&
        ((plot_type == 2) || (plot_type == 3))) {
        if (length(models_name) == 1) {
          print("skipping comparison plots because
                only one dataset was requested")
        }
        next
      }

      # setting up path and parameters
      exp <- models_name[model_idx]
      year1 <- models_start_year[model_idx]
      year2 <- models_end_year[model_idx]

      # Years to be considered based on namelist and cfg_file
      years <- year1:year2
      if (ryearplot[1] == "ALL") {
        years <- year1:year2
      } else if (ryearplot[1] == "FIRST") {
        years <- year1
      } else {
        years <- years[match(ryearplot, years)]
        years <- years[!is.na(years)]
      }
      nyears <- length(years)

      # Remove deserts if required
      if (removedesert) {
        filename <- getfilename_indices(
          work_dir_exp, diag_base, model_idx,
          season
        )
        if ((rgrid == F) & ((plot_type == 2) | (plot_type == 3))) {
          # regrid when comparing
          pry <-
            ncdf_opener(
              filename,
              "pry",
              "lon",
              "lat",
              rotate = "no",
              interp2grid = T,
              grid = ref_filename
            )
        } else {
          pry <- ncdf_opener(filename, "pry", "lon", "lat", rotate = "no")
        }
        retdes <- which(pry < 0.5)
        pry[retdes] <- NA
        # create mask with NAs for deserts and 1's for non-desert
        exp_retdes2D <- apply(pry * 0, c(1, 2), sum) + 1
        exp_retdes3D <-
          replicate(dim(pry)[length(dim(pry))], exp_retdes2D)
      }

      #-----------------Loading data-----------------------#
      # open experiment field
      for (field in field_names) {
        infile <- getfilename_indices(
          work_dir_exp, diag_base, model_idx,
          season
        )
        print(paste("Reading ", field, " from experiment ", infile))
        if ((rgrid == F) & ((plot_type == 2) | (plot_type == 3))) {
          # regrid when comparing
          field_exp <-
            ncdf_opener(
              infile,
              field,
              "lon",
              "lat",
              rotate = "no",
              interp2grid = T,
              grid = ref_filename
            )
        } else {
          field_exp <- ncdf_opener(infile, field, "lon", "lat", rotate = "no")
        }
        if (removedesert) {
          field_exp <- field_exp * exp_retdes3D
        }
        if (masksealand) {
          field_exp <- apply_elevation_mask(
            field_exp, relevation,
            sealandelevation
          )
        }
        # if requested calculate multiyear average and store it at time=1
        if (rmultiyear_mean) {
          years <- year1:year2
          retyears <- seq_along(years)
          skipyears <- which(as.logical(match(
            years,
            norm_years[1]:norm_years[2]
          )))
          retyears[skipyears] <- NA
          retyears <- retyears[which(is.finite(retyears))]
          field_exp[, , 1] <- apply(field_exp[, , retyears],
            c(1, 2), mean,
            na.rm = T
          )
        }
        if (highreselevation_only) {
          field_exp[] <- NA
        }
        assign(paste(field, "_exp", sep = ""), field_exp)
      }

      #---------------Multiyear mean-----#
      if (rmultiyear_mean) {
        nyears <- 1
      }

      #-----------------Producing figures------------------------#

      print(paste0(diag_base, ": starting figures"))

      # Set figure dimensions
      plot_size <-
        scale_figure(
          plot_type,
          diag_script_cfg,
          length(selfields),
          npancol,
          npanrow
        )
      if (boxregion != 0) {
        # boxregion will plot region boxes over a global map of selected field
        nregions <- 1
      }

      # LOOP over selected regions
      for (iselregion in 1:nregions) {
        iregion <- selregions[iselregion]
        print(paste("region: ", region_names[iregion]))

        # Startup graphics for multiple years in one figure
        if (plot_type == 4) {
          field_label <- "multiindex"
          figname <- getfilename_figure(
            plot_dir_exp,
            field_label,
            year1,
            year2,
            model_idx,
            season,
            "multiyear",
            region_codes[iregion],
            label,
            "map",
            output_file_type
          )
          graphics_startup(figname, output_file_type, plot_size)
          par(
            mfrow = c(nyears, nfields),
            cex.main = 1.3,
            cex.axis = 1.2,
            cex.lab = 1.2,
            mar = c(2, 2, 2, 2),
            oma = c(1, 1, 1, 1)
          )
        }
        # LOOP over years defined in parameter file
        for (iyear in c(1:nyears)) {
          if (ryearplot_ref[1] == "EXP") {
            iyear_ref <- iyear
          } else {
            iyear_ref <- match(ryearplot_ref, years_ref)
          }
          time_label <- years[iyear]
          time_label_ref <- years[iyear_ref]
          time_label_fig <- time_label
          if (rmultiyear_mean) {
            time_label <- paste(year1, year2, sep = "-")
            time_label_ref <- paste(year1_ref, year2_ref, sep = "-")
            time_label_fig <- "myearmean"
          }
          print(paste0(
            diag_base,
            ": plotting data for  ",
            region_names[iregion],
            "-",
            time_label
          ))

          # standard properties
          info_exp <- paste(exp, time_label) # ,season)
          info_ref <- paste(dataset_ref, time_label_ref) # ,season)

          #  Startup graphics for multiple fields/quantities in one figure
          if (plot_type == 3) {
            field_label <- "multiindex"
            figname <- getfilename_figure(
              plot_dir_exp,
              field_label,
              year1,
              year2,
              model_idx,
              season,
              time_label_fig,
              region_codes[iregion],
              label,
              "map",
              output_file_type
            )
            graphics_startup(figname, output_file_type, plot_size)
            par(
              mfrow = c(nfields, 3),
              cex.main = 1.3,
              cex.axis = 1.2,
              cex.lab = 1.2,
              mar = c(2, 2, 2, 6),
              oma = c(1, 1, 1, 1)
            )
          }
          # LOOP over fields
          for (field in field_names) {
            ifield <- which(field == field_names)
            if (anyNA(title_unit_m[ifield, 1:3])) {
              title_unit_m[ifield, 1:3] <- field
              title_unit_m[ifield, 4] <- ""
            }

            # get fields
            field_ref <- get(paste(field, "_ref", sep = ""))
            field_exp <- get(paste(field, "_exp", sep = ""))

            # MAPS: select required year (if requested, multiyear average
            #       is stored at iyear=1)
            field_ref <- field_ref[, , iyear]
            field_exp <- field_exp[, , iyear_ref]
            tmp_field <- field_exp

            # define quantity-dependent properties (exp, ref, exp-ref)
            tmp.colorbar <- c(F, T, T)
            if (plot_type == 1) {
              tmp.colorbar <- T
            }
            tmp.palette <- palette_giorgi2011
            if (is.na(levels_m[ifield, 1]) |
              is.na(levels_m[ifield, 2])) {
              print("No value for range: assigning min and max")
              tmp.levels <- seq(min(field_ref, na.rm = T),
                max(field_ref, na.rm = T),
                len = nlev
              )
            } else {
              tmp.levels <- seq(levels_m[ifield, 1], levels_m[ifield, 2],
                len = nlev
              )
            }
            if (highreselevation_only) {
              title_unit_m[ifield, 1] <- "Elevation"
            }
            tmp.titles <- paste0(
              title_unit_m[ifield, 1],
              ": ",
              region_names[iregion],
              "-",
              c(info_exp, info_ref, "Difference")
            )
            if (plot_type == 4) {
              tmp.titles <- paste(title_unit_m[ifield, 1], time_label)
            }

            # Startup graphics for individual fields and multi
            # quantities in each figure
            if (plot_type == 2) {
              figname <- getfilename_figure(
                plot_dir_exp,
                field,
                year1,
                year2,
                model_idx,
                season,
                time_label_fig,
                region_codes[iregion],
                label,
                "comp_map",
                output_file_type
              )
              graphics_startup(figname, output_file_type, plot_size)
              par(
                mfrow = c(3, 1),
                cex.main = 2,
                cex.axis = 1.5,
                cex.lab = 1.5,
                mar = c(5, 5, 4, 8),
                oma = c(1, 1, 1, 1)
              )
            }

            # --- MAPS ----
            # LOOP over quantity (exp,ref,exp-ref difference) to be plotted
            for (iquantity in c(1:nquantity[plot_type])) {
              if (iquantity == 2) {
                tmp_field <- field_ref
                ipsilon <- ipsilon_ref
                ics <- ics_ref
              }
              if (iquantity == 3) {
                tmp.palette <- palette2
                tmp_field <- field_exp - field_ref
                if (is.na(levels_m[ifield, 3]) |
                  is.na(levels_m[ifield, 4])) {
                  tmp_field_max <- max(abs(tmp_field), na.rm = T)
                  tmp.levels <- seq(-tmp_field_max, tmp_field_max,
                    len = nlev
                  )
                } else {
                  tmp.levels <- seq(levels_m[ifield, 3], levels_m[ifield, 4],
                    len = nlev
                  )
                }
              }
              # Startup graphics for individual field in each figure
              if (plot_type == 1) {
                figname <- getfilename_figure(
                  plot_dir_exp,
                  field,
                  year1,
                  year2,
                  model_idx,
                  season,
                  time_label_fig,
                  region_codes[iregion],
                  label,
                  "map",
                  output_file_type
                )
                graphics_startup(figname, output_file_type, plot_size)
                lonlat_aratio <- (max(ics) - min(ics)) /
                  (max(ipsilon) - min(ipsilon))
                par(
                  mfrow = c(1, 1),
                  cex.main = 2,
                  cex.axis = 1.5,
                  cex.lab = 1.5,
                  mar = c(5, 5, 4, 8),
                  oma = c(1, 1, 1, 1)
                )
                # mar = c(3, 3, 4, 8), oma = c(1, 1, 1, 1))
              }

              # set active panel
              if (plot_type == 3) {
                par(mfg = c(ifield, iquantity, nfields, 3))
              }
              if (plot_type == 4) {
                par(mfg = c(iyear, ifield, nyears, nfields))
              }

              # scale autolevels if required
              if (autolevels && (autolevels_scale != 1)) {
                autorange <- max(tmp.levels) - min(tmp.levels)
                meanrange <- mean(tmp.levels)
                tmp.levels <-
                  seq(
                    meanrange - autorange * autolevels_scale,
                    meanrange + autorange * autolevels_scale,
                    len = nlev
                  )
              }

              cex_main <- 1.4
              if (plot_type == 1) {
                cex_main <- 1.3
              }

              # drop data outside region limits
              retlon <- which(ics < regions[iregion, 1]
              | ics > regions[iregion, 2])
              retlat <- which(ipsilon < regions[iregion, 3]
              | ipsilon > regions[iregion, 4])
              mask_field <- tmp_field
              mask_field[retlon, ] <- NA
              mask_field[, retlat] <- NA
              tmp_field <- mask_field

              # contours
              filled_contour3(
                ics,
                ipsilon,
                tmp_field,
                xlab = "Longitude",
                ylab = "Latitude",
                main = tmp.titles[iquantity],
                levels = tmp.levels,
                color.palette = tmp.palette,
                xlim = c(regions[iregion, 1], regions[iregion, 2]),
                ylim = c(regions[iregion, 3], regions[iregion, 4]),
                axes = F,
                asp = 1,
                cex.main = cex_main
              )
              # continents
              continents_col <- "white"
              if (map_continents <= 0) {
                continents_col <- "gray30"
              }
              map(
                "world",
                regions = ".",
                interior = map_continents_regions,
                exact = F,
                boundary = T,
                add = T,
                col = continents_col,
                lwd = abs(map_continents)
              )
              # rect(regions[iregion, 1], regions[iregion, 3],
              #     regions[iregion, 2], regions[iregion, 4],
              #     border = "grey90", lwd = 3)
              # grid points
              if (oplot_grid) {
                # build up grid if needed
                ics2 <- replicate(length(ipsilon), ics)
                ipsilon2 <- t(replicate(length(ics), ipsilon))
                points(
                  ics2,
                  ipsilon2,
                  pch = 1,
                  col = "grey40",
                  cex = oplot_grid
                )
              }
              # add highres elevation contours
              if (highreselevation) {
                palette(terrain.colors(10))
                contour(
                  highresel$lon_el,
                  highresel$lat_el,
                  highresel$elevation,
                  levels = seq(500, 5000, length.out = 10),
                  col = 1:10,
                  add = T
                )
              }
              # boxes
              box(col = "grey60")
              if (boxregion != 0) {
                box_col <- "white"
                if (boxregion <= 0) {
                  box_col <- "grey30"
                }
                for (ireg in 2:length(selregions)) {
                  iselreg <- selregions[ireg]
                  rect(
                    regions[iselreg, 1],
                    regions[iselreg, 3],
                    regions[iselreg, 2],
                    regions[iselreg, 4],
                    border = box_col,
                    lwd = abs(boxregion)
                  )
                  text(
                    regions[iselreg, 1],
                    regions[iselreg, 3],
                    paste0("         ", region_codes[iselreg]),
                    col = box_col,
                    pos = 3,
                    offset = 0.5
                  )
                }
              }
              # axis
              if (plot_type <= 2) {
                if ((regions[iregion, 2] - regions[iregion, 1] > 90)
                |
                  (regions[iregion, 4] - regions[iregion, 3] > 90)) {
                  axis(1,
                    col = "grey40",
                    at = seq(-180, 180, 45)
                  )
                  axis(2, col = "grey40", at = seq(-90, 90, 30))
                } else {
                  axis(1, col = "grey40")
                  axis(2, col = "grey40")
                }
              } else if (plot_type == 3) {
                if (iquantity == 1) {
                  if ((regions[iregion, 2] - regions[iregion, 1] > 90)
                  |
                    (regions[iregion, 4] - regions[iregion, 3] > 90)) {
                    axis(2,
                      col = "grey40",
                      at = seq(-90, 90, 30)
                    )
                  } else {
                    axis(2, col = "grey40")
                  }
                }
                if (ifield == length(field_names)) {
                  if ((regions[iregion, 2] - regions[iregion, 1] > 90)
                  |
                    (regions[iregion, 4] - regions[iregion, 3] > 90)) {
                    axis(1,
                      col = "grey40",
                      at = seq(-180, 180, 45)
                    )
                  } else {
                    axis(1, col = "grey40")
                  }
                }
              } else if (plot_type == 4) {
                if (iyear == nyears) {
                  if ((regions[iregion, 2] - regions[iregion, 1] > 90)
                  |
                    (regions[iregion, 4] - regions[iregion, 3] > 90)) {
                    axis(1,
                      col = "grey40",
                      at = seq(-180, 180, 45)
                    )
                  } else {
                    axis(1, col = "grey40")
                  }
                }
                if (field == "int_norm") {
                  if ((regions[iregion, 2] - regions[iregion, 1] > 90)
                  |
                    (regions[iregion, 4] - regions[iregion, 3] > 90)) {
                    axis(2,
                      col = "grey40",
                      at = seq(-90, 90, 30)
                    )
                  } else {
                    axis(2, col = "grey40")
                  }
                }
              }

              # colorbar
              new_fig_scale <- c(-0.11, -0.04, 0.1, -0.1)
              line_label <- 2.7
              cex_label <- 1.2
              cex_colorbar <- 1
              if (plot_type == 2) {
                new_fig_scale <- c(-0.07, -0.02, 0.1, -0.1)
                line_label <- 2.7
                cex_label <- 1
                cex_colorbar <- 1.5
              }
              if (plot_type == 3) {
                new_fig_scale <- c(-0.11, -0.03, 0.1, -0.1)
                line_label <- 3
                cex_label <- 1
                cex_colorbar <- 1.2
              }
              if ((tmp.colorbar[iquantity]) & add_colorbar) {
                image_scale3(
                  volcano,
                  levels = tmp.levels,
                  new_fig_scale = new_fig_scale,
                  color.palette = tmp.palette,
                  colorbar.label =
                    paste(
                      title_unit_m[ifield, 1],
                      title_unit_m[ifield, 4]
                    ),
                  cex.colorbar = cex_colorbar,
                  cex.label = cex_label,
                  colorbar.width = 1,
                  line.label = line_label,
                  line.colorbar = 1.0
                )
              }
            } # close loop over quantity
            if (plot_type == 1) {
              graphics_close(figname)
              # Store data for provenance
              caption <-
                paste0(
                  "Map for index ",
                  field,
                  " over region ",
                  region_codes[iregion],
                  " according to ",
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
            if (plot_type == 2) {
              graphics_close(figname)
              # Store data for provenance
              caption <-
                paste0(
                  "Map for index ",
                  field,
                  " over region ",
                  region_codes[iregion],
                  " according to ",
                  models_name[model_idx],
                  " in comparison to reference dataset"
                )
              anc_list <- flatten_lists(c(prov_info[[infile]]$ancestors,
                  prov_info[[ref_filename]]$ancestors))
              prov_fig_now <- list(
                figname = figname,
                caption = caption,
                model_idx = list(model_idx, ref_idx),
                ancestors = anc_list
              )
              prov_info[[figname]] <- prov_fig_now
            }
          } # close loop over field
          if (plot_type == 3) {
            graphics_close(figname)
            # Store data for provenance
            caption <- paste0(
              "Comparison maps for multiple indices",
              " over region ",
              region_codes[iregion]
            )

            anc_list <- flatten_lists(c(prov_info[[infile]]$ancesstors,
                prov_info[[ref_filename]]$ancestors))
            prov_fig_now <- list(
              figname = figname,
              caption = caption,
              model_idx = list(model_idx, ref_idx),
              ancestors = anc_list
            )
            prov_info[[figname]] <- prov_fig_now
          }
        } # close loop over years
        if (plot_type == 4) {
          graphics_close(figname)
          # Store data for provenance
          caption <-
            paste0("Maps for multiple indices over selected years")
          anc_list <- flatten_lists(prov_info[[infile]]$ancestors)
          prov_fig_now <- list(
            figname = figname,
            caption = caption,
            model_idx = list(model_idx),
            ancestors = anc_list
          )
          prov_info[[figname]] <- prov_fig_now
        }
      } # close loop over regions
    } # close loop over models
    return(prov_info)
  } # close function
