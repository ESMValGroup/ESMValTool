; CLOUDS
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates annual/seasonal means of 2-d (cloud) parameters for comparison
;   with a reference data set. Optionally, differences to the reference data
;   set are also plotted.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   embracesetup:       True = 2 plots per line, False = 4 plots per line
;                       (default)
;   explicit_cn_levels: explicit contour levels (array)
;   extralegend:        plot legend(s) to extra file(s)
;   filename_add:       optionally add this string to plot filesnames
;   panel_labels:       label individual panels (true, false)
;   PanelTop:           manual override for "@gnsPanelTop" used by panel
;                       plot(s)
;   projection:         map projection for plotting (default =
;                       "CylindricalEquidistant")
;   showdiff            calculate and plot differences (default = False)
;   rel_diff:           if showdiff = True, then plot relative differences (%)
;                       (default = False)
;   ref_diff_min:       lower cutoff value in case of calculating relative
;                       differences
;                       (in units of input variable)
;   region:             show only selected geographic region given as latmin,
;                       latmax, lonmin, lonmax
;   timemean:           time averaging - "seasonal" = DJF, MAM, JJA, SON),
;                                        "annual" = annual mean
;   treat_var_as_error: treat variable as error when averaging (true, false)
;                       true:  avg = sqrt(mean(var*var))
;                       false: avg = mean(var)
;
; Required variable attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   long_name:         variable description
;   reference_dataset: reference dataset; REQUIRED when calculating
;                      differences (showdiff = True)
;   units:             variable units (for labeling plot only)
;
; Caveats
;   none
;
; Modification history
;   20190220-lauer_axel: added output of provenance (v2.0)
;   20181119-lauer_axel: adapted code to multi-variable capable framework
;   20180923-lauer_axel: added writing of results to netcdf
;   20180518-lauer_axel: code rewritten for ESMValTool v2.0
;   20170621-lauer_axel: reworked code to add tags for reporting
;   20160901-lauer_axel: added regridding option 1 deg x 1 deg
;   20151027-lauer_axel: moved call to 'write_references' to the beginning
;                        of the code
;   20150415-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "embrace_setup", False)
  set_default_att(diag_script_info, "extralegend", False)
  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "panel_labels", True)
  set_default_att(diag_script_info, "rel_diff", False)
  set_default_att(diag_script_info, "rel_diff_min", -1.0e19)
  set_default_att(diag_script_info, "showdiff", False)
  set_default_att(diag_script_info, "timemean", "annualclim")
  set_default_att(diag_script_info, "treat_var_as_error", False)

  flag_diff = diag_script_info@showdiff
  flag_rel_diff = diag_script_info@rel_diff
  flag_rel_diff_min = diag_script_info@rel_diff_min

  if (.not.flag_diff .and. flag_rel_diff) then
    log_info("flag_rel_diff = True has no effect until flag_diff is also " \
             + "set to True")
  end if

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  embracesetup = diag_script_info@embrace_setup

  if (isatt(diag_script_info, "projection")) then
    projection = diag_script_info@projection
    perim = False
  else
    projection = "CylindricalEquidistant"
    perim = True
  end if

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  timemean = diag_script_info@timemean
  numseas = 1          ; default
  season = (/"annual"/)

  if (timemean.eq."seasonalclim") then
    numseas = 4
    delete(season)
    season = (/"DJF", "MAM", "JJA", "SON"/)
  end if

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  panel_labels = diag_script_info@panel_labels

  treat_var_as_error = diag_script_info@treat_var_as_error

  extralegend = diag_script_info@extralegend

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  if (config_user_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if

  ref_ind = -1  ; set to invalid value

  ; if attribute is present, use it so correlations can be calculated
  if (isvar("refname")) then
    ; set reference model
    ref_ind = ind(names .eq. refname)
    if (ismissing(ref_ind)) then
      log_info("warning: reference dataset (" + refname + ") not found.")
      ref_ind = -1
    end if
  end if

  climofiles = metadata_att_as_array(info0, "filename")

  outfile = new(numseas, string)
  outfile(:) = ""

  if (flag_diff) then
    outfile_d = new(numseas, string)
    outfile_d(:) = ""

    ; check for reference model definition
    if (.not.isvar("refname")) then
      error_msg("f", DIAG_SCRIPT, "", \
                "no reference dataset defined in recipe")
    end if

    ; set reference model

    ref_ind = ind(names .eq. refname)
    if (ismissing(ref_ind)) then
      error_msg("f", DIAG_SCRIPT, "", "reference dataset (" \
                + refname + ") is missing")
    end if
  end if

end

begin
  ; ###########################################
  ; # get data and average time               #
  ; ###########################################

  maps = new((/dim_MOD, 4/), graphic)
  maps_d = new((/dim_MOD, 4/), graphic)

  ind_all_sorted = ispan(0, dim_MOD - 1, 1)  ; create array

  if (ref_ind .ge. 0) then
    ind_wo_ref = ind(names .ne. refname)
    ind_all_sorted(0) = ref_ind
    ind_all_sorted(1:dim_MOD - 1) = ind_wo_ref
  end if

  corr = new((/numseas/), float)
  gavg = new((/numseas/), float)
  rmsd = new((/numseas/), float)
  bias = new((/numseas/), float)

  ; filenames for netcdf output

  nc_filename_bias = work_dir + "clouds_" + var0 + "_bias.nc"
  nc_filename_bias@existing = "append"
  nc_filename_mean = work_dir + "clouds_" + var0 + "_mean.nc"
  nc_filename_mean@existing = "append"

  do ii = 0, dim_MOD - 1

    imod = ind_all_sorted(ii)
    log_info("processing " + names(imod))

    if (isvar("data1")) then
      delete(data1)
    end if

    if (isvar("A0")) then
      delete(A0)
    end if

    A0 = read_data(info0[imod])

    ; check dimensions

    dims = getvardims(A0)
    if (dimsizes(dims) .lt. 2) then
      error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + \
                " dimensions, need 2 or 3")
    end if
    idx = ind(dims .eq. "lat")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
    end if
    idx = ind(dims .eq. "lon")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
    end if

    ; average over time

    ; if variable is an error variable, we have to square it before
    ; averaging and then calculate the square-root afterwards

    if (treat_var_as_error) then
      log_info(" ++++++++++++++ Treating variable as error " + \
               "variable when averaging ")
      A0 = A0 * A0
    end if

    data1 = time_operations(A0, -1, -1, "average", timemean, True)

    if (treat_var_as_error) then
      data1 = sqrt(data1)
    end if

    delete(A0)

    ; if requested, select geographical region

    if (isatt(diag_script_info, "region")) then
      region = diag_script_info@region
      data1 := area_operations(data1, region(0), region(1), region(2), \
                               region(3), "extract", False)
      if (region(2).eq.0. .and. region(3).eq.360.) then
      else
        data1@res_gsnAddCyclic = False
      end if
      data1@res_mpMinLatF = region(0)    ; range to zoom in on
      data1@res_mpMaxLatF = region(1)
      data1@res_mpMinLonF = region(2)
      data1@res_mpMaxLonF = region(3)
      data1@res_mpCenterLonF = 0.5 * (region(2) + region(3))
      delete(region)
    end if

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/shared/plot/style.ncl

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    data1@res_cnFillOn       = True      ; color plot desired
    data1@res_cnLineLabelsOn = False     ; contour lines

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    ; annotation

    data1@res_tiMainOn             = False
    data1@res_cnLevelSelectionMode = "ExplicitLevels"
    data1@res_cnLinesOn            = False

    data1@res_mpOutlineOn          = True
    data1@res_mpFillOn             = False

    ; variable specific plotting settings

    if (var0.eq."pr") then
      data1@res_cnLevels = fspan(0.5, 10, 20)
      ; convert from kg m-2 s-1 to mm day-1
      data1 = data1 * 86400.0
      data1@units = "mm day-1"
    end if

    if (var0.eq."lwp") then
      data1@res_cnLevels            = ispan(10, 200, 10) * 0.001
      data1@res_mpOutlineOn         = False
      data1@res_mpFillOn            = True
      data1@res_mpLandFillColor     = "Black"
      pal = read_colormap_file("$diag_scripts/shared/plot/rgb/qcm3.rgb")
      data1@res_cnFillColors        = pal
    end if

    if (var0.eq."tas") then
      data1@res_cnLevels            = ispan(-30, 30, 3)
      pal = read_colormap_file("$diag_scripts/shared/plot/rgb/ipcc-tas.rgb")
      data1@res_cnFillColors        = pal
      ; convert from K to degC
      data1 = data1 - 273.15
      data1@units = "degC"
    end if

    if (var0.eq."clt") then
      data1@res_cnLevels            = fspan(5, 100, 20)
    end if

    if (var0.eq."clivi") then
      data1@res_cnLevels            = ispan(10, 200, 10) * 0.001
    end if

    if (var0.eq."clwvi") then
      data1@res_cnLevels            = ispan(10, 300, 10) * 0.001
    end if

    if (var0.eq."swcre") then
      data1@res_cnLevels            = ispan(-100, 0, 10)
    end if

    if (var0.eq."lwcre") then
      data1@res_cnLevels            = ispan(0, 100, 10)
    end if

    if (var0.eq."netcre") then
      data1@res_cnLevels            = ispan(-70, 70, 10)
    end if

    data1@res_lbLabelBarOn         = False
    data1@res_gsnRightString       = ""

    data1@res_mpFillDrawOrder       = "PostDraw"    ; draw map last
    data1@res_cnMissingValFillColor = "Gray"

    ; no tickmarks and no labels

    data1@res_tmYLLabelsOn       = False
    data1@res_tmYLOn             = False
    data1@res_tmYRLabelsOn       = False
    data1@res_tmYROn             = False
    data1@res_tmXBLabelsOn       = False
    data1@res_tmXBOn             = False
    data1@res_tmXTLabelsOn       = False
    data1@res_tmXTOn             = False
    data1@res_cnInfoLabelOn      = False    ; turn off cn info label
    data1@res_mpPerimOn          = perim    ; draw line around map

    ; specified in namelist

    data1@res_mpProjection       = projection

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      data1@res_cnLevelSelectionMode = "ExplicitLevels"
      data1@res_cnLevels = diag_script_info@explicit_cn_levels
    end if

    if (.not. isatt(data1, "res_cnLevels")) then
      log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
      log_info("info: using default contour levels")
      data1@res_cnLevels = fspan(min(data1), max(data1), 20)
    end if

    ; ###########################################
    ; # other Metadata: diag_script, var        #
    ; ###########################################
    ; add to data1 as attributes without prefix

    if (isatt(data1, "diag_script")) then  ; add to existing entries
      temp = data1@diag_script
      delete(data1@diag_script)
      data1@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
      delete(temp)
    else  ; add as new attribute
      data1@diag_script = (/DIAG_SCRIPT/)
    end if

    if (isatt(variable_info[0], "long_name")) then
      data1@var_long_name = variable_info[0]@long_name
    end if

    data1@var = var0

    if (isatt(variable_info[0], "units")) then
      data1@var_units = variable_info[0]@units
    else
      data1@var_units = ""
    end if

    if (.not. isvar("ref_data")) then
      ref_data = data1
    end if

    ; check if data are on same grid (for calculating difference, RMSD,
    ; correlation)

    same_grid = False

    if (all(dimsizes(ref_data) .eq. dimsizes(data1))) then
      if (max(abs(ref_data&lat - data1&lat)) .le. 1.0e-6) then
        if (max(abs(ref_data&lon - data1&lon)) .le. 1.0e-6) then
          same_grid = True
        end if
      end if
    end if

    if (flag_diff .and. .not.same_grid) then
      flag_diff = False
      error_msg("f", DIAG_SCRIPT, "", \
                "Data are not on same grid, cannot calculate differences. " \
                + "Set showdiff to False in namelist or regrid data to " \
                + "common grid (check/adjust " \
                + "preprocessor settings in recipe).")
    end if

    corr = corr@_FillValue
    gavg = gavg@_FillValue

    if (.not.all(ismissing(data1))) then
      if (numseas.gt.1) then
        do is = 0, numseas - 1
          if (same_grid .and. (ref_ind .ge. 0)) then
            corr(is) = calculate_metric(ref_data(is, :, :), data1(is, :, :), \
                                        "correlation")
          end if
          gavg(is) = area_operations(data1(is, :, :), -90., 90., 0., 360., \
                                     "average", True)
        end do
      else
        if (same_grid .and. (ref_ind .ge. 0)) then
          corr(0) = calculate_metric(ref_data, data1, "correlation")
        end if
        gavg(0) = area_operations(data1, -90., 90., 0., 360., "average", True)
      end if
    end if

    data1@res_gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
    data1@res_gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    data1@res_gsnDraw        = False  ; do not draw yet
    data1@res_gsnFrame       = False  ; don't advance frame

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      ; note: an array of workspaces (i.e. wks(numseas)) does not work as
      ;       attributes cannot be assigned to each array element
      ;       individually
      wks0 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                     "_" + season(0) + filename_add)
      ; difference plots will be saved to a different file
      if (flag_diff) then
        wks0d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                        "_bias_" + season(0) + filename_add)
      end if
      if (numseas.gt.1) then
        wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                       "_" + season(1) + filename_add)
        wks2 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                       "_" + season(2) + filename_add)
        wks3 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                       "_" + season(3) + filename_add)
        ; difference plots will be saved to a different files
        if (flag_diff) then
          wks1d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                          "_bias_" + season(1) + filename_add)
          wks2d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                          "_bias_" + season(2) + filename_add)
          wks3d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                          "_bias_" + season(3) + filename_add)
        end if
      end if
    end if

    if (numseas.gt.1) then
      do is = 0, numseas - 1
        if (.not.ismissing(corr(is))) then
          data1@res_gsnRightString = "corr = " + sprintf("%6.3f", corr(is))
        else
          data1@res_gsnRightString = ""
        end if
        if (.not.ismissing(gavg(is))) then
          data1@res_gsnLeftString = "mean = " + sprintf("%6.3f", gavg(is))
        else
          data1@res_gsnLeftString = ""
        end if

        if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
          data1@res_gsnRightString = ""
        end if

        if (is.eq.0) then
          maps(imod, is) = contour_map(wks0, data1(is, :, :), var0)
        end if
        if (is.eq.1) then
          maps(imod, is) = contour_map(wks1, data1(is, :, :), var0)
        end if
        if (is.eq.2) then
          maps(imod, is) = contour_map(wks2, data1(is, :, :), var0)
        end if
        if (is.eq.3) then
          maps(imod, is) = contour_map(wks3, data1(is, :, :), var0)
        end if
      end do
    else
      if (.not.ismissing(corr(0))) then
        data1@res_gsnRightString = "corr = " + sprintf("%6.3f", corr(0))
      else
        data1@res_gsnRightString = ""
      end if
      if (.not.ismissing(gavg(0))) then
        data1@res_gsnLeftString = "mean = " + sprintf("%6.3f", gavg(0))
      else
        data1@res_gsnLeftString = ""
      end if

      if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
        data1@res_gsnRightString = ""
      end if

      maps(imod, 0) = contour_map(wks0, data1, var0)
    end if

    ; mandatory netcdf output

    data1@var = var0 + "_mean_" + names(imod)
    nc_outfile_mean = ncdf_write(data1, nc_filename_mean)

    ; =======================================================================
    ; Create difference plots (if requested)
    ; =======================================================================

    if (flag_diff .and. (imod .ne. ref_ind)) then

      diff = data1
      if (flag_rel_diff) then
        diff = (diff - ref_data) / ref_data * 100.0
        diff = where(ref_data .le. rel_diff_min, diff@_FillValue, diff)
      else
        diff = diff - ref_data
      end if

      diff@res_gsnLeftString  = ""
      diff@res_gsnRightString = ""

      rmsd = rmsd@_FillValue
      bias = bias@_FillValue

      if (numseas.gt.1) then
        do is = 0, numseas - 1
          if (.not. flag_rel_diff) then
            if (same_grid) then
              rmsd(is) = calculate_metric(ref_data(is, :, :), \
                                          data1(is, :, :), "RMSD")
            end if
            bias(is) = area_operations(diff(is, :, :), -90., 90., 0., 360., \
                                       "average", True)
          end if
        end do
      else
        if (.not. flag_rel_diff) then
          if (same_grid) then
            rmsd(0) = calculate_metric(ref_data, data1, "RMSD")
          end if
          bias(0) = area_operations(diff, -90., 90., 0., 360., "average", \
                                    True)
        end if
      end if

      ; ----------------------------------------------------------------------

      ; ###########################################
      ; # plot ressources                         #
      ; ###########################################

      diff@res_gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
                                               / tofloat((dim_MOD + 1) / 2)/))
      diff@res_gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
                                               / tofloat((dim_MOD + 1) / 2)/))

      diff@res_tiMainOn       = False

      diff@res_cnFillOn       = True      ; color plot desired
      diff@res_cnLineLabelsOn = False     ; contour lines
      diff@res_cnLinesOn      = False

      ; colors
      ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

      ; annotation

      diff@res_cnLevelSelectionMode = "ExplicitLevels"
      diff@res_mpOutlineOn          = True
      diff@res_mpFillOn             = False

      ; variable specific plotting settings

      ; set contour levels / colors

      if (.not.isvar("cnLevels")) then

        if (isatt(diff, "res_cnLevels")) then
          delete(diff@res_cnLevels)
        end if
        if (isatt(diff, "res_cnFillColors")) then
          delete(diff@res_cnFillColors)
        end if
        if (isvar("pal")) then
          delete(pal)
        end if

        if (var0.eq."pr") then
          diff@res_cnLevels = ispan(-30, 30, 5) * 0.1
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "ipcc-precip-delta.rgb")
          diff@res_cnFillColors = pal
          diff@res_lbOrientation = "horizontal"
        end if

        if ((var0.eq."tas") .or. (var0.eq."ts")) then
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "ipcc-tas-delta.rgb")
          diff@res_cnFillPalette = pal
          if (var0.eq."ts") then
            diff@res_cnLevels = ispan(-5, 5, 1) * 0.5
          end if
        end if

        if (var0.eq."lwp") then
          diff@res_cnLevels           = ispan(-50, 50, 10) * 0.001
          diff@res_mpOutlineOn        = False
          diff@res_mpFillOn           = True
          diff@res_mpLandFillColor    = "Black"
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/qcm3.rgb")
          diff@res_cnFillColors       = pal
        end if

        if (var0.eq."clt") then
          diff@res_cnLevels           = fspan(-25, 25, 11)
        end if

        if (var0.eq."clivi") then
          diff@res_cnLevels           = ispan(-70, 70, 10) * 0.001
        end if

        if (var0.eq."clwvi") then
          diff@res_cnLevels           = ispan(-50, 50, 10) * 0.001
        end if

        if (var0.eq."swcre") then
          data1@res_cnLevels          = ispan(-30, 30, 5)
        end if

        if (var0.eq."lwcre") then
          data1@res_cnLevels          = ispan(-30, 30, 5)
        end if

        if (var0.eq."netcre") then
          data1@res_cnLevels          = ispan(-30, 30, 5)
        end if

        ; ******************************************************
        ; *** relative differences: use specific color table ***
        ; ******************************************************

        if (flag_rel_diff) then
          if (isatt(diff, "res_cnLevels")) then
            delete(diff@res_cnLevels)
          end if
          if (isatt(diff, "res_cnFillColors")) then
            delete(diff@res_cnFillColors)
          end if
          diff@res_cnLevels = fspan(-100, 100, 21)
          if (isvar("pal")) then
            delete(pal)
          end if
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "percent100.rgb")
          diff@res_cnFillColors = pal
        end if

        ; ******************************************************

        if (.not. isatt(diff, "res_cnLevels")) then
          log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
          log_info("info: using default contour levels")
          diff@res_cnLevels = fspan(min(diff), max(diff), 20)
        end if

        cnLevels = diff@res_cnLevels
        if (isatt(diff, "res_cnFillColors")) then
          cnFillColors = diff@res_cnFillColors
        end if

      else  ; use previously defined colors and contour intervals

        if (isatt(diff, "res_cnLevels")) then
          delete(diff@res_cnLevels)
        end if
        if (isatt(diff, "res_cnFillColors")) then
          delete(diff@res_cnFillColors)
        end if

        diff@res_cnLevels = cnLevels

        if (isvar("cnFillColors")) then
          diff@res_cnFillColors = cnFillColors
        end if

      end if  ; if .not.isvar("cnLevels")

      if (imod.eq.ref_ind) then
        diff@res_lbLabelBarOn = True
      else
        diff@res_lbLabelBarOn = False
      end if

      ; map attributes

      diff@res_mpFillDrawOrder       = "PostDraw"    ; draw map last
      diff@res_cnMissingValFillColor = "Gray"

      ; no tickmarks and no labels

      diff@res_tmYLLabelsOn       = False
      diff@res_tmYLOn             = False
      diff@res_tmYRLabelsOn       = False
      diff@res_tmYROn             = False
      diff@res_tmXBLabelsOn       = False
      diff@res_tmXBOn             = False
      diff@res_tmXTLabelsOn       = False
      diff@res_tmXTOn             = False
      diff@res_cnInfoLabelOn      = False    ; turn off cn info label

      ; specified in namelist

      diff@res_mpProjection       = projection

      ; set explicit contour levels

      if (isatt(diag_script_info, "explicit_cn_levels")) then
        diff@res_cnLevelSelectionMode = "ExplicitLevels"
        if (isatt(diff, "res_cnLevels")) then
          delete(diff@res_cnLevels)
        end if
        diff@res_cnLevels = diag_script_info@explicit_cn_levels
      end if

      ; ###########################################
      ; # other Metadata: diag_script, var        #
      ; ###########################################
      ; add to diff as attributes without prefix

      if (isatt(variable_info, "long_name")) then
        diff@var_long_name = variable_info@long_name
      end if
      if (isatt(variable_info, "units")) then
        diff@var_units = variable_info@units
      else
        diff@var_units = ""
      end if

      ; ###########################################
      ; # create the plot                         #
      ; ###########################################

      diff@res_gsnDraw        = False  ; do not draw yet
      diff@res_gsnFrame       = False  ; don't advance frame

      ; ----------------------------------------------------------------------

      if (numseas.gt.1) then
        do is = 0, numseas - 1
          if (.not.ismissing(rmsd(is))) then
            diff@res_gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(is))
          else
            diff@res_gsnRightString = ""
          end if
          if (.not.ismissing(bias(is))) then
            diff@res_gsnLeftString = "bias = " + sprintf("%6.3f", bias(is))
          else
            diff@res_gsnLeftString = ""
          end if

          if (is.eq.0) then
            maps_d(imod, is) = contour_map(wks0d, diff(is, :, :), var0)
          end if
          if (is.eq.1) then
            maps_d(imod, is) = contour_map(wks1d, diff(is, :, :), var0)
          end if
          if (is.eq.2) then
            maps_d(imod, is) = contour_map(wks2d, diff(is, :, :), var0)
          end if
          if (is.eq.3) then
            maps_d(imod, is) = contour_map(wks3d, diff(is, :, :), var0)
          end if
        end do
      else
        if (.not.ismissing(rmsd(0))) then
          diff@res_gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(0))
        else
          diff@res_gsnRightString = ""
        end if
        if (.not.ismissing(bias(0))) then
          diff@res_gsnLeftString = "bias = " + sprintf("%6.3f", bias(0))
        else
          diff@res_gsnLeftString = ""
        end if
        maps_d(imod, 0) = contour_map(wks0d, diff, var0)
      end if

      ; mandatory netcdf output

      diff@var = var0 + "_bias_" + names(imod)
      nc_outfile_bias = ncdf_write(diff, nc_filename_bias)

    end if  ; if flag_diff

    ; =======================================================================

  end do  ; ii-loop (models)

  if (write_plots) then
    ; save default color map in case it is needed later for optionally
    ; plotting color bar to a separate file

    tmp_colors = gsn_retrieve_colormap(wks0)
    cdims = dimsizes(tmp_colors)
    nboxes = dimsizes(data1@res_cnLevels)
    clen = cdims(0)
    stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
    fill_colors = ispan(2, clen(0) - 1, stride)
    mean_colors = tmp_colors(fill_colors, :)
    delete(tmp_colors)
    delete(fill_colors)
    delete(cdims)

    ; sort plots if needed (observations go first)

    plottmp = ispan(0, dim_MOD - 1, 1)
    plotind = plottmp

    ; move plots of observational datasets (if present) into the first line(s)
    ; of the panel plot

    j = 0
    do i = 0, dimsizes(plottmp) - 1
      if (i.eq.ref_ind) then
        plotind(j) = plottmp(i)
        j = j + 1
      else if (plottmp(i) .lt. dimsizes(projects)) then
        if (isStrSubset(str_lower(projects(plottmp(i))), \
                        "obs")) then
          plotind(j) = plottmp(i)
          j = j + 1
        end if
      end if
      end if
    end do

    do i = 0, dimsizes(plottmp) - 1
      if ((isStrSubset(str_lower(projects(plottmp(i))), \
                       "obs")).or.(i.eq.ref_ind)) then
      else
        plotind(j) = plottmp(i)
        j = j + 1
      end if
    end do

    pres                      = True    ; needed to override
                                        ; panelling defaults
    pres@gsnPanelLabelBar     = True    ; add common colorbar
    if (panel_labels) then
      ; print dataset name on each panel
      pres@gsnPanelFigureStrings = names(plotind)
    end if
    pres@gsnPanelFigureStringsFontHeightF = min((/0.01, 0.01 * 6.0 \
                                                / tofloat((dim_MOD + 1) / 2)/))
    pres@lbLabelFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                                / tofloat((dim_MOD + 1) / 2)/))
    pres@lbAutoManage                     = False
    pres@lbTopMarginF                     = 0.1
    pres@lbTitleOn                        = True
    pres@lbTitleFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                                / tofloat((dim_MOD + 1) / 2)/))
    pres@lbTitlePosition                  = "Bottom"
    pres@lbTitleString                    = data1@long_name + " (" \
      + data1@units + ")"
    pres@lbPerimOn                        = False   ; draw line around label
                                                    ; bar area
    pres@gsnPanelCenter                   = False
    if (dim_MOD.le.8) then
      pres@pmLabelBarOrthogonalPosF       = -0.03
    else
      pres@pmLabelBarOrthogonalPosF       = -0.01   ; shift label bar a bit to
                                                    ; the bottom
    end if

    if (embracesetup) then
      if (numseas.gt.1) then
        pres@txString = season(0)
        outfile(0) = panelling(wks0, maps(plotind, 0), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(1)
        outfile(1) = panelling(wks1, maps(plotind, 1), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(2)
        outfile(2) = panelling(wks2, maps(plotind, 2), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(3)
        outfile(3) = panelling(wks3, maps(plotind, 3), (dim_MOD + 3) / 4, \
                               4, pres)
        log_info(" Wrote " + outfile)
      else
        pres@gsnPanelRowSpec = True             ; tell panel what order to plt
        pres@gsnPanelYWhiteSpacePercent = 5
        pres@gsnPanelXWhiteSpacePercent = 5
        if (isatt(diag_script_info, "PanelTop")) then
          top = tofloat(diag_script_info@PanelTop)
        else
          top = 0.99  ; default
        end if
        pres@gsnPanelTop = top

        if (isvar("plotsperline")) then
          delete(plotsperline)
        end if

        plotsperline = new((dim_MOD + 1) / 2, integer)
        plotsperline = 2

        if ((isStrSubset(str_lower(projects(plotind(0))), \
                         "obs")).and. \
           .not.(isStrSubset(str_lower(projects(plotind(1))), \
                             "obs"))) then
          plotsperline(0) = 1
        end if

        if (sum(plotsperline).gt.dimsizes(plotind)) then
          plotsperline(dimsizes(plotsperline) - 1) = 1
        end if

        if (sum(plotsperline).lt.dimsizes(plotind)) then
          xadd = 1
          xtmp = array_append_record(plotsperline, xadd, 0)
          delete(plotsperline)
          plotsperline = xtmp
          delete(xtmp)
        end if

        gsn_panel(wks0, maps(plotind, 0), plotsperline, pres)
        outfile(0) = wks0@fullname
      end if
    else  ; if embracesetup
      if (numseas.gt.1) then
        pres@txString = season(0)
        outfile(0) = panelling(wks0, maps(plotind, 0), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(1)
        outfile(1) = panelling(wks1, maps(plotind, 1), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(2)
        outfile(2) = panelling(wks2, maps(plotind, 2), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(3)
        outfile(3) = panelling(wks3, maps(plotind, 3), (dim_MOD + 3) / 4, \
                               4, pres)
      else
        outfile(0) = panelling(wks0, maps(plotind, 0), (dim_MOD + 3) / 4, \
                               4, pres)
      end if
    end if  ; if embracesetup

    do is = 0, numseas - 1
      log_info("Wrote " + outfile(is))
    end do

    ; ------------------------------------------------------------------------
    ; write provenance to netcdf output and plot file(s) (mean)
    ; ------------------------------------------------------------------------

    statistics = (/"clim", "mean"/)
    if (isatt(diag_script_info, "region")) then
      domain = "reg"
    else
      domain = "global"
    end if
    plottype = "geo"

    do is = 0, numseas - 1
      caption = "Mean values for variable " + var0 \
                + " (" + allseas + ")."
      log_provenance(nc_outfile_mean, outfile(is), caption, statistics, \
                     domain, plottype, "", "", climofiles)
    end do

    ; ========================================================================

    if (flag_diff) then
      pres@lbTitleString = "~F33~D~F21~" + diff@long_name + " (" + \
                           diff@units + ")"

      ; save default color map in case it is needed later for optionally
      ; plotting color bar to a separate file

      if (isvar("nboxes")) then
        delete(nboxes)
      end if

      tmp_colors = gsn_retrieve_colormap(wks0d)
      cdims = dimsizes(tmp_colors)
      nboxes = dimsizes(diff@res_cnLevels)
      clen = cdims(0)
      stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
      fill_colors = ispan(2, clen(0) - 1, stride)
      diff_colors = tmp_colors(fill_colors, :)
      delete(tmp_colors)
      delete(fill_colors)
      delete(cdims)

      if (isvar("plottmp")) then
        delete(plottmp)
      end if

      if (isvar("plotind")) then
        delete(plotind)
      end if

      plottmp = ind(ispan(0, dim_MOD - 1, 1).ne.ref_ind)
      plotind = plottmp

      ; if there is a second observational dataset, move the corresponding
      ; plot to the first line of the panel plot

      j = 0
      do i = 0, dimsizes(plottmp) - 1
        if (isStrSubset(str_lower(projects(plottmp(i))), "obs")) then
          plotind(j) = plottmp(i)
          j = j + 1
        end if
      end do
      do i = 0, dimsizes(plottmp) - 1
        if (isStrSubset(str_lower(projects(plottmp(i))), "obs")) then
        else
          plotind(j) = plottmp(i)
          j = j + 1
        end if
      end do

      if (isatt(pres, "gsnPanelFigureStrings")) then
        delete(pres@gsnPanelFigureStrings)
      end if
      if (panel_labels) then
        pres@gsnPanelFigureStrings = names(plotind)
      end if

      if (dimsizes(plotind).eq.1) then
        pres@gsnPanelRight = 0.5
      end if

      if (embracesetup) then
        if (numseas.gt.1) then
          pres@txString = season(0)
          outfile_d(0) = panelling(wks0d, maps_d(plotind, 0), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString = season(1)
          outfile_d(1) = panelling(wks1d, maps_d(plotind, 1), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString = season(2)
          outfile_d(2) = panelling(wks2d, maps_d(plotind, 2), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString = season(3)
          outfile_d(3) = panelling(wks3, maps_d(plotind, 3), \
                                   (dim_MOD + 3) / 4, 4, pres)
        else
          pres@gsnPanelRowSpec = True           ; tell panel what order to plt
          pres@gsnPanelYWhiteSpacePercent = 5
          pres@gsnPanelXWhiteSpacePercent = 5
          pres@gsnPanelTop = tofloat(diag_script_info@PanelTop)

          if (isvar("plotsperline")) then
            delete(plotsperline)
          end if

          plotsperline = new(max((/1, dim_MOD / 2/)), integer)
          plotsperline = 2

          if (dimsizes(plotind).gt.1) then
            if ((isStrSubset(str_lower(projects(plotind(0))), "obs")).and. \
               .not. \
                (isStrSubset(str_lower(projects(plotind(1))), "obs"))) then
              plotsperline(0) = 1
            end if
          end if

          if (sum(plotsperline).gt.dimsizes(plotind)) then
            plotsperline(dimsizes(plotsperline) - 1) = 1
          end if

          if (sum(plotsperline).lt.dimsizes(plotind)) then
            xadd = 1
            xtmp = array_append_record(plotsperline, xadd, 0)
            delete(plotsperline)
            plotsperline = xtmp
            delete(xtmp)
          end if

          gsn_panel(wks0d, maps_d(plotind, 0), plotsperline, pres)
          outfile_d(0) = wks0d@fullname
        end if
      else  ; embracesetup = False
        if (numseas.gt.1) then
          pres@txString  = season(0)
          outfile_d(0) = panelling(wks0d, maps_d(plotind, 0), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString  = season(1)
          outfile_d(1) = panelling(wks1d, maps_d(plotind, 1), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString  = season(2)
          outfile_d(2) = panelling(wks2d, maps_d(plotind, 2), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString  = season(3)
          outfile_d(3) = panelling(wks3d, maps_d(plotind, 3), \
                                   (dim_MOD + 3) / 4, 4, pres)
        else
          outfile_d(0) = panelling(wks0d, maps_d(plotind, 0), \
                                   (dim_MOD + 3) / 4, 4, pres)
        end if
      end if  ; end if embracesetup

      do is = 0, numseas - 1
        log_info(" Wrote " + outfile(is))

        ; --------------------------------------------------------------------
        ; write provenance to netcdf output and plot file(s) (bias)
        ; --------------------------------------------------------------------

        statistics = (/"clim", "diff"/)
        if (isatt(diag_script_info, "region")) then
          domain = "reg"
        else
          domain = "global"
        end if
        plottype = "geo"

        ; note: because function log_provenance does not yet support to attach
        ;       different captions to netcdf (contains all seasons) and plots
        ;       (contain one season each), the caption cannot specifiy the
        ;       season plotted; using "annual" or "DJF/MAM/JJA/SON" instead.

        caption = "Differences for variable " + var0 \
                  + " (" + allseas + "), reference = " + refname + "."
        log_provenance(nc_outfile_bias, outfile_d(is), caption, statistics, \
                       domain, plottype, "", "", climofiles)
      end do

    end if  ; if flag_diff

    ; optionally save legend(s) to extra file(s)

    if (extralegend) then
      nboxes = dimsizes(data1@res_cnLevels) + 1
      wksleg = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 \
                       + "_legend")
      pres@lbMonoFillPattern = True
      pres@lbOrientation = "Horizontal"
      pres@vpWidthF = 0.7
      pres@vpHeightF = 0.1
      pres@lbLabelFontHeightF = 0.015
      pres@lbLabelAlignment = "InteriorEdges"
      pres@lbTitleFontHeightF = 0.015
      pres@lbTitleString = data1@long_name + " (" + data1@units + ")"

      labels = tostring(data1@res_cnLevels)

      ; remove trailing zeros from strings

      do i = 0, dimsizes(labels) - 1
        i1  = str_index_of_substr(labels(i), ".", -1)
        if (.not.ismissing(i1)) then
          tmp = stringtochar(labels(i))
          do j = dimsizes(tmp) - 2, i1, 1
            if ((tmp(j).ne.".").and.(tmp(j).ne."0")) then
              break
            end if
          end do
          labels(i) = chartostring(tmp(0:j))
          delete(tmp)
        end if
      end do

      if (isatt(data1, "res_cnFillColors")) then
        pres@lbFillColors = data1@res_cnFillColors
      else if (isatt(data1, "res_cnFillPalette")) then
        pres@lbFillColors = data1@res_cnFillPalette
      else
        pres@lbFillColors = mean_colors  ; default colors
      end if
      end if

      gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)

      delete(wksleg)
      delete(labels)
      delete(pres@lbFillColors)

      if (flag_diff) then
        nboxes = dimsizes(diff@res_cnLevels) + 1
        wksleg = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 \
                         + "_diff_legend")

        labels = tostring(diff@res_cnLevels)

        ; remove trailing zeros from strings

        do i = 0, dimsizes(labels) - 1
          i1  = str_index_of_substr(labels(i), ".", -1)
          if (.not.ismissing(i1)) then
            tmp = stringtochar(labels(i))
            do j = dimsizes(tmp) - 2, i1, 1
              if ((tmp(j).ne.".").and.(tmp(j).ne."0")) then
                break
              end if
            end do
            labels(i) = chartostring(tmp(0:j))
            delete(tmp)
          end if
        end do

        if (flag_rel_diff) then
          pres@lbTitleString = "~F33~D~F21~" + data1@long_name + " (%)"
        else
          pres@lbTitleString = "~F33~D~F21~" + data1@long_name + " (" + \
                        data1@units + ")"
        end if

        if (isatt(diff, "res_cnFillColors")) then
          pres@lbFillColors = diff@res_cnFillColors
        else if (isatt(diff, "res_cnFillPalette")) then
          pres@lbFillColors = diff@res_cnFillPalette
        else
          pres@lbFillColors = diff_colors  ; default colors
        end if
        end if

        gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)
      end if  ; if (flag_diff)
    end if  ; if (extralegend)
  end if  ; if write_plots

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
