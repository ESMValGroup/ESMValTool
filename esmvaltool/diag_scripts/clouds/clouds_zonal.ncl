; CLOUDS_ZONAL
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates annual/seasonal means of zonally averaged 3-d (cloud)
;   parameters for comparison with a reference data set. Optionally,
;   differences to the reference data set are also plotted.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   embracesetup:       True = 2 plots per line, False = 4 plots per line
;                       (default)
;   explicit_cn_levels: explicit contour levels for mean values (array)
;   explicit_cn_dlevels: explicit contour levels for differences (array)
;   extralegend:        plot legend(s) to extra file(s)
;   filename_add:       optionally add this string to plot filesnames
;   panel_labels:       label individual panels (true, false)
;   PanelTop:           manual override for "@gnsPanelTop" used by panel
;                       plot(s)
;   showdiff:           calculate and plot differences (default = False)
;   showyears:          add start and end years to the plot titles
;                       (default = False)
;   rel_diff:           if showdiff = True, then plot relative differences (%)
;                       (default = False)
;   rel_diff_min:       lower cutoff value in case of calculating relative
;                       differences
;                       (in units of input variable)
;   t_test:             perform t-test when calculating differences
;                       (default = False)
;   timemean:           time averaging - "seasonal" = DJF, MAM, JJA, SON),
;                                        "annual" = annual mean
;   units_to:           target units (automatic conversion)
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
;   20230117-lauer_axel: added support for ICON (code from Manuel)
;   20200211-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"

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
  set_default_att(diag_script_info, "showyears", False)
  set_default_att(diag_script_info, "t_test", False)
  set_default_att(diag_script_info, "timemean", "annualclim")
  set_default_att(diag_script_info, "units_to", "")

  flag_diff = diag_script_info@showdiff
  flag_rel_diff = diag_script_info@rel_diff
  rel_diff_min = diag_script_info@rel_diff_min
  t_test = diag_script_info@t_test

  if (.not.flag_diff) then
    if (flag_rel_diff) then
      log_info("rel_diff = True has no effect until showdiff is also " \
               + "set to True")
    end if
    if (t_test) then
      log_info("t_test = True has no effect until showdiff is also " \
               + "set to True")
      t_test = False
    end if
  end if

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  embracesetup = diag_script_info@embrace_setup

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

  units_to = diag_script_info@units_to

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  panel_labels = diag_script_info@panel_labels

  extralegend = diag_script_info@extralegend

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

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
                " dimensions, need 3")
    end if
    idx = ind(dims .eq. "lat")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
    end if
    idx = ind(dims .eq. "time")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no time dimension")
    end if

;    ; if coordinate variables do not have the attribute "long_name",
;    ; try to use attribute "standard_name" as a substitute
;
;    do n = 0, dimsizes(dims) - 1
;      if (.not.isatt(A0&$dims(n)$, "long_name")) then
;        if (isatt(A0&$dims(n)$, "standard_name")) then
;          A0&$dims(n)$@long_name = A0&$dims(n)$@standard_name
;        end if
;      end if
;    end do

    vcoord = dims(1)  ; save name of vertical coordinate variable
    if (ii .eq. 0) then
      vcoord0_var = A0&$vcoord$
      vcoord0 = vcoord
      if (isatt(vcoord0_var, "units")) then
        vcoord0_units = vcoord0_var@units
      else
        vcoord0_units = ""
      end if
    end if

    delete(dims)

    ; vertical coordinate is assumed to be the dimension not being
    ; "time" and "lat"

    ; average over time

    data1 = time_operations(A0, -1, -1, "average", timemean, True)

    if (t_test) then
      start_year = info0[imod]@start_year
      end_year = info0[imod]@end_year
      nyears = end_year - start_year + 1

      if (nyears .lt. 3) then
        log_info("warning: cannot calculate t-test for dataset " \
                 + names(imod) + "; need at least 3 years, dataset " \
                 + "length = " + tostring(nyears) + " years; disabling t-test")
        t_test = False
      else
        if (isvar("data1_t_stddev")) then
          delete(data1_t_stddev)
        end if
        if (isvar("data1_t_mean")) then
          delete(data1_t_mean)
        end if
        data1_t_stddev = interannual_variability(A0, -1, -1, timemean, "None")
        data1_t_mean = data1
      end if
    end if

    delete(A0)

    if (units_to .ne. "") then
      data0 = convert_units(data1, units_to)
      delete(data1)
      data1 = data0
      delete(data0)
      ; if attribute is present, overwrite with user specified units
      if (isatt(variable_info[0], "units")) then
        variable_info[0]@units = units_to
      end if
    end if

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/shared/plot/style.ncl

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    res = True

    res@cnFillOn       = True      ; color plot desired
    res@cnLineLabelsOn = False     ; contour lines

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    ; annotation

    ; if desired, add years to plot title
    years_str = ""
    if (diag_script_info@showyears) then
      years_str = " (" + variable_info[0]@start_year
      if (variable_info[0]@start_year .ne. variable_info[0]@end_year) then
        years_str = years_str + "-" + variable_info[0]@end_year
      end if
      years_str = years_str + ")"
    end if

;    res@tiMainOn             = False
    res@tiMainString = names(imod) + years_str
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn            = False

;    res@lbLabelBarOn         = False
    res@gsnRightString       = ""

    res@cnMissingValFillColor = "Gray"

    res@cnInfoLabelOn      = False    ; turn off cn info label

    if (isatt(data1&$vcoord$, "standard_name")) then
      res@tiYAxisString = data1&$vcoord$@standard_name
      if (isatt(data1&$vcoord$, "units")) then
        res@tiYAxisString = res@tiYAxisString + " (" \
          + data1&$vcoord$@units + ")"
      end if
    end if

    if (isatt(data1&lat, "standard_name")) then
      res@tiXAxisString = data1&lat@standard_name
    end if

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@cnLevels = diag_script_info@explicit_cn_levels
    end if

    if (.not.isatt(res, "cnLevels")) then
      if (var0.eq."clcalipso") then
        res@cnLevels            = fspan(5, 50, 10)
      else
        log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
        log_info("info: using default contour levels")
        res@cnLevels = fspan(min(data1), max(data1), 20)
      end if
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

    res@lbTitleOn          = True
    res@lbTitleString      = data1@var_units
    res@lbTitlePosition    = "Bottom"
    res@lbTitleFontHeightF = 0.015

    if (.not. isvar("ref_data")) then
      ref_data = data1
      vcoord_ref = vcoord
      if (t_test) then  ; save mean in same units as stddev
        ref_data_t_stddev = data1_t_stddev
        ref_data_t_mean = data1_t_mean
        nyears_ref = nyears
      end if
    end if

    ; check if data are on same grid (for calculating difference, RMSD,
    ; correlation)

    same_grid = False

    if (all(dimsizes(ref_data) .eq. dimsizes(data1))) then
      if (max(abs(ref_data&lat - data1&lat)) .le. 1.0e-6) then
        if (max(abs(ref_data&$vcoord_ref$ - data1&$vcoord$)) .le. 1.0e-6) then
          same_grid = True
        end if
      end if
    end if

    if (flag_diff .and. .not.same_grid) then
      flag_diff = False
      error_msg("f", DIAG_SCRIPT, "", \
                "Data are not on same grid (horizontal and / or vertical), " \
                + "cannot calculate differences. " \
                + "Set showdiff to False in recipe or regrid data to " \
                + "common grid (check/adjust " \
                + "preprocessor settings in recipe).")
    end if

    corr = corr@_FillValue
    gavg = gavg@_FillValue

;    if (.not.all(ismissing(data1))) then
;      if (numseas.gt.1) then
;        do is = 0, numseas - 1
;          if (same_grid .and. (ref_ind .ge. 0)) then
;            corr(is) = calculate_metric(ref_data(is, :, :), data1(is, :, :), \
;                                        "correlation")
;          end if
;          gavg(is) = dim_avg_n_Wrap(data1(is, :, :), (/0, 1/))
;        end do
;      else
;        if (same_grid .and. (ref_ind .ge. 0)) then
;          corr(0) = calculate_metric(ref_data, data1, "correlation")
;        end if
;        gavg(0) = dim_avg_n_Wrap(data1, (/0, 1/))
;      end if
;    end if

    res@gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
                                         / tofloat((dim_MOD + 1) / 2)/))
    res@gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
                                         / tofloat((dim_MOD + 1) / 2)/))

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      ; note: an array of workspaces (i.e. wks(numseas)) does not work as
      ;       attributes cannot be assigned to each array element
      ;       individually
      wks0 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + var0 + \
                     "_" + season(0) + filename_add)
      ; difference plots will be saved to a different file
      if (flag_diff) then
        wks0d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + \
                        var0 + "_bias_" + season(0) + filename_add)
        ndframe = 0
      end if
      if (numseas.gt.1) then
        wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + var0 + \
                       "_" + season(1) + filename_add)
        wks2 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + var0 + \
                       "_" + season(2) + filename_add)
        wks3 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + var0 + \
                       "_" + season(3) + filename_add)
        ; difference plots will be saved to a different files
        if (flag_diff) then
          wks1d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + \
                          var0 + "_bias_" + season(1) + filename_add)
          wks2d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + \
                          var0 + "_bias_" + season(2) + filename_add)
          wks3d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_zonal_" + \
                          var0 + "_bias_" + season(3) + filename_add)
        end if
      end if
    end if

    if (numseas.gt.1) then
      do is = 0, numseas - 1
        if (.not.ismissing(corr(is))) then
          res@gsnRightString = "corr = " + sprintf("%6.3f", corr(is))
        else
          res@gsnRightString = ""
        end if
        if (.not.ismissing(gavg(is))) then
          res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg(is))
        else
          res@gsnLeftString = ""
        end if

        if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
          res@gsnRightString = ""
        end if

        if (vcoord0_units .eq. "Pa") then
          if (is.eq.0) then
            maps(imod, is) = gsn_csm_pres_hgt(wks0, data1(is, :, :), res)
          end if
          if (is.eq.1) then
            maps(imod, is) = gsn_csm_pres_hgt(wks1, data1(is, :, :), res)
          end if
          if (is.eq.2) then
            maps(imod, is) = gsn_csm_pres_hgt(wks2, data1(is, :, :), res)
          end if
          if (is.eq.3) then
            maps(imod, is) = gsn_csm_pres_hgt(wks3, data1(is, :, :), res)
          end if
        else
          if (is.eq.0) then
            maps(imod, is) = gsn_csm_contour(wks0, data1(is, :, :), res)
          end if
          if (is.eq.1) then
            maps(imod, is) = gsn_csm_contour(wks1, data1(is, :, :), res)
          end if
          if (is.eq.2) then
            maps(imod, is) = gsn_csm_contour(wks2, data1(is, :, :), res)
          end if
          if (is.eq.3) then
            maps(imod, is) = gsn_csm_contour(wks3, data1(is, :, :), res)
          end if
        end if
      end do
    else
      if (.not.ismissing(corr(0))) then
        res@gsnRightString = "corr = " + sprintf("%6.3f", corr(0))
      else
        res@gsnRightString = ""
      end if
      if (.not.ismissing(gavg(0))) then
        res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg(0))
      else
        res@gsnLeftString = ""
      end if

      if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
        res@gsnRightString = ""
      end if
      if (vcoord0_units .eq. "Pa") then
        maps(imod, 0) = gsn_csm_pres_hgt(wks0, data1, res)
      else
        maps(imod, 0) = gsn_csm_contour(wks0, data1, res)
      end if
    end if

    ; mandatory netcdf output

    data1@var = var0 + "_mean_" + names(imod)

    ; vertical coordinate might have a different name, which prevents
    ; all data to be written to the same netcdf file
    ; --> rename vertical coordinate to match the one of the first dataset
    ;     written to the netcdf file
    delete(data1&$vcoord$)
    data1!0 = vcoord0
    data1&$vcoord0$ = vcoord0_var

    nc_outfile_mean = ncdf_write(data1, nc_filename_mean)

    ; =======================================================================
    ; Create difference plots (if requested)
    ; =======================================================================

    if (flag_diff .and. (imod .ne. ref_ind)) then

      diff = data1
      if (flag_rel_diff) then
        ref_data = where(ref_data .le. 1.0e-19, ref_data@_FillValue, \
                         ref_data)
        diff = (diff - ref_data) / ref_data * 100.0
        diff = where(ref_data .le. rel_diff_min, diff@_FillValue, diff)
      else
        diff = diff - ref_data
      end if

      dres = res

      dres@gsnLeftString  = ""
      dres@gsnRightString = ""

      rmsd = rmsd@_FillValue
      bias = bias@_FillValue

;      if (numseas.gt.1) then
;        do is = 0, numseas - 1
;          if (.not. flag_rel_diff) then
;            if (same_grid) then
;              rmsd(is) = calculate_metric(ref_data(is, :, :), \
;                                          data1(is, :, :), "RMSD")
;            end if
;            bias(is) = dim_avg_n_Wrap(diff(is, :, :), (/0, 1/))
;          end if
;        end do
;      else
;        if (.not. flag_rel_diff) then
;          if (same_grid) then
;            rmsd(0) = calculate_metric(ref_data, data1, "RMSD")
;          end if
;          bias(0) = dim_avg_n_Wrap(diff, (/0, 1/))
;        end if
;      end if

      ; ----------------------------------------------------------------------

      ; ###########################################
      ; # plot ressources                         #
      ; ###########################################

      dres@gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
                                            / tofloat((dim_MOD + 1) / 2)/))
      dres@gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
                                            / tofloat((dim_MOD + 1) / 2)/))

;      dres@tiMainOn       = False
      dres@tiMainString = names(imod) + " - " + refname + years_str

      dres@cnFillOn       = True      ; color plot desired
      dres@cnLineLabelsOn = False     ; contour lines
      dres@cnLinesOn      = False

      ; colors
      ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

      ; annotation

      dres@cnLevelSelectionMode = "ExplicitLevels"

      ; variable specific plotting settings

      ; set contour levels / colors

      if (.not.isvar("cnLevels")) then

        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        if (isatt(dres, "cnFillColors")) then
          delete(dres@cnFillColors)
        end if
        if (isvar("pal")) then
          delete(pal)
        end if

        if (var0.eq."clcalipso") then
          dres@cnLevels           = fspan(-25, 25, 11)
        end if

        ; ******************************************************
        ; *** relative differences: use specific color table ***
        ; ******************************************************

        if (flag_rel_diff) then
          if (isatt(dres, "cnLevels")) then
            delete(dres@cnLevels)
          end if
          if (isatt(dres, "cnFillColors")) then
            delete(dres@cnFillColors)
          end if
          dres@cnLevels = fspan(-100, 100, 21)
          if (isvar("pal")) then
            delete(pal)
          end if
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "percent100.rgb")
          dres@cnFillColors = pal
        end if

        ; ******************************************************

        if (.not. isatt(dres, "cnLevels")) then
          log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
          log_info("info: using default contour levels")
          dres@cnLevels = fspan(min(diff), max(diff), 20)
        end if

        cnLevels = dres@cnLevels
        if (isatt(dres, "cnFillColors")) then
          cnFillColors = dres@cnFillColors
        end if

      else  ; use previously defined colors and contour intervals

        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        if (isatt(dres, "cnFillColors")) then
          delete(dres@cnFillColors)
        end if

        dres@cnLevels = cnLevels

        if (isvar("cnFillColors")) then
          dres@cnFillColors = cnFillColors
        end if

      end if  ; if .not.isvar("cnLevels")

;      if (imod.eq.ref_ind) then
;        dres@lbLabelBarOn = True
;      else
;        dres@lbLabelBarOn = False
;      end if

      ; map attributes

      dres@cnMissingValFillColor = "Gray"

      dres@cnInfoLabelOn      = False    ; turn off cn info label

      ; set explicit contour levels

      if (isatt(diag_script_info, "explicit_cn_dlevels")) then
        dres@cnLevelSelectionMode = "ExplicitLevels"
        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        dres@cnLevels = diag_script_info@explicit_cn_dlevels
      end if

      ; ###########################################
      ; # other Metadata: diag_script, var        #
      ; ###########################################
      ; add to diff as attributes without prefix

      if (isatt(variable_info[0], "long_name")) then
        diff@var_long_name = variable_info[0]@long_name
      end if
      if (isatt(variable_info[0], "units")) then
        diff@var_units = variable_info[0]@units
      else
        diff@var_units = ""
      end if

      ; ###########################################
      ; # create the plot                         #
      ; ###########################################

      if (t_test) then
        dres@gsnDraw        = False  ; do not draw yet
        dres@gsnFrame       = False  ; don't advance frame
      end if

      ; ----------------------------------------------------------------------

      if (numseas.gt.1) then
        do is = 0, numseas - 1
          if (.not.ismissing(rmsd(is))) then
            dres@gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(is))
          else
            dres@gsnRightString = ""
          end if
          if (.not.ismissing(bias(is))) then
            dres@gsnLeftString = "bias = " + sprintf("%6.3f", bias(is))
          else
            dres@gsnLeftString = ""
          end if

          if (vcoord0_units .eq. "Pa") then
            if (is.eq.0) then
              maps_d(imod, is) = gsn_csm_pres_hgt(wks0d, diff(is, :, :), dres)
            end if
            if (is.eq.1) then
              maps_d(imod, is) = gsn_csm_pres_hgt(wks1d, diff(is, :, :), dres)
            end if
            if (is.eq.2) then
              maps_d(imod, is) = gsn_csm_pres_hgt(wks2d, diff(is, :, :), dres)
            end if
            if (is.eq.3) then
              maps_d(imod, is) = gsn_csm_pres_hgt(wks3d, diff(is, :, :), dres)
            end if
          else
            if (is.eq.0) then
              maps_d(imod, is) = gsn_csm_contour(wks0d, diff(is, :, :), dres)
            end if
            if (is.eq.1) then
              maps_d(imod, is) = gsn_csm_contour(wks1d, diff(is, :, :), dres)
            end if
            if (is.eq.2) then
              maps_d(imod, is) = gsn_csm_contour(wks2d, diff(is, :, :), dres)
            end if
            if (is.eq.3) then
              maps_d(imod, is) = gsn_csm_contour(wks3d, diff(is, :, :), dres)
            end if
          end if
        end do
      else
        if (.not.ismissing(rmsd(0))) then
          dres@gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(0))
        else
          dres@gsnRightString = ""
        end if
        if (.not.ismissing(bias(0))) then
          dres@gsnLeftString = "bias = " + sprintf("%6.3f", bias(0))
        else
          dres@gsnLeftString = ""
        end if
        if (vcoord0_units .eq. "Pa") then
          maps_d(imod, 0) = gsn_csm_pres_hgt(wks0d, diff, dres)
        else
          maps_d(imod, 0) = gsn_csm_contour(wks0d, diff, dres)
        end if
      end if
      ndframe = ndframe + 1

      ; mandatory netcdf output

      diff@var = var0 + "_bias_" + names(imod)
      nc_outfile_bias = ncdf_write(diff, nc_filename_bias)

      ; ---------------------------------------------------------------------
      ; optionally mask non-significant grid cells
      ; ---------------------------------------------------------------------

      if (t_test) then
        n1 = nyears
        n2 = nyears_ref

        tres = True

        tres@gsnDraw        = False  ; do not draw yet
        tres@gsnFrame       = False  ; don't advance frame
        tres@cnMissingValFillColor = -1
        tres@cnLevelSelectionMode = "ExplicitLevels"
        tres@cnLevels = 0.95
;        tres@cnFillColors = (/"gray70", "transparent"/)
        tres@cnFillColors = (/"black", "transparent"/)
        tres@cnFillPattern = 17
        tres@cnFillOn       = True      ; color plot desired
        tres@cnInfoLabelOn = False
        tres@cnLinesOn = False
        tres@cnLineLabelsOn = False
        tres@lbLabelBarOn = False
        tres@gsnRightString   = ""
        tres@gsnLeftString    = ""
        tres@gsnCenterString  = ""
        tres@tiYAxisOn = False
        tres@tmXBBorderOn = False
        tres@tmXTBorderOn = False
        tres@tmYLBorderOn = False
        tres@tmYRBorderOn = False
        tres@tmXBLabelsOn = False
        tres@tmYLLabelsOn = False
        tres@tmXBOn = False
        tres@tmXTOn = False
        tres@tmYLOn = False
        tres@tmYROn = False

        if (numseas.gt.1) then
          do is = 0, numseas - 1
            x1 = data1_t_mean(is, :, :)
            x2 = ref_data_t_mean(is, :, :)
            s1 = data1_t_stddev(is, :, :) ^ 2
            s2 = ref_data_t_stddev(is, :, :) ^ 2

            prob = ttest(x1, s1, n1, x2, s2, n2, True, False)
            plot_var = 1. - prob
            copy_VarCoords(data1_t_mean, plot_var)

            if (vcoord0_units .eq. "Pa") then
              if (is.eq.0) then
                tmask = gsn_csm_pres_hgt(wks0d, plot_var, tres)
              end if
              if (is.eq.1) then
                tmask = gsn_csm_pres_hgt(wks1d, plot_var, tres)
              end if
              if (is.eq.2) then
                tmask = gsn_csm_pres_hgt(wks2d, plot_var, tres)
              end if
              if (is.eq.3) then
                tmask = gsn_csm_pres_hgt(wks3d, plot_var, tres)
              end if
            else
              if (is.eq.0) then
                tmask = gsn_csm_contour(wks0d, plot_var, tres)
              end if
              if (is.eq.1) then
                tmask = gsn_csm_contour(wks1d, plot_var, tres)
              end if
              if (is.eq.2) then
                tmask = gsn_csm_contour(wks2d, plot_var, tres)
              end if
              if (is.eq.3) then
                tmask = gsn_csm_contour(wks3d, plot_var, tres)
              end if
            end if

            overlay(maps_d(imod, is), tmask)
            delete(tmask)
            draw(maps_d(imod, is))

            if (is.eq.0) then
              frame(wks0d)
            end if
            if (is.eq.1) then
              frame(wks1d)
            end if
            if (is.eq.2) then
              frame(wks2d)
            end if
            if (is.eq.3) then
              frame(wks3d)
            end if

          end do
        else
          x1 = data1_t_mean
          x2 = ref_data_t_mean
          s1 = data1_t_stddev ^ 2
          s2 = ref_data_t_stddev ^ 2

          prob = ttest(x1, s1, n1, x2, s2, n2, True, False)
          plot_var = 1. - prob
          copy_VarCoords(data1_t_mean, plot_var)

          if (vcoord0_units .eq. "Pa") then
            tmask = gsn_csm_pres_hgt(wks0d, plot_var, tres)
          else
            tmask = gsn_csm_contour(wks0d, plot_var, tres)
          end if

          overlay(maps_d(imod, 0), tmask)
          delete(tmask)
          draw(maps_d(imod, 0))
          frame(wks0d)
        end if

        delete([/x1, x2, s1, s2, prob/])
      end if  ; if t_test

      ; ---------------------------------------------------------------------

    end if  ; if flag_diff

    ; =======================================================================

  end do  ; ii-loop (models)

  ; save default color map in case it is needed later for optionally
  ; plotting color bar to a separate file

  tmp_colors = gsn_retrieve_colormap(wks0)
  cdims = dimsizes(tmp_colors)
  nboxes = dimsizes(res@cnLevels)
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
;  pres@gsnPanelLabelBar     = True    ; add common colorbar
  if (panel_labels) then
    ; print dataset name on each panel
    pres@gsnPanelFigureStrings = names(plotind)
  end if
  pres@gsnPanelFigureStringsJust = "TopRight"
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
  domain = "global"
  plottype = "zonal"

  do is = 0, numseas - 1
    caption = "Zonal mean values for variable " + var0 \
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
    nboxes = dimsizes(dres@cnLevels)
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
        outfile_d(3) = panelling(wks3d, maps_d(plotind, 3), \
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
    ndframe = ndframe + 1

    do is = 0, numseas - 1
      log_info(" Wrote " + outfile(is))

      ; --------------------------------------------------------------------
      ; write provenance to netcdf output and plot file(s) (bias)
      ; --------------------------------------------------------------------

      statistics = (/"clim", "diff"/)
      domain = "global"
      plottype = "zonal"

      ; note: because function log_provenance does not yet support to attach
      ;       different captions to netcdf (contains all seasons) and plots
      ;       (contain one season each), the caption cannot specifiy the
      ;       season plotted; using "annual" or "DJF/MAM/JJA/SON" instead.

      caption = "Differences for zonally averaged variable " + var0 \
                + " (" + allseas + "), reference = " + refname + "."

      if (t_test) then
        caption = caption \
          + " Non-significant grid cells are hatched."
      end if

      log_provenance(nc_outfile_bias, outfile_d(is), caption, statistics, \
                     domain, plottype, "", "", climofiles)
    end do

  end if  ; if flag_diff

  ; optionally save legend(s) to extra file(s)

  if (extralegend) then
    nboxes = dimsizes(res@cnLevels) + 1
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

    labels = tostring(res@cnLevels)

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
      pres@lbFillColors = res@cnFillColors
    else if (isatt(data1, "res_cnFillPalette")) then
      pres@lbFillColors = res@cnFillPalette
    else
      pres@lbFillColors = mean_colors  ; default colors
    end if
    end if

    gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)

    delete(wksleg)
    delete(labels)
    delete(pres@lbFillColors)

    if (flag_diff) then
      nboxes = dimsizes(dres@cnLevels) + 1
      wksleg = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 \
                       + "_diff_legend")

      labels = tostring(dres@cnLevels)

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
        pres@lbFillColors = dres@cnFillColors
      else if (isatt(diff, "res_cnFillPalette")) then
        pres@lbFillColors = dres@cnFillPalette
      else
        pres@lbFillColors = diff_colors  ; default colors
      end if
      end if

      gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)
    end if  ; if (flag_diff)
  end if  ; if (extralegend)

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
