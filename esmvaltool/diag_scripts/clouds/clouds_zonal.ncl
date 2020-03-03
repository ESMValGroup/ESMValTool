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
;   explicit_cn_levels: explicit contour levels (array)
;   extralegend:        plot legend(s) to extra file(s)
;   filename_add:       optionally add this string to plot filesnames
;   panel_labels:       label individual panels (true, false)
;   PanelTop:           manual override for "@gnsPanelTop" used by panel
;                       plot(s)
;   showdiff            calculate and plot differences (default = False)
;   rel_diff:           if showdiff = True, then plot relative differences (%)
;                       (default = False)
;   ref_diff_min:       lower cutoff value in case of calculating relative
;                       differences
;                       (in units of input variable)
;   timemean:           time averaging - "seasonal" = DJF, MAM, JJA, SON),
;                                        "annual" = annual mean
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
;   20200211-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
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
  set_default_att(diag_script_info, "timemean", "annualclim")

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

    ; average over time

    data1 = time_operations(A0, -1, -1, "average", timemean, True)
    delete(A0)

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

    res@tiMainOn             = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn            = False

    res@lbLabelBarOn         = False
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

    if (.not. isvar("ref_data")) then
      ref_data = data1
    end if

    ; check if data are on same grid (for calculating difference, RMSD,
    ; correlation)

    same_grid = False

    if (all(dimsizes(ref_data) .eq. dimsizes(data1))) then
      if (max(abs(ref_data&lat - data1&lat)) .le. 1.0e-6) then
        if (max(abs(ref_data&$vcoord$ - data1&$vcoord$)) .le. 1.0e-6) then
          same_grid = True
        end if
      end if
    end if

    if (flag_diff .and. .not.same_grid) then
      flag_diff = False
      error_msg("f", DIAG_SCRIPT, "", \
                "Data are not on same grid (horizontal and / or vertical), " \
                + "cannot calculate differences. " \
                + "Set showdiff to False in namelist or regrid data to " \
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

    res@gsnDraw        = False  ; do not draw yet
    res@gsnFrame       = False  ; don't advance frame

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
      maps(imod, 0) = gsn_csm_contour(wks0, data1, res)
    end if

    ; mandatory netcdf output

    data1@var = var0 + "_mean_" + names(imod)
;;;    nc_outfile_mean = ncdf_write(data1, nc_filename_mean)

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

      dres@tiMainOn       = False

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

      if (imod.eq.ref_ind) then
        dres@lbLabelBarOn = True
      else
        dres@lbLabelBarOn = False
      end if

      ; map attributes

      dres@cnMissingValFillColor = "Gray"

      dres@cnInfoLabelOn      = False    ; turn off cn info label

      ; set explicit contour levels

      if (isatt(diag_script_info, "explicit_cn_levels")) then
        dres@cnLevelSelectionMode = "ExplicitLevels"
        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        dres@cnLevels = diag_script_info@explicit_cn_levels
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

      dres@gsnDraw        = False  ; do not draw yet
      dres@gsnFrame       = False  ; don't advance frame

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
        maps_d(imod, 0) = gsn_csm_contour(wks0d, diff, dres)
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
;;      log_provenance(nc_outfile_mean, outfile(is), caption, statistics, \
;;                     domain, plottype, "", "", climofiles)
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
  end if  ; if write_plots

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
