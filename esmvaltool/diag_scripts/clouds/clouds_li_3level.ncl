; CLOUDS_LI_3LEVEL
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates vertical integrals of cloud liquid water / cloud ice in the
;   3 vertical ranges "low", "mid" and "high".
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   explicit_cn_levels: explicit contour levels (array)
;   explicit_cn_dlevels: explicit contour levels (array) for difference plots
;   filename_add:       optionally add this string to plot filesnames
;   showdiff            calculate and plot differences (default = False)
;   rel_diff:           if showdiff = True, then plot relative differences (%)
;                       (default = False)
;   ref_diff_min:       lower cutoff value in case of calculating relative
;                       differences
;                       (in units of input variable)
;   var:                short_name of variable to process (default = "" - use
;                       first variable in variable list)
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
;   20210416-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; define vertical level boundaries for vertical integration

  nlev = 3
  levi = (/70000., 45000., 5000./)     ; level intersections given in Pa
  levlabel = (/"low", "mid", "high"/)  ; level labels for plotting

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "rel_diff", False)
  set_default_att(diag_script_info, "rel_diff_min", -1.0e19)
  set_default_att(diag_script_info, "showdiff", False)
  set_default_att(diag_script_info, "var", "")

  if (diag_script_info@var .eq. "") then
    var0 = variable_info[0]@short_name
  else
    var0 = diag_script_info@var
  end if

  variables = metadata_att_as_array(variable_info, "short_name")
  varidx = ind(variables .eq. var0)
  if (ismissing(varidx)) then
    errstr = "diagnostic " + diag + " requires the following variable: var0"
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)

  if (isatt(variable_info[varidx], "reference_dataset")) then
    refname = variable_info[varidx]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  flag_diff = diag_script_info@showdiff
  flag_rel_diff = diag_script_info@rel_diff
  flag_rel_diff_min = diag_script_info@rel_diff_min

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  if (isatt(diag_script_info, "projection")) then
    projection = diag_script_info@projection
    perim = False
  else
    projection = "CylindricalEquidistant"
    perim = True
  end if

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
    end if
  end if

  climofiles = metadata_att_as_array(info0, "filename")

  if (flag_diff) then
    ; check for reference model definition
    if (.not.isvar("refname")) then
      error_msg("f", DIAG_SCRIPT, "", \
                "no reference dataset defined in recipe")
    end if

    ; set reference model
    if (ref_ind .lt. 0) then
      error_msg("f", DIAG_SCRIPT, "", "cannot calculate differences as " \
                + "reference dataset (" + refname + ") is missing")
    end if
  end if

end

begin
  ; ###########################################
  ; # get data and average time               #
  ; ###########################################

  nplots = dim_MOD

  maps = new((/nplots, nlev/), graphic)
  maps_d = new((/nplots, nlev/), graphic)

  ind_all_sorted = ispan(0, nplots, 1)  ; create array

  if (ref_ind .ge. 0) then
    ind_wo_ref = ind(names .ne. refname)
    ind_all_sorted(0) = ref_ind
    n = dimsizes(names)
    ind_all_sorted(1:n - 1) = ind_wo_ref
  end if

  ; filenames for netcdf output

  nc_filename_bias = work_dir + "clouds_li_3level_" + var0 + "_bias.nc"
  nc_filename_bias@existing = "append"
  nc_filename_mean = work_dir + "clouds_li_3level_" + var0 + "_mean.nc"
  nc_filename_mean@existing = "append"

  res = True

  corr = new(nlev, float)
  gavg = new(nlev, float)
  if (flag_diff) then
    rmsd = new(nlev, float)
    bias = new(nlev, float)
  end if

  do ii = 0, nplots - 1

    imod = ind_all_sorted(ii)
    log_info("processing " + names(imod) + " ***")

    if (isvar("A0")) then
      delete(A0)
    end if

    A0 = read_data(info0[imod])
    ; check dimensions

    dims = getvardims(A0)
    if (dimsizes(dims) .lt. 3) then
      error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + \
                " dimensions, need 3 or 4")
    end if
    idx = ind(dims .eq. "lat")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
    end if
    idx = ind(dims .eq. "lon")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
    end if
    idx = str_match_ind(dims, "lev")
    if (ismissing(idx)) then
      idx = str_match_ind(dims, "air_pressure")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lev dimension")
      end if
    end if

    ; check units of vertical levels

    levs = A0&$dims(idx)$

    if (levs@units .ne. "Pa") then
      error_msg("f", DIAG_SCRIPT, "", "vertical levels have to be in " \
                + "units Pa (currently " + levs@units + "). Use " \
                + "preprocessor to convert units to Pa.")
    else
      ; make sure levs is float
      levs := tofloat(levs)
    end if

    ; make sure levels are monotonic and sorted from high values (near surface)
    ; to low values (top)

    if (.not.(isMonotonic(levs) .eq. -1)) then
      error_msg("f", DIAG_SCRIPT, "", "vertical levels have to decreasing " \
                + "monotonically (from surface to top). Use " \
                + "preprocessor to sort/convert levels.")
    end if

    ; average over time

    data0 = time_operations(A0, -1, -1, "average", "annualclim", True)

    delete(A0)

    if (isvar("data1")) then
      delete(data1)
    end if

    ; calculate vertical integrals over 3 levels "low", "mid", "high"

    data1 = new((/nlev, dimsizes(data0&lat), dimsizes(data0&lon)/), float)
    data1!0 = "xlev"
    data1&xlev = tofloat(ispan(0, nlev - 1, 1))

    nl = dimsizes(levs)
    wgt = new((/nlev, nl/), float)
    wgt = 0.

    do n = 0, nl - 1
      if (n .gt. 0) then
        pbot = 0.5 * (levs(n - 1) + levs(n))
      else
        pbot = levs(n) + 0.5 * (levs(n) - levs(n + 1))
      end if
      if (n .lt. nlev - 1) then
        ptop = 0.5 * (levs(n) + levs(n + 1))
      else
        ptop = levs(n) - 0.5 * (levs(n - 1)  - levs(n))
      end if
      if (pbot .ge. levi(0)) then
        wgt(0, n) = pbot - max((/ptop, levi(0)/))
      end if
      do l = 1, nlev - 1
        if ((ptop .le. levi(l - 1)) .and. (pbot .ge. levi(l))) then
          wgt(l, n) = min((/pbot, levi(l - 1)/)) - max((/ptop, levi(l)/))
        end if
      end do
    end do

    wgt = wgt / 9.806  ; wgt contains layer thickness (in Pa) divided by
                       ; gravity --> corresponds to air mass in layer (kg/m2)

    ; calculate vertical sums (level = low, mid, high)

    ; filter possibly present NaNs
    data0 = where(isnan_ieee(data0), data0@_FillValue, data0)

    do l = 0, nlev - 1
      data1(l, :, :) = dim_sum_wgt_n_Wrap(data0, wgt(l, :), 1, 0)
    end do

    ; clean up

    delete(data0)
    delete(levs)
    delete(wgt)

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/shared/plot/style.ncl

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    res@gsnDraw      = False    ; Don't draw individual plot.
    res@gsnFrame     = False    ; Don't advance frame.

    res@cnFillOn       = True      ; color plot desired
    res@cnLineLabelsOn = False     ; contour lines

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    ; annotation

;    res@tiMainOn             = False

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn            = False

    res@mpOutlineOn          = True
    res@mpFillOn             = False

    res@lbLabelBarOn         = False
    res@gsnRightString       = ""

    res@mpFillDrawOrder       = "PostDraw"    ; draw map last
    res@cnMissingValFillColor = "Gray"

    ; no tickmarks and no labels

    res@tmYLLabelsOn       = False
    res@tmYLOn             = False
    res@tmYRLabelsOn       = False
    res@tmYROn             = False
    res@tmXBLabelsOn       = False
    res@tmXBOn             = False
    res@tmXTLabelsOn       = False
    res@tmXTOn             = False
    res@cnInfoLabelOn      = False    ; turn off cn info label
    res@mpPerimOn          = perim    ; draw line around map

    ; specified in namelist

    res@mpProjection       = projection

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@cnLevels = diag_script_info@explicit_cn_levels
    end if

    if (.not. isatt(res, "cnLevels")) then
      log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
      log_info("info: using default contour levels")
      res@cnLevels = fspan(min(data1), max(data1), 20)
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

    if (isatt(variable_info[varidx], "long_name")) then
      data1@long_name = variable_info[varidx]@long_name
    end if

    data1@var = var0

    if (isatt(variable_info[varidx], "units")) then
      data1@units = variable_info[varidx]@units
    else
      data1@units = ""
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
                + "common grid (check/adjust preprocessor settings in " \
                + "recipe).")
    end if

    corr = corr@_FillValue
    gavg = gavg@_FillValue

    if (.not.all(ismissing(data1))) then
      do il = 0, nlev - 1
        if (same_grid .and. (ref_ind .ge. 0)) then
          corr(il) = calculate_metric(ref_data(il, :, :), data1(il, :, :), \
                                      "correlation")
        end if
        gavg(il) = area_operations(data1(il, :, :), -90., 90., 0., 360., \
                                   "average", True)
      end do
    end if

    res@lbTitleString      = data1@units
    res@lbTitlePosition    = "Bottom"
    res@lbTitleFontHeightF = 0.015

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      wks0 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_li_3level_" + \
                     var0 + filename_add)
      ; difference plots will be saved to a different file
      if (flag_diff) then
        wks0d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_li_3level_" + \
                        var0 + "_bias" + filename_add)
      end if
    end if

    do il = 0, nlev - 1
      if (.not.ismissing(corr(il))) then
        res@gsnRightString = "corr = " + sprintf("%6.3f", corr(il))
      else
        res@gsnRightString = ""
      end if
      if (.not.ismissing(gavg(il))) then
        res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg(il))
      else
        res@gsnLeftString = ""
      end if

      if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
        res@gsnRightString = ""
      end if

      if (il.eq.0) then
        maps(imod, il) = gsn_csm_contour_map(wks0, data1(il, :, :), res)
      end if
      if (il.eq.1) then
        maps(imod, il) = gsn_csm_contour_map(wks0, data1(il, :, :), res)
      end if
      if (il.eq.2) then
        maps(imod, il) = gsn_csm_contour_map(wks0, data1(il, :, :), res)
      end if
    end do

    pres                      = True    ; needed to override
                                        ; panelling defaults
    ; print dataset name on each panel
    pres@gsnPanelLabelBar = True
    pres@gsnPanelFigureStrings = levlabel
    pres@gsnPanelFigureStringsFontHeightF = 0.015

    pres@gsnFrame           = False                ; don't advance panel plot
    pres@gsnPanelLabelBar   = True                 ; add common colorbar
    pres@gsnPanelBottom     = 0.05
    pres@gsnPanelMainString = names(imod)

    gsn_panel(wks0, maps(imod, :), (/3,1/), pres)

    frame(wks0)

    ; mandatory netcdf output

    data1@var = var0 + "_mean_" + names(imod)
    nc_outfile_mean = ncdf_write(data1, nc_filename_mean)

    ; =======================================================================
    ; Create difference plots (if requested)
    ; =======================================================================

    if (flag_diff) then
      dres = True

      if (imod .ne. ref_ind) then
        diff = data1
        if (flag_rel_diff) then
          diff = (diff - ref_data) / ref_data * 100.0
          diff = where(ref_data .le. rel_diff_min, diff@_FillValue, diff)
          diff@units = "%"
        else
          diff = diff - ref_data
        end if

        dres@gsnLeftString  = ""
        dres@gsnRightString = ""

        rmsd = rmsd@_FillValue
        bias = bias@_FillValue

        do il = 0, nlev - 1
          if (.not. flag_rel_diff) then
            if (same_grid) then
              rmsd(il) = calculate_metric(ref_data(il, :, :), \
                                          data1(il, :, :), "RMSD")
            end if
            bias(il) = area_operations(diff(il, :, :), -90., 90., 0., 360., \
                                       "average", True)
          end if
        end do
      else
        continue
      end if

      ; ----------------------------------------------------------------------

      ; ###########################################
      ; # plot ressources                         #
      ; ###########################################

;      dres@gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
;                                           / tofloat((dim_MOD + 1) / 2)/))
;      dres@gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
;                                           / tofloat((dim_MOD + 1) / 2)/))

;      dres@tiMainOn       = False

      dres@cnFillOn       = True      ; color plot desired
      dres@cnLineLabelsOn = False     ; contour lines
      dres@cnLinesOn      = False

      dres@lbLabelBarOn         = False
      dres@lbTitleString      = diff@units
      dres@lbTitlePosition    = "Bottom"
      dres@lbTitleFontHeightF = 0.015

      ; colors
      ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

      ; annotation

      dres@cnLevelSelectionMode = "ExplicitLevels"
      dres@mpOutlineOn          = True
      dres@mpFillOn             = False

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

        if (var0.eq."clw") then
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/qcm3.rgb")
          dres@cnFillColors       = pal
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

      ; map attributes

      dres@mpFillDrawOrder       = "PostDraw"    ; draw map last
      dres@cnMissingValFillColor = "Gray"

      ; no tickmarks and no labels

      dres@tmYLLabelsOn       = False
      dres@tmYLOn             = False
      dres@tmYRLabelsOn       = False
      dres@tmYROn             = False
      dres@tmXBLabelsOn       = False
      dres@tmXBOn             = False
      dres@tmXTLabelsOn       = False
      dres@tmXTOn             = False
      dres@cnInfoLabelOn      = False    ; turn off cn info label

      ; specified in namelist

      dres@mpProjection       = projection

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

      if (isatt(variable_info, "long_name")) then
        diff@long_name = variable_info@long_name
      end if
      if (isatt(variable_info, "units")) then
        diff@units = variable_info@units
      else
        diff@units = ""
      end if

      ; ###########################################
      ; # create the plot                         #
      ; ###########################################

      dres@gsnDraw        = False  ; do not draw yet
      dres@gsnFrame       = False  ; don't advance frame

      ; ----------------------------------------------------------------------

      do il = 0, nlev - 1
        if (.not.ismissing(rmsd(il))) then
          dres@gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(il))
        else
          dres@gsnRightString = ""
        end if
        if (.not.ismissing(bias(il))) then
          dres@gsnLeftString = "bias = " + sprintf("%6.3f", bias(il))
        else
          dres@gsnLeftString = ""
        end if

        if (il.eq.0) then
          maps_d(imod, il) = gsn_csm_contour_map(wks0d, diff(il, :, :), dres)
        end if
        if (il.eq.1) then
          maps_d(imod, il) = gsn_csm_contour_map(wks0d, diff(il, :, :), dres)
        end if
        if (il.eq.2) then
          maps_d(imod, il) = gsn_csm_contour_map(wks0d, diff(il, :, :), dres)
        end if
      end do

      pres@gsnPanelMainString = names(imod) + " - " + refname

      gsn_panel(wks0d, maps_d(imod, :), (/3, 1/), pres)
      frame(wks0d)

      ; mandatory netcdf output

      diff@var = var0 + "_bias_" + names(imod)
      nc_outfile_bias = ncdf_write(diff, nc_filename_bias)

    end if  ; if flag_diff

    ; =======================================================================

  end do  ; ii-loop (models)

  outfile = wks0@fullname

  log_info("Wrote " + outfile)

    ; ------------------------------------------------------------------------
    ; write provenance to netcdf output and plot file(s) (mean)
    ; ------------------------------------------------------------------------

    statistics = (/"clim", "mean"/)
    domain = "global"
    plottype = "geo"

    caption = "Mean values for variable " + var0 \
              + " (" + str_join(levlabel, ",") + ")."
    log_provenance(nc_outfile_mean, outfile, caption, statistics, \
                   domain, plottype, "", "", climofiles)

    ; ========================================================================

    if (flag_diff) then

      outfile_d = wks0d@fullname

      log_info(" Wrote " + outfile_d)

      ; --------------------------------------------------------------------
      ; write provenance to netcdf output and plot file(s) (bias)
      ; --------------------------------------------------------------------

      statistics = (/"clim", "diff"/)
      domain = "global"
      plottype = "geo"

      caption = "Differences for variable " + var0 \
                + " (" + str_join(levlabel, ",") + "), reference = " \
                + refname + "."
      log_provenance(nc_outfile_bias, outfile_d, caption, statistics, \
                     domain, plottype, "", "", climofiles)

    end if  ; if flag_diff

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
