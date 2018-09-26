; CLOUDS_BIAS
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates the multi-model mean bias, absolute difference and relative
;   difference of annual mean 2-d cloud variables compared with a
;   reference dataset (observations).
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   plot_abs_diff: additionally also plot absolute differences (True, False)
;   plot_rel_diff: additionally also plot relative differences (True, False)
;   projection:    map projection, e.g., Mollweide, Mercator
;   timemean:      time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
;                  "annualclim" (annual mean)
;
; Required variable_info attributes (variable specific)
;   reference_dataset: name of reference datatset
;
; Optional variable_info attributes (variable specific)
;   long_name: description of variable
;
; Caveats
;   TO DO (ESMValTool v2.0.0)
;   1) get tags from namelist and write meta data to plots
;
; Modification history
;   20180923-A_laue_ax: added writing of results to netcdf
;   20180914-A_laue_ax: code rewritten for ESMValTool v2.0
;   20170620-A_laue_ax: added tags for reporting
;   20160901-A_laue_ax: added regridding option 1 deg x 1 deg
;   20151027-A_laue_ax: moved call to 'write_references' to the beginning
;                       of the code
;   20150428-A-laue_ax: written.
;
; ############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/contour_maps.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; ========================================================================
  ; ========================= write references =============================
  ; ========================================================================

  write_references(diag_script,  \   ; script name
                   "A_laue_ax",  \   ; authors
                   "",  \            ; contributors
                   "",  \            ; diag_references
                   "",  \            ; obs_references
                   (/"P_embrace"/))  ; proj_references

  dim_MOD = dimsizes(input_file_info@dataset)
  dim_VAR = 1
  var0 = variables(0)
  field_type0 = field_types(0)

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "plot_abs_diff", False)
  set_default_att(diag_script_info, "plot_rel_diff", False)
  set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
  set_default_att(diag_script_info, "timemean", "annualclim")

  timemean = diag_script_info@timemean

  if (timemean.eq."seasonalclim") then
    numseas = 4
    season = (/"DJF", "MAM", "JJA", "SON"/)
  else
    numseas = 1  ; default
    season = (/"annual"/)
  end if

  ; make sure path for (optional) netcdf output exists

  if (config_user_info@write_netcdf.eq."True") then
    write_nc = True
    work_dir = config_user_info@work_dir + "/"
    ; Create work dir
    system("mkdir -p " + work_dir)
  else
    write_nc = False
  end if

  if (config_user_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if

end

begin
  ; ========================================================================
  ; ========================== initialization ==============================
  ; ========================================================================

  ; check for reference dataset definition

  if (.not.(isatt(variable_info, "reference_dataset"))) then
    error_msg("f", diag_script, "", "no reference dataset is selected")
  end if

  plot_abs_diff = diag_script_info@plot_abs_diff
  plot_rel_diff = diag_script_info@plot_rel_diff

  ; get reference dataset

  ref_ind = ind(input_file_info@dataset .eq. variable_info@reference_dataset)
  ref_model = variable_info@reference_dataset
  if (ismissing(ref_ind)) then
    error_msg("f", diag_script, "", "reference model (" \
              + variable_info@reference_dataset + ") is missing")
  end if

  ; get multi-model mean index

  mm_ind = ind(input_file_info@dataset .eq. "MultiModelMean")

  if (ismissing(mm_ind)) then
    error_msg("f", diag_script, "", "multi-model mean is missing (required)")
  end if

  ; basename of diag_script

  diag_script_base = basename(diag_script)

  ; ========================================================================
  ; =========================== calculations ===============================
  ; ========================================================================

  ; note: 1) masking is handled by the backend
  ;       2) multi-model mean is calculated by the backend

  ; read data

  A0 = read_data(mm_ind, variables(0), field_types(0))
  mmdata = time_operations(A0, -1, -1, "average", timemean, True)
  delete(A0)

  A0 = read_data(ref_ind, variables(0), field_types(0))
  refdata = time_operations(A0, -1, -1, "average", timemean, True)
  delete(A0)

  ; differences between multi-model mean and reference data set

  diff = mmdata - refdata
  mmdata@diag_script = diag_script
  copy_VarMeta(mmdata, diff)

  delete(refdata)

  ; convert units for plotting (precipitation: kg m-2 s-1 --> mm day-1)

  if (var0.eq."pr") then
    mmdata = convert_units(mmdata, "mm/day")
    diff = convert_units(diff, "mm/day")
    ; for plotting IPCC AR5 fig 9.4
    diff@res_cnLevels = ispan(-30, 30, 5) * 0.1
  end if

  if (var0.eq."tas") then
    mmdata = convert_units(mmdata, "degC")
    diff@units = "degC"
    ; for plotting IPCC AR5 fig 9.2
    diff@res_cnLevels = ispan(-5, 5, 1)
  end if

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  climofiles = new(2, string)
  climofiles(0) = input_file_info@filename(mm_ind)
  climofiles(1) = input_file_info@filename(ref_ind)

;  alltags = array_append_record(tags, (/"PT_geo", "ST_clim", "ST_diff", \
;                                        "DM_global"/), 0)

  diff@res_gsnMaximize  = True  ; use full page for the plot
  diff@res_cnFillOn  = True  ; color plot desired
  diff@res_cnLineLabelsOn  = False  ; contour lines
  diff@res_cnLinesOn  = False
  diff@res_tiMainOn  = True
  diff@res_gsnLeftStringFontHeightF = 0.015
  diff@res_cnLevelSelectionMode  = "ExplicitLevels"
  diff@res_mpOutlineOn  = True
  if (.not.isatt(diff, "res_cnLevels")) then
    diff@res_cnLevels = fspan(min(diff), max(diff), 20)
  end if
  diff@res_mpFillOn  = False
  diff@res_lbLabelBarOn  = True
  diff@res_gsnRightString  = ""
  diff@res_mpFillDrawOrder  = "PostDraw"  ; draw map fill last
  diff@res_cnMissingValFillColor  = "Gray"
  diff@res_tmYLLabelsOn  = False
  diff@res_tmYLOn  = False
  diff@res_tmYRLabelsOn  = False
  diff@res_tmYROn  = False
  diff@res_tmXBLabelsOn  = False
  diff@res_tmXBOn  = False
  diff@res_tmXTLabelsOn  = False
  diff@res_tmXTOn  = False
;  diff@res_gsnRightString  = variable_info@long_name + " in " \
;                             + diff@units
  diff@res_cnInfoLabelOn  = False  ; turn off cn info label
  diff@res_mpProjection  = diag_script_info@projection

  diff@var = var0  ; Overwrite existing entry
  if (isatt(variable_info, "long_name")) then
  diff@var_long_name = variable_info@long_name
  end if
  diff@var_units  = diff@units

  plots = new((/4, numseas/), graphic)

  ; --------------------------------------------------------------------
  ; plot contour map

  diff@res_gsnDraw  = False  ; Do not draw yet
  diff@res_gsnFrame  = False  ; Don't advance frame.
  diff@res_mpPerimOn  = False

;  diff@res_lbTitleString  = "~F8~D~F21~" + variable_info@long_name + \
;                            " (" + diff@units + ")"
  diff@res_lbTitleString  = "(" + diff@units + ")"
  diff@res_lbTitlePosition  = "Bottom"

  diff@res_lbLabelFontHeightF  = 0.014
;  diff@res_lbAutoManage  = False
  diff@res_lbTopMarginF  = 0.1
  diff@res_lbTitleFontHeightF  = 0.014
;  diff@res_pmLabelBarWidthF  = 0.8

  diff@res_tiMainFontHeightF  = 0.016

  diff@res_tiMainString = "Multi Model Mean Bias"

  copy_VarMeta(diff, mmdata)
  delete(mmdata@res_cnLevels)

  if (var0.eq."tas") then
    mmdata@res_cnLevels = ispan(-30, 30, 3)

    pal = read_colormap_file("diag_scripts/shared/plot/rgb/" \
                             + "ipcc-tas.rgb")
    mmdata@res_cnFillColors = pal
    pal2 = read_colormap_file("diag_scripts/shared/plot/rgb/" \
                              + "ipcc-tas-delta.rgb")
    diff@res_cnFillColors = pal2
    if (plot_abs_diff) then
      pal3 = read_colormap_file("diag_scripts/shared/plot/rgb/" \
                                + "ipcc-tas-absdelta.rgb")
    end if
    else if (var0.eq."pr") then
      mmdata@res_cnLevels = ispan(10, 100, 15) * 0.1

      pal = read_colormap_file("diag_scripts/shared/plot/rgb/ipcc-precip.rgb")
      mmdata@res_cnFillColors = pal
      pal2 = read_colormap_file("diag_scripts/shared/plot/rgb/" \
                                + "ipcc-precip-delta.rgb")
      diff@res_cnFillColors = pal2
      if (plot_abs_diff) then
        pal3 = read_colormap_file("diag_scripts/shared/plot/rgb/" \
                                  + "ipcc-precip-absdelta.rgb")
      end if
      if (plot_rel_diff) then
        pal4 = read_colormap_file("diag_scripts/shared/plot/rgb/" \
                                  + "ipcc-precip-reldelta.rgb")
      end if
    else
      mmdata@res_cnLevels = fspan(min(mmdata), max(mmdata), 20)
    end if
  end if

  mmdata@res_tiMainString = "Multi Model Mean"

  plotsperline = (/2, 0/)
  plotind = (/0, 1/)  ; mmm and mean bias are always plotted

  ; absolute differences

  if (plot_abs_diff) then
    absdiff = abs(diff)
    copy_VarMeta(diff, absdiff)
    if (isvar("pal3")) then
      delete(absdiff@res_cnFillColors)
      absdiff@res_cnFillColors = pal3
      nt = dimsizes(pal3)
      n = nt(0)
    else
      n = dimsizes(diff@res_cnLevels)
    end if
    mx = max(diff@res_cnLevels)
    delete(absdiff@res_cnLevels)

    tmp = fspan(0.0, mx, n)

    absdiff@res_cnLevels = tmp(1:dimsizes(tmp)-1)
    delete(tmp)

    absdiff@res_tiMainString = "Multi Model Mean of Absolute Error"

    iadd = 2
    itmp = array_append_record(plotind, iadd, 0)
    delete(plotind)
    plotind = itmp
    delete(itmp)
    plotsperline(1) = plotsperline(1) + 1
  end if

  ; relative differences

  if (plot_rel_diff) then
    ; replace "epsilon" values with missing value
    mm = where(abs(mmdata).lt.1.0e-6, mmdata@_FillValue, mmdata)

    reldiff = 100.0 * diff / mm

    copy_VarMeta(diff, reldiff)
    delete(reldiff@res_cnLevels)
    reldiff@res_cnLevels = fspan(-90.0, 90.0, 13)
    reldiff@res_tiMainString = "Multi Model Mean of Relative Error"
    reldiff@units = "%"
    reldiff@res_lbTitleString = "(" + reldiff@units + ")"
    if (isvar("pal4")) then
      delete(reldiff@res_cnFillColors)
      reldiff@res_cnFillColors = pal4
    end if

    iadd = 3
    itmp = array_append_record(plotind, iadd, 0)
    delete(plotind)
    plotind = itmp
    delete(itmp)
    plotsperline(1) = plotsperline(1) + 1
  end if

  ; panelling resources

  pres  = True
  pres@gsnPanelCenter  = False
  pres@gsnPanelRowSpec = True  ; tell panel what order to plot
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5

  if (write_plots) then
    do is = 0, numseas - 1
      ; --------------------------------------------------------------------
      ; create workspace

      if (isvar("wks")) then
        delete(wks)
      end if

      wks = get_wks("dummy_for_wks", diag_script, "clouds_bias_" + var0 \
                    + "_" + season(is))

      if (numseas.gt.1) then
        pres@txString = season(is)
        plots(0, is) = contour_map(wks, mmdata(is, :, :), var0)
        plots(1, is) = contour_map(wks, diff(is, :, :), var0)
        if (plot_abs_diff) then
          plots(2, is) = contour_map(wks, absdiff(is, :, :), var0)
        end if
        if (plot_rel_diff) then
          plots(3, is) = contour_map(wks, reldiff(is, :, :), var0)
        end if
        gsn_panel(wks, plots(plotind, is), plotsperline, pres)
      else
        plots(0, 0) = contour_map(wks, mmdata, var0)
        plots(1, 0) = contour_map(wks, diff, var0)
        if (plot_abs_diff) then
          plots(2, 0) = contour_map(wks, absdiff, var0)
        end if
        if (plot_rel_diff) then
          plots(3, 0) = contour_map(wks, reldiff, var0)
        end if
        gsn_panel(wks, plots(plotind, 0), plotsperline, pres)
      end if

      ; add meta data to plot (for reporting)

      caption = "Multi model values, from top left to bottom right: " \
                + "mean, bias"
      if (plot_abs_diff) then
        caption = caption + ", absolute error"
      end if
      if (plot_rel_diff) then
        caption = caption + ", relative error"
      end if
      caption = caption + " for variable " + variables(0) \
        + " (" + season(is) + "), reference = " + \
        input_file_info@dataset(ref_ind) + "."
      id = diag_script + "_" + variables(0) + "_" + season(is)

      contrib_authors = "A_laue_ax"

;      ESMValMD(wks@fullname, alltags, caption, id, variables(0), \
;               input_file_info@dataset, climofiles, diag_script, \
;               contrib_authors)

    end do  ; is-loop (seasons)
  end if  ; if write_plots

  ; ###########################################
  ; # Optional output to netCDF  #
  ; ###########################################

  if (write_nc) then
    nc_filename = work_dir + "clouds_bias_mmm_" + var0 + ".nc"
    nc_outfile = ncdf_write(mmdata, nc_filename)
    nc_filename = work_dir + "clouds_bias_diff_" + var0 + ".nc"
    nc_outfile = ncdf_write(diff, nc_filename)
    if (isvar("absdiff")) then
      nc_filename = work_dir + "clouds_bias_abs_diff_" + var0 + ".nc"
      nc_outfile = ncdf_write(absdiff, nc_filename)
    end if
    if (isvar("reldiff")) then
      nc_filename = work_dir + "clouds_bias_rel_diff_" + var0 + ".nc"
      nc_outfile = ncdf_write(reldiff, nc_filename)
    end if
  end if

  leave_msg(diag_script, "")

end
