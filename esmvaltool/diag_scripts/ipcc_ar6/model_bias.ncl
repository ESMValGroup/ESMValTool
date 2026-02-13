; MODEL_BIAS
; ############################################################################
; Author: Lisa Bock (DLR, Germany), Axel Lauer (DLR, Germany)
; PROJECT-NAME IPCC-AR6
; ############################################################################
; Description
;   Calculates the multi-model mean bias, absolute difference, relative
;   difference and root mean square difference of annual mean 2-d cloud
;   variables compared with a reference dataset (observations).
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   plot_abs_diff: additionally also plot absolute differences (True, False)
;   plot_rel_diff: additionally also plot relative differences (True, False)
;   plot_rms_diff: additionally also plot root mean square differences
;                  (True, False)
;   projection:    map projection, e.g., Mollweide, Mercator
;
; Required variable_info attributes (variable specific)
;   reference_dataset: name of reference datatset
;
; Optional variable_info attributes (variable specific)
;   long_name: description of variable
;
; Caveats
;   none
;
; Modification history
;   20190312-bock_lisa: added calculation of RMSD and IPCC AR6 color scheme
;   20190312-bock_lisa: adapted code (based on clouds/clouds_bias.ncl)
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/scaling.ncl"
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
  infiles = metadata_att_as_array(info0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "plot_abs_diff", False)
  set_default_att(diag_script_info, "plot_rms_diff", False)
  set_default_att(diag_script_info, "plot_rel_diff", False)
  set_default_att(diag_script_info, "stippandhatch", False)
  set_default_att(diag_script_info, "projection", "Robinson")

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  input_dir = diag_script_info@input_files(1)

end

begin
  ; ========================================================================
  ; ========================== initialization ==============================
  ; ========================================================================

  ; check for reference dataset definition

  if (.not.(isvar("refname"))) then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
  end if

  plot_abs_diff = diag_script_info@plot_abs_diff
  plot_rel_diff = diag_script_info@plot_rel_diff
  plot_rms_diff = diag_script_info@plot_rms_diff

  stippandhatch = diag_script_info@stippandhatch

  ; get reference dataset

  ref_ind = ind(names .eq. refname)
  if (ismissing(ref_ind)) then
    error_msg("f", DIAG_SCRIPT, "", "reference dataset (" \
              + refname + ") is missing")
  end if

  ; get multi-model mean index

  mm_ind = ind(names .eq. "MultiModelMean")

  if (ismissing(mm_ind)) then
    error_msg("f", DIAG_SCRIPT, "", "multi-model mean is missing (required)")
  end if

  ; basename of diag_script

  diag_script_base = basename(DIAG_SCRIPT)

  ; ========================================================================
  ; =========================== Read data ==================================
  ; ========================================================================

  ; note: 1) masking is handled by the backend
  ;       2) multi-model mean is calculated by the backend

  ; read data

  A0 = read_data(info0[mm_ind])
  mmdata = time_operations(A0, -1, -1, "average", "annualclim", True)
  if (dimsizes(A0&time).gt.12) then
    mm_mon = time_operations(A0, -1, -1, "average", "monthlyclim", True)
  else
    mm_mon = A0
  end if
  delete(A0)

  A0 = read_data(info0[ref_ind])
  refdata = time_operations(A0, -1, -1, "average", "annualclim", True)
  if (dimsizes(A0&time).gt.12) then
    ref_mon = time_operations(A0, -1, -1, "average", "monthlyclim", True)
  else
    ref_mon = A0
  end if
  delete(A0)

  dim = dimsizes(mmdata)
  data = new((/dim_MOD-2, dim(0), dim(1)/), float)
  data_mon = new((/dim_MOD-2, 12, dim(0), dim(1)/), float)

  ii = 0

  do imod = 0, dim_MOD - 1

    if(imod .ne. mm_ind .and. imod .ne. ref_ind) then
      log_info("processing " + names(imod))

      if (isvar("data1")) then
        delete(data1)
      end if

      if (isvar("A0")) then
        delete(A0)
      end if

      A0 = read_data(info0[imod])

      data(ii, :, :) = time_operations(A0, -1, -1, "average", "annualclim", \
                                       True)
      data_mon(ii, :, :, :) = time_operations(A0, -1, -1, "average", \
                                              "monthlyclim", True)

      ii = ii + 1

      delete(A0)

    end if

  end do

  ; -----------------------------------------------------------
  ; read precalculated natural variability for multi-model-mean
  ; -----------------------------------------------------------

  if (stippandhatch) then
    log_info("Read precalculated natural variability for multi-model")
    natvar_file = input_dir + "/" + "IAV_mmm_piControl_" + \
                              var0 + "_annualclim.nc"
    natvar = addfile(natvar_file, "r")
    natvar_mmm = natvar->iav
    tmp = natvar_mmm
    delete(natvar_mmm)
    natvar_mmm = rm_single_dims(tmp)
    delete(tmp)
  end if

  ; -----------------------------------------------------------

  ; convert units for plotting
  if (var0.eq."pr") then
    ; kg m-2 s-1 --> mm day-1
    mmdata = convert_units(mmdata, "mm/day")
    refdata = convert_units(refdata, "mm/day")
    data = convert_units(data, "mm/day")
    ref_mon = convert_units(ref_mon, "mm/day")
    data_mon = convert_units(data_mon, "mm/day")
  elseif (var0.eq."tas") then
    mmdata = convert_units(mmdata, "degC")
    data = convert_units(data, "degC")
    refdata = convert_units(refdata, "degC")
  elseif (var0.eq."tos" .or. var0.eq."thetao") then
    if (mmdata@units .eq. "K") then
      mmdata = convert_units(mmdata, "degC")
      data_mon = convert_units(data_mon, "degC")
    end if
    if (refdata@units .eq. "K") then
      refdata = convert_units(refdata, "degC")
      ref_mon = convert_units(ref_mon, "degC")
    end if
  elseif (var0.eq."sos") then
    refdata = refdata * 1000.
  end if

  ; ========================================================================
  ; =========================== Calculations ===============================
  ; ========================================================================

  ; differences between multi-model mean and reference data set
  ; (multi-model bias)

  diff = mmdata - refdata
  mmdata@diag_script = DIAG_SCRIPT
  copy_VarMeta(mmdata, diff)

  if any(var0 .eq. (/"tas", "tos", "thetao"/)) then
    diff@units = "~F34~0~F~ C"
    diff@res_cnLevels = ispan(-6, 6, 1)
  elseif (var0 .eq. "pr") then
    diff@res_cnLevels = ispan(-30, 30, 5) * 0.1
  elseif any(var0 .eq. (/"sos", "so"/)) then
    diff@units = "PSS-78"
    diff@res_cnLevels = ispan(-30, 30, 5) * 0.1
  end if

  diff_i    = new((/dim_MOD-2, dim(0), dim(1)/), float)
  absdiff_i = new((/dim_MOD-2, dim(0), dim(1)/), float)
  rmsdiff_i = new((/dim_MOD-2, dim(0), dim(1)/), float)
  reldiff_i = new((/dim_MOD-2, dim(0), dim(1)/), float)

  if (stippandhatch .or. plot_rms_diff .or. plot_abs_diff \
      .or. plot_rel_diff) then

    ii = 0

    do imod = 0, dim_MOD - 1

      if(imod .ne. mm_ind .and. imod .ne. ref_ind) then
        log_info("processing " + names(imod))

        ; difference
        diff_i(ii, :, :) = data(ii, :, :) - refdata

        ; absolute difference
        absdiff_i(ii, :, :) = abs(diff_i(ii, :, :))

        ; root mean square error
        if (plot_rms_diff) then
          rmsdiff_i(ii, :, :) = calculate_metric(data_mon(ii, :, :, :), \
                                                 ref_mon, "RMSDxy")
        end if

        ; relative difference
        ; replace "epsilon" values with missing value
        datam = where(abs(data) .lt. 1.e-6, data@_FillValue, data)
        reldiff_i(ii, :, :) = 100.0 * (diff_i(ii, :, :) / datam(ii, :, :))

        ii = ii + 1

      end if

    end do

  end if

  ; multi-model mean of absolute difference between models and reference
  ; data set

  absdiff = dim_avg_n_Wrap(absdiff_i, 0)
  copy_VarMeta(diff, absdiff)

  rmsdiff = dim_avg_n_Wrap(rmsdiff_i, 0)
  copy_VarMeta(diff, rmsdiff)

  reldiff = dim_avg_n_Wrap(reldiff_i, 0)
  copy_VarMeta(diff, reldiff)

  rmsd = calculate_metric(refdata, mmdata, "RMSD")
  bias = area_operations(diff, -90., 90., 0., 360., "average", True)
  mean = area_operations(mmdata, -90., 90., 0., 360., "average", True)

  ; ----------------------------------------------
  ; Determine sigma and sign for significance
  ; ----------------------------------------------
  if (stippandhatch) then
    var_pos_signif = dim_num_n(where(diff_i .gt. 0., 1., \
                                     diff@_FillValue) .eq. 1., 0)
    var_neg_signif = dim_num_n(where(diff_i .lt. 0., 1., \
                                     diff@_FillValue) .eq. 1., 0)
    var_pos_signif_tmp = \
      var_pos_signif / where((var_pos_signif + var_neg_signif) \
                             .ne. 0., var_pos_signif + \
                             var_neg_signif, diff@_FillValue)
    var_neg_signif_tmp = \
      var_neg_signif / where((var_pos_signif + var_neg_signif) \
                             .ne. 0., var_pos_signif + \
                             var_neg_signif, diff@_FillValue)
    signif_and_sign = where(var_pos_signif_tmp .ge. 0.8 .or. \
                            var_neg_signif_tmp .ge. 0.8, 1., 0.)
    signif_and_sign = where(ismissing(signif_and_sign), 0.,\
                            signif_and_sign)
    ; variablity threshold: sqrt(2)*1.645=~2.3264
    var_thres = abs(2.3264 * natvar_mmm)
    tmp := new(dimsizes(diff_i), "float")
    tmp = 0
    do imod = 0, dim_MOD - 3
      tmp(imod, :, :) = where(abs(diff_i(imod, :, :)) .gt. var_thres(:, :), \
                              1, 0)
    end do
    var_pos_signif := dim_sum_n(tmp, 0)
    delete(tmp)
    var_pos_signif_tmp := var_pos_signif / (dim_MOD - 2)
    sigma_and_sign = where(var_pos_signif_tmp .ge. 0.66, 1., 0.)

    to_plot_not_robust = where(sigma_and_sign .eq. 0, 1., 0.)

    to_plot_confl_sign = where(signif_and_sign .eq. 0, sigma_and_sign, 0.)

    copy_VarMeta(diff, to_plot_confl_sign)
    copy_VarMeta(diff, to_plot_not_robust)

    if (plot_rms_diff) then
      to_plot_not_signif_rmsd = where(rmsdiff(:, :) .lt. \
                                      abs(natvar_mmm), 1., 0.)
      copy_VarMeta(diff, to_plot_not_signif_rmsd)
    end if
  end if

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  diff@res_gsnMaximize  = True  ; use full page for the plot
  diff@res_cnFillOn  = True  ; color plot desired
  diff@res_cnLineLabelsOn  = False  ; contour lines
  diff@res_cnLinesOn  = False
  diff@res_tiMainOn  = True
  diff@res_gsnLeftStringFontHeightF = 0.015
  diff@res_gsnRightStringFontHeightF = 0.015
  diff@res_cnLevelSelectionMode  = "ExplicitLevels"
  diff@res_mpOutlineOn  = True
  if (.not.isatt(diff, "res_cnLevels")) then
    diff@res_cnLevels = fspan(min(diff), max(diff), 20)
  end if
  diff@res_mpFillOn  = False
  diff@res_lbLabelBarOn  = True
  diff@res_gsnRightString  = ""
  diff@res_gsnLeftString   = ""
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
  diff@res_cnInfoLabelOn  = False  ; turn off cn info label
  diff@res_mpProjection  = diag_script_info@projection

  if (var0 .eq. "tos" .or. var0 .eq. "sos") then
    diff@res_mpCenterLonF = 180
  end if

  diff@var = var0  ; Overwrite existing entry
  if (isatt(variable_info[0], "long_name")) then
    diff@var_long_name = variable_info[0]@long_name
  end if
  diff@var_units  = diff@units

  plots = new((/5/), graphic)

  ; --------------------------------------------------------------------
  ; plot contour map

  diff@res_gsnDraw  = False  ; Do not draw yet
  diff@res_gsnFrame  = False  ; Don't advance frame.
  diff@res_mpPerimOn  = False

  diff@res_mpGeophysicalLineColor = "gray42"
  diff@res_mpPerimOn    = False
  diff@res_mpGridLineColor   = "transparent"
  diff@res_mpGridAndLimbOn   = True
  diff@res_mpOutlineOn = True

  diff@res_lbTitleString  = "(" + diff@units + ")"
  diff@res_lbTitlePosition  = "Bottom"

  diff@res_lbLabelFontHeightF  = 0.014
  diff@res_lbTopMarginF  = 0.1
  diff@res_lbTitleFontHeightF  = 0.014

  diff@res_lbBoxEndCapStyle = "TriangleBothEnds"
  diff@res_lbBoxLineColor = "gray10"

  diff@res_tiMainFontHeightF  = 0.016

  diff@res_gsnRightStringParallelPosF = 0.8

  if(isatt(diag_script_info, "bias_caption")) then
    diff@res_tiMainString = diag_script_info@bias_caption
  else
    diff@res_tiMainString = "Multi-Model Mean Bias"
  end if

  copy_VarMeta(diff, mmdata)
  delete(mmdata@res_cnLevels)

  if(isatt(diag_script_info, "mean_caption")) then
    mmdata@res_tiMainString = diag_script_info@mean_caption
  else
    mmdata@res_tiMainString = "Multi-Model Mean"
  end if

  mmdata@res_gsnLeftStringFontHeightF = 0.015
  mmdata@res_gsnRightStringFontHeightF = 0.015
  mmdata@res_gsnRightString  = " "
  mmdata@res_gsnLeftString   = " "

  if (var0 .eq. "tas") then
    mmdata@res_cnLevels = ispan(-35, 35, 5)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_temperature_div.rgb")
    mmdata@res_cnFillPalette = pal
    diff@res_cnFillPalette = pal
    if (plot_abs_diff) then
      pal3 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_temperature_div.rgb")
    end if
    if (plot_rel_diff) then
      pal4 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_temperature_10.rgb")
    end if
    if (plot_rms_diff) then
      pal5 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_temperature_seq.rgb")
    end if
  elseif (var0.eq."pr") then
    mmdata@res_cnLevels = ispan(1, 10, 1)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_precipitation_seq.rgb")
    mmdata@res_cnFillPalette = pal
    mmdata@res_lbBoxEndCapStyle = "TriangleHighEnd"
    pal2 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                              + "ipcc-ar6_precipitation_div.rgb")
    diff@res_cnFillPalette = pal2
    if (plot_abs_diff) then
      pal3 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_precipitation_10.rgb")
    end if
    if (plot_rel_diff) then
      pal4 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_precipitation_10.rgb")
    end if
    if (plot_rms_diff) then
      pal5 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_precipitation_seq.rgb")
    end if
  elseif (var0.eq."tos" .or. var0.eq."thetao") then
    mmdata@res_cnLevels = ispan(-30, 30, 5)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_temperature_div.rgb")
    mmdata@res_cnFillPalette = pal
    diff@res_cnFillPalette = pal
    if (plot_abs_diff) then
      pal3 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_temperature_10.rgb")
    end if
    if (plot_rel_diff) then
      pal4 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_temperature_10.rgb")
    end if
    if (plot_rms_diff) then
      pal5 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_temperature_seq.rgb")
    end if
  elseif (var0.eq."sos" .or. var0 .eq. "so") then
    mmdata@res_cnLevels = ispan(26, 37, 1)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_misc_seq_1.rgb")
    mmdata@res_cnFillPalette = pal
    pal2 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                              + "ipcc-ar6_misc_div.rgb")
    diff@res_cnFillPalette = pal2
  else
    mmdata@res_cnLevels = fspan(min(mmdata), max(mmdata), 20)
  end if

  plotsperline = (/2, 0/)
  plotind = (/0, 1/)  ; mmm and mean bias are always plotted

  ; absolute differences

  if (plot_abs_diff) then
    copy_VarMeta(diff, absdiff)

    absdiff@res_tiMainString = "Multi Model Mean of Absolute Error"

    iadd = 2
    itmp = array_append_record(plotind, iadd, 0)
    delete(plotind)
    plotind = itmp
    delete(itmp)
    plotsperline(1) = plotsperline(1) + 1
  end if

  ; root mean square differences

  if (plot_rms_diff) then
    copy_VarMeta(diff, rmsdiff)
    if (isvar("pal5")) then
      rmsdiff@res_cnFillPalette = pal5
      rmsdiff@res_lbBoxEndCapStyle = "TriangleHighEnd"
    end if
    if (var0 .eq. "tas") then
      delete(rmsdiff@res_cnLevels)
      rmsdiff@res_cnLevels = ispan(1, 12, 1) * 0.5
    elseif (var0 .eq. "pr")
      delete(rmsdiff@res_cnLevels)
      rmsdiff@res_cnLevels = ispan(1, 8, 1) * 0.5
    end if

    if(isatt(diag_script_info, "rmsd_caption")) then
      rmsdiff@res_tiMainString = diag_script_info@rmsd_caption
    else
      rmsdiff@res_tiMainString = "Multi Model Mean of Root Mean Square Error"
    end if

    iadd = 3
    itmp = array_append_record(plotind, iadd, 0)
    delete(plotind)
    plotind = itmp
    delete(itmp)
    plotsperline(1) = plotsperline(1) + 1
  end if

  ; relative differences

  if (plot_rel_diff) then
    copy_VarMeta(diff, reldiff)
    delete(reldiff@res_cnLevels)
    reldiff@res_cnLevels = fspan(-80.0, 80.0, 9)
    reldiff@res_tiMainString = "Multi-Model Mean of Relative Error"
    reldiff@units = "%"
    reldiff@res_lbTitleString = "(" + reldiff@units + ")"
    if (isvar("pal4")) then
      delete(reldiff@res_cnFillColors)
      reldiff@res_cnFillColors = pal4
    end if

    iadd = 4
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
  if (isatt(diag_script_info, "panel_num")) then
    pres@gsnPanelFigureStrings = diag_script_info@panel_num
  else
    pres@gsnPanelFigureStrings = (/"a)", "b)", "c)", "d)"/)
  end if
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  pres@gsnPanelFigureStringsJust = "TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.016

  if (isatt(diag_script_info, "panel_title")) then
    pres@gsnPanelMainString = diag_script_info@panel_title
  end if

  ; --------------------------------------------------------------------
  ; create workspace

  if (isvar("wks")) then
    delete(wks)
  end if

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "model_bias_" + var0 \
                + "_annualclim_" + input_file_info[1]@project)

  plotfile = wks@fullname

  plots(0) = contour_map(wks, mmdata, var0)
  plots(1) = contour_map(wks, diff, var0)
  if (plot_abs_diff) then
    plots(2) = contour_map(wks, absdiff, var0)
  end if
  if (plot_rms_diff) then
    plots(3) = contour_map(wks, rmsdiff, var0)
  end if
  if (plot_rel_diff) then
    plots(4) = contour_map(wks, reldiff, var0)
  end if

  ; add significance
  if (stippandhatch) then
    ; resources cross-hatching
    ress = True
    ress@gsnDraw = False
    ress@gsnFrame = False
    ress@gsnAddCyclic = True
    ress@cnLinesOn = False
    ress@cnLineLabelsOn = False
    ress@cnLevelSelectionMode = "ExplicitLevels"
    ress@cnLevels = (/.3, 1.3/)
    ress@cnInfoLabelOn = False
    ress@tiMainOn = False
    ress@cnFillOn = True
    ress@cnFillMode = "AreaFill"
    ress@cnMonoFillPattern = False
    ress@cnFillPatterns = (/-1, 6, 6/)
    ress@cnMonoFillColor = True
    ress@lbLabelBarOn = False
    ress@gsnLeftString = ""
    ress@gsnRightString = ""
    ress@cnFillDotSizeF       = 0.001
    ress@cnFillScaleF        = 1.        ; add extra density

    ; plot cross-hatching
    if (max(to_plot_confl_sign(:, :)) .gt. \
        min(to_plot_confl_sign(:, :))) then
      plot1 = gsn_csm_contour(wks, to_plot_confl_sign(:, :), ress)
      overlay(plots(1), plot1)
    end if

    ; resources hatching
    resb = True
    resb@gsnDraw = False
    resb@gsnFrame = False
    resb@gsnAddCyclic = True
    resb@cnLinesOn = False
    resb@cnLineLabelsOn = False
    resb@cnLevelSelectionMode = "ExplicitLevels"
    resb@cnLevels = (/.3, 1.3/)
    resb@cnInfoLabelOn = False
    resb@tiMainOn = False
    resb@cnFillOn = True
    resb@cnFillMode = "AreaFill"
    resb@cnRasterSmoothingOn = False
    resb@cnMonoFillPattern = False
    resb@cnFillPatterns = (/-1, 4, 4/)
    resb@cnMonoFillColor = True
    resb@lbLabelBarOn = False
    resb@gsnLeftString = ""
    resb@gsnRightString = ""
    resb@cnFillScaleF = 0.5

    ; plot hatching
    if (max(to_plot_not_robust(:, :)) .gt. \
        min(to_plot_not_robust(:, :))) then
      plot2 = gsn_csm_contour(wks, to_plot_not_robust(:, :), resb)
      overlay(plots(1), plot2)
    end if

  end if

  gsn_panel(wks, plots(plotind), plotsperline, pres)

  ; add meta data to plot (for reporting)

  caption = "Multi model values, from top left to bottom right: " \
            + "mean, bias"
  if (plot_abs_diff) then
    caption = caption + ", absolute error"
  end if
  if (plot_rms_diff) then
    caption = caption + ", root mean square error"
  end if
  if (plot_rel_diff) then
    caption = caption + ", relative error"
  end if

  ; ###########################################
  ;  output to netCDF
  ; ###########################################

  nc_filename = work_dir + "model_bias_" + var0 + "_annualclim_" \
                         + input_file_info[1]@project + ".nc"

  nc_filename@existing = "overwrite"

  mmdata@var = var0 + "_mean"
  mmdata@long_name = var0 + " (multi-model mean)"
  nc_outfile = ncdf_write(mmdata, nc_filename)

  nc_filename@existing = "append"

  refdata@var = var0 + "_ref"
  refdata@long_name = var0 + " (reference data)"
  nc_outfile = ncdf_write(refdata, nc_filename)

  diff@var = var0 + "_bias"
  diff@long_name = var0 + " (multi-model bias)"
  nc_outfile = ncdf_write(diff, nc_filename)

  if (isvar("absdiff")) then
    absdiff@var = var0 + "_abs_bias"
    absdiff@long_name = var0 + " (multi-model absolute bias)"
    nc_outfile = ncdf_write(absdiff, nc_filename)
  end if

  if (isvar("rmsdiff")) then
    rmsdiff@var = var0 + "_rms_bias"
    rmsdiff@long_name = var0 + " (multi-model root mean square difference)"
    nc_outfile = ncdf_write(rmsdiff, nc_filename)
  end if

  if (isvar("reldiff")) then
    reldiff@var = var0 + "_rel_bias"
    reldiff@long_name = var0 + " (multi-model relative bias)"
    reldiff@units = reldiff@units
    nc_outfile = ncdf_write(reldiff, nc_filename)
  end if

  if (stippandhatch) then
    to_plot_confl_sign@diag_script = (/DIAG_SCRIPT/)
    to_plot_confl_sign@var = "confl_sign"
    to_plot_confl_sign@var_long_name = "more than 66% of models show change " \
                                       + "greater than variability threshold" \
                                       + " and <80% agree on sign of change"
    to_plot_confl_sign@units     = 1
    ncdf_outfile = ncdf_write(to_plot_confl_sign, nc_filename)

    to_plot_not_robust@diag_script = (/DIAG_SCRIPT/)
    to_plot_not_robust@var = "not_robust"
    to_plot_not_robust@var_long_name = "more than 66% of models show change " \
                                       + "greater than variability threshold"
    to_plot_not_robust@units     = 1
    ncdf_outfile = ncdf_write(to_plot_not_robust, nc_filename)
  end if

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s) (mean)
  ; ------------------------------------------------------------------------

  authors = (/"bock_lisa"/)
  statistics = (/"clim", "diff"/)
  domain = ("global")
  plottype = ("geo")
  references = (/"eyring21ipcc"/)
  prov_caption = caption + " for variable " + var0 \
    + ", reference = " + names(ref_ind) + "."
  infiles := array_append_record(infiles, natvar_file, 0)

  log_provenance(nc_filename, plotfile, prov_caption, statistics, \
                 domain, plottype, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")

end
