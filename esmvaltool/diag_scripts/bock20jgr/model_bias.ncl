; MODEL_BIAS
; ############################################################################
; Author: Axel Lauer (DLR, Germany), Lisa Bock (DLR, Germany)
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
;   none
;
; Modification history
;   20190312-A_bock_ls: added calculation of RMSD and IPCC AR6 color scheme
;   20190222-A_laue_ax: added output of provenance (v2.0)
;   20181119-A_laue_ax: adapted code to multi-variable capable framework
;   20180923-A_laue_ax: added writing of results to netcdf
;   20180914-A_laue_ax: code rewritten for ESMValTool v2.0
;   20170620-A_laue_ax: added tags for reporting
;   20160901-A_laue_ax: added regridding option 1 deg x 1 deg
;   20151027-A_laue_ax: moved call to 'write_references' to the beginning
;                       of the code
;   20150428-A-laue_ax: written.
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

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "plot_abs_diff", False)
  set_default_att(diag_script_info, "plot_rms_diff", False)
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

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

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

  if (.not.(isvar("refname"))) then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
  end if

  plot_abs_diff = diag_script_info@plot_abs_diff
  plot_rel_diff = diag_script_info@plot_rel_diff
  plot_rms_diff = diag_script_info@plot_rms_diff

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
  ; =========================== calculations ===============================
  ; ========================================================================

  ; note: 1) masking is handled by the backend
  ;       2) multi-model mean is calculated by the backend

  ; read data

  A0 = read_data(info0[mm_ind])
  mmdata = time_operations(A0, -1, -1, "average", timemean, True)
  if (dimsizes(A0&time).gt.12) then
    mm_mon = time_operations(A0, -1, -1, "average", "monthlyclim", True)
  else
    mm_mon = A0
  end if
  delete(A0)

  A0 = read_data(info0[ref_ind])
  refdata = time_operations(A0, -1, -1, "average", timemean, True)
  if (dimsizes(A0&time).gt.12) then
    ref_mon = time_operations(A0, -1, -1, "average", "monthlyclim", True)
  else
    ref_mon = A0
  end if
  delete(A0)

  dim = dimsizes(mmdata)
  data = new((/dim_MOD-2, dim(0), dim(1)/), float)
  data_mon = new((/dim_MOD - 2, 12, dim(0), dim(1)/), float)

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

      data(ii, :, :) = time_operations(A0, -1, -1, "average", timemean, True)
      data_mon(ii, :, :, :) = time_operations(A0, -1, -1, "average", \
                                              "monthlyclim", True)

      ii = ii + 1

      delete(A0)

    end if

  end do

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
  elseif (var0.eq."so") then
    refdata = refdata * 1000.
  end if

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
  elseif (var0 .eq. "so") then
    diff@units = "0.001"
    diff@res_cnLevels = ispan(-30, 30, 5) * 0.1
  end if

  diff_i    = new((/dim_MOD-2, dim(0), dim(1)/), float)
  absdiff_i = new((/dim_MOD-2, dim(0), dim(1)/), float)
  rmsdiff_i = new((/dim_MOD-2, dim(0), dim(1)/), float)
  reldiff_i = new((/dim_MOD-2, dim(0), dim(1)/), float)

  if (plot_rms_diff .or. plot_abs_diff .or. plot_rel_diff) then

    ii = 0

    do imod = 0, dim_MOD - 1

      if(imod .ne. mm_ind .and. imod .ne. ref_ind) then
        log_info("processing " + names(imod))

        ; difference
        diff_i(ii, :, :) = data(ii, :, :) - refdata

        ; absolute difference
        absdiff_i(ii, :, :) = abs(diff_i(ii, :, :))

        ; root mean square error
        rmsdiff_i(ii, :, :) = calculate_metric(data_mon(ii, :, :, :), ref_mon,\
                                               "RMSDxy")

        ; relative difference
        ; replace "epsilon" values with missing value
        datam = where(abs(data).lt.1.e-6, data@_FillValue, data)
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

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  climofiles = new(2, string)
  climofiles(0) = infiles(mm_ind)
  climofiles(1) = infiles(ref_ind)

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

  diff@var = var0  ; Overwrite existing entry
  if (isatt(variable_info[0], "long_name")) then
    diff@var_long_name = variable_info[0]@long_name
  end if
  diff@var_units  = diff@units

  plots = new((/5, numseas/), graphic)

  ; --------------------------------------------------------------------
  ; plot contour map

  diff@res_gsnDraw  = False  ; Do not draw yet
  diff@res_gsnFrame  = False  ; Don't advance frame.
  diff@res_mpPerimOn  = False

  diff@res_lbTitleString  = "(" + diff@units + ")"
  diff@res_lbTitlePosition  = "Bottom"

  diff@res_lbLabelFontHeightF  = 0.014
  diff@res_lbTopMarginF  = 0.1
  diff@res_lbTitleFontHeightF  = 0.014

  diff@res_tiMainFontHeightF  = 0.016

  diff@res_tiMainString = "Multi-Model Mean Bias"

  copy_VarMeta(diff, mmdata)
  delete(mmdata@res_cnLevels)

  mmdata@res_tiMainString = "Multi-Model Mean"

  mmdata@res_gsnLeftStringFontHeightF = 0.015
  mmdata@res_gsnRightStringFontHeightF = 0.015
  mmdata@res_gsnRightString  = " "
  mmdata@res_gsnLeftString   = " "

  if (var0.eq."tas") then
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
                                + "ipcc-ar6_temperature_div.rgb")
    end if
  elseif (var0.eq."pr") then
    mmdata@res_cnLevels = ispan(1, 10, 1)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_precipitation_seq.rgb")
    mmdata@res_cnFillPalette = pal
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
                                + "ipcc-ar6_precipitation_div.rgb")
    end if
  elseif (var0.eq."tos" .or. var0.eq."thetao") then
    mmdata@res_tiMainString = "Sea Surface Temperature ~C~Multi Model Mean"
    mmdata@res_cnLevels = ispan(-30, 30, 5)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_temperature_div.rgb")
    mmdata@res_cnFillPalette = pal
    diff@res_tiMainString = "Sea Surface Temperature ~C~Multi Model Mean Bias"
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
                                + "ipcc-ar6_temperature_10.rgb")
    end if
  elseif (var0.eq."so") then
    mmdata@res_tiMainString = "Sea Surface Salinity ~C~Multi Model Mean"
    mmdata@res_cnLevels = ispan(310, 370, 5) * 0.1

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_misc_seq_1.rgb")
    mmdata@res_cnFillPalette = pal
    diff@res_tiMainString = "Sea Surface Salinity ~C~Multi Model Mean Bias"
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
;    if (isvar("pal3")) then
;      delete(absdiff@res_cnFillColors)
;      absdiff@res_cnFillColors = pal3(4:0,:)
;      nt = dimsizes(pal3(4:0,:))
;      n = nt(0)
;    else
;      n = dimsizes(diff@res_cnLevels)
;    end if
;
;    mx = max(diff@res_cnLevels)
;    delete(absdiff@res_cnLevels)
;
;    tmp = fspan(0.0, mx, n)
;
;    absdiff@res_cnLevels = tmp(1:dimsizes(tmp)-1)
;    delete(tmp)

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
;    if (isvar("pal5")) then
;      delete(rmsdiff@res_cnFillColors)
;      if (var0 .eq. "tas") then
;        rmsdiff@res_cnFillColors = pal5(4:0,:)
;        nt = dimsizes(pal5(4:0,:))
;      elseif (var0 .eq. "pr") then
;        rmsdiff@res_cnFillColors = pal5(5:10,:)
;        nt = dimsizes(pal5(5:10,:))
;      elseif (var0 .eq. "tos") then
;        rmsdiff@res_cnFillColors = pal5(4:0,:)
;        nt = dimsizes(pal5(4:0,:))
;      end if
;      n = nt(0)
;    else
;      n = dimsizes(diff@res_cnLevels)
;    end if
;    mx = max(diff@res_cnLevels)
;    delete(rmsdiff@res_cnLevels)
;
;    tmp = fspan(0.0, mx, n)
;
;    rmsdiff@res_cnLevels = tmp(1:dimsizes(tmp)-1)
;    delete(tmp)
;
;    if (var0 .eq. "pr") then
;      delete(rmsdiff@res_cnLevels)
;      rmsdiff@res_cnLevels = (/0.5,1.,1.5,2.,2.5/)
;    end if

    rmsdiff@res_tiMainString = "Multi Model Mean of Root Mean Square Error"

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

  ; add global bias and rmsd value
  diff@res_gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd)
  diff@res_gsnLeftString  = "bias = " + sprintf("%6.3f", bias)

  ; panelling resources
  pres  = True
  pres@gsnPanelCenter  = False
  pres@gsnPanelRowSpec = True  ; tell panel what order to plot
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  pres@gsnPanelFigureStrings = (/"a)", "b)", "c)", "d)"/)
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsJust = "TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.016

  plotfile = new(numseas, string)
  plotfile(:) = ""

  if (write_plots) then
    do is = 0, numseas - 1
      ; --------------------------------------------------------------------
      ; create workspace

      if (isvar("wks")) then
        delete(wks)
      end if

      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "model_bias_" + var0 \
                    + "_" + season(is))

      plotfile(is) = wks@fullname

      if (numseas.gt.1) then
        pres@txString = season(is)
        plots(0, is) = contour_map(wks, mmdata(is, :, :), var0)
        plots(1, is) = contour_map(wks, diff(is, :, :), var0)
        if (plot_abs_diff) then
          plots(2, is) = contour_map(wks, absdiff(is, :, :), var0)
        end if
        if (plot_rms_diff) then
          plots(3, is) = contour_map(wks, rmsdiff(is, :, :), var0)
        end if
        if (plot_rel_diff) then
          plots(4, is) = contour_map(wks, reldiff(is, :, :), var0)
        end if
        gsn_panel(wks, plots(plotind, is), plotsperline, pres)
      else
        plots(0, 0) = contour_map(wks, mmdata, var0)
        plots(1, 0) = contour_map(wks, diff, var0)
        if (plot_abs_diff) then
          plots(2, 0) = contour_map(wks, absdiff, var0)
        end if
        if (plot_rms_diff) then
          plots(3, 0) = contour_map(wks, rmsdiff, var0)
        end if
        if (plot_rel_diff) then
          plots(4, 0) = contour_map(wks, reldiff, var0)
        end if
        gsn_panel(wks, plots(plotind, 0), plotsperline, pres)
      end if

      ; add meta data to plot (for reporting)

      caption = "Multi-model values, from top left to bottom right: " \
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

    end do  ; is-loop (seasons)
  end if  ; if write_plots

  ; ###########################################
  ; # output to netCDF                        #
  ; ###########################################

  nc_filename = work_dir + "model_bias_" + var0 + ".nc"
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

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s) (mean)
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "diff"/)
  domain = ("global")
  plottype = ("geo")
  prov_caption = caption + " for variable " + var0 \
    + " (" + allseas + "), reference = " + names(ref_ind) + "."

  do is = 0, numseas - 1
    log_provenance(nc_outfile, plotfile(is), prov_caption, statistics, \
                   domain, plottype, "", "", climofiles)
  end do

  leave_msg(DIAG_SCRIPT, "")

end
