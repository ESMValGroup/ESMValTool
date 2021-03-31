; MODEL_BIAS
; ############################################################################
; Author: Lisa Bock (DLR, Germany), Axel Lauer (DLR, Germany)
; Crescendo, IPCC-AR6
; ############################################################################
; Description
;   Calculates the multi-model mean bias of annual mean
;   variables compared with a reference dataset (observations).
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
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
;   20190312-A_bock_ls: adapted code (based on clouds/clouds_bias.ncl)
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
  project = input_file_info[1]@project
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

  ; read data

  A0 = read_data(info0[mm_ind])
  mmdata = time_operations(A0, -1, -1, "average", timemean, True)
  delete(A0)

  A0 = read_data(info0[ref_ind])
  refdata = time_operations(A0, -1, -1, "average", timemean, True)
  delete(A0)

  ; convert units for plotting
  if (var0.eq."pr") then
    ; kg m-2 s-1 --> mm day-1
    mmdata = convert_units(mmdata, "mm/day")
    refdata = convert_units(refdata, "mm/day")
  elseif (var0.eq."tas") then
    mmdata = convert_units(mmdata, "degC")
    refdata = convert_units(refdata, "degC")
  elseif (var0.eq."tos" .or. var0.eq."thetao") then
    if (mmdata@units .eq. "K") then
      mmdata = convert_units(mmdata, "degC")
    end if
    if (refdata@units .eq. "K") then
      refdata = convert_units(refdata, "degC")
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

  plots = new((/2, numseas/), graphic)

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

  diff@res_tiMainString = project + " Multi-Model Mean Bias"

  copy_VarMeta(diff, mmdata)
  delete(mmdata@res_cnLevels)

  mmdata@res_tiMainString = project + " Multi-Model Mean"

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
  elseif (var0.eq."pr") then
    mmdata@res_cnLevels = ispan(1, 10, 1)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_precipitation_seq.rgb")
    mmdata@res_cnFillPalette = pal
    pal2 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                              + "ipcc-ar6_precipitation_div.rgb")
    diff@res_cnFillPalette = pal2
  elseif (var0.eq."tos" .or. var0.eq."thetao") then
    mmdata@res_tiMainString = "Sea Surface Temperature ~C~Multi Model Mean"
    mmdata@res_cnLevels = ispan(-30, 30, 5)

    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-ar6_temperature_div.rgb")
    mmdata@res_cnFillPalette = pal
    diff@res_tiMainString = "Sea Surface Temperature ~C~Multi Model Mean Bias"
    diff@res_cnFillPalette = pal
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

  ; root mean square differences

  ; relative differences

  ; add global bias and rmsd value
  diff@res_gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd)
  diff@res_gsnLeftString  = "bias = " + sprintf("%6.3f", bias)

  ; panelling resources
  pres  = True
  pres@gsnPanelCenter  = False
  pres@gsnPanelRowSpec = True  ; tell panel what order to plot
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  pres@gsnPanelFigureStrings = (/"a)", "b)"/)
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
                    + "_" + season(is) + "_" + project)

      plotfile(is) = wks@fullname

      if (numseas.gt.1) then
        pres@txString = season(is)
        plots(0, is) = contour_map(wks, mmdata(is, :, :), var0)
        plots(1, is) = contour_map(wks, diff(is, :, :), var0)
        gsn_panel(wks, plots(plotind, is), plotsperline, pres)
      else
        plots(0, 0) = contour_map(wks, mmdata, var0)
        plots(1, 0) = contour_map(wks, diff, var0)
        gsn_panel(wks, plots(plotind, 0), plotsperline, pres)
      end if

      ; add meta data to plot (for reporting)

      caption = "Multi-model values, from top left to bottom right: " \
                + "mean, bias"

    end do  ; is-loop (seasons)
  end if  ; if write_plots

  ; ###########################################
  ; # output to netCDF                        #
  ; ###########################################

  nc_filename = work_dir + "model_bias_" + var0 + "_" + project + ".nc"
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
