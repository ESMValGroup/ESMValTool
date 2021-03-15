; ch09_fig09-3
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME CRESCENDO
; ############################################################################
; Description
;    Calculates the multi-model mean seasonality (DJF - JJA), multi-model mean
;    of absolute seasonality, multi-model mean bias in seasonality, and multi-
;    model mean bias in absolute seasonality (similar to IPCC AR5 ch. 9
;    fig. 9.3).
;
; Required diag_script_info attributes (diagnostics specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   projection: map projection, e.g., Mollweide, Mercator (default = Robinson)

; Required variable_info attributes (variable specific)
;   reference_dataset: name of reference datatset
;
; Optional variable_info attributes (variable specific)
;   map_diff_levels: explicit contour levels for plotting
;
; Caveats
;   none
;
; Modification history
;   20200714-lauer_axel: code rewritten for ESMValTool v2.0, renamed to
;                        ch09_fig09-3.ncl (was seasonality_mm.ncl)
;   20170622-lauer_axel: added tags for reporting
;   20170320-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
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

  ; check for reference dataset definition
  if (.not.(isvar("refname"))) then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
  end if

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "projection", "Robinson")

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  if (diag_script_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if
end

begin
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

  climofiles = new(2, string)
  climofiles(0) = infiles(mm_ind)
  climofiles(1) = infiles(ref_ind)

  ; ========================================================================
  ; =========================== calculations ===============================
  ; ========================================================================

  ; note: 1) masking is handled by the backend
  ;       2) multi-model mean is calculated by the backend

  ; read data and calculate seasonal means

  A0 = read_data(info0[mm_ind])
  mmdata_djf = time_operations(A0, -1, -1, "average", "DJF", True)
  mmdata_jja = time_operations(A0, -1, -1, "average", "JJA", True)
  delete(A0)

  A0 = read_data(info0[ref_ind])
  refdata_djf = time_operations(A0, -1, -1, "average", "DJF", True)
  refdata_jja = time_operations(A0, -1, -1, "average", "JJA", True)
  delete(A0)

  ; seasonality = DJF - JJA

  seas_ref = refdata_djf - refdata_jja
  seas_mm = mmdata_djf - mmdata_jja
  seas_mm@diag_script = (/DIAG_SCRIPT/)
  copy_VarMeta(refdata_djf, seas_ref)
  copy_VarMeta(seas_ref, seas_mm)
  delete(refdata_djf)
  delete(refdata_jja)
  delete(mmdata_djf)
  delete(mmdata_jja)

  abs_seas = abs(seas_mm)
  bias_seas = seas_mm - seas_ref
  bias_abs_seas = abs_seas - abs(seas_ref)

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  if (write_plots) then
    seas_mm@res_gsnMaximize       = True         ; use full page for the plot
    seas_mm@res_cnFillOn          = True         ; color plot desired
    seas_mm@res_cnLineLabelsOn    = False        ; contour lines
    seas_mm@res_cnLinesOn         = False
    seas_mm@res_tiMainOn                 = True
    seas_mm@res_gsnLeftStringFontHeightF = 0.015
    seas_mm@res_cnLevelSelectionMode     = "ExplicitLevels"
    seas_mm@res_mpOutlineOn              = True
    seas_mm@res_mpFillOn                 = False
    seas_mm@res_lbLabelBarOn             = True
    seas_mm@res_gsnRightString           = ""
    seas_mm@res_mpFillDrawOrder          = "PostDraw"     ; draw map fill last
    seas_mm@res_cnMissingValFillColor    = "Gray"
    seas_mm@res_tmYLLabelsOn             = False
    seas_mm@res_tmYLOn                   = False
    seas_mm@res_tmYRLabelsOn             = False
    seas_mm@res_tmYROn                   = False
    seas_mm@res_tmXBLabelsOn             = False
    seas_mm@res_tmXBOn                   = False
    seas_mm@res_tmXTLabelsOn             = False
    seas_mm@res_tmXTOn                   = False
    seas_mm@res_cnInfoLabelOn            = False      ; turn off cn info label
    seas_mm@res_mpProjection             = diag_script_info@projection

    ; --------------------------------------------------------------------
    ; create workspace

    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "seasonality_" + var0)
;    drawNDCGrid(wks) ; debugging option

    plots = new((/4/), graphic)

    ; --------------------------------------------------------------------
    ; plot contour map

    seas_mm@res_gsnDraw         = False  ; Do not draw yet
    seas_mm@res_gsnFrame        = False  ; Don't advance frame.
    seas_mm@res_mpPerimOn       = False

    seas_mm@res_lbTitleString   = "(" + seas_mm@units + ")"
    seas_mm@res_lbTitlePosition      = "Bottom"

    seas_mm@res_lbLabelFontHeightF   = 0.014
    seas_mm@res_lbTopMarginF         = 0.1
    seas_mm@res_lbTitleFontHeightF   = 0.014

    seas_mm@res_tiMainFontHeightF    = 0.016

    seas_mm@res_tiMainString = "Multi Model Mean Bias"

    if (isatt(seas_mm, "res_cnLevels")) then
      delete(seas_mm@res_cnLevels)
    end if

    if (isatt(variable_info[0], "map_diff_levels")) then
      seas_mm@res_cnLevels = variable_info@map_diff_levels
    end if

    copy_VarMeta(seas_mm, abs_seas)
    copy_VarMeta(seas_mm, bias_seas)
    copy_VarMeta(seas_mm, bias_abs_seas)

    if (var0.eq."tas") then
      pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                               + "ipcc-tas.rgb")
      seas_mm@res_cnFillColors = pal
      pal2 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-tas-absdelta.rgb")
      abs_seas@res_cnFillColors = pal2
      pal3 = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-tas-seasdelta.rgb")
      bias_seas@res_cnFillColors = pal3
      bias_abs_seas@res_cnFillColors = pal3

      if (isatt(seas_mm, "res_cnLevels")) then
        delete(seas_mm@res_cnLevels)
      end if
      if (isatt(abs_seas, "res_cnLevels")) then
        delete(abs_seas@res_cnLevels)
      end if
      if (isatt(bias_seas, "res_cnLevels")) then
        delete(bias_seas@res_cnLevels)
      end if
      if (isatt(bias_abs_seas, "res_cnLevels")) then
        delete(bias_abs_seas@res_cnLevels)
      end if

      ; IPCC ch.9 fig. 9.3 scales for tas
      seas_mm@res_cnLevels = fspan(-50, 50, 21)
      abs_seas@res_cnLevels = fspan(5, 50, 10)
      bias_seas@res_cnLevels = fspan(-8, 8, 17)
      bias_abs_seas@res_cnLevels = fspan(-8, 8, 17)
    end if

    seas_mm@res_tiMainString = "Multi Model Mean"
    abs_seas@res_tiMainString = "Multi Model Mean of Absolute Seasonality"
    bias_seas@res_tiMainString = "Multi Model Mean Bias in Seasonality"
    bias_abs_seas@res_tiMainString = "Multi Model Mean Bias in Absolute " \
                                     + "Seasonality"

    ; panelling resources

    pres                 = True
    pres@gsnPanelCenter  = False
    pres@gsnPanelRowSpec = True
    pres@gsnPanelYWhiteSpacePercent = 5
    pres@gsnPanelXWhiteSpacePercent = 5

    plots(0) = contour_map(wks, seas_mm, var0)
    plots(1) = contour_map(wks, abs_seas, var0)
    plots(2) = contour_map(wks, bias_seas, var0)
    plots(3) = contour_map(wks, bias_abs_seas, var0)

    gsn_panel(wks, plots, (/2, 2/), pres)

    plotfile = wks@fullname
  else
    plotfile = ""
  end if

  ; ###########################################
  ; # output to netCDF                        #
  ; ###########################################

  seas_mm@var = "seas_" + var0
  seas_mm@var_long_name = "Multi-model mean seasonality"
  seas_mm@var_units = seas_mm@units

  abs_seas@var = "abs_seas_" + var0
  abs_seas@var_long_name = "Multi-model mean of absolute seasonality"
  abs_seas@var_units = abs_seas@units

  bias_seas@var = "bias_seas_" + var0
  bias_seas@var_long_name = "Multi-model mean bias in seasonality"
  bias_seas@var_units = bias_seas@units

  bias_abs_seas@var = "bias_abs_seas_" + var0
  bias_abs_seas@var_long_name = "Multi-model mean bias in absolute seasonality"
  bias_abs_seas@var_units = bias_abs_seas@units

  nc_filename = work_dir + "seasonality_" + var0 + ".nc"
  nc_filename@existing = "append"

  nc_outfile = ncdf_write(seas_mm, nc_filename)
  nc_outfile = ncdf_write(abs_seas, nc_filename)
  nc_outfile = ncdf_write(bias_seas, nc_filename)
  nc_outfile = ncdf_write(bias_abs_seas, nc_filename)

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s) (mean)
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "diff", "mean"/)
  domain = "global"
  plottype = "geo"

  caption = "Multi model values for variable " + var0 + ", from top left " \
    + "to bottom right: mean, mean of absolute seasonality, mean bias in " \
    + "seasonality, mean bias in absolute seasonality. Reference dataset: " \
    + names(ref_ind) + ". Similar to IPCC AR5, fig. 9.3."

  log_provenance(nc_outfile, plotfile, caption, statistics, \
                 domain, plottype, "", "", climofiles)

  leave_msg(DIAG_SCRIPT, "")

end
