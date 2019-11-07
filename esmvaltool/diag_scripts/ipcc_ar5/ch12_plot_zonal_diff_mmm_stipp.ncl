; #############################################################################
; Multi-model mean change zonal plot (with significance and non-significance)
; Author: Ruth Lorenz (ETH, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Description
;    Annual mean or seasonal change in multi-model mean with signficance
;    and non-significance
;    lev-lon field precalculated, plotting only
;
; Caveats
;    Only handles one variable at a time
;
; Modification history
;    20181105-A_lorenz_ruth: ported to v2
;    20171002-A_lorenz_ruth: Added tags for searching
;    20170523-A_lorenz_ruth: written for ESMValTool
;                        based on IPCC AR5 scripts.
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/zonalmean_profile.ncl"
load "$diag_scripts/shared/plot/contourplot.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(DIAG_SCRIPT, "")

  ; 'datasets', 'variables' are fetched from the
  ; above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  ; Check required diag_script_info attributes
  req_atts = (/"scenarios", "periods"/)
  exit_if_missing_atts(diag_script_info, req_atts)
end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if
  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  netcdf_dir = diag_script_info@input_files(1)

  dim_scen = dimsizes(diag_script_info@scenarios)
  dim_per = dimsizes(diag_script_info@periods)
  if (isatt(diag_script_info, "seasons")) then
    dim_seas = dimsizes(diag_script_info@seasons)
  else
    dim_seas = 1
  end if
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data                                #
  ; ###########################################
  datapath = netcdf_dir + "/" + variable_info[0]@diagnostic + ".nc"
  to_plot = ncdf_read(datapath, var0)

  ; Check dimensions
  dims = getvardims(to_plot)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
  end if
  idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
  end if
  lev_name = dims(idx)

  if (isatt(diag_script_info, "sig") .and. \
      diag_script_info@sig .eq. True) then
    ; read already calculated significance
    to_plot_signif = ncdf_read(datapath, "signif")
  end if
  if (isatt(diag_script_info, "not_sig") .and. \
      diag_script_info@not_sig .eq. True) then
    ; read already calculated significance
    to_plot_not_signif = ncdf_read(datapath, "not_signif")
  end if
  if (isatt(diag_script_info, "model_nr") .and. \
      diag_script_info@model_nr .eq. True) then
    ; read already calculated number of models included in plot
    model_number = ncdf_read(datapath, "model_nr")
  end if
  if (isatt(diag_script_info, "base_cn")) then
    to_plot_base = ncdf_read(datapath, var0 + "_base")
  end if
  ; ######################################################
  ; # Separate output to a graphics file, Panelling      #
  ; ######################################################
  ; function in aux_plotting.ncl
  wks = get_wks(0, variable_info[0]@diagnostic, \
                var0 + "_" + diag_script_info@time_avg)

  if (isatt(diag_script_info, "max_vert")) then
    ; Maximum allowed number of plots per page (vertical)
    nvert = toint(diag_script_info@max_vert)
  else
    nvert = dim_scen
  end if
  if (isatt(diag_script_info, "max_hori")) then
    ; Maximum allowed number of plots per page (horizontal)
    nhori = toint(diag_script_info@max_hori)
  else
    nhori = 1
  end if

  to_plot@res_gsnDraw           = False  ; Do not draw yet
  to_plot@res_gsnFrame          = False  ; Don't advance frame.
  to_plot@res_cnFillOn          = True            ; Color plot desired
  to_plot@res_cnLineLabelsOn    = False           ; Contour lines
  to_plot@res_cnLinesOn         = False
  to_plot@res_cnInfoLabelOn     = False
  to_plot@res_tiMainOn = False
  if(isatt(diag_script_info, "colormap")) then
    colormap = read_colormap_file(diag_script_info@colormap)
  else
    colormap = "BlueYellowRed"
  end if
  to_plot@res_cnFillPalette = colormap
  if (isatt(diag_script_info, "span")) then
    to_plot@res_cnSpanFillPalette     = diag_script_info@span
  else
    to_plot@res_cnSpanFillPalette     = True   ; use full colormap
  end if
  to_plot@res_lbLabelBarOn = False
  if (lev_name .eq. "lev") then
    to_plot@res_trYReverse            = True
    to_plot@res_gsnYAxisIrregular2Linear = True
  else
    to_plot@res_trYMinF = 50.
    to_plot@res_trYMaxF = 1000.
  end if
  if (isatt(diag_script_info, "base_cn")) then
    to_plot_base@res_gsnFrame = False
    to_plot_base@res_gsnDraw = False
    to_plot_base@res_gsnStringFontHeightF = 0.019
    to_plot_base@res_cnInfoLabelOn             = True
    to_plot_base@res_cnFillOn                  = False
    to_plot_base@res_cnLinesOn                 = True
    to_plot_base@res_cnLineLabelsOn            = False
    to_plot_base@res_lbLabelBarOn             = False
    to_plot_base@res_trYReverse               = True
    to_plot_base@res_trYMinF = 50.
    to_plot_base@res_trYMaxF = 1000.
    to_plot_base@res_gsnContourNegLineDashPattern = 1
    to_plot_base@res_gsnContourZeroLineThicknessF = 2.

    to_plot_base@res_gsnCenterString = ""
    to_plot_base@res_gsnLeftString = ""
    to_plot_base@res_gsnLeftStringOrthogonalPosF = -.003
    to_plot_base@res_gsnRightString = ""  ; avg_period
    if (isatt(diag_script_info, "base_cnMaxLevel")) then
      to_plot_base@res_cnMaxLevelValF   = diag_script_info@base_cnMaxLevel
    else
      to_plot_base@res_cnMaxLevelValF   = 40
    end if
    if (isatt(diag_script_info, "base_cnMinLevel")) then
      to_plot_base@res_cnMinLevelValF   = diag_script_info@base_cnMinLevel
    else
      to_plot_base@res_cnMinLevelValF   = - (to_plot_base@res_cnMaxLevelValF)
    end if
    if (isatt(diag_script_info, "base_cnLevelSpacing")) then
      to_plot_base@res_cnLevelSpacingF = diag_script_info@base_cnLevelSpacing
    else
      to_plot_base@res_cnLevelSpacingF = 10
    end if
  end if
  if (isatt(diag_script_info, "sig") .and. \
      diag_script_info@sig .eq. True) then
    ; resources dots
    to_plot_signif@res_gsnDraw = False
    to_plot_signif@res_gsnFrame = False
    to_plot_signif@res_gsnAddCyclic = True
    to_plot_signif@res_cnLinesOn = False
    to_plot_signif@res_cnLineLabelsOn = False
    to_plot_signif@res_cnLevelSelectionMode = "ExplicitLevels"
    to_plot_signif@res_cnLevels = (/.3, 1.3/)
    to_plot_signif@res_cnInfoLabelOn = False
    to_plot_signif@res_tiMainOn = False
    to_plot_signif@res_cnFillOn = True
    to_plot_signif@res_cnFillMode = "AreaFill"
    to_plot_signif@res_cnMonoFillPattern = False
    to_plot_signif@res_cnFillPatterns = (/-1, 17, 17/)
    to_plot_signif@res_cnMonoFillColor = True
    to_plot_signif@res_lbLabelBarOn = False
    to_plot_signif@res_gsnLeftString = ""
    to_plot_signif@res_gsnRightString = ""
    to_plot_signif@res_cnFillDotSizeF       = 0.001
    to_plot_signif@res_cnFillScaleF        = .6
  end if
  if (isatt(diag_script_info, "not_sig") .and. \
      diag_script_info@not_sig .eq. True) then
    ; resources hatch
    to_plot_not_signif@res_gsnDraw = False
    to_plot_not_signif@res_gsnFrame = False
    to_plot_not_signif@res_gsnAddCyclic = True
    to_plot_not_signif@res_cnLinesOn = False
    to_plot_not_signif@res_cnLineLabelsOn = False
    to_plot_not_signif@res_cnLevelSelectionMode = "ExplicitLevels"
    to_plot_not_signif@res_cnLevels = (/.3, 1.3/)
    to_plot_not_signif@res_cnInfoLabelOn = False
    to_plot_not_signif@res_tiMainOn = False
    to_plot_not_signif@res_cnFillOn = True
    to_plot_not_signif@res_cnFillMode = "AreaFill"
    to_plot_not_signif@res_cnRasterSmoothingOn = False
    to_plot_not_signif@res_cnMonoFillPattern = False
    to_plot_not_signif@res_cnFillPatterns = (/-1, 10, 10/)
    to_plot_not_signif@res_cnMonoFillColor = True
    to_plot_not_signif@res_lbLabelBarOn = False
    to_plot_not_signif@res_gsnLeftString = ""
    to_plot_not_signif@res_gsnRightString = ""
  end if

  nmap = dim_scen * dim_per * dim_seas

  maps = new(nmap, graphic)  ; collect individual "map"s in a graphic array
  do i = 0, nmap - 1
    to_plot@res_gsnLeftString = ""
    to_plot@res_gsnRightString = ""
    if(isatt(diag_script_info, "label")) then
      to_plot@res_gsnCenterString = diag_script_info@label(i)
    end if
    ; Set explicit contour levels
    if(isatt(diag_script_info, "diff_levs")) then
      to_plot@res_cnLevelSelectionMode = "ExplicitLevels"
      to_plot@res_cnLevels = diag_script_info@diff_levs
    end if
    if (lev_name .eq. "plev") then
      maps(i) = zonalmean_profile(wks, to_plot(i, :, :), var0)
    else
      to_plot@res_tiYAxisString = "Depth [m]"
      maps(i) = contourplot(wks, to_plot(i, :, :), var0)
    end if
    if (isatt(diag_script_info, "sig") .and. \
        diag_script_info@sig .eq. True) then
      ; plot dots
      if all(to_plot_signif(i, :, :) .eq. 1) then
        to_plot_signif(i, 0, 0) = 0
      end if
      if (max(to_plot_signif(i, :, :)) .gt. \
          min(to_plot_signif(i, :, :))) then
        to_plot_signif@res_gsnAddCyclic = False
        if (lev_name .eq. "plev") then
          plot1 = zonalmean_profile(wks, to_plot_signif(i, :, :), var0)
        else
          plot1 = contourplot(wks, to_plot_signif(i, :, :), var0)
        end if
        overlay(maps(i), plot1)
      end if
    end if
    if (isatt(diag_script_info, "not_sig") .and. \
        diag_script_info@not_sig .eq. True) then
      ; plot hatch
      if all(to_plot_not_signif(i, :, :) .eq. 1) then
        to_plot_not_signif(i, 0, 0) = 0
      end if
      if (max(to_plot_not_signif(i, :, :)) .gt. \
          min(to_plot_not_signif(i, :, :))) then
        to_plot_not_signif@res_gsnAddCyclic = False
        to_plot_not_signif@res_cnFillScaleF = 0.45
        if (lev_name .eq. "plev") then
          plot2 = zonalmean_profile(wks, to_plot_not_signif(i, :, :), var0)
        else
          plot2 = contourplot(wks, to_plot_not_signif(i, :, :), var0)
        end if
        overlay(maps(i), plot2)
      end if
    end if

    if (isatt(diag_script_info, "base_cn")) then
      if (lev_name .eq. "plev") then
        cont1 = zonalmean_profile(wks, to_plot_base(i, :, :), var0)
      else
        cont1 = contourplot(wks, to_plot_base(i, :, :), var0)
      end if
      overlay(maps(i), cont1)
    end if
  end do

  txres                       = True
  txres@txFontHeightF         = 0.03
  amres = True
  amres@amParallelPosF   = 0.48    ; This is the right edge of the plot.
  amres@amOrthogonalPosF = -0.48    ; This is the bottom edge of the plot.
  amres@amJust           = "TopRight"
  if ((isatt(diag_script_info, "model_nr")) .and. \
      (diag_script_info@model_nr .eq. True)) then
    do i = 0, dimsizes(maps) - 1
      txid1 = gsn_create_text(wks, "" + model_number(i), txres)
      annoid1 = gsn_add_annotation(maps(i), txid1, amres)
    end do
  end if
  psles                  = True   ; needed to override panelling defaults
  psles@gsnPanelLabelBar = True  ; general label bar desired here

  a4_height = 29.7  ; in centimeters, if my
  a4_width = 23.0   ; reference is correct
  cm_per_inch = 2.54

  psles@gsnPaperWidth = a4_width / cm_per_inch
  psles@gsnPaperHeight = a4_height / cm_per_inch
  psles@gsnPaperOrientation = "portrait"
  psles@gsnPanelTop = 0.96

  if (isatt(diag_script_info, "units")) then
    psles@lbTitleString = diag_script_info@units
  else
    unit_string = format_units(to_plot@units)
    psles@lbTitleString = "(" + unit_string + ")"
  end if
  psles@lbTitleFontHeightF = 0.017
  psles@lbTitleDirection = "across"
  psles@lbTitlePosition = "Right"
  psles@lbTitleJust = "CenterLeft"
  psles@lbLabelFontHeightF = 0.014
  psles@lbLabelJust = "CenterCenter"
  psles@lbLabelAutoStride = False

  psles@txFontHeightF = 0.02
  if(isatt(diag_script_info, "title")) then
    psles@txString  = diag_script_info@title
  end if
  psles@pmLabelBarParallelPosF = 0.06
  psles@pmLabelBarOrthogonalPosF = -0.01

  outfile = panelling(wks, maps, nvert, nhori, psles)

  log_info(" Wrote " + outfile)

  ; collect meta-data
  nc_file = datapath
  plot_file = outfile
  caption = "Multi-model changes in zonal mean " + \
            variable_info[0]@long_name + \
            " relative to " + info_items[0]@start_year + "-" +  \
            info_items[0]@end_year + " for 20 year periods starting in " + \
            str_join(diag_script_info@periods, ", ") + " under the " + \
            str_join(diag_script_info@scenarios, ", ") + \
            " forcing scenarios. Hatching indicates regions where the " + \
            "multi-model mean change is less than one standard " + \
            "deviation of internal variability. Stippling indicates " + \
            "regions  where the multi-model change mean is greater " + \
            "than two standard deviations of internal variability and " + \
            "where at least 90% of the models agree on the sign of " + \
            "change. Similar to Figure 12.12" + \
            " of Collins et al. (2013)."
  statistics = ("mean")
  domains = ("global")
  plot_types = ("zonal")
  authors = (/"lorenz_ruth"/)
  references = (/"collins13ipcc"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, plot_file, caption, statistics, domains, \
                 plot_types, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")
end
