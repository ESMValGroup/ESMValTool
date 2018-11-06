; #############################################################################
; Multi-model mean change map with significance
; Author: Ruth Lorenz (ETH, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Caveats
;    Loop to plot is variables, seasons, scenarios, periods
;    Can either plot seasonal means or annual means in one panel, but not both
;
; Modification history
;    20181307-A_lore_ru: Ported to version 2
;    20171002-A_lore_ru: Added tags for searching
;    20170530-A_lore_ru: multiple variables possible
;    20170120-A_lore_ru: adapted towritten for ESMValTool
;                        based on IPCC AR5 scripts.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/contour_maps.ncl"
load "./diag_scripts/shared/plot/aux_plotting.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ; 'diag_script' and 'variables' are fetched from the
  ; above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; 'models', 'variables' and 'field_types' are fetched from the
  ; above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)

  ; check if diagnostic is change scaled by global temperature
  ; in this case even if two variables present only one to be plotted
  tf = isStrSubset(variable_info@diagnostic, "scaleT")
  if (tf) then
    dim_VAR = 1
  end if
  var0 = variables(0)
  field_type0 = field_types(0)

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
  if(isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data                                #
  ; ###########################################
  if (dim_VAR .eq. 1) then
    datapath = netcdf_dir + "/" + variable_info@diagnostic + ".nc"
  else
    datapath = netcdf_dir + "/" + var0 + "_" + variable_info@diagnostic + ".nc"
  end if
  tmp_plot = ncdf_read(datapath, var0)

  ; Check dimensions
  dims = getvardims(tmp_plot)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  nlat = dimsizes(tmp_plot&$dims(idx)$)
  idx = ind(dims .eq. "lon")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lon dimension")
  end if
  nlon = dimsizes(tmp_plot&$dims(idx)$)

  ; predefine data arrays
  to_plot = new((/dim_VAR * dim_scen * dim_per * dim_seas, nlat, nlon/), \
                float)
  if (isatt(diag_script_info, "sig") .and. \
      diag_script_info@sig .eq. True) then
    to_plot_signif = new((/dim_VAR * dim_scen * dim_per * dim_seas, nlat, \
                          nlon/), float)
  end if
  if (isatt(diag_script_info, "not_sig") .and. \
      diag_script_info@not_sig .eq. True) then
    to_plot_not_signif = new((/dim_VAR * dim_scen * dim_per * dim_seas, \
                              nlat, nlon/), float)
  end if
  if (isatt(diag_script_info, "model_nr") .and. \
      diag_script_info@model_nr .eq. True) then
    model_number = new((/dim_VAR * dim_scen * dim_per * dim_seas/), float)
  end if
  if (tf) then
    pan_var = dim_per * dim_seas
  else
    pan_var = dim_scen * dim_per * dim_seas
  end if

  do v = 0, dim_VAR - 1
    if (v .eq. 0) then
      ind1 = 0
      ind2 = pan_var - 1
    else
      datapath = netcdf_dir + "/" + variables(v) + "_" + \
                 variable_info@diagnostic + ".nc"
      tmp_plot = ncdf_read(datapath, variables(v))
      ind1 = v * pan_var
      ind2 = ind1 + pan_var - 1
    end if
    to_plot(ind1 : ind2, :, :) = tmp_plot

    if (isatt(diag_script_info, "sig") .and. \
        diag_script_info@sig .eq. True) then
      ; read already calculated significance
      to_plot_signif(ind1 : ind2, :, :) = ncdf_read(datapath, "signif")
    end if
    if (isatt(diag_script_info, "not_sig") .and. \
        diag_script_info@not_sig .eq. True) then
      ; read already calculated significance
      to_plot_not_signif(ind1 : ind2, :, :) = ncdf_read(datapath, "not_signif")
    end if
    if (isatt(diag_script_info, "model_nr") .and. \
        diag_script_info@model_nr .eq. True) then
      ; read already calculated number of models included in plot
      model_number(ind1 : ind2) = ncdf_read(datapath, "model_nr")
    end if

  end do  ; loop over variables
  ; ######################################################
  ; # Separate output to a graphics file, Panelling      #
  ; ######################################################
  to_plot@res_gsnDraw           = False  ; Do not draw yet
  to_plot@res_gsnFrame          = False  ; Don't advance frame.
  to_plot@res_cnFillOn          = True            ; Color plot desired
  to_plot@res_cnLineLabelsOn    = False           ; Contour lines
  to_plot@res_cnLinesOn         = False
  if (isatt(diag_script_info, "colormap")) then
    col = read_colormap_file(diag_script_info@colormap)
  else
    col = "BlueYellowRed"
  end if
  to_plot@res_cnFillPalette     = col
  if (isatt(diag_script_info, "span")) then
    to_plot@res_cnSpanFillPalette     = diag_script_info@span
  else
    to_plot@res_cnSpanFillPalette     = True   ; use full colormap
  end if
  to_plot@res_lbLabelBarOn = False

  if (isatt(diag_script_info, "sig") .and. \
      diag_script_info@sig .eq. True) then
    ; resources dots
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
    ress@cnFillPatterns = (/-1, 17, 17/)
    ress@cnMonoFillColor = True
    ress@lbLabelBarOn = False
    ress@gsnLeftString = ""
    ress@gsnRightString = ""
    ress@cnFillDotSizeF       = 0.001
    ress@cnFillScaleF        = .6        ; add extra density
  end if
  if (isatt(diag_script_info, "not_sig") .and. \
      diag_script_info@not_sig .eq. True) then
    ; resources hatch
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
    resb@cnFillPatterns = (/-1, 10, 10/)
    resb@cnMonoFillColor = True
    resb@lbLabelBarOn = False
    resb@gsnLeftString = ""
    resb@gsnRightString = ""
    resb@cnFillScaleF = 0.25
  end if
  if (isatt(diag_script_info, "seasons")) then
    dim_seas = dimsizes(diag_script_info@seasons)
  else
    dim_seas = 1
  end if
  if ((isatt(diag_script_info, "label")) .and. \
      (dimsizes(diag_script_info@label) .ne. \
      (dim_scen * dim_per * dim_seas))) then
    dim_scen = 1
  end if

  nmap = dim_VAR * dim_scen * dim_per * dim_seas

  ; function in aux_plotting.ncl
  if (isatt(diag_script_info, "pltname")) then
    wks = get_wks(0, variable_info@diagnostic, \
                  diag_script_info@pltname + "_" + \
                  diag_script_info@time_avg)
  else
    wks = get_wks(0, variable_info@diagnostic, \
                  var0 + "_" + diag_script_info@time_avg)
  end if
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
    nhori = dim_per
  end if
  if ((tf) .and. (nhori * nvert .lt. nmap)) then
    nhori = nmap
  end if
  maps = new(nmap, graphic)  ; collect individual maps in a graphic array
  do i = 0, nmap - 1  ; this stupid loop creates nmap times the same plot
    if(isatt(diag_script_info, "label")) then
      to_plot@res_tiMainString = diag_script_info@label(i)
    end if
    to_plot@res_gsnRightString = ""
    to_plot@res_gsnLeftString = ""
    if(isatt(diag_script_info, "projection")) then
      to_plot@res_mpProjection      = diag_script_info@projection
    else
      to_plot@res_mpProjection    = "Robinson"
    end if
    ; Set explicit contour levels
    if (isatt(diag_script_info, "diff_levs")) then
      to_plot@res_cnLevelSelectionMode = "ExplicitLevels"
      to_plot@res_cnLevels = diag_script_info@diff_levs
    end if

    to_plot@res_mpPerimOn         = False

    maps(i) = contour_map(wks, to_plot(i, :, :), var0)

    if (isatt(diag_script_info, "sig") .and. \
        diag_script_info@sig .eq. True) then
      ; plot dots
      if all(to_plot_signif(i, :, :) .eq. 1) then
        to_plot_signif(i, 0, 0) = 0
      end if
      if (max(to_plot_signif(i, :, :)) .gt. \
          min(to_plot_signif(i, :, :))) then
        plot1 = gsn_csm_contour(wks, to_plot_signif(i, :, :), ress)
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
        plot2 = gsn_csm_contour(wks, to_plot_not_signif(i, :, :), resb)
        overlay(maps(i), plot2)
      end if
    end if
  end do

  txres                       = True
  txres@txFontHeightF         = 0.02
  amres = True
  amres@amParallelPosF   = 0.48    ; This is the right edge of the plot.
  amres@amOrthogonalPosF = -0.48    ; This is the bottom edge of the plot.
  amres@amJust           = "TopRight"
  if (isatt(diag_script_info, "model_nr") .and. \
      diag_script_info@model_nr .eq. True) then
    do i = 0, dimsizes(maps) - 1
      txid1 = gsn_create_text(wks, "" + model_number(i), txres)
      annoid1 = gsn_add_annotation(maps(i), txid1, amres)
    end do
  end if
  pres                  = True   ; needed to override panelling defaults
  pres@gsnPanelLabelBar = True  ; no general label bar desired here

  a4_height = 29.7  ; in centimeters
  a4_width = 23.0   ; reference is correct
  cm_per_inch = 2.54

  pres@gsnPaperWidth = a4_width / cm_per_inch
  pres@gsnPaperHeight = a4_height / cm_per_inch
  pres@gsnPaperOrientation = "portrait"
  pres@gsnPanelTop = 0.96
  if (isatt(diag_script_info, "units")) then
    pres@lbTitleString = diag_script_info@units
  else
    if (isatt(diag_script_info, "percent")) then
      if (diag_script_info@percent .eq. 1) then
        unit_string = to_plot@var_units
      else
        unit_string = format_units(to_plot@var_units)
      end if
    else
      unit_string = format_units(to_plot@var_units)
    end if
    pres@lbTitleString = "(" + unit_string + ")"
  end if
  pres@lbTitleFontHeightF = 0.017
  pres@lbTitleDirection = "across"
  pres@lbTitlePosition = "Right"
  pres@lbTitleJust = "CenterLeft"
  pres@lbLabelFontHeightF = 0.014
  pres@lbLabelJust = "CenterCenter"
  pres@lbLabelAutoStride = True

  pres@txFontHeightF = 0.02
  if (isatt(diag_script_info, "title")) then
    pres@txString  = diag_script_info@title
  end if
  pres@pmLabelBarParallelPosF = 0.06
  pres@pmLabelBarWidthF = .6

  outfile = panelling(wks, maps, nvert, nhori, pres)

  log_info("Wrote " + outfile)
end
