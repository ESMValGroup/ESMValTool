; #############################################################################
; Plot global (or regional) timeseries for historical and RCPs 1850-2?00
; Author: Ruth Lorenz (ETH, Switzerland)
; PROJECT-NAME CRESCENDO
; #############################################################################
;
; Description
;    Plot global or regional timeseries in historical and RCPs
;    1850-2?00 for all available models precalculated by
;    ts_line_collins_ipcc_ch12.ncl
;
; Caveats
;    Variable number of scenarios? So far hardcoded 4 scenarios (rcp2.6,
;    rcp4.5, rcp6.0, rcp8.5)
;
; Modification history
;    20171002-A_lore_ru: Added tags for searching
;    20170523-A_lore_ru: separated plotting from calculations
;    20160625-A_lore_ru: adapted to ESMValTool in ts_line_collins_ipcc_ch12
;    20130506-A_sedl_ja: written for IPCC AR5.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
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
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
  write_references(diag_script, \  ; script name
                   "A_lore_ru", \  ; authors
                   "", \  ; contributors
                   "D_collins13ipcc",    \  ; diag_references
                   "",          \  ; obs_references
                   "P_crescendo")   ; proj_references
  ; Check required diag_script_info attributes
  req_atts = (/"syears", "eyears"/)
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

  ; Output netcdf directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
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
    log_debug("MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; climofiles = new(dim_MOD, string)
  ; do ii = 0, dim_MOD - 1
  ;   climofiles(ii) = \
  ;     interface_get_inpaths(ii) + "/" \
  ;     + interface_get_infile(variables(0), field_types(0), ii)
  ; end do

  ; ###########################################
  ; # Get data                                #
  ; ###########################################
  datapath = netcdf_dir + "/" + variable_info@diagnostic + ".nc"
  to_plot = ncdf_read(datapath, var0)
  if(isatt(diag_script_info, "model_nr") .and. \
     diag_script_info@model_nr .eq. True) then
    datapath2 = netcdf_dir + "/nr_runs_" + variable_info@diagnostic + ".nc"
    runs = ncdf_read(datapath2, "nr")
  end if
  ; start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)

  dims_data = (dimsizes(to_plot&metric) - 1) / 3

  ; ###########################################
  ; *********** PLOTTING ************
  ; Create plot variables
  outfile = var0 + "_ts_line_" + syears(0) + "_" + eyears(nperiods - 1)
  wks = get_wks(0, diag_script, outfile)
  if(isatt(diag_script_info, "colormap")) then
    colormap = RGBtoCmap(diag_script_info@colormap)
  else
    colormap = "StepSeq25"
  end if
  gsn_define_colormap(wks, colormap)

  colors_main = (/1, 3, 4, 5, 6/)
  colors_spread = (/2, 3, 4, 5, 6/)

  res                   = True

  res@gsnDraw = False
  res@gsnFrame = False
  res@txFont   = 25
  res@vpHeightF = 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.8
  res@vpXF     = 0.15                   ; start plot at x ndc coord

  if (isatt(diag_script_info, "title")) then
    res@tiMainString = diag_script_info@title
  end if
  if(isatt(diag_script_info, "yaxis")) then
    res@tiYAxisString = diag_script_info@yaxis
  end if
  res@tiXAxisString = "Year"

  res@xyDashPatterns     = 0
  res@trXMinF = toint(syears(0))
  res@trXMaxF = toint(eyears(nperiods - 1))
  if(isatt(diag_script_info, "ymin")) then
    res@trYMinF = diag_script_info@ymin
  end if
  if(isatt(diag_script_info, "ymax")) then
    res@trYMaxF = diag_script_info@ymax
  end if

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 50)
  res@tmXBLabels = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 50)
  res@tmXBMinorValues  = ispan(toint(syears(0)), \
                               toint(eyears(nperiods - 1)), 10)

  res@xyLineThicknesses = 3.5
  res@xyLineColor      = colors_main(0)

  ; create plot
  plot  = gsn_csm_xy(wks, to_plot&year, to_plot(0, :), res)

  res@gsnXYFillColors = colors_spread(0)
  res@xyLineColor     = -1
  spread_plot  = gsn_csm_xy(wks, to_plot&year, to_plot(1 : 2, :), res)
  overlay(plot, spread_plot)

  delete(res@gsnXYFillColors)

  do j = 1, dims_data
    if(all(ismissing(to_plot(j * 3, :)))) then
      error_msg("w", diag_script, "", "all missing values for " + \
                "dimension " + j * 3 + " in to_plot, need at least 5" + \
                " models, go to next")
      continue
    end if
    res@xyLineColor      = colors_main(j)
    plot1  = gsn_csm_xy(wks, to_plot&year, to_plot(j * 3, :), res)
    overlay(plot, plot1)

    res@gsnXYFillColors = colors_spread(j)
    res@gsnXYFillOpacities = (/0.2, 0.2, 0.2/)
    res@xyLineColor     = -1
    spread_plot = gsn_csm_xy(wks, to_plot&year,\
                             to_plot((j * 3) + 1 : (j * 3) + 2, :), res)
    overlay(plot, spread_plot)
    delete(res@gsnXYFillColors)
    delete(res@gsnXYFillOpacities)
  end do

  polyres                  = True
  polyres@gsLineThicknessF = 1.5
  ; add polylines
  dum = new(nperiods, graphic)
  do n = 0, nperiods - 1
    dum(n) = gsn_add_polyline(wks, plot, (/toint(eyears(n)), \
                              toint(eyears(n))/), (/-20, 20/), polyres)
  end do
  delete(res@xyLineColor)

  ; Attach a legend
  lgres                    = True
  lgres@lgLineColors       = colors_main  ; (/1, 3, 4, 5, 6/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.15           ; width of legend (NDC)
  lgres@vpHeightF          = 0.16           ; height of legend (NDC)
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgLineThicknessF        = 4
  lgres@lgPerimOn = False
  lgres@lgItemOrder = (/4, 3, 2, 1, 0/)
  if (isatt(diag_script_info, "label")) then
    labels = "  " + diag_script_info@label
    legend = gsn_create_legend(wks, 5, labels, lgres)
    amres                  = True
    amres@amJust           = "TopLeft"    ; Use bottom right corner of box
                                              ; for determining its location.
    amres@amParallelPosF   = -0.46              ; Move legend to right
    if (var0 .ne. "snw") then
      amres@amParallelPosF   = -0.45              ; Move legend to right
      amres@amOrthogonalPosF = -0.45              ; Move legend up.
    else
      amres@amParallelPosF   = -0.47              ; Move legend to right
      amres@amOrthogonalPosF = 0.12               ; Move legend down.
    end if
    annoid = gsn_add_annotation(plot, legend, amres)  ; add legend to plot
  end if
  if(isatt(diag_script_info, "model_nr")) then
    getvalues plot                     
      "tmYLLabelFontHeightF"   : fheight
      "trYMinF"                : trYMinF 
      "trYMaxF"                : trYMaxF
    end getvalues
    text = new(1 + (dims_data - 1) * nperiods, graphic)
    txres               = True
    txres@txFontHeightF = fheight
    label = tostring(runs(0, 0))
    xloc = (toint(syears(0)) + (toint(eyears(0)) - toint(syears(0))) / 2)
    yloc = (trYMaxF + trYMinF) / 2
    ydiff = (trYMaxF + trYMinF) / 7
    if yloc .le. 0 then
      ydiff = - ydiff
    end if
    text(0) = gsn_add_text(wks, plot, label, xloc, yloc, txres)
    i = 1
    do nscen = 1, dims_data
      do nper = 1, nperiods - 1
        txres@txFontColor = colors_main(nscen)
        label = tostring(runs(nscen, nper))
        xloc = (toint(eyears(nper)) - toint(syears(nper))) / 2
        text(i) = gsn_add_text(wks, plot, label, toint(syears(nper)) + xloc, \
                               yloc - nscen * ydiff, txres)
        i = i + 1
      end do
    end do
  end if
  draw(wks)
  frame(wks)
  ; alltags = array_append_record(tags,\
  ;   (/"DM_glob", "PT_times", "ST_mean", "ST_stddev", "ST_anomaly"/), 0)
  ; string1 = "Time series of global annual mean  " + variable_info@long_name +\
  ;           " anomalies (relative to " + dataset_info@start_year(0) + "-" + \
  ;           dataset_info@end_year(0) + ") from CMIP5 concentration drive " + \
  ;           "experiments. Projections are shown for each RCP for the " + \
  ;           "multi-model mean (solid lines) and the 5 to 95% range " + \
  ;           "(+-1.64 standar deviation) across the distribution of " + \
  ;           "individual models (shading). "
  ; if (nperiods .gt. 1) then
  ;   string2 = "Discontinuities at 2100 are due to different numbers of" + \
  ;     "models performing the extension runs beyond the 21st " + \
  ;     "century and have no physical meaning. Only one ensemble " + \
  ;     "member is used from each model and numbers in the " + \
  ;     "figure indicate the number of different models " + \
  ;     "contributing to the different time periods. No ranges " + \
  ;     "are given for the RCP6.0 projections beyond 2100 as only" + \
  ;     " 2 models are available."
  ; else
  ;   string2 = ""
  ; end if
  ; caption = string1 + string2
  ; id = diag_script + "_" + variables(0)

  ; contrib_authors = (/"A_lore_ru", "A_sedl_ja"/)
  ; ESMValMD(wks@fullname, alltags, \
  ;          caption, id, variables(0), dataset_info@dataset, \
  ;          climofiles, diag_script, contrib_authors)

  ;delete([/alltags, caption, id, climofiles/])

  log_info(" Wrote " + plot_dir + "/" + outfile + "." + file_type)
  leave_msg(diag_script, "")
end
