; #############################################################################
; Plot global timeseries for historical and SSP with constrained projections
; in 2100
; Author: Kevin Debeire (DLR, Germany)
; #############################################################################
;
; Description
; Plot global timeseries of TAS and bar plots of constrained projections from
; predefined csv files.
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
begin
  enter_msg(DIAG_SCRIPT, "")

  ; 'models', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

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
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data                                #
  ; ###########################################
  datapath = netcdf_dir + "/" + variable_info[0]@diagnostic + ".nc"
  to_plot = ncdf_read(datapath, var0)
  if(isatt(diag_script_info, "model_nr") .and. \
     diag_script_info@model_nr .eq. True) then
    datapath2 = netcdf_dir + "/nr_runs_" + variable_info[0]@diagnostic + ".nc"
    runs = ncdf_read(datapath2, "nr")
  end if
  ; start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)

  dims_data = (dimsizes(to_plot&metric) - 1) / 3

  ; Constrained projections lower bound, mean and upper bound
  constrained_projections_2100 = new((/dims_data, 4/), float)

  lower_bound = diag_script_info@lower_constrained_projections
  upper_bound = diag_script_info@upper_constrained_projections
  mean_bound = diag_script_info@mean_constrained_projections

  do nscen = 0, dims_data - 1
    constrained_projections_2100(nscen, :) = (/2100, lower_bound(nscen), \
                                               mean_bound(nscen), \
                                               upper_bound(nscen) /)
  end do

  ; ###########################################
  ; *********** PLOTTING ************
  ; Create plot variables
  outfile = var0 + "_ts_line_" + syears(0) + "_" + eyears(nperiods - 1)
  wks = get_wks(0, DIAG_SCRIPT, outfile)
  if(isatt(diag_script_info, "colormap")) then
    colormap = RGBtoCmap(diag_script_info@colormap)
  else
    colormap = "StepSeq25"
  end if
  gsn_define_colormap(wks, colormap)

  colors_main = (/1, 4, 5, 7, 8/)
  colors_spread = (/2, 4, 5, 7, 8/)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@txFont = 22
  res@txFontHeightF = 1.8
  res@vpHeightF = 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.7
  res@vpXF     = 0.15                   ; start plot at x ndc coord
  caption = ""
  if (isatt(diag_script_info, "title")) then
    res@tiMainString = diag_script_info@title
    caption = diag_script_info@title
    res@tiMainFontHeightF = 0.04
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
  res@xyLineColor  = colors_main(0)
  res@tmYROn= False
  res@tmYRBorderOn = False
  ; create second axis
  diff_avg = diag_script_info@baseline_offset
  res2 = True
  res2@xyLineThicknesses = 0.001
  res2@xyLineColor  = -1
  res2@tiYAxisString = "Relative to 1850-1900 (~S~o~N~ C)"
  if(isatt(diag_script_info, "ymin")) then
    res2@trYMinF = diff_avg + res@trYMinF
  end if
  if(isatt(diag_script_info, "ymax")) then
    res2@trYMaxF = res@trYMaxF+diff_avg
  end if

  ; create plot
  plot  = gsn_csm_xy(wks, to_plot&year, to_plot(0, :), res)

  res@gsnXYFillColors = colors_spread(0)
  res@xyLineColor     = -1
  spread_plot  = gsn_csm_xy(wks, to_plot&year, to_plot(1 : 2, :), res)
  overlay(plot, spread_plot)

  delete(res@gsnXYFillColors)
  delete(res@tiMainString)
  do j = 1, dims_data
    if(all(ismissing(to_plot(j * 3, :)))) then
      error_msg("w", DIAG_SCRIPT, "", "all missing values for " + \
                "dimension " + j * 3 + " in to_plot, need at least 5" + \
                " models, go to next")
      continue
    end if
    res@xyLineOpacities = 1
    res@xyLineColor      = colors_main(j)
    plot1  = gsn_csm_xy(wks, to_plot&year, to_plot(j * 3, :), res)
    overlay(plot, plot1)

    res@gsnXYFillColors = colors_spread(j)
    res@gsnXYFillOpacities = (/0.1, 0.1, 0.1/)
    res@xyLineColor     = -1
    if (j .eq. 1) then
        spread_plot = gsn_csm_xy2(wks, to_plot&year, to_plot((j * 3) + 1 \
                                  : (j * 3) + 2, :), to_plot((j * 3) + 1 \
                                  : (j * 3) + 2, :) * 0, res, res2)
    else
        spread_plot = gsn_csm_xy(wks, to_plot&year, \
                                 to_plot((j * 3) + 1 : (j * 3) + 2, :), res)
    end if
    overlay(plot, spread_plot)
    delete(res@gsnXYFillColors)
    delete(res@gsnXYFillOpacities)
  end do
  ; plot a plyline at the start of scenario
  polyres                  = True
  polyres@gsLineThicknessF = 1.5
  ; add polylines
  dum = new(nperiods, graphic)
  do n = 0, nperiods - 1
    gsn_polyline(wks, plot, (/toint(eyears(n)), \
                 toint(eyears(n))/), (/-20, 20/), polyres)
  end do
  delete(res@xyLineColor)
  delete(dum)
  ; plot a polyline indicating constrained projection
  polyres2 = True
  polyres2@gsLineThicknessF = 7
  polyres2@gsLineOpacityF = 0.9
  polyres2@gsLineColor      = colors_main(3)
  dumm = new(dims_data*4, graphic)

  do idx = 0, dims_data - 1
    polyres2@gsLineColor      = colors_main(idx+1)
    dumm(idx*4) = gsn_add_polyline(wks, plot, (/toint(eyears(nperiods - 1)-10)\
                                   , toint(eyears(nperiods - 1))-10/)\
                                   , (/constrained_projections_2100(idx, 1)\
                                   , constrained_projections_2100(idx, 3)/)\
                                   , polyres2)
    dumm(idx*4+1) = gsn_add_polyline(wks, plot, (/toint(eyears(nperiods - 1))\
                                     , toint(eyears(nperiods - 1)-19)/)\
                                     , (/constrained_projections_2100(idx, 2)\
                                     , constrained_projections_2100(idx, 2)/)\
                                     , polyres2)
    dumm(idx*4+2) = gsn_add_polyline(wks, plot, (/toint(eyears(nperiods - 1))\
                                     , toint(eyears(nperiods - 1)-19)/)\
                                     , (/constrained_projections_2100(idx, 3)\
                                     , constrained_projections_2100(idx, 3)/)\
                                     , polyres2)
    dumm(idx*4+3) = gsn_add_polyline(wks, plot, (/toint(eyears(nperiods - 1))\
                                     , toint(eyears(nperiods - 1)-19)/)\
                                     , (/constrained_projections_2100(idx, 1)\
                                     , constrained_projections_2100(idx, 1)/)\
                                     , polyres2)
  end do

  ; Attach a legend
  lgres                    = True
  lgres@lgLineColors       = colors_main  ; (/1, 3, 4, 5/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.2           ; width of legend (NDC)
  lgres@vpHeightF          = 0.2           ; height of legend (NDC)
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgLineThicknessF        = 5
  lgres@lgPerimOn = False
  lgres@lgItemOrder = (/3, 2, 1, 0/)
  if (isatt(diag_script_info, "label")) then
    labels = "  " + diag_script_info@label
    legend = gsn_create_legend(wks, 4, labels, lgres)
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
    txres@txFontHeightF = fheight * 0.8
    label = tostring(runs(0, 0))
    xloc = 1938
    yloc = trYMaxF - (trYMaxF-trYMinF)/9
    ydiff = (trYMaxF - trYMinF) / 10.9
    if yloc .le. 0 then
      ydiff = - ydiff
    end if
    text(0) = gsn_add_text(wks, plot, label, xloc, yloc, txres)
    i = 1
    do nscen = 1, dims_data
      do nper = 1, nperiods - 1
        txres@txFontColor = colors_main(nscen)
        label = tostring(runs(nscen, nper))
        text(i) = gsn_add_text(wks, plot, label, xloc, \
                               yloc - nscen * ydiff, txres)
        i = i + 1
      end do
    end do
  end if
  draw(wks)
  frame(wks)

  log_info(" Wrote " + plot_dir + "/" + outfile + "." + file_type)

  ; collect meta-data
  nc_file = datapath
  plot_file = outfile

  if (isatt(diag_script_info, "spread")) then
    spread_str = diag_script_info@spread
  else
    spread_str = "1.0"
  end if
  statistics = ("mean")
  domains = ("global")
  plot_types = ("times")
  authors = (/"debeire_kevin"/)
  references = (/"tebaldi21esd"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, plot_file, caption, statistics, domains, \
                 plot_types, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")
end
