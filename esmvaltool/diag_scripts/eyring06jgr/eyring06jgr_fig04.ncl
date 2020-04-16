; ############################################################################
; EYRING06JGR_FIG04
; Author:  Irene Cionni (ENEA, Italy) , Veronika Eyring (DLR, Germany)
; PROJECT-NAME CCMI
; ############################################################################
; Description
; Time series of monthly mean anomalies respect
; cliamtological mean for selected month and latitudinal region
;
; Required diag_script_info attributes (diagnostics specific)
;    latmin: minimum latitude for region selection
;    latmax: maximum latitude for region selection
;    season: season  when variable is averaged (i.e. "DJF","MAM","JJA","SON")
;    multimean: (i.e. False/True)
;    reference_yearmax: min year for calculating the mean for the anomalies
;    reference_yearmax: max year for calculating  the mean for the anomalies
;    show_trend: (i.e. True/False) ; write the  trend in the legend
;    trend_yearmin: min year for the trend
;    trend_yearmax: max year for the trend
;
; Optional diag_script_info attributes (diagnostic specific)
;  XMin : X axis min
;  XMax : X axis max
;  YMin : Y Axis min
;  YMax : Y Axis max

; ############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  ; Get variables and datasets
  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD1 = ListCount(info0)
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")
  y1 = metadata_att_as_array(info0, "start_year")
  y2 = metadata_att_as_array(info0, "end_year")

  if(isatt(variable_info, "units")) then
    units = variable_info@units
  else
    units = "K"
  end if

  if(isatt(variable_info, "long_name")) then
    long_name = variable_info@long_name
  else
    long_name = var0
  end if

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "season", "ANN")
  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "latmin", -90.)
  set_default_att(diag_script_info, "latmax", 90.)
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "reference_yearmin", max(y1))
  set_default_att(diag_script_info, "reference_yearmax", min(y2))
  set_default_att(diag_script_info, "show_trend", "False")
  set_default_att(diag_script_info, "xMin", 1960)
  set_default_att(diag_script_info, "xMax", 2001)
  set_default_att(diag_script_info, "yMin", -20)
  set_default_att(diag_script_info, "yMax", 20)

  year_ref_max = toint(diag_script_info@reference_yearmax)
  year_ref_min = toint(diag_script_info@reference_yearmin)
  seasons = diag_script_info@season
  latrange = new((/dimsizes(seasons), 2/), "float")
  latrange(:, 0) = diag_script_info@latmin
  latrange(:, 1) = diag_script_info@latmax
  seasons = diag_script_info@season
  trend = diag_script_info@show_trend

  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

  flag_mod = where(projects.ne."OBS" .and. projects.ne."OBS6" \
                   .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  XMax = diag_script_info@xMax
  XMin = diag_script_info@xMin
  YMax = diag_script_info@yMax
  YMin = diag_script_info@yMin

  multimean = diag_script_info@multimean
  if (multimean)
    dimTot = dim_MOD+dim_OBS+1
  else
    dimTot = dim_MOD+dim_OBS
  end if

  ; Create output directories
  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

end


begin
  plotpath = config_user_info@plot_dir + "eyring06f04_" + var0
  wks = gsn_open_wks(file_type, plotpath)
  plotpath1 = config_user_info@plot_dir + "eyring06f04_legend_" + var0
  wks1 = gsn_open_wks(file_type, plotpath1)
  output_filename = wks@fullname
  nPanel = dimsizes(seasons)
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  min_y1 = min(y1)
  max_y2 = max(y2)
  years = ispan(min_y1, max_y2, 1)
  years!0 = "years"
  years&years = years
  nyears = dimsizes(years)

  tmpData = new((/dim_MOD1, nPanel, nyears, 2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "years"
  tmpData!3 = "stat"
  tmpData&years = years

  Xdata = new((/dimTot, nPanel, nyears, 2/), "double")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "years"
  Xdata!3 = "stat"
  Xdata&years = years

  Xref = new((/ nPanel, nyears, 2/), "double")
  Xref!0 = "season"
  Xref!1 = "years"
  Xref!2 = "stat"
  Xref&years = years

  plotMod = new(dimTot, typeof(names))
  plotColor = new(dimTot, typeof(colors))
  plotDash = new(dimTot, typeof(dashes))
  plotThicks = new(dimTot, "float")
  plotScatters = new(dimTot, typeof(markers))
  LineScatters = new(dimTot, typeof(names))
  ik = 0

  ; Start datasets loop
  do imod = 0, dim_MOD1 - 1

    ; Start seasons loop
    do ipan = 0, dimsizes(seasons) - 1

      ; Extract data
      A0 = read_data(info0[imod])

      ; Select latitude range
      coslat = cos(A0&lat * rad)
      coslat!0 = "lat"
      coslat&lat = A0&lat
      dummy1 = dim_avg_wgt_Wrap(A0(:, {latrange(ipan, 0) : latrange(ipan, \
                                1)}), coslat({latrange(ipan, 0) : \
                                              latrange(ipan, 1)}), 1)
      delete(coslat)
      delete(A0)
      if(seasons(ipan).ne."ANN") then
        dummy2a = month_to_season(dummy1, seasons(ipan))
        ; 3 year smoothing
        dummy2 = runave(dummy2a, 3, 0)
        dummy2_avg = time_operations(dummy1, year_ref_min, year_ref_max, \
                                     "average", seasons(ipan), True)
        delete(dummy2a)
        delete(dummy1)
      else
        dummy2 = month_to_annual(dummy1, 1)
        dummy2_avg = time_operations(dummy1, year_ref_min, year_ref_max, \
                                     "average", "annualclim", True)
        delete(dummy1)
      end if
      tmp = dummy2
      ano_tmp = (/dummy2-dummy2_avg/)
      delete(dummy2)
      delete(dummy2_avg)
      tmpData(imod, ipan, {y1(imod):y2(imod)}, 0) = (/ano_tmp/)
      tmpData(imod, ipan, {y1(imod):y2(imod)}, 1) = (/tmp/)
      delete(tmp)
      delete(ano_tmp)
    end do      ; end seasons loop

    Xdata(imod, :, :, :) = (/tmpData(imod, :, :, :)/)
    plotMod(imod) = names(imod)
    plotColor(imod) = colors(imod)
    plotDash(imod) = dashes(imod)
    plotThicks(imod) = thicks(imod)
  end do    ; end loop model

  plotScatters(index_mod) = (/0/)
  plotScatters(index_obs) = (/1/)
  LineScatters(index_obs) = (/"Markers"/)
  LineScatters(index_mod) = (/"Lines"/)

  ; Calculate multimodel mean
  if (multimean) then
    Xdata(dimTot - 1, :, :, :) = (/dim_avg_n(Xdata(index_mod, :, :, :), 0)/)
    plotMod(dimTot - 1) = "MMM"
    plotColor(dimTot - 1) = "Black"
    plotDash(dimTot - 1) = 0
    LineScatters(dimTot - 1) = "Lines"
    plotThicks(dimTot - 1) = 4.
  end if
  X_MMM_std = (/dim_stddev_n(Xdata(index_mod, :, :, :), 0)/)

  ; Calculate trend
  if (trend.eq."True") then
    year_Min_trend = diag_script_info@trend_yearmin
    year_Max_trend = diag_script_info@trend_yearmax
    trend_Xdata = regCoef(years({year_Max_trend : year_Min_trend}), \
                          Xdata(line|:, season|:, {years|year_Max_trend : \
                                                   year_Min_trend}, stat|0))
    plotMod1 = plotMod + " trend"
  end if

  ; Calculate multimodel mean errors
  xp_M = new((/2 * dimsizes(Xdata&years), dimsizes(seasons)/), typeof(Xdata))
  yp_M = new((/2 * dimsizes(Xdata&years), dimsizes(seasons)/), typeof(Xdata))
  if (multimean.eq."True") then
    do j = 0, dimsizes(years) - 1
      yp_M(j, :) = Xdata(dimTot - 1, :, j, 0) - ((X_MMM_std(:, j, 0)))
      xp_M(j, :) = tofloat(years(j))
      xp_M(dimsizes(years)+j, :) = tofloat(years(dimsizes(years) - 1 - j))
      yp_M(dimsizes(years)+j, :) = Xdata(dimTot-1, :, \
                                         (dimsizes(years) - 1 - j), 0) + \
                                        (X_MMM_std(:, (dimsizes(years)\
                                         - 1 - j), 0))
    end do
  end if

  ; Plotting
  plot_d = new(dimsizes(seasons), "graphic")
  plot_xp_R_d = new(dimsizes(seasons), "graphic")
  plot_xp_M_d = new(dimsizes(seasons), "graphic")

  do ipanel = 0, dimsizes(seasons) - 1
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tfPolyDrawOrder = "Predraw"
    res@trXMaxF = XMax
    res@trXMinF = XMin
    res@trYMaxF = YMax(ipanel)
    res@trYMinF = YMin(ipanel)
    res@tmLabelAutoStride = True
    res@vpHeightF = 0.3
    res@vpWidthF = 0.9
    res@tiYAxisString = "Anomaly " + "(" + units + ")"
    res@tiMainFontHeightF = 0.023
    res@tiMainString = seasons(ipanel) + " " + \
      lat_names(latrange(ipanel, 0)) + "-" + \
      lat_names(latrange(ipanel, 1))
    res@xyMarkLineModes = LineScatters
    res@xyLineThicknessF = 3.0
    res@xyLineDashSegLenF = 0.13
    res@xyLineColors = plotColor
    res@xyMarkerColors = plotColor
    res@xyDashPatterns = plotDash
    res@xyMarkers = plotDash
    res@xyMarkerSizeF = 0.004
    res@xyMarkerThicknessF = 1.8
    plot_d(ipanel) = gsn_csm_xy(wks, years, Xdata(line|:, season|ipanel, \
                                                  years|:, stat|0), res)
    add_line(wks, plot_d(ipanel), (/years(0), years(dimsizes(years)-1)/), \
                                  (/0., 0./), "gray", 1)
    gsres = True
    gsres@gsFillColor = "Snow3"
    if (multimean.eq."True") then
      gsres@gsFillColor = "cyan"
      plot_xp_M_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel),\
                                            xp_M(:, ipanel), yp_M(:, ipanel),\
                                            gsres)
    end if
  end do   ; end loop season

  ; Call legend-plotting routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@gsnPanelFigureStrings = plot_num(0 : dimsizes(seasons) - 1)
  Pres@gsnPanelFigureStringsFontHeightF = 0.025
  Pres@gsnPanelFigureStringsPerimOn = "False"
  Pres@amJust = "TopLeft"

  ; Call panel routine
  if (dimsizes(seasons).gt.1) then
    gsn_panel(wks, plot_d, (/dimsizes(seasons), 1/), Pres)
  else
    gsn_panel(wks, plot_d, (/1, 1/), Pres)
  end if

  plotThicks = plotThicks+1
  frame(wks)
  if (trend.eq."True") then
    txid = new(dimsizes(seasons), "graphic")
    do ipanel = 0, dimsizes(seasons) - 1
      just_strs = seasons(ipanel) + " " + lat_names(latrange(ipanel, 0))\
                                + "-" + lat_names(latrange(ipanel, 1))
      xp = 0.5
      yp = 0.9 - ipanel * 0.3
      txres = True
      txres@gsnDraw = False
      txres@gsnFrame = False
      txres@txFontHeightF = 0.015
      gsn_text_ndc(wks1, just_strs, xp, yp, txres)
      plotMod1 = plotMod + " =  " + sprintf("%5.2f", 10 * trend_Xdata(:,\
                                            ipanel))
      pos = yp - 0.1
      add_legend(wks1, plotMod1, plotColor, tofloat(plotDash), LineScatters,\
                 plotThicks, tostring(pos))

    end do
      frame(wks1)
  else
    pos = 0.9

    add_legend(wks1, plotMod1, plotColor, tofloat(plotDash), LineScatters,\
               plotThicks, tostring(pos))
    frame(wks1)
  end if

  ; File output nc
  ; Time series
  nc_filename = config_user_info@work_dir + "eyring06f04_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f04"
  Xdata&stat = (/"anomalies", "raw"/)
  Xdata&season = seasons
  Xdata&line = plotMod
  Xdata&years = years
  nc_outfile = ncdf_write(Xdata, nc_filename)
  ; Trends
  if (trend.eq."True") then
    nc_filename1 = config_user_info@work_dir + \
                   "eyring06f04_trend" + var0 + ".nc"
    trend_Xdata@var = var0
    trend_Xdata@diag_script = "eyring06f04"
    trend_Xdata!0 = "line"

    trend_Xdata!1 = "season"
    trend_Xdata&line = plotMod
    trend_Xdata&season = seasons
    nc_outfile1 = ncdf_write(trend_Xdata, nc_filename1)
  end if

  ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Time series of " + var0 + " at select level", \
                 (/"anomaly", "other"/), \
                 (/"nhpolar", "shpolar", "global"/), \
                 "times", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))

  log_provenance(nc_filename1, \
                 plotpath + "." + file_type, \
                 "Trend of " + var0 + " at select level", \
                 "trend", \
                 (/"nhpolar", "shpolar", "global"/), \
                 "bar",\
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))
end
