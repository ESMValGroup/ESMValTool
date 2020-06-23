; ############################################################################
; SCRIPT TO REPRODUCE FIG. 15 OF EYRING ET AL. (JGR, 2006)
; Authors:  Irene Cionni (ENEA, Italy) and Veronika Eyring (DLR, Germany)
; CCMI Project
; ############################################################################
;
; Description:
;    Time series of seasonal  mean anomalie of total column of ozones
;    calculated with a methodology described in Appendix A of Eyring06
;    Annual cycle for the reference years are also shown
;   Ad more bullets if required
;
; Required diag_script_info attributes: (diagnostics specific)
;    latrange:  minimum and maximum latitude for regional average,
;               e.g. (/60, 90/))
;    season: season for temporal average (e.g, "DJF", "MAM", "JJA", "SON")
;    multimean: calculate multi-model mean (i.e. False/True)
;    refModel:  name of the reference model or observation for
;                         bias calculation (e.g. (/"ERA-Interim"/))
;
; Optional diag_script_info attributes:
;    latmin: minimum latitude for region selection
;    latmax: maximum latitude for region selection
;    Xmin: lower limit for the X axis
;    Xmax: upper limit for the X axis
;    Ymin: lower limit Y axis
;    Ymax: upper limit Y axis
;    start_year: start year for the climatology calculation
;    end_year: end year for the climatology calculation
;    Cycle_Ymin: lower limit Y axis (for right figures)
;    Cycle_Ymax: upper limit Y axis (for right figures)
;
; Required variable_info attributes:
;
; Optional variable_info attributes:
;
; Modification history
;    20200323-hassler_birgit: cleaned code
;    20180320-lauer_axel: added tags for reporting and call to write_references
;    20170315-cionni_irene: written
;
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

  ; Check variable attributes
  if (isatt(variable_info, "units")) then
    UNITS = variable_info@units
  else
    UNITS = "DU"
  end if

  if (isatt(variable_info, "long_name")) then
    LONG_NAME = variable_info@long_name
  else
    LONG_NAME = var0
  end if

  ; Set default values for optional diag_script_info attributes
  set_default_att(diag_script_info, "season", "ANN")
  set_default_att(diag_script_info, "latmin", -90.)
  set_default_att(diag_script_info, "latmax", 90.)
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "cycle_yearmin", max(y1))
  set_default_att(diag_script_info, "cycle_yearmax", min(y2))
  set_default_att(diag_script_info, "Xmin", 1960)
  set_default_att(diag_script_info, "Xmax", 2005)
  set_default_att(diag_script_info, "Ymin", -10)
  set_default_att(diag_script_info, "Ymax", 10)
  set_default_att(diag_script_info, "cycle_Ymin", 260.)
  set_default_att(diag_script_info, "cycle_Ymax", 380.)

  year_cycle_max = toint(diag_script_info@cycle_yearmax)
  year_cycle_min = toint(diag_script_info@cycle_yearmin)

  seasons = diag_script_info@season
  latrange = new((/dimsizes(seasons), 2/), "float")
  latrange(:, 0) = diag_script_info@latmin
  latrange(:, 1) = diag_script_info@latmax

  ; Styles
  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

  flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  XMax = diag_script_info@Xmax
  XMin = diag_script_info@Xmin
  YMax = diag_script_info@Ymax
  YMin = diag_script_info@Ymin
  cycle_max = toint(diag_script_info@cycle_Ymax)
  cycle_min = toint(diag_script_info@cycle_Ymin)

  multimean = diag_script_info@multimean
  if (multimean) then
    dimTot = dim_MOD + dim_OBS + 1
  else
    dimTot = dim_MOD + dim_OBS
  end if

  ; make sure path for (mandatory) netcdf output exists
  work_dir = config_user_info@work_dir + "/"

  ; Create output directory
  system("mkdir -p " + work_dir)

    ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

end

begin

  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f15_" + var0)
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

  tmpData = new((/dim_MOD1, nPanel, nyears, 2/), "float")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "years"
  tmpData!3 = "stat"
  tmpData&years = years

  Xdata = new((/dimTot, nPanel, nyears, 2/), "float")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "years"
  Xdata!3 = "stat"
  Xdata&years = years
  Xref = new((/ nPanel, nyears, 2/), "float")

  Xref!0 = "season"
  Xref!1 = "years"
  Xref!2 = "stat"
  Xref&years = years

  cytmpData = new((/dim_MOD + dim_OBS, nPanel, 366, 2/), "float")
  cytmpData!0 = "line"
  cytmpData!1 = "season"
  cytmpData!2 = "day"
  cytmpData!3 = "stat"
  cytmpData&day = ispan(1, 366, 1)

  cyXdata = new((/dimTot, nPanel, 366, 2/), "float")
  cyXdata!0 = "line"
  cyXdata!1 = "season"
  cyXdata!2 = "day"
  cyXdata!3 = "stat"
  cyXdata&day = ispan(1, 366, 1)

  plotMod = new(dimTot, typeof(names))
  plotColor = new(dimTot, typeof(names))
  plotDash = new(dimTot, typeof(names))
  plotThicks = new(dimTot, "float")
  plotScatters = new(dimTot, typeof(names))
  LineScatters = new(dimTot, typeof(names))
  ik = 0

  ; Start datasets loop
  do imod = 0, dim_MOD1 - 1

    ; Start seasons loop
    do ipan = 0, dimsizes(seasons) - 1

      ; Extract data
      A0 = read_data(info0[imod])

      if (all(isdim(A0, (/"time", "lat", "lon"/)))) then
        dummy1 = area_operations(A0(time|:, lat|:, lon|:), \
                                 latrange(ipan, 0), \
                                 latrange(ipan, 1), 0., 360., \
                                 "average", True)
        coslat = cos(A0&lat * rad)
        coslat!0 = "lat"
        coslat&lat = A0&lat
        delete(coslat)
        delete(A0)
      else
        if (all(isdim(A0, (/"time", "lat"/)))) then
          coslat = cos(A0&lat * rad)
          coslat!0 = "lat"
          coslat&lat = A0&lat
          dummy1 = \
            dim_avg_wgt_Wrap(A0(:, \
                             {latrange(ipan, 0):latrange(ipan, 1)}), \
                             coslat({latrange(ipan, 0):latrange(ipan, 1)}), 1)
          delete(coslat)
          delete(A0)
        end if
      end if

      ; Appendix A methodology
      ; (function calcRegCoeffs;
      ; function calcMeanAnnCycleAnnual
      ; procedure rmMeanAnnCycle)
      ; These functions and procedures are stored in aux_plotting.ncl
      tim = dummy1&time
      tim_min = doubletointeger(ut_inv_calendar(y1(imod), \
                                1, 1, 0, 0, 0, tim@units, 0))
      tim_max = doubletointeger(ut_inv_calendar(y2(imod), \
                                12, 31, 59, 0, 0, tim@units, 0))
      delete(tim)

      coeff_ab = calcRegCoeffs(dummy1({tim_min:tim_max}))
      coeff_a = coeff_ab(0:8)
      delete(coeff_ab)
      printVarSummary(coeff_a)
      FMeanAnnCycleDay = calcMeanAnnCycleAnnual(tofloat(coeff_a))   ; daily
      FMeanAnnCycleDay@units = dummy1@units
      FMeanAnnCycleDay@_FillValue = 10e21

      ; Anomalies
      coeff_a1 = coeff_a(0)
      FMeanAnnCycleMon = calcMeanAnnCycleMonthly(tofloat(coeff_a))  ; monthly
      FMeanAnnCycleMon@_FillValue = 10e21

      ; remove mean annual cycle and calculate anomalies
      ; set input for rmMeanAnnCycle
      nmon = 12
      ; set outputs for rmMeanAnnCycle
      dimx = dimsizes(dummy1)
      ntim = dimx(0)
      nyr = ntim/nmon

      FMonthlyAnomalies = new((/ntim/), typeof(dummy1))  ; Anomalies montlhy
      FAnnualAnomalies = new((/nyr/), typeof(dummy1))    ; Anomalies Anuual

      rmMeanAnnCycle(dummy1, FMeanAnnCycleMon, coeff_a1, \
                     FMonthlyAnomalies, FAnnualAnomalies)

      FMonthlyAnomalies@_FillValue = dummy1@_FillValue
      FAnnualAnomalies@_FillValue = dummy1@_FillValue
      FMonthlyAnomalies!0 = "time"
      FMonthlyAnomalies&time = dummy1&time
      delete(dummy1)

      if (seasons(ipan).ne."ANN") then
        FSeasAnomalies = month_to_season(FMonthlyAnomalies, seasons(ipan))
        date = ut_calendar(FSeasAnomalies&time, 0)
        year = date(:, 0)
        FSeasAnomalies!0 = "year"
        FSeasAnomalies&year = year
      else
        FSeasAnomalies = month_to_annual(FMonthlyAnomalies, 1)
        date = ut_calendar(FMonthlyAnomalies&time, 0)
        year = date(::12, 0)
        FSeasAnomalies&year = year
      end if   ; end of Appendix A methodology

      tmpData(imod, ipan, {y1(imod):y2(imod)}, 0) = (/tofloat(FSeasAnomalies)/)
      Xdata(imod, ipan, {y1(imod):y2(imod)}, 0) = (/tofloat(FSeasAnomalies)/)

      cytmpData(imod, ipan, :, 0) = (/tofloat(FMeanAnnCycleDay)/)
      cyXdata(imod, ipan, :, 0) = (/tofloat(FMeanAnnCycleDay)/)

      delete(coeff_a)
      delete(FMeanAnnCycleDay)
      delete(FSeasAnomalies)
      delete(date)
      delete(year)
      delete(FMonthlyAnomalies)
      delete(FAnnualAnomalies)
      delete(tim_min)
      delete(tim_max)

    end do  ; ipan

    plotMod(imod) = names(imod)
    plotColor(imod) = colors(imod)
    plotDash(imod) = dashes(imod)
    plotThicks(imod) = thicks(imod)
    plotScatters(imod) = markers(imod)

  end do  ; imod

  LineScatters(index_obs) = "Markers"
  ind_mon = (/0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 304, 335, 365/)
  save_obs = cyXdata(index_obs, :, ind_mon, 0)
  cyXdata(index_obs, :, :, 0) = cyXdata@_FillValue
  cyXdata(index_obs, :, ind_mon, 0) = save_obs
  save_obs = cyXdata(index_obs, :, ind_mon, 1)
  cyXdata(index_obs, :, :, 1) = cyXdata@_FillValue
  cyXdata(index_obs, :, ind_mon, 1) = save_obs
  delete(save_obs)

  ; (obs dataset is reduced to monthly data because is plotted as scatters )
  LineScatters(index_mod) = "Lines"

  ; Calculate multimodel mean
  if (multimean) then
    Xdata(dimTot-1, :, :, 0) = \
      (/dim_avg_n(Xdata(index_mod, :, :, 0), 0)/)
    Xdata(dimTot-1, :, :, 1) = \
      (/dim_stddev_n(Xdata(index_mod, :, :, 0), 0)/)
    cyXdata(dimTot-1, :, :, 0) = \
      (/dim_avg_n_Wrap(cyXdata(index_mod, :, :, 0), 0)/)
    cyXdata(dimTot-1, :, :, 1) = \
      (/dim_stddev_n_Wrap(cyXdata(index_mod, :, :, 0), 0)/)
    plotMod(dimTot-1) = "MMM"
    plotColor(dimTot-1) = "Black"
    plotDash(dimTot-1) = 0
    plotThicks(dimTot-1) = 4.
    LineScatters(dimTot-1) = "Lines"
    plotScatters(dimTot-1) = 0
  end if

  ; Plotting
  plot_d = new(dimTot, "graphic")

  plot_diff = new((/dimsizes(seasons)/), "graphic")
  plot_ann = new((/dimsizes(seasons)/), "graphic")
  plot_xp_M_d = new((/dimsizes(seasons)/), "graphic")
  plot_Mult = new((/dimsizes(seasons)/), "graphic")
  plot_ann_Mult = new((/dimsizes(seasons)/), "graphic")
  baseplot = new((/dimsizes(seasons) * 2/), "graphic")

  plot_xp_r_d = new(dimsizes(seasons), "graphic")
  plot_xp_m_d = new(dimsizes(seasons), "graphic")

  if (multimean) then
    xp_M = new((/2 * dimsizes(Xdata&years), \
                dimsizes(seasons)/), typeof(Xdata))
    yp_M = new((/2 * dimsizes(Xdata&years), \
                dimsizes(seasons)/), typeof(Xdata))
    do j = 0, dimsizes(years) - 1
      yp_M(j, :) = Xdata(dimTot-1, :, j, 0) - Xdata(dimTot-1, :, j, 1)
      xp_M(j, :) = tofloat(years(j))
      xp_M(dimsizes(years) + j, :) = tofloat(years(dimsizes(years) - 1 - j))
      yp_M(dimsizes(years) + j, :) = \
        Xdata(dimTot - 1, :, (dimsizes(years) - 1 - j), 0) + \
        Xdata(dimTot - 1, :, (dimsizes(years) - 1 - j), 1)
    end do
    xp_M2 = new((/2 * dimsizes(cyXdata&day), \
                 dimsizes(seasons)/), typeof(cyXdata))
    yp_M2 = new((/2 * dimsizes(cyXdata&day), \
                 dimsizes(seasons)/), typeof(cyXdata))
    do j = 0, 366 - 1
      yp_M2(j, :) = cyXdata(dimTot - 1, :, j, 0) - \
                    cyXdata(dimTot - 1, :, j, 1)
      xp_M2(j, :) = j + 1
      xp_M2(dimsizes(cyXdata&day) + j, :) = 366 - j
      yp_M2(dimsizes(cyXdata&day) + j, :) = \
        cyXdata(dimTot - 1, :, (dimsizes(cyXdata&day) - 1 - j), 0) + \
        cyXdata(dimTot - 1, :, (dimsizes(cyXdata&day) - 1 - j), 1)
    end do
  end if

  ; Plotting
  plotpath_anom = config_user_info@plot_dir + "eyring06f15_Anomalies_" + var0
  wks = gsn_open_wks(file_type, plotpath_anom)
  plotpath1_anom = config_user_info@plot_dir + "eyring06f01_legend_" + var0
  wks1 = gsn_open_wks(file_type, plotpath1_anom)

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
    res@vpHeightF = 0.3       ; change aspect ratio of plot
    res@vpWidthF = 0.7
    res@tiYAxisString = "Anomaly(%) "  ; y-axis label
    res@tiMainFontHeightF = 0.023
    res@tiMainString = seasons(ipanel) + " " + \
      lat_names(latrange(ipanel, 0)) + "-" + \
      lat_names(latrange(ipanel, 1))
    res@xyMarkLineModes = LineScatters
    res@xyLineThicknessF = 3.0
    res@xyLineDashSegLenF = 0.13
    res@xyLineColors = plotColor   ; line color
    res@xyDashPatterns = plotDash
    res@xyMarkers = plotScatters
    res@xyMarkerSizeF = 0.004
    res@xyMarkerThicknessF = 1.8
    plot_d(ipanel) = \
      gsn_csm_xy(wks, Xdata&years, \
                 Xdata(line|:, season|ipanel, years|:, stat|0), \
                 res)

    add_line(wks, plot_d(ipanel), \
             (/years(0), years(dimsizes(years) - 1)/), \
             (/0., 0./), "gray", 1)

    if (multimean) then
      gsres = True
      gsres@gsFillColor = "cyan"   ; color chosen
      plot_xp_M_d(ipanel) = \
        gsn_add_polygon(wks, plot_d(ipanel), \
                        xp_M(:, ipanel), \
                        yp_M(:, ipanel), gsres)
    end if

    res2 = True
    res2@gsnDraw = False
    res2@gsnFrame = False
    res2@tfPolyDrawOrder = "Predraw"
    res2@trXMaxF = 366
    res2@trXMinF = 1
    res2@trYMaxF = cycle_max(ipanel)
    res2@trYMinF = cycle_min(ipanel)
    res2@tmXBMode = "Explicit"
    res2@tmXBValues = (/1, 32, 61, 92, 122, 153, \
                       183, 214, 245, 275, 305, 336/)
    res2@tmXBLabels = (/"J", "F", "M", "A", "M", "J", \
                       "J", "A", "S", "O", "N", "D"/)

    res2@vpHeightF = 0.3    ; change aspect ratio of plot
    res2@vpWidthF = 0.3
    res2@xyMarkerSizeF = 0.004
    res2@xyMarkerThicknessF = 1.8
    ; y-axis label
    res2@tiYAxisString = "Detrended mean annual cycle (" + UNITS + ")"

    res2@xyMarkLineModes = LineScatters(:)
    res2@xyLineThicknessF = 3.0
    res2@xyLineDashSegLenF = 0.13
    res2@xyLineColors = plotColor(:)   ; line color
    res2@xyDashPatterns = plotDash(:)
    res2@xyMarkers = plotScatters(:)
    res2@xyMarkerSizeF = 0.004
    res2@xyMarkerThicknessF = 1.8
    plot_ann(ipanel) = \
      gsn_csm_xy(wks, cyXdata&day, \
                 cyXdata(line|:, season|ipanel, day|:, stat|0), \
                 res2)

    if (multimean) then
      gsres2 = True
      gsres2@gsFillColor = "cyan"
      plot_ann_Mult(ipanel) = \
        gsn_add_polygon(wks, plot_ann(ipanel), \
                        xp_M2(:, ipanel), \
                        yp_M2(:, ipanel), \
                        gsres2)
    end if

  end do   ; end ipanel

  do j = 0, dimsizes(seasons) - 1
    baseplot(2 * j) = plot_d(j)
    baseplot(2 * j + 1) = plot_ann(j)
  end do

  ; Call legend - plotting routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelXWhiteSpacePercent = 0
  Pres@gsnPanelBottom = 0.15
  Pres@gsnPanelLeft = 0.15
  gsn_panel(wks, baseplot, (/dimsizes(seasons), 2/), Pres)

  plotThicks = plotThicks + 1
  plotDash(index_obs) = plotScatters(index_obs)
  add_legend(wks1, plotMod, plotColor, tofloat(plotDash), \
             LineScatters, plotThicks, "0.9")
  frame(wks)
  frame(wks1)

  ; Anomalies output
  nc_filename = work_dir + "eyring06f15_Anomalies_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f15"
  Xdata&stat = (/"mean", "avg"/)
  Xdata&season = seasons
  Xdata&line = plotMod
  nc_outfile = ncdf_write(Xdata, nc_filename)

  ; AnnualCycle output
  nc_filename1 = work_dir + "eyring06f15_AnnualCycle_" + var0 + ".nc"
  cyXdata@var = var0
  cyXdata@diag_script = "eyring06f15"
  cyXdata&stat = (/"mean", "avg"/)
  cyXdata&season = seasons
  cyXdata&line = plotMod
  nc_outfile1 = ncdf_write(cyXdata, nc_filename1)

  ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath_anom + "." + file_type, \
                 "Anomaly time series and seasonal anomalies of " + var0, \
                 (/"anomaly", "mean"/), \
                 (/"nhpolar", "shpolar", "global"/), \
                 "pro", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
