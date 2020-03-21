; #############################################################################
; SCRIPT TO REPRODUCE FIG. 2 OF EYRING ET AL. (JGR, 2006)
; Authors: Irene Cionni (ENEA, Italy) and Veronika Eyring (DLR, Germany)
; CCMI Project
; #############################################################################
;
; Description
;    Descent of the zero zonal mean wind lines at 60°S based on the
;    climatological mean annual cycle calculated from the monthly mean zonal
;    mean winds, climatological mean for selected seasons and latitudinal
;    region
;
; Required diag_script_info attributes
;    latmin: minimum latitude for region selection
;    latmax: maximum latitude for region selection
;    monthstart: start month when variable is plotted (i.e. "OCT", "DEC")
;    monthend: end month when variable is plotted (i.e. "FEB", "MAR")
;    multimean: calculate multi-model mean (i.e. False/True)
;    refModel: name of the reference model or observation for
;              bias calculation (e.g. (/"ERA-Interim"/))
;
; Optional diag_script_info attributes
;    levmax:   max limit Y axis
;    levmin:   min limit Y axis
;    start_year: year when start the climatology calculation
;    start_year: year when end the climatology calculation
;
; Required variable_info attributes
;
; Optional variable_info attributes:
;
; Modification history
;    20200321-righi_mattia: cleaned code
;    20180320-lauer-axel: added tags for reporting and call to write_references
;    20170315-cionni_irene: written
;
; #############################################################################
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

  ; Check for reference dataset definition
  if (variable_info[0]@reference_dataset.eq."None") then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset is specified")
  end if
  iref = ind(names.eq.variable_info[0]@reference_dataset)

  ; Check variable attributes
  if(isatt(variable_info, "units")) then
    units = variable_info@units
  else
    units = "m/s"
  end if

  if(isatt(variable_info, "long_name")) then
    long_name = variable_info@long_name
  else
    long_name = var0
  end if

  ; Set default values for optional diag_script_info attributes
  set_default_att(diag_script_info, "season", "ONDJ")
  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "levmin", 1.)
  set_default_att(diag_script_info, "levmax", 1000.)
  set_default_att(diag_script_info, "latrange", "60S")
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "start_year", max(y1))
  set_default_att(diag_script_info, "end_year", min(y2))

  seasons = diag_script_info@season
  latrange = diag_script_info@latrange
  ind_month = calc_season_index(seasons)
  ind_start = ind_month(0)
  ind_end = ind_month(dimsizes(ind_month)-1)
  if (ind_end.le.ind_start) then
    ind_end = ind_end + 12
  end if

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

  Max_lev = new((/dimsizes(seasons)/), "float")
  Min_lev = new((/dimsizes(seasons)/), "float")
  XMax = new((/dimsizes(seasons)/), "float")
  XMin = new((/dimsizes(seasons)/), "float")
  Max_lev = diag_script_info@levmax
  Min_lev = diag_script_info@levmin

  multimean = diag_script_info@multimean
  if (multimean)
    dimTot = dim_MOD + dim_OBS
  else
    dimTot = dim_MOD + dim_OBS - 1
  end if

  year_start = toint(diag_script_info@start_year)
  year_end = toint(diag_script_info@end_year)

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

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f02_" + var0)
  output_filename = wks@fullname
  nPanel = dimsizes(seasons)
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = acos(-1.) / 180.0
  A0 = read_data(info0[iref])
  lev0 = A0&plev
  dimLev = dimsizes(lev0)
  delete(A0)

  tmpData = new((/dim_MOD + dim_OBS, nPanel, 12, dimsizes(lev0), 2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "month"
  tmpData!3 = "plev"
  tmpData!4 = "stat"

  XXdata = new((/dimTot, nPanel, 12, 2/), "double")
  XXdata!0 = "line"
  XXdata!1 = "season"
  XXdata!2 = "month"
  XXdata!3 = "stat"

  XXref = new((/nPanel, 12, 3/), "double")
  XXref!0 = "season"
  XXref!1 = "month"
  XXdata!2 = "stat"

  Xdata = new((/dimTot, nPanel, 730, 2/), "double")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "day"
  Xdata!3 = "stat"

  Xref = new((/nPanel, 730, 3/), "double")
  Xref!0 = "season"
  Xref!1 = "day"
  Xref!2 = "stat"

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
      Aref = read_data(info0[iref])
      dummy1 = time_operations(A0, year_start, year_end, "average", \
                               "monthlyclim", True)
      dummy2_b = time_operations(A0, year_start, year_end, "extract", \
                                 "monthly", True)
      date = yyyymm_time(year_start, year_end, "integer")
      yyyy = date / 100
      mm = date - yyyy * 100
      weights = days_in_month(yyyy, mm)
      dummy2 = dummy2_b(0:11, :)  ; Copy metadata
      do im = 0, 11
        do ilev = 0, dimsizes(lev0) - 1
          dummy2(im, ilev) = \
            dim_stddev_wgt_Wrap(dummy2_b(im::12, ilev), weights(im::12), 1)
        end do
      end do
      delete(weights)
      delete(date)
      delete(yyyy)
      delete(mm)
      delete(dummy2_b)
      dummy2!0 = dummy1!0
      dummy2!1 = dummy1!1
      tmpData(imod, ipan, :, :, 0) = (/dummy1/)
      tmpData(imod, ipan, :, :, 1) = (/dummy2/)
      delete(dummy1)
      delete(dummy2)
    end do
  end do

  do imonth = 0, 11
    do ilat = 0, nPanel - 1
      ikk = 0
      do imod = 0, dim_MOD1 - 1
        id0 = tmpData(imod, ilat, imonth, :, 0)
        id1 = tmpData(imod, ilat, imonth, :, 1)
        if (imod.ne.iref) then
          if (.not.ismissing(any(ind(id0.ge.0.)))) then
            index_lev = ind(id0.ge.0.)
            XXdata(ikk, ilat, imonth, 0) = (/min(lev0(index_lev))/)
            delete(index_lev)
            plotMod(ikk) = names(imod)
            plotColor(ikk) = colors(imod)
            plotDash(ikk) = dashes(imod)
            plotThicks(ikk) = thicks(imod)
            plotScatters(ikk) = 0.
            LineScatters(ikk) = "Lines"
          end if
          ikk = ikk + 1
        else
          if (.not.ismissing(any(ind(id0.ge.0.)))) then
            index_lev = ind(id0.ge.0.)
            index_lev_p = ind((id0 + id1).ge.0.)
            index_lev_m = ind((id0 - id1).ge.0.)
            XXref(ilat, imonth, 0) = (/min(lev0(index_lev))/)
            XXref(ilat, imonth, 1) = (/min(lev0((index_lev_m)))/)
            XXref(ilat, imonth, 2) = (/min(lev0((index_lev_p)))/)
            delete(index_lev_p)
            delete(index_lev_m)
            delete(index_lev)
          end if
        end if
      end do
    end do
  end do

  ; Interpolate to daily data
  mid_month = (/15.5, 45, 74.5, 105, 135.5, 166, \
               196.5, 227.5, 258, 288.5, 319, 349.5/)
  XXdata!0 = "line"
  XXdata!1 = "season"
  XXdata!2 = "month"
  XXdata!3 = "stat"

  a = linint1_n_Wrap(mid_month, XXdata, False, ispan(1, 365, 1), 0, 2)
  a!2 = "day"
  Xdata(:, :, :364, :) = (/a(line|:, season|:, day|:, stat|:)/)
  Xdata(:, :, 365:, :) = (/a(line|:, season|:, day|:, stat|:)/)
  delete(a)
  XXref!0 = "season"
  XXref!1 = "month"
  XXref!2 = "stat"

  b = linint1_n_Wrap(mid_month, XXref, False, ispan(1, 365, 1), 0, 1)
  b!1 = "day"
  Xref(:, :364, :) = (/b(season|:, day|:, stat|:)/)
  Xref(:, 365:, :) = (/b(season|:, day|:, stat|:)/)
  delete(b)

  if (multimean) then
    Xdata(dimTot - 1, :, :, 0) = (/dim_avg_n(Xdata(index_mod, :, :, 0), 0)/)
    plotMod(dimTot - 1) = "MMM"
    plotColor(dimTot - 1) = "Black"
    plotDash(dimTot - 1) = 0
    LineScatters(dimTot - 1) = "Lines"
  end if

  X_MMM_std = (/dim_stddev_n(Xdata(index_mod, :, :, :), 0)/)
  YData = ispan(1, 730, 1)
  YData!0 = "day"
  YData&day = YData
  plot_xp_M_d = new((nPanel), "graphic")
  xp_M = new((/2 * dimsizes(YData), nPanel/), typeof(Xdata))
  yp_M = new((/2 * dimsizes(YData), nPanel/), "double")

  if (multimean.eq."True") then
    do j = 0, dimsizes(YData) - 1
      xp_M(j, :) = Xdata(dimTot - 1, :, j, 0) - ((X_MMM_std(:, j, 0)))
      yp_M(j, :) = YData(j)
      yp_M(dimsizes(YData) + j, :) = YData(dimsizes(YData) - 1 - j)
      xp_M(dimsizes(YData) + j, :) = \
        Xdata(dimTot - 1, :, (dimsizes(YData) - 1 - j), 0) + \
        X_MMM_std(:, (dimsizes(YData) - 1 - j), 0)
    end do
  end if
  xp_R = new((/2 * dimsizes(YData), nPanel/), typeof(Xdata))
  yp_R = new((/2 * dimsizes(YData), nPanel/), "double")

  if (isatt(variable_info[0],  "reference_dataset")) then
    do j = 0, dimsizes(YData) - 1
      xp_R(j, :) = Xref(:, j, 1) / 100.
      yp_R(j, :) = YData(j)
      yp_R(dimsizes(YData) + j, :) = YData(dimsizes(YData) - 1 - j)
      xp_R(dimsizes(YData) + j, :) = \
        Xref(:, (dimsizes(YData) - 1 - j), 2) / 100.
    end do
  end if

  plot_d = new(dimsizes(seasons), "graphic")
  plot_xp_R_d = new(dimsizes(seasons), "graphic")
  plot_xp_M_d = new(dimsizes(seasons), "graphic")
  plot_R = new(dimsizes(seasons), "graphic")

  ; Plotting
  mid_month2 = array_append_record(mid_month, mid_month + 365, 0)
  month_names = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", \
                 "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"/)
  month_names := array_append_record(month_names, month_names, 0)

  Xdata = (/Xdata / 100./)
  Xref = (/Xref / 100./)

  do ipanel = 0, dimsizes(seasons) - 1
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tfPolyDrawOrder = "Predraw"
    res@trXMinF = mid_month2(ind_start) - 15
    res@trXMaxF = mid_month2(ind_end) + 15
    res@tmXBMode = "Explicit"
    res@tmXBValues = mid_month2
    res@tmXBLabels = month_names
    res@trYMaxF = Max_lev
    res@trYMinF = Min_lev
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/1, 3, 10, 30, 100, 300, 1000/)
    res@tmYLLabels = (/1, 3, 10, 30, 100, 300, 1000/)
    res@tmLabelAutoStride = True
    res@tiYAxisString = "Pressure (hPa)"
    res@tiXAxisString = "Date"
    res@xyMarkLineModes = LineScatters
    res@trYAxisType = "LogAxis"
    res@trYReverse = True
    res@xyLineThicknessF = 2.0
    res@xyLineDashSegLenF = 0.13
    res@xyLineColors = plotColor
    res@xyDashPatterns = plotDash
    res@xyMarkers = plotDash
    res@xyMarkerSizeF = 0.004
    res@xyMarkerThicknessF = 1.8
    res@tmXMajorGrid = True
    res@tmXMajorGridThicknessF = 1.
    res@tmXMajorGridLineDashPattern = 2

    Xdata!0 = "line"
    Xdata!1 = "season"
    Xdata!2 = "day"
    Xdata!3 = "stat"
    plot_d(ipanel) = gsn_csm_xy(wks, YData, \
                                Xdata(line|:, season|ipanel, day|:, stat|0), \
                                res)

    if (isatt(variable_info[0], "reference_dataset")) then
      gsres = True
      gsres@gsFillOpacityF = 0.5
      gsres@gsFillColor = "blue"
      resR = True

      resR@gsLineColor = colors(iref)
      resR@gsLineDashPattern = dashes(iref)
      resR@gsFontThicknessF = thicks(iref)
      plot_xp_R_d(ipanel) = gsn_add_polygon(wks, \
                                            plot_d(ipanel), \
                                            yp_R(:, ipanel), \
                                            xp_R(:, ipanel), \
                                            gsres)
    end if

    if (multimean) then
      gsres1 = True
      gsres1@gsFillColor = "cyan"
      gsres1@gsFillOpacityF = 0.2
      plot_xp_M_d(ipanel) = gsn_add_polygon(wks, \
                                            plot_d(ipanel), \
                                            yp_M(:, ipanel), \
                                            xp_M(:, ipanel), \
                                            gsres1)
    end if
  end do

  ; Call legend-plotting routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@gsnPanelFigureStrings = plot_num(0:dimsizes(seasons) - 1)
  Pres@gsnPanelFigureStringsFontHeightF = 0.025
  Pres@gsnPanelFigureStringsPerimOn = "False"
  Pres@amJust = "TopLeft"

  if (dimsizes(seasons).gt.1) then
    gsn_panel(wks, plot_d, (/2, dimsizes(seasons) / 2/), Pres)
  else
    gsn_panel(wks, plot_d, (/1, 1/), Pres)
  end if

  if (dimsizes(plotMod).gt.21) then
    frame(wks)
  end if
  Lin = new(dimsizes(names), "string")
  Lin = "Lines"
  add_legend(wks, names, colors, tofloat(dashes), Lin, thicks, "0.2")

  frame(wks)

  nc_filename = work_dir + "eyring06f02_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f02"
  Xdata&stat = (/"mean", "avg"/)
  Xdata&season = seasons
  Xdata&line = plotMod
  Xdata&day = YData
  nc_outfile = ncdf_write(Xdata, nc_filename)

end
