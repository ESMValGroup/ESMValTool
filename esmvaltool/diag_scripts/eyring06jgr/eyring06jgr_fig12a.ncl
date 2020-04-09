; #############################################################################
; SCRIPT TO REPRODUCE FIG.12A OF EYRING ET AL. (JGR, 2006)
; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
; PROJECT-NAME CCMI
; #############################################################################
; Description
;    Vertical profile climatological mean ,
;    climatological mean for selected month and latitudinal region
;
; Required diag_script_info attributes (diagnostics specific)
;    latmin:  min lat where variable is averaged
;    latmax:  max lat where variable is averaged
;    month:   month when variable is averaged
;                         (i.e. "3", "4", "9")
;    _multimean: calculate multi-model mean (i.e. False/True)
;
; Optional diag_script_info attributes (diagnostic specific)
;    XMax:     max limit X axis
;    XMin:     min limit X axis
;    levmax:   max limit Y axis in Km
;    levmin:   min limit Y axis in Km
;    start_year: year when start the climatology calculation
;    end_year: year when end the climatology calculation
; Required variable_info attributes (variable specific)
;
; Caveats
;
; Modification history
;    20200325-A_cion_ir: cleaned
;    20180320-A_laue_ax: added tags for reporting and call to write_references
;    20170315-A_cion_ir: written
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
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
    do imod = 0, dimsizes(names) - 1
      if (names(imod) .eq. refname) then
        iref = imod
      end if
    end do
  else
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset  " + \
              "(variable_info[0]@reference_dataset) needs to be defined.")
  end if
    ; Check variable attributes
  if(isatt(variable_info, "units")) then
    units = variable_info@units
  else
    units = "10e9"
  end if

  if(isatt(variable_info, "long_name")) then
    long_name = variable_info@long_name
  else
    long_name = var0
  end if
  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "month", 03)
  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "latmin", 75.)
  set_default_att(diag_script_info, "latmax", 85.)
  set_default_att(diag_script_info, "XMax", 10.)
  set_default_att(diag_script_info, "XMin", 0.)
  set_default_att(diag_script_info, "start_year", max(y1))
  set_default_att(diag_script_info, "end_year", min(y2))
  set_default_att(diag_script_info, "multimean", False)

  seasons = diag_script_info@month
  month_names = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG",\
                  "SEP", "OCT", "NOV", "DEC"/)
  year_start = toint(diag_script_info@start_year)
  year_end = toint(diag_script_info@end_year)
  latrange = new((/dimsizes(seasons), 2/), "float")
  latrange(:, 0) = diag_script_info@latmin
  latrange(:, 1) = diag_script_info@latmax
  seasons = diag_script_info@month

  ; Styles
  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

  flag_mod = where(projects.ne."OBS6" .and. projects.ne."OBS" \
                   .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  Max_lev = diag_script_info@levmax
  Min_lev = diag_script_info@levmin
  XMax = diag_script_info@XMax
  XMin = diag_script_info@XMin

  multimean = diag_script_info@multimean
  if (multimean)
    dimTot = dim_MOD + dim_OBS + 1
  else
    dimTot = dim_MOD + dim_OBS
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
  nPanel = dimsizes(seasons)
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  A0 = read_data(info0[iref])
  lev0 = A0&plev
  a = stdatmus_p2tdz(lev0 / 100) / 1000.
  zlev = (- 7 * log(lev0 / 100000))
  delete(lev0)
  lev0 = a(2, :)
  dimLev = dimsizes(lev0)
  delete(A0)
  delete(zlev)
  delete(a)

  tmpData = new((/dim_MOD + dim_OBS, nPanel, dimsizes(lev0), 2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "plev"
  tmpData!3 = "stat"

  Xdata = new((/dimTot, nPanel, dimsizes(lev0), 2/), "double")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "plev"
  Xdata!3 = "stat"
  Xdata&plev = lev0

  Xref = new((/nPanel, dimsizes(lev0), 2/), "double")
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

      coslat   = cos(A0&lat * rad)
      coslat!0 = "lat"
      coslat&lat = A0&lat

      dummy1 = dim_avg_wgt_Wrap(A0(:, :, {latrange(ipan, 0): \
                                   latrange(ipan, 1)}),\
                                coslat({latrange(ipan, 0): \
                                       latrange(ipan, 1)}), 1)
      delete(coslat)
      delete(A0)

      if (seasons(ipan) .ne. "ANN") then
        dummy2 = time_operations(dummy1, year_start, year_end, "extract",\
                                 seasons(ipan), True)
      else
        time = dummy1&time
        years1 = cd_calendar(time, 0)
        dummy2a = month_to_annual(dummy1(:, :), 1)
        dummy2a&year = years1(:: 12, 0)
        dummy2 = dummy2a({year_start : year_end}, :)
        delete(time)
        delete(years1)
        delete(dummy2a)
      end if
      delete(dummy1)
      tmp = dim_avg_n_Wrap(dummy2, 0)
      tmp_st = dim_stddev_n_Wrap(dummy2, 0)
      delete(dummy2)
      tmpData(imod, ipan, :, 0) = (/tmp/)
      tmpData(imod, ipan, :, 1) = (/tmp_st/)
      delete(tmp)
      delete(tmp_st)
    end do  ; end seasons loop

    Xdata(imod, :, :, :) = (/tmpData(imod, :, :, :)/)
    plotMod(imod) = names(imod)
    plotColor(imod) = colors(imod)
    plotDash(imod) = dashes(imod)
    plotThicks(imod) = thicks(imod)
    plotScatters(imod) = 0.
    LineScatters(imod) = "Lines"

    if (imod .eq. iref) then
      Xref = tmpData(imod, :, :, :)
      plotDash(imod) = 16.
      plotScatters(imod) = 16.
      LineScatters(imod) = "Markers"
    end if

  end do
  ; Calculate multimodel mean
  if (multimean) then
    Xdata(dimTot - 1, :, :, :) = (/dim_avg_n(Xdata(index_mod, :, :, :), 0)/)
    plotMod(dimTot - 1) = "MMM"
    plotColor(dimTot - 1) = "Black"
    plotDash(dimTot - 1) = 0.
    plotThicks(dimTot - 1) = 4.
    LineScatters(dimTot - 1) = "Lines"
  end if
  X_MMM_std = (/dim_stddev_n(Xdata(index_mod, :, :, :), 0)/)
; UNITS CHANGE TO PPMV it would be changed in the preprocessor
  if ((units .eq. "10e6") .or. (units .eq. "ppmv")) then
    Xdata = (/ Xdata * 10 ^ 3/)
    tmpData = (/ tmpData * 10 ^ 3/)
    Xref = (/ Xref * 10 ^ 3/)
    units = "ppbv"
  end if
  Xdata&plev = lev0
  tmpData&plev = lev0

  Ydata = lev0

  plot_d = new(dimsizes(seasons), "graphic")

  plot_xp_R_d = new(dimsizes(seasons), "graphic")
  plot_xp_M_d = new(dimsizes(seasons), "graphic")

  xp_R = new((/2 * dimsizes(lev0), dimsizes(seasons)/), typeof(Xref))
  yp_R = new((/2 * dimsizes(lev0), dimsizes(seasons)/), typeof(Xref))

  xp_M = new((/2 * dimsizes(lev0), dimsizes(seasons)/), typeof(Xref))
  yp_M = new((/2 * dimsizes(lev0), dimsizes(seasons)/), typeof(Xref))

  do j = 0, dimsizes(lev0) - 1
    xp_R(j, :) = (Xref(:, j, 0) - Xref(:, j, 1))
    yp_R(j, :) = tofloat(lev0(j) / 100.)
    yp_R(dimsizes(lev0) + j, :) = tofloat(lev0(dimsizes(lev0) - 1 - j) \
                                          / 100)
    xp_R(dimsizes(lev0) + j, :) = (Xref(:, (dimsizes(lev0) - 1 - j), 0) \
                                   + Xref(:, (dimsizes(lev0) - 1 - j), 1))
  end do
  if (multimean) then
    do j = 0, dimsizes(lev0) - 1
      xp_M(j, :) = Xdata(dimTot - 1, :, j, 0) - X_MMM_std(:, j, 0)
      yp_M(j, :) = tofloat(Ydata(j))
      yp_M(dimsizes(lev0) + j, :) = tofloat(Ydata(dimsizes(Ydata) - 1 - j))
      xp_M(dimsizes(lev0) + j, :) = Xdata(dimTot-1, :,\
                                          (dimsizes(lev0) - 1 - j), 0) + \
        X_MMM_std(:, (dimsizes(lev0) - 1 - j), 0)
    end do
  end if
  ; Plotting
  plotpath = config_user_info@plot_dir + "eyring06f12a_" + var0
  wks = gsn_open_wks(file_type, plotpath)
  plotpath1 = config_user_info@plot_dir + "eyring06f12a_legend" + var0
  wks1 = gsn_open_wks(file_type, plotpath1)

  do ipanel = 0, dimsizes(seasons) - 1
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tfPolyDrawOrder = "Predraw"
    res@trXMaxF = XMax(ipanel)
    res@trXMinF = XMin(ipanel)
    res@trYMaxF = Max_lev(ipanel)
    res@trYMinF = Min_lev(ipanel)
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/10, 20, 30, 40, 50/)
    res@tmYLLabels = (/10, 20, 30, 40, 50/)
    res@tmLabelAutoStride = True
    res@vpHeightF = 0.6
    res@vpWidthF = 0.45
    res@tiYAxisString = "Z (Km)"
    res@tiXAxisString = long_name + " " + "(" + units + ")"
    if (seasons(ipanel) .ne. "ANN") then
      seasons_name = month_names(toint(seasons) - 1)
    else
      seasons_name = "ANN"
    end if
    res@tiMainString = long_name + " " + lat_names(latrange(ipanel, 0)) + \
      " - " + lat_names(latrange(ipanel, 1)) + \
      " " + seasons_name

    res@xyMarkLineModes = LineScatters
    res@xyLineThicknessF = 4.0
    res@xyLineDashSegLenF = 0.13
    res@xyLineColors = plotColor(:)
    res@xyDashPatterns = plotDash(:)
    res@xyMarkers = plotDash(:)
    res@xyMarkerSizeF = 0.01
    res@xyMarkerThicknessF = 1.8

    Xdata&plev = lev0
    Ydata!0 = "plev"
    Ydata&plev = lev0

    plot_d(ipanel) = gsn_csm_xy(wks, Xdata(line|:, season|ipanel, plev|:, \
                                stat|0), Ydata(plev|:), res)

    gsres = True
    gsres@gsFillColor = "Snow3"
    if (isatt(variable_info[0], "reference_dataset")) then
      plot_xp_R_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel), \
                                            xp_R(:, ipanel), \
                                            yp_R(:, ipanel), \
                                            gsres)
    end if
    if (multimean) then
      gsres1 = True
      gsres1@gsFillColor = "cyan"
      gsres1@gsFillOpacityF = 0.2
      plot_xp_M_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel),\
                                            xp_M(:, ipanel), \
                                            yp_M(:, ipanel), \
                                            gsres1)
    end if
  end do

; Call legend-plotting routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@gsnPanelFigureStrings = plot_num(0 : dimsizes(seasons) - 1)
  Pres@gsnPanelFigureStringsFontHeightF = 0.025
  Pres@gsnPanelFigureStringsPerimOn = "False"
  Pres@amJust = "TopLeft"

  ; Call panel routine
  gsn_panel(wks, plot_d, (/1, dimsizes(seasons)/), Pres)

  frame(wks)

  add_legend(wks1, plotMod, plotColor, tofloat(plotDash), LineScatters, \
             plotThicks + 2., "0.9")
  frame(wks1)
  ; Common output file
  nc_filename = config_user_info@work_dir + "eyring06f12a_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f12a"
  Xdata&stat = (/"mean", "stdev"/)
  Xdata&season = seasons
  Xdata&line = plotMod
  Xdata&plev = lev0
  nc_outfile = ncdf_write(Xdata, nc_filename)
  ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Seasonal vertical profiles of " + var0, \
                 (/"mean", "stddev"/), \
                 (/"shpolar"/), \
                 "vert", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
