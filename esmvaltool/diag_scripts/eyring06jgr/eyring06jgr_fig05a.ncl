; ############################################################################
; SCRIPT TO REPRODUCE FIG. 5a OF EYRING ET AL. (JGR, 2006)
; Authors: Irene Cionni (ENEA, Italy) and Veronika Eyring (DLR, Germany)
; CCMI PROJECT
; ############################################################################

; Description
;    Vertical profile climatological mean ,
;    climatological mean for selected month and latitudinal region
;
; Required diag_script_info attributes (diagnostics specific)
;    latrange:  min lat and max lat where variable is averaged
;               (i.e. (/60, 90/))
;    season:   month when variable is averaged
;              (i.e. "03", "04", "09")
;    multimean: calculate multi-model mean (i.e. False/True)
;    refModel:  name of the reference model or observation
;               (e.g. (/"ERA-Interim"/))
;
; Optional diag_script_info attributes (diagnostic specific)
;    XMax:     max limit X axis
;    XMin:     min limit X axis
;    levmax:   max limit Y axis
;    levmin:   min limit Y axis
;    start_year: year when start the climatology calculation
;    start_year: year when end the climatology calculation
;
; Required variable_info attributes (variable specific)
;
; Caveats
;    units for the variable are now handled in the pre-processor which
;    means that the preprocessor must be adjusted if the variable is
;    changed in the recipe from ozone to something else
;
; Modification history
;    20200403-hassler_birgit: cleaned code
;    20180320-lauer_axel: added tags for reporting and call to write_references
;    20170315-cionni_irene: written
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

  if (isatt(variable_info, "long_name")) then
    LONG_NAME = variable_info@long_name
  else
    LONG_NAME = var0
  end if

  ; Set default values for optional diag_script_info attributes
  set_default_att(diag_script_info, "month", 03)
  set_default_att(diag_script_info, "latmin", 75.)
  set_default_att(diag_script_info, "latmax", 85.)
  set_default_att(diag_script_info, "Xmax", 10.)
  set_default_att(diag_script_info, "Xmin", 0.)
  set_default_att(diag_script_info, "start_year", max(y1))
  set_default_att(diag_script_info, "end_year", min(y2))
  set_default_att(diag_script_info, "multimean", False)

  seasons = diag_script_info@month
  month_names = \
    (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", \
      "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"/)
  seasons_name = month_names(toint(seasons) - 1)

  latrange = new((/dimsizes(seasons), 2/), "float")
  latrange(:, 0) = diag_script_info@latmin
  latrange(:, 1) = diag_script_info@latmax
  seasons = diag_script_info@month
  flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  ; Styles
  styleset = diag_script_info@styleset
  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  Max_lev = diag_script_info@levmax
  Min_lev = diag_script_info@levmin
  XMax = diag_script_info@Xmax
  XMin = diag_script_info@Xmin

  multimean = diag_script_info@multimean
  if (multimean) then
    dimTot = dim_MOD + dim_OBS + 1
  else
    dimTot = dim_MOD + dim_OBS
  end if

  year_start = toint(diag_script_info@start_year)
  year_end = toint(diag_script_info@end_year)

  ; make sure path for (mandatory) netcdf output exists
  work_dir = config_user_info@work_dir + "/"

  ; Create work dir
  system("mkdir -p " + work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

end

begin

  aux_title_info = ""
  plotpath = config_user_info@plot_dir + "eyring06f05a_" + var0
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f05a_" + var0)
  wks = gsn_open_wks(file_type, plotpath)
  plotpath1 = config_user_info@plot_dir + "eyring06f05a_legend_" + var0
  wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f05a_legend_" + var0)
  wks1 = gsn_open_wks(file_type, plotpath1)
  nPanel = dimsizes(seasons)
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  A0 = read_data(info0[iref])
  lev0 = A0&plev
  dimLev = dimsizes(lev0)
  delete(A0)

  tmpData = new((/dim_MOD + dim_OBS, nPanel, dimsizes(lev0), 2/), "float")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "plev"
  tmpData!3 = "stat"

  Xdata = new((/dimTot, nPanel, dimsizes(lev0), 2/), "float")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "plev"
  Xdata!3 = "stat"
  Xdata&plev = lev0

  Xref = new((/ nPanel, dimsizes(lev0), 2/), "float")
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

      printVarSummary(A0)
      coslat   = cos(A0&lat * rad)
      coslat!0 = "lat"
      coslat&lat = A0&lat
      dummy1 = \
        dim_avg_wgt_Wrap(A0(:, :, {latrange(ipan, 0):latrange(ipan, 1)}), \
                         coslat({latrange(ipan, 0):latrange(ipan, 1)}), 1)
      delete(coslat)
      delete(A0)
      unit_string = dummy1@units

      if (seasons(ipan).ne."ANN") then
        dummy2 = time_operations(dummy1, year_start, year_end, \
                                 "extract", seasons(ipan), True)
      else
        time = dummy1&time
        years1 = cd_calendar(time, 0)
        dummy2a = month_to_annual(dummy1(:, :), 1)
        dummy2a&year = years1(::12, 0)
        dummy2 = dummy2a({year_start:year_end}, :)
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

    if (imod.eq.iref) then
      Xref = tmpData(imod, :, :, :)
      plotDash(imod) = 16.
      plotScatters(imod) = 16.
      LineScatters(imod) = "Markers"
    end if

  end do  ; end model loop

  if (multimean) then
    Xdata(dimTot - 1, :, :, :) = \
      (/dim_avg_n(Xdata(index_mod, :, :, :), 0)/)
    plotMod(dimTot - 1) = "MMM"
    plotColor(dimTot - 1) = "Black"
    plotDash(dimTot - 1) = 0.
    plotThicks(dimTot - 1) = 4.
    LineScatters(dimTot - 1) = "Lines"
  end if

  ; Calculate multimodel mean
  X_MMM_std = (/dim_stddev_n(Xdata(index_mod, :, :, :), 0)/)

  Ydata = lev0 / 100.

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
    yp_R(dimsizes(lev0) + j, :) = \
      tofloat(lev0(dimsizes(lev0) - 1 - j) / 100)
    xp_R(dimsizes(lev0) + j, :) = \
      Xref(:, (dimsizes(lev0) - 1 - j), 0) + \
      Xref(:, (dimsizes(lev0) - 1 - j), 1)
  end do

  if (multimean) then
    do j = 0, dimsizes(lev0) - 1
      xp_M(j, :) = Xdata(dimTot - 1, :, j, 0) - X_MMM_std(:, j, 0)
      yp_M(j, :) = tofloat(lev0(j) / 100.)
      yp_M(dimsizes(lev0) + j, :) = \
        tofloat(lev0(dimsizes(lev0) - 1 - j) / 100.)
      xp_M(dimsizes(lev0) + j, :) = \
        Xdata(dimTot - 1, :, dimsizes(lev0) - 1 - j, 0) + \
        X_MMM_std(:, dimsizes(lev0) - 1 - j, 0)
    end do
  end if

  ; Plotting
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
    res@tmYLValues = (/1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 200/)
    res@tmYLLabels = \
      (/"1", "2", "3", "5", "7", "10", \
        "20", "30", "50", "70", "100", "200"/)
    res@tmLabelAutoStride = True

    res@vpHeightF = 0.6                   ; change aspect ratio of plot
    res@vpWidthF = 0.45

    ; y-axis label
    res@tiYAxisString = "Pressure (hPa)"
    ; x-axis label
    res@tiXAxisString = LONG_NAME + " " + "(" + unit_string + ")"

    res@tiMainString = \
      LONG_NAME + " " + \
      lat_names(latrange(ipanel, 0)) + " - " + \
      lat_names(latrange(ipanel, 1)) + " " + \
      seasons_name(ipanel)

    res@xyMarkLineModes = LineScatters

    res@trYAxisType = "LogAxis"
    res@trYReverse = True
    res@xyLineThicknessF = 4.0
    res@xyLineDashSegLenF = 0.13
    res@xyLineColors = plotColor                 ; line color
    res@xyDashPatterns = plotDash
    res@xyMarkers = plotDash
    res@xyMarkerSizeF = 0.01
    res@xyMarkerThicknessF = 1.8

    Xdata&plev = lev0
    Ydata!0 = "plev"
    Ydata&plev = lev0
    plot_d(ipanel) = \
      gsn_csm_xy(wks, \
                 Xdata(line|:, season|ipanel, plev|:, stat|0), \
                 Ydata(plev|:), \
                 res)

    gsres = True

    gsres@gsFillColor = "Snow3"
    plot_xp_R_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel), \
                                          xp_R(:, ipanel), \
                                          yp_R(:, ipanel), gsres)
    if (multimean) then
      gsres1 = True
      gsres1@gsFillColor = "cyan"
      gsres1@gsFillOpacityF = 0.2
      plot_xp_M_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel), \
                                            xp_M(:, ipanel), \
                                            yp_M(:, ipanel), gsres1)
    end if
  end do

  ; Call legend - plotting routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@gsnPanelFigureStrings = plot_num(0:dimsizes(seasons) - 1)
  Pres@gsnPanelFigureStringsFontHeightF = 0.025
  Pres@gsnPanelFigureStringsPerimOn = "False"
  Pres@amJust = "TopLeft"

  gsn_panel(wks, plot_d, (/1, dimsizes(seasons)/), Pres)
  frame(wks)

  add_legend(wks1, plotMod, plotColor, tofloat(plotDash), \
             LineScatters, plotThicks + 2., "0.9")
  frame(wks1)

  nc_filename = work_dir + "eyring06f05a_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f05a"
  printVarSummary(Xdata)
  Xdata&stat = (/"mean", "avg"/)
  Xdata&season = seasons
  Xdata&line = plotMod
  Xdata&plev = lev0
  printVarSummary(Xdata)
  nc_outfile = ncdf_write(Xdata, nc_filename)

  ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Vertical profile climatological means of " + var0, \
                 (/"mean"/), \
                 (/"nhpolar", "shpolar", "trop"/), \
                 "vert", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))


end
