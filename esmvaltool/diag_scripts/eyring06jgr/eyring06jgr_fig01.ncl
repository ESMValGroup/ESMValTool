; #############################################################################
; SCRIPT TO REPRODUCE FIG. 1 OF EYRING ET AL. (JGR, 2006)
; Authors: Irene Cionni (ENEA, Italy) and Veronika Eyring (DLR, Germany)
; CCMI Project
; #############################################################################
;
; Description:
;    Vertical profile climatological mean bias and climatological mean for
;    selected seasons and latitudinal regions.
;
; Required diag_script_info attributes:
;    latrange:  minimum and maximum latitude for regional average,
;               e.g. (/60, 90/))
;    season: seasonfor temporal average (e.g, "DJF", "MAM", "JJA", "SON")
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
;
; Required variable_info attributes:
;
; Optional variable_info attributes:
;
; Modification history
;    20200320-righi_mattia: cleaned code
;    20180320-lauer_axel: added tags for reporting and call to write_references
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
  if (isatt(variable_info, "units")) then
    units = variable_info@units
  else
    units = "K"
  end if

  if (isatt(variable_info, "long_name")) then
    long_name = variable_info@long_name
  else
    long_name = var0
  end if

  ; Set default values for optional diag_script_info attributes
  set_default_att(diag_script_info, "season", "ANN")
  set_default_att(diag_script_info, "latmin", -90.)
  set_default_att(diag_script_info, "latmax", 90.)
  set_default_att(diag_script_info, "Xmin", -30.)
  set_default_att(diag_script_info, "Xmax", 20.)
  set_default_att(diag_script_info, "Ymin", 1.)
  set_default_att(diag_script_info, "Ymax", 1000.)
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "start_year", max(y1))
  set_default_att(diag_script_info, "end_year", min(y2))

  seasons = diag_script_info@season
  styleset = diag_script_info@styleset
  nseasons = dimsizes(seasons)
  latrange = new((/nseasons, 2/), "float")
  latrange(:, 0) = diag_script_info@latmin
  latrange(:, 1) = diag_script_info@latmax
  flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  multimean = diag_script_info@multimean
  if (multimean) then
    dimTot = dim_MOD + dim_OBS
  else
    dimTot = dim_MOD + dim_OBS - 1
  end if

  year_start = toint(diag_script_info@start_year)
  year_end = toint(diag_script_info@end_year)

  ; Styles
  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

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

  nPanel = nseasons
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = acos(-1.) / 180.0
  A0 = read_data(info0[iref])
  lev0 = A0&plev
  nlev0 = dimsizes(lev0)
  delete(A0)

  tmpData = new((/dim_MOD + dim_OBS, nPanel, nlev0, 2/), "float")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "plev"
  tmpData!3 = "stat"

  Xdata = new((/dimTot, nPanel, nlev0, 2/), "float")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "plev"
  Xdata!3 = "stat"
  Xdata&plev = lev0

  Xref = new((/nPanel, nlev0, 2/), "float")
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
    do ipan = 0, nseasons - 1

      ; Extract data
      A0 = read_data(info0[imod])

      ; Select latitude range
      coslat = cos(A0&lat * rad)
      coslat!0 = "lat"
      coslat&lat = A0&lat
      A0 := \
        dim_avg_wgt_Wrap(A0(:, :, {latrange(ipan, 0):latrange(ipan, 1)}), \
                         coslat({latrange(ipan, 0):latrange(ipan, 1)}), 1)
      delete(coslat)

      if (seasons(ipan).ne."ANN") then
        tmp = time_operations(A0, year_start, year_end, \
                              "extract", seasons(ipan), True)
        dummy2 = time_operations(tmp, year_start, year_end, \
                                 "average", "yearly", True)
        delete(tmp)
      else
        dummy2 = time_operations(A0, year_start, year_end, \
                                 "average", "yearly", True)
      end if
      delete(A0)
      tmp = dim_avg_n_Wrap(dummy2, 0)
      tmp_st = dim_stddev_n_Wrap(dummy2, 0)
      delete(dummy2)
      tmpData(imod, ipan, :, 0) = (/tmp/)
      tmpData(imod, ipan, :, 1) = (/tmp_st/)
      delete(tmp)
      delete(tmp_st)
    end do  ; end seasons loop

    if (imod.ne.iref) then
      Xdata(ik, :, :, :) = (/tmpData(imod, :, :, :)/)
      plotMod(ik) = names(imod)
      plotColor(ik) = colors(imod)
      plotDash(ik) = dashes(imod)
      plotThicks(ik) = thicks(imod)
      plotScatters(ik) = 0.
      LineScatters(ik) = "Lines"
      ik = ik + 1
    else
      Xref = tmpData(imod, :, :, :)
    end if

  end do

  ; Calculate multimodel mean
  if (multimean) then
    Xdata(dimTot - 1, :, :, :) = (/dim_avg_n(Xdata(index_mod, :, :, :), 0)/)
    plotMod(dimTot - 1) = "MMM"
    plotColor(dimTot - 1) = "black"
    plotDash(dimTot - 1) = 0
    LineScatters(dimTot - 1) = "Lines"
    plotThicks(dimTot - 1) = 4.
  end if

  X_MMM_std = (/dim_stddev_n(Xdata(index_mod, :, :, :), 0)/)

  Xdata = (/Xdata - conform_dims(dimsizes(Xdata), Xref, (/1, 2, 3/))/)
  Ydata = lev0 / 100.

  plot_d = new(nseasons, "graphic")
  plot_xp_R_d = new(nseasons, "graphic")
  plot_xp_M_d = new(nseasons, "graphic")
  xp_R = new((/2 * nlev0, nseasons/), typeof(Xref))
  yp_R = new((/2 * nlev0, nseasons/), typeof(Xref))
  xp_M = new((/2 * nlev0, nseasons/), typeof(Xref))
  yp_M = new((/2 * nlev0, nseasons/), typeof(Xref))
  do j = 0, nlev0 - 1
    xp_R(j, :) = -1. * Xref(:, j, 1)
    yp_R(j, :) = tofloat(lev0(j) / 100.)
    yp_R(nlev0 + j, :) = tofloat(lev0(nlev0 - 1 - j) / 100)
    xp_R(nlev0 + j, :) = (Xref(:, (nlev0 - 1 - j), 1))
  end do

  if (multimean) then
    do j = 0, nlev0 - 1
      xp_M(j, :) = Xdata(dimTot - 1, :, j, 0) - X_MMM_std(:, j, 0)
      yp_M(j, :) = tofloat(lev0(j) / 100.)
      yp_M(nlev0 + j, :) = tofloat(lev0(nlev0 - 1 - j) / 100.)
      xp_M(nlev0 + j, :) = Xdata(dimTot - 1, :, (nlev0 - 1 - j), 0) + \
        X_MMM_std(:, (nlev0 - 1 - j), 0)
    end do
  end if

  ; Plotting
  plotpath = config_user_info@plot_dir + "eyring06f01_" + var0
  wks = gsn_open_wks(file_type, plotpath)
  plotpath1 = config_user_info@plot_dir + "eyring06f01_legend_" + var0
  wks1 = gsn_open_wks(file_type, plotpath1)

  do ipanel = 0, nseasons - 1
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tfPolyDrawOrder = "Predraw"
    res@trXMinF = diag_script_info@Xmin(ipanel)
    res@trXMaxF = diag_script_info@Xmax(ipanel)
    res@trYMinF = diag_script_info@Ymin(ipanel)
    res@trYMaxF = diag_script_info@Ymax(ipanel)
    res@tmYLMode = "Explicit"
    res@tmYLValues = (/1, 3, 10, 30, 100, 300, 1000/)
    res@tmYLLabels = (/1, 3, 10, 30, 100, 300, 1000/)
    res@tmLabelAutoStride = True

    res@tiXAxisString = "Model " + long_name + " bias " + "(" + units + ")"
    res@tiYAxisString = "Pressure (hPa)"
    res@tiMainString = \
      seasons(ipanel) + " " + \
      lat_names(latrange(ipanel, 0)) + " - " + \
      lat_names(latrange(ipanel, 1))
    res@xyMarkLineModes = LineScatters
    res@trYAxisType = "LogAxis"
    res@trYReverse = True
    res@xyLineThicknessF = 2.
    res@xyLineDashSegLenF = 0.13
    res@xyLineColors = plotColor
    res@xyDashPatterns = plotDash
    res@xyMarkers = plotDash
    res@xyMarkerSizeF = 0.004
    res@xyMarkerThicknessF = 1.8
    res@tmXMajorGrid = True
    res@tmXMajorGridThicknessF = 1.
    res@tmXMajorGridLineDashPattern = 2

    Xdata&plev = lev0
    Ydata!0 = "plev"
    Ydata&plev = lev0
    plot_d(ipanel) = \
      gsn_csm_xy(wks, \
                 Xdata(line|:, season|ipanel, {plev|80000:}, stat|0), \
                 Ydata({plev|80000:}), \
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
  Pres@gsnPanelFigureStrings = plot_num(0:nseasons - 1)
  Pres@gsnPanelFigureStringsFontHeightF = 0.025
  Pres@gsnPanelFigureStringsPerimOn = "False"
  Pres@amJust = "TopLeft"

  ; Call panel routine
  if (nseasons.gt.1) then
    gsn_panel(wks, plot_d, (/2, nseasons/2/), Pres)
  else
    gsn_panel(wks, plot_d, (/1, 1/), Pres)
  end if

  pos_leg = "0.9"
  add_legend(wks1, plotMod, plotColor, tofloat(plotDash), \
             LineScatters, plotThicks, pos_leg)
  frame(wks1)
  frame(wks)

  nc_filename = config_user_info@work_dir + "eyring06f01_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f01"
  Xdata&stat = (/"mean", "stddev"/)
  Xdata&season = seasons
  Xdata&line = plotMod
  nc_outfile = ncdf_write(Xdata,  nc_filename)

  ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Seasonal vertical profiles of " + var0 + " bias", \
                 (/"mean", "stddev"/), \
                 (/"nhpolar", "shpolar"/), \
                 "vert", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
