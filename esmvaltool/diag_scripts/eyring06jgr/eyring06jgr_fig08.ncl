; #############################################################################
; SCRIPT TO REPRODUCE FIG. 8 OF EYRING ET AL. (JGR, 2006)
; Authors: Irene Cionni (ENEA, Italy) and Veronika Eyring (DLR, Germany)
; CCMI Project
; #############################################################################
; Description
;     Time-height sections of water vapor mixing ratio shown as the
;     deviation (in parts per million by volume) from the time mean profile,
;     averaged between min latitude (10°S) and max latitude (10°)
;
; Required diag_script_info attributes (diagnostics specific)
;     @multimean: calculate multi-model mean (i.e. False/True)
;                   bias calculation (e.g. (/"ERA-Interim"/))
;
; Optional diag_script_info attributes (diagnostic specific)
;    cnMax: max contour
;    cnMmin: min contour
;    cnSpacing: contour interval
;    latrange: name of the lat range as in the "Equator"
;    start_year: year when start the climatology calculation
;    start_year: year when end the climatology calculation
;
; Required variable_info attributes (variable specific)
;
; Caveat
;
; Modification history
;    20200227-cionni_irene: written
;    20170325-cionni_irene: written
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
  y1 = (/ metadata_att_as_array(info0, "start_year")/)
  y2 = (/ metadata_att_as_array(info0, "end_year") /)

  ; Check for reference dataset definition
  iref = 0
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
    do imod = 0, dimsizes(names) - 1
      if (names(imod) .eq. refname) then
        iref = imod
      end if
    end do
  end if

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "cnMin", -1.)
  set_default_att(diag_script_info, "cnMax", 1.)
  set_default_att(diag_script_info, "cnSpacing", 0.2)
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "start_year", max(y1))
  set_default_att(diag_script_info, "end_year", min(y2))
  set_default_att(diag_script_info, "latrange", "10S-10N")

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

  cnMax = diag_script_info@cnMax
  cnMin = diag_script_info@cnMin
  cnSpacing =  diag_script_info@cnSpacing
  latrange = diag_script_info@latrange
  multimean = diag_script_info@multimean
  if (multimean)
    dimTot = dim_MOD + dim_OBS + 1
  else
    dimTot = dim_MOD + dim_OBS
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
  nPanel = dimTot
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0

  A0 = read_data(info0[iref])
  lev0 = A0&plev
  dimLev = dimsizes(lev0)
  delete(A0)

  tmpData = new((/dim_MOD + dim_OBS, 25, dimsizes(lev0), 2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "month"
  tmpData!2 = "lev"
  tmpData!3 = "stat"

  Xdata = new((/dimTot, 25, dimsizes(lev0), 2/), "double")
  Xdata!0 = "line"
  Xdata!1 = "month"
  Xdata!2 = "lev"
  Xdata!3 = "stat"

  plotMod = new(dimTot, typeof(names))

  ik = 0

  ; Start datasets loop
  do imod = 0, dim_MOD1 - 1
    ; Extract data
    A0 = read_data(info0[imod])
    if (var0.eq."hus") then
      A0 = (/A0/0.62/)
      UNITS = "ppmv"
      A0 = (/A0 * 10 ^ 6/)
    end if

    dummy1 = time_operations(A0, year_start, year_end, "average", \
                             "monthlyclim", True)
    dummy2_b = tofloat(A0)
    date = yyyymm_time(year_start, year_end, "integer")
    yyyy = date / 100
    mm = date-yyyy * 100
    weights = days_in_month(yyyy, mm)
    dummy2 = dummy2_b(0:11, :)
    do ilev = 0, dimsizes(lev0) - 1
      do im = 0, 11
        dummy2(im, ilev) = dim_stddev_Wrap(dummy2_b(im :: 12, ilev))
      end do
    end do

    dummy3 = dummy1
    dummy3 = (/dummy1 - conform(dummy1, dim_avg_n(dummy1, 0), 1)/)
    tmpData(imod, 0:11, :, 0) = (/dummy3/)
    tmpData(imod, 0:11, :, 1) = (/dummy2/)
    tmpData(imod, 12:23, :, 0) = (/dummy3/)
    tmpData(imod, 12:23, :, 1) = (/dummy2/)
    tmpData(imod, 24, :, 0) = (/dummy3(0, :)/)
    tmpData(imod, 24, :, 1) = (/dummy2(0, :)/)
    delete(A0)
    delete(weights)
    delete(date)
    delete(yyyy)
    delete(mm)
    delete(dummy2_b)
    delete(dummy1)
    delete(dummy2)
    delete(dummy3)
  end do   ; end loop Dataset

  Xdata(: dim_MOD1 - 1, :, :, :) = (/tmpData/)
  plotMod(: dim_MOD1 - 1) = names

  ; Calculate multimodel mean
  if (multimean) then
    Xdata(dimTot - 1, :, :, 0) = (/dim_avg_n(Xdata(index_mod, :, :, 0), 0)/)
    plotMod(dimTot - 1) = "MMM"
    X_MMM_std = (/dim_stddev_n(Xdata(index_mod, :, :, :), 0)/)
  end if

  XXData = ispan(1, 25, 1)
  XXData!0 = "month"
  XXData&month = XXData

  Ydata = lev0 / 100.
  Ydata!0 = "lev"
  Ydata&lev = lev0 / 100.
  Ydata@units = "hPa"

  ; Plotting
  plotpath = config_user_info@plot_dir + "eyring06f08_" + var0
  wks = gsn_open_wks(file_type, plotpath)
  if (dimTot.gt.21) then
    plotpath1 = config_user_info@plot_dir + "eyring06f08_2_" + var0
    wks1 = gsn_open_wks(file_type, plotpath1)
  end if

  plot_d = new(dimTot, "graphic")

  do imod = 0, dimTot - 1
    res = True
    res@tiMainString = plotMod(imod)
    res@gsnDraw = False
    res@gsnFrame = False
    res@lbLabelBarOn = False
    res@cnLineLabelsOn = False
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnLevelSpacingF = cnSpacing
    res@cnMinLevelValF = cnMin
    res@cnMaxLevelValF = cnMax
    res@lbLabelStrings = (/"-1", "-0.8", "-0.6", "-0.4", "-0.2", "0.", \
                           "0.2", "0.4", "0.6", "0.8", "1."/)
    res@cnFillOn = True
    res@cnInfoLabelOn = False
    res@gsnSpreadColors = True
    res@gsnSpreadColorStart = 10
    res@lbLabelAutoStride = True
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@gsnLeftStringFontHeightF = 0.015
    res@tmXBMode = "Explicit"
    res@gsnAddCyclic = True
    res@tmXBValues = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,\
                      16, 17, 18, 19, 20, 21, 22, 23, 24, 25/)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N",\
                       "D", "J", "F", "M", "A", "M", "J", "J", "A", "S", "O",\
                       "N", "D", "J"/)
    res@trYMaxF = 100
    res@trYMinF = 10
    res@trXMinF = 0
    res@trXMaxF = 25
    res@vpHeightF = 0.4
    res@vpWidthF = 0.75
    res@lbLabelFontHeightF = 0.011
    res@tiYAxisString = "Pressure [hPa]"
    Xdata!2 = "lev"
    Xdata&lev = Ydata
    plot_d(imod) = gsn_csm_pres_hgt(wks, Xdata(line|imod, lev|:, month|:,\
                                    stat|0), res)
    delete(res)
  end do

  ; Call panel routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@gsnPanelFigureStrings = plot_num(0 : dimTot - 1)
  Pres@gsnPanelFigureStringsFontHeightF = 0.015
  Pres@gsnPanelFigureStringsPerimOn = "False"
  Pres@amJust = "TopLeft"
  Pres@lbLabelAutoStride = True
  Pres@gsnPanelLabelBar = True
  if (dimsizes(plotMod).gt.21) then
    gsn_panel(wks, plot_d(0 : 20), (/7, 3/), Pres)
    gsn_panel(wks1, plot_d(21:), (/1 + (dimTot-21) / 3, 3/), Pres)
    frame(wks)
    frame(wks1)
  else
    gsn_panel(wks, plot_d, (/1 + dimTot / 3, 3/), Pres)
    frame(wks)
  end if
  ; Common output file
  nc_filename = config_user_info@work_dir + "eyring06f08_" + var0 + ".nc"
  Xdata@var = var0
  Xdata@diag_script = "eyring06f08"
  Xdata&stat = (/"mean", "stddev"/)
  Xdata&line = plotMod
  month_names = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", \
                  "D", "J", "F", "M", "A", "M", "J", "J", "A", "S", "O", \
                  "N", "D", "J"/)
  Xdata&month = month_names
  nc_outfile = ncdf_write(Xdata, nc_filename)
  ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Annual cycle vertical profiles of " + var0 + " anomalies", \
                 (/"mean", "stddev"/), \
                 (/"trop"/), \
                 "vert", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))
end
