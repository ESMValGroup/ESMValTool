; Based on EYRING13JGR_FIG12
; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany),
;         Lisa Bock (DLR, Germany), Birgit Hassler (DLR, Germany)
; #############################################################################
; Description
;  Calculation of the long-term mean and trend of the zonal wind
;
;
; Optional diag_script_info attributes (diagnostic specific)
;    @e13fig12_start_year: year when the climatology calculation starts
;    @e13fig12_end_year: year when the climatology calculation ends
;    @e13fig12_multimean: multimodel mean
;    @e13fig12_season: season
;    @e13fig12_exp_MMM: name of the experiments dor the MMM
; Caveats
;
; Modification history
; #############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD1 = ListCount(info0)
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")
  y1 = (/ metadata_att_as_array(info0, "start_year")/)
  y2 = (/ metadata_att_as_array(info0, "end_year") /)
  work_dir = config_user_info@work_dir + "/"
  ; Create work dir

  exps_name = where(ismissing(exps), "OBS", exps)
  system("mkdir -p " + work_dir)

  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
    do imod = 0, dimsizes(names) - 1
      if (names(imod) .eq. refname) then
        iref = imod
      end if
    end do
  end if

  if(isatt(variable_info, "units")) then
    UNITS = variable_info@units
  else
    UNITS = "m s~S~-1~N~"
  end if

  if(isatt(variable_info, "long_name")) then
    LONG_NAME = variable_info@long_name
  else
    LONG_NAME = var0
  end if

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "e13fig12_multimean", False)
  set_default_att(diag_script_info, "e13fig12_start_year", max(y1))
  set_default_att(diag_script_info, "e13fig12_end_year", min(y2))
  set_default_att(diag_script_info, "e13fig12_season", "ANN")
  flag_mod = where(projects.ne."OBS6" .and. projects.ne."OBS" .and.\
                   projects.ne."obs4mips", 1, 0)

  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  multimean = diag_script_info@e13fig12_multimean
  season = diag_script_info@e13fig12_season
  if (multimean) then
    exp_MMM = diag_script_info@e13fig12_exp_MMM
    flag_exp = new((/dimsizes(projects), dimsizes(exp_MMM)/), "integer")
    do id = 0, dimsizes(exp_MMM) - 1
      flag_exp(:, id) = where(projects.eq.exp_MMM(id), 0, 1)
    end do
    dimTot = dim_MOD + dim_OBS + 1
    if (dimsizes(exp_MMM).gt.1) then
      dimTot = dim_MOD + dim_OBS + dimsizes(exp_MMM)
    end if
  else
    dimTot = dim_MOD + dim_OBS
  end if

  ; make sure path for (mandatory) netcdf output exists
  year_start = toint(diag_script_info@e13fig12_start_year)
  year_end = toint(diag_script_info@e13fig12_end_year)

end

begin

  aux_title_info = ""
  ; wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f12_" + var0 +\
  ;              "_" + plotMod)
  A0 = read_data(info0[iref])
  lat0 = A0&lat
  plev0 = A0&plev
  delete(A0)

  tmpData = new((/dim_MOD + dim_OBS, dimsizes(plev0), dimsizes(lat0), 3/),\
                "double")
  tmpData!0 = "line"
  tmpData!1 = "plev"
  tmpData!2 = "lat"
  tmpData!3 = "stat"

  XXdata = new((/dimTot, dimsizes(plev0), dimsizes(lat0), 3/), "double")
  XXdata!0 = "line"
  XXdata!1 = "plev"
  XXdata!2 = "lat"
  XXdata!3 = "stat"

  XXref = new((/dimsizes(plev0), dimsizes(lat0), 3/), "double")
  XXref!0 = "plev"
  XXref!1 = "lat"
  XXref!2 = "stat"

  years1 = new(dimTot, "integer")
  years2 = new(dimTot, "integer")
  plot = new(dimTot, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))

  plotMod = new(dimTot, typeof(names))

  ik = 0

; ******************************************************************
  do imod = 0, dim_MOD1 - 1
    years1(imod) = y1(imod)
    years2(imod) = y2(imod)
    years = tofloat(ispan(years1(imod), years2(imod), 1))
    plotMod(imod) = names(imod)
    ; Extract data
    A0 = read_data(info0[imod])
    if (season.eq."ANN") then
      A1 = time_operations(A0, year_start, year_end, "average",\
                           "annualclim", True)
      A2 = time_operations(A0, -1, -1, "extract", "yearly", True)
      A2trend = regCoef(years, A2(plev|:, lat|:, time|:))
      tval_A2 = A2trend@tval
      b_A2 = 0.5
      df_A2 = A2trend@nptxy - 2
      prob_A2 = (/1 - betainc(df_A2 / (df_A2 + tval_A2 ^ 2),\
                 df_A2 / 2.0, b_A2)/)
    else
      A1 = time_operations(A0, -1, -1, "average", season, True)
      A2 = month_to_season(A0, season)
      A2trend = regCoef(years, A2(plev|:, lat|:, time|:))
      tval_A2 = A2trend@tval
      b_A2 = tval_A2
      b_A2 = (/0.5/)
      df_A2 = A2trend@nptxy - 2
      prob_A2 = A2(0, :, :)
      prob_A2 = onedtond((/1 - betainc(df_A2 / (df_A2 + tval_A2 ^ 2), \
                          df_A2 / 2.0, b_A2)/), dimsizes(A2(0, :, :)))
      delete(df_A2)
      delete(tval_A2)
      delete(b_A2)
    end if
    delete(A0)
    tmpData(imod, :, :, 0) = (/A1/)
    tmpData(imod, :, :, 1) = (/A2trend * 10/)
    tmpData(imod, :, :, 2) = (/prob_A2/)
  end do
  XXdata(:dim_MOD1 - 1, :, :, :) = (/tmpData/)

  XXref = (/tmpData(iref, :, :, :)/)

  XXdata!0 = "line"
  XXdata!1 = "plev"
  XXdata!2 = "lat"
  XXdata!3 = "stat"
  XXdata&lat = lat0
  plev0@long_name = "Pressure [hPa]"
  XXdata&plev = plev0

  ; *****************************************************
  if (multimean) then
    do id = 0, dimsizes(exp_MMM) - 1
      flag_exp(:, id) = where(projects.eq.exp_MMM(id), 0, 1)
      index_exp = ind(flag_exp(:, id).gt.0)
      XXdata(dimTot - 1, :, :, :) = \
        (/dim_avg_n(XXdata(index_exp, :, :, :), 0)/)
      ; plotMod(dimTot - 1) = "MMM_" + exps_name(id)
      plotMod(dimTot - 1) = projects(flag_mod(0))
      years1(dimTot - 1) = max(y1)
      years2(dimTot - 1) = min(y2)
      delete(index_exp)
    end do
  end if

  X_MMM_std = (/dim_stddev_n(XXdata(index_mod, :, :, :), 0)/)
  plot_d = new(2, "graphic")

  ; *************************plotting******************************************
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "zonal_westerly_winds")
  res = True                  ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False

  res@lbLabelBarOn = False

  pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                           + "ipcc-ar6_wind_div.rgb")
  res@cnFillPalette = pal

  res@tiYAxisFont = "helvetica"
  res@tiMainFont = "helvetica"
  res@tmXBLabelFont = "helvetica"
  res@tmYLLabelFont = "helvetica"
  res@lbLabelFont = "helvetica"
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = (/-5., -2., -1., -0.5, -0.2, -0.1, 0., 0.1, 0.2, 0.5, 1., \
                  2., 5./)
  res@tmXBMode = "Explicit"
  res@tmXBValues = (/-90., -60., -30., 0./)
  res@tmXBLabels = (/"90~F34~0~F~S", "60~F34~0~F~S", "30~F34~0~F~S", \
                    "0~F34~0~F~"/)
  res@trXMinF = -90.
  res@trXMaxF = 0.
  res@tmYRMode = "Automatic"
  res@tmYRLabelsOn = False
  res@tmYLMode = "Explicit"
  res@tmYLValues = (/1000., 500., 100., 50., 10., 5., 1./)
  res@tmYLLabels = (/"1000", "500", "100", "50", "10", "5", "1"/)
  res@trYMinF = 1.
  res@trYMaxF = 1000.
  res@cnLineLabelsOn = False           ; turn off line labels
  res@cnFillOn = True                  ; turn on color fill
  res@cnLinesOn = False
  res@cnInfoLabelOn = False            ; turn off contour label
  res@gsnSpreadColors = True           ; use full range of colors
  res@lbLabelAutoStride = True         ; optimal labels
  res@gsnLeftString = ""

  res@gsnLeftStringFontHeightF = 0.017
  res@gsnRightStringFontHeightF = 0.017
  res@tmXBMode = "Explicit"            ; explicit labels

  res@gsnRightString = ""

  res@tiMainString = plotMod(iref)    ; title

  plot_d(0) = gsn_csm_pres_hgt(wks, XXdata(iref, :, :, 1), res)

  res@tiMainString = plotMod(dimTot - 1)    ; title
  res@tiYAxisString = ""

  plot_d(1) = gsn_csm_pres_hgt(wks, XXdata(dimTot - 1, :, :, 1), res)

  res2 = True
  res2@gsnDraw = False
  res2@gsnFrame = False
  res2@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res2@cnMinLevelValF = -20.              ; set min contour level
  res2@cnMaxLevelValF =  40.              ; set max contour level
  res2@cnLevelSpacingF =  10.              ; set contour spacing
  res2@cnLineLabelsOn = True
  res2@gsnContourZeroLineThicknessF = 2. 	; doubles thickness of zero contour
  res2@gsnContourNegLineDashPattern = 1
  res2@cnFillOn = False
  res2@cnInfoLabelOn = False
  res2@cnLinesOn = True
  res2@cnLineColor = "Black"             ; color of second contours
  res2@cnLineThicknessF = 2.             ; line thickness
  res2@tmYRMode = "Automatic"
  res2@tmXBMode = "Explicit"
  res2@tmXBValues = (/-90., -60., -30., 0./)
  res2@tmYLMode = "Explicit"
  res2@tmYLValues = (/1000., 500., 100., 50., 10., 5., 1./)
  res2@lbTitleFontHeightF = 0.01  ; default: 0.025
  res2@lbTitleOffsetF = -0.25  ; default: 0.03 (negative values move it down)
  res2@lbTopMarginF = 0.1  ; default: 0.05 (negative values move it up)
  res2@lbTitleJust = "BottomCenter"  ; default: CenterCenter
  res2@lbLabelFontHeightF = 0.025
  ; model specific ressources
  res2@tiMainFontHeightF = 0.03

  res2@tiMainString = plotMod(iref)
  plotv_0 = gsn_csm_pres_hgt(wks, XXdata(iref, :, :, 0), res2)

  res2@tiMainString = plotMod(dimTot - 1)
  plotv_1 = gsn_csm_pres_hgt(wks, XXdata(dimTot - 1, :, :, 0), res2)

  delete(res)
  res3 = True
  res3@gsnDraw = False
  res3@gsnFrame = False
  res3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
  res3@cnMinLevelValF = 0.              ; set min contour level
  res3@cnMaxLevelValF =  1.              ; set max contour level
  res3@cnLevelSpacingF =  0.005              ; set contour spacing
  res3@cnLineLabelsOn = False
  res3@cnInfoLabelOn = False
  res3@cnFillOn = False
  res3@cnLinesOn = False
  res3@tmYRMode = "Automatic"

  plotz_0 = gsn_csm_pres_hgt(wks, XXdata(iref, :, :, 2), res3)

  plotz_1 = gsn_csm_pres_hgt(wks, XXdata(dimTot - 1, :, :, 2), res3)

  opt = True
  opt@gsnShadeFillType = "pattern"      ; pattern fill
  opt@gsnShadeHigh = 17
  plotz_0 = gsn_contour_shade(plotz_0, 1., 0.95, opt)
  plotz_1 = gsn_contour_shade(plotz_1, 1., 0.95, opt)

  overlay(plot_d(0), plotv_0)

  overlay(plot_d(1), plotv_1)

  ; Call legend-plotting routine
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@amJust = "TopLeft"

  Pres@gsnPanelMainString = "Speed-up of zonal-mean wind (" + years1(iref) \
                            + "-" + years2(iref) + ")"
  Pres@gsnPanelXWhiteSpacePercent = 0.2

  Pres@gsnPanelFigureStrings = (/"a)", "b)"/)
  Pres@gsnPanelFigureStringsPerimOn = False
  Pres@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  Pres@gsnPanelFigureStringsJust = "TopLeft"
  Pres@gsnPanelFigureStringsFontHeightF = 0.016

  Pres@gsnPanelLabelBar    = True                ; add common colorbar
  Pres@lbLabelFontHeightF  = 0.007               ; make labels smaller
  Pres@lbBoxEndCapStyle = "TriangleBothEnds"
  Pres@lbBoxLineColor = "gray10"
  Pres@lbTitleString  = "zonal wind (" + UNITS + "/decade)"
  Pres@lbTitlePosition  = "Bottom"
  Pres@lbLabelFontHeightF  = 0.012
  Pres@pmLabelBarWidthF      = 0.7
  Pres@lbTitleFontHeightF  = 0.016

  gsn_panel(wks, plot_d, (/1, 2/), Pres)

  frame(wks)

  ; *******************OUTPUT FILE******************************************
  nc_filename = work_dir + "eyring13f12_" + var0 + ".nc"
  Xdata = XXdata(line|:, plev|:, lat|:, stat|:)
  Xdata@var = var0
  Xdata@diag_script = "eyring13f12"
  Xdata&stat = (/"mean", "trend", "prob"/)
  Xdata&line = plotMod
  Xdata&plev = plev0
  Xdata&lat = lat0
  nc_outfile = ncdf_write(Xdata, nc_filename)

  ; Call provenance logger
  plotpath = config_user_info@plot_dir + "zonal_westerly_winds"
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Long-term mean and trend of " + var0, \
                 (/"mean", "trend"/), \
                 (/"global"/), \
                 "pro", \
                 (/"bock_lisa", "cionni_irene", "eyring_veronika", \
                  "hassler_birgit"/), \
                 (/"eyring13jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
