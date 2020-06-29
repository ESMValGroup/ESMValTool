; #######################################################################
; carbon_co2_cycle.ncl
; Author: Sabrina Zechlau (DLR, Germany)
; #######################################################################
; Description:
; Creates Figure 1 of Wenzel et al. 2016 Nature:
;         Panel: a) correlation of the [co2] seasonal cycle amplitude and
;                   [co2] increase
;                b) barchart of the trend of the upper correlation
; and Figure 3 of Wenzel et al. 2016 Nature:
;         Panel: a) Emergent Constraint between beta and the trend of the
;                   [co2] seasonal cycle amplitude
;                b) conditional PDF of the EC
;
; Required info attributes (scripts):
; - styleset       project for line, color, symbol styles
; - nc_infile      path of netCDF file containing beta
;                  (output from carbon_beta.ncl)
;
; Optional info attributes (scripts):
; - bc_xmax_year   end year (default = last year of all model datasets
;                  available)
; - bc_xmin_year   start year (default = first year of all model datasets
;                  available)
;
; Required variable attributes (variables):
; - reference_dataset: name of reference datatset (observations)
;
; Modification history
;    * 20200406-zechlau_sabrina: code rewritten for ESMValTool v2.0
;    * 2015xxxx-wenzel_sabrina: written
; ########################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/carbon_ec/carbon_aux.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"
load "$diag_scripts/shared/plot/legends.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/carbon_plots.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  tmp = metadata_att_as_array(input_file_info, "variable_group")
  co2s_ann = ind(tmp.eq."co2s")
  co2s_amp = ind(tmp.eq."co2s_amp")
  co2s_ann_obs = ind(tmp.eq."co2s_obs")
  co2s_amp_obs = ind(tmp.eq."co2s_amp_obs")
  delete(tmp)

  ; Load var_info
  var = variable_info[0]@short_name
  info = select_metadata_by_name(input_file_info, var)
  info0 = info[co2s_ann]
  ListPush(info0, info[co2s_ann_obs])
  info1 = info[co2s_amp]
  ListPush(info1, info[co2s_amp_obs])

  datasetnames = metadata_att_as_array(info0, "dataset")
  dim_MOD = ListCount(info0)

  reference_dataset = variable_info[0]@reference_dataset
  iref = ind(reference_dataset .eq. datasetnames)
  idat = ind(reference_dataset .ne. datasetnames)

  log_info("+++++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var + ")")
  log_info("+++++++++++++++++++++++++++++++++++++++++++++")

  ; ******************************************************************
  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  system("mkdir -p " + config_user_info@work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Check for required settings
  exit_if_missing_atts(diag_script_info, "styleset")

  ; optional input parameters
  if(isatt(diag_script_info, "bc_xmax_year")) then
    xMax_year = toint(diag_script_info@bc_xmax_year)
  else
    xMax_year = max(metadata_att_as_array(info, "end_year"))
  end if
  if(isatt(diag_script_info, "bc_xmin_year")) then
    xMin_year = toint(diag_script_info@bc_xmin_year)
  else
    xMin_year = min(metadata_att_as_array(info, "start_year"))
  end if

  ; Call plot scripts
  plot_file = "amplitude_" + var + "_" + xMin_year + "-" + xMax_year
  wks = gsn_open_wks(file_type, plot_dir + plot_file)
  colors  = project_style(info0, diag_script_info, "colors")
  markers = project_style(info0, diag_script_info, "markers")
  thicks  = project_style(info0, diag_script_info, "avgstd")
  lgLabels = datasetnames

  ; path to beta file
  temp_dir = config_user_info@work_dir + "/" + diag_script_info@nc_infile
  temp_dir = temp_dir + "beta_" + xMax_year + "-" + xMin_year + ".nc"
  ; ==================================================================
  ; setup data holders for plots
  nyMax = max(metadata_att_as_array(info0, "end_year"))
  nyMin = min(metadata_att_as_array(info0, "start_year"))
  nMax = nyMax - nyMin + 1

  Ydata  = new((/dim_MOD, nMax/), float)
  Yfitd = new((/dim_MOD, nMax/), float)
  Xdata = new((/dim_MOD, nMax/), float)
  rcc   = new((/dim_MOD, 2/), float)
  stdrc = new((/dim_MOD, 2/), float)
  Xdata!1 = "year"
  Xdata&year = ispan(nyMin+1, nyMax + 1, 1)

  ; ##############################################################
  ; loop for models
  do imod = 0, dim_MOD-1

    ; Read data
    A0 = read_data(info0[imod])
    A1 = read_data(info1[imod])

    ; Save input data with one year leap
    Xdata(imod, :dimsizes(A0)-2) = tofloat(A0(1:))
    Ydata(imod, :dimsizes(A1)-1) = tofloat(A1)

    ; Regression between [CO2] Amplitude (x) and [CO2] ann mean (Y)
    rc = regline_stats(Xdata(imod, :), Ydata(imod, :))
    rcc(imod, 1) = rc
    rcc(imod, 0) = rc@b(0)
    stdrc(imod, :) = rc@stderr
    Yfitd(imod, :dimsizes(rc@Yest)-1) = rc@Yest

    delete([/rc,  A1, A0/])
  end do     ; model loop

  ; read file with beta values
  beta = ncdf_read(temp_dir, "beta")

  ; Save in array w/o OBS
  betaGPP = new((/2, dim_MOD-1/), float)
  rcGPP   = new((/(dim_MOD-1)*2/), float)

  betaGPP(:, :) = beta
  rcGPP(1::2) = betaGPP@rcgpp

  rcGPP!0 = "models"
  tmpm = new((/(dim_MOD-1)*2/), string)
  tmpm(::2) = datasetnames(idat)
  tmpm(1::2) = datasetnames(idat)
  rcGPP&models = tmpm

  delete(beta)
  delete(tmpm)

  ; Calculate constraint on BETA_GPP and
  nfit = fspan(0, 50, 50)

  ; amplitude trends
  minx = min(rcc(idat, 1))-1.1*(max(rcc(idat, 1)) - min(rcc(idat, 1)))
  maxx = max(rcc(idat, 1))+1.1*(max(rcc(idat, 1)) - min(rcc(idat, 1)))

  rc = regline_stats(rcc(idat, 1), betaGPP(0, :))
  arc_bco2  = rc@r
  prc_bco2  = rc@F_pval
  xfit_bco2 = minx + (maxx - minx) * nfit/dimsizes(nfit)
  yfit_bco2 = rc * xfit_bco2 + rc@yintercept
  beta_mco2 = rc * rcc(iref, 1) + rc@yintercept
  beta_eco2 = sqrt(stdrc(iref, 1) ^ 2 + rc@rstd ^ 2)

  ; Calculate prior distribution
  PDF = Cond_PDF(xfit_bco2, yfit_bco2, rcc(idat, 1), \
                 betaGPP(0, :), rc@MSE, \
                 (/rcc(iref, 1), stdrc(iref, 1)/))

  if (stdrc(iref, 1).ne.0) then
    co2_Py  = PDF@Py
  else
    co2_Py  = 0.
  end if
  co2_Py_pr = PDF@Py_pr
  co2_y     = PDF@y

  ; this is for plotting the regression line with confidence intervalls
  yfit_bco2_std = new((/2, dimsizes(yfit_bco2)/), typeof(yfit_bco2))
  yfit_bco2_std(0, :) = yfit_bco2 + PDF@sigma
  yfit_bco2_std(1, :) = yfit_bco2 - PDF@sigma

  yBETA = new((/2, dim_MOD+1/), float)
  yBETA(0, 0)  = avg(betaGPP(0, :))
  yBETA(0, 1)  = beta_mco2
  yBETA(0, 2:) = betaGPP(0, :)
  yBETA(1, 0)  = stddev(betaGPP(0, :))
  yBETA(1, 1)  = beta_eco2
  yBETA(1, 2:) = betaGPP(1, :)

  xCO2 =  new((/2, dim_MOD+1/), float)
  xCO2(0, 1)  = rcc(iref, 1)
  xCO2(0, 2:) = rcc(idat, 1)
  xCO2(1, 1)  = stdrc(iref, 1)
  xCO2(1, 2:) = stdrc(idat, 1)

  delete(rc)
  delete(PDF)
  delete(minx)
  delete(maxx)

  ; Plot section
  ; -----------------------------------------------------------

  marker_thicks = new(dimsizes(datasetnames), float)
  marker_sizes  = new(dimsizes(datasetnames), float)
  marker_thicks = 3
  marker_sizes  = 0.02
  marker   = new(dim_MOD, "graphic")

  ; Create Plots 1a:
  ; -----------------------------------------------------------
  XStg    = "Annual CO~B~2~N~ [ppmv]"
  YStg    = "CO~B~2~N~ Amplitude [ppmv]"
  data_arr = (/Xdata(0, :), Ydata(0, :)/)
  data_arr!0 = "datasets"
  data_arr&datasets = (/XStg, YStg/)
  data_arr@res_tiMainString      = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString   = ""
  data_arr@res_gsnMaximize       = True
  data_arr@res_xyMarkerSizeF     = 1.
  data_arr@res_trXMinF           = min(Xdata(:, :))
  data_arr@res_trXMaxF           = max(Xdata(:, :))
  data_arr@res_trYMinF           = 0
  data_arr@res_trYMaxF = max(Ydata(:, :)) + 0.25 * max(Ydata(:, :))
  data_arr@res_vpWidthF          = 0.5
  data_arr@res_vpHeightF         = 0.4
  data_arr@diag_script = DIAG_SCRIPT
  diag_script_info@scatter_log = False

  data_arr@res_xyMarkLineMode = "Markers"

  plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

  txres = True
  txres@gsMarkerSizeF      = 10.
  txres@gsMarkerThicknessF = 3

  ; Add year strings as markers
  lineres = True
  do imod  = 0, dim_MOD - 1
    idx := ind(.not.ismissing(Xdata(imod, :)))
    idy := ind(.not.ismissing(Ydata(imod, :)))
    if (datasetnames(imod).eq."OBS") then
      lineres@gsLineThicknessF  = 4
    else
      lineres@gsLineThicknessF  = 3
    end if
    txres@gsMarkerColor      = colors(imod)
    txres@gsMarkerIndex      = markers(imod)
    add_markers(wks, plot, txres, Xdata(imod, idx), Ydata(imod, idy))

    lineres@gsLineDashPattern = 0
    lineres@gsLineColor = colors(imod)
    marker(imod) = gsn_add_polyline(wks, plot, Xdata(imod, :), \
                                    Yfitd(imod, :), lineres)
  end do

  draw(plot)
  frame(wks)
  delete([/plot, XStg, data_arr/])

  ; Create Plot 1b:
  ; ----------------------------------------------------------
  rcc@info = info0

  YStg    = "Sensitivity of CO~B~2~N~ Amplitude to CO~B~2~N~"

  data_arr = new((/dim_MOD, 1/), float)
  data_arr(:, 0) = rcc(:, 1)
  data_arr!0 = "datasets"
  data_arr&datasets(0) = datasetnames
  data_arr!1 = "yaxis"
  data_arr&yaxis = YStg
  data_arr@res_tiMainString    = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString = ""

  plot = barchart(wks, data_arr, var, info0)
  draw(plot)
  frame(wks)

  delete([/plot, data_arr/])

  ; Create Plots 3a:
  ; ----------------------------------------------------------
  tmpm = datasetnames(idat)
  delete(datasetnames)
  datasetnames = new(dim_MOD+1, typeof(tmpm))
  datasetnames(0) = "OBS"
  datasetnames(1) = "constraint mean"
  datasetnames(2:) = tmpm
  delete(tmpm)

  xmarkers = new(dim_MOD+1, typeof(markers))
  xmarkers(0) = 16
  xmarkers(1) = 3
  xmarkers(2:) = markers(idat)

  xcolors = new(dim_MOD+1, typeof(colors))
  xcolors(0) = "black"
  xcolors(1) = "red"
  xcolors(2:) = colors(idat)

  MainStg = ""
  XStg    = "Sensitivity of CO~B~2~N~ Amplitude to CO~B~2~N~"
  YStg    = "GPP(2xCO~B~2~N~)/GPP(1xCO~B~2~N~)"

  data_arr = (/xCO2(0, :), yBETA(0, :)/)
  data_arr!0 = "datasets"
  data_arr&datasets = (/XStg, YStg/)
  data_arr@res_tiMainString      = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString   = ""
  data_arr@res_gsnMaximize       = True
  data_arr@res_xyMarkerSizeF     = 1.
  data_arr@res_trXMinF           = min(xCO2(0, :))
  data_arr@res_trXMaxF           = 0.1  ; max(xCO2(0, :))
  data_arr@res_trYMinF           = 1.1  ; 0
  data_arr@res_trYMaxF = max(yBETA(0, :)) + 0.25 * max(yBETA(0, :))
  data_arr@res_vpWidthF          = 0.5
  data_arr@res_vpHeightF         = 0.4
  data_arr@diag_script = DIAG_SCRIPT
  diag_script_info@scatter_log = False

  plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

  ; add confidence intervall of regression line
  peres = True
  peres@tfPolyDrawOrder = "PreDraw"
  peres@xyLineColor     = "orange"
  peres@gsFillColor     = "orange"
  xin = array_append_record(xfit_bco2, xfit_bco2(::-1), 0)
  yin = array_append_record(yfit_bco2_std(0, :), yfit_bco2_std(1, ::-1), 0)
  add_prediction_error(wks, plot, peres, xin, yin)
  delete([/xin, yin/])

  if (rcc(iref, 1).ne.0) then
    ores = True
    ores@tfPolyDrawOrder = "PreDraw"
    ores@xyLineColor     = "grey"
    ores@gsFillColor     = "grey"
    xinobs = (/rcc(iref, 1) + stdrc(iref, 1), \
              rcc(iref, 1) + stdrc(iref, 1), \
              rcc(iref, 1) - stdrc(iref, 1), \
              rcc(iref, 1) - stdrc(iref, 1)/)
    yinobs = (/500, -500, -500, 500/)
    add_prediction_error(wks, plot, ores, xinobs, yinobs)
    delete([/xinobs, yinobs, ores/])
  end if

  ; add regression line
  lres = True
  lres@gsLineColor       = "red"
  lres@gsLineThicknessF  = 1.5
  fit_line = gsn_add_polyline(wks, plot, xfit_bco2, yfit_bco2, lres)

  ; add colored markers
  mres = True
  mres@gsMarkerSizeF = 20.
  mres@gsMarkerThicknessF = 3.
  do imod = 1, dim_MOD - 1
    mres@gsMarkerColor = xcolors(imod)
    mres@gsMarkerIndex = xmarkers(imod)
    add_markers(wks, plot, mres, xCO2(0, imod), yBETA(0, imod))
  end do

  ; add xy error bars
  eres = True
  eres@gsLineColor = colors(idat)
  add_errorbar(wks, plot, eres, xCO2(:, 2:), yBETA(:, 2:))

  draw(plot)
  frame(wks)

  delete([/plot, peres, eres, mres, lres, data_arr/])

  ; Create Plots 3b:
  ; -----------------------------------------------
  plot_pdf = new(1, graphic)

  bres  = True
  bres@gsnDraw          = False
  bres@gsnFrame         = False
  bres@xyLineThicknessF = 2.
  bres@xyLineColor      = (/"black"/)
  bres@tiXAxisString    = "GPP(2xCO~B~2~N~)/GPP(1xCO~B~2~N~)"
  bres@tiYAxisString    = "Probability Density"
  bres@trXMinF          = 0.6
  bres@trXMaxF          = 2.0
  bres@trYMinF          = 0.
  bres@trYMaxF          = max(co2_Py)*1.2
  bres@gsnXYBarChart    = True
  bres@tiMainString     = ""
  bres@vpWidthF         = 0.5
  bres@vpHeightF        = 0.4
  bres@tmXTLabelFontHeightF = 0.025
  bres@tmYLLabelFontHeightF = 0.025
  bres@tiXAxisFontHeightF   = 0.025
  bres@tiYAxisFontHeightF   = 0.022

  ; create line plot with PDFs
  diag_script_info@multi_model_mean = False
  diag_script_info@scatter_log = False
  diag_script_info@xy_line_legend = False

  pres = True
  pres@tiMainString       = ""
  pres@tiXAxisString      = XStg
  pres@tiYAxisString      = "Probability Density"
  pres@vpWidthF           = 0.4
  if (rcc(iref, 1).eq.0) then
    pres@xyLineColors   = "black"
    pres@xyDashPatterns = 2
    pres@trYMinF        = 0
      xy_line(wks, co2_Py_pr, co2_y, 0., pres, input_file_info)
  else
    pres@xyLineColors   = (/"black", "red"/)
    pres@xyDashPatterns = (/2, 0/)
    pres@trYMinF        = 0
      xy_line(wks, (/co2_Py_pr, co2_Py/), (/co2_y, co2_y/), \
              0., pres, input_file_info)
  end if

  obsres = True
  obsres@gsLineColor = "red"
  obsres@gsLineThicknessF  = 2
  PDF3 = gsn_add_polyline(wks, plot_pdf, co2_y, co2_Py, obsres)

  draw(plot_pdf)
  frame(wks)
  delete(plot_pdf)

  ; create separate legend

  leg = True
  leg@txFontQuality = "High"
  leg@txFont        = 25
  leg@txFontHeightF = 0.02
  leg@diag_script   = DIAG_SCRIPT
  leg@annots        = lgLabels(::-1)
  leg@markers       = markers(::-1)
  leg@thicks        = marker_thicks
  leg@sizes         = marker_sizes
  leg@ncols         = 1
  leg@colors        = colors(::-1)  ; rgbcolors

  create_legend_lines(leg@annots, leg, plot_dir + \
                      DIAG_SCRIPT + "_legend", "markers")

  ; ---------------------------------------------------------------
  ; Output to NetCDF
  new_path = config_user_info@work_dir
  new_path = new_path + "amplitude_" + xMax_year + "-" + xMin_year + ".nc"

  ; Attach attributes to the results
  CO2var = new((/2, dim_MOD+1/), float)
  CO2var(0, :) = xCO2(0, :)
  CO2var(1, :) = xCO2(1, :)
  CO2var!0        = "case"
  CO2var&case     = (/"mean", "stddev"/)
  CO2var!1        = "model"
  CO2var&model    = datasetnames
  CO2var@ncdf     = new_path
  CO2var@var      = var

  CO2var@diag_script = (/DIAG_SCRIPT/)

  ; Write NetCDF output
  ncdf_outfile = ncdf_write(CO2var, new_path)

  plotname = plot_dir + plot_file + "." + file_type

  ; -----------------------------------------------------------
  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 plotname, \
                 XStg + " vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"global"/),\
                 (/""/), \
                 (/"zechlau_sabrina"/), \
                 (/"wenzel16nat"/), \
                 metadata_att_as_array(info, "filename"))

end
; ###############################################################
