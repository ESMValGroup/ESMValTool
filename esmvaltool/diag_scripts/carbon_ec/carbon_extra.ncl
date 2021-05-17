; #######################################################################
; carbon_beta.ncl
; Author: Sabrina Zechlau (DLR, Germany)
; #######################################################################
; Description:
; creates two panels to diagnos beta GPP at 2x[CO2] increase
; Panel: a) correlation between GPP and [co2] increase,
;        b) barchart of the trend of the upper correlation
;
; Required info attributes (scripts):
; - styleset       project for line, color, symbol styles
;
; Optional info attributes (scripts):
; - bc_xmax_year   end year (default = last year of all model datasets
;                  available)
; - bc_xmin_year   start year (default = first year of all model datasets
;                  available)
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
  tmpv = get_unique_values(tmp)
  var = "gpp"
  gpp_bgc = ind(tmp.eq."gpp_bgc")
  gpp_hist = ind(tmp.eq."gpp_hist")
  co2s = ind(tmp.eq."co2s")
  ;co2s_amp = ind(tmp.eq."co2s_amp")

  delete(tmp)

  ; Load var_info
  infoa = select_metadata_by_name(input_file_info, "gpp")
  info0 = infoa[gpp_bgc]
  info1 = infoa[gpp_hist]

  datasetnames = metadata_att_as_array(info0, "dataset")
  dim_MOD = ListCount(info0)

  infob = select_metadata_by_name(input_file_info, "co2s")
  info2 = infob[co2s-(2*dim_MOD)]
  ;info3 = infob[co2s_amp-(2*dim_MOD)]

  delete(infoa)
  delete(infob)

  log_info("+++++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + tmpv(0) + " " + tmpv(1) + " " + tmpv(2) + ")")
  log_info("+++++++++++++++++++++++++++++++++++++++++++++")
  delete(tmpv)

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
  if (isatt(diag_script_info, "bc_xmax_year")) then
    xMax_year = toint(diag_script_info@bc_xmax_year)
  else
    xMax_year = max(metadata_att_as_array(info, "end_year"))
  end if
  if (isatt(diag_script_info, "bc_xmin_year")) then
    xMin_year = toint(diag_script_info@bc_xmin_year)
  else
    xMin_year = max(metadata_att_as_array(info, "start_year"))
  end if

  ; Call plot scripts
  plot_file = DIAG_SCRIPT + "_" + xMin_year + "-" + xMax_year
  wks = gsn_open_wks(file_type, plot_dir + plot_file)
  colors  = project_style(info0, diag_script_info, "colors")
  markers = project_style(info0, diag_script_info, "markers")
  thicks  = project_style(info0, diag_script_info, "avgstd")
  lgLabels = datasetnames

  ; ===============================================================
  ; setup data holders for plots
  nyMax = max(metadata_att_as_array(info0, "end_year"))
  nyMin = max(metadata_att_as_array(info0, "start_year"))
  nMax = nyMax - nyMin + 11

  pctco = new((/nMax/), double)
  pctco(0) = 285
  do tt = 1, nMax-1
    pctco(tt) = pctco(tt-1) + 0.01 * pctco(tt-1)
  end do
  pctco!0 = "year"
  pctco&year = ispan(nyMin - 10, nyMax, 1)
  ; --------------------------------------------------------------
  ; define arrey
  nfit = fspan(0, 50, 50)
  rcgpp_bgc  = new((/dim_MOD/), float)
  rcgpp_hist = new((/dim_MOD/), float)
  aY0data = new((/dim_MOD, nMax/), double)
  aY1data = new((/dim_MOD, nMax+6/), double)
  aY2data = new((/dim_MOD, nMax+6/), double)
  aY3data = new((/dim_MOD, nMax+6/), double)
  mean_a1 = new((/dim_MOD, nMax+6/), double)
  Yfitda  = new((/dim_MOD, nMax/), double)
  Yfit_hist   = new((/dim_MOD, nMax+6/), double)
  Yfit_change = new((/dim_MOD, dimsizes(nfit)/), double)
  Xfit_change = new((/dim_MOD, dimsizes(nfit)/), double)

  ; loop for models
  do imod = 0, dim_MOD-1

    ; Read data
    A0 = read_data(info0[imod])
    A0&time = pctco&year(10:)
    A1 = read_data(info1[imod])
    A2 = read_data(info2[imod])
    ;A3 = read_data(info3[imod])

    ; calculate annual mean for VAR:
    ymin = 0
    ymax = toint(dimsizes(A0)-1)
    tmp0 = (A0 * 3600. * 24. * 365) / 1e12
    copy_VarMeta(A0, tmp0)
    aY0data(imod, ymin:ymax) = (/tmp0/)

    rc_bgc = regline_stats(pctco(ymin:ymax), aY0data(imod, ymin:ymax))
    Yfitda(imod, ymin:ymax) = rc_bgc@Yest

    rcgpp_bgc(imod) = tofloat(rc_bgc)

    delete([/tmp0, ymin, ymax, A0, rc_bgc/])

    ; calculate annual mean for VAR:
    ymin = 0
    ymax = toint(dimsizes(A1)-1)
    tmp0 = (A1 * 3600. * 24. * 365) / 1e12
    copy_VarMeta(A1, tmp0)
    aY1data(imod, ymin:ymax) = (/tmp0/)
    aY2data(imod, ymin:ymax) = (/A2/)
    ;aY3data(imod, ymin:ymax) = A3 - dim_avg(A3)

    rc_hist= regline_stats(aY2data(imod, ymin:ymax), aY1data(imod, ymin:ymax))
    Yfit_hist(imod, ymin:ymax) = rc_hist@Yest
    rcgpp_hist(imod) = tofloat(rc_hist)

    ; calculate regression linse for Fig 2
    ;mean_a1(imod, ymin:ymax) = aY1data(imod, ymin:ymax) - dim_avg(aY1data(imod, ymin:ymax))
    ;rc_change= regline_stats(aY3data(imod, ymin:ymax), mean_a1(imod, ymin:ymax))
    ;Yfit_change(imod, ymin:ymax) = rc_change@Yest

    ;minx = min(aY3data(imod, ymin:ymax))-1.1*(max(aY3data(imod, ymin:ymax)) - min(aY3data(imod, ymin:ymax)))
    ;maxx = max(aY3data(imod, ymin:ymax))+1.1*(max(aY3data(imod, ymin:ymax)) + max(aY3data(imod, ymin:ymax))) 

    ;Xfit_change(imod, 0:dimsizes(nfit)-1) = minx + (maxx - minx) * nfit/dimsizes(nfit)
    ;Yfit_change(imod, 0:dimsizes(nfit)-1) = rc_change * Xfit_change(imod, 0:dimsizes(nfit)-1) + \
    ;                                        rc_change@yintercept

    delete([/tmp0, ymin, ymax, A1, A2, rc_hist/])

  end do     ; model loop

  rcgppco2 = regline_stats(rcgpp_bgc, rcgpp_hist)
  Yfit_total = rcgppco2@Yest

  ; -----------------------------------------------------------
  ; creat plots

  ; PLOT - Ext. Fig 2
  ;XStg    = "Change in CO~B~2~N~ amplitude, BRW [ppmv]"
  ;YStg    = "Change in Mean GPP [GtC/yr], 60N-90N"
  ;data_arr = (/aY3data(0, :), mean_a1(0, :)/)
  ;data_arr!0 = "datasets"
  ;data_arr&datasets = (/XStg, YStg/)
  ;data_arr@res_tiMainString    = ""
  ;data_arr@res_tiMainFontHeightF = 0.025
  ;data_arr@res_gsnCenterString = ""
  ;data_arr@res_gsnMaximize     = True
  ;data_arr@res_xyMonoMarker    = False
  ;data_arr@res_xyMarker        = 1
  ;data_arr@res_xyMarkerColors  = "white"
  ;data_arr@res_trXMinF         = -10 
  ;data_arr@res_trXMaxF         = 14
  ;data_arr@res_trYMinF         = -2
  ;data_arr@res_trYMaxF         = 4
  ;data_arr@res_vpWidthF        = 0.5
  ;data_arr@res_vpHeightF       = 0.4

  ;data_arr@diag_script = DIAG_SCRIPT
  ;diag_script_info@scatter_log = False

  ;plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

  ;txres = True
  ;MarkerSizeF = where(datasetnames.eq."CESM1-BGC", 20., 17.)
  ;txres@gsMarkerThicknessF = 3

  ;lineres = True
  ;lineres@gsLineDashPattern = 0
  ;lineres@gsLineThicknessF  = 3
  ;linehist = lineres

  ; Add markers
  ;do imod  = 0, dim_MOD-1
  ;  marker = unique_string("marker")
  ;  txres@gsMarkerColor = colors(imod)
  ;  txres@gsMarkerSizeF = 10. ;MarkerSizeF(imod)
  ;  txres@gsMarkerIndex = 5 ;circl
  ;  add_markers(wks, plot, txres, aY3data(imod, :), mean_a1(imod, :))

  ;  lineres@gsLineColor = colors(imod)
  ;  plot@$marker$ = gsn_add_polyline(wks, plot, Xfit_change(imod, :), Yfit_change(imod, :), lineres)

  ;end do
  ;draw(plot)
  ;frame(wks)
  ;delete([/plot, XStg, YStg, data_arr, marker, txres/])

  ; PLOT - Panel A
  XStg    = "CO~B~2~N~ at BRW [ppmv]"
  YStg    = "GPP [GtC/yr], 60N-90N"
  data_arr = (/pctco, aY0data(0, :)/)
  data_arr!0 = "datasets"
  data_arr&datasets = (/XStg, YStg/)
  data_arr@res_tiMainString    = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString = ""
  data_arr@res_gsnMaximize     = True
  data_arr@res_xyMonoMarker    = False
  data_arr@res_xyMarker        = 1
  data_arr@res_xyMarkerColors  = "white"
  data_arr@res_trXMinF         = 280
  data_arr@res_trXMaxF         = 400
  data_arr@res_trYMinF         = 0
  data_arr@res_trYMaxF         = 15
  data_arr@res_vpWidthF        = 0.5
  data_arr@res_vpHeightF       = 0.4

  data_arr@diag_script = DIAG_SCRIPT
  diag_script_info@scatter_log = False

  plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

  txres = True
  MarkerSizeF = where(datasetnames.eq."CESM1-BGC", 15., 10.)
  txres@gsMarkerThicknessF = 3

  lineres = True
  lineres@gsLineDashPattern = 0
  lineres@gsLineThicknessF  = 3
  linehist = lineres

  ; Add markers
  do imod  = 0, dim_MOD-1
    marker1 = unique_string("marker1")
    marker2 = unique_string("marker2")
    txres@gsMarkerColor = colors(imod)
    txres@gsMarkerSizeF = MarkerSizeF(imod)
    txres@gsMarkerIndex = 4 ;circl
    add_markers(wks, plot, txres, pctco, aY0data(imod, :))

    lineres@gsLineColor = colors(imod)
    plot@$marker1$ = gsn_add_polyline(wks, plot, pctco, Yfitda(imod, :), lineres)

    txres@gsMarkerIndex = 5 ;asterisks
    linehist@gsLineColor = colors(imod)
    add_markers(wks, plot, txres, aY2data(imod, :), aY1data(imod, :))
    plot@$marker2$ = gsn_add_polyline(wks, plot, aY2data(imod, :), Yfit_hist(imod, :), linehist)

  end do
  draw(plot)
  frame(wks)
  delete([/plot, XStg, YStg, data_arr, marker1, marker2, txres/])


  ; PLOT - Panel B
  XStg    = "dGPP/dCO~B~2~N~ at BRW [GtC/yr/ppmv], BGC"
  YStg    = "dGPP/dCO~B~2~N~ at BRW [GtC/yr/ppmv], Hist"

  data_arr = new((/2, dimsizes(rcgpp_bgc)/), float)
  data_arr!0      = "datasets"
  data_arr!1      = "model"
  data_arr&datasets = (/XStg, YStg/)
  data_arr&model = datasetnames
  data_arr(0, :) = tofloat(rcgpp_bgc)
  data_arr(1, :) = tofloat(rcgpp_hist)
  data_arr@res_tiMainString    = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString = ""
  data_arr@res_gsnMaximize     = True
  data_arr@res_xyMonoMarker    = False
  data_arr@res_xyMarkers       = markers
  data_arr@res_xyMarkerColors  = "white";colors
  data_arr@res_trXMinF         = -0.0005 ;min(rcgpp_bgc)- 0.25 * min(rcgpp_bgc) ;
  data_arr@res_trXMaxF         = 0.035;max(rcgpp_bgc)+ 0.25 * max(rcgpp_bgc) ;0.035
  data_arr@res_trYMinF         = -0.0005 ;min(rcgpp_hist)- 0.25 * min(rcgpp_hist) ;data_arr@res_trXMinF
  data_arr@res_trYMaxF         = 0.035 ;max(rcgpp_hist)+ 0.25 * max(rcgpp_hist) ;data_arr@res_trXMaxF
  data_arr@res_vpWidthF        = 0.5
  data_arr@res_vpHeightF       = 0.4
  data_arr@diag_script = DIAG_SCRIPT
  diag_script_info@scatter_log = False

  plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

  ;Add markers
  txres = True
  MarkerSizeF = where(datasetnames.eq."CESM1-BGC", 20., 17.)
  ncyc = (/"ACCESS-ESM1-5", "CESM2", "MIROC-ES2L", "MPI-ESM1-2-LR", "NorESM2-LM", "UKESM1-0-LL"/)
  txres@gsMarkerThicknessF = 4

  do imod  = 0, dim_MOD-1
    txres@gsMarkerColor = colors(imod)
    txres@gsMarkerSizeF = MarkerSizeF(imod)
    if (any(datasetnames(imod) .eq. ncyc)) then
      txres@gsMarkerIndex = 5 ;markers(imod)
    else
      txres@gsMarkerIndex = 4
    end if
    markers(imod) = txres@gsMarkerIndex
      
    add_markers(wks, plot, txres,  rcgpp_bgc(imod), rcgpp_hist(imod))
  end do

  ; add one-to-one line
  o5lres = True
  o5lres@gsLineColor       = "black"
  o5lres@gsLineDashPattern = 2
  o5lres@gsLineThicknessF  = 2
  CMIP5_fit = gsn_add_polyline(wks, plot, (/-20, 20/), (/-20, 20/), o5lres)

  ;add total reg line
  tres = True
  tres@gsLineColor       = "red"
  tres@gsLineDashPattern = 1
  tres@gsLineThicknessF  = 2
  totalline = gsn_add_polyline(wks, plot, rcgpp_bgc, Yfit_total, tres)

  draw(plot)
  frame(wks)
  delete([/plot, XStg, YStg, data_arr/])

  ; create separate legend
  marker_thicks = datasetnames
  marker_sizes  = datasetnames

  marker_thicks = 3
  marker_sizes  = 0.02

  leg = True
  leg@txFontQuality = "High"
  leg@txFont        = 25
  leg@txFontHeightF = 0.02
  leg@diag_script   = DIAG_SCRIPT
  leg@annots        = datasetnames(::-1)
  leg@colors        = colors(::-1)
  leg@markers       = markers(::-1)
  leg@thicks        = marker_thicks
  leg@sizes         = marker_sizes
  leg@ncols         = 1

  create_legend_lines(leg@annots, leg, plot_dir + \
                      DIAG_SCRIPT + "_legend", "markers")

  plotname = plot_dir + plot_file + "." + file_type

  ; Call provenance logger
;  log_provenance(ncdf_outfile, \
;                 plotname, \
;                 "Climate models vs " + YStg, \
;                 (/"anomaly", "corr", "stddev"/), \
;                 (/"global"/),\
;                 (/""/), \
;                 (/"zechlau_sabrina"/), \
;                 (/"wenzel16nat"/), \
;,                 metadata_att_as_array(info, "filename"))

end
