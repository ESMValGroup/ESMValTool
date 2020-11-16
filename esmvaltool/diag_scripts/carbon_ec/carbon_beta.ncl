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

  tmp = metadata_att_as_array(input_file_info, "short_name")
  variables = get_unique_values(tmp)
  delete(tmp)

  ; Load var_info
  var = variable_info[0]@short_name
  info = select_metadata_by_name(input_file_info, var)
  datasetnames = metadata_att_as_array(info, "dataset")
  dim_MOD = ListCount(info)

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
  plot_file = "beta_" + xMin_year + "-" + xMax_year
  wks = gsn_open_wks(file_type, plot_dir + plot_file)
  colors  = project_style(info, diag_script_info, "colors")
  markers = project_style(info, diag_script_info, "markers")
  thicks  = project_style(info, diag_script_info, "avgstd")
  lgLabels = datasetnames

  ; ===============================================================
  ; setup data holders for plots
  nyMax = max(metadata_att_as_array(info, "end_year"))
  nyMin = max(metadata_att_as_array(info, "start_year"))
  nMax = nyMax - nyMin + 1

  pctco = new((/nMax+10/), double)
  pctco(0) = 285
  do tt = 1, nMax+10-1
    pctco(tt) = pctco(tt-1) + 0.01 * pctco(tt-1)
  end do
  pctco!0 = "year"
  pctco&year = ispan(nyMin-10, nyMax, 1)

  ; --------------------------------------------------------------
  ; define arrey
  betaGPP = new((/dim_MOD/), float)
  delta0  = new((/dim_MOD/), float)
  delta0_std = new((/dim_MOD/), float)
  rcgpp   = new((/dim_MOD/), float)
  aY0data = new((/dim_MOD, nMax/), double)
  Yfitda  = new((/dim_MOD, nMax/), double)

  ; loop for models
  do imod = 0, dim_MOD-1

    ; Read data
    A0 = read_data(info[imod])
    A0&time = pctco&year(10:)

    ; calculate annual mean for VAR:
    ymin = 0
    ymax = toint(dimsizes(A0)-1)
    tmp0 = (A0 * 3600. * 24. * 365) / 1e12
    copy_VarMeta(A0, tmp0)
    aY0data(imod, ymin:ymax) = (/tmp0/)

    rc = regline_stats(aY0data(imod, ymin:ymax), \
                       pctco({year|nyMin:nyMax}))
    Yfitda(imod, ymin:ymax) = rc@Yest

    rcgpp(imod) = tofloat(rc)

    if (xMin_year.ge.info[imod]@start_year) then
      delta0(imod) = (avg(tofloat(tmp0({time|xMax_year-4:xMax_year})))/ \
                      avg(tofloat(tmp0({time|xMin_year:xMin_year+4}))))

      delta0_std(imod) = sqrt( \
                         (stddev(tofloat(tmp0( \
                             {time|xMax_year-4:xMax_year})))/ \
                          avg(tofloat(tmp0( \
                              {time|xMin_year:xMin_year+4})))) ^ 2 + \
                         (stddev(tofloat(tmp0( \
                             {time|xMin_year:xMin_year+4}))) * \
                          avg(tofloat(tmp0( \
                              {time|xMax_year-4:xMax_year})))/ \
                          avg(tofloat(tmp0( \
                              {time|xMin_year:xMin_year+4}))) ^ 2) ^ 2)
    else
      delta0(imod) = (avg(tmp0({time|xMax_year-4:xMax_year}))/ \
                      avg(tmp0(0:4)))

      delta0_std(imod) = sqrt( \
        (stddev(tmp0({time|xMax_year-4:xMax_year}))/ \
          avg(tmp0(0:4))) ^ 2 + (stddev(tmp0(0:4)) * \
                                 avg(tmp0({time|xMax_year-4:xMax_year}))/ \
                                 avg(tmp0(0:4)) ^ 2) ^ 2)

    end if
    betaGPP(imod) = delta0(imod) * 1/0.9

    delete([/tmp0, ymin, ymax, A0, rc/])
  end do     ; model loop

  ; ---------------------------------------------------------------------
  ; Output to NetCDF
  new_path = config_user_info@work_dir + "/"
  new_path = new_path + "beta_" + xMax_year + "-" + xMin_year + ".nc"

  ; Attach attributes to the results
  betavar = (/betaGPP, delta0_std/)
  betavar!0        = "case"
  betavar&case     = (/"mean", "stddev"/)
  betavar!1        = "model"
  betavar&model    = datasetnames
  betavar@ncdf     = new_path
  betavar@var      = "beta"
  betavar@diag_script = (/DIAG_SCRIPT/)
  betavar@rcgpp    = rcgpp
  betavar@betaGPP  = betaGPP
  betavar@delta0_std = delta0_std

  ; Write NetCDF output
  ncdf_outfile = ncdf_write(betavar, new_path)

  ; -----------------------------------------------------------
  ; creat plots

  ; PLOT - Panel A
  XStg    = "Annual CO~B~2~N~ [ppmv]"
  YStg    = "Annual GPP [GtC/yr]"
  data_arr = (/pctco({year|nyMin:nyMax}), aY0data(0, :)/)
  data_arr!0 = "datasets"
  data_arr&datasets = (/XStg, YStg/)
  data_arr@res_tiMainString    = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString = ""
  data_arr@res_gsnMaximize     = True
  data_arr@res_xyMonoMarker    = False
  data_arr@res_xyMarkers       = markers
  data_arr@res_xyMarkerColors  = "white"
  data_arr@res_trXMinF         = min(pctco({year|nyMin:nyMax}))
  data_arr@res_trXMaxF         = max(pctco({year|nyMin:nyMax}))
  data_arr@res_trYMinF         = 0
  data_arr@res_trYMaxF         = max(aY0data) + 0.5*max(aY0data)
  data_arr@res_vpWidthF        = 0.5
  data_arr@res_vpHeightF       = 0.4

  data_arr@diag_script = DIAG_SCRIPT
  diag_script_info@scatter_log = False

  plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

  txres = True
  MarkerSizeF = where(datasetnames.eq."CESM1-BGC", 20., 17.)
  txres@gsMarkerThicknessF = 3

  lineres = True
  lineres@gsLineDashPattern = 0
  lineres@gsLineThicknessF  = 3

  ; Add year strings as markers
  do imod  = 0, dim_MOD-1
    marker = unique_string("marker")
    txres@gsMarkerColor = colors(imod)
    txres@gsMarkerSizeF = MarkerSizeF(imod)
    txres@gsMarkerIndex = markers(imod)
    add_markers(wks, plot, txres, pctco({year|nyMin:nyMax}), aY0data(imod, :))

    lineres@gsLineColor = colors(imod)
    plot@$marker$ = gsn_add_polyline(wks, plot, Yfitda(imod, :), \
                                     aY0data(imod, :), lineres)
  end do
  draw(plot)
  frame(wks)
  delete([/plot, XStg, YStg, data_arr, marker/])

  ; PLOT - Panel B
  print(datasetnames + " " + betaGPP)

  YStg = str_upper(var) + "(2xCO~B~2~N~)/" + str_upper(var) + "(1xCO~B~2~N~)"

  data_arr = new((/dimsizes(betaGPP), 1/), float)
  data_arr(:, 0) = betaGPP
  data_arr!0 = "datasets"
  data_arr&datasets(0) = datasetnames
  data_arr!1 = "yaxis"
  data_arr&yaxis = YStg
  data_arr@res_tiMainString    = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString = ""

  plot = barchart(wks, data_arr, var, info)
  draw(plot)
  frame(wks)

  delete([/betaGPP, plot, delta0, plot, Yfitda/])

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
  log_provenance(ncdf_outfile, \
                 plotname, \
                 "Climate models vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"global"/),\
                 (/""/), \
                 (/"zechlau_sabrina"/), \
                 (/"wenzel16nat"/), \
                 metadata_att_as_array(info, "filename"))

end
