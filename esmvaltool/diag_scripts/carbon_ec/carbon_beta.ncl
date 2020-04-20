; #######################################################################
; carbon_beta.ncl
; Author: Sabrina Zechlau (DLR, Germany)
; ####################################################################### 
; Description:
; creates two panels to diagnos beta GPP at 2x[CO2] increase
; Panel: a) correlation between GPP and [co2] increase,
;        b) barchart of the trend of the upper correlation
;
; Required info attributes (diag_script_info):
; - styleset       project for line styles
; - bc_xmax_year   start year for time series plots
; - bc_xmin_year   end year for time series plots
;
; Modification history
;    * 20200406-A_wenz_sa: code rewritten for ESMValTool v2.0
;    * 2015xxxx-A_wenz_sa: written
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

  ;******************************************************************  
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

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "cl_output", False)
 
  ; optional input parameters 
  if(isatt(diag_script_info,"bc_xmax_year")) then 
    xMax_year = toint(diag_script_info@bc_xmax_year)
  else
    xMax_year = max(metadata_att_as_array(info, "end_year"))
  end if 
  if(isatt(diag_script_info,"bc_xmin_year")) then 
    xMin_year = toint(diag_script_info@bc_xmin_year)
  else
    xMin_year = max(metadata_att_as_array(info, "start_year"))
  end if 

  if(isatt(diag_script_info,"cl_output")) then 
    Output = diag_script_info@cl_output
  end if 

  ; Call plot scripts
  plot_file = "beta_" + var + "_" + xMin_year + "-" + xMax_year
  wks = gsn_open_wks(file_type, plot_dir + plot_file)
  colors  = project_style(info, diag_script_info, "colors")
  markers = project_style(info, diag_script_info, "markers")
  thicks  = project_style(info, diag_script_info, "avgstd")
  lgLabels = datasetnames

  ; ================================================================== 
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
  ; -----------------------------------------------------------------
  ; define arrey
  betaGPP = new((/dim_MOD/), float)
  delta0  = new((/dim_MOD/), float)
  delta0_std = new((/dim_MOD/), float)
  rcgpp   = new((/dim_MOD/), float)
  aY0data = new((/dim_MOD,nMax/), double)
  Yfitda  = new((/dim_MOD,nMax/), double)

  ; loop for models
  do imod = 0, dim_MOD-1  
    
    ; Read data
    A0 = read_data(info[imod])

    ; use common time units
    utc_date = cd_calendar(A0&time, 0)
    utc_year = tointeger(utc_date(:, 0))
    if (max(tointeger(utc_date(:, 0))).lt.200) then

      utc_year = tointeger(utc_date(:, 0)) + 1849
      month  = tointeger(utc_date(:, 1))
      day    = tointeger(utc_date(:, 2))
      hour   = tointeger(utc_date(:, 3))
      minute = tointeger(utc_date(:, 4))
      second = utc_date(:, 5)

      if (A0&time@units.eq."days since 0000-01-01 00:00:00") then
        units ="days since 1849-01-01 00:00:00"
      else
        units  = "days since 1850-01-01 00:00:00"
      end if

      A0&time = cd_inv_calendar(utc_year, month, day, hour, \
                                minute, second, units, 0)
      delete([/month, day, hour, minute, second/])
    end if
    delete(utc_date)

    ; calculate annual mean for VAR0 and VAR1:
    tmp0 = time_operations(A0, -1, -1, "average", "yearly", False)
    ymin = 0
    ymax = toint(dimsizes(tmp0)-1)
    tmp0 = (tmp0*3600.*24.*365)/1e12
    aY0data(imod, ymin:ymax) = (/tmp0/)

    ; calculate amplitude of gpp
    ti0 := cd_calendar(A0&time, 0) 
    ttime := ti0(:, 0)+(ti0(:, 1)-1)/12.
    seas_cyc := seas_cycle(ttime, A0, "amplitude")

    ; cange units of VAR1 to GtC/yr
    seas_cyc@Yamp = (seas_cyc@Yamp*3600.*24.*365)/1e12

    rc = regline_stats(aY0data(imod, ymin:ymax), \
                       pctco({year|min(utc_year): \ 
                                   max(utc_year)}))
    Yfitda(imod, ymin:ymax) = rc@Yest

    rcgpp(imod) = tofloat(rc)
    delete(seas_cyc)

     if (xMin_year.ge.info[imod]@start_year) then
       delta0(imod) = (avg(tofloat(tmp0({year|xMax_year-4:xMax_year})))/ \
                       avg(tofloat(tmp0({year|xMin_year:xMin_year+4}))))
       delta0_std(imod) = sqrt((stddev(tofloat(tmp0({year|xMax_year-4:xMax_year})))/ \
                                       avg(tofloat(tmp0({year|xMin_year:xMin_year+4}))))^2+ \
                               (stddev(tofloat(tmp0({year|xMin_year:xMin_year+4})))* \
                                       avg(tofloat(tmp0({year|xMax_year-4:xMax_year})))/ \
                                       avg(tofloat(tmp0({year|xMin_year:xMin_year+4})))^2)^2)
     else
       delta0(imod) = (avg(tmp0({year|xMax_year-4:xMax_year}))/ \
                       avg(tmp0(0:4)))
       delta0_std(imod) = sqrt((stddev(tmp0({year|xMax_year-4:xMax_year}))/ \
                          avg(tmp0(0:4)))^2+(stddev(tmp0(0:4))* \
                          avg(tmp0({year|xMax_year-4:xMax_year}))/ \
                          avg(tmp0(0:4))^2)^2)
     end if

     if (xMax_year.eq."1989".and.datasetnames(imod).ne."GFDL-ESM2M") then
       dco2 = 285*3
     else if (xMax_year.eq."1920".or.datasetnames(imod).eq."GFDL-ESM2M") then
       dco2 = 285
     else
       dco2 = avg(tmp1({year|xMax_year-9}))-avg(tmp1({year|xMin_year+9}))
     end if
     end if

     betaGPP(imod) = delta0(imod)*1/0.9

     delete([/tmp0, ymin, ymax, A0, rc/])
   end do     ;model loop 

; ---------------------------------------------------------------------
  ; Optional output to NetCDF
  if (config_user_info@write_netcdf) then
    new_path = config_user_info@run_dir
    new_path = new_path + "beta_" \
                        + diag_script_info@bc_xmax_year + "-" \
                        + diag_script_info@bc_xmin_year + ".nc"

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

   end if
;-----------------------------------------------------------
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
       add_markers(wks, plot, txres, pctco({year|nyMin:nyMax}), aY0data(imod,:))
       
       lineres@gsLineColor = colors(imod)
       plot@$marker$ = gsn_add_polyline(wks, plot, Yfitda(imod, :), \
                                        aY0data(imod, :), lineres)
   end do
   draw(plot)
   frame(wks)
   delete([/plot, XStg, YStg, data_arr, marker/])

   ; PLOT - Panel B
   print(datasetnames + " " + betaGPP)
   tmpb = betaGPP
   qsort(betaGPP)
   tmpm = datasetnames
   tmp_orig = datasetnames

   do imod = 0, dim_MOD-1
       id         = ind(tmpb.eq.betaGPP(imod))
       tmpm(imod) = datasetnames(id)
       delete(id)
   end do
   datasetnames = tmpm

   MainStg = ""
   XStg = datasetnames
   YStg = str_upper(var) + "(2xCO~B~2~N~)/" + \
          str_upper(var) + "(1xCO~B~2~N~)"
   betaGPP@info = info

   plot = barchart(betaGPP, MainStg, XStg, YStg, wks, "draw")
   info@name = tmp_orig

   delete([/betaGPP, delta0, plot, XStg, Yfitda, tmpm/])

  ; create separate legend
  marker_thicks = datasetnames
  marker_sizes  = datasetnames

  marker_thicks = 3 
  marker_sizes  = 0.02

  leg = True
  leg@txFontQuality = "High"
  leg@txFont        = 25
  leg@txFontHeightF = 0.02  ; used in
  leg@diag_script   = DIAG_SCRIPT
  leg@annots        = datasetnames(::-1)
  leg@colors        = colors(::-1)
  leg@markers       = markers(::-1)
  leg@thicks        = marker_thicks
  leg@sizes         = marker_sizes
  leg@ncols         = 1

  create_legend_lines(leg@annots, leg, plot_dir + DIAG_SCRIPT + "_legend", "markers")

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 config_user_info@plot_dir + "." + file_type, \
                 "Climate models vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"global"/),\
                 (/""/), \
                 (/"zechlau_sabrina"/), \
                 (/"wenzel16nature"/), \
                 metadata_att_as_array(info, "filename"))

end 
;#################################################################################

