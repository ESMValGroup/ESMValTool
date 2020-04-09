; #######################################################################
; carbon_cycle_analysis.ncl
; Author: Sabrina Wenzel (DLR, Germany)
; #######################################################################
;                                                       
; Description:
; creates several panels to diagnos beta from CO2 and compares with
; station data:
; Panel: a) [co2],
;        b) detrendet [co2],
;        c) rel. seasonal amplitude of [co2],
;        d) zero crossing day of seasonal [co2]
;
; Required info attributes (from variable_info):
;   - units
;
; Required info attributes (diag_script_info):
; - styleset       project for line styles
; - station        name of station
; - cl_lonrange    array of p-levels for time series plots if 2D or 3D
; - cl_latrange    array of latitudes for time series plots if 2D or 3D
; - cl_levrange    array of pressure levels for time series plots if 3D
; - bc_xmax_year   start year for time series plots
; - bc_xmin_year   end year for time series plots
; - Obs            (optional) reads OBS and compares with models

;Note: cd_calendar requires NCL v6.0
; A) This script is based on the tsline.ncl
; B) The station data file must be of type T0M

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
;load "$diag_scripts/shared/plot/carbon_plots.ncl"

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
  exit_if_missing_atts(diag_script_info, "bc_latmin")
  exit_if_missing_atts(diag_script_info, "bc_latmax")
  exit_if_missing_atts(diag_script_info, "bc_lonmin")
  exit_if_missing_atts(diag_script_info, "bc_lonmax")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "cl_output", False)
  set_default_att(diag_script_info, "station", 1)
    
  if(isatt(diag_script_info,"bc_latmin").and.\
     isatt(diag_script_info,"bc_latmax")) then 
    latrange = (/diag_script_info@bc_latmin, \
                 diag_script_info@bc_latmax/)
  end if 
  if(isatt(diag_script_info,"bc_lonmin").and.\
     isatt(diag_script_info,"bc_lonmax")) then 
    lonrange = (/diag_script_info@bc_lonmin, \
                 diag_script_info@bc_lonmax/)
  end if 

  if(isatt(diag_script_info,"con_units")) then 
    UNITS = variable_info@units 
  end if
 
   ; optional input parameters 
  if(isatt(diag_script_info,"bc_xmax_year")) then 
    xMax_year = toint(diag_script_info@bc_xmax_year)
  end if 
  if(isatt(diag_script_info,"bc_xmin_year")) then 
    xMin_year = toint(diag_script_info@bc_xmin_year)
  end if 

  if(isatt(diag_script_info,"cl_output")) then 
    Output = diag_script_info@cl_output
  end if 

  ; Call plot scripts
  plot_file = "beta_"+var+"_"+xMin_year+"-"+xMax_year
  wks = gsn_open_wks(file_type, plot_dir + plot_file)
  colors  = project_style(info, diag_script_info, "colors")
  markers = project_style(info, diag_script_info, "markers")
  thicks  = project_style(info, diag_script_info, "avgstd")
  lgLabels = datasetnames

  ;;================================================================== 
  ;; setup data holders for plots 
  nyMax = max(metadata_att_as_array(info, "end_year"))
  nyMin = min(metadata_att_as_array(info, "start_year"))
  nMax = nyMax - nyMin + 1

  pctco = new((/nMax+10/), double) 
  pctco(0) = 285
  do tt = 1, nMax+10-1
    pctco(tt) = pctco(tt-1) + 0.01 * pctco(tt-1)
  end do
  pctco!0 = "year"
  pctco&year = ispan(nyMin-10, nyMax, 1)
;------------------------------------------------------------------
  ; define arrey
    betaGPP = new((/dim_MOD/), float)
    betaNBP = new((/dim_MOD/), float)
    delta0  = new((/dim_MOD/), float)
    delta0_std = new((/dim_MOD/), float)
    delta1  = new((/dim_MOD/), float)
    rcgpp   = new((/dim_MOD/), float)
    aXdata  = new((/dim_MOD,nMax/), double)
    aY0data = new((/dim_MOD,nMax/), double)
    aY1data = new((/dim_MOD,nMax/), double)
    aY2data = new((/dim_MOD,nMax/), double)
    Yfitd   = new((/dim_MOD,nMax/), double)
    Yfitda  = new((/dim_MOD,nMax/), double)

  ;;loop for models
  do imod = 0, dim_MOD-1  
    
    ; Read data
    A0 = read_data(info[imod])

    ; x array (one value per year) 
    tmp = ispan(nyMin, nyMax,1) * 1.
    nyr = dimsizes(tmp) 
    delete(tmp) 
    tmp = ispan(1,nyr,1)*1.
    aXdata(imod, 0:nyr(0) - 1) = tmp
    delete(tmp)

    ; calculate area sum for VAR0 and VAR1:
    tmpa0 = area_operations(A0, latrange(0), latrange(1), \
                            lonrange(0), lonrange(1), "sum", False)

    ; calculate annual mean for VAR0 and VAR1:
    tmp0 = time_operations(tmpa0, -1, -1, "average", "yearly", False)
    ymin = 0
    ymax = toint(dimsizes(tmp0)-1)
    tmp0 = (tmp0*3600.*24.*365)/1e12
    aY0data(imod, ymin:ymax) = (/tmp0/)

    ; calculate amplitude of gpp
    ti0 := cd_calendar(A0&time, 0) 
    ttime := ti0(:, 0)+(ti0(:, 1)-1)/12.
    seas_cyc := seas_cycle(ttime, tmpa0, "amplitude")

    ; cange units of VAR1 to GtC/yr
    seas_cyc@Yamp = (seas_cyc@Yamp*3600.*24.*365)/1e12

    rc = reg_multlin_stats(pctco({year|toint(info@start_year(imod)): \
                                       toint(info@end_year(imod))}), \
                                       aY0data(imod, ymin:ymax), False)
    Yfitd(imod,ymin:ymax) = rc@Yest

    idx := ind(.not.ismissing(seas_cyc@Yamp))
    aY2data(imod, :dimsizes(idx)-1) = (/seas_cyc@Yamp(idx)/)
    delete(rc)

    rc = regline_stats(aY0data(imod, ymin:ymax), \
                       pctco({year|toint(info@start_year(imod)): \
                                   toint(info@end_year(imod))}))
    Yfitda(imod, ymin:ymax) = rc@Yest

    rcgpp(imod) = tofloat(rc)
    delete([/seas_cyc, idx, tmpa0/])

     if (xMin_year.ge.info@start_year(imod)) then
       delta0(imod) = (avg(tmp0({year|xMax_year-4:xMax_year}))/ \
                       avg(tmp0({year|xMin_year:xMin_year+4})))
       delta0_std(imod) = sqrt((stddev(tmp0({year|xMax_year-4:xMax_year}))/ \
                                       avg(tmp0({year|xMin_year:xMin_year+4})))^2+ \
                                       (stddev(tmp0({year|xMin_year:xMin_year+4}))* \
                                       avg(tmp0({year|xMax_year-4:xMax_year}))/ \
                                       avg(tmp0({year|xMin_year:xMin_year+4}))^2)^2)
     else
       delta0(imod) = (avg(tmp0({year|xMax_year-4:xMax_year}))/ \
                       avg(tmp0(0:4)))
       delta0_std(imod) = sqrt((stddev(tmp0({year|xMax_year-4:xMax_year}))/ \
                          avg(tmp0(0:4)))^2+(stddev(tmp0(0:4))* \
                          avg(tmp0({year|xMax_year-4:xMax_year}))/ \
                          avg(tmp0(0:4))^2)^2)
     end if

     if (xMax_year.eq."1989".and.info@name(imod).ne."GFDL-ESM2M") then
       dco2 = 285*3
     else if (xMax_year.eq."1920".or.info@name(imod).eq."GFDL-ESM2M") then
       dco2 = 285
     else
       dco2 = avg(tmp1({year|xMax_year-9}))-avg(tmp1({year|xMin_year+9}))
     end if
     end if

     betaGPP(imod) = delta0(imod)*1/0.9

     delete(tmp0)
     delete(ymin)
     delete(ymax)
     delete(A0)
     delete(rc)
   end do     ;model loop 

; ---------------------------------------------------------------------
; write arrays into table:
; ------------------------
   f_txt = config_user_info@plot_dir+"/" + diag_script + \
   "_fractional"+str_upper(var) + "_" + \
   diag_script_info@station + ".txt"
     
   header = datasetnames(0)
   alist = aY0data(0,:)
   do modi = 1, dimsizes(datasetnames)-1
     system("echo " + datasetnames(modi) + ">>" + header)
     system("echo " + aY0data(modi,:) + ">>" + alist)
   end do

   hlist = [/header/]
   write_table(f_txt, "w", hlist, "%s ")
   write_table(f_txt, "a", alist, "%i %f %f %f %f %f %f %f")

   delete(f_txt)
   delete(hlist)
   delete(alist)
;------------------------------------------------------------------------------
   ;;Optional output to NetCDF
  if (getenv("ESMValTool_write_netcdf").eq."True") then

      f_out1 = plot_file

      ;;Attach attributes to the results
      var = (/betaGPP,delta0_std/)
      var!0        = "case"
      var&case     = (/"mean", "stddev"/)
      var!1        = "model"
      var&model    = info@name
      var@ncdf     = f_out1
      var@var      = "beta"
      var@rcgpp    = rcgpp     
      var@diag_script = (/diag_script/)

      ;;Write NetCDF output
      ncdf_outfile = ncdf_write(var, f_out1)
      delete([/f_out1, ncdf_outfile, var/])

   end if
;-----------------------------------------------------------
  ;;creat plots

   if (latrange(0).eq.-90.and.latrange(1).eq.90) then
     latname = "Global"
   else
     latname = latrange(0) + "~S~o~N~N-" + \
               latrange(1) + "~S~o~N~N"
   end if

   ;;PLOT - Panel A
   XStg    = "Annual CO~B~2~N~ [ppmv]"
   YStg    = "Annual GPP [GtC/yr], " + latname
   data_arr = (/pctco({year|nyMin:nyMax}), aY0data(0, :)/)
   data_arr!0 = "models"
   data_arr&info = (/XStg, YStg/)
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

   data_arr@diag_script = diag_script
   diag_script_info@scatter_log = False

   plot = scatterplot(wks, data_arr, var, False, False, input_file_info)

   txres = True
   MarkerSizeF = where(datasetnames.eq."CESM1-BGC", 20., 17.)
   txres@gsMarkerThicknessF = 3

   lineres = True
   lineres@gsLineDashPattern = 0
   lineres@gsLineThicknessF  = 3

   ;; Add year strings as markers
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

   ;;PLOT - Panel B
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
          str_upper(var) + "(1xCO~B~2~N~), " + latname

   plot = barchart(betaGPP, MainStg, XStg, YStg, wks, "draw")
   info@name = tmp_orig

   delete([/aXdata, aY1data, betaGPP, delta0, plot, XStg, Yfitd, tmpm/])

  ;;create separate legend
  marker_thicks = datasetnames
  marker_sizes  = datasetnames

  marker_thicks = 3 
  marker_sizes  = 0.02

  leg = True
  leg@txFontQuality = "High"
  leg@txFont        = 25
  leg@txFontHeightF = 0.02  ; used in
  leg@diag_script   = diag_script
  leg@annots        = datasetnames(::-1)
  leg@colors        = colors(::-1)
  leg@markers       = markers(::-1)
  leg@thicks        = marker_thicks
  leg@sizes         = marker_sizes
  leg@ncols         = 1

  create_legend_lines(leg@annots, leg, plot_dir + diag_script + "_legend", "markers")

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 output_dir + "." + file_type, \
                 XStg + " vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"trop", "global"/),\
                 (/""/), \
                 (/"zechlau_sabrina"/), \
                 (/"wenzel16nature"/), \
                 metadata_att_as_array(info, "filename"))

end 
;#################################################################################

