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
 
  if(isatt(diag_script_info,"cl_latrange")) then 
    latrange = diag_script_info@cl_latrange 
  end if
   if(isatt(diag_script_info,"cl_lonrange")) then 
    lonrange = diag_script_info@cl_lonrange 
  else 
    info_output("plot requires variable attribute 'cl_lonrange'",verbosity,1) 
    exit 
  end if 

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "cl_output", False)
 
   ; optional input parameters 
  if(isatt(diag_script_info,"cl_xmax_year")) then 
    xMax_year = diag_script_info@cl_xmax_year
  else
    xMax_year = max(metadata_att_as_array(info, "end_year")) 
  end if 
  if(isatt(diag_script_info,"cl_xmin_year")) then 
    xMin_year = diag_script_info@cl_xmin_year 
  else
    xMin_year = max(metadata_att_as_array(info, "start_year"))
  end if 
  if(isatt(diag_script_info,"cl_period")) then 
    period = diag_script_info@cl_period
  end if
  if(isatt(diag_script_info,"cycle")) then 
    cycle = diag_script_info@cycle
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

   Ydata  = new((/dim_MOD+1, nMax/), float) 
   newYdata = new((/dim_MOD+1, nMax/), float)
   Yfitd = new((/dim_MOD+1, nMax/), float)
   XPdata = new((/dim_MOD+1, nMax/), float)
   A0data = new((/dim_MOD+11, nMax*12/), float) 
   rcc   = new((/dim_MOD+1, 2/), float)
   stdrc = new((/dim_MOD+1, 2/), float)  
   arc   = new((/dim_MOD+1/), float)
   yfit  = new((/dim_MOD+1, nMax-2/), float) 
   A0data!1 = "time"
   A0data&time = fspan(nyMin,nyMax,nMax*12)
   newYdata!1 = "time" 
   newYdata&time = ispan(nyMin+1,nyMax+1,1)
   XPdata!1 = "year" 
   XPdata&year = ispan(nyMin+1,nyMax+1,1)

  ; ##############################################################
  ; loop for models
  do imod = 0, dim_MOD-1  
    
    ; Read data
    A0 = read_data(info[imod])

    A0data(imod,:dimsizes(tmp0)-1) = A0
    tmpd0 = time_operations(tmp0, -1, -1, "average", "yearly", False)

    mod_max = toint(models@end_year(imod)) - toint(models@start_year(imod)) +1
    XPdata(imod, :dimsizes(tmpd0)-2) = (/tmpd0(1:)/)

    ti0 = cd_calendar(tmp0&time, 0) 
    tt  = ti0(:, 0)+(ti0(:, 1)-1)/12.

    seas_cyc = seas_cycle(tt, tmp0, cycle)
    tmp = seas_cyc@$cycle$
    idx = ind(.not.ismissing(seas_cyc@$cycle$))
    tmpA2 = tmp(idx)
    idx3 = ind(.not.ismissing(tmpA2))
    Ydata(imod, :dimsizes(idx3)-1) = (/tmpA2(idx3)/)
    newYdata(imod, :dimsizes(idx)-1) = (/tmp(idx)/)

    delete([/seas_cyc, tmp, tt, ti0, idx/])

    rc = regline_stats(XPdata(imod,:),Ydata(imod,:))
    rcc(imod,1) = rc
    rcc(imod,0) = rc@b(0)
    stdrc(imod,:) = rc@stderr
    arc(imod) = rc@r
    Yfitd(imod,:dimsizes(rc@Yest)-1) = rc@Yest

    delete([/rc, tmpd0, tmpA2, tmp0, A0, idx3, mod_max/])
   end do     ;model loop

   betaGPP = new((/2,dim_MOD+1/),float)
   rcGPP   = new((/dim_MOD*2/),float)

   fils = systemfunc("ls "+plot_dir+"/climo/CMIP5_beta_gpp_60N-90N_1860-1920.nc")
   tmp = addfile(fils, "r")
   betaGPP(:,:dim_MOD-1) = tmp->beta

   ;;asign gradients
   rcGPP(::2)  = (/0.02249557,0.01586392,0.01731257,0.004872571,0.02033394,0.006016535,0.005555378/)
   tmpg = betaGPP@rcgpp
   rcGPP(1::2) = tmpg
   delete(tmpg)

   rcGPP!0 = "models"
   tmpm = new((/dim_MOD*2/),string)
   tmpm(::2) = models@name(:dim_MOD-1)
   tmpm(1::2) = models@name(:dim_MOD-1);"";
   rcGPP&models = tmpm

   delete(fils)
   delete(tmp)
   delete(tmpm)
   nfit = fspan(0,50,50)

   ;;Calculate constraint on BETA_GPP
   ;;amplitude trends
   minx = min(rcc(:dim_MOD-1,1))-1.1*(max(rcc(:dim_MOD-1,1))-min(rcc(:dim_MOD-1,1)))
   maxx = max(rcc(:dim_MOD-1,1))+1.1*(max(rcc(:dim_MOD-1,1))-min(rcc(:dim_MOD-1,1)))

   rc = regline_stats(rcc(:dim_MOD-1,1),betaGPP(0,:dim_MOD-1))
   arc_bco2  = rc@r
   prc_bco2  = rc@F_pval
   xfit_bco2 = minx+(maxx-minx)*nfit/dimsizes(nfit)
   yfit_bco2 = rc*xfit_bco2 + rc@yintercept
   beta_mco2 = rc*rcc(dim_MOD,1)+rc@yintercept
   beta_eco2 = sqrt(stdrc(dim_MOD,1)^2+rc@rstd^2);beta_mco2*(stdrc(dim_MOD,1)/rcc(dim_MOD,1))

   ;;Calculating the prediction error
   D    = new((/dim_MOD,2/),float)   
   xobs = new((/1,2/),float)
   I    = new((/dim_MOD,dim_MOD/),float)
   I      = 1.
   D      = 1
   D(:,1) = rcc(:dim_MOD-1,1) 
   xobs   = 1.
   xobs(0,1) = rcc(dim_MOD,1)
   df = dim_MOD-1

   cfint  = cdft_t(0.05/2., df)^2*(1.+xobs#inverse_matrix(transpose(D)#D)#transpose(xobs))
   SSE1   = betaGPP(0,:dim_MOD-1)#(I-D#inverse_matrix(transpose(D)#D)#transpose(D))
   sigma  = (1./df)*SSE1
   tmp_cf = cfint(0,0)*sigma
   cf     = (sqrt(tmp_cf#betaGPP(0,:dim_MOD-1)))/100 ;;prediction error

   info_output("beta gpp-co2 "+beta_mco2+" +/- "+cf, verbosity,1)
   info_output("mean gpp-co2 "+avg(betaGPP(0,:dim_MOD-1))+" +/- "+stddev(betaGPP(0,:dim_MOD-1)), verbosity,1)

   PDF = Cond_PDF(xfit_bco2,yfit_bco2,rcc(:dim_MOD-1,1),betaGPP(0,:dim_MOD-1),\
                               rc@MSE,(/rcc(dim_MOD,1),stdrc(dim_MOD,1)/))
 
   co2_Py    = PDF@Py
   co2_Py_pr = PDF@Py_pr
   co2_y     = PDF@y

   ;;Calculate prior distribution
   mn_pr  = avg(betaGPP(0,:dim_MOD-1))
   std_pr = stddev(betaGPP(0,:dim_MOD-1))

   Py_pr := (1/(std_pr*sqrt(2*PI))) * exp(-0.5*(((betaGPP(0,:dim_MOD-1)-mn_pr)/std_pr)^2))

   ;;Calculate equal distand Histogramm
   bin_Py  = new(3,float)
   bin_Py@bin_center = fspan(min(betaGPP(0,:)),max(betaGPP(0,:)),3)
   mn_norm = dim_MOD/mn_pr/dim_MOD

   if (diag_script_info@station.eq."BRW") then
       bin_Py(0) = 4*mn_norm
       bin_Py(1) = 2*mn_norm
       bin_Py(2) = 1*mn_norm
   else if (any(diag_script_info@station.eq.(/"MLO","KMK"/))) then
       bin_Py(0) = 2*mn_norm
       bin_Py(1) = 4*mn_norm
       bin_Py(2) = 2*mn_norm
   end if
   end if

   ;;this is for plotting the regression line with confidence intervalls
   yfit_bco2_std = new((/2,dimsizes(yfit_bco2)/),typeof(yfit_bco2))
   yfit_bco2_std(0,:) = yfit_bco2+PDF@sigma
   yfit_bco2_std(1,:) = yfit_bco2-PDF@sigma

   yBETA = new((/2,dim_MOD+2/),float)
   yBETA(0,0)  = avg(betaGPP(0,:dim_MOD-1))
   yBETA(0,1)  = beta_mco2
   yBETA(0,2:) = betaGPP(0,:dim_MOD-1)
   yBETA(1,0)  = stddev(betaGPP(0,:dim_MOD-1))
   yBETA(1,1)  = beta_eco2
   yBETA(1,2:) = betaGPP(1,:dim_MOD-1)
   
   xCO2 =  new((/2,dim_MOD+2/),float)
   xCO2(0,1)  = rcc(dim_MOD,1)
   xCO2(0,2:) = rcc(:dim_MOD-1,1)
   xCO2(1,1)  = stdrc(dim_MOD,1)
   xCO2(1,2:) = stdrc(:dim_MOD-1,1)
   
   fit_reg = (/xfit_bco2,yfit_bco2/)

   delete(rc)
   delete(PDF)
   delete(minx)
   delete(maxx)
   info_output("r = "+arc_bco2+"; P = "+prc_bco2, verbosity,1)
;-----------------------------------------------------------
   colors = project_style(diag_script_info, "colors")
   thicks = project_style(diag_script_info, "markers")
   plot_num = integertochar(ispan(97,122,1))
   lgLabels = models@name
   marker   = new(dim_MOD+1,"graphic")

   if (diag_script_info@cl_latrange_t(0).eq.-90.and.diag_script_info@cl_latrange_t(1).eq.90) then
     latname = "Global"
   else
     latname = diag_script_info@cl_latrange_t(0)+"~S~o~N~N-"+diag_script_info@cl_latrange_t(1)+"~S~o~N~N"
   end if   

   ;;create separate legend
   marker_thicks = datasetnames
   marker_sizes  = datasetnames
 
   marker_thicks = 3 
   marker_sizes  = 0.02
 
   leg = True
   leg@txFontQuality = "High"
   leg@txFont        = 25
   leg@txFontHeightF = 0.02
   leg@diag_script   = DIAG_SCRIPT
   leg@annots        = datasetnames
   leg@colors        = colors
   leg@markers       = thicks
   leg@thicks        = marker_thicks
   leg@sizes         = marker_sizes
   leg@ncols         = 1

   create_legend_lines(leg@annots, leg, diri_out + DIAG_SCRIPT + \
                       "_legend", "markers")

   ;;Call plot scripts
   plot_file  = var0 + "-" + cycle + "_" + \
                diag_script_info@station + "_" + nyMin + "-" + nyMax
   wks = gsn_open_wks(file_type, output_dir + plot_file)

   ; Create Plots 1:
   ; -----------------------------------------------
   ;Loop over models
   MainStg = ""
   XStg    = "Annual CO~B~2~N~ [ppmv]"
   YStg    = "CO~B~2~N~ Amplitude [ppmv]"

   data_arr = (/XPdata(0, :),Ydata(0, :)/)
   data_arr!0 = "models"
   data_arr&models = (/XStg, YStg/)
   data_arr@res_tiMainString      = MainStg
   data_arr@res_tiMainFontHeightF = 0.025
   data_arr@res_gsnCenterString   = ""
   data_arr@res_gsnMaximize       = True
   data_arr@res_xyMonoMarker      = False
   data_arr@res_xyMarkers         = thicks 
   data_arr@res_xyMarkerSizeF     = 10.
   data_arr@res_trXMinF           = min(XPdata(:, :))
   data_arr@res_trXMaxF           = max(XPdata(:, :))
   data_arr@res_trYMinF           = 0;min(Ydata(:, :))
   data_arr@res_trYMaxF           = max(Ydata(:, :)) + 0.25 * max(Ydata(:, :))
   data_arr@res_vpWidthF          = 0.5
   data_arr@res_vpHeightF         = 0.4

   data_arr@diag_script = DIAG_SCRIPT
   diag_script_info@scatter_log = False

   plot = scatterplot(wks, data_arr, var0, False, False)

   txres = True
   txres@gsMarkerSizeF      = 10.
   txres@gsMarkerThicknessF = 3

   lineres = True
   do imod  = 0, dim_MOD
  ;; Add year strings as markers
       idx := ind(.not.ismissing(XPdata(imod, :)))
       idy := ind(.not.ismissing(Ydata(imod, :)))
       if (models@name(imod).eq."OBS") then
           lineres@gsLineThicknessF  = 4
       else
           lineres@gsLineThicknessF  = 3
       end if
       txres@gsMarkerColor      = colors(imod)
       txres@gsMarkerIndex      = thicks(imod)
       add_markers(wks, plot, txres, XPdata(imod, idx), Ydata(imod, idy))
       
       lineres@gsLineDashPattern = 0
       lineres@gsLineColor = colors(imod)
       marker(imod) = gsn_add_polyline(wks, plot, XPdata(imod, :), \
                                       Yfitd(imod, :), lineres)
   end do

   draw(plot)
   frame(wks)
   delete([/plot, XStg/])

   tmp_orig = models@name
   tmpb = rcc(:dim_MOD-1,1)
   tmpb_orig = rcc(:,1)
   
   rcc(:dim_MOD-1, 1) = tmpb(::-1)
   tmpm = models@name(:dim_MOD-1)
   models@name(:dim_MOD-1) = tmpm(::-1)

   MainStg = ""
   XStg    = models@name
   YStg    = "Sensitivity of CO~B~2~N~ Amplitude to CO~B~2~N~"
   plot = barchart(rcc(:, 1), MainStg, XStg, YStg, wks, "draw")

   models@name = tmp_orig 
   delete([/plot, XStg, tmpm, tmpb, tmp_orig/])

   ; Set resources
   res = True
   res@trXMinF  = 1974
   res@trXMaxF  = 2005+0.25 * (2005-1974)
   res@tmXBMode = "Explicit"
   res@trXMinF  = nyMin
   res@trXMaxF  = nyMax + 0.25 * (nyMax - nyMin)
   res@tiYAxisString = "CO~B~2~N~ [ppmv]"
   xy_line(wks, A0data, A0data&time, 0, res, False)

   tmpm = models@name(:dim_MOD-1)
   delete(models@name)
   models@name = new(dim_MOD+2,typeof(tmpm))
   models@name(0) = "OBS"
   models@name(1) = "constraint mean"
   models@name(2:) = tmpm
   delete(tmpm)

   ; Create Plots 2:
   ; -----------------------------------------------
   MainStg = ""
   XStg    = "Sensitivity of CO~B~2~N~ Amplitude to CO~B~2~N~"
   YStg    = "GPP(2xCO~B~2~N~)/GPP(1xCO~B~2~N~)"

   plot = mder_scatter_plot(xCO2(0,:),xCO2(1,:),yBETA(0,:),yBETA(1,:),\
                            fit_reg,yfit_bco2_std,\
                            MainStg,XStg,YStg,wks,"draw")
   delete(plot)

   bres  = True
   bres@gsnDraw          = False
   bres@gsnFrame         = False
   bres@xyLineThicknessF = 2.        
   bres@xyLineColor      = (/"black"/)  
   bres@tiXAxisString    = "GPP(2xCO~B~2~N~)/GPP(1xCO~B~2~N~)"
   bres@tiYAxisString    = "Probability Density"
   bres@trXMinF          = 0.6;min(co2_y)
   bres@trXMaxF          = 2.0;max(co2_y)
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
  
   plot = gsn_csm_xy(wks, bin_Py@bin_center, bin_Py, bres)

   obsres = True
   obsres@gsLineColor = "red"
   obsres@gsLineThicknessF  = 2
   PDF3 = gsn_add_polyline(wks,plot,co2_y,co2_Py,obsres)
   
   draw(plot)
   frame(wks)
   delete([/plot/])

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 config_user_info@plot_dir + "." + file_type, \
                 XStg + " vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"global"/),\
                 (/""/), \
                 (/"zechlau_sabrina"/), \
                 (/"wenzel16nature"/), \
                 metadata_att_as_array(info, "filename"))
 
end 
; ##############################################################################
