; #############################################################################
; carbon_constraint.ncl
; Author: Sabrina Wenzel (DLR, Germany)    
; #############################################################################
;                                                                            
; Description:
;    Calculates the long term sensitivity GAMMA_LT according to Wenzel et al. 
;    2013 with observations from GCP.
;    This plotscript calls a file wich contains the results from calculating 
;    GAMMA_IAV (carbon_corr_2var.ncl)
;
; Required diag_script_info attributes (diagnostics specific)
;    * con_latrange   array of latitudes for time series plots if 2D or 3D 
;    * con_lonrange   array of pressure levels for time series plots if 2D or 3D
;    * gIAV_file      file path for gamma_IAV values calculated with 
;                     carbon_corr_2var.ncl
;
; Modification history
;    202001-A_wenz_sa: code rewritten for ESMValTool v2.0
;    * 20170822-A_bock_ls: added tags for reporting
;    * 20151027-A_laue_ax: added call to 'write_references'
;    * 201509-A_Sabrina Wenzel: updated to new ESMValTool structure
;    * 201303-A_Sabrina Wenzel: written 
; #############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/carbon_ec/carbon_aux.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  tmp = metadata_att_as_array(input_file_info, "short_name")
  variables = get_unique_values(tmp)
  delete(tmp) 
print(variables)
  ; get indices associated with specific variables
  tas_idx = 0;ind(variables.eq."tas")
  nbp_idx = 1;ind(variables.eq."nbp")
print(input_file_info)
  ; Load var_info
  var0 = variable_info[tas_idx]@short_name;variables(tas_idx)
  info0 = select_metadata_by_name(input_file_info, var0)

  var_1pct = variable_info[nbp_idx]@short_name;variables(nbp_idx)
  info1 = select_metadata_by_name(input_file_info, var_1pct)

  var_esmFix = variable_info[nbp_idx]@short_name;variables(nbp_idx)
  info2 = select_metadata_by_name(input_file_info, var_esmFix)

  log_info("+++++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + " + " + var1 + ")")
  log_info("+++++++++++++++++++++++++++++++++++++++++++++")

  dim_MOD = ListCount(info0)
  datasetnames0 = metadata_att_as_array(info0, "dataset")
  datasetnames1 = metadata_att_as_array(info1, "dataset")
  datasetnames2 = metadata_att_as_array(info2, "dataset")

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
  set_default_att(diag_script_info, "ts_anom", False)
  set_default_att(diag_script_info, "volcanoes", False)
    
  if(isatt(diag_script_info, "ec_anom")) then
    anom = diag_script_info@ec_anom
  end if

  if(isatt(diag_script_info,"con_units")) then 
    UNITS = diag_script_info@con_units 
  end if 

  file_list = 

  ;;================================================================== 
  ;; setup data holders for plots 
  nMax = max(toint(models@end_year)) - min(toint(models@start_year)) + 1
  xMax_year = max(toint(models@end_year))
  xMin_year = min(toint(models@start_year))

  ;;carbon cycle climate feedback
  gamma_LT = new((/dim_MOD/),"float") 
  gamma_LT_std = new((/dim_MOD/),"float") 
  ;;trend lines 
  Yreg= new((/dim_MOD-1/),"float")

  ;;slope 
  arc=new((/dim_MOD/),"float") 
    
  ;;read uncoupled VAR0
  climo_dir = getenv("ESMValTool_climo_dir")
  ;;================================================================== 
  ;;loop for models
  do imod = 0, dim_MOD-1  
    
       ; Read data
    A0 = read_data(info0[imod])
    A1 = read_data(info1[imod])
    A2 = read_data(info2[imod])
    
    idx_mod := imod

    ; Convert units if appropriate
    if (isatt(variable_info[tas_idx], "plot_units")) then
      A0 := convert_units(A0,variable_info[tas_idx]@plot_units)
    end if
    if (isatt(variable_info[nbp_idx], "plot_units")) then
      A1 := convert_units(A1_co2flux,variable_info[nbp_idx]@plot_units)
    end if
    if (isatt(variable_info[nbp_idx], "plot_units")) then
      A2 := convert_units(A2_co2flux,variable_info[nbp_idx]@plot_units)
    end if

    tmpcA0 = new((/dimsizes(tmpA0)+1/),"float")
    tmpcA2 = new((/dimsizes(tmpA2)+1/),"float")
    tmpcA0(0) = 0.0
    tmpcA2(0) = 0.0
    
    do y = 0,nMax-1
      tmpcA0(y+1) = tmpcA0(y)+A0(y)
      tmpcA2(y+1) = tmpcA2(y)+A2(y)
    end do
    delete(y)

    ; calculate integral line:
    tidx := dimsizes(tmpA0)-1
    Y1data := avg(tmpcA0(tidx-9:tidx))-avg(tmpcA0(:9))
    Y2data := avg(tmpcA2(tidx-9:tidx))-avg(tmpcA2(:9))
    Y3data := avg(tmpA1(tidx-9:tidx))-avg(tmpA1(:9))

    ; calculate gamma_LT:
    gamma_LT(imod) = ((Y1data-Y2data)/(Y3data))
    gamma_LT_std(imod) = sqrt((stddev(tmpA0(tidx-9:tidx))/Y3data)^2\
    +(1/stddev(tmpA1(tidx-9:tidx)))^2\
    +(stddev(tmpA2(tidx-9:tidx))/Y3data)^2)
    
    delete([/A0,A1,A2,tidx/])
    delete([/tmpA0,tmpA1,tmpA2/])
    delete([/tmpaA0,tmpaA1,tmpaA2/])
    delete([/Y1data,Y2data,Y3data/])
    delete([/tmpcA0,tmpcA2/])
  end do     ;model loop

   ;;read file with GAMMA_IAV values
   IAV = ncdf_read(file_list, "gIAV")

   if (dimsizes(IAV(0,:)).gt.dim_MOD) then
       Y6data  = IAV(:,1:dim_MOD)         ;;Gamma_IAV from the models
       obs_std = IAV(:,0)                 ;;Gamma_IAV from the OBS
   else
       Y6data = IAV
       obs_std = (/0,0/)
   end if

   ;;exclude models from regression if appropriate
   if(isatt(diag_script_info,"con_units")) then 
     ind_MOD = diag_script_info@reg_models
     aX = where(datasetnames0.eq.diag_script_info@reg_models,Y6data,Y6data@_FillValue)
     aY = where(datasetnames0.eq.diag_script_info@reg_models,gamma_LT,gamma_LT@_FillValue)
   end if

   ;;Calculate confidence limits on fit (see Wikipedia page on "Simple Linear Regression")
   minx = min(aX(0,:))-1.1*(max(aX(0,:))-min(aX(0,:)))
   maxx = max(aX(0,:))+1.1*(max(aX(0,:))-min(aX(0,:)))
   nfit = fspan(0,50,50)

   ;;calculate regression line:
   rc = regline_stats(aX(0,:),aY(:)) 
   arc  = rc@r
   prc  = rc@F_pval
   xfit = minx+(maxx-minx)*nfit/dimsizes(nfit)
   yfit = rc*xfit + rc@yintercept

   PDF = Cond_PDF(xfit,yfit,aX(0,:),aY(:),rc@MSE,obs_std)

   if (obs_std(0).ne.0) then
       Py    = PDF@Py
       gamma_mco2 = rc*obs_std(0)+rc@yintercept
       gamma_eco2 = 0.5*cdft_t(0.05/2., 5)*sqrt(rc@rstd^2+obs_std(1)^2);gamma_mco2*(obs_std(1)/obs_std(0))
       info_output("constraint gamma_LT = "+gamma_mco2+" +/- "+gamma_eco2, verbosity,1)
   end if
   Py_pr = PDF@Py_pr
   y     = PDF@y

   ;;this is for plotting the regression line with confidence intervalls
   yfit_std = new((/2,dimsizes(yfit)/),typeof(yfit))
   yfit_std(0,:) = yfit+PDF@sigma
   yfit_std(1,:) = yfit-PDF@sigma

  ;;================================================================== 
  ;;PLOTTING

   ;; Get environment variables
   plot_dir = getenv("ESMValTool_plot_dir")
   diag_script_base = DIAG_SCRIPT
   output_dir = config_user_info@plot_dir

   plot_file = "constr_"+var0+"-"+var1+"_"+xMin_year+"-"+xMax_year
   wks = gsn_open_wks(file_type, output_dir + plot_file)

   colors = project_style(info0,diag_script_info, "colors")
   thicks = project_style(info0,diag_script_info, "thicks")
   plot_num = integertochar(ispan(97,122,1))

   XStg    = "Sensitivity of land+ocean CO~B~2~N~ flux, ~F8~g~F21~~B~LT~N~ (GtC/yr/K)"
   YStg    = "Carbon Cycle Climate Feedback, ~F8~g~F21~~B~LT~N~ (GtC/K)"

   data_arr = (/Y6data(0,:), gamma_LT/)
   data_arr!0 = "models"
   data_arr&models = (/XStg,YStg/)
   data_arr@diag_script = DIAG_SCRIPT

   data_arr@res_tiMainString      = plot_num(0)+") "
   data_arr@res_tiMainFontHeightF = 0.025
   data_arr@res_gsnCenterString   = ""
   data_arr@res_trXMinF           = min(Y6data(0,:))-0.2*abs(min(Y6data(0,:)))
   data_arr@res_trXMaxF           = max(Y6data(0,:))+0.2*abs(max(Y6data(0,:)))
   data_arr@res_trYMinF           = min(gamma_LT)-0.2*abs(min(gamma_LT))
   data_arr@res_trYMaxF           = max(gamma_LT)+0.2*abs(max(gamma_LT))
   plot = scatterplot(wks, data_arr, var0, False, False)
   
   peres = True
   peres@tfPolyDrawOrder = "PreDraw"
   peres@xyLineColor     = "orange"
   peres@gsFillColor     = "orange"
   xin = array_append_record(xfit,xfit(::-1),0)
   yin = array_append_record(yfit_std(0,:),yfit_std(1,::-1),0)
   add_prediction_error(wks, plot, peres, xin, yin)
   delete([/xin,yin/])

   if (obs_std(0).ne.0) then
       xin = (/obs_std(0)+obs_std(1), obs_std(0)+obs_std(1), \
               obs_std(0)-obs_std(1), obs_std(0)-obs_std(1)/)
       yin = (/100,-150,-150,100/)
       add_prediction_error(wks, plot, False, xin, yin)
       delete([/xin,yin/])   

;;        ecres = True 
;;        ecres@gsLineColor       = "slateblue3"
;;        ecres@gsLineDashPattern = 2
;;        ecres@gsLineThicknessF  = 1.5
;;        ec_line = gsn_add_polyline(wks, plot, (/data_arr@res_trXMinF, obs_std(0), obs_std(0)/), \
;;                                       (/gamma_mco2, gamma_mco2, data_arr@res_trYMinF/), ecres)
   end if

   lres = True
   lres@gsLineColor       = "red"
   lres@gsLineThicknessF  = 1.5
   fit_line = gsn_add_polyline(wks, plot, xfit, yfit, lres)

   ;;add colored markers
   mres = True
   mres@xyMarker      = 0
   do imod = 0, dim_MOD-1
       mres@xyMarkerColor = colors(imod)
       add_markers(wks, plot, mres, Y6data(0,imod), gamma_LT(imod))
   end do

   ;;add xy error bars
   eres = True
   eres@gsLineColor      = colors
   add_errorbar(wks, plot, eres, Y6data, (/gamma_LT,gamma_LT_std/))
   draw (plot) 
   frame (wks) 

   ;;create line plot with PDFs
   diag_script_info@multi_model_mean = "n"
   pres = True
   pres@tiMainString       = plot_num(1)+") "
   pres@tiXAxisString      = XStg
   pres@tiYAxisString      = "Probability Density"
   pres@vpWidthF           = 0.4
   if (obs_std(0).eq.0) then
       pres@xyLineColors   = "black"
       pres@xyDashPatterns = 2   
       pres@trYMinF        = 0.
       xy_line(wks, Py_pr, y, 0., pres, False)
   else
       pres@xyLineColors   = (/"black","red"/)
       pres@xyDashPatterns = (/2, 0/) 
       pres@trYMinF        = 0.
       xy_line(wks, (/Py_pr,Py/), (/y,y/), 0., pres, False)
   end if

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 output_dir + "." + file_type, \
                 XStg + " vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"trop", "global"/),\
                 (/""/), \
                 (/"wenzel_sabrina"/), \
                 (/"wenzel14jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end 
;#################################################################################
