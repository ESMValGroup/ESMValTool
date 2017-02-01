;;#############################################################################
;; carbon_constraint.ncl
;; Author: Sabrina Wenzel (DLR, Germany)    
;; EMBRACE project 
;;#############################################################################
;;                                                                            
;; Description:
;;    Calculates the long term sensitivity GAMMA_LT according to Wenzel et al. 
;;    2013 with observations from GCP.
;;    This plotscript calls a file wich contains the results from calculating 
;;    GAMMA_IAV (carbon_corr_2var.ncl)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    * con_latrange   array of latitudes for time series plots if 2D or 3D 
;;    * con_lonrange   array of pressure levels for time series plots if 2D or 3D
;;    * gIAV_file      file path for gamma_IAV values calculated with 
;;                     carbon_corr_2var.ncl
;;
;; Modification history
;;    * 20151027-A_laue_ax: added call to 'write_references'
;;    * 201509-A_Sabrina Wenzel: updated to new ESMValTool structure
;;    * 201303-A_Sabrina Wenzel: written 
;; #############################################################################
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"

load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/carbon/statistics.ncl"
load "./diag_scripts/uniform_coordinates.ncl"

load "./plot_scripts/ncl/scatterplot.ncl"
load "./plot_scripts/ncl/xy_line.ncl"
load "./plot_scripts/ncl/legends.ncl"

begin
  diag_script = "carbon_constraint"
  
  verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
  info_output("<<<<<<<< Entering " + diag_script + ".ncl", verbosity, 4)
  info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
  info_output( diag_script + ".ncl (var: " + variables(0) + \
                              " vs. " + variables(1) + ")", verbosity, 1)
  info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

; A-laue_ax+
    ; ###########################################
    ; # references                              #
    ; ###########################################

    write_references(diag_script,      \  ; script name
                     "A_wenz_sa",      \  ; authors
                     "",               \  ; contributors
                     "D_wenzel14jgr",  \  ; diag_references
                     "",               \  ; obs_references
                     (/""/))              ; proj_references
; A-laue_ax-

  ;; fields supported
  field_number = field_types
  field = field_number(0)

  var0 = variables(0)
  var1 = variables(1)

  dim_MOD = dimsizes(models@name)
  dim_VAR = dimsizes(variables)

  data_0 = True
  data_1 = True
  do idx_mod = 0, dim_MOD - 1
    data_0 = read_data(idx_mod, var0, field)
    data_1 = read_data(idx_mod, var1, field)
  end do
end

; Load var_info
begin
    var0deffile = "./interface_data/" + var0 + "_info.tmp"
    loadscript(var0deffile)
    var1deffile = "./interface_data/" + var1 + "_info.tmp"
    loadscript(var1deffile)
end

begin

  ;****************************************************************** 
  ;check the input parameters 

  if(isvar("latrange")) then 
    delete(latrange) 
  end if 
  if(isatt(diag_script_info,"con_latrange")) then 
    latrange = diag_script_info@con_latrange 
  else 
    print("plot requires variable attribute 'con_latrange'") 
    exit 
  end if 
  
  if(isvar("lonrange")) then 
    delete(lonrange) 
  end if 
  if(isatt(diag_script_info,"con_lonrange")) then 
    lonrange = diag_script_info@con_lonrange 
  end if 
  
  if(isatt(diag_script_info,"con_units")) then 
    UNITS = diag_script_info@con_units 
  end if 

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
  ;; loop for models
  do imod = 0,dim_MOD-1 
    info_output(models@name(imod),verbosity,1)
    
    ;;get data 
    A0 := extract_data_extended(imod,var0,data_0,-1,0,0)
    uniform_coordinates(A0,field)

    A1 := extract_data_extended(imod,var1,data_1,-1,0,0)
    uniform_coordinates(A1,field)

    file2 := systemfunc("ls " + climo_dir + "/"+ models@project(imod)+ "/*" \
                              + diag_script_info@pcBGC + "_" + models@name(imod) \
                              + "_" + models@ensemble(imod) + "_" + field \
                              + "_" + var0 + "*.nc")

    f2 := addfile(file2,"r")
    A2 := f2->$var0$
    uniform_coordinates(A2,field)

    ;;average yearly
    tmpaA0 := time_operations(A0,-1,-1,"average","yearly",True)
    tmpaA1 := time_operations(A1,-1,-1,"average","yearly",True)
    tmpaA2 := time_operations(A2,toint(models@start_year(imod)),\
                                 toint(models@end_year(imod)),\
                                      "average","yearly",True)

    ;;sum VAR0 and VAR2 over spatial region
    tmpA0 := area_operations(tmpaA0,latrange(0),latrange(1),lonrange(0),lonrange(1),"sum",True)
    tmpA2 := area_operations(tmpaA2,latrange(0),latrange(1),lonrange(0),lonrange(1),"sum",True)

    ;;average VAR1 over spatial region
    tmpA1 := area_operations(tmpaA1,latrange(0),latrange(1),lonrange(0),lonrange(1),"average",True)

    ;; Convert units if appropriate
    if (isatt(variable_info, "plot_units")) then
        tmpA0@units = A0@units
        tmpA0 := convert_units(tmpA0,variable_info@plot_units)
        tmpA2@units = A2@units
        tmpA2 := convert_units(tmpA2,variable_info@plot_units)
    end if

    tmpcA0 = new((/dimsizes(tmpA0)+1/),"float")
    tmpcA2 = new((/dimsizes(tmpA2)+1/),"float")
    tmpcA0(0) = 0.0
    tmpcA2(0) = 0.0

    do y = 0,nMax-1
        tmpcA0(y+1) = tmpcA0(y)+tmpA0(y)
        tmpcA2(y+1) = tmpcA2(y)+tmpA2(y)
    end do
    delete(y)

    ;;calculate integral line:
    tidx := dimsizes(tmpA0)-1
    Y1data := avg(tmpcA0(tidx-9:tidx))-avg(tmpcA0(:9))
    Y2data := avg(tmpcA2(tidx-9:tidx))-avg(tmpcA2(:9))
    Y3data := avg(tmpA1(tidx-9:tidx))-avg(tmpA1(:9))

    ;;calculate gamma_LT:
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

   ;; Read gIAV values
   xml_name = getenv("ESMValTool_xml_name")
   suffix = get_file_suffix(xml_name, 0)
   name = str_sub_str(xml_name, suffix, "")
   temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
                name + "_gIAV.nc"

   if (isfilepresent_esmval(temp_dir)) then
       ;; Read grading list
       temp_file = addfile(temp_dir, "r")
       temp_list = temp_file->temp_list
       file_list = tostring(temp_list)
       file_list@_FillValue = default_fillvalue(typeof(file_list))
   else
       error_msg("f", diag_script, "", "cannot find gIAV file")
   end if

   ;;read file with GAMMA_IAV values
   IAV = ncdf_read(file_list, "gIAV")

   if (dimsizes(IAV(0,:)).gt.dim_MOD) then
       Y6data  = IAV(:,1:dim_MOD)         ;;Gamma_IAV from the models
       obs_std = IAV(:,0)                 ;;Gamma_IAV from the OBS
   else
       Y6data = IAV
       obs_std = (/0,0/)
   end if

   ;; Remove temporary list
   system("rm -f " + temp_dir)
   delete([/IAV,xml_name,suffix,name,temp_dir,temp_file,temp_list/])

   ;;exclude models from regression if appropriate
   ind_MOD = diag_script_info@reg_models
   if (all(ind_MOD.eq.(/"all"/))) then
       aX = Y6data
       aY = gamma_LT
   else
       aX = Y6data(:,ind_MOD)
       aY = gamma_LT(ind_MOD)
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
   diag_script_base = basename(diag_script)
   output_dir = get_output_dir(plot_dir, diag_script_base)

   ;;Define path and names of plot files  
   file_type = getenv("ESMValTool_output_file_type")
   if(ismissing(file_type)) then
       file_type = "PS"
   end if

   plot_file = "constr_"+var0+"-"+var1+"_"+xMin_year+"-"+xMax_year
   wks = gsn_open_wks(file_type,output_dir+plot_file)

   colors = project_style(diag_script_info, "colors")
   thicks = project_style(diag_script_info, "thicks")
   plot_num = integertochar(ispan(97,122,1))

   XStg    = "Sensitivity of land+ocean CO~B~2~N~ flux, ~F8~g~F21~~B~LT~N~ (GtC/yr/K)"
   YStg    = "Carbon Cycle Climate Feedback, ~F8~g~F21~~B~LT~N~ (GtC/K)"

   data_arr = (/Y6data(0,:), gamma_LT/)
   data_arr!0 = "models"
   data_arr&models = (/XStg,YStg/)
   data_arr@diag_script = diag_script

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


  info_output("gv " + output_dir + plot_file + "." + file_type, verbosity, 1)

end 
;#################################################################################
