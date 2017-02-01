;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 15 of                                                      #
;# Eyring et al., Assessment of temperature, trace species and ozone in           #
;# chemistry-climate model simulations of the recent past, J. Geophys. Res., 111, #
;# D22308, doi:10.1029/2006JD007327, 2006.                                        #
;# (hereinafter "E06")                                                            #
;##################################################################################
;# Authors and Contact:                                                           #
;#                                                                                #
;# Irene Cionni (irene.cionni@enea.it) & Veronika Eyring (veronika.eyring@dlr.de) #
;##################################################################################
;# The NCL script has been developed for use within the CCMVal diagnostic tool    #
;# (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)                 #
;##################################################################################
;# Description: This script is designed to plot seasonal total column ozone       #
;# anomalies time series and associated seasonal cycles for a selected number     #
;# of latitude bands and seasons. Variables plotted can also be output to         #
;# netcdf files.                                                                  #
;# The number of panels in the plot is determined by the number of latitude       #
;# bands selected in the info array (info@fig15_lat_min ...).                     #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/toz_att.ncl'   (${WKDIR} specified in main namelist)  #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG15.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06FIG15/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the toz_att.ncl file
; info attributes required:
; fig15_lat_min         array of latitude (start) values for area averaging
; fig15_lat_max         array of latitude (end) values for area averaging
; fig15_season          array of months for averaging. Can use (/"ANN","ANN","ANN"/)
;                       for annual average. Note: this array must have the same
;                       dimensions as fig15_lat_min and fig15_lat_max
; fig15_xmin            plot range (x start) for anomaly time series
; fig15_xmax            plot range (x end) for anomaly time series
;
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG15.ncl"
load "./plot_type/get_field_number.ncl"
begin
   print("Begin ./plot_type/E06FIG15.ncl")
   ;fields supported T2Ms T2Ds
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************
   ;check the input parameters & translate them
   ; year_min and year_max = range of years to use for climatological analysis

   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   ; lat_min, lat_max = latitude ranges for spatial averaging
   if(isvar("lat_min")) then
      delete(lat_min)
   end if   
   if(isvar("lat_max")) then
      delete(lat_max)
   end if     
   if(isvar("lat_name_min")) then
      delete(lat_name_min)
   end if  
   if(isvar("lat_name_max")) then
      delete(lat_name_max)
   end if    
   if(isatt(info,"fig15_lat_min") .and. isatt(info,"fig15_lat_max")) then
      lat_min = info@fig15_lat_min
      lat_max = info@fig15_lat_max
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
   else
      print("fatal:plot requires variable attributes 'fig15_lat_min' and 'fig15_lat_max'")
      status_exit(1)
   end if

    ; xmin/xmax = plot ranges (min and max years for plot)
   if(isatt(info,"fig15_xmin") .and. isatt(info,"fig15_xmax")) then
      min_x = info@fig15_xmin
      max_x = info@fig15_xmax
   else
      print("fatal:plot requires variable attributes 'fig15_xmin' and 'fig15_xmax'")
      status_exit(1)
   end if

   ; season = season (three months) for corresponding spatial area. Can use (/"ANN", "ANN", "ANN"/)
   if(isvar("season")) then
      delete(season)
   end if    
   if(isatt(info,"fig15_season")) then
      season = info@fig15_season
   else
      print("fatal:plot requires variable attribute 'info@fig15_season'")
      status_exit(1)
   end if

   ; optional input parameters
   if(isatt(info,"fig15_obs_name") .and. isatt(info,"fig15_obs_file")) then
      OBS = info@fig15_obs_name
      obs_filename = info@fig15_obs_file
   else
      OBS = (/"NaN"/)
   end if

   ;******************************************************************
   dim_OBS = dimsizes(OBS)
   dim_MOD = dimsizes(models@name)
   if (OBS.ne."NaN") then
      Obs_mod = new(dim_OBS+dim_MOD,"string")
      Obs_mod(0) = OBS
      Obs_mod(1:) = models@name
   else
      dim_OBS = 0
      Obs_mod = new(dim_MOD,"string")
      Obs_mod = models@name
   end if

   ;******************************************************************
   mi_x = 1900
   ma_x = 2150
   t_array = new((/12*(ma_x-mi_x+1)/),"integer")

   do n_year=mi_x,ma_x
      do n_month=0,11
         t_array((n_year-mi_x)*12+n_month) = n_year*100+(n_month+1)
      end do
   end do

   y_array = ispan(mi_x,ma_x,1)
   
   v_wgt_sel1 = new((/(dim_MOD+dim_OBS),dimsizes(y_array),dimsizes(lat_min)/),"float")
   v_wgt_sel1!0 = "model"
   v_wgt_sel1!1 = "time"
   v_wgt_sel1!2 = "season"
   v_wgt_sel1&time = y_array
   v_wgt_sel1&model = Obs_mod
   v_wgt_sel1@units = "DU"

   v_ann_sel = new((/(dim_MOD+dim_OBS),366,dimsizes(lat_min)/),"float")
   v_ann_sel!0 = "model"
   v_ann_sel!1 = "time"
   v_ann_sel!2 = "season"
   days = ispan(1,366,1)
   v_ann_sel&time = days
   v_ann_sel&model = Obs_mod
  
   ;//////////////////////////here////////////////////////////////////

   do k=0,dimsizes(lat_min)-1

      print("******************************************************************")
      print("*** processing lat range: "+lat_name_max(k)+"-"+lat_name_min(k)+" (set "+(k+1)+" of " +dimsizes(lat_min)+") ***")

      if (dim_OBS.ne.0) then

         ;A0 = extract_data(obs,data,0)
         f_obs = addfile(obs_filename,"r")
         A0 = f_obs->$variable$
         field = get_field_number(A0)
         delete(f_obs)
         if (typeof(A0).eq."double") then
           A00=doubletofloat(A0)
           copy_VarCoords(A0,A00)
           delete (A0)
           A0=A00
           delete(A00)
         end if

         tim = A0&time
         tim = coord_list("time",tim)
         tim_mini = cd_calendar(tim(0),0)
         tim_maxi = cd_calendar(tim(dimsizes(tim)-1),0)

         year_mini = floattointeger(tim_mini(0,0))
         year_maxi = floattointeger(tim_maxi(0,0))

         ; call function to calculate annual cycle
         v_wgt_sel1(model|0,{time|year_mini:year_maxi},season|k) = (/figure15_anom(A0,field,lat_min(k),lat_max(k),season(k,:),info@fig15_yr_min,info@fig15_yr_max)/)
         ; call function to calculate anomalies
         v_ann_sel(0,:,k) = (/figure15(A0,field,lat_min(k),lat_max(k),info@fig15_yr_min,info@fig15_yr_max)/)

         delete(A0)
         delete(tim)
         delete(tim_mini)
         delete(tim_maxi)
         delete(year_mini)
         delete(year_maxi)
      end if

      do imod = 0,dim_MOD-1
         A0 = extract_data(imod,data,-1,0,0)
         if (typeof(A0).eq."double") then
           A00=doubletofloat(A0)
           copy_VarCoords(A0,A00)
           delete (A0)
           A0=A00
           
           delete(A00)
         end if
         ; tim = time span of this particular model run
         tim = A0&time
         tim = coord_list("time",tim)
         tim_mini = cd_calendar(tim(0),0)
         tim_maxi = cd_calendar(tim(dimsizes(tim)-1),0)

         year_mini = floattointeger(tim_mini(0,0))
         year_maxi = floattointeger(tim_maxi(0,0))

         print("processing "+models@name(imod)+" ("+(imod+1)+" of " +dim_MOD+")")
	 v_wgt_sel1(model|imod+dim_OBS,{time|year_mini:year_maxi},season|k) = (/figure15_anom(A0,field_number,lat_min(k),lat_max(k),season(k,:),year_min(imod),year_max(imod))/)
         v_ann_sel(model|imod+dim_OBS,time|:,season|k) = (/figure15(A0,field_number,lat_min(k),lat_max(k),year_min(imod),year_max(imod))/)

         year = v_wgt_sel1&time

         delete(A0)
         delete(tim)
         delete(tim_mini)
         delete(tim_maxi)
         delete(year_mini)
         delete(year_maxi)
      end do
   end do

   ;******************************************************************
   if (project_name.eq."CCMVal2") then
       colors = CCMVal2_ColorTable(Obs_mod)
       dashs = CCMVal2_DashTable(Obs_mod)
   else
       colors = CCMVal1_ColorTable(Obs_mod)
       dashs = CCMVal1_DashTable(Obs_mod)
   end if

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG15_"+variable)

   ;***************  plot  *******************************************
   plot_diff = new((/dimsizes(lat_max)/),"graphic")
   plot_ann = new((/dimsizes(lat_max)/),"graphic")
   baseplot = new((/dimsizes(lat_max)*2/),"graphic")

   do i=0,dimsizes(lat_max)-1

      ;***************************************************************
      res                  = True
      res@gsnDraw          = False
      res@gsnFrame         = False
      res@tfPolyDrawOrder  = "Predraw"
      res@trXMaxF          = max_x
      res@trXMinF          = min_x

      res@vpHeightF        = 0.3                   ; change aspect ratio of plot
      res@vpWidthF         = 0.6

      res@tiYAxisString    = "Anomaly(%) "         ; y-axis label
      res@tiXAxisString    = "Year"                ; y-axis label
      if (season(i,0).eq."ANN") then
         res@tiMainString  = season(i,0)+" "+lat_name_max(i)+"-"+lat_name_min(i)
      else
         res@tiMainString  = season(i,0)+"-"+season(i,1)+"-"+season(i,2)+" "+lat_name_max(i)+"-"+lat_name_min(i)
      end if

      res@xyLineColors     = colors                ; line color
      res@xyDashPatterns   = dashs

      plot_diff(i) = gsn_csm_xy (wks,v_wgt_sel1&time,v_wgt_sel1(model|:,season|i,time|:),res)

      res2                 = True
      res2@gsnDraw         = False
      res2@gsnFrame        = False

      res2@trXMaxF         = 366
      res2@trXMinF         = 1
      res2@tmXBMode        = "Explicit"
      res2@tmXBValues      = (/1,32,61,92,122,153,183,214,245,275,305,336/)
      res2@tmXBLabels      = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)
      ;res@tmLabelAutoStride = True

      res2@vpHeightF       = 0.3                   ; change aspect ratio of plot
      res2@vpWidthF        = 0.2

      res2@tiYAxisString   = "Detrended mean annual cycle ("+v_wgt_sel1@units+")"    ; y-axis label
      res2@tiXAxisString   = "Year"     ; y-axis label
      if (season(i,0).eq."ANN") then
         res2@tiMainString = season(i,0)+" "+lat_name_max(i)+"-"+lat_name_min(i)
      else
         res2@tiMainString = season(i,0)+"-"+season(i,1)+"-"+season(i,2)+" "+lat_name_max(i)+"-"+lat_name_min(i)
      end if

      res2@xyLineColors    = colors                ; line color
      res2@xyDashPatterns  = dashs

      plot_ann(i) = gsn_csm_xy (wks,v_ann_sel&time,v_ann_sel(model|:,season|i,time|:),res2)

   end do ; end loop season

   do j=0,2
      baseplot(2*j) = plot_diff(j)
      baseplot(2*j+1) = plot_ann(j)
   end do

   Pres          = True
   Pres@gsnFrame = False
   Pres@gsnPanelBottom = 0.33
   gsn_panel(wks,baseplot,(/3,2/),Pres)

   add_legenda(wks,Obs_mod,colors,dashs,"bottom")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG15.ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

   ;*******************OUTPUT FILES*******************************************
   if (OPTIONS_FILE.eq."True") then
      do k = 0,dimsizes(lat_max(k))-1
 
         fili = project_name+"_"+field_number+"_"+variable+".nc"
         f_out = diri_out+"/E06FIG15/E06FIG15_"+lat_max(k)+lat_min(k)+"_"+fili
         delete(fili)

         system("/bin/rm -f "+f_out)
         fout = addfile(f_out,"c")
         setfileoption(fout,"DefineMode",True)

         ;===================================================================
         ; create global attributes of the file
         ;===================================================================
         fAtt               	= True            ; assign file attributes
         fAtt@title 		= "Climatological-mean data (month)"
         ;fAtt@source_file	= f
         fAtt@creation_date	= systemfunc("date")
         fAtt@lat_max		= lat_max(k)
         fAtt@lat_min		= lat_min(k)
         fAtt@season		= season(k,0)+season(k,1)+season(k,2)
         do im=0,dimsizes(Obs_mod)-1                
	    str=unique_string("model_")
	    fAtt@$str$= Obs_mod(im)
         end do
         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes

         dimNames = (/"model","year","month"/)
         dimSizes = (/dimsizes(Obs_mod),dimsizes(v_wgt_sel1&time),dimsizes(v_ann_sel&time)/)
         dimUnlim = (/False,False,False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;===================================================================
         var_name1 = "Anomaly"
         var_name2 = "Annual_cycle"
         filevardef(fout, "year", typeof(v_wgt_sel1&time), "year")
         filevardef(fout, "month", typeof(v_ann_sel&time), "month")
         filevardef(fout, var_name2, typeof(v_ann_sel), (/"model","month"/))
         filevardef(fout, var_name1, typeof(v_wgt_sel1) ,(/"model","year"/))
         ;===================================================================

         setfileoption(fout,"DefineMode",False)

         fout->year = (/v_wgt_sel1&time/)
         fout->month= (/v_ann_sel&time/)
         fout->$var_name1$ = (/v_wgt_sel1(model|:,time|:,season|k)/)
         fout->$var_name2$ = (/v_ann_sel(model|:,season|k,time|:)/)
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("!")
         print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG15/  !!")
	 print("!")
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      end do
   else
      print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
   end if
   delete(lat_min)
   delete(t_array)
   delete(dim_MOD)
   delete(lat_name_max)
   delete(lat_name_min)
   delete(n_month)
   delete(v_ann_sel)
   delete(dashs)
   delete(days)
   delete(min_x)
   delete(colors)
   delete(season)
   delete(res)
   delete(k)
   delete(res2)
   delete(j)
   delete(Pres)
   delete(i)
   delete(max_x)
   delete(lat_max)
   delete(imod)
   delete(ma_x)
   delete(v_wgt_sel1)
   delete(year)
   delete(models)
   delete(y_array)
   delete(mi_x)
   delete(Obs_mod)
   delete(n_year)
   print ("=================== FINISH (E06FIG15) ====================")

end
