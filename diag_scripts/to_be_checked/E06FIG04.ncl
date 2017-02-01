;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 4 of                                                       #
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
; Description: This script plots time series of temperature anomalies (area	  #
; averaged) and calculates temperature trends (area weighted). Variables 	  #
; plotted can also be output to netcdf files.	                                  #
; The number of panels in the plot is determined by the number of latitude        #
; bands selected in the info array (info@fig04_lat_min ...)                       #
;;#################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/ta_att.ncl'   (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG04.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06FIG04/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the /var_att/ta_att.ncl file
; info attributes required:
; fig04_trend_year_min  start year of trend calculation
; fig04_trend_year_max  end year of trend calculation
;                       Note:if fig04_trend_year_min and fig04_trend_year_max are
;                       not specified then all years of the run are used in the
;                       trend calculation
; fig04_lat_min         array of latitude (start) values for area averaging
; fig04_lat_max         array of latitude (end) values for area averaging
; fig04_xmin            start year for timeseries plot
; fig04_xmax            start year for timeseries plot
;                       Note:if not specified then 1960 and 2100 are used as
;                       as max and min of the X-Axis.
; fig04_lev_sel         selected pressure level (hPa)
; fig04_season          array of months for averaging. Can use (/"ANN","ANN","ANN"/)
;                       for annual average. Note: this array must have the same
;                       dimensions as fig04_lat_min and fig04_lat_max

;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG04.ncl"
load "./plot_type/get_field_number.ncl"
load "./plot_type/uniform_coordinates.ncl"

begin
   ;fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG04.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************
   ;check the input parameters & translate theprint(OPTIONS_FILE)
    
   ; year_min and year_max = range of years to use as baseline for anomaly calculation
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)
   if(isatt(info,"fig04_trend_year_min") .and. isatt(info,"fig04_trend_year_max")) then
      trend_year_min = info@fig04_trend_year_min
      trend_year_max = info@fig04_trend_year_max
   else
      print("fatal:plot requires variable attributes 'fig04_trend_year_min' and 'fig04_trend_year_max'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if
   if(isatt(info,"fig04_detrend_year_min") .and. isatt(info,"fig04_detrend_year_max")) then
     detrend_year_min = info@fig04_detrend_year_min
     detrend_year_max = info@fig04_detrend_year_max
   else
     print("fatal:plot requires variable attributes 'fig04_detrend_year_min' and 'fig04_detrend_year_max'")
     print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
     print("fatal:!No files or plots are been produced!")
     print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
     status_exit(1)
  end if

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
      
   if(isatt(info,"fig04_lat_min") .and. isatt(info,"fig04_lat_max")) then
      lat_min = info@fig04_lat_min
      lat_max = info@fig04_lat_max
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
   else
      print("fatal:plot requires variable attributes 'fig04_lat_min' and 'fig04_lat_max'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   ; xmin/xmax = plot ranges (min and max years for plot)
   if(isatt(info,"fig04_xmin") .and. isatt(info,"fig04_xmax")) then
      min_x = info@fig04_xmin
      max_x = info@fig04_xmax
   else
      print("plot requires variable attributes 'fig04_xmin' and 'fig04_xmax' if not specified then 1960 and 2100 are used!")
      min_x=1960
      max_x=2100
   end if

   ; lev_sel = pressure level selected
   if(isvar("lev_sel")) then
       delete(lev_sel)
   end if  
 
   if(isatt(info,"fig04_lev_sel")) then
      lev_sel = info@fig04_lev_sel
   else
      print("fatal:plot requires variable attribute 'fig04_lev_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   ; season = season (three months) for corresponding spatial area. Can use (/"ANN", "ANN", "ANN"/)
   if(isvar("season")) then
      delete(season)
   end if   

   if(isatt(info,"fig04_season")) then
      season = info@fig04_season
   else
      print("fatal:plot requires variable attribute 'info@fig04_season'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   dim_MOD = dimsizes(models@name)

   mi_x = 1900
   ma_x = 2200

   yea = ispan(mi_x,ma_x,1)

   v_wgt_sel1 = new((/dim_MOD,dimsizes(yea),dimsizes(lat_min)/),"double")
   v_wgt_sel1!0 = "model"
   v_wgt_sel1!1 = "time"
   v_wgt_sel1!2 = "season"
   v_wgt_sel1&time = yea

   trend = new((/dim_MOD,dimsizes(lat_min)/),"double")

   do k = 0,dimsizes(lat_min)-1

      print("*************************************************************")
      print("*** processing lat range: "+lat_name_max(k)+"-"+lat_name_min(k)+" (set "+(k+1)+" of " +dimsizes(lat_min)+") ***")

      ; MODELS
      do imod = 0,dim_MOD-1
         A0 = extract_data(imod,data,-1,0,0)
         ; tim = time span of this particular model run
         tim = A0&time
         tim = coord_list("time",tim)
         tim_mini = cd_calendar(tim(0),0)
         tim_maxi = cd_calendar(tim(dimsizes(tim)-1),0)
         
         year_mini = floattointeger(tim_mini(0,0))
         year_maxi = floattointeger(tim_maxi(0,0))
         year_m = ispan(year_mini,year_maxi,1)

         print("processing "+models@name(imod)+" ("+(imod+1)+" of " +dim_MOD+")")
         ; call function to calculate temperature anomalies
         v_wgt_sel1(model|imod,{time|year_m},season|k)=figure4(A0,field_number,lev_sel,lat_min(k),lat_max(k),season(k,:),detrend_year_min,detrend_year_max)
         ; call function to calculate temperature trends
         trend(imod,k) = figure4_trend(A0,field_number,lev_sel,lat_min(k),lat_max(k),season(k,:),trend_year_min,trend_year_max)

         year = v_wgt_sel1&time

         delete(A0)
         delete(tim)
         delete(tim_mini)
         delete(tim_maxi)
         delete(year_mini)
         delete(year_maxi)
         delete(year_m)
      end do

   end do

   v_wgt_sel1!0 = "model"
   v_wgt_sel1!1 = "year"
   v_wgt_sel1!2 = "season"

   ;v_wgt_sel1&model = Obs_mod
   v_wgt_sel1&model = models@name
   v_wgt_sel1&year = year

   trend!0="model"
   trend!1="season"

   ;**************************************************************************
   ;******************plot****************************************************

   ;******************************************************************
   ; The function CCMVal_Color and CCMVal_Dash  in misc_function.ncl
   ; associate at any model a color and dash pattern
   ;******************************************************************
   if (project_name.eq."CCMVal2") then
      colors = CCMVal2_ColorTable(models@name)
      dashs = CCMVal2_DashTable(models@name)
   else
      colors = CCMVal1_ColorTable(models@name)
      dashs = CCMVal1_DashTable(models@name)
   end if

   plot_diff = new((/dimsizes(lat_max)/),"graphic")
   lbid = new((/dimsizes(lat_max)/),"graphic")

   ; open worksheet
   wks = gsn_open_wks ("ps",diri_out+"/E06FIG04")

   do i = 0,dimsizes(lat_max)-1

      ;***************************************************************
      res                 = True
      res@gsnDraw         = False
      res@gsnFrame        = False
      res@tfPolyDrawOrder = "Predraw"
      res@trXMaxF         = max_x
      res@trXMinF         = min_x

      res@vpHeightF       = 0.4                    ; change aspect ratio of plot
      res@vpWidthF        = 1.

      res@tiYAxisString = "Anomaly ("+v_wgt_sel1@units+")"    ; y-axis label
      res@tiXAxisString = "Year"     ; y-axis label
      if (season(i,0).eq."ANN") then
         res@tiMainString = v_wgt_sel1@long_name+" "+season(i,0)+" "+lat_name_max(i)+"-"+lat_name_min(i)
      else
         res@tiMainString = v_wgt_sel1@long_name+" "+season(i,0)+"-"+season(i,1)+"-"+season(i,2)+" "+lat_name_max(i)+"-"+lat_name_min(i)
      end if

      res@xyLineColors  = colors                    ; line color
      res@xyDashPatterns = dashs

      plot_diff(i) = gsn_csm_xy (wks,year,v_wgt_sel1(model|:,season|i,year|:),res)

      lgres            = True
      lgres@lgItemType ="Lines"
      lgres@lgPerimOn  = False
      lgres@lgLineColors = colors(::-1)
      lgres@lgDashIndexes = dashs(::-1)
      lgres@vpWidthF   = 0.28
      lgres@vpHeightF  = 0.45
      lgres@lgLabelFontHeightF=0.14
      lgres@lgLineDashSegLenF=0.07

      labels = new(dim_MOD,"string")
      labels = models@name(::-1)+" "+sprintf("%2.2f",trend(::-1,i))
      nitems = dimsizes(labels)
      lbid(i) = gsn_create_legend(wks,nitems,labels,lgres)

      amres                  = True
      amres@amParallelPosF   = 0.67
      amres@amOrthogonalPosF = 0.0

      annoid = gsn_add_annotation(plot_diff(i),lbid(i),amres)

   end do

   Pres        = True
   gsn_panel(wks,plot_diff,(/dimsizes(lat_max),1/),Pres)
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG04.ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   ;**************************************************************************
   ;*******************OUTPUT FILES*******************************************
   if (OPTIONS_FILE.eq."True") then

      do k = 0,dimsizes(lat_name_max)-1

         fili = project_name+"_"+variable+".nc"
         f_out = diri_out+"/E06FIG04/E06FIG04_"+lat_name_max(k)+"_"+lat_name_min(k)+"_trend_"+trend_year_min+"_"+trend_year_max+"_"+fili
         delete(fili)

         system("/bin/rm -f "+f_out)
         fout = addfile(f_out,"c")
         setfileoption(fout,"DefineMode",True)

         ;===================================================================
         ; create global attributes of the file
         ;===================================================================
         fAtt               = True            ; assign file attributes
         fAtt@title         = "Climatological-mean data (month)"
         ;fAtt@source_file = f
         fAtt@creation_date= systemfunc("date")
         fAtt@lat_max      = lat_name_max(k)
         fAtt@lat_min      = lat_name_min(k)
         fAtt@season       = season(k,0)+season(k,1)+season(k,2)
         do im=0,dimsizes(models@name)-1
            str=unique_string("model_")
            fAtt@$str$= models@name(im)
         end do
         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes

         dimNames = (/"model","year"/)
         dimSizes = (/dimsizes(models@name),dimsizes(v_wgt_sel1&year)/)
         dimUnlim = (/False,False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;===================================================================
         var_name = "Anomaly"

         filevardef(fout, "year", typeof(v_wgt_sel1&year), "year")
         filevardef(fout, "trend", typeof(trend), "model")
         filevardef(fout, var_name, typeof(v_wgt_sel1) ,(/"model","year"/))
         ;===================================================================

         setfileoption(fout,"DefineMode",False)

         fout->year = (/v_wgt_sel1&year/)
         fout->$var_name$ = (/v_wgt_sel1(model|:,year|:,season|k)/)
         fout->trend =(/trend(model|:,season|k)/)
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("!")
         print("!  Files are been produced in "+diri_out+"/E06FIG04/  !!")
         print("!")
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

      end do
   else
      print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
   end if

   delete(lat_min)
   delete(lev_sel)
   delete(yea)
   delete(nitems)
   delete(dim_MOD)
   delete(lat_name_max)
   delete(lat_name_min)
   delete(year_max)
   delete(dashs)
   delete(trend)
   delete(colors)
   delete(season)
   delete(plot_types)
   delete(res)
   delete(k)
   delete(Pres)
   delete(i)
   delete(labels)
   delete(lgres)
   delete(lat_max)
   delete(amres)
   delete(v_wgt_sel1)
   delete(year)
   delete(year_min)
   delete(trend_year_max)
   delete(trend_year_min)
   ;delete(Obs_mod)
   ;delete(OBS)

   print ("=================== FINISH (E06FIG04) ====================")

end
;**************************************************************************
