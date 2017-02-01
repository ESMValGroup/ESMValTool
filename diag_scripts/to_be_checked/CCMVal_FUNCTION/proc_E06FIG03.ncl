;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 3 of                                            #
;# Eyring et al., Assessment of temperature, trace species and ozone in           #
;# chemistry-climate model simulations of the recent past, J. Geophys. Res., 111, #
;# D22308, doi:10.1029/2006JD007327, 2006.                                        #
;# (hereinafter "E06")                                                            #
;##################################################################################
;##################################################################################
;# Authors and Contact:                                                           #
;#                                                                                #
;# Irene Cionni (i.cionni@niwa.co.nz) and Veronika Eyring (veronika.eyring@dlr.de)#
;##################################################################################
;# The NCL script is called by the NCL shell script E06FIG03.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure3")
function figure3(A:numeric,field:string,plev_sel:numeric,lat_min:numeric,lat_max:numeric,month_sel:string,model:string,year_min:numeric,year_max:numeric)

local ta,field,vName,field,plev,plev_sel,lat_min,lat_max,month_sel,model,year_min,year_max,\
tim,tim_min,tim_max,lat,lon,va, va_mean,v_wgtmean,times,month_sel_n,v_selected,v_sel_wgt,timi,\
tim_min2,tim_max2
begin
   tim = A&time
   tim = coord_list("time",tim)
   option=0
   option@calendar=tim@calendar

   tim_min = doubletofloat(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,option))
   tim_max = doubletofloat(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,option))
      
   if ((tim_min+31).lt.tim(0)).or.((tim_max-31).gt.tim(dimsizes(tim)-1)) then
      dum_t=(/tim(0),tim(dimsizes(tim)-1)/)
      copy_VarAtts(tim,dum_t)
      dum_ut = cd_calendar(dum_t,-1)
      print(model+" ["+dum_ut(0)+"-"+dum_ut(1)+"]  WARNING: Year range specified "+year_min+"-"+year_max+" is not included in the time series!!! ")
      delete(dum_t)
      delete(dum_ut)
   end if

   ;*********************************************************
   if((field.eq."T3D").or.(field.eq."T3M")) then

      ta = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)

      uniform_coordinates(ta,field)     
      lon = ta&lon
      lon = coord_list("lon",lon)
      plev = ta&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then   	
         plev=plev/100.  
         print("Warning: level units was not hPa! automatic plev/100!")
         ta&plev=plev
      end if
   end if

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then
      
      ta = A({time|tim_min:tim_max},plev|:,lat|:)
      uniform_coordinates(ta,field) 
      plev = ta&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then   	
         plev=plev/100.  
         print("Warning: level units was not hPa! automatic plev/100!")
         ta&plev=plev
      end if
   end if
   if((field.eq."T1Ds").or.(field.eq."T1Ms")) then
      ta = A({time|tim_min:tim_max},lat|:)
      uniform_coordinates(ta,field) 
   end if
   ;printVarSummary(ta)
   lat = ta&lat
   lat = coord_list("lat",lat)
  

   time = ta&time
   time = coord_list("time",time)

   do i=0,dimsizes(lat)-1
      if (lat(i).eq.-9999) then
         ;print("error reading latitude, some missing values")
         lat(i) = 0.
      end if
   end do
   ta&lat = lat
   if((field.eq."T1Ds").or.(field.eq."T1Ms")) then
         
         ta = A({time|tim_min:tim_max},lat|:)
         uniform_coordinates(ta,field) 
   end if
   
   lat = ta&lat
   lat = coord_list("lat",lat)
   
   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.
   if (field.eq."T3M") then
      ta_mean = ta
   end if

   if (field.eq."T2Mz") then
      ta_mean = ta     
      ;print(ta_mean)
   end if
   if (field.eq."T1Ms") then
      ta_mean = ta     
   end if
   if (field.eq."T3D") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ta_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat),dimsizes(lon)/),typeof(ta))
      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            do k_plev=0,dimsizes(plev)-1
               ta_mean(:,k_plev,j_lat,i_lon) = monthly_mean((/ta(time|:,plev|k_plev,lat|j_lat,lon|i_lon)/),time)
            end do
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ta_mean!1 = "plev"
      ta_mean!2 = "lat"
      ta_mean!3 = "lon"
      ta_mean&plev = ta@plev
      ta_mean&lat = ta@lat
      ta_mean&lon = ta@lon

      field = "T3M"
   end if

   if (field.eq."T2Dz") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ta_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat)/),typeof(ta))
      do i_lat=0,dimsizes(lat)-1
         do j_lev=0,dimsizes(plev)-1
            ta_mean(:,j_lev,i_lat) = monthly_mean((/ta(time|:,plev|j_lev,lat|i_lat)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ta_mean!1 = "plev"
      ta_mean!2 = "lat"
      ta_mean&plev = ta@plev
      ta_mean&lat = ta@lat

      field = "T2Mz"
   end if

   if (field.eq."T1Dz") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ta_mean = new((/n_dimt,dimsizes(lat)/),typeof(ta))
      do i_lat=0,dimsizes(lat)-1
            ta_mean(:,i_lat) = monthly_mean((/ta(time|:,lat|i_lat)/),time)
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ta_mean!1 = "lat"
      ta_mean&lat = ta@lat

      field = "T1Mz"
   end if
   v_wgtmean=new(dimsizes(ta_mean&time),typeof(ta_mean))
;*******************************************************************
;###################################################################

;    LAT AND LON WEIGHT MEAN  for the 3D field 
;*******************************************************************

   if (field.eq."T3D").or.(field.eq."T3M") then

      do ntime=0,dimsizes(ta_mean&time)-1
         v_wgtmean(ntime)=weight_area_latlon(ta_mean(time|ntime,{plev|plev_sel},lat|:,lon|:),lat,lon,lat_min,lat_max)
      end do
   end if
;###################################################################
;      LAT WEIGHT  MEAN     for the zonal field
;*******************************************************************

   if (field.eq."T2Dz").or.(field.eq."T2Mz") then

      do ntime=0,dimsizes(ta_mean&time)-1
         v_wgtmean(ntime)=weight_area_lat(ta_mean(time|ntime,{plev|plev_sel},lat|:),lat,lat_min,lat_max)         
      end do  
   end if

   if (field.eq."T1Ds").or.(field.eq."T1Ms") then
      do ntime=0,dimsizes(ta_mean&time)-1
         ind_sel=ind((lat.ge.lat_min).and.(lat.le.lat_max))         
         v_wgtmean(ntime)=weight_area_lat(ta_mean(time|ntime,lat|:),lat,lat_min,lat_max)   
      end do
   end if
   times=ta_mean&time
;###########################################################
   month_sel_n=new(dimsizes(month_sel),"integer")

   if (.not.ismissing(ind(month_sel.eq."JAN"))) then
      month_sel_n(ind(month_sel.eq."JAN"))=1
   end if 
   if (.not.ismissing(ind(month_sel.eq."FEB"))) then
      month_sel_n(ind(month_sel.eq."FEB"))=2
   end if 
   if (.not.ismissing(ind(month_sel.eq."MAR"))) then
      month_sel_n(ind(month_sel.eq."MAR"))=3
   end if
   if (.not.ismissing(ind(month_sel.eq."APR"))) then
      month_sel_n(ind(month_sel.eq."APR"))=4
   end if
   if (.not.ismissing(ind(month_sel.eq."MAY"))) then
      month_sel_n(ind(month_sel.eq."MAY"))=5
   end if 
   if (.not.ismissing(ind(month_sel.eq."JUN"))) then
      month_sel_n(ind(month_sel.eq."JUN"))=6
   end if 
   if (.not.ismissing(ind(month_sel.eq."JUL"))) then
      month_sel_n(ind(month_sel.eq."JUL"))=7
   end if
   if (.not.ismissing(ind(month_sel.eq."AUG"))) then
      month_sel_n(ind(month_sel.eq."AUG"))=8
   end if
   if (.not.ismissing(ind(month_sel.eq."SEP"))) then
      month_sel_n(ind(month_sel.eq."SEP"))=9
   end if 
   if (.not.ismissing(ind(month_sel.eq."OCT"))) then
      month_sel_n(ind(month_sel.eq."OCT"))=10
   end if 
   if (.not.ismissing(ind(month_sel.eq."NOV"))) then
      month_sel_n(ind(month_sel.eq."NOV"))=11
   end if
   if (.not.ismissing(ind(month_sel.eq."DEC"))) then
      month_sel_n(ind(month_sel.eq."DEC"))=12
   end if

;print(month_sel)
;print(month_sel_n)

;######################################################
;SELECT month_selS and MEAN

   v_wgtmean!0="time"
   times@units=tim@units
   timi=(cd_calendar(times,-1))

   v_wgtmean&time=timi
   ;printVarSummary(v_wgtmean)
   n_year=year_max-year_min
   ;printVarSummary(n_year)
   ;printVarSummary(v_wgtmean)
   year=new((/n_year+1/),"integer")

   v_selected=new((/n_year+1/),typeof(v_wgtmean))
   v_sel_wgt=new((/n_year+1,dimsizes(month_sel_n)/),typeof(v_wgtmean))	
   do k=year_min,year_max
       
      year(k-year_min)=k
      ;print(year)
      do imon=0,dimsizes(month_sel_n)-1
      a=k*100+month_sel_n(imon)
      yyyy_mm=ind(v_wgtmean&time.eq.(k*100+month_sel_n(imon)))
      ;print(month_sel_n)
         ;printVarSummary(v_wgtmean)
         if(.not.any(ismissing(yyyy_mm))) then
            v_sel_wgt(k-year_min,imon)=dim_avg(v_wgtmean(time|yyyy_mm))
         end if
        
         delete (yyyy_mm)
      end do
   end do
   ;print(v_sel_wgt)

   v_selected=dim_avg(v_sel_wgt)
   v_selected!0="time"
   v_selected&time=year
   v_selected@long_name=ta@long_name
   ;printVarSummary(v_selected)

   return(v_selected)
end
;##################################################################################
