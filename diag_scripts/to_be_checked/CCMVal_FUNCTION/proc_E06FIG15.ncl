;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 15 of                                           #
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
;# The NCL script is called by the NCL shell script E06FIG15.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure15")
function figure15(A:numeric,field:string,lat_min:numeric,lat_max:numeric,year_min:numeric,year_max:numeric)

local ua,ua1,tim,tim_min,tim_max,lat,lon,lat1,time,ua_mean,v_wgtmean,lev,ntime,lat_min,lat_max,month_sele,\
v_wgtmean1,FMeanAnnCycleMon,FMeanAnnCycleDay,times,coeff_ab,coeff_a,v_season_anom

begin
   tim = A&time

   tim = coord_list("time",tim)
   option=0
   option@calendar=tim@calendar
   tim_min = doubletointeger(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,option))
   tim_max = doubletointeger(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,option))

   ;*********************************************************
   if((field.eq."T2Ds").or.(field.eq."T2Ms")) then

      if (dimsizes(dimsizes(A)).gt.3) then
         ua = A(:,0,:,:)
      else
         ua = A
      end if
      
      uniform_coordinates(ua,field)        
      lat = ua&lat
      lat = coord_list("lat",lat)
      lon = ua&lon
      lon = coord_list("lon",lon)
; printVarSummary(A)
      ; sort latitudes if required
      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i) = 0.
         end if
      end do
      time = ua&time
      time = coord_list("time",time)
   else
      print ("function figure15 only expect field number = T2Ds/T2Ms")
      exit
   end if

   ;##########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do seasonal mean.

   if (field.eq."T2Ds") then
      Timi = cd_calendar(time, 0)

      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))

      n_dimt = (ye_max-ye_min+1)*12

      ua_mean = new((/n_dimt,dimsizes(lat),dimsizes(lon)/),typeof(ua))

      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            ua_mean(:,j_lat,i_lon) = monthly_mean((/ua(:,j_lat,i_lon)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ua_mean!1 = "lat"
      ua_mean!2 = "lon"
      ua_mean&lat = ua@lat
      ua_mean&lon = ua@lon
   else
      ua_mean = ua
   end if

   ua_mean@units = ua@units
   ua_mean@_FillValue = 1e20
   delete(ua)

   ;##################################################################
   ;      LAT WEIGHT  MEAN     for the zonal field
   ;******************************************************************
   v_wgtmean = new((/dimsizes(ua_mean&time)/),typeof(ua_mean),ua_mean@FillValue)
   ;******************************************************************
;    printVarSummary(ua_mean)
;    printVarSummary(v_wgtmean)

   do ntime=0,dimsizes(ua_mean&time)-1    
      v_wgtmean(ntime) = weight_area_latlon(ua_mean(ntime,:,:),lat,lon,lat_min,lat_max)
   end do

   v_wgtmean@units = ua_mean@units
   times = ua_mean&time
   times = coord_list("time",times)
   delete(ua_mean)

   v_wgtmean!0 = "time"
   v_wgtmean&time = times
 
   coeff_ab = calcRegCoeffs(v_wgtmean({tim_min:tim_max}))
   coeff_a = coeff_ab(0:8)
   delete(coeff_ab)

   ;**********************************************************
   ; calculate mean annual cycle
   ; FMeanAnnCycleMon = calcMeanAnnCycleMonthly(coeff_a)  ;monthly
   FMeanAnnCycleDay = calcMeanAnnCycleAnnual(coeff_a)   ;daily

   FMeanAnnCycleDay@units = v_wgtmean@units
   FMeanAnnCycleDay@_FillValue = v_wgtmean@_FillValue
; printVarSummary(FMeanAnnCycleDay)
   return(FMeanAnnCycleDay)

end
;-----------------------------------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;ANOMALIES;;;;;;;;;;ANOMALIES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef("figure15_anom")
;print("procedure")
function figure15_anom(A:numeric,field:string,lat_min:numeric,lat_max:numeric,month_sele:string,year_min:numeric,year_max:numeric)

local ua,ua1,tim,tim_min,tim_max,lat,lon,lat1,time,ua_mean,v_wgtmean,lev,ntime,lat_min,lat_max,\
v_wgtmean1,FMeanAnnCycleMon,FMeanAnnCycleDay,times,coeff_ab,coeff_a,dimx,ntim,nyr,FMonthlyAnomalies,\
FAnnualAnomalies,month_sel_n,month_sele,year,v_season_anom,yyyy_mm,v_selected,ye_max,ye_min,n_dimt,FMeanAnnCycleMon,	year_mini,year_maxi

begin
   tim = A&time
   tim = coord_list("time",tim)
;   printVarSummary(A)
   tim_min = doubletointeger(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,0))
   tim_max = doubletointeger(cd_inv_calendar(year_max,12,31,0,0,0,tim@units,0))
;   print(tim_min)
;   print(tim_max)
   ;*********************************************************
   if((field.eq."T2Ds").or.(field.eq."T2Ms")) then
      if(dimsizes(dimsizes(A)).gt.3) then
         ua=A(:,0,:,:)
      else
         ua = A
      end if 
      uniform_coordinates(ua,field)        
      lat = ua&lat
      lat = coord_list("lat",lat)
      lon = ua&lon
      lon = coord_list("lon",lon)

      ; sort latitudes if required
      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i) = 0.
         end if
      end do
      time = ua&time
      time = coord_list("time",time)
   else
      print ("function figure15 only expect field number = T2Ds/T2Ms")
      exit
   end if
   

   ;##########################################################
   ;                      MONTHLY MEAN
   ;##########################################################
   ; the field must to be montlhy to do seasonal mean.

   if ((field.eq."T2Ds")) then
      Timi = cd_calendar(time, 0)

      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))

      n_dimt = (ye_max-ye_min+1)*12

      ua_mean = new((/n_dimt,dimsizes(lat),dimsizes(lon)/),typeof(ua))

      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            ua_mean(:,j_lat,i_lon) = monthly_mean((/ua(:,j_lat,i_lon)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ua_mean!1 = "lat"
      ua_mean!2 = "lon"
      ua_mean&lat = ua@lat
      ua_mean&lon = ua@lon
   else
      ua_mean = ua
   end if

   ua_mean@_FillValue = ua@_FillValue
   delete(ua)
;    printVarSummary(ua_mean)
   ;##################################################################
   ;      LAT WEIGHT  MEAN     for the zonal field
   ;******************************************************************
   v_wgtmean = new((/dimsizes(ua_mean&time)/),typeof(ua_mean))
   ;*******************************************************************

   do ntime=0,dimsizes(ua_mean&time)-1
      v_wgtmean(ntime) = weight_area_latlon(ua_mean(ntime,:,:),lat,lon,lat_min,lat_max)
   end do

   v_wgtmean@_FillValue = ua_mean@_FillValue
   times = ua_mean&time
   times = coord_list("time",times)
   delete(ua_mean)

   v_wgtmean!0 = "time"
   v_wgtmean&time = times
;   printVarSummary(v_wgtmean)
   ti=v_wgtmean&time
   ;print(ti(0:300))
   ;print(tim_min)
   ;print(tim_max)
;   printVarSummary(v_wgtmean({tim_min:tim_max}))
   coeff_ab = calcRegCoeffs(v_wgtmean({tim_min:tim_max}))
   coeff_a = coeff_ab(0:8)
   delete(coeff_ab)

   ;**********************************************************
   ; calculate mean annual cycle
   FMeanAnnCycleMon = calcMeanAnnCycleMonthly(coeff_a)  ;monthly
   ; FMeanAnnCycleDay = calcMeanAnnCycleAnnual(coeff_a)   ;daily

   FMeanAnnCycleMon@_FillValue = v_wgtmean@_FillValue

   ; remove mean annual cycle and calculate anomalies
   ; set input for rmMeanAnnCycle
   coeff_a1 = coeff_a(0)

   nmon = 12
   ; set outputs for rmMeanAnnCycle
   dimx = dimsizes(v_wgtmean)
   ntim = dimx(0)
   nyr = ntim/nmon

   FMonthlyAnomalies = new((/ntim/),typeof(v_wgtmean))  ;Anomalies montlhy
   FAnnualAnomalies = new((/nyr/),typeof(v_wgtmean))    ;Anomalies Anuual

   rmMeanAnnCycle(v_wgtmean,FMeanAnnCycleMon,coeff_a1,FMonthlyAnomalies,FAnnualAnomalies)

   FMonthlyAnomalies@_FillValue = v_wgtmean@_FillValue
   FAnnualAnomalies@_FillValue=  v_wgtmean@_FillValue
   delete(v_wgtmean)

   ;------------------------------------------------------------------
   ;###########################################################

   if (any(month_sele.eq."ANN")) then
      month_sel_n = (/1,2,3,4,5,6,7,8,9,10,11,12/)
   else
      month_sel_n = new(dimsizes(month_sele),"integer")
   end if

   if (.not.ismissing(ind(month_sele.eq."JAN"))) then
      month_sel_n(ind(month_sele.eq."JAN")) = 1
   end if
   if (.not.ismissing(ind(month_sele.eq."FEB"))) then
      month_sel_n(ind(month_sele.eq."FEB")) = 2
   end if
   if (.not.ismissing(ind(month_sele.eq."MAR"))) then
      month_sel_n(ind(month_sele.eq."MAR")) = 3
   end if
   if (.not.ismissing(ind(month_sele.eq."APR"))) then
      month_sel_n(ind(month_sele.eq."APR")) = 4
   end if
   if (.not.ismissing(ind(month_sele.eq."MAY"))) then
      month_sel_n(ind(month_sele.eq."MAY")) = 5
   end if
   if (.not.ismissing(ind(month_sele.eq."JUN"))) then
      month_sel_n(ind(month_sele.eq."JUN")) = 6
   end if
   if (.not.ismissing(ind(month_sele.eq."JUL"))) then
      month_sel_n(ind(month_sele.eq."JUL")) = 7
   end if
   if (.not.ismissing(ind(month_sele.eq."AUG"))) then
      month_sel_n(ind(month_sele.eq."AUG")) = 8
   end if
   if (.not.ismissing(ind(month_sele.eq."SEP"))) then
      month_sel_n(ind(month_sele.eq."SEP")) = 9
   end if
   if (.not.ismissing(ind(month_sele.eq."OCT"))) then
      month_sel_n(ind(month_sele.eq."OCT")) = 10
   end if
   if (.not.ismissing(ind(month_sele.eq."NOV"))) then
      month_sel_n(ind(month_sele.eq."NOV")) = 11
   end if
   if (.not.ismissing(ind(month_sele.eq."DEC"))) then
      month_sel_n(ind(month_sele.eq."DEC")) = 12
   end if

   ;######################################################
   ;SELECT month_selS and MEAN

   timi = cd_calendar(times,1)
   tim_mini = cd_calendar(tim(0),0)
   tim_maxi = cd_calendar(tim(dimsizes(tim)-1),0)

   year_mini = floattointeger(tim_mini(0,0))
   year_maxi = floattointeger(tim_maxi(0,0))

   FMonthlyAnomalies!0 = "time"
   FMonthlyAnomalies&time = timi

   n_year = ((year_maxi))-((year_mini))
   year = new((/n_year+1/),"integer")

   v_season_anom = new((/n_year+1,dimsizes(month_sel_n)/),typeof(FMonthlyAnomalies))

   do k=year_mini,year_maxi
      year(k-year_mini) = k

      do i=0,dimsizes(month_sel_n)-1
         yyyy_mm = k*100+month_sel_n(i)
         v_season_anom(k-year_mini,i) = (/dim_avg(FMonthlyAnomalies({yyyy_mm}))/)
      end do
   end do

   delete(FMonthlyAnomalies)
   delete(FAnnualAnomalies)
   v_selected = (/dim_avg(v_season_anom)/)

   v_selected!0 = "year"
   v_selected&year = year
   smoothAnomalies(v_selected,30)
;  printVarSummary(v_selected)

   return(v_selected)

end
;##################################################################################
