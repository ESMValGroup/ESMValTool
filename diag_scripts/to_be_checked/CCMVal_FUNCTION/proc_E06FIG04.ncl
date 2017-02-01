;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 4 of                                            #
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
;# The NCL script is called by the NCL shell script E06FIG04.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"

undef("figure4")
function figure4(A:numeric,field:string,plev_sel:numeric,lat_min:numeric,lat_max:numeric,month_sel:string,year_min:numeric,year_max:numeric)

local field,lat_max,lat_min,ua,time,lat,\
lon,plev,ua_mean,v_Clim,v_Clim_sdt,v_lev,v_lev_sdt,nlev,ntime,\
year_maxi,year_mini,year_min,year_max,tim_min,tim_max

begin

   ; tim = time span of the data
   tim = A&time
   tim = coord_list("time",tim)
   ; tim_min/tim_max = time range for baseline calculation (used for calculating anomaly)
   option=0
   option@calendar=tim@calendar
   tim_min = doubletofloat(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,option))
   tim_max = doubletofloat(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,option))
   if ((tim_min+31).lt.tim(0)).or.((tim_max-31).gt.tim(dimsizes(tim)-1)) then
      dum_t=(/tim(0),tim(dimsizes(tim)-1)/)
      copy_VarAtts(tim,dum_t)
      dum_ut = cd_calendar(dum_t,-1)
      print("WARNING: Year range specified "+year_min+"-"+year_max+" is not included in the time series ("+dum_ut(0)+"-"+dum_ut(1)+")!!! ")
      delete(dum_t)
      delete(dum_ut)
   end if

   ;*********************************************************
   if((field.eq."T3D").or.(field.eq."T3M")) then
      ua = A
      uniform_coordinates(ua,field)

      ; sort latitudes if required

      lat = ua&lat
      lat = coord_list("lat",lat)
      lon = ua&lon
      lon = coord_list("lon",lon)
      plev = ua&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then   	
         plev=plev/100.  
         print("Warning: level units was not hPa! automatic plev/100!")
         ua&plev=plev
      end if
      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do
   end if

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then
      ua = A
      uniform_coordinates(ua,field)
      ; sort latitudes if required
      lat = ua&lat
      lat = coord_list("lat",lat)
      plev = ua&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
         ua&plev=plev
      end if
      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do
   end if
   
   time = ua&time
   time = coord_list("time",time)

   ;********************************************************************
   ;                      MONTHLY MEAN
   ;********************************************************************
   ; the field must to be montlhy to do season mean.

   if ((field.eq."T3M") .or. (field.eq."T2Mz")) then
      ua_mean = ua
   end if

   if (field.eq."T3D") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ua_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat),dimsizes(lon)/),typeof(ua))
      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            do k_plev=0,dimsizes(plev)-1
               ua_mean(:,k_plev,j_lat,i_lon) = monthly_mean((/ua(time|:,plev|k_plev,lat|j_lat,lon|i_lon)/),time)
            end do
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ua_mean!1 = "plev"
      ua_mean!2 = "lat"
      ua_mean!3 = "lon"
      ua_mean&plev= ua&plev
      ua_mean&lat = ua&lat
      ua_mean&lon = ua&lon

      field = "T3M"
   end if

   if (field.eq."T2Dz") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ua_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat)/),typeof(ua))
      do i_lat=0,dimsizes(lat)-1
         do j_lev=0,dimsizes(plev)-1
            ua_mean(:,j_lev,i_lat) = monthly_mean((/ua(time|:,plev|j_lev,lat|i_lat)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ua_mean!1 = "plev"
      ua_mean!2 = "lat"
      ua_mean&plev= ua&plev
      ua_mean&lat = ua&lat

      field = "T2Mz"
    end if

   ;********************************************************************
   v_wgtmean1 = new((/dimsizes(ua_mean&time)/),typeof(ua_mean))
   v_wgtmean = new((/dimsizes(ua_mean&time)/),typeof(ua_mean))

   ;*******************************************************************
   ; LAT AND LON WEIGHT MEAN  for the 3D field
   ;*******************************************************************
   if (field.eq."T3M") then
      do ntime = 0,dimsizes(ua_mean&time)-1
         v_wgtmean(ntime) = weight_area_latlon(ua_mean(time|ntime,{plev|plev_sel},lat|:,lon|:),lat,lon,lat_min,lat_max)
      end do
   end if
    
   ;*******************************************************************
   ; LAT WEIGHT MEAN for the zonal field
   ;*******************************************************************
   if (field.eq."T2Mz") then
      do ntime = 0,dimsizes(ua_mean&time)-1
         v_wgtmean(ntime) = weight_area_lat(ua_mean(time|ntime,{plev|plev_sel},lat|:),lat,lat_min,lat_max)
      end do
   end if
    ;print(v_wgtmean)
   ; times = times over which baseline is calculated
   times = ua_mean&time
   times = coord_list("time",times)
   ;********************************************************************
   ;print(month_sel)
   ;do i = 0,dimsizes(month_sel)-1
      if (any(month_sel.eq."ANN")) then
         month_sel_n = (/1,2,3,4,5,6,7,8,9,10,11,12/)
      else
         month_sel_n = new(dimsizes(month_sel),"integer")
      end if
   ;end do

   if (.not.ismissing(ind(month_sel.eq."JAN"))) then
      month_sel_n(ind(month_sel.eq."JAN")) = 1
   end if
   if (.not.ismissing(ind(month_sel.eq."FEB"))) then
      month_sel_n(ind(month_sel.eq."FEB")) = 2
   end if
   if (.not.ismissing(ind(month_sel.eq."MAR"))) then
      month_sel_n(ind(month_sel.eq."MAR")) = 3
   end if
   if (.not.ismissing(ind(month_sel.eq."APR"))) then
      month_sel_n(ind(month_sel.eq."APR")) = 4
   end if
   if (.not.ismissing(ind(month_sel.eq."MAY"))) then
      month_sel_n(ind(month_sel.eq."MAY")) = 5
   end if
   if (.not.ismissing(ind(month_sel.eq."JUN"))) then
      month_sel_n(ind(month_sel.eq."JUN")) = 6
   end if
   if (.not.ismissing(ind(month_sel.eq."JUL"))) then
      month_sel_n(ind(month_sel.eq."JUL")) = 7
   end if
   if (.not.ismissing(ind(month_sel.eq."AUG"))) then
      month_sel_n(ind(month_sel.eq."AUG")) = 8
   end if
   if (.not.ismissing(ind(month_sel.eq."SEP"))) then
      month_sel_n(ind(month_sel.eq."SEP")) = 9
   end if
   if (.not.ismissing(ind(month_sel.eq."OCT"))) then
      month_sel_n(ind(month_sel.eq."OCT")) = 10
   end if
   if (.not.ismissing(ind(month_sel.eq."NOV"))) then
      month_sel_n(ind(month_sel.eq."NOV")) = 11
   end if
   if (.not.ismissing(ind(month_sel.eq."DEC"))) then
      month_sel_n(ind(month_sel.eq."DEC")) = 12
   end if
   
   ;********************************************************************
   ;SELECT month_selS and MEAN

   ; maximum and minimum of data
   tim_mini = cd_calendar(times(0),0)
   tim_maxi = cd_calendar(times(dimsizes(times)-1),0)
   year_mini = floattointeger(tim_mini(0,0))
   year_maxi = floattointeger(tim_maxi(0,0))
   n_year = year_maxi - year_mini
   year = new((/n_year+1/),"integer")
   ; timi = integer representation of time spanning baseline calculation
   timi = cd_calendar(times,1)
   v_wgtmean!0 = "time"
   v_wgtmean&time = timi

   v_selected = new((/n_year+1/),typeof(v_wgtmean))
   v_sel_wgt = new((/n_year+1,dimsizes(month_sel_n)/),typeof(v_wgtmean))
   ; loop over baseline - calculate average
   do k = year_mini,year_maxi
      year(k-year_mini) = k
      do i = 0,dimsizes(month_sel_n)-1
         yyyy_mm=new(1,typeof(timi))
         yyyy_mm = (/k*100+month_sel_n(i)/)
         ind_yymm=ind(v_wgtmean&time.eq.yyyy_mm)
         ;print(ind_yymm)
         if (.not.ismissing(ind_yymm(0))) then
            v_sel_wgt(k-year_mini,i) = (/dim_avg(v_wgtmean(time|ind_yymm(0):ind_yymm(0)))/)
         end if
         delete(ind_yymm)
      end do
   end do
  ; print(v_sel_wgt)
   v_selected = dim_avg(v_sel_wgt)
  ; print(v_selected)
   v_selected!0 = "time"
   v_selected&time = year
;   print(year_min)
;   print(year(0))
;   print(year_max)
;   print(year(dimsizes(year)-1))

   if ((year_min.lt.year(0)).or.(year_max.gt.year(dimsizes(year)-1))) then
      print("ERROR: The year range selected ("+year_min+"-"+year_max+") is outside of the time series ("+year(0)+"-"+year(dimsizes(year)-1)+")")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("!No files or plots are been produced!")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      exit
   end if
   
   v_selected_mean = dim_avg(v_selected({time|year_min:year_max}))

   v_selected = v_selected - v_selected_mean
   if ((lat_min.eq.-90.).and.(lat_max.le.-60)) then
      v_selected = runave_Wrap(v_selected,3,0)
   end if

   if ((lat_max.eq.90.).and.(lat_min.ge.60)) then
      v_selected = runave_Wrap(v_selected,3,0)
   end if

   v_selected@units = ua@units
   v_selected@long_name = ua@long_name
   
   return(v_selected)

end

;***********************************************************************

undef("figure4_trend")
function figure4_trend(A:numeric,field:string,plev_sel:numeric,lat_min:numeric,lat_max:numeric,month_sel:string,trend_year_min:numeric,trend_year_max:numeric)

local field,lat_max,lat_min,ua,time,lat,\
lon,plev,ua_mean,v_Clim,v_Clim_sdt,v_lev,v_lev_sdt,nlev,ntime,\
year_maxi,year_mini,year_min,year_max,tim_min,tim_max,trend,year_trend,tim,\
Timi,ye_max,ye_min,n_dimt,ua_mean, v_wgtmean,v_wgtmean1,times,month_sel_n,tim_mini,\
tim_maxi,year_mini,year_maxi,year_trend, v_selected,v_sel_wgt,year_trend,trend,year_min_utc,year_max_utc
begin

   ; tim = time span of the data
   tim = A&time
   tim = coord_list("time",tim)
    
   ; tim_min/tim_max = time range for trend calculation
   tim_min = doubletofloat(cd_inv_calendar(trend_year_min,1,1,0,0,0,tim@units,0))
   tim_max = doubletofloat(cd_inv_calendar(trend_year_max,12,31,0,0,0,tim@units,0))

   if (((tim_min-31).gt.tim(dimsizes(tim)-1)).or.((tim_max+31).lt.tim(0))) then
      print("Warning!!!!Trend bounds outside the range of the run, all years of the run are selected for the trend!!!!")
      if (typeof(tim).eq."double") then
         tim_min=doubletofloat(tim(0))
         tim_max=doubletofloat(tim(dimsizes(tim)-1))
      else
         tim_min=tim(0)
         tim_max=tim(dimsizes(tim)-1)
      end if
   end if

   if (((tim_min+31).lt.tim(0))) then
      print("Warning!!!!Minimum Trend bound outside the range of the run, the minimum bound of the trend is the first value of the run!!!!")
      if (typeof(tim).eq."double") then
         tim_min=doubletofloat(tim(0))
      else
         tim_min=tim(0)
      end if
   end if

    if (((tim_max-31).gt.tim(dimsizes(tim)-1))) then
       print("Warning!!!!Maximus Trend bound outside the range of the run, the maximus bound of the trend is the last value of the run!!!!")
   
       if (typeof(tim).eq."double") then
          tim_min=doubletofloat(tim(dimsizes(tim)-1))
       else
          tim_max=tim(dimsizes(tim)-1)
       end if
    end if

    tim_min@units=tim@units
    tim_max@units=tim@units

    year_min_utc=cd_calendar(tim_min,0)
    year_max_utc=cd_calendar(tim_max,0)
    year_min=floattointeger(year_min_utc(0,0))
    year_max=floattointeger(year_max_utc(0,0))

   ;***********************************************************************
   if((field.eq."T3D").or.(field.eq."T3M")) then

      ua = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)

      uniform_coordinates(ua,field)
      ; sort latitudes if required
      lat  = ua&lat
      lat  = coord_list("lat",lat)
      lon  = ua&lon
      lon  = coord_list("lon",lon)
      plev = ua&plev
      plev = coord_list("plev",plev)

      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do
      time = ua&time

   end if

   ;***********************************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then

      ua = A({time|tim_min:tim_max},plev|:,lat|:)

      uniform_coordinates(ua,field)
      ; sort latitudes if required
      lat  = ua&lat
      lat  = coord_list("lat",lat)
      plev = ua&plev
      plev = coord_list("plev",plev)

      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do

      time = ua&time
   end if
   time = coord_list("time",time)
;print(ua)

   ;********************************************************************
   ;                      MONTHLY MEAN
   ;********************************************************************
   ; the field must to be montlhy to do season mean.

   if ((field.eq."T3M") .or. (field.eq."T2Mz")) then
      ua_mean = ua
   end if

   if (field.eq."T3D") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ua_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat),dimsizes(lon)/),typeof(ua))
      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            do k_plev=0,dimsizes(plev)-1
               ua_mean(:,k_plev,j_lat,i_lon) = monthly_mean((/ua(time|:,plev|k_plev,lat|j_lat,lon|i_lon)/),time)
            end do
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ua_mean!1 = "plev"
      ua_mean!2 = "lat"
      ua_mean!3 = "lon"
      ua_mean&plev = ua&plev
      ua_mean&lat  = ua&lat
      ua_mean&lon  = ua&lon

      field = "T3M"
   end if

   if (field.eq."T2Dz") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ua_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat)/),typeof(ua))
      do i_lat=0,dimsizes(lat)-1
         do j_lev=0,dimsizes(plev)-1
            ua_mean(:,j_lev,i_lat) = monthly_mean((/ua(time|:,plev|j_lev,lat|i_lat)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ua_mean!1 = "plev"
      ua_mean!2 = "lat"
      ua_mean&plev = ua&plev
      ua_mean&lat  = ua&lat

      field = "T2Mz"
   end if

   ;********************************************************************
   v_wgtmean1 = new((/dimsizes(ua_mean&time)/),typeof(ua_mean))
   v_wgtmean = new((/dimsizes(ua_mean&time)/),typeof(ua_mean))

   ;********************************************************************
   ; LAT AND LON WEIGHT MEAN  for the 3D field
   ;********************************************************************
   if (field.eq."T3M") then
      do ntime = 0,dimsizes(ua_mean&time)-1
         v_wgtmean(ntime) = weight_area_latlon(ua_mean(time|ntime,{plev|plev_sel},lat|:,lon|:),lat,lon,lat_min,lat_max)
      end do
   end if

   ;********************************************************************
   ; LAT WEIGHT  MEAN     for the zonal field
   ;********************************************************************
   if (field.eq."T2Mz") then
      do ntime = 0,dimsizes(ua_mean&time)-1
         v_wgtmean(ntime) = weight_area_lat(ua_mean(time|ntime,{plev|plev_sel},lat|:),lat,lat_min,lat_max)
      end do
   end if

   ; times = times over which trend is calculated
   times = ua_mean&time
   times = coord_list("time",times)

   ;********************************************************************
   ;do i = 0,dimsizes(month_sel)-1
      if (any(month_sel.eq."ANN")) then
         month_sel_n = (/1,2,3,4,5,6,7,8,9,10,11,12/)
      else
         month_sel_n = new(dimsizes(month_sel),"integer")
      end if
   ;end do

   if (.not.ismissing(ind(month_sel.eq."JAN"))) then
      month_sel_n(ind(month_sel.eq."JAN")) = 1
   end if
   if (.not.ismissing(ind(month_sel.eq."FEB"))) then
      month_sel_n(ind(month_sel.eq."FEB")) = 2
   end if
   if (.not.ismissing(ind(month_sel.eq."MAR"))) then
      month_sel_n(ind(month_sel.eq."MAR")) = 3
   end if
   if (.not.ismissing(ind(month_sel.eq."APR"))) then
      month_sel_n(ind(month_sel.eq."APR")) = 4
   end if
   if (.not.ismissing(ind(month_sel.eq."MAY"))) then
      month_sel_n(ind(month_sel.eq."MAY")) = 5
   end if
   if (.not.ismissing(ind(month_sel.eq."JUN"))) then
      month_sel_n(ind(month_sel.eq."JUN")) = 6
   end if
   if (.not.ismissing(ind(month_sel.eq."JUL"))) then
      month_sel_n(ind(month_sel.eq."JUL")) = 7
   end if
   if (.not.ismissing(ind(month_sel.eq."AUG"))) then
      month_sel_n(ind(month_sel.eq."AUG")) = 8
   end if
   if (.not.ismissing(ind(month_sel.eq."SEP"))) then
      month_sel_n(ind(month_sel.eq."SEP")) = 9
   end if
   if (.not.ismissing(ind(month_sel.eq."OCT"))) then
      month_sel_n(ind(month_sel.eq."OCT")) = 10
   end if
   if (.not.ismissing(ind(month_sel.eq."NOV"))) then
      month_sel_n(ind(month_sel.eq."NOV")) = 11
   end if
   if (.not.ismissing(ind(month_sel.eq."DEC"))) then
      month_sel_n(ind(month_sel.eq."DEC")) = 12
   end if
   ;********************************************************************
   ;SELECT month_selS and MEAN
   ; maximum and minimum of data
   tim_mini = cd_calendar(tim(0),0)
   tim_maxi = cd_calendar(tim(dimsizes(tim)-1),0)
   year_mini = floattointeger(tim_mini(0,0))
   year_maxi = floattointeger(tim_maxi(0,0))
   n_year = year_max - year_min
   year = new((/n_year+1/),"integer")

   ; timi = integer representation of time for trend calculation
   timi = cd_calendar(times,1)
   v_wgtmean!0 = "time"
   v_wgtmean&time = timi

   v_selected = new((/n_year+1/),typeof(v_wgtmean))
   v_sel_wgt = new((/n_year+1,dimsizes(month_sel_n)/),typeof(v_wgtmean))

   ; loop over trend period - calculate average
   do k = year_min,year_max
      year(k-year_min) = k
      do i = 0,dimsizes(month_sel_n)-1
         yyyy_mm = k*100+month_sel_n(i)
         ;v_sel_wgt(k-year_min,i) = dim_avg(v_wgtmean({time|yyyy_mm}))
         ind_yymm=ind(v_wgtmean&time.eq.yyyy_mm)
         if (.not.any(ismissing(ind_yymm))) then
            ;v_sel_wgt(k-year_min,i) = (/dim_avg(v_wgtmean(time|ind_yymm))/)
            v_sel_wgt(k-year_min,i) = (/dim_avg(v_wgtmean(time|ind_yymm(0):ind_yymm(dimsizes(ind_yymm)-1)))/)
         end if
         delete(ind_yymm)
      end do
   end do

   v_selected = dim_avg(v_sel_wgt)
   v_selected!0 = "time"
   v_selected&time = year

   ; calculate trend
   year_trend = ispan(year_min,year_max,1)
   trend = dtrend_msg(year_trend,v_selected({time|year_min:year_max}),False,True)
   slope = trend@slope*10

   delete(tim_min)
   delete(tim_max)

   return(slope)

end
;********************************************************************

