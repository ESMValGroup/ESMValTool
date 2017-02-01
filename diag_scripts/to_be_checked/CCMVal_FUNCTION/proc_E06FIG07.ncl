;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 7 of                                            #
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
;# The NCL script is called by the NCL shell script E06FIG07.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure7")
function figure7(A:numeric,field:string,lat_max:numeric,lat_min:numeric,lev_sel,year_min:numeric,year_max:numeric,opt:string)
;==========================
; NOTE: opt = "stdv"/"mean"
;==========================
local field,plev,plev_sel,lat_min,lat_max,month_sel,year_min,year_max,\
tim,tim_min,tim_max,lat,lon,va,va_mean,v_wgtmean,times,month_sel_n,v_selected,v_sel_wgt,timi,opt,\
va1,lat1,ye_max,ye_min,ndimt,v_Clim_zon,v_Clim_sdt_zon,v_Clim,v_Clim_sdt,ind_lat,v_Clim_sel,v_Clim_sdt_sel

begin

   ; tim = time span of the data
   tim = A&time
   tim = coord_list("time",tim)

   ;tim_min/tim_max = time range for averaging
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
   ;*****************************************************************

   if((field.eq."T3D").or.(field.eq."T3M")) then

      va = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)
      ;sort latitudes if required
      uniform_coordinates(va,field)

      lat = va&lat
      lat = coord_list("lat",lat)

      lon = va&lon
      lon = coord_list("lon",lon)

      plev = va&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
         va&plev=plev
      end if
      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do

      time = va&time
   end if

   ;*****************************************************************
;   print(field)
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then
      ; sort latitudes if required
      va = A({time|tim_min:tim_max},plev|:,lat|:)
      uniform_coordinates(va,field)
;      printVarSummary(va)
      
      lat = va&lat
      lat = coord_list("lat",lat)

      plev = va&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
         va&plev=plev
      end if

      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do

      time = va&time

   end if

   ;*****************************************************************
   if((field.eq."T2Ds").or.(field.eq."T2Ms")) then

      va = A({time|tim_min:tim_max},lat|:,lon|:)

      ; sort latitudes if required
      ; sort latitudes if required
      uniform_coordinates(va,field)
      lat = va&lat
      lat = coord_list("lat",lat)
      lon = va&lon
      lon = coord_list("lon",lon)
      do i = 0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            ;print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do
      time = va&time
   end if

   ;********************************************************************
      if((field.eq."T1Ds").or.(field.eq."T1Ms")) then

         va = A({time|tim_min:tim_max},lat|:)

         ; sort latitudes if required
         ; sort latitudes if required
         uniform_coordinates(va,field)
         lat = va&lat
         lat = coord_list("lat",lat)
         do i = 0,dimsizes(lat)-1
            if (lat(i).eq.-9999) then
               ;print("error reading latitude, some missing values")
               lat(i)=0.
            end if
         end do
         time = va&time
      end if
   time = coord_list("time",time)

   ;*****************************************************************
   ;                     MONTHLY MEAN
   ;*****************************************************************
   ; the field must to be monthly to do season mean.

   if((field.eq."T3M") .or. (field.eq."T2Mz") .or. (field.eq."T2Ms").or.(field.eq."T1Ms"))then
      va_mean = va
;      printVarSummary(va_mean)
   end if

   if (field.eq."T3D") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      va_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat),dimsizes(lon)/),typeof(va))
      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            do k_plev=0,dimsizes(plev)-1
               va_mean(:,k_plev,j_lat,i_lon) = monthly_mean((/va(time|:,plev|k_plev,lat|j_lat,lon|i_lon)/),time)
            end do
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      va_mean!1 = "plev"
      va_mean!2 = "lat"
      va_mean!3 = "lon"
      va_mean&plev = va@plev
      va_mean&lat = va@lat
      va_mean&lon = va@lon

      field = "T3M"
   end if

   if (field.eq."T2Dz") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      va_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat)/),typeof(va))
      do i_lat=0,dimsizes(lat)-1
         do j_lev=0,dimsizes(plev)-1
            va_mean(:,j_lev,i_lat) = monthly_mean((/va(time|:,plev|j_lev,lat|i_lat)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      va_mean!1 = "plev"
      va_mean!2 = "lat"
      va_mean&plev = va@plev
      va_mean&lat = va@lat

      field = "T2Mz"
   end if

   if (field.eq."T1Ds") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      va_mean = new((/n_dimt,dimsizes(lat)/),typeof(va))
      do i_lat=0,dimsizes(lat)-1
         va_mean(:,i_lat) = monthly_mean((/va(time|:,lat|i_lat)/),time)
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      va_mean!1 = "lat"
      va_mean&lat = va@lat

      field = "T1Ms"
   end if

   if (field.eq."T2Ds") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      va_mean = new((/n_dimt,dimsizes(lat),dimsizes(lon)/),typeof(va))
      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            va_mean(:,j_lat,i_lon) = monthly_mean((/va(time|:,lat|j_lat,lon|i_lon)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      va_mean!1 = "lat"
      va_mean!2 = "lon"
      va_mean&lat = va@lat
      va_mean&lon = va@lon

      field = "T2Ms"
   end if

   ;*****************************************************************
   ;  CLIMATOLOGY
   ;*****************************************************************
   ; va_mean is the monthly mean
   v_Clim = Climatology_CCMVal_mean(va_mean,time,field)
   v_Clim_sdt = Climatology_CCMVal_std(va_mean,time,field)
   if((field.eq."T3M") .or. (field.eq."T2Ms")) then
      v_Clim_zon = dim_avg_Wrap(v_Clim)
      v_Clim_sdt_zon = dim_avg_Wrap(v_Clim_sdt)
      delete(v_Clim)
      delete(v_Clim_sdt)
      v_Clim = v_Clim_zon
      v_Clim_sdt = v_Clim_sdt_zon
   end if

   ;*****************************************************************
   ; LEVEL and LATITUDE selection
   ;*****************************************************************
   if ((lat_max.eq.0.).and.(lat_min.eq.0.)) then
      ind_lat = ind((lat.ge.(-5.)).and.(lat.le.(5.)))
   else
      ind_lat = ind((lat.ge.(lat_min)).and.(lat.le.lat_max))
   end if

   if((field.eq."T3M") .or. (field.eq."T2Mz")) then
      v_Clim_sel = dim_avg(v_Clim(month|:,{plev|lev_sel},lat|ind_lat))
      v_Clim_sdt_sel = dim_avg(v_Clim_sdt(month|:,{plev|lev_sel},lat|ind_lat))
   end if

   if (field.eq."T2Ms") then
      v_Clim_sel = dim_avg(v_Clim(month|:,lat|ind_lat))
      v_Clim_sdt_sel = dim_avg(v_Clim_sdt(month|:,lat|ind_lat))
   end if

   if (field.eq."T1Ms") then
      v_Clim_sel = dim_avg(v_Clim(month|:,lat|ind_lat))
      v_Clim_sdt_sel = dim_avg(v_Clim_sdt(month|:,lat|ind_lat))
   end if

   v_Clim_sel@units = va@units
   v_Clim_sel@long_name = va@long_name
   v_Clim_sdt_sel@units=va@units
   v_Clim_sdt_sel@long_name=va@long_name

   if (opt.eq."stdv") then
      return(v_Clim_sdt_sel)
   else
      return(v_Clim_sel)
   end if

end
