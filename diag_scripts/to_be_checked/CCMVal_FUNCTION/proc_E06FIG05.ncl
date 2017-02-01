;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 5 of                                            #
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
;# The NCL script is called by the NCL shell script E06FIG05.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"

undef("figure5")
function figure5(A:numeric,field:string,MODEL:string,year_min:numeric,year_max:numeric)

local field,season,MODEL,lat_max,lat_min,ta,time,season_n,lat,\
lon,plev,ta_mean,v_Clim,v_Clim_sdt,nlev,ntime,v_wgtmean,\
v_wgtmean_sdt,year_min,year_max,tim_min,tim_max

begin

   tim = A&time
   tim = coord_list("time",tim)
   opt=0
   opt@calendar=tim@calendar

; following options for calendars were added by David W. Pierce, the developer of ncview
;   if((tim@calendar).eq."360_day") then
;      opt@calendar="360"
;   end if
;   if((tim@calendar).eq."365_day") then
;      opt@calendar="365"
;   end if
;   if((tim@calendar).eq."noleap") then
;      opt@calendar="365"
;   end if

   tim_min = doubletofloat(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,opt))
   tim_max = doubletofloat(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,opt))
   delete(opt)

   if ((tim_min+31).lt.tim(0)).or.((tim_max-31).gt.tim(dimsizes(tim)-1)) then
       dum_t=(/tim(0),tim(dimsizes(tim)-1)/)
       copy_VarAtts(tim,dum_t)
       dum_ut = cd_calendar(dum_t,-1)
       print(MODEL+" ["+dum_ut(0)+"-"+dum_ut(1)+"]  WARNING: Year range specified "+year_min+"-"+year_max+" is not included in the time series!!! ")
       delete(dum_t)
       delete(dum_ut)
   end if

   ;*********************************************************
   if((field.eq."T3D").or.(field.eq."T3M")) then

      ta = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)

      uniform_coordinates(ta,field)
      lon = ta&lon
      lon = coord_list("lon",lon)

   end if

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then

      ta = A({time|tim_min:tim_max},plev|:,lat|:)
      uniform_coordinates(ta,field)
;      printVarSummary(ta)
   end if

   lat = ta&lat
   lat = coord_list("lat",lat)
   plev = ta&plev
   plev = coord_list("plev",plev)
   if(max(plev).gt.2000.) then
      plev=plev/100.
      print("Warning: level units was not hPa! automatic plev/100!")
      ta&plev=plev
   end if
   if(min(ta).gt.100000.) then
      ta=ta/1000000.
      print("Warning: wrong units")
      ta&plev=plev
   end if

   time = ta&time
  ; time = coord_list("time",time)

   do i=0,dimsizes(lat)-1
      if (lat(i).eq.-9999) then
         ;print("error reading latitude, some missing values")
         lat(i) = 0.
      end if
   end do
   ta&lat = lat

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.
   if (field.eq."T3M") then
      ta_mean = ta
   end if

   if (field.eq."T2Mz") then
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

   ;#########################################################
   ;                      CLIMATOLOGY
   ;#########################################################
   ; ta_mean is the monthly mean
;   printVarSummary(ta_mean)
   v_Clim = Climatology_CCMVal_mean(ta_mean,time,field)
   if (field.eq."T3M") then
      v_Clim_zon = dim_avg_Wrap(v_Clim)
      delete(v_Clim)
      v_Clim = v_Clim_zon
   end if
   ;v_Clim_sdt = Climatology_CCMVal_std(ta_mean,time,field)
   ;con questa funzione qualunque sia il field ottengo la climatologia cioe# un file medio di 12 mesi se l#input e# mensile o di 365 se l#input e# dayli.
   ;##################################################################

   v_Clim@long_name = MODEL+" Climatological_mean"

   return(v_Clim)
end
;##################################################################################

undef("figure5_std")
function figure5_std(A:numeric,field:string,MODEL:string,year_min:numeric,year_max:numeric)

local field,season,MODEL,lat_max,lat_min,ta,time,season_n,lat,\
lon,plev,ta_mean,v_Clim,v_Clim_sdt,nlev,ntime,v_wgtmean,\
v_wgtmean_sdt,year_min,year_max,tim_min,tim_max

begin
   tim = A&time
   tim = coord_list("time",tim)
   tim_min = doubletofloat(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,0))
   tim_max = doubletofloat(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,0))

   if ((tim_min+31).lt.tim(0)).or.((tim_max-31).gt.tim(dimsizes(tim)-1)) then
      dum_t=(/tim(0),tim(dimsizes(tim)-1)/)
      copy_VarAtts(tim,dum_t)
      dum_ut = cd_calendar(dum_t,-1)
      print(MODEL+" ["+dum_ut(0)+"-"+dum_ut(1)+"]  WARNING: Year range specified "+year_min+"-"+year_max+" is not included in the time series!!! ")
      delete(dum_t)
      delete(dum_ut)
   end if
   ;*********************************************************
   if((field.eq."T3D").or.(field.eq."T3M")) then

      ta = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)

      uniform_coordinates(ta,field)
      lon = ta&lon
      lon = coord_list("lon",lon)

   end if

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then
      ta = A({time|tim_min:tim_max},plev|:,lat|:)
      uniform_coordinates(ta,field)
   end if

   lat = ta&lat
   lat = coord_list("lat",lat)
   plev = ta&plev
   plev = coord_list("plev",plev)
   if(max(plev).gt.2000.) then
      plev=plev/100.
      print("Warning: level units was not hPa! automatic plev/100!")
      va&plev=plev
   end if

   time = ta&time
   time = coord_list("time",time)

   do i=0,dimsizes(lat)-1
      if (lat(i).eq.-9999) then
         ;print("error reading latitude, some missing values")
         lat(i) = 0.
      end if
   end do
   ta&lat = lat

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.
   if (field.eq."T3M") then
      ta_mean = ta
   end if

   if (field.eq."T2Mz") then
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

   ;#########################################################
   ;                      CLIMATOLOGY
   ;#########################################################
   ; ta_mean is the monthly mean
   v_Clim_std = Climatology_CCMVal_std(ta_mean,time,field)
   if (field.eq."T3M") then
      v_Clim_zon = dim_avg_Wrap(v_Clim_std)
      delete(v_Clim)
      v_Clim_std = v_Clim_zon
   end if
   ;v_Clim_sdt = Climatology_CCMVal_std(ta_mean,time,field)
   ;con questa funzione qualunque sia il field ottengo la climatologia cioe# un file medio di 12 mesi se l#input e# mensile o di 365 se l#input e# dayli.
   ;##################################################################

   v_Clim_std@long_name = MODEL+" Climatological_mean"

   return(v_Clim_std)
end
;##################################################################################
