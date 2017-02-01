;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 8 of                                            #
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
;# The NCL script is called by the NCL shell script E06FIG08.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure8")
function figure8(A:numeric,field:string,lat_max:numeric,lat_min:numeric,year_max:numeric,year_min:numeric)

local field,season,lat_max,lat_min,year_max,year_min,tim_min,tim_max,ta,time,season_n,lat,\
lon,plev,ta_mean,v_Clim,v_Clim_sdt,nlev,ntime,v_wgtmean,\
v_wgtmean_sdt,year_min,year_max,tim_min,tim_max,ta_std
begin
   if (field.eq."C2Mz") then

      ta = A(case|0,month|:,plev|:,lat|:)
      ta_std =A(case|1,month|:,plev|:,lat|:)
      uniform_coordinates(ta,field)
      uniform_coordinates(ta_std,field)
      lat = ta&lat

      lat = coord_list("lat",lat)
      plev = ta&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
         ta&plev=plev
         ta_std&plev=plev
      end if
;      printVarSummary(ta)
;      printVarSummary(ta_std)
      v_Clim=dim_rmvmean_Wrap(ta(lat|:,plev|:,month|:))
      v_Clim_std=(ta_std(lat|:,plev|:,month|:))
   else

   tim = A&time
   tim = coord_list("time",tim)
   option =0
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

      ta = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)

      uniform_coordinates(ta,field)
      lon = ta&lon
      lon = coord_list("lon",lon)
      lat = ta&lat
      lat = coord_list("lat",lat)
      plev = ta&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
         ta&plev=plev
      end if

      do i=0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            print("error reading latitude, some missing values")
            lat(i) = 0.
         end if
      end do

      an_ta = dim_rmvmean_Wrap(ta(lat|:,lon|:,plev|:,time|:))
      ta = an_ta(time|:,plev|:,lat|:,lon|:)
      time = ta&time

   end if

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then

      ta = A({time|tim_min:tim_max},plev|:,lat|:)

      uniform_coordinates(ta,field)

      lat = ta&lat
      lat = coord_list("lat",lat)
      plev = ta&plev
      plev = coord_list("plev",plev)
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
         ta&plev=plev
      end if

      do i=0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            print("error reading latitude, some missing values")
            lat(i) = 0.
         end if
      end do
      lat=where(ismissing(lat),0.,lat)
      lat=where(lat.eq.-9999,0.,lat)
      ;  print(lat)
      an_ta = dim_rmvmean_Wrap(ta(lat|:,plev|:,time|:))
      ta = an_ta(time|:,plev|:,lat|:)
      time = ta&time
      ta&lat=lat
   end if

   time = coord_list("time",time)

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.

   if ((field.eq."T3M") .or. (field.eq."T2Mz")) then
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
      ta_mean&lat = lat
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
   v_Clim = Climatology_CCMVal_mean(ta_mean,time,field)
   v_Clim_std = Climatology_CCMVal_std(ta_mean,time,field)
   ;con questa funzione qualunque sia il field ottengo la climatologia cioe# un file medio di 12 mesi se l#input e# mensile o di 365 se l#input e# dayli.
   ;#########################################################
   end if
   ;*******************************************************************
   ;#########################################################
   ;    LAT AND LON WEIGHT MEAN  for the 3D field
   v_wgt = new((/2,24,dimsizes(plev)/),typeof(v_Clim))
   ;v_wgt_sdt = new((/24,dimsizes(plev)/),typeof(v_Clim))

   if (field.eq."T3M") then

      do nlev=0,dimsizes(plev)-1
         do ntime=0,12-1
            lat=v_Clim&lat
            v_wgt(0,ntime,nlev) = weight_area_latlon(v_Clim(month|ntime,plev|nlev,lat|:,lon|:),lat,lon,lat_min,lat_max)
            v_wgt(0,ntime+12,nlev) = weight_area_latlon(v_Clim(month|ntime,plev|nlev,lat|:,lon|:),lat,lon,lat_min,lat_max)
            v_wgt(1,ntime,nlev) = weight_area_latlon(v_Clim_sdt(month|ntime,plev|nlev,lat|:,lon|:),lat,lon,lat_min,lat_max)
            v_wgt(1,ntime+12,nlev) = weight_area_latlon(v_Clim_std(month|ntime,plev|nlev,lat|:,lon|:),lat,lon,lat_min,lat_max)
         end do
      end do
   end if

   ;#########################################################
   ;      LAT WEIGHT  MEAN     for the zonal field

   if (field.eq."T2Mz").or.(field.eq."C2Mz") then

      do nlev=0,dimsizes(plev)-1
         do ntime=0,12-1
            lat = v_Clim&lat
      ;      print(lat)
      ;      print(lat_min)
      ;      print(lat_max)
            v_wgt(0,ntime,nlev) = weight_area_lat(v_Clim(month|ntime,plev|nlev,lat|:),lat,lat_min,lat_max)
            v_wgt(0,ntime+12,nlev) = weight_area_lat(v_Clim(month|ntime,plev|nlev,lat|:),lat,lat_min,lat_max)
            v_wgt(1,ntime,nlev) = weight_area_lat(v_Clim_std(month|ntime,plev|nlev,lat|:),lat,lat_min,lat_max)
            v_wgt(1,ntime+12,nlev) = weight_area_lat(v_Clim_std(month|ntime,plev|nlev,lat|:),lat,lat_min,lat_max)
         end do
      end do
   end if
   v_wgt!0="case"        ;0 mean 1 standard deviation
   v_wgt!1 = "month"
   v_wgt!2 = "plev"
   month2 = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24/)
   plev@units = "hPa"
   v_wgt&month = month2
   v_wgt&plev = plev
   v_wgt@long_name = ta@long_name+" Climatological_mean"
   v_wgt@units = ta@units

   return(v_wgt)

end
;##################################################################################
