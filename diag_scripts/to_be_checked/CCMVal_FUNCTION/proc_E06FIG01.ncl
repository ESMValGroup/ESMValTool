;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 1 of                                            #
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
;# The NCL script is called by the NCL shell script E06FIG01.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure1")
function figure1(A:numeric,field:string,season:string,MODEL:string,lat_max:numeric,lat_min:numeric,year_min:numeric,year_max:numeric)

local field,season,MODEL,lat_max,lat_min,ta,time,season_n,lat,\
lon,plev,ta_mean,v_Clim,v_Clim_sdt,v_seas,v_seas_sdt,nlev,ntime,v_wgtmean,\
v_wgtmean_sdt,year_min,year_max,tim_min,tim_max
begin

   if ((field.ne."C2Mz").and.(field.ne."C3M").and.(field.ne."C1M")) then
;      print(field)
      tim = A&time
      tim = coord_list("time",tim)
      option=0
      option@calendar=tim@calendar
      tim_min = doubletofloat(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,option))
      tim_max = doubletofloat(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,option))

      if ((tim_min+31).lt.tim(0)).or.((tim_max-31).gt.tim(dimsizes(tim)-1)) then
         if ((MODEL.eq."ERA_40").or.(MODEL.eq."ERA40").or.(MODEL.eq."OBS").or.(MODEL.eq."UKMO")) then
            print("WARNING:  the observation years inside the file are outside the specified climatological bounds then ALL data in the file is used to construct the climatology. ")
         else
            print("WARNING: the model years inside the file are outside the specified climatological bounds then ALL data in the file is used to construct the climatology. Select years in the file namelist in the range of the run.")
         end if
         if (typeof(tim) .eq. "double") then
            tim_min = doubletofloat(tim(0))
            tim_max = doubletofloat(tim(dimsizes(tim)-1))
         else
            tim_min = tim(0)
            tim_max = tim(dimsizes(tim)-1)
         end if
      end if
   end if

   season_n = new(dimsizes(season),"integer")

   if (.not.ismissing(ind(season.eq."DJF"))) then
      season_n(ind(season.eq."DJF")) = 120102
   end if
   if (.not.ismissing(ind(season.eq."MAM"))) then
      season_n(ind(season.eq."MAM")) = 030405
   end if
   if (.not.ismissing(ind(season.eq."JJA"))) then
      season_n(ind(season.eq."JJA")) = 060708
   end if
   if (.not.ismissing(ind(season.eq."SON"))) then
      season_n(ind(season.eq."SON")) = 091011
   end if

   ;*********************************************************
   if ((field.ne."C2Mz").and.(field.ne."C3M").and.(field.ne."C1M")) then
      if((field.eq."T3D").or.(field.eq."T3M")) then
;         if (isdim(A,("plev"))) then
;            dim="plev"
;         end if
;         if (isdim(A,("lev"))) then
;            dim="lev"
;         end if
;         if (isdim(A,("p"))) then
;            dim="p"
;         end if
;         ta = A({time|tim_min:tim_max},$dim$|:,lat|:,lon|:)
;         ta!1="plev"
         ta = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)
         uniform_coordinates(ta,field)
         lon = ta&lon
         lon = coord_list("lon",lon)
      end if

   ;*********************************************************
      if((field.eq."T2Dz").or.(field.eq."T2Mz")) then
;         if (isdim(A,("plev"))) then
;            dim="plev"
;         end if
;         if (isdim(A,("lev"))) then
;            dim="lev"
;         end if
;         if (isdim(A,("p"))) then
;            dim="p"
;         end if

         if (iscoord(A,("lon")).and.(dimsizes(A&lon).lt.2)) then
            ta = A({time|tim_min:tim_max},plev|:,lat|:,lon|0)
;            ta = A({time|tim_min:tim_max},$dim$|:,lat|:,lon|0)
;            ta!1="plev"
         else
            ta = A({time|tim_min:tim_max},plev|:,lat|:)
;            ta = A({time|tim_min:tim_max},$dim$|:,lat|:)
;            ta!1="plev"
            uniform_coordinates(ta,field)
         end if
      end if
      plev = ta&plev
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
      end if
      plev = coord_list("plev",plev)
      ta&plev=plev
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

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.

      if ((field.eq."T3M")) then
         ta_mean = ta
      end if

      if ((field.eq."T2Mz")) then
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
      ;                      SEASONAL MEAN
      ;#########################################################
      ;calcolo le media stagionali
      ;if monthly
      v_seas = month_to_seasonN(ta_mean, season)

      ;                      CLIMATOLOGY
      ;#########################################################
      ; ta_mean is the monthly mean

      v_Clim = Climatology_CCMVal_mean(v_seas,v_seas&time,field)
      v_Clim_sdt = Climatology_CCMVal_std(v_seas,v_seas&time,field)
   else
      v_seas=month_to_seasonN(A,season)
      v_Clim=dim_avg_Wrap(v_seas(season|:,plev|:,lat|:,time|:))
      v_Clim_sdt=v_Clim
      v_Clim_sdt=v_Clim-v_Clim
      ;printVarSummary(v_Clim)
      ;v_Clim!0="month"
      ;v_Clim&month=(/1,2,3,4,5,6,7,8,9,10,11,12/)
      ;v_Clim_sdt!0="month"
      ;v_Clim_sdt&month=(/1,2,3,4,5,6,7,8,9,10,11,12/)
      plev = v_Clim&plev
      if(max(plev).gt.2000.) then
         plev=plev/100.
         print("Warning: level units was not hPa! automatic plev/100!")
      end if
      uniform_coordinates(v_Clim,field)
      uniform_coordinates(v_Clim_sdt,field)
      plev=v_Clim&plev
      plev = coord_list("plev",plev)
      v_Clim&plev=plev
      v_Clim_sdt&plev=plev
      lat = v_Clim&lat
      lat = coord_list("lat",lat)
      v_Clim&lat=lat
      v_Clim_sdt&lat=lat
   end if
   ;con questa funzione qualunque sia il field ottengo la climatologia cioe# un file medio di 12 mesi se l#input e# mensile o di 365 se l#input e# dayli.
   ;printVarSummary(v_Clim)

   ;******************************************************************
   v_wgt = new((/2,dimsizes(season),dimsizes(plev)/),typeof(v_Clim))
   ;################################################################
   ;    LAT AND LON WEIGHT MEAN  for the 3D field
   if (field.eq."T3M") then
      do nlev=0,dimsizes(plev)-1
         do ntime=0,dimsizes(season)-1
            v_wgt(0,ntime,nlev) = weight_area_latlon(v_Clim(season|ntime,plev|nlev,lat|:,lon|:),lat,lon,lat_min,lat_max)
            v_wgt(1,ntime,nlev) = weight_area_latlon(v_Clim_sdt(season|ntime,plev|nlev,lat|:,lon|:),lat,lon,lat_min,lat_max)
         end do
      end do
   end if

   ;##################################################################
   ;      LAT WEIGHT  MEAN     for the zonal field
   if ((field.eq."T2Mz").or.(field.eq."C2Mz")) then
      do nlev=0,dimsizes(plev)-1
         do ntime=0,dimsizes(season)-1
            ;printVarSummary(v_Clim)
            v_wgt(0,ntime,nlev) = weight_area_lat(v_Clim(season|ntime,plev|nlev,lat|:),lat,lat_min,lat_max)
            v_wgt(1,ntime,nlev) = weight_area_lat(v_Clim_sdt(season|ntime,plev|nlev,lat|:),lat,lat_min,lat_max)
         end do
      end do
   end if

   v_wgt!0 = "case"
   v_wgt!1 = "season"
   v_wgt!2 = "plev"

   v_wgt&case = (/"mean","sdt"/)
   v_wgt&season = season
   v_wgt&plev = plev

   v_wgt@long_name=MODEL+" Climatological_mean"

   return(v_wgt)
end
