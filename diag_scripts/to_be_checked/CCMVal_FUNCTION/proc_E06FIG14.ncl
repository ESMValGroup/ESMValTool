;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 14 of                                           #
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
;# The NCL script is called by the NCL shell script E06FIG14.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure14")
function figure14(A:numeric,field:string,MODEL:string,year_min:numeric,year_max:numeric)

local field,season,MODEL,lat_max,lat_min,ta,time,season_n,lat,\
lon,plev,ta_mean,v_Clim,v_Clim_sdt,nlev,ntime,v_wgtmean,\
v_wgtmean_sdt,year_min,year_max,tim_min,tim_max
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
      print("WARNING: Year range specified "+year_min+"-"+year_max+" is not included in the time series ("+dum_ut(0)+"-"+dum_ut(1)+")!!! ")
      delete(dum_t)
      delete(dum_ut)   
   end if

   ;*********************************************************
   if((field.eq."T2Ds").or.(field.eq."T2Ms")) then
;      print(A&time)
;      print(tim_min)
;      print(tim_max)
;      printVarSummary(A)
      if (dimsizes(dimsizes(A)).gt.3) then
        ta = A({time|tim_min:tim_max},plev|0,lat|:,lon|:)
      else
        ta = A({time|tim_min:tim_max},lat|:,lon|:)
      end if
      uniform_coordinates(ta,field)
      lon = ta&lon
      lon = coord_list("lon",lon)
      lat = ta&lat
      lat = coord_list("lat",lat)

      do i=0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            print("error reading latitude, some missing values")
            lat(i)=0.
         end if
      end do

      time = ta&time
      time = coord_list("time",time)

   end if

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.

   if ((field.eq."T2Ds")) then
      Timi = cd_calendar(time, 0)

      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))

      n_dimt = (ye_max-ye_min+1)*12

      ta_mean = new((/n_dimt,dimsizes(lat),dimsizes(lon)/),typeof(ta1))

      do i_lon=0,dimsizes(lon)-1
         do j_lat=0,dimsizes(lat)-1
            ta_mean(:,j_lat,i_lon) = monthly_mean((/ta1(time|:,lat|j_lat,lon|i_lon)/),time)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      ta_mean!1 = "lat"
      ta_mean!2 = "lon"
      ta_mean&lat = ta1@lat
      ta_mean&lon = ta1@lon
   else
      ta_mean = ta
   end if

   ;##################################################################
   ;                      CLIMATOLOGY
   ;##################################################################
   ; ta_mean is the monthly mean
   v_Clim = Climatology_CCMVal_mean(ta_mean,time,field)
   ;v_Clim_sdt = Climatology_CCMVal_std(ta_mean,time,field)
   ;con questa funzione qualunque sia il field ottengo la climatologia cioe# un file medio di 12 mesi se l#input e# mensile o di 365 se l#input e# dayli.
   ;###############################################################

   ;******************************************************************
   ;##################################################################
   ;      ZONAL MEAN     for the zonal field

   v_wgt = dim_avg_Wrap(v_Clim(month|:,lat|:,lon|:))

   v_wgt!0 = "month"
   v_wgt!1 = "lat"
   month = (/1,2,3,4,5,6,7,8,9,10,11,12/)

   v_wgt&month = month
   v_wgt&lat = ta&lat

   v_wgt@long_name = MODEL+" Climatological_mean"

   return(v_wgt)
end
;##################################################################################
