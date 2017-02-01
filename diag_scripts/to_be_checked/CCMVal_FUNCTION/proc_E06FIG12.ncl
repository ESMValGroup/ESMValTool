;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 12 of                                           #
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
;# The NCL script is called by the NCL shell script E06FIG12.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"
undef("figure12")
function figure12(A:numeric,field:string,sel_month:numeric,MODEL:string)

local field,season,MODEL,lat_max,lat_min,ta,time,season_n,lat,\
lon,plev,ta_mean,v_Clim,v_Clim_sdt,nlev,ntime,v_wgtmean,sel_month

begin
   tim = A&time
   tim = coord_list("time",tim)

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then

      ta = A(time|:,plev|:,lat|:)
      if ((dimsizes(A&lat)).gt.1) then
         uniform_coordinates(ta,field)
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
  
      
      do i=0,dimsizes(lat)-1
         if (lat(i).eq.-9999) then
            print("error reading latitude, some missing values")
            lat(i) = 0.
         end if
      end do

      time = ta&time
      time = coord_list("time",time)
   end if

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.

   if ((field.eq."T2Mz")) then
      ta_mean = ta
   end if

   if (field.eq."T2Dz") then
      Timi = cd_calendar(time, 0)
      ye_max = floattointeger(max(Timi(:,0)))
      ye_min = floattointeger(min(Timi(:,0)))
      n_dimt = (ye_max-ye_min+1)*12

      ta_mean = new((/n_dimt,dimsizes(plev),dimsizes(lat)/),typeof(ta))
      do i_lat=0,dimsizes(lat)-1
         do j_lev=0,dimsizes(plev)-1
            ta_mean(:,j_lev,i_lat) = (/monthly_mean((/ta(time|:,plev|j_lev,lat|i_lat)/),time)/)
         end do
      end do

      delete(Timi)
      delete(ye_max)
      delete(ye_min)
      delete(n_dimt)

      field = "T2Mz"
    end if

   ;##################################################################
   ;  Month series
   ;##################################################################

   timee = ta_mean&time
   timee = coord_list("time",timee)
   dates = cd_calendar(timee,0)

   yyyy = dates(:,0)
   mm = floattointeger(dates(:,1))
   dd = dates(:,2)

   ind_sel_month = ind(mm.eq.sel_month)
   v_sel = (ta_mean(time|ind_sel_month,plev|:,lat|:))
   v_sel!0 = "year"
   v_sel&year = yyyy(ind_sel_month)

   ;##################################################################
   ;##################################################################

   v_sel@long_name = MODEL+ta@long_name

   return(v_sel)
end
;##################################################################################
