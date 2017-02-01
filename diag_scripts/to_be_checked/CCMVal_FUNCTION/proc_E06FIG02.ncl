;##################################################################################
;# This NCL script processes CF compliant netCDF CCMVal model output as specified #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and is used to produces FIGURE 2 of                                            #
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
;# The NCL script is called by the NCL shell script E06FIG02.ncl                  #
;##################################################################################
load "./reformat_scripts/coord_function.ncl"
load "./plot_type/uniform_coordinates.ncl"

undef("figure2")
function figure2(A:numeric,field:string,MODEL:string,lat_sel:numeric,year_min:numeric,year_max:numeric)

local field,MODEL,lat_sel,ua,time,lat,\
lon,plev,ua_mean,v_Clim,v_Clim_sdt,v_lev,v_lev_sdt,nlev,ntime,k,\
year_min,year_max,tim_min,tim_max,v_sel,v_sel_sdt,indice_lev,indice_lev_p,indice_lev_m,v_sel_p,v_sel_m,d_lat,ind_lat

begin
   tim = A&time
   tim = coord_list("time",tim)
   option=0
   option@calendar=tim@calendar

   tim_min = doubletofloat(cd_inv_calendar(year_min,1,1,0,0,0,tim@units,option))
   tim_max = doubletofloat(cd_inv_calendar(year_max,12,31,23,59,59,tim@units,option))

   if ((tim_min+31).lt.(tim(0))).or.((tim_max-31).gt.tim(dimsizes(tim)-1)) then
      dum_t=(/tim(0),tim(dimsizes(tim)-1)/)
      copy_VarAtts(tim,dum_t)
      dum_ut = cd_calendar(dum_t,-1)
      print(MODEL+" ["+dum_ut(0)+"-"+dum_ut(1)+"]  WARNING: Year range specified "+year_min+"-"+year_max+" is not included in the time series!!! ")
      delete(dum_t)
      delete(dum_ut)
   end if

   ;*********************************************************
   if((field.eq."T3D").or.(field.eq."T3M")) then

      ua = A({time|tim_min:tim_max},plev|:,lat|:,lon|:)
      uniform_coordinates(ua,field)
      lon = ua&lon
      lon = coord_list("lon",lon)

   end if

   ;*********************************************************
   if((field.eq."T2Dz").or.(field.eq."T2Mz")) then
      ua = A({time|tim_min:tim_max},plev|:,lat|:)
      uniform_coordinates(ua,field)
   end if

   lat = ua&lat
   lat = coord_list("lat",lat)
   plev = ua&plev
   plev = coord_list("plev",plev)
   if(max(plev).gt.2000.) then   	
      plev=plev/100.  
      print("Warning: level units was not hPa! automatic plev/100!")
      ua&plev=plev
   end if
   
   time = ua&time
   time = coord_list("time",time)

   do i=0,dimsizes(lat)-1
      if (lat(i).eq.-9999) then
         ;print("error reading latitude, some missing values")
         lat(i) = 0.
      end if
   end do
   ua&lat = lat

   ;#########################################################
   ;                      MONTHLY MEAN
   ;#########################################################
   ; the field must to be montlhy to do season mean.

   if ((field.eq."T3M")) then
      ua_mean = ua
   end if

   if ((field.eq."T2Mz")) then
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
      ua_mean&plev = ua@plev
      ua_mean&lat = ua@lat
      ua_mean&lon = ua@lon

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
      ua_mean&plev = ua@plev
      ua_mean&lat = ua@lat

      field = "T2Mz"
    end if

   ;#########################################################
   ;                      CLIMATOLOGY
   ;#########################################################
   ; ua_mean is the monthly mean
   v_Clim = Climatology_CCMVal_mean(ua_mean,time,field)
   v_Clim_sdt = Climatology_CCMVal_std(ua_mean,time,field)
   ;con quesua funzione qualunque sia il field ottengo la climatologia cioe# un file medio di 12 mesi se l#input e# mensile o di 365 se l#input e# dayli.
   ;###############################################################

   ;****************************************************************
   ;################################################################
   ; LAT SEL
   if (field.eq."T3M") then
      v_Clim_zon = dim_avg_Wrap(v_Clim)
      v_Clim_sdt_zon = dim_avg_Wrap(v_Clim_sdt)
      delete(v_Clim)
      delete(v_Clim_sdt)
      v_Clim = v_Clim_zon
      v_Clim_sdt = v_Clim_sdt_zon
      field = "T2Mz"
   end if

   ; LAT SEL for the zonal field
   if (field.eq."T2Mz") then
      lat = v_Clim&lat
      d_lat = abs(lat(0)-lat(1))

      if (lat_sel.eq.(90.)) then
         ind_lat = ind((lat.ge.(lat_sel-d_lat)).and.(lat.le.lat_sel))
         if(dimsizes(ind_lat).eq.1) then
            v_sel = v_Clim(month|:,plev|:,lat|ind_lat)
            v_sel_sdt = v_Clim_sdt(month|:,plev|:,lat|ind_lat)
         else
            v_sel = dim_avg_Wrap(v_Clim(month|:,plev|:,lat|ind_lat))
            v_sel_sdt = dim_avg_Wrap(v_Clim_sdt(month|:,plev|:,lat|ind_lat))
         end if
      end if

      if (lat_sel.eq.(-90.)) then
         ind_lat = ind((lat.le.(lat_sel+d_lat)).and.(lat.ge.lat_sel))
         if(dimsizes(ind_lat).eq.1) then
            v_sel = v_Clim(month|:,plev|:,lat|ind_lat)
            v_sel_sdt = v_Clim_sdt(month|:,plev|:,lat|ind_lat)
         else
            v_sel = dim_avg_Wrap(v_Clim(month|:,plev|:,lat|ind_lat))
            v_sel_sdt = dim_avg_Wrap(v_Clim_sdt(month|:,plev|:,lat|ind_lat))
         end if
      end if

      if ((lat_sel.lt.90).and.(lat_sel.ge.0)) then
         ind_lat = ind((lat.le.(lat_sel+d_lat)).and.(lat.ge.(lat_sel-d_lat)))
         if(dimsizes(ind_lat).eq.1) then
            v_sel = v_Clim(month|:,plev|:,lat|ind_lat)
            v_sel_sdt = v_Clim_sdt(month|:,plev|:,lat|ind_lat)
         else
            v_sel = dim_avg_Wrap(v_Clim(month|:,plev|:,lat|ind_lat))
            v_sel_sdt = dim_avg_Wrap(v_Clim_sdt(month|:,plev|:,lat|ind_lat))
         end if
      end if

      if ((lat_sel.gt.-90).and.(lat_sel.le.0)) then
         ind_lat = ind((lat.le.(lat_sel+d_lat)).and.(lat.ge.(lat_sel-d_lat)))
         if(dimsizes(ind_lat).eq.1) then
            v_sel = v_Clim(month|:,plev|:,lat|ind_lat)
            v_sel_sdt = v_Clim_sdt(month|:,plev|:,lat|ind_lat)
         else
            v_sel = dim_avg_Wrap(v_Clim(month|:,plev|:,lat|ind_lat))
            v_sel_sdt = dim_avg_Wrap(v_Clim_sdt(month|:,plev|:,lat|ind_lat))
         end if
      end if
   end if

   ;##################################################################
   ; FIND the zero
   v_sel_p = v_sel
   v_sel_p = v_sel+v_sel_sdt

   v_sel_m = v_sel
   v_sel_m = v_sel-v_sel_sdt

   v_lev = new((/3,12/),typeof(plev))

   do k=0,11
      ; max_lev=ind(max(plev))
      ; min_lev=ind(min(plev))

      indice_lev = ind((v_sel(month|k,plev|:).ge.0.))
      indice_lev_p = ind((v_sel_p(month|k,plev|:).ge.0.))
      indice_lev_m = ind((v_sel_m(month|k,plev|:).ge.0.))
      
      if (plev(0).le.plev(1)) then
      
         if (.not.all(ismissing(indice_lev))) then
            v_lev(0,k) = min(plev((indice_lev)))
         end if
         if (.not.all(ismissing(indice_lev_m))) then
            v_lev(1,k) = min(plev((indice_lev_m)))
         end if
         if (.not.all(ismissing(indice_lev_p))) then
            v_lev(2,k) = min(plev((indice_lev_p)))
         end if
      else

         if (.not.all(ismissing(indice_lev))) then
            v_lev(0,k) = min(plev((indice_lev)))
         end if
         if (.not.all(ismissing(indice_lev_m))) then
            v_lev(1,k) = min(plev((indice_lev_m)))
         end if
         if (.not.all(ismissing(indice_lev_p))) then
            v_lev(2,k) = min(plev((indice_lev_p)))
         end if
      end if

      delete(indice_lev)
      delete(indice_lev_p)
      delete(indice_lev_m)
   end do

   v_lev!0 = "case"
   v_lev!1 = "month"
   v_lev&case = (/"climatological_mean","minus_sdt","plus_sdt"/)
   v_lev&month = v_sel&month
   
   return(v_lev)

end
;##################################################################################
