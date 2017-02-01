;##################################################################################
;# This NCL script Calculates long term monthly means (monthly mean climatology)  #    
;##################################################################################
;##################################################################################
;# Author:                                                                        #
;#                                                                                #
;# Irene Cionni (i.cionni@niwa.co.nz)                                             #
;##################################################################################
;# The NCL script is called by the NCL shell scripts E06FIGXX.ncl                 #
;##################################################################################
; At  the begin of the main programm must be called
;   load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;   load "/cydrive/f/Programmi/ncl_program/FIGURE1.ncl/new_function.ncl" 
 
                
undef("Climatology_CCMVal_mean")
function  Climatology_CCMVal_mean (v:numeric,time:numeric,field:string)
local v,field,time,v_Clim,v_Clim_ord,TIME, year, month,yyyymm,yyyyddd, time1,time2,last_month,v1,lat,lon,plev,day
begin
   ; time@calendar="standard"

   ;#############################################################
   ;##################montly fields##############################
   TIME    = cd_calendar(time(:), 0)          ; type float

   year    = floattointeger( TIME(:,0) )
   month   = floattointeger( TIME(:,1) )
   yyyymm  = year*100+month

   if (all(ismissing(TIME)))then
      print("Error!! Please check the units of the time coordinates")
      exit
   end if

   if (v!0.eq."season") then

      if ((field.eq."T2Ms").or.(field.eq."T2Ds")) then
         v_Clim = dim_avg_Wrap(v(season|:,lat|:,lon|:,time|:))
         v_Clim_ord = v_Clim(season|:,lat|:,lon|:)
      end if
      if ((field.eq."T3M").or.(field.eq."T3Ds")) then
         v_Clim = dim_avg_Wrap(v(season|:,plev|:,lat|:,lon|:,time|:))
         v_Clim_ord = v_Clim(season|:,plev|:,lat|:,lon|:)
      end if
      if ((field.eq."T2Mz").or.(field.eq."T2Dz")) then
         v_Clim = dim_avg_Wrap(v(season|:,plev|:,lat|:,time|:))
         v_Clim_ord = v_Clim(season|:,plev|:,lat|:)
      end if
      if (field.eq."T1Ms") then
         v_Clim = dim_avg_Wrap(v(season|:,lat|:,time|:))
         v_Clim_ord = v_Clim(season|:,lat|:)
      end if

   else

      ;#############################################################

      ;series must begin on January and end on December
      if (month(0).gt.1) then                    ;if the first months is not january
         time1 = new((dimsizes(year)+month(0)-1),"integer")   ;    new time variable dimension as if it begin in January
         do k=1,month(0)-1                                                ;the first months that in the series are missing
            time1(k-1) = year(0)*100+k
         end do
         time1(month(0):dimsizes(year)+month(0)-1) = year(:)*100+month(:)
      else
         time1 = new((dimsizes(year)+month(0)-1),"integer")   ;if  the is January time1 is as time
         time1 = year(:)*100+month(:)
      end if

      last_month = time1(dimsizes(time1)-1)-(year(dimsizes(year)-1)*100) ;last months of the series

      if (last_month.lt.12) then      ; if the last month is not december then
         time2 = new((dimsizes(time1)+12-last_month),"integer")
         do k=1,12-last_month(0)
            time2(dimsizes(time1)-1+k) = time1(dimsizes(time1)-1)+k
         end do
         time2(0:dimsizes(time1)-1) = time1
      else                                                  ;if the last month is December
         time2 = time1
      end if

      ;******************************************************************
      ;*******************************************************************
      ; extend the variable to a series that begin in January and end at December, using Fill_Value

      if (field.eq."T2Ms") then                        ; That field is (month,lat,lon) I need to now to do the correct
         v1 = new((/dimsizes(time2),dimsizes(v&lat),dimsizes(v&lon)/),typeof(v))
         v1!0 = "time"
         v1!1 = "lat"
         v1!2 = "lon"
         v1&time = time2

         ;;printVarSummary(v)

         v1({time|yyyymm},lat|:,lon|:) = (/v(time|:,lat|:,lon|:)/)
         v1&lat = v&lat
         v1&lon = v&lon

         v_Clim= clmMonLLT(v1(lat|:,lon|:,time|:))
         v_Clim_ord=v_Clim(month|:,lat|:,lon|:)
      end if

      if (field.eq."T3M") then
         v1 = new((/dimsizes(time2),dimsizes(v&plev),dimsizes(v&lat),dimsizes(v&lon)/),typeof(v))
         v1!1 = "plev"
         v1!2 = "lat"
         v1!3 = "lon"
         v1!0 = "time"
         v1&time = time2

         v1({time|yyyymm},plev|:,lat|:,lon|:) = (/v(time|:,plev|:,lat|:,lon|:)/)
         v1!1 = "plev"
         v1!2 = "lat"
         v1!3 = "lon"
         v1&lat = v&lat
         v1&lon = v&lon
         v1&plev = v&plev

         v_Clim = clmMonLLLT(v1(plev,lat,lon,time))
         v_Clim_ord = v_Clim(month|:,plev|:,lat|:,lon|:)
      end if

      if (field.eq."T2Mz") then
         v1 = new((/dimsizes(time2),dimsizes(v&plev),dimsizes(v&lat)/),typeof(v),v@_FillValue)
         v1!1 = "plev"
         v1!2 = "lat"
         v1!0 = "time"
         v1&time = time2

         v1({time|yyyymm},plev|:,lat|:) = (/v(time|:,plev|:,lat|:)/)
         v1!1 = "plev"
         v1!2 = "lat"
         v1&lat = v&lat
         v1&plev = v&plev

         v_Clim = clmMonLevLatT(v1(plev|:,lat|:,time|:))
         v_Clim_ord = v_Clim(month|:,plev|:,lat|:)
      end if

      if (field.eq."T1Ms") then
         v1 = new((/dimsizes(time2),dimsizes(v&lat)/),typeof(v))
         v1!0 = "time"
         v1!1 = "lat"
         v1&time = time2

         v1({time|yyyymm},lat|:) = (/v(time|:,lat|:)/)
         v1&lat = v&lat

         v_Clim = clmMonLatT(v1(lat|:,time|:))
         v_Clim_ord = v_Clim(month|:,lat|:)
      end if

      if (field.eq."T0") then
         v1 = new((/dimsizes(time2)/),typeof(v))
         v1!0="time"
         v1&time = time2

         v1({time|yyyymm})= (/v(time|:)/)

         v_Clim = clmMonT(v1(time|:))
         v_Clim_ord = v_Clim(month|:)
      end if

      ;##############################################################
      ;##################daily fields################################
      TIME    = cd_calendar(time, 0)          ; type float
      year    = floattointeger( TIME(:,0) )
      month   = floattointeger( TIME(:,1) )
      day     = floattointeger( TIME(:,2) )
      ddd     = day_of_year(year, month, day)
      yyyyddd = year*1000 + ddd
      ;#############################################################

      if ((field.eq."T0I")) then
         v&time = yyyyddd
         v_new = new((/dimsizes(v),1,1/),typeof(v))
         v_new(:,:,:) = (/v(:)/)
         v_new!0 = "time"
         v_new&time = v&time
         v_new!1 = "lat"
         v_new!2 = "lon"

         v_Clim = clmDayTLL(v_new(time|:,lat|:,lon|:), yyyyddd)
         ;printVarSummary(v_Clim)
         v_Clim_ord = (/v_Clim(:,0,0)/)
         v_Clim_ord!0 = "year_day"
         v_Clim_ord&year_day = v_Clim&year_day
      end if

      if ((field.eq."T2Ds")) then
         v&time = yyyyddd
         v_Clim = clmDayTLL(v(time|:,lat|:,lon|:),yyyyddd)
         v_Clim_ord = v_Clim(day|:,lat|:,lon|:)
      end if

      if (field.eq."T3D") then
         v&time = yyyyddd
         v_Clim = clmDayTLevLatLon(v(time|:,plev|:,lat|:,lon|:), yyyyddd)
         v_Clim_ord = v_Clim(day|:,plev|:,lat|:,lon|:)
      end if

      if ((field.eq."T2Dz")) then
         v&time = yyyyddd
         v_Clim = clmDayTLevLat(v(time|:,plev|:,lat|:), yyyyddd)
         v_Clim_ord = v_Clim(day|:,plev|:,lat|:)
      end if
   end if
   return (v_Clim_ord)
end


;######################################################################
;**********************************************************************
;######################################################################
undef("Climatology_CCMVal_std")
function  Climatology_CCMVal_std(v:numeric,time:numeric,field:string)
local v,field,v_Clim_std, v_Clim_std, TIME, year, month,yyyymm,yyyyddd, time1,time2,last_month,v1,lat,lon,plev,day
begin
   ;;printVarSummary(v)

   ;#############################################################
   ;##################montly fields##############################
   TIME    = cd_calendar(time(:), 0)          ; type float
   year    = floattointeger( TIME(:,0) )
   month   = floattointeger( TIME(:,1) )
   yyyymm  = year*100+month
   ;;printVarSummary(year)
   ;#############################################################

   if (v!0.eq."season") then

      if ((field.eq."T2Ms").or.(field.eq."T2Ds")) then
         v_Clim_std = dim_stddev_Wrap(v(season|:,lat|:,lon|:,time|:))
         v_Clim_std_ord=v_Clim_std(season|:,lat|:,lon|:)
      end if
      if ((field.eq."T3M").or.(field.eq."T3Ds")) then
         v_Clim_std = dim_stddev_Wrap(v(season|:,plev|:,lat|:,lon|:,time|:))
         v_Clim_std_ord=v_Clim_std(season|:,plev|:,lat|:,lon|:)
      end if
      if ((field.eq."T2Mz").or.(field.eq."T2Dz")) then
         v_Clim_std = dim_stddev_Wrap(v(season|:,plev|:,lat|:,time|:))
         v_Clim_std_ord=v_Clim_std(season|:,plev|:,lat|:)
      end if
      if (field.eq."T1Ms") then
         v_Clim_std = dim_stddev_Wrap(v(season|:,lat|:,time|:))
         v_Clim_std_ord=v_Clim_std(season|:,time|:,lat|:)
      end if

   else
      ;series must begin on January and end on December

      if (month(0).gt.1) then                    ;if the first months is not january
         time1 = new((dimsizes(year)+month(0)-1),"integer")   ;    new time variable dimension as if it begin in January
         do k=1,month(0)-1                                                ;the first months that in the series are missing
            time1(k) = year(0)*100+k
         end do
         time1(month(0):dimsizes(year)+month(0)-1) = year(:)*100+month(:)
      else
         time1 = new((dimsizes(year)+month(0)-1),"integer")   ;if  the is January time1 is as time
         time1 = year(:)*100+month(:)
      end if

      last_month=time1(dimsizes(time1)-1)-(year(dimsizes(year)-1)*100) ;last months of the series
      ;;print(last_month)
      if (last_month.lt.12) then      ; if the last month is not december then
         time2 = new((dimsizes(time1)+12-last_month),"integer")
         do k=1,12-last_month(0)
            time2(dimsizes(time1)-1+k) = time1+k
         end do
         time2(0:dimsizes(time1)-1) = time1
      else                                                  ;if the last month is December
         time2 = time1
      end if

      ;;printVarSummary(time2)

      ;******************************************************************+
      ;*******************************************************************
      ; extend the variable to a series that begin in January and end at December, using Fill_Value

      if (field.eq."T2Ms") then                        ; That field is (month,lat,lon) I need to now to do the correct
         v1 = new((/dimsizes(time2),dimsizes(v&lat),dimsizes(v&lon)/),typeof(v))
         v1!0 = "time"
         v1!1 = "lat"
         v1!2 = "lon"
         v1&time = time2

         v1({time|yyyymm},lat|:,lon|:) = v(time|:,lat|:,lon|:)
         v1&lat = v&lat
         v1&lon = v&lon

         v_Clim_std = stdMonLLT(v1(lat|:,lon|:,time|:))
         v_Clim_std_ord = v_Clim_std(month|:,lat|:,lon|:)
      end if

      if (field.eq."T3M") then
         v1 = new((/dimsizes(time2),dimsizes(v&plev),dimsizes(v&lat),dimsizes(v&lon)/),typeof(v))
         v1!0 = "time"
         v1!1 = "plev"
         v1!2 = "lat"
         v1!3 = "lon"
         v1&time = time2

         v1({time|yyyymm},plev|:,lat|:,lon|:) = v(time|:,plev|:,lat|:,lon|:)
         v1&lat = v&lat
         v1&lon = v&lon
         v1&plev = v&plev

         v_Clim_std = stdMonLLLT(v1(plev,lat,lon,time))
         v_Clim_std_ord = v_Clim_std(month|:,plev|:,lat|:,lon|:)
      end if

      if (field.eq."T2Mz") then
         v1 = new((/dimsizes(time2),dimsizes(v&plev),dimsizes(v&lat)/),typeof(v))
         v1!0 = "time"
         v1!1 = "plev"
         v1!2 = "lat"
         v1&time = time2

         v1({time|yyyymm},plev|:,lat|:) = (/v(time|:,plev|:,lat|:)/)
         v1&lat = v&lat
         v1&plev = v&plev

         v_Clim_std = stdMonLevLatT(v1(plev|:,lat|:,time|:))
         v_Clim_std_ord = v_Clim_std(month|:,plev|:,lat|:)
         ;;printVarSummary(v_Clim_std)
         ;;printVarSummary(v_Clim_std_ord)
      end if

      if (field.eq."T1Ms") then
         v1 = new((/dimsizes(time2),dimsizes(v&lat)/),typeof(v))
         v1!0 = "time"
         v1!1 = "lat"
         v1&time = time2

         v1({time|yyyymm},lat|:) = v(time|:,lat|:)
         v1&lat = v&lat

         v_Clim_std = stdMonLatT(v1(lat|:,time|:))
         v_Clim_std_ord = v_Clim_std(month|:,lat|:)
      end if

      ;##############################################################
      ;##################daily fields################################
      TIME    = cd_calendar(time, 0)          ; type float
      year    = floattointeger( TIME(:,0) )
      month   = floattointeger( TIME(:,1) )
      day     = floattointeger( TIME(:,2) )
      ddd     = day_of_year(year, month, day)
      yyyyddd = year*1000 + ddd
     ;#############################################################

      if ((field.eq."T2Ds")) then
          v&time = yyyyddd
          v_Clim_std = stdDayTLL(v(day|:,lat|:,lon|:), yyyyddd)
          v_Clim_std_ord = v_Clim(day|:,lat|:,lon|:)
      end if

      if (field.eq."T3D") then
         v&time = yyyyddd
         v_Clim_std = stdDayTLevLatLon(v(time|:,plev|:,lat|:,lon|:), yyyyddd)
         v_Clim_std_ord = v_Clim(day|:,plev|:,lat|:,lon|:)
      end if

      if ((field.eq."T2Dz")) then
         v&time = yyyyddd
         v_Clim_std = sdtDayTLevLat(v(time|:,plev|:,lat|:), yyyyddd)
         v_Clim_std_ord = v_Clim(day|:,plev|:,lat|:)
      end if
   end if

   return (v_Clim_std_ord)

end
;##################################################################################
