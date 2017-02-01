;##################################################################################
;# This NCL script contains the new functions for the CCMVal Diagnostic Tool      #
;##################################################################################
;# Authors and Contact:                                                           #
;#                                                                                #
;# I. Cionni and H.Shiona NIWA                                                    #
;##################################################################################
;# The NCL script is called by the NCL shell scripts E06FIGXX.ncl                 #
;##################################################################################
;********************************1************************************
; Calculate long term monthly means (monthly climatology)

; The time dimension must be a multiple of 12
;
;    x(lat,lev,time)  <==== INPUT DIMENSION ORDER
;    x!0 = "lat"
;    x!1 = "lev"
;    x!2 = "time"
;
; Usage:  moClm = clmMonLevLatT (x)

undef("clmMonLevLatT")
function clmMonLevLatT (x[*][*][*]:numeric)
local dimx, nlat, nlev, ntim, i, nmo, nmos, aveMonth
begin

   dimx = dimsizes (x)
   nlat = dimx(0)
   nlev = dimx(1)
   ntim = dimx(2)

   nmos = 12
   modCheck ("clmMonLevLatT", ntim, nmos)   ; error check
;
; Compute all 12 monthly averages first.
;
   aveMonth = new((/nlat,nlev,nmos/),typeof(x),getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(:,:,nmo) = dim_avg( x(:,:,nmo:ntim-1:nmos) )
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,aveMonth)
   aveMonth@time_op_ncl = "Climatology: "+ (ntim/nmos) +" years"
   aveMonth@info = "function clmMonLevLatT: new_function.ncl"

; copy spatial (eg, lat/lon) coordinate variables

   do i=0,1
      if (.not.ismissing(x!i)) then
         aveMonth!i = x!i
         if (iscoord(x,x!i)) then
            aveMonth&$aveMonth!i$ = x&$x!i$
         end if
      end if
   end do

   aveMonth!2     = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var

   return (aveMonth)
end

;********************************2************************************
;*********************************************************************
; I. Cionni
; Calculate long term monthly means (monthly climatology)
;
; The time dimension must be a multiple of 12
;
;    x(lat,time)  <==== INPUT DIMENSION ORDER
;    x!0 = "lat"
;    x!1 = "time"
;
; Usage:  moClm = clmMonLatT (x)

undef("clmMonLatT")
function clmMonLatT (x[*][*]:numeric)
local dimx, nlat, ntim, i, nmo, nmos, aveMonth
begin

   dimx = dimsizes (x)
   nlat = dimx(0)
   ntim = dimx(1)

   nmos = 12
   modCheck ("clmMonLatT", ntim, nmos)   ; error check
;
; Compute all 12 monthly averages first.
;
   aveMonth = new((/nlat,nmos/),typeof(x),getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(:,nmo) = dim_avg( x(:,nmo:ntim-1:nmos) )
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,aveMonth)
   aveMonth@time_op_ncl = "Climatology: "+ (ntim/nmos) +" years"
   aveMonth@info = "function clmMonLatT: new_function.ncl"

; copy spatial (eg, lat/lon) coordinate variables

   do i=0,1
      if (.not.ismissing(x!i)) then
         aveMonth!i = x!i
         if (iscoord(x,x!i)) then
            aveMonth&$aveMonth!i$ = x&$x!i$
         end if
      end if
   end do

   aveMonth!2     = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var

   return (aveMonth)
end

;********************************3************************************
;*********************************************************************
undef("clmMonT")
function clmMonT (x[*]:numeric)
local dimx,  ntim, i, nmo, nmos, aveMonth
begin

   dimx = dimsizes (x)
   ntim = dimx(0)

   nmos = 12
   modCheck ("clmMonT", ntim, nmos)   ; error check
;
; Compute all 12 monthly averages first. 
;
   aveMonth = new((/nmos/),typeof(x),getFillValue(x))

   do nmo=0,nmos-1
      aveMonth(nmo) = dim_avg( x(nmo:ntim-1:nmos) )
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,aveMonth)
   aveMonth@time_op_ncl = "Climatology: "+ (ntim/nmos) +" years"
   aveMonth@info = "function clmMonT: new_function.ncl"

   aveMonth!0     = "month"           ; create a "month" named dim
   aveMonth&month = ispan(0,nmos-1,1) ; create a month coord var

   return (aveMonth)
end

;********************************4************************************
;*********************************************************************
undef("clmDayTLevLatLon")
function clmDayTLevLatLon (x[*][*][*][*]:numeric, yyyyddd:integer)   
;
; calculate the mean Annual Cycle from daily data. 
; The return array will gave the raw climatology at each grid point
;
;              x(time,plev,lat,lon)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "plev"
;              x!2 = "lat"
;              x!2 = "lon"
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
;    Leap
;              1908001 =>  Jan  1, 1908]
;              1908032 =>  Feb  1, 1908]
;              1908059 =>  Feb 28, 1908]
;              1908060 =>  Feb 29, 1908]
;              1908061 =>  Mar  1, 1908]
;              1908366 =>  Dec 31, 1908]
;
; Usage:  xClmDay = clmDAY_TLL (x, yyyyddd)
; -------

local dimx, ntim, nlat, mlon,nlev, ndys, days, clmDay, ndy, indx, year_day, nFill
begin

   dimx = dimsizes (x)

   ntim = dimx(0)
   nlev = dimx(1)
   nlat = dimx(2)
   mlon = dimx(3)
   ndys = 366                                     ; allow for leap year

   days = yyyyddd - (yyyyddd/1000)*1000           ; strip year info [yyyy]

   clmDay = new((/ndys,nlev,nlat,mlon/),typeof(x), getFillValue(x) ) ; daily climatology
;
; Compute averages for each sequebtial day of the year.
; This uses dimension swapping.
;
   do ndy=0,ndys-2                                 ; ndy=0->364 ==> day of year 1->365
      indx = ind( days.eq.(ndy+1) )
      clmDay(ndy,:,:,:) = dim_avg(x(lev|:,lat|:,lon|:,time|indx))
      delete(indx)
   end do
                                                   ; nominal day 366
                                                   ; ave(31 Dec + 1 Jan)=leap
   clmDay(ndys-1,:,:,:) = (clmDay(0,:,:,:) + clmDay(ndys-2,:,:,:))*0.5

   nFill = num(ismissing(clmDay))
   if (nFill.eq.0) then
      delete(clmDay@_FillValue)
   end if

   clmDay@long_name   = "Daily Climatology"
   if (isatt(x,"long_name")) then
      clmDay@long_name = clmDay@long_name +": "+x@long_name
   end if
   if (isatt(x,"units")) then
      clmDay@units    = x@units
   end if
   clmDay@information = "Raw daily averages across all years"
   clmDay@smoothing   = "None"

   year_day           = ispan(1,ndys,1)
   year_day@long_name = "day of year"
   year_day@units     = "ddd"

   clmDay!0           = "year_day"
   clmDay&year_day    =  year_day

   copy_VarCoords(x(0,:,:,:), clmDay(0,:,:,:))   ; trick
   delete(clmDay@year_day)                   ; clean up

   return (clmDay)
end

;********************************5************************************
;*********************************************************************
undef("clmDayTLevLat")
function clmDayTLevLat (x[*][*][*]:numeric, yyyyddd:integer)
;
; calculate the mean Annual Cycle from daily data.
; The return array will gave the raw climatology at each grid point
;
;              x(time,plev,lat,lon)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "plev"
;              x!2 = "lat"
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
;    Leap
;              1908001 =>  Jan  1, 1908]
;              1908032 =>  Feb  1, 1908]
;              1908059 =>  Feb 28, 1908]
;              1908060 =>  Feb 29, 1908]
;              1908061 =>  Mar  1, 1908]
;              1908366 =>  Dec 31, 1908]
;
; Usage:  xClmDay = clmDAY_TLL (x, yyyyddd)
; -------

local dimx, ntim, nlat,nlev, ndys, days, clmDay, ndy, indx, year_day, nFill
begin

   dimx = dimsizes (x)

   ntim = dimx(0)
   nlev = dimx(1)
   nlat = dimx(2)
   ndys = 366                                     ; allow for leap year

   days = yyyyddd - (yyyyddd/1000)*1000           ; strip year info [yyyy]

   clmDay = new((/ndys,nlev,nlat,mlon/),typeof(x), getFillValue(x) ) ; daily climatology
;
; Compute averages for each sequebtial day of the year.
; This uses dimension swapping.
;
   do ndy=0,ndys-2                                 ; ndy=0->364 ==> day of year 1->365
      indx = ind( days.eq.(ndy+1) )
      clmDay(ndy,:,:,:) = dim_avg(x(lev|:,lat|:,lon|:,time|indx))
      delete(indx)
   end do
                                                   ; nominal day 366
                                                   ; ave(31 Dec + 1 Jan)=leap
   clmDay(ndys-1,:,:,:) = (clmDay(0,:,:,:) + clmDay(ndys-2,:,:,:))*0.5

   nFill = num(ismissing(clmDay))
   if (nFill.eq.0) then
      delete(clmDay@_FillValue)
   end if

   clmDay@long_name   = "Daily Climatology"
   if (isatt(x,"long_name")) then
      clmDay@long_name = clmDay@long_name +": "+x@long_name
   end if
   if (isatt(x,"units")) then
      clmDay@units    = x@units
   end if
   clmDay@information = "Raw daily averages across all years"
   clmDay@smoothing   = "None"

   year_day           = ispan(1,ndys,1)
   year_day@long_name = "day of year"
   year_day@units     = "ddd"

   clmDay!0           = "year_day"
   clmDay&year_day    =  year_day

   copy_VarCoords(x(0,:,:,:), clmDay(0,:,:,:))   ; trick
   delete(clmDay@year_day)                   ; clean up

   return (clmDay)
end

;********************************6************************************
;*********************************************************************
; Modified from a D Shea function
; compute a user specified seasonal mean [all are three-month means]
; DJF,JFM,FMA,MAM,AMJ,MJJ,JJA,JAS,ASO,SON,OND,NDJ

; first (DJF=JF) /last (NDJ=ND) seasons are 2-month averages
;
; x(time,lat,lon),  x(time,lev,lat,lon)
;   ^^^^^^^^^^^^      ^^^^^^^^^^^^^^^^
; must have named dim BUT can be ANY names
;
; The input "x" are assumed to contain monthly mean data
; The size of "time" MUST be divisible by 12.
; Also, it is assumed the "Jan" is the 1st month.
;
; xMon(time)  or  xMon(time,lat,lon)   or  xMon(time,lev,lat,lon)
; USAGE     xJJA = month_to_season (xMon, "JJA")
;
; RESULT    xJJA(time/12,lev,lat,lon)    xJJA(time/12,lat,lon)
;
; Note: this returns (7/2003) NMO as an attribute

undef ("month_to_season1")
function month_to_season1 (xMon:numeric, SEASON:string)

local season,NMO,dimx,rank,ntim,nlat,mlon,nmos,nyrs,con \
    , nyrStrt,nyrLast,nyr,n,xSea, klev, dName,cv,xSea
begin

   season  =  (/"DJF","JFM","FMA","MAM","AMJ","MJJ" \
               ,"JJA","JAS","ASO","SON","OND","NDJ" /)

   NMO     = ind(season.eq.SEASON)  ; index corresponding to season
   if (ismissing(NMO)) then
      ;;;print ("contributed: month_to_season: bad season: SEASON="+SEASON)
   end if

   dimx    = dimsizes(xMon)
   rank    = dimsizes(dimx)
   if (rank.eq.2 .or. rank.ge.5) then
      ;;;print ("contributed: month_to_season: rank="+rank)
      ;;;print ("----- rank currently not handled -----")
   end if

   nmos    = 12
   ntim    = dimx(0)
   modCheck ("month_to_season", ntim, nmos)

   if (rank.ge.3) then
      nlat = dimx(rank-2)
      mlon = dimx(rank-1)
   end if
   nyrs    = ntim/nmos
   con     = 1./3.

   nyrStrt = 0
   nyrLast = nyrs-1
   if (NMO.eq.0) then
      nyrStrt = 1
   end if
   if (NMO.eq.nmos-1) then
      nyrLast = nyrs-2
   end if

   if (rank.eq.1) then
      xSea = new ( nyrs, typeof(xMon), getFillValue(xMon))
      do nyr=nyrStrt,nyrLast
         n = nyr*nmos + NMO
         xSea(nyr) = (xMon(n-1) + xMon(n) + xMon(n+1))*con
      end do
                                        ; special for beginning/end points
      if (NMO.eq.0) then
         n = 0
         xSea(0) = (xMon(n) + xMon(n+1))*0.5
      end if
      if (NMO.eq.nmos-1) then
         n = (nyrs-1)*nmos + NMO
         xSea(nyrs-1) = (xMon(n) + xMon(n-1))*0.5
      end if
   end if

   if (rank.eq.3) then
      xSea = new ( (/nyrs,nlat,mlon/), typeof(xMon), getFillValue(xMon))
      do nyr=nyrStrt,nyrLast
         n = nyr*nmos + NMO
         xSea(nyr,:,:) = (xMon(n-1,:,:) + xMon(n,:,:) + xMon(n+1,:,:))*con
      end do
                                        ; special for beginning/end points
      if (NMO.eq.0) then
         n = 0
         xSea(0,:,:) = (xMon(n,:,:) + xMon(n+1,:,:))*0.5
      end if
      if (NMO.eq.nmos-1) then
         n = (nyrs-1)*nmos + NMO
         xSea(nyrs-1,:,:) = (xMon(n,:,:) + xMon(n-1,:,:))*0.5
      end if
   end if

   if (rank.eq.4) then
      klev = dimx(1)
      xSea = new ( (/nyrs,klev,nlat,mlon/), typeof(xMon), getFillValue(xMon))
      do nyr=nyrStrt,nyrLast
         n = nyr*nmos + NMO
         xSea(nyr,:,:,:) = (xMon(n-1,:,:,:) + xMon( n ,:,:,:) \
                                            + xMon(n+1,:,:,:))*0.33333
      end do

      if (NMO.eq.0) then
         n = 0
         xSea(0,:,:,:) = (xMon(n,:,:,:) + xMon(n+1,:,:,:))*0.5
      end if
      if (NMO.eq.nmos-1) then
         n = (nyrs-1)*nmos + NMO
         xSea(nyrs-1,:,:,:) = (xMon(n,:,:,:) + xMon(n-1,:,:,:))*0.5
      end if
   end if

   copy_VarAtts (xMon, xSea)
   if (isatt(xMon,"long_name") .or. isatt(xMon,"description") .or. \
      isatt(xMon,"standard_name") ) then
      xSea@long_name = SEASON+": "+getLongName(xMon)
   end if

   do n=1,rank-1                  ; copy spatial coordinates
      if (.not.ismissing(xMon!n)) then
         xSea!n = xMon!n
         if(iscoord(xMon,xMon!n))
            xSea&$xSea!n$ = xMon&$xMon!n$
         end if
      end if
   end do

   dName        = xMon!0
   xSea!0       = dName

   if(iscoord(xMon,dName)) then
      cv = xMon&$dName$(NMO:ntim-1:nmos)
      xSea&$dName$ = cv
   end if

   xSea@NMO = NMO   ; for possible use in subscripting
                    ; eg: nStrt= xSea@NMO ; time(nStrt:ntim-1,12)
   return (xSea)
end

;********************************7************************************
;*********************************************************************
; Modified from  D. Shea function
; Calculate standard deviations of monthly means 
;
; The time dimension must be a multiple of 12
;
;    x(lat,lon,time)  <==== INPUT DIMENSION ORDER
;    x!0 = "lat"   
;    x!1 = "lon"
;    x!2 = "time"
;
; Usage:  moStd = stdMonLatLonT (x)

undef("stdMonLatLonT")
function stdMonLatLonT (x[*][*][*]:numeric)
local dimx, nlat, mlon, ntim, i, nmo, nmos, stdMonth
begin

   dimx  = dimsizes (x)
   nlat  = dimx(0)
   mlon  = dimx(1)
   ntim  = dimx(2)

   nmos  = 12
   modCheck ("stdMonLatLonT", ntim, nmos)   ; error check
;
; Compute all 12 monthly standard deviations first.
;
   stdMonth = new((/nlat,mlon,nmos/),typeof(x) \
                                    ,getFillValue(x))

   do nmo=0,nmos-1
      stdMonth(:,:,nmo) = dim_stddev( x(:,:,nmo:ntim-1:nmos) )
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,stdMonth)

   stdMonth@time_op_ncl = " Monthly Standard Deviation: "+ (ntim/nmos) +" years"
   stdMonth@info    = "function stdMonLatLonT: contributed.ncl"

; copy spatial (eg, lat/lon) coordinate variables

   do i=0,1
      if (.not.ismissing(x!i)) then
          stdMonth!i = x!i
          if (iscoord(x,x!i)) then
              stdMonth&$stdMonth!i$ = x&$x!i$
          end if
      end if
   end  do

   stdMonth!2     = "month"           ; create a "month" named dim
   stdMonth&month = ispan(0,nmos-1,1) ; create a month coord var

   return (stdMonth)
end

;********************************8************************************
;*********************************************************************
; Modified from a D. Shea
; Calculate standard deviations of monthly means
;
; The time dimension must be a multiple of 12
;
;    x(lev,lat,time)  <==== INPUT DIMENSION ORDER
;
; Usage:  moStd = stdMonLevLatT (x)
;         moStd(lev,lat,lon,12)

undef("stdMonLevLatT")
function stdMonLevLatT (x[*][*][*]:numeric)
local dimx, klvl, nlat, ntim, i, nmo, nmos, stdMonth
begin

   dimx  = dimsizes (x)
   klvl  = dimx(0)
   nlat  = dimx(1)
   ntim  = dimx(2)

   nmos  = 12
   modCheck ("stdMonLevLatT", ntim, nmos)   ; error check
;
; Compute all 12 monthly standard deviations first.
;
   stdMonth = new((/klvl,nlat,nmos/),typeof(x) \
                                         ,getFillValue(x))

   do nmo=0,nmos-1
      stdMonth(:,:,nmo) = dim_stddev( x(:,:,nmo:ntim-1:nmos) )
   end do

; copy attributes and add a couple of informational attributes
; only need to do this if plotting or writing to netCDF file

   copy_VarAtts (x,stdMonth)

   stdMonth@time_op_ncl = " Monthly Standard Deviation: "+ (ntim/nmos) +" years"
   stdMonth@info    = "function stdMonLevLatT: new_function.ncl"

; copy spatial (eg, lat/lon) coordinate variables

   do i=0,1
      if (.not.ismissing(x!i)) then
          stdMonth!i = x!i
          if (iscoord(x,x!i)) then
              stdMonth&$stdMonth!i$ = x&$x!i$
          end if
      end if
   end  do

   stdMonth!2     = "month"           ; create a "month" named dim
   stdMonth&month = ispan(0,nmos-1,1) ; create a month coord var

   return (stdMonth)
end

;********************************9************************************
;*********************************************************************
undef("clmDayTLevLat1")
function clmDayTLevLat1 (x[*][*][*]:numeric, yyyyddd:integer)
;
; calculate the mean Annual Cycle from daily data.
; The return array will gave the raw climatology at each grid point
;
;              x(time,lev,lat)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "lev"
;              x!2 = "lat"
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
;    Leap
;              1908001 =>  Jan  1, 1908]
;              1908032 =>  Feb  1, 1908]
;              1908059 =>  Feb 28, 1908]
;              1908060 =>  Feb 29, 1908]
;              1908061 =>  Mar  1, 1908]
;              1908366 =>  Dec 31, 1908]
;
; Usage:  xstdDay = clmDAY_TLL (x, yyyyddd)
; -------

local x, dimx, ntim, nlev, mlat, ndys, days, clmDay, ndy, indx, year_day, nFill
begin

   dimx  = dimsizes (x)

   ntim  = dimx(0)
   nlev  = dimx(1)
   mlat  = dimx(2)
   ndys  = 366                                     ; allow for leap year

   days  = yyyyddd - (yyyyddd/1000)*1000           ; strip year info [yyyy]

   clmDay= new((/ndys,nlev,mlat/),typeof(x), getFillValue(x) ) ; daily climatology
;
; Compute averages for each sequebtial day of the year.
; This uses dimension swapping.
;
   do ndy=0,ndys-2                                 ; ndy=0->364 ==> day of year 1->365
      indx = ind( days.eq.(ndy+1) )
      stdDay(ndy,:,:) = dim_stddev(x(lev|:,lat|:,time|indx))
      delete(indx)
   end do
                                                   ; nominal day 366
                                                   ; ave(31 Dec + 1 Jan)=leap
   stdDay(ndys-1,:,:) = (stdDay(0,:,:) + stdDay(ndys-2,:,:))*0.5

   nFill = num(ismissing(stdDay))
   if (nFill.eq.0) then
       delete(stdDay@_FillValue)
   end if

   stdDay@long_name   = "Daily Standard Deviation"
   if (isatt(x,"long_name")) then
       stdDay@long_name = stdDay@long_name +": "+x@long_name
   end if
   if (isatt(x,"units")) then
       stdDay@units     = x@units
   end if
   stdDay@information = "Raw daily averages across all years"
   stdDay@smoothing   = "None"

   year_day           = ispan(1,ndys,1)
   year_day@long_name = "day of year"
   year_day@units     = "ddd"

   stdDay!0           = "year_day"
   stdDay&year_day    =  year_day

   copy_VarCoords(x(0,:,:), stdDay(0,:,:))   ; trick
   delete(stdDay@year_day)                   ; clean up

   return (stdDay)
end

;*******************************10************************************
;*********************************************************************
undef("stdDayTLevLatLon")
function stdDayTLevLatLon (x[*][*][*][*]:numeric, yyyyddd:integer)
;
; calculate the mean Annual Cycle from daily data.
; The return array will gave the raw climatology at each grid point
;
;              x(time,plev,lat,lon)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "plev"
;              x!2 = "lat"
;              x!2 = "lon"
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
;    Leap
;              1908001 =>  Jan  1, 1908]
;              1908032 =>  Feb  1, 1908]
;              1908059 =>  Feb 28, 1908]
;              1908060 =>  Feb 29, 1908]
;              1908061 =>  Mar  1, 1908]
;              1908366 =>  Dec 31, 1908]
;
; Usage:  xStdDay = stdDAY_TLL (x, yyyyddd)
; -------

local dimx, ntim, nlat, mlon,nlev, ndys, days, clmDay, ndy, indx, year_day, nFill
begin

   dimx  = dimsizes (x)

   ntim  = dimx(0)
   nlev  = dimx(1)
   nlat  = dimx(2)
   mlon  = dimx(3)
   ndys  = 366                                     ; allow for leap year

   days  = yyyyddd - (yyyyddd/1000)*1000           ; strip year info [yyyy]

   stdDay= new((/ndys,nlev,nlat,mlon/),typeof(x), getFillValue(x) ) ; daily climatology
;
; Compute averages for each sequebtial day of the year. 
; This uses dimension swapping.
;
   do ndy=0,ndys-2                                 ; ndy=0->364 ==> day of year 1->365       
      indx = ind( days.eq.(ndy+1) )       
      stdDay(ndy,:,:,:) = dim_stddev(x(lev|:,lat|:,lon|:,time|indx))
      delete(indx)                  
   end do
                                                   ; nominal day 366 
                                                   ; ave(31 Dec + 1 Jan)=leap
   stdDay(ndys-1,:,:,:) = (stdDay(0,:,:,:) + clmDay(ndys-2,:,:,:))*0.5

   nFill = num(ismissing(clmDay))
   if (nFill.eq.0) then
       delete(stdDay@_FillValue)
   end if

   stdDay@long_name   = "Daily Standard Deviation"
   if (isatt(x,"long_name")) then
       stdDay@long_name = stdDay@long_name +": "+x@long_name
   end if
   if (isatt(x,"units")) then
       stdDay@units     = x@units
   end if
   stdDay@information = "Raw daily standard deviation across all years"
   stdDay@smoothing   = "None"

   year_day           = ispan(1,ndys,1)
   year_day@long_name = "day of year"
   year_day@units     = "ddd"

   clmDay!0           = "year_day"
   clmDay&year_day    =  year_day
   
   copy_VarCoords(x(0,:,:,:), stdDay(0,:,:,:))   ; trick
   delete(stdDay@year_day)                   ; clean up
 
   return (stdDay)
end

;*******************************11************************************
;*********************************************************************
;/////////////////////////////////////////////////////////////////////////////////////
;author  I. Cionni Niwa  i.cionnii@niwa.co.nz
;////////////////////////////////////////////////////////////////////////////////////
undef("area_cell")
function  area_cell(lata[*]:numeric,lona[*]:numeric,R:numeric)

local lata,lona,lat1,lat2,lat3,lat4,lon1,lon2,lon3,lon4,Area_h,nlat,mlon
begin

   nlat = dimsizes(lata)
   mlon = dimsizes(lona)

   lat1 = new(nlat,typeof(lata))
   lat2 = new(nlat,typeof(lata))
   lat3 = new(nlat,typeof(lata))
   lat4 = new(nlat,typeof(lata))

   lon1 = new(mlon,typeof(lona))
   lon2 = new(mlon,typeof(lona))
   lon3 = new(mlon,typeof(lona))
   lon4 = new(mlon,typeof(lona))

   lat1(0) = -90.
   lat2(0) = -90.
   lat3(0) = -90.+abs((lata(1)-(-90.))/2)
   lat4(0) = -90.+abs((lata(1)-(-90.))/2)

   lat1(nlat-1) = 90.-abs((90.-lata(nlat-2))/2)
   lat2(nlat-1) = 90.-abs((90.-lata(nlat-2))/2)
   lat3(nlat-1) = 90.
   lat4(nlat-1) = 90.

   lata(0) = -90
   lata(nlat-1) = 90

   lat1(1:nlat-2) = lata(1:nlat-2)-abs((lata(1:nlat-2)-lata(0:nlat-3))/2)
   lat2(1:nlat-2) = lata(1:nlat-2)-abs((lata(1:nlat-2)-lata(0:nlat-3))/2)
   lat3(1:nlat-2) = lata(1:nlat-2)+abs((lata(2:nlat-1)-lata(1:nlat-2))/2)
   lat4(1:nlat-2) = lata(1:nlat-2)+abs((lata(2:nlat-1)-lata(1:nlat-2))/2)

   lon1(0) = 0.
   lon2(0) = 0+abs((lona(1)-0)/2)
   lon3(0) = 0+abs((lona(1)-0)/2)
   lon4(0) = 0.

   lon1(mlon-1) = 360.-abs((360.-lona(mlon-2))/2)
   lon2(mlon-1) = 360.
   lon3(mlon-1) = 360.
   lon4(mlon-1) = 360.-abs((360.-lona(mlon-2))/2)

   lona(0) = 0.
   lona(mlon-1) = 360.

   lon1(1:mlon-2) = lona(1:mlon-2)-abs((lona(1:mlon-2)-lona(0:mlon-3))/2)
   lon2(1:mlon-2) = lona(1:mlon-2)+abs((lona(2:mlon-1)-lona(1:mlon-2))/2)
   lon3(1:mlon-2) = lona(1:mlon-2)+abs((lona(2:mlon-1)-lona(1:mlon-2))/2)
   lon4(1:mlon-2) = lona(1:mlon-2)-abs((lona(1:mlon-2)-lona(0:mlon-3))/2)

   Area_h = new((/nlat,mlon/),typeof(lata))

   do i=0,nlat-1
      do k=0,mlon-1
         Area_h(i,k) = gc_qarea((/lat1(i),lat2(i),lat3(i),lat4(i)/),(/lon1(k),lon2(k),lon3(k),lon4(k)/))*R^2
      end do
   end do

   Area_h!0 = "lat"
   Area_h!1 = "lon"
   Area_h&lat = lata
   Area_h&lon = lona

   return(Area_h)
end

;*******************************12************************************
;*********************************************************************
;Function to calculate the weight area using the function area_cell
undef("weight_area_latlon_old")
function weight_area_latlon_old (var[*][*]:numeric,lat:numeric,lon:numeric,lat_min:numeric,lat_max:numeric)
local lat,lon,var,lat_min,lat_max,area_cella,R,Area_Tot,lat1var1,var_wgt_area

begin

   R = 6378   ;this is the value of the earth ray.
   Area_Tot = 4*acos(-1.0)*R^2

   if (lat(0).lt.lat(1)) then
      lat1 = lat
      var1 = var
   else
      lat1 = lat(::-1)
      var1 = var(::-1,:)
   end if
   var1!0 = "lat"
   var1&lat = lat1

   lon1 = lon
   lat_max = lat_max
   lat_min = lat_min

   ;if the coordinates variables are defined as double the computing required can give a segmentation fault
   ;then in the following part variables are defined as float

   if(typeof(lat1).eq."double") then
      lat1_flt = doubletofloat(lat1)
      delete(lat1)
      lat1 = lat1_flt
   end if

   if(typeof(lon1).eq."double") then
      lon_flt = doubletofloat(lon1)
      delete(lon1)
      lon1 = lon_flt
   end if

   area_cella = area_cell(lat1,lon1,R)

   ;somma = sum(area_cella)

   area_cella = area_cella/Area_Tot   ; area of cella are normalized using total area

   var_wgt_area = wgt_areaave2(var1({lat_min:lat_max},:), area_cella({lat_min:lat_max},:), 0)

   ;the option =0: the area sum is calculated using avaible non-missing data
   return(var_wgt_area)
end
;*******************************12vers2************************************
;*********************************************************************
;Function to calculate the weight area using the function area_cell
undef("weight_area_latlon")
function weight_area_latlon (var[*][*]:numeric,lat:numeric,lon:numeric,lat_min:numeric,lat_max:numeric)
local lat,lon,var,lat_min,lat_max,rad,con,dx,dy,dxdy,wgt

begin
   if (typeof(lat).ne."double") then
        R=6.37122e06
   else
    R = 6.37122e06 ;this is the value of the earth ray.
   end if
     
   rad = 4*atan(1.0)/180.

   con=R*rad

   if (lat(0).lt.lat(1)) then
      lat1 = lat
      var1 = var
   else
      lat1 = lat(::-1)
      var1 = var(::-1,:)
   end if
   var1!0 = "lat"
   var1&lat = lat1
 clat=cos(lat1*rad)

   dlon=abs(lon(2)-lon(1))
   dlat=abs(lat(2)-lat(1))
   dx=dlon*con*clat
   ;print(dlon)
   ;print(dlat)
   ;print(dx)

   dy=con*dlat

   dxdy=dy*dx

    wgt=new((/dimsizes(lat1),dimsizes(lon)/),typeof(dxdy))

   wgt=conform(wgt,dxdy,0)
   wgt!0="lat"
   wgt&lat=lat1
   wgt!1="lon"
   wgt&lon=lon

 ;  lon1 = lon
   lat_max = lat_max
   lat_min = lat_min

   ;if the coordinates variables are defined as double the computing required can give a segmentation fault
   ;then in the following part variables are defined as float

     ; area of cella are normalized using total area
   ;print(var1({lat_min:lat_max},:))
   ;print(wgt({lat_min:lat_max},:))

   var_wgt_area = wgt_areaave2(var1({lat_min:lat_max},:), wgt({lat_min:lat_max},:), 0)
   ;print(var_wgt_area)

   ;the option =0: the area sum is calculated using avaible non-missing data
   return(var_wgt_area)
end

;*******************************13************************************
;*********************************************************************
undef("weight_area_lat_old")
function weight_area_lat_old(var[*]:numeric,lat:numeric,lat_min:numeric,lat_max:numeric)

;this function calculate the latitude average using weight calculated with the function area

local R,Area_Tot,lat1,var1,lat,lon,var,lat_min,lat_max,area_cella,area_lat,var_wgt_area,Area_Tot,var_wgt_area

begin
   R = 6378
   Area_Tot = 4*acos(-1.0)*R^2

   if (lat(0).lt.lat(1)) then
      lat1 = lat
      var1 = var
   else
      lat1 = lat(::-1)
      var1 = var(::-1)
   end if
   var1!0 = "lat"
   var1&lat = lat1

   ;This function calculates the weight only for latitude coordinates.
   lon = (/0.,90.,270.,360./)
   area_cella = area_cell(lat1,lon,R)

   somma = sum(area_cella)
   ;somma_pesi = sum(area_cella)
   ;print(somma)
   ;print(Area_Tot)

   area_lat = dim_sum_Wrap(area_cella)
   ;print (area_lat)
   ;option =1 means compute the weighted mean of all non-missing values
   var_wgt_area = dim_avg_wgt_Wrap(var1({lat_min:lat_max}), area_lat({lat_min:lat_max}), 1)
   ;var_wgt_area = dim_avg_Wrap(var1({lat_min:lat_max}))
   return(var_wgt_area)

end

;*******************************14************************************
;*********************************************************************
undef("weight_area_lat")
function weight_area_lat(var[*]:numeric,lat:numeric,lat_min:numeric,lat_max:numeric)

;this function calculate the latitude average using weight calculated with the function area

local R,Area_Tot,lat1,var1,lat,lon,var,lat_min,lat_max,area_cella,area_lat,var_wgt_area,Area_Tot,var_wgt_area

begin
   R = 6.371e06
   rad = 4*atan(-1.0)/180

   con=R*rad

   if (lat(0).lt.lat(1)) then
      lat1 = lat
      var1 = var
   else
      lat1 = lat(::-1)
      var1 = var(::-1)
   end if
   var1!0 = "lat"
   var1&lat = lat1

    clat=cos(lat1*rad)
   ;This function calculates the weight only for latitude coordinates.

   dlon=360.
   dlat=abs(lat1(2)-lat1(1))
   dx=dlon*con*clat
   dy=con*dlat
   dxdy=dy*dx

   wgt=new((/dimsizes(lat1)/),typeof(dxdy))
   wgt=dxdy
   wgt!0="lat"
   wgt&lat=lat1

   var_wgt_area = dim_avg_wgt_Wrap(var1({lat_min:lat_max}), wgt({lat_min:lat_max}), 1)
   ;var_wgt_area = dim_avg_Wrap(var1({lat_min:lat_max}))
   return(var_wgt_area)

end

;*******************************14************************************
;*********************************************************************
undef ("month_selly_mean")
function month_selly_mean(x:numeric, time:numeric)
local x, time, yy,mm,dd,date1,date2,date3,n_yy,n_mm,n_tim,x_monavg,indice, date_p,nt,n_mo,n_yr,yy_d,mm_d,dd_d
begin
   date1 = cd_calendar(time,2)
   x!0 = "date1"
   x&date1 = date1

   yy = doubletointeger(date1/10000)
   mm = doubletointeger((date1-(yy*10000))/100)
   dd = doubletointeger((date1-(yy*10000))-(mm*100))

   date2 = yy*100+mm
   x!0 = "date2"
   x&date2 = date2

   n_yy = max(yy)-min(yy)+1
   n_mm = 12

   n_tim = n_yy*n_mm

   x_monavg = new((/n_tim/),typeof(x))
   date3 = new((/n_tim/),typeof(date1))

   nt = -1
   do n_yr=0,n_yy-1
      do n_mo=0,n_mm-1
         date_p = (min(yy)+n_yr)*100+(n_mo+1)

         indice = ind(date2.eq.date_p)

         nt = nt+1

         x_monavg(nt) = (/dim_avg(x(indice))/)
         date3(nt) = date_p
         delete(date_p)
         delete(indice)
      end do
   end do

   x_monavg!0 = "time"
   bo = date3/100.

   yy_d = doubletointeger(date3/100.)
   mm_d = doubletointeger(date3-yy_d*100.)
   dd_d = yy_d
   mim_d = yy_d
   hh_d = yy_d
   min_d = yy_d
   sec_d = yy_d

   dd_d(:) = 15
   hh_d(:) = 0
   min_d(:) = 0
   sec_d(:) = 0

   time1 = cd_inv_calendar(yy_d,mm_d,dd_d,hh_d,min_d,sec_d,time@units,0)
   x_monavg&time = time1

   return (x_monavg)

end

;*******************************15************************************
;*********************************************************************
;_______________________________________________________________________________
;      function  calcRegCoeffs(y:numeric); y(time,lat)
;_______________________________________________________________________________
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
undef("calcRegCoeffs")
function calcRegCoeffs(y:numeric)
local y,d_y,d_t,ma,x,funcx,pi,t,a
begin
   d_y = dimsizes(y)
   d_t = d_y(0)

   ma = 14
   x = fspan(1,d_t,d_t)
   a = new((/ma/),typeof(y))

   funcx = new((/ma,d_t/),"float")
   pi = 4.0 * atan(1.0)
   t = 2*pi*(x-0.5)/12

;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   funcx(0,:) = 1
   funcx(1,:) = sin(t)
   funcx(2,:) = cos(t)
   funcx(3,:) = sin(2*t)
   funcx(4,:) = cos(2*t)
   funcx(5,:) = sin(3*t)
   funcx(6,:) = cos(3*t)
   funcx(7,:) = sin(4*t)
   funcx(8,:) = cos(4*t)
   funcx(9,:) = x
   funcx(10,:) = sin(t)*x
   funcx(11,:) = cos(t)*x
   funcx(12,:) = sin(2*t)*x
   funcx(13,:) = cos(2*t)*x

;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;::::::::::::::::LEAST SQUARE FIT::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

   a(:) = reg_multlin(y(:),funcx,False)
   a@_FillValue =1e20
   delete_VarAtts(a, "constant")

   a@long_name = "Coefficients A and B"
   a!0 = "number of coefficients"

   return(a)
end

;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


;; -------------------------------------------------------
; CCMOzoneModelData.ncl
; delphi unit gb_CCMColumnOzoneModelData
;; -------------------------------------------------------
;*******************************16************************************
;*********************************************************************
undef("genZonalMeans")
function genZonalMeans (x:numeric, lat:numeric, startLat:numeric, endLat:numeric)
local dimx, ntim, nlat, nzone, z, t, i, j, cnt, idx_st, idx_ed, xZonMean
begin
   dimx = dimsizes(x)
   ntim = dimx(0)
   nlat = dimx(1)
   nlon = dimx(2)
   nzone = dimsizes(startLat)

   nmos = 12
   modCheck ("genZonalMeans", ntim, nmos)   ; error check

   xZonMean = new((/ntim,nzone/),typeof(x))

   pi  = 4.*atan(1.0)
   rad = (pi/180.)

   ; loop through lat_band
   do t=0,(ntim-1)
      do z=0,(nzone-1)
         idx_st = 0
         do while (lat(idx_st).lt.startLat(z))
            idx_st = idx_st + 1
         end do
         idx_ed = nlat-1;
         do while (lat(idx_ed).gt.endLat(z))
            idx_ed = idx_ed - 1
         end do

         zoneTotalX = 0.0
         ZoneTotalArea = 0.0

         do i=idx_st,idx_ed
            if (i.eq.idx_st) then
               Lat1 = startLat(z)
            else
               Lat1 = (lat(i)+lat(i-1))/2
            end if
            if (i.eq.idx_ed) then
               Lat2 = endLat(z)
            else
               Lat2 = (lat(i)+lat(i+1))/2
            end if
            lat_slice = x(t,i,:)
            idx_lon = ind(.not.ismissing(lat_slice))
            if (.not.all(ismissing(idx_lon))) then
               CellArea = abs(sin(Lat1*rad)-sin(Lat2*rad))
               zoneTotalX = zoneTotalX + sum(lat_slice) *  CellArea
               ZoneTotalArea = ZoneTotalArea + dimsizes(idx_lon) * CellArea
            end if
            delete(idx_lon)
         end do  ; loop(i)

         if (ZoneTotalArea.ne.0) then
            xZonMean(t,z) = zoneTotalX / ZoneTotalArea
         end if
      end do  ; loop(z)
   end do  ; loop(t)

   ; Create an informational attribute:
   xZonMean@info = "zonal mean: CCMOzoneModelData.ncl"

   return (xZonMean)
end

;*******************************17************************************
;*********************************************************************
undef("calcMeanAnnCycleMonthly")
function calcMeanAnnCycleMonthly (RegCoeffAs:numeric)
local dimx, z, month, t, xMeanAnnCycleMon
begin
   dimx = dimsizes(RegCoeffAs)

   pi  = 4.*atan(1.0)

   month = fspan(1,12,12)
   t = 2*pi*(month-0.5)/12.0

   xMeanAnnCycleMon = new((/12/),"float")

   xMeanAnnCycleMon(:) = RegCoeffAs(0)+ \
      RegCoeffAs(1)*sin(t)+RegCoeffAs(2)*cos(t)+ \
      RegCoeffAs(3)*sin(2*t)+RegCoeffAs(4)*cos(2*t)+ \
      RegCoeffAs(5)*sin(3*t)+RegCoeffAs(6)*cos(3*t)+ \
      RegCoeffAs(7)*sin(4*t)+RegCoeffAs(8)*cos(4*t)

   ; Create an informational attribute:
   xMeanAnnCycleMon@info = "mean annual cycle monthly: CCMOzoneModelData.ncl"

   return (xMeanAnnCycleMon)
end

;*******************************18************************************
;*********************************************************************
undef("calcMeanAnnCycleAnnual")
function calcMeanAnnCycleAnnual(RegCoeffAs:numeric)
local dimx, day, t, xMeanAnnCycleDay
begin
   dimx = dimsizes(RegCoeffAs)

   pi  = 4.*atan(1.0)

   day = fspan(1,366,366)
   t = 2*pi*(day-0.5)/366.0

   xMeanAnnCycleDay = new((/366/),"float")

   xMeanAnnCycleDay(:) = RegCoeffAs(0)+ \
       RegCoeffAs(1)*sin(t)+RegCoeffAs(2)*cos(t)+ \
       RegCoeffAs(3)*sin(2*t)+RegCoeffAs(4)*cos(2*t)+ \
       RegCoeffAs(5)*sin(3*t)+RegCoeffAs(6)*cos(3*t)+ \
       RegCoeffAs(7)*sin(4*t)+RegCoeffAs(8)*cos(4*t)

   ; Create an informational attribute:
   xMeanAnnCycleDay@info = "mean annual cycle daily: CCMOzoneModelData.ncl"

   return (xMeanAnnCycleDay)
end

;*******************************19************************************
;*********************************************************************
undef("rmMeanAnnCycle")
procedure rmMeanAnnCycle (x:numeric, xMeanAnnCycle:numeric, RegCoeffA1:numeric, xMonAnom:numeric, xAnnAnom:numeric)
local dimx, ntim,t, i, idx_st, idx_ed, nyr, xdata
begin
   dimx = dimsizes(x)
   ntim = dimx(0)

   nmos  = 12
   modCheck ("rmMeanAnnCycle", ntim, nmos)   ; error check

   nyr = ntim/nmos

   ; loop through lat_band
   do t=0,(ntim-1)
      mon = mod(t,12)
      xMonAnom(t) = (x(t)-xMeanAnnCycle(mon))*100.0/RegCoeffA1
   end do  ; loop(t)

   do t=0,(nyr-1)
      idx_st = t * 12
      idx_ed = (t + 1) * 12 - 1
      xdata = xMonAnom(idx_st:idx_ed)
      if(.not.any(ismissing(xdata))) then
         xAnnAnom(t) = avg(xdata)
      end if
   end do  ; loop(t)

   ; Create an informational attribute:
   xMonAnom@info = "monthly anomalies: CCMOzoneModelData.ncl"
   xAnnAnom@info = "annual anomalies: CCMOzoneModelData.ncl"

end

;*******************************20************************************
;*********************************************************************
undef("apfiltersmooth")
function apfiltersmooth(xMonAnom:numeric, filter:numeric, iterations:integer)
local dimx, dimf, ntim, nwin, z, t, i, j, hfw, idx_ed, avgwidth,\
      totalwgt, total, num, minval,finished, data, hold, xSmthMonAnom
begin
   dimx = dimsizes(xMonAnom)
   ntim = dimx(0)

   nmos = 12
   modCheck ("apfiltersmooth", ntim, nmos)   ; error check

   dimf = dimsizes(filter)
   nwin = dimf(0)
   hfw = nwin/2

   data = new(nwin,typeof(xMonAnom))
   hold = new(ntim,typeof(xMonAnom))
   xSmthMonAnom = new((/ntim/),typeof(xMonAnom))

   do i=1,iterations
      ;{Transfer the data to the storage array.}
      if (i.eq.1) then
         hold(:) = xMonAnom(:,z)
      else
         hold(:) = xSmthMonAnom(:,z)
      end if

      ;{Fill the data array for the first time.}
      do j=1,hfw
         data(j-1) = -999.99
      end do
      do j=hfw,(nwin-1)
         data(j) = hold(j-hfw)
      end do
      idx_ed = nwin-hfw

      ;{Process the data array.}
      do t=0,(ntim-1)
         if (data(hfw).lt.-999.0) then
            xSmthMonAnom(t,z) = -999.99
         else
            avgwidth = 0
            totalwgt = filter(hfw)
            total = filter(hfw) * data(hfw)
            num = 1
            finished = False
            do while (.not.finished)
               avgwidth = avgwidth+1
               if ((data(hfw-avgwidth).gt.-999.0).and.(data(hfw+avgwidth).gt.-999.0)) then
                  totalwgt = totalwgt + filter(hfw-avgwidth) + filter(hfw+avgwidth)
                  total = total+(filter(hfw-avgwidth)*data(hfw-avgwidth))+(filter(hfw+avgwidth)*data(hfw+avgwidth))
                  num = num+2
               else
                  finished = True
               end if
               if (avgwidth.eq.hfw) then
                  finished = True
               end if
            end do
            if (i.eq.iterations) then
               minval = hfw
            else
               minval = 0
            end if
            if (num.gt.minval) then
               xSmthMonAnom(t) = total/totalwgt
            else
               xSmthMonAnom(t) = -999.99
            end if
         end if

         do j=1,(nwin-1)
            data(j-1) = data(j)
         end do
         idx_ed = idx_ed+1
         if (idx_ed.gt.ntim) then
            data(nwin-1) = -999.99
         else
            data(nwin-1) = hold(idx_ed-1)
         end if
      end do  ; loop(t)
    end do  ; loop(i)

   ; Create an informational attribute:
   xSmthMonAnom@info = "smoothed monthly anomalies: CCMOzoneModelData.ncl"
   xSmthMonAnom@_FillValue = -999.99

   return (xSmthMonAnom)

end

;*******************************21************************************
;*********************************************************************
undef("smoothAnomalies")
procedure smoothAnomalies(xAnnAnom:numeric, iterations:integer)
local dimx, dimf, ntim, nset, nwin, s, t, i, j, hfw, idx_ed, avgwidth,\
      totalwgt, total, num, minval,finished, smthdata, hold, buffer,\
      filter, method
begin
   filter = (/1.,2.,1./)
   method = "EdgeP"

   dimx = dimsizes(xAnnAnom)
   ntim = dimx(0)
   ;nset = dimx(1)

   dimf = dimsizes(filter)
   nwin = dimf(0)
   hfw = nwin/2

   buffer = new(nwin,typeof(xAnnAnom))
   smthdata = new(ntim,typeof(xAnnAnom))
   hold = new(ntim,typeof(xAnnAnom))

   ;do s=0,(nset-1)
   do i=1,iterations
      ;{Transfer the data to the storage array.}
      if (i.eq.1) then
         hold = (/xAnnAnom/)
      else
         hold = (/smthdata/)
      end if

      ;{Fill the data array for the first time.}
      do j=1,hfw
         buffer(j-1) = -999.99
      end do
      do j=hfw,(nwin-1)
         buffer(j) =(/hold(j-hfw)/)
      end do
      idx_ed = nwin-hfw

      ;{Process the data array.}
      do t=0,(ntim-1)
         if (buffer(hfw).lt.-999.0) then
            smthdata(t) = -999.99
         else
            avgwidth = 0
            totalwgt = filter(hfw)
            total = filter(hfw) * buffer(hfw)
            num = 1
            finished = False
            do while (.not.finished)
               avgwidth = avgwidth+1
               if ((buffer(hfw-avgwidth).gt.-999.0).and.(buffer(hfw+avgwidth).gt.-999.0)) then
                  totalwgt = totalwgt + filter(hfw-avgwidth) + filter(hfw+avgwidth)
                  total = total+(filter(hfw-avgwidth)*buffer(hfw-avgwidth))+(filter(hfw+avgwidth)*buffer(hfw+avgwidth))
                  num = num+2
               else
                  finished = True
               end if
               if (avgwidth.eq.hfw) then
                  finished = True
               end if
            end do

            if (method.eq."EdgeP") then
               if (num.eq.1) then
                  if (buffer(hfw+1).lt.-999.0) then
                     smthdata(t) = (buffer(hfw) + buffer(hfw-1))/2.0
                  else
                     smthdata(t) = (buffer(hfw) + buffer(hfw+1))/2.0
                  end if
               else
                  smthdata(t) = total/totalwgt
               end if
            else
               if (i.eq.iterations) then
                  minval = hfw
               else
                  minval = 0
               end if
               if (num.gt.minval) then
                  smthdata(t) = total/totalwgt
               else
                  smthdata(t) = -999.99
               end if
            end if ; if(method)
         end if
         do j=1,(nwin-1)
            buffer(j-1) = buffer(j)
         end do
         idx_ed = idx_ed+1
         if (idx_ed.gt.ntim) then
            buffer(nwin-1) = -999.99
         else
            buffer(nwin-1) = hold(idx_ed-1)
         end if
      end do  ; loop(t)
   end do  ; loop(i)
   xAnnAnom(:) = (/smthdata(:)/)
   ;end do  ; loop(s)

   ; Create an informational attribute:
   xAnnAnom@info = "smoothed annual anomalies: CCMOzoneModelData.ncl"
   xAnnAnom@_FillValue = -999.99

end

;*******************************22************************************
;*********************************************************************
undef("Z")
function Z()
local j
begin
   j = ind(mod(t,12).eq.0)
end

;*******************************23************************************
;*********************************************************************
;From lat_band string to inf and sup value
;************************************************************************
undef("read_band")
procedure read_band(lat_band:string,inf:numeric,sup:numeric)

local lat_band,inf,sup,e,numei,numeis

begin
   inf = new((/dimsizes(lat_band)/),"float")
   sup = new((/dimsizes(lat_band)/),"float")

   do i=0,(dimsizes(lat_band)-1)
      e = stringtocharacter(lat_band(i))
      numeii = stringtofloat(charactertostring(e(0:1)))

      if (charactertointeger(e(2)).eq.78)then
         inf(i) = numeii
      else
         inf(i) = -numeii
      end if
      numeis = stringtofloat(charactertostring(e(4:5)))

      if (charactertointeger(e(6)).eq.78) then
         sup(i) = numeis
      else
         sup(i) = -numeis
      end if
   end do

   return (inf)
   return (sup)

end

;*******************************24************************************
;*********************************************************************
;_______________________________________________________________________________

;      function  least_square_fourier(y:numeric,ma:numeric,d_t:numeric)
;_______________________________________________________________________________
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
undef("least_square_fourier")
procedure least_square_fourier(y:numeric,ma:numeric,d_t:numeric,a:numeric,det_y:numeric)
local dt,x,funcx,pi,t,a,det_y
begin
   d_t = dimsizes(y)
   x = fspan(1,d_t,d_t)
                                         ;number coefficients
   funcx = new((/ma,d_t/),"float")
   pi = 3.141592654
   t = 2*pi*(x-0.5)/12
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   funcx(0,:) = 1
   funcx(1,:) = sin(t)
   funcx(2,:) = cos(t)
   funcx(3,:) = sin(2*t)
   funcx(4,:) = cos(2*t)
   funcx(5,:) = sin(3*t)
   funcx(6,:) = cos(3*t)
   funcx(7,:) = sin(4*t)
   funcx(8,:) = cos(4*t)
   funcx(9,:) = x
   funcx(10,:) = sin(t)*x
   funcx(11,:) = cos(t)*x
   funcx(12,:) = sin(2*t)*x
   funcx(13,:) = cos(2*t)*x
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;::::::::::::::::LEAST SQUARE FIT::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   a = reg_multlin (y,funcx,False)
   a@_FillValue = 1e20
   delete_VarAtts(a, "constant")

   det_y = a(0)+a(1)*sin(t)+a(2)*cos(t)+a(3)*sin(2*t)+a(4)*cos(2*t)+a(5)*sin(3*t)+a(6)*cos(3*t)+a(7)*sin(4*t)+a(8)*cos(4*t)

   return(a)
   return(det_y)

end
;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


;*******************************25************************************
;*********************************************************************
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;I.Cionni. NIWA.
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
undef ("monthly_mean")
function monthly_mean(x:numeric, time:numeric)

local x, time, yy,mm,dd,date1,date2,date3,n_yy,n_mm,n_tim,x_monavg,indice, date_p,nt,n_mo,n_yr,yy_d,mm_d,dd_d,sec_d,min_d,time1

begin
   date1 = cd_calendar(time,2)   ; time as yyyymmdd

   yy = doubletointeger(date1/10000)
   mm = doubletointeger((date1-(yy*10000))/100)
   dd = doubletointeger((date1-(yy*10000))-(mm*100))

   date2 = yy*100+mm              ;time as yyyymm
   x!0 = "date2"
   x&date2 = date2
   n_yy = max(yy)-min(yy)+1
   n_mm = 12

   n_tim = n_yy*n_mm

   x_monavg = new((/n_tim/),typeof(x))
   date3 = new((/n_tim/),typeof(date1))

   nt = -1
   do n_yr=0,n_yy-1
      do n_mo=0,n_mm-1
         date_p = (min(yy)+n_yr)*100+(n_mo+1)
         indice = ind(date2.eq.date_p)

         nt = nt+1

         x_monavg(nt) = (/dim_avg(x(date2|indice))/)
         date3(nt) = date_p
         delete(date_p)
         delete(indice)
      end do
   end do

   x_monavg!0 = "time"
   bo = date3/100.

   yy_d = doubletointeger(date3/100.)
   mm_d = doubletointeger(date3-yy_d*100.)
   dd_d = yy_d
   mim_d = yy_d
   hh_d = yy_d
   min_d = yy_d
   sec_d = yy_d

   dd_d(:) = 15
   hh_d(:) = 0
   min_d(:) = 0
   sec_d(:) = 0

   time1 = cd_inv_calendar(yy_d,mm_d,dd_d,hh_d,min_d,sec_d,time@units,0)
   x_monavg&time = time1

   return (x_monavg)

end
;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::