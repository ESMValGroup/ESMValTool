; file: functions_contrib.ncl
; some contributed functions

;************************************************************
; Mark Stevens
; compute area of seaice in km^2 x 10^6 from seaice fraction 
; for each hemisphere of regular lat/lon grid
;************************************************************
undef("iceareaFIX")
function iceareaFIX (icefrac[*][*]:numeric,hem:integer) 
local lat, nlat, lon, nlon, RE, pi, area, d2r, wgt, tsum, \
      nwgt, boxarea, hemarea, j, icearea
begin
; icearea is the fractional sea ice area (0.0-1.0)
; icefrac first dimension is latitude, second is longitude
; hem = 0 (SH), 1 (NH)

  if (.not.isatt(icefrac,"_FillValue")) then
    icefrac@_FillValue = getFillValue(icefrac)
  end if 
  lat = icefrac&lat
  nlat = dimsizes(lat)
  lon = icefrac&lon
  nlon = dimsizes(lon)

  RE = 6.37122e3             ; radius of earth in km
  pi = acos(-1.0)
  area = 4.*pi*RE^2          ; surface area of Earth 
  d2r = pi/180.0             ; convert degrees to radians

  if (typeof(lat) .eq. "double") then
    wgt = doubletofloat(NormCosWgtGlobe(lat))
  else
    wgt = NormCosWgtGlobe(lat)
  end if
  tsum = sum(wgt)            ; sum of all weights
  nwgt = wgt/tsum            ; frac of sphere of each lat band
  boxarea = area*nwgt/nlon   ; area of each grid box (lat) in km^2
  hemarea = new(nlat,float)

  if (hem .eq. 0) then       ; Southern Hemisphere
    do j = 0, nlat/2-1 
      hemarea(j) =  sum(boxarea(j)*icefrac(j,:))
    end do
  else                       ; Northern Hemisphere
    do j = nlat/2, nlat-1 
      hemarea(j) =  sum(boxarea(j)*icefrac(j,:))
    end do
  end if
  icearea = sum(hemarea)/1.e6

  return(icearea)      ; return area of ice km^2 x 10^6 
end

;************************************************************
; Mark Stevens
; compute area of seaice in km^2 x 10^6 from seaice fraction 
; for each hemisphere of gaussian grid
;************************************************************
undef("iceareaGAU")
function iceareaGAU (icefrac[*][*]:numeric,hem:integer)
local lat, nlat, lon, nlon, RE, pi, area, gw, tsum, \
      nwgt, boxarea, hemarea, j, icearea 
begin
; icearea is the fractional sea ice area (0.0-1.0)
; icefrac first dimension is latitude, second is longitude
; hem = 0 (SH), 1 (NH)

  if (.not.isatt(icefrac,"_FillValue")) then
    icefrac@_FillValue = getFillValue(icefrac)
  end if 
  lat = icefrac&lat
  nlat = dimsizes(lat)
  lon = icefrac&lon
  nlon = dimsizes(lon)

  RE = 6.37122e3            ; radius of earth in km
  pi = acos(-1.0)
  area = 4.*pi*RE^2         ; surface area of Earth 

  gw = latGauWgt(nlat,"lat","gaussian weights","none")
  tsum = sum(gw)            ; sum of all weights
  nwgt = gw/tsum            ; frac of sphere of each lat band
  boxarea = area*nwgt/nlon  ; area of each grid box (lat) in km^2
  hemarea = new(nlat,float)

  if (hem .eq. 0) then         ; Southern Hemisphere
    do j = 0, nlat/2-1 
      hemarea(j) =  sum(boxarea(j)*icefrac(j,:))
    end do
  else
    do j = nlat/2, nlat-1 
      hemarea(j) =  sum(boxarea(j)*icefrac(j,:))
    end do
  end if
  icearea = sum(hemarea)/1.e6

  return(icearea)      ; return area of ice km^2 x 10^6 
end

;************************************************************
undef("infoTimeStamp")
procedure infoTimeStamp (wks:graphic, chSize:float, FileInfo:string)

; Place text on plot to indicate the data source and time created
; This procedure must be invoked prior to any advance
; frame call.
; examples:
;          (0) wks =  gsn_open_wks("ncgm","bogus") 
;          (1) infoTimeStamp (wks, 0.01, "Bogus File")
;              [plot]
;
;          (2) mssPath = "/SHEA/ECMWF/"
;              mssName = "sample"
;              size    = 0.015
;              infoTimeStamp (wks, size, mssPath+mssName )
;              [plot]

local yBot, xLeft, xRight, txres, TimeStamp
begin
  yBot                = 0.001
  xLeft               = 0.001
  xRight              = 0.999
                                        ; Sylvia Murphy mods
  if (NhlClassName(wks).eq."psWorkstationClass") then
    getvalues wks
      "wkDeviceLowerX"   : ps_xLeft
      "wkDeviceLowerY"   : ps_yLeft
      "wkDeviceUpperX"   : ps_xRight
    end getvalues
      
    if (ps_xLeft .lt. 0) then         ; 36 is default, 0 is the beginning
      xoff_set = fabs(ps_xLeft/612.)  ; 612 is 8.5 * 72 in-1
      xLeft = 0.02941 + xoff_set      ; .02941 is 1/4 of an inch
    else                              ; which is the margin required
      xLeft = 0.02941                 ; when printing postscript 
    end if

    if (ps_xRight .gt. 612) then      ; 576 is default,612 is end
      xoff_set = fabs(1 - (ps_xRight/612.)) ; 612 is 8.5 * 72 in-1
      xRight= 0.97059 - xoff_set
    else
      xRight = 0.97059
    end if

    if (ps_yLeft .lt. 0) then         ; 126 is default, 0 is the beginning
      yoff_set = fabs(ps_yLeft/792.)  ; 792 is 11  * 72 in-1
      yBot  = 0.02941  + yoff_set
    else
      yBot  = 0.02941
    end if
  end if

  txres               = True          ; additional info
  txres@txFontHeightF = chSize        ; size of meta data
  txres@txJust        = "BottomRight"
  gsn_text_ndc (wks, FileInfo  , xRight , yBot , txres)

  TimeStamp           = systemfunc( "date" )
  txres@txJust        = "BottomLeft"
  gsn_text_ndc (wks,"Created: "+TimeStamp , xLeft, yBot, txres)
end

;*****************************************************************
undef("FixZeroContour")
function FixZeroContour (CNLVLS[*]:float, label:string)
; called internally
local eps, indEps
begin
  cnlvls = CNLVLS                          ; historical
  if (dimsizes(cnlvls).gt.1) then
    eps    = 1.e-09                          ; arbitrary
    indEps=ind(fabs(cnlvls).le.eps) 
    if (.not.ismissing(indEps)) then
      cnlvls(indEps) = 0.0            ; the "zero" line ==>-0.8e-09
    ;else                            ; debug print
    ;    print (label+": no zero contours")
    end if
  end if
  return (cnlvls)                                      
end

;*******************************************************************
undef("get_cnLevels")
function get_cnLevels (plot:graphic) 
local cnlvls
begin
  if (isatt(plot,"contour")) then
    getvalues plot@contour       
      "cnLevels"   : cnlvls
    end getvalues
  else
    getvalues plot
      "cnLevels"   : cnlvls
    end getvalues
  end if
  return(cnlvls)
end

;*********************************************************************
undef("ZeroNegDashLineContour")
function ZeroNegDashLineContour (plot:graphic) 

; operates on a plot object created by "gsn_csm.ncl" code
; Make zero line twice as thick and set neg contourlines to dash
; Dash line patterns: http://ngwww.ucar.edu/ngdoc/ng/ref/dashpatterns.html

local cnlvls, cnlinepat, cnlinethk, n, N
begin
  cnlvls = get_cnLevels (plot)

  N = dimsizes(cnlvls)
  if (ismissing(N) .or. N.le.0) then
    print("ZeroNegDashLineContour: dimsizes(cnlvls)="+N+" return (non-fatal)")
    return (plot)
  else 
    cnlvls = FixZeroContour (cnlvls, "ZeroNegDashLineContour")
  end if

  if (any(cnlvls.le.0.)) then
    cnlinepat = new (dimsizes(cnlvls), integer) ; line pattern vector
    cnlinepat = 0                               ; default is solid (=0)
    cnlinethk = new (dimsizes(cnlvls), integer) ; line thick vector
    cnlinethk = 1                               ; default

    do n=0,N-1
      if (cnlvls(n).lt.0.) then
        cnlinepat(n) = 5               ; simple dash line pattern
      end if
      if (cnlvls(n).eq.0.) then
        cnlinethk(n) = 2               ; make the zero contour thicker
      end if
    end do
    
    if (isatt(plot,"contour")) then
      setvalues plot@contour
        "cnMonoLineDashPattern"  : False
        "cnLineDashPatterns"     : cnlinepat
        "cnMonoLineThickness"    : False
        "cnLineThicknesses"      : cnlinethk
      end setvalues
    else
      setvalues plot
       "cnMonoLineDashPattern"  : False
       "cnLineDashPatterns"     : cnlinepat
       "cnMonoLineThickness"    : False
       "cnLineThicknesses"      : cnlinethk
      end setvalues
    end if
  end if                                ; any
  return (plot)
end

;************************************************************
undef("SigPattern")
function SigPattern (plot:graphic, CnLt:float, FillPatLt:integer ) 

; variation of Dennis Shea's ShadeLtContour
; operates on a plot object created by "gsn_csm.ncl" code
; Fill all contours less than or equal to "CnLt" to the
;      pattern corresponding to "FillPatLt" 
;      http://ngwww.ucar.edu/ngdoc/ng/ref/fillpatterns.html

local cnlvls, patterns, i, N
begin
                                        ; Retrieve contour levels.
   cnlvls = get_cnLevels (plot)

   N = dimsizes(cnlvls)
   if (ismissing(N) .or. N.le.0) then
     print ("SigPattern: dimsizes(cnlvls)=" \
            +N+" return (non-fatal)")
     return (plot)
   end if

   if (any(cnlvls.le.CnLt)) then
     patterns  = new(dimsizes(cnlvls)+1,integer) ; Create array for fill 
     patterns(:)  = -1                           ; patterns and initialize
                                                   ; it to transparent.
     do i=0,N-1                        ; Fill contour levels depending on 
       if(cnlvls(i).le.CnLt) then      ; different criteria.
         patterns(i) = FillPatLt       ; see above URL
       end if
     end do
     
     if (isatt(plot,"contour")) then
       setvalues plot@contour           ; Retrieve contour levels.
         "cnFillOn"          : True
         "cnMonoFillColor"   : True
         "cnMonoFillPattern" : False
         "cnFillPatterns"    : patterns
       end setvalues
     else
       setvalues plot
         "cnFillOn"          : True
         "cnMonoFillColor"   : True
         "cnMonoFillPattern" : False
         "cnFillPatterns"    : patterns
       end setvalues
     end if
   end if           
 
   return (plot)
end
