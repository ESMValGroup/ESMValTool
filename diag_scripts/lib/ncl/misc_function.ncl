
;##################################################################################
;# This NCL script contains the new functions for the CCMVal Diagnostic Tool      #
;##################################################################################
;# Author:                                                                        #
;# function month_sel(month_names:string)
;# function lat_names(lat_val:numeric) 
;# procedure add_line(wks,plot,x,yval,line_color,line_dash)    
;# procedure add_scatt(wks,plot,x,yval,mark_color,mark_ind)
;# procedure add_legenda(wks,MODEL:string,colors:string,dashs:numeric,scatters:string,ticknesses:numeric,place:string)   
;# Irene Cionni (i.cionni@niwa.co.nz)                                             #
;##################################################################################
;# The NCL script is called by the NCL shell scripts E06FIGXX.ncl                 #
;##################################################################################

;Function to calculate the weight area using the functionarea_cell
undef("month_sel")
function month_sel(month_names:string)
local month_names,month_n
begin
   if (.not.any(ismissing(ind(month_names.eq."ANN")))) then
      month_n = new(12,"integer")
      month_n = (/1,2,3,4,5,6,7,8,9,10,11,12/)
   else
      month_n = new(dimsizes(month_names),"integer")
      do i=0,dimsizes(month_names)-1
         if (.not.ismissing(ind(month_names(i).eq."JAN"))) then
            month_n(i)=1
         end if
         if (.not.ismissing(ind(month_names(i).eq."FEB"))) then
            month_n(i)=2
         end if
         if (.not.ismissing(ind(month_names(i).eq."MAR"))) then
            month_n(i)=3
         end if
         if (.not.ismissing(ind(month_names(i).eq."APR"))) then
            month_n(i)=4
         end if
         if (.not.ismissing(ind(month_names(i).eq."MAY"))) then
            month_n(i)=5
         end if
         if (.not.ismissing(ind(month_names(i).eq."JUN"))) then
            month_n(i)=6
         end if
         if (.not.ismissing(ind(month_names(i).eq."JUL"))) then
            month_n(i)=7
         end if
         if (.not.ismissing(ind(month_names(i).eq."AUG"))) then
            month_n(i)=8
         end if
         if (.not.ismissing(ind(month_names(i).eq."SEP"))) then
            month_n(i)=9
         end if
         if (.not.ismissing(ind(month_names(i).eq."OCT"))) then
            month_n(i)=10
         end if
         if (.not.ismissing(ind(month_names(i).eq."NOV"))) then
            month_n(i)=11
         end if
         if (.not.ismissing(ind(month_names(i).eq."DEC"))) then
            month_n(i)=12
         end if
      end do
   end if

   return(month_n)

end

;*************************************************
undef("lat_names")
function lat_names(lat_val:numeric)
local lat_name,lat_val
begin
   lat_name = new((/dimsizes(lat_val)/),"string")
   do i=0,dimsizes(lat_val)-1
      if (lat_val(i).gt.0.) then
         lat_name(i) = lat_val(i)+"N"
      end if
      if (lat_val(i).eq.0.) then
         lat_name(i) = "Equator"
      end if
      if (lat_val(i).lt.0.) then
         lat_name(i) = -lat_val(i)+"S"
      end if
   end do

   return (lat_name)

end


;*******ADDLINE PROCEDURE***************************************
undef("add_line")
procedure add_line(wks,plot,x,yval,line_color,line_dash)
local plres,str,y,plot,x,yval,line_color,line_dash
begin
; A-laue_ax+
   ; gsn_add_polyline crashes if the arrays x and/or yval contain
   ; only missing values
   if (all(ismissing(x)).or.all(ismissing(yval))) then
       return
   end if
; A-laue_ax-
   plres = True
   plres@gsLineColor = line_color
   plres@gsLineDashPattern = line_dash
   plres@gsLineThicknessF=2.5
   ;y=(/yval,yval/)
   str = unique_string("polyline")
   plot@$str$ = gsn_add_polyline(wks,plot,x,yval,plres)
end

;*******ADDSCATT PROCEDURE***************************************
undef("add_scatt")
procedure add_scatt(wks,plot,x,yval,mark_color,mark_ind)
local plres,str,y,plot,x,yval,line_color,line_dash
begin
   plres = True
   plres@gsMarkerColor = mark_color
   plres@gsMarkerIndex = mark_ind
   ;y=(/yval,yval/)
   str = unique_string("polymarker")
   plot@$str$ = gsn_add_polymarker(wks,plot,x,yval,plres)
end

;********LEGEND PROCEDURE*************************************
undef("add_legenda")
procedure add_legenda(wks,MODEL:string,colors:string,dashs:numeric,scatters:string,ticknesses:numeric,place:string)
local MODEL,dim_mod,colors,dashs,edge_x,edge_y,edge_res,tx1res,dim_lb,xleg,xl,yleg,yl,xtxt,ytxt,place,dim_lines,space_lines,lin_sp,scatters,ticknesses
begin
   ;***************Create Legend***********************************
   dim_mod = dimsizes(MODEL)
   dim_lines = floattoint(dim_mod/3)+1
   if ((mod(dim_mod,3)).eq.0 ) then
      dim_lines = floattoint(dim_mod/3)+1
   else
      dim_lines = floattoint(dim_mod/3)+2
   end if
   lin_sp = 0.025
   space_lines = (dim_lines+1)*lin_sp

   if (place.eq."bottom")
      top_line = 0.20;8
   end if
   if (place.eq."middle")
      top_line = 0.42
   end if
   if (place.eq."top")
      top_line = 0.99
   end if
   if ((place.ne."bottom").and.(place.ne."middle").and.(place.ne."top")) then
        top_line = tofloat(place)
   end if
   bot_line = top_line-space_lines
   edge_x = (/ 0.12,0.92,0.92,0.12,0.12/)
   
;   edge_x = (/ 0.22,0.76,0.76,0.22,0.22/)
   

   edge_y = (/ space_lines+bot_line,space_lines+bot_line,bot_line,bot_line,space_lines+bot_line/)
   edge_res = True             ; Indicate you want to set some resources.
   edge_res@gsLineColor="white"
   gsn_polyline_ndc(wks,edge_x,edge_y,edge_res)

   ;****************************************************************

   dim_lb = dimsizes(MODEL)

   gs1res               = True
   tx1res               = True
   tx1res@txFontHeightF = lin_sp/5*2

   xleg = new((/3*dim_lines,4/),"float")
   xl = new((/3,4/),"float")

   do i=0,2
      xl(0,:) = fspan((edge_x(0)+0.02),(edge_x(0)+0.07),4)
      xl(1,:) = fspan(((edge_x(2)+edge_x(0))/2-0.09),((edge_x(2)+edge_x(0))/2-0.04),4)
      xl(2,:) = fspan((edge_x(2)-0.20),(edge_x(2)-0.15),4)
   end do

   nt = 0
   do  j=0,dim_lines-1
      do i=0,2
         xleg(nt,:) = (/xl(i,:)/)
         nt = nt+1
      end do
   end do

   yleg = new((/3*dim_lines,4/),"float")
   yl = new(dim_lines,"float")

   do i=0,dim_lines-1
      yl(dim_lines-1-i) = edge_y(3)+(((edge_y(0) - (edge_y(3)))/dim_lines)*(i+1))-(lin_sp)
   end do

   nt = 0
   do j=0,dim_lines-1
      do i=0,2
         yleg(nt,:) = (/yl(j)/)
         nt = nt+1
      end do
   end do

   xtxt = new((/dim_lines*3/),"float")
   ytxt = new((/dim_lines*3/),"float")

   nt = 0
   do j=0,dim_lines-1
      do i=0,2
         xtxt(nt) = xl(i,3)+lin_sp/5.
         nt = nt+1
      end do
   end do
   nt = 0
   do i=0,dim_lines-1
      do j=0,2
         ytxt(nt) = (/yl(i)/)
         nt = nt+1
      end do
   end do

   do i=0,(dimsizes(MODEL)-1)
      gs1res@gsLineColor = colors(i)
      gs1res@gsLineDashPattern = dashs(i)
      gs1res@gsLineThicknessF=ticknesses(i)
      tx1res@txFontColor = colors(i)
      tx1res@txJust      = "CenterLeft"             ; change justification
      ;if (((MODEL(i).eq."HALOE").and.(dashs(i).eq.16)).or.((MODEL(i).eq."NCEP").and.(dashs(i).eq.16)).or.((MODEL(i).eq."UKMO").and.(dashs(i).eq.2))) then
     if (scatters(i).eq."Markers")       
         gs1res@gsMarkerColor = colors(i)
         gs1res@gsMarkerIndex = dashs(i)
         gs1res@gsMarkerSizeF = 0.01
         gs1res@gsMarkerThicknessF=1.5 
         gsn_polymarker_ndc(wks, xleg(i,2),yleg(i,2),gs1res)
      else

         gsn_polyline_ndc(wks, xleg(i,:),yleg(i,:),gs1res)
      end if

      gsn_text_ndc(wks,MODEL(i),xtxt(i),ytxt(i),tx1res)
   end do

end

;##################################################################################

;::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;##################################################################################
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

;:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
;****************************************************
; D. Shea
; Take a monthly climatology and make a daily climatology
; Current for gregorian / standard year.
;
; Supported: leftmost dimension must be 12
;            x(12), x(12,N), x(12,N1,N2), x(12,N1,N2,N3)
; x must have named dimensions on entry
; opt - not used set to zero [0]
;
undef("clmMon2clmDayn")
function clmMon2clmDayn( x:numeric, retOrder:integer, opt:integer )
local dNames, dimx, rank, X, midMon, day
begin
  if (.not.(retOrder.eq.0 .or. retOrder.eq.1)) then
      print("clmMon2clmDay: retOrder must be 0 or 1, retOrder=" +retOrder)
      exit
  end if

  dNames = getvardims( x )
  if (any(ismissing(dNames))) then
      print("clmMon2clmDay: named dimensions required:" +dNames)
      exit
  end if

  dimx   = dimsizes(x)
  if (dimx(0).ne.12) then
      print("clmMon2clmDay: leftmost dimension must be size=12: SIZE="+dimx(0))
      exit
  end if

  rank   = dimsizes( dimx )
  if (rank.gt.4) then
      print("clmMon2clmDay: currently max of 4 dimension supported: rank="+rank)
      exit
  end if

 ;if (isatt(x,"_FillValue")) then
 ;    nFill = num(ismissing(x))
 ;    if (nFill.gt.0) then
 ;        print("clmMon2clmDay: input is assumed to have no missing values, nFill="+nFill)
 ;        exit
 ;    end if
 ;end if
                          ; transfer to work arrsy,if necessary, reorder array
  if (rank.eq.1) then
      X =  x
  end if
  if (rank.eq.2) then
      X =  x($dNames(1)$|:, $dNames(0)$|:)
  end if
  if (rank.eq.3)
      X =  x($dNames(1)$|:, $dNames(2)$|:, $dNames(0)$|:)
  end if
  if (rank.eq.4)
      X =  x($dNames(1)$|:, $dNames(2)$|:, $dNames(3)$|:, $dNames(0)$|:)
  end if
                          ; mid day of each month
  if (isatt(opt,"midmon")) then
      if (dimsizes(opt@midMon).eq.12) then
          midMon = opt@midMon
      else
          print("clmMon2clmDay: midMon required to be size 12: size="+dimsizes(opt@midMon))
          exit
      end if
  else
      midMon = (/ 15.5, 45  , 74.5,105  ,135.5,166  \
                ,196.5,227.5,258  ,288.5,319  ,349.5/)
  end if
  midMon@long_name = "middle of month"

  day    = ispan(0,364,1)    ; use 0 => 364 for interpolation
  day!0  = "day"

  Z      = linint1_Wrap (midMon, X, True, day, 0)
  Z@info = "NCL: clmMon2clmDay"

  day    = ispan(1,365,1)    ; use 1 => 365 for coord variable
  day@long_name = "day of year: no leap"
  day@units     = "1=Jan 1, 32=Feb 1, ..., 365-Dec 31"
  Z!(rank-1) = "day"
  Z&day      =  day

  if (retOrder.eq.1) then
      return( Z )
  end if

  if (retOrder.eq.0) then
      if (rank.eq.1) then
          return(Z)
      end if
      if (rank.eq.2) then
          return( Z(day|:, $dNames(1)$|:) )
      end if
      if (rank.eq.3) then
          return( Z(day|:, $dNames(1)$|:, $dNames(2)$|:) )
      end if
      if (rank.eq.4) then
          return( Z(day|:, $dNames(1)$|:, $dNames(2)$|:, $dNames(3)$|:) )
      end if
  end if
end
