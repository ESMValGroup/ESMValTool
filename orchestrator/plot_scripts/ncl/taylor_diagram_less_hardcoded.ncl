; This is a donated script. It has no guarnteed support from the NCL project.
; See further comments below. For more info, and any improvements or bugs, 
; please contact Jatin Kala (Jatin.Kala.JK  *AT* gmail.com)
;===================================================================================
function taylor_diagram (wks:graphic, legendwks:graphic, RATIO[*][*]:numeric, CC[*][*]:numeric \
                                     ,rOpts:logical)

local dimR, nCase, nVar, X, Y, nc, angle, xyMin, xyOne, ceil_max_val, xyMax, xyMax_Panel, \
      FontHeightF, rxy, tempo_values, tempo_values1, temp_string, temp1_string, temp_string_y, \
      temp_string_y1, npts, xx, yy, sLabels, cLabels, rad, angC, taylor, rsrRes, dum0, dum1, \
      dum2, nStnRad, rr, radC, xC, yC, txRes, dum4, plRes, tmEnd, radTM, xTM, yTM, dum5, dum6, \
      i, mTM, angmTM, radmTM, dum7, angRL, rlRes, dum8, xRL, yRL, respl, dx, ncon, npts, ang, dum9, \
      dist_n, Markers, Colors, gsRes, ptRes, markerTxYOffset, dum10, dum11,n, lgres, nModel, lbid, \
      lbid2, amres, annoid1, nVar, varLabelsFontHeightF, txres, delta_y, ys, dum12 
;--------------------------------------------------------------------
; This version of taylor_diagram supports "paneling"
; It requires NCL version 4.2.0.a034 because it uses "gsn_create_legend"
; This version of the taylor_diagram function differs from the default provided 
; in that all variables are now passed to local (hopefully), and various
; bits and pieces of hard-coding of axis limits etc have been removed. 
; See further comments below. For more info, and any improvements or bugs, 
; please contact Jatin Kala (Jatin.Kala.JK@gmail.com)
;--------------------------------------------------------------------
;
; Generate a Taylor Diagram:
; Generate Multiple Aspects of Model Performance in a Single Diagram
; Taylor, K. E., J. Geophys. Res., 106, D7, 7183-7192, 2001
;
; An example:
; http://www.grida.no/climate/ipcc_tar/wg1/fig8-4.htm
;
; This expects one or more datasets. The left dimension 
; is the number of datasets. The rightmost is the number of pts.
;
; Markers are at: 
; http://www.ncl.ucar.edu/Document/Graphics/Resources/gs.shtml
;
; By default, the function can handle up to 10 variable comparisons..
; To expand ...  modify the 'Colors' and 'Markers' attributes.
; The user can change / add some default settings.
;
; The defaults that the user can modify:
;
; rOpts                 = True 
;                                  ; 'made-up' resources
; rOpts@Colors          =  (/ "blue" , "red", "green", "cyan", "black" \
;                           , "torquoise", "brown", "yellow"/)
; rOpts@Markers         =  (/ 2, 3, 6, 14, 9, 12, 7, 4/) ; Marker Indices
; rOpts@markerTxOffset  = 0.0175   ; offset for text above marker
; rOpts@stnRad          = (/ 1. /) ;  (/ 0.50, 0.75, 1.5 /) 
; rOpts@centerDiffRMS   = False    ;  True mean draw additional radii from REF 
; rOpts@caseLabelsFontHeightF = 0.05
; rOpts@varLabelsFontHeightF  = 0.013
; rOpts@varLabelsYloc         = 0.65
; rOpts@legendWidth           = 0.015
; rOpts@legendHeight          = 0.030*nCase
; rOpts@taylorDraw            = True
; rOpts@taylorFrame           = True
;
;                                  ; standard NCL resources
; rOpts@tiMainString    = "Taylor" ; not using title makes plot bigger
; rOpts@gsMarkerSizeF   = 0.0085   ; marker size   
; rOpts@gsMarkerThicknessF = 1.0
; rOpts@txFontHeightF   = 0.0125   ; text size 
; rOpts@tiMainFontHeightF = 0.0225 ; tiMainString size
;
; It returns to the user a graphic object containing the 
; Taylor background and plotted x/y pts.
; This graphic object contains a simple Taylor background appropriate
; for standardized data and the markers for the datasets.
; ==================================================================
; This version allows paneling:
;      The 'cumbersome' "dum" variables were added by 
;      Adam Phillips to allow paneling via "gsn_add_?".
;      Jatin Kala removed some hard-coding of axis-limits etc
;      1. xyMax = 1.65, 
;         changed to nearest x.5:
;         ceil_max_val          = ceil(max(RATIO))
;         if ceil_max_val - max(RATIO) .ge. 0.5
;             xyMax = ceil_max_val - 0.5
;          else
;             xyMax = ceil_max_val
;          end if
;      2. rxy@tmXBValues        = (/0.0,0.25,0.50,0.75,1.00,1.25,1.5/)    ; major tm
;         changed to: 
;         tempo_values        = decimalPlaces(fspan(0.0,ceil(max(RATIO)),toint(ceil(max(RATIO))/0.25)+1),2,True)
;         tempo_values1 = tempo_values(ind(tempo_values .lt. xyMax))
;         delete(tempo_values)
;         rxy@tmXBValues = tempo_values1
;         delete(tempo_values1)
;      3. rxy@tmXBLabels        = (/"    ","0.25","0.50","0.75","REF" ,"1.25","1.50"/)
;          changed to:
;          temp_string = tostring(sprintf("%4.2f",rxy@tmXBValues))
;          temp_string(0) = "    "
;          temp1_string = where(temp_string .eq. "1.00","REF",temp_string)
;          rxy@tmXBLabels = temp1_string
;          if (rOpts .and. isatt(rOpts,"OneX") )  then                     ; eg: rOpts@OneX="1.00" 
;             ;rxy@tmXBLabels        = (/"0.00","0.25","0.50","0.75",rOpts@OneX,"1.25","1.50"/)
;             ; rxy@tmXBLabels        = (/"    ","0.25","0.50","0.75",rOpts@OneX,"1.25","1.50"/)
;             rxy@tmXBLabels = temp_string
;          end if
;          delete(temp_string)
;          delete(temp1_string)
;      4.  rxy@tmYLValues        = (/0.0, .25,0.50, 0.75, 1.00, 1.25, 1.5/) ; major tm
;          changed to:
;          rxy@tmYLValues = rxy@tmXBValues
;      5.  rxy@tmYLLabels        = (/"0.00","0.25","0.50","0.75","1.00","1.25","1.50"/)
;          changed to:
;          temp_string_y = rxy@tmXBLabels
;          temp_string_y(0) = "0.00"
;          temp_string_y1 = where(temp_string_y .eq. "REF","1.00",temp_string_y)
;          delete(temp_string_y)
;          rxy@tmYLLabels = temp_string_y1
;          delete(temp_string_y1)
;       6. dum4 = gsn_add_text(wks,taylor,"Correlation",1.30,1.30,txRes)
;          changed to:
;          dum4 = gsn_add_text(wks,taylor,"Correlation",max(rxy@tmXBValues),max(rxy@tmXBValues),txRres)
;       7. Drawing of concentric RMS circles has been chaged, so circles are only drawn if there is enough space for a "full" semi-circle. Need more work here      
;          dist_n = dx*n + 1.0
;          if dist_n .gt. max(rxy@tmXBValues) + 0.25
;             ;  print("cric too big, not drawing")
;             ; need to find point of intersection between RMS circule and outer arc
;             ;  dum9(n-1) = gsn_add_polyline(wks,taylor,xx,yy,respl)
;           else
;               ;print("circ Ok")
;               dum9(n-1) = gsn_add_polyline(wks,taylor,xx,yy,respl)
;           end if
;           delete(dist_n)
; ==================================================================
begin
  dimR                  = dimsizes(RATIO)
  nCase                 = dimR(0)    ; # of cases [models] 
  nVar                  = dimR(1)    ; # of variables

                                     ; x/y coordinates for plotting
  X    = new ( (/nCase,nVar/) , typeof(RATIO) )
  Y    = new ( (/nCase,nVar/) , typeof(RATIO) )

  do nc=0,nCase-1
     angle      = acos( CC(nc,:) )   ; array operation                                    
     X(nc,:)    = RATIO(nc,:)*cos( angle )     
     Y(nc,:)    = RATIO(nc,:)*sin( angle )    
  end do

  ; fix to nearest x.5
  xyMin                 = 0.  
  xyOne                 = 1.00
; A_laue_ax+
; al  ceil_max_val          = ceil(max(RATIO))
; al  if ceil_max_val - max(RATIO) .ge. 0.5
; al        xyMax = ceil_max_val - 0.5
; al   else
; al        xyMax = ceil_max_val
; al   end if
  rmax = ceil(max(RATIO))
  if (rmax.lt.1.5) then
     rmax = 1.5
  end if
  if (rmax.gt.3.0) then
     rmax = 3.0
  end if
  xyMax = rmax + 0.1
; A_laue_ax-
  xyMax_Panel           = xyMax+ 0.10            ; paneling purposes
 
  if (rOpts .and. isatt(rOpts,"txFontHeightF"))  then 
      FontHeightF       = rOpts@txFontHeightF    ; user wants to specify size
  else
      FontHeightF       = 0.0175
  end if
 
; ----------------------------------------------------------------
; Part 1:
; base plot: Based upon request of Mark Stevens
; basic x-y and draw the 1.0 observed and the outer curve at 1.65
; ----------------------------------------------------------------

  rxy                   = True       
  rxy@gsnDraw           = False
  rxy@gsnFrame          = False
  rxy@vpHeightF         = 0.65
  rxy@vpWidthF          = 0.65
; A_laue_ax+
  if (rOpts .and. isatt(rOpts,"plotSize") ) then
     rxy@vpHeightF         = rOpts@plotSize
     rxy@vpWidthF          = rOpts@plotSize
     rxy@vpXF              = 0.125 ; left position of y-axis
     rxy@vpYF              = 0.825 ; top position of y-axis
  end if
; A_laue_ax-
  rxy@tmYLBorderOn      = False
  rxy@tmXBBorderOn      = False

  rxy@tiYAxisString     = "Standardized Deviations (Normalized)"
  rxy@tiYAxisFontHeightF= FontHeightF                        ; default=0.025 
  
  rxy@tmXBMode          = "Explicit" 
  ;rxy@tmXBValues        = (/0.0,0.25,0.50,0.75,1.00,1.25,1.5/)    ; major tm
                                                                  ; default  "OBS" or "REF"
; A_laue_ax+
  delta = 0.5
  if (rmax.le.1.5) then
     delta = 0.25
  end if
; al  tempo_values        = decimalPlaces(fspan(0.0,ceil(max(RATIO)),toint(ceil(max(RATIO))/0.25)+1),2,True)
  tempo_values        = decimalPlaces(fspan(0.0,rmax,toint(rmax/delta)+1),2,True)
; A_laue_ax-
  tempo_values1 = tempo_values(ind(tempo_values .lt. xyMax))
  delete(tempo_values)
  rxy@tmXBValues = tempo_values1
  delete(tempo_values1)
  ;rxy@tmXBLabels        = (/"0.00","0.25","0.50","0.75","REF" ,"1.25","1.50"/)
  ;rxy@tmXBLabels        = (/"    ","0.25","0.50","0.75","REF" ,"1.25","1.50"/)
  temp_string = tostring(sprintf("%4.2f",rxy@tmXBValues))
  temp_string(0) = "    "
  temp1_string = where(temp_string .eq. "1.00","REF",temp_string)
  rxy@tmXBLabels = temp1_string
  if (rOpts .and. isatt(rOpts,"OneX") )  then                     ; eg: rOpts@OneX="1.00" 
     ;rxy@tmXBLabels        = (/"0.00","0.25","0.50","0.75",rOpts@OneX,"1.25","1.50"/)
     ; rxy@tmXBLabels        = (/"    ","0.25","0.50","0.75",rOpts@OneX,"1.25","1.50"/)
       rxy@tmXBLabels = temp_string
  end if
  delete(temp_string)
  delete(temp1_string)

  rxy@tmXBMajorLengthF  = 0.015      ; default=0.02 for a vpHeightF=0.6
  rxy@tmXBLabelFontHeightF = FontHeightF
  rxy@tmXBMinorOn       = False
  rxy@trXMaxF           = xyMax_Panel

  rxy@tmYLMode          = "Manual"
  rxy@tmYLMinorOn       = False
  rxy@tmYLMajorLengthF  = rxy@tmXBMajorLengthF
  rxy@tmYLLabelFontHeightF = FontHeightF
  rxy@tmYLMode          = "Explicit" 
;  rxy@tmYLValues        = (/0.0, .25,0.50, 0.75, 1.00, 1.25, 1.5/) ; major tm
  rxy@tmYLValues = rxy@tmXBValues
;  rxy@tmYLLabels        = (/"0.00","0.25","0.50","0.75","1.00","1.25","1.50"/)
 ;rxy@tmYLLabels        = (/"    ","0.25","0.50","0.75","1.00","1.25","1.50"/)
  temp_string_y = rxy@tmXBLabels
  temp_string_y(0) = "0.00"
  temp_string_y1 = where(temp_string_y .eq. "REF","1.00",temp_string_y)
  delete(temp_string_y)
  rxy@tmYLLabels = temp_string_y1
  delete(temp_string_y1)
  rxy@trYMaxF           = xyMax_Panel

  rxy@tmYRBorderOn      = False
  rxy@tmYROn            = False      ; Turn off right tick marks.

  rxy@tmXTBorderOn      = False
  rxy@tmXTOn            = False      ; Turn off right tick marks.

  rxy@xyDashPatterns    = (/ 0 /)    ; line characteristics (dash,solid)
  rxy@xyLineThicknesses = (/ 2./)    ; choose line thickness

  rxy@gsnFrame          = False      ; Don't advance the frame.

                                            ; create outer 'correlation axis'
  npts    = 100                        ; arbitrary
  xx      = fspan(xyMin,xyMax,npts) 
  yy      = sqrt(xyMax^2 - xx^2    )   ; outer correlation line (xyMax)


  sLabels = (/"0.0","0.1","0.2","0.3","0.4","0.5","0.6" \ ; correlation labels
             ,"0.7","0.8","0.9","0.95","0.99","1.0"     /); also, major tm
  cLabels = stringtofloat(sLabels)
  rad     = 4.*atan(1.0)/180.
  angC    = acos(cLabels)/rad                     ; angles: correlation labels
                                                                       
  if (rOpts .and. isatt(rOpts,"tiMainString")) then
      rxy@tiMainString      = rOpts@tiMainString
     ;rxy@tiMainOffsetYF    = 0.015               ; default  0.0
      if (isatt(rOpts,"tiMainFontHeightF")) then
           rxy@tiMainFontHeightF = rOpts@tiMainFontHeightF
      else
           rxy@tiMainFontHeightF = 0.0225         ; default  0.025              
      end if
  end if
;;if (rOpts .and. isatt(rOpts,"gsnCenterString")) then
;;    rxy@gsnCenterString  = rOpts@gsnCenterString      ; only gsn_csm_xy
;;end if

  taylor  = gsn_xy(wks,xx,yy,rxy)                 ; Create and draw XY plot.

; A_laue_ax+
; --------------------------------------------------------------
; draw observational uncertainties given as normalized RMS error
; around the observational reference point at (1, 0)

  if (rOpts .and. isatt(rOpts, "rmsobs")) then
      npts = 100                                     ; arbitrary
      rmsobj = new(2, graphic)

      respl0                    = True
      respl0@gsLineThicknessF   = 4.0
      respl0@gsLineDashPattern  = 0
      respl0@gsLineColor        = (/0.50, 0.95, 0.47/) ;"green2"
      respl0@gsFillColor        = (/0.85, 1.00, 0.80/) ;"PaleGreen1"
      rr  = rOpts@rmsobs

      ; uncertainty circle covers whole plotting area
      if (rr.ge.xyMax) then
          ang  = fspan(180, 360, npts) * rad
          xx  = 0.995 * xyMax * cos(ang)
          yy  = fabs(0.995 * xyMax * sin(ang))
          rmsobj(0) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
          rmsobj(1) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
      else  ; uncertainty circle covers only part of the plotting area
          dist = rr + 1.0
          if (dist .gt. max(rxy@tmXBValues)) then ; uncertainty circle has to be cropped
              ; need to find point of intersection between uncertainty circule and outer arc
              xtmp = (xyMax^2 - rr^2 + 1.0) / 2.0
              ytmp = sqrt(xyMax^2 - xtmp^2)
              xtmp = xtmp - 1.0
              if (xtmp.eq.0.0) then
                  phitmp = 3.1415927
              else if (xtmp.gt.0.0) then
                  phitmp = 2.0 * 3.1415927 - atan(ytmp / xtmp)
              else
                  phitmp = 3.1415927 + atan(ytmp / xtmp)
              end if
              end if

              ang  = fspan(3.1415927, phitmp, npts)
              xx1 = 1.0 + rr * cos(ang)
              yy1 = fabs(rr * sin(ang))

              xtmp = xx1(npts - 1)
              ytmp = yy1(npts - 1)
              if (xtmp.eq.0.0) then
                  phitmp = 3.1415927
              else if (xtmp.gt.0.0) then
                  phitmp = 2.0 * 3.1415927 - atan(ytmp / xtmp)
              else
                  phitmp = 3.1415927 + atan(ytmp / xtmp)
              end if
              end if

              ang2 = fspan(phitmp, 0.0, npts)

              xx2 = 0.995 * xyMax * cos(ang2)
              yy2 = 0.995 * fabs(xyMax * sin(ang2))

              if (isvar("xx")) then
                  delete(xx)
              end if
              if (isvar("yy")) then
                  delete(yy)
              end if
              xx  = array_append_record(xx1, xx2, 0)
              yy  = array_append_record(yy1, yy2, 0)
              rmsobj(0) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
              rmsobj(1) = gsn_add_polyline(wks, taylor, xx1, yy1, respl0)
              delete(xx1)
              delete(xx2)
              delete(yy1)
              delete(yy2)
              delete(ang2)
          else  ; uncertainty circle does not need any adjustments
              ang  = fspan(180, 360, npts) * rad
              xx  = 1.0 + rr * cos(ang)
              yy  = fabs(rr * sin(ang))
              rmsobj(0) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
              rmsobj(1) = gsn_add_polyline(wks, taylor, xx, yy, respl0)
          end if  ; if uncertainty circle has to be cropped

          delete(ang)
          delete(xx)
          delete(yy)

      end if  ; if uncertainty circle covers whole plotting area
  end if  ; if an uncertainty estimate is available
; --------------------------------------------------------------
; A_laue_ax-

  rsrRes  = True
  rsrRes@gsLineThicknessF  = rxy@xyLineThicknesses(0)  ; line thickness
  rsrRes@gsLineDashPattern = 0                    ; solid line pattern
                                                  ; draw x and y to xyMax
  dum0 = gsn_add_polyline(wks,taylor,(/0.,  0. /),(/0.,xyMax/), rsrRes)
  dum1 = gsn_add_polyline(wks,taylor,(/0.,xyMax/),(/0.,  0. /), rsrRes)

  xx   = fspan(xyMin, xyOne ,npts)                ; draw 1.0 standard radius
  yy   = sqrt(xyOne - xx^2)   
  rsrRes@gsLineDashPattern = 1                    ; dashed line pattern
  rsrRes@gsLineThicknessF  = rxy@xyLineThicknesses(0)  ; line thickness
  dum2 = gsn_add_polyline(wks,taylor,xx,yy, rsrRes)
  delete(xx)
  delete(yy)
                                                  
  if (rOpts .and. isatt(rOpts,"stnRad") ) then
      rsrRes@gsLineThicknessF  = 1   ; rxy@xyLineThicknesses(0)  
      nStnRad = dimsizes(rOpts@stnRad)

      dum3  = new(nStnRad,graphic)
      do n=0,nStnRad-1
         rr = rOpts@stnRad(n)
; A_laue_ax+
         if (rr.gt.rmax) then
            continue
         end if
; A_laue_ax-
         xx = fspan(xyMin, rr ,npts) 
         yy = sqrt(rr^2   - xx^2)   
         dum3(n) = gsn_add_polyline(wks,taylor,xx,yy, rsrRes)
      end do
      taylor@$unique_string("dum")$ = dum3

      delete(xx)
      delete(yy)
  end if

  getvalues taylor                                ; get style info from taylor
    "tmYLLabelFont"        : tmYLLabelFont        ; use for correlation axis
    "tmYLLabelFontHeightF" : tmYLLabelFontHeightF
  end getvalues

; ----------------------------------------------------------------
; Part 2:
; Correlation labels
; ----------------------------------------------------------------
; A_laue_ax+
; al  radC    = xyMax                                  ; for correlation labels
  radC    = 1.02*xyMax                                 ; for correlation labels
; A_laue_ax-
  xC      = radC*cos(angC*rad)
  yC      = radC*sin(angC*rad)
; added to get some separation
; A_laue_ax+
; al  xC      = xC + 0.020*cos(rad*angC)
; al  yC      = yC + 0.060*sin(rad*angC)
; A_laue_ax-

  txRes               = True                      ; text mods desired
  txRes@txFontHeightF = FontHeightF               ; match YL 
  txRes@tmYLLabelFont = tmYLLabelFont             ; match YL
  txRes@txAngleF      = -45.
  if (.not.isatt(rOpts,"drawCorLabel") .or. rOpts@drawCorLabel) then 
      ;dum4 = gsn_add_text(wks,taylor,"Correlation",1.30,1.30,txRes)
; A_laue_ax+
; al       dum4 = gsn_add_text(wks,taylor,"Correlation",max(rxy@tmXBValues),max(rxy@tmXBValues),txRes)
       xdum = 1.15 * xyMax * cos(45.0*rad)
       ydum = 1.15 * xyMax * sin(45.0*rad)
       dum4 = gsn_add_text(wks,taylor,"Correlation",xdum,ydum,txRes)
; A_laue_ax-
	 taylor@$unique_string("dum")$ = dum4
  end if
  txRes@txAngleF      = 0.0 
  txRes@txFontHeightF = FontHeightF*0.50          ; bit smaller

;;dum0 = gsn_add_text(wks,taylor,"OBSERVED",1.00,0.075,txRes)

  plRes               = True
  plRes@gsLineThicknessF = 2.

; A_laue_ax+  
; al  txRes@txJust        = "CenterLeft"              ; Default="CenterCenter".
; A_laue_ax-
  txRes@txFontHeightF = FontHeightF               ; match YL 
 ;txRes@txBackgroundFillColor = "white"

  tmEnd = 0.975
  radTM = xyMax*tmEnd                             ; radius end: major TM 
  xTM   = new( 2 , "float")
  yTM   = new( 2 , "float")

  dum5 = new(dimsizes(sLabels),graphic)
  dum6 = dum5

  do i=0,dimsizes(sLabels)-1                      ; Loop to draw strings
    txRes@txAngleF = angC(i)
; A_laue_ax+
    if (angC(i).le.45.0) then
       txRes@txJust        = "CenterLeft"              ; Default="CenterCenter".
    else
       txRes@txJust        = "BottomCenter"              ; Default="CenterCenter".
    end if
; A_laue_ax-
    dum5(i) = gsn_add_text(wks, taylor, sLabels(i),xC(i),yC(i),txRes) ; cor label
    xTM(0)   = xyMax*cos(angC(i)*rad)             ; major tickmarks at
    yTM(0)   = xyMax*sin(angC(i)*rad)             ; correlation labels
    xTM(1)   = radTM*cos(angC(i)*rad)             
    yTM(1)   = radTM*sin(angC(i)*rad)
    dum6(i) = gsn_add_polyline(wks,taylor,xTM,yTM,plRes)
  end do
                                                  ; minor tm locations
  mTM     = (/0.05,0.15,0.25,0.35,0.45,0.55,0.65 \ 
             ,0.75,0.85,0.91,0.92,0.93,0.94,0.96,0.97,0.98  /)
  angmTM  = acos(mTM)/rad                         ; angles: correlation labels
  radmTM  = xyMax*(1.-(1.-tmEnd)*0.5)             ; radius end: minor TM 

  dum7 = new(dimsizes(mTM),graphic)

  do i=0,dimsizes(mTM)-1                          ; manually add tm
    xTM(0)   = xyMax*cos(angmTM(i)*rad)           ; minor tickmarks
    yTM(0)   = xyMax*sin(angmTM(i)*rad)
    xTM(1)   = radmTM*cos(angmTM(i)*rad)          
    yTM(1)   = radmTM*sin(angmTM(i)*rad)
    dum7(i)  = gsn_add_polyline(wks,taylor,xTM,yTM,plRes)
  end do
                                                  ; added for Wanli
  if (rOpts .and. isatt(rOpts,"ccRays") ) then
      angRL = acos(rOpts@ccRays)/rad             ; angles: radial lines

      rlRes = True
      rlRes@gsLineDashPattern= 2  ; line pattern
      rlRes@gsLineThicknessF = 1  ; choose line thickness
      if (isatt(rOpts,"ccRays_color")) then
; A_laue_ax+
; al          rlRes@gsLineColor    =  "LightGray"
          rlRes@gsLineColor    =  rOpts@ccRays_color
; A_laue_ax-
      end if

      dum8 = new(dimsizes(angRL),graphic)
      do i=0,dimsizes(angRL)-1
         xRL     = xyMax*cos(angRL(i)*rad)
         yRL     = xyMax*sin(angRL(i)*rad)
         dum8(i) = gsn_add_polyline(wks,taylor,(/0, xRL /),(/0,  yRL  /),rlRes)
      end do
      taylor@$unique_string("dum")$ = dum8
  end if
  
; ----------------------------------------------------------------
; Part 3:
; Concentric about 1.0 on XB axis
; ----------------------------------------------------------------
  if (rOpts .and. isatt(rOpts,"centerDiffRMS") \
            .and. rOpts@centerDiffRMS) then
      respl                    = True                ; polyline mods desired
;al      respl@xyLineThicknessF   = 1.0                 ; line thickness
;al      respl@xyLineDashPattern  = 2                   ; short dash lines
      respl@gsLineThicknessF   = 2.0                 ; line thickness
;al      respl@gsLineDashPattern  = 2                   ; short dash lines
      respl@gsLineColor        = "Black"             ; line color     
      if (isatt(rOpts,"centerDiffRMS_color")) then
; A_laue_ax+
; al          respl@gsLineColor    =  "LightGray"
          respl@gsLineColor    =  rOpts@centerDiffRMS_color
; A_laue_ax-
      end if

; A_laue_ax+
; al      dx   = 0.25
      dx   = delta
      ncon = 4                                       ; 0.75, 0.50, 0.25, 0.0
; A_laue_ax-
      npts = 100                                     ; arbitrary
      ang  = fspan(180,360,npts)*rad

      dum9 = new(ncon,graphic)

      do n=1,ncon 
         rr  = n*dx            ; radius from 1.0 [OBS] abscissa
         xx  = 1. + rr*cos(ang)
         yy  = fabs( rr*sin(ang) )
          dist_n = dx*n + 1.0
          if dist_n .gt. max(rxy@tmXBValues) ; + 0.25
             ;  print("cric too big, not drawing")
             ; need to find point of intersection between RMS circule and outer arc
; A_laue_ax+
             xtmp = (xyMax^2 - rr^2 + 1.0) / 2.0
             ytmp = sqrt(xyMax^2 - xtmp^2)
             xtmp = xtmp - 1.0
             if (xtmp.eq.0.0) then
                 phitmp = 3.1415927
             else if (xtmp.gt.0.0) then
                 phitmp = 2.0 * 3.1415927 - atan(ytmp / xtmp)
             else
                 phitmp = 3.1415927 + atan(ytmp / xtmp)
             end if
             end if
             ang2  = fspan(3.1415927, phitmp, npts)
             xx2  = 1.0 + rr * cos(ang2)
             yy2  = fabs(rr * sin(ang2))
             dum9(n - 1) = gsn_add_polyline(wks, taylor, xx2, yy2, respl)
             delete(ang2)
             delete(xx2)
             delete(yy2)
; A_laue_ax-
           ;    dum9(n-1) = gsn_add_polyline(wks,taylor,xx,yy,respl)
           else
               ;print("circ Ok")
               dum9(n-1) = gsn_add_polyline(wks,taylor,xx,yy,respl)
           end if 
           delete(dist_n)

     ;    if (n.le.2) then
     ;        dum9(n-1) = gsn_add_polyline(wks,taylor,xx,yy,respl)
     ;    end if
     ;    if (n.eq.3) then
     ;        n3 = floattointeger( 0.77*npts ) 
     ;        dum9(n-1) = gsn_add_polyline(wks,taylor,xx(0:n3),yy(0:n3),respl)
     ;    end if
     ;    if (n.eq.4) then
     ;        n4 = floattointeger( 0.61*npts ) 
     ;        dum9(n-1) = gsn_add_polyline(wks,taylor,xx(0:n4),yy(0:n4),respl)
     ;    end if
      end do
      delete(ang)
      delete(xx)
      delete(yy)
      taylor@$unique_string("dum")$ = dum9

  end if

; A_laue_ax+
  ; ------------------------------------------------------------------------
  ; draw filled black circle at (1, 0) to mark the reference point
  if (rOpts .and. isatt(rOpts,"gsMarkerSizeF")) then
      gsMarkerSizeF = rOpts@gsMarkerSizeF
  else
      gsMarkerSizeF = 0.0085
  end if
  gsRes = True
  gsRes@gsMarkerIndex = 16
  gsRes@gsMarkerColor = "Black"
  gsRes@gsMarkerSizeF = gsMarkerSizeF
  getvalues taylor
    "vpXF"      : vpx
    "vpYF"      : vpy
    "vpWidthF"  : vpw
    "vpHeightF" : vph
  end getvalues
  dotx = vpx + 1.0 / xyMax_Panel * vpw
  doty = vpy - vph
  gsn_polymarker_ndc(wks, dotx, doty, gsRes)  ; half circle below the x-axis
  dotobj = new(1, graphic)
  dotobj = gsn_add_polymarker(wks, taylor, 1.0, 0.0, gsRes)  ; half circle above the x-axis
  ; ------------------------------------------------------------------------
; A_laue_ax-

; ---------------------------------------------------------------
; Part 4:
; generic resources that will be applied to all users data points
; of course, these can be changed 
; http://www.ncl.ucar.edu/Document/Graphics/Resources/gs.shtml
; ---------------------------------------------------------------
  if (rOpts .and. isatt(rOpts,"Markers")) then
      Markers = rOpts@Markers
  else
      Markers = (/ 4, 6, 8,  0, 9, 12, 7, 2, 11, 16/) ; Marker Indices
  end if

  if (rOpts .and. isatt(rOpts,"Colors")) then
      Colors  = rOpts@Colors
  else
      Colors  = (/ "red", "blue", "green", "cyan", "orange" \
                 , "torquoise", "brown", "yellow", "purple", "black"/)
  end if

  if (rOpts .and. isatt(rOpts,"gsMarkerThicknessF")) then
      gsMarkerThicknessF = rOpts@gsMarkerThicknessF
  else
      gsMarkerThicknessF = 2.0
  end if

  if (rOpts .and. isatt(rOpts,"gsMarkerSizeF")) then
      gsMarkerSizeF      = rOpts@gsMarkerSizeF
  else
      gsMarkerSizeF      = 0.0085                  ; Default: 0.007
  end if

  gsRes = True
  gsRes@gsMarkerThicknessF = gsMarkerThicknessF      ; default=1.0
  gsRes@gsMarkerSizeF      = gsMarkerSizeF           ; Default: 0.007 

  ptRes = True                        ; text options for points
  ptRes@txJust             = "BottomCenter"; Default="CenterCenter".
  ptRes@txFontThicknessF   = 1.2      ; default=1.00
  ptRes@txFontHeightF      = 0.0125   ; default=0.05
  if (rOpts .and. isatt(rOpts,"txFontHeightF")) then
      ptRes@txFontHeightF  = rOpts@txFontHeightF  
  end if

  markerTxYOffset          = 0.0175   ; default
  if (rOpts .and. isatt(rOpts,"markerTxYOffset")) then
      markerTxYOffset = rOpts@markerTxYOffset             ; user defined offset
  end if

  dum10 = new((nCase*nVar),graphic)
  dum11 = dum10

  do n=0,nCase-1
     gsRes@gsMarkerIndex   = Markers(n)             ; marker style (+)
     gsRes@gsMarkerColor   = Colors(n)              ; marker color
     ptRes@txFontColor     = gsRes@gsMarkerColor
    do i=0,nVar-1
       dum10(n*nVar+i) = gsn_add_polymarker(wks,taylor,X(n,i),Y(n,i),gsRes) 
; A_laue_ax+
;al       print(rOpts@caseLabels(n))
;al       print("X(n,i) = " + X(n,i) + " " + "Y(n,i) = " + Y(n,i)) 
;al       print("CC(n,i) = " + CC(n,i) + " " + "RATIO(n,i) = " + RATIO(n,i))
;al       print("---------------------------------------------------------------------------")
; al       dum11(n*nVar+i) = gsn_add_text(wks,taylor,(i+1),X(n,i),Y(n,i)+markerTxYOffset,ptRes)
       printdl = True
       if (rOpts .and. isatt(rOpts,"printDataLabels")) then
          if (.not.rOpts@printDataLabels) then
              printdl = False
          end if
       end if
       if (printdl) then
          dum11(n*nVar+i) = gsn_add_text(wks,taylor,(i+1),X(n,i),Y(n,i)+markerTxYOffset,ptRes)
       end if
; A_laue_ax-
    end do
  end do

; ---------------------------------------------------------------
; Part 5: ; add case legend and variable labels
; ---------------------------------------------------------------

  extrafile = False

  if (rOpts .and. isatt(rOpts, "legendExtraFile")) then
      if (rOpts@legendExtraFile) then
          extrafile = True
      end if
  end if

  if (rOpts .and. isatt(rOpts,"caseLabels")) then 

      if (isatt(rOpts,"caseLabelsFontHeightF")) then
          caseLabelsFontHeightF = rOpts@caseLabelsFontHeightF
      else
          caseLabelsFontHeightF = 0.05  
      end if

      lgres                    = True
      lgres@lgMarkerColors     = Colors        ; colors of markers
      lgres@lgMarkerIndexes    = Markers       ; Markers 
      lgres@lgMarkerSizeF      = gsMarkerSizeF ; Marker size
      lgres@lgItemType         = "Markers"     ; draw markers only
      lgres@lgLabelFontHeightF = caseLabelsFontHeightF  ; font height of legend case labels

      nModel                   = dimsizes( rOpts@caseLabels )

      if (isatt(rOpts,"legendWidth")) then
          lgres@vpWidthF       = rOpts@legendWidth
      else
          lgres@vpWidthF       = 0.15           ; width of legend (NDC)
      end if

      if (nModel.gt.20) then
          lgres@vpWidthF = lgres@vpWidthF * 0.5
          lgres@lgMarkerSizeF = lgres@lgMarkerSizeF * 0.75
      end if

      lgres@lgPerimOn          = False         ; turn off perimeter
; A_laue_ax+
      if (isatt(rOpts,"legendBox")) then
          if (rOpts@legendBox) then
             lgres@lgPerimOn = True
             lgres@lgRightMarginF = 0.1
          end if
      end if
; A_laue_ax-
;      lbid = gsn_create_legend(legendwks,nModel,rOpts@caseLabels,lgres)

    if (nModel.le.20) then
        if (isatt(rOpts,"reverseOrder")) then
            if (rOpts@reverseOrder) then
                lgres@lgItemOrder = ispan(nModel-1,0,1)
            end if
        end if
        if (isatt(rOpts,"legendHeight")) then
            lgres@vpHeightF      = rOpts@legendHeight
        else   
            lgres@vpHeightF      = 0.030*nCase
        end if
        lbid = gsn_create_legend(legendwks,nModel,rOpts@caseLabels,lgres)
        amres = True
        if (isatt(rOpts,"legendXpos")) then
            amres@amParallelPosF = rOpts@legendXpos
            amres@amJust = "TopLeft"
        else
            amres@amParallelPosF = 0.35
        end if
        if (isatt(rOpts,"legendYpos")) then
            amres@amOrthogonalPosF = rOpts@legendYpos
            amres@amJust = "TopLeft"
        else
            amres@amOrthogonalPosF = -0.35
        end if
        if (extrafile) then
            gsn_legend_ndc(legendwks, nModel, rOpts@caseLabels(0:nModel-1), 0.78, 0.75, lgres)
        else
            annoid1 = gsn_add_annotation(taylor,lbid,amres)       ; add legend to plot
        end if
    else
        if (isatt(rOpts,"modelsperrow")) then
            modelsperrow = rOpts@modelsperrow
        else
            modelsperrow = (nModel + 1 ) / 2
        end if
        do lr = 0, nModel / modelsperrow
            ix0 = lr * modelsperrow
            if (ix0.ge.nCase) then
                break
            end if
            ix1 = min((/nCase - 1, ix0 + modelsperrow - 1/))
            dumdelta = ix1 - ix0 + 1

            if (isatt(rOpts,"reverseOrder")) then
                if (rOpts@reverseOrder) then
                    if (isatt(lgres, "lgItemOrder")) then
                        delete(lgres@lgItemOrder)
                    end if
                    lgres@lgItemOrder = ispan(dumdelta - 1, 0, 1)
                end if
            end if

            lgres@vpHeightF      = 0.026 * dumdelta

            delete(lgres@lgMarkerColors)
            delete(lgres@lgMarkerIndexes)
;            delete(lgres@lgMarkerSizeF)
            lgres@lgMarkerColors  = Colors(ix0:ix1)
            lgres@lgMarkerIndexes = Markers(ix0:ix1)
            lgres@lgPerimOn = False ; no line around the legend(s)

            if (extrafile) then
                legend_x = 0.05 + 0.11 * lr
            else
                legend_x = 0.78 + 0.11 * lr
            end if

            gsn_legend_ndc(legendwks, dumdelta, rOpts@caseLabels(ix0:ix1), \
                           legend_x, 0.75, lgres)
        end do
    end if
; A_laue_ax-

;      amres = True
; A_laue_ax+
; al      amres@amParallelPosF     =  0.35           
; al      amres@amOrthogonalPosF   = -0.35             
;      if (isatt(rOpts,"legendXpos")) then
;         amres@amParallelPosF = rOpts@legendXpos
;         amres@amJust = "TopLeft"
;      else
;         amres@amParallelPosF = 0.35
;      end if
;      if (isatt(rOpts,"legendYpos")) then
;         amres@amOrthogonalPosF = rOpts@legendYpos
;         amres@amJust = "TopLeft"
;      else
;         amres@amOrthogonalPosF = -0.35
;      end if
; A_laue_ax-
;      annoid1 = gsn_add_annotation(taylor,lbid,amres)       ; add legend to plot
;      if (nModel.gt.20) then
;          amres@amParallelPosF = amres@amParallelPosF + lgres@vpWidthF * 1
;          annoid1 = gsn_add_annotation(taylor,lbid2,amres)  ; add second legend to plot
;      end if
  end if

  if (rOpts .and. isatt(rOpts,"varLabels")) then 
      nVar    = dimsizes(rOpts@varLabels)

      if (isatt(rOpts,"varLabelsFontHeightF")) then
          varLabelsFontHeightF = rOpts@varLabelsFontHeightF
      else
          varLabelsFontHeightF = 0.013
      end if

      txres = True
      txres@txFontHeightF = varLabelsFontHeightF
      txres@txJust = "CenterLeft"              ; justify to the center left

     ;delta_y = 0.02       
      delta_y = 0.06   
      if (rOpts .and. isatt(rOpts,"varLabelsYloc")) then
          ys  = rOpts@varLabelsYloc            ; user specified
      else
          ys  = max( (/nVar*delta_y , 0.30/) )
      end if

      
      do i = 1,nVar     
         if (i.eq.1) then
             dum12 = new(nVar,graphic)
	 end if

         dum12(i-1) = gsn_add_text(wks,taylor,i+" - "+rOpts@varLabels(i-1), .125,ys,txres)
         ys = ys- delta_y
      end do

      taylor@$unique_string("dum")$ = dum12
  end if

  taylor@$unique_string("dum")$ = dum0   ; x-axis
  taylor@$unique_string("dum")$ = dum1   ; y-axis
  taylor@$unique_string("dum")$ = dum2   ; 1.0 std curve
  taylor@$unique_string("dum")$ = dum5   ; labels [COR]
  taylor@$unique_string("dum")$ = dum6   ; major tm [COR]
  taylor@$unique_string("dum")$ = dum7   ; minor tm
  taylor@$unique_string("dum")$ = dum10  ; markers
  taylor@$unique_string("dum")$ = dum11  ; text

  if (.not.isatt(rOpts,"taylorDraw") .or. \
     (isatt(rOpts,"taylorDraw") .and. rOpts@taylorDraw)) then 
	draw(taylor)
  end if
  if (.not.isatt(rOpts,"taylorFrame") .or. \
     (isatt(rOpts,"taylorFrame") .and. rOpts@taylorFrame)) then 
	frame(wks)
  end if

  return(taylor)
end
