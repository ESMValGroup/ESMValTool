;*******************************************************************
; procedure to draw the Taylor Diagram Metrics Table
;
; AWWG METRICS 0.3 ;;;
; 27 Jun 06   ASP
;  5 Jul 06   DJS 
;
; Arguments:
;     mfname     : name of the output table
;     varNames   : variable names [metrics]
;     cases      : case (model) names
;     seasons    : season names
;     values     : array containing the values to be plotted
;     opt        : used to pass optional arguments

procedure taylor_metrics_table(mfname[1]:string        \       ; plot name
                              ,varNames[*]:string      \       ; variables
                              ,cases[*]:string         \
                              ,seasons[*]:string       \
                              ,values[*][*][*]:numeric  \
                              ,opt:logical              )


;;;
;;; Output of metrics to a table [grid]

begin  
    if (.not.isatt(opt,"pltType") .or. \
       (opt@pltType.eq."eps" .or. opt@pltType.eq."png" .or. opt@pltType.eq."gif")) then
        wks = gsn_open_wks("eps",mfname)    ; default 
    else
        wks = gsn_open_wks(opt@pltType,mfname)   
    end if
    
    xbox0  = (/0.0,1.0,1.0,0.0,0.0/) ; Box template.
    ybox0  = (/0.0,0.0,1.0,1.0,0.0/)
    
    nvar   = dimsizes(varNames)
    ncases = dimsizes(cases)
    nseas  = dimsizes(seasons)
    
    dimt   = dimsizes(values)        ; table dimensions
    
    if (ncases.ne.dimt(0)) then
       print("fatal:taylor_metrics_table fatal error: Number of case names ("+ncases+ \
             ") does not match the number of cases passed in ("+dimt(0)+")")
       status_exit(1)
    end if
    if (nseas.ne.dimt(1)) then
       print("fatal:taylor_metrics_table fatal error: Number of season names ("+nseas+ \
             ") does not match the number of seasons passed in ("+dimt(1)+")")
       status_exit(1)
    end if
    if (dimsizes(varNames).ne.dimt(2)) then
       print("fatal:taylor_metrics_table fatal error: Number of metric variables names ("+dimsizes(varNames)+ \
             ") does not match the number of metric variables passed in ("+dimt(2)+")")
       status_exit(1)
    end if
    
    ncolm1 = dimt(0)*dimt(1)
    if (ncolm1.le.12) then
       if (ncolm1.le.6) then
          tt_width = .4
       end if
       if (ncolm1.ge.7.and.ncolm1.le.12) then
          tt_width = .2
       end if
       tt_height     = 0.1
       tt_theight    = 0.02
       mn_height     = 0.03  
       mn_theight    = 0.012   
    else
       if (ncolm1.le.20) then
          tt_width   = .15
          tt_height  = 0.06
          tt_theight = 0.012
          mn_height  = 0.018  
          mn_theight = 0.0072
       else
          tt_width   = .1
          tt_height  = 0.05
          tt_theight = 0.01
          mn_height  = 0.015  
          mn_theight = 0.006
       end if
    end if
    if (ncolm1.gt.26.or.dimt(2).gt.63) then
       print("warning: Recommended maximium size of array passed into taylor_metrics_table is "+ \
             "26 columns (#cases * #seasons) by 63 rows (#variables), continuing")
    end if

    mn_width   = tt_width
    mv_height  = mn_height
    mv_width   = 0.1
    mv_theight = mn_theight
    
    tt_pRes = True
    tt_pRes@gsFillColor   = "CornflowerBlue"      ; background color for + values
    
    tt_tRes  = True
    tt_tRes@txFontHeightF = tt_theight
    
    mn_pRes = True
    mn_pRes@gsFillColor   = "Gray70"              ; background color for variable names
    
    mn_tRes = True
    mn_tRes@txFontHeightF = mn_theight
    
    mv_pRes = True
    mv_pRes@gsFillColor   = "White"               ; background color for reference values
    
    mv_tRes = True
    mv_tRes@txFontHeightF = mn_theight     
    
; Title
    
    xbox   = tt_width*xbox0
    ybox   = 1.-tt_height*ybox0
    
    ixtbox = xbox(0)+0.5*(xbox(1)-xbox(0))
    iytbox = ybox(0)+0.5*(ybox(2)-ybox(0))
    
    gsn_polygon_ndc(wks,xbox,ybox,tt_pRes)
    if (isatt(opt,"tableTitle") ) then
        gsn_text_ndc(wks, opt@tableTitle ,ixtbox,iytbox, tt_tRes)
    else
        gsn_text_ndc(wks,"CAM METRICS",ixtbox,iytbox, tt_tRes)
    end if
    gsn_polyline_ndc(wks,xbox,ybox,False)
    
    do im = 0,nvar-1
      ybox   = min(ybox)-(ybox0*mn_height)
      ixtbox = xbox(0)+0.5*(xbox(1)-xbox(0))
      iytbox = ybox(0)+0.5*(ybox(2)-ybox(0))
      
      gsn_polygon_ndc(wks,xbox,ybox,mn_pRes)
      gsn_polyline_ndc(wks,xbox,ybox,False)
      gsn_text_ndc(wks,varNames(im),ixtbox,iytbox, mn_tRes)
    end do
    
    do icase = 0, ncases-1
      ybox   = 1.-0.5*tt_height*ybox0
      xbox   = ((1.-tt_width)/ncases)*(xbox0+icase)+tt_width
      xboxi  = min(xbox)
      ixtbox = xbox(0)+0.5*(xbox(1)-xbox(0))
      iytbox = ybox(0)+0.5*(ybox(2)-ybox(0))
      gsn_polygon_ndc(wks,xbox,ybox,mn_pRes)
      gsn_polyline_ndc(wks,xbox,ybox,False)
      gsn_text_ndc(wks,cases(icase),ixtbox,iytbox, mn_tRes)
       
      do iseas = 0, nseas-1
      ybox     = 1.-0.5*tt_height-0.5*tt_height*ybox0
      xbox     = xboxi +  ((1.-tt_width)/(ncases*nseas))*(xbox0+iseas) 
      ixtbox   = xbox(0)+0.5*(xbox(1)-xbox(0))
      iytbox   = ybox(0)+0.5*(ybox(2)-ybox(0))
      gsn_polygon_ndc(wks,xbox,ybox,mn_pRes)
      gsn_polyline_ndc(wks,xbox,ybox,False)
      gsn_text_ndc(wks,seasons(iseas),ixtbox,iytbox, mn_tRes) 
      
      do im=0,nvar-1
          if (im .eq. 0) then
            ybox = 1.-tt_height-mv_height*ybox0
          else
            ybox = ybox - mv_height
          end if
          iytbox = ybox(0)+0.5*(ybox(2)-ybox(0))
          mv_pRes@gsFillColor = "White" 
          if (icase .gt. 0) then
            if (ismissing(values(icase,iseas,im)) .or. \    ; ???
                ismissing(values(icase,iseas,im)))then 
            else
              if (values(icase,iseas,im).le.values(0,iseas,im)) then
    		  if (isatt(opt,"color0")) then
                   mv_pRes@gsFillColor = opt@color0
                else
                   mv_pRes@gsFillColor = "DarkOliveGreen3"
    		  end if
              else
    		  if (isatt(opt,"color1")) then
                   mv_pRes@gsFillColor = opt@color1
    		  else	
    		   mv_pRes@gsFillColor = "IndianRed1"
    		  end if
              end if
            end if
          end if
          gsn_polygon_ndc(wks,xbox,ybox,mv_pRes)
          gsn_polyline_ndc(wks,xbox,ybox,False)
          gsn_text_ndc(wks,sprintf("%4.3f",values(icase,iseas,im)),ixtbox,iytbox, mv_tRes)     
        end do
      end do
    end do
    
    draw(wks)
  ;; activate the ;; lines if this is made a function
  ;;if (.not.isatt(opt,"gsnFrame") .or. opt@gsnFrame) then
        frame(wks)
                             ; if png or gif then use convert
        if (isatt(opt,"pltType")    .and. \ 
           (opt@pltType.eq."png" .or. opt@pltType.eq."gif")) then
            system("convert -trim +repage "+mfname+".eps "+mfname+"."+opt@pltType)
        end if
  ;;end if
    
end
