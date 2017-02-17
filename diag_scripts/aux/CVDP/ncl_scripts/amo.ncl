; Calculates the AMO pattern, timeseries, and spectra.
;
; Variables used: TS
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: amo.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")
  
  nsim = numAsciiRow("namelist_byvar/namelist_ts")
  na = asciiread("namelist_byvar/namelist_ts",(/nsim/),"string")
  names = new(nsim,"string")
  paths = new(nsim,"string")
  syear = new(nsim,"integer",-999)
  eyear = new(nsim,"integer",-999)
  delim = "|"

  do gg = 0,nsim-1
     names(gg) = str_strip(str_get_field(na(gg),1,delim))
     paths(gg) = str_strip(str_get_field(na(gg),2,delim))
     syear(gg) = stringtointeger(str_strip(str_get_field(na(gg),3,delim)))
     eyear(gg) = stringtointeger(str_strip(str_get_field(na(gg),4,delim)))
  end do
  nyr = eyear-syear+1
  nyr_max = max(nyr)

  pi=4.*atan(1.0)
  rad=(pi/180.)
         
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks = gsn_open_wks(wks_type,getenv("OUTDIR")+"amo")
  wks2 = gsn_open_wks(wks_type,getenv("OUTDIR")+"amo.powspec")
  wks3 = gsn_open_wks(wks_type,getenv("OUTDIR")+"amo.timeseries")  
  if (COLORMAP.eq."0") then
     gsn_define_colormap(wks,"ncl_default")   
     gsn_define_colormap(wks2,"cb_9step")
     gsn_define_colormap(wks3,"ncl_default")
  end if
  if (COLORMAP.eq."1") then
     gsn_define_colormap(wks,"BlueDarkRed18")   
     gsn_define_colormap(wks2,"cb_9step")
     gsn_define_colormap(wks3,"ncl_default")
  end if
  
  map = new(nsim,"graphic")  
  pspec = new(nsim,"graphic")
  xyplot = new(nsim,"graphic")
  xyplot2 = new(nsim,"graphic")
  if (isfilepresent2("obs_ts")) then
     pspec_obs = new(nsim,"graphic")
  end if
 
  do ee = 0,nsim-1
     sstT = data_read_in(paths(ee),"TS",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(sstT,"is_all_missing")) then
        delete(sstT)
        continue
     end if  
     if (isfilepresent2("obs_ts").and.ee.eq.0) then   ; mask out land in array going into EOF functions
     else
        sstT = where(sstT.le.-1.8,-1.8,sstT)
        d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
        basemap = d->LSMASK
        lsm = landsea_mask(basemap,sstT&lat,sstT&lon)
        sstT = mask(sstT,conform(sstT,lsm,(/1,2/)).ge.1,False)
        delete([/lsm,basemap/])
        delete(d)
     end if    
     
     sstT = lonFlip(sstT)   ; orient longitudes from -180:180 (set to 0:360 in data_read_in function)
     if (OPT_CLIMO.eq."Full") then
        sstT = rmMonAnnCycTLL(sstT)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = sstT
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(sstT&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        sstT   = calcMonAnomTLL(sstT,climo) 
        delete(climo)
     end if
          
     coswgt=cos(rad*sstT&lat)
     coswgt!0 = "lat"
     coswgt&lat= sstT&lat
     natl_aa =   wgt_areaave(sstT(:,{0:60},{-80:0}),coswgt({0.:60.}),1.0,0)
     global_aa = wgt_areaave(sstT(:,{-60:60},:),coswgt({-60.:60.}),1.0,0)
     
     finarr = new((/2,dimsizes(natl_aa)/),"float",-999.)    ; timeseries plot
     finarr!1 = "time"
     finarr&time = sstT&time
     finarr(0,:) = (/ natl_aa - global_aa /)
     finarr(1,:) = (/ runave(finarr(0,:),61,0) /)
     delete(natl_aa)
     delete(global_aa)  

     finreg = sstT(0,:,:)
     finreg = (/ regCoef(finarr(0,:),sstT(lat|:,lon|:,time|:))  /) 
     delete([/sstT/])     
     do gg = 0,2
        finreg = (/ smth9(finreg,0.5,0.25,True) /)
     end do
     delete([/coswgt/])
     finreg@syear = syear(ee)
     finreg@eyear = eyear(ee) 
;---------------------------------------------------------------------------------------------     
     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        fn = getenv("OUTDIR")+modname+".cvdp_data."+syear(ee)+"-"+eyear(ee)+".nc"
        if (.not.isfilepresent2(fn)) then
           z = addfile(fn,"c")
           z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
           z@notes = "Data from "+names(ee)+" from "+syear(ee)+"-"+eyear(ee)
           if (OPT_CLIMO.eq."Full") then
              z@climatology = syear(ee)+"-"+eyear(ee)+" climatology removed prior to all calculations (other than means)"
           else
              z@climatology = CLIMO_SYEAR+"-"+CLIMO_EYEAR+" climatology removed prior to all calculations (other than means)"
           end if
        else
           z = addfile(fn,"w")
        end if
        z->amo_pattern_mon = finreg
        amo_ts = finarr(0,:)
        amo_ts@units = "C"
        z->amo_timeseries_mon = amo_ts
        delete([/modname,fn,amo_ts/])
     end if     
;------------------------------------------------------------------------
     iopt = 0
     jave = (7*nyr(ee))/100
     val1 = .95
     val2 = .99
     if (jave.eq.0) then
       jave = 1
     end if
     pct = 0.1     
     spectra_mvf = False        ; missing value flag for nino3.4
     if (any(ismissing(finarr(0,:)))) then
        print("Missing data detected for "+names(ee)+", power spectra function does not allow missing data, not creating AMO spectra")
        spectra_mvf = True
        if (isfilepresent2("obs_ts").and.ee.eq.0) then
           spectra_mvf_obs = True    ; missing value flag for obs nino3.4
        end if
     else
        if (isfilepresent2("obs_ts").and.ee.eq.0) then
           spectra_mvf_obs = False   ; missing value flag for obs nino3.4
        end if
        sdof = specx_anal(dim_standardize(finarr(0,:),0),iopt,jave,pct)
        splt1 = specx_ci(sdof,val1,val2)
        if (OUTPUT_DATA.eq."True") then
           splt1!0 = "ncurves"
           splt1&ncurves = ispan(0,3,1)
           splt1!1 = "frequency"
           splt1&frequency = sdof@frq
           splt1@units_info = "df refers to frequency interval; data are standardized so there are no physical units"
           splt1@units = "1/df"
           splt1@info = "(0,:)=spectrum,(1,:)=Markov red noise spectrum, (2,:)="+val1+"% confidence bound for Markhov, (3,:)="+val2+"% confidence bound for Markhov"
           z->amo_spectra = splt1
        end if
        if (isfilepresent2("obs_ts").and.ee.eq.0) then
           sdof_obs = sdof
        end if
        delete([/iopt,jave,pct/])
     end if
     if (isvar("z")) then
        delete(z)
     end if
;=========================================================================          
     res = True
     res@mpProjection = "WinkelTripel"
     res@mpGeophysicalLineColor = "gray42"
     res@mpPerimOn    = False
     res@mpGridLatSpacingF =  90            ; change latitude  line spacing
     res@mpGridLonSpacingF = 180.           ; change longitude line spacing
     res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
     res@mpFillOn = False
     res@mpCenterLonF = 0.
     res@mpOutlineOn = True  
     res@mpDataSetName = "Earth..4"
     res@gsnDraw      = False
     res@gsnFrame     = False
     res@vpYF = 0.95
     res@vpHeightF = 0.3
     res@vpXF = 0.2
     res@vpWidthF = 0.6
  
;     res@cnFillMode = "RasterFill"
     res@cnLevelSelectionMode = "ExplicitLevels"
     if (COLORMAP.eq."0") then
        res@cnLevels = fspan(-4.,4.,21)
     end if
     if (COLORMAP.eq."1") then
        res@cnLevels = fspan(-3.2,3.2,17)
     end if

     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     res@gsnRightStringOrthogonalPosF = -0.05
     res@gsnRightStringParallelPosF = 0.96
     res@gsnLeftStringOrthogonalPosF = -0.05
     res@gsnLeftStringParallelPosF = 0.005
     res@gsnLeftStringFontHeightF = 0.014
     res@gsnCenterStringFontHeightF = 0.018
     res@gsnRightStringFontHeightF = 0.014   
     res@gsnLeftString = syear(ee)+"-"+eyear(ee)
     res@gsnCenterString = names(ee)
     res@gsnRightString = ""
     
     if (isfilepresent2("obs_ts").and.ee.eq.0) then    ; for pattern correlation table
        patcor = new((/nsim,dimsizes(finreg&lat),dimsizes(finreg&lon)/),typeof(finreg))
        patcor!1 = "lat"
        patcor&lat = finreg&lat
        patcor!2 = "lon"
        patcor&lon = finreg&lon
        patcor(ee,:,:) = (/ finreg /)
     end if
     if (isfilepresent2("obs_ts").and.ee.ge.1) then
        patcor(ee,:,:) = (/ linint2(finreg&lon,finreg&lat,finreg,True,patcor&lon,patcor&lat,0) /)
     end if

     map(ee) = gsn_csm_contour_map(wks,finreg,res)
     delete([/finreg/])
     
     pres = True
     pres@vpXF = 0.07
     pres@trYMinF = 0.
     pres@trXMinF = 0.0
;     pres@trYMaxF = 82.
     pres@trXMaxF = 0.0832
     pres@tiYAxisString = "Power"              ; yaxis
     pres@xyLineColor   =  "black"         
     pres@gsnFrame      = False
     pres@gsnDraw       = False
  
     pres@tmXBLabelDeltaF = -.8
     pres@tmXTLabelDeltaF = -.8
     pres@pmLegendDisplayMode    = "Never"
     pres@xyLineThicknesses   = (/3.5,2.,1.,1./) 
     pres@xyDashPatterns      = (/0,0,0,0/)
     pres@xyLineColors        = (/"foreground","red","blue","green"/)      
     pres@xyLabelMode = "custom"
     pres@xyLineLabelFontColors = pres@xyLineColors 
     pres@xyExplicitLabels = (/"","",val1*100+"%",val2*100+"%"/)
     pres@tmXTOn = True
     pres@tmYROn = False
     pres@tmXTLabelsOn = True
     pres@tmXUseBottom = False
     pres@tmXTMode   = "Explicit"  
     pres@tmXBMode   = "Explicit"            
     pres@tmXTValues = (/".00167",".00833",".01667",".02778",".0416",".0556",".0832"/)
     pres@tmXTLabels = (/"50","10","5","3","2","1.5","1"/)           
     pres@tmXBValues = (/".0",".01",".02",".03",".042",".056",".083"/)
     pres@tmXBLabels =     pres@tmXBValues
     pres@tmXTLabelFontHeightF = 0.018
     pres@tmXBLabelFontHeightF = 0.018
     pres@tmYLLabelFontHeightF = 0.018
     pres@tiYAxisString = "Variance" ;"Power (~S~o~N~C~S~2~N~ / cycles mo~S~-1~N~)"              ; yaxis
     pres@tiXAxisString = "Frequency (cycles mo~S~-1~N~)"
     pres@tiMainString = ""
     pres@txFontHeightF = 0.015
     pres@xyLineLabelFontHeightF = 0.022
     pres@tiXAxisFontHeightF = 0.025
     pres@tiYAxisFontHeightF = 0.025
     pres@tiMainFontHeightF = 0.03
     
     pres@tiMainOn = False
     pres@gsnCenterString = "Period (years)"
     pres@gsnCenterStringFontHeightF =    pres@tiYAxisFontHeightF
     pres@gsnRightString = ""
     pres@gsnLeftString = "" 
     if (wks_type.eq."png") then
        pres@xyLineThicknessF   = 3.5
        res@mpGeophysicalLineThicknessF = 2.  
     else
        pres@xyLineThicknessF   = 1.5
        res@mpGeophysicalLineThicknessF = 1.  
     end if
     pres@gsnCenterString = names(ee)
     if (spectra_mvf.eq.False) then
        pspec(ee) = gsn_csm_xy(wks2,sdof@frq,splt1,pres) 
     
        if (isfilepresent2("obs_ts").and.ee.ge.1.and.spectra_mvf_obs.eq.False) then
           pres@xyLineColors        = (/"gray70","black","black","black"/)
           pres@gsnCenterString = ""
           pspec_obs(ee) = gsn_csm_xy(wks2,sdof_obs@frq,sdof_obs@spcx,pres)
           overlay(pspec(ee),pspec_obs(ee)) 
        end if 
        delete([/sdof,splt1/])
     end if 
     
     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnRightString = ""
     xyres@gsnLeftString = ""
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"
     xyres@gsnXYBarChart = False
     xyres@gsnAboveYRefLineColor = 185
     xyres@gsnBelowYRefLineColor = 35
     xyres@xyLineThicknessF = 0.1
     xyres@xyLineColor = "gray70"
;     xyres@xyLineColors = (/   xyres@gsnAboveYRefLineColor,   xyres@gsnBelowYRefLineColor/)
     xyres@tiYAxisString = ""
     if (nsim.le.5) then
        xyres@tmXBLabelFontHeightF = 0.0125
        xyres@tmYLLabelFontHeightF = 0.0125
        xyres@gsnStringFontHeightF = 0.017     
     else
        xyres@tmXBLabelFontHeightF = 0.018
        xyres@tmYLLabelFontHeightF = 0.018
        xyres@gsnStringFontHeightF = 0.024
     end if

     xyres@vpXF = 0.05
     xyres@vpHeightF = 0.15
     if (SCALE_TIMESERIES.eq."True") then
        xyres@vpWidthF = 0.9*((nyr(ee)*1.)/nyr_max)
     else
        xyres@vpWidthF = 0.9
     end if    
     xyres@gsnCenterString = ""
     xyres@trXMinF = syear(ee)-.5
     xyres@trXMaxF = eyear(ee)+1.5
     
     xyres2 =    xyres
     delete(xyres2@gsnXYBarChart)
     delete(xyres2@gsnAboveYRefLineColor)
     delete(xyres2@gsnBelowYRefLineColor)
;     delete(xyres2@xyLineColors)
     xyres2@xyLineColor = "black"
     if (wks_type.eq."png") then
        xyres2@xyLineThicknessF = 3.5
     else
        xyres2@xyLineThicknessF = 2.5
     end if
     
     xyres@gsnCenterString = names(ee)
     xyplot(ee) = gsn_csm_xy(wks3,fspan(syear(ee),eyear(ee)+.91667,dimsizes(finarr&time)),finarr(0,:),xyres)   ; use standardized timeseries  
     xyplot2(ee) = gsn_csm_xy(wks3,fspan(syear(ee),eyear(ee)+.91667,dimsizes(finarr&time)),finarr(1,:),xyres2)
     overlay(xyplot(ee),xyplot2(ee))
     delete([/val1,val2,finarr/])
  end do

  if (isfilepresent2("obs_ts")) then    ; for pattern correlation table  
     clat = cos(0.01745329*patcor&lat)
;     finpaco = "AMO (Monthly)     "    ; Must be 18 characters long
;     finrms  = finpaco
     finpr   = "AMO (Monthly)     "    ; Must be 18 characters long
     line3   = "                  "    ; Must be 18 characters long
     line4   = line3
     do hh = 1,nsim-1
        dimY = dimsizes(tochar(names(hh)))
        nchar = dimY
        nchar = where(nchar.le.10,10,nchar)
        if (dimY.lt.10) then
           ntb = ""
           do ii = 0,10-dimY-1
              ntb = ntb+" "
           end do
           ntb = ntb+names(hh)
        else
           ntb = names(hh)
        end if
        
        ntc = ""
        do ii = 0,nchar-1
           ntc = ntc+"-"
        end do
        format2 = "%"+(nchar-5+1)+".2f"
        format3 = "%4.2f"
        line3 = line3+" "+ntb   
        line4 = line4+" "+ntc 
        if (all(ismissing(patcor(hh,:,:)))) then
           finpr = finpr+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr = finpr+sprintf(format2,(pattern_cor(patcor(0,:,:),patcor(hh,:,:),clat,0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(NewCosWeight(patcor(0,:,:))),ndtooned(NewCosWeight(patcor(hh,:,:))))))
        end if
     end do
     if (dimsizes(tochar(line4)).ge.8190) then   ; system or fortran compiler limit
        print("Metrics table warning: Not creating metrics table as size of comparison results in a invalid ascii row size.")   
     else
        if (.not.isfilepresent2(getenv("OUTDIR")+"metrics_orig.txt")) then
           header = (/"","Pattern Correlations/RMS Differences   Observations vs. Model(s)",""/)
           write_table(getenv("OUTDIR")+"metrics_orig.txt","w",[/header/],"%s")
           write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/line3/],"%s")
           write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/line4/],"%s")
           delete(header)
        end if
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr/],"%s")
     end if
     delete([/finpr,line3,line4,format2,format3,nchar,ntc,clat,patcor,ntb,dimY/])
  end if
  
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPanelBottom = 0.05
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.55
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  if (nsim.le.5) then
     panres@txFontHeightF = 0.024
  else
     panres@txFontHeightF = 0.016
  end if
  
  panres@txString = "AMO (Monthly)"
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel(wks,map,(/nrow,ncol/),panres)
  delete(wks)

  delete(panres@gsnPanelLabelBar)
  panres@txString = "AMO (Monthly)"
  gsn_panel(wks2,pspec,(/nrow,ncol/),panres)  
  delete(wks2)
  
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres@gsnPanelScalePlotIndex = tt(0)
     delete(tt)
  end if
  panres@txString = "AMO (Monthly)"
  if (nsim.le.12) then
     lp = (/nsim,1/)
  else
     lp = (/nrow,ncol/)   ;(/nsim/2+1,nsim/8+1/)  
  end if
  gsn_panel(wks3,xyplot,lp,panres)  
  delete(wks3)
  delete([/map,pspec,syear,eyear,nyr,nyr_max,SCALE_TIMESERIES,lp/])
  print("Finished: amo.ncl")
end
    
