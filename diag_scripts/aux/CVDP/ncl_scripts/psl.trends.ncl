; Calculates PSL global trends
;
; Variables used: psl
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: psl.trends.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")
  
  nsim = numAsciiRow("namelist_byvar/namelist_psl")
  na = asciiread("namelist_byvar/namelist_psl",(/nsim/),"string")
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
  wks_trends_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.trends.mon")

  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_trends_djf,"ncl_default")   
     gsn_define_colormap(wks_trends_mam,"ncl_default")  
     gsn_define_colormap(wks_trends_jja,"ncl_default") 
     gsn_define_colormap(wks_trends_son,"ncl_default") 
     gsn_define_colormap(wks_trends_ann,"ncl_default") 
     gsn_define_colormap(wks_trends_mon,"ncl_default") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_trends_djf,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_mam,"BlueDarkRed18")  
     gsn_define_colormap(wks_trends_jja,"BlueDarkRed18")  
     gsn_define_colormap(wks_trends_son,"BlueDarkRed18")  
     gsn_define_colormap(wks_trends_ann,"BlueDarkRed18")  
     gsn_define_colormap(wks_trends_mon,"BlueDarkRed18")  
  end if

  map_djf = new(nsim,"graphic")  
  map_mam = new(nsim,"graphic")  
  map_jja = new(nsim,"graphic")  
  map_son = new(nsim,"graphic")  
  map_ann = new(nsim,"graphic")  
  map_mon = new(nsim,"graphic")  

  do ee = 0,nsim-1
     psl = data_read_in(paths(ee),"PSL",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(psl,"is_all_missing")) then
        delete(psl)
        continue
     end if 
     if (OPT_CLIMO.eq."Full") then
        psl = rmMonAnnCycTLL(psl)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = psl
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(psl&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        psl   = calcMonAnomTLL(psl,climo) 
        delete(climo)
     end if
 
     tttt = dtrend_msg_n(ispan(0,dimsizes(psl&time)-1,1),psl,False,True,0)
     psl_trends_mon = psl(0,:,:)
     psl_trends_mon = (/ onedtond(tttt@slope, (/dimsizes(psl&lat),dimsizes(psl&lon)/) ) /)
     psl_trends_mon = psl_trends_mon*dimsizes(psl&time)
     psl_trends_mon@units = psl@units+" "+nyr(ee)+"yr~S~-1~N~"
     delete(tttt)
     
     psl_seas = runave_n_Wrap(psl,3,0,0)
     psl_seas(0,:,:) = (/ dim_avg_n(psl(:1,:,:),0) /)
     psl_seas(dimsizes(psl&time)-1,:,:) = (/ dim_avg_n(psl(dimsizes(psl&time)-2:,:,:),0) /)
     psl_ann = runave_n_Wrap(psl,12,0,0)
     delete(psl)
     
     psl_trends_seas = psl_seas(:3,:,:)
     psl_trends_seas = psl_trends_seas@_FillValue
     psl_trends_ann  = psl_trends_seas(0,:,:)
     do ff = 0,4
        if (ff.le.3) then
           tarr = psl_seas(ff*3::12,:,:)     
        end if  
        if (ff.eq.4) then
           tarr = psl_ann(5::12,:,:)
        end if
        tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
        if (ff.le.3) then
           psl_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
        end if
        if (ff.eq.4) then
           psl_trends_ann = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
        end if
        delete([/tarr,tttt/])        
     end do
     psl_trends_seas = psl_trends_seas*nyr(ee)
     psl_trends_seas@units = psl_seas@units+" "+nyr(ee)+"yr~S~-1~N~"
     psl_trends_ann = psl_trends_ann*nyr(ee)
     psl_trends_ann@units = psl_ann@units+" "+nyr(ee)+"yr~S~-1~N~"         
     delete([/psl_seas,psl_ann/])    
     
     
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
        z->psl_trends_djf     = psl_trends_seas(0,:,:)
        z->psl_trends_mam     = psl_trends_seas(1,:,:)
        z->psl_trends_jja     = psl_trends_seas(2,:,:)
        z->psl_trends_son     = psl_trends_seas(3,:,:)
        z->psl_trends_ann     = psl_trends_ann
        z->psl_trends_mon     = psl_trends_mon
        delete(z)
     end if
     
;========================================================================
     res = True
     res@mpProjection = "WinkelTripel"
     res@mpGeophysicalLineColor = "gray42"
     res@mpGeophysicalLineThicknessF = 2.          
     res@mpPerimOn    = False
     res@mpGridLatSpacingF =  90            ; change latitude  line spacing
     res@mpGridLonSpacingF = 180.           ; change longitude line spacing
     res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
     res@mpFillOn = False
     res@mpCenterLonF = 210.
     res@mpOutlineOn = True  
     res@mpDataSetName = "Earth..4"
     res@gsnDraw      = False
     res@gsnFrame     = False
  
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLevels = ispan(-8,8,1)
     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     res@gsnLeftStringOrthogonalPosF = -0.05
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.05
     res@gsnRightStringParallelPosF = 0.96
     res@gsnRightString = ""
     res@gsnLeftString = ""
     res@gsnLeftStringFontHeightF = 0.014
     res@gsnCenterStringFontHeightF = 0.018
     res@gsnRightStringFontHeightF = 0.014
     res@gsnLeftString = syear(ee)+"-"+eyear(ee)
 
     res@gsnRightString = psl_trends_seas@units
     res@gsnCenterString = names(ee)
     map_djf(ee) = gsn_csm_contour_map(wks_trends_djf,psl_trends_seas(0,:,:),res)
     map_mam(ee) = gsn_csm_contour_map(wks_trends_mam,psl_trends_seas(1,:,:),res)
     map_jja(ee) = gsn_csm_contour_map(wks_trends_jja,psl_trends_seas(2,:,:),res)
     map_son(ee) = gsn_csm_contour_map(wks_trends_son,psl_trends_seas(3,:,:),res)
     map_ann(ee) = gsn_csm_contour_map(wks_trends_ann,psl_trends_ann,res)
     map_mon(ee) = gsn_csm_contour_map(wks_trends_mon,psl_trends_mon,res)
   
     delete([/psl_trends_seas,psl_trends_ann,psl_trends_mon,res/])  
  end do  
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPanelBottom = 0.05
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.65
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  if (nsim.le.5) then
     panres@txFontHeightF = 0.024
  else
     panres@txFontHeightF = 0.016
  end if
  panres@lbLabelStride = 1
  
  panres@txString = "PSL Trends (DJF)"
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel(wks_trends_djf,map_djf,(/nrow,ncol/),panres)
  delete(wks_trends_djf)
  
  panres@txString = "PSL Trends (MAM)"
  gsn_panel(wks_trends_mam,map_mam,(/nrow,ncol/),panres)
  delete(wks_trends_mam)
  
  panres@txString = "PSL Trends (JJA)"
  gsn_panel(wks_trends_jja,map_jja,(/nrow,ncol/),panres)
  delete(wks_trends_jja)
  
  panres@txString = "PSL Trends (SON)"
  gsn_panel(wks_trends_son,map_son,(/nrow,ncol/),panres)
  delete(wks_trends_son)
  
  panres@txString = "PSL Trends (Annual)"
  gsn_panel(wks_trends_ann,map_ann,(/nrow,ncol/),panres)
  delete(wks_trends_ann)
  
  panres@txString = "PSL Trends (Monthly)"
  gsn_panel(wks_trends_mon,map_mon,(/nrow,ncol/),panres)
  delete(wks_trends_mon)  
  delete([/nrow,ncol,map_djf,map_mam,map_jja,map_son,map_ann,map_mon,panres/])
  print("Finished: psl.trends.ncl")
end
