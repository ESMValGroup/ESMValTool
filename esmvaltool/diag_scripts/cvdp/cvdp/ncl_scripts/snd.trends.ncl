; Calculates snow depth global trends
;
; Variables used: snd
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: snd.trends.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")
  
  nsim = numAsciiRow("namelist_byvar/namelist_snowdp")
  na = asciiread("namelist_byvar/namelist_snowdp",(/nsim/),"string")
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
  wks_trends_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.trends.mon")

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
  cmap = gsn_retrieve_colormap(wks_trends_djf)

  map_djf = new(nsim,"graphic")  
  map_mam = new(nsim,"graphic")  
  map_jja = new(nsim,"graphic")  
  map_son = new(nsim,"graphic")  
  map_ann = new(nsim,"graphic")  
  map_mon = new(nsim,"graphic")  

  do ee = 0,nsim-1
     snd = data_read_in(paths(ee),"SNOWDP",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up    
     if (isatt(snd&lat,"_FillValue")) then   ; required in v6.2.0-beta to reset _FillValue to avoid error message
        snd&lat@_FillValue = 1.e20
        snd&lat@missing_value = snd&lat@_FillValue
     end if
     if (isatt(snd&lon,"_FillValue")) then
        snd&lon@_FillValue = 1.e20
        snd&lon@missing_value = snd&lon@_FillValue
     end if
     
     if (isatt(snd,"is_all_missing")) then
        delete(snd)
        continue
     end if 
     if (OPT_CLIMO.eq."Full") then
        snd = rmMonAnnCycTLL(snd)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = snd
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(snd&time,-1)
        if (CLIMO_SYEAR.lt.0) then
           climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
        else
           climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        end if   
        delete(temp_arr)
        snd   = calcMonAnomTLL(snd,climo) 
        delete(climo)
     end if
 
     tttt = dtrend_msg_n(ispan(0,dimsizes(snd&time)-1,1),snd,False,True,0)
     snd_trends_mon = snd(0,:,:)
     snd_trends_mon = (/ onedtond(tttt@slope, (/dimsizes(snd&lat),dimsizes(snd&lon)/) ) /)
     snd_trends_mon = snd_trends_mon*dimsizes(snd&time)
     snd_trends_mon@units = snd@units+" "+nyr(ee)+"yr~S~-1~N~"
     delete(tttt)
     
     snd_seas = runave_n_Wrap(snd,3,0,0)
     snd_seas(0,:,:) = (/ dim_avg_n(snd(:1,:,:),0) /)
     snd_seas(dimsizes(snd&time)-1,:,:) = (/ dim_avg_n(snd(dimsizes(snd&time)-2:,:,:),0) /)
     snd_ann = runave_n_Wrap(snd,12,0,0)
     delete(snd)
     
     snd_trends_seas = snd_seas(:3,:,:)
     snd_trends_seas = snd_trends_seas@_FillValue
     snd_trends_ann  = snd_trends_seas(0,:,:)
     do ff = 0,4
        if (ff.le.3) then
           tarr = snd_seas(ff*3::12,:,:)     
        end if  
        if (ff.eq.4) then
           tarr = snd_ann(5::12,:,:)
        end if
        tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
        if (ff.le.3) then
           snd_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
        end if
        if (ff.eq.4) then
           snd_trends_ann = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
        end if
        delete([/tarr,tttt/])        
     end do
     snd_trends_seas = snd_trends_seas*nyr(ee)
     snd_trends_seas@units = snd_seas@units+" "+nyr(ee)+"yr~S~-1~N~"
     snd_trends_ann = snd_trends_ann*nyr(ee)
     snd_trends_ann@units = snd_ann@units+" "+nyr(ee)+"yr~S~-1~N~"         
     delete([/snd_seas,snd_ann/])    
     
     
     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.snd.trends."+syear(ee)+"-"+eyear(ee)+".nc"
        if (.not.isfilepresent2(fn)) then
           z = addfile(fn,"c")
           z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
           z@notes = "Data from "+names(ee)+" from "+syear(ee)+"-"+eyear(ee)
           if (OPT_CLIMO.eq."Full") then
              z@climatology = syear(ee)+"-"+eyear(ee)+" climatology removed prior to all calculations (other than means)"
           else
              if (CLIMO_SYEAR.lt.0) then
                 z@climatology = (eyear(ee)+CLIMO_SYEAR)+"-"+(eyear(ee)+CLIMO_EYEAR)+" climatology removed prior to all calculations (other than means)"
              else
                 z@climatology = CLIMO_SYEAR+"-"+CLIMO_EYEAR+" climatology removed prior to all calculations (other than means)"
              end if
           end if
           z@Conventions = "CF-1.6"
        else
           z = addfile(fn,"w")
        end if
        snd_seas = (/ snd_trends_seas /)
        snd_seas!1 = "LAT"
        snd_seas&LAT = snd_trends_seas&lat
        snd_seas!2 = "LON"
        snd_seas&LON = snd_trends_seas&lon
        copy_VarAtts(snd_trends_seas,snd_seas)
        
        snd_ann = (/ snd_trends_ann /)
        snd_ann!0 = "LAT"
        snd_ann&LAT = snd_trends_ann&lat
        snd_ann!1 = "LON"
        snd_ann&LON = snd_trends_ann&lon
        copy_VarAtts(snd_trends_ann,snd_ann)
        
        snd_mon = (/ snd_trends_mon /)
        snd_mon!0 = "LAT"
        snd_mon&LAT = snd_trends_mon&lat
        snd_mon!1 = "LON"
        snd_mon&LON = snd_trends_mon&lon
        copy_VarAtts(snd_trends_mon,snd_mon)
        
        z->snd_trends_djf     = set_varAtts(snd_seas(0,:,:),"snd linear trends (DJF)","","")
        z->snd_trends_mam     = set_varAtts(snd_seas(1,:,:),"snd linear trends (MAM)","","")
        z->snd_trends_jja     = set_varAtts(snd_seas(2,:,:),"snd linear trends (JJA)","","")
        z->snd_trends_son     = set_varAtts(snd_seas(3,:,:),"snd linear trends (SON)","","")
        z->snd_trends_ann     = set_varAtts(snd_ann,"snd linear trends (annual)","","")
        z->snd_trends_mon     = set_varAtts(snd_mon,"snd linear trends (monthly)","","")
        delete(z)
        delete([/snd_seas,snd_ann,snd_mon/])
     end if
     
     snd_trends_seas = where(abs(snd_trends_seas).le..005,snd_trends_seas@_FillValue,snd_trends_seas)  ; .005m = arbitrary # to white out
     snd_trends_ann = where(abs(snd_trends_ann).le..005,snd_trends_ann@_FillValue,snd_trends_ann)      ; areas w/very very small trends..
     snd_trends_mon = where(abs(snd_trends_mon).le..005,snd_trends_mon@_FillValue,snd_trends_mon)
;========================================================================
;     cmap = read_colormap_file("ncl_default")
     
     
     res = True
     res@mpProjection = "WinkelTripel"
     res@mpGeophysicalLineColor = "gray42"
     if (wks_type.eq."png") then
        res@mpGeophysicalLineThicknessF = 2.  
     else
        res@mpGeophysicalLineThicknessF = 1.  
     end if     
     res@mpPerimOn    = False
     res@mpGridLatSpacingF =  90            ; change latitude  line spacing
     res@mpGridLonSpacingF = 180.           ; change longitude line spacing
     res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
     res@mpFillOn = False
     res@mpCenterLonF = 0.
     res@mpOutlineOn = True  
     res@gsnDraw      = False
     res@gsnFrame     = False
  
     res@cnFillPalette = cmap(2::-1,:)
     res@cnFillMode = "RasterFill"
     res@cnLevelSelectionMode = "ExplicitLevels"
     if (COLORMAP.eq.0) then
        res@cnLevels = fspan(-.5,.5,21)
     end if
     if (COLORMAP.eq.1) then
        res@cnLevels = fspan(-.8,.8,17)
     end if
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
 
     res@gsnRightString = snd_trends_seas@units
     res@gsnCenterString = names(ee)
     map_djf(ee) = gsn_csm_contour_map(wks_trends_djf,snd_trends_seas(0,:,:),res)
     map_mam(ee) = gsn_csm_contour_map(wks_trends_mam,snd_trends_seas(1,:,:),res)
     map_jja(ee) = gsn_csm_contour_map(wks_trends_jja,snd_trends_seas(2,:,:),res)
     map_son(ee) = gsn_csm_contour_map(wks_trends_son,snd_trends_seas(3,:,:),res)
     map_ann(ee) = gsn_csm_contour_map(wks_trends_ann,snd_trends_ann,res)
     map_mon(ee) = gsn_csm_contour_map(wks_trends_mon,snd_trends_mon,res)
   
     delete([/snd_trends_seas,snd_trends_ann,snd_trends_mon/])
     delete(res)  
  end do  
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.65
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  panres@lbLabelFontHeightF = 0.013
  if (nsim.le.4) then
     if (nsim.eq.1) then
        panres@txFontHeightF = 0.022
        panres@gsnPanelBottom = 0.50
     else
        panres@txFontHeightF = 0.0145
        panres@gsnPanelBottom = 0.50
     end if
  else
     panres@txFontHeightF = 0.016
     panres@gsnPanelBottom = 0.05
  end if
  panres@lbLabelStride = 1
  
  panres@txString = "SND Trends (DJF)"
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel2(wks_trends_djf,map_djf,(/nrow,ncol/),panres)
  delete(wks_trends_djf)
  
  panres@txString = "SND Trends (MAM)"
  gsn_panel2(wks_trends_mam,map_mam,(/nrow,ncol/),panres)
  delete(wks_trends_mam)
  
  panres@txString = "SND Trends (JJA)"
  gsn_panel2(wks_trends_jja,map_jja,(/nrow,ncol/),panres)
  delete(wks_trends_jja)
  
  panres@txString = "SND Trends (SON)"
  gsn_panel2(wks_trends_son,map_son,(/nrow,ncol/),panres)
  delete(wks_trends_son)
  
  panres@txString = "SND Trends (Annual)"
  gsn_panel2(wks_trends_ann,map_ann,(/nrow,ncol/),panres)
  delete(wks_trends_ann)
  
  panres@txString = "SND Trends (Monthly)"
  gsn_panel2(wks_trends_mon,map_mon,(/nrow,ncol/),panres)
  delete(wks_trends_mon)  
  delete([/nrow,ncol,map_djf,map_mam,map_jja,map_son,map_ann,map_mon,panres,cmap/])
  print("Finished: snd.trends.ncl")
end
