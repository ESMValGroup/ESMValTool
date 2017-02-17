; Calculates precipitation global trends and timeseries
;
; Variables used: pr
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: pr.trends_timeseries.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")  
  
  nsim = numAsciiRow("namelist_byvar/namelist_prect")
  na = asciiread("namelist_byvar/namelist_prect",(/nsim/),"string")
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
  wks_trends_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.trends.mon")
    
  wks_aa_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.timeseries.djf")
  wks_aa_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.timeseries.mam")
  wks_aa_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.timeseries.jja")
  wks_aa_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.timeseries.son")
  wks_aa_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.timeseries.ann")
  wks_aa_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.timeseries.mon")
  
  wks_rt_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.runtrend.mon")

  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_trends_djf,"precip_diff_12lev")   
     gsn_define_colormap(wks_trends_mam,"precip_diff_12lev")  
     gsn_define_colormap(wks_trends_jja,"precip_diff_12lev") 
     gsn_define_colormap(wks_trends_son,"precip_diff_12lev") 
     gsn_define_colormap(wks_trends_ann,"precip_diff_12lev") 
     gsn_define_colormap(wks_trends_mon,"precip_diff_12lev") 
     gsn_define_colormap(wks_aa_djf,"ncl_default")   
     gsn_define_colormap(wks_aa_mam,"ncl_default")  
     gsn_define_colormap(wks_aa_jja,"ncl_default") 
     gsn_define_colormap(wks_aa_son,"ncl_default") 
     gsn_define_colormap(wks_aa_ann,"ncl_default") 
     gsn_define_colormap(wks_aa_mon,"ncl_default") 
     gsn_define_colormap(wks_rt_mon,"ncl_default") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_trends_djf,"BrownBlue12")     
     gsn_define_colormap(wks_trends_mam,"BrownBlue12")    
     gsn_define_colormap(wks_trends_jja,"BrownBlue12")   
     gsn_define_colormap(wks_trends_son,"BrownBlue12")   
     gsn_define_colormap(wks_trends_ann,"BrownBlue12")   
     gsn_define_colormap(wks_trends_mon,"BrownBlue12")   
     gsn_define_colormap(wks_aa_djf,"ncl_default")   
     gsn_define_colormap(wks_aa_mam,"ncl_default")  
     gsn_define_colormap(wks_aa_jja,"ncl_default") 
     gsn_define_colormap(wks_aa_son,"ncl_default") 
     gsn_define_colormap(wks_aa_ann,"ncl_default") 
     gsn_define_colormap(wks_aa_mon,"ncl_default") 
     gsn_define_colormap(wks_rt_mon,"ncl_default") 
  end if
  
  map_djf = new(nsim,"graphic")  
  map_mam = new(nsim,"graphic")  
  map_jja = new(nsim,"graphic")  
  map_son = new(nsim,"graphic")  
  map_ann = new(nsim,"graphic")  
  map_mon = new(nsim,"graphic")  
  xy_djf = new(nsim,"graphic")  
  xy_mam = new(nsim,"graphic")  
  xy_jja = new(nsim,"graphic")  
  xy_son = new(nsim,"graphic")  
  xy_ann = new(nsim,"graphic")  
  xy_mon = new(nsim,"graphic")  
  
  xy_rt_mon_8 = new(nsim,"graphic")  
  xy_rt_mon_10 = new(nsim,"graphic")  
  xy_rt_mon_12 = new(nsim,"graphic")  
  xy_rt_mon_14 = new(nsim,"graphic")  
  xy_rt_mon_16 = new(nsim,"graphic")  
  
  if (isfilepresent2("obs_prect")) then
     xy_obs_djf = new(nsim,"graphic")  
     xy_obs_mam = new(nsim,"graphic")  
     xy_obs_jja = new(nsim,"graphic")  
     xy_obs_son = new(nsim,"graphic")  
     xy_obs_ann = new(nsim,"graphic")  
     xy_obs_mon = new(nsim,"graphic")  
  end if
  do ee = 0,nsim-1
     ppt = data_read_in(paths(ee),"PRECT",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(ppt,"is_all_missing")) then
        delete(ppt)
        continue
     end if 

     if (OPT_CLIMO.eq."Full") then
        ppt = rmMonAnnCycTLL(ppt)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = ppt
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(ppt&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        ppt   = calcMonAnomTLL(ppt,climo) 
        delete(climo)
     end if

     coswgt=cos(rad*ppt&lat)
     coswgt!0 = "lat"
     coswgt&lat= ppt&lat
     
     ppt_aa_mon = wgt_areaave_Wrap(ppt,coswgt,1.0,0)
     tttt = dtrend_msg_n(ispan(0,dimsizes(ppt&time)-1,1),ppt,False,True,0)
     ppt_trends_mon = ppt(0,:,:)
     ppt_trends_mon = (/ onedtond(tttt@slope, (/dimsizes(ppt&lat),dimsizes(ppt&lon)/) ) /)
     ppt_trends_mon = ppt_trends_mon*dimsizes(ppt&time)
     ppt_trends_mon@units = ppt@units+" "+nyr(ee)+"yr~S~-1~N~"
     delete(tttt)

     ppt_seas = runave_n_Wrap(ppt,3,0,0)
     ppt_seas(0,:,:) = (/ dim_avg_n(ppt(:1,:,:),0) /)
     ppt_seas(dimsizes(ppt&time)-1,:,:) = (/ dim_avg_n(ppt(dimsizes(ppt&time)-2:,:,:),0) /)
     ppt_ann = runave_n_Wrap(ppt,12,0,0)
     delete(ppt)
     
     ppt_trends_seas = ppt_seas(:3,:,:)
     ppt_trends_seas = ppt_trends_seas@_FillValue
     ppt_trends_ann  = ppt_trends_seas(0,:,:)
     ppt_aa_seas = new((/4,nyr(ee)/),typeof(ppt_seas))
     ppt_aa_seas!1 = "time"
     ppt_aa_seas&time = ispan(syear(ee),eyear(ee),1)
     ppt_aa_seas&time@units = "YYYY"
     ppt_aa_seas&time@long_name = "time"
     ppt_aa_ann = ppt_aa_seas(0,:)
     do ff = 0,4
        if (ff.le.3) then
           tarr = ppt_seas(ff*3::12,:,:)     
        end if  
        if (ff.eq.4) then
           tarr = ppt_ann(5::12,:,:)
        end if
        tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
        if (ff.le.3) then
           ppt_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
           ppt_aa_seas(ff,:) = (/ wgt_areaave(tarr,coswgt,1.0,0) /)
        end if
        if (ff.eq.4) then
           ppt_trends_ann = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
           ppt_aa_ann = (/ wgt_areaave(tarr,coswgt,1.0,0) /)
        end if
        delete([/tarr,tttt/])        
     end do
     ppt_trends_seas = ppt_trends_seas*nyr(ee)
     ppt_trends_seas@units = ppt_seas@units+" "+nyr(ee)+"yr~S~-1~N~"
     ppt_trends_ann = ppt_trends_ann*nyr(ee)
     ppt_trends_ann@units = ppt_ann@units+" "+nyr(ee)+"yr~S~-1~N~"         
     delete([/ppt_seas,ppt_ann,coswgt/])    

     if (isfilepresent2("obs_prect").and.ee.eq.0) then
        ppt_aa_seas@syear = syear(ee)
        ppt_aa_seas@eyear = eyear(ee)
        ppt_aa_mon@syear = syear(ee)
        ppt_aa_mon@eyear = eyear(ee)
        ppt_aa_ann@syear = syear(ee)
        ppt_aa_ann@eyear = eyear(ee)
        ppt_aa_seas_obs = ppt_aa_seas
        ppt_aa_mon_obs  = ppt_aa_mon
        ppt_aa_ann_obs  = ppt_aa_ann
     end if

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
        ppt_aa_seas2 = ppt_aa_seas
        ppt_aa_seas2!1 = "TIME"
        ppt_aa_seas2&TIME = ispan(syear(ee),eyear(ee),1)
        ppt_aa_seas2&TIME@units = "YYYY"
        ppt_aa_seas2&TIME@long_name = "time"
        ppt_aa_ann2 = ppt_aa_ann
        ppt_aa_ann2!0 = "TIME"
        ppt_aa_ann2&TIME = ispan(syear(ee),eyear(ee),1)
        ppt_aa_ann2&TIME@units = "YYYY"
        ppt_aa_ann2&TIME@long_name = "time"
        z->ppt_global_avg_mon = ppt_aa_mon
        z->ppt_global_avg_djf = ppt_aa_seas2(0,:)
        z->ppt_global_avg_mam = ppt_aa_seas2(1,:)
        z->ppt_global_avg_jja = ppt_aa_seas2(2,:)
        z->ppt_global_avg_son = ppt_aa_seas2(3,:)
        z->ppt_global_avg_ann = ppt_aa_ann2
        z->ppt_trends_djf     = ppt_trends_seas(0,:,:)
        z->ppt_trends_mam     = ppt_trends_seas(1,:,:)
        z->ppt_trends_jja     = ppt_trends_seas(2,:,:)
        z->ppt_trends_son     = ppt_trends_seas(3,:,:)
        z->ppt_trends_ann     = ppt_trends_ann
        z->ppt_trends_mon     = ppt_trends_mon
        delete(z)
        delete(ppt_aa_seas2)
        delete(ppt_aa_ann2)
     end if    
;========================================================================
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
     res@mpCenterLonF = 210.
     res@mpOutlineOn = True  
     res@mpDataSetName = "Earth..4"
     res@gsnDraw      = False
     res@gsnFrame     = False
  
     res@cnLevelSelectionMode = "ExplicitLevels"
     if (COLORMAP.eq.0) then
        res@cnLevels = (/-6,-4,-2,-1,-0.5,-0.2,0,0.2,0.5,1,2,4,6/)
        res@cnFillColors = (/2,3,4,5,6,7,8,8,9,10,11,12,13,14/)
     end if
     if (COLORMAP.eq.1) then
         res@cnLevels = (/-4,-2,-1,-0.5,-0.2,0,0.2,0.5,1,2,4/)
     end if
     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     res@gsnLeftStringOrthogonalPosF = -0.05
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.05
     res@gsnRightStringParallelPosF = 0.975
     res@gsnRightString = ""
     res@gsnLeftString = ""
     res@gsnLeftStringFontHeightF = 0.014
     res@gsnCenterStringFontHeightF = 0.018
     res@gsnRightStringFontHeightF = 0.014
     res@gsnLeftString = syear(ee)+"-"+eyear(ee)
 
     res@gsnRightString = ppt_trends_seas@units
     res@gsnCenterString = names(ee)
     map_djf(ee) = gsn_csm_contour_map(wks_trends_djf,ppt_trends_seas(0,:,:),res)
     map_mam(ee) = gsn_csm_contour_map(wks_trends_mam,ppt_trends_seas(1,:,:),res)
     map_jja(ee) = gsn_csm_contour_map(wks_trends_jja,ppt_trends_seas(2,:,:),res)
     map_son(ee) = gsn_csm_contour_map(wks_trends_son,ppt_trends_seas(3,:,:),res)
     map_ann(ee) = gsn_csm_contour_map(wks_trends_ann,ppt_trends_ann,res)
     map_mon(ee) = gsn_csm_contour_map(wks_trends_mon,ppt_trends_mon,res)
     
     
     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnFrame = False
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"
     
     if (wks_type.eq."png") then
        xyres@xyLineThicknessF = 4.
     else
        xyres@xyLineThicknessF = 2.0
     end if
     if (isfilepresent2("obs_prect").and.ee.eq.0) then
        xyres@xyLineColor = "black"
     else
        xyres@xyLineColor = "royalblue"
     end if
     xyres@tiYAxisString = ""
     if (nsim.le.5) then
        xyres@tmXBLabelFontHeightF = 0.0125
        xyres@tmYLLabelFontHeightF = 0.0125
        xyres@gsnLeftStringFontHeightF = 0.017     
        xyres@gsnRightStringFontHeightF = 0.013     
     else
        xyres@tmXBLabelFontHeightF = 0.018
        xyres@tmYLLabelFontHeightF = 0.018
        xyres@gsnLeftStringFontHeightF = 0.024
        xyres@gsnRightStringFontHeightF = 0.020     
     end if
     xyres@vpXF = 0.05
     xyres@vpHeightF = 0.15
     if (SCALE_TIMESERIES.eq."True") then
        xyres@vpWidthF = 0.9*((nyr(ee)*1.)/nyr_max)
     else
        xyres@vpWidthF = 0.9
     end if
     xyres@gsnLeftString = ""     
     xyres@gsnCenterString = ""
     xyres@gsnRightString = ""
     
     xyres@trXMinF = syear(ee)-.5
     xyres@trXMaxF = eyear(ee)+0.5
     
     xyres2 = xyres
     xyres2@xyLineColor = "gray60"
     xyres2@xyCurveDrawOrder = "PreDraw"
     
     xyres@gsnLeftString = names(ee)
     tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_seas&time)-1,1),ppt_aa_seas(0,:),False,True)
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xyres@trYMinF = min((/min(ppt_aa_seas(0,:)),min(ppt_aa_seas_obs(0,:))/))-.005
        xyres@trYMaxF = max((/max(ppt_aa_seas(0,:)),max(ppt_aa_seas_obs(0,:))/))+.005
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+ppt_trends_seas@units
     xy_djf(ee)     = gsn_csm_xy(wks_aa_djf,ispan(syear(ee),eyear(ee),1),ppt_aa_seas(0,:),xyres)
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xy_obs_djf(ee) = gsn_csm_xy(wks_aa_djf,ispan(ppt_aa_seas_obs@syear,ppt_aa_seas_obs@eyear,1),ppt_aa_seas_obs(0,:),xyres2)
        overlay(xy_djf(ee),xy_obs_djf(ee))
     end if
     
     tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_seas&time)-1,1),ppt_aa_seas(1,:),False,True)  
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xyres@trYMinF = min((/min(ppt_aa_seas(1,:)),min(ppt_aa_seas_obs(1,:))/))-.005
        xyres@trYMaxF = max((/max(ppt_aa_seas(1,:)),max(ppt_aa_seas_obs(1,:))/))+.005
     end if 
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+ppt_trends_seas@units
     xy_mam(ee)     = gsn_csm_xy(wks_aa_mam,ispan(syear(ee),eyear(ee),1),ppt_aa_seas(1,:),xyres)  
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xy_obs_mam(ee) = gsn_csm_xy(wks_aa_mam,ispan(ppt_aa_seas_obs@syear,ppt_aa_seas_obs@eyear,1),ppt_aa_seas_obs(1,:),xyres2)
        overlay(xy_mam(ee),xy_obs_mam(ee))
     end if
     
     tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_seas&time)-1,1),ppt_aa_seas(2,:),False,True)   
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xyres@trYMinF = min((/min(ppt_aa_seas(2,:)),min(ppt_aa_seas_obs(2,:))/))-.005
        xyres@trYMaxF = max((/max(ppt_aa_seas(2,:)),max(ppt_aa_seas_obs(2,:))/))+.005
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+ppt_trends_seas@units
     xy_jja(ee)     = gsn_csm_xy(wks_aa_jja,ispan(syear(ee),eyear(ee),1),ppt_aa_seas(2,:),xyres)  
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xy_obs_jja(ee) = gsn_csm_xy(wks_aa_jja,ispan(ppt_aa_seas_obs@syear,ppt_aa_seas_obs@eyear,1),ppt_aa_seas_obs(2,:),xyres2)
        overlay(xy_jja(ee),xy_obs_jja(ee))
     end if
     
     tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_seas&time)-1,1),ppt_aa_seas(3,:),False,True) 
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xyres@trYMinF = min((/min(ppt_aa_seas(3,:)),min(ppt_aa_seas_obs(3,:))/))-.005
        xyres@trYMaxF = max((/max(ppt_aa_seas(3,:)),max(ppt_aa_seas_obs(3,:))/))+.005
     end if  
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+ppt_trends_seas@units
     xy_son(ee)     = gsn_csm_xy(wks_aa_son,ispan(syear(ee),eyear(ee),1),ppt_aa_seas(3,:),xyres)   
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xy_obs_son(ee) = gsn_csm_xy(wks_aa_son,ispan(ppt_aa_seas_obs@syear,ppt_aa_seas_obs@eyear,1),ppt_aa_seas_obs(3,:),xyres2)
        overlay(xy_son(ee),xy_obs_son(ee))
     end if
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_mon&time)-1,1),ppt_aa_mon,False,True)   
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xyres@trYMinF = min((/min(ppt_aa_mon),min(ppt_aa_mon_obs)/))-.005
        xyres@trYMaxF = max((/max(ppt_aa_mon),max(ppt_aa_mon_obs)/))+.005
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(ppt_aa_mon&time),2,True)+ppt_trends_mon@units
     xy_mon(ee)     = gsn_csm_xy(wks_aa_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(ppt_aa_mon)),ppt_aa_mon,xyres)  
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xy_obs_mon(ee) = gsn_csm_xy(wks_aa_mon,fspan(ppt_aa_seas_obs@syear,ppt_aa_seas_obs@eyear+.91667,dimsizes(ppt_aa_mon_obs)),ppt_aa_mon_obs,xyres2)
        overlay(xy_mon(ee),xy_obs_mon(ee))
     end if
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_ann&time)-1,1),ppt_aa_ann,False,True)   
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xyres@trYMinF = min((/min(ppt_aa_ann),min(ppt_aa_ann_obs)/))-.005
        xyres@trYMaxF = max((/max(ppt_aa_ann),max(ppt_aa_ann_obs)/))+.005
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+ppt_trends_ann@units
     xy_ann(ee)     = gsn_csm_xy(wks_aa_ann,ispan(syear(ee),eyear(ee),1),ppt_aa_ann,xyres)   
     if (isfilepresent2("obs_prect").and.ee.ge.1) then
        xy_obs_ann(ee) = gsn_csm_xy(wks_aa_ann,ispan(ppt_aa_seas_obs@syear,ppt_aa_seas_obs@eyear,1),ppt_aa_ann_obs,xyres2)
        overlay(xy_ann(ee),xy_obs_ann(ee))
        delete(xyres@trYMinF)
        delete(xyres@trYMaxF)
     end if
     
     delete([/ppt_trends_seas,ppt_trends_ann,ppt_trends_mon/])  
     delete([/ppt_aa_seas,ppt_aa_mon,ppt_aa_ann,xyres,xyres2,res,tttt/])
  end do
  if (isfilepresent2("obs_prect")) then
     delete([/ppt_aa_seas_obs,ppt_aa_mon_obs,ppt_aa_ann_obs/])
  end if 
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
  
  
  panres@txString = "PPT Trends (DJF)"
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel(wks_trends_djf,map_djf,(/nrow,ncol/),panres)
  delete(wks_trends_djf)
  
  panres@txString = "PPT Trends (MAM)"
  gsn_panel(wks_trends_mam,map_mam,(/nrow,ncol/),panres)
  delete(wks_trends_mam)

  panres@txString = "PPT Trends (JJA)"
  gsn_panel(wks_trends_jja,map_jja,(/nrow,ncol/),panres)
  delete(wks_trends_jja)

  panres@txString = "PPT Trends (SON)"
  gsn_panel(wks_trends_son,map_son,(/nrow,ncol/),panres)
  delete(wks_trends_son)

  panres@txString = "PPT Trends (Annual)"
  gsn_panel(wks_trends_ann,map_ann,(/nrow,ncol/),panres)
  delete(wks_trends_ann)

  panres@txString = "PPT Trends (Monthly)"
  gsn_panel(wks_trends_mon,map_mon,(/nrow,ncol/),panres)
  delete(wks_trends_mon)  
  
  panres2 = True
  panres2@gsnMaximize = True
  panres2@gsnPaperOrientation = "portrait"
  panres2@gsnPanelYWhiteSpacePercent = 3.0  
  if (nsim.le.5) then
     panres2@txFontHeightF = 0.024
  else
     panres2@txFontHeightF = 0.016
  end if
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres2@gsnPanelScalePlotIndex = tt(0)
     delete(tt)
  end if
  if (nsim.le.12) then
     lp = (/nsim,1/)
  else
     lp = (/nrow,ncol/)  ;(/nsim/2+1,nsim/8+1/)  
  end if
  panres2@txString = "PR Global Average (DJF)"
  gsn_panel(wks_aa_djf,xy_djf,lp,panres2)  
  delete(wks_aa_djf)
  
  panres2@txString = "PR Global Average (MAM)"
  gsn_panel(wks_aa_mam,xy_mam,lp,panres2)  
  delete(wks_aa_mam)
  
  panres2@txString = "PR Global Average (JJA)"
  gsn_panel(wks_aa_jja,xy_jja,lp,panres2)  
  delete(wks_aa_jja)
  
  panres2@txString = "PR Global Average (SON)"
  gsn_panel(wks_aa_son,xy_son,lp,panres2)  
  delete(wks_aa_son)
  
  panres2@txString = "PR Global Average (Annual)"
  gsn_panel(wks_aa_ann,xy_ann,lp,panres2)  
  delete(wks_aa_ann)
  
  panres2@txString = "PR Global Average (Monthly)"
  gsn_panel(wks_aa_mon,xy_mon,lp,panres2)  
  delete(wks_aa_mon)
  
  delete([/nrow,ncol,lp,map_djf,map_mam,map_jja,map_son,map_ann,map_mon,xy_djf,xy_mam,xy_jja,xy_son,xy_ann,xy_mon/])
  delete(panres2)
  if (isfilepresent2("obs_prect")) then
     delete([/xy_obs_djf,xy_obs_mam,xy_obs_jja,xy_obs_son,xy_obs_ann,xy_obs_mon/])
  end if
  print("Finished: pr.trends_timeseries.ncl")
end
