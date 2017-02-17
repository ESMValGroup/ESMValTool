; Calculates 2m air temperature global trends, running global trends and timeseries
;
; Variables used: TREFHT (tas)
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: tas.trends_timeseries.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")  
  
  nsim = numAsciiRow("namelist_byvar/namelist_trefht")
  na = asciiread("namelist_byvar/namelist_trefht",(/nsim/),"string")
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
  wks_trends_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.trends.mon")
    
  wks_aa_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.timeseries.djf")
  wks_aa_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.timeseries.mam")
  wks_aa_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.timeseries.jja")
  wks_aa_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.timeseries.son")
  wks_aa_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.timeseries.ann")
  wks_aa_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.timeseries.mon")
  
  wks_rt_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"air2m.runtrend.mon")

  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_trends_djf,"ncl_default")   
     gsn_define_colormap(wks_trends_mam,"ncl_default")  
     gsn_define_colormap(wks_trends_jja,"ncl_default") 
     gsn_define_colormap(wks_trends_son,"ncl_default") 
     gsn_define_colormap(wks_trends_ann,"ncl_default") 
     gsn_define_colormap(wks_trends_mon,"ncl_default") 
     gsn_define_colormap(wks_aa_djf,"ncl_default")   
     gsn_define_colormap(wks_aa_mam,"ncl_default")  
     gsn_define_colormap(wks_aa_jja,"ncl_default") 
     gsn_define_colormap(wks_aa_son,"ncl_default") 
     gsn_define_colormap(wks_aa_ann,"ncl_default") 
     gsn_define_colormap(wks_aa_mon,"ncl_default") 
     gsn_define_colormap(wks_rt_mon,"ncl_default") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_trends_djf,"BlueDarkRed18")     
     gsn_define_colormap(wks_trends_mam,"BlueDarkRed18")    
     gsn_define_colormap(wks_trends_jja,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_son,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_ann,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_mon,"BlueDarkRed18")   
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
  
  if (isfilepresent2("obs_trefht")) then
     xy_obs_djf = new(nsim,"graphic")  
     xy_obs_mam = new(nsim,"graphic")  
     xy_obs_jja = new(nsim,"graphic")  
     xy_obs_son = new(nsim,"graphic")  
     xy_obs_ann = new(nsim,"graphic")  
     xy_obs_mon = new(nsim,"graphic")  
  end if
  do ee = 0,nsim-1
     tas = data_read_in(paths(ee),"TREFHT",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(tas,"is_all_missing")) then
        delete(tas)
        continue
     end if 
     if (OPT_CLIMO.eq."Full") then
        tas = rmMonAnnCycTLL(tas)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = tas
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(tas&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        tas   = calcMonAnomTLL(tas,climo) 
        delete(climo)
     end if
  
     coswgt=cos(rad*tas&lat)
     coswgt!0 = "lat"
     coswgt&lat= tas&lat
     
     tas_aa_mon = wgt_areaave_Wrap(tas,coswgt,1.0,0)
     tttt = dtrend_msg_n(ispan(0,dimsizes(tas&time)-1,1),tas,False,True,0)
     tas_trends_mon = tas(0,:,:)
     tas_trends_mon = (/ onedtond(tttt@slope, (/dimsizes(tas&lat),dimsizes(tas&lon)/) ) /)
     tas_trends_mon = tas_trends_mon*dimsizes(tas&time)
     tas_trends_mon@units = tas@units+" "+nyr(ee)+"yr~S~-1~N~"
     delete(tttt)
     
     tas_seas = runave_n_Wrap(tas,3,0,0)
     tas_seas(0,:,:) = (/ dim_avg_n(tas(:1,:,:),0) /)
     tas_seas(dimsizes(tas&time)-1,:,:) = (/ dim_avg_n(tas(dimsizes(tas&time)-2:,:,:),0) /)
     tas_ann = runave_n_Wrap(tas,12,0,0)
     delete(tas)
     
     tas_trends_seas = tas_seas(:3,:,:)
     tas_trends_seas = tas_trends_seas@_FillValue
     tas_trends_ann  = tas_trends_seas(0,:,:)
     tas_aa_seas = new((/4,nyr(ee)/),typeof(tas_seas))
     tas_aa_seas!1 = "time"
     tas_aa_seas&time = ispan(syear(ee),eyear(ee),1)
     tas_aa_seas&time@units = "YYYY"
     tas_aa_seas&time@long_name = "time"
     tas_aa_ann = tas_aa_seas(0,:)
     do ff = 0,4
        if (ff.le.3) then
           tarr = tas_seas(ff*3::12,:,:)     
        end if  
        if (ff.eq.4) then
           tarr = tas_ann(5::12,:,:)
        end if
        tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
        if (ff.le.3) then
           tas_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
           tas_aa_seas(ff,:) = (/ wgt_areaave(tarr,coswgt,1.0,0) /)
        end if
        if (ff.eq.4) then
           tas_trends_ann = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
           tas_aa_ann = (/ wgt_areaave(tarr,coswgt,1.0,0) /)
        end if
        delete([/tarr,tttt/])        
     end do
     tas_trends_seas = tas_trends_seas*nyr(ee)
     tas_trends_seas@units = tas_seas@units+" "+nyr(ee)+"yr~S~-1~N~"
     tas_trends_ann = tas_trends_ann*nyr(ee)
     tas_trends_ann@units = tas_ann@units+" "+nyr(ee)+"yr~S~-1~N~"         
     delete([/tas_seas,tas_ann,coswgt/])
     
     if (isfilepresent2("obs_trefht").and.ee.eq.0) then
        tas_aa_seas@syear = syear(ee)
        tas_aa_seas@eyear = eyear(ee)
        tas_aa_mon@syear = syear(ee)
        tas_aa_mon@eyear = eyear(ee)
        tas_aa_ann@syear = syear(ee)
        tas_aa_ann@eyear = eyear(ee)
        tas_aa_seas_obs = tas_aa_seas
        tas_aa_mon_obs  = tas_aa_mon
        tas_aa_ann_obs  = tas_aa_ann
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
        tas_aa_seas2 = tas_aa_seas
        tas_aa_seas2!1 = "TIME"
        tas_aa_seas2&TIME = ispan(syear(ee),eyear(ee),1)
        tas_aa_seas2&TIME@units = "YYYY"
        tas_aa_seas2&TIME@long_name = "time"
        tas_aa_ann2 = tas_aa_ann
        tas_aa_ann2!0 = "TIME"
        tas_aa_ann2&TIME = ispan(syear(ee),eyear(ee),1)
        tas_aa_ann2&TIME@units = "YYYY"
        tas_aa_ann2&TIME@long_name = "time"
        z->tas_global_avg_mon = tas_aa_mon
        z->tas_global_avg_djf = tas_aa_seas2(0,:)
        z->tas_global_avg_mam = tas_aa_seas2(1,:)
        z->tas_global_avg_jja = tas_aa_seas2(2,:)
        z->tas_global_avg_son = tas_aa_seas2(3,:)
        z->tas_global_avg_ann = tas_aa_ann2
        z->tas_trends_djf     = tas_trends_seas(0,:,:)
        z->tas_trends_mam     = tas_trends_seas(1,:,:)
        z->tas_trends_jja     = tas_trends_seas(2,:,:)
        z->tas_trends_son     = tas_trends_seas(3,:,:)
        z->tas_trends_ann     = tas_trends_ann
        z->tas_trends_mon     = tas_trends_mon
        delete(z)
        delete([/tas_aa_seas2,tas_aa_ann2/])
     end if
     
     dimT = dimsizes(tas_aa_mon)      ; calculate running trends from the monthly data
     tas_rt_mon = new((/5,dimT/),typeof(tas_aa_mon))
     tas_rt_mon!1 = "time"
     tas_rt_mon&time = tas_aa_mon&time     
     rt_nyr = (/8,10,12,14,16/)
     do ff = 0,dimsizes(rt_nyr)-1
        incr = rt_nyr(ff)*12
        do gg = 0,dimT-incr-1
           tttt = dtrend_msg(ispan(0,incr-1,1),tas_aa_mon(gg:gg+incr-1),False,True)
           tas_rt_mon(ff,gg) = (/ tttt@slope*incr /)
           delete(tttt) 
        end do
     end do
     delete(dimT)      
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
        res@cnLevels = (/-8,-6,-5,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,5,6,8/)
     end if
     if (COLORMAP.eq.1) then
        res@cnLevels = (/-6,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,6/)
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
 
     res@gsnRightString = tas_trends_seas@units
     res@gsnCenterString = names(ee)
     map_djf(ee) = gsn_csm_contour_map(wks_trends_djf,tas_trends_seas(0,:,:),res)
     map_mam(ee) = gsn_csm_contour_map(wks_trends_mam,tas_trends_seas(1,:,:),res)
     map_jja(ee) = gsn_csm_contour_map(wks_trends_jja,tas_trends_seas(2,:,:),res)
     map_son(ee) = gsn_csm_contour_map(wks_trends_son,tas_trends_seas(3,:,:),res)
     map_ann(ee) = gsn_csm_contour_map(wks_trends_ann,tas_trends_ann,res)
     map_mon(ee) = gsn_csm_contour_map(wks_trends_mon,tas_trends_mon,res)

     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     

     xyres@gsnFrame = False
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"
     
     if (wks_type.eq."png") then
        xyres@xyLineThicknessF = 4.
     else
        xyres@xyLineThicknessF = 2.
     end if
     if (isfilepresent2("obs_trefht").and.ee.eq.0) then
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
     tttt = dtrend_msg(ispan(0,dimsizes(tas_aa_seas&time)-1,1),tas_aa_seas(0,:),False,True)   
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(tas_aa_seas(0,:)),min(tas_aa_seas_obs(0,:))/))-.01
        xyres@trYMaxF = max((/max(tas_aa_seas(0,:)),max(tas_aa_seas_obs(0,:))/))+.01
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+tas_trends_seas@units
     xy_djf(ee)     = gsn_csm_xy(wks_aa_djf,ispan(syear(ee),eyear(ee),1),tas_aa_seas(0,:),xyres)
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xy_obs_djf(ee) = gsn_csm_xy(wks_aa_djf,ispan(tas_aa_seas_obs@syear,tas_aa_seas_obs@eyear,1),tas_aa_seas_obs(0,:),xyres2)
        overlay(xy_djf(ee),xy_obs_djf(ee))
     end if
     
     tttt = dtrend_msg(ispan(0,dimsizes(tas_aa_seas&time)-1,1),tas_aa_seas(1,:),False,True)   
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(tas_aa_seas(1,:)),min(tas_aa_seas_obs(1,:))/))-.01
        xyres@trYMaxF = max((/max(tas_aa_seas(1,:)),max(tas_aa_seas_obs(1,:))/))+.01
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+tas_trends_seas@units
     xy_mam(ee)     = gsn_csm_xy(wks_aa_mam,ispan(syear(ee),eyear(ee),1),tas_aa_seas(1,:),xyres)  
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xy_obs_mam(ee) = gsn_csm_xy(wks_aa_mam,ispan(tas_aa_seas_obs@syear,tas_aa_seas_obs@eyear,1),tas_aa_seas_obs(1,:),xyres2)
        overlay(xy_mam(ee),xy_obs_mam(ee))
     end if
     
     tttt = dtrend_msg(ispan(0,dimsizes(tas_aa_seas&time)-1,1),tas_aa_seas(2,:),False,True)   
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(tas_aa_seas(2,:)),min(tas_aa_seas_obs(2,:))/))-.01
        xyres@trYMaxF = max((/max(tas_aa_seas(2,:)),max(tas_aa_seas_obs(2,:))/))+.01
     end if
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+tas_trends_seas@units
     xy_jja(ee)     = gsn_csm_xy(wks_aa_jja,ispan(syear(ee),eyear(ee),1),tas_aa_seas(2,:),xyres)  
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xy_obs_jja(ee) = gsn_csm_xy(wks_aa_jja,ispan(tas_aa_seas_obs@syear,tas_aa_seas_obs@eyear,1),tas_aa_seas_obs(2,:),xyres2)
        overlay(xy_jja(ee),xy_obs_jja(ee))
     end if
     
     tttt = dtrend_msg(ispan(0,dimsizes(tas_aa_seas&time)-1,1),tas_aa_seas(3,:),False,True)
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(tas_aa_seas(3,:)),min(tas_aa_seas_obs(3,:))/))-.01
        xyres@trYMaxF = max((/max(tas_aa_seas(3,:)),max(tas_aa_seas_obs(3,:))/))+.01
     end if   
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+tas_trends_seas@units
     xy_son(ee)     = gsn_csm_xy(wks_aa_son,ispan(syear(ee),eyear(ee),1),tas_aa_seas(3,:),xyres)   
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xy_obs_son(ee) = gsn_csm_xy(wks_aa_son,ispan(tas_aa_seas_obs@syear,tas_aa_seas_obs@eyear,1),tas_aa_seas_obs(3,:),xyres2)
        overlay(xy_son(ee),xy_obs_son(ee))
     end if
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(tas_aa_mon&time)-1,1),tas_aa_mon,False,True)
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(tas_aa_mon),min(tas_aa_mon_obs)/))-.01
        xyres@trYMaxF = max((/max(tas_aa_mon),max(tas_aa_mon_obs)/))+.01
     end if   
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(tas_aa_mon&time),2,True)+tas_trends_mon@units
     xy_mon(ee)     = gsn_csm_xy(wks_aa_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tas_aa_mon)),tas_aa_mon,xyres)  
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xy_obs_mon(ee) = gsn_csm_xy(wks_aa_mon,fspan(tas_aa_seas_obs@syear,tas_aa_seas_obs@eyear+.91667,dimsizes(tas_aa_mon_obs)),tas_aa_mon_obs,xyres2)
        overlay(xy_mon(ee),xy_obs_mon(ee))
     end if
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(tas_aa_ann&time)-1,1),tas_aa_ann,False,True)   
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(tas_aa_ann),min(tas_aa_ann_obs)/))-.01
        xyres@trYMaxF = max((/max(tas_aa_ann),max(tas_aa_ann_obs)/))+.01
     end if   
     xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+tas_trends_ann@units
     xy_ann(ee)     = gsn_csm_xy(wks_aa_ann,ispan(syear(ee),eyear(ee),1),tas_aa_ann,xyres)   
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xy_obs_ann(ee) = gsn_csm_xy(wks_aa_ann,ispan(tas_aa_seas_obs@syear,tas_aa_seas_obs@eyear,1),tas_aa_ann_obs,xyres2)
        overlay(xy_ann(ee),xy_obs_ann(ee))
        delete(xyres@trYMinF)
        delete(xyres@trYMaxF)
     end if
     
     xyres@gsnRightString = ""     ; this section could be made to be more generic, in case rt_nyr gets changed. 1 graphical array would be needed..
     xy_rt_mon_8(ee)  = gsn_csm_xy(wks_rt_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tas_aa_mon&time)),tas_rt_mon(0,:),xyres)
     xy_rt_mon_10(ee) = gsn_csm_xy(wks_rt_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tas_aa_mon&time)),tas_rt_mon(1,:),xyres)
     xy_rt_mon_12(ee) = gsn_csm_xy(wks_rt_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tas_aa_mon&time)),tas_rt_mon(2,:),xyres)
     xy_rt_mon_14(ee) = gsn_csm_xy(wks_rt_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tas_aa_mon&time)),tas_rt_mon(3,:),xyres)
     xy_rt_mon_16(ee) = gsn_csm_xy(wks_rt_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tas_aa_mon&time)),tas_rt_mon(4,:),xyres)
         
     delete([/tas_trends_seas,tas_trends_ann,tas_trends_mon/])  
     delete([/tas_aa_seas,tas_aa_mon,tas_aa_ann,xyres,xyres2,res,tttt,tas_rt_mon/])
  end do
  if (isfilepresent2("obs_trefht")) then
     delete([/tas_aa_seas_obs,tas_aa_mon_obs,tas_aa_ann_obs/])
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
  
  panres@txString = "TAS Trends (DJF)"
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel(wks_trends_djf,map_djf,(/nrow,ncol/),panres)
  delete(wks_trends_djf)
  
  panres@txString = "TAS Trends (MAM)"
  gsn_panel(wks_trends_mam,map_mam,(/nrow,ncol/),panres)
  delete(wks_trends_mam)
  
  panres@txString = "TAS Trends (JJA)"
  gsn_panel(wks_trends_jja,map_jja,(/nrow,ncol/),panres)
  delete(wks_trends_jja)
  
  panres@txString = "TAS Trends (SON)"
  gsn_panel(wks_trends_son,map_son,(/nrow,ncol/),panres)
  delete(wks_trends_son)
  
  panres@txString = "TAS Trends (Annual)"
  gsn_panel(wks_trends_ann,map_ann,(/nrow,ncol/),panres)
  delete(wks_trends_ann)
  
  panres@txString = "TAS Trends (Monthly)"
  gsn_panel(wks_trends_mon,map_mon,(/nrow,ncol/),panres)
  delete(wks_trends_mon)
  delete(panres)
  
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
     lp = (/nrow,ncol/)   ;(/nsim/2+1,nsim/8+1/)  
  end if
  panres2@txString = "TAS Global Average (DJF)"
  gsn_panel(wks_aa_djf,xy_djf,lp,panres2)  
  delete(wks_aa_djf)
  
  panres2@txString = "TAS Global Average (MAM)"
  gsn_panel(wks_aa_mam,xy_mam,lp,panres2)  
  delete(wks_aa_mam)
  
  panres2@txString = "TAS Global Average (JJA)"
  gsn_panel(wks_aa_jja,xy_jja,lp,panres2)  
  delete(wks_aa_jja)
  
  panres2@txString = "TAS Global Average (SON)"
  gsn_panel(wks_aa_son,xy_son,lp,panres2)  
  delete(wks_aa_son)
  
  panres2@txString = "TAS Global Average (Annual)"
  gsn_panel(wks_aa_ann,xy_ann,lp,panres2)  
  delete(wks_aa_ann)
  
  panres2@txString = "TAS Global Average (Monthly)"
  gsn_panel(wks_aa_mon,xy_mon,lp,panres2)  
  delete(wks_aa_mon)
  
  panres2@txString = "TAS Running 8yr Trend (Monthly)"
  gsn_panel(wks_rt_mon,xy_rt_mon_8,lp,panres2)  
  
  panres2@txString = "TAS Running 10yr Trend (Monthly)"
  gsn_panel(wks_rt_mon,xy_rt_mon_10,lp,panres2)  
  
  panres2@txString = "TAS Running 12yr Trend (Monthly)"
  gsn_panel(wks_rt_mon,xy_rt_mon_12,lp,panres2)  
  
  panres2@txString = "TAS Running 14yr Trend (Monthly)"
  gsn_panel(wks_rt_mon,xy_rt_mon_14,lp,panres2)  
  
  panres2@txString = "TAS Running 16yr Trend (Monthly)"
  gsn_panel(wks_rt_mon,xy_rt_mon_16,lp,panres2)  
  delete(wks_rt_mon)
  
  delete([/nrow,ncol,lp,map_djf,map_mam,map_jja,map_son,map_ann,map_mon,xy_djf,xy_mam,xy_jja,xy_son,xy_ann,xy_mon/])
  delete([/xy_rt_mon_8,xy_rt_mon_10,xy_rt_mon_12,xy_rt_mon_14,xy_rt_mon_16/])
  delete(panres2)
  if (isfilepresent2("obs_trefht")) then
     delete([/xy_obs_djf,xy_obs_mam,xy_obs_jja,xy_obs_son,xy_obs_ann,xy_obs_mon/])
  end if
  OUTDIR = getenv("OUTDIR")
  if (wks_type.eq."png") then  
     system("mv "+OUTDIR+"air2m.runtrend.mon.000001.png "+OUTDIR+"air2m."+rt_nyr(0)+"yr_runtrend.mon.png")
     system("mv "+OUTDIR+"air2m.runtrend.mon.000002.png "+OUTDIR+"air2m."+rt_nyr(1)+"yr_runtrend.mon.png")
     system("mv "+OUTDIR+"air2m.runtrend.mon.000003.png "+OUTDIR+"air2m."+rt_nyr(2)+"yr_runtrend.mon.png")
     system("mv "+OUTDIR+"air2m.runtrend.mon.000004.png "+OUTDIR+"air2m."+rt_nyr(3)+"yr_runtrend.mon.png")
     system("mv "+OUTDIR+"air2m.runtrend.mon.000005.png "+OUTDIR+"air2m."+rt_nyr(4)+"yr_runtrend.mon.png")
  else
     system("psplit "+OUTDIR+"air2m.runtrend.mon.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"air2m."+rt_nyr(0)+"yr_runtrend.mon.ps")
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"air2m."+rt_nyr(1)+"yr_runtrend.mon.ps")
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"air2m."+rt_nyr(2)+"yr_runtrend.mon.ps")
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"air2m."+rt_nyr(3)+"yr_runtrend.mon.ps")
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"air2m."+rt_nyr(4)+"yr_runtrend.mon.ps")
  end if
  delete(OUTDIR)
  print("Finished: tas.trends_timeseries.ncl")
end
