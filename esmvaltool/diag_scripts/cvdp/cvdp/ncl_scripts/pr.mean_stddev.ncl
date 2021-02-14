; Calculates precipitation global means, zonal means, and standard deviations
;
; Variables used: pr
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: pr.mean_stddev.ncl")
  
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
  
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_stddev_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.stddev.djf")
  wks_stddev_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.stddev.mam")
  wks_stddev_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.stddev.jja")
  wks_stddev_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.stddev.son")
  wks_stddev_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.stddev.ann")
  wks_mean = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.mean")
  wks_za_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.za.djf")
  wks_za_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.za.mam")
  wks_za_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.za.jja")
  wks_za_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.za.son")
  wks_za_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"pr.za.ann")
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_stddev_djf,"precip3_16lev") 
     gsn_define_colormap(wks_stddev_mam,"precip3_16lev")  
     gsn_define_colormap(wks_stddev_jja,"precip3_16lev") 
     gsn_define_colormap(wks_stddev_son,"precip3_16lev")  
     gsn_define_colormap(wks_stddev_ann,"precip3_16lev")   
     gsn_define_colormap(wks_mean,"precip3_16lev") 
     gsn_define_colormap(wks_za_djf,"cb_9step")  
     gsn_define_colormap(wks_za_mam,"cb_9step")  
     gsn_define_colormap(wks_za_jja,"cb_9step")  
     gsn_define_colormap(wks_za_son,"cb_9step")  
     gsn_define_colormap(wks_za_ann,"cb_9step")  
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_stddev_djf,"cb_rainbow") 
     gsn_define_colormap(wks_stddev_mam,"cb_rainbow")  
     gsn_define_colormap(wks_stddev_jja,"cb_rainbow") 
     gsn_define_colormap(wks_stddev_son,"cb_rainbow")  
     gsn_define_colormap(wks_stddev_ann,"cb_rainbow")  
     gsn_define_colormap(wks_mean,"BlueDarkRed18")      
     gsn_define_colormap(wks_za_djf,"cb_9step")  
     gsn_define_colormap(wks_za_mam,"cb_9step")  
     gsn_define_colormap(wks_za_jja,"cb_9step")  
     gsn_define_colormap(wks_za_son,"cb_9step")  
     gsn_define_colormap(wks_za_ann,"cb_9step")  
  end if

  plot_mean_djf = new(nsim,"graphic")  
  plot_mean_mam = new(nsim,"graphic")  
  plot_mean_jja = new(nsim,"graphic")  
  plot_mean_son = new(nsim,"graphic")   
  plot_mean_ann = new(nsim,"graphic")  
  plot_stddev_djf = new(nsim,"graphic")  
  plot_stddev_mam = new(nsim,"graphic")  
  plot_stddev_jja = new(nsim,"graphic")  
  plot_stddev_son = new(nsim,"graphic")   
  plot_stddev_ann = new(nsim,"graphic")   
  
  plot_za_djf = new(nsim,"graphic")  
  plot_za_mam = new(nsim,"graphic")  
  plot_za_jja = new(nsim,"graphic")  
  plot_za_son = new(nsim,"graphic")   
  plot_za_ann = new(nsim,"graphic")  
  
  if (isfilepresent2("obs_pr")) then
     c1 = 1
  else
     c1 = 76
  end if  
;  color = (/c1,2,6,11,5,3,7,15,23,31,39,47,55,63,71,79,c1,2,6,11,5,3,7,15,23,31,39,47,55,63,71,79,c1,2,6,11,5,3,7,15,23,31,39,47,55,63,71,79,c1,2/)   
;  dash   = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3/)   

  if (nsim.le.15) then
     color = (/c1,2,6,11,5,3,7,15,23,31,39,47,55,63,71,79/)   
     dash   = (/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)   
  else
     zt = (nsim/16)+1
     color = new((/zt*16/),integer)
     dash = color
     eind = 0
     do dd = 0,zt-1
        color(eind:eind+15) = (/c1,2,6,11,5,3,7,15,23,31,39,47,55,63,71,79/)  
        if (dd.le.16) then 
           dash(eind:eind+15)  = dd
        else
           dash(eind:eind+15)  = mod(dd,16)
        end if
        eind = eind+16
     end do 
     delete([/zt,eind/])
  end if

  do ee = 0,nsim-1
     ppt = data_read_in(paths(ee),"PRECT",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up   
     if (isatt(ppt,"is_all_missing")) then
        delete(ppt)
        continue
     end if 
     do ff = 0,1     
        pptT = ppt
        if (ff.eq.1) then
           if (OPT_CLIMO.eq."Full") then
              pptT = rmMonAnnCycTLL(pptT)
           else
              check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = pptT
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(pptT&time,-1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if    
              delete(temp_arr)
              pptT   = calcMonAnomTLL(pptT,climo) 
              delete(climo)
           end if
        end if     
        ppt_seas = runave_n_Wrap(pptT,3,0,0)
        ppt_seas(0,:,:) = (/ dim_avg_n(pptT(:1,:,:),0) /)
        ppt_seas(dimsizes(pptT&time)-1,:,:) = (/ dim_avg_n(pptT(dimsizes(pptT&time)-2:,:,:),0) /)
        ppt_ann = runave_n_Wrap(pptT,12,0,0)
        delete(pptT)
        
        if (ff.eq.0) then
           ppt_mean_djf = dim_avg_n_Wrap(ppt_seas(0::12,:,:),0)
           ppt_mean_mam = dim_avg_n_Wrap(ppt_seas(3::12,:,:),0)
           ppt_mean_jja = dim_avg_n_Wrap(ppt_seas(6::12,:,:),0)
           ppt_mean_son = dim_avg_n_Wrap(ppt_seas(9::12,:,:),0)
           ppt_mean_ann = dim_avg_n_Wrap(ppt_ann(5::12,:,:),0)
  
           ppt_zamean_djf = dim_avg_n_Wrap(ppt_mean_djf,1)
           ppt_zamean_mam = dim_avg_n_Wrap(ppt_mean_mam,1)
           ppt_zamean_jja = dim_avg_n_Wrap(ppt_mean_jja,1)
           ppt_zamean_son = dim_avg_n_Wrap(ppt_mean_son,1)
           ppt_zamean_ann = dim_avg_n_Wrap(ppt_mean_ann,1)
        end if
        if (ff.eq.1) then
           ppt_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),ppt_seas(0::12,:,:),False,False,0),0)
           ppt_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),ppt_seas(3::12,:,:),False,False,0),0)
           ppt_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),ppt_seas(6::12,:,:),False,False,0),0)
           ppt_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),ppt_seas(9::12,:,:),False,False,0),0)
           ppt_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),ppt_ann(5::12,:,:),False,False,0),0)
       end if
       delete([/ppt_seas,ppt_ann/])
     end do
     delete(ppt)
     copy_VarMeta(ppt_mean_djf,ppt_sd_djf)
     copy_VarMeta(ppt_mean_mam,ppt_sd_mam)
     copy_VarMeta(ppt_mean_jja,ppt_sd_jja)
     copy_VarMeta(ppt_mean_son,ppt_sd_son)
     copy_VarMeta(ppt_mean_ann,ppt_sd_ann)     
     
     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.pr.mean_stddev."+syear(ee)+"-"+eyear(ee)+".nc"
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

        z->pr_spatialmean_djf     = set_varAtts(ppt_mean_djf,"pr mean (DJF)","","")
        z->pr_spatialmean_mam     = set_varAtts(ppt_mean_mam,"pr mean (MAM)","","")
        z->pr_spatialmean_jja     = set_varAtts(ppt_mean_jja,"pr mean (JJA)","","")
        z->pr_spatialmean_son     = set_varAtts(ppt_mean_son,"pr mean (SON)","","")
        z->pr_spatialmean_ann     = set_varAtts(ppt_mean_ann,"pr mean (annual)","","")
        
        z->pr_spatialstddev_djf     = set_varAtts(ppt_sd_djf,"pr standard deviation (DJF)","","")
        z->pr_spatialstddev_mam     = set_varAtts(ppt_sd_mam,"pr standard deviation (MAM)","","")
        z->pr_spatialstddev_jja     = set_varAtts(ppt_sd_jja,"pr standard deviation (JJA)","","")
        z->pr_spatialstddev_son     = set_varAtts(ppt_sd_son,"pr standard deviation (SON)","","")
        z->pr_spatialstddev_ann     = set_varAtts(ppt_sd_ann,"pr standard deviation (annual)","","")
        delete(z)
        delete([/modname,fn/])
     end if
;==========================================================================================
     res = True
     res@mpProjection = "WinkelTripel"
     res@mpGeophysicalLineColor = "gray42"
     res@mpPerimOn    = False
     res@mpGridLatSpacingF =  90            ; change latitude  line spacing
     res@mpGridLonSpacingF = 180.           ; change longitude line spacing
     res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
     res@mpFillOn = False
     res@mpCenterLonF = 210.
     res@mpOutlineOn = True  
     if (wks_type.eq."png") then
        res@mpGeophysicalLineThicknessF = 2.  
     else
        res@mpGeophysicalLineThicknessF = 1.  
     end if
     res@gsnDraw      = False
     res@gsnFrame     = False
  
     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     res@cnLevelSelectionMode = "ExplicitLevels"
     

     res@gsnLeftStringOrthogonalPosF = -0.05
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.05
     res@gsnRightStringParallelPosF = 0.96
     res@gsnRightString = ""
     res@gsnLeftString = ""
     res@gsnLeftStringFontHeightF = 0.014
     res@gsnCenterStringFontHeightF = 0.018
     res@gsnRightStringFontHeightF = 0.014

     sres = res
     sres@cnLevels  = (/0.5,1,2,3,4,5,6,7,8,9,10,12,14,16,18/)
     res@cnLevels = (/.2,.4,.6,1.0,1.5,2.0,2.5,3.5/)
     if (COLORMAP.eq.0) then
        res@cnFillColors = (/2,4,6,8,10,12,14,16,18/)
     end if
     if (COLORMAP.eq.1) then
        res@cnFillColors = (/35,47,63,79,95,111,124,155,175/)
     end if


     
     if (isfilepresent2("obs_prect").and.ee.eq.0) then    ; for pattern correlation table
        patcor = new((/nsim,dimsizes(ppt_sd_ann&lat),dimsizes(ppt_sd_ann&lon)/),typeof(ppt_sd_ann))
        patcor!1 = "lat"
        patcor&lat = ppt_sd_ann&lat
        patcor!2 = "lon"
        patcor&lon = ppt_sd_ann&lon
        patcor(ee,:,:) = (/ ppt_sd_ann /)
     end if
     if (isfilepresent2("obs_prect").and.ee.ge.1.and.isvar("patcor")) then
        patcor(ee,:,:) = (/ totype(linint2(ppt_sd_ann&lon,ppt_sd_ann&lat,ppt_sd_ann,True,patcor&lon,patcor&lat,0),typeof(patcor)) /)
     end if
     
     res@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     res@gsnRightString = ppt_mean_djf@units
     res@gsnCenterString = names(ee)     
     plot_stddev_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,ppt_sd_djf,res)
     plot_stddev_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,ppt_sd_mam,res)
     plot_stddev_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,ppt_sd_jja,res)
     plot_stddev_son(ee) = gsn_csm_contour_map(wks_stddev_son,ppt_sd_son,res)
     plot_stddev_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,ppt_sd_ann,res)
     
     sres@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     sres@gsnRightString = ppt_mean_djf@units
     sres@gsnCenterString = names(ee)
     plot_mean_djf(ee) = gsn_csm_contour_map(wks_mean,ppt_mean_djf,sres)
     plot_mean_mam(ee) = gsn_csm_contour_map(wks_mean,ppt_mean_mam,sres)
     plot_mean_jja(ee) = gsn_csm_contour_map(wks_mean,ppt_mean_jja,sres)
     plot_mean_son(ee) = gsn_csm_contour_map(wks_mean,ppt_mean_son,sres)
     plot_mean_ann(ee) = gsn_csm_contour_map(wks_mean,ppt_mean_ann,sres)
     delete([/ppt_sd_djf,ppt_sd_mam,ppt_sd_jja,ppt_sd_son,ppt_sd_ann,ppt_mean_djf,ppt_mean_mam,ppt_mean_jja,ppt_mean_son,ppt_mean_ann,res,sres/])
    
     zres = True
     zres@vpYF = 0.8
     zres@vpXF = 0.14
     zres@vpWidthF = 0.55
     zres@vpHeightF = 0.55
     zres@trYMinF = 0.
     zres@trYMaxF = 11.0
     zres@gsnDraw = False
     zres@gsnFrame = False
 
     zres@tmXTLabelFontHeightF = 0.018
     zres@tmXBLabelFontHeightF = 0.018
     zres@tmYLLabelFontHeightF = 0.018
     zres@tiMainString = ""
     zres@txFontHeightF = 0.015
     zres@xyLineLabelFontHeightF = 0.016
     zres@tiXAxisFontHeightF = 0.019
     zres@tiYAxisFontHeightF = 0.019
     zres@tiMainFontHeightF = 0.03
  
     zres@pmLegendDisplayMode    = "Never"
     zres@tiYAxisString = "mm day~S~-1~N~"
  
     zres@xyLineColor   =  "black"         
     zres@xyDashPattern = 0
     if (wks_type.eq."png") then
        zres@xyLineThicknessF = 3.5
        if (isfilepresent2("obs_prect").and.ee.eq.0) then
           zres@xyLineThicknessF   = 7.
        end if
     else
        zres@xyLineThicknessF = 2.
        if (isfilepresent2("obs_prect").and.ee.eq.0) then
           zres@xyLineThicknessF   = 4.
        end if
     end if
     
     zres@xyDashPattern = dash(ee)   ;dash(mod(ee,50))
     zres@xyLineColor = color(ee)    ;color(mod(ee,50))
     zres@tiMainFont = "helvetica"
     
     polyres = True
     polyres@gsLineColor = color(mod(ee,50))
     polyres@gsLineThicknessF = zres@xyLineThicknessF
     polyres@gsLineDashPattern = dash(mod(ee,50))
     
     txres = True
     if (nsim.le.15) then
        txres@txFontHeightF = 0.012
        yeval = .02
     end if
     if (nsim.ge.16.and.nsim.le.45) then
        txres@txFontHeightF = 0.009
        yeval = .0175
     end if
     if (nsim.ge.46.and.nsim.le.72) then
        txres@txFontHeightF = 0.006
        yeval = .011
     end if
     if (nsim.ge.73.and.nsim.le.106) then
        txres@txFontHeightF = 0.004
        yeval = .0075
     end if
     if (nsim.ge.107.and.nsim.le.228) then
        txres@txFontHeightF = 0.002
        yeval = .0035
     end if
     if (nsim.ge.229) then
        txres@txFontHeightF = 0.001
        yeval = .002
     end if
          
     txres@txJust = "CenterLeft"
     
     zres@tiMainString = "PR Zonal Average (DJF)"
     zres@gsnRightString = "mm/day"

     plot_za_djf(ee) = gsn_csm_xy(wks_za_djf,ppt_zamean_djf&lat,ppt_zamean_djf,zres) 
     if (ee.ne.0) then
        overlay(plot_za_djf(0),plot_za_djf(ee))
     end if     
     gsn_text_ndc(wks_za_djf,names(ee),0.765,0.8-(ee*yeval),txres)
     gsn_polyline_ndc(wks_za_djf,(/0.72,.75/),(/0.8-(ee*yeval),0.8-(ee*yeval)/),polyres)
     
     zres@tiMainString = "PR Zonal Average (MAM)"
     plot_za_mam(ee) = gsn_csm_xy(wks_za_mam,ppt_zamean_mam&lat,ppt_zamean_mam,zres) 
     if (ee.ne.0) then
        overlay(plot_za_mam(0),plot_za_mam(ee))
     end if     
     gsn_text_ndc(wks_za_mam,names(ee),0.765,0.8-(ee*yeval),txres)
     gsn_polyline_ndc(wks_za_mam,(/0.72,.75/),(/0.8-(ee*yeval),0.8-(ee*yeval)/),polyres)
     
     zres@tiMainString = "PR Zonal Average (JJA)"
     plot_za_jja(ee) = gsn_csm_xy(wks_za_jja,ppt_zamean_jja&lat,ppt_zamean_jja,zres) 
     if (ee.ne.0) then
        overlay(plot_za_jja(0),plot_za_jja(ee))
     end if     
     gsn_text_ndc(wks_za_jja,names(ee),0.765,0.8-(ee*yeval),txres)
     gsn_polyline_ndc(wks_za_jja,(/0.72,.75/),(/0.8-(ee*yeval),0.8-(ee*yeval)/),polyres)
     
     zres@tiMainString = "PR Zonal Average (SON)"
     plot_za_son(ee) = gsn_csm_xy(wks_za_son,ppt_zamean_son&lat,ppt_zamean_son,zres) 
     if (ee.ne.0) then
        overlay(plot_za_son(0),plot_za_son(ee))
     end if     
     gsn_text_ndc(wks_za_son,names(ee),0.765,0.8-(ee*yeval),txres)
     gsn_polyline_ndc(wks_za_son,(/0.72,.75/),(/0.8-(ee*yeval),0.8-(ee*yeval)/),polyres)
     
     zres@tiMainString = "PR Zonal Average (ANN)"
     plot_za_ann(ee) = gsn_csm_xy(wks_za_ann,ppt_zamean_ann&lat,ppt_zamean_ann,zres) 
     if (ee.ne.0) then
        overlay(plot_za_ann(0),plot_za_ann(ee))
     end if     
     gsn_text_ndc(wks_za_ann,names(ee),0.765,0.8-(ee*yeval),txres)
     gsn_polyline_ndc(wks_za_ann,(/0.72,.75/),(/0.8-(ee*yeval),0.8-(ee*yeval)/),polyres)
     delete([/zres,polyres,txres,ppt_zamean_djf,ppt_zamean_mam,ppt_zamean_jja,ppt_zamean_son,ppt_zamean_ann/])
  end do    
  
   if (isvar("patcor")) then    ; for pattern correlation table  
     clat = cos(0.01745329*patcor&lat)
     finpr   = "pr Std Dev (Ann)  "    ; Must be 18 characters long
     line3   = "                  "    ; Must be 18 characters long
     line4   = line3
     header = (/"","Pattern Correlations/RMS Differences   Observations vs. Model(s)",""/)
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
           finpr = finpr+sprintf(format2,(pattern_cor(patcor(0,:,:),patcor(hh,:,:),clat,0)))+"/"+sprintf(format3,(wgt_arearmse(patcor(0,:,:),patcor(hh,:,:),clat,1.0,0)))
        end if
     end do
     if (dimsizes(tochar(line4)).ge.8190) then   ; system or fortran compiler limit
        print("Metrics table warning: Not creating metrics table as size of comparison results in a invalid ascii row size.")   
     else          
        write_table(getenv("OUTDIR")+"metrics.pr.mean_stddev.txt","w",[/header/],"%s")
        write_table(getenv("OUTDIR")+"metrics.pr.mean_stddev.txt","a",[/line3/],"%s")
        write_table(getenv("OUTDIR")+"metrics.pr.mean_stddev.txt","a",[/line4/],"%s")
        write_table(getenv("OUTDIR")+"metrics.pr.mean_stddev.txt","a",[/finpr/],"%s")
     end if
     delete([/finpr,line3,line4,format2,format3,nchar,ntc,clat,patcor,dimY,ntb,header/])
  end if
     
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.65
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
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
  panres@lbLabelFontHeightF = 0.013
  panres@lbLabelStride = 1
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  
  panres@txString = "PR Standard Deviations (DJF)"
  gsn_panel2(wks_stddev_djf,plot_stddev_djf,(/nrow,ncol/),panres)
  delete(wks_stddev_djf)
  
  panres@txString = "PR Standard Deviations (MAM)"
  gsn_panel2(wks_stddev_mam,plot_stddev_mam,(/nrow,ncol/),panres)
  delete(wks_stddev_mam)
  
  panres@txString = "PR Standard Deviations (JJA)"
  gsn_panel2(wks_stddev_jja,plot_stddev_jja,(/nrow,ncol/),panres)
  delete(wks_stddev_jja)
  
  panres@txString = "PR Standard Deviations (SON)"
  gsn_panel2(wks_stddev_son,plot_stddev_son,(/nrow,ncol/),panres)
  delete(wks_stddev_son)
  
  panres@txString = "PR Standard Deviations (Annual)"
  gsn_panel2(wks_stddev_ann,plot_stddev_ann,(/nrow,ncol/),panres)
  delete(wks_stddev_ann)
  
  panres@txString = "PR Means (DJF)"
  gsn_panel2(wks_mean,plot_mean_djf,(/nrow,ncol/),panres)
  
  panres@txString = "PR Means (MAM)"
  gsn_panel2(wks_mean,plot_mean_mam,(/nrow,ncol/),panres)
  
  panres@txString = "PR Means (JJA)"
  gsn_panel2(wks_mean,plot_mean_jja,(/nrow,ncol/),panres)
  
  panres@txString = "PR Means (SON)"
  gsn_panel2(wks_mean,plot_mean_son,(/nrow,ncol/),panres)
  
  panres@txString = "PR Means (Annual)"
  gsn_panel2(wks_mean,plot_mean_ann,(/nrow,ncol/),panres)
  delete(wks_mean)
  delete(panres)
  
  draw(plot_za_djf(0))
  frame(wks_za_djf)
  delete(wks_za_djf)
  
  draw(plot_za_mam(0))
  frame(wks_za_mam)
  delete(wks_za_mam)
  
  draw(plot_za_jja(0))
  frame(wks_za_jja)
  delete(wks_za_jja)
  
  draw(plot_za_son(0))
  frame(wks_za_son)
  delete(wks_za_son)
 
  draw(plot_za_ann(0))
  frame(wks_za_ann)
  delete(wks_za_ann)
;--------------------------------------------------------------------------------------------------------------------------------------------  
  OUTDIR = getenv("OUTDIR") 
  if (wks_type.eq."png") then  
     system("mv "+OUTDIR+"pr.mean.000001.png "+OUTDIR+"pr.mean.djf.png") 
     system("mv "+OUTDIR+"pr.mean.000002.png "+OUTDIR+"pr.mean.mam.png") 
     system("mv "+OUTDIR+"pr.mean.000003.png "+OUTDIR+"pr.mean.jja.png") 
     system("mv "+OUTDIR+"pr.mean.000004.png "+OUTDIR+"pr.mean.son.png") 
     system("mv "+OUTDIR+"pr.mean.000005.png "+OUTDIR+"pr.mean.ann.png")
  else
     system("psplit "+OUTDIR+"pr.mean.ps "+OUTDIR+"pr_m")
     system("mv "+OUTDIR+"pr_m0001.ps "+OUTDIR+"pr.mean.djf.ps") 
     system("mv "+OUTDIR+"pr_m0002.ps "+OUTDIR+"pr.mean.mam.ps") 
     system("mv "+OUTDIR+"pr_m0003.ps "+OUTDIR+"pr.mean.jja.ps") 
     system("mv "+OUTDIR+"pr_m0004.ps "+OUTDIR+"pr.mean.son.ps") 
     system("mv "+OUTDIR+"pr_m0005.ps "+OUTDIR+"pr.mean.ann.ps") 
     system("rm "+OUTDIR+"pr.mean.ps")
  end if
  print("Finished: pr.mean_stddev.ncl")
end
