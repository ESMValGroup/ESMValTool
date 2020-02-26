; Calculates 2m air temperature global means and standard deviations
;
; Variables used: snd
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: snd.mean_stddev.ncl")
  
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
  
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_stddev_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.stddev.djf")
  wks_stddev_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.stddev.mam")
  wks_stddev_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.stddev.jja")
  wks_stddev_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.stddev.son")
  wks_stddev_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.stddev.ann")
  wks_mean_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.mean.djf")
  wks_mean_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.mean.mam")
  wks_mean_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.mean.jja")
  wks_mean_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.mean.son")
  wks_mean_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"snd.mean.ann")
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_stddev_djf,"rainbow+white") 
     gsn_define_colormap(wks_stddev_mam,"rainbow+white")  
     gsn_define_colormap(wks_stddev_jja,"rainbow+white") 
     gsn_define_colormap(wks_stddev_son,"rainbow+white")  
     gsn_define_colormap(wks_stddev_ann,"rainbow+white")
     gsn_define_colormap(wks_mean_djf,"ncl_default") 
     gsn_define_colormap(wks_mean_mam,"ncl_default")  
     gsn_define_colormap(wks_mean_jja,"ncl_default") 
     gsn_define_colormap(wks_mean_son,"ncl_default")  
     gsn_define_colormap(wks_mean_ann,"ncl_default")
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_stddev_djf,"cb_rainbow") 
     gsn_define_colormap(wks_stddev_mam,"cb_rainbow")  
     gsn_define_colormap(wks_stddev_jja,"cb_rainbow") 
     gsn_define_colormap(wks_stddev_son,"cb_rainbow")  
     gsn_define_colormap(wks_stddev_ann,"cb_rainbow")
     gsn_define_colormap(wks_mean_djf,"BlueDarkRed18") 
     gsn_define_colormap(wks_mean_mam,"BlueDarkRed18")  
     gsn_define_colormap(wks_mean_jja,"BlueDarkRed18") 
     gsn_define_colormap(wks_mean_son,"BlueDarkRed18")  
     gsn_define_colormap(wks_mean_ann,"BlueDarkRed18")     
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
     do ff = 0,1     
        sndT = snd
        if (ff.eq.1) then
           if (OPT_CLIMO.eq."Full") then
              sndT = rmMonAnnCycTLL(sndT)
           else
              check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = sndT
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(sndT&time,-1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if   
              delete(temp_arr)
              sndT   = calcMonAnomTLL(sndT,climo) 
              delete(climo)
           end if
        end if         
        snd_seas = runave_n_Wrap(sndT,3,0,0)
        snd_seas(0,:,:) = (/ dim_avg_n(sndT(:1,:,:),0) /)
        snd_seas(dimsizes(sndT&time)-1,:,:) = (/ dim_avg_n(sndT(dimsizes(sndT&time)-2:,:,:),0) /)
        snd_ann = runave_n_Wrap(sndT,12,0,0)
        delete(sndT)
        
        if (ff.eq.0) then
           snd_mean_djf = dim_avg_n_Wrap(snd_seas(0::12,:,:),0)
           snd_mean_mam = dim_avg_n_Wrap(snd_seas(3::12,:,:),0)
           snd_mean_jja = dim_avg_n_Wrap(snd_seas(6::12,:,:),0)
           snd_mean_son = dim_avg_n_Wrap(snd_seas(9::12,:,:),0)
           snd_mean_ann = dim_avg_n_Wrap(snd_ann(5::12,:,:),0)
        end if
        if (ff.eq.1) then     
           snd_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),snd_seas(0::12,:,:),False,False,0),0)
           snd_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),snd_seas(3::12,:,:),False,False,0),0)
           snd_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),snd_seas(6::12,:,:),False,False,0),0)
           snd_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),snd_seas(9::12,:,:),False,False,0),0)
           snd_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),snd_ann(5::12,:,:),False,False,0),0)
        end if
        delete([/snd_seas,snd_ann/])
     end do
     delete(snd)
     copy_VarMeta(snd_mean_djf,snd_sd_djf)
     copy_VarMeta(snd_mean_mam,snd_sd_mam)
     copy_VarMeta(snd_mean_jja,snd_sd_jja)
     copy_VarMeta(snd_mean_son,snd_sd_son)
     copy_VarMeta(snd_mean_ann,snd_sd_ann)
     
     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.snd.mean_stddev."+syear(ee)+"-"+eyear(ee)+".nc"
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

        mean_djf = (/ snd_mean_djf /)
        mean_djf!0 = "LAT"
        mean_djf&LAT = snd_mean_djf&lat
        mean_djf!1 = "LON"
        mean_djf&LON = snd_mean_djf&lon
        copy_VarAtts(snd_mean_djf,mean_djf)
        mean_mam = (/ snd_mean_mam /)
        copy_VarMeta(mean_djf,mean_mam)
        mean_jja = (/ snd_mean_jja /)
        copy_VarMeta(mean_djf,mean_jja)
        mean_son = (/ snd_mean_son /)
        copy_VarMeta(mean_djf,mean_son)
        mean_ann = (/ snd_mean_ann /)
        copy_VarMeta(mean_djf,mean_ann)
        
        sd_djf = (/ snd_sd_djf /)
        sd_djf!0 = "LAT"
        sd_djf&LAT = snd_sd_djf&lat
        sd_djf!1 = "LON"
        sd_djf&LON = snd_sd_djf&lon
        copy_VarAtts(snd_sd_djf,sd_djf)
        sd_mam = (/ snd_sd_mam /)
        copy_VarMeta(sd_djf,sd_mam)
        sd_jja = (/ snd_sd_jja /)
        copy_VarMeta(sd_djf,sd_jja)
        sd_son = (/ snd_sd_son /)
        copy_VarMeta(sd_djf,sd_son)
        sd_ann = (/ snd_sd_ann /)
        copy_VarMeta(sd_djf,sd_ann)

        z->snd_spatialmean_djf     = set_varAtts(mean_djf,"snd mean (DJF)","","")
        z->snd_spatialmean_mam     = set_varAtts(mean_mam,"snd mean (MAM)","","")
        z->snd_spatialmean_jja     = set_varAtts(mean_jja,"snd mean (JJA)","","")
        z->snd_spatialmean_son     = set_varAtts(mean_son,"snd mean (SON)","","")
        z->snd_spatialmean_ann     = set_varAtts(mean_ann,"snd mean (annual)","","")
        
        z->snd_spatialstddev_djf     = set_varAtts(sd_djf,"snd standard deviation (DJF)","","")
        z->snd_spatialstddev_mam     = set_varAtts(sd_mam,"snd standard deviation (MAM)","","")
        z->snd_spatialstddev_jja     = set_varAtts(sd_jja,"snd standard deviation (JJA)","","")
        z->snd_spatialstddev_son     = set_varAtts(sd_son,"snd standard deviation (SON)","","")
        z->snd_spatialstddev_ann     = set_varAtts(sd_ann,"snd standard deviation (annual)","","")
        delete(z)
        delete(modname)
        delete([/mean_djf,mean_mam,mean_jja,mean_son,mean_ann,sd_djf,sd_mam,sd_jja,sd_son,sd_ann/])
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
     res@cnFillMode = "RasterFill"
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
     
     res@cnLevels = fspan(.05,.45,9)    
     if (COLORMAP.eq.0) then
        res@cnFillColors = (/0,54,80,95,125,175,185,195,205,236/)
        sres@cnLevels = fspan(0.05,1.5,30)
        sres@cnFillColors = ispan(8,248,8)
        sres@cnFillColors(0) = 0
     end if
     if (COLORMAP.eq.1) then
        res@cnFillColors = (/0,35,47,63,79,95,111,124,155,175/)
        sres@cnLevels = fspan(0.05,1.45,15)
        sres@cnFillColors = (/0,4,5,6,7,8,9,10, 11,12,13,14,15,16,17,18/)
     end if
     
     res@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     res@gsnRightString = snd_sd_djf@units
     res@gsnCenterString = names(ee)         

     plot_stddev_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,snd_sd_djf,res)
     plot_stddev_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,snd_sd_mam,res)
     plot_stddev_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,snd_sd_jja,res)
     plot_stddev_son(ee) = gsn_csm_contour_map(wks_stddev_son,snd_sd_son,res)
     plot_stddev_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,snd_sd_ann,res)
     
     sres@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     sres@gsnRightString = snd_mean_djf@units
     sres@gsnCenterString = names(ee)
     plot_mean_djf(ee) = gsn_csm_contour_map(wks_mean_djf,snd_mean_djf,sres)
     plot_mean_mam(ee) = gsn_csm_contour_map(wks_mean_mam,snd_mean_mam,sres)
     plot_mean_jja(ee) = gsn_csm_contour_map(wks_mean_jja,snd_mean_jja,sres)
     plot_mean_son(ee) = gsn_csm_contour_map(wks_mean_son,snd_mean_son,sres)
     plot_mean_ann(ee) = gsn_csm_contour_map(wks_mean_ann,snd_mean_ann,sres)
     delete([/snd_sd_djf,snd_sd_mam,snd_sd_jja,snd_sd_son,snd_sd_ann,snd_mean_djf,snd_mean_mam,snd_mean_jja,snd_mean_son,snd_mean_ann,res,sres/])
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
  
  panres@txString = "SND Standard Deviations (DJF)"
  gsn_panel2(wks_stddev_djf,plot_stddev_djf,(/nrow,ncol/),panres)
  delete(wks_stddev_djf)
  
  panres@txString = "SND Standard Deviations (MAM)"
  gsn_panel2(wks_stddev_mam,plot_stddev_mam,(/nrow,ncol/),panres)
  delete(wks_stddev_mam)
  
  panres@txString = "SND Standard Deviations (JJA)"
  gsn_panel2(wks_stddev_jja,plot_stddev_jja,(/nrow,ncol/),panres)
  delete(wks_stddev_jja)
  
  panres@txString = "SND Standard Deviations (SON)"
  gsn_panel2(wks_stddev_son,plot_stddev_son,(/nrow,ncol/),panres)
  delete(wks_stddev_son)
  
  panres@txString = "SND Standard Deviations (Annual)"
  gsn_panel2(wks_stddev_ann,plot_stddev_ann,(/nrow,ncol/),panres)
  delete(wks_stddev_ann)
  
  panres@txString = "SND Means (DJF)"
  gsn_panel2(wks_mean_djf,plot_mean_djf,(/nrow,ncol/),panres)
  delete(wks_mean_djf)
  
  panres@txString = "SND Means (MAM)"
  gsn_panel2(wks_mean_mam,plot_mean_mam,(/nrow,ncol/),panres)
  delete(wks_mean_mam)
  
  panres@txString = "SND Means (JJA)"
  gsn_panel2(wks_mean_jja,plot_mean_jja,(/nrow,ncol/),panres)
  delete(wks_mean_jja)
  
  panres@txString = "SND Means (SON)"
  gsn_panel2(wks_mean_son,plot_mean_son,(/nrow,ncol/),panres)
  delete(wks_mean_son)
  
  panres@txString = "SND Means (Annual)"
  gsn_panel2(wks_mean_ann,plot_mean_ann,(/nrow,ncol/),panres)
  delete(wks_mean_ann)
  delete(panres)
  print("Finished: snd.mean_stddev.ncl")
end
