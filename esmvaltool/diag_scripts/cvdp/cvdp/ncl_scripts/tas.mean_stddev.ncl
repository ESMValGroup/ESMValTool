; Calculates 2m air temperature global means and standard deviations
;
; Variables used: tas
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: tas.mean_stddev.ncl")
  
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
  
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_stddev_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.stddev.djf")
  wks_stddev_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.stddev.mam")
  wks_stddev_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.stddev.jja")
  wks_stddev_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.stddev.son")
  wks_stddev_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.stddev.ann")
  wks_mean_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.mean.djf")
  wks_mean_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.mean.mam")
  wks_mean_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.mean.jja")
  wks_mean_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.mean.son")
  wks_mean_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.mean.ann")
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_stddev_djf,"precip3_16lev") 
     gsn_define_colormap(wks_stddev_mam,"precip3_16lev")  
     gsn_define_colormap(wks_stddev_jja,"precip3_16lev") 
     gsn_define_colormap(wks_stddev_son,"precip3_16lev")  
     gsn_define_colormap(wks_stddev_ann,"precip3_16lev")
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
     tas = data_read_in(paths(ee),"TREFHT",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up   
     if (isatt(tas,"is_all_missing")) then
        delete(tas)
        continue
     end if 
     do ff = 0,1     
        tasT = tas
        if (ff.eq.1) then
           if (OPT_CLIMO.eq."Full") then
              tasT = rmMonAnnCycTLL(tasT)
           else
              check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = tasT
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(tasT&time,-1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if
              delete(temp_arr)
              tasT   = calcMonAnomTLL(tasT,climo) 
              delete(climo)
           end if
        end if        
        tas_seas = runave_n_Wrap(tasT,3,0,0)
        tas_seas(0,:,:) = (/ dim_avg_n(tasT(:1,:,:),0) /)
        tas_seas(dimsizes(tasT&time)-1,:,:) = (/ dim_avg_n(tasT(dimsizes(tasT&time)-2:,:,:),0) /)
        tas_ann = runave_n_Wrap(tasT,12,0,0)
        delete(tasT)
        
        if (ff.eq.0) then
           tas_mean_djf = dim_avg_n_Wrap(tas_seas(0::12,:,:),0)
           tas_mean_mam = dim_avg_n_Wrap(tas_seas(3::12,:,:),0)
           tas_mean_jja = dim_avg_n_Wrap(tas_seas(6::12,:,:),0)
           tas_mean_son = dim_avg_n_Wrap(tas_seas(9::12,:,:),0)
           tas_mean_ann = dim_avg_n_Wrap(tas_ann(5::12,:,:),0)
        end if
        if (ff.eq.1) then     
           tas_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),tas_seas(0::12,:,:),False,False,0),0)
           tas_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),tas_seas(3::12,:,:),False,False,0),0)
           tas_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),tas_seas(6::12,:,:),False,False,0),0)
           tas_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),tas_seas(9::12,:,:),False,False,0),0)
           tas_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),tas_ann(5::12,:,:),False,False,0),0)
        end if
        delete([/tas_seas,tas_ann/])
     end do
     delete(tas)
     copy_VarMeta(tas_mean_djf,tas_sd_djf)
     copy_VarMeta(tas_mean_mam,tas_sd_mam)
     copy_VarMeta(tas_mean_jja,tas_sd_jja)
     copy_VarMeta(tas_mean_son,tas_sd_son)
     copy_VarMeta(tas_mean_ann,tas_sd_ann)
     
     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.tas.mean_stddev."+syear(ee)+"-"+eyear(ee)+".nc"
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

        z->tas_spatialmean_djf     = set_varAtts(tas_mean_djf,"tas mean (DJF)","","")
        z->tas_spatialmean_mam     = set_varAtts(tas_mean_mam,"tas mean (MAM)","","")
        z->tas_spatialmean_jja     = set_varAtts(tas_mean_jja,"tas mean (JJA)","","")
        z->tas_spatialmean_son     = set_varAtts(tas_mean_son,"tas mean (SON)","","")
        z->tas_spatialmean_ann     = set_varAtts(tas_mean_ann,"tas mean (annual)","","")
        
        z->tas_spatialstddev_djf     = set_varAtts(tas_sd_djf,"tas standard deviation (DJF)","","")
        z->tas_spatialstddev_mam     = set_varAtts(tas_sd_mam,"tas standard deviation (MAM)","","")
        z->tas_spatialstddev_jja     = set_varAtts(tas_sd_jja,"tas standard deviation (JJA)","","")
        z->tas_spatialstddev_son     = set_varAtts(tas_sd_son,"tas standard deviation (SON)","","")
        z->tas_spatialstddev_ann     = set_varAtts(tas_sd_ann,"tas standard deviation (annual)","","")
        delete(z)
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
     
     res@cnLevels = fspan(.4,3.2,8)    
     if (COLORMAP.eq.0) then
        res@cnFillColors = (/2,4,6,8,10,12,14,16,18/)
        sres@cnLevels = ispan(-40,40,2)
     end if
     if (COLORMAP.eq.1) then
        res@cnFillColors = (/35,47,63,79,95,111,124,155,175/)
        sres@cnLevels = ispan(-20,40,4)
     end if
     
     res@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     res@gsnRightString = tas_sd_djf@units
     res@gsnCenterString = names(ee)         
     plot_stddev_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,tas_sd_djf,res)
     plot_stddev_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,tas_sd_mam,res)
     plot_stddev_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,tas_sd_jja,res)
     plot_stddev_son(ee) = gsn_csm_contour_map(wks_stddev_son,tas_sd_son,res)
     plot_stddev_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,tas_sd_ann,res)
     
     sres@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     sres@gsnRightString = tas_mean_djf@units
     sres@gsnCenterString = names(ee)
     plot_mean_djf(ee) = gsn_csm_contour_map(wks_mean_djf,tas_mean_djf,sres)
     plot_mean_mam(ee) = gsn_csm_contour_map(wks_mean_mam,tas_mean_mam,sres)
     plot_mean_jja(ee) = gsn_csm_contour_map(wks_mean_jja,tas_mean_jja,sres)
     plot_mean_son(ee) = gsn_csm_contour_map(wks_mean_son,tas_mean_son,sres)
     plot_mean_ann(ee) = gsn_csm_contour_map(wks_mean_ann,tas_mean_ann,sres)
     delete([/tas_sd_djf,tas_sd_mam,tas_sd_jja,tas_sd_son,tas_sd_ann,tas_mean_djf,tas_mean_mam,tas_mean_jja,tas_mean_son,tas_mean_ann,res,sres/])
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
  
  panres@txString = "TAS Standard Deviations (DJF)"
  gsn_panel2(wks_stddev_djf,plot_stddev_djf,(/nrow,ncol/),panres)
  delete(wks_stddev_djf)
  
  panres@txString = "TAS Standard Deviations (MAM)"
  gsn_panel2(wks_stddev_mam,plot_stddev_mam,(/nrow,ncol/),panres)
  delete(wks_stddev_mam)
  
  panres@txString = "TAS Standard Deviations (JJA)"
  gsn_panel2(wks_stddev_jja,plot_stddev_jja,(/nrow,ncol/),panres)
  delete(wks_stddev_jja)
  
  panres@txString = "TAS Standard Deviations (SON)"
  gsn_panel2(wks_stddev_son,plot_stddev_son,(/nrow,ncol/),panres)
  delete(wks_stddev_son)
  
  panres@txString = "TAS Standard Deviations (Annual)"
  gsn_panel2(wks_stddev_ann,plot_stddev_ann,(/nrow,ncol/),panres)
  delete(wks_stddev_ann)
  
  panres@txString = "TAS Means (DJF)"
  gsn_panel2(wks_mean_djf,plot_mean_djf,(/nrow,ncol/),panres)
  delete(wks_mean_djf)
  
  panres@txString = "TAS Means (MAM)"
  gsn_panel2(wks_mean_mam,plot_mean_mam,(/nrow,ncol/),panres)
  delete(wks_mean_mam)
  
  panres@txString = "TAS Means (JJA)"
  gsn_panel2(wks_mean_jja,plot_mean_jja,(/nrow,ncol/),panres)
  delete(wks_mean_jja)
  
  panres@txString = "TAS Means (SON)"
  gsn_panel2(wks_mean_son,plot_mean_son,(/nrow,ncol/),panres)
  delete(wks_mean_son)
  
  panres@txString = "TAS Means (Annual)"
  gsn_panel2(wks_mean_ann,plot_mean_ann,(/nrow,ncol/),panres)
  delete(wks_mean_ann)
  delete(panres)
  print("Finished: tas.mean_stddev.ncl")
end
