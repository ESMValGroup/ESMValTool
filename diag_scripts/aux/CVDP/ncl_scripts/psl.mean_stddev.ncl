; Calculates PSL global means and standard deviations
;
; Variables used: psl
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: psl.mean_stddev.ncl")
  
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
  
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_stddev_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.stddev.djf")
  wks_stddev_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.stddev.mam")
  wks_stddev_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.stddev.jja")
  wks_stddev_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.stddev.son")
  wks_stddev_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.stddev.ann")
  wks_mean_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.mean.djf")
  wks_mean_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.mean.mam")
  wks_mean_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.mean.jja")
  wks_mean_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.mean.son")
  wks_mean_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"psl.mean.ann")

  
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
     psl = data_read_in(paths(ee),"PSL",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up   
     if (isatt(psl,"is_all_missing")) then
        delete(psl)
        continue
     end if 
     psl_seas = runave_n_Wrap(psl,3,0,0)
     psl_seas(0,:,:) = (/ dim_avg_n(psl(:1,:,:),0) /)
     psl_seas(dimsizes(psl&time)-1,:,:) = (/ dim_avg_n(psl(dimsizes(psl&time)-2:,:,:),0) /)
     psl_ann = runave_n_Wrap(psl,12,0,0)
     delete(psl)
     
     psl_mean_djf = dim_avg_n_Wrap(psl_seas(0::12,:,:),0)
     psl_mean_mam = dim_avg_n_Wrap(psl_seas(3::12,:,:),0)
     psl_mean_jja = dim_avg_n_Wrap(psl_seas(6::12,:,:),0)
     psl_mean_son = dim_avg_n_Wrap(psl_seas(9::12,:,:),0)
     psl_mean_ann = dim_avg_n_Wrap(psl_ann(5::12,:,:),0)
     
     psl_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),psl_seas(0::12,:,:),False,False,0),0)
     psl_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),psl_seas(3::12,:,:),False,False,0),0)
     psl_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),psl_seas(6::12,:,:),False,False,0),0)
     psl_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),psl_seas(9::12,:,:),False,False,0),0)
     psl_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),psl_ann(5::12,:,:),False,False,0),0)

     copy_VarMeta(psl_mean_djf,psl_sd_djf)
     copy_VarMeta(psl_mean_mam,psl_sd_mam)
     copy_VarMeta(psl_mean_jja,psl_sd_jja)
     copy_VarMeta(psl_mean_son,psl_sd_son)
     copy_VarMeta(psl_mean_ann,psl_sd_ann)

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

        z->psl_spatialmean_djf     = psl_mean_djf
        z->psl_spatialmean_mam     = psl_mean_mam
        z->psl_spatialmean_jja     = psl_mean_jja
        z->psl_spatialmean_son     = psl_mean_son
        z->psl_spatialmean_ann     = psl_mean_ann
        
        z->psl_spatialstddev_djf     = psl_sd_djf
        z->psl_spatialstddev_mam     = psl_sd_mam
        z->psl_spatialstddev_jja     = psl_sd_jja
        z->psl_spatialstddev_son     = psl_sd_son
        z->psl_spatialstddev_ann     = psl_sd_ann
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
     res@mpDataSetName = "Earth..4"
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
     
     res@cnLevels = fspan(.4,6.0,8)
     if (COLORMAP.eq.0) then
        res@cnFillColors = (/54,80,95,125,175,185,195,205,236/)
        sres@cnLevels  = ispan(972,1044,4)
     end if
     if (COLORMAP.eq.1) then
        res@cnFillColors = (/35,47,63,79,95,111,124,155,175/)
        sres@cnLevels  = ispan(980,1036,4)
     end if


     
     if (isfilepresent2("obs_psl").and.ee.eq.0) then    ; for pattern correlation table
        patcor = new((/nsim,dimsizes(psl_sd_ann&lat),dimsizes(psl_sd_ann&lon)/),typeof(psl_sd_ann))
        patcor!1 = "lat"
        patcor&lat = psl_sd_ann&lat
        patcor!2 = "lon"
        patcor&lon = psl_sd_ann&lon
        patcor(ee,:,:) = (/ psl_sd_ann /)
     end if
     if (isfilepresent2("obs_psl").and.ee.ge.1) then
        patcor(ee,:,:) = (/ linint2(psl_sd_ann&lon,psl_sd_ann&lat,psl_sd_ann,True,patcor&lon,patcor&lat,0) /)
     end if

     res@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     res@gsnRightString = psl_mean_djf@units
     res@gsnCenterString = names(ee)     
     plot_stddev_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,psl_sd_djf,res)
     plot_stddev_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,psl_sd_mam,res)
     plot_stddev_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,psl_sd_jja,res)
     plot_stddev_son(ee) = gsn_csm_contour_map(wks_stddev_son,psl_sd_son,res)
     plot_stddev_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,psl_sd_ann,res)
     
     sres@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     sres@gsnRightString = psl_mean_djf@units
     sres@gsnCenterString = names(ee)
     plot_mean_djf(ee) = gsn_csm_contour_map(wks_mean_djf,psl_mean_djf,sres)
     plot_mean_mam(ee) = gsn_csm_contour_map(wks_mean_mam,psl_mean_mam,sres)
     plot_mean_jja(ee) = gsn_csm_contour_map(wks_mean_jja,psl_mean_jja,sres)
     plot_mean_son(ee) = gsn_csm_contour_map(wks_mean_son,psl_mean_son,sres)
     plot_mean_ann(ee) = gsn_csm_contour_map(wks_mean_ann,psl_mean_ann,sres)
     delete([/psl_sd_djf,psl_sd_mam,psl_sd_jja,psl_sd_son,psl_sd_ann,psl_mean_djf,psl_mean_mam,psl_mean_jja,psl_mean_son,psl_mean_ann,res,sres/])
  end do     
  
  if (isfilepresent2("obs_psl")) then    ; for pattern correlation table  
     clat = cos(0.01745329*patcor&lat)
     finpr   = "PSL Std Dev (Ann) "    ; Must be 18 characters long
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
     end if
     write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr/],"%s")
     delete([/finpr,line3,line4,format2,format3,nchar,ntc,clat,patcor,dimY,ntb/])
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
  panres@lbLabelFontHeightF = 0.013
  panres@lbLabelStride = 1
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  
  panres@txString = "PSL Standard Deviations (DJF)"
  gsn_panel(wks_stddev_djf,plot_stddev_djf,(/nrow,ncol/),panres)
  delete(wks_stddev_djf)
  
  panres@txString = "PSL Standard Deviations (MAM)"
  gsn_panel(wks_stddev_mam,plot_stddev_mam,(/nrow,ncol/),panres)
  delete(wks_stddev_mam)
  
  panres@txString = "PSL Standard Deviations (JJA)"
  gsn_panel(wks_stddev_jja,plot_stddev_jja,(/nrow,ncol/),panres)
  delete(wks_stddev_jja)
  
  panres@txString = "PSL Standard Deviations (SON)"
  gsn_panel(wks_stddev_son,plot_stddev_son,(/nrow,ncol/),panres)
  delete(wks_stddev_son)
  
  panres@txString = "PSL Standard Deviations (Annual)"
  gsn_panel(wks_stddev_ann,plot_stddev_ann,(/nrow,ncol/),panres)
  delete(wks_stddev_ann)
  
  panres@txString = "PSL Means (DJF)"
  gsn_panel(wks_mean_djf,plot_mean_djf,(/nrow,ncol/),panres)
  delete(wks_mean_djf)
  
  panres@txString = "PSL Means (MAM)"
  gsn_panel(wks_mean_mam,plot_mean_mam,(/nrow,ncol/),panres)
  delete(wks_mean_mam)
  
  panres@txString = "PSL Means (JJA)"
  gsn_panel(wks_mean_jja,plot_mean_jja,(/nrow,ncol/),panres)
  delete(wks_mean_jja)
  
  panres@txString = "PSL Means (SON)"
  gsn_panel(wks_mean_son,plot_mean_son,(/nrow,ncol/),panres)
  delete(wks_mean_son)
  
  panres@txString = "PSL Means (Annual)"
  gsn_panel(wks_mean_ann,plot_mean_ann,(/nrow,ncol/),panres)
  delete(wks_mean_ann)
  delete(panres)
  print("Finished: psl.mean_stddev.ncl")
end
