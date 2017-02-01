; Calculates SST global means and standard deviations
;
; Variables used: TS
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: sst.mean_stddev.ncl")
  
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
  
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_stddev_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.stddev.djf")
  wks_stddev_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.stddev.mam")
  wks_stddev_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.stddev.jja")
  wks_stddev_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.stddev.son")
  wks_stddev_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.stddev.ann")
  wks_mean_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.mean.djf")
  wks_mean_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.mean.mam")
  wks_mean_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.mean.jja")
  wks_mean_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.mean.son")
  wks_mean_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst.mean.ann")
  
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
     sst = data_read_in(paths(ee),"TS",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(sst,"is_all_missing")) then
        delete(sst)
        continue
     end if 
     if (isfilepresent2("obs_ts").and.ee.eq.0) then   ; mask out land in array going into EOF functions
     else
        sst = where(sst.le.-1.8,-1.8,sst)
        d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
        basemap = d->LSMASK
        lsm = landsea_mask(basemap,sst&lat,sst&lon)
        sst = mask(sst,conform(sst,lsm,(/1,2/)).ge.1,False)
;        printStatInfo(sst)
        delete([/lsm,basemap/])
        delete(d)
     end if    
     sst_seas = runave_n_Wrap(sst,3,0,0)
     sst_seas(0,:,:) = (/ dim_avg_n(sst(:1,:,:),0) /)
     sst_seas(dimsizes(sst&time)-1,:,:) = (/ dim_avg_n(sst(dimsizes(sst&time)-2:,:,:),0) /)
     sst_ann = runave_n_Wrap(sst,12,0,0)
     delete(sst)
     
     sst_mean_djf = dim_avg_n_Wrap(sst_seas(0::12,:,:),0)
     sst_mean_mam = dim_avg_n_Wrap(sst_seas(3::12,:,:),0)
     sst_mean_jja = dim_avg_n_Wrap(sst_seas(6::12,:,:),0)
     sst_mean_son = dim_avg_n_Wrap(sst_seas(9::12,:,:),0)
     sst_mean_ann = dim_avg_n_Wrap(sst_ann(5::12,:,:),0)
     
     sst_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),sst_seas(0::12,:,:),False,False,0),0)
     sst_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),sst_seas(3::12,:,:),False,False,0),0)
     sst_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),sst_seas(6::12,:,:),False,False,0),0)
     sst_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),sst_seas(9::12,:,:),False,False,0),0)
     sst_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),sst_ann(5::12,:,:),False,False,0),0)

     copy_VarMeta(sst_mean_djf,sst_sd_djf)
     copy_VarMeta(sst_mean_mam,sst_sd_mam)
     copy_VarMeta(sst_mean_jja,sst_sd_jja)
     copy_VarMeta(sst_mean_son,sst_sd_son)
     copy_VarMeta(sst_mean_ann,sst_sd_ann)

     delete([/sst_seas,sst_ann/])
     
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

        z->sst_spatialmean_djf     = sst_mean_djf
        z->sst_spatialmean_mam     = sst_mean_mam
        z->sst_spatialmean_jja     = sst_mean_jja
        z->sst_spatialmean_son     = sst_mean_son
        z->sst_spatialmean_ann     = sst_mean_ann
        
        z->sst_spatialstddev_djf     = sst_sd_djf
        z->sst_spatialstddev_mam     = sst_sd_mam
        z->sst_spatialstddev_jja     = sst_sd_jja
        z->sst_spatialstddev_son     = sst_sd_son
        z->sst_spatialstddev_ann     = sst_sd_ann
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
     
     res@cnLevels = fspan(.2,1.6,8)
     if (COLORMAP.eq.0) then
        res@cnFillColors = (/54,80,95,125,175,185,195,205,236/)
        res@mpLandFillColor = "gray75"
        sres@cnLevels = ispan(0,36,2)
     end if
     if (COLORMAP.eq.1) then
        res@cnFillColors = (/35,47,63,79,95,111,124,155,175/)
        res@mpLandFillColor = "gray30"
        sres@cnLevels = ispan(4,32,2)
     end if

     if (isfilepresent2("obs_ts").and.ee.eq.0) then    ; for pattern correlation table
        patcor = new((/nsim,dimsizes(sst_sd_ann&lat),dimsizes(sst_sd_ann&lon)/),typeof(sst_sd_ann))
        patcor!1 = "lat"
        patcor&lat = sst_sd_ann&lat
        patcor!2 = "lon"
        patcor&lon = sst_sd_ann&lon
        patcor(ee,:,:) = (/ sst_sd_ann /)
     end if
     if (isfilepresent2("obs_ts").and.ee.ge.1) then
        patcor(ee,:,:) = (/ linint2(sst_sd_ann&lon,sst_sd_ann&lat,sst_sd_ann,True,patcor&lon,patcor&lat,0) /)
     end if
     
     res@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     res@gsnRightString = sst_mean_djf@units
     res@gsnCenterString = names(ee)     
     plot_stddev_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,sst_sd_djf,res)
     plot_stddev_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,sst_sd_mam,res)
     plot_stddev_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,sst_sd_jja,res)
     plot_stddev_son(ee) = gsn_csm_contour_map(wks_stddev_son,sst_sd_son,res)
     plot_stddev_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,sst_sd_ann,res)
     
     sres@gsnLeftString = syear(ee)+"-"+eyear(ee) 
     sres@gsnRightString = sst_mean_djf@units
     sres@gsnCenterString = names(ee)
     plot_mean_djf(ee) = gsn_csm_contour_map(wks_mean_djf,sst_mean_djf,sres)
     plot_mean_mam(ee) = gsn_csm_contour_map(wks_mean_mam,sst_mean_mam,sres)
     plot_mean_jja(ee) = gsn_csm_contour_map(wks_mean_jja,sst_mean_jja,sres)
     plot_mean_son(ee) = gsn_csm_contour_map(wks_mean_son,sst_mean_son,sres)
     plot_mean_ann(ee) = gsn_csm_contour_map(wks_mean_ann,sst_mean_ann,sres)
     delete([/sst_sd_djf,sst_sd_mam,sst_sd_jja,sst_sd_son,sst_sd_ann,sst_mean_djf,sst_mean_mam,sst_mean_jja,sst_mean_son,sst_mean_ann,res,sres/])
  end do     
  
  if (isfilepresent2("obs_ts")) then    ; for pattern correlation table  
     clat = cos(0.01745329*patcor&lat)
     finpr   = "SST Std Dev (Ann) "    ; Must be 18 characters long
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
  
  panres@txString = "TS Standard Deviations (DJF)"
  gsn_panel(wks_stddev_djf,plot_stddev_djf,(/nrow,ncol/),panres)
  delete(wks_stddev_djf)
  
  panres@txString = "TS Standard Deviations (MAM)"
  gsn_panel(wks_stddev_mam,plot_stddev_mam,(/nrow,ncol/),panres)
  delete(wks_stddev_mam)
  
  panres@txString = "TS Standard Deviations (JJA)"
  gsn_panel(wks_stddev_jja,plot_stddev_jja,(/nrow,ncol/),panres)
  delete(wks_stddev_jja)
  
  panres@txString = "TS Standard Deviations (SON)"
  gsn_panel(wks_stddev_son,plot_stddev_son,(/nrow,ncol/),panres)
  delete(wks_stddev_son)
  
  panres@txString = "TS Standard Deviations (Annual)"
  gsn_panel(wks_stddev_ann,plot_stddev_ann,(/nrow,ncol/),panres)
  delete(wks_stddev_ann)
  
  panres@txString = "TS Means (DJF)"
  gsn_panel(wks_mean_djf,plot_mean_djf,(/nrow,ncol/),panres)
  delete(wks_mean_djf)
  
  panres@txString = "TS Means (MAM)"
  gsn_panel(wks_mean_mam,plot_mean_mam,(/nrow,ncol/),panres)
  delete(wks_mean_mam)
  
  panres@txString = "TS Means (JJA)"
  gsn_panel(wks_mean_jja,plot_mean_jja,(/nrow,ncol/),panres)
  delete(wks_mean_jja)
  
  panres@txString = "TS Means (SON)"
  gsn_panel(wks_mean_son,plot_mean_son,(/nrow,ncol/),panres)
  delete(wks_mean_son)
  
  panres@txString = "TS Means (Annual)"
  gsn_panel(wks_mean_ann,plot_mean_ann,(/nrow,ncol/),panres)
  delete(wks_mean_ann)
  delete(panres)
  print("Finished: sst.mean_stddev.ncl")
end
