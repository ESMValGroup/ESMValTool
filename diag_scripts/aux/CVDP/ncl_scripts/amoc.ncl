; Calculates MOC means/standard deviations, AMOC EOF1/PC1, 
; regressions onto AMOC PC1, and lag correlations vs. AMO
;
; Variables used: MOC, TS, TREFHT
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: amoc.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")
  
  nsim = numAsciiRow("namelist_byvar/namelist_moc")
  na = asciiread("namelist_byvar/namelist_moc",(/nsim/),"string")
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
 
;-----------TS/TAS read in for AMOC regressions------------------------------------------ 
  
  nsim_trefht = numAsciiRow("namelist_byvar/namelist_trefht")
  na_trefht = asciiread("namelist_byvar/namelist_trefht",(/nsim_trefht/),"string")
  names_trefht = new(nsim_trefht,"string")
  paths_trefht = new(nsim_trefht,"string")
  syear_trefht = new(nsim_trefht,"integer",-999)
  eyear_trefht = new(nsim_trefht,"integer",-999)

  do gg = 0,nsim_trefht-1
     names_trefht(gg) = str_strip(str_get_field(na_trefht(gg),1,delim))
     paths_trefht(gg) = str_strip(str_get_field(na_trefht(gg),2,delim))
     syear_trefht(gg) = stringtointeger(str_strip(str_get_field(na_trefht(gg),3,delim)))
     eyear_trefht(gg) = stringtointeger(str_strip(str_get_field(na_trefht(gg),4,delim)))
  end do
  delete(na_trefht)
  nyr_trefht = eyear_trefht-syear_trefht+1   
  
  nsim_ts = numAsciiRow("namelist_byvar/namelist_ts")
  na_ts = asciiread("namelist_byvar/namelist_ts",(/nsim_ts/),"string")
  names_ts = new(nsim_ts,"string")
  paths_ts = new(nsim_ts,"string")
  syear_ts = new(nsim_ts,"integer",-999)
  eyear_ts = new(nsim_ts,"integer",-999)

  do gg = 0,nsim_ts-1
     names_ts(gg) = str_strip(str_get_field(na_ts(gg),1,delim))
     paths_ts(gg) = str_strip(str_get_field(na_ts(gg),2,delim))
     syear_ts(gg) = stringtointeger(str_strip(str_get_field(na_ts(gg),3,delim)))
     eyear_ts(gg) = stringtointeger(str_strip(str_get_field(na_ts(gg),4,delim)))
  end do
  delete(na_ts)
  nyr_ts = eyear_ts-syear_ts+1   

  pi=4.*atan(1.0)
  rad=(pi/180.)

  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_mean = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.mean.ann")  
  wks_stddev = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.stddev.ann") 
  wks_amoc = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.ann") 
  wks_amoc_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.timeseries.ann")
  wks_amoc_powspec = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.powspec.ann")
  wks_amoc_sstreg  = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.sstreg.ann")
  wks_amoc_tasreg  = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc.tasreg.ann")
  wks_amoc_amo     = gsn_open_wks(wks_type,getenv("OUTDIR")+"amoc_amo.leadlag.ann")
 
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_mean,"rainbow+white")
     gsn_define_colormap(wks_stddev,"rainbow+white")
     gsn_define_colormap(wks_amoc,"ncl_default")
     gsn_define_colormap(wks_amoc_ts,"ncl_default")
     gsn_define_colormap(wks_amoc_powspec,"ncl_default")
     gsn_define_colormap(wks_amoc_sstreg,"ncl_default")
     gsn_define_colormap(wks_amoc_tasreg,"ncl_default")
     gsn_define_colormap(wks_amoc_amo,"ncl_default")
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_mean,"BlueDarkRed18")
     gsn_define_colormap(wks_stddev,"cb_rainbow")
     gsn_define_colormap(wks_amoc,"BlueDarkRed18")
     gsn_define_colormap(wks_amoc_ts,"ncl_default")
     gsn_define_colormap(wks_amoc_powspec,"cb_9step")
     gsn_define_colormap(wks_amoc_sstreg,"BlueDarkRed18")
     gsn_define_colormap(wks_amoc_tasreg,"BlueDarkRed18")
     gsn_define_colormap(wks_amoc_amo,"ncl_default")
  end if 
  plot_mean_ann    = new(nsim,"graphic")   
  plot_stddev_ann  = new(nsim,"graphic") 
  plot_amoc_ann    = new(nsim,"graphic") 
  plot_amoc_ts_ann = new(nsim,"graphic") 
  plot_amoc_powspec_ann = new(nsim,"graphic") 
  plot_amoc_sstreg_ann  = new(nsim,"graphic") 
  plot_amoc_tasreg_ann  = new(nsim,"graphic") 
  plot_amoc_amo_ann  = new(nsim,"graphic") 

  do ee = 0,nsim-1
     mocT = data_read_in_ocean_MOC(paths(ee),"MOC",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up        
     if (isatt(mocT,"is_all_missing")) then
        delete(mocT)
        continue
     end if 
     lat = tofloat(mocT&lat)
     lev = tofloat(mocT&lev)
     ny = dimsizes(lat)
     nz = dimsizes(lev)
         
     moc_ann = runave_n_Wrap(mocT,12,0,0)                             ; compute annual means and standard deviations
     moc_mean_ann = dim_avg_n_Wrap(moc_ann(5::12,:,:),0)
     moc_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),moc_ann(5::12,:,:),False,False,0),0)         
     copy_VarMeta(moc_mean_ann,moc_sd_ann)     
     delete(moc_ann)
     moc_sd_ann = where(moc_sd_ann.lt.0.001,moc_sd_ann@_FillValue,moc_sd_ann)
     moc_mean_ann = where(ismissing(moc_sd_ann),moc_mean_ann@_FillValue,moc_mean_ann) 
     
;-----------------------------------------------------------------------------------
;-----Calculate AMOC EOF1
;-----------------------------------------------------------------------------------     
     if (OPT_CLIMO.eq."Full") then
        mocT = rmMonAnnCycTLL(mocT)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = mocT
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(mocT&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        mocT   = calcMonAnomTLL(mocT,climo) 
        delete(climo)
     end if

     temp = runave_n_Wrap(mocT,12,0,0)   ; form annual averages
     amoc = temp(5::12,:,:)
     delete([/temp,mocT/])
     
     atl_begin = ind_nearest_coord (-33.0, lat, 0)          ; set missing values based on variance, and mask Atlantic south of 33S
     vareps = 1.e-6
     amocvar = conform(amoc,dim_variance_n_Wrap(amoc,0),(/1,2/))
     amoc@_FillValue = 1.e20
     amoc(:,:,0:atl_begin-1) = amoc@_FillValue
     amoc = where(amocvar.lt.vareps,amoc@_FillValue,amoc)
     delete(amocvar)     

     dlat = lat           ; Compute latitudinal weights (in meters)
     rEarth = 6.37122e8 		; Earth radius in m
     do iy=0,ny-1
        if (iy.gt.0.and.iy.lt.ny-1) then
           dy0 = (lat(iy)-lat(iy-1))/2.0
           dy1 = (lat(iy+1)-lat(iy))/2.0
           dlat(iy) = (dy0+dy1)*rEarth
           delete(dy0)
           delete(dy1)
        end if
        if (iy.eq.0) then
           dy1 = (lat(iy+1)-lat(iy))/2.0
           dlat(iy) = (2.*dy1)*rEarth
           delete(dy1)
        end if
        if (iy.eq.ny-1) then
           dy0 = (lat(iy)-lat(iy-1))/2.0
           dlat(iy) = (2.*dy0)*rEarth
           delete(dy0)
        end if
     end do
    
     dz = lev        ; compute vertical weights (in meters)
     do iz=0,nz-1
        if (iz.gt.0.and.iz.lt.nz-1) then
           dz(iz) = (lev(iz)-lev(iz-1))/2.0 + (lev(iz+1)-lev(iz))/2.0
        end if
        if (iz.eq.0) then
           dz(iz) = (lev(iz+1)-lev(iz))/2.0
        end if
        if (iz.eq.nz-1) then
           dz(iz) = (lev(iz)-lev(iz-1))/2.0
        end if
     end do
     
     test = dlat(0)*dz(0)
     wgt    = new((/nz,ny/),typeof(test))      ; weight the data
     delete(test)
     do iz=0,nz-1
        do iy=0,ny-1
           wgt(iz,iy) = dlat(iy)*dz(iz)
        end do
     end do
     amocW   = amoc*conform(amoc, wgt, (/1,2/))  ; same units as "amoc"
     delete(wgt)
     copy_VarMeta(amoc,amocW)
     amocW@long_name = "area weighted "+amoc@long_name

     workeof    = eofunc_Wrap(amocW(lev|:,lat|:,time|:), 3, 75)
     workeof_ts = eofunc_ts_Wrap (amocW(lev|:,lat|:,time|:), workeof, False)
     delete(amocW)
     amoc_pc_ann = dim_standardize(workeof_ts(0,:),0)
     moc_reg_ann = amoc(0,:,:)
     moc_reg_ann = (/ regCoef(amoc_pc_ann,amoc(lev|:,lat|:,time|:)) /)
     moc_reg_ann@pcvari = tofloat(sprintf("%4.1f", workeof@pcvar(0)))
     delete([/atl_begin,lat,lev,dz,dlat,amoc,workeof,workeof_ts,ny,nz/])
     if (max(moc_reg_ann&lev).ge.2000) then
        if (.not.ismissing(moc_reg_ann({2000.},{38}))) then
           if (moc_reg_ann({2000.},{38}).lt.0) then   ; arbitrary attempt to make all plots have the same sign..
              moc_reg_ann = moc_reg_ann*-1.
              amoc_pc_ann = amoc_pc_ann*-1.
           end if
        end if
     end if
;----------------------------------------------------------------------------------------     
     iopt = 0                         ; calculate spectra of AMOC PC1
     jave = (1*nyr(ee))/100
     if (jave.eq.0) then
       jave = 1
     end if
     val1 = .95
     val2 = .99
     pct = 0.1     
     spectra_mvf = False        ; missing value flag for nino3.4
     if (any(ismissing(amoc_pc_ann))) then
        print("Missing data detected for "+names(ee)+", power spectra function does not allow missing data, not creating AMOC spectra")
        spectra_mvf = True
;        if (isfilepresent2("obs_ts").and.ee.eq.0) then
;           spectra_mvf_obs = True    ; missing value flag for obs nino3.4
;        end if
     else
 ;       if (isfilepresent2("obs_ts").and.ee.eq.0) then
 ;          spectra_mvf_obs = False   ; missing value flag for obs nino3.4
 ;       end if
        sdof = specx_anal(amoc_pc_ann,iopt,jave,pct)
        splt1 = specx_ci(sdof,val1,val2)
        splt1!0 = "ncurves"
        splt1&ncurves = ispan(0,3,1)
        splt1!1 = "frequency2"
        splt1&frequency2 = sdof@frq
        splt1@units_info = "df refers to frequency interval; data are standardized so there are no physical units"
        splt1@units = "1/df"
        splt1@info = "(0,:)=spectrum,(1,:)=Markov red noise spectrum, (2,:)="+val1+"% confidence bound for Markhov, (3,:)="+val2+"% confidence bound for Markhov"

        delete([/iopt,jave,pct/])
     end if
;-------------------Read in TS and TAS for regressions onto PC1 and for AMO calculation----------------------------------------     
     trefht_offset = 0
     ts_offset = 0
     if (isfilepresent2("obs_ts")) then
        ts_offset = 1
     end if
     if (isfilepresent2("obs_trefht")) then
        trefht_offset = 1
     end if
     tas = data_read_in(paths_trefht(ee+trefht_offset),"TREFHT",syear_trefht(ee+trefht_offset),eyear_trefht(ee+trefht_offset))
     sst = data_read_in(paths_ts(ee+ts_offset),"TS",syear_ts(ee+ts_offset),eyear_ts(ee+ts_offset))
     delete([/ts_offset,trefht_offset/])
     
     sst = where(sst.le.-1.8,-1.8,sst)
     d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
     basemap = d->LSMASK
     lsm = landsea_mask(basemap,sst&lat,sst&lon)
     sst = mask(sst,conform(sst,lsm,(/1,2/)).ge.1,False)
     delete([/lsm,basemap/])
     delete(d)
     
     if (OPT_CLIMO.eq."Full") then
        sst = rmMonAnnCycTLL(sst)
        tas = rmMonAnnCycTLL(tas)
     else
        check_custom_climo(names_ts(ee),syear_ts(ee),eyear_ts(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = sst
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(sst&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        sst   = calcMonAnomTLL(sst,climo) 
        delete(climo)
        
        check_custom_climo(names_trefht(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = tas
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(tas&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        tas   = calcMonAnomTLL(tas,climo) 
        delete(climo)
     end if
     
     temp = runave_n_Wrap(tas,12,0,0)   ; form annual averages
     tas_ann = temp(5::12,:,:)
     delete([/temp,tas/])

     temp = runave_n_Wrap(sst,12,0,0)   ; form annual averages
     sst_ann = temp(5::12,:,:)
     delete([/temp,sst/])

     sst_reg_ann = sst_ann(0,:,:)                                                 ; SST regression onto AMOC PC1
     sst_reg_ann = (/ regCoef(amoc_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
     
     tas_reg_ann = tas_ann(0,:,:)                                                 ; TAS regression onto AMOC PC1
     tas_reg_ann = (/ regCoef(amoc_pc_ann,tas_ann(lat|:,lon|:,time|:)) /)
     delete(tas_ann)     
;-------------Compute AMO timeseries----------------------------------------------------------------------------------
     sst_ann = lonFlip(sst_ann)   ; orient longitudes from -180:180 (set to 0:360 in data_read_in function) for AMO calculation
     coswgt=cos(rad*sst_ann&lat)
     coswgt!0 = "lat"
     coswgt&lat= sst_ann&lat
     natl_aa =   wgt_areaave(sst_ann(:,{0:60},{-80:0}),coswgt({0.:60.}),1.0,0)
     global_aa = wgt_areaave(sst_ann(:,{-60:60},:),coswgt({-60.:60.}),1.0,0)
     
     AMO = new((/dimsizes(natl_aa)/),"float",-999.)    ; timeseries plot
     AMO!0 = "time"
     AMO&time = sst_ann&time
     AMO = (/ natl_aa - global_aa /)
     delete([/coswgt,natl_aa,global_aa,sst_ann/])
;---------------compute AMOC/AMO lead/lags----------------------------------------------------------------------------- 
     if (nyr(ee).ge.90) then    ; need a minimum number of years to compute lead/lag correlations
        nwt = 51
        pda = 15            ; longest period
        pdb = 1            ; shortest period
        fca = 1./pda       ;  ==> lowest allowed frequency
        fcb = 1./pdb       ;  ==> highest allowed frequency
        ihp = 0            ;  0 ==> low pass filter, fcb ignored
        nsigma = 1.
        twgt = filwgts_lanczos (nwt, ihp, fca, fcb, nsigma)
        AMO_wgt = wgt_runave_Wrap(AMO,twgt,0)
        amoc_pc_ann_wgt = wgt_runave_Wrap(amoc_pc_ann,twgt,0)

        mxlag = 15
        x_Lead_y = esccr(amoc_pc_ann_wgt,AMO_wgt,mxlag)
        y_Lead_x = esccr(AMO_wgt,amoc_pc_ann_wgt,mxlag)    ; switch the order of the series
     
        ccr = new ( 2*mxlag+1, float)    
        ccr(0:mxlag-1) = y_Lead_x(1:mxlag:-1)  ; "negative lag", -1 reverses order
        ccr(mxlag:)    = x_Lead_y(0:mxlag)     ; "positive lag"
        delete([/x_Lead_y,y_Lead_x,AMO_wgt,amoc_pc_ann_wgt/])
     end if
     delete(AMO)
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
     
        if (spectra_mvf.eq.False) then
           z->amoc_spectra_ann = splt1
        end if
        z->amoc_sst_regression_ann = sst_reg_ann
        z->amoc_tas_regression_ann = tas_reg_ann
        if (nyr(ee).ge.90) then
           time_lag_cor = ispan(mxlag*-1,mxlag,1)
           time_lag_cor@units = "months since 0000-01-01 00:00:00"
           time_lag_cor@long_name = "Time"
           time_lag_cor@standard_name = "time"
           time_lag_cor@calendar = "standard"
           time_lag_cor!0 = "time_lag_cor"
           time_lag_cor&time_lag_cor = time_lag_cor
           ccr!0 = "time_lag_cor"
           ccr&time_lag_cor = time_lag_cor
           z->amoc_amo_lag_cor = ccr
           delete(time_lag_cor)
        end if
        TIME = ispan(0,dimsizes(amoc_pc_ann)-1,1)
        TIME@units = "years since "+syear(ee)+"-01-15 00:00:00"
        TIME@long_name = "Time"
        TIME@standard_name = "time"
        TIME@calendar = "standard"
        TIME!0 = "TIME"
        TIME&TIME = TIME
        amoc_pc_ann!0 = "TIME"
        amoc_pc_ann&TIME = TIME
        z->amoc_timeseries_ann = amoc_pc_ann
        lat_amoc = moc_reg_ann&lat
        lat_amoc!0 = "lat_amoc"
        lat_amoc&lat_amoc = lat_amoc
        delete(moc_reg_ann&lat)
        moc_reg_ann!1 = "lat_amoc"
        moc_reg_ann&lat_amoc = lat_amoc
        delete(moc_mean_ann&lat)
        moc_mean_ann!1 = "lat_amoc"
        moc_mean_ann&lat_amoc = lat_amoc
        delete(moc_sd_ann&lat)
        moc_sd_ann!1 = "lat_amoc"
        moc_sd_ann&lat_amoc = lat_amoc
        z->amoc_mean_ann   = moc_mean_ann
        z->amoc_stddev_ann = moc_sd_ann
        z->amoc_ann        = moc_reg_ann 
        delete([/modname,fn,TIME,lat_amoc/])
     end if     

     moc_mean_ann&lev = moc_mean_ann&lev/1000.
     moc_mean_ann&lev@units = "km"
     moc_sd_ann&lev = moc_sd_ann&lev/1000.
     moc_sd_ann&lev@units = "km"
     moc_reg_ann&lev = moc_reg_ann&lev/1000.
     moc_reg_ann&lev@units = "km"
;=======================================================================================================================     
     mocres                      = True                 ; plot mods desired
     mocres@gsnDraw		   = False		
     mocres@gsnFrame		   = False		
     mocres@cnFillOn             = True                 ; turn on color fill
     mocres@cnMissingValFillColor = "gray50"
     mocres@cnMissingValFillPattern = 0
     mocres@cnLinesOn            = True         
     mocres@cnLineLabelsOn 	   = False
;     mocres@cnLineLabelFontColor = "white"
;     mocres@cnLineLabelBackgroundColor = -1
     mocres@lbLabelBarOn        = False

     mocres@cnInfoLabelOn       = False       ; Turn off informational label
     mocres@cnLevelSelectionMode = "ExplicitLevels"	; manually set the contour levels
     if (nsim.le.7) then
        mocres@tmXBLabelFontHeightF = 0.01
        mocres@tmYLLabelFontHeightF = 0.01
        mocres@gsnLeftStringFontHeightF   = 0.0125
        mocres@gsnCenterStringFontHeightF = 0.0125
        mocres@gsnRightStringFontHeightF  = 0.011    
     else
        mocres@tmXBLabelFontHeightF = 0.014
        mocres@tmYLLabelFontHeightF = 0.014
        mocres@gsnLeftStringFontHeightF   = 0.017
        mocres@gsnCenterStringFontHeightF = 0.017
        mocres@gsnRightStringFontHeightF  = 0.0155
     end if
     mocres@tiYAxisFontHeightF = mocres@tmXBLabelFontHeightF 
     mocres@gsnCenterStringOrthogonalPosF = -0.96
     mocres@gsnCenterStringParallelPosF = 0.80
     
     mocres@tmXBLabelsOn          = True
     mocres@tmXTLabelsOn          = False
     mocres@tmXTOn          = False
     mocres@tmYRLabelsOn          = False
     mocres@tmYROn          = False
     mocres@cnMonoLineThickness = False
     mocres@cnMonoLineDashPattern = False
     mocres@vpWidthF        = 0.375
     mocres@vpHeightF       = 0.28

     mocres@trYReverse           = True       ; reverses y-axis
     mocres@gsnYAxisIrregular2Linear = True
     mocres@gsnXAxisIrregular2Linear = True
     mocres@tiYAxisString= "depth (km)"
     mocres@tiXAxisString= ""
     mocres@gsnCenterString = ""
     mocres@trXMinF = 0.
     mocres@trXMaxF = 90.
     mocres@tmXBMode = "Explicit"
     mocres@tmXBValues = (/0.,30.,60.,90./)
     mocres@tmXBLabels = (/"0~S~o~N~N","30~S~o~N~N","60~S~o~N~N","90~S~o~N~N"/)

     mocres@gsnCenterString = syear(ee)+"-"+eyear(ee) 
     mocres@gsnRightString = moc_mean_ann@units
     mocres@gsnLeftString = names(ee)
 
 
     mocres@cnLevels = ispan(-4,28,2)
     mocres@cnLineThicknesses  = (/1,1,2,1,1,1,1,2,1,1,1,1,2,1,1,1,1/)
     mocres@cnLineDashPatterns = (/1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
     plot_mean_ann(ee) = gsn_csm_contour(wks_mean,moc_mean_ann,mocres)  
     delete(moc_mean_ann)  
     delete(mocres@cnLevels)
     delete(mocres@cnLineThicknesses)
     delete(mocres@cnLineDashPatterns)
     
     mocres@cnLevels = (/.1,.3,.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0/)  
     if (COLORMAP.eq.0) then     
        mocres@cnFillColors = (/20,38,54,80,95,125,175,185,195,205,236/)
     end if
     if (COLORMAP.eq.1) then     
        mocres@cnFillColors = (/14,23,35,47,63,79,95,111,124,155,175/)
     end if
     plot_stddev_ann(ee) = gsn_csm_contour(wks_stddev,moc_sd_ann,mocres)
     delete(moc_sd_ann)
     delete(mocres@cnLevels)
     delete(mocres@cnFillColors)

     mocres@cnLevels = fspan(-2,2,41)
     mocres@cnMonoLineThickness = True
     mocres@cnMonoLineDashPattern = True
     mocres@cnLineDashPattern = 0
     mocres@gsnCenterString = syear(ee)+"-"+eyear(ee) 
     mocres@gsnRightString = moc_reg_ann@pcvari+"%"
     mocres@gsnLeftString = names(ee)
     plot_amoc_ann(ee) = gsn_csm_contour(wks_amoc,moc_reg_ann,mocres)  
     delete(moc_reg_ann)  
     delete(mocres)    

     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnXYBarChart = False
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"
     xyres@gsnAboveYRefLineColor = 185
     xyres@gsnBelowYRefLineColor = 35
     if (wks_type.eq."png") then  
        xyres@xyLineThicknessF = .5
     else
        xyres@xyLineThicknessF = .2
     end if
     xyres@xyLineColor = "gray52"
     xyres@tiYAxisString = ""
     xyres@tiXAxisString = ""
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
     xyres@gsnLeftString = ""
     xyres@gsnRightString = ""
     xyres@trXMinF = syear(ee)-.5
     xyres@trXMaxF = eyear(ee)+1.5

     xyres@gsnCenterString = names(ee)
     plot_amoc_ts_ann(ee) = gsn_csm_xy(wks_amoc_ts,ispan(syear(ee),eyear(ee),1),amoc_pc_ann,xyres)  ; use standardized timeseries  
     delete(amoc_pc_ann)
     delete(xyres)
    
     pres = True
     pres@vpXF = 0.07
     pres@trYMinF = 0.
     pres@trXMinF = 0.0
     pres@trXMaxF = 0.5
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
     pres@tmXTValues = (/".02",".10",".20",".3333",".50"/)
     pres@tmXTLabels = (/"50","10","5","3","2"/)         

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
     pres@gsnCenterString = names(ee)
     if (spectra_mvf.eq.False) then
        plot_amoc_powspec_ann(ee) = gsn_csm_xy(wks_amoc_powspec,sdof@frq,splt1,pres) 
        delete([/splt1,sdof/])
     end if
     delete(pres)
     
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
     if (wks_type.eq."png") then
        res@mpGeophysicalLineThicknessF = 2.  
     else
        res@mpGeophysicalLineThicknessF = 1.  
     end if
     res@gsnDraw      = False
     res@gsnFrame     = False
     res@vpYF = 0.95
     res@vpHeightF = 0.3
     res@vpXF = 0.2
     res@vpWidthF = 0.6
  
;     res@cnFillMode = "RasterFill"
     res@cnLevelSelectionMode = "ExplicitLevels"
     if (COLORMAP.eq."0") then
        res@cnLevels = fspan(-.5,.5,21)
     end if
     if (COLORMAP.eq."1") then
        res@cnLevels = fspan(-.4,.4,17)
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

     plot_amoc_sstreg_ann(ee) = gsn_csm_contour_map(wks_amoc_sstreg,sst_reg_ann,res)
     delete(sst_reg_ann)
     plot_amoc_tasreg_ann(ee) = gsn_csm_contour_map(wks_amoc_tasreg,tas_reg_ann,res)
     delete(tas_reg_ann)
     delete(res)
     
     res2 = True
     res2@gsnDraw = False
     res2@gsnFrame = False
     res2@trYMinF = -1.  ; min((/-0.4,min(ccr)/))
     res2@trYMaxF = 1.   ; max((/0.6,max(ccr)/))
     res2@vpWidthF = 0.6
     res2@vpHeightF = 0.4
     res2@gsnYRefLine = 0.0
     res2@gsnYRefLineColor = "gray42"
     res2@gsnXRefLine = 0.0
     res2@gsnXRefLineColor = "gray42"
     res2@xyLineColor = "royalblue"
     if (wks_type.eq."png") then
        res2@xyLineThicknessF = 3.5
     else
        res2@xyLineThicknessF = 1.75
     end if
     if (nsim.le.5) then
        res2@tmXBLabelFontHeightF = 0.0125
        res2@tmYLLabelFontHeightF = 0.0125
        res2@gsnLeftStringFontHeightF = 0.013
        res2@gsnCenterStringFontHeightF = 0.017
        res2@gsnRightStringFontHeightF = 0.013         
     else
        res2@tmXBLabelFontHeightF = 0.018
        res2@tmYLLabelFontHeightF = 0.018
        res2@gsnLeftStringFontHeightF = 0.020
        res2@gsnCenterStringFontHeightF = 0.024
        res2@gsnRightStringFontHeightF = 0.020    
     end if
     res2@gsnLeftStringOrthogonalPosF = -1.01
     res2@gsnRightStringOrthogonalPosF = -1.01
     res2@gsnLeftStringParallelPosF = 0.01
     res2@gsnRightStringParallelPosF = 0.99
     
     res2@gsnLeftString = "AMO leads"
     res2@gsnCenterString = names(ee)
     res2@gsnRightString = "AMOC PC1 leads"
     if (nyr(ee).ge.90) then    ; need a minimum number of years to compute lead/lag correlations
        res2@trXMinF = mxlag*-1
        res2@trXMaxF = mxlag
        plot_amoc_amo_ann(ee) = gsn_csm_xy(wks_amoc_amo,ispan(mxlag*-1,mxlag,1),ccr,res2)
        delete([/mxlag,ccr/])
     end if
     delete(res2)
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
  panres@lbLabelFontHeightF = 0.013
  panres@lbLabelStride = 1
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  
  panres@txString = "AMOC Means (Annual)"
  gsn_panel(wks_mean,plot_mean_ann,(/nrow,ncol/),panres)
  delete(wks_mean)

  panres@txString = "AMOC Standard Deviations (Annual)"
  gsn_panel(wks_stddev,plot_stddev_ann,(/nrow,ncol/),panres)
  delete(wks_stddev)


  panres@txString = "AMOC (Annual)"
  gsn_panel(wks_amoc,plot_amoc_ann,(/nrow,ncol/),panres)
  delete(wks_amoc)

  panres@txString = "AMOC TAS Regressions (Annual)"
  gsn_panel(wks_amoc_tasreg,plot_amoc_tasreg_ann,(/nrow,ncol/),panres)
  delete(wks_amoc_tasreg)

  panres@txString = "AMOC SST Regressions (Annual)"
  gsn_panel(wks_amoc_sstreg,plot_amoc_sstreg_ann,(/nrow,ncol/),panres)
  delete(wks_amoc_sstreg)
  delete(panres)

  panres2 = True    
  if (nsim.le.5) then
     panres2@txFontHeightF = 0.024
  else
     panres2@txFontHeightF = 0.016
  end if  
  panres2@gsnMaximize = True
  panres2@gsnPaperOrientation = "portrait"
  if (nsim.le.12) then
     lp = (/nsim,1/)
  else
     lp = (/nrow,ncol/)  
  end if
  panres2@gsnPanelYWhiteSpacePercent = 3.0
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres2@gsnPanelScalePlotIndex = tt(0)
     delete(tt)
  end if
  panres2@txString = "AMOC PC1 (Annual)"
  gsn_panel(wks_amoc_ts,plot_amoc_ts_ann,(/nrow,ncol/),panres2)
  delete(wks_amoc_ts)

  if (isatt(panres2,"gsnPanelScalePlotIndex")) then
     delete(panres2@gsnPanelScalePlotIndex)
  end if
  panres2@txString = "AMOC PC1 (Annual, detrended)"
  gsn_panel(wks_amoc_powspec,plot_amoc_powspec_ann,(/nrow,ncol/),panres2)
  delete(wks_amoc_powspec)


  panres2@txString = "AMO / AMOC PC1 Lag Correlation (Annual)"
  gsn_panel(wks_amoc_amo,plot_amoc_amo_ann,(/nrow,ncol/),panres2)
  delete(wks_amoc_amo)
  delete(panres2)
  print("Finished: amoc.ncl")
end


