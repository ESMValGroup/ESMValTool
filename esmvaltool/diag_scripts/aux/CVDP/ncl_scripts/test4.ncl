; Calculates precipitation global trends and timeseries
;
; Variables used: pr
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions_dk.ncl"

begin
  print("Starting: test4.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")  
;;------------------PRECT---------------------------------------------  
  nsim_ppt = numAsciiRow("namelist_byvar/namelist_prect")
  na_ppt = asciiread("namelist_byvar/namelist_prect",(/nsim_ppt/),"string")
  names_ppt = new(nsim_ppt,"string")
  paths_ppt = new(nsim_ppt,"string")
  syear_ppt = new(nsim_ppt,"integer",-999)
  eyear_ppt = new(nsim_ppt,"integer",-999)
  delim = "|"

  do gg = 0,nsim_ppt-1
     names_ppt(gg) = str_strip(str_get_field(na_ppt(gg),1,delim))
     paths_ppt(gg) = str_strip(str_get_field(na_ppt(gg),2,delim))
     syear_ppt(gg) = stringtointeger(str_strip(str_get_field(na_ppt(gg),3,delim)))
     eyear_ppt(gg) = stringtointeger(str_strip(str_get_field(na_ppt(gg),4,delim)))
  end do
  nyr_ppt = eyear_ppt-syear_ppt+1
  delete(na_ppt) 
;;----------------TS-------------------------------------------------
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
 ;;------------------CLT--------------------------------------------
  nsim_clt = numAsciiRow("namelist_byvar/namelist_clt")
  na_clt = asciiread("namelist_byvar/namelist_clt",(/nsim_clt/),"string")
  names_clt = new(nsim_clt,"string")
  paths_clt = new(nsim_clt,"string")
  syear_clt = new(nsim_clt,"integer",-999)
  eyear_clt = new(nsim_clt,"integer",-999)
  delim = "|"

  do gg = 0,nsim_clt-1
     names_clt(gg) = str_strip(str_get_field(na_clt(gg),1,delim))
     paths_clt(gg) = str_strip(str_get_field(na_clt(gg),2,delim))
     syear_clt(gg) = stringtointeger(str_strip(str_get_field(na_clt(gg),3,delim)))
     eyear_clt(gg) = stringtointeger(str_strip(str_get_field(na_clt(gg),4,delim)))
  end do
  nyr_clt = eyear_clt-syear_clt+1
  delete(na_clt)
;;----------------------------------------------------------------------------------
  
  ;wgt = (/1.,2.,1./)   
  ;wgt = wgt/sum(wgt)
  pi=4.*atan(1.0)
  rad=(pi/180.)
  rad1 =flt2dble(rad)
  delete(rad)
  rad = rad1
  delete(rad1)

  n34sc_frame = 1     ; n34sc_frame = flag to create regressions .ps/.png files. Created/used instead of n34sc_plot_flag
                      ; so that if sst regressions are not created for the last simulation listed that .ps/png files are created

  if (isfilepresent("obs_ts")) then
    if (nsim.eq.nsim_ppt.and.nsim.eq.nsim_clt) then
                                ; print("Total number of simulations/obs are equal between psl/trefht/ts namelists, proceeding")
      n34sc_obs = 0
    else
      print("Total number of simulations/obs between ts/clt and ts/precip do not match, skipping")
      n34sc_obs = 1
    end if 
  else
    n34sc_obs = 1
  end if
  if (n34sc_obs.eq.0) then
    if (syear(0).eq.syear_ppt(0).and.syear(0).eq.syear_clt(0)) then
      if (eyear(0).eq.eyear_clt(0).and.eyear(0).eq.eyear_ppt(0)) then
        n34sc_obs = 0
      else
        print("Observational end years do not match, skipping")
        n34sc_obs = 1
      end if
    else
      print("Observational start years do not match, skipping nino3.4 spatial composite")
      n34sc_obs = 1
    end if
  end if
 
;-----------------------------------------------------------------------------------
   wks_type = OUTPUT_TYPE
   
   if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
   end if
   
   wks_aa_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst_precip.timeseries")
   wks_aa_mon1 = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst_cloud.timeseries")
   wks_corr_cloud = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34_cloud.corr")
   wks_corr_rain = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34_precip.corr")
   wks_sst_reg = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34_sst_reg")
   wks_n34sc_clt = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.spatialcomp.clt")
   wks_n34sc_ppt = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.spatialcomp.ppt")
   wks_lagppt = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst_precip.lag")
   wks_lagclt = gsn_open_wks(wks_type,getenv("OUTDIR")+"sst_clt.lag")
   
   if (COLORMAP.eq.0) then
      gsn_define_colormap(wks_corr_rain,"ncl_default")
      gsn_define_colormap(wks_corr_cloud,"ncl_default")
      gsn_define_colormap(wks_aa_mon,"ncl_default") 
      gsn_define_colormap(wks_aa_mon1,"ncl_default")
      gsn_define_colormap(wks_sst_reg,"ncl_default")
      gsn_define_colormap(wks_n34sc_clt,"NCV_blue_red")
      gsn_define_colormap(wks_n34sc_ppt,"MPL_PuOr")
      gsn_define_colormap(wks_lagppt,"ncl_default") 
      gsn_define_colormap(wks_lagclt,"ncl_default") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_corr_rain,"ncl_default")
     gsn_define_colormap(wks_corr_cloud,"ncl_default")
     gsn_define_colormap(wks_aa_mon,"ncl_default") 
     gsn_define_colormap(wks_aa_mon1,"ncl_default")  
     gsn_define_colormap(wks_sst_reg,"ncl_default")
     gsn_define_colormap(wks_n34sc_clt,"BlueDarkRed18")   
     gsn_define_colormap(wks_n34sc_ppt,"BlueDarkRed18")   
     gsn_define_colormap(wks_lagppt,"ncl_default")
     gsn_define_colormap(wks_lagclt,"ncl_default") 
 end if
   
   plot_corr_rain = new(nsim,"graphic")
   plot_corr_cloud = new(nsim,"graphic")
   plot_sst_reg = new(nsim, "graphic")
   
   lagppt_xy = new(nsim,"graphic")
   lagclt_xy = new(nsim,"graphic")
   
   xy_mon = new(nsim,"graphic")   
   xy_mon1 = new(nsim, "graphic")

   map_n34sc_jja0_clt = new(nsim,"graphic")  
   map_n34sc_son0_clt = new(nsim,"graphic")  
   map_n34sc_djf1_clt = new(nsim,"graphic")  
   map_n34sc_mam1_clt = new(nsim,"graphic")  

   map_n34sc_jja0_ppt = new(nsim,"graphic")  
   map_n34sc_son0_ppt = new(nsim,"graphic")  
   map_n34sc_djf1_ppt = new(nsim,"graphic")  
   map_n34sc_mam1_ppt = new(nsim,"graphic")
   
;  if (isfilepresent("obs_prect")) then
;    xy_obs_mon = new(nsim,"graphic")  
;  end if
;------------------ Starting Loop Models and Observations
  
   do ee = 0,nsim-1
    sst = data_read_in(paths(ee),"TS",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
    
    if (isatt(sst,"is_all_missing")) then
      delete(sst)
      continue
    end if
    
    if (isfilepresent("obs_ts").and.ee.eq.0) then   ; mask out land in array going into EOF functions
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
    
    if (OPT_CLIMO.eq."Full") then
      sst = rmMonAnnCycTLL(sst)
    else
      check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
      temp_arr = sst
      delete(temp_arr&time)
      temp_arr&time = cd_calendar(sst&time,-1)
      climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
      delete(temp_arr)
      sst   = calcMonAnomTLL(sst,climo) 
      delete(climo)
    end if
    
    coswgt=cos(rad*sst&lat)
    coswgt!0 = "lat"
    coswgt&lat= sst&lat
    
    llats = -5.     ; nino3.4
    llatn = 5.
    llonw = 190.
    llone = 240.
    nino34 = wgt_areaave_Wrap(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
    nino34@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
    nino34@units = sst@units
    nino34@long_name = "nino3.4"
    
    sst_nino34 = sst(0,:,:)
    sst_nino34 = (/ regCoef(nino34,sst(lat|:,lon|:,time|:)) /)

    plats = -35. ;; lat/lon composite sst vs. cloud/precip timeseries    
    platn = -15.
    plonw = 115.
    plone = 150.
    
    sst_aa_mon = wgt_areaave_Wrap(sst(:,{plats:platn},{plonw:plone}),coswgt({plats:platn}),1.0,0)  
    sst_aa_mon@area = plats+":"+platn+"N, "+plonw+":"+plone+"E"
    sst_aa_mon@units = sst@units
    sst_aa_mon@long_name = "sst anomalies "
    area = sst_aa_mon@area
    delete([/coswgt/])
;---------------------------------------------------- Nino3.4 Composites---------------------------------------------------
    wgt = (/1.,2.,1./)   
    wgt = wgt/sum(wgt)

    coswgt=cos(rad*sst&lat)
    coswgt!0 = "lat" 
    coswgt&lat= sst&lat 

    nino34T = wgt_runave(nino34,wgt,1)     ; for use in ENSO composites / hovmuellers / running standard deviations
    nino34_ndj = nino34T(11:dimsizes(nino34T)-13:12)   ; cannot count last 1yr as spatial composite uses +1yrs data betond NDJ..  
    nino34_ndj!0 = "time"
    nino34_ndj&time = ispan(syear(ee),eyear(ee)-1,1)    
    nino34_ndj = dtrend_msg(ispan(0,dimsizes(nino34_ndj&time)-1,1),nino34_ndj,True,False)
    nino34_ndj = dim_standardize(nino34_ndj,0)
    
    sst = (/ dtrend_msg_n(ispan(0,nyr(ee)*12-1,1),sst,False,False,0) /) ; detrend the sst array
    
    
    delete([/coswgt/])
;--------------------------------------------------------------------------------------------------------------------------
    n34sc_plot_flag = 0   
    if (ee.eq.0.and.n34sc_obs.eq.1) then    ; first time thru, and obs = False
      n34sc_plot_flag = 1                  ; do not plot
      if (.not.isfilepresent("obs_ts")) then  ; however, if there are no TS obs whatsoever, the 1st iteration of the
        n34sc_plot_flag = 0                  ; ee do loop will be the first model, so turn it on.
      end if
    end if
      
    if (n34sc_plot_flag.eq.0) then        
      ppt_offset    = 0
      clt_offset = 0
      if (ee.eq.0.and.n34sc_obs.eq.0) then
      else
        if (isfilepresent("obs_ts")) then
          if (.not.isfilepresent("obs_prect")) then
            ppt_offset = -1
          end if
          if (.not.isfilepresent("obs_clt")) then
            clt_offset = -1
          end if
        else
          if (isfilepresent("obs_prect")) then
            ppt_offset = 1
          end if
          if (isfilepresent("obs_clt")) then
            clt_offset = 1
          end if
        end if
      end if
      
      ppt = data_read_in(paths_ppt(ee+ppt_offset),"PRECT",syear_ppt(ee+ppt_offset),eyear_ppt(ee+ppt_offset))
      clt_double = data_read_in(paths_clt(ee+clt_offset),"CLT",syear_clt(ee+clt_offset),eyear_clt(ee+clt_offset))
      
      clt = dble2flt(clt_double)
      delete(clt_double)
 ;----------------------------------Nino3.4 composites section-----------------------------------------------------
        TIME = sst&time
        yyyymm = cd_calendar(sst&time,-1)  ; convert tas, ts, and sst from CF-conforming time to YYYYMM for coding below
        delete(sst&time)
        sst&time = yyyymm 
        delete(yyyymm)

        yyyymm = cd_calendar(clt&time,-1)  
        delete(clt&time)
        clt&time = yyyymm 
        delete(yyyymm)

        yyyymm = cd_calendar(ppt&time,-1) 
        delete(ppt&time)
        ppt&time = yyyymm 
        delete(yyyymm)
;----------------------------------------------------------------------------------------------------------------
        if (isatt(ppt,"is_all_missing").or.isatt(clt,"is_all_missing")) then
          n34sc_plot_flag = 1
          delete([/clt,ppt/])
        end if
        
        if (n34sc_plot_flag.eq.0) 
          if (OPT_CLIMO.eq."Full") then
            ppt = rmMonAnnCycTLL(ppt)
            clt = rmMonAnnCycTLL(clt)           
          else
            check_custom_climo(names_ppt(ee+ppt_offset),syear_ppt(ee+ppt_offset),eyear_ppt(ee+ppt_offset),CLIMO_SYEAR,CLIMO_EYEAR)
            climo = clmMonTLL(ppt({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
            ppt   = calcMonAnomTLL(ppt,climo) 
            delete(climo)
            
            check_custom_climo(names_clt(ee+clt_offset),syear_clt(ee+clt_offset),eyear_clt(ee+clt_offset),CLIMO_SYEAR,CLIMO_EYEAR)
            climo = clmMonTLL(clt({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
            clt   = calcMonAnomTLL(clt,climo) 
            delete(climo)
          end if
          
          coswgt=cos(rad*ppt&lat)
          coswgt!0 = "lat"
          coswgt&lat= ppt&lat
          
          ppt_aa_mon = wgt_areaave_Wrap(ppt(:,{plats:platn},{plonw:plone}),coswgt({plats:platn}),1.0,0)
          ppt_aa_mon@area = plats+":"+platn+"N, "+plonw+":"+plone+"E"
          ppt_aa_mon@long_name = "precipitation"
          delete([/coswgt/])
          
          coswgt=cos(rad*clt&lat)
          coswgt!0 = "lat"
          coswgt&lat= clt&lat
          
          clt_aa_mon = wgt_areaave_Wrap(clt(:,{plats:platn},{plonw:plone}),coswgt({plats:platn}),1.0,0)
          clt_aa_mon@area = plats+":"+platn+"N, "+plonw+":"+plone+"E"
          clt_aa_mon@long_name = "cloud cover"
          delete([/coswgt/])
          
        end if
      end if
      
;---------------------------------Nino3.4 composites section-------------------------------------
      ppt = (/ dtrend_msg_n(ispan(0,dimsizes(ppt&time)-1,1),ppt,False,False,0) /)     ; sst detrended up above
      clt = (/ dtrend_msg_n(ispan(0,dimsizes(clt&time)-1,1),clt,False,False,0) /)  
      
      ta = dim_avg_n(sst(:1,:,:),0)
      sst = runave_n_Wrap(sst,3,0,0)
      sst(0,:,:) = (/ ta /)
      delete(ta)
      
      ta = dim_avg_n(ppt(:1,:,:),0)
      ppt = runave_n_Wrap(ppt,3,0,0)
      ppt(0,:,:) = (/ ta /)
      delete(ta)
      
      ta = dim_avg_n(clt(:1,:,:),0)
      clt = runave_n_Wrap(clt,3,0,0)
      clt(0,:,:) = (/ ta /)
      delete(ta)
    

      hicntr = 0
      locntr = 0
      hiyr = new(dimsizes(nino34_ndj&time),integer)
      loyr = hiyr
      
;           print(nino34_ndj)
      
      do hh = 0,dimsizes(nino34_ndj)-1
        if (nino34_ndj(hh).ge.1) then
                                ; print("High year = "+nino34_ndj&time(hh))
          hiyr(hicntr) = nino34_ndj&time(hh)
          hicntr = hicntr+1
        end if
        if (nino34_ndj(hh).le.-1) then
                                ; print("Low year = "+nino34_ndj&time(hh))
          loyr(locntr) = nino34_ndj&time(hh)
          locntr = locntr+1
        end if
      end do
      
      if (hicntr.eq.0) then     ; for simulations with climatological SSTs
        highyr = hiyr(0)
      else
        highyr = hiyr(:hicntr-1)
      end if
      delete([/hiyr,hicntr/])
      if (locntr.eq.0) then
        lowyr = loyr(0)
      else
        lowyr = loyr(:locntr-1)
      end if
      delete([/loyr,locntr/])
      
      dimS = dimsizes(ppt&time)    ; change time from YYYYMM->YYYY.frac
      tmin = ppt&time(0)/100
      tmax = ppt&time(dimS-1)/100
      delete(ppt&time)
      ppt&time = fspan(tmin*1.,(tmax*1.)+(11/12.),dimS)
      dimS = dimsizes(clt&time)
      tmin = clt&time(0)/100
      tmax = clt&time(dimS-1)/100
      delete(clt&time)
      clt&time = fspan(tmin*1.,(tmax*1.)+(11/12.),dimS)
      dimS = dimsizes(sst&time)
      tmin = sst&time(0)/100
      tmax = sst&time(dimS-1)/100
      delete(sst&time)
      sst&time = fspan(tmin*1.,(tmax*1.)+(11/12.),dimS)
      delete([/dimS,tmin,tmax/])
                                ; print(sst&time)
        
      sc_clt_hi = clt(:23,:,:)
      sc_clt_lo = clt(:23,:,:)    
      sc_sst_hi = sst(:23,:,:)
      sc_sst_lo = sst(:23,:,:)  
      sc_ppt_hi = ppt(:23,:,:)
      sc_ppt_lo = ppt(:23,:,:)  
      
      sc_clt_hi = sc_clt_hi@_FillValue
      sc_clt_lo = sc_clt_lo@_FillValue
      sc_sst_hi = sc_sst_hi@_FillValue
      sc_sst_lo = sc_sst_lo@_FillValue
      sc_ppt_hi = sc_ppt_hi@_FillValue
      sc_ppt_lo = sc_ppt_lo@_FillValue
      
      if (dimsizes(highyr).le.1) then
        print("For "+names(ee)+", 1 or less (normalized) nino3.4 value greater than one standard deviation found, setting nino3.4 spatial composites to missing")  ; sc_*_hi arrays left to _FillValue
        sc_clt_hi@is_all_missing = True
      else           
        do gg = 0,23
          tt = gg/12.
          sc_ppt_hi(gg,:,:) = (/ dim_avg_n(ppt({highyr+tt},:,:),0) /)
          sc_sst_hi(gg,:,:) = (/ dim_avg_n(sst({highyr+tt},:,:),0) /)    
          sc_clt_hi(gg,:,:) = (/ dim_avg_n(clt({highyr+tt},:,:),0) /)
        end do
        delete([/tt,highyr/])
      end if
      if (dimsizes(lowyr).le.1) then
        print("For "+names(ee)+", 1 or less (normalized) nino3.4 value less than -1 standard deviation found, setting nino3.4 spatial composites to missing")   ; sc_*_lo arrays left to _FillValue
      else           
        do gg = 0,23
          tt = gg/12.
          sc_ppt_lo(gg,:,:) = (/ dim_avg_n(ppt({lowyr+tt},:,:),0) /)
          sc_sst_lo(gg,:,:) = (/ dim_avg_n(sst({lowyr+tt},:,:),0) /)     
          sc_clt_lo(gg,:,:) = (/ dim_avg_n(clt({lowyr+tt},:,:),0) /)
        end do
        delete([/tt,lowyr/])
      end if
      
      n34sc_ppt = sc_ppt_hi
      n34sc_ppt = (/ sc_ppt_hi - sc_ppt_lo /)
      n34sc_sst = sc_sst_hi
      n34sc_sst = (/ sc_sst_hi - sc_sst_lo /)
      n34sc_clt = sc_clt_hi
      n34sc_clt = (/ sc_clt_hi - sc_clt_lo /)
      delete([/sc_clt_hi,sc_clt_lo,sc_sst_hi,sc_sst_lo,sc_ppt_hi,sc_ppt_lo/])
      delete(sst&time)
      sst&time = TIME    
      delete(TIME)

      if (n34sc_frame.eq.1.and.n34sc_plot_flag.eq.0) then    ; n34sc_frame = flag to create composites .ps/.png files 
        n34sc_frame = 0
      end if   
      if (isvar("TIME")) then
        delete(TIME)
      end if
    
;-------------------------------Nino34 vs. PRECIP/CLT---------------------------------------------
      
     if (dimsizes(nino34).eq.dimsizes(ppt&time).or.dimsizes(nino34).eq.dimsizes(clt&time)) then
       rclt  = escorc(nino34,clt(lat|:,lon|:,time|:))
       ;rclt = pattern_cor(sst_nino34(lat|:,lon|:,time|:),clt(lat|:,lon|:,time|:), 1.0, 0)
       copy_VarCoords(clt(0,:,:), rclt)
       rclt@long_name = "Correlation: Nino34-CLOUD"
       
       ;rpr = pattern_cor(sst_nino34(lat|:,lon|:,time|:),ppt(lat|:,lon|:,time|:), 1.0, 0)
       rpr  = escorc(nino34,ppt(lat|:,lon|:,time|:))                            
       copy_VarCoords(ppt(0,:,:), rpr)                                                                                                                
       rpr@long_name = "Correlation: Nino34-RAIN"       
     else
       print("Obs do not match")
       continue
     end if
      
      delete ([/ppt,clt,sst, wgt/])
      delete(nino34_ndj)
;========================================================================   
      
      scres = True    ; scres = spatial composite res
      scres@mpProjection = "WinkelTripel"
      scres@mpGeophysicalLineColor = "gray42"
      
     scres@mpPerimOn    = False
     scres@mpGridLatSpacingF =  90            ; change latitude  line spacing
     scres@mpGridLonSpacingF = 180.           ; change longitude line spacing
     scres@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     scres@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
     scres@mpFillOn = False
     scres@mpCenterLonF = 210.
     scres@mpOutlineOn = True  
     scres@mpDataSetName = "Earth..4"
     scres@gsnDraw      = False
     scres@gsnFrame     = False
     
     scres@cnLineLabelsOn = False
     scres@cnFillOn        = True
     scres@cnLinesOn       = False
     scres@lbLabelBarOn    = False
     scres@cnInfoLabelOn = False
     scres@gsnAddCyclic = True
     

     scres@gsnLeftStringOrthogonalPosF = -0.05
     scres@gsnLeftStringParallelPosF = .005
     scres@gsnRightStringOrthogonalPosF = -0.05
     scres@gsnRightStringParallelPosF = 0.96
     scres@gsnRightString = ""
     scres@gsnLeftString = ""
     scres@gsnLeftStringFontHeightF = 0.014
     scres@gsnCenterStringFontHeightF = 0.018
     scres@gsnRightStringFontHeightF = 0.014
     scres@gsnLeftString = syear(ee)+"-"+eyear(ee)
 
     scres@gsnRightString = ""
     scres@gsnCenterString = names(ee)  
      
      scres2 = scres
      scres3 = scres
           
      scres@mpOutlineDrawOrder = "PostDraw"
      scres@cnFillMode = "AreaFill"
           
      if (wks_type.eq."png") then
        ;scres3@cnLineThicknessF = 3.
        scres@mpGeophysicalLineThicknessF = 2.  
      else
        ;scres3@cnLineThicknessF = 1.25
        scres@mpGeophysicalLineThicknessF = 1.  
      end if
      
      scres2 = scres
      scres3 = scres
      
      if (ee.eq.0.and.n34sc_obs.eq.0) then
        scres2@gsnCenterString = scres@gsnCenterString+" / "+names_ppt(ee)
        scres3@gsnCenterString = scres@gsnCenterString+" / "+names_clt(ee)
      end if
      
      scres3@cnLevelSelectionMode = "ExplicitLevels"
      scres3@cnLevels = (/-30,-25,-20,-17.5, -15,-12.5,-10,-7.5,-5,-2.5,-1.25,0,1.25,2.5,5,7.5,10,12.5,15,17.5,20,25,30/)
      
      scres2@cnLevelSelectionMode = "ExplicitLevels"
      scres2@cnLevels = (/-7,-6,-5,-4.5,-4,-3.5,-3,-2.5,-2,-1.75,-1.5,-1.25,-1,-0.75,-.5,-.25,0,.25,.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.5,5,6,7/)

      if (n34sc_plot_flag.eq.0) then
        if (n34sc_obs.eq.0.and.ee.eq.0) then    ; for metrics table
           patcor_clt = new((/nsim,dimsizes(n34sc_clt&lat),dimsizes(n34sc_clt&lon)/),typeof(n34sc_clt))
           patcor_clt!1 = "lat"
           patcor_clt&lat = n34sc_clt&lat
           patcor_clt!2 = "lon"
           patcor_clt&lon = n34sc_clt&lon
           patcor_ppt = new((/nsim,dimsizes(n34sc_ppt&lat),dimsizes(n34sc_ppt&lon)/),typeof(n34sc_ppt))
           patcor_ppt!1 = "lat"
           patcor_ppt&lat = n34sc_ppt&lat
           patcor_ppt!2 = "lon"
           patcor_ppt&lon = n34sc_ppt&lon
           patcor_clt(ee,:,:) = (/ n34sc_clt(12,:,:) /)
           patcor_ppt(ee,:,:) = (/ n34sc_ppt(12,:,:) /)
        end if
        if (n34sc_obs.eq.0.and.ee.ge.1) then
           patcor_clt(ee,:,:) = (/ linint2(n34sc_clt&lon,n34sc_clt&lat,n34sc_clt(12,:,:),True,patcor_clt&lon,patcor_clt&lat,0) /)
           patcor_ppt(ee,:,:) = (/ linint2(n34sc_ppt&lon,n34sc_ppt&lat,n34sc_ppt(12,:,:),True,patcor_ppt&lon,patcor_ppt&lat,0) /)
         end if
         
        map_n34sc_jja0_ppt(ee) = gsn_csm_contour_map(wks_n34sc_ppt,n34sc_ppt(6,:,:),scres2)    ; 6 = JJA 0    
        map_n34sc_son0_ppt(ee) = gsn_csm_contour_map(wks_n34sc_ppt,n34sc_ppt(9,:,:),scres2)    ; 9 = SON 0      
        map_n34sc_djf1_ppt(ee) = gsn_csm_contour_map(wks_n34sc_ppt,n34sc_ppt(12,:,:),scres2)    ; 12 = DJF+1             
        map_n34sc_mam1_ppt(ee) = gsn_csm_contour_map(wks_n34sc_ppt,n34sc_ppt(15,:,:),scres2)    ; 15 = MAM+1      
        
        map_n34sc_jja0_clt(ee) = gsn_csm_contour_map(wks_n34sc_clt,n34sc_clt(6,:,:),scres3)    ; 6 = JJA 0    
        map_n34sc_son0_clt(ee) = gsn_csm_contour_map(wks_n34sc_clt,n34sc_clt(9,:,:),scres3)    ; 9 = SON 0      
        map_n34sc_djf1_clt(ee) = gsn_csm_contour_map(wks_n34sc_clt,n34sc_clt(12,:,:),scres3)    ; 12 = DJF+1             
        map_n34sc_mam1_clt(ee) = gsn_csm_contour_map(wks_n34sc_clt,n34sc_clt(15,:,:),scres3)    ; 15 = MAM+1   
       
        delete([/n34sc_sst,n34sc_clt,n34sc_ppt,scres/])
     end if
 ;--------------------------------------------------------------------------------------------              
      xyres2 = True
      xyres2@gsnDraw = False
      xyres2@gsnFrame = False
      xyres2@gsnYRefLine = 0.0
      xyres2@gsnYRefLineColor = "gray42"
      xyres2@gsnAboveYRefLineColor = "red"
      xyres2@gsnBelowYRefLineColor = "blue"
      if (wks_type.eq."png") then
        xyres2@xyLineThicknessF = 4.
      else
        xyres2@xyLineThicknessF = 2.0
      end if
      xyres2@tiYAxisString = ""
      if (nsim.le.5) then
        xyres2@tmXBLabelFontHeightF = 0.0125
        xyres2@tmYLLabelFontHeightF = 0.0125
        xyres2@gsnLeftStringFontHeightF = 0.017     
        xyres2@gsnRightStringFontHeightF = 0.013     
      else
         xyres2@tmXBLabelFontHeightF = 0.018
         xyres2@tmYLLabelFontHeightF = 0.018
         xyres2@gsnLeftStringFontHeightF = 0.024
         xyres2@gsnRightStringFontHeightF = 0.020     
       end if
       xyres2@vpXF = 0.05
       xyres2@vpHeightF = 0.15
       if (SCALE_TIMESERIES.eq."True") then
         xyres2@vpWidthF = 0.9*((nyr(ee)*1.)/nyr_max)
       else
         xyres2@vpWidthF = 0.9
      end if
      xyres2@gsnLeftString = ""     
      xyres2@gsnCenterString = ""
      xyres2@gsnRightString = ""
      
      xyres2@trXMinF = syear(ee)-.5
      xyres2@trXMaxF = eyear(ee)+0.5
      
      xyres2@xyLineColor = "gray60"
      xyres2@xyCurveDrawOrder = "PreDraw"
      xyres2@trYMinF = -3.
      xyres2@trYMaxF = 3.
      
      xyres = True
      if (ee.eq.0.and.n34sc_obs.eq.0) then
        xyres@gsnLeftString = names(ee)+" / "+names_ppt(ee)
      else  
        xyres@gsnLeftString = names(ee) 
      end if  
      xyres@trYMinF = -1.
      xyres@trYMaxF = 1.
      xyres@xyLineColor = "black"
      xyres@xyLineThicknessF = 2
      
      xyres3 = True                    
      if (ee.eq.0.and.n34sc_obs.eq.0) then
        xyres3@gsnLeftString = names(ee)+" / "+names_clt(ee)
      else  
        xyres3@gsnLeftString = names(ee) 
      end if                            
      xyres3@trYMinF = -16.
      xyres3@trYMaxF = 16.
      xyres3@xyLineColor = "black"
      xyres3@xyLineThicknessF = 2
            
      ihp     = 0                             ; low pass
      sigma   = 1.0                           ; Lanczos sigma
      
      nWgt    = 13                           ; loose 24 months each end                            
      fca     = 1./6.                       ; start freq
      wgt     = filwgts_lanczos (nWgt, ihp, fca, -999., sigma )
      
      clt_fil    = wgt_runave ( clt_aa_mon, wgt, 0 )      ; 20-100 day
      ppt_fil    = wgt_runave ( ppt_aa_mon, wgt, 0 )      ; 20-100 day
      sst_fil    = wgt_runave ( sst_aa_mon, wgt, 0 )      ; 20-100 day
          
      corr1 = escorc(sst_fil,ppt_fil)
      corr2 = escorc(sst_fil,clt_fil)
      corr3 = escorc(sst_fil, nino34)

      xyres@gsnRightString = decimalPlaces(corr1,2,True)
      xyres3@gsnRightString = decimalPlaces(corr2,2,True)     
      
      ;--------------------------------- lag correlation -------------------------------------------------
      maxlag=15                                     ; set lag
      totlag=2*maxlag                               ; set total lag
      
      ccr1_ppt = esccr(nino34,ppt_fil,maxlag)                  ; calc pos lag cross cor
      ccr2_ppt = esccr(ppt_fil,nino34,maxlag)                  ; calc neg lag cross cor

      ccr1_clt = esccr(nino34,clt_fil,maxlag)                  ; calc pos lag cross cor
      ccr2_clt = esccr(clt_fil,nino34,maxlag)                  ; calc neg lag cross cor
      ;************************************************
      ;  combine pos and neg into one time series;;;
      ;************************************************
      ccrtot_ppt = new( (/totlag/),float)               ; allocate memory 
      ccrtot_ppt(0:maxlag-1) = ccr2_ppt(0:maxlag-1:-1)      
      ccrtot_ppt(maxlag:)    = ccr1_ppt(0:maxlag-1)
      
      ccrtot_clt = new( (/totlag/),float)               ; allocate memory 
      ccrtot_clt(0:maxlag-1) = ccr2_clt(0:maxlag-1:-1)      
      ccrtot_clt(maxlag:)    = ccr1_clt(0:maxlag-1)
 
      x = ispan(-maxlag,maxlag,1)    

      lagres1                 = True
      lagres1@gsnDraw              = False
      lagres1@gsnFrame             = False                           ; make plot mods
      ;lagres1@tiMainString    = "SST vs. PPT"  ; title
      lagres1@tiXAxisString   = "SST leads        PPT leads"                        ; x-axis label
      lagres1@gsnLeftString = names(ee)
      lagres1@gsnRightString = names_ppt(ee)
      
      lagres2                 = True 
      lagres2@gsnDraw              = False
      lagres2@gsnFrame             = False                   ; plot correlation
      ;lagres2@tiMainString    = "SST vs. CLT"  ; title
      lagres2@tiXAxisString   = "SST leads        CLT leads"                        ; x-axis label
      lagres2@gsnLeftString = names(ee)
      lagres2@gsnRightString = names_ppt(ee)

      lagppt_xy(ee) = gsn_csm_xy(wks_lagppt,x,ccrtot_ppt,lagres1)                       ; plot correlation ppt
      lagclt_xy(ee) = gsn_csm_xy(wks_lagclt,x,ccrtot_clt,lagres2)                       ; plot correlation clt
      
      delete([/ccrtot_ppt, ccrtot_clt, maxlag, totlag, ccr1_ppt, ccr2_ppt, ccr1_clt, ccr2_clt, x/])
      ;---------------------------------------------------------------------------------------------------
      cttt = dtrend_msg(ispan(0,dimsizes(clt_aa_mon&time)-1,1),clt_aa_mon,False,True) 
      tttt = dtrend_msg(ispan(0,dimsizes(ppt_aa_mon&time)-1,1),ppt_aa_mon,False,True) 
      sttt = dtrend_msg_n(ispan(0,dimsizes(nino34)-1,1),nino34,False,True,0)
      sttt1 = dtrend_msg_n(ispan(0,dimsizes(sst_aa_mon)-1,1),sst_aa_mon,False,True,0)
      
      
      arr = new((/2,dimsizes(nino34)/),typeof(nino34))
      arr(0,:) = (/ nino34 /)
      arr(1,:) = (/  (ispan(0,dimsizes(nino34)-1,1)*sttt@slope)+sttt@y_intercept /)
     
      xy_mon(ee) = gsn_csm_xy2(wks_aa_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nino34)),nino34,ppt_fil,xyres2, xyres)
      xy_mon1(ee) = gsn_csm_xy2(wks_aa_mon1,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nino34)),nino34,clt_fil,xyres2, xyres3)
      
      delete([/sttt, sttt1, cttt, tttt, arr, ppt_aa_mon, clt_aa_mon, sst_aa_mon, nino34,corr1, corr2, xyres, xyres2, xyres3,\
      ppt_fil, clt_fil, sst_fil, wgt/])

;------------------------------Correlation Maps--------------------------------------------------------------------------------------------
     rescn = True
     rescn@gsnDraw              = False
     rescn@gsnFrame             = False
     rescn@gsnSpreadColors      = True     ; use full range of color map
     rescn@cnFillOn             = True     ; turn on color fill
     rescn@cnLinesOn            = False    ; turn of contour lines
     rescn@cnLineLabelsOn       = False    ; turn of contour line labels
     rescn@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     rescn@cnMinLevelValF       = -0.9               ; set min contour level
     rescn@cnMaxLevelValF       =  0.9               ; set max contour level
     rescn@cnLevelSpacingF      =  0.1               ; set contour spacing
     rescn@lbLabelBarOn         = False              ; turn off individual cb's
     rescn@mpCenterLonF         = 180.
     rescn@mpFillOn             = False              ; do not grey fill land
     
     rescn@mpProjection         = "Robinson"
     rescn@mpPerimOn            = False     
     rescn@mpGridAndLimbOn      = True 
     rescn@mpGridLatSpacingF    =  90                ; change latitude  line spacing
     rescn@mpGridLonSpacingF    = 360.               ; change longitude line spacing
                                ;rescn@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
     rescn@mpGridAndLimbOn      = True               ; turn on lat/lon lines
     rescn1 = rescn
     rescn2 = rescn
     if (ee.eq.0.and.n34sc_obs.eq.0) then
       rescn1@gsnLeftString = names(ee)+" / "+names_clt(ee)
     else  
       rescn1@gsnLeftString = names(ee) 
     end if
     if (ee.eq.0.and.n34sc_obs.eq.0) then
       rescn2@gsnLeftString = names(ee)+" / "+names_ppt(ee)
     else  
       rescn2@gsnLeftString = names(ee)   
     end if                           ; eliminate really small values
     
     rescn3 = True
     rescn3@gsnDraw              = False
     rescn3@gsnFrame             = False
     rescn3@gsnSpreadColors      = True     ; use full range of color map
     rescn3@cnFillOn             = True     ; turn on color fill
     rescn3@cnLinesOn            = False    ; turn of contour lines
     rescn3@cnLineLabelsOn       = False    ; turn of contour line labels
     rescn3@cnLevelSelectionMode = "ManualLevels"     ; set manual contour levels
     rescn3@cnMinLevelValF       = -1.4               ; set min contour level
     rescn3@cnMaxLevelValF       =  1.4               ; set max contour level
     rescn3@cnLevelSpacingF      =  0.1               ; set contour spacing
     rescn3@lbLabelBarOn         = False              ; turn off individual cb's
     rescn3@mpCenterLonF         = 180.
     rescn3@mpFillOn             = False              ; do not grey fill land
     
     rescn3@mpProjection         = "Robinson"
     rescn3@mpPerimOn            = False     
     rescn3@mpGridAndLimbOn      = True 
     rescn3@mpGridLatSpacingF    =  90                ; change latitude  line spacing
     rescn3@mpGridLonSpacingF    = 360.               ; change longitude line spacing
                                ;rescn@mpGridLineColor      = "transparent"      ; trick ncl into drawing perimeter
     rescn3@mpGridAndLimbOn      = True               ; turn on lat/lon lines rescn3@gsnLeftString = names(ee)
     rescn3@gsnLeftString = names(ee)   
     
     rclt = where(abs(rclt).lt.0.10, rclt@_FillValue, rclt)
     rpr = where(abs(rpr).lt.0.10, rpr@_FillValue, rpr)

     rescn@gsnRightString = "test1"
     plot_corr_cloud(ee) = gsn_csm_contour_map(wks_corr_cloud, rclt, rescn1)  ; create plot 
     
     rescn@gsnRightString = "test2"
     plot_corr_rain(ee) = gsn_csm_contour_map(wks_corr_rain, rpr, rescn2)  ; create plot

     rescn@gsnRightString = "test3"
     plot_sst_reg(ee) = gsn_csm_contour_map(wks_sst_reg, sst_nino34 , rescn3)  ; create plot
     
     delete ([/rclt,rpr, rescn, sst_nino34/])
;------------------------------------------------------------------------------------------------------------------------------------         
   end do
 ;--------------------------------------------Panelling------------------------------------------------------------------------------
   ncol = floattointeger(sqrt(nsim))
   nrow = (nsim/ncol)+mod(nsim,ncol) 
   
   panres = True
   panres@gsnMaximize = True
   panres@gsnPaperOrientation = "portrait"
   panres@gsnPanelYWhiteSpacePercent = 3.0
   if (nsim.le.10) then
     panres@txFontHeightF = 0.016
  else
    panres@txFontHeightF = 0.012
  end if
  
  if (SCALE_TIMESERIES.eq."True") then
    tt = ind(nyr.eq.nyr_max)
    panres@gsnPanelScalePlotIndex = tt(0)
    delete(tt)
  end if
  if (nsim.le.12) then
    lp = (/nsim,1/)
  else
    lp = (/nrow,ncol/)   ;(/nsim/2+1,nsim/8+1/)  
  end if

  if(n34sc_frame.eq.0) then
    panres@txString = "sst vs. precip "+area
    gsn_panel(wks_aa_mon,xy_mon,lp,panres)   
    delete(wks_aa_mon)

    panres@txString = "sst vs. cloud "+area
    gsn_panel(wks_aa_mon1,xy_mon1,lp,panres)
    delete(wks_aa_mon1)

    panres@txString = "sst vs. precip lag"+area
    gsn_panel(wks_lagppt,lagppt_xy,lp,panres)   
    delete(wks_lagppt)

    panres@txString = "sst vs. clt lag"+area
    gsn_panel(wks_lagclt,lagclt_xy,lp,panres)   
    delete(wks_lagclt)
  end if
  delete(area)
;-------------------------------------------------------------------------
   resP                     = panres                ; modify the panel plot
   resP@txString            = "Nino34 Correlation: lag=0"
   resP@gsnPanelLabelBar    = True                ; add common colorbar
   resP@lbLabelAutoStride   = True
   resP@gsnMaximize         = True                ; add common colorbar

   if(n34sc_frame.eq.0) then
     resP@txString = "Correlation Nino34 - Cloud Cover"
     gsn_panel(wks_corr_cloud,plot_corr_cloud,lp,resP)               ; now draw as one plot
     delete(wks_corr_cloud)
     
     resP@txString = "Correlation Nino34 - Precipitation"
     gsn_panel(wks_corr_rain,plot_corr_rain,lp,resP)               ; now draw as one plot
     delete(wks_corr_rain)

     resP@txString = "Nino3.4 - SST Regression"
     gsn_panel(wks_sst_reg,plot_sst_reg,lp,resP)               ; now draw as one plot
     delete(wks_sst_reg)
   end if

;-------------------------------------------------------------------------------- 
   panres2 = True
   panres2@gsnMaximize = True
   panres2@gsnPaperOrientation = "portrait"
   panres2@gsnPanelLabelBar = True
   panres2@lbLabelStride = 1
   panres2@pmLabelBarWidthF = 0.4
   panres2@pmLabelBarHeightF = 0.06
   panres2@lbLabelFontHeightF = 0.013
   panres2@txString = ""   
   if (nsim.le.5) then
     panres2@txFontHeightF = 0.024
   else
     panres2@txFontHeightF = 0.016
   end if
   panres2@gsnPanelXWhiteSpacePercent = 8.5
   panres2@gsnPanelBottom = 0.05
   
   if (n34sc_frame.eq.0) then
     panres2@txString = "nino3.4 CLT Spatial Composite (JJA~S~0~N~)"
     gsn_panel(wks_n34sc_clt,map_n34sc_jja0_clt,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 CLT Spatial Composite (SON~S~0~N~)"
     gsn_panel(wks_n34sc_clt,map_n34sc_son0_clt,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 CLT Spatial Composite (DJF~S~+1~N~)"
     gsn_panel(wks_n34sc_clt,map_n34sc_djf1_clt,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 CLT Spatial Composite (MAM~S~+1~N~)"
     gsn_panel(wks_n34sc_clt,map_n34sc_mam1_clt,(/nrow,ncol/),panres2)   
     delete(wks_n34sc_clt)
     
     panres2@txString = "nino3.4 PPT Spatial Composite (JJA~S~0~N~)"
     gsn_panel(wks_n34sc_ppt,map_n34sc_jja0_ppt,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 PPT Spatial Composite (SON~S~0~N~)"
     gsn_panel(wks_n34sc_ppt,map_n34sc_son0_ppt,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 PPT Spatial Composite (DJF~S~+1~N~)"
     gsn_panel(wks_n34sc_ppt,map_n34sc_djf1_ppt,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 PPT Spatial Composite (MAM~S~+1~N~)"
     gsn_panel(wks_n34sc_ppt,map_n34sc_mam1_ppt,(/nrow,ncol/),panres2)   
     delete(wks_n34sc_ppt)
     delete([/map_n34sc_djf1_clt,map_n34sc_jja0_clt,map_n34sc_son0_clt,map_n34sc_mam1_clt,\
     map_n34sc_djf1_ppt,map_n34sc_jja0_ppt,map_n34sc_son0_ppt,map_n34sc_mam1_ppt/])
  end if     
       

  OUTDIR = getenv("OUTDIR")
  
  if (n34sc_frame.eq.0) then
     if (wks_type.eq."png") then  
        system("mv "+OUTDIR+"nino34.spatialcomp.clt.000001.png "+OUTDIR+"nino34.spatialcomp.clt.jja0.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.clt.000002.png "+OUTDIR+"nino34.spatialcomp.clt.son0.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.clt.000003.png "+OUTDIR+"nino34.spatialcomp.clt.djf1.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.clt.000004.png "+OUTDIR+"nino34.spatialcomp.clt.mam1.png")
        
        system("mv "+OUTDIR+"nino34.spatialcomp.ppt.000001.png "+OUTDIR+"nino34.spatialcomp.ppt.jja0.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.ppt.000002.png "+OUTDIR+"nino34.spatialcomp.ppt.son0.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.ppt.000003.png "+OUTDIR+"nino34.spatialcomp.ppt.djf1.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.ppt.000004.png "+OUTDIR+"nino34.spatialcomp.ppt.mam1.png")
        
     else
        system("psplit "+OUTDIR+"nino34.spatialcomp.clt.ps "+OUTDIR+"pict")
        system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nino34.spatialcomp.clt.jja0.ps")
        system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nino34.spatialcomp.clt.son0.ps")
        system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nino34.spatialcomp.clt.djf1.ps")
        system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nino34.spatialcomp.clt.mam1.ps")
        system("rm "+OUTDIR+"nino34.spatialcomp.clt.ps")

        system("psplit "+OUTDIR+"nino34.spatialcomp.ppt.ps "+OUTDIR+"pict")
        system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nino34.spatialcomp.ppt.jja0.ps")
        system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nino34.spatialcomp.ppt.son0.ps")
        system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nino34.spatialcomp.ppt.djf1.ps")
        system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nino34.spatialcomp.ppt.mam1.ps")
        system("rm "+OUTDIR+"nino34.spatialcomp.ppt.ps")
     end if
   end if
  
;---------------------------------------------------------------------------------
  delete([/nrow,ncol,lp,xy_mon,xy_mon1, panres,rad, pi, resP, plot_corr_rain, plot_corr_cloud/])

  print("Finished: test4.ncl")
end
