; Calculates a variety of oceanic indices, as well as hovmollers, spectra, 
; monthly standard deviations, running standard deviations, and spatial
; composites based on the nino3.4 index. 
;
; Variables used: TS, PSL, and TREFHT (tas)
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: sst.indices.ncl")
  
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
  delete(na)
  nyr = eyear-syear+1
  nyr_max = max(nyr)
;----------- nino3.4 spatial composite coding-------------  
  nsim_psl = numAsciiRow("namelist_byvar/namelist_psl")
  na_psl = asciiread("namelist_byvar/namelist_psl",(/nsim_psl/),"string")
  names_psl = new(nsim_psl,"string")
  paths_psl = new(nsim_psl,"string")
  syear_psl = new(nsim_psl,"integer",-999)
  eyear_psl = new(nsim_psl,"integer",-999)

  do gg = 0,nsim_psl-1
     names_psl(gg) = str_strip(str_get_field(na_psl(gg),1,delim))
     paths_psl(gg) = str_strip(str_get_field(na_psl(gg),2,delim))
     syear_psl(gg) = stringtointeger(str_strip(str_get_field(na_psl(gg),3,delim)))
     eyear_psl(gg) = stringtointeger(str_strip(str_get_field(na_psl(gg),4,delim)))
  end do
  delete(na_psl)
  nyr_psl = eyear_psl-syear_psl+1
    
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

  if (isfilepresent2("obs_ts")) then
     if (nsim.eq.nsim_psl.and.nsim.eq.nsim_trefht.and.nsim_psl.eq.nsim_trefht) then
        ; print("Total number of simulations/obs are equal between psl/trefht/ts namelists, proceeding")
        n34sc_obs = 0
     else
        print("Total number of simulations/obs between psl/trefht/ts do not match, skipping observational nino3.4 spatial composite.")
        n34sc_obs = 1
        ; delete([/nsim_psl,na_psl,names_psl,paths_psl,syear_psl,eyear_psl,nyr_psl,nsim_trefht,na_trefht,names_trefht,paths_trefht,syear_trefht,eyear_trefht,nyr_trefht/])
     end if  
  else
     n34sc_obs = 1
  end if
  if (n34sc_obs.eq.0) then
     if (syear(0).eq.syear_psl(0).and.syear(0).eq.syear_trefht(0).and.syear_psl(0).eq.syear_trefht(0)) then
        if (eyear(0).eq.eyear_psl(0).and.eyear(0).eq.eyear_trefht(0).and.eyear_psl(0).eq.eyear_trefht(0)) then
;           print("Obs years match")
           n34sc_obs = 0
        else
           print("Observational end years do not match, skipping nino3.4 spatial composite")
           n34sc_obs = 1
        end if
     else
        print("Observational start years do not match, skipping nino3.4 spatial composite")
        n34sc_obs = 1
     end if
  end if
;-------------------------------------------------------------------------------------------------
        
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_n34 = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.timeseries")  
  wks_n4 = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino4.timeseries")
  wks_n3 = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino3.timeseries")
  wks_n12 = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino12.timeseries")
  wks_tna = gsn_open_wks(wks_type,getenv("OUTDIR")+"tna.timeseries")
  wks_tsa = gsn_open_wks(wks_type,getenv("OUTDIR")+"tsa.timeseries")
  wks_tio = gsn_open_wks(wks_type,getenv("OUTDIR")+"tio.timeseries")

  wks_n34_tlon_hi = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.hov.elnino")
  wks_n34_tlon_lo = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.hov.lanina")

  wks_n34_p = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.powspec")
  
  wks_n34_rst = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.runstddev")

  wks_n34_mst = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.monstddev")
  
  wks_n34sc = gsn_open_wks(wks_type,getenv("OUTDIR")+"nino34.spatialcomp")
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_n34,"ncl_default")    
     gsn_define_colormap(wks_n4,"ncl_default")  
     gsn_define_colormap(wks_n3,"ncl_default")  
     gsn_define_colormap(wks_n12,"ncl_default")   
     gsn_define_colormap(wks_tna,"ncl_default")
     gsn_define_colormap(wks_tsa,"ncl_default")
     gsn_define_colormap(wks_tio,"ncl_default")  
     gsn_merge_colormaps(wks_n34_tlon_hi,"BlueDarkRed18",(/"gray30","gray50","gray70"/))
     gsn_merge_colormaps(wks_n34_tlon_lo,"BlueDarkRed18",(/"gray30","gray50","gray70"/))
     gsn_define_colormap(wks_n34_p,"cb_9step")
     gsn_define_colormap(wks_n34_rst,"ncl_default")   
     gsn_define_colormap(wks_n34_mst,"ncl_default")   
     gsn_define_colormap(wks_n34sc,"ncl_default")   
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_n34,"ncl_default")    
     gsn_define_colormap(wks_n4,"ncl_default")  
     gsn_define_colormap(wks_n3,"ncl_default")  
     gsn_define_colormap(wks_n12,"ncl_default")   
     gsn_define_colormap(wks_tna,"ncl_default")
     gsn_define_colormap(wks_tsa,"ncl_default")
     gsn_define_colormap(wks_tio,"ncl_default")  
     gsn_merge_colormaps(wks_n34_tlon_hi,"BlueDarkRed18",(/"gray30","gray50","gray70"/))
     gsn_merge_colormaps(wks_n34_tlon_lo,"BlueDarkRed18",(/"gray30","gray50","gray70"/))
     gsn_define_colormap(wks_n34_p,"cb_9step")
     gsn_define_colormap(wks_n34_rst,"ncl_default")   
     gsn_define_colormap(wks_n34_mst,"ncl_default")   
     gsn_define_colormap(wks_n34sc,"BlueDarkRed18")   
  end if

  xyn34 = new(nsim,"graphic")  
  xyn4  = new(nsim,"graphic")  
  xyn3  = new(nsim,"graphic")  
  xyn12 = new(nsim,"graphic")  
  xytna = new(nsim,"graphic")  
  xytsa = new(nsim,"graphic") 
  xytio = new(nsim,"graphic")  
  xyiod = new(nsim,"graphic")  
  
  plot_n34hi = new(nsim,"graphic")  
  plot_n34lo = new(nsim,"graphic")  
  
  map_n34sc_jja0 = new(nsim,"graphic")  
  map_n34sc_son0 = new(nsim,"graphic")  
  map_n34sc_djf1 = new(nsim,"graphic")  
  map_n34sc_mam1 = new(nsim,"graphic")  

  xyn34_rst = new(nsim,"graphic")  
  xyn34_mst = new(nsim,"graphic")  

  pspec = new(nsim,"graphic")
  if (isfilepresent2("obs_ts")) then
     pspec_obs = new(nsim,"graphic")
  end if

  wgt = (/1.,2.,1./)   
  wgt = wgt/sum(wgt)
  pi=4.*atan(1.0)
  rad=(pi/180.)
  
  n34sc_frame = 1     ; n34sc_frame = flag to create regressions .ps/.png files. Created/used instead of n34sc_plot_flag
                      ; so that if sst regressions are not created for the last simulation listed that .ps/png files are created
  do ee = 0,nsim-1
;     print(paths(ee)+" ee = "+ee)
     sst = data_read_in(paths(ee),"TS",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     
     if (isatt(sst,"is_all_missing")) then
        delete(sst)
        continue
     end if
       
     if (isfilepresent2("obs_ts").and.ee.eq.0) then   ; mask out land in array going into EOF functions
        n34sc_flag = 1     ; this flag MAY be overruled by n34sc_obs
     else
        n34sc_flag = 0     ; flag set to trigger nino3.4 composite
        
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
     nino34@long_name = "nino3.4 timeseries"

     llats = -5.    ; nino3
     llatn = 5.
     llonw = 210.
     llone = 270.
     nino3 = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino3@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,nino3)
     nino3@units = sst@units
     nino3@long_name = "nino3 timeseries"
     
     llats = -5.    ; nino4
     llatn = 5.
     llonw = 160.
     llone = 210.
     nino4 = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino4@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,nino4)
     nino4@units = sst@units
     nino4@long_name = "nino4 timeseries"
     
     llats = -10.    ; nino1+2
     llatn = 0.
     llonw = 270.
     llone = 280.
     nino12 = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino12@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,nino12)
     nino12@units = sst@units
     nino12@long_name = "nino1+2 timeseries"

     ssttemp = lonFlip(sst)
     llats = -20.    ;  Tropical Southern Atlantic Index
     llatn = 0.
     llonw = -30.
     llone = 10.
     tsa = wgt_areaave(ssttemp(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     tsa@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,tsa)
     tsa@units = sst@units
     tsa@long_name = "Tropical Southern Atlantic SST timeseries"
     delete(ssttemp)
     
     llats = 5.5    ;  Tropical Northern Atlantic Index
     llatn = 23.5
     llonw = 302.5
     llone = 345.
     tna = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     tna@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,tna)
     tna@units = sst@units
     tna@long_name = "Tropical Northern Atlantic SST timeseries"

     llats = -15.    ;  Indian Ocean SST index
     llatn = 15.
     llonw = 40.
     llone = 110.
     tio = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     tio@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,tio)
     tio@units = sst@units
     tio@long_name = "Tropical Indian Ocean SST timeseries"
     
                     ; Indian Ocean Dipole Index http://www.bom.gov.au/climate/IOD/about_IOD.shtml
     iod = wgt_areaave(sst(:,{-10.:10.},{50.:70.}),coswgt({-10.:10.}),1.0,0) - wgt_areaave(sst(:,{-10.:0.},{90.:110.}),coswgt({-10.:0.}),1.0,0)  
     iod@area = "(-10:10N, 50:70E) - (-10:0N, 90:110E)"
     copy_VarCoords(nino34,iod)
     iod@units = sst@units
     iod@long_name = "Indian Ocean Dipole Index"    
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
           date = cd_calendar(nino34&time,-1)
           date@long_name = "current date (YYYYMM)"
           delete(date@calendar)
           date!0 = "time"
           date&time = nino34&time
           z->date = date
           delete(date)
        else
           z = addfile(fn,"w")
        end if
        z->nino34 = nino34
        z->nino12 = nino12
        z->nino3  = nino3
        z->nino4  = nino4
        z->north_tropical_atlantic = tna
        z->south_tropical_atlantic = tsa
        z->tropical_indian_ocean = tio
        z->indian_ocean_dipole = iod
        delete([/modname,fn/])
     end if
;---------------------------------------------------------------------------------------------
     nino34T = wgt_runave(nino34,wgt,1)     ; for use in ENSO composites / hovmuellers / running standard deviations
     nino34_ndj = nino34T(11:dimsizes(nino34T)-13:12)   ; cannot count last 1yr as spatial composite uses +1yrs data betond NDJ..  
     nino34_ndj!0 = "time"
     nino34_ndj&time = ispan(syear(ee),eyear(ee)-1,1)    
     nino34_ndj = dtrend_msg(ispan(0,dimsizes(nino34_ndj&time)-1,1),nino34_ndj,True,False)
     nino34_ndj = dim_standardize(nino34_ndj,0)

     sst = (/ dtrend_msg_n(ispan(0,nyr(ee)*12-1,1),sst,False,False,0) /) ; detrend the sst array
     
     sstr = sst(:,{-3:3},{120:280})                        ; ENSO hovmuellers based on NDJ nino34
     ; delete(sst)
     finsst_hi = sstr(:60,0,:)    ; for Jan-2 -> Jan+3
     finsst_hi!0 = "time"
     finsst_hi&time = ispan(0,60,1)
     finsst_hi = 0.
     finsst_lo = finsst_hi
     finsst_mid = finsst_hi
     cntr_hi = 0
     cntr_lo = 0
     cntr_mid = 0
     cntr_lo@_FillValue = default_fillvalue(typeof(cntr_lo)) 
     cntr_mid@_FillValue = default_fillvalue(typeof(cntr_mid)) 
     cntr_hi@_FillValue = default_fillvalue(typeof(cntr_hi)) 
        
     mocntr = 24   ; note: if this is set at 24 gg should start at 2
     do gg = 2,dimsizes(nino34_ndj)-3   ; remember that Dec is month 11. End @ -3 because we need to grab + 3 yrs and 1 month from there (nino34_ndj already ends at eyear-1)
        if (.not.ismissing(nino34_ndj(gg))) then    ; note that finsst_* indices 24:52 (Jan+0 -> May +2) are all that is shown in the hovmoller plots 
           if (nino34_ndj(gg).ge.1.) then
              finsst_hi = (/ finsst_hi+dim_avg_n(sstr(mocntr-24:mocntr+36,:,:),1) /)    ; nino34_ndj value is at sstr index mocntr+11     
               cntr_hi = cntr_hi+1
           end if
           if (nino34_ndj(gg).ge.-0.5.and.nino34_ndj(gg).le.0.5) then
              finsst_mid = (/ finsst_mid+dim_avg_n(sstr(mocntr-24:mocntr+36,:,:),1) /)          
              cntr_mid = cntr_mid+1
           end if
           if (nino34_ndj(gg).le.-1.) then
              finsst_lo = (/ finsst_lo+dim_avg_n(sstr(mocntr-24:mocntr+36,:,:),1) /)            
              cntr_lo = cntr_lo+1
           end if
        end if
        mocntr = mocntr+12
     end do
     delete([/sstr,mocntr/])

     cntr_hi  = where(cntr_hi.eq.0, cntr_hi@_FillValue, cntr_hi)
     cntr_mid = where(cntr_mid.eq.0,cntr_mid@_FillValue,cntr_mid)
     cntr_lo  = where(cntr_lo.eq.0, cntr_lo@_FillValue, cntr_lo)
     finsst_hi  = (/ finsst_hi/cntr_hi /)
     finsst_mid = (/ finsst_mid/cntr_mid /)
     finsst_lo  = (/ finsst_lo/cntr_lo /)
     delete([/coswgt/])
     
     if (OUTPUT_DATA.eq."True") then
        hov_hi = (/ finsst_hi(24:52,:) /)  ; 24:52 runs from Jan+0->May+2 and matches range shown in plot
        time_mon1 = ispan(0,28,1)
        time_mon1@units = "months since 0000-01-01 00:00:00"
        time_mon1@long_name = "Time"
        time_mon1@standard_name = "time"
        time_mon1@calendar = "standard"
        time_mon1!0 = "time_mon1"
        time_mon1&time_mon1 = time_mon1
        hov_hi!0 = "time_mon1"
        hov_hi&time_mon1 = time_mon1
        longitude = finsst_hi&lon
        longitude@standard_name = "longitude"
        hov_hi!1 = "longitude"
        hov_hi&longitude = longitude
        delete([/time_mon1,longitude/])
        hov_lo = (/ finsst_lo(24:52,:) /)  ; 24:52 runs from Jan+0->May+2 and matches range shown in plot
        copy_VarCoords(hov_hi,hov_lo)
        hov_hi@number_of_events = cntr_hi
        hov_lo@number_of_events = cntr_lo
        hov_hi@units = "C"
        hov_lo@units = "C"
        z->n34hov_elnino = hov_hi
        z->n34hov_lanina = hov_lo
        delete([/hov_hi,hov_lo/])
     end if     
;- - - - - -nino3.4 spatial composite section- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     n34sc_plot_flag = 0   
     if (ee.eq.0.and.n34sc_obs.eq.1) then    ; first time thru, and obs = False
        n34sc_plot_flag = 1                  ; do not plot
        if (.not.isfilepresent2("obs_ts")) then  ; however, if there are no TS obs whatsoever, the 1st iteration of the
           n34sc_plot_flag = 0                  ; ee do loop will be the first model, so turn it on.
        end if
     end if
     
     if (n34sc_plot_flag.eq.0) then        
        psl_offset    = 0
        trefht_offset = 0
        if (ee.eq.0.and.n34sc_obs.eq.0) then
        else
           if (isfilepresent2("obs_ts")) then
              if (.not.isfilepresent2("obs_psl")) then
                 psl_offset = -1
              end if
              if (.not.isfilepresent2("obs_trefht")) then
                 trefht_offset = -1
              end if
           else
              if (isfilepresent2("obs_psl")) then
                 psl_offset = 1
              end if
              if (isfilepresent2("obs_trefht")) then
                 trefht_offset = 1
              end if
           end if
        end if
        
        ; print("Data to be read in: "+paths_trefht(ee+trefht_offset)+" from "+syear_trefht(ee+trefht_offset)+":"+eyear_trefht(ee+trefht_offset))
        tas = data_read_in(paths_trefht(ee+trefht_offset),"TREFHT",syear_trefht(ee+trefht_offset),eyear_trefht(ee+trefht_offset))
        ; print("Data to be read in: "+paths_psl(ee+psl_offset)+" from "+syear_psl(ee+psl_offset)+":"+eyear_psl(ee+psl_offset))
        psl = data_read_in(paths_psl(ee+psl_offset),"PSL",syear_psl(ee+psl_offset),eyear_psl(ee+psl_offset))

        TIME = sst&time
        yyyymm = cd_calendar(sst&time,-1)  ; convert tas, ts, and sst from CF-conforming time to YYYYMM for coding below
        delete(sst&time)
        sst&time = yyyymm 
        delete(yyyymm)

        yyyymm = cd_calendar(tas&time,-1)  
        delete(tas&time)
        tas&time = yyyymm 
        delete(yyyymm)

        yyyymm = cd_calendar(psl&time,-1) 
        delete(psl&time)
        psl&time = yyyymm 
        delete(yyyymm)

        if (isatt(tas,"is_all_missing").or.isatt(psl,"is_all_missing")) then
           n34sc_plot_flag = 1
           delete([/tas,psl/])
        end if
        
        if (nyr(ee).lt.15) then   ; 15+ years needed for composites
           n34sc_plot_flag = 1
        end if

        if (n34sc_plot_flag.eq.0) then     ; only continue if all 3 fields are present
         ;  d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")     ; mask ocean for TAS array
         ;  basemap = d->LSMASK                                              ; This is now done right before plotting.
         ;  lsm = landsea_mask(basemap,tas&lat,tas&lon)                      ; so that the entire TAS array is used
         ;  tas = mask(tas,conform(tas,lsm,(/1,2/)).eq.0,False)              ; in the nino3.4 pattern correlations
         ;  delete([/lsm,basemap/])                                          ; (Even if the land portion of TAS is the
         ;  delete(d)                                                        ; only portion plotted as SST shown over oceans.)

           if (OPT_CLIMO.eq."Full") then
              tas = rmMonAnnCycTLL(tas)
              psl = rmMonAnnCycTLL(psl)           
           else
              check_custom_climo(names_trefht(ee+trefht_offset),syear_trefht(ee+trefht_offset),eyear_trefht(ee+trefht_offset),CLIMO_SYEAR,CLIMO_EYEAR)
              climo = clmMonTLL(tas({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              tas   = calcMonAnomTLL(tas,climo) 
              delete(climo)
           
              check_custom_climo(names_psl(ee+psl_offset),syear_psl(ee+psl_offset),eyear_psl(ee+psl_offset),CLIMO_SYEAR,CLIMO_EYEAR)
              climo = clmMonTLL(psl({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              psl   = calcMonAnomTLL(psl,climo) 
              delete(climo)
           end if
           tas = (/ dtrend_msg_n(ispan(0,dimsizes(tas&time)-1,1),tas,False,False,0) /)     ; sst detrended up above
           psl = (/ dtrend_msg_n(ispan(0,dimsizes(psl&time)-1,1),psl,False,False,0) /)  
        
           ta = dim_avg_n(sst(:1,:,:),0)
           sst = runave_n_Wrap(sst,3,0,0)
           sst(0,:,:) = (/ ta /)
           delete(ta)
        
           ta = dim_avg_n(psl(:1,:,:),0)
           psl = runave_n_Wrap(psl,3,0,0)
           psl(0,:,:) = (/ ta /)
           delete(ta)
        
           ta = dim_avg_n(tas(:1,:,:),0)
           tas = runave_n_Wrap(tas,3,0,0)
           tas(0,:,:) = (/ ta /)
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
        
           dimS = dimsizes(psl&time)    ; change time from YYYYMM->YYYY.frac
           tmin = psl&time(0)/100
           tmax = psl&time(dimS-1)/100
           delete(psl&time)
           psl&time = fspan(tmin*1.,(tmax*1.)+(11/12.),dimS)
           dimS = dimsizes(tas&time)
           tmin = tas&time(0)/100
           tmax = tas&time(dimS-1)/100
           delete(tas&time)
           tas&time = fspan(tmin*1.,(tmax*1.)+(11/12.),dimS)
           dimS = dimsizes(sst&time)
           tmin = sst&time(0)/100
           tmax = sst&time(dimS-1)/100
           delete(sst&time)
           sst&time = fspan(tmin*1.,(tmax*1.)+(11/12.),dimS)
           delete([/dimS,tmin,tmax/])
           ; print(sst&time)
        
           sc_tas_hi = tas(:23,:,:)
           sc_tas_lo = tas(:23,:,:)    
           sc_sst_hi = sst(:23,:,:)
           sc_sst_lo = sst(:23,:,:)  
           sc_psl_hi = psl(:23,:,:)
           sc_psl_lo = psl(:23,:,:)  
           
           sc_tas_hi = sc_tas_hi@_FillValue
           sc_tas_lo = sc_tas_lo@_FillValue
           sc_sst_hi = sc_sst_hi@_FillValue
           sc_sst_lo = sc_sst_lo@_FillValue
           sc_psl_hi = sc_psl_hi@_FillValue
           sc_psl_lo = sc_psl_lo@_FillValue
           
           if (dimsizes(highyr).le.1) then
              print("For "+names(ee)+", 1 or less (normalized) nino3.4 value greater than one standard deviation found, setting nino3.4 spatial composites to missing")  ; sc_*_hi arrays left to _FillValue
              sc_tas_hi@is_all_missing = True
           else           
              do gg = 0,23
                 tt = gg/12.
                 sc_psl_hi(gg,:,:) = (/ dim_avg_n(psl({highyr+tt},:,:),0) /)
                 sc_sst_hi(gg,:,:) = (/ dim_avg_n(sst({highyr+tt},:,:),0) /)    
                 sc_tas_hi(gg,:,:) = (/ dim_avg_n(tas({highyr+tt},:,:),0) /)
              end do
              delete([/tt,highyr/])
           end if
           if (dimsizes(lowyr).le.1) then
              print("For "+names(ee)+", 1 or less (normalized) nino3.4 value less than -1 standard deviation found, setting nino3.4 spatial composites to missing")   ; sc_*_lo arrays left to _FillValue
           else           
              do gg = 0,23
                 tt = gg/12.
                 sc_psl_lo(gg,:,:) = (/ dim_avg_n(psl({lowyr+tt},:,:),0) /)
                 sc_sst_lo(gg,:,:) = (/ dim_avg_n(sst({lowyr+tt},:,:),0) /)     
                 sc_tas_lo(gg,:,:) = (/ dim_avg_n(tas({lowyr+tt},:,:),0) /)
              end do
              delete([/tt,tas,psl,lowyr/])
           end if
           
           n34sc_psl = sc_psl_hi
           n34sc_psl = (/ sc_psl_hi - sc_psl_lo /)
           n34sc_sst = sc_sst_hi
           n34sc_sst = (/ sc_sst_hi - sc_sst_lo /)
           n34sc_tas = sc_tas_hi
           n34sc_tas = (/ sc_tas_hi - sc_tas_lo /)
           delete([/sc_psl_hi,sc_psl_lo,sc_sst_hi,sc_sst_lo,sc_tas_hi,sc_tas_lo/])
           delete(sst&time)
           sst&time = TIME    
           delete(TIME)       
                      
           if (OUTPUT_DATA.eq."True") then
              n34sc_sst&lat@standard_name = "latitude"
              n34sc_sst&lon@standard_name = "longitude"
              z->nino34spacomp_sst_jja0 = n34sc_sst(6,:,:)    
              z->nino34spacomp_sst_son0 = n34sc_sst(9,:,:)
              z->nino34spacomp_sst_djf1 = n34sc_sst(12,:,:)
              z->nino34spacomp_sst_mam1 = n34sc_sst(15,:,:)
              if (isfilepresent2("obs_ts").and.ee.eq.0) then   ; will not write out PSL/TAS  to OBS SST files due to different
              else                                             ; variable going into SST OBS file..
                 z->nino34spacomp_psl_jja0 = n34sc_psl(6,:,:)    
                 z->nino34spacomp_psl_son0 = n34sc_psl(9,:,:)
                 z->nino34spacomp_psl_djf1 = n34sc_psl(12,:,:)
                 z->nino34spacomp_psl_mam1 = n34sc_psl(15,:,:)
                 z->nino34spacomp_tas_jja0 = n34sc_tas(6,:,:)    
                 z->nino34spacomp_tas_son0 = n34sc_tas(9,:,:)
                 z->nino34spacomp_tas_djf1 = n34sc_tas(12,:,:)
                 z->nino34spacomp_tas_mam1 = n34sc_tas(15,:,:)
              end if
           end if
        end if
     end if
     if (isvar("TIME")) then
        delete(TIME)
     end if
     delete([/sst,nino34_ndj/])

     if (n34sc_frame.eq.1.and.n34sc_plot_flag.eq.0) then    ; n34sc_frame = flag to create composites .ps/.png files 
        n34sc_frame = 0
     end if   
;-----------------------------------------------------------------------------------------     
     if (nyr(ee).ge.35) then    ; need a minimum number of years to compute running nino3.4 standard deviations
        nino34T = dtrend_msg(ispan(0,dimsizes(nino34T)-1,1),nino34T,True,False)
        nino34T!0 = "time"
        nino34T&time = nino34&time
        sd_run = nino34T
        sd_run = sd_run@_FillValue
        sd_run@units = nino34@units
        sd_run@long_name = "nino3.4 30yr running standard deviation"
        do gg = 179,dimsizes(nino34T)-181
           sd_run(gg) = (/ dim_stddev(nino34T(gg-179:gg+180)) /)
        end do
        if (OUTPUT_DATA.eq."True") then
           z->nino34_runstddev = sd_run
        end if
     end if
     delete(nino34T)
;-----------------------------------------------------------------------------------------
     iopt = 0                ; nino3.4 power spectra
     jave = (7*nyr(ee))/100
     val1 = .95
     val2 = .99
     pct = 0.1 
     spectra_mvf = False        ; missing value flag for nino3.4
     ;if (any(ismissing(nino34)).or.isatt(n34sc_tas,"is_all_missing")) then     ; n34sc_tas@is_all_missing results from there being no standardized nino34 values above 1 std dev, and that  <-- modified back to CVDP3.3.0 version, because got a error message      
     if (any(ismissing(nino34))) then
        print("Missing data detected for "+names(ee)+", power spectra function does not allow missing data, not creating spectra in sst.indices.ncl")  ; likely occurred due to climatological SSTs being present. Climatological SSTs may cause specx_ci to fail.
        spectra_mvf = True
        if (isfilepresent2("obs_ts").and.ee.eq.0) then
           spectra_mvf_obs = True    ; missing value flag for obs nino3.4
        end if
     else
        if (isfilepresent2("obs_ts").and.ee.eq.0) then
           spectra_mvf_obs = False   ; missing value flag for obs nino3.4
        end if
        nino34_dt = dtrend_msg(ispan(0,dimsizes(nino34)-1,1),nino34,True,False) 
    
        sdof = specx_anal(nino34_dt,iopt,jave,pct)
        splt1 = specx_ci(sdof,val1,val2)
        if (OUTPUT_DATA.eq."True") then
           splt1!0 = "ncurves"
           splt1&ncurves = ispan(0,3,1)
           splt1!1 = "frequency"
           splt1&frequency = sdof@frq
           splt1@units_info = "df refers to frequency interval"
           splt1@units = "C^2/df"
           splt1@info = "(0,:)=spectrum,(1,:)=Markov red noise spectrum, (2,:)="+val1+"% confidence bound for Markhov, (3,:)="+val2+"% confidence bound for Markhov"
           z->nino34_spectra = splt1
        end if
        if (isfilepresent2("obs_ts").and.ee.eq.0) then
           sdof_obs = sdof
        end if
        delete([/nino34_dt,iopt,jave,pct/])
     end if
;------------------------------------------------------------------------------------------
     nino34_dt = dtrend_msg(ispan(0,dimsizes(nino34&time)-1,1),nino34,True,False)
     nino34_mon_sd = new(12,typeof(nino34))
     
     do hh = 0,11
        nino34_mon_sd(hh) = (/ dim_stddev(nino34_dt(hh::12)) /)
     end do
     nino34_mon_sd@units = "C"
     delete(nino34_dt)    
     if (OUTPUT_DATA.eq."True") then
        time_mon2 = ispan(0,11,1)
        time_mon2@units = "months since 0000-01-01 00:00:00"
        time_mon2@long_name = "Time"
        time_mon2@standard_name = "time"
        time_mon2@calendar = "standard"
        time_mon2!0 = "time_mon2"
        time_mon2&time_mon2 = time_mon2
        nino34_mon_sd!0 = "time_mon2"
        nino34_mon_sd&time_mon2 = time_mon2
        z->nino34_monthly_stddev = nino34_mon_sd
        delete(time_mon2)
     end if
     if (isvar("z")) then
        delete(z)
     end if
;==========================================================================================
     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnRightString = ""
     xyres@gsnLeftString = ""
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"
     xyres@xyLineColor = "gray62"
     if (wks_type.eq."png") then
        xyres@xyLineThicknessF = .75  
     else
        xyres@xyLineThicknessF = .5  
     end if   
     xyres@tiYAxisString = ""
     if (nsim.le.5) then
        xyres@tmXBLabelFontHeightF = 0.0125
        xyres@tmYLLabelFontHeightF = 0.0125
        xyres@gsnLeftStringFontHeightF = 0.017
        xyres@gsnCenterStringFontHeightF = 0.017
        xyres@gsnRightStringFontHeightF = 0.013         
     else
        xyres@tmXBLabelFontHeightF = 0.018
        xyres@tmYLLabelFontHeightF = 0.018
        xyres@gsnLeftStringFontHeightF = 0.024
        xyres@gsnCenterStringFontHeightF = 0.024
        xyres@gsnRightStringFontHeightF = 0.020    
     end if
     xyres@vpXF = 0.05
     xyres@vpHeightF = 0.3
     if (SCALE_TIMESERIES.eq."True") then
        xyres@vpWidthF = 0.9*((nyr(ee)*1.)/nyr_max)
     else
        xyres@vpWidthF = 0.9
     end if    
     xyres@gsnCenterString = ""
     
     xyres@trXMinF = syear(ee)-.5
     xyres@trXMaxF = eyear(ee)+1.5
     xyres@tiMainOn = False
     
     xyres2 = xyres
     xyres2@vpHeightF = 0.15
     xyres2@gsnXYBarChart = False

     xyres2@xyLineColor = "royalblue"
     xyres2@trYMinF = 0.3    ; hard wire YMinF and YMaxF for running stddev plots
     xyres2@trYMaxF = 1.8       
     if (wks_type.eq."png") then
        xyres2@xyLineThicknessF = 3.5
     else
        xyres2@xyLineThicknessF = 1.75
     end if
     delete(xyres2@gsnYRefLine)
     xyres2@gsnYRefLine = (/.6,0.9,1.2,1.5/)
     xyres2@gsnYRefLineColor = "gray85"
     
     xyres3 = xyres            ; resource list for monthly nino3.4 standard deviations
     xyres3@trXMinF = 0.5
     xyres3@trXMaxF = 12.5
     xyres3@vpWidthF = 0.65
     xyres3@vpHeightF = 0.35
     xyres3@trYMinF = 0.2
     xyres3@trYMaxF = 2.0
     xyres3@gsnAboveYRefLineColor = "gray50"
     xyres3@xyLineColor = "black"
     if (wks_type.eq."png") then
        xyres3@xyLineThicknessF = 3.5
     else
        xyres3@xyLineThicknessF = 1.75
     end if
     xyres3@gsnXYBarChart = True
     xyres3@gsnXYBarChartBarWidth = 0.75
     xyres3@tmXBMode    = "Explicit"        ; explicit labels
     xyres3@tmXBValues     = ispan(1,12,1)
     xyres3@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
     xyres3@tmXTOn = False
     xyres3@gsnLeftStringOrthogonalPosF = 0.025
     xyres3@gsnCenterStringOrthogonalPosF = xyres3@gsnLeftStringOrthogonalPosF
     xyres3@gsnRightStringOrthogonalPosF = xyres3@gsnLeftStringOrthogonalPosF
     
     
     xyres@gsnXYAboveFillColors = "red"
     xyres@gsnXYBelowFillColors = "blue"
     xyres@gsnLeftString = names(ee)
     arr = new((/2,dimsizes(nino34)/),typeof(nino34))
   
     tttt = dtrend_msg(ispan(0,dimsizes(nino34)-1,1),nino34,False,True)
     arr(0,:) = (/ nino34 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino34)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(nino34),2,True)+nino34@units+" "+nyr(ee)+"yr~S~-1~N~"
     xyn34(ee) = gsn_csm_xy(wks_n34,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nino34)),arr,xyres) 
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(nino3)-1,1),nino3,False,True)
     arr(0,:) = (/ nino3 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino3)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(nino3),2,True)+nino3@units+" "+nyr(ee)+"yr~S~-1~N~"
     xyn3(ee) = gsn_csm_xy(wks_n3,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nino3)),arr,xyres) 
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(nino4)-1,1),nino4,False,True)
     arr(0,:) = (/ nino4 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino4)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(nino4),2,True)+nino4@units+" "+nyr(ee)+"yr~S~-1~N~"
     xyn4(ee) = gsn_csm_xy(wks_n4,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nino4)),arr,xyres)
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(nino12)-1,1),nino12,False,True)
     arr(0,:) = (/ nino12 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino12)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(nino12),2,True)+nino12@units+" "+nyr(ee)+"yr~S~-1~N~"
     xyn12(ee) = gsn_csm_xy(wks_n12,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nino12)),arr,xyres)
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(tna)-1,1),tna,False,True)
     arr(0,:) = (/ tna /)
     arr(1,:) = (/  (ispan(0,dimsizes(tna)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(tna),2,True)+tna@units+" "+nyr(ee)+"yr~S~-1~N~"
     xytna(ee) = gsn_csm_xy(wks_tna,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tna)),arr,xyres)
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(tsa)-1,1),tsa,False,True)
     arr(0,:) = (/ tsa /)
     arr(1,:) = (/  (ispan(0,dimsizes(tsa)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(tsa),2,True)+tsa@units+" "+nyr(ee)+"yr~S~-1~N~"
     xytsa(ee) = gsn_csm_xy(wks_tsa,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tsa)),arr,xyres)
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(tio)-1,1),tio,False,True)
     arr(0,:) = (/ tio /)
     arr(1,:) = (/  (ispan(0,dimsizes(tio)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(tio),2,True)+tio@units+" "+nyr(ee)+"yr~S~-1~N~"
     xytio(ee) = gsn_csm_xy(wks_tio,fspan(syear(ee),eyear(ee)+.91667,dimsizes(tio)),arr,xyres)
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(iod)-1,1),iod,False,True)
     arr(0,:) = (/ iod /)
     arr(1,:) = (/  (ispan(0,dimsizes(iod)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(iod),2,True)+iod@units+" "+nyr(ee)+"yr~S~-1~N~"
     xyiod(ee) = gsn_csm_xy(wks_tio,fspan(syear(ee),eyear(ee)+.91667,dimsizes(iod)),arr,xyres)
     delete([/arr,tttt/])
     
     xyres2@gsnLeftString = names(ee)
     if (nyr(ee).ge.35) then
        xyres2@gsnRightString = sprintf("%4.2f", min(sd_run))+" / "+sprintf("%4.2f", avg(sd_run))+" / "+sprintf("%4.2f", max(sd_run))
        xyn34_rst(ee) = gsn_csm_xy(wks_n34_rst,fspan(syear(ee),eyear(ee)+.91667,dimsizes(sd_run)),sd_run,xyres2)
     end if
     
     xyres3@gsnRightStringFontHeightF = xyres3@gsnCenterStringFontHeightF
     xyres3@gsnLeftString = syear(ee)+"-"+eyear(ee)
     xyres3@gsnCenterString = names(ee)
     xyres3@gsnRightString = "(~S~o~N~C)"
     xyn34_mst(ee) = gsn_csm_xy(wks_n34_mst,ispan(1,12,1),nino34_mon_sd,xyres3)
     delete([/xyres,xyres2,xyres3/])
;- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     res = True
     res@vpHeightF = 0.45
     res@vpWidthF = 0.35
     res@gsnFrame = False
     res@gsnDraw = False 
  

     res@tmYLMode = "Explicit"
;     res@tmYLValues = ispan(0,72,6)
;     res@tmYLLabels = (/"Jan~S~-2~N~","Jul~S~-2~N~","Jan~S~-1~N~","Jul~S~-1~N~", \
;                     "Jan~S~0~N~","Jul~S~0~N~","Jan~S~+1~N~","Jul~S~+1~N~", \
;                     "Jan~S~+2~N~","Jul~S~+2~N~","Jan~S~+3~N~","Jul~S~+3~N~","Jan~S~+4~N~"/)
     res@trYMinF = 24
     res@trYMaxF = 52  
     res@tmYLValues = ispan(24,52,4)
     res@tmYLLabels = (/"Jan~S~0~N~","May~S~0~N~","Sep~S~0~N~","Jan~S~+1~N~", \
                     "May~S~+1~N~","Sep~S~+1~N~","Jan~S~+2~N~","May~S~+2~N~"/)    
     res@tmYLMinorValues = ispan(24,52,2)               
     res@tmYLLabelJust = "CenterCenter"               
     res@tmYLLabelDeltaF = 1.3    ;0.05
     res@cnFillOn = True
     res@gsnSpreadColors = True  
     res@gsnSpreadColorEnd = 19
     
     res@lbLabelBarOn = False
  
     res@tiMainOn = False
     res@cnInfoLabelOn = False
     res@cnLinesOn = True
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-.75,-.5,-.25,0,.25,.5,.75,1,1.5,2,2.5,3/)   ;fspan(-2.,2.,17)
     carr = new(dimsizes(res@cnLevels),"string")
     carr = "transparent"
     carr(8) = "gray50"
     res@cnMonoLineColor = False
     res@cnLineColors = carr
     res@cnLineLabelsOn = False
     res@tmYLLabelFontHeightF = 0.014
     res@tmXBLabelFontHeightF = 0.014
     res@gsnMajorLonSpacing = 30.
     res@gsnMinorLonSpacing = 10.
     res@tiYAxisOn = False
     
     if (wks_type.eq."png") then
        res@cnLineThicknessF = 2.  
     else
        res@cnLineThicknessF = 1.  
     end if
  
     res@gsnCenterStringFontHeightF = 0.017
     res@gsnLeftStringFontHeightF = 0.017
     res@gsnRightStringFontHeightF = 0.017

     res@gsnLeftString = ""
     res@gsnCenterString= ""
     res@gsnRightString = ""

     if (isfilepresent2("obs_ts").and.ee.eq.0) then    ; for metrics table
        patcor_hov_hi = new((/nsim,dimsizes(finsst_hi&time),dimsizes(finsst_hi&lon)/),typeof(finsst_hi))
        patcor_hov_hi!1 = "time"
        patcor_hov_hi&time = finsst_hi&time
        patcor_hov_hi!2 = "lon"
        patcor_hov_hi&lon = finsst_hi&lon
           
        patcor_hov_lo = patcor_hov_hi
           
        patcor_hov_hi(ee,:,:) = (/ finsst_hi /)
        patcor_hov_lo(ee,:,:) = (/ finsst_lo /)
     end if
     if (isfilepresent2("obs_ts").and.ee.ge.1) then           
        dimT = dimsizes(finsst_hi&time)
        do hh = 0,dimT-1     ; need to loop over each timestep, using linint1 to interpolate to set longitudes. 
           patcor_hov_hi(ee,hh,:) = (/ linint1(finsst_hi&lon,finsst_hi(hh,:),False,patcor_hov_hi&lon,0) /)
           patcor_hov_lo(ee,hh,:) = (/ linint1(finsst_lo&lon,finsst_lo(hh,:),False,patcor_hov_lo&lon,0) /)
        end do
     end if

     res@gsnCenterString = names(ee)   ;"El Nin~H-13V2F35~D~FV-2H3F21~o"
     res@gsnRightString = cntr_hi
     plot_n34hi(ee) = gsn_csm_hov(wks_n34_tlon_hi,finsst_hi,res)
     
     res@gsnRightString = cntr_lo
     plot_n34lo(ee) = gsn_csm_hov(wks_n34_tlon_lo,finsst_lo,res)
     delete([/finsst_hi,finsst_lo,finsst_mid,cntr_hi,cntr_lo,cntr_mid/])
     delete(res)
;- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     pres = True
     pres@vpXF = 0.07
     pres@trYMinF = 0.
     pres@trXMinF = 0.0
;     pres@trYMaxF = 82.
     pres@trXMaxF = 0.0832
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
     pres@tmXBMode   = "Explicit"            
     pres@tmXTValues = (/".00167",".00833",".01667",".02778",".0416",".0556",".0832"/)
     pres@tmXTLabels = (/"50","10","5","3","2","1.5","1"/)           
     pres@tmXBValues = (/".0",".01",".02",".03",".042",".056",".083"/)
     pres@tmXBLabels =     pres@tmXBValues
     pres@tmXTLabelFontHeightF = 0.018
     pres@tmXBLabelFontHeightF = 0.018
     pres@tmYLLabelFontHeightF = 0.018
     pres@tiYAxisString = "Power (~S~o~N~C~S~2~N~ / cycles mo~S~-1~N~)"              ; yaxis
     pres@tiXAxisString = "Frequency (cycles mo~S~-1~N~)"
     pres@tiMainString = ""
     pres@txFontHeightF = 0.015
     pres@xyLineLabelFontHeightF = 0.022
     pres@tiXAxisFontHeightF = 0.025
     pres@tiYAxisFontHeightF = 0.025
     pres@tiMainFontHeightF = 0.03
     
     if (spectra_mvf.eq.False) then
        if (isfilepresent2("obs_ts").and.ee.ge.1.and.spectra_mvf_obs.eq.False) then
           val = new(2,typeof(sdof_obs@spcx))
           val(0) = max(sdof_obs@spcx)
           val(1) = max(splt1(0,:))
           pres@trYMaxF = max(val)+15.   ; 15 = arbitrary
           delete(val)
        else
           pres@trYMaxF = max(splt1(0,:))+15
        end if
     end if
     
     pres@tiMainOn = False
     pres@gsnCenterString = "Period (years)"
     pres@gsnCenterStringFontHeightF =    pres@tiYAxisFontHeightF
     pres@gsnRightString = ""
     pres@gsnLeftString = "" 
     if (wks_type.eq."png") then
        pres@xyLineThicknessF   = 3.5    
     else
        pres@xyLineThicknessF   = 1.5
     end if
     pres@gsnCenterString = names(ee)
     if (spectra_mvf.eq.False) then
        pspec(ee) = gsn_csm_xy(wks_n34_p,sdof@frq,splt1,pres) 
        if (isfilepresent2("obs_ts").and.ee.ge.1.and.spectra_mvf_obs.eq.False) then
           pres@xyLineColors        = (/"gray70","black","black","black"/)
           pres@gsnCenterString = ""
           pspec_obs(ee) = gsn_csm_xy(wks_n34_p,sdof_obs@frq,sdof_obs@spcx,pres)
           overlay(pspec(ee),pspec_obs(ee)) 
        end if
        delete([/sdof,splt1/])
     end if     
;- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
;     if (ee.eq.0) then
;     if (.not.isvar("title_n34")) then
        title_n34 = nino34@area
        title_n4  = nino4@area
        title_n3  = nino3@area
        title_n12 = nino12@area

        title_tna = tna@area
        title_tsa = tsa@area
        
        title_tio = tio@area
        title_iod = iod@area
;     end if
     delete([/nino34,nino3,nino4,nino12,tsa,tna,tio,iod,val1,val2,pres/])
     if (nyr(ee).ge.35) then
        delete(sd_run)
     end if
     
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
     
     scres@cnLevelSelectionMode = "ExplicitLevels"
     scres@cnLevels = (/-4,-3,-2,-1.5,-1,-.5,-.25,0,.25,.5,1,1.5,2,3.,4/)

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
     
     scres2 = True
     scres2@gsnDraw      = False
     scres2@gsnFrame     = False
     scres2@cnLevelSelectionMode = "ExplicitLevels"
     scres2@cnLevels = (/-4,-3,-2,-1.5,-1,-.5,-.25,0,.25,.5,1,1.5,2,3.,4/)

     scres2@cnLineLabelsOn = False
     scres2@cnFillOn        = True
     scres2@cnLinesOn       = False
     scres2@lbLabelBarOn    = False
     scres2@cnInfoLabelOn = False
     scres2@gsnRightString = ""
     scres2@gsnLeftString = "" 
     scres2@gsnCenterString = ""   
     scres2@gsnAddCyclic = True
     
     scres3 = True          ; PSL resources
     scres3@cnLineColor = "black"
     scres3@cnLineLabelsOn = False  
     scres3@cnLevelSelectionMode = "ExplicitLevels"
;     scres3@cnInfoLabelOn = False
     scres3@tiMainOn = False      
     new_index = NhlNewDashPattern(wks_n34sc,"$_$_$_$_$_$_$_$_$_")
     scres3@gsnContourNegLineDashPattern  = new_index
     scres3@cnLineDashSegLenF = 0.08
     scres3@gsnDraw = False
     scres3@gsnFrame = False
     scres3@gsnLeftString = ""
     scres3@gsnRightString = ""
     scres3@gsnCenterString = ""
     scres3@cnLevels = ispan(-16,16,2)
     
     if (ee.eq.0.and.n34sc_obs.eq.0) then
        scres@mpOutlineDrawOrder = "PostDraw"
        scres2@cnFillMode = "AreaFill"
        scres@cnFillMode = "RasterFill"
;        scres@cnMissingValFillColor = "white"   ; partial fix for rare weird contouring effect where some areas are not drawn?
        scres@gsnCenterString = scres@gsnCenterString+" / "+names_trefht(ee)+" / "+names_psl(ee)
     else
        scres@mpOutlineDrawOrder = "Draw"
        scres2@cnFillMode = "AreaFill"
        scres@cnFillMode = "AreaFill"    
     end if
     
     if (wks_type.eq."png") then
        scres3@cnLineThicknessF = 3.
        scres@mpGeophysicalLineThicknessF = 2.  
     else
        scres3@cnLineThicknessF = 1.25
        scres@mpGeophysicalLineThicknessF = 1.  
     end if
     
     if (n34sc_plot_flag.eq.0) then
        if (n34sc_obs.eq.0.and.ee.eq.0) then    ; for metrics table
           patcor_tas = new((/nsim,dimsizes(n34sc_tas&lat),dimsizes(n34sc_tas&lon)/),typeof(n34sc_tas))
           patcor_tas!1 = "lat"
           patcor_tas&lat = n34sc_tas&lat
           patcor_tas!2 = "lon"
           patcor_tas&lon = n34sc_tas&lon
           patcor_psl = new((/nsim,dimsizes(n34sc_psl&lat),dimsizes(n34sc_psl&lon)/),typeof(n34sc_psl))
           patcor_psl!1 = "lat"
           patcor_psl&lat = n34sc_psl&lat
           patcor_psl!2 = "lon"
           patcor_psl&lon = n34sc_psl&lon
           patcor_tas(ee,:,:) = (/ n34sc_tas(12,:,:) /)
           patcor_psl(ee,:,:) = (/ n34sc_psl(12,:,:) /)
        end if
        if (n34sc_obs.eq.0.and.ee.ge.1) then
           patcor_tas(ee,:,:) = (/ linint2(n34sc_tas&lon,n34sc_tas&lat,n34sc_tas(12,:,:),True,patcor_tas&lon,patcor_tas&lat,0) /)
           patcor_psl(ee,:,:) = (/ linint2(n34sc_psl&lon,n34sc_psl&lat,n34sc_psl(12,:,:),True,patcor_psl&lon,patcor_psl&lat,0) /)
        end if
     
        d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")     ; mask ocean for TAS array
        basemap = d->LSMASK
        lsm = landsea_mask(basemap,n34sc_tas&lat,n34sc_tas&lon)
        n34sc_tas = mask(n34sc_tas,conform(n34sc_tas,lsm,(/1,2/)).eq.0,False)
        delete([/lsm,basemap/])
        delete(d)
     
        map_n34sc_jja0(ee) = gsn_csm_contour_map(wks_n34sc,n34sc_sst(6,:,:),scres)    ; 6 = JJA 0    
        o1 = gsn_csm_contour(wks_n34sc,n34sc_tas(6,:,:),scres2)  
        o2 = gsn_csm_contour(wks_n34sc,n34sc_psl(6,:,:),scres3) 
        overlay(map_n34sc_jja0(ee),o1)
        overlay(map_n34sc_jja0(ee),o2)
        delete([/o1,o2/])
        
        map_n34sc_son0(ee) = gsn_csm_contour_map(wks_n34sc,n34sc_sst(9,:,:),scres)    ; 9 = SON 0      
        o3 = gsn_csm_contour(wks_n34sc,n34sc_tas(9,:,:),scres2)  
        o4 = gsn_csm_contour(wks_n34sc,n34sc_psl(9,:,:),scres3) 
        overlay(map_n34sc_son0(ee),o3)
        overlay(map_n34sc_son0(ee),o4)
        delete([/o3,o4/])
        
        
        map_n34sc_djf1(ee) = gsn_csm_contour_map(wks_n34sc,n34sc_sst(12,:,:),scres)    ; 12 = DJF+1      
        o5 = gsn_csm_contour(wks_n34sc,n34sc_tas(12,:,:),scres2)  
        o6 = gsn_csm_contour(wks_n34sc,n34sc_psl(12,:,:),scres3) 
        overlay(map_n34sc_djf1(ee),o5)
        overlay(map_n34sc_djf1(ee),o6)
        delete([/o5,o6/])
        
        map_n34sc_mam1(ee) = gsn_csm_contour_map(wks_n34sc,n34sc_sst(15,:,:),scres)    ; 15 = MAM+1      
        o7 = gsn_csm_contour(wks_n34sc,n34sc_tas(15,:,:),scres2)  
        o8 = gsn_csm_contour(wks_n34sc,n34sc_psl(15,:,:),scres3) 
        overlay(map_n34sc_mam1(ee),o7)
        overlay(map_n34sc_mam1(ee),o8)
        delete([/o7,o8/])       
        delete([/n34sc_sst,n34sc_tas,n34sc_psl/])
     end if
  end do
  
  if (n34sc_plot_flag.eq.0.and.n34sc_obs.eq.0.and.ee.ge.1) then    ; for pattern correlation table  
     clat_sst = cos(0.01745329*patcor_tas&lat)
     clat_psl = cos(0.01745329*patcor_psl&lat)
     finpr_sst = "ENSO TAS (DJF+1)  "    ; Must be 18 characters long
     finpr_psl = "ENSO PSL (DJF+1)  " 
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
        if (all(ismissing(patcor_tas(hh,:,:)))) then
           finpr_sst = finpr_sst+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr_sst = finpr_sst+sprintf(format2,(pattern_cor(patcor_tas(0,:,:),patcor_tas(hh,:,:),clat_sst,0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(NewCosWeight(patcor_tas(0,:,:))),ndtooned(NewCosWeight(patcor_tas(hh,:,:))))))
        end if
        if (all(ismissing(patcor_psl(hh,:,:)))) then
           finpr_psl = finpr_psl+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr_psl = finpr_psl+sprintf(format2,(pattern_cor(patcor_psl(0,:,:),patcor_psl(hh,:,:),clat_psl,0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(NewCosWeight(patcor_psl(0,:,:))),ndtooned(NewCosWeight(patcor_psl(hh,:,:))))))
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
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr_sst/],"%s")
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr_psl/],"%s")
     end if  
     delete([/finpr_sst,finpr_psl,line3,line4,format2,format3,nchar,ntc,clat_sst,clat_psl,patcor_tas,patcor_psl,dimY,ntb/])
  end if
  if (isfilepresent2("obs_ts").and.ee.ge.1) then           ; for pattern correlation table  
     finpr_hi = "El Nino Hovmoller "    ; Must be 18 characters long
     finpr_lo = "La Nina Hovmoller " 
     line3   = "                  "    ; Must be 18 characters long                                 patcor_hov_hi
     line4   = line3
     do hh = 1,nsim-1
        nchar = dimsizes(tochar(names(hh)))
        nchar = where(nchar.le.10,10,nchar)
        ntc = ""
        do ii = 0,nchar-1
           ntc = ntc+"-"
        end do
        format2 = "%"+(nchar-5+1)+".2f"
        format3 = "%4.2f"
        line3 = line3+" "+names(hh)   ;
        line4 = line4+" "+ntc 
        if (all(ismissing(patcor_hov_hi(hh,24:52,:)))) then                ; 24:52 refers to Jan+0->May+2, which is the range shown in the hovmoller plots.
           finpr_hi = finpr_hi+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr_hi = finpr_hi+sprintf(format2,(pattern_cor(patcor_hov_hi(0,24:52,:),patcor_hov_hi(hh,24:52,:),1.0,0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(patcor_hov_hi(0,24:52,:)),ndtooned(patcor_hov_hi(hh,24:52,:)))))
        end if
        if (all(ismissing(patcor_hov_lo(hh,24:52,:)))) then
           finpr_lo = finpr_lo+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr_lo = finpr_lo+sprintf(format2,(pattern_cor(patcor_hov_lo(0,24:52,:),patcor_hov_lo(hh,24:52,:),1.0,0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(patcor_hov_lo(0,24:52,:)),ndtooned(patcor_hov_lo(hh,24:52,:)))))
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
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr_hi/],"%s")
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr_lo/],"%s")
     end if
     delete([/finpr_hi,finpr_lo,line3,line4,format2,format3,nchar,ntc,patcor_hov_hi,patcor_hov_lo/])  
  end if
  
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

  panres@txString = "Nino3.4 (Monthly, "+title_n34+")"
  gsn_panel(wks_n34,xyn34,lp,panres)    

  panres@txString = "Nino4 (Monthly, "+title_n4+")"
  gsn_panel(wks_n4,xyn4,lp,panres)  

  panres@txString = "Nino3 (Monthly, "+title_n3+")"
  gsn_panel(wks_n3,xyn3,lp,panres)       

  panres@txString = "Nino1+2 (Monthly, "+title_n12+")"
  gsn_panel(wks_n12,xyn12,lp,panres)       
   
  panres@txString = "Tropical North Atlantic (Monthly, "+title_tna+")"
  gsn_panel(wks_tna,xytna,lp,panres)
  
  panres@txString = "Tropical South Atlantic (Monthly, "+title_tsa+")"
  gsn_panel(wks_tsa,xytsa,lp,panres)
     
  panres@txString = "Tropical Indian Ocean (Monthly, "+title_tio+")"
  gsn_panel(wks_tio,xytio,lp,panres)   
  
  panres@txString = "Indian Ocean Dipole (Monthly, "+title_iod+")"
  gsn_panel(wks_tio,xyiod,lp,panres)  
  delete(wks_tio) 

  if (all(ismissing(xyn34_rst))) then
     print("No valid running standard deviation plots, skipping")
  else
     panres@txString = "Running nino3.4 standard deviation (Monthly)"
     gsn_panel(wks_n34_rst,xyn34_rst,lp,panres) 
     delete(xyn34_rst)  
  end if
  panres@gsnPanelYWhiteSpacePercent = 0.5
  panres@txString = "nino3.4 standard deviation (Monthly)"
  gsn_panel(wks_n34_mst,xyn34_mst,(/nrow,ncol/),panres)   
  
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

  
  panres2@txString = "El Nin~H-13V2F35~D~FV-2H3F21~o Composite (3~S~o~N~S:3~S~o~N~N)"
  gsn_panel(wks_n34_tlon_hi,plot_n34hi,(/nrow,ncol/),panres2)   
  panres2@txString = "La Nin~H-13V2F35~D~FV-2H3F21~a Composite (3~S~o~N~S:3~S~o~N~N)"
  gsn_panel(wks_n34_tlon_lo,plot_n34lo,(/nrow,ncol/),panres2)   
  
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@txString = "nino3.4 (Monthly, detrended)"
  gsn_panel(wks_n34_p,pspec,(/nrow,ncol/),panres)  
  delete(wks_n34_p)
  
  if (n34sc_frame.eq.0) then
     panres2@txString = "nino3.4 TS,TAS,PSL Spatial Composite (JJA~S~0~N~)"
     gsn_panel(wks_n34sc,map_n34sc_jja0,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 TS,TAS,PSL Spatial Composite (SON~S~0~N~)"
     gsn_panel(wks_n34sc,map_n34sc_son0,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 TS,TAS,PSL Spatial Composite (DJF~S~+1~N~)"
     gsn_panel(wks_n34sc,map_n34sc_djf1,(/nrow,ncol/),panres2)   
  
     panres2@txString = "nino3.4 TS,TAS,PSL Spatial Composite (MAM~S~+1~N~)"
     gsn_panel(wks_n34sc,map_n34sc_mam1,(/nrow,ncol/),panres2)   
     delete(wks_n34sc)
     
     delete([/map_n34sc_djf1,map_n34sc_jja0,map_n34sc_son0,map_n34sc_mam1/])
  end if     
  delete([/xyn34,xyn4,xyn3,xyn12,xytna,xytsa,xytio,xyiod,plot_n34hi,plot_n34lo,pspec,pi,rad,wgt,lp/])     

  OUTDIR = getenv("OUTDIR")
  
  if (n34sc_frame.eq.0) then
     if (wks_type.eq."png") then  
        system("mv "+OUTDIR+"nino34.spatialcomp.000001.png "+OUTDIR+"nino34.spatialcomp.jja0.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.000002.png "+OUTDIR+"nino34.spatialcomp.son0.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.000003.png "+OUTDIR+"nino34.spatialcomp.djf1.png")
        system("mv "+OUTDIR+"nino34.spatialcomp.000004.png "+OUTDIR+"nino34.spatialcomp.mam1.png")
     else
        system("psplit "+OUTDIR+"nino34.spatialcomp.ps "+OUTDIR+"pict")
        system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nino34.spatialcomp.jja0.ps")
        system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nino34.spatialcomp.son0.ps")
        system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nino34.spatialcomp.djf1.ps")
        system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nino34.spatialcomp.mam1.ps")
        system("rm "+OUTDIR+"nino34.spatialcomp.ps")
     end if
  end if
  
  if (wks_type.eq."png") then  
     system("mv "+OUTDIR+"tio.timeseries.000001.png "+OUTDIR+"tio.timeseries.png")
     system("mv "+OUTDIR+"tio.timeseries.000002.png "+OUTDIR+"iod.timeseries.png")
  else
     system("psplit "+OUTDIR+"tio.timeseries.ps "+OUTDIR+"pict")
     system("rm "+OUTDIR+"tio.timeseries.ps")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"tio.timeseries.ps")
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"iod.timeseries.ps")
  end if
  print("Finished: sst.indices.ncl")
end
