; Calculates the coupled and atmospheric modes of variability (patterns and timeseries), 
; as well as regressions of those PC timeseries onto the land-masked TS field. 
;
; Variables used: PSL and TS
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: psl.modes_indices.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COMPUTE_MODES_MON = getenv("COMPUTE_MODES_MON")
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

  pi=4.*atan(1.0)
  rad=(pi/180.)
  
;---------SST Regressions coding-------------------------------------------------
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
    
  if (isfilepresent2("obs_psl")) then
     if (nsim.eq.nsim_ts) then
        ; print("Total number of simulations/obs are equal between psl/ts namelists")
        sstreg_obs = 0
     else
        print("Total number of simulations/obs between psl/ts do not match, skipping observational SST regressions in psl.modes.indices.ncl.")
        sstreg_obs = 1
        ; delete([/nsim_ts,na_ts,names_ts,paths_ts,syear_ts,eyear_ts,nyr_ts/])
     end if  
  else
     sstreg_obs = 1
  end if
  
  if (sstreg_obs.eq.0) then     ; further check for obs
     if (syear(0).eq.syear_ts(0)) then
        if (eyear(0).eq.eyear_ts(0)) then
;           print("Obs years match")
           sstreg_obs = 0
        else
           print("Obs end years do not match, skipping observed SST regressions in psl.modes.indices.ncl.")
           sstreg_obs = 1
        end if
     else
        print("Obs start years do not match, skipping observed SST regressions in psl.modes.indices.ncl.")
        sstreg_obs = 1
     end if
  end if
;-------------------------------------------------------------------------------------------------
    
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_nam = gsn_open_wks(wks_type,getenv("OUTDIR")+"nam")
  wks_nam_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"nam.timeseries")
  
  wks_sam = gsn_open_wks(wks_type,getenv("OUTDIR")+"sam")
  wks_sam_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"sam.timeseries")
  
  wks_nao = gsn_open_wks(wks_type,getenv("OUTDIR")+"nao")
  wks_nao_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"nao.timeseries")
  
  wks_pna = gsn_open_wks(wks_type,getenv("OUTDIR")+"pna")
  wks_pna_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"pna.timeseries")
  
  wks_npo = gsn_open_wks(wks_type,getenv("OUTDIR")+"npo")
  wks_npo_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"npo.timeseries")
  
  wks_psa1 = gsn_open_wks(wks_type,getenv("OUTDIR")+"psa1")
  wks_psa1_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"psa1.timeseries")

  wks_psa2 = gsn_open_wks(wks_type,getenv("OUTDIR")+"psa2")
  wks_psa2_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"psa2.timeseries")
  
;   Note: Currently near maximum number of open workstations: using 14 out of 15.  
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_nam,"ncl_default")   
     gsn_define_colormap(wks_nam_ts,"ncl_default")
     gsn_define_colormap(wks_sam,"ncl_default")   
     gsn_define_colormap(wks_sam_ts,"ncl_default")
     gsn_define_colormap(wks_nao,"ncl_default")   
     gsn_define_colormap(wks_nao_ts,"ncl_default")
     gsn_define_colormap(wks_pna,"ncl_default")   
     gsn_define_colormap(wks_pna_ts,"ncl_default")
     gsn_define_colormap(wks_npo,"ncl_default")   
     gsn_define_colormap(wks_npo_ts,"ncl_default")
     gsn_define_colormap(wks_psa1,"ncl_default")   
     gsn_define_colormap(wks_psa1_ts,"ncl_default")
     gsn_define_colormap(wks_psa2,"ncl_default")   
     gsn_define_colormap(wks_psa2_ts,"ncl_default")
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_nam,"BlueDarkRed18")   
     gsn_define_colormap(wks_nam_ts,"ncl_default")
     gsn_define_colormap(wks_sam,"BlueDarkRed18")   
     gsn_define_colormap(wks_sam_ts,"ncl_default")
     gsn_define_colormap(wks_nao,"BlueDarkRed18")   
     gsn_define_colormap(wks_nao_ts,"ncl_default")
     gsn_define_colormap(wks_pna,"BlueDarkRed18")   
     gsn_define_colormap(wks_pna_ts,"ncl_default")
     gsn_define_colormap(wks_npo,"BlueDarkRed18")   
     gsn_define_colormap(wks_npo_ts,"ncl_default")
     gsn_define_colormap(wks_psa1,"BlueDarkRed18")   
     gsn_define_colormap(wks_psa1_ts,"ncl_default")
     gsn_define_colormap(wks_psa2,"BlueDarkRed18")   
     gsn_define_colormap(wks_psa2_ts,"ncl_default")
  end if
  
  
  
  map_nam_djf = new(nsim,"graphic")  
  map_nam_mam = new(nsim,"graphic")  
  map_nam_jja = new(nsim,"graphic")  
  map_nam_son = new(nsim,"graphic")  
  map_nam_ann = new(nsim,"graphic")  
  map_nam_mon = new(nsim,"graphic")  
  xy_nam_djf = new(nsim,"graphic")
  xy_nam_mam = new(nsim,"graphic")
  xy_nam_jja = new(nsim,"graphic")
  xy_nam_son = new(nsim,"graphic")
  xy_nam_ann = new(nsim,"graphic")
  xy_nam_mon = new(nsim,"graphic")
  reg_nam_djf = new(nsim,"graphic")  
  reg_nam_mam = new(nsim,"graphic")  
  reg_nam_jja = new(nsim,"graphic")  
  reg_nam_son = new(nsim,"graphic")  
  reg_nam_ann = new(nsim,"graphic")  
  reg_nam_mon = new(nsim,"graphic")  
  
  map_sam_djf = new(nsim,"graphic")  
  map_sam_mam = new(nsim,"graphic")  
  map_sam_jja = new(nsim,"graphic")  
  map_sam_son = new(nsim,"graphic")  
  map_sam_ann = new(nsim,"graphic")  
  map_sam_mon = new(nsim,"graphic")  
  xy_sam_djf = new(nsim,"graphic")
  xy_sam_mam = new(nsim,"graphic")
  xy_sam_jja = new(nsim,"graphic")
  xy_sam_son = new(nsim,"graphic")
  xy_sam_ann = new(nsim,"graphic")
  xy_sam_mon = new(nsim,"graphic")
  reg_sam_djf = new(nsim,"graphic")  
  reg_sam_mam = new(nsim,"graphic")  
  reg_sam_jja = new(nsim,"graphic")  
  reg_sam_son = new(nsim,"graphic")  
  reg_sam_ann = new(nsim,"graphic")  
  reg_sam_mon = new(nsim,"graphic")  

  map_nao_djf = new(nsim,"graphic")  
  map_nao_mam = new(nsim,"graphic")  
  map_nao_jja = new(nsim,"graphic")  
  map_nao_son = new(nsim,"graphic")  
  map_nao_ann = new(nsim,"graphic")  
  map_nao_mon = new(nsim,"graphic")  
  xy_nao_djf = new(nsim,"graphic")
  xy_nao_mam = new(nsim,"graphic")
  xy_nao_jja = new(nsim,"graphic")
  xy_nao_son = new(nsim,"graphic")
  xy_nao_ann = new(nsim,"graphic")
  xy_nao_mon = new(nsim,"graphic")
  reg_nao_djf = new(nsim,"graphic")  
  reg_nao_mam = new(nsim,"graphic")  
  reg_nao_jja = new(nsim,"graphic")  
  reg_nao_son = new(nsim,"graphic")  
  reg_nao_ann = new(nsim,"graphic")  
  reg_nao_mon = new(nsim,"graphic")  

  map_pna_djf = new(nsim,"graphic")  
  map_pna_mam = new(nsim,"graphic")  
  map_pna_jja = new(nsim,"graphic")  
  map_pna_son = new(nsim,"graphic")  
  map_pna_ann = new(nsim,"graphic")  
  map_pna_mon = new(nsim,"graphic")  
  xy_pna_djf = new(nsim,"graphic")
  xy_pna_mam = new(nsim,"graphic")
  xy_pna_jja = new(nsim,"graphic")
  xy_pna_son = new(nsim,"graphic")
  xy_pna_ann = new(nsim,"graphic")
  xy_pna_mon = new(nsim,"graphic")
  reg_pna_djf = new(nsim,"graphic")  
  reg_pna_mam = new(nsim,"graphic")  
  reg_pna_jja = new(nsim,"graphic")  
  reg_pna_son = new(nsim,"graphic")  
  reg_pna_ann = new(nsim,"graphic")  
  reg_pna_mon = new(nsim,"graphic")  
  
  map_npo_djf = new(nsim,"graphic")  
  map_npo_mam = new(nsim,"graphic")  
  map_npo_jja = new(nsim,"graphic")  
  map_npo_son = new(nsim,"graphic")  
  map_npo_ann = new(nsim,"graphic")  
  map_npo_mon = new(nsim,"graphic")  
  xy_npo_djf = new(nsim,"graphic")
  xy_npo_mam = new(nsim,"graphic")
  xy_npo_jja = new(nsim,"graphic")
  xy_npo_son = new(nsim,"graphic")
  xy_npo_ann = new(nsim,"graphic")
  xy_npo_mon = new(nsim,"graphic")
  reg_npo_djf = new(nsim,"graphic")  
  reg_npo_mam = new(nsim,"graphic")  
  reg_npo_jja = new(nsim,"graphic")  
  reg_npo_son = new(nsim,"graphic")  
  reg_npo_ann = new(nsim,"graphic")  
  reg_npo_mon = new(nsim,"graphic")  
  
  map_psa1_djf = new(nsim,"graphic")  
  map_psa1_mam = new(nsim,"graphic")  
  map_psa1_jja = new(nsim,"graphic")  
  map_psa1_son = new(nsim,"graphic")  
  map_psa1_ann = new(nsim,"graphic")  
  map_psa1_mon = new(nsim,"graphic")  
  xy_psa1_djf = new(nsim,"graphic")
  xy_psa1_mam = new(nsim,"graphic")
  xy_psa1_jja = new(nsim,"graphic")
  xy_psa1_son = new(nsim,"graphic")
  xy_psa1_ann = new(nsim,"graphic")
  xy_psa1_mon = new(nsim,"graphic")
  reg_psa1_djf = new(nsim,"graphic")  
  reg_psa1_mam = new(nsim,"graphic")  
  reg_psa1_jja = new(nsim,"graphic")  
  reg_psa1_son = new(nsim,"graphic")  
  reg_psa1_ann = new(nsim,"graphic")  
  reg_psa1_mon = new(nsim,"graphic")  
  
  map_psa2_djf = new(nsim,"graphic")  
  map_psa2_mam = new(nsim,"graphic")  
  map_psa2_jja = new(nsim,"graphic")  
  map_psa2_son = new(nsim,"graphic")  
  map_psa2_ann = new(nsim,"graphic")  
  map_psa2_mon = new(nsim,"graphic")  
  xy_psa2_djf = new(nsim,"graphic")
  xy_psa2_mam = new(nsim,"graphic")
  xy_psa2_jja = new(nsim,"graphic")
  xy_psa2_son = new(nsim,"graphic")
  xy_psa2_ann = new(nsim,"graphic")
  xy_psa2_mon = new(nsim,"graphic")
  reg_psa2_djf = new(nsim,"graphic")  
  reg_psa2_mam = new(nsim,"graphic")  
  reg_psa2_jja = new(nsim,"graphic")  
  reg_psa2_son = new(nsim,"graphic")  
  reg_psa2_ann = new(nsim,"graphic")  
  reg_psa2_mon = new(nsim,"graphic")  
  
  xy_npi = new(nsim,"graphic")
  sstreg_frame = 1    ; sstreg_frame = flag to create regressions .ps/.png files. Created/used instead of sstreg_plot_flag
                      ; so that if sst regressions are not created for the last simulation listed that .ps/png files are created
  do ee = 0,nsim-1
;     print(paths(ee)+" "+syear(ee)+" "+eyear(ee)) 
     arr = data_read_in(paths(ee),"PSL",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(arr,"is_all_missing")) then
        delete(arr)
        continue
     end if 
     
     if (OPT_CLIMO.eq."Full") then
        arr = rmMonAnnCycTLL(arr)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = arr
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(arr&time,-1)
        climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        delete(temp_arr)
        arr   = calcMonAnomTLL(arr,climo) 
        delete(climo)
     end if

     arrT =  runave_n_Wrap(arr,3,0,0)               ; form DJF averages
     arrT(0,:,:) = (/ dim_avg_n(arr(:1,:,:),0) /)   
     arr_djf = arrT(0::12,:,:)
     arr_mam = arrT(3::12,:,:)
     arr_jja = arrT(6::12,:,:)                       ; form JJA averages
     arr_son = arrT(9::12,:,:)
     delete(arrT)
     
     arrU = runave_n_Wrap(arr,5,0,0)
     arrU(0,:,:) = (/ dim_avg_n(arr(:2,:,:),0) /)
     arr_ndjfm = arrU(0::12,:,:)
     delete(arrU)
     
     arrV = runave_n_Wrap(arr,12,0,0)
     arr_ann = arrV(5::12,:,:)
     delete(arrV)
;          
;     arr_djf = (/ dtrend_msg_n(ispan(0,dimsizes(arr_djf&time)-1,1),arr_djf,True,False,0) /)
;     arr_mam = (/ dtrend_msg_n(ispan(0,dimsizes(arr_mam&time)-1,1),arr_mam,True,False,0) /)
;     arr_jja = (/ dtrend_msg_n(ispan(0,dimsizes(arr_jja&time)-1,1),arr_jja,True,False,0) /)
;     arr_son = (/ dtrend_msg_n(ispan(0,dimsizes(arr_son&time)-1,1),arr_son,True,False,0) /)
;     
;     arr_ann = (/ dtrend_msg_n(ispan(0,dimsizes(arr_ann&time)-1,1),arr_ann,True,False,0) /)
;     
;     arr_ndjfm = (/ dtrend_msg_n(ispan(0,dimsizes(arr_ndjfm&time)-1,1),arr_ndjfm,True,False,0) /)
;     
;     arr = (/ dtrend_msg_n(ispan(0,dimsizes(arr&time)-1,1),arr,True,False,0) /)
;---------SST Regressions coding-------------------------------------------------
     sstreg_plot_flag = 0   
     if (ee.eq.0.and.sstreg_obs.eq.1) then    ; first time thru, and obs = False
        sstreg_plot_flag = 1                  ; do not plot
        if (.not.isfilepresent2("obs_psl")) then  ; however, if there are no PSL obs whatsoever, the 1st iteration of the
           sstreg_plot_flag = 0                  ; ee do loop will be the first model, so turn it on.
        end if
     end if

     if (sstreg_plot_flag.eq.0) then 
        ts_offset     = 0                           ; offsets necessary if OBS is present for one variable, but not another.
        if (ee.eq.0.and.sstreg_obs.eq.0) then
        else
           if (isfilepresent2("obs_psl")) then
              if (.not.isfilepresent2("obs_ts")) then
                 ts_offset = -1
              end if
           else
              if (isfilepresent2("obs_ts")) then
                 ts_offset = 1
              end if
           end if
        end if

        ; print("Data to be read in: "+paths_ts(ee+ts_offset)+" from "+syear_ts(ee+ts_offset)+":"+eyear_ts(ee+ts_offset))
        sst = data_read_in(paths_ts(ee+ts_offset),"TS",syear_ts(ee+ts_offset),eyear_ts(ee+ts_offset))
        if (isatt(sst,"is_all_missing")) then
           sstreg_plot_flag = 1
           delete(sst)
        end if
        
        if (sstreg_plot_flag.eq.0) then     ; only continue if both PSL/TS fields are present
           sst = where(sst.le.-1.8,-1.8,sst)
           d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
           basemap = d->LSMASK
           lsm = landsea_mask(basemap,sst&lat,sst&lon)
           sst = mask(sst,conform(sst,lsm,(/1,2/)).ge.1,False)
           delete(lsm)
           
           if (OPT_CLIMO.eq."Full") then
              sst = rmMonAnnCycTLL(sst)
           else
              check_custom_climo(names_ts(ee+ts_offset),syear_ts(ee+ts_offset),eyear_ts(ee+ts_offset),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = sst
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(sst&time,1)
              climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              delete(temp_arr)
              sst   = calcMonAnomTLL(sst,climo) 
              delete(climo)
           end if
;           sst = (/ dtrend_msg_n(ispan(0,dimsizes(sst&time)-1,1),sst,False,False,0) /)  
        
           sstT =  runave_n_Wrap(sst,3,0,0)               ; form DJF averages
           sstT(0,:,:) = (/ dim_avg_n(sst(:1,:,:),0) /)   
           sst_djf = sstT(0::12,:,:)
           sst_mam = sstT(3::12,:,:)
           sst_jja = sstT(6::12,:,:)                       ; form JJA averages
           sst_son = sstT(9::12,:,:)
           delete(sstT)
     
           sstV = runave_n_Wrap(sst,12,0,0)
           sst_ann = sstV(5::12,:,:)
           delete(sstV)
        end if
     end if
;----------------NPI calculation-----------------------------------     
     coswgt=cos(rad*arr_djf&lat)
     coswgt!0 = "lat"
     coswgt&lat = arr_djf&lat
     npi_ndjfm = wgt_areaave(arr_ndjfm(:,{30:65},{160:220}),coswgt({30.:65.}),1.0,0)
     npi_ndjfm!0 = "TIME"
     npi_ndjfm&TIME = ispan(syear(ee),eyear(ee),1)
     npi_ndjfm&TIME@units = "YYYY"
     npi_ndjfm&TIME@long_name = "time"
     
     npi_ndjfm@area = "30:65N, 160:220E"
     npi_ndjfm@units = arr_ndjfm@units
     npi_ndjfm@long_name = "North Pacific Index"
     delete([/coswgt,arr_ndjfm/])
;---------NAM calculations----------------------------------------------------------     
     arr_djf_CW = SqrtCosWeight(arr_djf)
     arr_mam_CW = SqrtCosWeight(arr_mam)
     arr_jja_CW = SqrtCosWeight(arr_jja)
     arr_son_CW = SqrtCosWeight(arr_son)
     arr_ann_CW = SqrtCosWeight(arr_ann)
     if (COMPUTE_MODES_MON.eq."True") then   
        arr_mon_CW = SqrtCosWeight(arr)
     else
        delete(arr)
        if (sstreg_plot_flag.eq.0) then
           delete(sst)
        end if
     end if
  
     evecv = eofunc(arr_djf_CW({lat|20:},lon|:,time|:),1,75)
     pcts = eofunc_ts(arr_djf_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_djf = dim_standardize(pcts(0,:),0)
     nam_djf = arr_djf(0,:,:)
     nam_djf = (/ regCoef(nam_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_djf = sst_djf(0,:,:)
        nam_sst_djf = (/ regCoef(nam_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_djf({85},{5}))) then
        if (nam_djf({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_djf = nam_djf*-1.
           nam_pc_djf = nam_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_djf = nam_sst_djf*-1.
           end if
        end if
     end if
     nam_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nam_pc_djf)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_mam_CW({lat|20:},lon|:,time|:),1,75)
     pcts = eofunc_ts(arr_mam_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_mam = dim_standardize(pcts(0,:),0)
     nam_mam = arr_mam(0,:,:)
     nam_mam = (/ regCoef(nam_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_mam = sst_mam(0,:,:)
        nam_sst_mam = (/ regCoef(nam_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_mam({85},{5}))) then
        if (nam_mam({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_mam = nam_mam*-1.
           nam_pc_mam = nam_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_mam = nam_sst_mam*-1.
           end if
        end if
     end if
     nam_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))     
     copy_VarCoords(npi_ndjfm,nam_pc_mam)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_jja_CW({lat|20:},lon|:,time|:),1,75)
     pcts = eofunc_ts(arr_jja_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_jja = dim_standardize(pcts(0,:),0)
     nam_jja = arr_jja(0,:,:)
     nam_jja = (/ regCoef(nam_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_jja = sst_jja(0,:,:)
        nam_sst_jja = (/ regCoef(nam_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_jja({85},{5}))) then
        if (nam_jja({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_jja = nam_jja*-1.
           nam_pc_jja = nam_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_jja = nam_sst_jja*-1.
           end if
        end if
     end if
     nam_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nam_pc_jja)  
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_son_CW({lat|20:},lon|:,time|:),1,75)
     pcts = eofunc_ts(arr_son_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_son = dim_standardize(pcts(0,:),0)
     nam_son = arr_son(0,:,:)
     nam_son = (/ regCoef(nam_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_son = sst_son(0,:,:)
        nam_sst_son = (/ regCoef(nam_pc_son,sst_son(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_son({85},{5}))) then
        if (nam_son({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_son = nam_son*-1.
           nam_pc_son = nam_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_son = nam_sst_son*-1.
           end if
        end if
     end if
     nam_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nam_pc_son)    
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_ann_CW({lat|20:},lon|:,time|:),1,75)
     pcts = eofunc_ts(arr_ann_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_ann = dim_standardize(pcts(0,:),0)
     nam_ann = arr_ann(0,:,:)
     nam_ann = (/ regCoef(nam_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_ann = sst_ann(0,:,:)
        nam_sst_ann = (/ regCoef(nam_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_ann({85},{5}))) then
        if (nam_ann({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_ann = nam_ann*-1.
           nam_pc_ann = nam_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_ann = nam_sst_ann*-1.
           end if
        end if
     end if
     nam_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))    
     copy_VarCoords(npi_ndjfm,nam_pc_ann) 
     delete([/evecv,pcts/])    
     
     if (COMPUTE_MODES_MON.eq."True") then
        evecv = eofunc(arr_mon_CW({lat|20:},lon|:,time|:),1,75)
        pcts = eofunc_ts(arr_mon_CW({lat|20:},lon|:,time|:),evecv,False)  
        nam_pc_mon = dim_standardize(pcts(0,:),0)
        nam_mon = arr(0,:,:)
        nam_mon = (/ regCoef(nam_pc_mon,arr(lat|:,lon|:,time|:)) /)
        if (sstreg_plot_flag.eq.0) then
           nam_sst_mon = sst(0,:,:)
           nam_sst_mon = (/ regCoef(nam_pc_mon,sst(lat|:,lon|:,time|:)) /)
        end if
        if (.not.ismissing(nam_mon({85},{5}))) then
           if (nam_mon({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              nam_mon = nam_mon*-1.
              nam_pc_mon = nam_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 nam_sst_mon = nam_sst_mon*-1.
              end if
           end if
        end if
        nam_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0))) 
        nam_pc_mon!0 = "time"
        nam_pc_mon&time = arr&time   
        delete([/evecv,pcts/])
     end if
;----------SAM/PSA1/PSA2 calculations----------------------------------------------------------------------     
     evecv = eofunc(arr_djf_CW({lat|:-20},lon|:,time|:),3,75)
     pcts = eofunc_ts(arr_djf_CW({lat|:-20},lon|:,time|:),evecv,False)  
     sam_pc_djf = dim_standardize(pcts(0,:),0)
     psa1_pc_djf = dim_standardize(pcts(1,:),0)
     psa2_pc_djf = dim_standardize(pcts(2,:),0)     
     sam_djf = arr_djf(0,:,:)
     sam_djf = (/ regCoef(sam_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     psa1_djf = arr_djf(0,:,:)
     psa1_djf = (/ regCoef(psa1_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     psa2_djf = arr_djf(0,:,:)
     psa2_djf = (/ regCoef(psa2_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        sam_sst_djf = sst_djf(0,:,:)
        sam_sst_djf = (/ regCoef(sam_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
        psa1_sst_djf = sst_djf(0,:,:)
        psa1_sst_djf = (/ regCoef(psa1_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
        psa2_sst_djf = sst_djf(0,:,:)
        psa2_sst_djf = (/ regCoef(psa2_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
     end if     
     
     if (.not.ismissing(sam_djf({-85},{5}))) then
        if (sam_djf({-85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           sam_djf = sam_djf*-1.
           sam_pc_djf = sam_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              sam_sst_djf = sam_sst_djf*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa1_djf({-62},{270}))) then
        if (psa1_djf({-62},{270}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa1_djf = psa1_djf*-1.
           psa1_pc_djf = psa1_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              psa1_sst_djf = psa1_sst_djf*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa2_djf({-52},{150}))) then
        if (psa2_djf({-52},{150}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa2_djf = psa2_djf*-1.
           psa2_pc_djf = psa2_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              psa2_sst_djf = psa2_sst_djf*-1.
           end if
        end if     
     end if
     sam_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     psa1_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     psa2_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(2)))
     copy_VarCoords(npi_ndjfm,sam_pc_djf) 
     copy_VarCoords(npi_ndjfm,psa1_pc_djf)
     copy_VarCoords(npi_ndjfm,psa2_pc_djf)  
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_mam_CW({lat|:-20},lon|:,time|:),3,75)
     pcts = eofunc_ts(arr_mam_CW({lat|:-20},lon|:,time|:),evecv,False)  
     sam_pc_mam = dim_standardize(pcts(0,:),0)
     psa1_pc_mam = dim_standardize(pcts(1,:),0)
     psa2_pc_mam = dim_standardize(pcts(2,:),0)     
     sam_mam = arr_mam(0,:,:)
     sam_mam = (/ regCoef(sam_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     psa1_mam = arr_mam(0,:,:)
     psa1_mam = (/ regCoef(psa1_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     psa2_mam = arr_mam(0,:,:)
     psa2_mam = (/ regCoef(psa2_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        sam_sst_mam = sst_mam(0,:,:)
        sam_sst_mam = (/ regCoef(sam_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
        psa1_sst_mam = sst_mam(0,:,:)
        psa1_sst_mam = (/ regCoef(psa1_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
        psa2_sst_mam = sst_mam(0,:,:)
        psa2_sst_mam = (/ regCoef(psa2_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
     end if  
     
     if (.not.ismissing(sam_mam({-85},{5}))) then
        if (sam_mam({-85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           sam_mam = sam_mam*-1.
           sam_pc_mam = sam_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              sam_sst_mam = sam_sst_mam*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa1_mam({-62},{270}))) then
        if (psa1_mam({-62},{270}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa1_mam = psa1_mam*-1.
           psa1_pc_mam = psa1_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              psa1_sst_mam = psa1_sst_mam*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa2_mam({-52},{150}))) then
        if (psa2_mam({-52},{150}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa2_mam = psa2_mam*-1.
           psa2_pc_mam = psa2_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              psa2_sst_mam = psa2_sst_mam*-1.
           end if
        end if     
     end if
     sam_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     psa1_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     psa2_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(2)))
     copy_VarCoords(npi_ndjfm,sam_pc_mam) 
     copy_VarCoords(npi_ndjfm,psa1_pc_mam)
     copy_VarCoords(npi_ndjfm,psa2_pc_mam)  
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_jja_CW({lat|:-20},lon|:,time|:),3,75)
     pcts = eofunc_ts(arr_jja_CW({lat|:-20},lon|:,time|:),evecv,False)  
     sam_pc_jja = dim_standardize(pcts(0,:),0)
     psa1_pc_jja = dim_standardize(pcts(1,:),0)
     psa2_pc_jja = dim_standardize(pcts(2,:),0)     
     sam_jja = arr_jja(0,:,:)
     sam_jja = (/ regCoef(sam_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     psa1_jja = arr_jja(0,:,:)
     psa1_jja = (/ regCoef(psa1_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     psa2_jja = arr_jja(0,:,:)
     psa2_jja = (/ regCoef(psa2_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        sam_sst_jja = sst_jja(0,:,:)
        sam_sst_jja = (/ regCoef(sam_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
        psa1_sst_jja = sst_jja(0,:,:)
        psa1_sst_jja = (/ regCoef(psa1_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
        psa2_sst_jja = sst_jja(0,:,:)
        psa2_sst_jja = (/ regCoef(psa2_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
     end if  
     if (.not.ismissing(sam_jja({-85},{5}))) then
        if (sam_jja({-85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           sam_jja = sam_jja*-1.
           sam_pc_jja = sam_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              sam_sst_jja = sam_sst_jja*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa1_jja({-62},{270}))) then
        if (psa1_jja({-62},{270}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa1_jja = psa1_jja*-1.
           psa1_pc_jja = psa1_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              psa1_sst_jja = psa1_sst_jja*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa2_jja({-52},{150}))) then
        if (psa2_jja({-52},{150}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa2_jja = psa2_jja*-1.
           psa2_pc_jja = psa2_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              psa2_sst_jja = psa2_sst_jja*-1.
           end if
        end if     
     end if
     sam_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     psa1_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     psa2_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(2)))
     copy_VarCoords(npi_ndjfm,sam_pc_jja) 
     copy_VarCoords(npi_ndjfm,psa1_pc_jja)
     copy_VarCoords(npi_ndjfm,psa2_pc_jja)  
     delete([/evecv,pcts/])
 
     evecv = eofunc(arr_son_CW({lat|:-20},lon|:,time|:),3,75)
     pcts = eofunc_ts(arr_son_CW({lat|:-20},lon|:,time|:),evecv,False)  
     sam_pc_son = dim_standardize(pcts(0,:),0)
     psa1_pc_son = dim_standardize(pcts(1,:),0)
     psa2_pc_son = dim_standardize(pcts(2,:),0)     
     sam_son = arr_son(0,:,:)
     sam_son = (/ regCoef(sam_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     psa1_son = arr_son(0,:,:)
     psa1_son = (/ regCoef(psa1_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     psa2_son = arr_son(0,:,:)
     psa2_son = (/ regCoef(psa2_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        sam_sst_son = sst_son(0,:,:)
        sam_sst_son = (/ regCoef(sam_pc_son,sst_son(lat|:,lon|:,time|:)) /)
        psa1_sst_son = sst_son(0,:,:)
        psa1_sst_son = (/ regCoef(psa1_pc_son,sst_son(lat|:,lon|:,time|:)) /)
        psa2_sst_son = sst_son(0,:,:)
        psa2_sst_son = (/ regCoef(psa2_pc_son,sst_son(lat|:,lon|:,time|:)) /)
     end if  
     if (.not.ismissing(sam_son({-85},{5}))) then
        if (sam_son({-85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           sam_son = sam_son*-1.
           sam_pc_son = sam_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              sam_sst_son = sam_sst_son*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa1_son({-62},{270}))) then
        if (psa1_son({-62},{270}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa1_son = psa1_son*-1.
           psa1_pc_son = psa1_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              psa1_sst_son = psa1_sst_son*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa2_son({-52},{150}))) then
        if (psa2_son({-52},{150}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa2_son = psa2_son*-1.
           psa2_pc_son = psa2_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              psa2_sst_son = psa2_sst_son*-1.
           end if
        end if     
     end if
     sam_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     psa1_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     psa2_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(2)))
     copy_VarCoords(npi_ndjfm,sam_pc_son) 
     copy_VarCoords(npi_ndjfm,psa1_pc_son)
     copy_VarCoords(npi_ndjfm,psa2_pc_son)  
     delete([/evecv,pcts/])

     evecv = eofunc(arr_ann_CW({lat|:-20},lon|:,time|:),3,75)
     pcts = eofunc_ts(arr_ann_CW({lat|:-20},lon|:,time|:),evecv,False)  
     sam_pc_ann = dim_standardize(pcts(0,:),0)
     psa1_pc_ann = dim_standardize(pcts(1,:),0)
     psa2_pc_ann = dim_standardize(pcts(2,:),0)     
     sam_ann = arr_ann(0,:,:)
     sam_ann = (/ regCoef(sam_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     psa1_ann = arr_ann(0,:,:)
     psa1_ann = (/ regCoef(psa1_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     psa2_ann = arr_ann(0,:,:)
     psa2_ann = (/ regCoef(psa2_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        sam_sst_ann = sst_ann(0,:,:)
        sam_sst_ann = (/ regCoef(sam_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
        psa1_sst_ann = sst_ann(0,:,:)
        psa1_sst_ann = (/ regCoef(psa1_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
        psa2_sst_ann = sst_ann(0,:,:)
        psa2_sst_ann = (/ regCoef(psa2_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
     end if  
     if (.not.ismissing(sam_ann({-85},{5}))) then
        if (sam_ann({-85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           sam_ann = sam_ann*-1.
           sam_pc_ann = sam_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              sam_sst_ann = sam_sst_ann*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa1_ann({-62},{270}))) then
        if (psa1_ann({-62},{270}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa1_ann = psa1_ann*-1.
           psa1_pc_ann = psa1_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              psa1_sst_ann = psa1_sst_ann*-1.
           end if
        end if
     end if
     if (.not.ismissing(psa2_ann({-52},{150}))) then
        if (psa2_ann({-52},{150}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           psa2_ann = psa2_ann*-1.
           psa2_pc_ann = psa2_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              psa2_sst_ann = psa2_sst_ann*-1.
           end if
        end if     
     end if
     sam_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     psa1_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     psa2_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(2)))
     copy_VarCoords(npi_ndjfm,sam_pc_ann) 
     copy_VarCoords(npi_ndjfm,psa1_pc_ann)
     copy_VarCoords(npi_ndjfm,psa2_pc_ann)  
     delete([/evecv,pcts/])
     
     if (COMPUTE_MODES_MON.eq."True") then     
        evecv = eofunc(arr_mon_CW({lat|:-20},lon|:,time|:),3,75)
        pcts = eofunc_ts(arr_mon_CW({lat|:-20},lon|:,time|:),evecv,False)  
        sam_pc_mon = dim_standardize(pcts(0,:),0)
        psa1_pc_mon = dim_standardize(pcts(1,:),0)
        psa2_pc_mon = dim_standardize(pcts(2,:),0)        
        sam_mon = arr(0,:,:)
        sam_mon = (/ regCoef(sam_pc_mon,arr(lat|:,lon|:,time|:)) /)
        psa1_mon = arr(0,:,:)
        psa1_mon = (/ regCoef(psa1_pc_mon,arr(lat|:,lon|:,time|:)) /)
        psa2_mon = arr(0,:,:)
        psa2_mon = (/ regCoef(psa2_pc_mon,arr(lat|:,lon|:,time|:)) /)
        if (sstreg_plot_flag.eq.0) then
           sam_sst_mon = sst(0,:,:)
           sam_sst_mon = (/ regCoef(sam_pc_mon,sst(lat|:,lon|:,time|:)) /)
           psa1_sst_mon = sst(0,:,:)
           psa1_sst_mon = (/ regCoef(psa1_pc_mon,sst(lat|:,lon|:,time|:)) /)
           psa2_sst_mon = sst(0,:,:)
           psa2_sst_mon = (/ regCoef(psa2_pc_mon,sst(lat|:,lon|:,time|:)) /)
        end if  
        if (.not.ismissing(sam_mon({-85},{5}))) then
           if (sam_mon({-85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              sam_mon = sam_mon*-1.
              sam_pc_mon = sam_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 sam_sst_mon = sam_sst_mon*-1.
              end if
           end if
        end if
        if (.not.ismissing(psa1_mon({-62},{270}))) then
           if (psa1_mon({-62},{270}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              psa1_mon = psa1_mon*-1.
              psa1_pc_mon = psa1_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 psa1_sst_mon = psa1_sst_mon*-1.
              end if
           end if
        end if
        if (.not.ismissing(psa2_mon({-52},{150}))) then
           if (psa2_mon({-52},{150}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              psa2_mon = psa2_mon*-1.
              psa2_pc_mon = psa2_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 psa2_sst_mon = psa2_sst_mon*-1.
              end if
           end if     
        end if
        sam_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
        psa1_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
        psa2_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(2)))
        sam_pc_mon!0 = "time"
        sam_pc_mon&time = arr&time   
        psa1_pc_mon!0 = "time"
        psa1_pc_mon&time = arr&time   
        psa2_pc_mon!0 = "time"
        psa2_pc_mon&time = arr&time   
        delete([/evecv,pcts/])
     end if
;----------PNA/NPO calculations (EOF1/2 of NP PSL)----------------------------------------------------------
     evecv = eofunc(arr_djf_CW({lat|20:85},{lon|120:240},time|:),2,75)
     pcts = eofunc_ts(arr_djf_CW({lat|20:85},{lon|120:240},time|:),evecv,False)  
     pna_pc_djf = dim_standardize(pcts(0,:),0)
     npo_pc_djf = dim_standardize(pcts(1,:),0)
     pna_djf = arr_djf(0,:,:)
     pna_djf = (/ regCoef(pna_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     npo_djf = arr_djf(0,:,:)
     npo_djf = (/ regCoef(npo_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        pna_sst_djf = sst_djf(0,:,:)
        pna_sst_djf = (/ regCoef(pna_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
        npo_sst_djf = sst_djf(0,:,:)
        npo_sst_djf = (/ regCoef(npo_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
     end if     
     if (.not.ismissing(pna_djf({50},{185}))) then
        if (pna_djf({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_djf = pna_djf*-1.
           pna_pc_djf = pna_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_djf = pna_sst_djf*-1.
           end if
        end if
     end if
     if (.not.ismissing(npo_djf({65},{185}))) then
        if (npo_djf({65},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           npo_djf = npo_djf*-1.
           npo_pc_djf = npo_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              npo_sst_djf = npo_sst_djf*-1.
           end if
        end if
     end if
     pna_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     npo_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     copy_VarCoords(npi_ndjfm,pna_pc_djf) 
     copy_VarCoords(npi_ndjfm,npo_pc_djf)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_mam_CW({lat|20:85},{lon|120:240},time|:),2,75)
     pcts = eofunc_ts(arr_mam_CW({lat|20:85},{lon|120:240},time|:),evecv,False)  
     pna_pc_mam = dim_standardize(pcts(0,:),0)
     npo_pc_mam = dim_standardize(pcts(1,:),0)
     pna_mam = arr_mam(0,:,:)
     pna_mam = (/ regCoef(pna_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     npo_mam = arr_mam(0,:,:)
     npo_mam = (/ regCoef(npo_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        pna_sst_mam = sst_mam(0,:,:)
        pna_sst_mam = (/ regCoef(pna_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
        npo_sst_mam = sst_mam(0,:,:)
        npo_sst_mam = (/ regCoef(npo_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
     end if   
     if (.not.ismissing(pna_mam({50},{185}))) then
        if (pna_mam({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_mam = pna_mam*-1.
           pna_pc_mam = pna_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_mam = pna_sst_mam*-1.
           end if
        end if
     end if
     if (.not.ismissing(npo_mam({65},{185}))) then
        if (npo_mam({65},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           npo_mam = npo_mam*-1.
           npo_pc_mam = npo_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              npo_sst_mam = npo_sst_mam*-1.
           end if
        end if
     end if
     pna_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     npo_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     copy_VarCoords(npi_ndjfm,pna_pc_mam) 
     copy_VarCoords(npi_ndjfm,npo_pc_mam)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_jja_CW({lat|20:85},{lon|120:240},time|:),2,75)
     pcts = eofunc_ts(arr_jja_CW({lat|20:85},{lon|120:240},time|:),evecv,False)  
     pna_pc_jja = dim_standardize(pcts(0,:),0)
     npo_pc_jja = dim_standardize(pcts(1,:),0)
     pna_jja = arr_jja(0,:,:)
     pna_jja = (/ regCoef(pna_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     npo_jja = arr_jja(0,:,:)
     npo_jja = (/ regCoef(npo_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        pna_sst_jja = sst_jja(0,:,:)
        pna_sst_jja = (/ regCoef(pna_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
        npo_sst_jja = sst_jja(0,:,:)
        npo_sst_jja = (/ regCoef(npo_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
     end if   
     if (.not.ismissing(pna_jja({50},{185}))) then
        if (pna_jja({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_jja = pna_jja*-1.
           pna_pc_jja = pna_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_jja = pna_sst_jja*-1.
           end if
        end if
     end if
     if (.not.ismissing(npo_jja({65},{185}))) then
        if (npo_jja({65},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           npo_jja = npo_jja*-1.
           npo_pc_jja = npo_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              npo_sst_jja = npo_sst_jja*-1.
           end if
        end if
     end if
     pna_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     npo_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     copy_VarCoords(npi_ndjfm,pna_pc_jja) 
     copy_VarCoords(npi_ndjfm,npo_pc_jja)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_son_CW({lat|20:85},{lon|120:240},time|:),2,75)
     pcts = eofunc_ts(arr_son_CW({lat|20:85},{lon|120:240},time|:),evecv,False)  
     pna_pc_son = dim_standardize(pcts(0,:),0)
     npo_pc_son = dim_standardize(pcts(1,:),0)
     pna_son = arr_son(0,:,:)
     pna_son = (/ regCoef(pna_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     npo_son = arr_son(0,:,:)
     npo_son = (/ regCoef(npo_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        pna_sst_son = sst_son(0,:,:)
        pna_sst_son = (/ regCoef(pna_pc_son,sst_son(lat|:,lon|:,time|:)) /)
        npo_sst_son = sst_son(0,:,:)
        npo_sst_son = (/ regCoef(npo_pc_son,sst_son(lat|:,lon|:,time|:)) /)
     end if      
     if (.not.ismissing(pna_son({50},{185}))) then
        if (pna_son({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_son = pna_son*-1.
           pna_pc_son = pna_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_son = pna_sst_son*-1.
           end if
        end if
     end if
     if (.not.ismissing(npo_son({65},{185}))) then
        if (npo_son({65},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           npo_son = npo_son*-1.
           npo_pc_son = npo_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              npo_sst_son = npo_sst_son*-1.
           end if
        end if
     end if
     pna_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     npo_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     copy_VarCoords(npi_ndjfm,pna_pc_son) 
     copy_VarCoords(npi_ndjfm,npo_pc_son)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_ann_CW({lat|20:85},{lon|120:240},time|:),2,75)
     pcts = eofunc_ts(arr_ann_CW({lat|20:85},{lon|120:240},time|:),evecv,False)  
     pna_pc_ann = dim_standardize(pcts(0,:),0)
     npo_pc_ann = dim_standardize(pcts(1,:),0)
     pna_ann = arr_ann(0,:,:)
     pna_ann = (/ regCoef(pna_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     npo_ann = arr_ann(0,:,:)
     npo_ann = (/ regCoef(npo_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        pna_sst_ann = sst_ann(0,:,:)
        pna_sst_ann = (/ regCoef(pna_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
        npo_sst_ann = sst_ann(0,:,:)
        npo_sst_ann = (/ regCoef(npo_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
     end if 
     if (.not.ismissing(pna_ann({50},{185}))) then
        if (pna_ann({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_ann = pna_ann*-1.
           pna_pc_ann = pna_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_ann = pna_sst_ann*-1.
           end if
        end if
     end if
     if (.not.ismissing(npo_ann({65},{185}))) then
        if (npo_ann({65},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           npo_ann = npo_ann*-1.
           npo_pc_ann = npo_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              npo_sst_ann = npo_sst_ann*-1.
           end if
        end if
     end if
     pna_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     npo_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
     copy_VarCoords(npi_ndjfm,pna_pc_ann) 
     copy_VarCoords(npi_ndjfm,npo_pc_ann)
     delete([/evecv,pcts/])
     
     if (COMPUTE_MODES_MON.eq."True") then          
        evecv = eofunc(arr_mon_CW({lat|20:85},{lon|120:240},time|:),2,75)
        pcts = eofunc_ts(arr_mon_CW({lat|20:85},{lon|120:240},time|:),evecv,False)  
        pna_pc_mon = dim_standardize(pcts(0,:),0)
        npo_pc_mon = dim_standardize(pcts(1,:),0)
        pna_mon = arr(0,:,:)
        pna_mon = (/ regCoef(pna_pc_mon,arr(lat|:,lon|:,time|:)) /)
        npo_mon = arr(0,:,:)
        npo_mon = (/ regCoef(npo_pc_mon,arr(lat|:,lon|:,time|:)) /)
        if (sstreg_plot_flag.eq.0) then
           pna_sst_mon = sst(0,:,:)
           pna_sst_mon = (/ regCoef(pna_pc_mon,sst(lat|:,lon|:,time|:)) /)
           npo_sst_mon = sst(0,:,:)
           npo_sst_mon = (/ regCoef(npo_pc_mon,sst(lat|:,lon|:,time|:)) /)
        end if 
        if (.not.ismissing(pna_mon({50},{185}))) then
           if (pna_mon({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              pna_mon = pna_mon*-1.
              pna_pc_mon = pna_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 pna_sst_mon = pna_sst_mon*-1.
              end if
           end if
        end if
        if (.not.ismissing(npo_mon({65},{185}))) then
           if (npo_mon({65},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              npo_mon = npo_mon*-1.
              npo_pc_mon = npo_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 npo_sst_mon = npo_sst_mon*-1.
              end if
           end if
        end if
        pna_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
        npo_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))
        pna_pc_mon!0 = "time"
        pna_pc_mon&time = arr&time   
        npo_pc_mon!0 = "time"
        npo_pc_mon&time = arr&time   
        delete([/evecv,pcts/])    
     end if
;----------NAO calculations--------------------------------------------------------------------------------     
     arr_djf_CW_LF = lonFlip(arr_djf_CW)
     arr_mam_CW_LF = lonFlip(arr_mam_CW)
     arr_jja_CW_LF = lonFlip(arr_jja_CW)
     arr_son_CW_LF = lonFlip(arr_son_CW)
     arr_ann_CW_LF = lonFlip(arr_ann_CW)
     if (COMPUTE_MODES_MON.eq."True") then   
        arr_mon_CW_LF = lonFlip(arr_mon_CW)
        delete(arr_mon_CW)
     end if
     delete([/arr_djf_CW,arr_mam_CW,arr_jja_CW,arr_son_CW,arr_ann_CW/])
         
     evecv = eofunc(arr_djf_CW_LF({lat|20:80},{lon|-90.:40},time|:),1,75)
     pcts = eofunc_ts(arr_djf_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_djf = dim_standardize(pcts(0,:),0)
     nao_djf = arr_djf(0,:,:)
     nao_djf = (/ regCoef(nao_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_djf = sst_djf(0,:,:)
        nao_sst_djf = (/ regCoef(nao_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
        delete(sst_djf)
     end if
     if (.not.ismissing(nao_djf({70},{350}))) then
        if (nao_djf({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_djf = nao_djf*-1.
           nao_pc_djf = nao_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_djf = nao_sst_djf*-1.
           end if
        end if
     end if
     nao_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nao_pc_djf) 
     delete([/evecv,pcts,arr_djf,arr_djf_CW_LF/])
     
     evecv = eofunc(arr_mam_CW_LF({lat|20:80},{lon|-90.:40},time|:),1,75)
     pcts = eofunc_ts(arr_mam_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_mam = dim_standardize(pcts(0,:),0)
     nao_mam = arr_mam(0,:,:)
     nao_mam = (/ regCoef(nao_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_mam = sst_mam(0,:,:)
        nao_sst_mam = (/ regCoef(nao_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
        delete(sst_mam)
     end if
     if (.not.ismissing(nao_mam({70},{350}))) then
        if (nao_mam({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_mam = nao_mam*-1.
           nao_pc_mam = nao_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_mam = nao_sst_mam*-1.
           end if
        end if
     end if
     nao_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nao_pc_mam) 
     delete([/evecv,pcts,arr_mam,arr_mam_CW_LF/])
     
     evecv = eofunc(arr_jja_CW_LF({lat|20:80},{lon|-90.:40},time|:),1,75)
     pcts = eofunc_ts(arr_jja_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_jja = dim_standardize(pcts(0,:),0)
     nao_jja = arr_jja(0,:,:)
     nao_jja = (/ regCoef(nao_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_jja = sst_jja(0,:,:)
        nao_sst_jja = (/ regCoef(nao_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
        delete(sst_jja)
     end if
     if (.not.ismissing(nao_jja({70},{350}))) then
        if (nao_jja({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_jja = nao_jja*-1.
           nao_pc_jja = nao_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_jja = nao_sst_jja*-1.
           end if
        end if
     end if
     nao_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nao_pc_jja) 
     delete([/evecv,pcts,arr_jja,arr_jja_CW_LF/])
     
     evecv = eofunc(arr_son_CW_LF({lat|20:80},{lon|-90.:40},time|:),1,75)
     pcts = eofunc_ts(arr_son_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_son = dim_standardize(pcts(0,:),0)
     nao_son = arr_son(0,:,:)
     nao_son = (/ regCoef(nao_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_son = sst_son(0,:,:)
        nao_sst_son = (/ regCoef(nao_pc_son,sst_son(lat|:,lon|:,time|:)) /)
        delete(sst_son)
     end if
     if (.not.ismissing(nao_son({70},{350}))) then
        if (nao_son({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_son = nao_son*-1.
           nao_pc_son = nao_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_son = nao_sst_son*-1.
           end if
        end if
     end if
     nao_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nao_pc_son) 
     delete([/evecv,pcts,arr_son,arr_son_CW_LF/])
     
     evecv = eofunc(arr_ann_CW_LF({lat|20:80},{lon|-90.:40},time|:),1,75)
     pcts = eofunc_ts(arr_ann_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_ann = dim_standardize(pcts(0,:),0)
     nao_ann = arr_ann(0,:,:)
     nao_ann = (/ regCoef(nao_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_ann = sst_ann(0,:,:)
        nao_sst_ann = (/ regCoef(nao_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
        delete(sst_ann)
     end if
     if (.not.ismissing(nao_ann({70},{350}))) then
        if (nao_ann({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_ann = nao_ann*-1.
           nao_pc_ann = nao_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_ann = nao_sst_ann*-1.
           end if
        end if
     end if
     nao_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
     copy_VarCoords(npi_ndjfm,nao_pc_ann) 
     delete([/evecv,pcts,arr_ann,arr_ann_CW_LF/])

     if (COMPUTE_MODES_MON.eq."True") then               
        evecv = eofunc(arr_mon_CW_LF({lat|20:80},{lon|-90.:40},time|:),1,75)
        pcts = eofunc_ts(arr_mon_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
        nao_pc_mon = dim_standardize(pcts(0,:),0)
        nao_mon = arr(0,:,:)
        nao_mon = (/ regCoef(nao_pc_mon,arr(lat|:,lon|:,time|:)) /)
        if (sstreg_plot_flag.eq.0) then
           nao_sst_mon = sst(0,:,:)
           nao_sst_mon = (/ regCoef(nao_pc_mon,sst(lat|:,lon|:,time|:)) /)
           delete(sst)
        end if
        if (.not.ismissing(nao_mon({70},{350}))) then
           if (nao_mon({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              nao_mon = nao_mon*-1.
              nao_pc_mon = nao_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 nao_sst_mon = nao_sst_mon*-1.
              end if
           end if
        end if
        nao_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))
        nao_pc_mon!0 = "time"
        nao_pc_mon&time = arr&time   
        delete([/evecv,pcts,arr,arr_mon_CW_LF/])
     end if
;------------------------------------------------------------------------------------------------------
     if (sstreg_frame.eq.1.and.sstreg_plot_flag.eq.0) then    ; sstreg_frame = flag to create regressions .ps/.png files 
        sstreg_frame = 0
     end if          
;------------------------------------------------------------------------------------------------------
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
        z->npi_ndjfm = npi_ndjfm
        
        z->nao_pc_djf = nao_pc_djf
        z->nao_pc_mam = nao_pc_mam
        z->nao_pc_jja = nao_pc_jja
        z->nao_pc_son = nao_pc_son
        z->nao_pc_ann = nao_pc_ann
        
        z->sam_pc_djf = sam_pc_djf
        z->sam_pc_mam = sam_pc_mam
        z->sam_pc_jja = sam_pc_jja
        z->sam_pc_son = sam_pc_son
        z->sam_pc_ann = sam_pc_ann
        
        z->psa1_pc_djf = psa1_pc_djf
        z->psa1_pc_mam = psa1_pc_mam
        z->psa1_pc_jja = psa1_pc_jja
        z->psa1_pc_son = psa1_pc_son
        z->psa1_pc_ann = psa1_pc_ann
        
        z->psa2_pc_djf = psa2_pc_djf
        z->psa2_pc_mam = psa2_pc_mam
        z->psa2_pc_jja = psa2_pc_jja
        z->psa2_pc_son = psa2_pc_son
        z->psa2_pc_ann = psa2_pc_ann
        
        z->nam_pc_djf = nam_pc_djf
        z->nam_pc_mam = nam_pc_mam
        z->nam_pc_jja = nam_pc_jja
        z->nam_pc_son = nam_pc_son
        z->nam_pc_ann = nam_pc_ann
        
        z->pna_pc_djf = pna_pc_djf
        z->pna_pc_mam = pna_pc_mam
        z->pna_pc_jja = pna_pc_jja
        z->pna_pc_son = pna_pc_son
        z->pna_pc_ann = pna_pc_ann
        
        z->npo_pc_djf = npo_pc_djf
        z->npo_pc_mam = npo_pc_mam
        z->npo_pc_jja = npo_pc_jja
        z->npo_pc_son = npo_pc_son
        z->npo_pc_ann = npo_pc_ann
        
        z->nao_djf = nao_djf
        z->nao_mam = nao_mam
        z->nao_jja = nao_jja
        z->nao_son = nao_son
        z->nao_ann = nao_ann
        
        z->sam_djf = sam_djf
        z->sam_mam = sam_mam
        z->sam_jja = sam_jja
        z->sam_son = sam_son
        z->sam_ann = sam_ann
        
        z->psa1_djf = psa1_djf
        z->psa1_mam = psa1_mam
        z->psa1_jja = psa1_jja
        z->psa1_son = psa1_son
        z->psa1_ann = psa1_ann
        
        z->psa2_djf = psa2_djf
        z->psa2_mam = psa2_mam
        z->psa2_jja = psa2_jja
        z->psa2_son = psa2_son
        z->psa2_ann = psa2_ann
        
        z->nam_djf = nam_djf
        z->nam_mam = nam_mam
        z->nam_jja = nam_jja
        z->nam_son = nam_son
        z->nam_ann = nam_ann
        
        z->pna_djf = pna_djf
        z->pna_mam = pna_mam
        z->pna_jja = pna_jja
        z->pna_son = pna_son
        z->pna_ann = pna_ann
        
        z->npo_djf = npo_djf
        z->npo_mam = npo_mam
        z->npo_jja = npo_jja
        z->npo_son = npo_son
        z->npo_ann = npo_ann

        if (isfilepresent2("obs_psl").and.ee.eq.0) then   ; will not write out SST regressions to OBS files due to different
        else                                             ; variable going into PSL OBS file..
           if (sstreg_plot_flag.eq.0) then
              z->nao_sst_regression_djf = nao_sst_djf
              z->nao_sst_regression_mam = nao_sst_mam
              z->nao_sst_regression_jja = nao_sst_jja
              z->nao_sst_regression_son = nao_sst_son
              z->nao_sst_regression_ann = nao_sst_ann
        
              z->sam_sst_regression_djf = sam_sst_djf
              z->sam_sst_regression_mam = sam_sst_mam
              z->sam_sst_regression_jja = sam_sst_jja
              z->sam_sst_regression_son = sam_sst_son
              z->sam_sst_regression_ann = sam_sst_ann
        
              z->psa1_sst_regression_djf = psa1_sst_djf
              z->psa1_sst_regression_mam = psa1_sst_mam
              z->psa1_sst_regression_jja = psa1_sst_jja
              z->psa1_sst_regression_son = psa1_sst_son
              z->psa1_sst_regression_ann = psa1_sst_ann
        
              z->psa2_sst_regression_djf = psa2_sst_djf
              z->psa2_sst_regression_mam = psa2_sst_mam
              z->psa2_sst_regression_jja = psa2_sst_jja
              z->psa2_sst_regression_son = psa2_sst_son
              z->psa2_sst_regression_ann = psa2_sst_ann
        
              z->nam_sst_regression_djf = nam_sst_djf
              z->nam_sst_regression_mam = nam_sst_mam
              z->nam_sst_regression_jja = nam_sst_jja
              z->nam_sst_regression_son = nam_sst_son
              z->nam_sst_regression_ann = nam_sst_ann
        
              z->pna_sst_regression_djf = pna_sst_djf
              z->pna_sst_regression_mam = pna_sst_mam
              z->pna_sst_regression_jja = pna_sst_jja
              z->pna_sst_regression_son = pna_sst_son
              z->pna_sst_regression_ann = pna_sst_ann
        
              z->npo_sst_regression_djf = npo_sst_djf
              z->npo_sst_regression_mam = npo_sst_mam
              z->npo_sst_regression_jja = npo_sst_jja
              z->npo_sst_regression_son = npo_sst_son
              z->npo_sst_regression_ann = npo_sst_ann
           end if
        end if
        
        if (COMPUTE_MODES_MON.eq."True") then   
           z->nao_pc_mon = nao_pc_mon
           z->sam_pc_mon = sam_pc_mon
           z->psa1_pc_mon = psa1_pc_mon
           z->psa2_pc_mon = psa2_pc_mon
           z->nam_pc_mon = nam_pc_mon
           z->pna_pc_mon = pna_pc_mon
           z->npo_pc_mon = npo_pc_mon
           z->nao_mon = nao_mon
           z->sam_mon = sam_mon
           z->psa1_mon = psa1_mon
           z->psa2_mon = psa2_mon
           z->nam_mon = nam_mon
           z->pna_mon = pna_mon
           z->npo_mon = npo_mon
           if (isfilepresent2("obs_psl").and.ee.eq.0) then   ; will not write out SST regressions to OBS files due to different
           else                                             ; variable going into PSL OBS file..
              if (sstreg_plot_flag.eq.0) then 
                 z->nao_sst_regression_mon = nao_sst_mon
                 z->sam_sst_regression_mon = sam_sst_mon
                 z->psa1_sst_regression_mon = psa1_sst_mon
                 z->psa2_sst_regression_mon = psa2_sst_mon
                 z->nam_sst_regression_mon = nam_sst_mon
                 z->pna_sst_regression_mon = pna_sst_mon
                 z->npo_sst_regression_mon = npo_sst_mon
              end if
           end if
        end if
        delete(z)
     end if     
     
;========================================================================
     res = True
     res@mpGeophysicalLineColor = "gray42"
     res@mpGeophysicalLineThicknessF = 2.     
     res@mpGridAndLimbOn = False
     res@mpFillOn = False
     res@mpOutlineOn = True  
     res@mpDataSetName = "Earth..4"
     res@gsnDraw      = False
     res@gsnFrame     = False
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     res@gsnLeftStringOrthogonalPosF = -0.03
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.03
     res@gsnRightStringParallelPosF = 0.96
     res@gsnRightString = ""
     res@gsnLeftString = ""
     if (nsim.le.5) then
        res@gsnLeftStringFontHeightF = 0.018
        res@gsnCenterStringFontHeightF = 0.022
        res@gsnRightStringFontHeightF = 0.018
     else
        res@gsnLeftStringFontHeightF = 0.024
        res@gsnCenterStringFontHeightF = 0.028
        res@gsnRightStringFontHeightF = 0.024     
     end if
     res@gsnPolar = "NH"     
     res@mpMinLatF    = 20.
     res@mpCenterLonF = 0.
     
     res@cnLevels = (/-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7./)
     res@gsnLeftString = syear(ee)+"-"+eyear(ee)
     res@gsnCenterString = names(ee)
     
     if (isfilepresent2("obs_psl").and.ee.eq.0) then    ; for pattern correlation table. Save entire lat/lon array
        patcor_nam_djf = new((/nsim,dimsizes(nam_djf&lat),dimsizes(nam_djf&lon)/),typeof(nam_djf))
        patcor_nam_djf!1 = "lat"
        patcor_nam_djf&lat = nam_djf&lat
        patcor_nam_djf!2 = "lon"
        patcor_nam_djf&lon = nam_djf&lon
        patcor_nam_jja = patcor_nam_djf
        patcor_nam_ann = patcor_nam_djf
        patcor_sam_djf = patcor_nam_djf
        patcor_sam_jja = patcor_nam_djf
        patcor_sam_ann = patcor_nam_djf
        patcor_pna_djf = patcor_nam_djf
        patcor_pna_jja = patcor_nam_djf
        patcor_pna_ann = patcor_nam_djf
        patcor_npo_djf = patcor_nam_djf
        patcor_npo_jja = patcor_nam_djf
        patcor_npo_ann = patcor_nam_djf
        patcor_nao_djf = patcor_nam_djf
        patcor_nao_jja = patcor_nam_djf
        patcor_nao_ann = patcor_nam_djf
        patcor_psa1_djf = patcor_nam_djf
        patcor_psa1_jja = patcor_nam_djf
        patcor_psa1_ann = patcor_nam_djf
        patcor_psa2_djf = patcor_nam_djf
        patcor_psa2_jja = patcor_nam_djf
        patcor_psa2_ann = patcor_nam_djf
        patcor_nam_djf(ee,:,:) = (/ nam_djf /)
        patcor_nam_jja(ee,:,:) = (/ nam_jja /)
        patcor_nam_ann(ee,:,:) = (/ nam_ann /)
        patcor_sam_djf(ee,:,:) = (/ sam_djf /)
        patcor_sam_jja(ee,:,:) = (/ sam_jja /)
        patcor_sam_ann(ee,:,:) = (/ sam_ann /)
        patcor_nao_djf(ee,:,:) = (/ nao_djf /)
        patcor_nao_jja(ee,:,:) = (/ nao_jja /)
        patcor_nao_ann(ee,:,:) = (/ nao_ann /)
        patcor_pna_djf(ee,:,:) = (/ pna_djf /)
        patcor_pna_jja(ee,:,:) = (/ pna_jja /)
        patcor_pna_ann(ee,:,:) = (/ pna_ann /)
        patcor_npo_djf(ee,:,:) = (/ npo_djf /)
        patcor_npo_jja(ee,:,:) = (/ npo_jja /)
        patcor_npo_ann(ee,:,:) = (/ npo_ann /)
        patcor_psa1_djf(ee,:,:) = (/ psa1_djf /)
        patcor_psa1_jja(ee,:,:) = (/ psa1_jja /)
        patcor_psa1_ann(ee,:,:) = (/ psa1_ann /)
        patcor_psa2_djf(ee,:,:) = (/ psa2_djf /)
        patcor_psa2_jja(ee,:,:) = (/ psa2_jja /)
        patcor_psa2_ann(ee,:,:) = (/ psa2_ann /)
     end if
     if (isfilepresent2("obs_psl").and.ee.ge.1) then
        patcor_nam_djf(ee,:,:) = (/ linint2(nam_djf&lon,nam_djf&lat,nam_djf,True,patcor_nam_djf&lon,patcor_nam_djf&lat,0) /)
        patcor_nam_jja(ee,:,:) = (/ linint2(nam_jja&lon,nam_jja&lat,nam_jja,True,patcor_nam_jja&lon,patcor_nam_jja&lat,0) /)
        patcor_nam_ann(ee,:,:) = (/ linint2(nam_ann&lon,nam_ann&lat,nam_ann,True,patcor_nam_ann&lon,patcor_nam_ann&lat,0) /)
        
        patcor_sam_djf(ee,:,:) = (/ linint2(sam_djf&lon,sam_djf&lat,sam_djf,True,patcor_sam_djf&lon,patcor_sam_djf&lat,0) /)
        patcor_sam_jja(ee,:,:) = (/ linint2(sam_jja&lon,sam_jja&lat,sam_jja,True,patcor_sam_jja&lon,patcor_sam_jja&lat,0) /)
        patcor_sam_ann(ee,:,:) = (/ linint2(sam_ann&lon,sam_ann&lat,sam_ann,True,patcor_sam_ann&lon,patcor_sam_ann&lat,0) /)
        
        patcor_nao_djf(ee,:,:) = (/ linint2(nao_djf&lon,nao_djf&lat,nao_djf,True,patcor_nao_djf&lon,patcor_nao_djf&lat,0) /)
        patcor_nao_jja(ee,:,:) = (/ linint2(nao_jja&lon,nao_jja&lat,nao_jja,True,patcor_nao_jja&lon,patcor_nao_jja&lat,0) /)
        patcor_nao_ann(ee,:,:) = (/ linint2(nao_ann&lon,nao_ann&lat,nao_ann,True,patcor_nao_ann&lon,patcor_nao_ann&lat,0) /)
        
        patcor_pna_djf(ee,:,:) = (/ linint2(pna_djf&lon,pna_djf&lat,pna_djf,True,patcor_pna_djf&lon,patcor_pna_djf&lat,0) /)
        patcor_pna_jja(ee,:,:) = (/ linint2(pna_jja&lon,pna_jja&lat,pna_jja,True,patcor_pna_jja&lon,patcor_pna_jja&lat,0) /)
        patcor_pna_ann(ee,:,:) = (/ linint2(pna_ann&lon,pna_ann&lat,pna_ann,True,patcor_pna_ann&lon,patcor_pna_ann&lat,0) /)
        
        patcor_npo_djf(ee,:,:) = (/ linint2(npo_djf&lon,npo_djf&lat,npo_djf,True,patcor_npo_djf&lon,patcor_npo_djf&lat,0) /)
        patcor_npo_jja(ee,:,:) = (/ linint2(npo_jja&lon,npo_jja&lat,npo_jja,True,patcor_npo_jja&lon,patcor_npo_jja&lat,0) /)
        patcor_npo_ann(ee,:,:) = (/ linint2(npo_ann&lon,npo_ann&lat,npo_ann,True,patcor_npo_ann&lon,patcor_npo_ann&lat,0) /)
        
        patcor_psa1_djf(ee,:,:) = (/ linint2(psa1_djf&lon,psa1_djf&lat,psa1_djf,True,patcor_psa1_djf&lon,patcor_psa1_djf&lat,0) /)
        patcor_psa1_jja(ee,:,:) = (/ linint2(psa1_jja&lon,psa1_jja&lat,psa1_jja,True,patcor_psa1_jja&lon,patcor_psa1_jja&lat,0) /)
        patcor_psa1_ann(ee,:,:) = (/ linint2(psa1_ann&lon,psa1_ann&lat,psa1_ann,True,patcor_psa1_ann&lon,patcor_psa1_ann&lat,0) /)
        
        patcor_psa2_djf(ee,:,:) = (/ linint2(psa2_djf&lon,psa2_djf&lat,psa2_djf,True,patcor_psa2_djf&lon,patcor_psa2_djf&lat,0) /)
        patcor_psa2_jja(ee,:,:) = (/ linint2(psa2_jja&lon,psa2_jja&lat,psa2_jja,True,patcor_psa2_jja&lon,patcor_psa2_jja&lat,0) /)
        patcor_psa2_ann(ee,:,:) = (/ linint2(psa2_ann&lon,psa2_ann&lat,psa2_ann,True,patcor_psa2_ann&lon,patcor_psa2_ann&lat,0) /)
     end if
               
     res@gsnRightString = nam_djf@pcvar+"%"
     map_nam_djf(ee) = gsn_csm_contour_map_polar(wks_nam,nam_djf,res)
     res@gsnRightString = nam_mam@pcvar+"%"
     map_nam_mam(ee) = gsn_csm_contour_map_polar(wks_nam,nam_mam,res)
     res@gsnRightString = nam_jja@pcvar+"%"
     map_nam_jja(ee) = gsn_csm_contour_map_polar(wks_nam,nam_jja,res)
     res@gsnRightString = nam_son@pcvar+"%"
     map_nam_son(ee) = gsn_csm_contour_map_polar(wks_nam,nam_son,res)
     res@gsnRightString = nam_ann@pcvar+"%"
     map_nam_ann(ee) = gsn_csm_contour_map_polar(wks_nam,nam_ann,res)     
     delete([/nam_djf,nam_mam,nam_jja,nam_son,nam_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = nam_mon@pcvar+"%"
        map_nam_mon(ee) = gsn_csm_contour_map_polar(wks_nam,nam_mon,res)
        delete([/nam_mon/])
     end if
     
     res@gsnRightString = nao_djf@pcvar+"%"
     map_nao_djf(ee) = gsn_csm_contour_map_polar(wks_nao,nao_djf,res)
     res@gsnRightString = nao_mam@pcvar+"%"
     map_nao_mam(ee) = gsn_csm_contour_map_polar(wks_nao,nao_mam,res)
     res@gsnRightString = nao_jja@pcvar+"%"
     map_nao_jja(ee) = gsn_csm_contour_map_polar(wks_nao,nao_jja,res)
     res@gsnRightString = nao_son@pcvar+"%"
     map_nao_son(ee) = gsn_csm_contour_map_polar(wks_nao,nao_son,res)
     res@gsnRightString = nao_ann@pcvar+"%"
     map_nao_ann(ee) = gsn_csm_contour_map_polar(wks_nao,nao_ann,res)     
     delete([/nao_djf,nao_mam,nao_jja,nao_son,nao_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = nao_mon@pcvar+"%"
        map_nao_mon(ee) = gsn_csm_contour_map_polar(wks_nao,nao_mon,res)
        delete([/nao_mon/])
     end if
     
     res@gsnRightString = pna_djf@pcvar+"%"
     map_pna_djf(ee) = gsn_csm_contour_map_polar(wks_pna,pna_djf,res)
     res@gsnRightString = pna_mam@pcvar+"%"
     map_pna_mam(ee) = gsn_csm_contour_map_polar(wks_pna,pna_mam,res)
     res@gsnRightString = pna_jja@pcvar+"%"
     map_pna_jja(ee) = gsn_csm_contour_map_polar(wks_pna,pna_jja,res)
     res@gsnRightString = pna_son@pcvar+"%"
     map_pna_son(ee) = gsn_csm_contour_map_polar(wks_pna,pna_son,res)
     res@gsnRightString = pna_ann@pcvar+"%"
     map_pna_ann(ee) = gsn_csm_contour_map_polar(wks_pna,pna_ann,res)     
     delete([/pna_djf,pna_mam,pna_jja,pna_son,pna_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = pna_mon@pcvar+"%"
        map_pna_mon(ee) = gsn_csm_contour_map_polar(wks_pna,pna_mon,res)
        delete([/pna_mon/])
     end if
     
     res@cnLevels = (/-4,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4/)
     res@gsnRightString = npo_djf@pcvar+"%"
     map_npo_djf(ee) = gsn_csm_contour_map_polar(wks_npo,npo_djf,res)
     res@gsnRightString = npo_mam@pcvar+"%"
     map_npo_mam(ee) = gsn_csm_contour_map_polar(wks_npo,npo_mam,res)
     res@gsnRightString = npo_jja@pcvar+"%"
     map_npo_jja(ee) = gsn_csm_contour_map_polar(wks_npo,npo_jja,res)
     res@gsnRightString = npo_son@pcvar+"%"
     map_npo_son(ee) = gsn_csm_contour_map_polar(wks_npo,npo_son,res)
     res@gsnRightString = npo_ann@pcvar+"%"
     map_npo_ann(ee) = gsn_csm_contour_map_polar(wks_npo,npo_ann,res)     
     delete([/npo_djf,npo_mam,npo_jja,npo_son,npo_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = npo_mon@pcvar+"%"
        map_npo_mon(ee) = gsn_csm_contour_map_polar(wks_npo,npo_mon,res)
        delete([/npo_mon/])
     end if
     
     if (sstreg_plot_flag.eq.0) then
        res@cnLevels = fspan(-1.05,1.05,15)
        res@gsnRightString = ""
        reg_nam_djf(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_djf,res)
        reg_nam_mam(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_mam,res)
        reg_nam_jja(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_jja,res)
        reg_nam_son(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_son,res)
        reg_nam_ann(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_ann,res)     
        delete([/nam_sst_djf,nam_sst_mam,nam_sst_jja,nam_sst_son,nam_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_nam_mon(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_mon,res)
           delete([/nam_sst_mon/])
        end if
     
        reg_nao_djf(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_djf,res)
        reg_nao_mam(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_mam,res)
        reg_nao_jja(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_jja,res)
        reg_nao_son(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_son,res)
        reg_nao_ann(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_ann,res)     
        delete([/nao_sst_djf,nao_sst_mam,nao_sst_jja,nao_sst_son,nao_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_nao_mon(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_mon,res)
           delete([/nao_sst_mon/])
        end if
     
        reg_pna_djf(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_djf,res)
        reg_pna_mam(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_mam,res)
        reg_pna_jja(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_jja,res)
        reg_pna_son(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_son,res)
        reg_pna_ann(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_ann,res)     
        delete([/pna_sst_djf,pna_sst_mam,pna_sst_jja,pna_sst_son,pna_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_pna_mon(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_mon,res)
           delete([/pna_sst_mon/])
        end if
     
        reg_npo_djf(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_djf,res)
        reg_npo_mam(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_mam,res)
        reg_npo_jja(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_jja,res)
        reg_npo_son(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_son,res)
        reg_npo_ann(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_ann,res)     
        delete([/npo_sst_djf,npo_sst_mam,npo_sst_jja,npo_sst_son,npo_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_npo_mon(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_mon,res)
           delete([/npo_sst_mon/])
        end if  
     end if
     
     res@cnLevels = (/-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7./)
     res@gsnPolar = "SH"     
     delete(res@mpMinLatF)
     res@mpMaxLatF    = -20.
     res@mpCenterLonF = 0.
     res@gsnLeftString = syear(ee)+"-"+eyear(ee)
     res@gsnCenterString = names(ee)
     res@gsnRightString = sam_djf@pcvar+"%"
     map_sam_djf(ee) = gsn_csm_contour_map_polar(wks_sam,sam_djf,res)
     res@gsnRightString = sam_mam@pcvar+"%"
     map_sam_mam(ee) = gsn_csm_contour_map_polar(wks_sam,sam_mam,res)
     res@gsnRightString = sam_jja@pcvar+"%"
     map_sam_jja(ee) = gsn_csm_contour_map_polar(wks_sam,sam_jja,res)
     res@gsnRightString = sam_son@pcvar+"%"
     map_sam_son(ee) = gsn_csm_contour_map_polar(wks_sam,sam_son,res)
     res@gsnRightString = sam_ann@pcvar+"%"
     map_sam_ann(ee) = gsn_csm_contour_map_polar(wks_sam,sam_ann,res)     
     delete([/sam_djf,sam_mam,sam_jja,sam_son,sam_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = sam_mon@pcvar+"%"
        map_sam_mon(ee) = gsn_csm_contour_map_polar(wks_sam,sam_mon,res)
        delete([/sam_mon/])
     end if
     
     res@cnLevels = (/-4,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4/)
     res@gsnRightString = psa1_djf@pcvar+"%"
     map_psa1_djf(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_djf,res)
     res@gsnRightString = psa1_mam@pcvar+"%"
     map_psa1_mam(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_mam,res)
     res@gsnRightString = psa1_jja@pcvar+"%"
     map_psa1_jja(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_jja,res)
     res@gsnRightString = psa1_son@pcvar+"%"
     map_psa1_son(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_son,res)
     res@gsnRightString = psa1_ann@pcvar+"%"
     map_psa1_ann(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_ann,res)     
     delete([/psa1_djf,psa1_mam,psa1_jja,psa1_son,psa1_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = psa1_mon@pcvar+"%"
        map_psa1_mon(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_mon,res)
        delete([/psa1_mon/])
     end if
     
     res@gsnRightString = psa2_djf@pcvar+"%"
     map_psa2_djf(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_djf,res)
     res@gsnRightString = psa2_mam@pcvar+"%"
     map_psa2_mam(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_mam,res)
     res@gsnRightString = psa2_jja@pcvar+"%"
     map_psa2_jja(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_jja,res)
     res@gsnRightString = psa2_son@pcvar+"%"
     map_psa2_son(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_son,res)
     res@gsnRightString = psa2_ann@pcvar+"%"
     map_psa2_ann(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_ann,res)     
     delete([/psa2_djf,psa2_mam,psa2_jja,psa2_son,psa2_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = psa2_mon@pcvar+"%"
        map_psa2_mon(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_mon,res)
        delete([/psa2_mon/])
     end if
     
     if (sstreg_plot_flag.eq.0) then
        res@cnLevels = fspan(-1.05,1.05,15)
        res@gsnRightString = ""
        reg_sam_djf(ee) = gsn_csm_contour_map_polar(wks_sam,sam_sst_djf,res)
        reg_sam_mam(ee) = gsn_csm_contour_map_polar(wks_sam,sam_sst_mam,res)
        reg_sam_jja(ee) = gsn_csm_contour_map_polar(wks_sam,sam_sst_jja,res)
        reg_sam_son(ee) = gsn_csm_contour_map_polar(wks_sam,sam_sst_son,res)
        reg_sam_ann(ee) = gsn_csm_contour_map_polar(wks_sam,sam_sst_ann,res)     
        delete([/sam_sst_djf,sam_sst_mam,sam_sst_jja,sam_sst_son,sam_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_sam_mon(ee) = gsn_csm_contour_map_polar(wks_sam,sam_sst_mon,res)
           delete([/sam_sst_mon/])
        end if
     
        reg_psa1_djf(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_sst_djf,res)
        reg_psa1_mam(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_sst_mam,res)
        reg_psa1_jja(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_sst_jja,res)
        reg_psa1_son(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_sst_son,res)
        reg_psa1_ann(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_sst_ann,res)     
        delete([/psa1_sst_djf,psa1_sst_mam,psa1_sst_jja,psa1_sst_son,psa1_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_psa1_mon(ee) = gsn_csm_contour_map_polar(wks_psa1,psa1_sst_mon,res)
           delete([/psa1_sst_mon/])
        end if
     
        reg_psa2_djf(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_sst_djf,res)
        reg_psa2_mam(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_sst_mam,res)
        reg_psa2_jja(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_sst_jja,res)
        reg_psa2_son(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_sst_son,res)
        reg_psa2_ann(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_sst_ann,res)     
        delete([/psa2_sst_djf,psa2_sst_mam,psa2_sst_jja,psa2_sst_son,psa2_sst_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_psa2_mon(ee) = gsn_csm_contour_map_polar(wks_psa2,psa2_sst_mon,res)
           delete([/psa2_sst_mon/])
        end if
     end if
     delete(res@mpMaxLatF)

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
     
     xyresmon = xyres
     xyresmon@gsnXYBarChart = False            
     xyresmon@xyLineThicknessF = .1 
     
     xy_nam_djf(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee),dimsizes(nam_pc_djf)),nam_pc_djf,xyres)   ; use standardized timeseries  
     xy_nam_mam(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee),dimsizes(nam_pc_mam)),nam_pc_mam,xyres)   ; use standardized timeseries  
     xy_nam_jja(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee),dimsizes(nam_pc_jja)),nam_pc_jja,xyres)   ; use standardized timeseries  
     xy_nam_son(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee),dimsizes(nam_pc_son)),nam_pc_son,xyres)   ; use standardized timeseries  
     xy_nam_ann(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee),dimsizes(nam_pc_ann)),nam_pc_ann,xyres)   ; use standardized timeseries  
     delete([/nam_pc_djf,nam_pc_mam,nam_pc_jja,nam_pc_son,nam_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then  
        xy_nam_mon(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee),dimsizes(nam_pc_mon)),nam_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/nam_pc_mon/])
     end if

     xy_nao_djf(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_djf)),nao_pc_djf,xyres)   ; use standardized timeseries  
     xy_nao_mam(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_mam)),nao_pc_mam,xyres)   ; use standardized timeseries  
     xy_nao_jja(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_jja)),nao_pc_jja,xyres)   ; use standardized timeseries  
     xy_nao_son(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_son)),nao_pc_son,xyres)   ; use standardized timeseries  
     xy_nao_ann(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_ann)),nao_pc_ann,xyres)   ; use standardized timeseries  
     delete([/nao_pc_djf,nao_pc_mam,nao_pc_jja,nao_pc_son,nao_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then     
        xy_nao_mon(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_mon)),nao_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/nao_pc_mon/])
     end if

     xy_pna_djf(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_djf)),pna_pc_djf,xyres)   ; use standardized timeseries  
     xy_pna_mam(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_mam)),pna_pc_mam,xyres)   ; use standardized timeseries  
     xy_pna_jja(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_jja)),pna_pc_jja,xyres)   ; use standardized timeseries  
     xy_pna_son(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_son)),pna_pc_son,xyres)   ; use standardized timeseries  
     xy_pna_ann(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_ann)),pna_pc_ann,xyres)   ; use standardized timeseries  
     delete([/pna_pc_djf,pna_pc_mam,pna_pc_jja,pna_pc_son,pna_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then
        xy_pna_mon(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_mon)),pna_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/pna_pc_mon/])
     end if

     xy_npo_djf(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_djf)),npo_pc_djf,xyres)   ; use standardized timeseries  
     xy_npo_mam(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_mam)),npo_pc_mam,xyres)   ; use standardized timeseries  
     xy_npo_jja(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_jja)),npo_pc_jja,xyres)   ; use standardized timeseries  
     xy_npo_son(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_son)),npo_pc_son,xyres)   ; use standardized timeseries  
     xy_npo_ann(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_ann)),npo_pc_ann,xyres)   ; use standardized timeseries  
     delete([/npo_pc_djf,npo_pc_mam,npo_pc_jja,npo_pc_son,npo_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then    
        xy_npo_mon(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_mon)),npo_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/npo_pc_mon/])
     end if

     xy_sam_djf(ee) = gsn_csm_xy(wks_sam_ts,fspan(syear(ee),eyear(ee),dimsizes(sam_pc_djf)),sam_pc_djf,xyres)   ; use standardized timeseries  
     xy_sam_mam(ee) = gsn_csm_xy(wks_sam_ts,fspan(syear(ee),eyear(ee),dimsizes(sam_pc_mam)),sam_pc_mam,xyres)   ; use standardized timeseries  
     xy_sam_jja(ee) = gsn_csm_xy(wks_sam_ts,fspan(syear(ee),eyear(ee),dimsizes(sam_pc_jja)),sam_pc_jja,xyres)   ; use standardized timeseries  
     xy_sam_son(ee) = gsn_csm_xy(wks_sam_ts,fspan(syear(ee),eyear(ee),dimsizes(sam_pc_son)),sam_pc_son,xyres)   ; use standardized timeseries  
     xy_sam_ann(ee) = gsn_csm_xy(wks_sam_ts,fspan(syear(ee),eyear(ee),dimsizes(sam_pc_ann)),sam_pc_ann,xyres)   ; use standardized timeseries  
     delete([/sam_pc_djf,sam_pc_mam,sam_pc_jja,sam_pc_son,sam_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then   
        xy_sam_mon(ee) = gsn_csm_xy(wks_sam_ts,fspan(syear(ee),eyear(ee),dimsizes(sam_pc_mon)),sam_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/sam_pc_mon/])
     end if

     xy_psa1_djf(ee) = gsn_csm_xy(wks_psa1_ts,fspan(syear(ee),eyear(ee),dimsizes(psa1_pc_djf)),psa1_pc_djf,xyres)   ; use standardized timeseries  
     xy_psa1_mam(ee) = gsn_csm_xy(wks_psa1_ts,fspan(syear(ee),eyear(ee),dimsizes(psa1_pc_mam)),psa1_pc_mam,xyres)   ; use standardized timeseries  
     xy_psa1_jja(ee) = gsn_csm_xy(wks_psa1_ts,fspan(syear(ee),eyear(ee),dimsizes(psa1_pc_jja)),psa1_pc_jja,xyres)   ; use standardized timeseries  
     xy_psa1_son(ee) = gsn_csm_xy(wks_psa1_ts,fspan(syear(ee),eyear(ee),dimsizes(psa1_pc_son)),psa1_pc_son,xyres)   ; use standardized timeseries  
     xy_psa1_ann(ee) = gsn_csm_xy(wks_psa1_ts,fspan(syear(ee),eyear(ee),dimsizes(psa1_pc_ann)),psa1_pc_ann,xyres)   ; use standardized timeseries  
     delete([/psa1_pc_djf,psa1_pc_mam,psa1_pc_jja,psa1_pc_son,psa1_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then   
        xy_psa1_mon(ee) = gsn_csm_xy(wks_psa1_ts,fspan(syear(ee),eyear(ee),dimsizes(psa1_pc_mon)),psa1_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/psa1_pc_mon/])
     end if

     xy_psa2_djf(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(psa2_pc_djf)),psa2_pc_djf,xyres)   ; use standardized timeseries  
     xy_psa2_mam(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(psa2_pc_mam)),psa2_pc_mam,xyres)   ; use standardized timeseries  
     xy_psa2_jja(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(psa2_pc_jja)),psa2_pc_jja,xyres)   ; use standardized timeseries  
     xy_psa2_son(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(psa2_pc_son)),psa2_pc_son,xyres)   ; use standardized timeseries  
     xy_psa2_ann(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(psa2_pc_ann)),psa2_pc_ann,xyres)   ; use standardized timeseries  
     delete([/psa2_pc_djf,psa2_pc_mam,psa2_pc_jja,psa2_pc_son,psa2_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then          
        xy_psa2_mon(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(psa2_pc_mon)),psa2_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/psa2_pc_mon/])
     end if 

     xy_npi(ee) = gsn_csm_xy(wks_psa2_ts,fspan(syear(ee),eyear(ee),dimsizes(npi_ndjfm)),npi_ndjfm,xyres)   ; throw NPI into wks_psa2_ts workstation
     delete(npi_ndjfm)   
     delete(sstreg_plot_flag)  
  end do
  
  if (isvar("clim_syear")) then
     delete(clim_syear)
  end if
  if (isvar("clim_eyear")) then
     delete(clim_eyear)
  end if
  
  if (isfilepresent2("obs_psl")) then    ; for pattern correlation table  
     clat = cos(0.01745329*patcor_nam_djf&lat)
     clat!0 = "lat"
     clat&lat = patcor_nam_djf&lat
     line3   = "                  "    ; Must be 18 characters long
     line4   = line3

     finpr_nam_djf  = "NAM (DJF)         "    ; Must be 18 characters long 
     finpr_sam_djf  = "SAM (DJF)         "     
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
        if (all(ismissing(patcor_nam_djf(hh,{20:},:)))) then
           finpr_nam_djf = finpr_nam_djf+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr_nam_djf = finpr_nam_djf+sprintf(format2,(pattern_cor(patcor_nam_djf(0,{20:},:),patcor_nam_djf(hh,{20:},:),clat({20:}),0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(NewCosWeight(patcor_nam_djf(0,{20:},:))),ndtooned(NewCosWeight(patcor_nam_djf(hh,{20:},:))))))
        end if
        if (all(ismissing(patcor_sam_djf(hh,{:-20},:)))) then
           finpr_sam_djf = finpr_sam_djf+sprintf(format2,9.99)+"/"+sprintf(format3,9.99)
        else
           finpr_sam_djf = finpr_sam_djf+sprintf(format2,(pattern_cor(patcor_sam_djf(0,{:-20},:),patcor_sam_djf(hh,{:-20},:),clat({:-20}),0)))+"/"+sprintf(format3,(dim_rmsd(ndtooned(NewCosWeight(patcor_sam_djf(0,{:-20},:))),ndtooned(NewCosWeight(patcor_sam_djf(hh,{:-20},:))))))
        end if        
     end do
;     
     if (dimsizes(tochar(line4)).ge.8190) then   ; system or fortran compiler limit
        print("Metrics table warning: Not creating metrics table as size of comparison results in a invalid ascii row size.")   
     else
        if (.not.isfilepresent2(getenv("OUTDIR")+"metrics_orig.txt")) then
           header = (/"","Pattern Correlations   Observations vs. Model(s)",""/)
           write_table(getenv("OUTDIR")+"metrics_orig.txt","w",[/header/],"%s")
           write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/line3/],"%s")
           write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/line4/],"%s")
           delete(header)
        end if    
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr_nam_djf/],"%s")
        write_table(getenv("OUTDIR")+"metrics_orig.txt","a",[/finpr_sam_djf/],"%s")
     end if
     delete([/line3,line4,format2,format3,nchar,ntc,clat,patcor_nam_djf,patcor_nam_jja,patcor_nam_ann,patcor_sam_djf,patcor_sam_jja,patcor_sam_ann/])
     delete([/patcor_nao_djf,patcor_nao_jja,patcor_nao_ann,patcor_npo_djf,patcor_npo_jja,patcor_npo_ann,patcor_pna_djf,patcor_pna_jja,patcor_pna_ann/])
     delete([/patcor_psa1_djf,patcor_psa1_jja,patcor_psa1_ann,patcor_psa2_djf,patcor_psa2_jja,patcor_psa2_ann,dimY,ntb/])
  end if 

  panres = True
  panres@gsnMaximize = True
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@gsnPanelBottom = 0.04
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.55
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  if (nsim.le.5) then
     panres@txFontHeightF = 0.024
  else
     panres@txFontHeightF = 0.016
  end if
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol) 
  
  panres@txString = "NAM (DJF)" 
  gsn_panel(wks_nam,map_nam_djf,(/nrow,ncol/),panres)
  delete(map_nam_djf)
  panres@txString = "NAM (MAM)"
  gsn_panel(wks_nam,map_nam_mam,(/nrow,ncol/),panres)
  delete(map_nam_mam)
  panres@txString = "NAM (JJA)"
  gsn_panel(wks_nam,map_nam_jja,(/nrow,ncol/),panres)
  delete(map_nam_jja)
  panres@txString = "NAM (SON)"
  gsn_panel(wks_nam,map_nam_son,(/nrow,ncol/),panres)
  delete(map_nam_son)
  panres@txString = "NAM (Annual)"
  gsn_panel(wks_nam,map_nam_ann,(/nrow,ncol/),panres)
  delete(map_nam_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "NAM (Monthly)"
     gsn_panel(wks_nam,map_nam_mon,(/nrow,ncol/),panres)
     delete(map_nam_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "NAM SST Regressions (DJF)" 
     gsn_panel(wks_nam,reg_nam_djf,(/nrow,ncol/),panres)
     delete(reg_nam_djf)
     panres@txString = "NAM SST Regressions (MAM)"
     gsn_panel(wks_nam,reg_nam_mam,(/nrow,ncol/),panres)
     delete(reg_nam_mam)
     panres@txString = "NAM SST Regressions (JJA)"
     gsn_panel(wks_nam,reg_nam_jja,(/nrow,ncol/),panres)
     delete(reg_nam_jja)
     panres@txString = "NAM SST Regressions (SON)"
     gsn_panel(wks_nam,reg_nam_son,(/nrow,ncol/),panres)
     delete(reg_nam_son)
     panres@txString = "NAM SST Regressions (Annual)"
     gsn_panel(wks_nam,reg_nam_ann,(/nrow,ncol/),panres)
     delete(reg_nam_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NAM SST Regressions (Monthly)"
        gsn_panel(wks_nam,reg_nam_mon,(/nrow,ncol/),panres)
        delete(reg_nam_mon)
     end if
     delete(wks_nam)  
  end if
  
  panres@txString = "NAO (DJF)" 
  gsn_panel(wks_nao,map_nao_djf,(/nrow,ncol/),panres)
  delete(map_nao_djf)
  panres@txString = "NAO (MAM)"
  gsn_panel(wks_nao,map_nao_mam,(/nrow,ncol/),panres)
  delete(map_nao_mam)
  panres@txString = "NAO (JJA)"
  gsn_panel(wks_nao,map_nao_jja,(/nrow,ncol/),panres)
  delete(map_nao_jja)
  panres@txString = "NAO (SON)"
  gsn_panel(wks_nao,map_nao_son,(/nrow,ncol/),panres)
  delete(map_nao_son)
  panres@txString = "NAO (Annual)"
  gsn_panel(wks_nao,map_nao_ann,(/nrow,ncol/),panres)
  delete(map_nao_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "NAO (Monthly)"
     gsn_panel(wks_nao,map_nao_mon,(/nrow,ncol/),panres)
     delete(map_nao_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "NAO SST Regressions (DJF)" 
     gsn_panel(wks_nao,reg_nao_djf,(/nrow,ncol/),panres)
     delete(reg_nao_djf)
     panres@txString = "NAO SST Regressions (MAM)"
     gsn_panel(wks_nao,reg_nao_mam,(/nrow,ncol/),panres)
     delete(reg_nao_mam)
     panres@txString = "NAO SST Regressions (JJA)"
     gsn_panel(wks_nao,reg_nao_jja,(/nrow,ncol/),panres)
     delete(reg_nao_jja)
     panres@txString = "NAO SST Regressions (SON)"
     gsn_panel(wks_nao,reg_nao_son,(/nrow,ncol/),panres)
     delete(reg_nao_son)
     panres@txString = "NAO SST Regressions (Annual)"
     gsn_panel(wks_nao,reg_nao_ann,(/nrow,ncol/),panres)
     delete(reg_nao_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NAO SST Regressions (Monthly)"
        gsn_panel(wks_nao,reg_nao_mon,(/nrow,ncol/),panres)
        delete(reg_nao_mon)
     end if
     delete(wks_nao)
  end if
  
  panres@txString = "PNA (DJF)" 
  gsn_panel(wks_pna,map_pna_djf,(/nrow,ncol/),panres)
  delete(map_pna_djf)
  panres@txString = "PNA (MAM)"
  gsn_panel(wks_pna,map_pna_mam,(/nrow,ncol/),panres)
  delete(map_pna_mam)
  panres@txString = "PNA (JJA)"
  gsn_panel(wks_pna,map_pna_jja,(/nrow,ncol/),panres)
  delete(map_pna_jja)
  panres@txString = "PNA (SON)"
  gsn_panel(wks_pna,map_pna_son,(/nrow,ncol/),panres)
  delete(map_pna_son)
  panres@txString = "PNA (Annual)"
  gsn_panel(wks_pna,map_pna_ann,(/nrow,ncol/),panres)
  delete(map_pna_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "PNA (Monthly)"
     gsn_panel(wks_pna,map_pna_mon,(/nrow,ncol/),panres)
     delete(map_pna_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "PNA SST Regressions (DJF)" 
     gsn_panel(wks_pna,reg_pna_djf,(/nrow,ncol/),panres)
     delete(reg_pna_djf)
     panres@txString = "PNA SST Regressions (MAM)"
     gsn_panel(wks_pna,reg_pna_mam,(/nrow,ncol/),panres)
     delete(reg_pna_mam)
     panres@txString = "PNA SST Regressions (JJA)"
     gsn_panel(wks_pna,reg_pna_jja,(/nrow,ncol/),panres)
     delete(reg_pna_jja)
     panres@txString = "PNA SST Regressions (SON)"
     gsn_panel(wks_pna,reg_pna_son,(/nrow,ncol/),panres)
     delete(reg_pna_son)
     panres@txString = "PNA SST Regressions (Annual)"
     gsn_panel(wks_pna,reg_pna_ann,(/nrow,ncol/),panres)
     delete(reg_pna_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "PNA SST Regressions (Monthly)"
        gsn_panel(wks_pna,reg_pna_mon,(/nrow,ncol/),panres)
        delete(reg_pna_mon)
     end if
     delete(wks_pna)
  end if
  
  panres@txString = "NPO (DJF)" 
  gsn_panel(wks_npo,map_npo_djf,(/nrow,ncol/),panres)
  delete(map_npo_djf)
  panres@txString = "NPO (MAM)"
  gsn_panel(wks_npo,map_npo_mam,(/nrow,ncol/),panres)
  delete(map_npo_mam)
  panres@txString = "NPO (JJA)"
  gsn_panel(wks_npo,map_npo_jja,(/nrow,ncol/),panres)
  delete(map_npo_jja)
  panres@txString = "NPO (SON)"
  gsn_panel(wks_npo,map_npo_son,(/nrow,ncol/),panres)
  delete(map_npo_son)
  panres@txString = "NPO (Annual)"
  gsn_panel(wks_npo,map_npo_ann,(/nrow,ncol/),panres)
  delete(map_npo_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "NPO (Monthly)"
     gsn_panel(wks_npo,map_npo_mon,(/nrow,ncol/),panres)
     delete(map_npo_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "NPO SST Regressions (DJF)" 
     gsn_panel(wks_npo,reg_npo_djf,(/nrow,ncol/),panres)
     delete(reg_npo_djf)
     panres@txString = "NPO SST Regressions (MAM)"
     gsn_panel(wks_npo,reg_npo_mam,(/nrow,ncol/),panres)
     delete(reg_npo_mam)
     panres@txString = "NPO SST Regressions (JJA)"
     gsn_panel(wks_npo,reg_npo_jja,(/nrow,ncol/),panres)
     delete(reg_npo_jja)
     panres@txString = "NPO SST Regressions (SON)"
     gsn_panel(wks_npo,reg_npo_son,(/nrow,ncol/),panres)
     delete(reg_npo_son)
     panres@txString = "NPO SST Regressions (Annual)"
     gsn_panel(wks_npo,reg_npo_ann,(/nrow,ncol/),panres)
     delete(reg_npo_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NPO SST Regressions (Monthly)"
        gsn_panel(wks_npo,reg_npo_mon,(/nrow,ncol/),panres)
        delete(reg_npo_mon)
     end if
     delete(wks_npo)
  end if
  
  panres@txString = "SAM (DJF)" 
  gsn_panel(wks_sam,map_sam_djf,(/nrow,ncol/),panres)
  delete(map_sam_djf)
  panres@txString = "SAM (MAM)"
  gsn_panel(wks_sam,map_sam_mam,(/nrow,ncol/),panres)
  delete(map_sam_mam)
  panres@txString = "SAM (JJA)"
  gsn_panel(wks_sam,map_sam_jja,(/nrow,ncol/),panres)
  delete(map_sam_jja)
  panres@txString = "SAM (SON)"
  gsn_panel(wks_sam,map_sam_son,(/nrow,ncol/),panres)
  delete(map_sam_son)
  panres@txString = "SAM (Annual)"
  gsn_panel(wks_sam,map_sam_ann,(/nrow,ncol/),panres)
  delete(map_sam_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "SAM (Monthly)"
     gsn_panel(wks_sam,map_sam_mon,(/nrow,ncol/),panres)
     delete(map_sam_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "SAM SST Regressions (DJF)" 
     gsn_panel(wks_sam,reg_sam_djf,(/nrow,ncol/),panres)
     delete(reg_sam_djf)
     panres@txString = "SAM SST Regressions (MAM)"
     gsn_panel(wks_sam,reg_sam_mam,(/nrow,ncol/),panres)
     delete(reg_sam_mam)
     panres@txString = "SAM SST Regressions (JJA)"
     gsn_panel(wks_sam,reg_sam_jja,(/nrow,ncol/),panres)
     delete(reg_sam_jja)
     panres@txString = "SAM SST Regressions (SON)"
     gsn_panel(wks_sam,reg_sam_son,(/nrow,ncol/),panres)
     delete(reg_sam_son)
     panres@txString = "SAM SST Regressions (Annual)"
     gsn_panel(wks_sam,reg_sam_ann,(/nrow,ncol/),panres)
     delete(reg_sam_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "SAM SST Regressions (Monthly)"
        gsn_panel(wks_sam,reg_sam_mon,(/nrow,ncol/),panres)
        delete(reg_sam_mon)
     end if
     delete(wks_sam)
  end if
  
  panres@txString = "PSA1 (DJF)" 
  gsn_panel(wks_psa1,map_psa1_djf,(/nrow,ncol/),panres)
  delete(map_psa1_djf)
  panres@txString = "PSA1 (MAM)"
  gsn_panel(wks_psa1,map_psa1_mam,(/nrow,ncol/),panres)
  delete(map_psa1_mam)
  panres@txString = "PSA1 (JJA)"
  gsn_panel(wks_psa1,map_psa1_jja,(/nrow,ncol/),panres)
  delete(map_psa1_jja)
  panres@txString = "PSA1 (SON)"
  gsn_panel(wks_psa1,map_psa1_son,(/nrow,ncol/),panres)
  delete(map_psa1_son)
  panres@txString = "PSA1 (Annual)"
  gsn_panel(wks_psa1,map_psa1_ann,(/nrow,ncol/),panres)
  delete(map_psa1_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "PSA1 (Monthly)"
     gsn_panel(wks_psa1,map_psa1_mon,(/nrow,ncol/),panres)
     delete(map_psa1_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "PSA1 SST Regressions (DJF)" 
     gsn_panel(wks_psa1,reg_psa1_djf,(/nrow,ncol/),panres)
     delete(reg_psa1_djf)
     panres@txString = "PSA1 SST Regressions (MAM)"
     gsn_panel(wks_psa1,reg_psa1_mam,(/nrow,ncol/),panres)
     delete(reg_psa1_mam)
     panres@txString = "PSA1 SST Regressions (JJA)"
     gsn_panel(wks_psa1,reg_psa1_jja,(/nrow,ncol/),panres)
     delete(reg_psa1_jja)
     panres@txString = "PSA1 SST Regressions (SON)"
     gsn_panel(wks_psa1,reg_psa1_son,(/nrow,ncol/),panres)
     delete(reg_psa1_son)
     panres@txString = "PSA1 SST Regressions (Annual)"
     gsn_panel(wks_psa1,reg_psa1_ann,(/nrow,ncol/),panres)
     delete(reg_psa1_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "PSA1 SST Regressions (Monthly)"
        gsn_panel(wks_psa1,reg_psa1_mon,(/nrow,ncol/),panres)
        delete(reg_psa1_mon)
     end if
     delete(wks_psa1)
  end if
  
  panres@txString = "PSA2 (DJF)" 
  gsn_panel(wks_psa2,map_psa2_djf,(/nrow,ncol/),panres)
  delete(map_psa2_djf)
  panres@txString = "PSA2 (MAM)"
  gsn_panel(wks_psa2,map_psa2_mam,(/nrow,ncol/),panres)
  delete(map_psa2_mam)
  panres@txString = "PSA2 (JJA)"
  gsn_panel(wks_psa2,map_psa2_jja,(/nrow,ncol/),panres)
  delete(map_psa2_jja)
  panres@txString = "PSA2 (SON)"
  gsn_panel(wks_psa2,map_psa2_son,(/nrow,ncol/),panres)
  delete(map_psa2_son)
  panres@txString = "PSA2 (Annual)"
  gsn_panel(wks_psa2,map_psa2_ann,(/nrow,ncol/),panres)
  delete(map_psa2_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "PSA2 (Monthly)"
     gsn_panel(wks_psa2,map_psa2_mon,(/nrow,ncol/),panres)
     delete(map_psa2_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     panres@txString = "PSA2 SST Regressions (DJF)" 
     gsn_panel(wks_psa2,reg_psa2_djf,(/nrow,ncol/),panres)
     delete(reg_psa2_djf)
     panres@txString = "PSA2 SST Regressions (MAM)"
     gsn_panel(wks_psa2,reg_psa2_mam,(/nrow,ncol/),panres)
     delete(reg_psa2_mam)
     panres@txString = "PSA2 SST Regressions (JJA)"
     gsn_panel(wks_psa2,reg_psa2_jja,(/nrow,ncol/),panres)
     delete(reg_psa2_jja)
     panres@txString = "PSA2 SST Regressions (SON)"
     gsn_panel(wks_psa2,reg_psa2_son,(/nrow,ncol/),panres)
     delete(reg_psa2_son)
     panres@txString = "PSA2 SST Regressions (Annual)"
     gsn_panel(wks_psa2,reg_psa2_ann,(/nrow,ncol/),panres)
     delete(reg_psa2_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "PSA2 SST Regressions (Monthly)"
        gsn_panel(wks_psa2,reg_psa2_mon,(/nrow,ncol/),panres)
        delete(reg_psa2_mon)
     end if
     delete(wks_psa2)
  end if

  panres2 = True    
  if (nsim.le.5) then
     panres2@txFontHeightF = 0.024
  else
     panres2@txFontHeightF = 0.016
  end if  
  panres2@gsnMaximize = True
  panres2@gsnPaperOrientation = "portrait"
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres2@gsnPanelScalePlotIndex = tt(0)
     delete(tt)  
  end if  
  if (nsim.le.12) then
     lp = (/nsim,1/)
  else
     lp = (/nrow,ncol/)  
  end if
  
  panres2@txString = "NAM (DJF)"
  gsn_panel(wks_nam_ts,xy_nam_djf,lp,panres2)  
  delete(xy_nam_djf)
  panres2@txString = "NAM (MAM)"
  gsn_panel(wks_nam_ts,xy_nam_mam,lp,panres2)  
  delete(xy_nam_mam)
  panres2@txString = "NAM (JJA)"
  gsn_panel(wks_nam_ts,xy_nam_jja,lp,panres2)  
  delete(xy_nam_jja)
  panres2@txString = "NAM (SON)"
  gsn_panel(wks_nam_ts,xy_nam_son,lp,panres2)  
  delete(xy_nam_son)
  panres2@txString = "NAM (Annual)"
  gsn_panel(wks_nam_ts,xy_nam_ann,lp,panres2)  
  delete(xy_nam_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "NAM (Monthly)"
     gsn_panel(wks_nam_ts,xy_nam_mon,lp,panres2)  
     delete(xy_nam_mon) 
  end if
  delete(wks_nam_ts)
    
  panres2@txString = "NAO (DJF)"
  gsn_panel(wks_nao_ts,xy_nao_djf,lp,panres2)  
  delete(xy_nao_djf)
  panres2@txString = "NAO (MAM)"
  gsn_panel(wks_nao_ts,xy_nao_mam,lp,panres2)  
  delete(xy_nao_mam)
  panres2@txString = "NAO (JJA)"
  gsn_panel(wks_nao_ts,xy_nao_jja,lp,panres2)  
  delete(xy_nao_jja)
  panres2@txString = "NAO (SON)"
  gsn_panel(wks_nao_ts,xy_nao_son,lp,panres2)  
  delete(xy_nao_son)
  panres2@txString = "NAO (Annual)"
  gsn_panel(wks_nao_ts,xy_nao_ann,lp,panres2)  
  delete(xy_nao_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "NAO (Monthly)"
     gsn_panel(wks_nao_ts,xy_nao_mon,lp,panres2)  
     delete(xy_nao_mon) 
  end if
  delete(wks_nao_ts)
  
  panres2@txString = "PNA (DJF)"
  gsn_panel(wks_pna_ts,xy_pna_djf,lp,panres2)  
  delete(xy_pna_djf)
  panres2@txString = "PNA (MAM)"
  gsn_panel(wks_pna_ts,xy_pna_mam,lp,panres2)  
  delete(xy_pna_mam)
  panres2@txString = "PNA (JJA)"
  gsn_panel(wks_pna_ts,xy_pna_jja,lp,panres2)  
  delete(xy_pna_jja)
  panres2@txString = "PNA (SON)"
  gsn_panel(wks_pna_ts,xy_pna_son,lp,panres2)  
  delete(xy_pna_son)
  panres2@txString = "PNA (Annual)"
  gsn_panel(wks_pna_ts,xy_pna_ann,lp,panres2)  
  delete(xy_pna_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "PNA (Monthly)"
     gsn_panel(wks_pna_ts,xy_pna_mon,lp,panres2)  
     delete(xy_pna_mon) 
  end if
  delete(wks_pna_ts)
  
  panres2@txString = "NPO (DJF)"
  gsn_panel(wks_npo_ts,xy_npo_djf,lp,panres2)  
  delete(xy_npo_djf)
  panres2@txString = "NPO (MAM)"
  gsn_panel(wks_npo_ts,xy_npo_mam,lp,panres2)  
  delete(xy_npo_mam)
  panres2@txString = "NPO (JJA)"
  gsn_panel(wks_npo_ts,xy_npo_jja,lp,panres2)  
  delete(xy_npo_jja)
  panres2@txString = "NPO (SON)"
  gsn_panel(wks_npo_ts,xy_npo_son,lp,panres2)  
  delete(xy_npo_son)
  panres2@txString = "NPO (Annual)"
  gsn_panel(wks_npo_ts,xy_npo_ann,lp,panres2)  
  delete(xy_npo_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "NPO (Monthly)"
     gsn_panel(wks_npo_ts,xy_npo_mon,lp,panres2)  
     delete(xy_npo_mon) 
  end if
  delete(wks_npo_ts)

  panres2@txString = "SAM (DJF)"
  gsn_panel(wks_sam_ts,xy_sam_djf,lp,panres2)  
  delete(xy_sam_djf)
  panres2@txString = "SAM (MAM)"
  gsn_panel(wks_sam_ts,xy_sam_mam,lp,panres2)  
  delete(xy_sam_mam)
  panres2@txString = "SAM (JJA)"
  gsn_panel(wks_sam_ts,xy_sam_jja,lp,panres2)  
  delete(xy_sam_jja)
  panres2@txString = "SAM (SON)"
  gsn_panel(wks_sam_ts,xy_sam_son,lp,panres2)  
  delete(xy_sam_son)
  panres2@txString = "SAM (Annual)"
  gsn_panel(wks_sam_ts,xy_sam_ann,lp,panres2)  
  delete(xy_sam_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "SAM (Monthly)"
     gsn_panel(wks_sam_ts,xy_sam_mon,lp,panres2)  
     delete(xy_sam_mon) 
  end if
  delete(wks_sam_ts)
    
  panres2@txString = "PSA1 (DJF)"
  gsn_panel(wks_psa1_ts,xy_psa1_djf,lp,panres2)  
  delete(xy_psa1_djf)
  panres2@txString = "PSA1 (MAM)"
  gsn_panel(wks_psa1_ts,xy_psa1_mam,lp,panres2)  
  delete(xy_psa1_mam)
  panres2@txString = "PSA1 (JJA)"
  gsn_panel(wks_psa1_ts,xy_psa1_jja,lp,panres2)  
  delete(xy_psa1_jja)
  panres2@txString = "PSA1 (SON)"
  gsn_panel(wks_psa1_ts,xy_psa1_son,lp,panres2)  
  delete(xy_psa1_son)
  panres2@txString = "PSA1 (Annual)"
  gsn_panel(wks_psa1_ts,xy_psa1_ann,lp,panres2)  
  delete(xy_psa1_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "PSA1 (Monthly)"
     gsn_panel(wks_psa1_ts,xy_psa1_mon,lp,panres2)  
     delete(xy_psa1_mon) 
  end if
  delete(wks_psa1_ts)
  
  panres2@txString = "PSA2 (DJF)"
  gsn_panel(wks_psa2_ts,xy_psa2_djf,lp,panres2)  
  delete(xy_psa2_djf)
  panres2@txString = "PSA2 (MAM)"
  gsn_panel(wks_psa2_ts,xy_psa2_mam,lp,panres2)  
  delete(xy_psa2_mam)
  panres2@txString = "PSA2 (JJA)"
  gsn_panel(wks_psa2_ts,xy_psa2_jja,lp,panres2)  
  delete(xy_psa2_jja)
  panres2@txString = "PSA2 (SON)"
  gsn_panel(wks_psa2_ts,xy_psa2_son,lp,panres2)  
  delete(xy_psa2_son)
  panres2@txString = "PSA2 (Annual)"
  gsn_panel(wks_psa2_ts,xy_psa2_ann,lp,panres2)  
  delete(xy_psa2_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "PSA2 (Monthly)"
     gsn_panel(wks_psa2_ts,xy_psa2_mon,lp,panres2)  
     delete(xy_psa2_mon) 
  end if
  
  panres2@txString = "NPI (NDJFM)"
  gsn_panel(wks_psa2_ts,xy_npi,lp,panres2)  
  delete(xy_npi) 
  delete(wks_psa2_ts)  
;--------------------------------------------------------------------------------------------------
  OUTDIR = getenv("OUTDIR") 
  if (wks_type.eq."png") then  
     system("mv "+OUTDIR+"nam.000001.png "+OUTDIR+"nam.djf.png") 
     system("mv "+OUTDIR+"nam.000002.png "+OUTDIR+"nam.mam.png") 
     system("mv "+OUTDIR+"nam.000003.png "+OUTDIR+"nam.jja.png") 
     system("mv "+OUTDIR+"nam.000004.png "+OUTDIR+"nam.son.png") 
     system("mv "+OUTDIR+"nam.000005.png "+OUTDIR+"nam.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"nam.000006.png "+OUTDIR+"nam.mon.png") 
        if (sstreg_frame.eq.0) then 
           system("mv "+OUTDIR+"nam.000007.png "+OUTDIR+"nam.sstreg.djf.png") 
           system("mv "+OUTDIR+"nam.000008.png "+OUTDIR+"nam.sstreg.mam.png") 
           system("mv "+OUTDIR+"nam.000009.png "+OUTDIR+"nam.sstreg.jja.png") 
           system("mv "+OUTDIR+"nam.000010.png "+OUTDIR+"nam.sstreg.son.png") 
           system("mv "+OUTDIR+"nam.000011.png "+OUTDIR+"nam.sstreg.ann.png")
           system("mv "+OUTDIR+"nam.000012.png "+OUTDIR+"nam.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"nam.000006.png "+OUTDIR+"nam.sstreg.djf.png") 
           system("mv "+OUTDIR+"nam.000007.png "+OUTDIR+"nam.sstreg.mam.png") 
           system("mv "+OUTDIR+"nam.000008.png "+OUTDIR+"nam.sstreg.jja.png") 
           system("mv "+OUTDIR+"nam.000009.png "+OUTDIR+"nam.sstreg.son.png") 
           system("mv "+OUTDIR+"nam.000010.png "+OUTDIR+"nam.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"nao.000001.png "+OUTDIR+"nao.djf.png") 
     system("mv "+OUTDIR+"nao.000002.png "+OUTDIR+"nao.mam.png") 
     system("mv "+OUTDIR+"nao.000003.png "+OUTDIR+"nao.jja.png") 
     system("mv "+OUTDIR+"nao.000004.png "+OUTDIR+"nao.son.png") 
     system("mv "+OUTDIR+"nao.000005.png "+OUTDIR+"nao.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"nao.000006.png "+OUTDIR+"nao.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"nao.000007.png "+OUTDIR+"nao.sstreg.djf.png") 
           system("mv "+OUTDIR+"nao.000008.png "+OUTDIR+"nao.sstreg.mam.png") 
           system("mv "+OUTDIR+"nao.000009.png "+OUTDIR+"nao.sstreg.jja.png") 
           system("mv "+OUTDIR+"nao.000010.png "+OUTDIR+"nao.sstreg.son.png") 
           system("mv "+OUTDIR+"nao.000011.png "+OUTDIR+"nao.sstreg.ann.png")
           system("mv "+OUTDIR+"nao.000012.png "+OUTDIR+"nao.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"nao.000006.png "+OUTDIR+"nao.sstreg.djf.png") 
           system("mv "+OUTDIR+"nao.000007.png "+OUTDIR+"nao.sstreg.mam.png") 
           system("mv "+OUTDIR+"nao.000008.png "+OUTDIR+"nao.sstreg.jja.png") 
           system("mv "+OUTDIR+"nao.000009.png "+OUTDIR+"nao.sstreg.son.png") 
           system("mv "+OUTDIR+"nao.000010.png "+OUTDIR+"nao.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"pna.000001.png "+OUTDIR+"pna.djf.png") 
     system("mv "+OUTDIR+"pna.000002.png "+OUTDIR+"pna.mam.png") 
     system("mv "+OUTDIR+"pna.000003.png "+OUTDIR+"pna.jja.png") 
     system("mv "+OUTDIR+"pna.000004.png "+OUTDIR+"pna.son.png") 
     system("mv "+OUTDIR+"pna.000005.png "+OUTDIR+"pna.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pna.000006.png "+OUTDIR+"pna.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pna.000007.png "+OUTDIR+"pna.sstreg.djf.png") 
           system("mv "+OUTDIR+"pna.000008.png "+OUTDIR+"pna.sstreg.mam.png") 
           system("mv "+OUTDIR+"pna.000009.png "+OUTDIR+"pna.sstreg.jja.png") 
           system("mv "+OUTDIR+"pna.000010.png "+OUTDIR+"pna.sstreg.son.png") 
           system("mv "+OUTDIR+"pna.000011.png "+OUTDIR+"pna.sstreg.ann.png")
           system("mv "+OUTDIR+"pna.000012.png "+OUTDIR+"pna.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pna.000006.png "+OUTDIR+"pna.sstreg.djf.png") 
           system("mv "+OUTDIR+"pna.000007.png "+OUTDIR+"pna.sstreg.mam.png") 
           system("mv "+OUTDIR+"pna.000008.png "+OUTDIR+"pna.sstreg.jja.png") 
           system("mv "+OUTDIR+"pna.000009.png "+OUTDIR+"pna.sstreg.son.png") 
           system("mv "+OUTDIR+"pna.000010.png "+OUTDIR+"pna.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"npo.000001.png "+OUTDIR+"npo.djf.png") 
     system("mv "+OUTDIR+"npo.000002.png "+OUTDIR+"npo.mam.png") 
     system("mv "+OUTDIR+"npo.000003.png "+OUTDIR+"npo.jja.png") 
     system("mv "+OUTDIR+"npo.000004.png "+OUTDIR+"npo.son.png") 
     system("mv "+OUTDIR+"npo.000005.png "+OUTDIR+"npo.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"npo.000006.png "+OUTDIR+"npo.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"npo.000007.png "+OUTDIR+"npo.sstreg.djf.png") 
           system("mv "+OUTDIR+"npo.000008.png "+OUTDIR+"npo.sstreg.mam.png") 
           system("mv "+OUTDIR+"npo.000009.png "+OUTDIR+"npo.sstreg.jja.png") 
           system("mv "+OUTDIR+"npo.000010.png "+OUTDIR+"npo.sstreg.son.png") 
           system("mv "+OUTDIR+"npo.000011.png "+OUTDIR+"npo.sstreg.ann.png")
           system("mv "+OUTDIR+"npo.000012.png "+OUTDIR+"npo.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"npo.000006.png "+OUTDIR+"npo.sstreg.djf.png") 
           system("mv "+OUTDIR+"npo.000007.png "+OUTDIR+"npo.sstreg.mam.png") 
           system("mv "+OUTDIR+"npo.000008.png "+OUTDIR+"npo.sstreg.jja.png") 
           system("mv "+OUTDIR+"npo.000009.png "+OUTDIR+"npo.sstreg.son.png") 
           system("mv "+OUTDIR+"npo.000010.png "+OUTDIR+"npo.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"sam.000001.png "+OUTDIR+"sam.djf.png") 
     system("mv "+OUTDIR+"sam.000002.png "+OUTDIR+"sam.mam.png") 
     system("mv "+OUTDIR+"sam.000003.png "+OUTDIR+"sam.jja.png") 
     system("mv "+OUTDIR+"sam.000004.png "+OUTDIR+"sam.son.png") 
     system("mv "+OUTDIR+"sam.000005.png "+OUTDIR+"sam.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"sam.000006.png "+OUTDIR+"sam.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"sam.000007.png "+OUTDIR+"sam.sstreg.djf.png") 
           system("mv "+OUTDIR+"sam.000008.png "+OUTDIR+"sam.sstreg.mam.png") 
           system("mv "+OUTDIR+"sam.000009.png "+OUTDIR+"sam.sstreg.jja.png") 
           system("mv "+OUTDIR+"sam.000010.png "+OUTDIR+"sam.sstreg.son.png") 
           system("mv "+OUTDIR+"sam.000011.png "+OUTDIR+"sam.sstreg.ann.png")
           system("mv "+OUTDIR+"sam.000012.png "+OUTDIR+"sam.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"sam.000006.png "+OUTDIR+"sam.sstreg.djf.png") 
           system("mv "+OUTDIR+"sam.000007.png "+OUTDIR+"sam.sstreg.mam.png") 
           system("mv "+OUTDIR+"sam.000008.png "+OUTDIR+"sam.sstreg.jja.png") 
           system("mv "+OUTDIR+"sam.000009.png "+OUTDIR+"sam.sstreg.son.png") 
           system("mv "+OUTDIR+"sam.000010.png "+OUTDIR+"sam.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"psa1.000001.png "+OUTDIR+"psa1.djf.png") 
     system("mv "+OUTDIR+"psa1.000002.png "+OUTDIR+"psa1.mam.png") 
     system("mv "+OUTDIR+"psa1.000003.png "+OUTDIR+"psa1.jja.png") 
     system("mv "+OUTDIR+"psa1.000004.png "+OUTDIR+"psa1.son.png") 
     system("mv "+OUTDIR+"psa1.000005.png "+OUTDIR+"psa1.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psa1.000006.png "+OUTDIR+"psa1.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psa1.000007.png "+OUTDIR+"psa1.sstreg.djf.png") 
           system("mv "+OUTDIR+"psa1.000008.png "+OUTDIR+"psa1.sstreg.mam.png") 
           system("mv "+OUTDIR+"psa1.000009.png "+OUTDIR+"psa1.sstreg.jja.png") 
           system("mv "+OUTDIR+"psa1.000010.png "+OUTDIR+"psa1.sstreg.son.png") 
           system("mv "+OUTDIR+"psa1.000011.png "+OUTDIR+"psa1.sstreg.ann.png")
           system("mv "+OUTDIR+"psa1.000012.png "+OUTDIR+"psa1.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psa1.000006.png "+OUTDIR+"psa1.sstreg.djf.png") 
           system("mv "+OUTDIR+"psa1.000007.png "+OUTDIR+"psa1.sstreg.mam.png") 
           system("mv "+OUTDIR+"psa1.000008.png "+OUTDIR+"psa1.sstreg.jja.png") 
           system("mv "+OUTDIR+"psa1.000009.png "+OUTDIR+"psa1.sstreg.son.png") 
           system("mv "+OUTDIR+"psa1.000010.png "+OUTDIR+"psa1.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"psa2.000001.png "+OUTDIR+"psa2.djf.png") 
     system("mv "+OUTDIR+"psa2.000002.png "+OUTDIR+"psa2.mam.png") 
     system("mv "+OUTDIR+"psa2.000003.png "+OUTDIR+"psa2.jja.png") 
     system("mv "+OUTDIR+"psa2.000004.png "+OUTDIR+"psa2.son.png") 
     system("mv "+OUTDIR+"psa2.000005.png "+OUTDIR+"psa2.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psa2.000006.png "+OUTDIR+"psa2.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psa2.000007.png "+OUTDIR+"psa2.sstreg.djf.png") 
           system("mv "+OUTDIR+"psa2.000008.png "+OUTDIR+"psa2.sstreg.mam.png") 
           system("mv "+OUTDIR+"psa2.000009.png "+OUTDIR+"psa2.sstreg.jja.png") 
           system("mv "+OUTDIR+"psa2.000010.png "+OUTDIR+"psa2.sstreg.son.png") 
           system("mv "+OUTDIR+"psa2.000011.png "+OUTDIR+"psa2.sstreg.ann.png")
           system("mv "+OUTDIR+"psa2.000012.png "+OUTDIR+"psa2.sstreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psa2.000006.png "+OUTDIR+"psa2.sstreg.djf.png") 
           system("mv "+OUTDIR+"psa2.000007.png "+OUTDIR+"psa2.sstreg.mam.png") 
           system("mv "+OUTDIR+"psa2.000008.png "+OUTDIR+"psa2.sstreg.jja.png") 
           system("mv "+OUTDIR+"psa2.000009.png "+OUTDIR+"psa2.sstreg.son.png") 
           system("mv "+OUTDIR+"psa2.000010.png "+OUTDIR+"psa2.sstreg.ann.png")
        end if
     end if
     
     system("mv "+OUTDIR+"nam.timeseries.000001.png "+OUTDIR+"nam.timeseries.djf.png") 
     system("mv "+OUTDIR+"nam.timeseries.000002.png "+OUTDIR+"nam.timeseries.mam.png") 
     system("mv "+OUTDIR+"nam.timeseries.000003.png "+OUTDIR+"nam.timeseries.jja.png") 
     system("mv "+OUTDIR+"nam.timeseries.000004.png "+OUTDIR+"nam.timeseries.son.png") 
     system("mv "+OUTDIR+"nam.timeseries.000005.png "+OUTDIR+"nam.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"nam.timeseries.000006.png "+OUTDIR+"nam.timeseries.mon.png")     
     end if
     
     system("mv "+OUTDIR+"nao.timeseries.000001.png "+OUTDIR+"nao.timeseries.djf.png") 
     system("mv "+OUTDIR+"nao.timeseries.000002.png "+OUTDIR+"nao.timeseries.mam.png") 
     system("mv "+OUTDIR+"nao.timeseries.000003.png "+OUTDIR+"nao.timeseries.jja.png") 
     system("mv "+OUTDIR+"nao.timeseries.000004.png "+OUTDIR+"nao.timeseries.son.png") 
     system("mv "+OUTDIR+"nao.timeseries.000005.png "+OUTDIR+"nao.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"nao.timeseries.000006.png "+OUTDIR+"nao.timeseries.mon.png")     
     end if
     
     system("mv "+OUTDIR+"pna.timeseries.000001.png "+OUTDIR+"pna.timeseries.djf.png") 
     system("mv "+OUTDIR+"pna.timeseries.000002.png "+OUTDIR+"pna.timeseries.mam.png") 
     system("mv "+OUTDIR+"pna.timeseries.000003.png "+OUTDIR+"pna.timeseries.jja.png") 
     system("mv "+OUTDIR+"pna.timeseries.000004.png "+OUTDIR+"pna.timeseries.son.png") 
     system("mv "+OUTDIR+"pna.timeseries.000005.png "+OUTDIR+"pna.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pna.timeseries.000006.png "+OUTDIR+"pna.timeseries.mon.png")     
     end if
     
     system("mv "+OUTDIR+"npo.timeseries.000001.png "+OUTDIR+"npo.timeseries.djf.png") 
     system("mv "+OUTDIR+"npo.timeseries.000002.png "+OUTDIR+"npo.timeseries.mam.png") 
     system("mv "+OUTDIR+"npo.timeseries.000003.png "+OUTDIR+"npo.timeseries.jja.png") 
     system("mv "+OUTDIR+"npo.timeseries.000004.png "+OUTDIR+"npo.timeseries.son.png") 
     system("mv "+OUTDIR+"npo.timeseries.000005.png "+OUTDIR+"npo.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"npo.timeseries.000006.png "+OUTDIR+"npo.timeseries.mon.png")     
     end if
     
     system("mv "+OUTDIR+"sam.timeseries.000001.png "+OUTDIR+"sam.timeseries.djf.png") 
     system("mv "+OUTDIR+"sam.timeseries.000002.png "+OUTDIR+"sam.timeseries.mam.png") 
     system("mv "+OUTDIR+"sam.timeseries.000003.png "+OUTDIR+"sam.timeseries.jja.png") 
     system("mv "+OUTDIR+"sam.timeseries.000004.png "+OUTDIR+"sam.timeseries.son.png") 
     system("mv "+OUTDIR+"sam.timeseries.000005.png "+OUTDIR+"sam.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"sam.timeseries.000006.png "+OUTDIR+"sam.timeseries.mon.png")     
     end if
     
     system("mv "+OUTDIR+"psa1.timeseries.000001.png "+OUTDIR+"psa1.timeseries.djf.png") 
     system("mv "+OUTDIR+"psa1.timeseries.000002.png "+OUTDIR+"psa1.timeseries.mam.png") 
     system("mv "+OUTDIR+"psa1.timeseries.000003.png "+OUTDIR+"psa1.timeseries.jja.png") 
     system("mv "+OUTDIR+"psa1.timeseries.000004.png "+OUTDIR+"psa1.timeseries.son.png") 
     system("mv "+OUTDIR+"psa1.timeseries.000005.png "+OUTDIR+"psa1.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psa1.timeseries.000006.png "+OUTDIR+"psa1.timeseries.mon.png")     
     end if
     
     system("mv "+OUTDIR+"psa2.timeseries.000001.png "+OUTDIR+"psa2.timeseries.djf.png") 
     system("mv "+OUTDIR+"psa2.timeseries.000002.png "+OUTDIR+"psa2.timeseries.mam.png") 
     system("mv "+OUTDIR+"psa2.timeseries.000003.png "+OUTDIR+"psa2.timeseries.jja.png") 
     system("mv "+OUTDIR+"psa2.timeseries.000004.png "+OUTDIR+"psa2.timeseries.son.png") 
     system("mv "+OUTDIR+"psa2.timeseries.000005.png "+OUTDIR+"psa2.timeseries.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psa2.timeseries.000006.png "+OUTDIR+"psa2.timeseries.mon.png")     
        system("mv "+OUTDIR+"psa2.timeseries.000007.png "+OUTDIR+"npi.timeseries.ndjfm.png")     
     else
        system("mv "+OUTDIR+"psa2.timeseries.000006.png "+OUTDIR+"npi.timeseries.ndjfm.png")     
     end if
  else
     system("psplit "+OUTDIR+"nam.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nam.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nam.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nam.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nam.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"nam.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"nam.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"nam.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"nam.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"nam.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"nam.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"nam.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"nam.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"nam.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"nam.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"nam.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"nam.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"nam.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"nao.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nao.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nao.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nao.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nao.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"nao.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"nao.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"nao.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"nao.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"nao.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"nao.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"nao.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"nao.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"nao.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"nao.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"nao.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"nao.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"nao.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"pna.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"pna.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"pna.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"pna.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"pna.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"pna.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"pna.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"pna.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"pna.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"pna.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"pna.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"pna.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"pna.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"pna.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"pna.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"pna.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"pna.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"pna.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"npo.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"npo.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"npo.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"npo.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"npo.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"npo.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"npo.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"npo.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"npo.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"npo.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"npo.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"npo.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"npo.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"npo.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"npo.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"npo.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"npo.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"npo.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"sam.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"sam.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"sam.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"sam.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"sam.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"sam.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"sam.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"sam.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"sam.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"sam.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"sam.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"sam.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"sam.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"sam.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"sam.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"sam.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"sam.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"sam.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"psa1.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"psa1.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"psa1.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"psa1.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"psa1.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"psa1.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"psa1.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"psa1.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"psa1.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"psa1.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"psa1.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"psa1.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"psa1.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"psa1.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"psa1.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"psa1.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"psa1.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"psa1.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"psa2.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"psa2.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"psa2.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"psa2.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"psa2.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"psa2.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"psa2.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"psa2.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"psa2.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"psa2.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"psa2.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0011.ps "+OUTDIR+"psa2.sstreg.ann.ps")
           system("mv "+OUTDIR+"pict0012.ps "+OUTDIR+"psa2.sstreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"psa2.sstreg.djf.ps") 
           system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"psa2.sstreg.mam.ps") 
           system("mv "+OUTDIR+"pict0008.ps "+OUTDIR+"psa2.sstreg.jja.ps") 
           system("mv "+OUTDIR+"pict0009.ps "+OUTDIR+"psa2.sstreg.son.ps") 
           system("mv "+OUTDIR+"pict0010.ps "+OUTDIR+"psa2.sstreg.ann.ps")     
        end if
     end if
     
     system("psplit "+OUTDIR+"nam.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nam.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nam.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nam.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nam.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"nam.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"nam.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"nao.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"nao.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"nao.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"nao.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"nao.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"nao.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"nao.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"pna.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"pna.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"pna.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"pna.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"pna.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"pna.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"pna.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"npo.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"npo.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"npo.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"npo.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"npo.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"npo.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"npo.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"sam.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"sam.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"sam.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"sam.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"sam.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"sam.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"sam.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"psa1.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"psa1.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"psa1.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"psa1.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"psa1.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"psa1.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"psa1.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"psa2.timeseries.ps "+OUTDIR+"pict")
     system("mv "+OUTDIR+"pict0001.ps "+OUTDIR+"psa2.timeseries.djf.ps") 
     system("mv "+OUTDIR+"pict0002.ps "+OUTDIR+"psa2.timeseries.mam.ps") 
     system("mv "+OUTDIR+"pict0003.ps "+OUTDIR+"psa2.timeseries.jja.ps") 
     system("mv "+OUTDIR+"pict0004.ps "+OUTDIR+"psa2.timeseries.son.ps") 
     system("mv "+OUTDIR+"pict0005.ps "+OUTDIR+"psa2.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"psa2.timeseries.mon.ps")    
        system("mv "+OUTDIR+"pict0007.ps "+OUTDIR+"npi.timeseries.ndjfm.ps")     
     else
        system("mv "+OUTDIR+"pict0006.ps "+OUTDIR+"npi.timeseries.ndjfm.ps")     
     end if
  end if
  print("Finished: psl.modes_indices.ncl")
end  
