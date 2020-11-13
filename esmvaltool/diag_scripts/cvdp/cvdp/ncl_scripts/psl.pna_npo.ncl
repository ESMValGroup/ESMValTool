; Calculates PNA and NPO (patterns and PC timeseries), as well as regressions 
; of those PC timeseries  onto ts, tas, and pr. Also calculates the NPI.
;
; Variables used: psl, ts, tas, and pr
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: psl.pna_npo.ncl")
  
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
;---------TAS Regressions coding-------------------------------------------------
  nsim_tas = numAsciiRow("namelist_byvar/namelist_trefht")
  na_tas = asciiread("namelist_byvar/namelist_trefht",(/nsim_tas/),"string")
  names_tas = new(nsim_tas,"string")
  paths_tas = new(nsim_tas,"string")
  syear_tas = new(nsim_tas,"integer",-999)
  eyear_tas = new(nsim_tas,"integer",-999)

  do gg = 0,nsim_tas-1
     names_tas(gg) = str_strip(str_get_field(na_tas(gg),1,delim))
     paths_tas(gg) = str_strip(str_get_field(na_tas(gg),2,delim))
     syear_tas(gg) = stringtointeger(str_strip(str_get_field(na_tas(gg),3,delim)))
     eyear_tas(gg) = stringtointeger(str_strip(str_get_field(na_tas(gg),4,delim)))
  end do
  delete(na_tas)
  nyr_tas = eyear_tas-syear_tas+1  
;---------PR Regressions coding-------------------------------------------------
  nsim_pr = numAsciiRow("namelist_byvar/namelist_prect")
  na_pr = asciiread("namelist_byvar/namelist_prect",(/nsim_pr/),"string")
  names_pr = new(nsim_pr,"string")
  paths_pr = new(nsim_pr,"string")
  syear_pr = new(nsim_pr,"integer",-999)
  eyear_pr = new(nsim_pr,"integer",-999)

  do gg = 0,nsim_pr-1
     names_pr(gg) = str_strip(str_get_field(na_pr(gg),1,delim))
     paths_pr(gg) = str_strip(str_get_field(na_pr(gg),2,delim))
     syear_pr(gg) = stringtointeger(str_strip(str_get_field(na_pr(gg),3,delim)))
     eyear_pr(gg) = stringtointeger(str_strip(str_get_field(na_pr(gg),4,delim)))
  end do
  delete(na_pr)
  nyr_pr = eyear_pr-syear_pr+1  
;-------------------------------------------------------------------------------------------------
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if

  wks_pna = gsn_open_wks(wks_type,getenv("OUTDIR")+"pna")
  wks_pna_pr = gsn_open_wks(wks_type,getenv("OUTDIR")+"pna.prreg")
  wks_pna_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"pna.timeseries")
  
  wks_npo = gsn_open_wks(wks_type,getenv("OUTDIR")+"npo")
  wks_npo_pr = gsn_open_wks(wks_type,getenv("OUTDIR")+"npo.prreg")
  wks_npo_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"npo.timeseries")
  
  wks_npi_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"npi.timeseries.ndjfm")
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_pna,"ncl_default")   
     gsn_define_colormap(wks_pna_ts,"ncl_default")
     gsn_define_colormap(wks_npo,"ncl_default")   
     gsn_define_colormap(wks_npo_ts,"ncl_default")
     gsn_define_colormap(wks_npi_ts,"ncl_default")
     gsn_define_colormap(wks_pna_pr,"MPL_BrBG")  
     gsn_define_colormap(wks_npo_pr,"MPL_BrBG")  
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_pna,"BlueDarkRed18")   
     gsn_define_colormap(wks_pna_ts,"ncl_default")
     gsn_define_colormap(wks_npo,"BlueDarkRed18")   
     gsn_define_colormap(wks_npo_ts,"ncl_default")
     gsn_define_colormap(wks_npi_ts,"ncl_default")
     gsn_define_colormap(wks_pna_pr,"MPL_BrBG")  
     gsn_define_colormap(wks_npo_pr,"MPL_BrBG")  
  end if

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
  reg_pna_pr_djf = new(nsim,"graphic")  
  reg_pna_pr_mam = new(nsim,"graphic")  
  reg_pna_pr_jja = new(nsim,"graphic")  
  reg_pna_pr_son = new(nsim,"graphic")  
  reg_pna_pr_ann = new(nsim,"graphic")  
  reg_pna_pr_mon = new(nsim,"graphic")   
  
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
  reg_npo_pr_djf = new(nsim,"graphic")  
  reg_npo_pr_mam = new(nsim,"graphic")  
  reg_npo_pr_jja = new(nsim,"graphic")  
  reg_npo_pr_son = new(nsim,"graphic")  
  reg_npo_pr_ann = new(nsim,"graphic")  
  reg_npo_pr_mon = new(nsim,"graphic")  
  
  xy_npi = new(nsim,"graphic")
  
  sstreg_frame = 1    ; sstreg_frame = flag to create regressions .ps/.png files. Created/used instead of sstreg_plot_flag
                      ; so that if sst regressions are not created for the last simulation listed that .ps/png files are created
  tasreg_frame = 1
  prreg_frame  = 1

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
        if (CLIMO_SYEAR.lt.0) then
           climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
        else
           climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        end if   
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
     if (any(ismissing((/syear(ee),syear_ts(ee),eyear(ee),eyear_ts(ee)/)))) then
        sstreg_plot_flag = 1
     else
        if (syear(ee).eq.syear_ts(ee)) then     ; check that the start and end years match for ts, trefht, and psl
           if (eyear(ee).eq.eyear_ts(ee)) then
              sstreg_plot_flag = 0
           else
              sstreg_plot_flag = 1
           end if
        else
           sstreg_plot_flag = 1
        end if
     end if 

     if (sstreg_plot_flag.eq.0) then 
        ; print("Data to be read in: "+paths_ts(ee)+" from "+syear_ts(ee)+":"+eyear_ts(ee))
        sst = data_read_in(paths_ts(ee),"TS",syear_ts(ee),eyear_ts(ee))
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
              check_custom_climo(names_ts(ee),syear_ts(ee),eyear_ts(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = sst
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(sst&time,1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if   
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
;---------TAS Regressions coding-------------------------------------------------
     if (any(ismissing((/syear(ee),syear_tas(ee),eyear(ee),eyear_tas(ee)/)))) then
        tasreg_plot_flag = 1
     else
        if (syear(ee).eq.syear_tas(ee)) then     ; check that the start and end years match for ts, tas, and psl
           if (eyear(ee).eq.eyear_tas(ee)) then
              tasreg_plot_flag = 0
           else
              tasreg_plot_flag = 1
           end if
        else
           tasreg_plot_flag = 1
        end if
        if (sstreg_plot_flag.eq.1) then          ; if the ts dataset is missing but the tas is not, do not
           tasreg_plot_flag = 1                  ; run through the tas calculations as both currently required
        end if
     end if 

     if (tasreg_plot_flag.eq.0) then 
        tas = data_read_in(paths_tas(ee),"TREFHT",syear_tas(ee),eyear_tas(ee))
        if (isatt(tas,"is_all_missing")) then
           tasreg_plot_flag = 1
           delete(tas)
        end if
        
        if (tasreg_plot_flag.eq.0) then     ; only continue if both PSL/TS fields are present
           d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
           basemap = d->LSMASK
           lsm = landsea_mask(basemap,tas&lat,tas&lon)
           tas = mask(tas,conform(tas,lsm,(/1,2/)).eq.0,False)
           delete(lsm)
           
           if (OPT_CLIMO.eq."Full") then
              tas = rmMonAnnCycTLL(tas)
           else
              check_custom_climo(names_tas(ee),syear_tas(ee),eyear_tas(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = tas
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(tas&time,1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if   
              delete(temp_arr)
              tas   = calcMonAnomTLL(tas,climo) 
              delete(climo)
           end if
;           tas = (/ dtrend_msg_n(ispan(0,dimsizes(tas&time)-1,1),tas,False,False,0) /)  
        
           tasT =  runave_n_Wrap(tas,3,0,0)               ; form DJF averages
           tasT(0,:,:) = (/ dim_avg_n(tas(:1,:,:),0) /)   
           tas_djf = tasT(0::12,:,:)
           tas_mam = tasT(3::12,:,:)
           tas_jja = tasT(6::12,:,:)                       ; form JJA averages
           tas_son = tasT(9::12,:,:)
           delete(tasT)
     
           tasV = runave_n_Wrap(tas,12,0,0)
           tas_ann = tasV(5::12,:,:)
           delete([/tasV/])
        end if
     end if
;---------PR Regressions coding-------------------------------------------------
     if (any(ismissing((/syear(ee),syear_pr(ee),eyear(ee),eyear_pr(ee)/)))) then
        prreg_plot_flag = 1
     else
        if (syear(ee).eq.syear_pr(ee)) then     ; check that the start and end years match for ts, tas, and psl
           if (eyear(ee).eq.eyear_pr(ee)) then
              prreg_plot_flag = 0
           else
              prreg_plot_flag = 1
           end if
        else
           prreg_plot_flag = 1
        end if
     end if 

     if (prreg_plot_flag.eq.0) then 
        pr = data_read_in(paths_pr(ee),"PRECT",syear_pr(ee),eyear_pr(ee))
        if (isatt(pr,"is_all_missing")) then
           prreg_plot_flag = 1
           delete(pr)
        end if
        
        if (prreg_plot_flag.eq.0) then     ; only continue if both PSL/TS fields are present
           if (OPT_CLIMO.eq."Full") then
              pr = rmMonAnnCycTLL(pr)
           else
              check_custom_climo(names_pr(ee),syear_pr(ee),eyear_pr(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = pr
              delete(temp_arr&time)
              temp_arr&time = cd_calendar(pr&time,1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if   
              delete(temp_arr)
              pr   = calcMonAnomTLL(pr,climo) 
              delete(climo)
           end if
;           pr = (/ dtrend_msg_n(ispan(0,dimsizes(pr&time)-1,1),pr,False,False,0) /)  
        
           prT =  runave_n_Wrap(pr,3,0,0)               ; form DJF averages
           prT(0,:,:) = (/ dim_avg_n(pr(:1,:,:),0) /)   
           pr_djf = prT(0::12,:,:)
           pr_mam = prT(3::12,:,:)
           pr_jja = prT(6::12,:,:)                       ; form JJA averages
           pr_son = prT(9::12,:,:)
           delete(prT)
     
           prV = runave_n_Wrap(pr,12,0,0)
           pr_ann = prV(5::12,:,:)
           delete([/prV/])
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

;------------------------------------------------------------------
     arr_djf_CW = SqrtCosWeight(arr_djf)
     arr_mam_CW = SqrtCosWeight(arr_mam)
     arr_jja_CW = SqrtCosWeight(arr_jja)
     arr_son_CW = SqrtCosWeight(arr_son)
     arr_ann_CW = SqrtCosWeight(arr_ann)
     if (COMPUTE_MODES_MON.eq."True") then   
        arr_mon_CW = SqrtCosWeight(arr)
     else
        if (isvar("arr")) then
           delete(arr)
        end if
        if (isvar("sst")) then
           delete(sst)
        end if 
        if (isvar("tas")) then
           delete(tas)
        end if 
        if (isvar("pr")) then
           delete(pr)
        end if 
     end if
;----------PNA/NPO calculations (EOF1/2 of NP PSL)----------------------------------------------------------
     evecv = eofunc(arr_djf_CW({lat|20:85},{lon|120:240},time|:),3,75)
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
        delete(sst_djf)
     end if   
     if (tasreg_plot_flag.eq.0) then
        pna_tas_djf = tas_djf(0,:,:)
        pna_tas_djf = (/ regCoef(pna_pc_djf,tas_djf(lat|:,lon|:,time|:)) /)
        npo_tas_djf = tas_djf(0,:,:)
        npo_tas_djf = (/ regCoef(npo_pc_djf,tas_djf(lat|:,lon|:,time|:)) /)
        delete(tas_djf)
     end if
     if (prreg_plot_flag.eq.0) then
        pna_pr_djf = pr_djf(0,:,:)
        pna_pr_djf = (/ regCoef(pna_pc_djf,pr_djf(lat|:,lon|:,time|:)) /)
        npo_pr_djf = pr_djf(0,:,:)
        npo_pr_djf = (/ regCoef(npo_pc_djf,pr_djf(lat|:,lon|:,time|:)) /)
        delete(pr_djf)
     end if  

     if (.not.ismissing(pna_djf({50},{185}))) then
        if (pna_djf({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_djf = pna_djf*-1.
           pna_pc_djf = pna_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_djf = pna_sst_djf*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              pna_tas_djf = pna_tas_djf*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              pna_pr_djf = pna_pr_djf*-1.
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
           if (tasreg_plot_flag.eq.0) then
              npo_tas_djf = npo_tas_djf*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              npo_pr_djf = npo_pr_djf*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(pna_pc_djf),False)
     if (sig_pcv(0)) then    ; if True then significant
        pna_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        pna_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     if (sig_pcv(1)) then    ; if True then significant
        npo_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%*"
     else
        npo_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(npi_ndjfm,pna_pc_djf) 
     copy_VarCoords(npi_ndjfm,npo_pc_djf)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_mam_CW({lat|20:85},{lon|120:240},time|:),3,75)
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
        delete(sst_mam)
     end if   
     if (tasreg_plot_flag.eq.0) then
        pna_tas_mam = tas_mam(0,:,:)
        pna_tas_mam = (/ regCoef(pna_pc_mam,tas_mam(lat|:,lon|:,time|:)) /)
        npo_tas_mam = tas_mam(0,:,:)
        npo_tas_mam = (/ regCoef(npo_pc_mam,tas_mam(lat|:,lon|:,time|:)) /)
        delete(tas_mam)
     end if
     if (prreg_plot_flag.eq.0) then
        pna_pr_mam = pr_mam(0,:,:)
        pna_pr_mam = (/ regCoef(pna_pc_mam,pr_mam(lat|:,lon|:,time|:)) /)
        npo_pr_mam = pr_mam(0,:,:)
        npo_pr_mam = (/ regCoef(npo_pc_mam,pr_mam(lat|:,lon|:,time|:)) /)
        delete(pr_mam)
     end if  

     if (.not.ismissing(pna_mam({50},{185}))) then
        if (pna_mam({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_mam = pna_mam*-1.
           pna_pc_mam = pna_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_mam = pna_sst_mam*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              pna_tas_mam = pna_tas_mam*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              pna_pr_mam = pna_pr_mam*-1.
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
           if (tasreg_plot_flag.eq.0) then
              npo_tas_mam = npo_tas_mam*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              npo_pr_mam = npo_pr_mam*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(pna_pc_mam),False)
     if (sig_pcv(0)) then    ; if True then significant
        pna_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        pna_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     if (sig_pcv(1)) then    ; if True then significant
        npo_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%*"
     else
        npo_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(npi_ndjfm,pna_pc_mam) 
     copy_VarCoords(npi_ndjfm,npo_pc_mam)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_jja_CW({lat|20:85},{lon|120:240},time|:),3,75)
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
        delete(sst_jja)
     end if   
     if (tasreg_plot_flag.eq.0) then
        pna_tas_jja = tas_jja(0,:,:)
        pna_tas_jja = (/ regCoef(pna_pc_jja,tas_jja(lat|:,lon|:,time|:)) /)
        npo_tas_jja = tas_jja(0,:,:)
        npo_tas_jja = (/ regCoef(npo_pc_jja,tas_jja(lat|:,lon|:,time|:)) /)
        delete(tas_jja)
     end if
     if (prreg_plot_flag.eq.0) then
        pna_pr_jja = pr_jja(0,:,:)
        pna_pr_jja = (/ regCoef(pna_pc_jja,pr_jja(lat|:,lon|:,time|:)) /)
        npo_pr_jja = pr_jja(0,:,:)
        npo_pr_jja = (/ regCoef(npo_pc_jja,pr_jja(lat|:,lon|:,time|:)) /)
        delete(pr_jja)
     end if  

     if (.not.ismissing(pna_jja({50},{185}))) then
        if (pna_jja({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_jja = pna_jja*-1.
           pna_pc_jja = pna_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_jja = pna_sst_jja*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              pna_tas_jja = pna_tas_jja*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              pna_pr_jja = pna_pr_jja*-1.
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
           if (tasreg_plot_flag.eq.0) then
              npo_tas_jja = npo_tas_jja*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              npo_pr_jja = npo_pr_jja*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(pna_pc_jja),False)
     if (sig_pcv(0)) then    ; if True then significant
        pna_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        pna_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     if (sig_pcv(1)) then    ; if True then significant
        npo_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%*"
     else
        npo_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(npi_ndjfm,pna_pc_jja) 
     copy_VarCoords(npi_ndjfm,npo_pc_jja)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_son_CW({lat|20:85},{lon|120:240},time|:),3,75)
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
        delete(sst_son)
     end if   
     if (tasreg_plot_flag.eq.0) then
        pna_tas_son = tas_son(0,:,:)
        pna_tas_son = (/ regCoef(pna_pc_son,tas_son(lat|:,lon|:,time|:)) /)
        npo_tas_son = tas_son(0,:,:)
        npo_tas_son = (/ regCoef(npo_pc_son,tas_son(lat|:,lon|:,time|:)) /)
        delete(tas_son)
     end if
     if (prreg_plot_flag.eq.0) then
        pna_pr_son = pr_son(0,:,:)
        pna_pr_son = (/ regCoef(pna_pc_son,pr_son(lat|:,lon|:,time|:)) /)
        npo_pr_son = pr_son(0,:,:)
        npo_pr_son = (/ regCoef(npo_pc_son,pr_son(lat|:,lon|:,time|:)) /)
        delete(pr_son)
     end if  

     if (.not.ismissing(pna_son({50},{185}))) then
        if (pna_son({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_son = pna_son*-1.
           pna_pc_son = pna_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_son = pna_sst_son*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              pna_tas_son = pna_tas_son*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              pna_pr_son = pna_pr_son*-1.
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
           if (tasreg_plot_flag.eq.0) then
              npo_tas_son = npo_tas_son*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              npo_pr_son = npo_pr_son*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(pna_pc_son),False)
     if (sig_pcv(0)) then    ; if True then significant
        pna_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        pna_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     if (sig_pcv(1)) then    ; if True then significant
        npo_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%*"
     else
        npo_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(npi_ndjfm,pna_pc_son) 
     copy_VarCoords(npi_ndjfm,npo_pc_son)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_ann_CW({lat|20:85},{lon|120:240},time|:),3,75)
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
        delete(sst_ann)
     end if   
     if (tasreg_plot_flag.eq.0) then
        pna_tas_ann = tas_ann(0,:,:)
        pna_tas_ann = (/ regCoef(pna_pc_ann,tas_ann(lat|:,lon|:,time|:)) /)
        npo_tas_ann = tas_ann(0,:,:)
        npo_tas_ann = (/ regCoef(npo_pc_ann,tas_ann(lat|:,lon|:,time|:)) /)
        delete(tas_ann)
     end if
     if (prreg_plot_flag.eq.0) then
        pna_pr_ann = pr_ann(0,:,:)
        pna_pr_ann = (/ regCoef(pna_pc_ann,pr_ann(lat|:,lon|:,time|:)) /)
        npo_pr_ann = pr_ann(0,:,:)
        npo_pr_ann = (/ regCoef(npo_pc_ann,pr_ann(lat|:,lon|:,time|:)) /)
        delete(pr_ann)
     end if  

     if (.not.ismissing(pna_ann({50},{185}))) then
        if (pna_ann({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           pna_ann = pna_ann*-1.
           pna_pc_ann = pna_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              pna_sst_ann = pna_sst_ann*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              pna_tas_ann = pna_tas_ann*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              pna_pr_ann = pna_pr_ann*-1.
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
           if (tasreg_plot_flag.eq.0) then
              npo_tas_ann = npo_tas_ann*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              npo_pr_ann = npo_pr_ann*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(pna_pc_ann),False)
     if (sig_pcv(0)) then    ; if True then significant
        pna_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        pna_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     if (sig_pcv(1)) then    ; if True then significant
        npo_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%*"
     else
        npo_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(npi_ndjfm,pna_pc_ann) 
     copy_VarCoords(npi_ndjfm,npo_pc_ann)
     delete([/evecv,pcts/])
     
     if (COMPUTE_MODES_MON.eq."True") then          
        evecv = eofunc(arr_mon_CW({lat|20:85},{lon|120:240},time|:),3,75)
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
           delete(sst)
        end if 
        if (tasreg_plot_flag.eq.0) then
           pna_tas_mon = tas(0,:,:)
           pna_tas_mon = (/ regCoef(pna_pc_mon,tas(lat|:,lon|:,time|:)) /)
           npo_tas_mon = tas(0,:,:)
           npo_tas_mon = (/ regCoef(npo_pc_mon,tas(lat|:,lon|:,time|:)) /)
           delete(tas)
        end if 
        if (prreg_plot_flag.eq.0) then
           pna_pr_mon = pr(0,:,:)
           pna_pr_mon = (/ regCoef(pna_pc_mon,pr(lat|:,lon|:,time|:)) /)
           npo_pr_mon = pr(0,:,:)
           npo_pr_mon = (/ regCoef(npo_pc_mon,pr(lat|:,lon|:,time|:)) /)
           delete(pr)
        end if 
        if (.not.ismissing(pna_mon({50},{185}))) then
           if (pna_mon({50},{185}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              pna_mon = pna_mon*-1.
              pna_pc_mon = pna_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 pna_sst_mon = pna_sst_mon*-1.
              end if
              if (tasreg_plot_flag.eq.0) then
                 pna_tas_mon = pna_tas_mon*-1.
              end if
              if (prreg_plot_flag.eq.0) then
                 pna_pr_mon = pna_pr_mon*-1.
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
              if (tasreg_plot_flag.eq.0) then
                 npo_tas_mon = npo_tas_mon*-1.
              end if
              if (prreg_plot_flag.eq.0) then
                 npo_pr_mon = npo_pr_mon*-1.
              end if
           end if
        end if
        sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(pna_pc_mon),False)
        if (sig_pcv(0)) then    ; if True then significant
           pna_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
        else
           pna_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
        end if
        if (sig_pcv(1)) then    ; if True then significant
           npo_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%*"
        else
           npo_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(1)))+"%"
        end if
        delete(sig_pcv)
        pna_pc_mon!0 = "time"
        pna_pc_mon&time = arr&time   
        npo_pc_mon!0 = "time"
        npo_pc_mon&time = arr&time   
        delete([/evecv,pcts,arr,arr_mon_CW/])    
     end if
     delete([/arr_djf_CW,arr_mam_CW,arr_jja_CW,arr_son_CW,arr_ann_CW/])
     delete([/arr_djf,arr_mam,arr_jja,arr_son,arr_ann/])
;------------------------------------------------------------------------------------------------------
     if (sstreg_frame.eq.1.and.sstreg_plot_flag.eq.0) then    ; sstreg_frame = flag to create regressions .ps/.png files 
        sstreg_frame = 0
     end if
     if (tasreg_frame.eq.1.and.tasreg_plot_flag.eq.0) then    ; tasreg_frame = flag to create regressions .ps/.png files 
        tasreg_frame = 0
     end if  
     if (prreg_frame.eq.1.and.prreg_plot_flag.eq.0) then      ; prreg_frame = flag to create regressions .ps/.png files 
        prreg_frame = 0
     end if   
;------------------------------------------------------------------------------------------------------
     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.psl.pna_npo."+syear(ee)+"-"+eyear(ee)+".nc"
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
        z->npi_ndjfm = npi_ndjfm

        z->pna_timeseries_djf = set_varAtts(pna_pc_djf,"PNA normalized principal component timeseries (DJF)","1","")
        z->pna_timeseries_mam = set_varAtts(pna_pc_mam,"PNA normalized principal component timeseries (MAM)","1","")
        z->pna_timeseries_jja = set_varAtts(pna_pc_jja,"PNA normalized principal component timeseries (JJA)","1","")
        z->pna_timeseries_son = set_varAtts(pna_pc_son,"PNA normalized principal component timeseries (SON)","1","")
        z->pna_timeseries_ann = set_varAtts(pna_pc_ann,"PNA normalized principal component timeseries (annual)","1","")
        
        z->npo_timeseries_djf = set_varAtts(npo_pc_djf,"NPO normalized principal component timeseries (DJF)","1","")
        z->npo_timeseries_mam = set_varAtts(npo_pc_mam,"NPO normalized principal component timeseries (MAM)","1","")
        z->npo_timeseries_jja = set_varAtts(npo_pc_jja,"NPO normalized principal component timeseries (JJA)","1","")
        z->npo_timeseries_son = set_varAtts(npo_pc_son,"NPO normalized principal component timeseries (SON)","1","")
        z->npo_timeseries_ann = set_varAtts(npo_pc_ann,"NPO normalized principal component timeseries (annual)","1","")
        
        z->pna_pattern_djf = set_varAtts(pna_djf,"PNA spatial pattern (DJF)","","")
        z->pna_pattern_mam = set_varAtts(pna_mam,"PNA spatial pattern (MAM)","","")
        z->pna_pattern_jja = set_varAtts(pna_jja,"PNA spatial pattern (JJA)","","")
        z->pna_pattern_son = set_varAtts(pna_son,"PNA spatial pattern (SON)","","")
        z->pna_pattern_ann = set_varAtts(pna_ann,"PNA spatial pattern (annual)","","")
        
        z->npo_pattern_djf = set_varAtts(npo_djf,"NPO spatial pattern (DJF)","","")
        z->npo_pattern_mam = set_varAtts(npo_mam,"NPO spatial pattern (MAM)","","")
        z->npo_pattern_jja = set_varAtts(npo_jja,"NPO spatial pattern (JJA)","","")
        z->npo_pattern_son = set_varAtts(npo_son,"NPO spatial pattern (SON)","","")
        z->npo_pattern_ann = set_varAtts(npo_ann,"NPO spatial pattern (annual)","","")

        if (COMPUTE_MODES_MON.eq."True") then   
           z->pna_timeseries_mon = set_varAtts(pna_pc_mon,"PNA principal component timeseries (monthly)","","")
           z->npo_timeseries_mon = set_varAtts(npo_pc_mon,"NPO principal component timeseries (monthly)","","")
           z->pna_pattern_mon = set_varAtts(pna_mon,"PNA spatial pattern (monthly)","","")
           z->npo_pattern_mon = set_varAtts(npo_mon,"NPO spatial pattern (monthly)","","")
        end if
        delete(z)
        delete([/modname,fn/])        

        if (sstreg_plot_flag.eq.0) then
           modname = str_sub_str(names_ts(ee)," ","_")
           bc = (/"/","'","(",")"/)
           do gg = 0,dimsizes(bc)-1
              modname = str_sub_str(modname,bc(gg),"_")
           end do
           fn = getenv("OUTDIR")+modname+".cvdp_data.psl.pna_npo.ts."+syear_ts(ee)+"-"+eyear_ts(ee)+".nc"
           if (.not.isfilepresent2(fn)) then
              z = addfile(fn,"c")
              z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
              z@notes = "Data from "+names_ts(ee)+" from "+syear_ts(ee)+"-"+eyear_ts(ee)
              if (OPT_CLIMO.eq."Full") then
                 z@climatology = syear_ts(ee)+"-"+eyear_ts(ee)+" climatology removed prior to all calculations (other than means)"
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
           z->pna_sst_regression_djf = set_varAtts(pna_sst_djf,"sst regression onto PNA principal component timeseries (DJF)","","")    
           z->pna_sst_regression_mam = set_varAtts(pna_sst_mam,"sst regression onto PNA principal component timeseries (MAM)","","")    
           z->pna_sst_regression_jja = set_varAtts(pna_sst_jja,"sst regression onto PNA principal component timeseries (JJA)","","")    
           z->pna_sst_regression_son = set_varAtts(pna_sst_son,"sst regression onto PNA principal component timeseries (SON)","","")    
           z->pna_sst_regression_ann = set_varAtts(pna_sst_ann,"sst regression onto PNA principal component timeseries (annual)","","")    
        
           z->npo_sst_regression_djf = set_varAtts(npo_sst_djf,"sst regression onto NPO principal component timeseries (DJF)","","")    
           z->npo_sst_regression_mam = set_varAtts(npo_sst_mam,"sst regression onto NPO principal component timeseries (MAM)","","")    
           z->npo_sst_regression_jja = set_varAtts(npo_sst_jja,"sst regression onto NPO principal component timeseries (JJA)","","")    
           z->npo_sst_regression_son = set_varAtts(npo_sst_son,"sst regression onto NPO principal component timeseries (SON)","","")    
           z->npo_sst_regression_ann = set_varAtts(npo_sst_ann,"sst regression onto NPO principal component timeseries (annual)","","")    
           if (COMPUTE_MODES_MON.eq."True") then     
              z->pna_sst_regression_mon = set_varAtts(pna_sst_mon,"sst regression onto PNA principal component timeseries (monthly)","","")    
              z->npo_sst_regression_mon = set_varAtts(npo_sst_mon,"sst regression onto NPO principal component timeseries (monthly)","","")    
           end if
           delete(z)
           delete([/modname,fn/])
        end if
        if (tasreg_plot_flag.eq.0) then
           modname = str_sub_str(names_tas(ee)," ","_")
           bc = (/"/","'","(",")"/)
           do gg = 0,dimsizes(bc)-1
              modname = str_sub_str(modname,bc(gg),"_")
           end do
           fn = getenv("OUTDIR")+modname+".cvdp_data.psl.pna_npo.tas."+syear_tas(ee)+"-"+eyear_tas(ee)+".nc"
           if (.not.isfilepresent2(fn)) then
              z = addfile(fn,"c")
              z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
              z@notes = "Data from "+names_tas(ee)+" from "+syear_tas(ee)+"-"+eyear_tas(ee)
              if (OPT_CLIMO.eq."Full") then
                 z@climatology = syear_tas(ee)+"-"+eyear_tas(ee)+" climatology removed prior to all calculations (other than means)"
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
           z->pna_tas_regression_djf = set_varAtts(pna_tas_djf,"tas regression onto PNA principal component timeseries (DJF)","","")    
           z->pna_tas_regression_mam = set_varAtts(pna_tas_mam,"tas regression onto PNA principal component timeseries (MAM)","","")    
           z->pna_tas_regression_jja = set_varAtts(pna_tas_jja,"tas regression onto PNA principal component timeseries (JJA)","","")    
           z->pna_tas_regression_son = set_varAtts(pna_tas_son,"tas regression onto PNA principal component timeseries (SON)","","")    
           z->pna_tas_regression_ann = set_varAtts(pna_tas_ann,"tas regression onto PNA principal component timeseries (annual)","","")    
        
           z->npo_tas_regression_djf = set_varAtts(npo_tas_djf,"tas regression onto NPO principal component timeseries (DJF)","","")    
           z->npo_tas_regression_mam = set_varAtts(npo_tas_mam,"tas regression onto NPO principal component timeseries (MAM)","","")    
           z->npo_tas_regression_jja = set_varAtts(npo_tas_jja,"tas regression onto NPO principal component timeseries (JJA)","","")    
           z->npo_tas_regression_son = set_varAtts(npo_tas_son,"tas regression onto NPO principal component timeseries (SON)","","")    
           z->npo_tas_regression_ann = set_varAtts(npo_tas_ann,"tas regression onto NPO principal component timeseries (annual)","","")    
           if (COMPUTE_MODES_MON.eq."True") then     
              z->pna_tas_regression_mon = set_varAtts(pna_tas_mon,"tas regression onto PNA principal component timeseries (monthly)","","")    
              z->npo_tas_regression_mon = set_varAtts(npo_tas_mon,"tas regression onto NPO principal component timeseries (monthly)","","")    
           end if
           delete(z)
           delete([/modname,fn/])
        end if
        if (prreg_plot_flag.eq.0) then
           modname = str_sub_str(names_pr(ee)," ","_")
           bc = (/"/","'","(",")"/)
           do gg = 0,dimsizes(bc)-1
              modname = str_sub_str(modname,bc(gg),"_")
           end do
           fn = getenv("OUTDIR")+modname+".cvdp_data.psl.pna_npo.pr."+syear_pr(ee)+"-"+eyear_pr(ee)+".nc"
           if (.not.isfilepresent2(fn)) then
              z = addfile(fn,"c")
              z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
              z@notes = "Data from "+names_pr(ee)+" from "+syear_pr(ee)+"-"+eyear_pr(ee)
              if (OPT_CLIMO.eq."Full") then
                 z@climatology = syear_pr(ee)+"-"+eyear_pr(ee)+" climatology removed prior to all calculations (other than means)"
              else
                 z@climatology = CLIMO_SYEAR+"-"+CLIMO_EYEAR+" climatology removed prior to all calculations (other than means)"
              end if
              z@Conventions = "CF-1.6"
           else
              z = addfile(fn,"w")
           end if
           z->pna_pr_regression_djf = set_varAtts(pna_pr_djf,"pr regression onto PNA principal component timeseries (DJF)","","")    
           z->pna_pr_regression_mam = set_varAtts(pna_pr_mam,"pr regression onto PNA principal component timeseries (MAM)","","")    
           z->pna_pr_regression_jja = set_varAtts(pna_pr_jja,"pr regression onto PNA principal component timeseries (JJA)","","")    
           z->pna_pr_regression_son = set_varAtts(pna_pr_son,"pr regression onto PNA principal component timeseries (SON)","","")    
           z->pna_pr_regression_ann = set_varAtts(pna_pr_ann,"pr regression onto PNA principal component timeseries (annual)","","")    
        
           z->npo_pr_regression_djf = set_varAtts(npo_pr_djf,"pr regression onto NPO principal component timeseries (DJF)","","")    
           z->npo_pr_regression_mam = set_varAtts(npo_pr_mam,"pr regression onto NPO principal component timeseries (MAM)","","")    
           z->npo_pr_regression_jja = set_varAtts(npo_pr_jja,"pr regression onto NPO principal component timeseries (JJA)","","")    
           z->npo_pr_regression_son = set_varAtts(npo_pr_son,"pr regression onto NPO principal component timeseries (SON)","","")    
           z->npo_pr_regression_ann = set_varAtts(npo_pr_ann,"pr regression onto NPO principal component timeseries (annual)","","")    
           if (COMPUTE_MODES_MON.eq."True") then     
              z->pna_pr_regression_mon = set_varAtts(pna_pr_mon,"pr regression onto PNA principal component timeseries (monthly)","","")    
              z->npo_pr_regression_mon = set_varAtts(npo_pr_mon,"pr regression onto NPO principal component timeseries (monthly)","","")    
           end if
           delete(z)
           delete([/modname,fn/])
        end if
     end if     
;========================================================================
     res = True
     res@mpGeophysicalLineColor = "gray42"
     res@mpGeophysicalLineThicknessF = 2.     
     res@mpGridAndLimbOn = False
     res@mpFillOn = False
     res@mpOutlineOn = True  
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

     res4 = res    ; res4 = pr regression resources     
     if (COLORMAP.eq.0) then
        res4@cnLevels := fspan(-.7,.7,15)  
     else
        res4@cnLevels := fspan(-.5,.5,11)     
     end if

     res2 = True
     res2@gsnDraw      = False
     res2@gsnFrame     = False
     res2@cnLevelSelectionMode = "ExplicitLevels"
     res2@cnLevels = res@cnLevels

     res2@cnLineLabelsOn = False
     res2@cnFillOn        = True
     res2@cnLinesOn       = False
     res2@cnFillMode = "AreaFill"
     res2@lbLabelBarOn    = False
     res2@cnInfoLabelOn = False
     res2@gsnRightString = ""
     res2@gsnLeftString = "" 
     res2@gsnCenterString = ""   
     res2@gsnAddCyclic = True
     
     res@gsnRightString = pna_djf@pcvar
     map_pna_djf(ee) = gsn_csm_contour_map_polar(wks_pna,pna_djf,res)
     res@gsnRightString = pna_mam@pcvar
     map_pna_mam(ee) = gsn_csm_contour_map_polar(wks_pna,pna_mam,res)
     res@gsnRightString = pna_jja@pcvar
     map_pna_jja(ee) = gsn_csm_contour_map_polar(wks_pna,pna_jja,res)
     res@gsnRightString = pna_son@pcvar
     map_pna_son(ee) = gsn_csm_contour_map_polar(wks_pna,pna_son,res)
     res@gsnRightString = pna_ann@pcvar
     map_pna_ann(ee) = gsn_csm_contour_map_polar(wks_pna,pna_ann,res)     
     delete([/pna_djf,pna_mam,pna_jja,pna_son,pna_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = pna_mon@pcvar
        map_pna_mon(ee) = gsn_csm_contour_map_polar(wks_pna,pna_mon,res)
        delete([/pna_mon/])
     end if
     
     res@cnLevels = (/-4,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4/)
     res@gsnRightString = npo_djf@pcvar
     map_npo_djf(ee) = gsn_csm_contour_map_polar(wks_npo,npo_djf,res)
     res@gsnRightString = npo_mam@pcvar
     map_npo_mam(ee) = gsn_csm_contour_map_polar(wks_npo,npo_mam,res)
     res@gsnRightString = npo_jja@pcvar
     map_npo_jja(ee) = gsn_csm_contour_map_polar(wks_npo,npo_jja,res)
     res@gsnRightString = npo_son@pcvar
     map_npo_son(ee) = gsn_csm_contour_map_polar(wks_npo,npo_son,res)
     res@gsnRightString = npo_ann@pcvar
     map_npo_ann(ee) = gsn_csm_contour_map_polar(wks_npo,npo_ann,res)     
     delete([/npo_djf,npo_mam,npo_jja,npo_son,npo_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = npo_mon@pcvar
        map_npo_mon(ee) = gsn_csm_contour_map_polar(wks_npo,npo_mon,res)
        delete([/npo_mon/])
     end if
     
     if (sstreg_plot_flag.eq.0) then
        res@cnLevels := fspan(-.7,.7,15)
        if (tasreg_plot_flag.eq.0) then
           if (names_ts(ee).eq.names_tas(ee)) then
              res@gsnCenterString = names_ts(ee)
           else
              res@gsnCenterString = names_ts(ee)+" / "+names_tas(ee)
           end if    
        else
           res@gsnCenterString = names_ts(ee)
        end if
        res@gsnRightString = ""
        reg_pna_djf(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_djf,res)
        reg_pna_mam(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_mam,res)
        reg_pna_jja(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_jja,res)
        reg_pna_son(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_son,res)
        reg_pna_ann(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_ann,res)     
        delete([/pna_sst_djf,pna_sst_mam,pna_sst_jja,pna_sst_son,pna_sst_ann/])
        if (tasreg_plot_flag.eq.0) then
           o_djf = gsn_csm_contour(wks_pna,pna_tas_djf,res2)
           o_mam = gsn_csm_contour(wks_pna,pna_tas_mam,res2)
           o_jja = gsn_csm_contour(wks_pna,pna_tas_jja,res2)
           o_son = gsn_csm_contour(wks_pna,pna_tas_son,res2)
           o_ann = gsn_csm_contour(wks_pna,pna_tas_ann,res2)
           delete([/pna_tas_djf,pna_tas_mam,pna_tas_jja,pna_tas_son,pna_tas_ann/])
           overlay(reg_pna_djf(ee),o_djf)
           overlay(reg_pna_mam(ee),o_mam)
           overlay(reg_pna_jja(ee),o_jja)
           overlay(reg_pna_son(ee),o_son)
           overlay(reg_pna_ann(ee),o_ann)
           delete([/o_djf,o_mam,o_jja,o_son,o_ann/])           
        end if
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_pna_mon(ee) = gsn_csm_contour_map_polar(wks_pna,pna_sst_mon,res)
           delete([/pna_sst_mon/])
           if (tasreg_plot_flag.eq.0) then
              o_mon = gsn_csm_contour(wks_pna,pna_tas_mon,res2)
              overlay(reg_pna_mon(ee),o_mon)
              delete([/o_mon,pna_tas_mon/])
           end if
        end if
     
        reg_npo_djf(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_djf,res)
        reg_npo_mam(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_mam,res)
        reg_npo_jja(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_jja,res)
        reg_npo_son(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_son,res)
        reg_npo_ann(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_ann,res)     
        delete([/npo_sst_djf,npo_sst_mam,npo_sst_jja,npo_sst_son,npo_sst_ann/])
        if (tasreg_plot_flag.eq.0) then
           o_djf = gsn_csm_contour(wks_npo,npo_tas_djf,res2)
           o_mam = gsn_csm_contour(wks_npo,npo_tas_mam,res2)
           o_jja = gsn_csm_contour(wks_npo,npo_tas_jja,res2)
           o_son = gsn_csm_contour(wks_npo,npo_tas_son,res2)
           o_ann = gsn_csm_contour(wks_npo,npo_tas_ann,res2)
           delete([/npo_tas_djf,npo_tas_mam,npo_tas_jja,npo_tas_son,npo_tas_ann/])
           overlay(reg_npo_djf(ee),o_djf)
           overlay(reg_npo_mam(ee),o_mam)
           overlay(reg_npo_jja(ee),o_jja)
           overlay(reg_npo_son(ee),o_son)
           overlay(reg_npo_ann(ee),o_ann)
           delete([/o_djf,o_mam,o_jja,o_son,o_ann/])           
        end if
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_npo_mon(ee) = gsn_csm_contour_map_polar(wks_npo,npo_sst_mon,res)
           delete([/npo_sst_mon/])
           if (tasreg_plot_flag.eq.0) then
              o_mon = gsn_csm_contour(wks_npo,npo_tas_mon,res2)
              overlay(reg_npo_mon(ee),o_mon)
              delete([/o_mon,npo_tas_mon/])
           end if
        end if 
     end if

     if (prreg_plot_flag.eq.0) then     ; PR regressions
        res4@gsnRightString = ""
        res4@gsnCenterString = names_pr(ee)
        reg_pna_pr_djf(ee) = gsn_csm_contour_map_polar(wks_pna_pr,pna_pr_djf,res4)
        reg_pna_pr_mam(ee) = gsn_csm_contour_map_polar(wks_pna_pr,pna_pr_mam,res4)
        reg_pna_pr_jja(ee) = gsn_csm_contour_map_polar(wks_pna_pr,pna_pr_jja,res4)
        reg_pna_pr_son(ee) = gsn_csm_contour_map_polar(wks_pna_pr,pna_pr_son,res4)
        reg_pna_pr_ann(ee) = gsn_csm_contour_map_polar(wks_pna_pr,pna_pr_ann,res4)     
        delete([/pna_pr_djf,pna_pr_mam,pna_pr_jja,pna_pr_son,pna_pr_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_pna_pr_mon(ee) = gsn_csm_contour_map_polar(wks_pna_pr,pna_pr_mon,res4)
           delete([/pna_pr_mon/])
        end if
     
        reg_npo_pr_djf(ee) = gsn_csm_contour_map_polar(wks_npo_pr,npo_pr_djf,res4)
        reg_npo_pr_mam(ee) = gsn_csm_contour_map_polar(wks_npo_pr,npo_pr_mam,res4)
        reg_npo_pr_jja(ee) = gsn_csm_contour_map_polar(wks_npo_pr,npo_pr_jja,res4)
        reg_npo_pr_son(ee) = gsn_csm_contour_map_polar(wks_npo_pr,npo_pr_son,res4)
        reg_npo_pr_ann(ee) = gsn_csm_contour_map_polar(wks_npo_pr,npo_pr_ann,res4)     
        delete([/npo_pr_djf,npo_pr_mam,npo_pr_jja,npo_pr_son,npo_pr_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_npo_pr_mon(ee) = gsn_csm_contour_map_polar(wks_npo_pr,npo_pr_mon,res4)
           delete([/npo_pr_mon/])
        end if
     end if


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
     xyres@gsnCenterStringOrthogonalPosF = 0.025
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
     xyres@gsnCenterStringOrthogonalPosF = 0.025

     xyres@gsnCenterString = names(ee)
     
     xyresmon = xyres
     xyresmon@gsnXYBarChart = False            
     xyresmon@xyLineThicknessF = .1 

     xy_pna_djf(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_djf)),pna_pc_djf,xyres)   ; use standardized timeseries  
     xy_pna_mam(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_mam)),pna_pc_mam,xyres)   ; use standardized timeseries  
     xy_pna_jja(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_jja)),pna_pc_jja,xyres)   ; use standardized timeseries  
     xy_pna_son(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_son)),pna_pc_son,xyres)   ; use standardized timeseries  
     xy_pna_ann(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee),dimsizes(pna_pc_ann)),pna_pc_ann,xyres)   ; use standardized timeseries  
     delete([/pna_pc_djf,pna_pc_mam,pna_pc_jja,pna_pc_son,pna_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then
        xy_pna_mon(ee) = gsn_csm_xy(wks_pna_ts,fspan(syear(ee),eyear(ee)+.91667,dimsizes(pna_pc_mon)),pna_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/pna_pc_mon/])
     end if

     xy_npo_djf(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_djf)),npo_pc_djf,xyres)   ; use standardized timeseries  
     xy_npo_mam(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_mam)),npo_pc_mam,xyres)   ; use standardized timeseries  
     xy_npo_jja(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_jja)),npo_pc_jja,xyres)   ; use standardized timeseries  
     xy_npo_son(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_son)),npo_pc_son,xyres)   ; use standardized timeseries  
     xy_npo_ann(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee),dimsizes(npo_pc_ann)),npo_pc_ann,xyres)   ; use standardized timeseries  
     delete([/npo_pc_djf,npo_pc_mam,npo_pc_jja,npo_pc_son,npo_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then    
        xy_npo_mon(ee) = gsn_csm_xy(wks_npo_ts,fspan(syear(ee),eyear(ee)+.91667,dimsizes(npo_pc_mon)),npo_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/npo_pc_mon/])
     end if

     xy_npi(ee) = gsn_csm_xy(wks_npi_ts,fspan(syear(ee),eyear(ee),dimsizes(npi_ndjfm)),npi_ndjfm,xyres)   ; throw NPI into wks_psa2_ts workstation
     delete(npi_ndjfm)   
     delete(sstreg_plot_flag)  
  end do
  
  if (isvar("clim_syear")) then
     delete(clim_syear)
  end if
  if (isvar("clim_eyear")) then
     delete(clim_eyear)
  end if

  panres = True
  panres@gsnMaximize = True
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.55
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
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol) 
  
  panres@txString = "PNA (DJF)" 
  gsn_panel2(wks_pna,map_pna_djf,(/nrow,ncol/),panres)
  delete(map_pna_djf)
  panres@txString = "PNA (MAM)"
  gsn_panel2(wks_pna,map_pna_mam,(/nrow,ncol/),panres)
  delete(map_pna_mam)
  panres@txString = "PNA (JJA)"
  gsn_panel2(wks_pna,map_pna_jja,(/nrow,ncol/),panres)
  delete(map_pna_jja)
  panres@txString = "PNA (SON)"
  gsn_panel2(wks_pna,map_pna_son,(/nrow,ncol/),panres)
  delete(map_pna_son)
  panres@txString = "PNA (Annual)"
  gsn_panel2(wks_pna,map_pna_ann,(/nrow,ncol/),panres)
  delete(map_pna_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "PNA (Monthly)"
     gsn_panel2(wks_pna,map_pna_mon,(/nrow,ncol/),panres)
     delete(map_pna_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     if (tasreg_frame.eq.0) then
        txt0 = "SST/TAS"
     else
        txt0 = "SST"
     end if
     panres@txString = "PNA "+txt0+" Regressions (DJF)" 
     gsn_panel2(wks_pna,reg_pna_djf,(/nrow,ncol/),panres)
     delete(reg_pna_djf)
     panres@txString = "PNA "+txt0+" Regressions (MAM)"
     gsn_panel2(wks_pna,reg_pna_mam,(/nrow,ncol/),panres)
     delete(reg_pna_mam)
     panres@txString = "PNA "+txt0+" Regressions (JJA)"
     gsn_panel2(wks_pna,reg_pna_jja,(/nrow,ncol/),panres)
     delete(reg_pna_jja)
     panres@txString = "PNA "+txt0+" Regressions (SON)"
     gsn_panel2(wks_pna,reg_pna_son,(/nrow,ncol/),panres)
     delete(reg_pna_son)
     panres@txString = "PNA "+txt0+" Regressions (Annual)"
     gsn_panel2(wks_pna,reg_pna_ann,(/nrow,ncol/),panres)
     delete(reg_pna_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "PNA "+txt0+" Regressions (Monthly)"
        gsn_panel2(wks_pna,reg_pna_mon,(/nrow,ncol/),panres)
        delete(reg_pna_mon)
     end if
     delete(wks_pna)
  end if
  if (prreg_frame.eq.0) then
     panres@txString = "PNA PR Regressions (DJF)" 
     gsn_panel2(wks_pna_pr,reg_pna_pr_djf,(/nrow,ncol/),panres)
     delete(reg_pna_pr_djf)
     panres@txString = "PNA PR Regressions (MAM)"
     gsn_panel2(wks_pna_pr,reg_pna_pr_mam,(/nrow,ncol/),panres)
     delete(reg_pna_pr_mam)
     panres@txString = "PNA PR Regressions (JJA)"
     gsn_panel2(wks_pna_pr,reg_pna_pr_jja,(/nrow,ncol/),panres)
     delete(reg_pna_pr_jja)
     panres@txString = "PNA PR Regressions (SON)"
     gsn_panel2(wks_pna_pr,reg_pna_pr_son,(/nrow,ncol/),panres)
     delete(reg_pna_pr_son)
     panres@txString = "PNA PR Regressions (Annual)"
     gsn_panel2(wks_pna_pr,reg_pna_pr_ann,(/nrow,ncol/),panres)
     delete(reg_pna_pr_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "PNA PR Regressions (Monthly)"
        gsn_panel2(wks_pna_pr,reg_pna_pr_mon,(/nrow,ncol/),panres)
        delete(reg_pna_pr_mon)
     end if
     delete(wks_pna_pr)  
  end if
  
  panres@txString = "NPO (DJF)" 
  gsn_panel2(wks_npo,map_npo_djf,(/nrow,ncol/),panres)
  delete(map_npo_djf)
  panres@txString = "NPO (MAM)"
  gsn_panel2(wks_npo,map_npo_mam,(/nrow,ncol/),panres)
  delete(map_npo_mam)
  panres@txString = "NPO (JJA)"
  gsn_panel2(wks_npo,map_npo_jja,(/nrow,ncol/),panres)
  delete(map_npo_jja)
  panres@txString = "NPO (SON)"
  gsn_panel2(wks_npo,map_npo_son,(/nrow,ncol/),panres)
  delete(map_npo_son)
  panres@txString = "NPO (Annual)"
  gsn_panel2(wks_npo,map_npo_ann,(/nrow,ncol/),panres)
  delete(map_npo_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "NPO (Monthly)"
     gsn_panel2(wks_npo,map_npo_mon,(/nrow,ncol/),panres)
     delete(map_npo_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     if (tasreg_frame.eq.0) then
        txt0 = "SST/TAS"
     else
        txt0 = "SST"
     end if
     panres@txString = "NPO "+txt0+" Regressions (DJF)" 
     gsn_panel2(wks_npo,reg_npo_djf,(/nrow,ncol/),panres)
     delete(reg_npo_djf)
     panres@txString = "NPO "+txt0+" Regressions (MAM)"
     gsn_panel2(wks_npo,reg_npo_mam,(/nrow,ncol/),panres)
     delete(reg_npo_mam)
     panres@txString = "NPO "+txt0+" Regressions (JJA)"
     gsn_panel2(wks_npo,reg_npo_jja,(/nrow,ncol/),panres)
     delete(reg_npo_jja)
     panres@txString = "NPO "+txt0+" Regressions (SON)"
     gsn_panel2(wks_npo,reg_npo_son,(/nrow,ncol/),panres)
     delete(reg_npo_son)
     panres@txString = "NPO "+txt0+" Regressions (Annual)"
     gsn_panel2(wks_npo,reg_npo_ann,(/nrow,ncol/),panres)
     delete(reg_npo_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NPO "+txt0+" Regressions (Monthly)"
        gsn_panel2(wks_npo,reg_npo_mon,(/nrow,ncol/),panres)
        delete(reg_npo_mon)
     end if
     delete(wks_npo)
  end if
  if (prreg_frame.eq.0) then
     panres@txString = "NPO PR Regressions (DJF)" 
     gsn_panel2(wks_npo_pr,reg_npo_pr_djf,(/nrow,ncol/),panres)
     delete(reg_npo_pr_djf)
     panres@txString = "NPO PR Regressions (MAM)"
     gsn_panel2(wks_npo_pr,reg_npo_pr_mam,(/nrow,ncol/),panres)
     delete(reg_npo_pr_mam)
     panres@txString = "NPO PR Regressions (JJA)"
     gsn_panel2(wks_npo_pr,reg_npo_pr_jja,(/nrow,ncol/),panres)
     delete(reg_npo_pr_jja)
     panres@txString = "NPO PR Regressions (SON)"
     gsn_panel2(wks_npo_pr,reg_npo_pr_son,(/nrow,ncol/),panres)
     delete(reg_npo_pr_son)
     panres@txString = "NPO PR Regressions (Annual)"
     gsn_panel2(wks_npo_pr,reg_npo_pr_ann,(/nrow,ncol/),panres)
     delete(reg_npo_pr_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NPO PR Regressions (Monthly)"
        gsn_panel2(wks_npo_pr,reg_npo_pr_mon,(/nrow,ncol/),panres)
        delete(reg_npo_pr_mon)
     end if
     delete(wks_npo_pr)  
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
  
  panres2@txString = "PNA (DJF)"
  gsn_panel2(wks_pna_ts,xy_pna_djf,lp,panres2)  
  delete(xy_pna_djf)
  panres2@txString = "PNA (MAM)"
  gsn_panel2(wks_pna_ts,xy_pna_mam,lp,panres2)  
  delete(xy_pna_mam)
  panres2@txString = "PNA (JJA)"
  gsn_panel2(wks_pna_ts,xy_pna_jja,lp,panres2)  
  delete(xy_pna_jja)
  panres2@txString = "PNA (SON)"
  gsn_panel2(wks_pna_ts,xy_pna_son,lp,panres2)  
  delete(xy_pna_son)
  panres2@txString = "PNA (Annual)"
  gsn_panel2(wks_pna_ts,xy_pna_ann,lp,panres2)  
  delete(xy_pna_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "PNA (Monthly)"
     gsn_panel2(wks_pna_ts,xy_pna_mon,lp,panres2)  
     delete(xy_pna_mon) 
  end if
  delete(wks_pna_ts)
  
  panres2@txString = "NPO (DJF)"
  gsn_panel2(wks_npo_ts,xy_npo_djf,lp,panres2)  
  delete(xy_npo_djf)
  panres2@txString = "NPO (MAM)"
  gsn_panel2(wks_npo_ts,xy_npo_mam,lp,panres2)  
  delete(xy_npo_mam)
  panres2@txString = "NPO (JJA)"
  gsn_panel2(wks_npo_ts,xy_npo_jja,lp,panres2)  
  delete(xy_npo_jja)
  panres2@txString = "NPO (SON)"
  gsn_panel2(wks_npo_ts,xy_npo_son,lp,panres2)  
  delete(xy_npo_son)
  panres2@txString = "NPO (Annual)"
  gsn_panel2(wks_npo_ts,xy_npo_ann,lp,panres2)  
  delete(xy_npo_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "NPO (Monthly)"
     gsn_panel2(wks_npo_ts,xy_npo_mon,lp,panres2)  
     delete(xy_npo_mon) 
  end if
  delete(wks_npo_ts)
  
  panres2@txString = "NPI (NDJFM)"
  gsn_panel2(wks_npi_ts,xy_npi,lp,panres2)  
  delete(xy_npi) 
  delete(wks_npi_ts)  
;--------------------------------------------------------------------------------------------------
  OUTDIR = getenv("OUTDIR") 
  if (wks_type.eq."png") then       
     system("mv "+OUTDIR+"pna.000001.png "+OUTDIR+"pna.djf.png") 
     system("mv "+OUTDIR+"pna.000002.png "+OUTDIR+"pna.mam.png") 
     system("mv "+OUTDIR+"pna.000003.png "+OUTDIR+"pna.jja.png") 
     system("mv "+OUTDIR+"pna.000004.png "+OUTDIR+"pna.son.png") 
     system("mv "+OUTDIR+"pna.000005.png "+OUTDIR+"pna.ann.png")
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"pna.000006.png "+OUTDIR+"pna.mon.png")  
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pna.000007.png "+OUTDIR+"pna.tempreg.djf.png") 
           system("mv "+OUTDIR+"pna.000008.png "+OUTDIR+"pna.tempreg.mam.png") 
           system("mv "+OUTDIR+"pna.000009.png "+OUTDIR+"pna.tempreg.jja.png") 
           system("mv "+OUTDIR+"pna.000010.png "+OUTDIR+"pna.tempreg.son.png") 
           system("mv "+OUTDIR+"pna.000011.png "+OUTDIR+"pna.tempreg.ann.png")
           system("mv "+OUTDIR+"pna.000012.png "+OUTDIR+"pna.tempreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"pna.000006.png "+OUTDIR+"pna.tempreg.djf.png") 
           system("mv "+OUTDIR+"pna.000007.png "+OUTDIR+"pna.tempreg.mam.png") 
           system("mv "+OUTDIR+"pna.000008.png "+OUTDIR+"pna.tempreg.jja.png") 
           system("mv "+OUTDIR+"pna.000009.png "+OUTDIR+"pna.tempreg.son.png") 
           system("mv "+OUTDIR+"pna.000010.png "+OUTDIR+"pna.tempreg.ann.png")
        end if
     end if

     if (prreg_frame.eq.0) then
        system("mv "+OUTDIR+"pna.prreg.000001.png "+OUTDIR+"pna.prreg.djf.png") 
        system("mv "+OUTDIR+"pna.prreg.000002.png "+OUTDIR+"pna.prreg.mam.png") 
        system("mv "+OUTDIR+"pna.prreg.000003.png "+OUTDIR+"pna.prreg.jja.png") 
        system("mv "+OUTDIR+"pna.prreg.000004.png "+OUTDIR+"pna.prreg.son.png") 
        system("mv "+OUTDIR+"pna.prreg.000005.png "+OUTDIR+"pna.prreg.ann.png")         
        if (COMPUTE_MODES_MON.eq."True") then     
           system("mv "+OUTDIR+"pna.prreg.000006.png "+OUTDIR+"pna.prreg.mon.png")         
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
           system("mv "+OUTDIR+"npo.000007.png "+OUTDIR+"npo.tempreg.djf.png") 
           system("mv "+OUTDIR+"npo.000008.png "+OUTDIR+"npo.tempreg.mam.png") 
           system("mv "+OUTDIR+"npo.000009.png "+OUTDIR+"npo.tempreg.jja.png") 
           system("mv "+OUTDIR+"npo.000010.png "+OUTDIR+"npo.tempreg.son.png") 
           system("mv "+OUTDIR+"npo.000011.png "+OUTDIR+"npo.tempreg.ann.png")
           system("mv "+OUTDIR+"npo.000012.png "+OUTDIR+"npo.tempreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"npo.000006.png "+OUTDIR+"npo.tempreg.djf.png") 
           system("mv "+OUTDIR+"npo.000007.png "+OUTDIR+"npo.tempreg.mam.png") 
           system("mv "+OUTDIR+"npo.000008.png "+OUTDIR+"npo.tempreg.jja.png") 
           system("mv "+OUTDIR+"npo.000009.png "+OUTDIR+"npo.tempreg.son.png") 
           system("mv "+OUTDIR+"npo.000010.png "+OUTDIR+"npo.tempreg.ann.png")
        end if
     end if

     if (prreg_frame.eq.0) then
        system("mv "+OUTDIR+"npo.prreg.000001.png "+OUTDIR+"npo.prreg.djf.png") 
        system("mv "+OUTDIR+"npo.prreg.000002.png "+OUTDIR+"npo.prreg.mam.png") 
        system("mv "+OUTDIR+"npo.prreg.000003.png "+OUTDIR+"npo.prreg.jja.png") 
        system("mv "+OUTDIR+"npo.prreg.000004.png "+OUTDIR+"npo.prreg.son.png") 
        system("mv "+OUTDIR+"npo.prreg.000005.png "+OUTDIR+"npo.prreg.ann.png")         
        if (COMPUTE_MODES_MON.eq."True") then     
           system("mv "+OUTDIR+"npo.prreg.000006.png "+OUTDIR+"npo.prreg.mon.png")         
        end if
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
  else
     system("psplit "+OUTDIR+"pna.ps "+OUTDIR+"psl_pn")
     system("mv "+OUTDIR+"psl_pn0001.ps "+OUTDIR+"pna.djf.ps") 
     system("mv "+OUTDIR+"psl_pn0002.ps "+OUTDIR+"pna.mam.ps") 
     system("mv "+OUTDIR+"psl_pn0003.ps "+OUTDIR+"pna.jja.ps") 
     system("mv "+OUTDIR+"psl_pn0004.ps "+OUTDIR+"pna.son.ps") 
     system("mv "+OUTDIR+"psl_pn0005.ps "+OUTDIR+"pna.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_pn0006.ps "+OUTDIR+"pna.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_pn0007.ps "+OUTDIR+"pna.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_pn0008.ps "+OUTDIR+"pna.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_pn0009.ps "+OUTDIR+"pna.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_pn0010.ps "+OUTDIR+"pna.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_pn0011.ps "+OUTDIR+"pna.tempreg.ann.ps")
           system("mv "+OUTDIR+"psl_pn0012.ps "+OUTDIR+"pna.tempreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_pn0006.ps "+OUTDIR+"pna.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_pn0007.ps "+OUTDIR+"pna.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_pn0008.ps "+OUTDIR+"pna.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_pn0009.ps "+OUTDIR+"pna.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_pn0010.ps "+OUTDIR+"pna.tempreg.ann.ps")     
        end if
     end if
     
     if (prreg_frame.eq.0) then
        system("psplit "+OUTDIR+"pna.prreg.ps "+OUTDIR+"pr_nn")
        system("mv "+OUTDIR+"pr_nn0001.ps "+OUTDIR+"pna.prreg.djf.ps") 
        system("mv "+OUTDIR+"pr_nn0002.ps "+OUTDIR+"pna.prreg.mam.ps") 
        system("mv "+OUTDIR+"pr_nn0003.ps "+OUTDIR+"pna.prreg.jja.ps") 
        system("mv "+OUTDIR+"pr_nn0004.ps "+OUTDIR+"pna.prreg.son.ps") 
        system("mv "+OUTDIR+"pr_nn0005.ps "+OUTDIR+"pna.prreg.ann.ps") 
        if (COMPUTE_MODES_MON.eq."True") then
           system("mv "+OUTDIR+"pr_nn0006.ps "+OUTDIR+"pna.prreg.mon.ps") 
        end if
     end if

     system("psplit "+OUTDIR+"npo.ps "+OUTDIR+"psl_pn")
     system("mv "+OUTDIR+"psl_pn0001.ps "+OUTDIR+"npo.djf.ps") 
     system("mv "+OUTDIR+"psl_pn0002.ps "+OUTDIR+"npo.mam.ps") 
     system("mv "+OUTDIR+"psl_pn0003.ps "+OUTDIR+"npo.jja.ps") 
     system("mv "+OUTDIR+"psl_pn0004.ps "+OUTDIR+"npo.son.ps") 
     system("mv "+OUTDIR+"psl_pn0005.ps "+OUTDIR+"npo.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_pn0006.ps "+OUTDIR+"npo.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_pn0007.ps "+OUTDIR+"npo.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_pn0008.ps "+OUTDIR+"npo.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_pn0009.ps "+OUTDIR+"npo.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_pn0010.ps "+OUTDIR+"npo.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_pn0011.ps "+OUTDIR+"npo.tempreg.ann.ps")
           system("mv "+OUTDIR+"psl_pn0012.ps "+OUTDIR+"npo.tempreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_pn0006.ps "+OUTDIR+"npo.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_pn0007.ps "+OUTDIR+"npo.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_pn0008.ps "+OUTDIR+"npo.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_pn0009.ps "+OUTDIR+"npo.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_pn0010.ps "+OUTDIR+"npo.tempreg.ann.ps")     
        end if
     end if
     if (prreg_frame.eq.0) then
        system("psplit "+OUTDIR+"npo.prreg.ps "+OUTDIR+"pr_nn")
        system("mv "+OUTDIR+"pr_nn0001.ps "+OUTDIR+"npo.prreg.djf.ps") 
        system("mv "+OUTDIR+"pr_nn0002.ps "+OUTDIR+"npo.prreg.mam.ps") 
        system("mv "+OUTDIR+"pr_nn0003.ps "+OUTDIR+"npo.prreg.jja.ps") 
        system("mv "+OUTDIR+"pr_nn0004.ps "+OUTDIR+"npo.prreg.son.ps") 
        system("mv "+OUTDIR+"pr_nn0005.ps "+OUTDIR+"npo.prreg.ann.ps") 
        if (COMPUTE_MODES_MON.eq."True") then
           system("mv "+OUTDIR+"pr_nn0006.ps "+OUTDIR+"npo.prreg.mon.ps") 
        end if
        system("rm "+OUTDIR+"pna.prreg.ps "+OUTDIR+"npo.prreg.ps")
     end if

     system("psplit "+OUTDIR+"pna.timeseries.ps "+OUTDIR+"psl_pn")
     system("mv "+OUTDIR+"psl_pn0001.ps "+OUTDIR+"pna.timeseries.djf.ps") 
     system("mv "+OUTDIR+"psl_pn0002.ps "+OUTDIR+"pna.timeseries.mam.ps") 
     system("mv "+OUTDIR+"psl_pn0003.ps "+OUTDIR+"pna.timeseries.jja.ps") 
     system("mv "+OUTDIR+"psl_pn0004.ps "+OUTDIR+"pna.timeseries.son.ps") 
     system("mv "+OUTDIR+"psl_pn0005.ps "+OUTDIR+"pna.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_pn0006.ps "+OUTDIR+"pna.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"npo.timeseries.ps "+OUTDIR+"psl_pn")
     system("mv "+OUTDIR+"psl_pn0001.ps "+OUTDIR+"npo.timeseries.djf.ps") 
     system("mv "+OUTDIR+"psl_pn0002.ps "+OUTDIR+"npo.timeseries.mam.ps") 
     system("mv "+OUTDIR+"psl_pn0003.ps "+OUTDIR+"npo.timeseries.jja.ps") 
     system("mv "+OUTDIR+"psl_pn0004.ps "+OUTDIR+"npo.timeseries.son.ps") 
     system("mv "+OUTDIR+"psl_pn0005.ps "+OUTDIR+"npo.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_pn0006.ps "+OUTDIR+"npo.timeseries.mon.ps") 
     end if
     system("rm "+OUTDIR+"npo.timeseries.ps "+OUTDIR+"pna.timeseries.ps "+OUTDIR+"npo.ps "+OUTDIR+"pna.ps")         
  end if
  print("Finished: psl.pna_npo.ncl")
end  
