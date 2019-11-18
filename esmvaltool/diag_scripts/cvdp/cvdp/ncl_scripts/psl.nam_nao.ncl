; Calculates NAM and NAO (patterns and PC timeseries), as well as 
; regressions of those PC timeseries  onto ts, tas, and pr.
;
; Variables used: psl, ts, tas, and pr
; 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: psl.nam_nao.ncl")
  
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
  wks_nam = gsn_open_wks(wks_type,getenv("OUTDIR")+"nam")
  wks_nam_pr = gsn_open_wks(wks_type,getenv("OUTDIR")+"nam.prreg")
  wks_nam_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"nam.timeseries")

  wks_nao = gsn_open_wks(wks_type,getenv("OUTDIR")+"nao")
  wks_nao_pr = gsn_open_wks(wks_type,getenv("OUTDIR")+"nao.prreg")
  wks_nao_ts = gsn_open_wks(wks_type,getenv("OUTDIR")+"nao.timeseries")
  
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_nam,"ncl_default")   
     gsn_define_colormap(wks_nam_ts,"ncl_default")
     gsn_define_colormap(wks_nao,"ncl_default")   
     gsn_define_colormap(wks_nao_ts,"ncl_default")
     gsn_define_colormap(wks_nam_pr,"MPL_BrBG")  
     gsn_define_colormap(wks_nao_pr,"MPL_BrBG")  
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_nam,"BlueDarkRed18")   
     gsn_define_colormap(wks_nam_ts,"ncl_default")
     gsn_define_colormap(wks_nao,"BlueDarkRed18")   
     gsn_define_colormap(wks_nao_ts,"ncl_default")
     gsn_define_colormap(wks_nam_pr,"BrownBlue12")    
     gsn_define_colormap(wks_nao_pr,"BrownBlue12")     
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
  reg_nam_pr_djf = new(nsim,"graphic")  
  reg_nam_pr_mam = new(nsim,"graphic")  
  reg_nam_pr_jja = new(nsim,"graphic")  
  reg_nam_pr_son = new(nsim,"graphic")  
  reg_nam_pr_ann = new(nsim,"graphic")  
  reg_nam_pr_mon = new(nsim,"graphic")  

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
  reg_nao_pr_djf = new(nsim,"graphic")  
  reg_nao_pr_mam = new(nsim,"graphic")  
  reg_nao_pr_jja = new(nsim,"graphic")  
  reg_nao_pr_son = new(nsim,"graphic")  
  reg_nao_pr_ann = new(nsim,"graphic")  
  reg_nao_pr_mon = new(nsim,"graphic")   
  
  xy_npi = new(nsim,"graphic")
  sstreg_frame = 1    ; *reg_frame = flag to create regressions .ps/.png files. Created/used instead of *reg_plot_flag
                      ; so that if {sst,tas,pr} regressions are not created for the last simulation listed that .ps/png files are created
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
        if (syear(ee).eq.syear_ts(ee)) then     ; check that the start and end years match for ts, tas, and psl
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
           delete([/sstV/])
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
        if (syear(ee).eq.syear_pr(ee)) then     ; check that the start and end years match for pr and psl
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
        
        if (prreg_plot_flag.eq.0) then     ; only continue if both PSL/PR fields are present
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

;-------------------------------------------------------------------     
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
;---------NAM calculations----------------------------------------------------------     
     evecv = eofunc(arr_djf_CW({lat|20:},lon|:,time|:),2,75)
     pcts = eofunc_ts(arr_djf_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_djf = dim_standardize(pcts(0,:),0)
     nam_djf = arr_djf(0,:,:)
     nam_djf = (/ regCoef(nam_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_djf = sst_djf(0,:,:)
        nam_sst_djf = (/ regCoef(nam_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
     end if
     if (tasreg_plot_flag.eq.0) then
        nam_tas_djf = tas_djf(0,:,:)
        nam_tas_djf = (/ regCoef(nam_pc_djf,tas_djf(lat|:,lon|:,time|:)) /)
     end if
     if (prreg_plot_flag.eq.0) then
        nam_pr_djf = pr_djf(0,:,:)
        nam_pr_djf = (/ regCoef(nam_pc_djf,pr_djf(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_djf({85},{5}))) then
        if (nam_djf({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_djf = nam_djf*-1.
           nam_pc_djf = nam_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_djf = nam_sst_djf*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nam_tas_djf = nam_tas_djf*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nam_pr_djf = nam_pr_djf*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nam_pc_djf),False)
     if (sig_pcv(0)) then    ; if True then significant
        nam_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nam_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     nam_pc_djf!0 = "TIME"
     nam_pc_djf&TIME = ispan(syear(ee),eyear(ee),1)
     nam_pc_djf&TIME@units = "YYYY"
     nam_pc_djf&TIME@long_name = "time"     
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_mam_CW({lat|20:},lon|:,time|:),2,75)
     pcts = eofunc_ts(arr_mam_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_mam = dim_standardize(pcts(0,:),0)
     nam_mam = arr_mam(0,:,:)
     nam_mam = (/ regCoef(nam_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_mam = sst_mam(0,:,:)
        nam_sst_mam = (/ regCoef(nam_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
     end if
     if (tasreg_plot_flag.eq.0) then
        nam_tas_mam = tas_mam(0,:,:)
        nam_tas_mam = (/ regCoef(nam_pc_mam,tas_mam(lat|:,lon|:,time|:)) /)
     end if
     if (prreg_plot_flag.eq.0) then
        nam_pr_mam = pr_mam(0,:,:)
        nam_pr_mam = (/ regCoef(nam_pc_mam,pr_mam(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_mam({85},{5}))) then
        if (nam_mam({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_mam = nam_mam*-1.
           nam_pc_mam = nam_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_mam = nam_sst_mam*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nam_tas_mam = nam_tas_mam*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nam_pr_mam = nam_pr_mam*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nam_pc_mam),False)
     if (sig_pcv(0)) then    ; if True then significant
        nam_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nam_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)    
     copy_VarCoords(nam_pc_djf,nam_pc_mam)
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_jja_CW({lat|20:},lon|:,time|:),2,75)
     pcts = eofunc_ts(arr_jja_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_jja = dim_standardize(pcts(0,:),0)
     nam_jja = arr_jja(0,:,:)
     nam_jja = (/ regCoef(nam_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_jja = sst_jja(0,:,:)
        nam_sst_jja = (/ regCoef(nam_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
     end if
     if (tasreg_plot_flag.eq.0) then
        nam_tas_jja = tas_jja(0,:,:)
        nam_tas_jja = (/ regCoef(nam_pc_jja,tas_jja(lat|:,lon|:,time|:)) /)
     end if
     if (prreg_plot_flag.eq.0) then
        nam_pr_jja = pr_jja(0,:,:)
        nam_pr_jja = (/ regCoef(nam_pc_jja,pr_jja(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_jja({85},{5}))) then
        if (nam_jja({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_jja = nam_jja*-1.
           nam_pc_jja = nam_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_jja = nam_sst_jja*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nam_tas_jja = nam_tas_jja*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nam_pr_jja = nam_pr_jja*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nam_pc_jja),False)
     if (sig_pcv(0)) then    ; if True then significant
        nam_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nam_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nam_pc_jja)  
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_son_CW({lat|20:},lon|:,time|:),2,75)
     pcts = eofunc_ts(arr_son_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_son = dim_standardize(pcts(0,:),0)
     nam_son = arr_son(0,:,:)
     nam_son = (/ regCoef(nam_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_son = sst_son(0,:,:)
        nam_sst_son = (/ regCoef(nam_pc_son,sst_son(lat|:,lon|:,time|:)) /)
     end if
     if (tasreg_plot_flag.eq.0) then
        nam_tas_son = tas_son(0,:,:)
        nam_tas_son = (/ regCoef(nam_pc_son,tas_son(lat|:,lon|:,time|:)) /)
     end if
     if (prreg_plot_flag.eq.0) then
        nam_pr_son = pr_son(0,:,:)
        nam_pr_son = (/ regCoef(nam_pc_son,pr_son(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_son({85},{5}))) then
        if (nam_son({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_son = nam_son*-1.
           nam_pc_son = nam_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_son = nam_sst_son*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nam_tas_son = nam_tas_son*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nam_pr_son = nam_pr_son*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nam_pc_son),False)
     if (sig_pcv(0)) then    ; if True then significant
        nam_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nam_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nam_pc_son)    
     delete([/evecv,pcts/])
     
     evecv = eofunc(arr_ann_CW({lat|20:},lon|:,time|:),2,75)
     pcts = eofunc_ts(arr_ann_CW({lat|20:},lon|:,time|:),evecv,False)  
     nam_pc_ann = dim_standardize(pcts(0,:),0)
     nam_ann = arr_ann(0,:,:)
     nam_ann = (/ regCoef(nam_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nam_sst_ann = sst_ann(0,:,:)
        nam_sst_ann = (/ regCoef(nam_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
     end if
     if (tasreg_plot_flag.eq.0) then
        nam_tas_ann = tas_ann(0,:,:)
        nam_tas_ann = (/ regCoef(nam_pc_ann,tas_ann(lat|:,lon|:,time|:)) /)
     end if
     if (prreg_plot_flag.eq.0) then
        nam_pr_ann = pr_ann(0,:,:)
        nam_pr_ann = (/ regCoef(nam_pc_ann,pr_ann(lat|:,lon|:,time|:)) /)
     end if
     if (.not.ismissing(nam_ann({85},{5}))) then
        if (nam_ann({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nam_ann = nam_ann*-1.
           nam_pc_ann = nam_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              nam_sst_ann = nam_sst_ann*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nam_tas_ann = nam_tas_ann*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nam_pr_ann = nam_pr_ann*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nam_pc_ann),False)
     if (sig_pcv(0)) then    ; if True then significant
        nam_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nam_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nam_pc_ann) 
     delete([/evecv,pcts/])    
     
     if (COMPUTE_MODES_MON.eq."True") then
        evecv = eofunc(arr_mon_CW({lat|20:},lon|:,time|:),2,75)
        pcts = eofunc_ts(arr_mon_CW({lat|20:},lon|:,time|:),evecv,False)  
        nam_pc_mon = dim_standardize(pcts(0,:),0)
        nam_mon = arr(0,:,:)
        nam_mon = (/ regCoef(nam_pc_mon,arr(lat|:,lon|:,time|:)) /)
        if (sstreg_plot_flag.eq.0) then
           nam_sst_mon = sst(0,:,:)
           nam_sst_mon = (/ regCoef(nam_pc_mon,sst(lat|:,lon|:,time|:)) /)
        end if
        if (tasreg_plot_flag.eq.0) then
           nam_tas_mon = tas(0,:,:)
           nam_tas_mon = (/ regCoef(nam_pc_mon,tas(lat|:,lon|:,time|:)) /)
        end if
        if (prreg_plot_flag.eq.0) then
           nam_pr_mon = pr(0,:,:)
           nam_pr_mon = (/ regCoef(nam_pc_mon,pr(lat|:,lon|:,time|:)) /)
        end if
        if (.not.ismissing(nam_mon({85},{5}))) then
           if (nam_mon({85},{5}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              nam_mon = nam_mon*-1.
              nam_pc_mon = nam_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 nam_sst_mon = nam_sst_mon*-1.
              end if
              if (tasreg_plot_flag.eq.0) then
                 nam_tas_mon = nam_tas_mon*-1.
              end if
              if (prreg_plot_flag.eq.0) then
                 nam_pr_mon = nam_pr_mon*-1.
              end if
           end if
        end if
        sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nam_pc_mon),False)
        if (sig_pcv(0)) then    ; if True then significant
           nam_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
        else
           nam_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
        end if
        delete(sig_pcv)
        nam_pc_mon!0 = "time"
        nam_pc_mon&time = arr&time   
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
         
     evecv = eofunc(arr_djf_CW_LF({lat|20:80},{lon|-90.:40},time|:),2,75)
     pcts = eofunc_ts(arr_djf_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_djf = dim_standardize(pcts(0,:),0)
     nao_djf = arr_djf(0,:,:)
     nao_djf = (/ regCoef(nao_pc_djf,arr_djf(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_djf = sst_djf(0,:,:)
        nao_sst_djf = (/ regCoef(nao_pc_djf,sst_djf(lat|:,lon|:,time|:)) /)
        delete(sst_djf)
     end if
     if (tasreg_plot_flag.eq.0) then
        nao_tas_djf = tas_djf(0,:,:)
        nao_tas_djf = (/ regCoef(nao_pc_djf,tas_djf(lat|:,lon|:,time|:)) /)
        delete(tas_djf)
     end if
     if (prreg_plot_flag.eq.0) then
        nao_pr_djf = pr_djf(0,:,:)
        nao_pr_djf = (/ regCoef(nao_pc_djf,pr_djf(lat|:,lon|:,time|:)) /)
        delete(pr_djf)
     end if
     if (.not.ismissing(nao_djf({70},{350}))) then
        if (nao_djf({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_djf = nao_djf*-1.
           nao_pc_djf = nao_pc_djf*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_djf = nao_sst_djf*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nao_tas_djf = nao_tas_djf*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nao_pr_djf = nao_pr_djf*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nao_pc_djf),False)
     if (sig_pcv(0)) then    ; if True then significant
        nao_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nao_djf@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nao_pc_djf) 
     delete([/evecv,pcts,arr_djf,arr_djf_CW_LF/])
     
     evecv = eofunc(arr_mam_CW_LF({lat|20:80},{lon|-90.:40},time|:),2,75)
     pcts = eofunc_ts(arr_mam_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_mam = dim_standardize(pcts(0,:),0)
     nao_mam = arr_mam(0,:,:)
     nao_mam = (/ regCoef(nao_pc_mam,arr_mam(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_mam = sst_mam(0,:,:)
        nao_sst_mam = (/ regCoef(nao_pc_mam,sst_mam(lat|:,lon|:,time|:)) /)
        delete(sst_mam)
     end if
     if (tasreg_plot_flag.eq.0) then
        nao_tas_mam = tas_mam(0,:,:)
        nao_tas_mam = (/ regCoef(nao_pc_mam,tas_mam(lat|:,lon|:,time|:)) /)
        delete(tas_mam)
     end if
     if (prreg_plot_flag.eq.0) then
        nao_pr_mam = pr_mam(0,:,:)
        nao_pr_mam = (/ regCoef(nao_pc_mam,pr_mam(lat|:,lon|:,time|:)) /)
        delete(pr_mam)
     end if
     if (.not.ismissing(nao_mam({70},{350}))) then
        if (nao_mam({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_mam = nao_mam*-1.
           nao_pc_mam = nao_pc_mam*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_mam = nao_sst_mam*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nao_tas_mam = nao_tas_mam*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nao_pr_mam = nao_pr_mam*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nao_pc_mam),False)
     if (sig_pcv(0)) then    ; if True then significant
        nao_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nao_mam@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nao_pc_mam) 
     delete([/evecv,pcts,arr_mam,arr_mam_CW_LF/])
     
     evecv = eofunc(arr_jja_CW_LF({lat|20:80},{lon|-90.:40},time|:),2,75)
     pcts = eofunc_ts(arr_jja_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_jja = dim_standardize(pcts(0,:),0)
     nao_jja = arr_jja(0,:,:)
     nao_jja = (/ regCoef(nao_pc_jja,arr_jja(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_jja = sst_jja(0,:,:)
        nao_sst_jja = (/ regCoef(nao_pc_jja,sst_jja(lat|:,lon|:,time|:)) /)
        delete(sst_jja)
     end if
     if (tasreg_plot_flag.eq.0) then
        nao_tas_jja = tas_jja(0,:,:)
        nao_tas_jja = (/ regCoef(nao_pc_jja,tas_jja(lat|:,lon|:,time|:)) /)
        delete(tas_jja)
     end if
     if (prreg_plot_flag.eq.0) then
        nao_pr_jja = pr_jja(0,:,:)
        nao_pr_jja = (/ regCoef(nao_pc_jja,pr_jja(lat|:,lon|:,time|:)) /)
        delete(pr_jja)
     end if
     if (.not.ismissing(nao_jja({70},{350}))) then
        if (nao_jja({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_jja = nao_jja*-1.
           nao_pc_jja = nao_pc_jja*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_jja = nao_sst_jja*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nao_tas_jja = nao_tas_jja*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nao_pr_jja = nao_pr_jja*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nao_pc_jja),False)
     if (sig_pcv(0)) then    ; if True then significant
        nao_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nao_jja@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nao_pc_jja) 
     delete([/evecv,pcts,arr_jja,arr_jja_CW_LF/])
     
     evecv = eofunc(arr_son_CW_LF({lat|20:80},{lon|-90.:40},time|:),2,75)
     pcts = eofunc_ts(arr_son_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_son = dim_standardize(pcts(0,:),0)
     nao_son = arr_son(0,:,:)
     nao_son = (/ regCoef(nao_pc_son,arr_son(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_son = sst_son(0,:,:)
        nao_sst_son = (/ regCoef(nao_pc_son,sst_son(lat|:,lon|:,time|:)) /)
        delete(sst_son)
     end if
     if (tasreg_plot_flag.eq.0) then
        nao_tas_son = tas_son(0,:,:)
        nao_tas_son = (/ regCoef(nao_pc_son,tas_son(lat|:,lon|:,time|:)) /)
        delete(tas_son)
     end if
     if (prreg_plot_flag.eq.0) then
        nao_pr_son = pr_son(0,:,:)
        nao_pr_son = (/ regCoef(nao_pc_son,pr_son(lat|:,lon|:,time|:)) /)
        delete(pr_son)
     end if
     if (.not.ismissing(nao_son({70},{350}))) then
        if (nao_son({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_son = nao_son*-1.
           nao_pc_son = nao_pc_son*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_son = nao_sst_son*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nao_tas_son = nao_tas_son*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nao_pr_son = nao_pr_son*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nao_pc_son),False)
     if (sig_pcv(0)) then    ; if True then significant
        nao_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nao_son@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nao_pc_son) 
     delete([/evecv,pcts,arr_son,arr_son_CW_LF/])
     
     evecv = eofunc(arr_ann_CW_LF({lat|20:80},{lon|-90.:40},time|:),2,75)
     pcts = eofunc_ts(arr_ann_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
     nao_pc_ann = dim_standardize(pcts(0,:),0)
     nao_ann = arr_ann(0,:,:)
     nao_ann = (/ regCoef(nao_pc_ann,arr_ann(lat|:,lon|:,time|:)) /)
     if (sstreg_plot_flag.eq.0) then
        nao_sst_ann = sst_ann(0,:,:)
        nao_sst_ann = (/ regCoef(nao_pc_ann,sst_ann(lat|:,lon|:,time|:)) /)
        delete(sst_ann)
     end if
     if (tasreg_plot_flag.eq.0) then
        nao_tas_ann = tas_ann(0,:,:)
        nao_tas_ann = (/ regCoef(nao_pc_ann,tas_ann(lat|:,lon|:,time|:)) /)
        delete(tas_ann)
     end if
     if (prreg_plot_flag.eq.0) then
        nao_pr_ann = pr_ann(0,:,:)
        nao_pr_ann = (/ regCoef(nao_pc_ann,pr_ann(lat|:,lon|:,time|:)) /)
        delete(pr_ann)
     end if
     if (.not.ismissing(nao_ann({70},{350}))) then
        if (nao_ann({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
           nao_ann = nao_ann*-1.
           nao_pc_ann = nao_pc_ann*-1.
           if (sstreg_plot_flag.eq.0) then
              nao_sst_ann = nao_sst_ann*-1.
           end if
           if (tasreg_plot_flag.eq.0) then
              nao_tas_ann = nao_tas_ann*-1.
           end if
           if (prreg_plot_flag.eq.0) then
              nao_pr_ann = nao_pr_ann*-1.
           end if
        end if
     end if
     sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nao_pc_ann),False)
     if (sig_pcv(0)) then    ; if True then significant
        nao_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
     else
        nao_ann@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
     end if
     delete(sig_pcv)
     copy_VarCoords(nam_pc_djf,nao_pc_ann) 
     delete([/evecv,pcts,arr_ann,arr_ann_CW_LF/])

     if (COMPUTE_MODES_MON.eq."True") then               
        evecv = eofunc(arr_mon_CW_LF({lat|20:80},{lon|-90.:40},time|:),2,75)
        pcts = eofunc_ts(arr_mon_CW_LF({lat|20:80},{lon|-90.:40},time|:),evecv,False)
        nao_pc_mon = dim_standardize(pcts(0,:),0)
        nao_mon = arr(0,:,:)
        nao_mon = (/ regCoef(nao_pc_mon,arr(lat|:,lon|:,time|:)) /)
        if (sstreg_plot_flag.eq.0) then
           nao_sst_mon = sst(0,:,:)
           nao_sst_mon = (/ regCoef(nao_pc_mon,sst(lat|:,lon|:,time|:)) /)
           delete(sst)
        end if
        if (tasreg_plot_flag.eq.0) then
           nao_tas_mon = tas(0,:,:)
           nao_tas_mon = (/ regCoef(nao_pc_mon,tas(lat|:,lon|:,time|:)) /)
           delete(tas)
        end if
        if (prreg_plot_flag.eq.0) then
           nao_pr_mon = pr(0,:,:)
           nao_pr_mon = (/ regCoef(nao_pc_mon,pr(lat|:,lon|:,time|:)) /)
           delete(pr)
        end if
        if (.not.ismissing(nao_mon({70},{350}))) then
           if (nao_mon({70},{350}).ge.0) then  ; arbitrary attempt to make all plots have the same sign..
              nao_mon = nao_mon*-1.
              nao_pc_mon = nao_pc_mon*-1.
              if (sstreg_plot_flag.eq.0) then
                 nao_sst_mon = nao_sst_mon*-1.
              end if
              if (tasreg_plot_flag.eq.0) then
                 nao_tas_mon = nao_tas_mon*-1.
              end if
              if (prreg_plot_flag.eq.0) then
                 nao_pr_mon = nao_pr_mon*-1.
              end if
           end if
        end if
        sig_pcv = eofunc_north2(evecv@pcvar,dimsizes(nao_pc_mon),False)
        if (sig_pcv(0)) then    ; if True then significant
           nao_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%*"
        else
           nao_mon@pcvar = tofloat(sprintf("%4.1f", evecv@pcvar(0)))+"%"
        end if
        delete(sig_pcv)
        nao_pc_mon!0 = "time"
        nao_pc_mon&time = arr&time   
        delete([/evecv,pcts,arr,arr_mon_CW_LF/])
     end if
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
        fn = getenv("OUTDIR")+modname+".cvdp_data.psl.nam_nao."+syear(ee)+"-"+eyear(ee)+".nc"
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
        
        z->nao_timeseries_djf = set_varAtts(nao_pc_djf,"NAO normalized principal component timeseries (DJF)","1","")
        z->nao_timeseries_mam = set_varAtts(nao_pc_mam,"NAO normalized principal component timeseries (MAM)","1","")
        z->nao_timeseries_jja = set_varAtts(nao_pc_jja,"NAO normalized principal component timeseries (JJA)","1","")
        z->nao_timeseries_son = set_varAtts(nao_pc_son,"NAO normalized principal component timeseries (SON)","1","")
        z->nao_timeseries_ann = set_varAtts(nao_pc_ann,"NAO normalized principal component timeseries (ANN)","1","")
        
        z->nam_timeseries_djf = set_varAtts(nam_pc_djf,"NAM normalized principal component timeseries (DJF)","1","")
        z->nam_timeseries_mam = set_varAtts(nam_pc_mam,"NAM normalized principal component timeseries (MAM)","1","")
        z->nam_timeseries_jja = set_varAtts(nam_pc_jja,"NAM normalized principal component timeseries (JJA)","1","")
        z->nam_timeseries_son = set_varAtts(nam_pc_son,"NAM normalized principal component timeseries (SON)","1","")
        z->nam_timeseries_ann = set_varAtts(nam_pc_ann,"NAM normalized principal component timeseries (ANN)","1","")
        
        z->nao_pattern_djf = set_varAtts(nao_djf,"NAO spatial pattern (DJF)","","")   
        z->nao_pattern_mam = set_varAtts(nao_mam,"NAO spatial pattern (MAM)","","")
        z->nao_pattern_jja = set_varAtts(nao_jja,"NAO spatial pattern (JJA)","","")
        z->nao_pattern_son = set_varAtts(nao_son,"NAO spatial pattern (SON)","","")
        z->nao_pattern_ann = set_varAtts(nao_ann,"NAO spatial pattern (annual)","","")
        
        z->nam_pattern_djf = set_varAtts(nam_djf,"NAM spatial pattern (DJF)","","")
        z->nam_pattern_mam = set_varAtts(nam_mam,"NAM spatial pattern (MAM)","","")
        z->nam_pattern_jja = set_varAtts(nam_jja,"NAM spatial pattern (JJA)","","")
        z->nam_pattern_son = set_varAtts(nam_son,"NAM spatial pattern (SON)","","")
        z->nam_pattern_ann = set_varAtts(nam_ann,"NAM spatial pattern (annual)","","")

        if (COMPUTE_MODES_MON.eq."True") then   
           z->nao_timeseries_mon = set_varAtts(nao_pc_mon,"NAO principal component timeseries (monthly)","","")
           z->nam_timeseries_mon = set_varAtts(nam_pc_mon,"NAM principal component timeseries (monthly)","","")
           z->nao_pattern_mon = set_varAtts(nao_mon,"NAO spatial pattern (monthly)","","")
           z->nam_pattern_mon = set_varAtts(nam_mon,"NAM spatial pattern (monthly)","","")
        end if
        delete(z)
        delete([/modname,fn/])

        if (sstreg_plot_flag.eq.0) then
           modname = str_sub_str(names_ts(ee)," ","_")
           bc = (/"/","'","(",")"/)
           do gg = 0,dimsizes(bc)-1
              modname = str_sub_str(modname,bc(gg),"_")
           end do
           fn = getenv("OUTDIR")+modname+".cvdp_data.psl.nam_nao.ts."+syear_ts(ee)+"-"+eyear_ts(ee)+".nc"
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
           z->nao_sst_regression_djf = set_varAtts(nao_sst_djf,"sst regression onto NAO principal component timeseries (DJF)","","")    
           z->nao_sst_regression_mam = set_varAtts(nao_sst_mam,"sst regression onto NAO principal component timeseries (MAM)","","")    
           z->nao_sst_regression_jja = set_varAtts(nao_sst_jja,"sst regression onto NAO principal component timeseries (JJA)","","")    
           z->nao_sst_regression_son = set_varAtts(nao_sst_son,"sst regression onto NAO principal component timeseries (SON)","","")    
           z->nao_sst_regression_ann = set_varAtts(nao_sst_ann,"sst regression onto NAO principal component timeseries (annual)","","")    
        
           z->nam_sst_regression_djf = set_varAtts(nam_sst_djf,"sst regression onto NAM principal component timeseries (DJF)","","")    
           z->nam_sst_regression_mam = set_varAtts(nam_sst_mam,"sst regression onto NAM principal component timeseries (MAM)","","")    
           z->nam_sst_regression_jja = set_varAtts(nam_sst_jja,"sst regression onto NAM principal component timeseries (JJA)","","")    
           z->nam_sst_regression_son = set_varAtts(nam_sst_son,"sst regression onto NAM principal component timeseries (SON)","","")    
           z->nam_sst_regression_ann = set_varAtts(nam_sst_ann,"sst regression onto NAM principal component timeseries (annual)","","")    
           if (COMPUTE_MODES_MON.eq."True") then   
              z->nao_sst_regression_mon = set_varAtts(nao_sst_mon,"sst regression onto NAO principal component timeseries (monthly)","","")    
              z->nam_sst_regression_mon = set_varAtts(nam_sst_mon,"sst regression onto NAM principal component timeseries (monthly)","","")    
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
           fn = getenv("OUTDIR")+modname+".cvdp_data.psl.nam_nao.tas."+syear_tas(ee)+"-"+eyear_tas(ee)+".nc"
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
           z->nao_tas_regression_djf = set_varAtts(nao_tas_djf,"tas regression onto NAO principal component timeseries (DJF)","","")    
           z->nao_tas_regression_mam = set_varAtts(nao_tas_mam,"tas regression onto NAO principal component timeseries (MAM)","","")    
           z->nao_tas_regression_jja = set_varAtts(nao_tas_jja,"tas regression onto NAO principal component timeseries (JJA)","","")    
           z->nao_tas_regression_son = set_varAtts(nao_tas_son,"tas regression onto NAO principal component timeseries (SON)","","")    
           z->nao_tas_regression_ann = set_varAtts(nao_tas_ann,"tas regression onto NAO principal component timeseries (annual)","","")    
        
           z->nam_tas_regression_djf = set_varAtts(nam_tas_djf,"tas regression onto NAM principal component timeseries (DJF)","","")    
           z->nam_tas_regression_mam = set_varAtts(nam_tas_mam,"tas regression onto NAM principal component timeseries (MAM)","","")    
           z->nam_tas_regression_jja = set_varAtts(nam_tas_jja,"tas regression onto NAM principal component timeseries (JJA)","","")    
           z->nam_tas_regression_son = set_varAtts(nam_tas_son,"tas regression onto NAM principal component timeseries (SON)","","")    
           z->nam_tas_regression_ann = set_varAtts(nam_tas_ann,"tas regression onto NAM principal component timeseries (annual)","","")    
           if (COMPUTE_MODES_MON.eq."True") then   
              z->nao_tas_regression_mon = set_varAtts(nao_tas_mon,"tas regression onto NAO principal component timeseries (monthly)","","")    
              z->nam_tas_regression_mon = set_varAtts(nam_tas_mon,"tas regression onto NAM principal component timeseries (monthly)","","")    
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
           fn = getenv("OUTDIR")+modname+".cvdp_data.psl.nam_nao.pr."+syear_pr(ee)+"-"+eyear_pr(ee)+".nc"
           if (.not.isfilepresent2(fn)) then
              z = addfile(fn,"c")
              z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
              z@notes = "Data from "+names_pr(ee)+" from "+syear_pr(ee)+"-"+eyear_pr(ee)
              if (OPT_CLIMO.eq."Full") then
                 z@climatology = syear_pr(ee)+"-"+eyear_pr(ee)+" climatology removed prior to all calculations (other than means)"
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
           z->nao_pr_regression_djf = set_varAtts(nao_pr_djf,"pr regression onto NAO principal component timeseries (DJF)","","")    
           z->nao_pr_regression_mam = set_varAtts(nao_pr_mam,"pr regression onto NAO principal component timeseries (MAM)","","")    
           z->nao_pr_regression_jja = set_varAtts(nao_pr_jja,"pr regression onto NAO principal component timeseries (JJA)","","")    
           z->nao_pr_regression_son = set_varAtts(nao_pr_son,"pr regression onto NAO principal component timeseries (SON)","","")    
           z->nao_pr_regression_ann = set_varAtts(nao_pr_ann,"pr regression onto NAO principal component timeseries (annual)","","")    
        
           z->nam_pr_regression_djf = set_varAtts(nam_pr_djf,"pr regression onto NAM principal component timeseries (DJF)","","")    
           z->nam_pr_regression_mam = set_varAtts(nam_pr_mam,"pr regression onto NAM principal component timeseries (MAM)","","")    
           z->nam_pr_regression_jja = set_varAtts(nam_pr_jja,"pr regression onto NAM principal component timeseries (JJA)","","")    
           z->nam_pr_regression_son = set_varAtts(nam_pr_son,"pr regression onto NAM principal component timeseries (SON)","","")    
           z->nam_pr_regression_ann = set_varAtts(nam_pr_ann,"pr regression onto NAM principal component timeseries (annual)","","")    
           if (COMPUTE_MODES_MON.eq."True") then   
              z->nao_pr_regression_mon = set_varAtts(nao_pr_mon,"pr regression onto NAO principal component timeseries (monthly)","","")    
              z->nam_pr_regression_mon = set_varAtts(nam_pr_mon,"pr regression onto NAM principal component timeseries (monthly)","","")    
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

     
     if (isfilepresent2("obs_psl").and.ee.eq.0) then    ; for pattern correlation table. Save entire lat/lon array
        patcor_nam_djf = new((/nsim,dimsizes(nam_djf&lat),dimsizes(nam_djf&lon)/),typeof(nam_djf))
        patcor_nam_djf!1 = "lat"
        patcor_nam_djf&lat = nam_djf&lat
        patcor_nam_djf!2 = "lon"
        patcor_nam_djf&lon = nam_djf&lon
        patcor_nam_jja = patcor_nam_djf
        patcor_nam_ann = patcor_nam_djf
        patcor_nao_djf = patcor_nam_djf
        patcor_nao_jja = patcor_nam_djf
        patcor_nao_ann = patcor_nam_djf
        patcor_nam_djf(ee,:,:) = (/ nam_djf /)
        patcor_nam_jja(ee,:,:) = (/ nam_jja /)
        patcor_nam_ann(ee,:,:) = (/ nam_ann /)
        patcor_nao_djf(ee,:,:) = (/ nao_djf /)
        patcor_nao_jja(ee,:,:) = (/ nao_jja /)
        patcor_nao_ann(ee,:,:) = (/ nao_ann /)
     end if
     if (isfilepresent2("obs_psl").and.ee.ge.1.and.isvar("patcor_nam_djf")) then
        patcor_nam_djf(ee,:,:) = (/ totype(linint2(nam_djf&lon,nam_djf&lat,nam_djf,True,patcor_nam_djf&lon,patcor_nam_djf&lat,0),typeof(patcor_nam_djf)) /)
        patcor_nam_jja(ee,:,:) = (/ totype(linint2(nam_jja&lon,nam_jja&lat,nam_jja,True,patcor_nam_jja&lon,patcor_nam_jja&lat,0),typeof(patcor_nam_jja)) /)
        patcor_nam_ann(ee,:,:) = (/ totype(linint2(nam_ann&lon,nam_ann&lat,nam_ann,True,patcor_nam_ann&lon,patcor_nam_ann&lat,0),typeof(patcor_nam_ann)) /)
        
        patcor_nao_djf(ee,:,:) = (/ totype(linint2(nao_djf&lon,nao_djf&lat,nao_djf,True,patcor_nao_djf&lon,patcor_nao_djf&lat,0),typeof(patcor_nao_djf)) /)
        patcor_nao_jja(ee,:,:) = (/ totype(linint2(nao_jja&lon,nao_jja&lat,nao_jja,True,patcor_nao_jja&lon,patcor_nao_jja&lat,0),typeof(patcor_nao_jja)) /)
        patcor_nao_ann(ee,:,:) = (/ totype(linint2(nao_ann&lon,nao_ann&lat,nao_ann,True,patcor_nao_ann&lon,patcor_nao_ann&lat,0),typeof(patcor_nao_ann)) /)
     end if
               
     res@gsnRightString = nam_djf@pcvar
     map_nam_djf(ee) = gsn_csm_contour_map_polar(wks_nam,nam_djf,res)
     res@gsnRightString = nam_mam@pcvar
     map_nam_mam(ee) = gsn_csm_contour_map_polar(wks_nam,nam_mam,res)
     res@gsnRightString = nam_jja@pcvar
     map_nam_jja(ee) = gsn_csm_contour_map_polar(wks_nam,nam_jja,res)
     res@gsnRightString = nam_son@pcvar
     map_nam_son(ee) = gsn_csm_contour_map_polar(wks_nam,nam_son,res)
     res@gsnRightString = nam_ann@pcvar
     map_nam_ann(ee) = gsn_csm_contour_map_polar(wks_nam,nam_ann,res)     
     delete([/nam_djf,nam_mam,nam_jja,nam_son,nam_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = nam_mon@pcvar
        map_nam_mon(ee) = gsn_csm_contour_map_polar(wks_nam,nam_mon,res)
        delete([/nam_mon/])
     end if
     
     res@gsnRightString = nao_djf@pcvar
     map_nao_djf(ee) = gsn_csm_contour_map_polar(wks_nao,nao_djf,res)
     res@gsnRightString = nao_mam@pcvar
     map_nao_mam(ee) = gsn_csm_contour_map_polar(wks_nao,nao_mam,res)
     res@gsnRightString = nao_jja@pcvar
     map_nao_jja(ee) = gsn_csm_contour_map_polar(wks_nao,nao_jja,res)
     res@gsnRightString = nao_son@pcvar
     map_nao_son(ee) = gsn_csm_contour_map_polar(wks_nao,nao_son,res)
     res@gsnRightString = nao_ann@pcvar
     map_nao_ann(ee) = gsn_csm_contour_map_polar(wks_nao,nao_ann,res)     
     delete([/nao_djf,nao_mam,nao_jja,nao_son,nao_ann/])
     if (COMPUTE_MODES_MON.eq."True") then                    
        res@gsnRightString = nao_mon@pcvar
        map_nao_mon(ee) = gsn_csm_contour_map_polar(wks_nao,nao_mon,res)
        delete([/nao_mon/])
     end if
     
     if (sstreg_plot_flag.eq.0) then     ; SSTs have to be present for regressions. TASs are optional
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
        reg_nam_djf(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_djf,res)
        reg_nam_mam(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_mam,res)
        reg_nam_jja(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_jja,res)
        reg_nam_son(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_son,res)
        reg_nam_ann(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_ann,res)     
        delete([/nam_sst_djf,nam_sst_mam,nam_sst_jja,nam_sst_son,nam_sst_ann/])
        if (tasreg_plot_flag.eq.0) then
           o_djf = gsn_csm_contour(wks_nam,nam_tas_djf,res2)
           o_mam = gsn_csm_contour(wks_nam,nam_tas_mam,res2)
           o_jja = gsn_csm_contour(wks_nam,nam_tas_jja,res2)
           o_son = gsn_csm_contour(wks_nam,nam_tas_son,res2)
           o_ann = gsn_csm_contour(wks_nam,nam_tas_ann,res2)
           delete([/nam_tas_djf,nam_tas_mam,nam_tas_jja,nam_tas_son,nam_tas_ann/])
           overlay(reg_nam_djf(ee),o_djf)
           overlay(reg_nam_mam(ee),o_mam)
           overlay(reg_nam_jja(ee),o_jja)
           overlay(reg_nam_son(ee),o_son)
           overlay(reg_nam_ann(ee),o_ann)
           delete([/o_djf,o_mam,o_jja,o_son,o_ann/])           
        end if
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_nam_mon(ee) = gsn_csm_contour_map_polar(wks_nam,nam_sst_mon,res)
           delete([/nam_sst_mon/])
           if (tasreg_plot_flag.eq.0) then
              o_mon = gsn_csm_contour(wks_nam,nam_tas_mon,res2)
              overlay(reg_nam_mon(ee),o_mon)
              delete([/o_mon,nam_tas_mon/])
           end if
        end if
     
        reg_nao_djf(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_djf,res)
        reg_nao_mam(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_mam,res)
        reg_nao_jja(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_jja,res)
        reg_nao_son(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_son,res)
        reg_nao_ann(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_ann,res)     
        delete([/nao_sst_djf,nao_sst_mam,nao_sst_jja,nao_sst_son,nao_sst_ann/])
        if (tasreg_plot_flag.eq.0) then
           o_djf = gsn_csm_contour(wks_nao,nao_tas_djf,res2)
           o_mam = gsn_csm_contour(wks_nao,nao_tas_mam,res2)
           o_jja = gsn_csm_contour(wks_nao,nao_tas_jja,res2)
           o_son = gsn_csm_contour(wks_nao,nao_tas_son,res2)
           o_ann = gsn_csm_contour(wks_nao,nao_tas_ann,res2)
           delete([/nao_tas_djf,nao_tas_mam,nao_tas_jja,nao_tas_son,nao_tas_ann/])
           overlay(reg_nao_djf(ee),o_djf)
           overlay(reg_nao_mam(ee),o_mam)
           overlay(reg_nao_jja(ee),o_jja)
           overlay(reg_nao_son(ee),o_son)
           overlay(reg_nao_ann(ee),o_ann)
           delete([/o_djf,o_mam,o_jja,o_son,o_ann/])           
        end if
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_nao_mon(ee) = gsn_csm_contour_map_polar(wks_nao,nao_sst_mon,res)
           delete([/nao_sst_mon/])
           if (tasreg_plot_flag.eq.0) then
              o_mon = gsn_csm_contour(wks_nao,nao_tas_mon,res2)
              overlay(reg_nao_mon(ee),o_mon)
              delete([/o_mon,nao_tas_mon/])
           end if
        end if
     end if

     if (prreg_plot_flag.eq.0) then     ; PR regressions
        res4@gsnRightString = ""
        res4@gsnCenterString = names_pr(ee)
        reg_nam_pr_djf(ee) = gsn_csm_contour_map_polar(wks_nam_pr,nam_pr_djf,res4)
        reg_nam_pr_mam(ee) = gsn_csm_contour_map_polar(wks_nam_pr,nam_pr_mam,res4)
        reg_nam_pr_jja(ee) = gsn_csm_contour_map_polar(wks_nam_pr,nam_pr_jja,res4)
        reg_nam_pr_son(ee) = gsn_csm_contour_map_polar(wks_nam_pr,nam_pr_son,res4)
        reg_nam_pr_ann(ee) = gsn_csm_contour_map_polar(wks_nam_pr,nam_pr_ann,res4)     
        delete([/nam_pr_djf,nam_pr_mam,nam_pr_jja,nam_pr_son,nam_pr_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_nam_pr_mon(ee) = gsn_csm_contour_map_polar(wks_nam_pr,nam_pr_mon,res4)
           delete([/nam_pr_mon/])
        end if
     
        reg_nao_pr_djf(ee) = gsn_csm_contour_map_polar(wks_nao_pr,nao_pr_djf,res4)
        reg_nao_pr_mam(ee) = gsn_csm_contour_map_polar(wks_nao_pr,nao_pr_mam,res4)
        reg_nao_pr_jja(ee) = gsn_csm_contour_map_polar(wks_nao_pr,nao_pr_jja,res4)
        reg_nao_pr_son(ee) = gsn_csm_contour_map_polar(wks_nao_pr,nao_pr_son,res4)
        reg_nao_pr_ann(ee) = gsn_csm_contour_map_polar(wks_nao_pr,nao_pr_ann,res4)     
        delete([/nao_pr_djf,nao_pr_mam,nao_pr_jja,nao_pr_son,nao_pr_ann/])
        if (COMPUTE_MODES_MON.eq."True") then                    
           reg_nao_pr_mon(ee) = gsn_csm_contour_map_polar(wks_nao_pr,nao_pr_mon,res4)
           delete([/nao_pr_mon/])
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
        xy_nam_mon(ee) = gsn_csm_xy(wks_nam_ts,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nam_pc_mon)),nam_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/nam_pc_mon/])
     end if

     xy_nao_djf(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_djf)),nao_pc_djf,xyres)   ; use standardized timeseries  
     xy_nao_mam(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_mam)),nao_pc_mam,xyres)   ; use standardized timeseries  
     xy_nao_jja(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_jja)),nao_pc_jja,xyres)   ; use standardized timeseries  
     xy_nao_son(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_son)),nao_pc_son,xyres)   ; use standardized timeseries  
     xy_nao_ann(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee),dimsizes(nao_pc_ann)),nao_pc_ann,xyres)   ; use standardized timeseries  
     delete([/nao_pc_djf,nao_pc_mam,nao_pc_jja,nao_pc_son,nao_pc_ann/])
     if (COMPUTE_MODES_MON.eq."True") then     
        xy_nao_mon(ee) = gsn_csm_xy(wks_nao_ts,fspan(syear(ee),eyear(ee)+.91667,dimsizes(nao_pc_mon)),nao_pc_mon,xyresmon)   ; use standardized timeseries  
        delete([/nao_pc_mon/])
     end if
     delete(sstreg_plot_flag)  
  end do
  
  if (isvar("clim_syear")) then
     delete(clim_syear)
  end if
  if (isvar("clim_eyear")) then
     delete(clim_eyear)
  end if
  
  if (isvar("patcor_nam_djf")) then    ; for pattern correlation table  
     clat = cos(0.01745329*patcor_nam_djf&lat)
     clat!0 = "lat"
     clat&lat = patcor_nam_djf&lat
     line3   = "                  "    ; Must be 18 characters long
     line4   = line3
     header = (/"","Pattern Correlations   Observations vs. Model(s)",""/)
     finpr_nam_djf  = "NAM (DJF)         "    ; Must be 18 characters long 
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
           finpr_nam_djf = finpr_nam_djf+sprintf(format2,(pattern_cor(patcor_nam_djf(0,{20:},:),patcor_nam_djf(hh,{20:},:),clat({20:}),0)))+"/"+sprintf(format3,(wgt_arearmse(patcor_nam_djf(0,{20:},:),patcor_nam_djf(hh,{20:},:),clat({20:}),1.0,0)))
        end if  
     end do
;     
     if (dimsizes(tochar(line4)).ge.8190) then   ; system or fortran compiler limit
        print("Metrics table warning: Not creating metrics table as size of comparison results in a invalid ascii row size.")   
     else           
        write_table(getenv("OUTDIR")+"metrics.psl.nam_nao.txt","w",[/header/],"%s")
        write_table(getenv("OUTDIR")+"metrics.psl.nam_nao.txt","a",[/line3/],"%s")
        write_table(getenv("OUTDIR")+"metrics.psl.nam_nao.txt","a",[/line4/],"%s")  
        write_table(getenv("OUTDIR")+"metrics.psl.nam_nao.txt","a",[/finpr_nam_djf/],"%s")
     end if
     delete([/line3,line4,format2,format3,nchar,ntc,clat,patcor_nam_djf,patcor_nam_jja,patcor_nam_ann/])
     delete([/patcor_nao_djf,patcor_nao_jja,patcor_nao_ann/])
     delete([/dimY,ntb,header/])
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
  
  panres@txString = "NAM (DJF)" 
  gsn_panel2(wks_nam,map_nam_djf,(/nrow,ncol/),panres)
  delete(map_nam_djf)
  panres@txString = "NAM (MAM)"
  gsn_panel2(wks_nam,map_nam_mam,(/nrow,ncol/),panres)
  delete(map_nam_mam)
  panres@txString = "NAM (JJA)"
  gsn_panel2(wks_nam,map_nam_jja,(/nrow,ncol/),panres)
  delete(map_nam_jja)
  panres@txString = "NAM (SON)"
  gsn_panel2(wks_nam,map_nam_son,(/nrow,ncol/),panres)
  delete(map_nam_son)
  panres@txString = "NAM (Annual)"
  gsn_panel2(wks_nam,map_nam_ann,(/nrow,ncol/),panres)
  delete(map_nam_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "NAM (Monthly)"
     gsn_panel2(wks_nam,map_nam_mon,(/nrow,ncol/),panres)
     delete(map_nam_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     if (tasreg_frame.eq.0) then
        txt0 = "SST/TAS"
     else
        txt0 = "SST"
     end if
     panres@txString = "NAM "+txt0+" Regressions (DJF)" 
     gsn_panel2(wks_nam,reg_nam_djf,(/nrow,ncol/),panres)
     delete(reg_nam_djf)
     panres@txString = "NAM "+txt0+" Regressions (MAM)"
     gsn_panel2(wks_nam,reg_nam_mam,(/nrow,ncol/),panres)
     delete(reg_nam_mam)
     panres@txString = "NAM "+txt0+" Regressions (JJA)"
     gsn_panel2(wks_nam,reg_nam_jja,(/nrow,ncol/),panres)
     delete(reg_nam_jja)
     panres@txString = "NAM "+txt0+" Regressions (SON)"
     gsn_panel2(wks_nam,reg_nam_son,(/nrow,ncol/),panres)
     delete(reg_nam_son)
     panres@txString = "NAM "+txt0+" Regressions (Annual)"
     gsn_panel2(wks_nam,reg_nam_ann,(/nrow,ncol/),panres)
     delete(reg_nam_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NAM "+txt0+" Regressions (Monthly)"
        gsn_panel2(wks_nam,reg_nam_mon,(/nrow,ncol/),panres)
        delete(reg_nam_mon)
     end if
     delete(wks_nam)  
  end if
  if (prreg_frame.eq.0) then
     panres@txString = "NAM PR Regressions (DJF)" 
     gsn_panel2(wks_nam_pr,reg_nam_pr_djf,(/nrow,ncol/),panres)
     delete(reg_nam_pr_djf)
     panres@txString = "NAM PR Regressions (MAM)"
     gsn_panel2(wks_nam_pr,reg_nam_pr_mam,(/nrow,ncol/),panres)
     delete(reg_nam_pr_mam)
     panres@txString = "NAM PR Regressions (JJA)"
     gsn_panel2(wks_nam_pr,reg_nam_pr_jja,(/nrow,ncol/),panres)
     delete(reg_nam_pr_jja)
     panres@txString = "NAM PR Regressions (SON)"
     gsn_panel2(wks_nam_pr,reg_nam_pr_son,(/nrow,ncol/),panres)
     delete(reg_nam_pr_son)
     panres@txString = "NAM PR Regressions (Annual)"
     gsn_panel2(wks_nam_pr,reg_nam_pr_ann,(/nrow,ncol/),panres)
     delete(reg_nam_pr_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NAM PR Regressions (Monthly)"
        gsn_panel2(wks_nam_pr,reg_nam_pr_mon,(/nrow,ncol/),panres)
        delete(reg_nam_pr_mon)
     end if
     delete(wks_nam_pr)  
  end if
  
  panres@txString = "NAO (DJF)" 
  gsn_panel2(wks_nao,map_nao_djf,(/nrow,ncol/),panres)
  delete(map_nao_djf)
  panres@txString = "NAO (MAM)"
  gsn_panel2(wks_nao,map_nao_mam,(/nrow,ncol/),panres)
  delete(map_nao_mam)
  panres@txString = "NAO (JJA)"
  gsn_panel2(wks_nao,map_nao_jja,(/nrow,ncol/),panres)
  delete(map_nao_jja)
  panres@txString = "NAO (SON)"
  gsn_panel2(wks_nao,map_nao_son,(/nrow,ncol/),panres)
  delete(map_nao_son)
  panres@txString = "NAO (Annual)"
  gsn_panel2(wks_nao,map_nao_ann,(/nrow,ncol/),panres)
  delete(map_nao_ann)
  if (COMPUTE_MODES_MON.eq."True") then    
     panres@txString = "NAO (Monthly)"
     gsn_panel2(wks_nao,map_nao_mon,(/nrow,ncol/),panres)
     delete(map_nao_mon)
  end if
  
  if (sstreg_frame.eq.0) then
     if (tasreg_frame.eq.0) then
        txt0 = "SST/TAS"
     else
        txt0 = "SST"
     end if
     panres@txString = "NAO "+txt0+" Regressions (DJF)" 
     gsn_panel2(wks_nao,reg_nao_djf,(/nrow,ncol/),panres)
     delete(reg_nao_djf)
     panres@txString = "NAO "+txt0+" Regressions (MAM)"
     gsn_panel2(wks_nao,reg_nao_mam,(/nrow,ncol/),panres)
     delete(reg_nao_mam)
     panres@txString = "NAO "+txt0+" Regressions (JJA)"
     gsn_panel2(wks_nao,reg_nao_jja,(/nrow,ncol/),panres)
     delete(reg_nao_jja)
     panres@txString = "NAO "+txt0+" Regressions (SON)"
     gsn_panel2(wks_nao,reg_nao_son,(/nrow,ncol/),panres)
     delete(reg_nao_son)
     panres@txString = "NAO "+txt0+" Regressions (Annual)"
     gsn_panel2(wks_nao,reg_nao_ann,(/nrow,ncol/),panres)
     delete(reg_nao_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NAO "+txt0+" Regressions (Monthly)"
        gsn_panel2(wks_nao,reg_nao_mon,(/nrow,ncol/),panres)
        delete(reg_nao_mon)
     end if
     delete(wks_nao)  
  end if
  if (prreg_frame.eq.0) then
     panres@txString = "NAO PR Regressions (DJF)" 
     gsn_panel2(wks_nao_pr,reg_nao_pr_djf,(/nrow,ncol/),panres)
     delete(reg_nao_pr_djf)
     panres@txString = "NAO PR Regressions (MAM)"
     gsn_panel2(wks_nao_pr,reg_nao_pr_mam,(/nrow,ncol/),panres)
     delete(reg_nao_pr_mam)
     panres@txString = "NAO PR Regressions (JJA)"
     gsn_panel2(wks_nao_pr,reg_nao_pr_jja,(/nrow,ncol/),panres)
     delete(reg_nao_pr_jja)
     panres@txString = "NAO PR Regressions (SON)"
     gsn_panel2(wks_nao_pr,reg_nao_pr_son,(/nrow,ncol/),panres)
     delete(reg_nao_pr_son)
     panres@txString = "NAO PR Regressions (Annual)"
     gsn_panel2(wks_nao_pr,reg_nao_pr_ann,(/nrow,ncol/),panres)
     delete(reg_nao_pr_ann)
     if (COMPUTE_MODES_MON.eq."True") then    
        panres@txString = "NAO PR Regressions (Monthly)"
        gsn_panel2(wks_nao_pr,reg_nao_pr_mon,(/nrow,ncol/),panres)
        delete(reg_nao_pr_mon)
     end if
     delete(wks_nao_pr)  
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
  gsn_panel2(wks_nam_ts,xy_nam_djf,lp,panres2)  
  delete(xy_nam_djf)
  panres2@txString = "NAM (MAM)"
  gsn_panel2(wks_nam_ts,xy_nam_mam,lp,panres2)  
  delete(xy_nam_mam)
  panres2@txString = "NAM (JJA)"
  gsn_panel2(wks_nam_ts,xy_nam_jja,lp,panres2)  
  delete(xy_nam_jja)
  panres2@txString = "NAM (SON)"
  gsn_panel2(wks_nam_ts,xy_nam_son,lp,panres2)  
  delete(xy_nam_son)
  panres2@txString = "NAM (Annual)"
  gsn_panel2(wks_nam_ts,xy_nam_ann,lp,panres2)  
  delete(xy_nam_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "NAM (Monthly)"
     gsn_panel2(wks_nam_ts,xy_nam_mon,lp,panres2)  
     delete(xy_nam_mon) 
  end if
  delete(wks_nam_ts)
    
  panres2@txString = "NAO (DJF)"
  gsn_panel2(wks_nao_ts,xy_nao_djf,lp,panres2)  
  delete(xy_nao_djf)
  panres2@txString = "NAO (MAM)"
  gsn_panel2(wks_nao_ts,xy_nao_mam,lp,panres2)  
  delete(xy_nao_mam)
  panres2@txString = "NAO (JJA)"
  gsn_panel2(wks_nao_ts,xy_nao_jja,lp,panres2)  
  delete(xy_nao_jja)
  panres2@txString = "NAO (SON)"
  gsn_panel2(wks_nao_ts,xy_nao_son,lp,panres2)  
  delete(xy_nao_son)
  panres2@txString = "NAO (Annual)"
  gsn_panel2(wks_nao_ts,xy_nao_ann,lp,panres2)  
  delete(xy_nao_ann)
  if (COMPUTE_MODES_MON.eq."True") then  
     panres2@txString = "NAO (Monthly)"
     gsn_panel2(wks_nao_ts,xy_nao_mon,lp,panres2)  
     delete(xy_nao_mon) 
  end if
  delete(wks_nao_ts)
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
           system("mv "+OUTDIR+"nam.000007.png "+OUTDIR+"nam.tempreg.djf.png") 
           system("mv "+OUTDIR+"nam.000008.png "+OUTDIR+"nam.tempreg.mam.png") 
           system("mv "+OUTDIR+"nam.000009.png "+OUTDIR+"nam.tempreg.jja.png") 
           system("mv "+OUTDIR+"nam.000010.png "+OUTDIR+"nam.tempreg.son.png") 
           system("mv "+OUTDIR+"nam.000011.png "+OUTDIR+"nam.tempreg.ann.png")
           system("mv "+OUTDIR+"nam.000012.png "+OUTDIR+"nam.tempreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"nam.000006.png "+OUTDIR+"nam.tempreg.djf.png") 
           system("mv "+OUTDIR+"nam.000007.png "+OUTDIR+"nam.tempreg.mam.png") 
           system("mv "+OUTDIR+"nam.000008.png "+OUTDIR+"nam.tempreg.jja.png") 
           system("mv "+OUTDIR+"nam.000009.png "+OUTDIR+"nam.tempreg.son.png") 
           system("mv "+OUTDIR+"nam.000010.png "+OUTDIR+"nam.tempreg.ann.png")
        end if
     end if

     if (prreg_frame.eq.0) then
        system("mv "+OUTDIR+"nam.prreg.000001.png "+OUTDIR+"nam.prreg.djf.png") 
        system("mv "+OUTDIR+"nam.prreg.000002.png "+OUTDIR+"nam.prreg.mam.png") 
        system("mv "+OUTDIR+"nam.prreg.000003.png "+OUTDIR+"nam.prreg.jja.png") 
        system("mv "+OUTDIR+"nam.prreg.000004.png "+OUTDIR+"nam.prreg.son.png") 
        system("mv "+OUTDIR+"nam.prreg.000005.png "+OUTDIR+"nam.prreg.ann.png")         
        if (COMPUTE_MODES_MON.eq."True") then     
           system("mv "+OUTDIR+"nam.prreg.000006.png "+OUTDIR+"nam.prreg.mon.png")         
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
           system("mv "+OUTDIR+"nao.000007.png "+OUTDIR+"nao.tempreg.djf.png") 
           system("mv "+OUTDIR+"nao.000008.png "+OUTDIR+"nao.tempreg.mam.png") 
           system("mv "+OUTDIR+"nao.000009.png "+OUTDIR+"nao.tempreg.jja.png") 
           system("mv "+OUTDIR+"nao.000010.png "+OUTDIR+"nao.tempreg.son.png") 
           system("mv "+OUTDIR+"nao.000011.png "+OUTDIR+"nao.tempreg.ann.png")
           system("mv "+OUTDIR+"nao.000012.png "+OUTDIR+"nao.tempreg.mon.png")     
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"nao.000006.png "+OUTDIR+"nao.tempreg.djf.png") 
           system("mv "+OUTDIR+"nao.000007.png "+OUTDIR+"nao.tempreg.mam.png") 
           system("mv "+OUTDIR+"nao.000008.png "+OUTDIR+"nao.tempreg.jja.png") 
           system("mv "+OUTDIR+"nao.000009.png "+OUTDIR+"nao.tempreg.son.png") 
           system("mv "+OUTDIR+"nao.000010.png "+OUTDIR+"nao.tempreg.ann.png")
        end if
     end if
     
     if (prreg_frame.eq.0) then
        system("mv "+OUTDIR+"nao.prreg.000001.png "+OUTDIR+"nao.prreg.djf.png") 
        system("mv "+OUTDIR+"nao.prreg.000002.png "+OUTDIR+"nao.prreg.mam.png") 
        system("mv "+OUTDIR+"nao.prreg.000003.png "+OUTDIR+"nao.prreg.jja.png") 
        system("mv "+OUTDIR+"nao.prreg.000004.png "+OUTDIR+"nao.prreg.son.png") 
        system("mv "+OUTDIR+"nao.prreg.000005.png "+OUTDIR+"nao.prreg.ann.png")         
        if (COMPUTE_MODES_MON.eq."True") then     
           system("mv "+OUTDIR+"nao.prreg.000006.png "+OUTDIR+"nao.prreg.mon.png")         
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
     
     
  else
     system("psplit "+OUTDIR+"nam.ps "+OUTDIR+"psl_nn")
     system("mv "+OUTDIR+"psl_nn0001.ps "+OUTDIR+"nam.djf.ps") 
     system("mv "+OUTDIR+"psl_nn0002.ps "+OUTDIR+"nam.mam.ps") 
     system("mv "+OUTDIR+"psl_nn0003.ps "+OUTDIR+"nam.jja.ps") 
     system("mv "+OUTDIR+"psl_nn0004.ps "+OUTDIR+"nam.son.ps") 
     system("mv "+OUTDIR+"psl_nn0005.ps "+OUTDIR+"nam.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_nn0006.ps "+OUTDIR+"nam.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_nn0007.ps "+OUTDIR+"nam.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_nn0008.ps "+OUTDIR+"nam.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_nn0009.ps "+OUTDIR+"nam.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_nn0010.ps "+OUTDIR+"nam.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_nn0011.ps "+OUTDIR+"nam.tempreg.ann.ps")
           system("mv "+OUTDIR+"psl_nn0012.ps "+OUTDIR+"nam.tempreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_nn0006.ps "+OUTDIR+"nam.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_nn0007.ps "+OUTDIR+"nam.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_nn0008.ps "+OUTDIR+"nam.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_nn0009.ps "+OUTDIR+"nam.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_nn0010.ps "+OUTDIR+"nam.tempreg.ann.ps")     
        end if
     end if

     if (prreg_frame.eq.0) then
        system("psplit "+OUTDIR+"nam.prreg.ps "+OUTDIR+"pr_nn")
        system("mv "+OUTDIR+"pr_nn0001.ps "+OUTDIR+"nam.prreg.djf.ps") 
        system("mv "+OUTDIR+"pr_nn0002.ps "+OUTDIR+"nam.prreg.mam.ps") 
        system("mv "+OUTDIR+"pr_nn0003.ps "+OUTDIR+"nam.prreg.jja.ps") 
        system("mv "+OUTDIR+"pr_nn0004.ps "+OUTDIR+"nam.prreg.son.ps") 
        system("mv "+OUTDIR+"pr_nn0005.ps "+OUTDIR+"nam.prreg.ann.ps") 
        if (COMPUTE_MODES_MON.eq."True") then
           system("mv "+OUTDIR+"pr_nn0006.ps "+OUTDIR+"nam.prreg.mon.ps") 
        end if
     end if

     system("psplit "+OUTDIR+"nao.ps "+OUTDIR+"psl_nn")
     system("mv "+OUTDIR+"psl_nn0001.ps "+OUTDIR+"nao.djf.ps") 
     system("mv "+OUTDIR+"psl_nn0002.ps "+OUTDIR+"nao.mam.ps") 
     system("mv "+OUTDIR+"psl_nn0003.ps "+OUTDIR+"nao.jja.ps") 
     system("mv "+OUTDIR+"psl_nn0004.ps "+OUTDIR+"nao.son.ps") 
     system("mv "+OUTDIR+"psl_nn0005.ps "+OUTDIR+"nao.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_nn0006.ps "+OUTDIR+"nao.mon.ps") 
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_nn0007.ps "+OUTDIR+"nao.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_nn0008.ps "+OUTDIR+"nao.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_nn0009.ps "+OUTDIR+"nao.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_nn0010.ps "+OUTDIR+"nao.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_nn0011.ps "+OUTDIR+"nao.tempreg.ann.ps")
           system("mv "+OUTDIR+"psl_nn0012.ps "+OUTDIR+"nao.tempreg.mon.ps") 
        end if
     else
        if (sstreg_frame.eq.0) then
           system("mv "+OUTDIR+"psl_nn0006.ps "+OUTDIR+"nao.tempreg.djf.ps") 
           system("mv "+OUTDIR+"psl_nn0007.ps "+OUTDIR+"nao.tempreg.mam.ps") 
           system("mv "+OUTDIR+"psl_nn0008.ps "+OUTDIR+"nao.tempreg.jja.ps") 
           system("mv "+OUTDIR+"psl_nn0009.ps "+OUTDIR+"nao.tempreg.son.ps") 
           system("mv "+OUTDIR+"psl_nn0010.ps "+OUTDIR+"nao.tempreg.ann.ps")     
        end if
     end if
     
     if (prreg_frame.eq.0) then
        system("psplit "+OUTDIR+"nao.prreg.ps "+OUTDIR+"pr_nn")
        system("mv "+OUTDIR+"pr_nn0001.ps "+OUTDIR+"nao.prreg.djf.ps") 
        system("mv "+OUTDIR+"pr_nn0002.ps "+OUTDIR+"nao.prreg.mam.ps") 
        system("mv "+OUTDIR+"pr_nn0003.ps "+OUTDIR+"nao.prreg.jja.ps") 
        system("mv "+OUTDIR+"pr_nn0004.ps "+OUTDIR+"nao.prreg.son.ps") 
        system("mv "+OUTDIR+"pr_nn0005.ps "+OUTDIR+"nao.prreg.ann.ps") 
        if (COMPUTE_MODES_MON.eq."True") then
           system("mv "+OUTDIR+"pr_nn0006.ps "+OUTDIR+"nao.prreg.mon.ps") 
        end if
        system("rm "+OUTDIR+"nam.prreg.ps "+OUTDIR+"nao.prreg.ps")
     end if

     system("psplit "+OUTDIR+"nam.timeseries.ps "+OUTDIR+"psl_nn")
     system("mv "+OUTDIR+"psl_nn0001.ps "+OUTDIR+"nam.timeseries.djf.ps") 
     system("mv "+OUTDIR+"psl_nn0002.ps "+OUTDIR+"nam.timeseries.mam.ps") 
     system("mv "+OUTDIR+"psl_nn0003.ps "+OUTDIR+"nam.timeseries.jja.ps") 
     system("mv "+OUTDIR+"psl_nn0004.ps "+OUTDIR+"nam.timeseries.son.ps") 
     system("mv "+OUTDIR+"psl_nn0005.ps "+OUTDIR+"nam.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_nn0006.ps "+OUTDIR+"nam.timeseries.mon.ps") 
     end if
     
     system("psplit "+OUTDIR+"nao.timeseries.ps "+OUTDIR+"psl_nn")
     system("mv "+OUTDIR+"psl_nn0001.ps "+OUTDIR+"nao.timeseries.djf.ps") 
     system("mv "+OUTDIR+"psl_nn0002.ps "+OUTDIR+"nao.timeseries.mam.ps") 
     system("mv "+OUTDIR+"psl_nn0003.ps "+OUTDIR+"nao.timeseries.jja.ps") 
     system("mv "+OUTDIR+"psl_nn0004.ps "+OUTDIR+"nao.timeseries.son.ps") 
     system("mv "+OUTDIR+"psl_nn0005.ps "+OUTDIR+"nao.timeseries.ann.ps") 
     if (COMPUTE_MODES_MON.eq."True") then     
        system("mv "+OUTDIR+"psl_nn0006.ps "+OUTDIR+"nao.timeseries.mon.ps") 
     end if
     system("rm "+OUTDIR+"nao.timeseries.ps "+OUTDIR+"nam.timeseries.ps "+OUTDIR+"nao.ps "+OUTDIR+"nam.ps")     
  end if
  print("Finished: psl.nam_nao.ncl")
end  
