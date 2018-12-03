; use the user specified namelist / namelist_obs files to locate
; the files to be used, and write those file paths to namelist_byvar/namelist_*
; for use by other CVDP scripts.
;
; Note: ".nc" is removed from the paths given in namelist. 
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: namelist.ncl")
  o = getenv("OBS")
  case_sens = getenv("MACHINE")
  
  if (o.eq."True") then
     obsflag = True
  else
     obsflag = False
  end if
  
  if (.not.isfilepresent2("namelist_byvar/")) then
     system("mkdir namelist_byvar")
  end if
  
  nsim = numAsciiRow("namelist")
  na = asciiread("namelist",(/nsim/),"string")
  
  blankrow = ind(na.eq."")
  if (.not.any(ismissing(blankrow))) then
     goodrows = ind(na.ne."")
     na2 = na(goodrows)
     delete(na)
     na = na2
     delete(na2)
     nsim = dimsizes(na)
  end if     

  system("   export NSIM="+nsim)
  
  nentry = numAsciiCol("namelist")
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
  
  delete([/na,delim/])
  
    
  do gg = 0,nsim-1         ; if path ends in .nc remove it.
;     print(paths(gg))      ;  (It will get appended to the end of the path automatically when searching below.)
     paths(gg) = str_sub_str(paths(gg),".nc","")
;     print(paths(gg))
  end do

;----- Read in namelist_obs, and check number of supplied Observational datasets ------

  maxnumobs = 0    ; set maximum number of obs datasets per variable. if(obsflag).eq.True, this will likely get altered.

  if (obsflag) then
     nobs    = numAsciiRow("namelist_obs")
     nentryB = numAsciiCol("namelist_obs")
     vnamesB  = new(nobs,"string")
     namesB  = new(nobs,"string")
     pathsB  = new(nobs,"string")
     syearBT  = new(nobs,"string")
     eyearBT  = new(nobs,"string")
  
     na = asciiread("namelist_obs",(/nobs/),"string")
;     print(na)
  
     delim = "|"

     do gg = 0,nobs-1
        vnamesB(gg) = str_strip(str_get_field(na(gg),1,delim))
        namesB(gg)  = str_strip(str_get_field(na(gg),2,delim))
        pathsB(gg)  = str_strip(str_get_field(na(gg),3,delim))
        syearBT(gg) = str_strip(str_get_field(na(gg),4,delim))
        eyearBT(gg) = str_strip(str_get_field(na(gg),5,delim))
     end do
     namesB = where(namesB.eq."",namesB@_FillValue,namesB)
     pathsB = where(pathsB.eq."",pathsB@_FillValue,pathsB)
     syearBT = where(syearBT.eq."",syearBT@_FillValue,syearBT)
     eyearBT = where(eyearBT.eq."",eyearBT@_FillValue,eyearBT)
     
     maxnumobs = max((/dimsizes(ind(vnamesB.eq."TS")),dimsizes(ind(vnamesB.eq."PSL")),dimsizes(ind(vnamesB.eq."TREFHT")), \
                       dimsizes(ind(vnamesB.eq."PRECT")),dimsizes(ind(vnamesB.eq."MOC")),dimsizes(ind(vnamesB.eq."SNOWDP")), \
                      dimsizes(ind(vnamesB.eq."aice_nh")),dimsizes(ind(vnamesB.eq."aice_sh"))/))
     syearB = stringtointeger(syearBT)
     eyearB = stringtointeger(eyearBT)
     
     
     do gg = 0,nobs-1      ; check to see if any names are duplicated. If they are, add a "_2", "_3" to the name
        dupn = ind(namesB.eq.namesB(gg))            ; this is needed so that each output .nc file has a different name
        if (dimsizes(dupn).ge.2) then
           do hh = 1,dimsizes(dupn)-1
              namesB(dupn(hh)) = namesB(dupn(hh))+"_"+hh
           end do
        end if
        delete(dupn)
     end do
     delete([/eyearBT,syearBT/])
     delete([/na,delim,nentryB,nobs/])
     asciiwrite("obs_maxnum",maxnumobs)
  end if
;  print(vnamesB+" "+namesB+" "+pathsB+" "+syearB+" "+eyearB)
;exit
;----- TS section---------------
  namelist_ts = new(nsim+maxnumobs,string)
  if (obsflag) then
     ts_i = ind(vnamesB.eq."TS")
     if (.not.ismissing(ts_i(0))) then
        incr = dimsizes(ts_i)
        do gg = 0,incr-1
           namelist_ts(gg) = namesB(ts_i(gg))+" | "+pathsB(ts_i(gg))+" | "+syearB(ts_i(gg))+" | "+eyearB(ts_i(gg))  
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_ts(0))) then
              nmiss = ind(ismissing(namelist_ts(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_ts(hh) = namesB(ts_i(0))+"_"+hh+" | "+pathsB(ts_i(0))+" | "+syearB(ts_i(0))+" | "+eyearB(ts_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_ts",namelist_ts(0)) 
        delete(incr)
     end if
     delete(ts_i)
  end if  
  
  if (case_sens.eq."True") then
     tstring0 = "TS_,.TS.,ts_,.ts.,t_surf_,t_surf.,sst.,sst_"
     tstring1 = "TS,ts,t_surf,sst"
  else
     tstring0 = "TS_,.TS.,t_surf_,t_surf.,sst.,sst_"
     tstring1 = "TS,t_surf,sst"
  end if
  
  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{"+tstring0+"}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_ts(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{"+tstring1+"}/")     ; explicitly specify TS,ts in directory structure to eliminate "/tsmin/" being used
        namelist_ts(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{"+tstring0+"}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  delete([/tstring0,tstring1/])
  asciiwrite("namelist_byvar/namelist_ts",namelist_ts)
;------- PSL section----------------------------  
  namelist_psl = new(nsim+maxnumobs,string)
  if (obsflag) then
     psl_i = ind(vnamesB.eq."PSL")
     if (.not.ismissing(psl_i(0))) then
        incr = dimsizes(psl_i)
        do gg = 0,incr-1
           namelist_psl(gg) = namesB(psl_i(gg))+" | "+pathsB(psl_i(gg))+" | "+syearB(psl_i(gg))+" | "+eyearB(psl_i(gg))  
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_psl(0))) then
              nmiss = ind(ismissing(namelist_psl(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_psl(hh) = namesB(psl_i(0))+"_"+hh+" | "+pathsB(psl_i(0))+" | "+syearB(psl_i(0))+" | "+eyearB(psl_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_psl",namelist_psl(0)) 
        delete(incr)
     end if
     delete(psl_i)
  end if  
  
  if (case_sens.eq."True") then
     tstring0 = "PSL_,PSL.,psl_,psl.,slp.,slp_"
     tstring1 = "PSL,psl,SLP,slp"
  else
     tstring0 = "PSL_,PSL.,slp.,slp_"
     tstring1 = "PSL,slp"
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{"+tstring0+"}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_psl(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{"+tstring1+"}/")
        namelist_psl(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{"+tstring0+"}*.nc | "+syear(gg)+" | "+eyear(gg)  
        delete(tpath)
     end if
     delete(fsst)
  end do  
  delete([/tstring0,tstring1/])
  asciiwrite("namelist_byvar/namelist_psl",namelist_psl)
;------- TREFHT section----------------------------  
  namelist_trefht = new(nsim+maxnumobs,string)
  if (obsflag) then
     trefht_i = ind(vnamesB.eq."TREFHT")
     if (.not.ismissing(trefht_i(0))) then
        incr = dimsizes(trefht_i)
        do gg = 0,incr-1
           namelist_trefht(gg) = namesB(trefht_i(gg))+" | "+pathsB(trefht_i(gg))+" | "+syearB(trefht_i(gg))+" | "+eyearB(trefht_i(gg)) 
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_trefht(0))) then
              nmiss = ind(ismissing(namelist_trefht(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_trefht(hh) = namesB(trefht_i(0))+"_"+hh+" | "+pathsB(trefht_i(0))+" | "+syearB(trefht_i(0))+" | "+eyearB(trefht_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_trefht",namelist_trefht(0)) 
        delete(incr)
     end if
     delete(trefht_i)
  end if 
  
  if (case_sens.eq."True") then
     tstring0 = "TREFHT_,TREFHT.,tas.,tas_,t_ref.,t_ref_,T2.,T2_"
     tstring1 = "TREFHT,tas,t_ref,T2"
  else
     tstring0 = "TREFHT_,TREFHT.,tas.,tas_,t_ref.,t_ref_,T2.,T2_"
     tstring1 = "TREFHT,tas,t_ref,T2"
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{"+tstring0+"}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_trefht(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{"+tstring1+"}/")
        namelist_trefht(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{"+tstring0+"}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_trefht",namelist_trefht)  
;------- PRECT section--(more complicated due to PRECC+PRECL, + pr being a common 2 letter combination)------  
  namelist_prect = new(nsim+maxnumobs,string)
  if (obsflag) then
     prect_i = ind(vnamesB.eq."PRECT")
     if (.not.ismissing(prect_i(0))) then
        incr = dimsizes(prect_i)
        do gg = 0,incr-1
           namelist_prect(gg) = namesB(prect_i(gg))+" | "+pathsB(prect_i(gg))+" | "+syearB(prect_i(gg))+" | "+eyearB(prect_i(gg)) 
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-ouprect
           if (.not.ismissing(namelist_prect(0))) then
              nmiss = ind(ismissing(namelist_prect(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_prect(hh) = namesB(prect_i(0))+"_"+hh+" | "+pathsB(prect_i(0))+" | "+syearB(prect_i(0))+" | "+eyearB(prect_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_prect",namelist_prect(0)) 
        delete(incr)
     end if
     delete(prect_i)
  end if 
  do gg = 0,nsim-1
     pstring = (/"pr_*","pr.*","_pr_*",".pr.*","PRECT.*","PRECT_*","PRECC.*","PRECC_*","precip_*","precip.*","prcp_*","prcp.*","prate_*","prate.*"/)
     do hh = 0,dimsizes(pstring)-1
        if (isvar((/"fsst"/))) then
           delete(fsst)
        end if
        fsst = systemfunc("bash -c 'ls "+paths(gg)+"*"+pstring(hh)+".nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
        if (dimsizes(fsst).eq.1) then
           if (.not.ismissing(fsst)) then
              if (pstring(hh).eq."PRECC.*".or.pstring(hh).eq."PRECC_*") then
                 tpath = str_sub_str(paths(gg),"/*/","/{PRECC,PRECL}/")
                 namelist_prect(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{PRECC,PRECL}*.nc | "+syear(gg)+" | "+eyear(gg)
                 delete(tpath)
              else
                 namelist_prect(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file         
              end if
              break        
           else
               namelist_prect(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; file is missing..          
           end if
        else
           if (pstring(hh).eq."PRECC.*".or.pstring(hh).eq."PRECC_*") then
              tpath = str_sub_str(paths(gg),"/*/","/{PRECC,PRECL}/")
              namelist_prect(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{PRECC,PRECL}*.nc | "+syear(gg)+" | "+eyear(gg)
           else
              if (pstring(hh).eq."pr_*".or.pstring(hh).eq."pr.*".or.pstring(hh).eq."*_pr_*".or.pstring(hh).eq."*.pr.*") then
                 tpath = str_sub_str(paths(gg),"/*/","/pr/")
              end if
              if (pstring(hh).eq."PRECT.*".or.pstring(hh).eq."PRECT_*") then
                 tpath = str_sub_str(paths(gg),"/*/","/PRECC/")
              end if
              if (isvar((/"tpath"/))) then
                 namelist_prect(gg+maxnumobs) =  names(gg)+" | "+tpath+"*"+pstring(hh)+"*.nc | "+syear(gg)+" | "+eyear(gg)
                 delete(tpath)
              else
                 namelist_prect(gg+maxnumobs) =  names(gg)+" | missing"
              end if
           end if
           break
        end if
     end do
  end do  
  delete(pstring)
  if (isvar((/"fsst"/))) then
     delete(fsst)
  end if
  asciiwrite("namelist_byvar/namelist_prect",namelist_prect)    
;----- SNOWDP section---------------
  namelist_snowdp = new(nsim+maxnumobs,string)
  if (obsflag) then
     snowdp_i = ind(vnamesB.eq."SNOWDP")
     if (.not.ismissing(snowdp_i(0))) then
        incr = dimsizes(snowdp_i)
        do gg = 0,incr-1
           namelist_snowdp(gg) = namesB(snowdp_i(gg))+" | "+pathsB(snowdp_i(gg))+" | "+syearB(snowdp_i(gg))+" | "+eyearB(snowdp_i(gg)) 
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_snowdp(0))) then
              nmiss = ind(ismissing(namelist_snowdp(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_snowdp(hh) = namesB(snowdp_i(0))+"_"+hh+" | "+pathsB(snowdp_i(0))+" | "+syearB(snowdp_i(0))+" | "+eyearB(snowdp_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_snowdp",namelist_snowdp(0)) 
        delete(incr)
     end if
     delete(snowdp_i)
  end if 
  sn_string = "SNOWDP_,SNOWDP.,snd_,snd."
  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{"+sn_string+"}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_snowdp(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{SNOWDP,snd}/")     ; explicitly specify SNOWDP/snd in directory structure to eliminate "/sndmin/" being used
        namelist_snowdp(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{SNOWDP_,SNOWDP.,snd_,snd.}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_snowdp",namelist_snowdp)
;------- MOC section----------------------------  
  namelist_moc = new(nsim+maxnumobs,string)
  if (obsflag) then
     moc_i = ind(vnamesB.eq."MOC")
     if (.not.ismissing(moc_i(0))) then
        incr = dimsizes(moc_i)
        do gg = 0,incr-1
           namelist_moc(gg) = namesB(moc_i(gg))+" | "+pathsB(moc_i(gg))+" | "+syearB(moc_i(gg))+" | "+eyearB(moc_i(gg)) 
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_moc(0))) then
              nmiss = ind(ismissing(namelist_moc(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_moc(hh) = namesB(moc_i(0))+"_"+hh+" | "+pathsB(moc_i(0))+" | "+syearB(moc_i(0))+" | "+eyearB(moc_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_moc",namelist_moc(0)) 
        delete(incr)
     end if
     delete(moc_i)
  end if 

  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{MOC_,MOC.,msftmyz.,msftmyz_,stfmmc.,stfmmc_}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_moc(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{MOC,msftmyz,stfmmc}/")
        namelist_moc(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{MOC_,MOC.,msftmyz.,msftmyz_,stfmmc.,stfmmc_}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_moc",namelist_moc) 

;------- aice_nh section----------------------------  
  namelist_aice_nh = new(nsim+maxnumobs,string)
  if (obsflag) then
     aice_nh_i = ind(vnamesB.eq."aice_nh".or.vnamesB.eq."AICE_NH")
     if (.not.ismissing(aice_nh_i(0))) then
        incr = dimsizes(aice_nh_i)
        do gg = 0,incr-1
           namelist_aice_nh(gg) = namesB(aice_nh_i(gg))+" | "+pathsB(aice_nh_i(gg))+" | "+syearB(aice_nh_i(gg))+" | "+eyearB(aice_nh_i(gg)) 
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_aice_nh(0))) then
              nmiss = ind(ismissing(namelist_aice_nh(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_aice_nh(hh) = namesB(aice_nh_i(0))+"_"+hh+" | "+pathsB(aice_nh_i(0))+" | "+syearB(aice_nh_i(0))+" | "+eyearB(aice_nh_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_aice_nh",namelist_aice_nh(0)) 
        delete(incr)
     end if
     delete(aice_nh_i)
  end if 

  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{aice_nh.,aice.,sic_,sic.,.CN.,_CN_}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_aice_nh(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{aice,sic,aice_nh,CN}/")
        namelist_aice_nh(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{aice_nh.,aice.,sic_,sic.,.CN.,_CN_}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_aice_nh",namelist_aice_nh) 
;------- aice_sh section----------------------------  
  namelist_aice_sh = new(nsim+maxnumobs,string)
  if (obsflag) then
     aice_sh_i = ind(vnamesB.eq."aice_sh".or.vnamesB.eq."AICE_SH")
     if (.not.ismissing(aice_sh_i(0))) then
        incr = dimsizes(aice_sh_i)
        do gg = 0,incr-1
           namelist_aice_sh(gg) = namesB(aice_sh_i(gg))+" | "+pathsB(aice_sh_i(gg))+" | "+syearB(aice_sh_i(gg))+" | "+eyearB(aice_sh_i(gg)) 
        end do
        if (incr.lt.maxnumobs) then     ; fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
           if (.not.ismissing(namelist_aice_sh(0))) then
              nmiss = ind(ismissing(namelist_aice_sh(:maxnumobs-1)))
              do hh = nmiss(0),nmiss(dimsizes(nmiss)-1)
                 namelist_aice_sh(hh) = namesB(aice_sh_i(0))+"_"+hh+" | "+pathsB(aice_sh_i(0))+" | "+syearB(aice_sh_i(0))+" | "+eyearB(aice_sh_i(0))
              end do
              delete(nmiss)
           end if
        end if
        asciiwrite("obs_aice_sh",namelist_aice_sh(0)) 
        delete(incr)
     end if
     delete(aice_sh_i)
  end if 

  do gg = 0,nsim-1
     fsst = systemfunc("bash -c 'ls "+paths(gg)+"*{aice_sh.,aice.,sic_,sic.,.CN.,_CN_}*.nc 2> /dev/null'")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_aice_sh(gg+maxnumobs) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{aice,sic,aice_sh,CN}/")
        namelist_aice_sh(gg+maxnumobs) =  names(gg)+" | "+tpath+"*{aice_sh.,aice.,sic_,sic.,.CN.,_CN_}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_aice_sh",namelist_aice_sh) 
;----------------------------------------------------------------------------
  if (obsflag) then
     delete([/vnamesB,namesB,pathsB,syearB,eyearB/])
  end if
  print("Finished: namelist.ncl")
end
