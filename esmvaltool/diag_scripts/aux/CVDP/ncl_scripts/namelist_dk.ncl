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
  print("Starting: Climate Variability Diagnostics Package")
  print("Starting: namelist.ncl")
  o = getenv("OBS")
  
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

     syearB = stringtointeger(syearBT)
     eyearB = stringtointeger(eyearBT)
     delete([/eyearBT,syearBT/])
     delete([/na,delim,nentryB,nobs/])
  end if
;  print(vnamesB+" "+namesB+" "+pathsB+" "+syearB+" "+eyearB)
;----- TS section---------------
  if (obsflag) then
     if (.not.ismissing(pathsB(0))) then
        namelist_ts = new(nsim+1,string)
        incr = 1
        namelist_ts(0) = namesB(0)+" | "+pathsB(0)+" | "+syearB(0)+" | "+eyearB(0)   ; grab first OBS file, which = TS
        asciiwrite("obs_ts",namelist_ts(0))    ; write file signaling that observational TS file is specified
     else
        namelist_ts = new(nsim,string)
        incr = 0
     end if
  else
     namelist_ts = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{TS_,TS.,ts_,ts.}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_ts(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{TS,ts}/")     ; explicitly specify TS,ts in directory structure to eliminate "/tsmin/" being used
        namelist_ts(gg+incr) =  names(gg)+" | "+tpath+"*{TS_,TS.,ts_,ts.}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_ts",namelist_ts)
;------- PSL section----------------------------  
  if (obsflag) then
     if (.not.ismissing(pathsB(1))) then
        namelist_psl = new(nsim+1,string)
        incr = 1
        namelist_psl(0) = namesB(1)+" | "+pathsB(1)+" | "+syearB(1)+" | "+eyearB(1)   ; grab 2nd OBS file, which = psl
        asciiwrite("obs_psl",namelist_psl(0))    ; write file signaling that observational psl file is specified
     else
        namelist_psl = new(nsim,string)
        incr = 0
     end if
  else
     namelist_psl = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{PSL_,PSL.,psl_,psl.,slp.,slp_}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_psl(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{PSL,psl}/")
        namelist_psl(gg+incr) =  names(gg)+" | "+tpath+"*{PSL_,PSL.,psl_,psl.,slp.,slp_}*.nc | "+syear(gg)+" | "+eyear(gg)  
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_psl",namelist_psl)
;------- TREFHT section----------------------------  
  if (obsflag) then
     if (.not.ismissing(pathsB(2))) then
        namelist_trefht = new(nsim+1,string)
        incr = 1
        namelist_trefht(0) = namesB(2)+" | "+pathsB(2)+" | "+syearB(2)+" | "+eyearB(2)   ; grab 2nd OBS file, which = trefht
        asciiwrite("obs_trefht",namelist_trefht(0))    ; write file signaling that observational trefht file is specified
     else
        namelist_trefht = new(nsim,string)
        incr = 0
     end if
  else
     namelist_trefht = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{TREFHT_,TREFHT.,tas.,tas_}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_trefht(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{TREFHT,tas}/")
        namelist_trefht(gg+incr) =  names(gg)+" | "+tpath+"*{TREFHT_,TREFHT.,tas.,tas_}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_trefht",namelist_trefht)  
;------- PRECT section--(more complicated due to PRECC+PRECL, + pr being a common 2 letter combination)------  
  if (obsflag) then
     if (.not.ismissing(pathsB(3))) then
        namelist_prect = new(nsim+1,string)
        incr = 1
        namelist_prect(0) = namesB(3)+" | "+pathsB(3)+" | "+syearB(3)+" | "+eyearB(3)   ; grab 4th OBS file, which = prect
        asciiwrite("obs_prect",namelist_prect(0))    ; write file signaling that observational prect file is specified
     else
        namelist_prect = new(nsim,string)
        incr = 0
     end if
  else
     namelist_prect = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     pstring = (/"pr_*","pr.*","_pr_*",".pr.*","PRECT.*","PRECT_*","PRECC.*","PRECC_*"/)
     do hh = 0,dimsizes(pstring)-1
        if (isvar((/"fsst"/))) then
           delete(fsst)
        end if
        fsst = systemfunc("ls "+paths(gg)+"*"+pstring(hh)+".nc 2> /dev/null")   ; /dev/null suppresses all standard error output
        if (dimsizes(fsst).eq.1) then
           if (.not.ismissing(fsst)) then
              if (pstring(hh).eq."PRECC.*".or.pstring(hh).eq."PRECC_*") then
                 tpath = str_sub_str(paths(gg),"/*/","/{PRECC,PRECL}/")
                 namelist_prect(gg+incr) =  names(gg)+" | "+tpath+"*{PRECC,PRECL}*.nc | "+syear(gg)+" | "+eyear(gg)
                 delete(tpath)
              else
                 namelist_prect(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file         
              end if
              break        
           else
               namelist_prect(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; file is missing..          
           end if
        else
           if (pstring(hh).eq."PRECC.*".or.pstring(hh).eq."PRECC_*") then
              tpath = str_sub_str(paths(gg),"/*/","/{PRECC,PRECL}/")
              namelist_prect(gg+incr) =  names(gg)+" | "+tpath+"*{PRECC,PRECL}*.nc | "+syear(gg)+" | "+eyear(gg)
           else
              if (pstring(hh).eq."pr_*".or.pstring(hh).eq."pr.*".or.pstring(hh).eq."*_pr_*".or.pstring(hh).eq."*.pr.*") then
                 tpath = str_sub_str(paths(gg),"/*/","/pr/")
              end if
              if (pstring(hh).eq."PRECT.*".or.pstring(hh).eq."PRECT_*") then
                 tpath = str_sub_str(paths(gg),"/*/","/PRECC/")
              end if
              if (isvar((/"tpath"/))) then
                 namelist_prect(gg+incr) =  names(gg)+" | "+tpath+"*"+pstring(hh)+"*.nc | "+syear(gg)+" | "+eyear(gg)
                 delete(tpath)
              else
                 namelist_prect(gg+incr) =  names(gg)+" | missing"
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
  if (obsflag) then
     if (.not.ismissing(pathsB(4))) then
        namelist_snowdp = new(nsim+1,string)
        incr = 1
        namelist_snowdp(0) = namesB(4)+" | "+pathsB(4)+" | "+syearB(4)+" | "+eyearB(4)   ; grab fifth OBS file, which = SNOWDP
        asciiwrite("obs_snowdp",namelist_snowdp(0))    ; write file signaling that observational TS file is specified
     else
        namelist_snowdp = new(nsim,string)
        incr = 0
     end if
  else
     namelist_snowdp = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{SNOWDP_,SNOWDP.,snd_,snd.}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_snowdp(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{SNOWDP,snd}/")     ; explicitly specify SNOWDP/snd in directory structure to eliminate "/sndmin/" being used
        namelist_snowdp(gg+incr) =  names(gg)+" | "+tpath+"*{SNOWDP_,SNOWDP.,snd_,snd.}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_snowdp",namelist_snowdp)
;------- MOC section----------------------------  
  if (obsflag) then
     if (.not.ismissing(pathsB(5))) then
        namelist_moc = new(nsim+1,string)
        incr = 1
        namelist_moc(0) = namesB(5)+" | "+pathsB(5)+" | "+syearB(5)+" | "+eyearB(5)   ; grab 2nd OBS file, which = moc
        asciiwrite("obs_moc",namelist_moc(0))    ; write file signaling that observational moc file is specified
     else
        namelist_moc = new(nsim,string)
        incr = 0
     end if
  else
     namelist_moc = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{MOC_,MOC.,msftmyz.,msftmyz_,stfmmc.,stfmmc_}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_moc(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","/{MOC,msftmyz,stfmmc}/")
        namelist_moc(gg+incr) =  names(gg)+" | "+tpath+"*{MOC_,MOC.,msftmyz.,msftmyz_,stfmmc.,stfmmc_}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_moc",namelist_moc) 
  
;----- CLT section---------------
  if (obsflag) then
     if (.not.ismissing(pathsB(6))) then
        namelist_clt = new(nsim+1,string)
        incr = 1
        namelist_clt(0) = namesB(6)+" | "+pathsB(6)+" | "+syearB(6)+" | "+eyearB(6)   ; grab sixth OBS file, which = CLT
        asciiwrite("obs_clt",namelist_clt(0))    ; write file signaling that observational TS file is specified
     else
        namelist_clt = new(nsim,string)
        incr = 0
     end if
  else
     namelist_clt = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{CLT_,CLT.,clt_,clt.}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_clt(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","{/CLT,clt/}")     ; explicitly specify CLT/clt in directory structure to eliminate "/sndmin/" being used
        namelist_clt(gg+incr) =  names(gg)+" | "+tpath+"*{CLT_,CLT.,clt_,clt.}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_clt",namelist_clt)
  
;----- CLT section---------------
  if (obsflag) then
     if (.not.ismissing(pathsB(6))) then
        namelist_clt = new(nsim+1,string)
        incr = 1
        namelist_clt(0) = namesB(6)+" | "+pathsB(6)+" | "+syearB(6)+" | "+eyearB(6)   ; grab sixth OBS file, which = CLT
        asciiwrite("obs_clt",namelist_clt(0))    ; write file signaling that observational TS file is specified
     else
        namelist_clt = new(nsim,string)
        incr = 0
     end if
  else
     namelist_clt = new(nsim,string)
     incr = 0
  end if
  do gg = 0,nsim-1
     fsst = systemfunc("ls "+paths(gg)+"*{CLT_,CLT.,clt_,clt.}*.nc 2> /dev/null")   ; /dev/null suppresses all standard error output
     if (dimsizes(fsst).eq.1) then
        namelist_clt(gg+incr) =  names(gg)+" | "+fsst(0)+" | "+syear(gg)+" | "+eyear(gg)   ; grab first file
     else
        tpath = str_sub_str(paths(gg),"/*/","{/CLT,clt/}")     ; explicitly specify CLT/clt in directory structure to eliminate "/sndmin/" being used
        namelist_clt(gg+incr) =  names(gg)+" | "+tpath+"*{CLT_,CLT.,clt_,clt.}*.nc | "+syear(gg)+" | "+eyear(gg)
        delete(tpath)
     end if
     delete(fsst)
  end do  
  asciiwrite("namelist_byvar/namelist_clt",namelist_clt)
   
  if (obsflag) then
     delete([/vnamesB,namesB,pathsB,syearB,eyearB/])
  end if
  print("Finished: namelist.ncl")
end
