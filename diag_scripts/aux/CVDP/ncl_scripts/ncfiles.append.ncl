; Concatenate all .nc files from same model/observational dataset 
; into a single .nc file.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
;load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: ncfiles.append.ncl")

  OUTDIR      = getenv("OUTDIR")
  o = getenv("OBS")
;  
  if (o.eq."True") then
     obsflag = True
  else
     obsflag = False
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
  
  nentry = numAsciiCol("namelist")
  names = new(nsim,"string")
  syear = new(nsim,"integer",-999)
  eyear = new(nsim,"integer",-999)
    
  delim = "|"
  do gg = 0,nsim-1
     names(gg) = str_sub_str(str_sub_str(str_sub_str(str_sub_str(str_sub_str(str_strip(str_get_field(na(gg),1,delim))," ","_"),"/","_"),"'","_"),"(","_"),")","_")
     syear(gg) = stringtointeger(str_strip(str_get_field(na(gg),3,delim)))
     eyear(gg) = stringtointeger(str_strip(str_get_field(na(gg),4,delim)))
  end do

  do gg = 0,nsim-1
     fils = systemfunc("ls "+OUTDIR+names(gg)+".*.nc 2> /dev/null")
     if (.not.ismissing(fils(0))) then
        dimf = dimsizes(fils)
        ofile = OUTDIR+names(gg)+".cvdp_data."+syear(gg)+"-"+eyear(gg)+".nc"
        if (dimf.eq.1) then
           system("mv "+fils(0)+" "+ofile)
        else
           if (isfilepresent(ofile)) then                                    ; if file master is present append individual data files to file master.
              do hh = 0,dimf-1
                 if (fils(hh).ne.ofile) then
                    system("ncks -A -h "+fils(hh)+" "+ofile)
                 end if
              end do
           else                                                              ; if file master is not present, append individual data files to last file in list,
              do hh = 0,dimf-2                                               ; and when done move the last file to be the master file
                 system("ncks -A -h "+fils(hh)+" "+fils(dimf-1))             
              end do
              system("mv "+fils(dimf-1)+" "+ofile)
           end if
           if (dimsizes(fils(:dimf-2)).ge.2) then
              system("rm "+str_sub_str(str_join(fils(:dimf-2)," "),ofile,""))   ; remove each script's file, but do not remove the master file (if present)
           end if
        end if
        system("ncks -O "+ofile+" "+ofile)   ; done to alphabetize output variable
        delete([/dimf,ofile/])
     else
;        print("NetCDF files not found for "+names+", skipping appending")
     end if
     delete(fils)
  end do
  delete([/nsim,na,blankrow,nentry,names,syear,eyear/])
;------------------------------------------------
  if (obsflag) then
     maxnumobs  = asciiread("obs_maxnum",(/1/),"integer")  

     namelist_files = (/"psl","prect","trefht","ts","snowdp","moc","aice_nh","aice_sh"/)
     delim = "|"
     cntr = 0
     namesB = new(maxnumobs*dimsizes(namelist_files),string)
     do gg = 0,dimsizes(namelist_files)-1                    ; grab all observational dataset names from namelist_$var files
        na = asciiread("namelist_byvar/namelist_"+namelist_files(gg),(/maxnumobs/),"string")
        namesB(cntr:cntr+maxnumobs-1) = str_sub_str(str_sub_str(str_sub_str(str_sub_str(str_sub_str(str_strip(str_get_field(na,1,delim))," ","_"),"/","_"),"'","_"),"(","_"),")","_")
        cntr = cntr+maxnumobs
        delete(na)
     end do

     namesB = where(namesB.eq."",namesB@_FillValue,namesB)     ; for blank names set them to _FillValue
     if (any(namesB.eq."missing")) then
        namesB(str_match_ind(namesB,"missing")) = namesB@_FillValue ; check for any names containing "missing", set to _FillValue
     end if
     delete([/delim,cntr,namelist_files/])

     do gg = 0,dimsizes(namesB)-1
        if (.not.ismissing(namesB(gg))) then
           fils = systemfunc("ls "+OUTDIR+namesB(gg)+".cvdp_data.*.nc 2> /dev/null")
           if (.not.ismissing(fils(0))) then
              dimf = dimsizes(fils)
              fil0 = tochar(fils(0))
              suffix = tostring(fil0(dimsizes(fil0)-12:dimsizes(fil0)-1))
              delete(fil0)
              ofi = OUTDIR+namesB(gg)+".cvdp_data."+suffix
              if (dimf.ge.2) then
                 if (isfilepresent(ofi)) then                   ; if file master is present append individual data files to file master.
                    do hh = 0,dimf-1
                       if (fils(hh).ne.ofi) then
                          system("ncks -A -h "+fils(hh)+" "+ofi)
                       end if
                    end do
                 else                                                        ; if file master is not present, append individual data files to last file in list,
                    do hh = 0,dimf-2                                         ; and when done move the last file to be the master file
                       system("ncks -A -h "+fils(hh)+" "+fils(dimf-1))
                    end do
                    system("mv "+fils(dimf-1)+" "+ofi)
                 end if

                 if (dimsizes(fils(:dimf-2)).ge.2) then
                    system("rm "+str_sub_str(str_join(fils(:dimf-2)," "),ofi,""))   ; remove each script's file, but do not remove the master file (if present)
                 end if
              else
                 if (fils(0).ne.ofi) then
                    system("mv "+fils(0)+" "+ofi)
                 end if
              end if
              system("ncks -O "+ofi+" "+ofi)   ; done to alphabetize output variable
              delete([/dimf,ofi/])
           else
;              print("NetCDF files not found for "+namesB(gg)+", skipping appending")
           end if
           delete(fils)
        end if
     end do
     delete([/namesB/])
  end if
  print("Finished: ncfiles.append.ncl")
end
