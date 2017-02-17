;=================================================================================================
; create blank array for use when something may be/is wrong.
;
function create_empty_array(yS:numeric,yE:numeric,mS:numeric,mE:numeric,opttype:string)
local yS,yE,mS,mE,opttype
begin
  timeT =  yyyymm_time(yS, yE, "integer")
  time = timeT({yS*100+mS:yE*100+mE})
  if (opttype.eq."time_lat_lon") then  
     blankarr = new((/dimsizes(time),90,180/),"float",1.e20)
     blankarr!0 = "time"    ; time coordinate variable assigned below  
     blankarr&time = time
     blankarr!1 = "lat"
     blankarr&lat = fspan(-89,89,90)
     blankarr!2 = "lon"
     blankarr&lon = fspan(0,358,180)
     blankarr&lat@units = "degrees_north"
     blankarr&lon@units = "degrees_east"
  end if
  if (opttype.eq."time_lev_lat") then  
     blankarr = new((/dimsizes(time),41,90/),"float",1.e20)
     blankarr!0 = "time"    ; time coordinate variable assigned below  
     blankarr&time = time
     blankarr!1 = "lev"
     blankarr&lev =fspan(0,5000,41)
     blankarr!2 = "lat"
     blankarr&lat = fspan(-89,89,90)
     blankarr&lat@units = "degrees_north"
     blankarr&lev@units = "m"
     blankarr&lev@positive = "down"
  end if
  blankarr@units = ""
  blankarr@is_all_missing = True
  return(blankarr)
  delete([/yS,yE,mS,mE,opttype,blankarr,timeT,time/])
end
;===================================================================================================
; read in atmospheric/land data from selected files
; assign time coordinate variables, check for issues with the array, assign _FillValue (if needed)
; assign dimension names (for ease-of-use), check and modify units
;
; vname settings at top of this script can be modified if a different variable name is
; encountered. For instance, if a TS data file has the TS array named as "sfc_t", one
; could add "sfc_t" to the vname TS coding as follows:
;   if (vn.eq."TS") then
;      vname = (/"TS","ts","sst","sfc_t"/)
;   end if
;
function data_read_in(zpath:string,vn:string,yearS:integer,yearE:integer)
; path for TS file(s), variable name, start year, and end year are read in.
local zpath,vn,cpath0,ta,tfiles,c,arr,farr,yearS,yearE,mocheck,fils_precc,fils_precl
begin
;  setfileoption("nc","Format","NetCDF4Classic")
;  setfileoption("nc","Compression",1)

  if (vn.eq."TS") then
     vname = (/"TS","ts","sst"/)
  end if
  if (vn.eq."PSL") then
     vname = (/"PSL","psl","slp","SLP","prmsl"/)
  end if
  if (vn.eq."TREFHT") then
     vname = (/"TREFHT","tas","temp","air","temperature_anomaly"/)
  end if
  if (vn.eq."PRECT") then
     vname = (/"PRECC","PRECL","PRECT","pr","PPT","ppt","p","P","precip","PRECIP"/)
  end if
  if (vn.eq."SNOWDP") then
     vname = (/"SNOWDP","snd"/)
  end if
    
  if (ismissing(zpath) ) then
     print("File missing, creating blank array of data. View "+vn+" namelist for details.")
     arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12
  else
     cpath0 = stringtochar(zpath)
  
     ta = stringtochar("*")
     if (any(cpath0.eq.ta(0)).or.any(cpath0.eq."{")) then   ; check for "*" and "{" denoting multiple files
        tfiles = systemfunc("ls "+zpath+" 2> /dev/null")   ; /dev/null suppresses all standard error output
        if (vn.eq."PRECT") then   ; special section for precip, as might need to do PRECC+PRECL
           b = addfile(tfiles(0),"r")   ; check for PRECC
           if (isfilevar(b,"PRECC").or.isfilevar(b,"PRECL")) then  ; PRECC/PRECL section
              fils_precc = str_match(tfiles,"PRECC") 
              fils_precl = str_match(tfiles,"PRECL") 
              if (any(ismissing(fils_precc)).or.any(ismissing(fils_precl))) then
                 print("Fatal: Need both PRECC and PRECL file(s), creating blank array")
                 print(fils_precc)
                 print(fils_precl)
                 arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
                 sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
                 eydata = yearE     ; as data array is totally missing..  
                 smdata = 1
                 emdata = 12
                 break
              end if
              c = addfiles(fils_precc,"r")
              arr = c[:]->PRECC
              c2  = addfiles(fils_precl,"r")
              arr = (/ arr+c2[:]->PRECL  /)
              arr@long_name = "Large-scale (stable) + convective precipitation rate (liq + ice)"
              delete([/c2,fils_precc,fils_precl/])            
           else    ; pr, ppt, PPT, PRECT multiple/single file read-in here..
              c = addfiles(tfiles,"r")
              do ii=0,dimsizes(vname)-1
                 if (isfilevar(c[0],vname(ii))) then
                    arr = c[:]->$vname(ii)$
                    break
                 end if
              end do
           end if
           delete(b)
        else
           c = addfiles(tfiles,"r")
           do ii=0,dimsizes(vname)-1
              if (isfilevar(c[0],vname(ii))) then
                 arr = c[:]->$vname(ii)$
                 break
              end if
           end do
        end if     
        nfil = dimsizes(tfiles)
        cpathS = stringtochar(tfiles(0))      ; this section will work for PRECC/PRECL, as it will read the first
        cpathE = stringtochar(tfiles(nfil-1)) ; PRECC file and the last PRECL file.  
        ncharS = dimsizes(cpathS)
        ncharE = dimsizes(cpathE)
        sydata = stringtointeger(charactertostring(cpathS(ncharS-17:ncharS-14)))
        smdata = stringtointeger(charactertostring(cpathS(ncharS-13:ncharS-12)))
        eydata = stringtointeger(charactertostring(cpathE(ncharE-10:ncharE-7)))
        emdata = stringtointeger(charactertostring(cpathE(ncharE-6:ncharE-5)))      
        delete([/cpathS,cpathE,ncharS,ncharE,nfil/])
;        delete(c)
     else
        c = addfile(zpath,"r")
        do i=0,dimsizes(vname)-1
           if (isfilevar(c,vname(i))) then
              arr = c->$vname(i)$
              break
           end if
        end do
        cpath = stringtochar(zpath)     
        nchar = dimsizes(cpath)
        sydata = stringtointeger(charactertostring(cpath(nchar-17:nchar-14)))
        smdata = stringtointeger(charactertostring(cpath(nchar-13:nchar-12)))
        eydata = stringtointeger(charactertostring(cpath(nchar-10:nchar-7)))
        emdata = stringtointeger(charactertostring(cpath(nchar-6:nchar-5))) 
        delete([/cpath,nchar/])
;        delete(c)     
     end if
     delete([/ta,cpath0/])
  end if
  
  if (isvar("arr").eq.False) then
     print("Variable ("+vn+") not found. Examine input file "+zpath+". Creating empty array and continuing")
     arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
  end if

  if (isshort(arr)) then
     arrT = short2flt(arr)
     delete(arr)
     arr = arrT
     delete(arrT)
  end if
  
  if (.not.isatt(arr,"_FillValue")) then     ; assign _FillValue if one is not present
     if (isatt(arr,"missing_value")) then
        arr@_FillValue = arr@missing_value
     else
        arr@_FillValue = default_fillvalue(typeof(arr)) 
     end if
  end if
  
  dimz = dimsizes(arr)
  if (any(dimz.eq.1)) then
     arrT = rm_single_dims(arr)
     delete(arr)
     arr = arrT
     delete(arrT)
  end if
  delete(dimz)
  
  arr!0 = "time"
  arr!1 = "lat"
  arr!2 = "lon"
  
  if (isatt(arr,"valid_range")) then    ; check to make sure data is in valid range. Reset to stay within the valid range if needed.
;     print("Data outside valid_range in "+zpath+", resetting data to be within valid_range")
     arr = where(arr.lt.arr@valid_range(0),arr@valid_range(0),arr)
     arr = where(arr.gt.arr@valid_range(1),arr@valid_range(1),arr)
  end if
  
  if (any(abs(arr).ge.1.e20)) then   ; check for inf values or values way out of range, reset to _FillValue.
     print("Values greater than 1.e20 or less than -1.e20 detected in "+zpath+", resetting to _FillValue")
     arr = where(abs(arr).ge.1.e20,arr@_FillValue,arr)
  end if
       
;  if (.not.iscoord(arr,"lat")) then
;     print("The data might have 2D latitudes. The functions.ncl coding is not set up to handle this, exiting from data_read_in.ncl")
;     print("zpath")
;     exit
;  end if
;  if (.not.iscoord(arr,"lon")) then
;     print("The data might have 2D longitudes. The functions.ncl coding is not set up to handle this, exiting from data_read_in.ncl")
;     print("zpath")
;  exit
;  end if
     
  if (yearS.lt.sydata.or.yearE.gt.eydata) then
     print("Requested "+yearS+"-"+yearE+" time span is outside the input file "+zpath+" time span of "+sydata+"-"+eydata+"")
     print("script: data_read_in.ncl")
     print("Setting array to all missing")
     delete(arr)
     arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12          
  else
     timeT =  yyyymm_time(sydata, eydata, "integer")
     time = timeT({sydata*100+smdata:eydata*100+emdata})  
     if (iscoord(arr,"time")) then
        delete(arr&time)
     end if
     dimz = dimsizes(arr)
     if (dimz(0).eq.dimsizes(time)) then
        arr&time = time
     else
        print("Possible mismatch detected between time specified in file name and file variables, setting array to missing")
        print("File = "+zpath)
        print("Read from file name: "+min(time)+"-"+max(time))        
        delete(arr)
        arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
        sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
        eydata = yearE     ; as data array is totally missing..  
        smdata = 1
        emdata = 12 
     end if
     delete(dimz)
     delete([/time,timeT/])
  end if
  delete([/sydata,smdata,eydata,emdata/])
    
;  printVarSummary(arr)
;  printVarSummary(arr({sy*100+1:ey*100+12},:,:))
  if (arr&lat(0).ge.0) then
     farr = arr({yearS*100+1:yearE*100+12},::-1,:)   ; flip the latitudes
  else
     farr = arr({yearS*100+1:yearE*100+12},:,:)
  end if
;  printVarSummary(farr)
  delete(arr)     
  
  mocheck = (/(yearS*100+1)-min(farr&time),(yearE*100+12) - max(farr&time)/)
  if (any(mocheck.ne.0)) then    ; previously: if (mod(dimsizes(farr&time),12).ne.0) then
     if (mocheck(0).ne.0) then
        print("First requested year is incomplete")
     end if
     if (mocheck(1).ne.0) then
        print("Last requested year is incomplete")
     end if
     print("Incomplete data year(s) requested for file "+zpath+", printing out time and creating blank array")
     print("Time requested: "+yearS+"-"+yearE)
     print(farr&time)
     delete(farr)
     farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
  end if
  delete(mocheck)
  
  if (farr&lon(0).lt.0) then
     farr = lonFlip(farr)     ; lon flip
  end if     
  if (min(farr&lon).lt.0.or.max(farr&lon).gt.360) then
     print(farr&lon)
     print("Model = "+names(gg)+", path = "+paths(gg))
     print("Fatal: Longitudes not in expected 0-360E range, creating blank array")
     delete(farr)
     farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
  end if
  
  if (vn.eq."TREFHT".or.vn.eq."TS") then      ; units check
     if (farr@units.eq."K".or.farr@units.eq."Kelvin") then   
        farr = farr-273.15
        farr@units = "C"
     end if
     if (farr@units.eq."degrees_C".or.farr@units.eq."degrees C") then
        farr@units = "C"
     end if
  end if
  if (vn.eq."PSL") then
     if (farr@units.eq."Pa".or.farr@units.eq."Pascals".or.farr@units.eq."Pascal") then
        farr = farr/100.
        farr@units = "hPa"
     end if
  end if
  if (vn.eq."PRECT") then    ; convert (if necessary) to mm/day
     if (farr@units.eq."m/s".or.farr@units.eq."m s-1") then
        farr = farr*86400000.
     end if
     if (farr@units.eq."kg m-2 s-1") then
        farr = farr*86400.
     end if
     if (farr@units.eq."cm".or.farr@units.eq."cm/month".or.farr@units.eq."mm".or.farr@units.eq."mm/month") then    
        yr = floattointeger(farr&time)/100
        mo = floattointeger(farr&time - (yr*100))
        days = days_in_month(yr,mo)
        do gg = 0,dimsizes(farr&time)-1
		         farr(gg,:,:) = (/ farr(gg,:,:) / days(gg) /)
	       end do
        if (farr@units.eq."cm".or.farr@units.eq."cm/month") then
           farr = farr*10.   ; convert from cm/day to mm/day
        end if        
     end if  
     farr@units = "mm/day" 
  end if    
  if (vn.eq."SNOWDP") then
     if (.not.isatt(farr,"is_all_missing")) then
        if (farr@units.ne."m".and.farr@units.ne."meters") then
           print("Warning: SNOWDP/snd units may not be in meters. listed units = "+farr@units)
        end if
     end if
  end if
  
  date = farr&time         ; switch time to be CF-conforming
  delete(farr&time)
  yyyy = date/100
  mm = date-(yyyy*100)
  days = (days_in_month(yyyy,mm))/2
  hms  = days
  hms = 0   ; hours, minutes, seconds all the same (=0)
  time = cd_inv_calendar(yyyy,mm,days,hms,hms,hms,"months since "+min(yyyy)+"-01-15 00:00:00",0)
  time@long_name = "Time"
  time@standard_name = "time"
  time!0 = "time"
  time&time = time
  farr&time = time
  delete([/time,yyyy,mm,days,hms,date/])
  return(farr)
end
;=================================================================================================
; read in ocean data from given files
;
; assign time coordinate variables, check for issues with the array, assign _FillValue (if needed)
; assign dimension names (for ease-of-use), check and modify units
;
; vname settings at top of this script can be modified if a different variable name is
; encountered. For instance, if a SSH data file has the SSH array named as "seasfc_height", one
; could add "seasfc_height" to the vname SSH coding as follows:
;   if (vn.eq."SSH") then
;      vname = (/"SSH","zos","seasfc_height"/)
;   end if
;
function data_read_in_ocean(zpath:string,vn:string,level:numeric,yearS:integer,yearE:integer)
; path for TS file(s), variable name, start year, and end year are read in.
local zpath,vn,cpath0,ta,tfiles,c,arr,farr,yearS,yearE,mocheck,fils_precc,fils_precl
begin
  if (vn.eq."SALT") then
     vname = (/"SALT","so"/)
  end if
  if (vn.eq."TEMP") then
     vname = (/"TEMP","thetao"/)
  end if
  if (vn.eq."SSH") then
     vname = (/"SSH","zos"/)
  end if

  if (ismissing(zpath) ) then
     print("File missing, creating blank array of data. View "+vn+" namelist for details.")
     arr = new((/(yearE-yearS+1)*12,90,180/),"float",1.e20)
     arr!0 = "time"    ; time coordinate variable assigned below
     arr!1 = "lat"
     arr&lat = fspan(-89,89,90)
     arr!2 = "lon"
     arr&lon = fspan(0,358,180)
     arr&lat@units = "degrees_north"
     arr&lon@units = "degrees_east"
     arr@units = ""
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12
     arr@is_all_missing = True
  else
     cpath0 = stringtochar(zpath) 
     ta = stringtochar("*")
     if (any(cpath0.eq.ta(0)).or.any(cpath0.eq."{")) then   ; check for "*" and "{" denoting multiple files
        tfiles = systemfunc("ls "+zpath+" 2> /dev/null")   ; /dev/null suppresses all standard error output
        c = addfiles(tfiles,"r")
        do ii=0,dimsizes(vname)-1
           if (isfilevar(c[0],vname(ii))) then
              if (vname(ii).eq."SALT".or.vname(ii).eq."TEMP") then
                 print("1")
                     
                 arr = c[0]->$vname(ii)$(:,{level*100.},:,:)  ; level in meters, CESM levels in cm
                 printVarSummary(arr)
                 exit
                 level@units = "m"
                 delete(arr&z_t)
                 print("2a")
                 arr!1 = "lev"
                 arr&lev = level
                 print("2")
              end if
              if (vname(ii).eq."so".or.vname(ii).eq."thetao") then
                 arr = c[:]->$vname(ii)$(:,{level},:,:)    ; level in meters
              end if
              if (vname(ii).eq."SSH".or.vname(ii).eq."zos") then
                 arr = c[:]->$vname(ii)$
              end if
              break
           end if
        end do   
        nfil = dimsizes(tfiles)
        cpathS = stringtochar(tfiles(0))      ; this section will work for PRECC/PRECL, as it will read the first
        cpathE = stringtochar(tfiles(nfil-1)) ; PRECC file and the last PRECL file.  
        ncharS = dimsizes(cpathS)
        ncharE = dimsizes(cpathE)
        sydata = stringtointeger(charactertostring(cpathS(ncharS-17:ncharS-14)))
        smdata = stringtointeger(charactertostring(cpathS(ncharS-13:ncharS-12)))
        eydata = stringtointeger(charactertostring(cpathE(ncharE-10:ncharE-7)))
        emdata = stringtointeger(charactertostring(cpathE(ncharE-6:ncharE-5)))      
        delete([/cpathS,cpathE,ncharS,ncharE,nfil/])
        delete(c)
     else
        c = addfile(zpath,"r")
        do i=0,dimsizes(vname)-1
           if (isfilevar(c,vname(i))) then
              if (vname(ii).eq."SALT".or.vname(ii).eq."TEMP") then
                 if (dimsizes(level).eq.1) then
                    arr = c->$vname(ii)$(:,{level*100.:level*100.},:,:)  ; level in meters, CESM levels in cm
                 else   
                    arr = c->$vname(ii)$(:,{level*100.},:,:)  ; level in meters, CESM levels in cm
                 end if
                 level@units = "m"
                 delete(arr&z_t)
                 arr!1 = "lev"
                 arr&lev = level
              end if
              if (vname(ii).eq."so".or.vname(ii).eq."thetao") then
                 if (dimsizes(level).eq.1) then
                    arr = c->$vname(ii)$(:,{level:level},:,:)    ; level in meters
                 else
                    arr = c->$vname(ii)$(:,{level},:,:)    ; level in meters
                 end if
              end if
              if (vname(ii).eq."SSH".or.vname(ii).eq."zos") then
                 arr = c->$vname(ii)$
              end if
              break
           end if
        end do
        cpath = stringtochar(zpath)     
        nchar = dimsizes(cpath)
        sydata = stringtointeger(charactertostring(cpath(nchar-17:nchar-14)))
        smdata = stringtointeger(charactertostring(cpath(nchar-13:nchar-12)))
        eydata = stringtointeger(charactertostring(cpath(nchar-10:nchar-7)))
        emdata = stringtointeger(charactertostring(cpath(nchar-6:nchar-5))) 
        delete([/cpath,nchar/])
;        delete(c)     
     end if
     delete([/ta,cpath0/])
  end if
  
  if (isshort(arr)) then
     arrT = short2flt(arr)
     delete(arr)
     arr = arrT
     delete(arrT)
  end if
  
  if (.not.isatt(arr,"_FillValue")) then     ; assign _FillValue if one is not present
     if (isatt(arr,"missing_value")) then
        arr@_FillValue = arr@missing_value
     else
        arr@_FillValue = default_fillvalue(typeof(arr)) 
     end if
  end if
  
  
  if (isatt(arr,"coordinates")) then
     strarr = str_splt(arr@coordinates)
     if (any(strarr.eq."TLONG")) then    ; CESM longitude 2D coordinate
        arr@lon2d = c->TLONG
     end if
     if (any(strarr.eq."TLAT")) then     ; CESM latitude 2D coordinate
        arr@lat2d = c->TLAT
     end if
     if (any(strarr.eq."lon")) then      ; IPCC longitude 2D coordinate
        arr@lon2d = c->lon
     end if
     if (any(strarr.eq."lat")) then      ; IPCC latitude 2D coordinate
        arr@lat2d = c->lat
     end if
  else
     print("2D coordinates for ocean data are not detected")
  end if

;  print("huh")

;  if (dimsizes(dimsizes(arr)).eq.3) then
;     arr!0 = "time"
;     arr!1 = "nlat"
;     arr!2 = "nlon"
;  end if
;  if (dimsizes(dimsizes(arr)).eq.4) then
     arr!0 = "time"
     arr!1 = "lev"
     arr!2 = "nlat"
     arr!3 = "nlon"
;  end if

;print("A")

  if (yearS.lt.sydata.or.yearE.gt.eydata) then
     print("Requested "+yearS+"-"+yearE+" time span is outside the input file "+zpath+" time span of "+sydata+"-"+eydata+"")
     print("script: data_read_in.ncl")
     print("Setting array to all missing")
     delete(arr)
     arr = new((/(yearE-yearS+1)*12,90,180/),"float",1.e20)
     arr!0 = "time"    ; time coordinate variable assigned below
     timeT =  yyyymm_time(yearS, yearE, "integer")
     time = timeT({yearS*100+smdata:yearE*100+emdata})  
     arr&time = time
     arr!1 = "lat"
     arr&lat = fspan(-89,89,90)
     arr!2 = "lon"
     arr&lon = fspan(0,358,180)
     arr&lat@units = "degrees_north"
     arr&lon@units = "degrees_east"
     arr@units = ""
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12
     arr@is_all_missing = True     
  else
     timeT =  yyyymm_time(sydata, eydata, "integer")
     time = timeT({sydata*100+smdata:eydata*100+emdata})  
     if (iscoord(arr,"time")) then
        delete(arr&time)
     end if
     arr&time = time
  end if
  delete([/timeT,time,sydata,smdata,eydata,emdata/])  
    
  printVarSummary(arr)
  printVarSummary(arr({sy*100+1:ey*100+12},:,:,:))
  farr = arr({yearS*100+1:yearE*100+12},:,:,:)  
  printVarSummary(farr)
  delete(arr)     
  
  mocheck = (/(yearS*100+1)-min(farr&time),(yearE*100+12) - max(farr&time)/)
  if (any(mocheck.ne.0)) then    ; previously: if (mod(dimsizes(farr&time),12).ne.0) then
     if (mocheck(0).ne.0) then
        print("First requested year is incomplete")
     end if
     if (mocheck(1).ne.0) then
        print("Last requested year is incomplete")
     end if
     print("Incomplete data year(s) requested for file "+zpath+", printing out time and exiting from data_read_in.ncl")
     print("Time requested: "+yearS+"-"+yearE)
     print(farr&time)
     exit
  end if
  delete(mocheck)

  if (vn.eq."SSH") then
     if (farr@units.eq."cm".or.farr@units.eq."CM") then
        farr = farr/100.
        farr@units = "m"
     end if
  end if
  if (vn.eq."TEMP") then      ; units check
     if (farr@units.eq."K".or.farr@units.eq."Kelvin") then   
        farr = farr-273.15
        farr@units = "C"
     end if
  end if  
  if (vn.eq."SALT") then
     if (farr@units.eq."gram/kilogram") then
        farr = (/ farr/1000. /)   ; convert to units of "1"
     end if
  end if
  
  date = farr&time         ; switch time to be CF-conforming
  delete(farr&time)
  yyyy = date/100
  mm = date-(yyyy*100)
  days = (days_in_month(yyyy,mm))/2
  hms  = days
  hms = 0   ; hours, minutes, seconds all the same (=0)
  time = cd_inv_calendar(yyyy,mm,days,hms,hms,hms,"months since "+min(yyyy)+"-01-15 00:00:00",0)
  time@long_name = "Time"
  time@standard_name = "time"
  time!0 = "time"
  time&time = time
  farr&time = time
  delete([/time,yyyy,mm,days,hms,date/])
  
  return(farr)
end
;=================================================================================================
; read in MOC ocean data from given files
;
; assign time coordinate variables, check for issues with the array, assign _FillValue (if needed)
; assign dimension names (for ease-of-use), check and modify units
;
function data_read_in_ocean_MOC(zpath:string,vn:string,yearS:integer,yearE:integer)
; path for MOC file(s), variable name, start year, and end year are read in.
local zpath,vn,cpath0,ta,tfiles,c,arr,farr,yearS,yearE,mocheck,dimC,lev
begin
  if (vn.eq."MOC") then
     vname = (/"MOC","msftmyz","stfmmc"/)
  end if

  if (ismissing(zpath) ) then
     print("File missing, creating blank array of data. View "+vn+" namelist for details.")
     arr = create_empty_array(yearS,yearE,1,12,"time_lev_lat")
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12
  else
     cpath0 = stringtochar(zpath)
  
     ta = stringtochar("*")
     if (any(cpath0.eq.ta(0)).or.any(cpath0.eq."{")) then   ; check for "*" and "{" denoting multiple files
        tfiles = systemfunc("ls "+zpath+" 2> /dev/null")   ; /dev/null suppresses all standard error output
        c = addfiles(tfiles,"r")
        do ii=0,dimsizes(vname)-1
           if (isfilevar(c[0],vname(ii))) then
              dimC = filevardimsizes(c[0],"MOC")
              if (vname(ii).eq."MOC") then    ; CCSM/CESM file     
                 if (dimC(2).ge.2) then
                    arr = dim_sum_n_Wrap(c[:]->$vname(ii)$(:,1,:,:,:),1)  ; select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region and sum over moc_comp
                 else
                    arr = c[:]->$vname(ii)$(:,1,0,:,:)        ; select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region and the only moc_comp dimension        
                 end if
              else                           ; CMIP file
                 arr = c[:]->$vname(ii)$(:,0,:,:)             ; CMIP file: 0th basin/region = atlantic_ocean (CMIP3) or atlantic_arctic_ocean (CMIP5)         
              end if   
              delete(dimC)  
              break
           end if
        end do     
        nfil = dimsizes(tfiles)
        cpathS = stringtochar(tfiles(0))      
        cpathE = stringtochar(tfiles(nfil-1)) 
        ncharS = dimsizes(cpathS)
        ncharE = dimsizes(cpathE)
        sydata = stringtointeger(charactertostring(cpathS(ncharS-17:ncharS-14)))
        smdata = stringtointeger(charactertostring(cpathS(ncharS-13:ncharS-12)))
        eydata = stringtointeger(charactertostring(cpathE(ncharE-10:ncharE-7)))
        emdata = stringtointeger(charactertostring(cpathE(ncharE-6:ncharE-5)))      
        delete([/cpathS,cpathE,ncharS,ncharE,nfil/])
;        delete(c)
     else
        c = addfile(zpath,"r")
        do i=0,dimsizes(vname)-1
           if (isfilevar(c,vname(i))) then
              dimC = filevardimsizes(c,"MOC")
              if (vname(i).eq."MOC") then    ; CCSM/CESM file     
                 if (dimC(2).ge.2) then
                    arr = dim_sum_n_Wrap(c->$vname(i)$(:,1,:,:,:),1)  ; select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region and sum over moc_comp
                 else
                    arr = c->$vname(i)$(:,1,0,:,:)        ; select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region         
                 end if
              else                           ; CMIP file
                 arr = c->$vname(i)$(:,0,:,:)             ; CMIP file: 0th basin/region = atlantic_ocean (CMIP3) or atlantic_arctic_ocean (CMIP5)         
              end if   
              delete(dimC)                       
              break
           end if                             
        end do                           
        cpath = stringtochar(zpath)     
        nchar = dimsizes(cpath)
        sydata = stringtointeger(charactertostring(cpath(nchar-17:nchar-14)))
        smdata = stringtointeger(charactertostring(cpath(nchar-13:nchar-12)))
        eydata = stringtointeger(charactertostring(cpath(nchar-10:nchar-7)))
        emdata = stringtointeger(charactertostring(cpath(nchar-6:nchar-5))) 
        delete([/cpath,nchar/])
;        delete(c)     
     end if
     delete([/ta,cpath0/])
  end if

  if (isvar("arr").eq.False) then
     print("Variable ("+vn+") not found. Examine input file "+zpath+". Creating empty array and continuing")
     arr = create_empty_array(yearS,yearE,1,12,"time_lev_lat")
  end if
  if (isshort(arr)) then
     arrT = short2flt(arr)
     delete(arr)
     arr = arrT
     delete(arrT)
  end if
  if (.not.isatt(arr,"_FillValue")) then     ; assign _FillValue if one is not present
     if (isatt(arr,"missing_value")) then
        arr@_FillValue = arr@missing_value
     else
        arr@_FillValue = default_fillvalue(typeof(arr)) 
     end if
  end if
  arr!0 = "time"
  arr!1 = "lev"
  arr!2 = "lat"
  
  if (arr&lev@units.eq."centimeters".or.arr&lev@units.eq."cm") then
     lev = arr&lev
     lev@units = "m"
     lev = lev/100.
     lev&lev = lev
     delete(arr&lev)
     arr&lev = lev
;    print("Level converted to m from cm")
;     printVarSummary(lev)
;     print(lev)
     delete(lev)
  end if
  
  if (arr&lev(2).lt.0) then   ; check for negative levels
     lev = arr&lev
     lev = lev*-1.
     if (any(lev.lt.0)) then
        print("Error detected in MOC level sign conversion")
        print(lev)
     end if
     lev@positive = "down"
     lev&lev = lev
     delete(arr&lev)
     arr&lev = lev
;     print("Levels converted from negative downwards to positive downwards")
;     printVarSummary(lev)
;     print(lev)
     delete(lev)
  end if
  
  
  if (isatt(arr,"valid_range")) then    ; check to make sure data is in valid range. Reset to stay within the valid range if needed.
;     print("Data outside valid_range in "+zpath+", resetting data to be within valid_range")
     arr = where(arr.lt.arr@valid_range(0),arr@valid_range(0),arr)
     arr = where(arr.gt.arr@valid_range(1),arr@valid_range(1),arr)
  end if
  
  if (any(abs(arr).ge.1.e20)) then   ; check for inf values or values way out of range, reset to _FillValue.
     print("Values greater than 1.e20 or less than -1.e20 detected in "+zpath+", resetting to _FillValue")
     arr = where(abs(arr).ge.1.e20,arr@_FillValue,arr)
  end if
     
  if (yearS.lt.sydata.or.yearE.gt.eydata) then
     print("Requested "+yearS+"-"+yearE+" time span is outside the input file "+zpath+" time span of "+sydata+"-"+eydata+"")
     print("script: data_read_in.ncl")
     print("Setting array to all missing")
     delete(arr)
     arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lev_lat")
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12          
  else
     timeT =  yyyymm_time(sydata, eydata, "integer")
     time = timeT({sydata*100+smdata:eydata*100+emdata})  
     if (iscoord(arr,"time")) then
        delete(arr&time)
     end if
     dimz = dimsizes(arr)
     if (dimz(0).eq.dimsizes(time)) then
        arr&time = time
     else
        print("Possible mismatch detected between time specified in file name and file variables, setting array to missing")
        print("File = "+zpath)
        print("Read from file name: "+min(time)+"-"+max(time))        
        delete(arr)
        arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lev_lat")
        sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
        eydata = yearE     ; as data array is totally missing..  
        smdata = 1
        emdata = 12 
     end if
     delete(dimz)
     delete([/time,timeT/])
  end if
  delete([/sydata,smdata,eydata,emdata/])
    
;  printVarSummary(arr)
;  printVarSummary(arr({sy*100+1:ey*100+12},:,:))
  if (arr&lat(0).ge.0) then
     farr = arr({yearS*100+1:yearE*100+12},:,::-1)   ; flip the latitudes
  else
     farr = arr({yearS*100+1:yearE*100+12},:,:)
  end if
;  printVarSummary(farr)
  delete(arr)     
  
  mocheck = (/(yearS*100+1)-min(farr&time),(yearE*100+12) - max(farr&time)/)
  if (any(mocheck.ne.0)) then    ; previously: if (mod(dimsizes(farr&time),12).ne.0) then
     if (mocheck(0).ne.0) then
        print("First requested year is incomplete")
     end if
     if (mocheck(1).ne.0) then
        print("Last requested year is incomplete")
     end if
     print("Incomplete data year(s) requested for file "+zpath+", printing out time and creating blank array")
     print("Time requested: "+yearS+"-"+yearE)
     print(farr&time)
     delete(farr)
     farr = create_empty_array(yearS,yearE,1,12,"time_lev_lat")
  end if
  delete(mocheck)
  
  ; check units for MOC array. CMIP5 = "kg s-1"  CMIP3 = "m3 s-1" CCSM3 = "Sverdrups" CCSM4 = "Sverdrups"

  if (farr@units.eq."Sverdrups") then
     farr@units = "Sv"
  end if
  if (farr@units.eq."kg s-1".or.farr@units.eq."KG S-1".or.farr@units.eq."kg/s".or.farr@units.eq."KG/S") then     ; 1 Sv = 1.e9 kg/s
     farr = (/ farr/1.e9 /)
     farr@units = "Sv"
  end if
  if (farr@units.eq."m3 s-1".or.farr@units.eq."M3 S-1".or.farr@units.eq."m3/s".or.farr@units.eq."M3/S") then     ; 1 Sv = 1.e6 m3/s
     farr = (/ farr/1.e6 /)
     farr@units = "Sv"
  end if
  
;  printVarSummary(farr)
  
  date = farr&time         ; switch time to be CF-conforming
  delete(farr&time)
  yyyy = date/100
  mm = date-(yyyy*100)
  days = (days_in_month(yyyy,mm))/2
  hms  = days
  hms = 0   ; hours, minutes, seconds all the same (=0)
  time = cd_inv_calendar(yyyy,mm,days,hms,hms,hms,"months since "+min(yyyy)+"-01-15 00:00:00",0)
  time@long_name = "Time"
  time@standard_name = "time"
  time!0 = "time"
  time&time = time
  farr&time = time
  delete([/time,yyyy,mm,days,hms,date/])
  return(farr)
end

;=================================================================================================
;  alters the formatting of the Y-axis
;
;  not currently used
;
function y_axis_check(temparr:numeric,tempres:logical)
local temparr,tempres,minval,maxval
begin
  minval = min(temparr)
  maxval = max(temparr)
  if (minval.gt.-1.and.minval.lt.0.and.maxval.lt.1.and.maxval.gt.0) then
     tempres@tmYLFormat = "0@;*.2f"
  else
     tempres@tmYLFormat = "0@*+^sg"
  end if
  return(tempres)
  delete([/tempres,temparr,minval,maxval/])
end
;=================================================================================================
; Check that the user-specified climatological period is within the time range of the data
;
procedure check_custom_climo(mn:string,startyear:numeric,endyear:numeric,climo_startyear:numeric,climo_endyear:numeric)
local startyear,endyear,climo_startyear,climo_endyear,mn
begin  
  do gg = 0,dimsizes(startyear)-1
     if (climo_startyear.ge.startyear(gg).and.climo_endyear.le.endyear(gg)) then
     else
        print("check_custom_climo: Warning! Beginning and/or ending of climatological period is outside time range of data.")
        print("Dataset: "+mn+", years = "+startyear(gg)+":"+endyear(gg)+", set climatological period = "+climo_startyear+":"+climo_endyear)
        print("The diagnostics package will proceed, but one or more dataset(s) will not have the full climatological period removed and/or tha package may fail with the following message: fatal:NclOneDValGetRangeIndex: start coordinate index out of range.")
     end if
  end do
end
;==================================================================================================
; In version 6.2.1 the behavior of isfilepresent switched, where only files readable by NCL return
; True. Previously if a file (or directory) simply existed, isfilepresent would return True. A new
; function has been created in v6.2.1, fileexists, that acts like the previous version of isfilepresent
; did. To compensate for this, check the NCL version number, and use isfilepresent/fileexists when 
; appropriate.
;
function isfilepresent2(fdpres:string)
local nclver, num0, num1, ra
begin
  nclver = stringtochar(get_ncl_version())
  
  num0 = toint(tostring(nclver(0)))
  num1 = toint(tostring(nclver(2)))
  num2 = toint(tostring(nclver(4)))
  if (num0.le.5) then
     ra = isfilepresent(fdpres)
  end if
  if (num0.eq.6) then
     if (num1.le.1) then
        ra = isfilepresent(fdpres)
     end if
     if (num1.eq.2) then
        if (num2.eq.0) then
           ra = isfilepresent(fdpres)
        else
           ra = fileexists(fdpres)
        end if
     end if
     if (num1.ge.3) then
        ra = fileexists(fdpres)
     end if
  end if  
  if (num0.ge.7) then
     ra = fileexists(fdpres)
  end if
  return(ra)
  delete([/nclver,num0,num1,ra/])
end
;=================================================================================================
;
function table_link_setup(ipath:string,iname:string,ltxt:string)
; image name, along with link text
local ipath, iname, ltxt, otxt, quote
begin
  quote = str_get_dq()
  if (isfilepresent2(ipath+iname)) then
     otxt = "<a href="+quote+iname+quote+">"+ltxt+"</a>" 
  else
     otxt = ltxt
  end if
  return(otxt)
  delete([/ipath,iname,ltxt,otxt,quote/])
end
