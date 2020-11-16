;=================================================================================================
; create blank array for use when something may be/is wrong.
;
undef("create_empty_array")
function create_empty_array(yS:numeric,yE:numeric,mS:numeric,mE:numeric,opttype:string)
local yS,yE,mS,mE,opttype
begin
  if (ismissing(yS).or.ismissing(yE)) then
     yS = 1
     yE = 50
  end if
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
undef("data_read_in")
function data_read_in(zpath:string,vn:string,yearS:integer,yearE:integer)
; path for TS file(s), variable name, start year, and end year are read in.
local zpath,vn,cpath0,ta,tfiles,c,arr,farr,yearS,yearE,mocheck,fils_precc,fils_precl
begin
  if (vn.eq."TS") then
     vname = (/"TS","ts","sst","t_surf","skt"/)
  end if
  if (vn.eq."PSL") then
     vname = (/"PSL","psl","slp","SLP","prmsl","msl","slp_dyn"/)
  end if
  if (vn.eq."TREFHT") then
     vname = (/"TREFHT","tas","temp","air","temperature_anomaly","temperature","t2m","t_ref","T2","tempanomaly"/)
  end if
  if (vn.eq."PRECT") then
     vname = (/"PRECC","PRECL","PRECT","pr","PPT","ppt","p","P","precip","PRECIP","tp","prcp","prate"/)
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
  if (dimsizes(dimz).le.2) then
     print("Possible curvilinear (or unstructured) grid detected. The CVDP cannot analyze curvilinear data. Please regrid to a rectilinear grid for inclusion in CVDP comparisons.") 
     print("Input file: "+zpath)
     print("Setting array to all missing")
     delete(arr)
     arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
     sydata = yearS     ; assign these variables based on yearS/yearE provided in namelist. Doesn't matter 
     eydata = yearE     ; as data array is totally missing..  
     smdata = 1
     emdata = 12 
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
     print("path = "+zpath)
     print("Fatal: Longitudes not in expected 0-360E range, creating blank array")
     delete(farr)
     farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
  end if
  
  if (vn.eq."TREFHT".or.vn.eq."TS") then      ; units check
     if (farr@units.eq."K".or.farr@units.eq."Kelvin".or.farr@units.eq."deg_k".or.farr@units.eq."deg_K") then   
        if (max(farr).ge.100) then    ; data sets can be anomalies with units of K, so check for range before subtracting
           farr = farr-273.15
        end if
        farr@units = "C"
     end if
     if (farr@units.eq."degrees_C".or.farr@units.eq."degrees C".or.farr@units.eq."degree_C".or.farr@units.eq."degree C") then
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
     if (farr@units.eq."kg m-2 s-1".or.farr@units.eq."kg/m2/s".or.farr@units.eq."kg/m^2/s".or.farr@units.eq."kg/(s*m2)".or.farr@units.eq."mm/s") then
        farr = farr*86400.
     end if
     if (farr@units.eq."m".or.farr@units.eq."m/month".or.farr@units.eq."cm".or.farr@units.eq."cm/month".or.farr@units.eq."mm".or.farr@units.eq."mm/month") then    
        yr = toint(farr&time)/100
        mo = toint(farr&time - (yr*100))
        days = days_in_month(yr,mo)
        do gg = 0,dimsizes(farr&time)-1
		         farr(gg,:,:) = (/ farr(gg,:,:) / days(gg) /)
	       end do
        if (farr@units.eq."cm".or.farr@units.eq."cm/month") then
           farr = farr*10.   ; convert from cm/day to mm/day
        end if        
        if (farr@units.eq."m".or.farr@units.eq."m/month") then
           farr = farr*1000.   ; convert from m/day to mm/day
        end if   
     end if  
     if (farr@units.eq."m/day".or.farr@units.eq."m day-1") then
        farr = farr*1000.
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
  time@actual_range = (/min(time),max(time)/)
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
undef("data_read_in_ocean_MOC")
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

  if (isatt(arr,"coordinates")) then
     delete(arr@coordinates)
  end if
  
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
; read in ice data from given files
;
; assign time coordinate variables, check for issues with the array, assign _FillValue (if needed)
; assign dimension names (for ease-of-use), check and modify units
;
undef("data_read_in_ice")
function data_read_in_ice(zpath:string,vn:string,yearS:integer,yearE:integer)
; path for ice file(s), variable name, start year, and end year are read in.
local zpath,vn,cpath0,ta,tfiles,c,arr,farr,yearS,yearE,mocheck,dimC,lev
begin
  if (vn.eq."aice_nh") then
     vname = (/"aice_nh","aice","sic","SIC","CN","ice","icec"/)
  end if
  if (vn.eq."aice_sh") then
     vname = (/"aice_sh","aice","sic","SIC","CN","ice","icec"/)
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
     tfiles = systemfunc("ls "+zpath+" 2> /dev/null")   ; /dev/null suppresses all standard error output
     c = addfiles(tfiles,"r")
     do ii=0,dimsizes(vname)-1
        if (isfilevar(c[0],vname(ii))) then
           if (vname(ii).eq."aice_nh".or.vname(ii).eq."aice_sh".or.vname(ii).eq."aice") then    ; CCSM/CESM file     
              arr = c[:]->$vname(ii)$
              if (isatt(arr,"coordinates")) then
                 strarr = str_split(arr@coordinates," ")
                 if (any(strarr.eq."TLON")) then    ; CESM longitude 2D coordinate
                    dimZ = dimsizes(c[0]->TLON)
                    if (dimsizes(dimZ).eq.3) then
                       arr@lon2d = c[0]->TLON(0,:,:)
                    else
                       arr@lon2d = c[0]->TLON
                    end if
                    delete(dimZ)
                 end if
                 if (any(strarr.eq."TLAT")) then     ; CESM latitude 2D coordinate
                    dimZ = dimsizes(c[0]->TLAT)
                    if (dimsizes(dimZ).eq.3) then
                       arr@lat2d = c[0]->TLAT(0,:,:)
                    else
                       arr@lat2d = c[0]->TLAT
                    end if
                    delete(dimZ)
                 end if
                 delete(strarr)
;              else
;                 print("2D coordinates for ice data are not detected")
              end if
              if (isatt(arr,"cell_measures").and.isfilevar(c[0],"tarea")) then    ; if an attribute named cell_measures exists, and tarea is on file(0)
                 if (arr@cell_measures.eq."area: tarea") then
                    arr@area = totype(c[0]->tarea,typeof(arr))      ; in units of m^2
                 end if
              end if
           else                           ; CMIP or other file
              if (vname(ii).eq."CN") then     ; GFDL file
                 arrT = c[:]->$vname(ii)$
                 arr = dim_sum_n_Wrap(arrT,1)
                 delete(arrT)
                 arr = where(arr.ge.1,1,arr) ; optional
              else
                 arr = c[:]->$vname(ii)$
              end if
              if (isatt(arr,"coordinates")) then
                 strarr = str_split(arr@coordinates," ")
                 if (any(strarr.eq."lon")) then      ; IPCC longitude 2D coordinate
                    arr@lon2d = c[0]->lon
                 end if
                 if (any(strarr.eq."lat")) then      ; IPCC latitude 2D coordinate
                    arr@lat2d = c[0]->lat
                 end if
                 if (any(strarr.eq."longitude")) then      ; NSIDC longitude 2D coordinate
                    arr@lon2d = c[0]->longitude
                 end if
                 if (any(strarr.eq."latitude")) then      ; NSIDC latitude 2D coordinate
                    arr@lat2d = c[0]->latitude
                 end if
                 delete(strarr)
;              else
;                 print("2D coordinates for ice data are not detected")
              end if  
              dir_name = str_split(tfiles(0),"/")
              if (dimsizes(dir_name).ge.8) then
                 dir_name_new = "/"+str_join(dir_name(:4),"/")+"/fx/areacello/"+dir_name(7)+"/r0i0p0/*.nc"
                 ufile = systemfunc("ls "+dir_name_new+" 2> /dev/null")   ; /dev/null suppresses all standard error output
                 delete(dir_name_new)
              else
                 ufile = new(1,string)
              end if
              if (.not.ismissing(ufile)) then
                 d = addfile(ufile,"r")
                 arr@area = totype(d->areacello,typeof(arr))
                 dimQ = dimsizes(arr)
                 if (dimsizes(ndtooned(arr@area)).ne.(dimQ(1)*dimQ(2))) then    ; the dimension sizes of areacello 
                    delete(arr@area)                                  ; do not match sizes of area j,i dimensions 
                 end if
                 delete(dimQ)
              end if

              if (isfilevar(c[0],"AREA")) then    ; check to see if there is an AREA array present and if so use it
                 areaT = c[0]->AREA
                 if (areaT@units.eq."km^2") then
                    area_unit_km2_to_m2 = True
                    areaT = areaT*1000000.
                    areaT@units = "m^2"
                 end if
                 areaT@_FillValue = 1.e20
                 arr@area = totype(areaT,typeof(arr))
                 if (isatt(areaT,"pole_hole_area")) then    ; format of ystart, yend, value, ystart, yend, value
                    if (isvar("area_unit_km2_to_m2")) then
                       extra_area = tofloat(areaT@pole_hole_area)
                       extra_area(2::3) = extra_area(2::3)*1000000.   ; convert pole hole area from km^2->m^2
                       arr@pole_hole_area = totype(extra_area,typeof(arr))
                       delete(extra_area)
                    else
                       arr@pole_hole_area = totype(areaT@pole_hole_area,typeof(arr))
                    end if
                 end if
                 delete(areaT)
              end if
              delete([/dir_name,ufile/]) 
           end if   
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
     delete([/cpathS,cpathE,ncharS,ncharE,nfil,ta,cpath0/])
  end if

  if (isvar("arr").eq.False) then
     print("Variable ("+vn+") not found. Examine input file "+zpath+". Creating empty array and continuing")
     arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
  end if

  if (.not.isatt(arr,"area")) then   ; calculate grid cell areas manually (not implemented)
;     print("Grid cell areas not found.")
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
  arr!1 = "j"
  arr!2 = "i"

  if (.not.isatt(arr,"lat2d")) then     ; if latitudes are 1D, make sure latitudes run from south to north +
     if (arr&j(0).ge.0) then            ; calculate area of 1D lat/lon arrays
        tarr = arr(:,::-1,:)
        delete(arr)
        arr = tarr
        delete(tarr)
     end if

     if (min(arr&i).ge.0.and.max(arr&i).le.360) then
        fctr = 111120   ; how many meters per degree of latitude (approximate)
	       pi=4.*atan(1.0)
        rad=(pi/180.)
        lat = tofloat(arr&j)
	       dimlat = dimsizes(lat)
	       latr = new(dimlat,typeof(lat))
	       do gg = 0,dimlat-1
	          if (gg.eq.0) then
		            latr(gg) = abs(-90-(lat(1)+lat(0))/2.)
		         end if
	          if (gg.ge.1.and.gg.lt.dimlat-1) then
	             latr(gg) = abs((lat(gg-1)+lat(gg))/2. - (lat(gg)+lat(gg+1))/2.)
		         end if
		         if (gg.eq.dimlat-1) then
		            latr(gg) = abs(90 - (lat(dimlat-2)+lat(dimlat-1))/2.)
		         end if
	       end do
        lon = tofloat(arr&i)
	       dimlon = dimsizes(lon)
	       lonr = new(dimlon,typeof(lon))
        do gg = 0,dimlon-1
	          if (gg.eq.0) then
		            lonr(gg) = abs( (lon(1)+lon(0))/2. - (((lon(dimlon-1)+(lon(0)+360))/2.)-360) )
		         end if
		         if (gg.ge.1.and.gg.lt.dimlon-1) then
			       	   lonr(gg) = abs((lon(gg)+lon(gg+1))/2. - (lon(gg-1)+lon(gg))/2.)
		         end if
		         if (gg.eq.dimlon-1) then
		            lonr(gg) = abs(((lon(dimlon-1)+(lon(0)+360))/2.) - (lon(gg-1)+lon(gg))/2.)
		         end if		
	       end do
        area = tofloat(arr(0,:,:))
	       area = area@_FillValue
	       area@long_name = "Area of grid box"
	       area@units = "m2"

;        printVarSummary(area)
	       do ff = 0,dimlat-1
	          do gg = 0,dimlon-1
		            area(ff,gg) = (/ (fctr*latr(ff))*(cos(rad*lat(ff))*lonr(gg)*fctr) /)     ; cosine weighting
		         end do
	       end do
;        print("Total area = "+sum(area))
        arr@area = totype(area,typeof(arr))
        delete([/lat,lon,latr,lonr,area,fctr,pi,rad,dimlat,dimlon/])		
     end if
  end if

  if (.not.isatt(arr,"is_all_missing")) then    ; erase data in hemisphere not specified via vn
     if (isatt(arr,"lat2d")) then
        tlat2 = conform(arr,arr@lat2d,(/1,2/))
        tlon2 = conform(arr,arr@lon2d,(/1,2/))
        if (vn.eq."aice_nh") then
           arr = where(tlat2.ge.0,arr,arr@_FillValue)
        end if
        if (vn.eq."aice_sh") then
           arr = where(tlat2.lt.0,arr,arr@_FillValue)
        end if
        delete([/tlat2,tlon2/])
     else
        if (vn.eq."aice_nh") then
           arr(:,{:-1.},:) = arr@_FillValue
        end if
        if (vn.eq."aice_sh") then
           arr(:,{0:},:) = arr@_FillValue
        end if
     end if
  end if
     
  if (yearS.lt.sydata.or.yearE.gt.eydata) then
     print("Requested "+yearS+"-"+yearE+" time span is outside the input file "+zpath+" time span of "+sydata+"-"+eydata+"")
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
  farr = arr({yearS*100+1:yearE*100+12},:,:)  
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
     print("From file: Times present from "+min(farr&time)+"-"+max(farr&time))
     delete(farr)
     farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
  end if
  delete(mocheck)

  if (farr@units.eq."0-1".or.farr@units.eq."1") then    ; GFDL units, NSIDC units
     farr = (/ farr*100. /)
     farr@units = "%"
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
;  alters the formatting of the Y-axis
;
;  not currently used
;
undef("y_axis_check")
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
undef("check_custom_climo")
procedure check_custom_climo(mn:string,startyear:numeric,endyear:numeric,climo_startyear:numeric,climo_endyear:numeric)
local startyear,endyear,climo_startyear,climo_endyear,mn
begin  
  do gg = 0,dimsizes(startyear)-1
     if (climo_startyear.ge.0) then    ; exact years specified for climatological period
        if (climo_startyear.ge.startyear(gg).and.climo_endyear.le.endyear(gg)) then
        else
           print("check_custom_climo: Warning! Beginning and/or ending of climatological period is outside time range of data.")
           print("Dataset: "+mn+", years = "+startyear(gg)+":"+endyear(gg)+", set climatological period = "+climo_startyear+":"+climo_endyear)
           print("The diagnostics package will proceed, but one or more dataset(s) will not have the full climatological period removed and/or the package may fail with the following message: fatal:NclOneDValGetRangeIndex: start coordinate index out of range.")
        end if
     else                              ; relative years specified for climatological period
        if ((endyear(gg)-startyear(gg)+1).lt.(climo_endyear-climo_startyear+1)) then
           print("check_custom_climo: Warning! Beginning and/or ending of climatological period is outside time range of data.")
           print("Dataset: "+mn+", years = "+startyear(gg)+":"+endyear(gg)+", set climatological period = "+(endyear(gg)+climo_startyear)+":"+(endyear(gg)+climo_endyear))
           print("The diagnostics package will proceed, but one or more dataset(s) will not have the full climatological period removed and/or the package may fail with the following message: fatal:NclOneDValGetRangeIndex: start coordinate index out of range.")
        end if
        if (abs(climo_startyear).ge.(endyear(gg)-startyear+1)) then
           print("check_custom_climo: Warning! Dataset: "+mn+", climatology start year "+(endyear(gg)+climo_startyear)+" is outside of analysis time period ("+startyear(gg)+":"+endyear(gg)+"), exiting script.")
           exit
        end if
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
undef("isfilepresent2")
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
undef("table_link_setup")
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
;=================================================================================================
undef ("gsn_panel2")
procedure gsn_panel2(wksp:graphic,plotp:graphic,lpl:numeric,panelres:logical)
; checks to make sure at least one image is present in plot before paneling,
; thereby eliminating error message: 
; Error: gsn_panel: all of the plots passed to gsn_panel appear to be invalid

local wksp, plotp, lpl, panelres
begin
  if (.not.all(ismissing(plotp))) then
     gsn_panel(wksp,plotp,lpl,panelres)
  end if  
end
;=================================================================================================
undef ("eofunc_north2")
function eofunc_north2(eval[*]:numeric, N[1]:integer, prinfo[1]:logical)
;
; North, G.R. et al (1982): Sampling Errors in the Estimation of Empirical Orthogonal Functions.
; Mon. Wea. Rev., 110, 699–706.
; doi: http://dx.doi.org/10.1175/1520-0493(1982)110<0699:SEITEO>2.0.CO;2
;
; Usage after 'eofunc'. Here ntim was used,  
;             prinfo = True
;             sig    = eval_north(eof@eval, ntim, prinfo)
;
; Copied directly from v6.3.0 contributed.ncl for use in the package regardless of NCL version.
;
local neval, dlam, low, high, sig, n
begin
  neval   = dimsizes(eval)
  if (neval.eq.1)
      print("eofunc_north: neval=1, no testing can be performed")
      sig = True 
      sig@long_name = "EOF separation is not testable N=1"
      sig@N         =  N
      return(sig)
  end if

  dlam    = eval * sqrt(2.0/N)   ; eq 24
  low     = eval-dlam
  high    = eval+dlam

  sig     = new(dimsizes(eval), logical)
  sig = False  ; default is not significantly separated

; first and last eigenvalues are special cales

  if (eval(0).gt.high(1)) then
      sig(0) = True
  end if
  if (eval(neval-1).lt.low(neval-2)) then
      sig(neval-1) = True
  end if

; loop over other eignvalues

  if (N.gt.2) then
      do n=1,neval-2
         if (eval(n).lt.low(n-1) .and. eval(n).gt.high(n+1)) then
             sig(n) = True
         end if
      end do
  end if

  if (prinfo) then
      print(dlam+"   "+low+"   "+eval+"   "+high+"  "+sig)
  end if

  sig@long_name = "EOF separation"
  sig@N         =  N
  return(sig)
end
;=================================================================================================
; Standardize and set attributes for array. Remove NCL-added and superfluous attributes, 
; set missing_value equal to _FillValue, provide options to set long_name, units, and comment_cvdp 
; attributes. For last three inputs "" means leave as set and "delete" means delete the attribute.
; This function will be used immediately prior to an array being written to a netCDF file.
;
undef("set_varAtts")
function set_varAtts(zarr:numeric,loname:string,uni:string,com_cvdp:string)

local zarr,loname,uni,com_cvdp
begin
  if (isatt(zarr,"anomaly_op_ncl")) then
     delete(zarr@anomaly_op_ncl)
  end if
  if (isatt(zarr,"average_op_ncl")) then
     delete(zarr@average_op_ncl)
  end if
  if (isatt(zarr,"cell_measures")) then
     delete(zarr@cell_measures)
  end if
  if (isatt(zarr,"cell_methods")) then
     delete(zarr@cell_methods)
  end if
  if (isatt(zarr,"lonFlip")) then
     delete(zarr@lonFlip)
  end if
  if (isatt(zarr,"runave_op_ncl")) then
     delete(zarr@runave_op_ncl)
  end if
  if (isatt(zarr,"stddev_op_ncl")) then
     delete(zarr@stddev_op_ncl)
  end if
  if (isatt(zarr,"sum_op_ncl")) then
     delete(zarr@sum_op_ncl)
  end if
  if (isatt(zarr,"time")) then
     delete(zarr@time)
  end if
  if (isatt(zarr,"wgt_areaave_op_ncl")) then
     delete(zarr@wgt_areaave_op_ncl)
  end if

  if (isatt(zarr,"_FillValue")) then      ; set missing_value = _FillValue
     zarr@missing_value = zarr@_FillValue
  end if
  
  if (loname.eq."delete") then
     delete(zarr@long_name)
  else
     if (loname.ne."") then    ; "" = leave as is
        zarr@long_name = loname
     end if
  end if
  if (uni.eq."delete") then
     delete(zarr@units)
  else
     if (uni.ne."") then    ; "" = leave as is
        zarr@units = uni
     end if
  end if
  if (com_cvdp.eq."delete") then
     delete(zarr@comment_cvdp)
  else
     if (com_cvdp.ne."") then    ; "" = leave as is
        zarr@comment_cvdp = com_cvdp
     end if
  end if

  return(zarr)
end

