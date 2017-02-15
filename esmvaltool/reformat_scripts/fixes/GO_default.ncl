
;; This script is called by ./reformat/reformatGO_main.ncl
;; For (switch .eq. "ocean-latlon") this script should write to fout:
;;  * variable(time,i,j) 
;;  * lat(i,j)
;;  * lon(i,j)
;;  * area(i,j)   -> grid cell areas in sqkm
;; 
;; It is only active specified in array "rewrite_variables" !!!!!!!!!!!!!!!!!!!!!!!!!!!
;; 
;; Note: Coordinates for variable(time,lat,lon) are set in attribute.ncl.
;;       Thus unusual original coordinate names must be defined there.
;; Here is taken care of renaming all coordinates according to the above 
;; specifications. This is supposed to work automatically. 
;; However, scaling of values (e.g. 0<sic<1, areacello) must be implemented 
;; here on a per-model basis (although something is auto already).
;; Creates 2D lat & 2D lon for Gaussian grids.
;;
;; Contact: K. Gottschaldt (klaus-dirk.gottschaldt@dlr.de)

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"

undef("fix_default")
function fix_default(dataX, switch)
local dNames,fn_out,lat_name,lon_name,lat,lon,lat_index,lon_index,fout,\
      areafile,g,fx,area,script_name,i,j,lat2D,lon2D
begin
    
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output(" In fix_default from GO_default.ncl ", verbosity, 3)
    script_name = getenv("ESMValTool_project")+"_"+getenv("ESMValTool_model")+".ncl"
 
    if (switch .eq. "ocean-latlon" ) then
      rewrite_variables = (/"chl","fgco2","spco2","sosstsst","o2","talk"/)

      variable   = getenv("ESMValTool_variable")

      info_output("variable: " + variabl, verbosity, 3)
        if(inlist(variable, rewrite_variables)) then  
            ;; Fetch necessary variables
	    fin = dataX(0)
	    fout = dataX(1)
            lat = fin->lat   ;get 2D lat from file
            lon = fin->lon   ;get 2D lon from file
            var = fout->$variable$
            dNames = getfilevardims(fout, variable)
            lat_index = ind(dNames .eq. "lat")  ; var(time,lat,lon)
            lon_index = ind(dNames .eq. "lon")

            info_output("lat_index= " + lat_index, verbosity, 3)
           
           delete(fout)


            if ( dimsizes(dimsizes(var)) .eq. 4 ) then
              depth_name = asciiread("./interface_data/depth_name.txt", 1, "string")
            end if
            fn_out = asciiread("./interface_data/fn_out.txt", 1, "string")
            lat_name = asciiread("./interface_data/lat_name.txt", 1, "string")
            lon_name = asciiread("./interface_data/lon_name.txt", 1, "string")        
            if (isfilevar(fin,lat_name)) then  ;get original coordinates (not always in file!)
               info_output("lat_name from file " + lat_name, verbosity, 3)
              i = fin->$lon_name$  
              j = fin->$lat_name$
            else     ;fall back option
              info_output("lat_name from variable " + var&lat, verbosity, 3)
              i = var&lon
              j = var&lat
            end if 

            ;; fix coordinates of variable 
            var!lat_index = lat_name  
            var!lon_index = lon_name  
            delete(var&$lat_name$)
            delete(var&$lon_name$)
            var!lat_index = "j"
            var!lon_index = "i" 
            var&i = i
            var&j = j

            ;; get grid cell areas
            areafile = getenv("ESMValTool_areafile_path")  ;;defined in main nml (e.g. for project CMIP5_gridfile)
            info_output("areafile: " + areafile, verbosity, 3)
            
            info_output("reading areafile = "+areafile,verbosity,2)
            
            g = addfile(areafile,"r")
            delete(areafile)
            area = g->areacello

            ;; update dNames for 2D
            delete(dNames)
            dNames = getfilevardims(g, "areacello")
            lat_index = ind(dNames .eq. lat_name)  ; assumes that areacello has the same
            lon_index = ind(dNames .eq. lon_name)  ;  original coordinate names as variable 

            ;; fix grid cell areas
            area!lon_index = "i"
            area!lat_index = "j"
            area&i = i
            area&j = j

            ;; make lat & lon 2D, if needed 
            if (dimsizes(dimsizes(lat)).eq. 1) then   ;check rank       
                lat2D = conform_dims(dimsizes(area),lat,0) ;make 2D to be compatible with
                lon2D = conform_dims(dimsizes(area),lon,1) ; irregular grid handling
                delete([/lat,lon/])
                lat = lat2D
                lon = lon2D
            end if

            ;; fix coordinates of lat
            lat!lat_index = "j"
            lat!lon_index = "i"
            lat&i = i
            lat&j = j

            ;; fix coordinates of lon
            lon!lat_index = "j"
            lon!lon_index = "i"
            lon&i = i
            lon&j = j

            ;; scale whatever is necessary
            ;; Adjust on a per-model basis HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (variable .eq. "sic") then
              if (max(var).gt.95. .and. max(var).lt.105.) then 
                var = var/100.  ;sic is expected to be within [0,1]
                info_output("Auto scaling "+var+"/100 in "+script_name,verbosity,3)
              else if (max(var).gt.0.95 .and. max(var).lt.1.05)
                info_output("No scaling for "+var+" in "+script_name,verbosity,3) 
              else 
                print("Warning: check scaling of sic in "+script_name) 
              end if   
              end if    
            end if
            if (area@units .eq. "m2") then
              area = area / 1000000.  ;convert from sqm to sqkm
              info_output("Auto scaling areacello/1e6 in "+script_name,verbosity,3)
            else if (area@units .eq. "m-2") then
              area = area / 1000000.  ;convert from sqm to sqkm for GFDL model family
              info_output("Auto scaling areacello/1e6 in "+script_name,verbosity,3)
            else if (area@units .eq. "m 2") then
              area = area / 1000000.  ;convert from sqm to sqkm for GFDL model family
              info_output("Auto scaling areacello/1e6 in "+script_name,verbosity,3)
            end if
            end if  
            end if
            ;; END of adjustment section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ;; write to fout
            system("rm -f " + fn_out)            
            fout = addfile(fn_out, "c")
            fout->lat = lat(i|:,j|:)          
            fout->lon = lon(i|:,j|:)            
            if (dimsizes(dimsizes(var)) .eq. 4 ) then
              fout->$variable$ = var(time|:,depth|:,i|:,j|:)
            else 
              fout->$variable$ = var(time|:,i|:,j|:)          
            end if
            fout->area = area(i|:,j|:) 
            info_output("Written area to " + fn_out, verbosity, 3)
           
        else
            info_output("No operation for ocean_latlon in "+script_name, verbosity, 4)
        end if
    else
        info_output(" switch '" + switch + "' in fix_data-routine invalid",verbosity,1)
;        status_exit(1)
    end if
    return(dataX)
end
