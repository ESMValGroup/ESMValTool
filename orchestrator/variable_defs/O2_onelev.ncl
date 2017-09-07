;
;  Requires: O2:TO2M
;
variable_info = True
variable_info@derived = True
variable_info@long_name="o2"
variable_info@units = "mmol/m^3"
variable_info@MyDiag_title = "Oxygen"

load "interface_scripts/data_handling.ncl"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_number [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'temp_ncl_code/filesinfo'-file
;;    variable - logical with relevant variable as string attribute
;;    field_number  - string with field number classification
;; Return value:
;;    data_new - logical
local tmp, dum, dimension 
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering O2_onelev.ncl", verbosity, 2)

  data_new = True
  tmp = read_data(index,"O2","TO2M")
   
   dum = extract_data(index,tmp,-1,0,0)
;;   if (dum@units .eq. "mol/m3" .or. dum@units .eq. "mol m-3" ) then
;;        dum = dum * 1000.
;;        dum@units = "mmol/m^3"
;;   else if (dum@units .eq. "mol/L" .or. dum@units .eq. "mol L-1" ) then
;;         dum = dum * 1000000.
;;         dum@units = "mmol/m^3" 
      
;;   else if ( dum@units .eq. "mmol/m3" .or. dum@units .eq. "mmol m-3" ) then
            dum = dum  
;;   else
;;            print("fatal: O2_onelev.ncl: Only mol/L, mol/m3, mol/m-3, mmol/m3 or mmol m-3 are valid units for O2")
;;            status_exit(1)
;;    end if
;;    end if
;;    end if
    
;  dum = tmp
   dum@long_name=variable_info@long_name
; If input data on NEMO grid the dimesnions will be i and j, else lon and lat
; maybe better to switch on model type?
   dNames = getvardims(dum)
    if ( models@project(index) .eq. "GO_gridfile" ) then
     lon_name = "i"
     lat_name = "j"
    else 
     lon_name = "lon" 
     lat_name = "lat"
   end if 
  
           var_tmp = dum(time|:, $lat_name$|:,  $lon_name$|:, \
             {depth|diag_script_info@$"depthrange_basic"$(0):diag_script_info@$"depthrange_basic"$(1)} )
                dNames = getvardims(var_tmp)
                print("dNames is: " + dNames )
            delete(dum)
            dum = var_tmp(time|:,$lon_name$|:,$lat_name$|:,depth|0)
            delete(var_tmp)
           
       derived_field_type="TO2M"
    add_data_var(index, data_new, dum, variable)
    if ( models@project(index) .eq. "GO_gridfile" ) then
      in_file = interface_get_fullpath("O2", derived_field_type, index)
      out_file = interface_get_fullpath(variable, derived_field_type, index)
      fout = addfile(out_file, "w")
      fin = addfile(in_file,"r")
      fout->lat = fin->lat
      fout->lon = fin->lon
      fout->area = fin->area
    end if 

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Leaving O2.ncl", verbosity, 2)
  return(data_new)
end
