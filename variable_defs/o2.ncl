;
;  Requires: o2:TO2Ms
;
variable_info = True
variable_info@derived = False
variable_info@long_name="o2"
variable_info@units = "mol m-3"
variable_info@MyDiag_title = "Oxygen"

load "interface_scripts/data_handling.ncl"
load "reformat_scripts/constants.ncl"
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
    info_output("<<<<<<<< Entering o2.ncl", verbosity, 2)

  data_new = True
  tmp = read_data(index,"o2","TO3M")
   
   dum = extract_data(index,tmp,-1,0,0)
 
  
   dum@long_name=variable_info@long_name
   check_units(dum,variable)
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
    in_file = interface_get_fullpath("o2", derived_field_type, index)
    system("mv " + in_file + " ./work/temp.nc " )
    add_data_var(index, data_new, dum, variable)
    if ( models@project(index) .eq. "GO_gridfile" ) then
      in_file = "./work/temp.nc"
      out_file = interface_get_fullpath(variable, derived_field_type, index)
      fout = addfile(out_file, "w")
      fin = addfile(in_file,"r")
      fout->lat = fin->lat
      fout->lon = fin->lon
      fout->area = fin->area
    end if 

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Leaving o2.ncl", verbosity, 2)
  return(data_new)
end
