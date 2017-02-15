;
;  Requires: o2:TO3M
;
variable_info = True
variable_info@derived = True
variable_info@long_name="o2_onelev"
variable_info@MyDiag_title = "Oxygen"
variable_info@units = "mol m-3"
load "interface_scripts/data_handling.ncl"
load "reformat_scripts/constants.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
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
local tmp, dum, dimension, ind_r
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering O2_onelev.ncl", verbosity, 2)

    data_new = True
    tmp = read_data(index, "o2", "TO3M")

    dum = extract_data(index, tmp, -1, 0, 0)

    dum@long_name=variable_info@long_name

    dNames = getvardims(dum)

    i_j_indices = diag_script_info@models_using_i_j_indices
    ind_r = inlist(models@name(index), i_j_indices)
    if (any(ind_r .eq. True)) then
      lon_name = "i"
      lat_name = "j"
    else
      lon_name = "lon"
      lat_name = "lat"
    end if

    var_tmp = dum(time|:, $lat_name$|:,  $lon_name$|:, \
              {depth|diag_script_info@$"depthrange_basic"$(0):diag_script_info@$"depthrange_basic"$(1)} )
    dNames = getvardims(var_tmp)

    delete(dum)
    dum = var_tmp(time|:, $lon_name$|:, $lat_name$|:, depth|0)
    delete(var_tmp)

print("*********************************************************************")
print("*********************************************************************")
print("*********************************************************************")
printVarSummary(dum)
print("*********************************************************************")
print("*********************************************************************")
print("*********************************************************************")

; A-laue_ax+
    derived_field_type="TO2M"
; A-laue_ax-

    add_data_var(index, data_new, dum, variable)

    if (models@project(index) .eq. "GO_gridfile") then
      in_file = interface_get_fullpath("o2", "TO3M", index)
      out_file = interface_get_fullpath(variable, derived_field_type, index)
      fout = addfile(out_file, "w")
      fin = addfile(in_file,"r")
      fout->lat = fin->lat
      fout->lon = fin->lon
      fout->area = fin->area
    end if

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Leaving o2_onelev.ncl", verbosity, 2)
    return(data_new)
end
