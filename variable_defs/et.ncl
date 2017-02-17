;
; Requires: hfls:T2*s
;
load "./interface_scripts/messaging.ncl"
load "interface_scripts/data_handling.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "Evapotranspiration"
variable_info@units = "mm/d"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - string with current variable
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
;; Description:
;;    Converts 'Surface Upward Latent Heat Flux' from units
;;    [W m-2] to [mm day-1]
local result, T, verbosity, et_latent_heat_of_vaporisation
begin
    et_latent_heat_of_vaporisation = 2.465e6  ; [J kg-1]
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering calculate (et)", verbosity, 4)
    result = read_data(index, "hfls", "T2*s")
    T = extract_data(index, result, -1, 0, 0)
    T = T * 24 * 3600 / et_latent_heat_of_vaporisation
    T@units = variable_info@units
    T@_FillValue = 1e20

    data_new = True
    add_data_var(index, data_new, T, variable)
    info_output(">>>>>>>> Leaving calculate (et)", verbosity, 4)
    return(data_new)
end

