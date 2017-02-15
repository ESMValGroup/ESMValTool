;
; Requires: tas:T2*s
;
load "interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

variable_info = True
variable_info@derived = True
variable_info@units = "degC"

; Contour levels
variable_info@map_ref_Levels = ispan(-30, 30, 3)
variable_info@map_diff_Levels = ispan(-50, 50, 10) * 0.1

;;#############################################################################
load "interface_scripts/data_handling.ncl"

undef("calculate")
function calculate(index[1]:integer,
                   variable[1]:string,
                   field_type[1]:string)
;;
;; Arguments
;;    index: index to current infile defined in the 
;;           'interface_data/ncl.interface'-file
;;    variable: a string with current variable
;;    field_type: a string with field type classification
;;
;; Return value
;;    data_new: logical
;;
;; Description
;;    Converts temperature from degC to K.
;;
;; Caveats
;;
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;
;; Modification history
;;    20151106-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, result, T
begin

    funcname = "calculate"
    scriptname = "variable_defs/tas-degC.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    result = read_data(index, "tas", "T2*s")
    T = extract_data_extended(index, "tas", result, -1, 0, 0)
    T = T - 273.15
    T@units = variable_info@units

    data_new = True

    add_data_var(index, data_new, T, variable)
    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
