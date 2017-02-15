;
; Requires: vmrco:T2*s
;
variable_info          = True
variable_info@derived  = True
variable_info@short_name = "CO VMR"
variable_info@units = "mole mole-1"
variable_info@plot_units = "nmol/mol"

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "interface_scripts/constants.ncl"
load "diag_scripts/lib/ncl/latlon.ncl"

undef("calculate")
function calculate(index[1] : integer,
                   variable[1] : string,
                   field_type[1] : string)
;;
;; Arguments
;;    index: index to current infile defined in interface_data/ncl.interface
;;    variable; logical with relevant variable as string attribute
;;    field_type: string with field type classification
;;
;; Return value
;;    data_new: logical
;;
;; Description:
;;    Extracts time series of a single position. For use with GLOBALVIEW data.
;;
;; Caveats
;;
;; Modification history
;;    20140630-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, tmp1, vmrco, vmrco_t, time, dum
begin

    funcname = "calculate"
    scriptname = "variable_defs/vmrco_nwr.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Parameters
    location = "Niwot Ridge, Colorado, United States"
    lat = 40.05
    lon = 360. - 105.58
    altitude = "3526 masl"

    ;; Fetch data
    tmp1 = read_data(index, "vmrco", "T2*s")
    vmrco_t = extract_data(index, tmp1, -1, 0, 0)

    ;; Define target field
    time = vmrco_t&time
    dum = new((/dimsizes(time)/), typeof(vmrco_t))
    dum!0 = "time"
    dum&time = time

    ;; Calculations
    dum = area_operations(vmrco_t, lat, lat, lon, lon, "extract", True)

    ;; Set metadata and output derived variable
    dum@long_name = variable_info@long_name
    dum@standard_name = variable_info@standard_name
    dum@short_name = variable_info@short_name

    dum@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dum, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
