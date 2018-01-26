;
; Requires: tro3:T3*
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "Ozone Profile"
variable_info@short_name = "Ozone Profile"
variable_info@units = "1e-9"
variable_info@plot_units = "ppmv"

; Contour levels
variable_info@colormap = "BlAqGrYeOrReVi200"
variable_info@acp_cnlvl = fspan(150., 500., 15)
variable_info@map_ref_Levels = ispan(150, 500, 15)
variable_info@map_diff_Levels = ispan(-60, 60, 5)

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "interface_scripts/constants.ncl"

undef("calculate")
function calculate(index[1]:integer,
                   variable[1]:string,
                   field_type[1]:string)
;;
;; Arguments
;;    index: index to current infile defined in interface_data/ncl.interface
;;    variable; logical with relevant variable as string attribute
;;    field_type: string with field type classification
;;
;; Return value
;;    data_new: logical
;;
;; Description
;;    Calculates total column ozone in [DU] given tro3 and ps. Assumes a top
;;    pressure of 10 Pa.
;;
;; Caveats
;;    Disabled section to re-order plev can be removed when ozone column
;;    diagnostics are tested and working.
;;
;; Reference
;;
;; Modification history
;;    20140328-A_gott_kl: written based on total_column.ncl by Irene Cionni.
;;
local funcname, scriptname, verbosity, tmp1, tro3_profile, tro3_t
begin

    funcname = "calculate"
    scriptname = "variable_defs/tro3_profile.ncl"
    verbosity = config_user_info@verbosity
    enter_msg(scriptname, funcname, 4)

    ;; Fetch data
    tmp1 = read_data(index, "tro3", "T3*")
    tro3_t = extract_data(index, tmp1, -1, 0, 0)

    ;; Convert units [1e9 --> mol mol-1]
    tro3_t = 1.e-9 * tro3_t

    ;; Compute zonal mean
    tro3_profile = dim_avg_Wrap(tro3_t)

    ;; Set metadata and return
    tro3_profile@long_name = variable_info@long_name
    tro3_profile@units = variable_info@units

    data_new = True
    add_data_var(index, data_new, tro3_profile, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
