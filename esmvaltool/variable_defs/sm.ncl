;
;  Requires: mrsos:T2*s
;
variable_info = True
variable_info@derived = True 
variable_info@units = "m3 m-3"

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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
;;    Calculates degree of saturation as the ratio between mrso and mrsofc
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20160824-A_laue_ax: written.
;;
local funcname, scriptname, verbosity, result, tmp, f, a, mrsofc, dos
begin

    funcname = "calculate"
    scriptname = "variable_defs/sm.ncl"
    verbosity = diag_script_info@verbosity
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    result = read_data(index, "mrsos", "T2*s")
    tmp = extract_data_extended(index, "mrsos" , result, -1, 0, 0)

    f = addfile(interface_get_fullpath("mrsos", field_types(0), index), "r")
    if (isfilevar(f, "depth_bnds")) then
        depth_bnds = f->depth_bnds
    else
        error_msg("f", scriptname, funcname, "depth_bnds not found")
    end if

    ;; calculate layer thickness
    depth = tofloat(depth_bnds(1) - depth_bnds(0))
    if (depth.le.0) then
        error_msg("f", scriptname, funcname, "invalid layer thickness")
    end if

    ;; convert moisture content of soil layer (kg/m2) into volumetric
    ;; soil moisture (m3/m3)
    ;; assume density of water = 998.2 kg/m3 (T = 20 degC)
    sm = tmp / depth / 998.2
    copy_VarCoords(tmp, sm)
    sm@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, sm, "sm")

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
