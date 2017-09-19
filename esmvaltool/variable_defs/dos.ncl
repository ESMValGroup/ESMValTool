;
;  Requires: mrso:T2*s
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
;;    20160523-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, result, tmp, f, a, mrsofc, dos
begin

    funcname = "calculate"
    scriptname = "variable_defs/dos.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    result = read_data(index, "mrso", "T2*s")
    tmp = extract_data_extended(index, "mrso" , result, -1, 0, 0)

    f = addfile(interface_get_fullpath("mrso", field_types(0), index), "r")
    if (isfilevar(f, "mrsofc")) then
        mrsofc = f->mrsofc
    else
        error_msg("f", scriptname, funcname, "mrsofc fx file not found")
    end if

    ;; Check units
    if (tmp@units.ne.mrsofc@units) then
        error_msg("f", scriptname, funcname, \
            "mrso and mrsofc must have the same units")
    end if

    ;; Calculate ratio
    dos = tmp
    mrsofc = where(mrsofc.eq.0, mrsofc@_FillValue, mrsofc)
    dos = tmp / conform(tmp, mrsofc, (/1, 2/))
    dos@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dos, "dos")

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
