;
;  Requires: cSoil:*2*s
;
variable_info = True
variable_info@derived = True 
variable_info@units = "kg m-2"
variable_info@plot_units = "PgC"

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
;;    Corrects cSoil for coastal regions. If available, the land mask for the
;;    given model is used, otherwise the NCL mask is applied.
;;
;; Caveats
;;
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;
;; Modification history
;;    20151029-A_righ_ma: written based on a previous version by wenz_sa.
;;
local funcname, scriptname, verbosity, result, tmp, f, a, sftlf, CSOIL
begin

    funcname = "calculate"
    scriptname = "variable_defs/CSOIL.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    result = read_data(index, "cSoil", "*2*s")
    tmp = extract_data_extended(index, "cSoil" , result, -1, 0, 0)

    f = addfile(interface_get_fullpath("cSoil", field_types(0), index), "r")
    if (isfilevar(f, "sftlf")) then
        sftlf = f->sftlf
        sftlf = sftlf / 100.  ; [%] --> [1]
    else
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
        a = f->LSMASK
        sftlf = byte2flt(landsea_mask(a, tmp&lat, tmp&lon))
        sftlf = where(sftlf.gt.1., 1., sftlf)
        info_output("Land-mask file not available, using NCL default mask", \
                    verbosity, 1)
    end if
    sftlf = where(sftlf.eq.0, sftlf@_FillValue, sftlf)

    ;; Apply mask
    CSOIL = tmp
    CSOIL = CSOIL * conform(tmp, sftlf, (/1, 2/))

    data_new = True
    add_data_var(index, data_new, CSOIL, "CSOIL")

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
