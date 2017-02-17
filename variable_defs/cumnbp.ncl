;
;  Requires: nbp:*2*s
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "Cumulated NBP"
variable_info@units = "kg m-2 s-1"
variable_info@plot_units = "PgC y-1"

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
;;    Calculates cumulated NBP.
;;
;; Caveats
;;
;; References
;;    The default NCL mask is $NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc
;;
;; Modification history
;;    20151029-A_wenz_sa: written.
;;
local funcname, scriptname, verbosity, result, tmp, f, a, sftlf, NBP
begin

    funcname = "calculate"
    scriptname = "variable_defs/cumnbp.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    result = read_data(index, "nbp", "*2*s")
    tmp = extract_data_extended(index, "nbp", result, -1, 0, 0)

    f = addfile(interface_get_fullpath("nbp", field_types(0), index), "r")
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
    NBP = tmp
    NBP = NBP * conform(tmp, sftlf, (/1, 2/))

    ;; Integrate
    idx_new = dimsizes(NBP)
    NBP_tmp = new((/idx_new(0), idx_new(1), idx_new(2)/), typeof(NBP))
    NBP_tmp(0, :, :) = 0.0
    do y = 0, dimsizes(NBP(:, 0, 0)) - 2
        NBP_tmp(y + 1, :, :) = NBP_tmp(y, :, :) + NBP(y, :, :)
    end do
    delete(y)
    NBP = (/NBP_tmp/)

    NBP!0 = "time"
    NBP&time = tmp&time
    NBP!1 = "lat"
    NBP&lat = tmp&lat
    NBP!2 = "lon"
    NBP&lon = tmp&lon
    NBP@units = variable_info@units
    NBP@_FillValue = tmp@_FillValue
    delete([/tmp, NBP_tmp/])

    copy_VarAtts(variable_info, NBP)
    data_new = True
    add_data_var(index, data_new, NBP, "cumnbp")

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
