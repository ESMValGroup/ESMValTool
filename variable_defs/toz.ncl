;
; Requires: tro3:T3*,ps:T2*s
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "Total column ozone"
variable_info@short_name = "Tot. Col. Ozone"
variable_info@units = "DU"

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
local funcname, scriptname, verbosity, tmp1, tmp2, tro3, tro3_t, ps_t, time, \
    lat, lon, dum, ptop, plev, dp, o3nd
begin

    funcname = "calculate"
    scriptname = "variable_defs/toz.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Top pressure of the column [Pa]
    ptop = 10.

    ;; Fetch data
    tmp1 = read_data(index, "tro3", "T3*")
    tmp2 = read_data(index, "ps", "T2*s")
    tro3_t = extract_data_extended(index, "tro3", tmp1, -1, 0, 0)
    ps_t = extract_data_extended(index, "ps", tmp2, -1, 0, 0)

    ;; Convert units [1e9 --> mol mol-1]
    tro3_t = 1.e-9 * tro3_t

    ;; Define target field
    time = tro3_t&time
    lat = tro3_t&lat
    lon = tro3_t&lon
    dum = new((/dimsizes(time), dimsizes(lat), dimsizes(lon)/), typeof(tro3_t))
    dum!0 = "time"
    dum!1 = "lat"
    dum!2 = "lon"
    dum&time = time
    dum&lat = lat
    dum&lon = lon

    ;; Calculations
    do itim = 0, dimsizes(time) - 1

        tro3 = tro3_t(itim, :, :, :)
        ps = ps_t(itim, :, :)

        ;; Regrid ps to tro3, if necessary
        ;; Needed if this auxiliary field is taken from a different model.
        regridding = False
        if (dimsizes(tro3&lat).ne.dimsizes(ps&lat) .or. \
            dimsizes(tro3&lon).ne.dimsizes(ps&lon)) then
            regridding = True
        else
            if (any((tro3&lat - ps&lat).ne.0.) .or. \
                any((tro3&lon - ps&lon).ne.0.)) then
                regridding = True
            end if
        end if

        if(regridding) then
            temp = linint2(ps&lon, ps&lat, ps(lat|:, lon|:), \
                           True, tro3&lon, tro3&lat, 0)
            delete(ps)
            ps = temp
            ps!0 = "lat"
            ps&lat = tro3&lat
            ps!1 = "lon"
            ps&lon = tro3&lon
            delete(temp)
        end if

        ;; Pressure coordinate
        plev = tro3_t&plev
        if(plev@units.ne."Pa") then
            error_msg("f", scriptname, funcname, "pressure must be in Pa, " + \
                      "is " + plev@units)
        end if

        ;; Calculate total column
        dp = dpres_plevel(plev, ps, ptop, 0)
        delete(ps)
        o3nd = tro3 * tofloat(dp) / (Grav * MWair / MWo3)
        delete(tro3)
        delete(dp)
        dum(itim, :, :) = (/tofloat(dim_sum_n(o3nd, 0) / DU)/)
        delete(o3nd)

    end do

    ;; Set metadata and return
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dum, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
