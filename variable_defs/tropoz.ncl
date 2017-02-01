;
; Requires: tro3:T3*,ps:T2*s,ta:T3*
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "Tropospheric column ozone"
variable_info@short_name = "Trop. Col. Ozone "
variable_info@units = "DU"

; Contour levels
variable_info@acp_cnlvl = fspan(0., 70., 15)
variable_info@acp_dilvl = ispan(-8, 12, 2)
variable_info@map_ref_Levels = ispan(0, 70, 5)
variable_info@map_diff_Levels = ispan(-30, 30, 5)

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "interface_scripts/constants.ncl"
load "interface_scripts/auxiliary.ncl"


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
;;    Calculates tropospheric column ozone in [DU] given tro3, ps and ta.
;;    Assumes a top pressure of 10 Pa.
;;
;; Caveats
;;    Disabled section to re-order plev can be removed when ozone column
;;    diagnostics are tested and working.
;;
;; Reference
;;
;; Modification history
;;    20140328-A_gott_kl: written based on tropospheric_column.ncl by
;;                        Irene Cionni.
;;
local funcname, scriptname, verbosity, tmp1, tmp2, tro3, tro3_t, ps_t, time, \
    lat, lon, dum, ptop, plev, dp, o3nd, tmp3, ta_t, ta, trop_h, trop_hD, pm, \
    trop_o3col
begin

    funcname = "calculate"
    scriptname = "variable_defs/tropoz.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Top pressure of the column [Pa]
    ptop = 10.  ; top of the pressure column [Pa]

    ;; Flag for NCL routine trop_wmo to test lapse rate
    opt = True
    opt@lapsec = 2.0  ; default

    ;; Fetch data
    tmp1 = read_data(index, "tro3", "T3*")
    tmp2 = read_data(index, "ps", "T2*s")
    tmp3 = read_data(index, "ta", "T3*")
    tro3_t = extract_data(index, tmp1, -1, 0, 0)  ; tro3
    ps_t = extract_data(index, tmp2, -1, 0, 0)  ; ps
    ta_t = extract_data(index, tmp3, -1, 0, 0)  ; ta

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
        ta = ta_t(itim, :, :, :)

        ;; Pressure coordinate
        plev = tro3_t&plev
        if(plev@units.ne."Pa") then
            error_msg("f", scriptname, funcname, "pressure must be in Pa, " + \
                      "is " + plev@units)
        end if

        ;; trop_wmo requires monotonically increasing pressure levels
        if (plev(0).gt.plev(1)) then
            plev = plev(::-1)
            tro3 = tro3(::-1, :, :)
        end if
        plev_ta = ta&plev
        if (plev_ta(0).gt.plev_ta(1)) then
            ta = ta(::-1, :, :)
            plev_ta = plev_ta(::-1)
        end if

        ;; Regrid ps & ta to tro3, if necessary
        ;; Needed, if those auxiliary fields are taken from a different model.
        ;; If pressure levels of tro3 and ta do not match, it is assumed
        ;; that ps and ta stem from the same model and need to be regridded.
        regridding = False
        if(dimsizes(plev).ne.dimsizes(plev_ta)) then
            regridding = True
        else
            if(any((plev - plev_ta).ne.0.)) then
                regridding = True
            end if
        end if

        if(regridding) then
            info_output("Regridding ta and ps in tropoz.ncl as they seem " + \
                        "to stem from a different simulation than tro3", \
                        verbosity, 1)

            ;; ta
            temp = int2p_n(plev_ta, ta, plev, 1, 0)
            copy_VarCoords_l1(ta, temp)
            temp!0 = "plev"
            temp&plev = plev
            delete(ta)
            ta = linint2(temp&lon, temp&lat, temp(plev|:, lat|:, lon|:), \
                         True, tro3&lon, tro3&lat, 0)
            copy_VarCoords(tro3, ta)
            delete(temp)

            ;; ps
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

        ;; WMO tropopause
        trop_h = trop_wmo(plev, ta(lat|:, lon|:, plev|:), 1, opt)
        trop_h!0 = "lat"
        trop_h!1 = "lon"
        trop_h&lat = ta&lat
        trop_h&lon = ta&lon
        trop_h@_FillValue = 1.e20
        delete(ta)
        trop_h = where(trop_h.le.0., trop_h@_FillValue, trop_h)
        trop_hD = conform_dims((/dimsizes(plev), dimsizes(lat),\
                                dimsizes(lon)/), trop_h, (/1, 2/))
        delete(trop_h)

        ;; calculate tropospheric column
        dp = dpres_plevel(plev, ps, ptop, 0)
        delete(ps)
        o3nd = tro3 * tofloat(dp) / (Grav * MWair / MWo3)
        delete(tro3)
        delete(dp)
        pm = conform_dims((/dimsizes(plev), dimsizes(lat), dimsizes(lon)/),\
                          plev, (/0/))
        trop_o3col = where(pm.gt.trop_hD, o3nd, o3nd@_FillValue)
            ;; mask everything except the troposphere
        delete(trop_hD)
        delete(o3nd)
        delete(pm)
        dum(itim, :, :) = (/tofloat(dim_sum_n(trop_o3col, 0) / DU)/)
        delete(trop_o3col)

    end do

    ;; Set metadata and return
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dum, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
