;
; Requires: tro3:T3*
;
variable_info          = True
variable_info@derived  = True
variable_info@short_name = "Ozone VMR"
variable_info@units = "1e-9"
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
;; Description
;;    Extracts time series of a single position. For use with Tilmes data.
;;
;; Caveats
;;
;; Modification history:
;;    20140630-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, tmp1, data, tro3, tro3_t, time, plev, \
    dum, lon, lat, location
begin

    funcname = "calculate"
    scriptname = "variable_defs/tro3_NHext.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; (Hardcoded) Parameters
    ;; See global attributes of OBS_Tilmes_insitu_NHext_T1M_tro3.nc
    ;; and reformat_scripts/obs/reformat_obs_Tilmes.ncl
    lat = tofloat((/83, 59, 52, 53, 80, 53, 48, 22, 35, 32, 52, 60, 52, \
                   26, 79, 47, 50, 74, 71, 36, 51, 40/))
    lon = tofloat((/298, 266, 5, 246, 274, 300, 11, 114, 87, 131, 21, 359, \
                   14, 123, 12, 7, 15, 265, 338, 140, 4, 284/))
    location = "NH Extratropics"

    ;; Fetch data
    tmp1 = read_data(index, "tro3", "T3*")
    tro3_t = extract_data(index, tmp1, -1, 0, 0)  ; tro3

    ;; Target = input coordinates
    time = tro3_t&time
    plev = tro3_t&plev

    ;; Collect profiles @ sondes positions in auxiliary field
    data = new((/dimsizes(time), dimsizes(plev), dimsizes(lat)/), \
               typeof(tro3_t))
    do i = 0, dimsizes(lat) - 1
        data(:, :, i) = area_operations(tro3_t, lat(i), lat(i), \
                                        lon(i), lon(i), "extract", True)
    end do

    ;; Average (no weighting, as in reformat_obs_Tilmes.ncl)
    dum = dim_avg(data)
    delete(data)

    ;; Add coordinates
    dum!0 = "time"
    dum&time = time
    dum!1 = "plev"
    dum&plev = plev

    ;; Set metadata and output derived variable
    dum@long_name = variable_info@long_name
    dum@standard_name = variable_info@standard_name
    dum@short_name = variable_info@short_name
    dum@location = location

    dum@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dum, variable)  ; data_handling.ncl

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
