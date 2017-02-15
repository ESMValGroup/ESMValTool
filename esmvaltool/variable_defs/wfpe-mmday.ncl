;
; Requires: wfpe:T2*s
;
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "Water Flux From Precipitation And Evaporation"
variable_info@units = "mm day-1"


undef("calculate")
function calculate(index[1] : integer,
                   variable[1] : string,
                   field_type[1] : string)
;;                 return_val[1] : logical
;; Arguments:
;;    index    - index to current infile defined in the
;;               'interface_data/ncl.interface'-file
;;    variable - string with current variable
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
;; Description:
;;    Converts water flux from P+E from units "kg m-2 s-1" to "mm day-1"
local result, scriptname, T, verbosity
begin
    scriptname = "variable_defs/wfpe-mmday.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, "", 4)

    result = read_data(index, "wfpe", "T2*s")
    T = extract_data(index, result, -1, 0, 0)
    T = T * 24 * 3600
    T@units = variable_info@units

    data_new = True
    add_data_var(index, data_new, T, variable)

    ;; if lat and lon are not coordinate variables (irregular grids), copy
    ;; lon from tauv and lat from tauu input file, this should be ok for both
    ;; staggered and destaggered grids.
    if (.not. iscoord(T, "lon") .or. .not. iscoord(T, "lat")) then
        ;; interface_get_fullpath is defined in 'interface_data/ncl.interface'
        in_file = interface_get_fullpath("wfpe", field_type, index)
        out_file = interface_get_fullpath(variable, field_type, index)

        fin = addfile(in_file, "r")
        fout = addfile(out_file, "w")

        vars = getfilevarnames(fin)
        if (.not. any("lon" .eq. vars) .or. .not. any("lat" .eq. vars)) then
            error_msg("f", get_script_name(), "", \
                      "variables lat and/or lon not found in " + in_file)
        end if

        fout->lon = fin->lon
        fout->lat = fin->lat

        ;; copy areacello if available
        if (any("areacello" .eq. vars)) then
            fout->areacello = fin->areacello
            info_output("lon, lat and areacello copied from " + in_file, \
                        verbosity, 2)
        else
            info_output("lon and lat copied from " + in_file, verbosity, 2)
        end if
    end if

    leave_msg(scriptname, "", 4)

    return(data_new)
end
