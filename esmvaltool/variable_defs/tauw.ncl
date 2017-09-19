;
;  Requires: tauu:TO2Ms,tauv:TO2Ms
;
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

variable_info = True
variable_info@derived = True
variable_info@standard_name = "surface_downward_stress"
variable_info@long_name = "Surface Downward Wind Stress"
variable_info@units = "Pa"

undef("calculate")
function calculate(index[1] : integer,
                   variable[1] : string,
                   field_type[1] : string)
;;                 return_val[1] : logical
;; Arguments:
;;    index    - index to current infile defined in the
;;               'interface_data/ncl.interface'-file
;;    variable - logical with releveant variable as string attribute
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp_u, tmp_v, dum_u, dum_v, dum_w, in_file, out_file, fin, fout, \
    field_type, vars, verbosity, scriptname
begin
    scriptname = "variable_defs/tauw.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, "", 4)

    field_type = getenv("ESMValTool_derived_field_type")

    tmp_u = read_data(index, "tauu", "TO2Ms")
    tmp_v = read_data(index, "tauv", "TO2Ms")

    ; See ./interface_scripts/data_handling.ncl
    dum_u = extract_data(index, tmp_u, -1, 0, 0)
    dum_v = extract_data(index, tmp_v, -1, 0, 0)

    ;; calculate wind stress magnitude
    dum_w = dum_u
    dum_w = sqrt(dum_u^2 + dum_v^2)
    dum_w@standard_name = variable_info@standard_name
    dum_w@long_name = variable_info@long_name
    if (dum_u@units .ne. variable_info@units .or. \
        dum_v@units .ne. variable_info@units) then
        error_msg("f", scriptname, "", \
                  "units mismatch and no conversion implemented")
    end if
    dum_w@units = variable_info@units

    ;; See ./interface_scripts/data_handling.ncl
    data_new = True
    add_data_var(index, data_new, dum_w, variable)

    ;; if lat and lon are not coordinate variables (irregular grids), copy
    ;; lon from tauv and lat from tauu input file, this should be ok for both
    ;; staggered and destaggered grids.
    if (.not. iscoord(dum_w, "lon") .or. .not. iscoord(dum_w, "lat")) then
        ; 'interface_get_fullpath' is defined in 'interface_data/ncl.interface'
        inu_file = interface_get_fullpath("tauu", "T*2Ms", index)
        inu_file = systemfunc("ls " + inu_file)  ; Remove globbing

        inv_file = interface_get_fullpath("tauv", "T*2Ms", index)
        inv_file = systemfunc("ls " + inv_file)  ; Remove globbing

        out_file = interface_get_fullpath(variable, field_type, index)

        finu = addfile(inu_file, "r")
        finv = addfile(inv_file, "r")
        fout = addfile(out_file, "w")

        varsu = getfilevarnames(finu)
        varsv = getfilevarnames(finv)
        if (.not. any("lon" .eq. varsv) .or. .not. any("lat" .eq. varsu)) then
            error_msg("f", scriptname, "", \
                      "variables lat and/or lon not found in " + inv_file + \
                      " and " + inu_file)
        end if

        fout->lon = finv->lon
        fout->lat = finu->lat

        ;; copy areacello if available
        if (any("areacello" .eq. varsu)) then
            fout->areacello = finu->areacello
            info_output("lat and areacello copied from " + inu_file, \
                        verbosity, 2)
        else
            info_output("lat copied from " + inu_file, verbosity, 2)
        end if
        info_output("lon copied from " + inv_file, verbosity, 2)
    end if

    leave_msg(scriptname, "", 4)

    return(data_new)
end
