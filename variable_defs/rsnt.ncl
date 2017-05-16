;
;  Requires: rsdt:T2*s,rsut:T2*s
;
load "interface_scripts/data_handling.ncl"
load "interface_scripts/auxiliary.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "TOA Net downward Shortwave Radiation"
variable_info@units = "W m-2"

; Contour levels
variable_info@map_ref_Levels = ispan(0, 360, 20)
variable_info@map_diff_Levels = ispan(-100, 125, 25)


undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the
;;               'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp, tmp1, tmp2, dum1, dum2, dum, i, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering rsnt.ncl", verbosity, 4)
    data_new = True
    tmp1 = read_data(index, "rsdt", "T2Ms")
    tmp2 = read_data(index, "rsut", "T2Ms")
    dum1 = extract_data_extended(index, "rsdt", tmp1, -1, 0, 0)
    dum2 = extract_data_extended(index, "rsut", tmp2, -1, 0, 0)

    dum = dum1
    dum = dum1 - dum2
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving rsnt.ncl", verbosity, 4)
    return(data_new)
end
