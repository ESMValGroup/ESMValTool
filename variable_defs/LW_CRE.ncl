;
;  Requires: rlut:T2*s,rlutcs:T2*s
;
load "interface_scripts/data_handling.ncl"
load "interface_scripts/auxiliary.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "Longwave cloud radiative effect"
variable_info@units = "W m-2"

; Contour levels
variable_info@map_ref_Levels = ispan(0, 100, 10)
variable_info@map_diff_Levels = ispan(-30, 30, 5)

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined ini the
;;               'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp, dum, i, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering LW_CRE.ncl", verbosity, 4)
    data_new = True
    tmp1 = read_data(index, "rlut", "T2Ms")
    tmp2 = read_data(index, "rlutcs", "T2Ms")
    dum1 = extract_data_extended(index, "rlut", tmp1, -1, 0, 0)
    dum2 = extract_data_extended(index, "rlutcs", tmp2, -1, 0, 0)

    dum = dum1
    dum = dum2 - dum1
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving LW_CRE.ncl", verbosity, 4)
    return(data_new)
end
