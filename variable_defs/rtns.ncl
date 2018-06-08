;
;  Requires: rsds:T2*s,rsus:T2*s,rlds:T2*s,rlus:T2*s
;
load "interface_scripts/auxiliary.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "Surface Net downward Total Radiation"
variable_info@units = "W m-2"

; contour levels
variable_info@map_ref_Levels = ispan(-70, 70, 10)
variable_info@map_diff_Levels = ispan(-90, 90, 15)

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
    info_output("<<<<<<<< Entering rtns.ncl", verbosity, 4)
    data_new = True
    tmp1 = read_data(index, "rsus", "T2Ms")
    tmp2 = read_data(index, "rsds", "T2Ms")
    dum1 = extract_data(index, tmp1, -1, 0, 0)
    dum2 = extract_data(index, tmp2, -1, 0, 0)

    dum = dum1
    dum = dum2 - dum1

    tmp1 = read_data(index, "rlus", "T2Ms")
    tmp2 = read_data(index, "rlds", "T2Ms")
    dum1 = extract_data(index, tmp1, -1, 0, 0)
    dum2 = extract_data(index, tmp2, -1, 0, 0)

    dum = dum + dum2 - dum1

    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving rtns.ncl", verbosity, 4)
    return(data_new)
end
