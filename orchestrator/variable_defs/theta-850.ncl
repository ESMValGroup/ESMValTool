;
;  Requires: ta:*3*
;
load "interface_scripts/data_handling.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name="850 hPa Pot. Temperature"
variable_info@units="K"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp, dum, i, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering theta-850.ncl", verbosity, 4)
    tmp = read_data(index, "ta", "*3*")
    dum = extract_data(index, tmp, 1, 85000., 85000.)

    dum = dum*(1000./850.)^0.286

    dum@long_name = variable_info@long_name
    dum@units = variable_info@units

    data_new = True
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving theta-850.ncl", verbosity, 4)
    return(data_new)
end
