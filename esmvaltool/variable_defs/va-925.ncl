;
;  Requires: va:*3*
;
load "interface_scripts/data_handling.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name="925mb wind"
variable_info@units="m s-1"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - Index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - String with field type classification
;; Return value:
;;    data_new - logical
local tmp, dum, i, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering va-925.ncl", verbosity, 4)
    tmp = read_data(index, "va", "*3*")
    dum = extract_data(index, tmp, 1, 92500., 92500.)

    dum@long_name = variable_info@long_name
    dum@units = variable_info@units

    data_new = True
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving va-925.ncl", verbosity, 4)
    return(data_new)
end
