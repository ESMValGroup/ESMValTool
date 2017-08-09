;
;  Requires: ta:*3*
;
load "interface_scripts/data_handling.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "air temperature"
variable_info@units = "K"
variable_info@MyDiag_title = "200 hPa"

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the
;;               'interface_data/ncl.interface'-file
;;    variable - logical with releveant variable as string attribute
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp, dum, dimension, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering MyVar.ncl", verbosity, 4)
    tmp = read_data(index, "ta", "*3*")
    dimension = 1  ; Level

    ; See ./interface_scripts/data_handling.ncl
    dum = extract_data(index, tmp, dimension, 20000., 20000.)
    dum@long_name = variable_info@long_name

    data_new = True
    ; See ./interface_scripts/data_handling.ncl
    add_data_var(index, data_new, dum, variable)
    info_output(">>>>>>>> Leaving MyVar.ncl", verbosity, 4)
    return(data_new)
end
