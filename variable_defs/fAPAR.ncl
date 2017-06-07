;
;  Requires: lai:T2*s
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "fAPAR"
variable_info@units = "-"

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
local tmp, dum, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering fAPAR.ncl", verbosity, 4)
    data_new = True
    tmp = read_data(index, "lai", "T2Ms")
    dum = extract_data_extended(index, "lai", tmp, -1, 0, 0)

    dum = 1-exp(-0.5*dum)
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving fAPAR.ncl", verbosity, 4)
    return(data_new)
end
