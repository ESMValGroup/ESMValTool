;
;  Requires: rsus:T2*s,rsds:T2*s
;
variable_info = True
variable_info@derived = True
variable_info@long_name = "sw_albedo"
variable_info@units = "W m-2 / W m-2"

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
local tmp1, tmp2, dum1, dum2, dum, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering alb.ncl", verbosity, 4)
    data_new = True
    tmp1 = read_data(index, "rsus", "T2Ms")
    tmp2 = read_data(index, "rsds", "T2Ms")
    dum1 = extract_data_extended(index, "rsus", tmp1, -1, 0, 0)
    dum2 = extract_data_extended(index, "rsds", tmp2, -1, 0, 0)

    dum = dum1
    dum = dum1 / where(dum2.ne.0, dum2, dum1@_FillValue)
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, dum, variable)

    info_output(">>>>>>>> Leaving alb.ncl", verbosity, 4)
    return(data_new)
end
