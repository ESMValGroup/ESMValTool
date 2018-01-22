;
;  Requires: shrubFrac:T2*s,treeFrac:T2*s
;
load "interface_scripts/data_handling.ncl"
load "interface_scripts/auxiliary.ncl"

variable_info = True
variable_info@derived = True
variable_info@long_name = "Upper Mid to High Vegetation Fraction"
variable_info@units = "%"

; Contour levels
variable_info@map_ref_Levels = ispan(0, 100, 10)
variable_info@map_diff_Levels = ispan(-30, 30, 5)

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;; Arguments:
;;    index    - index to current infile defined ini the
;;               'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local funcname, scriptname, verbosity, variable, tmp1, tmp2, dum1, dum2, dum1a, dum2a, dum
begin
    
    funcname = "calculate"
    scriptname = "variable_defs/shrubNtreeFrac.ncl"
    verbosity = diag_script_info@verbosity
    enter_msg(scriptname, funcname, 4)

    data_new = True

    tmp1 = read_data(index, "shrubFrac", "T2Ms")
    tmp2 = read_data(index, "treeFrac", "T2Ms")

    dum1 = extract_data(index, tmp1, -1, 0, 0)
    dum2 = extract_data(index, tmp2, -1, 0, 0)

    dum1a = where(ismissing(dum1), 0.0, dum1)
    dum2a = where(ismissing(dum2), 0.0, dum2)

    dum = dum1a + dum2a

    delete(dum1a)
    delete(dum2a)

    dum1 = where(ismissing(dum1).and.ismissing(dum2), dum1@_FillValue, dum)

    delete(dum)
    delete(dum2)

    dum1@long_name = variable_info@long_name
    dum1@units = variable_info@units
    add_data_var(index, data_new, dum1, variable)

    delete(dum1)

    leave_msg(scriptname, funcname, 4)
    return(data_new)
end
