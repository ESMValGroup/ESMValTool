;
; Requires: clcci:T4*
;
variable_info = True
variable_info@derived = True
variable_info@long_name  = "CCI Middle Level Thick Cloud Area Fraction"
variable_info@short_name = "CCI Middle Level Thick Cloud Area Fraction"
variable_info@units      = "%"
variable_info@plot_units = "%"

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "interface_scripts/constants.ncl"

undef("calculate")
function calculate(index[1]:integer,
                   variable[1]:string,
                   field_type[1]:string)
;;
;; Arguments
;;    index: index to current infile defined in interface_data/ncl.interface
;;    variable; logical with relevant variable as string attribute
;;    field_type: string with field type classification
;;
;; Return value
;;    data_new: logical
;;
;; Description
;;    Calculates the CCI middle level thick cloud area fraction.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20170222-A_laue_ax: written.
;;
local funcname, scriptname, verbosity,  tmp1, clmtkcci
begin

    funcname = "calculate"
    scriptname = "variable_defs/clmtkcci.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Fetch data
    tmp  = read_data(index, "clcci", "T4*")
    tmp1 = extract_data_extended(index, "clcci", tmp, -1, 0, 0)

    dimNames = getvardims(tmp1)

    taus = ind(tmp1&$dimNames(1)$.gt.23)
    plevs = ind(tmp1&$dimNames(2)$.gt.44000 .and. tmp1&$dimNames(2)$.le.68000.)

    ; sum over tau and pressure level bins

    if ((dimsizes(taus).gt.1).and.(dimsizes(plevs).gt.1)) then
        avgdims = (/1, 2/)
    else if (dimsizes(taus).gt.1) then
        avgdims = (/1/)
    else
        avgdims = (/2/)
    end if
    end if

    clmtkcci = dim_sum_n_Wrap(tmp1(:,taus,plevs,:,:), avgdims)

    ;; Set metadata and return
    clmtkcci@long_name = variable_info@long_name
    clmtkcci@units = variable_info@units
    data_new = True
    add_data_var(index, data_new, clmtkcci, variable)

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
