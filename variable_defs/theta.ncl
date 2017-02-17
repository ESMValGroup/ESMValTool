;
; Requires: ta:*3*
;
load "interface_scripts/data_handling.ncl"

variable_info          = True
variable_info@derived  = True
variable_info@long_name="Potential Temperature"
variable_info@units="K"

; Countour levels for zonalmean_profile.ncl
;variable_info@zonal_ref_Levels = ispan(200, 300, 5)
;variable_info@zonal_diff_Levels = ispan(-10, 10, 2)

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - logical with releveant variable as string attribute
;;    field_type  - string with field type classification
;; Return value:
;;    data_new - logical
local tmp, theta, i, dum, verbosity, pp
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering theta.ncl", verbosity, 2)
    data_new = True
    tmp = read_data(index, "ta", "*3*")
    dum = extract_data(index, tmp, -1, 0, 0)
    p0 = 100000.
    plev = dum&plev

    pp = new( dimsizes(dum), float)
    
    nlev = dimsizes(plev)
    ilev = 0
    do while(ilev.lt.nlev)
      pp(:,ilev,:,:) = tofloat((p0/plev(ilev))^0.286)
      ilev=ilev+1
    end do
    dum = dum*pp
    dum@long_name = variable_info@long_name
    dum@units = variable_info@units
    add_data_var(index, data_new, dum, variable)
    info_output(">>>>>>>> Leaving theta.ncl", verbosity, 2)
    return(data_new)
end
