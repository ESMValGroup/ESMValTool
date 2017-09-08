;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: mmrbcfree
;; RECIPE:   sum over BC-free modes.
;; HISTORY:  20140912-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local modes, aero, spec, im, is, mw, yy
begin

    ;; Define all possible modes and species
    modes = get_aerosol_modes()
    aero = get_aerosol_spec()

    ;; Loop over species (excluding BC)
    do is = 0, dimsizes(aero) - 1
        if (isStrSubset(aero(is), "BC")) then
            continue
        end if
        ;; Loop over modes
        do im = 0, dimsizes(modes) - 1
            spec = aero(is) + modes(im)
            yy = find_and_read_var(spec, date, False)
            if (typeof(yy).ne."string") then
                info_output("Adding species " + spec + " (" + date + ")", \
                            verbosity, 2)
                mw = tofloat(yy@molarmass)
                if (isdefined("xx")) then
                    xx = xx + yy * mw  ; [mol/mol] --> [g/mol]
                else
                    xx = yy  ; save metadata
                    xx = yy * mw  ; [mol/mol] --> [g/mol]
                end if
            end if
            delete(yy)
        end do
    end do

    ;; Convert units [g/mol] --> [kg kg-1]
    if (xx@units.eq."mol/mol" .or. xx@units.eq."mol mol-1") then
        xx = xx / MWair
        xx@units = "kg kg-1"
    end if

    return(xx)

end
