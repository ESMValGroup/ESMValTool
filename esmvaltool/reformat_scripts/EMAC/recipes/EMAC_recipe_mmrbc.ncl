;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: mmrbc
;; RECIPE:   sum over BC modes
;; HISTORY:  20140912-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local modes, im, yy
begin

    ;; Selected species
    selspec = "BC"

    ;; Define all possible modes
    modes = get_aerosol_modes()
    aero = get_aerosol_spec()

    ;; Loop over species
    do is = 0, dimsizes(aero) - 1
        if (.not.isStrSubset(aero(is), selspec)) then
            continue
        end if
        ;; Loop over modes
        do im = 0, dimsizes(modes) - 1
            yy = find_and_read_var(aero(is) + modes(im), date, False)
            if (typeof(yy).ne."string") then
                info_output("Adding species " +  aero(is) + modes(im) + \
                            " (" + date + ")", verbosity, 2)
                if (isdefined("xx")) then
                    xx = xx + yy
                else
                    xx = yy
                end if
            end if
            delete(yy)
        end do
    end do

    ;; Convert units [mol/mol] --> [kg kg-1]
    if (xx@units.eq."mol/mol" .or. xx@units.eq."mol mol-1") then
        xx = xx * tofloat(xx@molarmass) / MWair
        xx@units = "kg kg-1"
    end if

    return(xx)

end
