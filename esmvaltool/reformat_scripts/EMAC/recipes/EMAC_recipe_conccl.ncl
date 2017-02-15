;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: conccl
;; RECIPE:   sum over Cl (SS) modes
;; HISTORY:  20150410-A_kais_ch: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local selspec, fac, modes, aero, repeat, is, im, yy, pp, tt
begin

    ;; Selected species
    selspec = "Cl"
    fac     = 1.0

    ;; Define all possible modes
    modes = get_aerosol_modes()
    aero = get_aerosol_spec()

    ;; If tracer "Cl" is not found, repeat procedure with `selname = "SS"' and
    ;; obtain Cl from SS by multiplication with `fac'
    repeat = True
    do while (repeat)
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
                        xx = xx + fac * yy
                    else
                        xx = yy
                    end if
                end if
                delete(yy)
            end do
        end do
        ;; Terminate while loop or repeat it with new variables
        if (selspec.eq."SS") then
            repeat = False
        else
            if (isdefined("xx")) then
                repeat = False
            else
                selspec = "SS"
                fac     = 0.5504
                info_output("Cl not found. Using " + \
                            tostring_with_format(fac, "%6.4f") + \
                            " * SS instead." , verbosity, 2)
            end if
        end if
    end do

    ;; Convert units [mol/mol] --> [kg m-3]
    if (xx@units.eq."mol/mol" .or. xx@units.eq."mol mol-1") then
        pp = find_and_read_var("press", date, True)
        tt = find_and_read_var("tm1", date, True)
        xx = xx * pp * tofloat(xx@molarmass) / Rgas / tt / 1000.
        xx@units = "kg m-3"
    end if

    return(xx)

end
