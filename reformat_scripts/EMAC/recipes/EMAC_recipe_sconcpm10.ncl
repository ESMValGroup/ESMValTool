;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: sconcpm10
;; RECIPE:   apply the error function for size selection to aerosol dry mass
;; HISTORY:  20141009-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local dthr, modes, aero, spec, im, is, mw, spmass, tmass, lnsig, xx, pp, tt
begin

    ;; Threshold diameter [m]
    dthr = 10.e-6

    ;; Define all possible modes and species
    modes = get_aerosol_modes()
    aero = get_aerosol_spec()

    ;; Loop over modes
    do im = 0, dimsizes(modes) - 1

        ;; Loop over species (excluding water)
        do is = 0, dimsizes(aero) - 1
            if (isStrSubset(aero(is), "H2O")) then
                continue
            end if
            spec = aero(is) + modes(im)
            spmass = find_and_read_var(spec, date, False)
            if (typeof(spmass).ne."string") then
                info_output("Adding species " + spec + " (" + date + ")", \
                            verbosity, 2)
                mw = tofloat(spmass@molarmass)
                if (isdefined("tmass")) then
                    tmass = tmass + spmass * mw  ; [mol/mol] --> [g/mol]
                else
                    tmass = spmass  ; save metadata
                    tmass = spmass * mw  ; [mol/mol] --> [g/mol]
                end if
            end if
            delete(spmass)
        end do

        ;; Apply error-function for size selection
        spec = "dryrad" + modes(im)
        dd = find_and_read_var(spec, date, False)
        if (typeof(dd).ne."string") then
            info_output("Calculating mode " + modes(im) + " (" + date + ")", \
                        verbosity, 2)
            lnsig = log(modes@sigma(im))
            dd = log(2. * dd) + 3. * lnsig^2  ; log(Dv)
            frac = 0.5 + 0.5 * erf((log(dthr) - dd) / (sqrt(2.) * lnsig))
            if (isdefined("xx")) then
                xx = xx + tmass * frac
            else
                xx = tmass  ; save metadata
                xx = tmass * frac
            end if
            delete(tmass)
        end if
        delete(dd)

    end do

    ;; Convert units [g/mol] --> [kg m-3]
    pp = find_and_read_var("press", date, True)
    tt = find_and_read_var("tm1", date, True)
    xx = xx * pp / Rgas / tt / 1000.
    xx@units = "kg m-3"

    return(xx)

end
