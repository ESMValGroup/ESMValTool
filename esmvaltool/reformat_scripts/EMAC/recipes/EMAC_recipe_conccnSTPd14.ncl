;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: conccnSTPd14
;; RECIPE:   number concentration of particles with dry diameter larger than
;;           14 nm in STP units
;; HISTORY:  20140507-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local modes, diams, sigma, im, nn, dd, yy
begin

    ;; Threshold diameter [m]
    dthr = 14.e-9

    ;; Define all possible modes and corresponding sigmas
    modes = get_aerosol_modes()
    sigma = modes@sigma

    ;; Number and dry diameters
    modes = "N" + modes
    diams = str_sub_str(modes, "N_", "dryrad_")

    ;; Loop over modes
    do im = 0, dimsizes(modes) - 1
        nn = find_and_read_var(modes(im), date, False)  ; number concentration
        dd = find_and_read_var(diams(im), date, False)  ; dry radius
        if (typeof(nn).ne."string" .and. typeof(dd).ne."string") then
            info_output("Adding mode " + modes(im) + \
                        ", " + diams(im) + ", " + sigma(im) + \
                        " (" + date + ")", verbosity, 2)
            dd = 2.*dd  ; radius to diameter
            yy = nn
            yy = yy/2. * erfc(log(dthr/dd) / (sqrt(2.)*log(sigma(im))))
            if (isdefined("xx")) then
                xx = xx + yy
            else
                xx = yy
            end if
        end if
        delete(nn)
        delete(dd)
    end do

    ;; Convert units [1/mol] --> [m-3-STP]
    if (xx@units.eq."1/mol" .or. xx@units.eq."mol-1") then
        xx = xx * p0stp  / Rgas / T0stp
        xx@units = "m-3-STP"
    end if

    return(xx)

end
