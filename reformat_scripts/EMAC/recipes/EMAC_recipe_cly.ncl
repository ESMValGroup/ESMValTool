;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: cly
;; RECIPE:   HCl + ClNO3 + HOCl + ClO + Cl + 2*Cl2O2 + 2*Cl2 + OClO + BrCl
;; HISTORY:  20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local specs, fac, sp, yy
begin

    ;; Define species
    specs = (/"HCl", "ClNO3", "HOCl", "ClO", "Cl", "Cl2O2", "Cl2", "OClO", \
              "BrCl"/)

    ;; Multiplication factor
    fac   = (/1, 1, 1, 1, 1, 2, 2, 1, 1/)

    ;; Loop over species
    do sp = 0, dimsizes(specs) - 1
        yy = find_and_read_var(specs(sp), date, False)
        if (typeof(yy).ne."string") then
            yy = fac(sp) * yy
            if (isdefined("xx")) then
                xx = xx + yy
            else
                xx = yy
            end if
        end if
    end do
    return(xx)

end
