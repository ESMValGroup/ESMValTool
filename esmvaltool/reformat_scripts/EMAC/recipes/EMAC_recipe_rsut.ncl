;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: rsut
;; RECIPE:   - srad0u
;; HISTORY:  20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local xx
begin

    xx = find_and_read_var("srad0u", date, True)
    xx = -1. * xx
    return(xx)

end
