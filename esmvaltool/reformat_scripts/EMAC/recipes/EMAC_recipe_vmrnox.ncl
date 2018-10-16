;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: vmrnox
;; RECIPE:   NO + NO2
;; HISTORY:  20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
begin

    xx = find_and_read_var("NO", date, True)
    xx = xx + find_and_read_var("NO2", date, True)
    return(xx)

end
