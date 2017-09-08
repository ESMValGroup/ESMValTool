;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: clwvi
;; RECIPE:   xlvi + xivi
;; HISTORY:  20150504-A_laue_ax: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local yy
begin

    xx = find_and_read_var("xlvi", date, True)
    yy = find_and_read_var("xivi", date, True)
    xx = xx + yy

    return(xx)

end
