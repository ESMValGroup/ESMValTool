;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: sftlf
;; RECIPE:   100. * slf
;; HISTORY:  20150824-A_laue_ax: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
begin

    xx = find_and_read_var("slf", date, True)
    xx = 100. * xx
    xx@units = "%"
    return(xx)

end
