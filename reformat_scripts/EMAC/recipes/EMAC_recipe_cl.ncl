;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: cl
;; RECIPE:   100. * aclcac
;; HISTORY:  20140929-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
begin

    xx = find_and_read_var("aclcac", date, True)
    xx = 100. * xx
    xx@units = "%"
    return(xx)

end
