;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: zg
;; RECIPE:   (geopot + geosp) / 9.81
;; HISTORY:  20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local yy
begin

    xx = find_and_read_var("geopot", date, True)
    yy = find_and_read_var("geosp", date, True)
    xx = xx + conform(xx, yy, (/0, ind("lat".eq.getvardims(xx)), \
                               ind("lon".eq.getvardims(xx))/))
    xx = xx / Grav  ; [m2/s2] --> [m]
    xx@units = "m"
    return(xx)

end
