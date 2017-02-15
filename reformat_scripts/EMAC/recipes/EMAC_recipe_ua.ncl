;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: ua
;; RECIPE:   um1 / cos(lat)
;; HISTORY:  20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local yy
begin

    xx = find_and_read_var("um1", date, True)
    yy = conform(xx, tofloat(xx&lat), ind(getvardims(xx).eq."lat"))
    yy = yy * PI / 180.  ; [deg] --> [rad]
    xx = xx / cos(yy)

    ;; Add 10% to account for time-filter error in the output
    ;; Only for MESSy versions prior to 2.51
    ;; (Markus Kunze, priv. comm.)
    if (VERSION.eq.1 .or. (VERSION.eq.2 .and. SUBVERSION.le.50)) then
        xx = 1.10 * xx
    end if

    return(xx)

end
