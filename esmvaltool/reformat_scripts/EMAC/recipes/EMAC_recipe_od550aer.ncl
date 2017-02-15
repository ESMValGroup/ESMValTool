;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: od550aer
;; RECIPE 1: aot550_total
;; RECIPE 2: dim_sum(aot_opt_TOT_550_total, 1)
;; HISTORY:  20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local xx
begin

    ;; MESSy1: simply read variable
    if (VERSION.eq.1) then
        yy = find_and_read_var("aot550_total", date, True)
    end if

    ;; MESSy2: sum over vertical coordinate
    if (VERSION.eq.2) then
        xx = find_and_read_var("aot_opt_TOT_550_total", date, True)
        yy = dim_sum_n_Wrap(xx, 1)
    end if

    yy@units = "1"
    return(yy)

end
