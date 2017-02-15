;;#############################################################################
;; Recipe to extract a complex variable from GFDL output
;;#############################################################################
;;
;; VARIABLE: snd
;; RECIPE: HS * dim_sum(CN, 1)
;; HISTORY: 20150710-maso_er: written.
;;
;;#############################################################################

load "./interface_scripts/constants.ncl"
undef("GFDL_recipe")
function GFDL_recipe(name: string)
local cn, hs, sic
begin

    ;; Read in the data
    cn = find_and_read_var("CN", True)
    hs = find_and_read_var("HS", True)

    ;; Integrate sea-ice concentration in the vertical
    sic = dim_sum_n_Wrap(cn, 1)

    ;; Calculate snd
    snd = hs  ; save metadata
    snd = hs * sic
    snd@units = "m"
    return(snd)

end
