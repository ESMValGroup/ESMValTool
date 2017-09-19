;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: rlutcs
;; RECIPE:   (1) - traf0
;;           (2) - top level of flxtf
;; HISTORY:  20150507-A_laue_ax: added option to use 3-d radiation field
;;                               if 2-d field is not available
;;                               (e.g. EMAC RC1/RC2 runs)
;;           20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local xx, yy
begin
    ; check for 2-d field first
    xx = find_and_read_var("traf0", date, False)
    if ((typeof(xx).eq."string").or.(abs(sum(xx)).lt.1.0e-6)) then
        ; 2-d field was not found or field contains only zeros
        ; ---> try using 3-d field instead
        delete(xx)
        ; flxtf = "longwave flux clearsky"
        xx = find_and_read_var("flxtf", date, True)
        ; extract top level
        yy = -1.0 * xx(:, 0, :, :)
        copy_VarMeta(xx(:, 0, :, :), yy)
        return(yy)
    else
        xx = -1. * xx
        return(xx)
    end if
end
