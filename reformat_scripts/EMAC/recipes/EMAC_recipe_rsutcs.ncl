;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: rsutcs
;; RECIPE:   (1) - (sraf0 - srad0d)
;;           (2) top level of flxs - srad0u - top level of flxsf
;; COMMENT:  Only the net outgoing sw radiation (sraf0) is available in EMAC.
;;           This is defined as sraf0 = sraf0d + sraf0u, where sraf0u is the
;;           clear-sky outgoing sw radiation. Therefore sraf0u = sraf0 - sraf0d
;;           and sraf0d = srad0d (incoming solar radiation). The minus sign is
;;           needed due to different conventions in EMAC and CMOR.
;; HISTORY:  20150507-A_laue_ax: added option to use 3-d radiation fields
;;                               if 2-d fields are not available
;;                               (e.g. EMAC RC1/RC2 runs)
;;           20140424-A_righ_ma: written.
;;
;;#############################################################################
load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local xx, yy
begin
    ; check for 2-d fields first
    xx = find_and_read_var("sraf0", date, False)
    if ((typeof(xx).eq."string").or.(abs(sum(xx)).lt.1.0e-6)) then
        ; 2-d field was not found or field contains only zeros
        ; ---> try using 3-d field instead
        delete(xx)
        ; flxs   = shortwave flux
        ; flxsf  = shortwave flux clearsky
        ; srad0u = top solar radiation upward
        xx = find_and_read_var("flxs", date, True)
        ; extract top level
        yy = xx(:, 0, :, :)
        copy_VarMeta(xx(:, 0, :, :), yy)
        delete(xx)
        xx = find_and_read_var("srad0u", date, True)
        yy = yy - xx
        delete(xx)
        xx = find_and_read_var("flxsf", date, True)
        yy = yy - xx(:, 0, :, :)
        return(yy)
    else
        yy = find_and_read_var("srad0d", date, True)
        xx = xx - yy
        xx = -1.0 * xx
        return(xx)
    end if
end
