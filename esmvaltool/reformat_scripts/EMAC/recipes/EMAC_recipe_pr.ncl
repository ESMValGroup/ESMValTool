;;#############################################################################
;; Recipe to extract a complex variable from EMAC output
;;#############################################################################
;;
;; VARIABLE: pr
;; RECIPE 1: rainflux + cv_precflx + snowflux + cv_snowflx
;; RECIPE 2: aprl + aprc
;; HISTORY:  20150507-A_kune_dm: written.
;;
;;#############################################################################
;load "./interface_scripts/constants.ncl"
;undef("EMAC_recipe")
;function EMAC_recipe(date: string)
;local xx, zz
;begin

 ;   xx = find_and_read_var("cv_precflx", date, True)
 ;   xx = xx + find_and_read_var("cv_snowflx", date, True)
 ;   xx = xx + find_and_read_var("rainflux", date, True)
 ;   xx = xx + find_and_read_var("snowflux", date, True)

    ;; Extract lowermost level
 ;   zz = getvardims(xx)
 ;   yy = xx(:, dimsizes(xx&$zz(1)$) - 1, :, :)
 ;   return(yy)

;end

load "./interface_scripts/constants.ncl"
undef("EMAC_recipe")
function EMAC_recipe(date: string)
local xx, yy, zz, qq
begin

    qq = find_and_read_var("aprl", date, False)

    if (typeof(qq).eq."string") then

        xx = find_and_read_var("cv_precflx", date, True)
        xx = xx + find_and_read_var("cv_snowflx", date, True)
        xx = xx + find_and_read_var("rainflux", date, True)
        xx = xx + find_and_read_var("snowflux", date, True)
    
        ;; Extract lowermost level
        zz = getvardims(xx)
        yy = xx(:, dimsizes(xx&$zz(1)$) - 1, :, :)
        return(yy)

    else

        ;; The contribution of snow is already included in aprl and aprc, for
        ;; large scale and convective precipitation, respectively.
        ;; There is no need to add aprs, otherwise snow will be double counted.
	xx = qq
    	xx = xx + find_and_read_var("aprc", date, True)
	return(xx)

    end if

end
