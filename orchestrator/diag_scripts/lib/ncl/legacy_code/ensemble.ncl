;#######################################################################
;# GENERAL ROUTINES FOR ENSEMBLE                                       #
;#######################################################################
;# Please consider using of extending existing routines before adding  #
;# new ones. Check the header of each routine for documentation.       #
;#                                                                     #
;# Contents:                                                           #
;#                                                                     #
;#    * function get_time_range: returns variable with time coordinate #
;#                                                                     #
;#######################################################################

;#######################################################################
undef("get_time_range")
function get_time_range(field : numeric)
;;
;; Description:
;;    * builds up time axis either consisting of an individual model's time
;;      range as specified in namelist, or only intersection of all models'
;;      time ranges, depending on range option set in cfg script
;;
;; Return value:
;;    field with time axis
;;
;; Arguments:
;;    field = variable to which time axis to be added
;;
;; Modification history:
;;    * 20140115 moved to legacy_code (use time_operations in statistics.ncl)
;;               by M. Righi (Mattia.Righi@dlr.de)
;;    * 20131217 adjusted from SeaIce_polcon.ncl by D. Senftleben 
;;               (Daniel.Senftleben@dlr.de)
;;    * ???????? written by K.-D. Gottschaldt (Klaus-Dirk.Gottschaldt@dlr.de)
;;
local field, start_year, end_year, var_time, var_ind, var_ann, years,\
    months, years_ind, months_ind, max_year, min_year, yyyy, yyyymm,\
    yyyymm_ind, imon, temp
begin
    if (diag_script_info@range_option .eq. 1) then
        start_year = max(stringtoint(models@start_year))
        end_year = min(stringtoint(models@end_year))
    else if (diag_script_info@range_option .eq. 0) then
        start_year = min(stringtoint(models@start_year))
        end_year = max(stringtoint(models@end_year))
    else
        info_output("fatal: " + diag_script + ": range_option not defined correctly", verbosity, 0)
        status_exit(1)
    end if
    end if

    ;; time axis
    var_time = cd_calendar(field&time, 0)
    years = floattoint(var_time(:, 0))
    max_year = max(years)
    max_year = min((/max_year, end_year/))
    min_year = min(years)
    min_year = max((/min_year, start_year/))
    months = floattoint(var_time(:, 1))
    yyyymm = years * 100 + months
    delete(var_time)

    ;; select period of years
    temp = ind((years .ge. min_year) .and. (years .le. max_year))
    yyyymm_ind = yyyymm(temp)
    delete(yyyymm)
    yyyymm = yyyymm_ind
    delete(yyyymm_ind)  ; ---------------------
    years_ind = years(temp)
    delete(years)
    years = years_ind
    delete(years_ind)  ; ----------------------
    months_ind = months(temp)
    delete(months)
    months = months_ind
    delete(months_ind)  ; ---------------------
    var_ind = field(temp, :, :)
    var = var_ind
    delete(var_ind)  ; ------------------------
    delete(temp)

    ;; select month or annual mean
    if diag_script_info@month .eq. "A" then  ; weighted annual mean
        yyyy = month_to_annual(int2flt(years), 1)
        var_ann = month_to_annual_weighted(yyyymm, var, 1)
    else  ; select a month
        imon = stringtofloat(diag_script_info@month)
        temp = ind(months .eq. imon)
        yyyy = years(temp)
        var_ann = var(temp, :, :)
        delete(temp)
        var_ann!0 = "year"
        var_ann&year = yyyy
    end if

    return(var_ann)
end
