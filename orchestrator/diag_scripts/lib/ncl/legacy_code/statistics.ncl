;#######################################################################
;# GENERAL ROUTINES FOR STATISTICS                                     #
;#######################################################################
;# Please consider using of extending existing routines before adding  #
;# new ones. Check the header of each routine for documentation.       #
;#                                                                     #
;# Contents:                                                           #
;#                                                                     #
;#    * function time_idx: returns the index of the time dimension     #
;#    * function time_weights: returns weights for time averaging      #
;#    * function time_avg: averages the time dimension                 #
;#                                                                     #
;#######################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

;#######################################################################
undef("time_idx")
function time_idx(field:numeric)
;;
;; Description:
;;    * loops over valid time dimensions until one is found in the dimensions
;;      of field
;;
;; Return value:
;;    * the index of the time (first) dimension;  -1, if none found
;;
;; Arguments:
;;    field = to be checked for a valid time dimension
;;
;; Modification history:
;;    * 20140115 moved to legacy_code (use time_operations in statistics.ncl)
;;               by M. Righi (Mattia.Righi@dlr.de)
;;    * 20131028 written by K.-D. Gottschaldt (Klaus-Dirk.Gottschaldt@dlr.de)
;;
local field, valid_time_dimensions, idx, dims
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering time_idx (statistics.ncl)", verbosity, 10)
    valid_time_dimensions = (/"time", "TIME", "Time"/)
    dims = getvardims(field)
    idx = -1  ; Flag
    do i = 0, dimsizes(valid_time_dimensions) - 1
        if(any(dims .eq. valid_time_dimensions(i))) then
            idx = ind(dims .eq. valid_time_dimensions(i))
            return(idx)
        end if
    end do
    info_output(">>>>>>>> Leaving time_idx (statistics.ncl)", verbosity, 10)
    return(idx)
end

;#######################################################################
undef("time_weights")
function time_weights(time)
;;
;; Description:
;;    * checks for valid time units (as recognized by cd_calendar),
;;      then converts to UTC and determines if time is monthly or annual
;;    * otherwise tries to guess if it is monthly or annual
;;    * considers leap years
;;
;; Return value:
;;    * weights for time (float, vector with dimension of time)
;;    * default -> all weights = 1
;;
;; Arguments:
;;    time = time coordinate variable
;;
;; Caveats:
;;    * no proper error handling for invalid time units (yet)
;;
;; Modification history:
;;    * 20140115 moved to legacy_code (use time_operations in statistics.ncl)
;;               by M. Righi (Mattia.Righi@dlr.de)
;;    * 20131028 written by K.-D. Gottschaldt (Klaus-Dirk.Gottschaldt@dlr.de)
;;
local time, utc_date, year, month, day, weights, uniq_y, uniq_m, uniq_d
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering time_weights (statistics.ncl)", \
                verbosity, 10)
    weights = new(dimsizes(time), float)
    weights(:) = 1.
    if (isatt(time, "units")) then
        info_output("function time_weights: time@units = " + time@units, verbosity, 5)
        ;; Proper units
        if (.not. ismissing(str_match_ic(time@units, "since"))) then
            ;; Convert to UTC
            utc_date = cd_calendar(time, 0)
            year   = tointeger(utc_date(:, 0))
            month  = tointeger(utc_date(:, 1))
            day    = tointeger(utc_date(:, 2))
            ;; Evaluate UTC
            uniq_y = dimsizes(UNIQ(year))
            uniq_m = dimsizes(UNIQ(month))
            uniq_d = dimsizes(UNIQ(day))
            ;; Cases
            if (uniq_m .eq. 1 .and. uniq_d .eq. 1) then
                weights = where(isleapyear(year), 366., 365.)  ; Annual values
                info_output("function time_weights: weighting annual values via UTC", verbosity, 5)
            end if
            if (uniq_m .ne. 1 .and. uniq_d .eq. 1) then
                weights = tofloat(days_in_month(year, month))  ; Monthly values
                info_output("function time_weights: weighting monthly values via UTC", verbosity, 5)
            end if
        end if
        ;; Units indicating "month"
        if (any((/"month", "months", "Month", "Months", "MONTH", "MONTHS"/) .eq. time@units)) then
            year = time
            year(:) = 2001  ; Just use some non-leap year
            weights = tofloat(days_in_month(year, time))  ; Monthly values
            info_output("function time_weights: weighting monthly values for non-leap years", verbosity, 5)
        end if
        ;; Units indicating "year"
        if (any((/"year", "years", "Year", "Years", "YEAR", "YEARS"/) .eq. time@units)) then
            weights = where(isleapyear(year), 366., 365.)  ; Annual values
            info_output("function time_weights: weighting annual values via years", verbosity, 5)
        end if
    end if
    if (all(weights .eq. 1.)) then
        print("warning: all weights for time average = 1 (function time_weights)")
    end if
    info_output(">>>>>>>> Leaving time_weights (statistics.ncl)", verbosity, 10)
    return(weights)
end


undef("time_avg")  ; *************************************************
function time_avg(field:numeric)
;;
;; Description:
;;    * checks if there is a time dimension in field
;;    * tries to determine calender for weighting
;;    * averages time dimension (with weighting, if possible)
;;
;; Return value:
;;    * field, with the time dimension averaged and removed
;;    * original field unchanged, if there is no valid time dimension
;;
;; Arguments:
;;    field = numeric field
;;
;; Caveats:
;;    * treatment of missing values not explicitely specified (yet).
;;
;; Modification history:
;;    * 20140115 moved to legacy_code (use time_operations in statistics.ncl)
;;               by M. Righi (Mattia.Righi@dlr.de)
;;    * 20131028 written by K.-D. Gottschaldt (Klaus-Dirk.Gottschaldt@dlr.de)
;;
local field, valid_time_dimensions, idx, dims
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering time_avg (statistics.ncl)", verbosity, 10)
    ;; Check if there is a time dimension in field
    idx = time_idx(field)  ; Function in time_operators.ncl
    if (idx .eq. -1) then
        print("warning: no valid time dimension (function time_avg)")
        return(field)  ; Return original field unchanged
    else
        ;; Weighting factors for time dimension
        dims = getvardims(field)
        time_dim = dims(idx)
        weights = time_weights(field&$time_dim$)  ; Function in time_operators.ncl
        ;; Average over time
        result = dim_avg_wgt_n_Wrap(field, weights, 1, idx)
    end if
    info_output(">>>>>>>> Leaving time_avg (statistics.ncl)", verbosity, 10)
    return(result)
end
