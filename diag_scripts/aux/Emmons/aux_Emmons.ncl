;#######################################################################
;# ROUTINES ONLY USED BY Emmons.ncl                                    #
;#######################################################################
;# Contents:                                                           #
;#                                                                     #
;#   * function toy_weights: creates  weights for the time axis of     #
;#                      model data, reflecting the distribution of obs #
;#   * function get_period: converts a date string specific to Emmons  #
;#                           data into an usable vector                #
;#   * function poi: checks how data period matches period of interest #
;#                                                                     #
;#######################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"


;#######################################################################
undef("toy_weights")
function toy_weights(start[2], final[2], time:numeric)
;;
;; Description:
;;  * toy_weights = time of year weights
;;  * creates a mask for the obs period
;;    with daily accuracy from multiple year data
;;  * if coarser data are given, weights reflect the number
;;    of matching obs days in each month or season
;;  * annual data are just weighted according to dpy
;;
;; Return value:
;;  * vector of dimsizes(t), with one weighting value for each element of t
;;  * type float, but values are essentially # of days, i.e. integers
;;
;; Arguments: (will be converted to long)
;;  * start:  (/start month, start day/); will be converted to integer
;;  * final:  (/final month, final day/); will be converted to integer
;;  * time: time axis of the model data
;;
;; Caveats:
;;  * processing of daily data not tested yet
;;
;; Modification history:
;;  * 20140211 written (Klaus-Dirk.Gottschaldt@dlr.de)
;;             for use with Emmons.ncl
;;
local m0, m1, d0, d1, wgt, ycal, mcal, dcal, doy, doy0, doy1,\
    wgt_366, wgt_365, months, wgt_mon, wgt_yea, wgt_sea, l_wgt, \
    t, dpm, utc
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering toy_weights (aux_Emmons.ncl)", \
                verbosity, 6)

    ;; convert period input to integer
    m0 = toint(start(0))
    d0 = toint(start(1))
    m1 = toint(final(0))
    d1 = toint(final(1))

    ;; mask observational period on a daily basis
    ;; leap year
    wgt_366 = new(366, float)
    doy0 = day_of_year(2000, m0, d0)
    doy1 = day_of_year(2000, m1, d1)
    info_output("d0 = " + d0 + "   m0 = " + m0 + "    doy0 = " + doy0 +\
                "    d1 = " + d1 + "   m1 = " + m1 + "    doy1 = " +\
                doy1, verbosity, 5)
    if(doy0.le.doy1) then  ; new year not included
        wgt_366(:) = 0.
        wgt_366(doy0 - 1:doy1 - 1) = 1.
    else  ; campaign covers new year's eve
        wgt_366(:) = 1.
        wgt_366(doy0 - 1:doy1 - 1) = 0.
    end if
    ;; non-leap year
    wgt_365 = new(365, float)
    wgt_365(0:58) = wgt_366(0:58)
    wgt_365(59:364) = wgt_366(60:365)

    ;; weights reflecting the monthly distribution of obs
    months = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    dpm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
          ;; neglect leap years here
    months = where(months.lt.m0, 0, months)
    months = where(months.gt.m1, 0, 1)
    wgt_mon = tofloat(dpm * months)
    wgt_mon(m0 - 1) = wgt_mon(m0 - 1) - tofloat(d0) + 1.
    wgt_mon(m1 - 1) = tofloat(d1)
    ;wgt_mon = wgt_mon / sum(wgt_mon)  ; check if needed

    ;; weights reflecting the seasonal distribution of obs,
    ;; for use with seasonal values indicated by an arbitrary
    ;; month of the corresponding standard seasonal range
    wgt_sea = wgt_mon
    wgt_sea((/ 0, 1, 11/)) = sum(wgt_mon((/ 0, 1, 11/)))  ; DJF
    wgt_sea((/ 2, 3,  4/)) = sum(wgt_mon((/ 2, 3,  4/)))  ; MAM
    wgt_sea((/ 5, 6,  7/)) = sum(wgt_mon((/ 5, 6,  7/)))  ; JJA
    wgt_sea((/ 8, 9, 10/)) = sum(wgt_mon((/ 8, 9, 10/)))  ; SON

    ;; translate model time axis to yyyy/mm/dd
    utc = calendar_decode2(time, -5)  ; utc = (/ntime/,/y,m,d,h,m,s/)
    ycal = utc(:, 0)  ; years
    mcal = utc(:, 1)  ; months
    dcal = utc(:, 2)  ; days

    ;; determine time step length in the data, then branch accordingly
    dt = tstep(utc)  ; interface_scripts/auxiliary.ncl

    ;; generate weights: selection of the data period / weighting for
    ;; time averaging will be done via a vector of weights that has one
    ;; one entry for each time step
    wgt = new(dimsizes(time), float)  ; tstep weighting
    if(dt.eq."yea") then  ; time average all data
        info_output("warning: Annual data do not allow accurate matching",\
                    verbosity, 2)
        info_output("         of the observational period. Continuing " +\
                    "anyway. ", verbosity, 2)
        info_output("         Please interpret results with caution.",\
                    verbosity, 2)
        l_wgt = isleapyear(ycal)  ; consider leap years
        wgt = where(l_wgt, 366., 365.)
    else if(dt.eq."sea") then
        info_output("warning: Seasonal data do not allow accurate matching ",\
                    verbosity, 2)
        info_output("         of the observational period. Continuing " +\
                    "anyway. ", verbosity, 2)
        info_output("         Please interpret results with caution.",\
                    verbosity, 2)
        do i = 0, dimsizes(mcal) - 1
            wgt(i) = wgt_sea(mcal(i) - 1)
        end do
    else if(dt.eq."mon") then
        info_output("warning: Monthly data do not allow accurate matching",\
                    verbosity, 2)
        info_output("         of the observational period. Continuing " +\
                    "anyway. ", verbosity, 2)
        info_output("         Please interpret results with caution.",\
                    verbosity, 2)
        do i = 0, dimsizes(mcal) - 1
            wgt(i) = wgt_mon(mcal(i) - 1)
        end do
    else  ; dt = d,h,m,s -> select period based on days
        ;; loop over data time steps & assign weighting based on day_of_year
        do i = 0, dimsizes(dcal) - 1
            doy = day_of_year(ycal(i), mcal(i), dcal(i))
            if(isleapyear(ycal(i))) then
                wgt(i) = wgt_366(doy - 1)
            else
                wgt(i) = wgt_365(doy - 1)
            end if
        end do
    end if
    end if
    end if

    info_output("Model data time step dt = " + dt, verbosity, 5)
    info_output(">>>>>>>> Leaving toy_weights (aux_Emmons.ncl)", \
                verbosity, 6)
    return(wgt)
end

;#######################################################################
undef("poi")
function poi(yyyy_dat[2], yyyy_poi[2])
;;
;; Description:
;;  * poi = period of interest
;;  * checks if the period described by yyyy_dat
;;      / is contained in /
;;      / has its center in /
;;      / overlaps with /
;;      / encloses /
;;      / is outside /
;;    the years described by yyyy_poi
;;
;; Return value:
;;  * string array containing all true statements from the following list:
;;    contained / center / overlap / encloses / outside
;;    (e.g. "contained" implies "center" and "overlap")
;;
;; Arguments: (will be converted to long)
;;    yyyy_poi = (/year min,year max/) of the period of interest
;;    yyyy_dat = (/year min,year max/) of the (data) to check
;;
;; Caveats:
;;
;; Modification history:
;;  * 20140210 written (Klaus-Dirk.Gottschaldt@dlr.de)
;;             for use with Emmons.ncl
;;
local result, result_all, datmin, datmax, poimin, poimax, poi_checks
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering poi (aux_Emmons.ncl)", verbosity, 6)

    ;; convert all input to integer
    datmin = tolong(yyyy_dat(0))
    datmax = tolong(yyyy_dat(1))
    poimin = tolong(yyyy_poi(0))
    poimax = tolong(yyyy_poi(1))

    ;; check periods (dat with respect to poi)
    poi_checks = new(4, logical)  ; holds results of the 4 checks
    poi_checks = (/datmin.le.poimin, datmin.le.poimax,\
                  datmax.le.poimin, datmax.le.poimax/)
    if (all(poi_checks.eq.(/False, True, False, True/))) then
        result = "contained"
    else if (all(poi_checks.eq.(/True, True, False, False/))) then
        result = "enclosing"
    else if (all(poi_checks.eq.True) .or. all(poi_checks.eq.False)) then
        result = "outside"
        if(datmin.eq.poimin .or. datmax.eq.poimax) then  ; catch datmin=datmax
            result = "overlap"  ; cautious choice
        end if
    else
        result = "overlap"
    end if
    end if
    end if

    ;; check "overlap" if the majority of poi_dat is within poi_roi
    if (result.eq."overlap") then
        a1 = datmax - datmin  ; length of poi_dat
        a3 = min((/datmax, poimax/)) - max((/datmin, poimin/))
             ;; overlap with poi_roi
        if (2. * a3.gt.a1) then result = "center" end if  ; compare
    end if

    ;; add other true, but less strong statements to the output vector
    if (result.eq."contained") then
        result_all = (/"contained", "center", "overlap"/)
    else if (result.eq."center") then
        result_all = (/"center", "overlap"/)
    else if (result.eq."encloses") then
        result_all = (/"encloses", "center", "overlap"/)
    else if (result.eq."overlap") then
        result_all = (/"overlap"/)
    else
        result_all = (/"outside"/)
    end if
    end if
    end if
    end if

    info_output("  poi = " + yyyy_poi(0) + "-" + yyyy_poi(1) +\
                " ; dat = " + yyyy_dat(0) + "-" + yyyy_dat(1), verbosity, 5)
    info_output(">>>>>>>> Leaving poi (aux_Emmons.ncl)", verbosity, 6)
    return(result_all)
end

;#######################################################################
undef("get_period")
function get_period(period_in[1]:string)
;;
;; Description:
;;  * splits the input string into yyyy, mm, dd
;;  * for start and end time
;;
;; Return value:
;;  * string array (3,2)
;;
;; Arguments:
;;  * period_in: string of the form "yyyy/mm/dd-yyyy/mm/dd"
;;
;; Caveats:
;;  * very specific
;;
;; Modification history:
;;  * 20140131 written (Klaus-Dirk.Gottschaldt@dlr.de)
;;             for use with Emmons.ncl
;;
local result
begin

    funcname = "get_period"

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_period (aux_Emmons.ncl)", \
                verbosity, 6)

    ;; test if period_in is in the expected format
    n1 = str_index_of_substr(period_in, "-", 0)
    n2 = str_index_of_substr(period_in, "/", 0)
    if (dimsizes(n1).ne.1 .or. dimsizes(n2).ne.4) then
        error_msg("f", "aux_Emmons.ncl", funcname, "period_in must be in " + \
                  "the form yyyy/mm/dd-yyyy/mm/dd, but is " + \
                  period_in)
    else

    ;; split period_in and put elements into result array
        start_final = str_split(period_in, "-")
        start = str_split(start_final(0), "/")
        final = str_split(start_final(1), "/")
        result = new((/3, 2/), string)
        result(0, 0) = start(0)  ; start year
        if(tolong(start(1)).le.12) then
            result(1, 0) = start(1)  ; start month
        else
            error_msg("f", "aux_Emmons.ncl", funcname, "start date must " + \
                      "have the form yyyy/mm/dd, but is " + start_final(0) + \
                      ": mm = " + start(1))
        end if
        if(tolong(start(2)).le.31) then
            result(2, 0) = start(2)  ; start day
        else
            error_msg("f", "aux_Emmons.ncl", funcname, "start date must " + \
                      "have the form yyyy/mm/dd, but is " + start_final(0) + \
                      ": mm = " + start(2))
        end if
        result(0, 1) = final(0)  ; final year
        if(tolong(final(1)).le.12) then
            result(1, 1) = final(1)  ; final month
        else
            error_msg("f", "aux_Emmons.ncl", funcname, "start date must " + \
                      "have the form yyyy/mm/dd, but is " + start_final(1) + \
                      ": mm = " + final(1))
        end if
        if(tolong(final(2)).le.31) then
            result(2, 1) = final(2)  ; final day
        else
            error_msg("f", "aux_Emmons.ncl", funcname, "start date must " + \
                      "have the form yyyy/mm/dd, but is " + start_final(1) + \
                      ": mm = " + final(2))
        end if
    end if
    info_output(">>>>>>>> Leaving get_period (aux_Emmons.ncl)", \
                verbosity, 6)
    return(result)
end
