;;#############################################################################
;; FUNCTIONS FOR THE reformat_obs_*.ncl SCRIPTS
;;#############################################################################
;; General purpose functions called by the reformat_obs_*.ncl scripts. This
;; version contains only "local" functions, i.e., functions that do not
;; depend on the ESMValTool framework.
;;
;; Contents
;;    function var_attrib
;;    function check_time_range
;;    function time_attrib
;;    function lat_attrib
;;    function lon_attrib
;;
;;#############################################################################

;; Global variables/constants
TUNITS = "days since 1950-01-01 00:00:00"
FILL = 1.e+20

;;#############################################################################

undef("var_attrib")
function var_attrib(var:numeric, attrs[1]:logical)
;;
;; Arguments
;;    var: the variable's field.
;;    attrs: logical with output variable attributes attached
;;
;; Return value
;;    An array of the same dimensionality of var.
;;
;; Description
;;    Sets standard variable attributes.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, tmp, p1, p2, natt, in_units, out_units, tmp, att
begin
    funcname = "var_attrib"

    ;; Set fill value first
    if(isatt(var, "_FillValue")) then
        var = where(var.eq.var@_FillValue, FILL, var)
    end if
    var@_FillValue = FILL

    ;; Remove attributes
    delete_VarAtts(var, -1)

    ;; Convert to float
    if (typeof(var).ne."float") then
        out = tofloat(var)
        copy_VarCoords(var, out)
    else
        out = var
    end if

    ;; Append attributes
    out@_FillValue = FILL
    att = getvaratts(attrs)
    do ii = 0, dimsizes(att) - 1
        out@$att(dimsizes(att) - 1 - ii)$ = attrs@$att(dimsizes(att) - 1 - ii)$
    end do

    return(out)
end

undef("check_time_range")
function check_time_range(time:numeric)
;;
;; Arguments
;;    time: input time coordinate.
;;
;; Return value
;;    An array of the same dimensionality of time.
;;
;; Description
;;    Checks the consistency of time coordinate and sets the day to 15th of the
;;    month (for monthly means) and the time to 12:00 (for daily means).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150119-A_gott_kl: defaults for YEAR1, YEAR2, FIELD.
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, date, ntime, tt, pos, y1, y2, opt, exp_ntime
begin

    funcname = "check_time_range"

    if(.not. isdefined("YEAR1")) then
        YEAR1 = time@YEAR1
    end if
    if(.not. isdefined("YEAR2")) then
        YEAR2 = time@YEAR2
    end if
    if(.not. isdefined("FIELD")) then
        FIELD = time@FIELD
    end if

    date = cd_calendar(time, 0)
    ntime = dimsizes(time)

    ;; Check monotonicity
    if (isMonotonic(time).ne.1) then
            print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                               + "non-monotonically-increasing time coordinate in input data")
            status_exit(1)
    end if

    ;; Boundaries
    if (date(0, 0) .gt. YEAR1) then
            print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                  + "selected period outside range (" + YEAR1\
                  + "<" + date(0, 0) + ")")
            status_exit(1)
    end if
    if (date(ntime - 1, 0).lt.YEAR2) then
            print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                               + "selected period outside range (" + YEAR2\
                               + ">" + date(ntime - 1, 0) + ")")
            status_exit(1)
    end if

    ;; Monthly mean case
    if (isStrSubset(FIELD, "M")) then

        ;; Check number of timesteps
        exp_ntime = 12 * (YEAR2 - YEAR1 + 1)
        if (ntime .ne. exp_ntime) then
            print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                               + "incorrect number of timesteps in input data ("\
                               + ntime + " != " + exp_ntime + "). Possibly due to "\
                               + "incomplete years in input.")
            status_exit(1)
        end if

        ;; Set time coordinate to 15th of each month at midnight
        do tt = 0, ntime - 1
            yy = YEAR1 + tt / 12
            mm = mod(tt, 12) + 1
            if (date(tt, 0).ne.yy .or. date(tt, 1).ne.mm) then
                print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                                   + "inconsistent time coordinate in input data"\
                                   + yy + sprinti("%0.2i", mm))
                status_exit(1)
            end if
        end do
        date(:, 2) = 15  ; 15th of the month
        date(:, 3:5) = 0  ; at midnight
    end if

    ;; Daily mean case
    if (isStrSubset(FIELD, "D")) then
        calendar = time@calendar

        ;; Special treatment for proleptic_gregorian calendars
        ;; (not supported by 'days_in_month(...)'
        if (calendar .eq. "proleptic_gregorian") then
            if (YEAR1 .lt. 1582) then
                print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                                   + "changing calendar from proleptic_gregorian "\
                                   + "to gregorian prior to 1582")
            status_exit(1)
            end if
            calendar = "gregorian"
        end if

        ;; Check number of timesteps
        y1 = YEAR1
        y1@calendar = calendar
        y2 = YEAR2
        y2@calendar = calendar
        opt = 0
        opt@calendar = calendar

        exp_ntime = cd_inv_calendar(y2, 12, days_in_month(y2, 12), \
                                    0, 0, 0, TUNITS, opt) - \
            cd_inv_calendar(y1, 1, 1, 0, 0, 0, TUNITS, opt) + 1
        if (ntime.ne.exp_ntime) then
            print("fatal: in " + funcname + "(reformat_obs_func-local.ncl), "\
                               + "incorrect number of timesteps in input data ("\
                               + ntime + " != " + exp_ntime + "). Possibly due to "\
                               + "incomplete years in input.")
            status_exit(1)
        end if

        ;; Set time coordinate to midnight
        ;; TO BE IMPLEMENTED

    end if

    ;; Instantaneous case
    ;; Nothing to check

    ;; Set calendar
    opt = 0
    if (isatt(time, "calendar")) then
        opt@calendar = time@calendar
    end if

    out = cd_inv_calendar(toint(date(:, 0)), toint(date(:, 1)), \
                          toint(date(:, 2)), toint(date(:, 3)), \
                          toint(date(:, 4)), toint(date(:, 5)), \
                          TUNITS, opt)

    return(out)

end

undef("time_attrib")
function time_attrib(coord:numeric)
;;
;; Arguments
;;    coord: input time coordinate.
;;
;; Return value
;;    An array of the same dimensionality of coord.
;;
;; Description
;;    Assigns attributes to the time coordinate (including range check).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, tmp, cal
begin

    funcname = "time_attrib"

    if (isatt(coord, "calendar")) then
        cal = coord@calendar
    else
        cal = "standard"
    end if

    tmp = check_time_range(coord)

    out_coord = todouble(tmp)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@calendar = cal
    out_coord@long_name = "time"
    out_coord@axis = "T"
    out_coord@units = TUNITS
    out_coord@standard_name = "time"
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    return(out_coord)

end

undef("lat_attrib")
function lat_attrib(coord)
;;
;; Arguments
;;    coord: input latitude coordinate.
;;
;; Return value
;;    An array of the same dimensionality of coord.
;;
;; Description
;;    Assigns attributes to the latitude coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, out_coord
begin

    funcname = "lat_attrib"

    out_coord = todouble(coord)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@long_name = "latitude"
    out_coord@axis = "Y"
    out_coord@units = "degrees_north"
    out_coord@standard_name = "latitude"
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    return(out_coord)

end

;;#############################################################################
undef("lon_attrib")
function lon_attrib(coord)
;;
;; Arguments
;;    coord: input longitude coordinate.
;;
;; Return value
;;    An array of the same dimensionality of coord.
;;
;; Description
;;    Assigns attributes to the longitude coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, out_coord
begin

    funcname = "lon_attrib"

    out_coord = todouble(coord)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@long_name = "longitude"
    out_coord@axis = "X"
    out_coord@units = "degrees_east"
    out_coord@standard_name = "longitude"
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    return(out_coord)

end
