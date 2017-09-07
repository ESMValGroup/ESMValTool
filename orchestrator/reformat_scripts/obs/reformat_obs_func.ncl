;;#############################################################################
;; FUNCTIONS FOR THE reformat_obs_*.ncl SCRIPTS
;;#############################################################################
;; General purpose functions called by the reformat_obs_*.ncl scripts.
;;
;; Contents
;;    function yyyyfrac_to_ut
;;    function create_timec
;;    procedure write_nc
;;    procedure write_nc_profile
;;    function process_EBAS_data
;;
;;#############################################################################

load "./interface_scripts/messaging.ncl"
load "./interface_scripts/auxiliary.ncl"

;;#############################################################################
undef("yyyyfrac_to_ut")
function yyyyfrac_to_ut(frac: numeric,
                        option: integer,
                        flag)
;;
;; Arguments
;;    frac: a one dimensional array of UTC decimal dates of the form yyyy.frac,
;;          i.e., year plus fraction of the year
;;    option: determines output format (see "option" of cd_calendar function)
;;    flag: not yet used
;;
;; Return value
;;    The format of the output is dependent upon the value of flag.
;;    In all cases the first dimension is the same size as frac.
;;    See http://www.ncl.ucar.edu/Document/Functions/Built-in/cd_calendar.shtml
;;
;; Description
;;    Converts UTC decimal dates to a form suitable for NCL calendar functions
;;
;; Caveats
;;    Considers leap years, but no leap seconds
;;    -> Not sure, how this is treated in the raw data ...
;;
;; Modification history
;;    20140624-A_gott_kl: written.
;;

local verbosity, funcname, years, YEAR, units, offset, sorted, spy, ileap, \
    iyears, iy, f, f_sec, time, ut
begin

    funcname = "yyyyfrac_to_ut"
    scriptname = "reformat_scripts/obs/reformat_obs_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Seconds per year
    years = toint(frac)
    spy = todouble(years)  ; just to create array
    spy(:) = 365. * 24. * 3600.
    ileap = ind(isleapyear(years))
    if(.not.all(ismissing(ileap))) then
        spy(ileap) = 366. * 24. * 3600.
    end if

    ;; Offset (i.e. previous years) in seconds
    offset = spy  ; just to create array
    iyears = UNIQ(years)
    sorted = new((/3, dimsizes(iyears)/), double)
    sorted(2, :) = todouble(years(iyears))
    sorted(1, :) = todouble(spy(iyears))
    isorted = dim_pqsort(sorted(2, :), 1)
    sorted(2, :) = sorted(2, isorted)
    sorted(1, :) = sorted(1, isorted)
    sorted(0, :) = 0.
    ;; calculate one offset per year
    sorted(0, 0) = 0.
    do i = 1, dimsizes(iyears) - 1
        sorted(0, i) = sorted(0, i - 1) + sorted(1, i - 1)
    end do
    ;; offsets for all time steps (independent of order)
    do i = 0, dimsizes(iyears) - 1
        iy = ind(years.eq.toint(sorted(2, i)))
        offset(iy) = sorted(0, i)
        delete(iy)
    end do

    ;; Determine units
    YEAR = tostring(sorted(2, 0))
    units = "seconds since " + YEAR + "-01-01 00:00:0.0"

    ;; Fractions
    f = todouble(frac) - todouble(years)  ; decimal fraction
    f_sec = todouble(spy) * f  ; fraction in seconds
    time = offset + f_sec

    ;; Add attributes
    time@units = units
    time@calendar = "standard"

    ;; Convert to UT
    ut = cd_calendar(time, flag)

    leave_msg(scriptname, funcname, 4)
    return(ut)
end

;;#############################################################################
undef("create_timec")
function create_timec(y1:integer,
                      y2:integer)
;;
;; Arguments
;;    y1: start year of the time range.
;;    y2: end year of the time range.
;;
;; Return value
;;    A one-dimensional array of size 12*(y2-y1+1).
;;
;; Description
;;    Create a monthly time coordinate for the given time range.
;;
;; Modification history
;;    20140124-A_righ_ma: written.
;;
local verbosity, funcname, yy, mm
begin

    funcname = "create_timec"
    scriptname = "reformat_scripts/obs/reformat_obs_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    out = new(12 * (y2 - y1 + 1), double)
    do yy = y1, y2
        do mm = 1, 12
            out(12 * (yy - y1) + mm - 1) = \
                cd_inv_calendar(yy, mm, 15, 0, 0, 0, TUNITS, 0)
        end do
    end do

    leave_msg(scriptname, funcname, 4)
    return(out)
end

;;#############################################################################
undef("write_nc")
procedure write_nc(outfile:string,
                   name:string,
                   var:numeric,
                   gAtt:logical)
;;
;; Arguments
;;    outfile: the name of the file to be written, including its path.
;;    name: the variable name.
;;    var: the variable field.
;;    gAtt: a logical variable, whose attributes are appended as file
;;          attributes.
;;
;; Description
;;    Writes the given variable to the given NetCDF file, appending also the
;;    provided global attributes.
;;    This is a modified version of the write_output function in
;;    ./reformat_scripts/constants.ncl, to account for more flexible input
;;    arguments.
;;
;; Modification history
;;    20140123-A_righ_ma: written.
;;
local verbosity, funcname, w, gAtt, dim_names, ndims, dim_sizes, dim_types, \
    dim_unlim, ii
begin

    funcname = "write_nc"
    scriptname = "reformat_scripts/obs/reformat_obs_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Open file
    if (isfilepresent(outfile)) then
        system("rm -f " + outfile)
    end if
    w = addfile(outfile, "c")
    setfileoption(w, "DefineMode", True)

    ; Get coordinates
    dim_names = getvardims(var)
    ndims = dimsizes(dim_names)
    dim_sizes = new(ndims, integer)
    dim_types = new(ndims, string)
    dim_unlim = new(ndims, logical)
    do ii = 0, ndims - 1
        dim_sizes(ii) = dimsizes(var&$dim_names(ii)$)
        dim_types(ii) = typeof(var&$dim_names(ii)$)
        dim_unlim(ii) = False
    end do

    ; Time coordinate must be unlimited
    if (any(dim_names.eq."time")) then
        dim_sizes(ind(dim_names.eq."time")) = -1
        dim_unlim(ind(dim_names.eq."time")) = True
    end if

    ; Define dimensions
    filedimdef(w, dim_names, dim_sizes, dim_unlim)
    do ii = 0, ndims - 1
        filevardef(w, dim_names(ii), dim_types(ii), dim_names(ii))
    end do
    filevardef(w, name, "float", dim_names)
    do ii = 0, ndims - 1
        filevarattdef(w, dim_names(ii), var&$dim_names(ii)$)
    end do
    filevarattdef(w, name, var)

    ;; Append global attributes
    fileattdef(w, gAtt)

    ; Write
    setfileoption(w, "DefineMode", False)
    do ii = 0, ndims - 1
        w->$dim_names(ii)$ = (/var&$dim_names(ii)$/)
    end do
    w->$name$ = (/var/)

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("write_nc_profile")
procedure write_nc_profile(outfile: string,
                           name: string,
                           var: numeric,
                           gAtt: logical)
;;
;; Arguments
;;    outfile: the name of the file to be written, including its path.
;;    name: the variable name.
;;    var: the variable field.
;;    gAtt: a logical variable, whose attributes are appended as file
;;          attributes.
;;
;; Description
;;    Writes the given variable to the given NetCDF file, appending also the
;;    provided global attributes.
;;    This is a modified version of the write_output function in
;;    ./reformat_scripts/constants.ncl, to account for more flexible input
;;    arguments.
;;    Designed to write multiple variables for the vertical profiles data.
;;
;; Modification history
;;    20140422-A_righ_ma: written.
;;
local verbosity, funcname, w, coords, cc, jj, locname, locvar, cname
begin

    funcname = "write_nc_profile"
    scriptname = "reformat_scripts/obs/reformat_obs_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Open file
    if (isfilepresent_esmval(outfile)) then
        system("rm -f " + outfile)
    end if
    w = addfile(outfile, "c")
    setfileoption(w, "DefineMode", True)

    ;; Attach global attributes
    fileattdef(w, gAtt)

    ;; Write dimensions
    coords = getvardims(var)
    do cc = 0, dimsizes(coords) - 2  ; skip column
        cname = coords(cc)
        filedimdef(w, cname, dimsizes(var&$cname$), False)
        filevardef(w, cname, typeof(var&$cname$), cname)
    end do

    ;; Write variable
    do jj = 0, dimsizes(var&column) - 1

        ;; Extract given column
        locname = str_sub_str(name + "_" + var&column(jj), "%", "")
        if (isdim(var, "case")) then
            locvar = var(:, :, jj)
        else
            locvar = var(:, jj)
        end if
        if (var&column(jj).eq."N") then
            locvar@units = "1"
        end if

        ;; Define dimensions
        filevardef(w, locname, "float", coords(0: dimsizes(coords) - 2))
        do cc = 0, dimsizes(coords) - 2
            cname = coords(cc)
            filevarattdef(w, cname, locvar&$cname$)
        end do
        filevarattdef(w, locname, locvar)

        ;; Write
        setfileoption(w, "DefineMode", False)
        do cc = 0, dimsizes(coords) - 2
            cname = coords(cc)
            w->$cname$ = (/locvar&$cname$/)
        end do
        w->$locname$ = (/locvar/)
        delete(locvar)
        delete(locname)

    end do

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("process_EBAS_data")
function process_EBAS_data(in_vars[*]: string,
                           in_units[*]:  string,
                           in_matrix[*]: string,
                           in_compon[*]: string,
                           in_column[*]: string,
                           indir[1]:string,
                           st_code:string,
                           y1:integer,
                           y2:integer)
;;
;; Arguments
;;    in_vars: variables standard name.
;;    in_units: variables units in the raw data.
;;    in_matrix: variables matrix in the raw data.
;;    in_compon: variables name in the raw data.
;;    in_column: variables name in the header.
;;    indir: the input directory for raw data.
;;    stcode: the code of the station to be processed (used for
;;            cross-checking)
;;    y1: start year of the considered time interval.
;;    y2: end year of the considered time interval.
;;
;; Return value
;;    A two-dimensional array (time, variable) with the monthly mean time
;;    series of for each of the processed variables.
;;
;; Description
;;    This function is used to process the data from the EBAS database (e.g.,
;;    EANET, EMEP).
;;
;; Caveats
;;    For the time coordinate in the input data, only units of days are
;;    currently accepted.
;;
;; Modification history
;;    20150413-A_righ_ma: improved time selection.
;;    20140124-A_righ_ma: written.
;;
local timec, datec, vID, fID, bn, en, head, hh, cline, syear, smonth, sday, \
    scode, comp, matr, unit, scale, fills, lline, cols, data_col, flag_col, \
    start_col, end_col, data, value, flag, start_time, end_time, mm, sidx, \
    monthind, stday, enday, nd, pt1, pt2
begin

    funcname = "process_EBAS_data"
    scriptname = "reformat_scripts/obs/reformat_obs_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; EBAS flags for valid measurements
    ;; (see http://www.nilu.no/projects/ccc/flags/index.html)
    validflags = (/798, 797, 782, 781, 780, 771, 770, 741, 740, 680, 679, \
                  678, 676, 675, 668, 665, 662, 660, 657, 656, 655, 654, \
                  653, 652, 651, 650, 649, 648, 644, 640, 559, 558, 557, \
                  556, 555, 532, 531, 521, 499, 498, 476, 475, 470, 458, \
                  457, 450, 440, 420, 410, 394, 392, 390, 382, 380, 370, \
                  299, 298, 276, 275, 258, 257, 250, 249, 248, 247, 220, \
                  211, 210, 191, 190, 189, 188, 187, 186, 185, 147, 120, \
                  111, 110, 103, 102, 101, 100, 000/)

    ;; Create time coordinate
    timec = create_timec(y1, y2)
    datec = cd_calendar(timec, -1)

    ;; Create output array
    data_arr = new((/dimsizes(timec), dimsizes(in_vars)/), float)
    data_arr!0 = "time"
    data_arr&time = timec
    data_arr@_FillValue = FILL

    ;; Create a temporary arrays for time averages and weights
    temp_arr = new(dimsizes(timec), float)
    temp_arr!0 = "time"
    temp_arr&time = timec
    temp_arr@_FillValue = FILL
    ndays_arr = new(dimsizes(timec), float)
    ndays_arr!0 = "time"
    ndays_arr&time = timec
    ndays_arr@_FillValue = FILL

    ;; Loop over variables
    do vID = 0, dimsizes(in_vars) - 1

        info_output("  Processing variable " + in_compon(vID), verbosity, 1)

        ;; Initialize
        temp_arr = 0.
        ndays_arr = 0.

        ;; Read file list
        cstr = "find " + indir + " -type f -name '" + \
            st_code + ".*." + in_compon(vID) + "." + in_matrix(vID) + "*.nas'"
        in_files = systemfunc(cstr)
        if (all(ismissing(in_files))) then
            delete(in_files)
            continue
        end if
        in_bnames = systemfunc(cstr + " -exec basename {} \;")
        sy = str_get_cols(str_get_field(in_bnames, 2, "."), 0, 5)
        delete(cstr)

        sqsort(in_files)
        sqsort(in_bnames)

        ;; Check for duplicates
        if (dimsizes(UNIQ(sy)).ne.dimsizes(sy)) then
            info_output("Duplicated data in input files", verbosity, 1)
            do fID = 0, dimsizes(in_files) - 1
                info_output("  " + in_files(fID), verbosity, 1)
            end do
            info_output("Remove duplicated files considering the " + \
                        "following criteria", verbosity, 1)
            info_output("  most recent revision date", verbosity, 1)
            info_output("  most complete time coverage", verbosity, 1)
            info_output("  same instrument in different years", verbosity, 1)
            error_msg("f", scriptname, funcname, \
                      "rerun this station after removing duplicates")
        end if
        delete(sy)

        ;; Loop over input files
        do fID = 0, dimsizes(in_files) - 1

            info_output("   Reading file " + in_bnames(fID), verbosity, 2)

            ;; Read header
            head = readAsciiHead(in_files(fID), "starttime")

            ;; Extract and check starting date
            syear = toint(str_get_field(head(6), 1, " "))
            smonth = toint(str_get_field(head(6), 2, " "))
            sday = toint(str_get_field(head(6), 3, " "))

            ;; Get time units
            utime = str_sub_str(head(8), "file reference point", "")
            if (.not.isStrSubset(utime, "days")) then
                error_msg("f", scriptname, funcname, \
                          "unexpected time units")
            end if
            utime = utime + syear + "-" + smonth + "-" + sday
            delete(syear)
            delete(smonth)
            delete(sday)

            ;; Use first file units as reference
            if (fID.eq.0) then
                ref_utime = utime
            end if

            ;; Check units consistency
            do hh = 0, dimsizes(head) - 1
                if (isStrSubset(head(hh), "Unit:")) then
                    unit = str_squeeze(str_get_field(head(hh), 2, ":"))
                    if (unit .ne. in_units(vID) .and. unit.ne."ug/m3") then
                        error_msg("f", scriptname, funcname, \
                                  "units in the file not as expected " + \
                                  "(" + unit + " vs. " + in_units(vID) + ")")
                    end if
                    delete(unit)
                end if

            end do

            ;; Get specific fill values and scale factors
            scale = tofloat(str_get_field(head(10), 2, " "))
            fills = tofloat(str_get_field(head(11), 2, " "))

            ;; Get column names
            lline = head(dimsizes(head) - 1)
            ncols = str_fields_count(lline, " ")
            cols = new(ncols, string)
            do cc = 0, ncols - 1
                cols(cc) = str_get_field(lline, cc + 1, " ")
            end do
            data_col = min(ind(cols.eq.in_column(vID))) + 1
            flag_col = (ind(cols.eq."flag_" + in_column(vID).or. \
                            cols.eq."flag")) + 1
            start_col = ind(cols.eq."starttime") + 1
            end_col = ind(cols.eq."endtime") + 1
            delete(cols)
            delete(ncols)

            ;; Read data
            data = readAsciiTable(in_files(fID), 1, "string", dimsizes(head))
            delete(head)

            ;; Read data (for the given month)
            value = tofloat(str_get_field(data(:, 0), data_col, " "))
            value@_FillValue = -999.

            ;; Read flag
            flag = toint(1000 * \
                         tofloat(str_get_field(data(:, 0), flag_col, " ")))

            ;; Filter for valid values
            value = where(value.eq.fills, value@_FillValue, value)
            value = where(value.lt.0, value@_FillValue, value)
            do jj = 0, dimsizes(value) - 1
                if (all(flag(jj).ne.validflags)) then
                    value(jj) = value@_FillValue
                end if
            end do
            delete(flag)
            delete(fills)

            ;; Apply scaling
            if (scale.ne.1) then
                value = value * scale
            end if
            delete(scale)

            ;; Get start and end time
            stt = todouble(str_get_field(data(:, 0), start_col, " "))
            stt@units = utime
            ent = todouble(str_get_field(data(:, 0), end_col, " "))
            ent@units = utime
            delete(data)

            ;; Convert to reference time units
            stt = cd_convert(stt, ref_utime)
            ent = cd_convert(ent, ref_utime)

            ;; Create time series
            if (fID.eq.0) then
                start_time = stt
                end_time = ent
                var = value
            else
                tmp = array_append_record(start_time, stt, 0)
                delete(start_time)
                start_time = tmp
                delete(tmp)
                tmp = array_append_record(end_time, ent, 0)
                delete(end_time)
                end_time = tmp
                delete(tmp)
                tmp = array_append_record(var, value, 0)
                delete(var)
                var = tmp
                delete(tmp)
            end if
            delete(stt)
            delete(ent)
            delete(value)

        end do
        delete(in_files)
        delete(in_bnames)

        ;; Check monotonicity
        if (isMonotonic(start_time).ne.1) then
            error_msg("f", scriptname, funcname, \
                      "non-monotonically increasing time-series, possible " + \
                      "duplicated data in input")
        end if

        ;; Calculate monthly mean from the time series
        do vv = 0, dimsizes(var) - 1

            if (ismissing(var(vv))) then
                continue
            end if

            pstart = start_time(vv)
            pend = -1000.d0

            do while (pend.lt.end_time(vv))

                wdate = cd_calendar(pstart, -5)
                wdatec = cd_calendar(pstart, -1)

                ;; Find beginning of next month
                if (wdate(0, 1).eq.12) then
                    wyear = wdate(0, 0) + 1
                    wmonth = 1
                else
                    wyear = wdate(0, 0)
                    wmonth = wdate(0, 1) + 1
                end if
                pend = cd_inv_calendar(wyear, wmonth, 1, 0, 0, 0, ref_utime, 0)

                if (pend.gt.end_time(vv)) then
                    pend = (/end_time(vv)/)
                end if

                didx = ind(wdatec.eq.datec)
                if (wdate(0, 0).lt.y1 .or. wdate(0, 0).gt.y2) then
                    pstart = pend
                    continue
                end if
                nd = tofloat(pend - pstart)
                temp_arr(didx) = temp_arr(didx) + var(vv) * nd
                ndays_arr(didx) = ndays_arr(didx) + nd

                ;; DEBUG+++
                ;;print(start_time(vv) +"  "+end_time(vv) + "  " + "(" + \
                ;;      cd_calendar(start_time(vv), -2) + "-" + \
                ;;      cd_calendar(end_time(vv), -2)+") " + datec(didx) + \
                ;;      " nd="+nd)
                ;; DEBUG---

                pstart = pend

            end do
            delete(pstart)
            delete(pend)

        end do

        delete(var)
        delete(start_time)
        delete(end_time)

        ;; Calculate monthly mean
        temp_arr = where(temp_arr.eq.0, temp_arr@_FillValue, temp_arr)
        ndays_arr = where(ndays_arr.eq.0, ndays_arr@_FillValue, ndays_arr)
        temp_arr = temp_arr / ndays_arr

        ;; Assign to global data array
        idx_nm = ind(.not.ismissing(temp_arr))
        if (all(ismissing(idx_nm))) then
            delete(idx_nm)
            continue
        end if
        do ii = 0, dimsizes(idx_nm) - 1
            data_arr(idx_nm(ii), vID) = temp_arr(idx_nm(ii))
        end do
        delete(idx_nm)

    end do  ; variables

    leave_msg(scriptname, funcname, 4)
    return(data_arr)

end
