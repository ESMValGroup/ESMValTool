;;#############################################################################
;; FUNCTIONS FOR THE GFDL REFORMAT TOOL
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function check_inpath
;;    function find_name
;;    function find_and_read_var
;;    function time_mean
;;    function extract_var
;;    function format_coord
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"

load "./reformat_scripts/constants.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"

;;#############################################################################
undef("check_inpath")
procedure check_inpath()
;;
;; Argument
;;
;; Description
;;    Checks the input path.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, subdir, ii, pp, len_inpath, last_char
begin

    funcname = "check_inpath"
    scriptname = "reformat_scripts/GFDL/reformat_GFDL_func.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Simulation path
    if (.not.isfilepresent_esmval(INPATH)) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, "model data " + \
                  "path (INPATH) not found")
    end if

    ;; Get rid of $ in path strings
    INPATH = systemfunc("echo " + INPATH)

    ;; Add "/" to paths if needed
    len_inpath = strlen(INPATH)
    last_char = str_get_cols(INPATH, len_inpath - 1, len_inpath - 1)
    if (last_char .ne. "/") then
        INPATH = INPATH + "/"
    end if

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("find_name")
function find_name(var:string)
;;
;; Arguments
;;    var: input variable.
;;
;; Return value
;;    The GFDL name of the selected variable.
;;
;; Description
;;    Finds the corresponding GFDL name of the input variable var.
;;
;; Caveats
;;    Needs an external input table "names_GFDL.dat"
;;
;; Reference
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, fname, data, std, emac, type, vert, id
begin

    funcname = "find_name"
    scriptname = "reformat_scripts/GFDL/reformat_GFDL_func.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    fname = "./reformat_scripts/GFDL/names_GFDL.dat"

    data = readAsciiTable(fname, 1, "string", \
                          dimsizes(readAsciiHead(fname, "STANDARD")) + 1)

    std  = str_squeeze(str_get_field(data(:, 0), 1, "|"))
    emac = str_squeeze(str_get_field(data(:, 0), 2, "|"))
    type = str_squeeze(str_get_field(data(:, 0), 3, "|"))
    vert = str_squeeze(str_get_field(data(:, 0), 4, "|"))

    id = ind(std.eq.var)

    if (ismissing(id)) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                  "cannot find GFDL name for " + var)
    else
        out = emac(id)
        out@type = type(id)
        out@vert = vert(id)
    end if

    leave_msg(scriptname, funcname, 4)
    return(out)

end

;;#############################################################################
undef("find_and_read_var")
function find_and_read_var(var:string,
                           err:logical)
;;
;; Arguments
;;    var: input variable.
;;    err: logical.
;;
;; Return value
;;    An array with the selected variable.
;;
;; Description
;;    Reads the variable after finding the GFDL file containing it.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20150708-A_righ_ma: written.
;;
local funcname,  scriptname, verbosity, flist, f
begin

    funcname = "find_and_read_var"
    scriptname = "reformat_scripts/GFDL/reformat_GFDL_func.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Define files list
    ;; Naming convention: <REALM>.<YYYYMM>-<YYYYMM>.<VAR>.nc
    flist = systemfunc("find " + INPATH + " -type f -name '" + REALM + "." + \
                       "*." + var + ".nc' 2>/dev/null")

    ;; If all missing: give up
    if (all(ismissing(flist))) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, "no GFDL files " + \
                  "found for variable " + var)
    end if

    ;; Multiple files will be in the same directory, grab
    ;; the first file from the list and check that the title
    ;; attribute matches the expected values

    f = addfile(flist(0), "r")

    ;; Check that experiment name matches
    if (.not.isatt(f, "title")) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, "attribute " + \
                  "'title' not found in " + flist(0))
    end if
    if (f@title.ne.EXPER) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, "data does not " + \
                  "correspond to the required experiment (" + EXPER + ")")
    end if

    undef("f")

; A-laue_ax+
    ; save list of files (to be) read to global array gFileList;
    ; gFileList and gNFiles are defined in reformat_GFDL_main.ncl
    nf = dimsizes(flist)
    if (nf.gt.dimsizes(gFileList)) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, "size of array " + \
              "gFileList too small. Increase size in reformat_GFDL_main.ncl")
    end if
    gFileList(0:nf-1) = flist(0:nf-1)
    gNFiles = nf
; A-laue_ax-

    ;; Read and return data
    f = addfiles(flist, "r")
    xx = f[:]->$var$
    leave_msg(scriptname, funcname, 4)
    return(xx)

end

;;#############################################################################
undef("time_mean")
function time_mean(var:numeric,
                   mean:string)
;;
;; Arguments
;;    var: input variable.
;;    date: type of mean.
;;
;; Return value
;;    The time average of var.
;;
;; Description
;;    Calculate monthly/daily mean of the given input variable if required.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, rank, date, opt, days, ndays, dd, \
    idx, dum
begin

    funcname = "time_mean"
    scriptname = "reformat_scripts/GFDL/reformat_GFDL_func.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    rank = dimsizes(getvardims(var))
    date = cd_calendar(var&time, 0)

    ;; Set calendar
    opt = 0
    opt@calendar = var&time@calendar

    ;; Monthly mean case
    if (mean.eq."M") then
        ;; Assume monthly fields assembled in a 1-year file for ntime=12.
        ;; Alternatively one could evaluate the time coordinate in detail,
        ;; which might be more flexible in the long run.
        if (dimsizes(var&time).gt.1 .and. dimsizes(var&time).ne.12) then
            if (rank.eq.5) then
                out = var(0:0, :, :, :, :)
                out(0, :, :, :, :) = dim_avg_n_Wrap(var, 0)
            else if (rank.eq.4) then
                out = var(0:0, :, :, :)
                out(0, :, :, :) = dim_avg_n_Wrap(var, 0)
            else if (rank.eq.3) then
                out = var(0:0, :, :)
                out(0, :, :) = dim_avg_n_Wrap(var, 0)
            else if (rank.eq.2) then
                out = var(0:0, :)
                out(0, :) = dim_avg_n_Wrap(var, 0)
            else if (rank.eq.1) then
                out = dim_avg_Wrap(var)
            else
                error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                          "invalid_rank")
            end if
            end if
            end if
            end if
            end if
            ;; Set time coordinate
            out&time = cd_inv_calendar(toint(date(0, 0)), toint(date(0, 1)), \
                                       15, 0, 0, 0, var&time@units, opt)
            leave_msg(scriptname, funcname, 4)
            return(out)
        else
            leave_msg(scriptname, funcname, 4)
            return(var)
        end if
    end if

    ;; Daily mean case
    if (mean.eq."D") then
        days = date(UNIQ(date(:, 2)), 2)
        ndays = dimsizes(days)
        if (dimsizes(var&time).gt.ndays) then
            if (rank.eq.5) then
                out = var(0:ndays - 1, :, :, :, :)
                do dd = 0, ndays - 1
                    idx = ind(date(:, 2).eq.days(dd))
                    out(dd, :, :, :, :) = \
                        dim_avg_n_Wrap(var(idx, :, :, :, :), 0)
                end do
            else if (rank.eq.4) then
                out = var(0:ndays - 1, :, :, :)
                do dd = 0, ndays - 1
                    idx = ind(date(:, 2).eq.days(dd))
                    out(dd, :, :, :) = dim_avg_n_Wrap(var(idx, :, :, :), 0)
                end do
            else if (rank.eq.3) then
                out = var(0:ndays - 1, :, :)
                do dd = 0, ndays - 1
                    idx = ind(date(:, 2).eq.days(dd))
                    out(dd, :, :) = dim_avg_n_Wrap(var(idx, :, :), 0)
                end do
            else if (rank.eq.2) then
                out = var(0:ndays - 1, :)
                do dd = 0, ndays - 1
                    idx = ind(date(:, 2).eq.days(dd))
                    out(dd, :, :) = dim_avg_n_Wrap(var(idx, :), 0)
                end do
            else if (rank.eq.1) then
                out = var(0:ndays - 1)
                do dd = 0, ndays - 1
                    idx = ind(date(:, 2).eq.days(dd))
                    out(dd) = dim_avg_Wrap(var(idx))
                end do
            else
                error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                          "invalid rank")
            end if
            end if
            end if
            end if
            end if

            ; Set time coordinate
            dum = new(ndays, integer)
            dum = 0
            out&time = cd_inv_calendar(toint(date(0:ndays - 1, 0)), \
                                       toint(date(0:ndays - 1, 1)), \
                                       toint(days), dum, dum, dum, \
                                       var&time@units, opt)
            leave_msg(scriptname, funcname, 4)
            return(out)
        else
            leave_msg(scriptname, funcname, 4)
            return(var)
        end if
    end if

    ;; Instantaneous output case (nothing to do)
    if (mean.eq."I") then
        leave_msg(scriptname, funcname, 4)
        return(var)
    end if

    error_msg("f", "reformat_GFDL_func.ncl", funcname, "unrecognized " + \
              "time mean option " + mean)

end

;;#############################################################################
undef("extract_var")
if (isfilepresent_esmval(RECIPEFILE)) then
    loadscript(RECIPEFILE)
end if
function extract_var(name:string,
                     vert:string,
                     year1:integer,
                     year2:integer,
                     stname:string)
;;
;; Arguments
;;    name: GFDL name of the variable.
;;    vert: flag for the vertical coordinate (not used).
;;    date: selected date.
;;    stname: standard CMOR name of the variable.
;;
;; Return value
;;    The selected variable reformatted according to vert and FIELD.
;;
;; Description
;;    Extracts the requested variable from the input
;;
;; Caveats
;;    Vertical integration not implemented.
;;
;; Reference
;;
;; Modification history
;;    20150708-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, xx, xz, rank, date, idx1, idx2
begin

    funcname = "extract_var"
    scriptname = "reformat_scripts/GFDL/reformat_GFDL_func.ncl"
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read variable
    if (name.ne."#RECIP") then  ; GFDL variable
        xx = find_and_read_var(name, True)
    else
        if (isfilepresent_esmval(RECIPEFILE)) then
            xx = GFDL_recipe(name)  ; complex variable
        else
            error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                      "Recipe file for " + stname + " not found")
        end if
    end if

    ;; Calculate zonal mean if required
    if (isStrSubset(FIELD, "z")) then
        xz = dim_avg_n_Wrap(xx, ind("lon".eq.getvardims(xp)))
    else
        xz = xx
    end if
    delete(xx)
    rank = dimsizes(dimsizes(xz))

    ;; Extract given time range
    date = cd_calendar(xz&time, -1)
    idx1 = ind(date.eq.(100 * year1 + 1))
    idx2 = ind(date.eq.(100 * year2 + 12))

    if (ismissing(idx1).or.ismissing(idx2)) then
        error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                  "selected time range (" + sprinti("%0.4i", year1) + \
                  "-" + sprinti("%0.4i", year2) + ") is " + \
                  "outside the model time range (" + \
                  sprinti("%0.4i", toint(min(date)/100)) + "-" + \
                  sprinti("%0.4i", toint(max(date)/100)) + ")")
    end if

    if (rank.eq.4) then
        xt = xz(idx1:idx2, :, :, :)
    else if (rank.eq.3) then
        xt = xz(idx1:idx2, :, :)
    else if (rank.eq.2) then
        xt = xz(idx1:idx2, :)
    else
        error_msg("f", "reformat_GFDL_func.ncl", funcname, "Invalid rank")
    end if
    end if
    end if

    leave_msg(scriptname, funcname, 4)
    return(xt)

end

;;#############################################################################
undef("format_coord")
procedure format_coord(var:numeric)
;;
;; Arguments
;;    var: input variable.
;;
;; Description
;;    Converts variable coordinates (time,plev,lat,lon) to the standard format.
;;
;; Caveats
;;
;; Reference
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, rank, ctime, cplev, clat, clon, tmp
begin

    funcname = "format_coord"
    scriptname = "reformat_scripts/GFDL/reformat_GFDL_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    rank = dimsizes(dimsizes(var))

    ;; Time
    if (iscoord(var, "time")) then

        ;; Check calendar and change case to be recognized by NCL
        if (isatt(var&time, "calendar")) then
            if (var&time@calendar.eq."NOLEAP") then
                var&time@calendar = "noleap"
            end if
        end if
        if (isatt(var&time, "calendar")) then
            if (var&time@calendar.ne."standard".and.\
                var&time@calendar.ne."gregorian".and.\
                var&time@calendar.ne."noleap") then
                error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                          "calendar type " + var&time@calendar + \
                          " not supported")
            end if
        end if
        ctime = time_attrib(var&time)
        delete(var&time)
        var&time = ctime

        ;; Add shift to the current time coordinate
        if (SHIFT_YEAR.ne.0) then
            stdate = cd_calendar(var&time(0), 0)
            stdate(0, 0) = stdate(0, 0) + SHIFT_YEAR
            opt = 0
            opt@calendar = stdate@calendar
            sttime = cd_inv_calendar(toint(stdate(0, 0)), \
                                     toint(stdate(0, 1)), \
                                     toint(stdate(0, 2)), \
                                     toint(stdate(0, 3)), \
                                     toint(stdate(0, 4)), \
                                     stdate(0, 5), TUNITS, opt)
            shift = sttime
            shift = sttime - var&time(0)
            var&time = var&time + shift
            delete(stdate)
            delete(sttime)
            delete(opt)
            delete(shift)
        end if

    end if

    ;; Level
    if (iscoord(var, "level")) then
        var!1 = "plev"
        var&plev = var&plev * 100  ; [hPa] --> [Pa]
        ;; Check ordering (must be monotonically decreasing)
        if (isMonotonic(var&plev) .eq. 0) then
            error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                      "non-monotonic vertical coordinate")
        end if
        if (isMonotonic(var&plev).eq.1) then
            pos = ind("plev".eq.getvardims(var))
            if (rank.eq.5 .and. pos.eq.1) then
                var = var(:, ::-1, :, :, :)
            else if (rank.eq.4 .and. pos.eq.1) then
                var = var(:, ::-1, :, :)
            else if (rank.eq.3 .and. pos.eq.1) then
                var = var(:, ::-1, :)
            else
                error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                          "invalid rank or plev position")
            end if
            end if
            end if
        end if
        cplev = plev_attrib(var&plev)
        delete(var&plev)
        var&plev = cplev
    end if

    ;; Latitude
    if (iscoord(var, "lat")) then
        ;; Reverse order [N-S] --> [S-N]
        if (var&lat(0).gt.var&lat(1)) then
            pos = ind("lat".eq.getvardims(var))
            if (rank.eq.4 .and. pos.eq.2) then
                var = var(:, :, ::-1, :)  ; time, lev, lat, lon
            else if (rank.eq.3 .and. pos.eq.1) then
                var = var(:, ::-1, :)  ; time, lat, lon
            else
                error_msg("f", "reformat_GFDL_func.ncl", funcname, \
                          "invalid rank or lat position")
            end if
            end if
        end if
        clat = lat_attrib(var&lat)
        delete(var&lat)
        var&lat = clat
    end if

    ;; Longitude
    if (iscoord(var, "lon")) then
        clon = lon_attrib(var&lon)
        delete(var&lon)
        var&lon = clon
    end if

    leave_msg(scriptname, funcname, 4)

end
