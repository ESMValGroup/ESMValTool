;;#############################################################################
;; FUNCTIONS FOR THE DEFAULT REFORMAT ROUTINES
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure check_inpath
;;    function find_alt_names
;;    function check_var
;;    procedure apply_fixes
;;    function read_var
;;    function check_rank
;;    function find_grid_type
;;    procedure check_fill
;;    procedure reformat_time
;;    procedure reformat_plev
;;    procedure reformat_lev
;;    procedure reformat_basin
;;    procedure reformat_lat
;;    function reformat_lon
;;    function reformat_coord
;;    procedure add_irregulargrid_info
;;    procedure add_lmask_info
;;    procedure add_omask_info
;;    procedure add_porosity_info
;;    procedure add_depth_info
;;    procedure add_regionbasin_info
;;
;;#############################################################################

load "./reformat_scripts/constants.ncl"
load "./interface_scripts/logging.ncl"
load "interface_scripts/interface.ncl"

begin
    if (fileexists(FIXFILE)) then
        loadscript(FIXFILE)
    end if
end

;;#############################################################################
undef("check_inpath")
procedure check_inpath()
;;
;; Arguments
;;
;; Description
;;    Checks the input path.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname
begin

    funcname = "check_inpath"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Simulation path
    if (.not. fileexists(INPATH)) then
        error_msg("f", scriptname, funcname, \
                  "simulation path INPATH = " + INPATH + ", not found")
    end if

    ;; Get rid of $ in path strings
    INPATH = systemfunc("echo " + INPATH)

    ;; Add "/" to paths if needed
    if (str_get_cols(INPATH, -1, -1).ne."/") then
        INPATH = INPATH + "/"
    end if
    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("find_alt_names")
function find_alt_names(name:string)
;;
;; Arguments
;;    name: the variable's standard name.
;;
;; Return value
;;    An array of strings with alternative names.
;;
;; Description
;;    Returns possible alternative names for the given variable.
;;    Needs an external input table F_VARIABLES.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, locname, data, col1, col2, std, alt, \
    id, tmp, ncols
begin

    funcname = "find_alt_names"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Create local copy
    locname = name

    ;; Read standard names file
    data = readAsciiTable(F_VARIABLES, 1, "string", 0)

    ;; Parse file
    col1 = str_squeeze(str_get_field(data(:, 0), 1, "="))
    col2 = str_squeeze(str_get_field(data(:, 0), 2, "="))
    std = col2(ind(col1.eq."std_name"))
    alt = col2(ind(col1.eq."alt_name"))
    id = ind(std.eq.name)
    if (ismissing(id)) then
        error_msg("f", scriptname, funcname, \
                  "variable " + name + " not in " + F_VARIABLES)
    else
        tmp = alt(id)
    end if

    ncols = str_fields_count(tmp, ", ")

    ;; Return alternative names
    if (ncols .ne. 0) then
        out = new(ncols, string)
        do ii = 0, ncols - 1
            out(ii) = str_squeeze(str_get_field(tmp, ii + 1, ", "))
        end do
        leave_msg(scriptname, funcname)
        return(out)
    else
        leave_msg(scriptname, funcname)
        return(locname)
    end if
    return(return_val)

end

;;#############################################################################
undef("check_var")
function check_var(fo:file, name:string)
;;
;; Arguments
;;    fo: an object for the input file.
;;    name: the variable's standard name.
;;
;; Return value
;;    A string indicating the name of the requested variable in the file.
;;
;; Description
;;    Checks if the given variable is in the file, check for alternative names
;;    and returns the actual variable name in the file.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: modified to return a missing if not found.
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, altnames, ii
begin

    funcname = "check_var"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Return value initialized to missing
    ret_name = new(1, string)

    if (isfilevar(fo, name)) then
        ret_name = name
    else
        ;; Try the alternatives
        altnames = find_alt_names(name)
        do ii = 0, dimsizes(altnames) - 1
            if (isfilevar(fo, altnames(ii)) .and. \
                ismissing(ret_name)) then
                log_info("  Variable " + name + \
                            " found with alternative name " + altnames(ii))
                ret_name = altnames(ii)
            end if
        end do
    end if

    leave_msg(scriptname, funcname)
    return(ret_name)

end

;;#############################################################################
undef("apply_fixes")
procedure apply_fixes(var:numeric, name:string, fixfile:string)
;;
;; Arguments
;;    var: the variable's field.
;;    name: the variable's standard name.
;;
;; Description
;;    Looks for a model-specific fixing routine and applies it to the variable.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20131021-A_righ_ma: written.
;;
local funcname, scriptname, ret
begin

    fixfile = ""  ; filename of "FIXFILE"

    funcname = "apply_fixes"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    if (fileexists(FIXFILE)) then
        ret = fixing_func(var, name)
        if (ret.eq.0) then
            log_info("  A data-fixing function (" + FIXFILE + ") " + \
                        "has been applied to this model/variable!")
            fixfile = FIXFILE
        end if
    end if
    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("read_var")
function read_var(name:string, filelist:string, fixfile:string)
;;
;; Arguments
;;    name: the variable's standard name.
;;    filelist: list of filenames processed to produce the file to be written
;;    fixfile: filename of a FIXFILE applied to process the data to be written
;;             (if applicable)
;;
;; Return value
;;    The variable's field.
;;
;; Description
;;    Reads the selected variable and extract the given time range
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150630-a_laue_ax: saving the following details to the output NetCDF
;;                        files (to be added to the ESMValTool log-files):
;;                        - ESMValTtool version
;;                        - list of input files processed
;;                        - name of the "FIXFILE" applied (if any)
;;    20150415-A_righ_ma: adapted to changes in check_var
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, files, ff, f, altname, tmp, rank, \
    tmp2, tname, opt, date, tcoord, lim_date1, lim_date2, curr_date1, \
    curr_date2, set_date1, set_date2, s1, s2, xx
begin

    funcname = "read_var"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Lists of files
    files = systemfunc("ls " + INFILE)

    ic = 0

    do ff = 0, dimsizes(files) - 1
        f = addfile(files(ff), "r")
        altname = check_var(f, name)
        if (ismissing(altname)) then
            error_msg("f", scriptname, funcname, \
                      "variable " + name + " is not in the file")
        end if
        tmp = f->$altname$

        ;; Check for model-specific fixes to this variable
        apply_fixes(tmp, name, fixfile)

        ;; Define rank
        rank = dimsizes(dimsizes(tmp))

        ;; Convert to float
        if (typeof(tmp) .ne. "float") then
            tmp2 = tofloat(tmp)
            copy_VarMeta(tmp, tmp2)
            delete(tmp)
            tmp = tmp2
            delete(tmp2)
        end if
        tname = tmp!0

        ; In order to support F2Ms fields (lon+lat, no time dimension),
        ; we need the possibility to skip the selection of the time rage
        ; ---> if (year1=0) .and. (year2=0) ---> skip selection of time range
        if ((YEAR1.eq.0).and.(YEAR2.eq.0)) then
            xx = tmp
        else

            ;; Set calendar for date conversion
            opt = 0
            if (isatt(tmp&$tname$, "calendar")) then
                opt@calendar = tmp&$tname$@calendar
            end if

            ;; Set standard units for time coordinate
            if (tmp&$tname$@units .ne. TUNITS) then
                date = cd_calendar(tmp&$tname$, 0)
                tcoord = cd_inv_calendar(toint(date(:, 0)), toint(date(:, 1)),\
                                         toint(date(:, 2)), toint(date(:, 3)),\
                                         toint(date(:, 4)), toint(date(:, 5)),\
                                         TUNITS, opt)
                delete(tmp&$tname$)
                tmp&$tname$ = tcoord
                delete(tcoord)
                delete(date)
            end if

            ;; Given time interval (from namelist)
            lim_date1 = cd_inv_calendar(YEAR1, 1, 1, 0, 0, 0,\
                                        tmp&$tname$@units, opt)
            lim_date2 = cd_inv_calendar(YEAR2, 12, 31, 23, 59, 59,\
                                        tmp&$tname$@units, opt)
            if (isatt(tmp&$tname$, "calendar")) then
                if (any(tmp&$tname$@calendar.eq.(/"360", "360_day"/))) then
                    lim_date2 = cd_inv_calendar(YEAR2, 12, 30, 23, 59, 59,\
                                                tmp&$tname$@units, opt)
                end if
            end if

            ;; Current time interval in this file
            curr_date1 = tmp&$tname$(0)
            curr_date2 = tmp&$tname$(dimsizes(tmp&$tname$) - 1)

            ;; Select time interval
            if (lim_date1 .gt. curr_date2 .or. lim_date2 .lt. curr_date1) then
                delete(tmp)
                continue
            end if
            set_date1 = max((/lim_date1, curr_date1/))
            set_date2 = min((/lim_date2, curr_date2/))

            ;; Set consistent time coordinate types
            if (typeof(tmp&$tname$).eq."float") then
                s1 = tofloat(set_date1)
                s2 = tofloat(set_date2)
                delete(set_date1)
                delete(set_date2)
                set_date1 = s1
                set_date2 = s2
                delete(s1)
                delete(s2)
            end if

            ;; Extract time interval
            if (rank.eq.4) then
                xx = tmp({set_date1:set_date2}, :, :, :)
            else if (rank.eq.3) then
                xx = tmp({set_date1:set_date2}, :, :)
            else if (rank.eq.2) then
                xx = tmp({set_date1:set_date2}, :)
            else if (rank.eq.1) then
                xx = tmp({set_date1:set_date2})
            end if
            end if
            end if
            end if

        end if
        delete(tmp)

        ;; Append to time series
        if (.not. isdefined("outvar")) then
            outvar = xx
        else
            tmp = array_append_record(outvar, xx, 0)
            delete(outvar)
            outvar = tmp
            delete(tmp)
        end if
        filelist(ic) = files(ff)
        ic = ic + 1
        delete(xx)
        delete(f)
        if (isvar("opt")) then
            delete(opt)
        end if

    end do

    if (isdefined("outvar")) then
        leave_msg(scriptname, funcname)
        return(outvar)
    else
        error_msg("f", scriptname, funcname, \
                  "cannot extract variable, check time range")
    end if

end

;;######################################################################
undef("check_rank")
function check_rank(var:numeric)
;;
;; Arguments
;;    var: the variable's field.
;;
;; Return value
;;    The variable's field with dummy coordinates removed.
;;
;; Description
;;    Checks consistency between variable rank and given field type
;;    and removes "dummy" coordinates
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, dims, rank, pos, lerr
begin

    funcname = "check_rank"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    dims = dimsizes(var)
    rank = dimsizes(dims)

    ;; Calculate zonal mean (T3? --> T2?z)
    if (rank.eq.4 .and. \
        isStrSubset(FIELD, "T2") .and. \
        isStrSubset(FIELD, "z")) then
        out = dim_avg_Wrap(var)
        delete(dims)
        delete(rank)
        dims = dimsizes(out)
        rank = dimsizes(dims)
    end if

    ;; Remove "dummy" coordinates (i.e. coordinates of size 1)
    pos = ind(dims.eq.1)
    if (dimsizes(pos) .gt. 1) then
        error_msg("f", scriptname, funcname, \
                  "multiple dummy coordinates not supported")
    end if

    if (.not. ismissing(pos)) then
        if (pos.eq.0) then
            continue  ; ignore time coordinate
        else
            ;; Rank 4
            if (rank.eq.4) then
                if (pos.eq.3) then
                    out = var(:, :, :, 0)
                else if (pos.eq.2) then
                    out = var(:, :, 0, :)
                else if (pos.eq.1) then
                    out = var(:, 0, :, :)
                else
                    error_msg("f", scriptname, funcname, \
                              "invalid dimension position")
                end if
                end if
                end if
            ;; Rank 3
            else if (rank.eq.3) then
                if (pos.eq.2) then
                    out = var(:, :, 0)
                else if (pos.eq.1) then
                    out = var(:, 0, :)
                else
                    error_msg("f", scriptname, funcname, \
                              "invalid dimension position")
                end if
                end if
            ;; Rank 2
            else if (rank.eq.2) then
                if (pos.eq.1) then
                    out = var(:, 0)
                else
                    error_msg("f", scriptname, funcname, \
                              "invalid dimension position")
                end if
            else
                error_msg("f", scriptname, funcname, \
                          "invalid rank")
            end if
            end if
            end if
        end if
    end if

    ;; Redefine rank
    if (isdefined("out")) then
        rank = dimsizes(dimsizes(out))
    else
        out = var
    end if

    lerr = False
    do ii = 0, 3
        if (isStrSubset(FIELD, "T" + ii) .and. rank .ne. ii + 1) then
            lerr = True
        end if
    end do

    if (lerr) then
        error_msg("f", scriptname, funcname, \
                  "variable rank (" + rank + ") " + \
                  "not consistent with given type (" + FIELD + ")")
    end if

    leave_msg(scriptname, funcname)
    return(out)


end

;;######################################################################
undef("find_grid_type")
function find_grid_type()
;;
;; Arguments
;;
;; Return value
;;    The grid type as string ("rectilinear" or "irregular").
;;
;; Description
;;    Identifies the grid for the given variable.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, ff, altname, lat, lon
begin

    funcname = "find_grid_type"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)
    files = systemfunc("ls " + INFILE)

    ff = addfile(files(0), "r")

    ;; Check if latitude is present, if not assume rectilinear grid
    altname = check_var(ff, "lat")
    if (ismissing(altname)) then
        out = "rectilinear"
        leave_msg(scriptname, funcname)
        return(out)
    end if
    lat = ff->$altname$
    delete(altname)

    ;; Check if longitude is present, if not assume rectilinear grid
    altname = check_var(ff, "lon")
    if (ismissing(altname)) then
        out = "rectilinear"
        leave_msg(scriptname, funcname)
        return(out)
    end if
    lon = ff->$altname$
    delete(altname)

    ;; If lat and lon are 2D arrays, assume irregular grid
    if (dimsizes(dimsizes(lat)).ge.2.and.dimsizes(dimsizes(lon)).ge.2) then
        out = "irregular"
        leave_msg(scriptname, funcname)
        return(out)
    end if

    ;; If lat and lon are 1D arrays, assume rectilinear grid
    if (dimsizes(dimsizes(lat)).ge.1.and.dimsizes(dimsizes(lon)).ge.1) then
        out = "rectilinear"
        leave_msg(scriptname, funcname)
        return(out)
    end if

    ;; Case not covered, issue an error
    error_msg("f", scriptname, funcname, \
              "cannot determine the grid type for this variable")

end

;;######################################################################
undef("check_fill")
procedure check_fill(var:numeric)
;;
;; Arguments
;;    var: the variable's field.
;;
;; Description
;;    Checks variable's missing values
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, altfill, ii
begin

    funcname = "check_fill"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    if (.not. isatt(var, "_FillValue")) then
        log_info("  No _FillValue found for this variable, " + \
                    "assigning default")
        var@_FillValue = FILL
    end if

    ;; Possible undefined missing values
    altfill = (/-999., -9999., -99999., -999999., 1.e15, -1.e34, 1.e36/)

    ;; Check for possibly undefined missing values
    do ii = 0, dimsizes(altfill) - 1
        if (any(var.eq.altfill(ii)) .and. \
            var@_FillValue .ne. altfill(ii)) then
            var = where(var.eq.altfill(ii), var@_FillValue, var)
            log_info("  Found probable missing values (" + altfill(ii) + \
                        ") not defined as _FillValue, fixing")
        end if
    end do

    ;; Set standard ESMVal-tool missing value (1.e20)
    var@_FillValue = FILL
    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("reformat_time")
procedure reformat_time(var:numeric,
                        pos:integer)
;;
;; Arguments
;;    var: the variable's field.
;;    pos: position in the array.
;;
;; Description
;;    Reformats the time coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: adapted to new reformat_coord procedure.
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, ctime
begin

    funcname = "reformat_time"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Check positioning
    if (pos.ne.0) then
        error_msg("f", scriptname, funcname, \
                  "time must be the dimension 0")
    end if

    ;; Assign standard name
    var!pos = "time"

    ;; Assign attributes (including range check)
    ctime = time_attrib(var&time)
    delete(var&time)
    var&time = ctime
    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("reformat_plev")
procedure reformat_plev(var:numeric,
                        pos:integer)
;;
;; Arguments
;;    var: the variable's field.
;;    pos: position in the array.
;;
;; Description
;;     Reformats the pressure level coordinate.
;;
;; Caveats
;;    plev is not a good name anymore as level can also be in m.
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: adapted to new reformat_coord procedure.
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, allowed_units, alt_units_hpa, \
    alt_units_bar, rank, cplev
begin

    funcname = "reformat_plev"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    allowed_units = (/"m", "M",\
                      "Pa", "pa", "PA", "Pascals", "pascals", "PASCALS"/)
    alt_units_hpa = (/"hpa", "hPa", "HPA", "hecto-pascals", "HECTO-PASCALS",\
                      "mbar", "mBar", "mb", "Mb", "MB", "millibar",\
                      "millibars", "MILLIBARS"/)
    alt_units_bar = (/"b", "B", "bar", "bars", "BARS"/)

    ;; Check positioning
    if (pos.ne.1) then
        error_msg("f", scriptname, funcname, \
                  "pressure level must be the dimension 1")
    end if

    ;; Assign standard name
    var!pos = "plev"

    ;; Check ordering (must be monotonically decreasing)
    rank = dimsizes(dimsizes(var))
    if (isMonotonic(var&plev).eq.0) then
        error_msg("f", scriptname, funcname, \
                  "non-monotonic pressure level coordinate")
    end if
    if (isMonotonic(var&plev).eq.1) then
        if (rank.eq.4) then
            var = var(:, ::-1, :, :)
        else if (rank.eq.3) then
            var = var(:, ::-1, :)
        else
            error_msg("f", scriptname, funcname, \
                      "invalid rank")
        end if
        end if
    end if

    ;; Check units
    if (all(var&plev@units.ne.allowed_units)) then
        if (any(var&plev@units.eq.alt_units_hpa)) then
            var&plev = var&plev * 100.  ; [hPa] --> [Pa]
        else if (any(var&plev@units.eq.alt_units_bar)) then
            var&plev = var&plev * 1.e5  ; [bar] --> [Pa]
        else
            error_msg("f", scriptname, funcname, \
                      "unrecognized units " + \
                      "(" + var&plev@units + ") in pressure level coordinate")
        end if
        end if
    end if

    ;; Assign attributes
    cplev = plev_attrib(var&plev)
    delete(var&plev)
    var&plev = cplev

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("reformat_lev")
procedure reformat_lev(var:numeric,
                       pos:integer)
;;
;; Arguments
;;    var: the variable's field.
;;    pos: position in the array.
;;
;; Description
;;    Reformats the ocean level coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: written.
;;
local funcname, scriptname, rank, clev
begin

    funcname = "reformat_lev"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Assign standard name
    var!pos = "lev"

    ;; Check standard name
    if (var&lev@standard_name.ne."depth") then
        error_msg("f", scriptname, funcname, \
                  "unexpected standard name (" + var&lev@standard_name + \
                  ") in ocean level coordinate")
    end if

    ;; Check units
    if (var&lev@units.ne."m") then
        error_msg("f", scriptname, funcname, \
                  "unrecognized units " + \
                  "(" + var&lev@units + ") in ocean level coordinate")
    end if

    ;; Check ordering (must be monotonically increasing)
    rank = dimsizes(dimsizes(var))
    if (isMonotonic(var&lev).eq.0) then
        error_msg("f", scriptname, funcname, \
                  "non-monotonic ocean level coordinate")
    end if
    if (isMonotonic(var&lev).eq.-1) then
        if (rank.eq.4) then
            if (pos.eq.2) then
                var = var(:, :, ::-1, :)
            else
                error_msg("f", scriptname, funcname, \
                          "invalid dimension position")
            end if
        else
            error_msg("f", scriptname, funcname, \
                      "invalid rank")
        end if
    end if

    ;; Assign attributes
    clev = lev_attrib(var&lev)
    delete(var&lev)
    var&lev = clev

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("reformat_basin")
procedure reformat_basin(var:numeric,
                         pos:integer)
;;
;; Arguments
;;    var: the variable's field.
;;    pos: position in the array.
;;
;; Description
;;    Reformats the basin coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: written.
;;
local funcname, scriptname, dsizes
begin

    funcname = "reformat_basin"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Assign standard name
    var!pos = "basin"

    ;; Define coordinate
    dsizes = dimsizes(var)
    if (.not.iscoord(var, "basin"))
        var&basin = ispan(1, dsizes(pos), 1)
        var&basin@long_name = "ocean basin"
        var&basin@standard_name = "region"
    end if

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("reformat_lat")
procedure reformat_lat(var:numeric,
                       pos:integer)
;;
;; Arguments
;;    var: the variable's field.
;;    pos: position in the array.
;;
;; Description
;;    Reformats the latitude coordinate
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: adapted to new reformat_coord procedure.
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, rank, clat
begin

    funcname = "reformat_lat"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Assign standard name
    var!pos = "lat"

    ;; Check ordering (must be S-N, -90:90)
    rank = dimsizes(dimsizes(var))
    if (isMonotonic(var&lat).eq.0) then
        error_msg("f", scriptname, funcname, \
                  "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(var&lat).eq.-1) then
        ;; Rank 4
        if (rank.eq.4) then
            if (pos.eq.2) then
                var = var(:, :, ::-1, :)
            else if (pos.eq.3)
                var = var(:, :, :, ::-1)
            else
                error_msg("f", scriptname, funcname, \
                          "invalid dimension position")
            end if
            end if
        ;; Rank 3
        else if (rank.eq.3)
            if (pos.eq.1) then
                var = var(:, ::-1, :)
            else
                error_msg("f", scriptname, funcname, \
                          "invalid dimension position")
            end if
        else
            error_msg("f", scriptname, funcname, \
                      "invalid rank")
        end if
        end if
    end if

    ;; Assign attributes
    clat = lat_attrib(var&lat)
    delete(var&lat)
    var&lat = clat

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("reformat_lon")
function reformat_lon(var:numeric,
                      pos:integer)
;;
;; Arguments
;;    var: the variable's field.
;;    pos: position in the array.
;;
;; Return value
;;    The input variable with reformatted longitude coordinate.
;;
;; Description
;;    Reformats the longitude coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: adapted to new reformat_coord procedure.
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, cvar, rank, tmp, ip, clon
begin

    funcname = "reformat_lon"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Copy to local var
    cvar = var

    ;; Check positioning
    rank = dimsizes(dimsizes(cvar))
    if ((rank - 1).ne.pos) then
        error_msg("f", scriptname, funcname, \
                  "longitude must be the dimension n-1")
    end if

    ;; Assign standard name
    cvar!pos = "lon"

    ;; Check double entries 0 = 360
    if (cvar&lon(0).eq.cvar&lon(dimsizes(cvar&lon) - 1) - 360.) then
        if (pos.eq.3) then
            tmp = cvar(:, :, :, 0:dimsizes(cvar&lon) - 2)
        else if (pos.eq.2) then
            tmp = cvar(:, :, 0:dimsizes(cvar&lon) - 2)
        else
            error_msg("f", scriptname, funcname, \
                      "cannot determine position of longitude coordinate")
        end if
        end if
        delete(cvar)
        cvar = tmp
        delete(tmp)
    end if

    ;; Check for values >= 360
    if (cvar&lon(dimsizes(cvar&lon) - 1) .ge. 360. .or. \
        cvar&lon(dimsizes(cvar&lon) - 1).eq.0.) then
        cvar&lon = where(cvar&lon .ge. 360., cvar&lon - 360., cvar&lon)
        ip = dim_pqsort(cvar&lon, 1)
        tmp = cvar
        if (pos.eq.3) then
            tmp = cvar(:, :, :, ip)
        else if (pos.eq.2) then
            tmp = cvar(:, :, ip)
        else
            error_msg("f", scriptname, funcname, \
                      "cannot determine position of longitude coordinate")
        end if
        end if
        tmp&lon = cvar&lon(ip)
        delete(cvar)
        cvar = tmp
        delete(tmp)
    end if

    ;; Convert -180:180 to 0:360
    if (any(cvar&lon .lt. 0.)) then
        cvar = lonFlip(cvar)
    end if

    ;; Check ordering
    if (isMonotonic(cvar&lon).eq.0) then
        error_msg("f", scriptname, funcname, \
                  "non-monotonic longitude coordinate")
    end if

    ;; Assign attributes
    clon = lon_attrib(cvar&lon)
    delete(cvar&lon)
    cvar&lon = clon
    out = cvar

    leave_msg(scriptname, funcname)
    return(out)

end

;;######################################################################
undef("reformat_coord")
function reformat_coord(var:numeric)
;;
;; Arguments
;;    var: the variable's field.
;;
;; Return value
;;    The input variable with all coordinates reformatted.
;;
;; Description
;;    Reformats all the coordinates by calling the corresponding function.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: written.
;;
local funcname, scriptname, locvar, dnames, ndnames, dd, tmp
begin

    funcname = "reformat_coord"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Create a local copy
    locvar = var

    ;; Get coordinate names
    dnames = getvardims(var)
    ndnames = dimsizes(dnames)

    ;; Loop over dimensions and call the corresponding function
    do dd = 0, ndnames - 1

        ;; Time
        if (dnames(dd).eq."time" .or. \
            any(dnames(dd).eq.find_alt_names("time"))) then
            reformat_time(locvar, dd)
            continue
        end if

        ;; Pressure level
        if (dnames(dd).eq."plev".or. \
            any(dnames(dd).eq.find_alt_names("plev"))) then
            reformat_plev(locvar, dd)
            continue
        end if

        ;; Ocean depth
        if (dnames(dd).eq."lev".or. \
            any(dnames(dd).eq.find_alt_names("lev"))) then
            reformat_lev(locvar, dd)
            continue
        end if

        ;; Basin
        if (dnames(dd).eq."basin") then
            reformat_basin(locvar, dd)
            continue
        end if

        ;; Latitude
        if (dnames(dd).eq."lat".or. \
            any(dnames(dd).eq.find_alt_names("lat"))) then
            if (GRID.eq."rectilinear") then
                reformat_lat(locvar, dd)
                continue
            else
                continue
            end if
        end if

        ;; Longitude
        if (dnames(dd).eq."lon".or. \
            any(dnames(dd).eq.find_alt_names("lon"))) then
            if (GRID.eq."rectilinear") then
                tmp = reformat_lon(locvar, dd)
                delete(locvar)
                locvar = tmp
                delete(tmp)
                continue
            else
                continue
            end if
        end if

        ;; Skip i and j (irregular grids)
        if (any(dnames(dd).eq.(/"i", "j"/))) then
            continue
        end if

        ;; Error if none of the above is called
        error_msg("f", scriptname, funcname, \
                  "a routine for reformatting coordinate " + dnames(dd) + \
                  " is not available")

    end do

    leave_msg(scriptname, funcname)
    return(locvar)

end

;;######################################################################
undef("add_irregulargrid_info")
procedure add_irregulargrid_info(var:numeric)
;;
;; Arguments
;;    var: the variable's field.
;;
;; Description
;;    Add irregular grid info to the output.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20161202-A_laue_ax: preserve attribute "bounds" if present
;;    20130528-A_righ_ma: written.
;;
local funcname, scriptname, vdims, vsize, vrank, files, ff, vars, \
    w, clat, altdname, ii, lat, ldims, lsize, clon, lon, area, tmp, latv, lonv
begin

    funcname = "add_irregulargrid_info"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Variable dimensions (must be consistent with areacello dimensions)
    vdims = getvardims(var)
    vsize = dimsizes(var)
    vrank = dimsizes(vsize)

    ;; Open input file
    files = systemfunc("ls " + INFILE)
    ff = addfile(files(0), "r")
    vars = getfilevarnames(ff)

    ;; Open output file
    w = addfile(OUTFILE, "w")
    setfileoption(w, "DefineMode", True)

    ;; Read, reformat and check latitude coordinate
    if (any("lat".eq.vars)) then
        clat = ff->lat
    else
        altdname = find_alt_names("lat")
        do ii = 0, dimsizes(altdname) - 1
            if (any(altdname(ii).eq.vars)) then
                clat = ff->$altdname(ii)$
                break
            end if
        end do
        delete(altdname)
    end if
    if (.not. isdefined("clat")) then
        error_msg("f", scriptname, funcname, \
                  "cannot find latitude coordinate for the irregular grid")
    end if
    lat = lat_attrib(clat)
; A_laue_ax+
    ; if present, preserve the attribute "bounds"
    if (isatt(clat, "bounds")) then
        lat@bounds = "lat_vertices"
    end if
; A_laue_ax-
    ldims = getvardims(clat)
    lsize = dimsizes(clat)
    if (ldims(1) .ne. vdims(vrank - 1) .or. \
        ldims(0) .ne. vdims(vrank - 2) .or. \
        lsize(1) .ne. vsize(vrank - 1) .or. \
        lsize(0) .ne. vsize(vrank - 2)) then
        error_msg("f", scriptname, funcname, \
                  "inconsistent latitude coordinate for the irregular grid")
    end if
    filevardef(w, "lat", typeof(lat), ldims)
    filevarattdef(w, "lat", lat)
    delete(ldims)
    delete(lsize)

    ;; Read, reformat and check longitude coordinate
    if (any("lon".eq.vars)) then
        clon = ff->lon
    else
        altdname = find_alt_names("lon")
        do ii = 0, dimsizes(altdname) - 1
            if (any(altdname(ii).eq.vars)) then
                clon = ff->$altdname(ii)$
                break
            end if
        end do
        delete(altdname)
    end if
    if (.not. isdefined("clon")) then
        error_msg("f", scriptname, funcname, \
                  "cannot find longitude coordinate for the irregular grid")
    end if
    lon = lon_attrib(clon)
; A_laue_ax+
    ; if present, preserve the attribute "bounds"
    if (isatt(clon, "bounds")) then
        lon@bounds = "lon_vertices"
    end if
; A_laue_ax-
    ldims = getvardims(clon)
    lsize = dimsizes(clon)
    if (ldims(1) .ne. vdims(vrank - 1) .or. \
        ldims(0) .ne. vdims(vrank - 2) .or. \
        lsize(1) .ne. vsize(vrank - 1) .or. \
        lsize(0) .ne. vsize(vrank - 2)) then
        error_msg("f", scriptname, funcname, \
                  "inconsistent longitude coordinate for the irregular grid")
    end if
    filevardef(w, "lon", typeof(lon), ldims)
    filevarattdef(w, "lon", lon)
    delete(ldims)
    delete(lsize)
    delete(vars)

    ;; Read, reformat and check the irregular grid area
    if (.not. fileexists(AREAFILE)) then
        print("")
        print("Cannot find the specified area file")
        print("      " + AREAFILE)
        print("")
        print("**************************************************************")
        print("* NOTE:                                                      *")
        print("**************************************************************")
        if (PROJECT.eq."CMIP5") then
            print("When using the CMIP5 project class in the namelist the")
            print("areafile location is given a default value in the areafile")
            print("routine in the CMIP5 class. See this class in the")
            print("'pyton_code/projects.py'-file for details. To fix this")
            print("error, move your areafile to the above location or change")
            print("to a project class that defines the areafile differently")
            print("e.g., CMIP5_gridfile or one of the CMIP5_????.")
        else if(PROJECT.eq."CMIP5_gridfile") then
            print("The explicit location for the areafile in your namelist")
            print("seems to be incorrect. Move your areafile to the above")
            print("location or update your namelist.")
        else if(isStrSubset(PROJECT, "CMIP5_")) then
            print("The path to the areafile as defined in the " + PROJECT)
            print("class seems incorrect. Check that the class defintion in")
            print("interface_scripts/projects.py for the areafile is correct.")
        end if
        end if
        end if
        print("See the pdf:s in the 'doc/'-folder for more details on the")
        print("handling of the various CMIP-classes")
        print("**************************************************************")
        print("")
        error_msg("f", scriptname, funcname, "")
    else
        ff = addfile(AREAFILE, "r")
        vars = getfilevarnames(ff)

        ;; Area
        if (any("areacello".eq.vars)) then
            area = ff->areacello
        else
            altdname = find_alt_names("areacello")
            do ii = 0, dimsizes(altdname) - 1
                if (any(altdname(ii).eq.vars)) then
                    area = ff->$altdname(ii)$
                    break
                end if
            end do
            delete(altdname)
        end if
        check_units(area, "areacello")
        if (.not. isdefined("area")) then
            error_msg("f", scriptname, funcname, \
                      "cannot find areacello for the irregular grid")
        end if
        tmp = var_attrib(area, "areacello")
        delete(area)
        area = tmp
        delete(tmp)
        ldims = getvardims(area)
        lsize = dimsizes(area)

        ;; Check consistency of dimension sizes
        if (lsize(1) .ne. vsize(vrank - 1) .or. \
            lsize(0) .ne. vsize(vrank - 2)) then
            error_msg("f", scriptname, funcname, \
                      "inconsistent areacello for the irregular grid")
        end if

        ;; Check consistency of dimension names
        if (ldims(1) .ne. vdims(vrank - 1) .or. \
            ldims(0) .ne. vdims(vrank - 2)) then
            area!1 = vdims(vrank - 1)
            area!0 = vdims(vrank - 2)
            ldims = getvardims(area)
        end if

        log_info("  Adding areacello to OUTFILE")
        filevardef(w, "areacello", typeof(area), ldims)
        filevarattdef(w, "areacello", area)
        delete(ldims)
        delete(lsize)

        ;; Latitude vertices
        if (any("lat_vertices".eq.vars)) then
            latv = ff->lat_vertices
        else
            altdname = find_alt_names("lat_vertices")
            do ii = 0, dimsizes(altdname) - 1
                if (any(altdname(ii).eq.vars)) then
                    latv = ff->$altdname(ii)$
                    break
                end if
            end do
            delete(altdname)
        end if

        if (isdefined("latv")) then
            ldims = getvardims(latv)
            lsize = dimsizes(latv)

            ;; Check consistency of dimension sizes
            if (lsize(1) .ne. vsize(vrank - 1) .or. \
                lsize(0) .ne. vsize(vrank - 2)) then
                error_msg("f", scriptname, funcname, \
                          "inconsistent lat_vertices for the irregular grid")
            end if

            ;; Check consistency of dimension names
            if (ldims(1) .ne. vdims(vrank - 1) .or. \
                ldims(0) .ne. vdims(vrank - 2)) then
                latv!1 = vdims(vrank - 1)
                latv!0 = vdims(vrank - 2)
                ldims = getvardims(latv)
            end if
            latv!2 = "vertices"

            log_info("  Adding lat_vertices to OUTFILE")
            filedimdef(w, (/"vertices"/), (/4/), (/False/))
            filevardef(w, "lat_vertices", typeof(latv), getvardims(latv))
            filevarattdef(w, "lat_vertices", latv)

            delete(ldims)
            delete(lsize)
        end if

        ;; Longitude vertices
        if (any("lon_vertices".eq.vars)) then
            lonv = ff->lon_vertices
        else
            altdname = find_alt_names("lon_vertices")
            do ii = 0, dimsizes(altdname) - 1
                if (any(altdname(ii).eq.vars)) then
                    lonv = ff->$altdname(ii)$
                    break
                end if
            end do
            delete(altdname)
        end if
        if (isdefined("lonv")) then
            ldims = getvardims(lonv)
            lsize = dimsizes(lonv)

            ;; Check consistency of dimension sizes
            if (lsize(1) .ne. vsize(vrank - 1) .or. \
                lsize(0) .ne. vsize(vrank - 2)) then
                error_msg("f", scriptname, funcname, \
                          "inconsistent lon_vertices for the irregular grid")
            end if

            ;; Check consistency of dimension names
            if (ldims(1) .ne. vdims(vrank - 1) .or. \
                ldims(0) .ne. vdims(vrank - 2)) then
                lonv!1 = vdims(vrank - 1)
                lonv!0 = vdims(vrank - 2)
                ldims = getvardims(lonv)
            end if
            lonv!2 = "vertices"

            log_info("  Adding lon_vertices to OUTFILE")
            filevardef(w, "lon_vertices", typeof(lonv), getvardims(lonv))
            filevarattdef(w, "lon_vertices", lonv)

            delete(ldims)
            delete(lsize)
        end if

    end if

    ;; Write output
    setfileoption(w, "DefineMode", False)
    w->lat = (/lat/)
    w->lon = (/lon/)
    w->areacello = (/area/)
    if (isdefined("latv")) then
        w->lat_vertices = (/latv/)
    end if
    if (isdefined("lonv")) then
        w->lon_vertices = (/lonv/)
    end if
    delete(w)

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("add_lmask_info")
procedure add_lmask_info(var:numeric)
;;
;; Arguments
;;
;; Description
;;    Add land mask to the output.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151029-A_righ_ma: written.
;;
local funcname, scriptname, lexist, ff, vname, altdname, maskl, \
    tmp, w, mdims, msize
begin

    funcname = "add_lmask_info"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Variable dimensions (must be consistent with sftlf dimensions)
    vdims = getvardims(var)
    vsize = dimsizes(var)
    vrank = dimsizes(vsize)

    lexist = False
    if (fileexists(LMASKFILE)) then
        ff = addfile(LMASKFILE, "r")
        if (isfilevar(ff, "sftlf")) then
            lexist = True
            vname = "sftlf"
        else
            altdname = find_alt_names("sftlf")
            do ii = 0, dimsizes(altdname) - 1
                if (isfilevar(ff, altdname(ii))) then
                    lexist = True
                    vname = altdname(ii)
                    break
                end if
            end do
            delete(altdname)
        end if
    end if

    if (lexist) then
        maskl = ff->$vname$

        ;; Reformat variable
        apply_fixes(maskl, "sftlf", "")
        check_units(maskl, "sftlf")
        tmp = var_attrib(maskl, "sftlf")
        delete(maskl)
        maskl = tmp
        delete(tmp)
        tmp = reformat_coord(maskl)
        delete(maskl)
        maskl = tmp
        delete(tmp)

        ;; Check latitude/longitude consistency
        lwrite = True
        mdims = getvardims(maskl)
        msize = dimsizes(maskl)
        mrank = dimsizes(msize)
        if (mdims(mrank - 2).ne.vdims(vrank - 2) .or. \
            msize(mrank - 2).ne.vsize(vrank - 2)) then
            lwrite = False
        end if
        if (mdims(mrank - 1).ne.vdims(vrank - 1) .or. \
            msize(mrank - 1).ne.vsize(vrank - 1)) then
            lwrite = False
        end if

        ;; Write output
        if (lwrite) then
            w = addfile(OUTFILE, "w")
            setfileoption(w, "DefineMode", True)
            filevardef(w, "sftlf", typeof(maskl), mdims)
            filevarattdef(w, "sftlf", maskl)
            setfileoption(w, "DefineMode", False)
            w->sftlf = (/maskl/)
            log_info("  Adding sftlf to OUTFILE")
        else
            log_debug("Dimensionality of sftlf not consistent with " + \
                        "dimensionality of " + VAR + ", skipping")
        end if
    end if

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("add_omask_info")
procedure add_omask_info(var:numeric)
;;
;; Arguments
;;
;; Description
;;    Add ocean mask to the output.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151029-A_righ_ma: written.
;;
local funcname, scriptname, lexist, ff, vname, altdname, maskl, \
    tmp, w, mdims, msize
begin

    funcname = "add_omask_info"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Variable dimensions (must be consistent with sftof dimensions)
    vdims = getvardims(var)
    vsize = dimsizes(var)
    vrank = dimsizes(vsize)

    lexist = False
    if (fileexists(OMASKFILE)) then
        ff = addfile(OMASKFILE, "r")
        if (isfilevar(ff, "sftof")) then
            lexist = True
            vname = "sftof"
        else
            altdname = find_alt_names("sftof")
            do ii = 0, dimsizes(altdname) - 1
                if (isfilevar(ff, altdname(ii))) then
                    lexist = True
                    vname = altdname(ii)
                    break
                end if
            end do
            delete(altdname)
        end if
    end if

    if (lexist) then
        masko = ff->$vname$

        ;; Reformat variable
        apply_fixes(masko, "sftof", "")
        check_units(masko, "sftof")
        tmp = var_attrib(masko, "sftof")
        delete(masko)
        masko = tmp
        delete(tmp)
        tmp = reformat_coord(masko)
        delete(masko)
        masko = tmp
        delete(tmp)

        ;; Check latitude/longitude consistency
        lwrite = True
        mdims = getvardims(masko)
        msize = dimsizes(masko)
        mrank = dimsizes(msize)
        if (mdims(mrank - 2).ne.vdims(vrank - 2) .or. \
            msize(mrank - 2).ne.vsize(vrank - 2)) then
            lwrite = False
        end if
        if (mdims(mrank - 1).ne.vdims(vrank - 1) .or. \
            msize(mrank - 1).ne.vsize(vrank - 1)) then
            lwrite = False
        end if

        ;; Write output
        if (lwrite) then
            w = addfile(OUTFILE, "w")
            setfileoption(w, "DefineMode", True)
            filevardef(w, "sftof", typeof(masko), mdims)
            filevarattdef(w, "sftof", masko)
            setfileoption(w, "DefineMode", False)
            w->sftof = (/masko/)
            log_info("  Adding sftof to OUTFILE")
        else
            log_debug("Dimensionality of sftof not consistent with " + \
                      "dimensionality of " + VAR + ", skipping")
        end if

    end if

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("add_porosity_info")
procedure add_porosity_info(var:numeric)
;;
;; Arguments
;;
;; Description
;;    Add porosity to the output.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20160823-A_righ_ma: written.
;;
local funcname, scriptname, vdims, vsize, vrank, lexist, ff, \
    vname, altdname, ii, maskl, tmp, lwrite, mdims, msize, mrank, w
begin

    funcname = "add_porosity_info"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Variable dimensions (must be consistent with mrsofc dimensions)
    vdims = getvardims(var)
    vsize = dimsizes(var)
    vrank = dimsizes(vsize)

    lexist = False
    if (fileexists(POROFILE)) then
        ff = addfile(POROFILE, "r")
        if (isfilevar(ff, "mrsofc")) then
            lexist = True
            vname = "mrsofc"
        else
            altdname = find_alt_names("mrsofc")
            do ii = 0, dimsizes(altdname) - 1
                if (isfilevar(ff, altdname(ii))) then
                    lexist = True
                    vname = altdname(ii)
                    break
                end if
            end do
            delete(altdname)
        end if
    end if

    if (lexist) then
        maskl = ff->$vname$

        ;; Reformat variable
        apply_fixes(maskl, "mrsofc", "")
        check_units(maskl, "mrsofc")
        tmp = var_attrib(maskl, "mrsofc")
        delete(maskl)
        maskl = tmp
        delete(tmp)
        tmp = reformat_coord(maskl)
        delete(maskl)
        maskl = tmp
        delete(tmp)

        ;; Check latitude/longitude consistency
        lwrite = True
        mdims = getvardims(maskl)
        msize = dimsizes(maskl)
        mrank = dimsizes(msize)
        if (mdims(mrank - 2).ne.vdims(vrank - 2) .or. \
            msize(mrank - 2).ne.vsize(vrank - 2)) then
            lwrite = False
        end if
        if (mdims(mrank - 1).ne.vdims(vrank - 1) .or. \
            msize(mrank - 1).ne.vsize(vrank - 1)) then
            lwrite = False
        end if

        ;; Write output
        if (lwrite) then
            w = addfile(OUTFILE, "w")
            setfileoption(w, "DefineMode", True)
            filevardef(w, "mrsofc", typeof(maskl), mdims)
            filevarattdef(w, "mrsofc", maskl)
            setfileoption(w, "DefineMode", False)
            w->mrsofc = (/maskl/)
            log_info("  Adding mrsofc to OUTFILE")
        else
            log_info("Dimensionality of mrsofc not consistent with " + \
                        "dimensionality of " + VAR + ", skipping")
        end if
    end if

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("add_depth_info")
procedure add_depth_info(var:numeric)
;;
;; Arguments
;;
;; Description
;;    Add depth to the output.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20160824-A_righ_ma: written.
;;
local funcname, scriptname, files, ff, depths, ii, lexist, vname, \
    altdname, w
begin

    funcname = "add_depth_info"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Any of the input files is ok, since depth/depth_bnds are constants
    files = systemfunc("ls " + INFILE)
    ff = addfile(files(0), "r")

    ;; Read and add depth variables
    depths = (/"depth", "depth_bnds"/)
    do ii = 0, dimsizes(depths) - 1
        lexist = False
        if (isfilevar(ff, depths(ii))) then
            lexist = True
            vname = depths(ii)
        else
            altdname = find_alt_names(depths(ii))
            do ii = 0, dimsizes(altdname) - 1
                if (isfilevar(ff, altdname(ii))) then
                    lexist = True
                    vname = altdname(ii)
                    break
                end if
            end do
            delete(altdname)
        end if

        if (lexist) then
            depth = ff->$vname$
            w = addfile(OUTFILE, "w")
            w->$depths(ii)$ = depth
            log_info("  Adding " + depths(ii) + " to OUTFILE")
            delete(depth)
        else
            error_msg("f", scriptname, funcname, \
                      "variable " + name + " not found in input file")
        end if
    end do

    leave_msg(scriptname, funcname)

end

;;######################################################################
undef("add_regionbasin_info")
procedure add_regionbasin_info()
;;
;; Arguments
;;
;; Description
;;    Add region/basin info to the output.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: written.
;;
local funcname, scriptname, files, ff, w, region, dims, size
begin

    funcname = "add_regionbasin_info"
    scriptname = "reformat_scripts/default/reformat_default_func.ncl"
    enter_msg(scriptname, funcname)

    ;; Open input file
    files = systemfunc("ls " + INFILE)
    ff = addfile(files(0), "r")
    vars = getfilevarnames(ff)

    ;; Open output file
    w = addfile(OUTFILE, "w")
    setfileoption(w, "DefineMode", True)

    ;; Read region coordinate
    if (isfilevar(ff, "region"))
        region = ff->region
    else
        error_msg("f", scriptname, funcname, \
                  "cannot find region coordinate for variable " + VAR)
    end if

    ;; Read attributes
    region&strlen = ispan(1, dimsizes(region(0, :)), 1)
    filedimdef(w, (/"strlen"/), (/21/), (/False/))
    dims = getvardims(region)
    size = dimsizes(region)

    filevardef(w, "region", typeof(region), dims)
    filevarattdef(w, "region", region)

    ;; Write output
    setfileoption(w, "DefineMode", False)
    w->region = (/region/)

    leave_msg(scriptname, funcname)

end
