;;#############################################################################
;; CONSTANTS AND GENERAL PURPOSE FUNCTIONS FOR THE REFORMAT ROUTINES
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function check_time_range
;;    function time_attrib
;;    function plev_attrib
;;    function lev_attrib
;;    function lat_attrib
;;    function lon_attrib
;;    function read_cmor
;;    procedure check_units
;;    function var_attrib
;;    function define_globatt
;;    procedure write_output
;;    function set_size_array
;;
;;#############################################################################

load "./interface_data/ncl.interface"

scriptname = "reformat_scripts/constants.ncl"

;;#############################################################################

;; Pressure levels for interpolation [hPa]
PLEVS = (/1000, 990, 980, 970, 960, 950, 900, 850, 700, 500, 400, 300, 250, \
         200, 170, 150, 130, 115, 100, 90, 80, 70, 50, 30, 20, 15, 10, \
         7, 5, 3, 2, 1.5, 1, 0.5, 0.3, 0.2, 0.1d/)
PLEVS@units = "hPa"

;; Time units
TUNITS = "days since 1950-01-01 00:00:00"

;; FillValue
FILL = 1.e+20

;; Data file for recognized units
F_UNITS = "./reformat_scripts/recognized_units.dat"

;; Data file for recognized variable names
F_VARIABLES = "./reformat_scripts/recognized_vars.dat"

;;#############################################################################
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
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    if(.not.isdefined("YEAR1")) then
        YEAR1 = time@YEAR1
    end if
    if(.not.isdefined("YEAR2")) then
        YEAR2 = time@YEAR2
    end if
    if(.not.isdefined("FIELD")) then
        FIELD = time@FIELD
    end if

    date = cd_calendar(time, 0)
    ntime = dimsizes(time)

    ;; Check monotonicity
    if (isMonotonic(time).ne.1) then
        error_msg("fatal", "constants.ncl", funcname, \
                  "non-monotonically-increasing time coordinate in input data")
    end if

    ;; Boundaries
    if (date(0, 0).gt.YEAR1) then
        error_msg("fatal", "constants.ncl", funcname, \
                  "selected period outside range (" + YEAR1 + \
                  "<" + date(0, 0) + ")")
    end if
    if (date(ntime - 1, 0).lt.YEAR2) then
        error_msg("fatal", "constants.ncl", funcname, \
                  "selected period outside range (" + YEAR2 + \
                  ">" + date(ntime - 1, 0) + ")")
    end if

    ;; Monthly mean case
    if (isStrSubset(FIELD, "M")) then

        ;; Check number of timesteps
        exp_ntime = 12 * (YEAR2 - YEAR1 + 1)
        if (ntime.ne.exp_ntime) then
            error_msg("fatal", "constants.ncl", funcname, \
                      "incorrect number of timesteps in input data (" + \
                      ntime + " != " + exp_ntime + "). Possibly due to " + \
                      "incomplete years in input.")
        end if

        ;; Set time coordinate to 15th of each month at midnight
        do tt = 0, ntime - 1
            yy = YEAR1 + tt / 12
            mm = mod(tt, 12) + 1
            if (date(tt, 0).ne.yy .or. date(tt, 1).ne.mm) then
                error_msg("fatal", "constants.ncl", funcname, \
                          "inconsistent time coordinate in input data" + \
                          yy + sprinti("%0.2i", mm))
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
                error_msg("w", "constants.ncl", funcname, \
                          "changing calendar from proleptic_gregorian " + \
                          "to gregorian prior to 1582")
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
            error_msg("fatal", "constants.ncl", funcname, \
                      "incorrect number of timesteps in input data (" + \
                      ntime + " != " + exp_ntime + "). Possibly due to " + \
                      "incomplete years in input.")
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

    leave_msg(scriptname, funcname, 4)
    return(out)

end

;;#############################################################################
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
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

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

    leave_msg(scriptname, funcname, 4)
    return(out_coord)

end

;;#############################################################################
undef("plev_attrib")
function plev_attrib(coord:numeric)
;;
;; Arguments
;;    coord: input pressure level coordinate.
;;
;; Return value
;;    An array of the same dimensionality of coord.
;;
;; Description
;;    Assigns attributes to the pressure level coordinate.
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

    funcname = "plev_attrib"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    out_coord = todouble(coord)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@positive = "down"
    out_coord@axis     = "Z"
    if (coord@units .eq. "m") then
        out_coord@long_name     = "ocean depth coordinate"
        out_coord@units         = "m"
        out_coord@standard_name = "depth"
    else
        out_coord@long_name     = "pressure"
        out_coord@units         = "Pa"
        out_coord@standard_name = "air_pressure"
    end if
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    leave_msg(scriptname, funcname, 4)
    return(out_coord)

end

;;#############################################################################
undef("lev_attrib")
function lev_attrib(coord:numeric)
;;
;; Arguments
;;    coord: input ocean level coordinate.
;;
;; Return value
;;    An array of the same dimensionality of coord.
;;
;; Description
;;    Assigns attributes to the ocean level coordinate.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150415-A_righ_ma: written.
;;
local funcname, verbosity, out_coord
begin

    funcname = "lev_attrib"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    out_coord = todouble(coord)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@long_name = "ocean depth coordinate"
    out_coord@positive = "down"
    out_coord@axis = "Z"
    out_coord@units = "m "
    out_coord@standard_name = "depth"
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    leave_msg(scriptname, funcname, 4)
    return(out_coord)

end

;;#############################################################################
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
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    out_coord = todouble(coord)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@long_name = "latitude"
    out_coord@axis = "Y"
    out_coord@units = "degrees_north"
    out_coord@standard_name = "latitude"
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    leave_msg(scriptname, funcname, 4)
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
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    out_coord = todouble(coord)  ; this also removes attributes
    copy_VarCoords(coord, out_coord)
    out_coord@long_name = "longitude"
    out_coord@axis = "X"
    out_coord@units = "degrees_east"
    out_coord@standard_name = "longitude"
    if (isatt(out_coord, "_FillValue")) then
        delete(out_coord@_FillValue)
    end if

    leave_msg(scriptname, funcname, 4)
    return(out_coord)

end

;;#############################################################################
undef("read_cmor")
function read_cmor(name:string)
;;
;; Arguments
;;    name: standard variable name.
;;
;; Return value
;;    A logical variable with the CMOR table attached as attributes.
;;
;; Description
;;    Read variable attributes from the CMOR tables (cmor/<name>.cmor).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, fn, data, idxu, idxd, attn, attv
begin

    funcname = "read_cmor"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read attributes from cmor table
    fn = "./reformat_scripts/cmor/CMOR_" + name + ".dat"
    if (.not.isfilepresent_esmval(fn)) then
        error_msg("fatal", "constants.ncl", funcname, \
                  "CMOR table for variable " + name + " not found")
    end if

    data = readAsciiTable(fn, 1, "string", 0)
    idxu = ind(str_get_field(data(:, 0), 1, ":").eq. \
               "! Variable attributes") + 2
    idxd = ind(str_get_field(data(:, 0), 1, ":").eq. \
               "! Additional variable information") - 2
    attn = str_squeeze(str_get_field(data(idxu:idxd, 0), 1, ":"))
    attv = str_squeeze(str_get_field(data(idxu:idxd, 0), 2, ":"))

    out = True
    do ii = 0, dimsizes(attn) - 1
        out@$attn(ii)$ = attv(ii)
    end do

    leave_msg(scriptname, funcname, 4)
    return(out)

end

;;#############################################################################
undef("check_units")
procedure check_units(var:numeric, name:string)
;;
;; Arguments
;;    var: the variable's field.
;;    name: the variable's standard name.
;;
;; Description
;;    Check and convert units.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, atts, data, col1, col2, std, alt, id, tmp1, tmp2, \
    ncols, alt_units, factors
begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    funcname = "check_units"
    info_output("In check_units for " + name +" units are " + var@units,\
                verbosity, 4)
    enter_msg(scriptname, funcname, 4)

    ;; CMOR units
    atts = read_cmor(name)
    if (isatt(atts, "units")) then
        cmor_units = atts@units
    else
        error_msg("fatal", "constants.ncl", funcname, \
                  "cannot find units attribute in the CMOR table for " + \
                  "variable " + name)
    end if
    delete(atts)

    ;; If units not defined assign CMOR units
    if (.not.isatt(var, "units")) then
        info_output("No units attribute found for this variable, " + \
                    "assigning CMOR units", verbosity, 1)
        var@units = cmor_units
    end if

    ;; Check units
    if (var@units.ne.cmor_units) then

        ;; Read alternative units and conversion factors
        data = readAsciiTable(F_UNITS, 1, "string", 0)
        col1 = str_squeeze(str_get_field(data(:, 0), 1, "="))
        col2 = str_squeeze(str_get_field(data(:, 0), 2, "="))
        std = col2(ind(col1.eq."std_unit"))
        alt = col2(ind(col1.eq."alt_unit"))
        id = ind(std.eq.cmor_units)
        if (ismissing(id)) then
            error_msg("fatal", "constants.ncl", funcname, \
                      "units " + var@units + " not in " + F_UNITS)
        else
            tmp1 = alt(id)
        end if
        ncols = str_fields_count(tmp1, ",")
        if (ncols.ne.0) then
            tmp2 = new(ncols, string)
            do ii = 0, ncols - 1
                tmp2(ii) = str_squeeze(str_get_field(tmp1, ii + 1, ","))
            end do
            alt_units = tmp2(::2)
            factors = tofloat(tmp2(1::2))
        else
            error_msg("fatal", "constants.ncl", funcname, \
                      "cannot convert units " + var@units + ", add an " + \
                      "appropriate unit conversion factor to " + F_UNITS)
        end if

        ;; Check consistency and convert
        if (any(var@units.eq.alt_units)) then
            id = ind(var@units.eq.alt_units)
            ;; Rename units
            var@units = cmor_units
            info_output("Renaming units to " + cmor_units, verbosity, 4)
            ;; Convert units
            if (factors(id).ne.1) then
                var = var * factors(id)
                convert_msg = "Converting units ["\
                            + alt_units(id)\
                            + "] to ["\
                            + cmor_units + "]"
                info_output(convert_msg, verbosity, 4)
                info_output("Conversion factor = " + factors(id), verbosity, 4)
            end if
        else
            error_msg("fatal", "constants.ncl", funcname, \
                      "cannot convert units " + var@units + ", add an " + \
                      "appropriate unit conversion factor to " + F_UNITS)
        end if

    end if

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("var_attrib")
function var_attrib(var:numeric, name:string)
;;
;; Arguments
;;    var: the variable's field.
;;    name: the variable's standard name.
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
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

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
    tmp = read_cmor(name)
    att = getvaratts(tmp)
    do ii = 0, dimsizes(att) - 1
        out@$att(dimsizes(att) - 1 - ii)$ = tmp@$att(dimsizes(att) - 1 - ii)$
    end do

    leave_msg(scriptname, funcname, 4)
    return(out)

end

;#######################################################################
undef("define_globatt")
function define_globatt()
;;
;; Arguments
;;
;; Return value
;;    A logical variable with  attached global attributes.
;;
;; Description
;;    Sets file global attributes.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, field, user
begin

    funcname = "define_globatt"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Define field name
    if (isStrSubset(FIELD, "M")) then
        field = FIELD + ": Monthly-mean "
    end if
    if (isStrSubset(FIELD, "D")) then
        field = FIELD + ": Daily-mean "
    end if
    if (isStrSubset(FIELD, "I")) then
        field = FIELD + ": Instantaneous "
    end if
    if (isStrSubset(FIELD, "T3")) then
      if (isStrSubset(FIELD, "h")) then
        field = field + "3D atmospheric data (time,lev,lat,lon)"
      else
        field = field + "3D atmospheric data (time,press,lat,lon)"
      end if
    end if
    if (isStrSubset(FIELD, "z")) then
      field = field + "2D atmospheric data (time,press,lat)"
    end if
    if (isStrSubset(FIELD, "s")) then
      field = field + "2D atmospheric data (time,lat,lon)"
    end if

    ;; Find user name
    user = getenv("USER")

    ;; Find host name
    global = True
    global@title = "Reformatted model output data for the ESMVal " + \
        "evaluation tool"
    global@source = "For the CMOR standard see " + \
        "http://www2-pcmdi.llnl.gov/cmor/tables/"
    global@project = PROJECT
    global@model = MODEL
    if (isdefined("EXPER")) then
        global@experiment = EXPER
    end if
    global@variable = VAR
    global@field = field
    global@period = sprinti("%0.4i",YEAR1) + "01-" + \
        sprinti("%0.4i",YEAR2) + "12"
    global@user = user
    global@host = systemfunc("echo $HOST")
    global@history = "Creation date " + systemfunc("date")

    leave_msg(scriptname, funcname, 4)
    return(global)

end

;;#############################################################################
undef("write_output")
; a_laue_ax+
;procedure write_output(var:numeric, name:string)
procedure write_output(var:numeric, name:string, filelist:string, \
                       fixfile:string, version:string)
; al_laue_ax-
;;
;; Arguments
;;    var: the variable's field.
;;    name: the variable's standard name.
;;    filelist: list of filenames processed to produce the file to be written
;;    fixfile: filename of a FIXFILE applied to process the data to be written
;;             (if applicable)
;;    version: version of the ESMValTool
;;
;; Description
;;    Writes the selected variable in the output file.
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
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, w, gAtt, dim_names, ndims, dim_sizes, dim_types, \
    dim_unlim, ii, lat, lon
begin

    funcname = "write_output"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Open file
    if (isfilepresent_esmval(OUTFILE)) then
        system("rm -f " + OUTFILE)
    end if
    w = addfile(OUTFILE, "c")
    setfileoption(w, "DefineMode", True)

    ;; Attach global attributes
    gAtt = define_globatt()
; a_laue_ax+
    ;; adding list of files processed, name of FIXFILE applied and ESMValTool
    ;; version to the global attributes

    do i = 0, dimsizes(filelist) - 1
        if (ismissing(filelist(i))) then
            break
        end if
    end do

    ;; we switch the order of the creation of the attributes "infile" so that
    ;; "infile_0000" will be on top of the list

    if (i.gt.0) then
        do j = i - 1, 0, 1
            gAtt@$sprinti("infile_%04d", j)$ = filelist(j)
        end do
    end if

    if (strlen(fixfile).gt.0) then
        gAtt@fixfile = fixfile
    else
        gAtt@fixfile = "none"
    end if
    gAtt@version = version
; a_laue_ax-
    fileattdef(w, gAtt)

    ;; Get coordinates
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

    ;; Time coordinate must be unlimited
    if (any(dim_names.eq."time")) then
        dim_sizes(ind(dim_names.eq."time")) = -1
        dim_unlim(ind(dim_names.eq."time")) = True
    end if

    ;; Define dimensions
    filedimdef(w, dim_names, dim_sizes, dim_unlim)
    do ii = 0, ndims - 1
        filevardef(w, dim_names(ii), dim_types(ii), dim_names(ii))
    end do
    filevardef(w, name, "float", dim_names)
    do ii = 0, ndims - 1
        filevarattdef(w, dim_names(ii), var&$dim_names(ii)$)
    end do
    filevarattdef(w, name, var)

    ;; Write
    setfileoption(w, "DefineMode", False)
    do ii = 0, ndims - 1
        w->$dim_names(ii)$ = (/var&$dim_names(ii)$/)
    end do
    w->$name$ = (/var/)

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("set_size_array")
function set_size_array()
;;
;; Arguments
;;
;; Return value
;;    An array of type double.
;;
;; Description
;;    Sets a logarithmic array of sizes to be used for particle size
;;    distribution calculations.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130528-A_righ_ma: written.
;;
local funcname, verbosity, minsize, maxsize, nbins, bin
begin

    funcname = "set_size_array"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Size range (0.5 nm - 10 um)
    minsize = 0.5e-9
    maxsize = 10.e-6
    nbins = 100

    ;; Generate array
    out = new(nbins, double)
    bin = 10.^(log10(maxsize / minsize) / (nbins - 1))
    out(0) = minsize
    do ii = 1, nbins - 1
        out(ii) = out(ii - 1) * bin
    end do

    leave_msg(scriptname, funcname, 4)
    return(out)

end
