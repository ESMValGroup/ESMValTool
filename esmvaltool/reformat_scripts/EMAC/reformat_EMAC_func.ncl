;;#############################################################################
;; FUNCTIONS FOR THE EMAC REFORMAT ROUTINES
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function fix_FillValue
;;    procedure check_inpath
;;    procedure find_messy_version
;;    function get_aerosol_modes
;;    function get_aerosol_spec
;;    function find_name
;;    function get_vfatmp
;;    procedure put_vfatmp
;;    function find_and_read_var
;;    function time_mean
;;    function extract_var
;;    procedure format_coord
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"

load "./reformat_scripts/constants.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"

;; MESSy version (initialize)
VERSION = 0
SUBVERSION = 0

;;#############################################################################
undef("fix_FillValue")
function fix_FillValue(var:numeric)
;;
;; Arguments
;;    var: input variable.
;;
;; Return value
;;    An array of the same dimensionality of var.
;;
;; Description
;;    Standard EMAC missing value is -1.e34, but interpolations in post-
;;    processing scripts might consider it as real. Set everything below a
;;    certain threshold to _FillValue.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140217-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, result, upper, lower, tol
begin

    funcname = "fix_FillValue"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    if (var@_FillValue.eq.-1.e34) then
        ;; specify threshold
        tol = 0.000001
        upper = var@_FillValue * tol
        error_msg("w", "reformat_EMAC_func.ncl", funcname, \
                  "original _FillValue = " + var@_FillValue + \
                  ". Considering val < " + upper + " as missing.")
        result = where(var.lt.upper, var@_FillValue, var)
        copy_VarMeta(var, result)
    end if
    leave_msg(scriptname, funcname, 4)
    return(result)

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
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, subdir, ii, pp, len_inpath, last_char
begin

    funcname = "check_inpath"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Simulation path
    if (.not.isfilepresent_esmval(INPATH)) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, "model data " + \
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
undef("find_messy_version")
procedure find_messy_version()
;;
;; Arguments
;;
;; Description
;;    Finds current MESSy version and subversion.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130505_A_gott_kl: modified do deal with leading letters.
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, found, ii, fn, f, tmp, jj, temp1
begin

    funcname = "find_messy_version"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    found = False
    fn = systemfunc("find " + INPATH + " -type f -name '" + EXPER + \
                    "*_" + YEAR1 + "01_*.nc' 2>/dev/null")

    if (all(ismissing(fn))) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, \
                  "no files found for experiment " + EXPER)
    end if

    ii = 0
    do while (.not.found)
        f = addfile(fn(ii), "r")
        if (isatt(f, "messy_version")) then
            tmp = str_squeeze(f@messy_version)
        end if
        if (isatt(f, "MESSy")) then
            tmp = str_sub_str(str_squeeze(\
                str_get_field(f@MESSy, 3, " ")), ",", "")
        end if
        if (isdefined("tmp")) then
            found = True
            SUBVERSION = toint(str_get_field(tmp, 2, "."))

            ;; Remove leading letters (e.g., the "d" in development versions)
            temp1 = tostring(ispan(0, 9, 1))  ; valid strings
            jj = 0
            do while (.not. any(temp1 .eq. str_get_cols( \
                                str_get_field(tmp, 1, "."), jj, jj)))
                jj = jj + 1
            end do
            VERSION = toint(str_get_cols(str_get_field(tmp, 1, "."), jj, -1))
            delete(temp1)
            delete(jj)
            delete(tmp)
        end if
        ii = ii + 1
    end do

    if (VERSION.eq.0 .and. SUBVERSION.eq.0) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, \
                  "cannot determine MESSy version of experiment " + EXPER)
    end if

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("get_aerosol_modes")
function get_aerosol_modes()
;;
;; Arguments
;;
;; Return value
;;    An array of string with aerosol modes, with corresponding sigmas as
;;    attribute.
;;
;; Description
;;    Returns all possible aerosol mode names.
;;
;; Caveats
;;    The current version referes to MADE, MADE-IN and MADE3.
;;
;; References
;;
;; Modification history
;;    20140915-A_righ_ma: written.
;;
local funcname, scriptname, verbosity
begin

    funcname = "get_aerosol_modes"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Possible modes names
    modes = (/"_i", "_km", "_ks", "_ki", "_is", "_sooti", \
             "_j", "_am", "_as", "_ai", "_js", "_sootj", \
             "_c", "_cm", "_cs", "_ci"/)

    ;; Corresponding sigmas
    modes@sigma = (/1.7, 1.7, 1.7, 1.7, 1.7, 1.7, \
                   2.0, 2.0, 2.0, 2.0, 2.0, 2.0, \
                   2.2, 2.2, 2.2, 2.2/)

    if (dimsizes(modes).ne.dimsizes(modes@sigma)) then
        error_msg("f", "reformat_EMAC_func", funcname, "Number of modes " + \
                  "and corresponding sigmas must be consistent")
    end if

    leave_msg(scriptname, funcname, 4)
    return(modes)

end

;;#############################################################################
undef("get_aerosol_spec")
function get_aerosol_spec()
;;
;; Arguments
;;
;; Return value
;;    An array of strings with aerosol species.
;;
;; Description
;;    Returns all possible aerosol species.
;;
;; Caveats
;;    The current version referes to MADE, MADE-IN and MADE3.
;;
;; References
;;
;; Modification history
;;    20140929-A_righ_ma: written.
;;
local funcname, scriptname, verbosity
begin

    funcname = "get_aerosol_spec"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Possible aerosol species
    aero = (/"SO4", "NH4", "NO3", "Na", "Cl", "POM", "POM2", "POMphil", \
            "POMphob", "BC", "BC2", "DU", "H2O", "SS"/)

    leave_msg(scriptname, funcname, 4)
    return(aero)

end

;;#############################################################################
undef("find_name")
function find_name(var:string)
;;
;; Arguments
;;    var: input variable.
;;
;; Return value
;;    The EMAC name of the selected variable.
;;
;; Description
;;    Finds the corresponding EMAC name of the input variable var.
;;
;; Caveats
;;    Needs an external input table "names_EMAC.dat"
;;
;; References
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, fname, data, std, emac, type, vert, id
begin

    funcname = "find_name"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    fname = "./reformat_scripts/EMAC/names_EMAC.dat"

    data = readAsciiTable(fname, 1, "string", \
                          dimsizes(readAsciiHead(fname, "STANDARD")) + 1)

    std  = str_squeeze(str_get_field(data(:, 0), 1, "|"))
    emac = str_squeeze(str_get_field(data(:, 0), 2, "|"))
    type = str_squeeze(str_get_field(data(:, 0), 3, "|"))
    vert = str_squeeze(str_get_field(data(:, 0), 4, "|"))

    id = ind(std.eq.var)

    if (ismissing(id)) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, \
                  "cannot find EMAC name for " + var)
    else
        out = emac(id)
        out@type = type(id)
        out@vert = vert(id)
    end if

    leave_msg(scriptname, funcname, 4)
    return(out)

end

;;#############################################################################
undef("put_vfatmp")
procedure put_vfatmp(var:string,
                     date:string,
                     fn:string)
;;
;; Arguments
;;    var: input variable.
;;    date: selected date (YYYYMM format).
;;    fn: name of the file containing the variable.
;;
;; Description
;;    Writes the variable-file association to the VFATMP temporary file.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150324-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, dn, bn, str
begin

    funcname = "put_vfatmp"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    dn = systemfunc("dirname " + fn)
    bn = systemfunc("basename " + fn)
    bn = str_sub_str(bn, date, "YYYYMM")

    str = var + "," + dn + "/" + bn

    system("echo " + str_get_dq + str + str_get_dq + " >> " + VFATMP)

    leave_msg(scriptname, funcname, 4)

end

;;#############################################################################
undef("get_vfatmp")
function get_vfatmp(var:string,
                    date:string)
;;
;; Arguments
;;    var: input variable.
;;    date: selected date (YYYYMM format).
;;
;; Return value
;;    The name of the file containing the variable.
;;
;; Description
;;    Reads the variable-file association from the VFATMP file.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150324-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, xx, vv, fn, id, dn, bn
begin

    funcname = "get_vfatmp"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    if (.not.isfilepresent_esmval(VFATMP)) then
        leave_msg(scriptname, funcname, 4)
        return("err")
    end if

    xx = readAsciiTable(VFATMP, 1, "string", 0)
    vv = str_get_field(xx(:, 0), 1, ",")
    fn = str_get_field(xx(:, 0), 2, ",")

    ;; Check if variable was read before (both var and var_ave)
    id = ind(var.eq.vv)
    if (.not.all(ismissing(id))) then
        dn = systemfunc("dirname " + fn(id))
        bn = systemfunc("basename " + fn(id))
        bn = str_sub_str(bn, "YYYYMM", date)
        str = dn + "/" + bn
        leave_msg(scriptname, funcname, 4)
        return(str)
    else
        id = ind(var + "_ave".eq.vv)
        if (.not.all(ismissing(id))) then
            dn = systemfunc("dirname " + fn(id))
            bn = systemfunc("basename " + fn(id))
            bn = str_sub_str(bn, "YYYYMM", date)
            str = dn + "/" + bn
            var = var + "_ave"  ; propagate updated name to find_and_read_var
            leave_msg(scriptname, funcname, 4)
            return(str)
        end if
    end if

    ;; Variable not read before
    leave_msg(scriptname, funcname, 4)
    return("err")

end

;;#############################################################################
undef("find_and_read_var")
function find_and_read_var(var:string,
                           date:string,
                           err:logical)
;;
;; Arguments
;;    var: input variable.
;;    date: selected date (YYYYMM format).
;;    err: logical.
;;
;; Return value
;;    An array with the selected variable.
;;
;; Description
;;    Reads the variable after finding the EMAC stream/channel containing it.
;;    When more than one stream/channel is found, the selection is made based
;;    on the following priority list:
;;        time coordinate matching the FIELD (monthly mean or daily output)
;;        data already on pressure (for T3M, T2Mz) or hybrid (for T3M) levels
;;        tracer_gp and tr_ channels
;;        first file in the list
;;    The err flag is to select whether the code should exit when the variable
;;    is not found or simply return an error string to the calling environment.
;;    To improve the performance, once a variable is found (usually in the
;;    first step of the time loop), the association to the file is written in
;;    a temporary file (procedure put_vfatmp) and can be read in the subsequent
;;    steps of the time loop (function get_vfatmp).
;;
;; Caveats
;;    Priority check hybrid vs. pressure to be implemented.
;;
;; References
;;
;; Modification history
;;    20150714-A_laue_ax: now saving all names of the files processed to the
;;                        global list gFileList (defined in reformat_EMAC_main.ncl)
;;                        in order to write the list of the files to the log-file
;;    20150322-A_righ_ma: a few modifications for processing REF-C1.
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, pfn, flist, lvar, not_ave, date, ii, \
    f, vlist, files, p, l_files, ii, tmp, dtime
begin

    funcname = "find_and_read_var"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Create a local copy
    lvar = var

    ;; Check if the variable was already found in a previous step
    pfn = get_vfatmp(lvar, date)
    if (pfn.ne."err") then
        f = addfile(pfn, "r")
        info_output("Reading " + lvar + " from " + pfn, verbosity, 2)
        xx = f->$lvar$
        leave_msg(scriptname, funcname, 4)
; a_laue_ax+
        gFileList(gNFiles) = pfn
        gNFiles = gNFiles + 1
        if (gNFiles.ge.dimsizes(gFileList)) then
            error_msg("f", "reformat_EMAC_func.ncl", funcname, "size of array " + \
                  "gFileList too small. Increase size in reformat_EMAC_main.ncl")
        end if
; a_laue_ax-
        return(xx)
    end if

    ;; Variables for which the _ave version (MESSy2) shall not be considered
    not_ave = (/"aps", "press", "tm1", "grvol", "gboxarea"/)

    ;; Define files list
    ;; Naming convention: <EXPER>_<YYYYMM>_<CHANNEL>[OTHER FLAGS].nc
    flist = systemfunc("find " + INPATH + " -type f -name '" + EXPER + "*_" + \
                       date + "_*.nc' 2>/dev/null")

    ;; If all missing: give up
    if (all(ismissing(flist))) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, "no EMAC files " + \
                  "found for date " + date)
    end if

    ;; Get variables list
    flist@_FillValue = "missing"
    do ii = 0, dimsizes(flist) - 1

        f = addfile(flist(ii), "r")
        vlist = getfilevarnames(f)

        ;; Consider var_ave as possible MESSy2 stream (not for all vars)
        if (all(var.ne.vlist)) then
            if (any(var.eq.not_ave)) then
                flist(ii) = flist@_FillValue
                delete(vlist)
                continue
            else
                if (all(var + "_ave".ne.vlist)) then
                    flist(ii) = flist@_FillValue
                    delete(vlist)
                    continue
                end if
            end if
        end if
        delete(vlist)
        delete(f)

    end do

    ;; No files containg the variable found
    if (all(ismissing(flist))) then
        if (err) then
            error_msg("f", "reformat_EMAC_func.ncl", funcname, "variable " + \
                      var + " (" + var + "_ave) not found")
        else
            leave_msg(scriptname, funcname, 4)
            return("err")
        end if
    end if

    ;; List of files containing var
    files = flist(ind(.not.ismissing(flist)))
    files = str_sub_str(files, INPATH, "")
    files@_FillValue = "missing"
    delete(flist)

    ;; For hybrid coordinates and any file is ok
    if (any(var.eq.(/"hyam", "hybm", "hyai", "hybi"/))) then
        sf = INPATH + files(0)
        f = addfile(sf, "r")
        put_vfatmp(var, date, sf)
        info_output("Reading " + var + " from " + sf, verbosity, 1)
        xx = f->$var$
        leave_msg(scriptname, funcname, 4)
; a_laue_ax+
        gFileList(gNFiles) = sf
        gNFiles = gNFiles + 1
; a_laue_ax-
        return(xx)
    end if

    ;; Only 1 file found: read var and return
    if (dimsizes(files).eq.1) then
        sf = INPATH + files
        f = addfile(sf, "r")
        lvar = var
        if (.not.isfilevar(f, var)) then
            lvar = lvar + "_ave"
        end if
        put_vfatmp(lvar, date, sf)
        info_output("Reading " + lvar + " from " + sf, verbosity, 1)
        xx = f->$lvar$
        leave_msg(scriptname, funcname, 4)
; a_laue_ax+
        gFileList(gNFiles) = sf
        gNFiles = gNFiles + 1
; a_laue_ax-
        return(xx)
    end if

    ;; More than 1 file found: sort out based on time coordinate
    lfiles = new(dimsizes(files), logical)
    lfiles = False
    do ii = 0, dimsizes(files) - 1
        f = addfile(INPATH + files(ii), "r")
        if (isStrSubset(FIELD, "M")) then  ; monthly mean
            if (dimsizes(f->time).eq.1) then
                lfiles(ii) = True
            end if
        else  ; daily mean or instantaneous output
            if (dimsizes(f->time).gt.1) then
                lfiles(ii) = True
            end if
        end if
    end do
    if (num(lfiles).eq.1) then
        sf = INPATH + files(ind(lfiles))
        f = addfile(sf, "r")
        lvar = var
        if (.not.isfilevar(f, var)) then
            lvar = lvar + "_ave"
        end if
        put_vfatmp(lvar, date, sf)
        info_output("Reading " + lvar + " from " + sf, verbosity, 1)
        xx = f->$lvar$
        leave_msg(scriptname, funcname, 4)
; a_laue_ax+
        gFileList(gNFiles) = sf
        gNFiles = gNFiles + 1
; a_laue_ax-
        return(xx)
    end if

    ;; More than 1 file found: priority to data on pressure/hybrid levels
    ;; TO DO
;    if (isStrSubset(FIELD, "T3")) then
;        PRIORITY TO PRESSURE
;    end if
;    if (isStrSubset(FIELD, "T2Ms")) then
;        PRIORITY TO HYBRID
;    end if

    ;; More than 1 file found: priority to tracer_gp and _tr_
    tmp = files(ind(lfiles))
    delete(files)
    files = tmp
    delete(tmp)
    delete(lfiles)
    lfiles = new(dimsizes(files), logical)
    lfiles = False
    do ii = 0, dimsizes(files) - 1
        if (isStrSubset(files(ii), "tracer_gp").or. \
            isStrSubset(files(ii), "_tr_")) then
            lfiles(ii) = True
        end if
    end do
    if (num(lfiles).eq.1) then
        sf = INPATH + files(ind(lfiles))
        f = addfile(sf, "r")
        lvar = var
        if (.not.isfilevar(f, var)) then
            lvar = lvar + "_ave"
        end if
        put_vfatmp(lvar, date, sf)
        info_output("Reading " + lvar + " from " + sf, verbosity, 1)
        xx = f->$lvar$
        leave_msg(scriptname, funcname, 4)
; a_laue_ax+
        gFileList(gNFiles) = sf
        gNFiles = gNFiles + 1
; a_laue_ax-
        return(xx)
    end if

    ;; More than 1 file found: pick the first one
    info_output("Variable " + var + " found in more than 1 file: ", \
                verbosity, 2)
    do ii = 0, dimsizes(files) - 1
        info_output("     " + files(ii), verbosity, 2)
    end do
    sf = INPATH + files(0)
    f = addfile(sf, "r")
    lvar = var
    if (.not.isfilevar(f, var)) then
        lvar = lvar + "_ave"
    end if
    put_vfatmp(lvar, date, sf)
    info_output("Reading " + lvar + " from " + sf, verbosity, 1)
    xx = f->$lvar$
    leave_msg(scriptname, funcname, 4)
; a_laue_ax+
    gFileList(gNFiles) = sf
    gNFiles = gNFiles + 1
; a_laue_ax-
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
;; References
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, rank, date, opt, days, ndays, dd, idx, \
    dum
begin

    funcname = "time_mean"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
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
                error_msg("f", "reformat_EMAC_func.ncl", funcname, \
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
                error_msg("f", "reformat_EMAC_func.ncl", funcname, \
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

    error_msg("f", "reformat_EMAC_func.ncl", funcname, "unrecognized " + \
              "time mean option " + mean)

end

;;#############################################################################
undef("extract_var")
if (isfilepresent_esmval(RECIPEFILE)) then
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    loadscript(RECIPEFILE)
end if
function extract_var(name:string,
                     vert:string,
                     date:string,
                     stname:string)
;;
;; Arguments
;;    name: EMAC name of the variable.
;;    vert: flag for the vertical coordinate.
;;    date: selected date.
;;    stname: standard CMOR name of the variable.
;;
;; Return value
;;    The selected variable reformatted according to vert and FIELD.
;;
;; Description
;;    Extracts the requested variable from the given stream and reformats it to
;;    the requested type (T3M, T2Mz etc.).
;;
;; Caveats
;;    Column integration currently implemented for ozone only and NOT tested.
;;
;; References
;;
;; Modification history
;;    20130506-A_gott_kl: changes for data already on pressure levels.
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, xx, recog_vert_coords, recog_press, \
    recog_units, l_inst, l_day, l_month, l_4dim, l_zonal, l_vert, l_press, \
    vert_coord, levels, l_press, l_colum, l_surf, ndims, hyam, hybm, aps, \
    vol, pp, tt, area, xt, maxplev, slev, reordxt, caps, ret, xp, xz, xc, \
    nlev, xs
begin

    funcname = "extract_var"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read variable
    if (name.ne."#RECIP") then  ; EMAC variable
        xx = find_and_read_var(name, date, True)
    else
        if (isfilepresent_esmval(RECIPEFILE)) then
            xx = EMAC_recipe(date)  ; complex variable
        else
            error_msg("f", "reformat_EMAC_func.ncl", funcname, \
                      "Recipe file for " + stname + " not found")
        end if
    end if

    ;; Recognized vertical coordinates, pressure names and units
    recog_vert_coords = (/"lev", "mlev", "plev", "lev_p", "lev_2"/)
    recog_press = (/"air_pressure", "pressure", "press", "PRESS", "PRESSURE"/)
    recog_units = (/"Pa", "hPa", "PA", "mbar", "MBAR"/)

    ;; Logical flags
    l_inst  = isStrSubset(FIELD, "I")  ; instantaneous output?
    l_day   = isStrSubset(FIELD, "D")  ; daily mean?
    l_month = isStrSubset(FIELD, "M")  ; monthly mean?
    l_4dim  = isStrSubset(FIELD, "T3")  ; (time,lev,lat,lon)?
    l_zonal = isStrSubset(FIELD, "z")  ; zonal mean?
    l_vert  = any(iscoord(xx, recog_vert_coords))  ; vertical coordinate?
    l_press = l_4dim.or.l_zonal  ; h2p interpolation?
    l_colum = isStrSubset(FIELD, "s").and.l_vert.and.(vert.eq."C")  ; column?
    l_surf  = isStrSubset(FIELD, "s").and.l_vert.and.(vert.eq."S")  ; surface?

    ;; Check if already on pressure levels?
    if (l_vert) then
        vert_coord = recog_vert_coords(ind(iscoord(xx, recog_vert_coords)))
        levels = xx&$vert_coord$
        if(inlist(levels@units, recog_units) .or. \
           inlist(levels@long_name, recog_press) .or. \
           inlist(levels@standard_name, recog_press)) then
            l_press = False
            info_output("reformat_EMAC_func: data are already on pressure " + \
                        "levels", verbosity, 2)
            xx!1 = "lev_p"
            if (levels@units.eq."Pa" .or. levels@units.eq."PA") then
                ;; Force units conversion to trigger standardization
                xx&lev_p = xx&lev_p * 0.01  ; [Pa] --> [hPa]
            end if
        end if
        delete([/recog_press, recog_units, levels, vert_coord/])
    end if
    delete(recog_vert_coords)

    ;; Number of dimensions
    ndims = dimsizes(dimsizes(xx))

    ;; Check consistency
    if (l_4dim.and.ndims.lt.4) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, "requested " + \
                  "variable has only " + ndims + " dimensions, cannot " + \
                  "apply a " + FIELD + " field")
    end if
    if (l_zonal.and..not.l_vert) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, "no vertical " + \
                  "coordinate exists, cannot apply a " + FIELD + " field")
    end if

    ;; Read surface pressure (for h2p interpolation)
    if (l_press) then
        hyam = find_and_read_var("hyam", date, True)
        hybm = find_and_read_var("hybm", date, True)
        if (l_month) then
            aps  = time_mean(find_and_read_var("aps", date, True), "M")
        end if
        if (l_day)
            aps = time_mean(find_and_read_var("aps", date, True), "D")
        end if
        if (l_inst)
            aps = find_and_read_var("aps", date, True)
        end if
    end if

    ;; Read volume, area, pressure and temperature
    if (l_colum) then
        if (l_month) then
            vol  = time_mean(find_and_read_var("grvol", date, True), "M")
            pp   = time_mean(find_and_read_var("press", date, True), "M")
            tt   = time_mean(find_and_read_var("tm1", date, True), "M")
            area = time_mean(find_and_read_var("gboxarea", date, True), "M")
        end if
        if (l_day) then
            vol  = time_mean(find_and_read_var("grvol", date, True), "D")
            pp   = time_mean(find_and_read_var("press", date, True), "D")
            tt   = time_mean(find_and_read_var("tm1", date, True), "D")
            area = time_mean(find_and_read_var("gboxarea", date, True), "D")
        end if
        if (l_inst) then
            vol  = find_and_read_var("grvol", date, True)
            pp   = find_and_read_var("press", date, True)
            tt   = find_and_read_var("tm1", date, True)
            area = find_and_read_var("gboxarea", date, True)
        end if
    end if

    ;; Start processing chain

    ;; ### TIME (xx --> xt)
    if (l_inst) then
        xt = xx
    end if
    if (l_day) then
        xt = time_mean(xx, "D")
    end if
    if (l_month) then
        xt = time_mean(xx, "M")
    end if
    delete(xx)

    ;; ### HYBRID-TO-PRESSURE (xt --> xp)
    if (l_press) then

        ;; Limit PLEVS range to the uppermost level of the input variable
        maxplev = hyam(0) / 100.
        slev = PLEVS(ind(PLEVS.ge.maxplev))

        ;; Check for missing values, as vinth2p cannot handle them
        if (num(ismissing(xt)).ne.0) then
            error_msg("w", "reformat_EMAC_func.ncl", funcname, "data to " + \
                      "be interpolated from hybrid to pressure levels " + \
                      "contains missing values.")
        end if

        ;; Interpolate
        if (any(getvardims(xt).eq."mode")) then  ; special case for size modes
            if (iscoord(xt, "mlev")) then
                reordxt = xt(mode|:, time|:, mlev|:, lat|:, lon|:)
            end if
            if (iscoord(xt, "lev")) then
                reordxt = xt(mode|:, time|:, lev|:, lat|:, lon|:)
            end if
            caps = conform(dim_avg_n_Wrap(reordxt, 2), aps, (/1, 2, 3/))
            ret = vinth2p(reordxt, hyam / 101325., hybm, slev, caps, 2, \
                          1013.25, 1, False)
            xp = ret(time|:, lev_p|:, lat|:, lon|:, mode|:)
            delete(reordxt)
            delete(caps)
            delete(ret)
        else  ; standard case
            xp = vinth2p(xt, hyam / 101325., hybm, slev, aps, 2, 1013.25, 1, \
                         False)
        end if
        xp&lev_p@units = "hPa"
        xp@units = xt@units

        ;; Print range of interpolated values and compare with _FillValue
        if (num(ismissing(xt)).ne.0) then
            error_msg("w", "reformat_EMAC_func.ncl", funcname, "if " + \
                      xp@_FillValue + " is close to " + min(xp) + " or " + \
                      max(xp) + ", consider to activate fix_FillValue.")
        end if

    else
        xp = xt
    end if
    delete(xt)

    ;; ### ZONAL MEAN (xp --> xz)
    if (l_zonal) then
        xz = dim_avg_n_Wrap(xp, ind("lon".eq.getvardims(xp)))
    else
        xz = xp
    end if
    delete(xp)

    ;; ### COLUMN INTEGRATION (xz --> xc)
    if (l_colum) then
        error_msg("f", "reformat_EMAC_func.ncl", funcname, "feature not " + \
                  "yet implemented")
    else
        xc = xz
    end if
    delete(xz)

    ;; ### SURFACE VALUE (xc --> xs)
    if (l_surf) then
        if (iscoord(xc, "lev")) then
            nlev = dimsizes(xc&lev)
        end if
        if (iscoord(xc, "mlev")) then
            nlev = dimsizes(xc&mlev)
        end if
        if (any(getvardims(xc).eq."mode")) then  ; special case for size modes
            xs = xc(:, nlev - 1, :, :, :)
        else
            xs = xc(:, nlev - 1, :, :)
        end if
    else
        xs = xc
    end if
    delete(xc)

    leave_msg(scriptname, funcname, 4)
    return(xs)

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
;; References
;;
;; Modification history
;;    20130502-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, rank, ctime, cplev, clat, clon, tmp
begin

    funcname = "format_coord"
    scriptname = "reformat_scripts/EMAC/reformat_EMAC_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    rank = dimsizes(dimsizes(var))

    ;; Time
    if (iscoord(var, "time")) then
        ;; Check calendar
        if (isatt(var&time, "calendar")) then
            if (var&time@calendar.ne."standard".and.\
                var&time@calendar.ne."gregorian") then
                error_msg("f", "reformat_EMAC_func.ncl", funcname, \
                          "calendar type " + var&time@calendar + \
                          " not supported")
            end if
        end if
        ctime = time_attrib(var&time)
        delete(var&time)
        var&time = ctime
    end if

    ;; Level
    if (iscoord(var, "lev_p")) then
        var!1 = "plev"
        var&plev = var&plev * 100  ; [hPa] --> [Pa]
        ;; Check ordering (must be monotonically decreasing)
        if (isMonotonic(var&plev) .eq. 0) then
            error_msg("f", "reformat_EMAC_func.ncl", funcname, \
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
                error_msg("f", "reformat_EMAC_func.ncl", funcname, \
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
            if (rank.eq.5 .and. pos.eq.2) then
                var = var(:, :, ::-1, :, :)  ; time, lev, lat, lon, mode
            else if (rank.eq.4 .and. pos.eq.2) then
                var = var(:, :, ::-1, :)  ; time, lev, lat, lon
            else if (rank.eq.4 .and. pos.eq.1) then
                var = var(:, ::-1, :, :)  ; time, lat, lon, mode
            else if (rank.eq.3 .and. pos.eq.1) then
                var = var(:, ::-1, :)  ; time, lat, lon
            else
                error_msg("f", "reformat_EMAC_func.ncl", funcname, \
                          "invalid rank or lat position")
            end if
            end if
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
