;;#############################################################################
;; GENERAL-PURPOSE AUXILIARY NCL ROUTINES
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure copy_VarCoords_l1
;;    procedure check_min_max_models
;;    function tstep
;;    function ncdf_varname
;;    function get_work_dir
;;    function get_ncdf_name
;;    function get_ncdf_dir
;;    function ncdf_read
;;    function ncdf_write
;;    function ncdf_att
;;    function ncdf_define
;;    function att2var_default
;;    function att2var
;;    function bname
;;    function PRIORITY
;;    function get_output_dir
;;    function basename
;;    function extract_years
;;    function extend_var_at
;;    function copy_CoordNames_n
;;    function empty_str
;;    function write_info
;;    function remove_index
;;    procedure tool_stop
;;    function filter_attrs
;;    procedure write_ignore_warnings
;;    function get_ref_model_idx
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_data/ncl.interface"
load "./interface_scripts/messaging.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"

;;#############################################################################
undef("copy_VarCoords_l1")
procedure copy_VarCoords_l1(var_from,
                            var_to)
;;
;; Arguments
;;    var_from: An array of any dimensionality.
;;              Must have named dimensions and coordinate variables.
;;    var_to: An array of any dimensionality.
;;
;; Description
;;    This procedure copies all named dimensions and coordinate variables
;;    from one variable to another except for the leftmost dimension.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140703-A_gott_kl: Written for use with mymm in
;;                        time_operations (statistics.ncl).
;;
local funcname, scriptname, verbosity, dims_from, dims, rank_from, rank, \
    dimdiff, idimdiff
begin

    funcname = "copy_VarCoords_l1"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    dims_from = getvardims(var_from)
    rank_from = dimsizes(dims_from)
    dims = getvardims(var_to)
    rank = dimsizes(dims)

    if(rank.gt.1) then

        ;; Check input
        if(rank_from.ne.rank) then  ; check rank
            error_msg("f", scriptname, funcname, \
                      "var_from and var_to must have same rank")
        else  ; check dimensions
            dimdiff = dimsizes(var_from) - dimsizes(var_to)
            idimdiff = ind(dimdiff.ne.0)
            if(.not.all(ismissing(idimdiff))) then
                if(sum(idimdiff).ne.0) then  ; more or other than leftmost dim
                    error_msg("f", scriptname, funcname, \
                              "var_from and var_to must have same " + \
                              "dimension sizes, except the leftmost")
                end if
            end if
        end if

        ;; Loop over dimensions
        do i = 1, rank - 1
            var_to!i = dims_from(i)
            var_to&$dims_from(i)$ = var_from&$dims_from(i)$
        end do

    else

        info_output("No coordinates to process in copy_VarCoords_l1", \
                    verbosity, 2)

    end if

    leave_msg(scriptname, funcname, 6)

end

;;#############################################################################
undef("check_min_max_models")
procedure check_min_max_models(no_models[1]:integer,
                               min_no_models[1]:integer,
                               max_no_models[1]:integer,
                               diag_script[1])
;;
;; Arguments
;;    no_models: number of models.
;;    min_no_models: minimum no of models required for this diagnostic script.
;;    max_no_models: maximum no of models required for this diagnostic script.
;;    diag_script: name of diagnostic script.
;;
;; Description
;;    Checks if sufficient number of models are present for the given
;;    diagnostic, tries to abort with a intelligible message.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity
begin

    funcname = "check_min_max_models"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    if (no_models .lt. min_no_models) then
        error_msg("f", scriptname, funcname,\
                  "Too few models (=" + no_models + ") for " + diag_script \
                  + ". Mininum is " + min_no_models)
    else if (no_models .gt. max_no_models) then
        error_msg("f", scriptname, funcname,\
                  "Too many models (=" + no_models + ") for " + diag_script \
                  + ". Maximum is " + max_no_models)
    else
        no_operation()
    end if
    end if

    leave_msg(scriptname, funcname, 6)

end

;;#############################################################################
undef("tstep")
function tstep(utc:numeric)
;;
;; Arguments
;;    utc: UTC time as returned by cd_calendar
;;         (dimsizes(time) x 6, with y,m,d,h,m,s).
;;
;; Return value
;;    "yea"=yearly, "sea"=seasonal, "mon"=monthly, "d"=daily, "h", "m", "s".
;;
;; Description
;;    Determines time step size in given UTC time axis.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140211-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, t, dt, nt, y, m, d, h, n, s, option, \
    units
begin

    funcname = "tstep"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    ;; extract vectors
    y = toint(utc(:, 0))
    m = toint(utc(:, 1))
    d = toint(utc(:, 2))
    h = toint(utc(:, 3))
    n = toint(utc(:, 4))  ; minutes
    s = toint(utc(:, 5))

    ;; convert to seconds after y(0)-01-01 00:00:00
    units = "seconds after " + tostring(y(0)) + "-01-01 00:00:00"
    option = 0
    option@calendar = "standard"
    t = cd_inv_calendar(y, m, d, h, n, s, units, option)

    ;; create vector of tstep lengths
    nt = dimsizes(t)
    dt = tolong(t(1:nt - 1)) - tolong(t(0:nt - 2))

    ;; check mean dt
    mdt = avg(dt)
    if (mdt.eq.1) then
        result = "s"
    else if (mdt.eq.60) then
        result = "m"
    else if (mdt.eq.3600) then
        result = "h"
    else if (mdt.eq.86400) then
        result = "d"
    else if (mdt.ge.2419200 .and. mdt.le.2678400) then
        result = "mon"  ; 28 ... 31 days
    else if (mdt.ge.7776000 .and. mdt.le.7948800) then
        result = "sea"  ; 90 ... 92 days
    else if (mdt.ge.31536000 .and. mdt.le.31622400) then
        result = "yea"  ; 365 ... 366 days
    else
        error_msg("f", scriptname, funcname, "could not determine " + \
                  "step length for mdt = " + mdt)
    end if
    end if
    end if
    end if
    end if
    end if
    end if

    ;; check if all dt match the result
    check_flag = False
    if (result.eq."s" .or. result.eq."m" .or. \
        result.eq."h" .or. result.eq."d") then
        if (max(dt) - min(dt) .ne. 0) then
            check_flag = True
        end if
    else if (result.eq."mon") then
        if (max(dt) - min(dt) .gt. 259200) then
            check_flag = True
        end if  ; 31-28 = 3 days
    else if (result.eq."sea") then
        if (max(dt) - min(dt) .gt. 172800) then
            check_flag = True
        end if  ; 92-90 = 2 days
    else if (result.eq."yea") then
        if (max(dt) - min(dt) .gt. 86400) then
            check_flag = True
        end if  ; 366-365 = 1 days
    end if
    end if
    end if
    end if
    if (check_flag) then
        error_msg("f", scriptname, funcname, "range of dt is not " + \
                  "consistent with " + result)
    end if

    leave_msg(scriptname, funcname, 6)
    return(result)

end

;;#############################################################################
undef("ncdf_varname")
function ncdf_varname(data)
;;
;; Arguments
;;    data: data with metadata (attributes).
;;
;; Return value
;;    A string with the variable name.
;;
;; Description
;;    Checks "data" for attribute "var" or "varname" or "name" or "long_name"
;;    (order = priority) and returns the name as a string variable.
;;    Returns the string "var" by default.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20131108-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, data, result, target_atts, att
begin

    funcname = "ncdf_varname"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    ;; Attributes names to check for, ordered by increasing priority
    target_atts = (/"long_name", "var_long_name", \
                    "var_name", "name", "varname", "var"/)
    result = "var"   ; Default
    do i = 0, dimsizes(target_atts) - 1
        if (isatt(data, target_atts(i))) then
            att = data@$target_atts(i)$
            if (typeof(att) .eq. "string" .and. dimsizes(att) .eq. 1) then
                result = att
            end if
            delete(att)
        end if
    end do

    leave_msg(scriptname, funcname, 6)
    return(result)

end

;;#############################################################################
undef("get_work_dir")
function get_work_dir()
;;
;; Arguments
;;
;; Return value
;;    A string with the full path of work_dir.
;;
;; Description
;;    Fetches $ESMValTool_wrk_dir and checks that it ends with a slash.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140414-A_righ_ma: bug fixed and unnecessary if case removed.
;;    20131108-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, curr_work_dir
begin

    funcname = "get_work_dir"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    curr_work_dir = getenv("ESMValTool_wrk_dir")

    ;; Make sure it ends with a slash
    if (str_get_cols(curr_work_dir, -1, -1).ne."/") then
        curr_work_dir = curr_work_dir + "/"
    end if

    leave_msg(scriptname, funcname, 10)
    return(curr_work_dir)

end

;;#############################################################################
undef("get_ncdf_name")
function get_ncdf_name(diag_script[1]:string,
                       var[1]:string)
;;
;; Arguments
;;    diag_script: name of diagnostic script(s).
;;    var: variable name.
;;
;; Return value
;;    A string with the full outfile name (incl. path)
;;    output_dir + diag_script_base + var + ".nc".
;;
;; Description
;;    Fetches work_dir, diag_script_base, output_dir via other scripts.
;;
;; Caveats
;;    diag_script may need to be extended by other contributing scripts.
;;
;; References
;;
;; Modification history
;;    20131108-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, work_dir, diag_script_base, \
    diag_script, output_dir, ncdf_dir
begin

    funcname = "get_ncdf_name"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    ;; Collect parts of outfile name
    work_dir = get_work_dir()
    diag_script_base = basename(diag_script)
    ncdf_dir = get_ncdf_dir(work_dir, diag_script_base)

    ;; Assemble final name
    outfile = ncdf_dir + diag_script_base + "_" + var + ".nc"

    leave_msg(scriptname, funcname, 10)
    return(outfile)

end

;;#############################################################################
undef("get_ncdf_dir")
function get_ncdf_dir(work_dir[1]:string,
                      subfolder[1]:string)
;;
;; Arguments
;;    work_dir: root path for intermediate files: expected to end with a slash.
;;    subfolder: subfolder for diag_script.
;;
;; Return value
;;    A string with the full path for netCDF output.
;;
;; Description
;;    Concatenates/generates work_dir + subfolder.
;;
;; Caveats
;;    diag_script may need to be extended by other contributing scripts.
;;
;; References
;;
;; Modification history
;;    20131108-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, work_dir, subfolder, ncdf_dir
begin

    funcname = "get_ncdf_dir"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    ncdf_dir = work_dir + subfolder + "/"
    if (.not. isfilepresent_esmval(ncdf_dir)) then
        system("mkdir -p " + ncdf_dir)
    end if

    leave_msg(scriptname, funcname, 10)
    return(ncdf_dir)

end

;;#############################################################################
undef("ncdf_read")
function ncdf_read(infile, varname)
;;
;; Arguments
;;    infile: full path & name of the netCDF file to be read.
;;
;; Return value
;;    Data with metadata:
;;        if varname="dummy" AND only one var in file: return that var
;;        if varname="default" AND only one var in file: return that var
;;        if varname="dummy" AND more than one var in file: exit
;;        if varname="default" AND more than one var in file: exit
;;        if varname=var AND varname exists in file: return var varname
;;        if varname=var AND varname does not exist: exit
;;
;; Description
;;    http://www.ncl.ucar.edu/Applications/r-support.shtml
;;
;; Caveats
;;    A variable is considered a coordinate variable, if it has the same name
;;    as first dimension (not sure if this is a satisfactory criterion).
;;    The netCDF format does not support data type "logical" and converts True
;;    (logical) to 1 (integer), False (logical) to 0 (integer), see
;;    http://www.ncl.ucar.edu/Support/talk_archives/2009/0367.html
;;    A workaround could be to convert logical data to string before writing to
;;    netCDF, and adjust reading routine accordingly.
;;    Currently ncdf_read interprets variable attributes of type "integer" and
;;    dimsizes = 1 & value 0 or 1 as logicals.
;;    If a named coordinate was converted to a character array in ncdf_write it
;;    is converted back to string and appended to the main variable. It only
;;    works if the variable has an attribute named_coordinate equal 1.
;;
;; References
;;
;; Modification history
;;    20131111-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, data, infile, ncdf, varname, i, ivar, \
    nvar, vdims, var, atts, temp, defaults
begin

    funcname = "ncdf_read"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    ;; Open file for reading
    if (isfilepresent_esmval(infile)) then
        ncdf = addfile(infile, "r")
    else
        error_msg("f", scriptname, funcname, infile + " does not exist")
    end if

    ;; Determine (number of) variables, exclusive coordinates
    vNames = getfilevarnames(ncdf)
    nvar = new(dimsizes(vNames), integer)
    nvar(:) = 1  ; used for flagging non-coordinate vars
    do ivar = 0, dimsizes(vNames) - 1
        vdims = getfilevardims(ncdf, vNames(ivar))
        if (vdims(0) .eq. vNames(ivar)) then
            nvar(ivar) = 0  ; coordinate variable
        end if
        delete(vdims)
    end do
    var = vNames(ind(nvar .eq. 1))  ; non-coordinate vars

    ;; Cases for reading
    defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
    if (inlist(varname, defaults)) then  ; allow synonyms
        varname = "dummy"  ; unify
    end if
    if (varname .eq. "dummy" .and. dimsizes(var) .eq. 1) then
        data = ncdf->$var$
    else if (varname .eq. "dummy" .and. dimsizes(var) .ne. 1) then
        error_msg("f", scriptname, funcname, "please select variable " + \
                  "from " + infile)
    else if (varname .ne. "dummy" .and. inlist(varname, vNames)) then
        i = ind(vNames .eq. varname)
        data = ncdf->$vNames(i)$
    else if (varname .ne. "dummy" .and. .not.inlist(varname, vNames)) then
        error_msg("f", scriptname, funcname, "variable " + varname + \
                  " not in " + infile)
    else
        error_msg("f", scriptname, funcname, "case not covered")
    end if
    end if
    end if
    end if

    ;; Convert attributes to type logical, if they meet above conditions
    atts = getvaratts(data)
    do i = 0, dimsizes(atts) - 1
        if (typeof(data@$atts(i)$) .eq. "integer" .and. \
            dimsizes(data@$atts(i)$) .eq. 1) then
            if (data@$atts(i)$ .eq. 1 .or. data@$atts(i)$ .eq. 0) then
                temp = new(1, logical)
                temp = data@$atts(i)$
                delete(data@$atts(i)$)
                data@$atts(i)$ = temp
                delete(temp)
            end if
        end if
    end do

    ;; Converts named coordinates back to strings and append them to the data
    coord_ind = ind(nvar.eq.0)
    do j = 0, dimsizes(vNames(coord_ind)) - 1
        new_coord = ncdf->$vNames(coord_ind(j))$
        if ((isatt(new_coord, "named_coordinate")).and. \
            (new_coord@named_coordinate.eq.1)) then
            data&$vNames(coord_ind(j))$ = chartostring(new_coord)
        end if
        delete(new_coord)
    end do

    leave_msg(scriptname, funcname, 10)
    return(data)

end

;;#############################################################################
undef("ncdf_write")
function ncdf_write(data,
                    outfile_in)
;;
;; Arguments
;;    data: data with metadata (attributes) to be written to a netCDF file.
;;    outfile_in: full path & name of the netCDF file to be written. It may
;;                may contain an attribute "existing" with the values
;;                "overwrite" (default) / "append" / "rename"
;;                default is used for "default", "dummy" or non-string values.
;; Return value
;;    A string with the outfile.
;;
;; Description
;;    Creates or modifies outfile, adding a single variable. All attributes of
;;    "data" are added to the variable in the netCDF file. It may rename an
;;    existing file to keep it.
;;    Two options for netCDF output, one disabled (see "Option"), change if
;;    needed.
;;    If data has coordinates which are named with strings, these names are
;;    converted to characeter arrays and the coordinates are appended as two
;;    dimensional variables.
;;
;; Caveats
;;    Appending not tested.
;;    netCDF does not support data type "logical" and converts True (logical)
;;    to 1 (integer), False (logical) to 0 (integer).
;;    A workaround could be to convert logical data to string before writing to
;;    netCDF, and adjust reading routine accordingly.
;;    Currently ncdf_read interprets variable attributes of type "integer" and
;;    dimsizes = 1 & value 0 or 1 as logicals.
;;
;; References
;;    http://www.ncl.ucar.edu/Document/Manuals/Ref_Manual/NclFormatSupport.shtml#NetCDF
;;    http://www.ncl.ucar.edu/Applications/o-netcdf.shtml
;;    http://www.ncl.ucar.edu/Support/talk_archives/2009/0367.html
;;
;; Modification history
;;    20131107-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, data, outfile, outfile_in, ncdf, varname
begin

    funcname = "ncdf_write"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    ;; Retrieve special attributes
    diag_script = att2var(data, "diag_script")
    varname = att2var(data, "var")

    ;; Determine outfile
    defaulting = (/"default", "dummy"/)
    if (typeof(outfile_in) .ne. "string" .or. \
        .not.ismissing(ind(defaulting .eq. outfile_in))) then
        diag_script = att2var(data, "diag_script")
        outfile = get_ncdf_name(diag_script, varname)
    else
        outfile = outfile_in
    end if

    ;; Determine mode for dealing with existing files
    if (isatt(outfile, "existing")) then
        existing = outfile@existing
    else
        existing = "overwrite"   ; Default
    end if

    ;; Open outfile (see NCL documentation for "addfile")
    if (isfilepresent_esmval(outfile)) then  ; Outfile already exists
        if (existing .eq. "rename") then
            sysdate = systemfunc("date +%Y%m%d%H%M%S")
            system("mv " + outfile + " " + outfile + "." + sysdate)
            ncdf = addfile(outfile, "c")
            info_output(" Keeping existing " + outfile + " with suffix ." + \
                        sysdate, verbosity, 2)
        else if (existing .eq. "append") then
            ncdf = addfile(outfile, "w")
        else
            system("rm -f " + outfile)
            ncdf = addfile(outfile, "c")
        end if
        end if
    else  ; Outfile does not exist yet
        ncdf = addfile(outfile, "c")
    end if

    ;; Write to netCDF
    ;; Option 1: q&d
;    ncdf->$varname$ = data
    ;; Option 2: write with explicit file definition
    dimNames = getvardims(data)
    ncdf = ncdf_define(ncdf, data)
    do i = 0, dimsizes(dimNames) - 1
        if (typeof(data&$dimNames(i)$).eq."string") then
            ;; Convert string array to two dimensional character array
            char_array = stringtochar(data&$dimNames(i)$)
            ncdf->$dimNames(i)$ = (/char_array/)
            delete(char_array)
        else
            ncdf->$dimNames(i)$ = (/data&$dimNames(i)$/)
        end if
    end do
    ncdf->$varname$ = (/data/)

    info_output("Wrote variable " + varname + " to " + outfile, verbosity, 1)

    leave_msg(scriptname, funcname, 10)
    return(outfile)

end

;;#############################################################################
undef("ncdf_att")
function ncdf_att(infile[1],
                  att[1]:string)
;;
;; Arguments
;;    infile: full path & name of the netCDF file to be read OR a file
;;            variable already created by addfile(...).
;;    att: attribute to be read.
;;
;; Return value
;;    Contents of attribute att (could be any type supported by netCDF).
;;
;; Description
;;    Checks existence of att in infile, then retrieves contents.
;;
;; Caveats
;;    str_infile does not resolve the file name for a file variable.
;;
;; References
;;
;; Modification history
;;    20140131-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, f, file_atts, i, result
begin

    funcname = "ncdf_att"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    if (typeof(infile).eq."string") then
        f = addfile(infile, "r")
        str_infile = infile
    else if (typeof(infile).eq."file") then
        f = infile
        str_infile = "file associated with file variable argument."
    else
        error_msg("f", scriptname, funcname, "wrong type of infile " + \
                  "argument")
    end if
    end if
    file_atts = getvaratts(f)  ; retrieve global attributes
    if (any(file_atts.eq.att)) then
        i = ind(file_atts.eq.att)
        result = f@$file_atts(i)$
    else
        error_msg("w", scriptname, funcname, att + " not found in " + \
                  str_infile)
        result = "missing"
    end if

    leave_msg(scriptname, funcname, 10)
    return(result)

end

;;#############################################################################
undef("ncdf_define")
function ncdf_define(ncdf, data)
;;
;; Arguments
;;    ncdf: file specified.
;;    data: data with metadata (attributes).
;;
;; Return value
;;    Modifies the file specifier "ncdf".
;;
;; Description
;;
;; Caveats
;;    No unlimited dimensions supported (yet).
;;
;; References
;;    http://www.ncl.ucar.edu/Applications/method_2.shtml
;;
;; Modification history
;;    20131112-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, data, diag_script, var, fAtt, \
    dimNames, dimSzs, dimUnlim
begin

    funcname = "ncdf_define"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    diag_script = att2var(data, "diag_script")
    var = att2var(data, "var")

    ;; Explicitly declare file definition mode. Improve efficiency.
    setfileoption(ncdf, "DefineMode", True)

    ;; Create global attributes of the file
    fAtt               = True  ; assign file attributes
    fAtt@creation_date = systemfunc("date")
    fAtt@created_by    = "ESMValTool, " + diag_script

    fileattdef(ncdf, fAtt)   ; copy file attributes

    ;; Predefine the coordinate variables and their dimensionality
    ;; Note: to get an UNLIMITED record dimension, we set the dimensionality
    ;; to -1 (or the actual size) and set the dimension name to True.
    dimNames = getvardims(data)
    dim_dim = dimsizes(dimNames)
    dimSzs = new(dimsizes(dimNames), long)
    dimUnlim = new(dimsizes(dimNames), logical)
    dimUnlim(:) = False

    do i = 0, dim_dim - 1
        dimSzs(i) = dimsizes(data&$dimNames(i)$)
    end do

    do ii = 0, dim_dim - 1
        if (typeof(data&$dimNames(ii)$).eq."string") then
            ;; Add necessary infomation to the dim* arrays
            dim_coord = dimsizes(stringtochar(data&$dimNames(ii)$))
            ;; If coordinate has length 1
            if (dimsizes(dim_coord).eq.1) then
                dimc = dim_coord(0)
            else
                dimc = dim_coord(1)
            end if
            new_dimSzs = array_append_record(dimSzs, tolong(dimc), 0)
            new_dimUnlim = array_append_record(dimUnlim, False, 0)
            new_dimNames = array_append_record(dimNames, \
                                               dimNames(ii) + "_char", 0)

            ;; Delete and reattach to previous variable
            delete([/dimSzs, dimUnlim, dimNames/])
            dimSzs = new_dimSzs
            dimUnlim = new_dimUnlim
            dimNames = new_dimNames
            delete([/new_dimSzs, new_dimUnlim, new_dimNames/])
            delete(dim_coord)
            jj = 1  ; counter initialization for later use
        end if
    end do

    filedimdef(ncdf, dimNames, dimSzs, dimUnlim)
    ;; Predefine the the dimensionality of the variables to be written out
    do i = 0, dim_dim - 1
        if (typeof(data&$dimNames(i)$).eq."string") then
            ;; Attach second dimension (character array) to variable
            filevardef(ncdf, dimNames(i), "character", \
                       (/getvardims(data&$dimNames(i)$), \
                       dimNames(dim_dim + jj - 1)/))
            data&$dimNames(i)$@named_coordinate = True  ; set flag attribute
            jj = jj + 1
        else
            filevardef(ncdf, dimNames(i), typeof(data&$dimNames(i)$),\
                       getvardims(data&$dimNames(i)$))
        end if
    end do
    filevardef(ncdf, var, typeof(data), getvardims(data))

    ;; Copy attributes associated with each variable to the file
    do i = 0, dim_dim - 1
        ;; No attribute assignment to variable if no attributes are available
        if (.not.(all(ismissing(getvaratts(data&$dimNames(i)$))))) then
            filevarattdef(ncdf, dimNames(i), data&$dimNames(i)$)
        end if
    end do
    filevarattdef(ncdf, var, data)

    ;; Explicitly exit file definition mode. **NOT REQUIRED**
    setfileoption(ncdf, "DefineMode", False)

    leave_msg(scriptname, funcname, 10)
    return(ncdf)

end

;;#############################################################################
undef("att2var_default")
function att2var_default(data,
                         meta[1]:string,
                         default,
                         option)
;;
;; Arguments
;;    data: variable of any type with plot specific metadata.
;;    meta: string specifying the requested attribute
;;          "var", "res_" and "diag_file" get special treatment.
;;    default: variable of any type.
;;    option: not used yet.
;;
;; Return value
;;    Requested attribute of data is returned as variable.
;;    See function  att2var for details on 'magic' strings.
;;    Those get special treatment and may trigger a crash, if missing.
;;    Default is returned if (non magic) data@meta is missing
;;
;; Description
;;    Puts attributes of data into individual variables,
;;    and returns default if the requested attribute is missing
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140403-A_GoKl: written.
;;
local funcname, scriptname, verbosity, data, meta, result, default, option
begin

    funcname = "att2var_default"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    result = att2var(data, meta)
    if(dimsizes(result).eq.1) then
        if(result.eq.default_fillvalue("string")) then
            info_output("Using default: " + meta + " = " + tostring(default), \
                        verbosity, 4)
            result = default
            copy_VarMeta(default, result)
        end if
    end if

    leave_msg(scriptname, funcname, 10)
    return(result)

end

;;#############################################################################
undef("att2var")
function att2var(data,
                 meta[1]:string)
;;
;; Arguments
;;    data: variable of any type with plot specific metadata.
;;    meta: string specifying the requested attribute, "var", "res_" and
;;          "diag_file" get special treatment.
;;
;; Return value
;;    Requested attribute of data is returned as variable:
;;      var (string): name string, with "var_" attached as attributes
;;                    (without prefix).
;;      diag_file (string): name bases of registred scripts, joined into one
;;                          string.
;;      res (logical): attributes needed for plotting.
;;      other attributes (any type): as requested by meta.
;;
;; Description
;;    Puts attributes of data into individual variables.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20131107-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, data, meta, result, temp, atts, att, \
    iatt
begin

    funcname = "att2var"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    ;; Attribute "var"
    if (meta .eq. "var")
        if (isatt(data, meta)) then
            result = data@$meta$
            atts = getvaratts(data)
            do iatt = 0, dimsizes(atts) - 1
                att = atts(iatt)
                if (isStrSubset(att, "var_")) then
                    ;; Remove prefix "var_"
                    temp = str_get_cols(att, 4, strlen(att))
                    ;; Assign contents to variable "var"
                    result@$temp$ = data@$att$
                end if
            end do
        else
            error_msg("f", scriptname, funcname, get_script_name() + \
                      " requires variable attribute " + meta)
        end if

    ;; Attribute "diag_script"
    else if (meta .eq. "diag_script")
        if (isatt(data, meta)) then
            base = bname(data@$meta$)
            ;; Join names of all contributing scripts
            result = str_join(base, "__")
        else
            error_msg("f", scriptname, funcname, get_script_name() + \
                      " requires variable attribute " + meta)
        end if

    ;; Attribute "res_"
    ;; DEBUG info: different attribute types are not anticipated
    ;; -> check this if it crashes
    else if (meta .eq. "res_")
        result = True
        atts = getvaratts(data)
        do iatt = 0, dimsizes(atts) - 1
            att = atts(iatt)
            if (isStrSubset(att, "res_")) then
                ;; Remove prefix "res_"
                temp = str_get_cols(att, 4, strlen(att))
                ;; Assign contents to variable "res"
                result@$temp$ = data@$att$
            end if
        end do

    ;; Optional/generic attributes
    else
        if (isatt(data, meta)) then
            result = data@$meta$
        else
            error_msg("w", scriptname, funcname, "attribute " + meta + \
                      " not found")
            result = default_fillvalue("string")
        end if
    end if
    end if
    end if

    leave_msg(scriptname, funcname, 10)
    return(result)

end

;;#############################################################################
undef("bname")
function bname(a:string)
;;
;; Arguments
;;    a: (array of) file name strings to analyse.
;;
;; Return value
;;    An array of strings, without paths and suffixes.
;;
;; Description
;;    Strips paths and suffixes from file names.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20131106-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, a, i, b, rm_prefix, rm_suffix
begin

    funcname = "bname"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    b = a
    do i = 0, dimsizes(a) - 1
        rm_prefix = systemfunc("basename " + a(i))
        rm_suffix = str_split(rm_prefix, ".")
        b(i) = rm_suffix(0)
        delete(rm_suffix)
    end do

    leave_msg(scriptname, funcname, 10)
    return(b)

end

;;#############################################################################
undef("PRIORITY")
function PRIORITY(a)
;;
;; Arguments
;;    a: array to analyse, could be numeric or string.
;;
;; Return value
;;    array of unique elements, ordered by commonness in "a", starting with the
;;    most common element, equally common elements are ordered alphabetically.
;;
;; Description
;;    Useful for assigning priorities to array elements.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130419-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, a, valid_types, u, idx, i, n, nu
begin

    funcname = "PRIORITY"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    valid_types = (/"integer", "float", "double", "string"/)
    if (.not. any(typeof(a) .eq. valid_types)) then
        error_msg("w", scriptname, funcname, "input is not the " + \
                  "correct type")
        leave_msg(scriptname, funcname, 10)
        return(new(1, float))    ; Return missing
    else
        if (typeof(a) .eq. "string") then  ; String
            sqsort(a)
        else  ; Numeric
            qsort(a)
        end if
        idx = UNIQ(a)
        ;; Unique elements in alphabetical/numerical order
        u = a(idx)
        ;; Shall collect number of occurences of each element
        nu = new(dimsizes(u), integer)
        delete(idx)
        do i = 0, dimsizes(u) - 1
            n = dimsizes(ind(a .eq. u(i)))
            nu(i) = n
        end do
        ;; -1 = return permutation resulting from sorting nu in
        ;; Decreasing order, do not sort nu
        idx = dim_pqsort(nu, -1)
        u = u(idx)  ; Sort u accordingly
    end if

    leave_msg(scriptname, funcname, 10)
    return(u)

end

;;#############################################################################
undef("get_output_dir")
function get_output_dir(plot_dir[1]:string,
                        subfolder[1]:string)
;;
;; Arguments
;;    plot_dir: root path for all plots.
;;    subfolder: subfolder for diag_script plots.
;;
;; Return value
;;    A string with the full path of the output directory for plots.
;;
;; Description
;;    Creates and return the path of the output directory for plots.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130419-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, output_dir
begin

    funcname = "get_output_dir"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    output_dir = plot_dir + "/" + subfolder + "/"
    if (.not. isfilepresent_esmval(output_dir)) then
        system("mkdir -p " + output_dir)
    end if

    leave_msg(scriptname, funcname, 10)
    return(output_dir)

end

;;#############################################################################
undef("basename")
function basename(name[1]:string)
;;
;; Arguments
;;    name: input string.
;;
;; Return value
;;    A string.
;;
;; Description
;;    Removes the suffix from the string 'name' and returns the string.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130419-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, suffix
begin

    funcname = "basename"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    suffix = get_file_suffix(name, 0)
    bsn = str_sub_str(name, suffix, "")

    leave_msg(scriptname, funcname, 10)
    return(bsn)

end

;;#############################################################################
undef("extract_years")
function extract_years(data,
                       startYr,
                       endYr)
;;
;; Arguments
;;
;; Return value
;;
;; Description
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, data_new, tim, t_utc, t_utc_yr, idx_t
begin

    funcname = "extract_years"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    ds = dimsizes(data)
    rank = dimsizes(ds)

    tim = data&time

    t_utc = cd_calendar(tim, 0)
    t_utc_yr = floattoint(t_utc(:, 0))
    idx_t = ind((t_utc_yr .ge. startYr) .and. (t_utc_yr .le. endYr))

    if (rank .eq. 2) then
        data_new = data(idx_t, :)
    end if
    if (rank .eq. 3) then
        data_new = data(idx_t, :, :)
    end if
    if (rank .eq. 4) then
        data_new = data(idx_t, :, :, :)
    end if

    leave_msg(scriptname, funcname, 10)
    return(data_new)

end

;;#############################################################################
undef("extend_var_at")
function extend_var_at(var,
                       coord:integer,
                       index:integer)
;;
;; Arguments
;;    var : array of variables to be extended .
;;    coord : index of coordinate to be extended.
;;    index : index where additional entry is created.
;;
;; Return value
;;    An array with an extra dimension.
;;
;; Description
;;    Adds to var in coord at index a new empty entry by retaining metadata.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130419-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, var, coord, index, var_new, var_dim, \
    var_dim_plus, rank
begin

    funcname = "extend_var_at"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    var_dim = dimsizes(var)
    var_dim_plus = var_dim
    var_dim_plus(coord) = var_dim(coord) + 1

    ;; Define new variables  with one more entry for diagnostics
    rank = dimsizes(var_dim)
    var_temp = new((/var_dim_plus/), typeof(var))

    ;; Name cooridnates
    var_new = copy_CoordNames_n(var, var_temp, coord)

    ;; Copy old entries
    if (rank.eq.1) then
        if (index - 1.ge.0) then
            var_new(:index - 1) = var(:index - 1)
        end if
        if (index.lt.var_dim(0)) then
            var_new(index + 1:) = var(index:)
        end if
    else if (rank.eq.2) then
        if (coord.eq.0) then
            if (index - 1.ge.0) then
                var_new(:index - 1, :) = var(:index - 1, :)
            end if
            if (index.lt.var_dim(0)) then
                var_new(index + 1:, :) = var(index:, :)
            end if
        else
            if (index - 1.ge.0) then
                var_new(:, :index - 1) = var(:, :index - 1)
            end if
            if (index.lt.var_dim(1)) then
                var_new(:, index + 1:) = var(:, index:)
            end if
        end if
    else if (rank.eq.3) then
        if (coord.eq.0) then
            if (index - 1.ge.0) then
                var_new(:index - 1, :, :) = var(:index - 1, :, :)
            end if
            if (index.lt.var_dim(0)) then
                var_new(index + 1:, :, :) = var(index:, :, :)
            end if
        else if (coord.eq.1) then
            if (index - 1.ge.0) then
                var_new(:, :index - 1, :) = var(:, :index - 1, :)
            end if
            if (index.lt.var_dim(1)) then
                var_new(:, index + 1:, :) = var(:, index:, :)
            end if
        else
            if (index - 1.ge.0) then
                var_new(:, :, :index - 1) = var(:, :, :index - 1)
            end if
            if (index.lt.var_dim(2)) then
                var_new(:, :, index + 1:) = var(:, :, index:)
            end if
        end if
        end if
    else
        error_msg("f", scriptname, funcname, "implemented only for " + \
                  "variables with dimension less or equal 3!")
    end if
    end if
    end if

    leave_msg(scriptname, funcname, 10)
    return(var_new)

end

;;#############################################################################
undef("copy_CoordNames_n")
function copy_CoordNames_n(var_from,
                           var_to,
                           n)
;;
;; Arguments
;;    var_from: an array of any dimensionality. Must have named dimensions and
;;              coordinate variables.
;;    var_to: an array of any dimensionality.
;;    n: index of coordinate which is not to be copied.
;;
;; Return value
;;    An array with the same dimensionality of var_from.
;;
;; Description
;:    Copies all named dimensions and corrdinate variables except the n-th one.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20130419-A_gott_kl: written.
;;
local funcname, scriptname, verbosity, var_from, var_to, n, rank
begin

    funcname = "copy_CoordNames_n"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    rank = dimsizes(dimsizes(var_from))

    if (rank.le.n) then
        error_msg("f", scriptname, funcname, "index is not in the " + \
                  "range of the variables")
    else
        names = getvardims(var_from)
        do ii = 0, rank - 1
            if (.not.ismissing(names(ii))) then
                var_to!ii = names(ii)
                if (ii.ne.n) then
                    var_to&$names(ii)$ = var_from&$names(ii)$
                end if
            end if
        end do
    end if

    leave_msg(scriptname, funcname, 10)
    return(var_to)

end

;;#############################################################################
undef("empty_str")
function empty_str(nn[1]:integer)
;;
;; Arguments
;;     nn: length of the returned string.
;;
;; Return value
;;     An empty string of length nn.
;;
;; Description
;;     Returns an empty string of length nn (useful for text formatting).
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;     20140228-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, ii
begin

    funcname = "empty_str"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 10)

    str = " "
    do ii = 1, nn - 1
        str = str + " "
    end do

    leave_msg(scriptname, funcname, 10)
    return(str)

end

;;#############################################################################
undef("write_info")
function write_info(display[1]:logical)
;;
;; Arguments
;;    display: logical with attributes listing the items to put in the infobox
;;
;; Return value
;;    A string to write to plot.
;;
;; Description
;;    Read the attributes in displayinfo and returns them in a formatted
;;    string.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, infobox, info_id
begin

    funcname = "write_info"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 8)

    infobox = ""
    infobox = infobox + "~F22~DEBUG INFO:~F21~ ~C~"
    display_varatts = getvaratts(display)
    exceptions = (/"description", "description_ycoord"/)
    do info_id = 0, dimsizes(display_varatts) - 1
        info_tag = display_varatts(info_id)
        if (any(info_tag .eq. exceptions)) then
            continue
        end if
        ;; ~C~ = new line
        infobox = infobox + info_tag + ": " + display@$info_tag$ + "~C~"
    end do

    leave_msg(scriptname, funcname, 8)
    return(infobox)

end

;;#############################################################################
undef("remove_index")
function remove_index(array,
                      idx:integer)
;;
;; Arguments
;;    array: an array of any type
;;    index: the index(es) to be removed.
;;
;; Description
;;    Remove the selected index(es) from the array's dimension 0.
;;
;; Caveats
;;    Only ranking upto 5 is supported.
;;
;; References
;;
;; Modification history
;;    20160704-A_righ_ma: written.
;;
local funcname, scriptname, verbosity, dsize, rank, iarr
begin

    funcname = "remove_index"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    dsize = dimsizes(array)
    rank = dimsizes(dsize)
    dsize := dsize(0)
    iarr = ispan(0, dsize - 1, 1)
    iarr@_FillValue = default_fillvalue(typeof(iarr))

    iarr(idx) = iarr@_FillValue
    iarr := iarr(ind(.not.ismissing(iarr)))

    if (rank.eq.1) then
        out = array(iarr)
    end if
    if (rank.eq.2) then
        out = array(iarr, :)
    end if
    if (rank.eq.3) then
        out = array(iarr, :, :)
    end if
    if (rank.eq.4) then
        out = array(iarr, :, :, :)
    end if
    if (rank.eq.5) then
        out = array(iarr, :, :, :, :)
    end if
    if (rank.gt.5) then
        error_msg("f", funcname, scriptname, \
                  "rank greater than 5 not supported")
    end if

    leave_msg(scriptname, funcname, 6)
    return(out)

end

;;#############################################################################
undef("tool_stop")
procedure tool_stop()
;;
;; Arguments
;;
;; Description
;;    Forces the tool execution to stop at any point of an ncl script.
;;    For dubugging purposes only.
;;
;; Caveats
;;    This should not be used for the standard exit-on-error. See the function
;;    error_msg in messaging.ncl instead.
;;
;; References
;;
;; Modification history
;;    20150521-A_righ_ma: written.
;;
begin

    print("fatal: tool stop forced by the user")
    status_exit(1)

end

;;#############################################################################
undef("filter_attrs")
function filter_attrs(source[1]:logical,\
                      target[1]:logical,\
                      lb_prefix[1]:string)
;;
;; Arguments
;;
;; Description
;;    Copy over certain attributes if they have a specified prefix,
;;    the prefix is removed on the target attribute
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150630-A_eval_ma: written.
;;
local source_attrs, att_idx, prefix_removed
begin

    funcname = "filter_attrs"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    source_attrs = getvaratts(source)
    do att_idx = 0, dimsizes(source_attrs) - 1
        if (.not. ismissing(str_match(source_attrs(att_idx), lb_prefix)))
            prefix_removed = str_get_cols(source_attrs(att_idx), \
                                          strlen(lb_prefix), -1)
            target@$prefix_removed$ = source@$source_attrs(att_idx)$
        end if
    end do

    leave_msg(scriptname, funcname, 6)
    return(target)

end

;;#############################################################################
undef("write_ignore_warnings")
procedure write_ignore_warnings(warnings[*]:string,
                                warn_file_str[1]:string)
;;
;; Arguments
;;
;; Description
;;    Write a file with warnings that will be disregarded by the Pyhton
;;    launcher.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150630-A_eval_ma: written.
;;
local funcname, scriptname, verbosity, idx_warn
begin

    funcname = "write_ignore_warnings"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    if (isfilepresent_esmval(warn_file_str)) then
        system("rm -f " + warn_file_str)
    end if
    do idx_warn = 0, dimsizes(warnings) - 1
        system("echo " + warnings(idx_warn) + " >> " + warn_file_str)
    end do

    leave_msg(scriptname, funcname, 6)

end

;;#############################################################################
undef("get_single_ref")
function get_single_ref(models[1]:logical,\
                        model_refs[1]:string)
;;
;; Arguments
;;    models: list of all models
;;    model_refs: reference model name
;;
;; Description
;;    Fetches the reference model index in models@name
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151109-A_righ_ma: moved from diag_scripts/aux/SAMonsson.ncl to here.
;;    2015????-A_eval_ma: written.
;;
local idx_single_ref
begin
    funcname = "get_single_ref"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    idx_single_ref = ind(model_attr_id .eq. model_refs)
    if (all(ismissing(idx_single_ref))) then
        delete(idx_single_ref)
        idx_single_ref = ind(models@name .eq. str_squeeze(model_refs))
    end if

    ;; Check if ref_model is in the current namelist
    if (all(ismissing(idx_single_ref))) then
        error_msg("f", scriptname, funcname, "missing reference model - " + \
                  "ref_model is not included in model list (namelist)")
    end if

    ;; In case we're comparing the same model (i.e., same models@name)
    idx_single_ref_sizes = dimsizes(idx_single_ref)
    if (idx_single_ref_sizes.gt.1) then
        idx_single_ref_tmp = idx_single_ref(0)
        delete(idx_single_ref)
        idx_single_ref = idx_single_ref_tmp
        delete(idx_single_ref_tmp)
    end if

    leave_msg(scriptname, funcname, 6)
    return idx_single_ref
end

;;#############################################################################
undef("get_ref_model_idx")
function get_ref_model_idx(models[1]:logical,\
                           model_ref_name[1]:string)
;;
;; Arguments
;;    models: list of all models
;;    model_ref_name: reference model name, possibly two, separated by comma
;;
;; Description
;;    Fetches the reference model index in models@name
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151109-A_righ_ma: moved from diag_scripts/aux/SAMonsson.ncl to here.
;;    2015????-A_eval_ma: written.
;;
local funcname, scriptname, verbosity, idx_ref, idx_ref_sizes, idx_ref_tmp, idx, idx_ref_size
begin

    funcname = "get_ref_model_idx"
    scriptname = "interface_scripts/auxiliary.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    model_refs = str_split(model_ref_name, ",")

    if (dimsizes(model_refs) .eq. 1) then
         idx_ref = get_single_ref(models, model_refs)
    else if (dimsizes(model_refs) .gt. 1) then
        idx_ref_size = dimsizes(model_refs)
        idx_ref = new(idx_ref_size, integer)
        do idx = 0, idx_ref_size - 1
            idx_ref(idx) = get_single_ref(models, model_refs(idx))
        end do
    end if
    end if

    leave_msg(scriptname, funcname, 6)
    return(idx_ref)
end
