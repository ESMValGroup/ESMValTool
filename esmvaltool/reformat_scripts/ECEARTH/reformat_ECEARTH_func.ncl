;;#############################################################################
;; FUNCTIONS FOR THE ECEARTH REFORMAT TOOL
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function find_name
;;    procedure add_ijcoordvars
;;    procedure add_ECEARTHgrid_info
;;    procedure convert_units
;;    procedure mask_land
;;    procedure sub_staggergrid
;;    procedure rename_lev
;;    procedure rename_time
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/messaging.ncl"

;;#############################################################################
undef("find_name")
function find_name(name:string, which:string)
;;
;; Arguments
;;    name: the variable's standard name.
;;    which: std: find the standard name given the EC-Earth name
;;         : ECEARTH (or anything else): find the EC-Earth name given the
;;           standard name
;;
;; Return value
;;    A string indicating the ECEARTH name of the requested variable.
;;
;; Description
;;    Finds the corresponding ECEARTH name of the input variable name.
;;    Needs an external input table "names_ECEARTH.dat".
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local col_in, col_out, data, ecearth, fname, funcname, grid, id, scriptname, \
    std, verbosity
begin

    funcname = "find_name"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; read text file with names
    fname = "./reformat_scripts/" + PROJECT + "/names_" + PROJECT + \
            ".dat"
    data = readAsciiTable(fname, 1, "string", \
           dimsizes(readAsciiHead(fname, "STANDARD")) + 1)
    std     = str_squeeze(str_get_field(data(:, 0), 1, "|"))
    ecearth = str_squeeze(str_get_field(data(:, 0), 2, "|"))
    grid    = str_squeeze(str_get_field(data(:, 0), 3, "|"))

    ;; std -> EC-Earth of vice versa
    if (which .eq. "std") then
        col_in = ecearth
        col_out    = std
    else
        col_in = std
        col_out    = ecearth
    end if

    ;; find our variable and store info
    id = ind(col_in .eq. name)
    if (dimsizes(id) .gt. 1) then
        error_msg("f", scriptname, funcname, \
                  "more than 1 entry matching name " + name + \
                  "found in " + fname + ". Can't decide which one " + \
                  "to use.")
    end if
    if (ismissing(id)) then
        error_msg("f", scriptname, funcname, \
                  "cannot find " + which + " name for " + name + ".")
    end if
    out = col_out(id)
    out@staggergrid = grid(id)

    leave_msg(scriptname, funcname, 4)

    return(out)

end

;#######################################################################
undef("add_ijcoordvars")
procedure add_ijcoordvars(var:numeric)
;;
;; Arguments
;;    var: the variable's field.
;;
;; Description
;;    Rename x and y dimensions to i and j and add i and j variables and
;;    assign them as coordinate variables.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local dim_sizes, funcname, i, ix, iy, j, scriptname, verbosity
begin

    funcname = "add_ijcoordvars"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; get dimension indices
    dim_sizes = dimsizes(var)
    ix = ind(getvardims(var).eq."x")
    iy = ind(getvardims(var).eq."y")
    if (ismissing(ix)) then
        error_msg("f", scriptname, funcname, \
                  "x dimension not found.")
    end if
    if (ismissing(iy)) then
        error_msg("f", scriptname, funcname, \
                  "y dimension not found.")
    end if

    ;; rename dimensions to i and j
    var!ix = "i"
    var!iy = "j"

    ;; create i and j variables
    i = ispan(1, dim_sizes(ix), 1)
    i@units = "1"
    i@long_name = "cell index along i dimension"
    j = ispan(1, dim_sizes(iy), 1)
    j@units = "1"
    j@long_name = "cell index along j dimension"

    ;; append i and j variables as coordinate variables
    var&i = i
    var&j = j

    leave_msg(scriptname, funcname, 4)

end

;#######################################################################
undef("add_ECEARTHgrid_info")
procedure add_ECEARTHgrid_info(var:numeric, name:string)
;;
;; Arguments
;;    var        : the variable's field.
;;    name       : the variable's ECEARTH name.
;;
;; Description
;;    Add ECEARTH grid info (lat, lon, areacello and grid sizes) to the
;;    output
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local area, ce1, ce2, ce3, clat, clon, e1_name, e2_name, e3_name, ff, \
    files, funcname, i, j, lat, lat_name, ldims, lgrid, lgrids, lon, \
    lon_name, lsize, nl, scriptname, tmp, vars, vdims, verbosity, vrank, \
    vsize, w
begin

    funcname = "add_ECEARTHgrid_info"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    nl = str_get_nl()

    ;; Variable dimensions (must be consistent with lat and lon dims)
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
    lat_name = find_name("lat", "ECEARTH")
    if (any(lat_name .eq. vars)) then
        clat = ff->$lat_name$
    else
        error_msg("f", scriptname, funcname, \
                  "cannot find latitude variable for the EC-Earth " + \
                  "grid.")
    end if
    lat = lat_attrib(clat)
    ldims = getvardims(clat)
    do i = 0, dimsizes(ldims) - 1, 1
        ldims(i) = find_name(ldims(i), "std")
    end do
    lsize = dimsizes(clat)
    if (ldims(1) .ne. vdims(vrank - 1) .or. \
        ldims(0) .ne. vdims(vrank - 2) .or. \
        lsize(1) .ne. vsize(vrank - 1) .or. \
        lsize(0) .ne. vsize(vrank - 2)) then
        error_msg("f", scriptname, funcname, \
                  "inconsistent latitude coordinate for the " + \
                  "EC-Earth grid.")
    end if
    filevardef(w, "lat", typeof(lat), ldims)
    filevarattdef(w, "lat", lat)
    delete(ldims)
    delete(lsize)

    ;; Read, reformat and check longitude coordinate
    lon_name = find_name("lon", "ECEARTH")
    if (any(lon_name .eq. vars)) then
        clon = ff->$lon_name$
    else
        error_msg("f", scriptname, funcname, \
                  "cannot find longitude variable for the EC-Earth " + \
                  "grid.")
    end if
    lon = lon_attrib(clon)
    ldims = getvardims(clon)
    do i = 0, dimsizes(ldims) - 1, 1
        ldims(i) = find_name(ldims(i), "std")
    end do
    lsize = dimsizes(clon)
    if (ldims(1) .ne. vdims(vrank - 1) .or. \
        ldims(0) .ne. vdims(vrank - 2) .or. \
        lsize(1) .ne. vsize(vrank - 1) .or. \
        lsize(0) .ne. vsize(vrank - 2)) then
        error_msg("f", scriptname, funcname, \
                  "inconsistent longitude coordinate for the " + \
                  "EC-Earth grid.")
    end if
    filevardef(w, "lon", typeof(lon), ldims)
    filevarattdef(w, "lon", lon)
    delete(ldims)
    delete(lsize)
    delete(vars)

    ;; Add detailed horizontal grid sizes (for grid on which variable is
    ;; known) and areacello (for T grid) to the output file
    if (.not. isfilepresent_esmval(HGRIDFILE)) then
        error_msg("f", scriptname, funcname, \
            nl + nl + \
            "cannot find the specified horizontal grid file " + \
            HGRIDFILE + nl + nl + \
            "********************************************************" + nl + \
            "* NOTE:                                                *" + nl + \
            "********************************************************" + nl + \
            "* The location of the hgrid file should be given       *" + nl + \
            "* inside the <AUXILIARIES>-tag by an 'fx_file'-tag     *" + nl + \
            "* entry with attribute:                                *" + nl + \
            "*      id=" + dq + "nemo_hgrid_file" + dq + "          *" + nl + \
            "* To fix this error, update your 'fx_file'-tag         *" + nl + \
            "********************************************************")
    else
        if (str_lower(name@staggergrid) .eq. "t") then
            ;; only T grid (detailed info and areacello)
            lgrids = (/ str_lower(name@staggergrid) /)
        else
            ;; T (areacello) and variable grid (detailed info), note
            ;; that order is important, T grid must be treated first
            lgrids = (/ "t", str_lower(name@staggergrid) /)
        end if

        ff = addfile(HGRIDFILE, "r")
        vars = getfilevarnames(ff)

        do j = 0, dimsizes(lgrids) - 1
            ;; get e1* and e2*
            lgrid = lgrids(j)
            e1_name = "e1" + lgrid
            e2_name = "e2" + lgrid
            if (any(e1_name .eq. vars)) then
                ce1 = ff->$e1_name$
            else
                error_msg("f", scriptname, funcname, \
                          "cannot find " + e1_name + " variable " + \
                          "for the EC-Earth grid in " + HGRIDFILE + ".")
            end if
            if (any(e2_name .eq. vars)) then
                ce2 = ff->$e2_name$
            else
                error_msg("f", scriptname, funcname, \
                          "cannot find " + e2_name + " variable " + \
                          "for the EC-Earth grid in " + HGRIDFILE + ".")
            end if

            ;; require 1 time step for variables
            if (ce1!0 .ne. "t" .or. dimsizes(ce1!0) .ne. 1) then
                error_msg("f", scriptname, funcname, \
                          "first dimension of " + e1_name + " must " + \
                          "be t and have length 1.")
            end if
            if (ce2!0 .ne. "t" .or. dimsizes(ce2!0) .ne. 1) then
                error_msg("f", scriptname, funcname, \
                          "first dimension of " + e2_name + " must " + \
                          "be t and have length 1.")
            end if

            ;; for T grid, calculate area, assume in metres and set dims
            ;; and atts
            if (lgrid .eq. "t") then
                area = ce1(0, :, :) * ce2(0, :, :)
                area!0 = ce1!1
                area!1 = ce1!2
                area@units = "m2"
                check_units(area, "areacello")
                tmp = var_attrib(area, "areacello")
                delete(area)
                area = tmp
                delete(tmp)
                ldims = getvardims(area)
                do i = 0, dimsizes(ldims) - 1, 1
                    ldims(i) = find_name(ldims(i), "std")
                end do
                lsize = dimsizes(area)

                ;; Check consistency of dimension sizes
                if (lsize(1) .ne. vsize(vrank - 1) .or. \
                    lsize(0) .ne. vsize(vrank - 2)) then
                    error_msg("f", scriptname, funcname, \
                              "inconsistent areacello for the " + \
                              "irregular grid.")
                end if

                ;; Check consistency of dimension names
                if (ldims(1) .ne. vdims(vrank - 1) .or. \
                    ldims(0) .ne. vdims(vrank - 2)) then
                    area!1 = vdims(vrank - 1)
                    area!0 = vdims(vrank - 2)
                    ldims = getvardims(area)
                end if

                filevardef(w, "areacello", typeof(area), ldims)
                filevarattdef(w, "areacello", area)
            end if
        end do

        ;; add e1* and e2* on grid on which variable is known
        filevardef(w, e1_name, typeof(ce1), ldims)
        filevardef(w, e2_name, typeof(ce2), ldims)

        ;; clean a bit
        delete(ldims)
        delete(lsize)
        delete(vars)
        delete(ff)
    end if

    ;; For 3D variables, add detailed vertical grid size (for grid on
    ;; which variable is known)
    if (isStrSubset(FIELD, "T3")) then
        if (.not. isfilepresent_esmval(ZGRIDFILE)) then
            error_msg("f", scriptname, funcname, \
                nl + nl + \
                "cannot find the specified vertical grid file " + \
                ZGRIDFILE + nl + nl + \
                "********************************************************" + nl + \
                "* NOTE:                                                *" + nl + \
                "********************************************************" + nl + \
                "* The location of the zgrid file should be given       *" + nl + \
                "* inside the <AUXILIARIES>-tag by an 'fx_file'-tag     *" + nl + \
                "* entry with attribute:                                *" + nl + \
                "*      id=" + dq + "nemo_zgrid_file" + dq + "          *" + nl + \
                "* To fix this error, update your 'fx_file'-tag         *" + nl + \
                "********************************************************")
        else
            lgrid = str_lower(name@staggergrid)

            ff = addfile(ZGRIDFILE, "r")
            vars = getfilevarnames(ff)

            ;; get e3*
            e3_name = "e3" + lgrid
            if (any(e3_name .eq. vars)) then
                ce3 = ff->$e3_name$
            else
                error_msg("f", scriptname, funcname, \
                          "cannot find " + e3_name + " variable " + \
                          "for the EC-Earth grid in " + ZGRIDFILE + ".")
            end if

            ;; require 1 time step for variable, then collapse it
            if (ce3!0 .ne. "t" .or. dimsizes(ce3!0) .ne. 1) then
                error_msg("f", scriptname, funcname, \
                          "first dimension of " + e3_name + " must " + \
                          "be t and have length 1.")
            end if
            tmp = ce3(0, :, :, :)
            delete(ce3)
            ce3 = tmp
            delete(tmp)

            ;; Check consistency of remaining dimension sizes and names
            ldims = getvardims(ce3)
            do i = 0, dimsizes(ldims) - 1, 1
                ldims(i) = find_name(ldims(i), "std")
            end do
            lsize = dimsizes(ce3)
            if (dimsizes(ldims) .ne. vrank - 1) then
                error_msg("f", scriptname, funcname, \
                          "variable " + name + \
                          " and vertical coordinate " + e3_name + \
                          " have an inconsistent number of " + \
                          "dimensions:" + nl + \
                          dimsizes(ldims) + " <-> " + vrank)
            end if
            do i = 0, dimsizes(lsize) - 1
                if (ldims(i) .ne. vdims(i + 1)) then
                    error_msg("f", scriptname, funcname, \
                              "variable " + name + \
                              " and vertical coordinate " + e3_name + \
                              " have different dimension names. " + \
                              nl + str_join(ldims, " | ") + \
                              nl + str_join(vdims, " | "))
                end if
                if (lsize(i) .ne. vsize(i + 1)) then
                    error_msg("f", scriptname, funcname, \
                              "variable " + name + \
                              " and vertical coordinate " + e3_name + \
                              " have different dimension sizes. " + \
                              nl + str_join(tostring(lsize), " | ") + \
                              nl + str_join(tostring(vsize), " | "))
                end if
            end do

            ;; add e3* on grid on which variable is known
            filevardef(w, e3_name, typeof(ce3), ldims)

            ;; clean a bit
            delete(ldims)
            delete(lsize)
            delete(vars)
            delete(ff)
        end if
    end if

    ;; Write output
    setfileoption(w, "DefineMode", False)
    w->lat = (/lat/)
    w->lon = (/lon/)
    w->areacello = (/area/)
    w->$e1_name$ = (/ce1/)
    w->$e2_name$ = (/ce2/)
    if (isdefined("ce3")) then
        w->$e3_name$ = (/ce3/)
    end if
    delete(w)

    leave_msg(scriptname, funcname, 4)

end

;#######################################################################
undef("convert_units")
procedure convert_units(var:numeric, name:string)
;;
;; Arguments
;;    var : the variable's field.
;;    name: the variable's standard name.
;;
;; Description
;;    Unit conversions that cannot be handled by check_units.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local atts, cmor_units, funcname, scriptname, verbosity
begin

    funcname = "convert_units"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; get CMOR units
    atts = read_cmor(name)
    if (isatt(atts, "units")) then
        cmor_units = atts@units
    else
        error_msg("f", scriptname, funcname, \
                  "cannot find units attribute in the CMOR table " + \
                  "for variable " + name + ".")
    end if
    delete(atts)

    ;: Degrees C to K
    if (var@units .eq. "C" .and. cmor_units .eq. "K") then
        var       = var + 273.15
        var@units = "K"
    end if

    ;; Go from Evaporation - precipitation to P+E
    if (name .eq. "wfpe") then
        var = -1 * var
    end if

    leave_msg(scriptname, funcname, 4)

end

;#######################################################################
undef("mask_land")
procedure mask_land(var:numeric)
;;
;; Arguments
;;    var: the variable's field.
;;
;; Description
;;    Land points have 0, not a fill value. The routine sets land points
;;    (as in the lsm file) to fill values.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local flsm, funcname, i, lsm, ldims, lrank, lsize, ocean_only, scriptname, \
    vdims, vrank, vsize, verbosity
begin

    funcname = "mask_land"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    nl = str_get_nl()
    dq = str_get_dq()

    ;; check the land mask file exists
    if (.not. isfilepresent_esmval(LSMFILE)) then
        error_msg("f", scriptname, funcname, \
            nl + nl + \
            "cannot find the specified land mask file " + \
            LSMFILE + nl + nl + \
            "********************************************************" + nl + \
            "* NOTE:                                                *" + nl + \
            "********************************************************" + nl + \
            "* The location of the land mask file should be given   *" + nl + \
            "* inside the <AUXILIARIES>-tag by an 'fx_file'-tag     *" + nl + \
            "* entry with attribute:                                *" + nl + \
            "*      id=" + dq + "nemo_lsm_file" + dq + ", and       *" + nl + \
            "*      id=" + dq + "nemo_lsm3d_file" + dq + "          *" + nl + \
            "* To fix this error, update your 'fx_file'-tags        *" + nl + \
            "********************************************************")
    end if

    ;; create local copy
    ocean_only = var

    ;; open landmask file
    flsm = addfile(LSMFILE, "r")
    if (isfilevar(flsm, "lsm3d")) then
        lsm = flsm->lsm3d
    else if (isfilevar(flsm, "lsm")) then
        lsm = flsm->lsm
    else
        error_msg("f", scriptname, funcname, \
                  "cannot find the lsm variable in land mask file " + \
                  LSMFILE + ". It must be named lsm or lsm3d.")
    end if
    end if

    ;; get dimension names and sizes
    ldims = getvardims(lsm)
    lsize = dimsizes(lsm)
    lrank = dimsizes(lsize)
    vdims = getvardims(var)
    vsize = dimsizes(var)
    vrank = dimsizes(vsize)

    ;; Check consistency of dimension sizes
    if (lrank .ne. vrank - 1) then
        error_msg("f", scriptname, funcname, \
                  "land mask has incorrect number of dimensions.")
    end if

    do i = 0, lrank - 1
        if (lsize(i) .ne. vsize(i + 1)) then
            error_msg("f", scriptname, funcname, \
                      "dimension sizes of land mask do not match " + \
                      "those of variable:" + nl + \
                      "dimension " + vdims(i + 1) + ": " + lsize(i) + \
                      " <-> " + vsize(i + 1))
        end if
    end do

    ;; Return the values of var where lsm.eq.0
    ocean_only = mask(var, lsm, 0)

    ;; Copy back
    var = ocean_only

    leave_msg(scriptname, funcname, 4)

end

;#######################################################################
undef("sub_staggergrid")
procedure sub_staggergrid(name:string)
;;
;; Arguments
;;    name: the variable's ECEARTH name.
;;
;; Description
;;    Substitute STAGGERGRID in INFILE with actual grid type.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local funcname, scriptname, verbosity
begin

    funcname = "sub_staggergrid"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; check INFILE format
    if (.not. isStrSubset(INFILE, "_grid_*")) then
        error_msg("f", scriptname, funcname, \
                  "string grid_* not found in INFILE format " + \
                  INFILE + ". Check get_cf_infile in " + \
                  "'interface_scripts/projects.py' for project " + \
                  PROJECT + ".")
    end if

    ;; Replace string
    INFILE = str_sub_str(INFILE, "_grid_*", "_grid_" + name@staggergrid)

    leave_msg(scriptname, funcname, 4)

end

;#######################################################################
undef("rename_lev")
procedure rename_lev(var:numeric, name:string)
;;
;; Arguments
;;    var : the variable's field.
;;    name: the variable's ECEARTH name.
;;
;; Description
;;    Vertical coordinate name depends on grid, rename it to lev and add
;;    standard name depth.
;;    Needs an external input table "names_ECEARTH.dat".
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local ecearth, fname, funcname, grid, i, id, levname, lfound, scriptname, \
    std, verbosity
begin

    funcname = "rename_lev"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; set name for level and file name with names
    levname = "lev"
    fname = "./reformat_scripts/" + PROJECT + "/names_" + PROJECT + \
            ".dat"

    ;; Applies to T3X type
    if (isStrSubset(FIELD, "T3")) then
        ;; read text file with names
        data    = readAsciiTable(fname, 1, "string", \
                  dimsizes(readAsciiHead(fname, "STANDARD")) + 1)
        std     = str_squeeze(str_get_field(data(:, 0), 1, "|"))
        ecearth = str_squeeze(str_get_field(data(:, 0), 2, "|"))
        grid    = str_squeeze(str_get_field(data(:, 0), 3, "|"))

        ;; find vertical coordinate name for the grid of this variable
        id = ind(std .eq. levname .and. grid .eq. name@staggergrid)
        if (dimsizes(id) .gt. 1) then
            error_msg("f", scriptname, funcname, \
                      "more than 1 entry matching name " + levname + \
                      " on grid " + name@staggergrid + " found in " + \
                      fname + ". Can't decide which one to use.")
        end if
        if (ismissing(id)) then
            error_msg("f", scriptname, funcname, "cannot find " + PROJECT + \
                      " name for " + levname + " on the " + \
                      name@staggergrid + " grid.")
        end if

        ;; rename vertical coordinate
        i = 0
        lfound = False
        do i = 0, dimsizes(dimsizes(var)) - 1
            if (var!i .eq. ecearth(id)) then
                var!i = std(id)
                lfound = True
            end if
        end do
        if (.not. lfound) then
            error_msg("f", scriptname, funcname, \
                      "cannot find vertical coordinate " + \
                      ecearth(id) + " for variable " + name + ".")
        end if

        ;; set standard name attribute, required for reformat_lev
        var&$std(id)$@standard_name = "depth"
    end if

    leave_msg(scriptname, funcname, 4)

end

;#######################################################################
undef("rename_time")
procedure rename_time(var:numeric, name:string)
;;
;; Arguments
;;    var : the variable's field.
;;    name: the variable's ECEARTH name.
;;
;; Description
;;    Rename time variable from EC-Earth name to standard name and
;;    remove_FillValue attribute.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20141205-A_vanu_be: written.
;;
local ecearth_name, funcname, i, lfound, scriptname, std_name, tmp, verbosity
begin

    funcname = "reformat_time_ecearth"
    scriptname = "reformat_scripts/ECEARTH/reformat_ECEARTH_func.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; get time name in EC-Earth
    std_name = "time"
    ecearth_name = find_name(std_name, "ECEARTH")

    ;; rename time dimension
    i = 0
    lfound = False
    do i = 0, dimsizes(dimsizes(var)) - 1
        if (var!i .eq. ecearth_name) then
            tmp = var&$ecearth_name$
            var!i = std_name
            var&$std_name$ = tmp
            delete(tmp)
            lfound = True
        end if
    end do
    if (.not. lfound) then
        error_msg("f", scriptname, funcname, \
                  "cannot find time dimension " + ecearth_name + \
                  " for variable " + name + ".")
    end if

    ;; delete time FillValue as it seems to confuse cd_calendar like
    ;; functions
    if (isatt(var&$std_name$, "_FillValue")) then
        delete(var&$std_name$@_FillValue)
    end if

    leave_msg(scriptname, funcname, 4)

end
