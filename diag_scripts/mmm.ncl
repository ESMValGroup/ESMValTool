;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME
;;#############################################################################
;; Description
;;    Regrids all models to specified grid (coarsest, finest, ref, 1x1) and
;;    then calculates the multi-model mean over all models. The results are
;;    saved to netCDF.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref, 1x1)
;;    ref_model:       name of reference data set (only needed if grid="ref")
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    stddev:          calculate inter-model standard deviation (True, False)
;;    lsmask:          0 (default): do nothing
;;                     1: mask land grid cells in final mmm result
;;                     2: mask ocean grid cells in final mmm result
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;    - missing values are ignored when averaging
;;    - no consistency checks (e.g. number of time steps, etc)
;;    - no handling of irregular grids
;;
;; Modification history
;;    20161028-A_laue_ax: added option to calculate inter-model standard deviation
;;    20161014-A_laue_ax: written.
;;
;;#############################################################################

; #####################################
; # load external NCL code, if needed #
; #####################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"


begin
    ; ###########################################################
    ; # Fetch general parameters, set in namelist_cloud.xml,    #
    ; # passed via environment variables by python code         #
    ; ###########################################################

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)

    ; 'diag_script' and 'variables' are fetched from the
    ; above 'ncl.interface' file

    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ; ###########################################
    ; # references                              #
    ; ###########################################

    write_references(diag_script,      \  ; script name
                     "A_laue_ax",      \  ; authors
                     "",               \  ; contributors
                     "", \  ; diag_references
                     "",               \  ; obs_references
                     (/""/))     ; proj_references

    ; 'models', 'variables' and 'field_types' are fetched from the
    ; above 'ncl.interface' file

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ; check for reference model definition
    if (.not.(isatt(diag_script_info, "ref_model"))) then
        error_msg("f", diag_script, "", "no reference model is selected")
    end if

    ; set regridding method

    if (.not.(isatt(diag_script_info, "grid"))) then
        error_msg("w", diag_script, "", \
                  "Regridding method not defined. " + "Using 1x1")
        selected_grid = "1x1"
    else
        selected_grid = diag_script_info@grid
    end if

    ; define the target grid for regridding

    models_ind = ispan(0, dim_MOD - 1, 1)
    info_output("Setting target grid:", verbosity, 2)
    if (selected_grid.eq."ref") then

        ; set reference model
        ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
        ref_model = diag_script_info@ref_model
        if (ismissing(ref_ind)) then
            if (diag_script_info@ref_model(0).eq."n.a.") then
                ref_ind = 0
            else
                error_msg("f", diag_script, "", \
                     "no adequate reference model provided")
            end if
        end if

        grid = find_destination_grid(ref_ind, \
                                    var0, field_type0, "coarsest")
    else if (selected_grid.eq."1x1") then
        grid     = new((/180, 360/), float)
        grid!0   = "lat"
        grid!1   = "lon"
        grid&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
        grid&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
        grid&lat@standard_name = "latitude"
        grid&lon@standard_name = "longitude"
        grid&lat@long_name     = "latitude"
        grid&lon@long_name     = "longitude"
        grid&lat@units         = "degrees_north"
        grid&lon@units         = "degrees_east"
        grid&lat@axis          = "Y"
        grid&lon@axis          = "X"
    else
        grid = find_destination_grid(models_ind, \
                                     var0, field_type0, selected_grid)
    end if
    end if
    delete(models_ind)

end

; load variable_info

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ; ############
    ; # get data #
    ; ############

    opt = True
    opt@plot_type = "regridded"
    opt@time_avg = selected_grid

    do imod = 0, dim_MOD - 1

        info_output("processing " + models@name(imod), verbosity, 0)

        data = read_data(imod, var0, field_type0)

        ; =====================================================================
        ; =============================== regridding ==========================
        ; =====================================================================

        info_output("Regridding " + models@name(imod), verbosity, 2)
        regr_path = clim_diag_file_path(var0, field_type0, imod, \
                                        diag_script, opt)

        ; Do not process data if already processed

        if (isfilepresent(regr_path) .and. \
            getenv("ESMValTool_force_calc").eq."False") then

            info_output("Reading in file = " + regr_path, verbosity, 1)
            var_hreg = ncdf_read(regr_path, var0)
        else
            A0 = extract_data(imod, data, -1, 0, 0)
            info_output(" ++++++++++++++ Variable summary of original data " \
                        + "(A0): ", verbosity, 5)
            if (verbosity .ge. 4) then
                printVarSummary(A0)
            end if

            ; Determine dimensions to regrid

            xregrid = True  ; lon
            yregrid = True  ; lat
            dnames = getVarDimNames(A0)

            lonfound = False
            latfound = False

            ; find longitude and latitude dimension

            do i = 0, dimsizes(dnames) - 1
                if (str_lower(dnames(i)) .eq. "lon") then
                    lonfound = True
                end if
                if (str_lower(dnames(i)) .eq. "lat") then
                    latfound = True
                end if
            end do

            ; alternatively, try if there are extra variables for lat/lon
            ; in the input file (should be the case for ocean grids)

            if (.not.latfound .or. .not.lonfound) then
                f = addfile(interface_get_fullpath(var0, field_type0, imod), "r")
            end if

            if (.not.lonfound) then
                if (isfilevar(f, "lon")) then
                   srclons = f->lon
                   lonfound = True
                end if
            else
                srclons = A0&lon
            end if

            if (.not.latfound) then
                if (isfilevar(f, "lat")) then
                   srclats = f->lat
                   latfound = True
                end if
            else
                srclats = A0&lat
            end if

            ; if lat/lon are still not found, issue error message

            if (.not.latfound .or. .not.lonfound) then
                error_msg("f", diag_script, "", \
                     "could not find latitude and/or longitude dimensions")
            end if

            ldims = dimsizes(srclons)

            if (dimsizes(ldims).eq.1) then
                if (ldims.eq.dimsizes(grid&lon)) then
                    if (max(abs(srclons - grid&lon)).lt.1.0d-5) then
                        xregrid = False
                    end if
                end if
            end if

            delete(ldims)
            ldims = dimsizes(srclats)

            if (dimsizes(ldims).eq.1) then
                if (ldims.eq.dimsizes(grid&lat)) then
                    if (max(abs(srclats - grid&lat)).lt.1.0d-5) then
                        yregrid = False
                    end if
                end if
            end if

            delete(ldims)

            ; regrid in the horizontal (if necessary)

            if (xregrid.or.yregrid) then
;                var_hreg = linint2_Wrap(srclons, srclats, A0, \
;                                        True, grid&lon, grid&lat, 0)
                Opt = True
                ;Opt@SrcMask2D  = where(ismissing(A0), 0, 1)
                Opt@SrcGridLat = srclats
                Opt@SrcGridLon = srclons
                Opt@DstGridLat = grid&lat
                Opt@DstGridLon = grid&lon
                Opt@ForceOverwrite = True
                Opt@PrintTimings   = True
                Opt@Debug          = True
                Opt@InterpMethod   = "patch"  ; "bilinear", "patch", "conserve"
                ; save regrid files for later use: careful with naming!
                ;Opt@WgtFileName = "ESMF_regrid_Wgt.nc"
                ;Opt@SrcFileName = "ESMF_regrid_Src.nc"
                ;Opt@DstFileName = "ESMF_regrid_Dst.nc"

                var_hreg = ESMF_regrid(A0, Opt)

                delete(Opt)
;            else if (xregrid.and..not.yregrid) then
;                dx = ind(dnames.eq."lon")
;                var_hreg = linint1_n_Wrap(srclons, A0, True, \
;                                          grid&lon, 0, dx)
;                info_output("Latitude regridding not required", \
;                            verbosity, 4)
;            else if (.not.xregrid.and.yregrid) then
;                dx = ind(dnames.eq."lat")
;                var_hreg = linint1_n_Wrap(srclats, A0, False, \
;                                          grid&lat, 0, dx)
;                info_output("Longitude regridding not required", \
;                            verbosity, 4)
;            else if (.not.xregrid .and. .not.yregrid) then
            else
                var_hreg = A0
                info_output("Horizontal regridding not required", \
                            verbosity, 4)
;            end if
;            end if
;            end if
            end if

            delete(A0)

            ; output to NetCDF

            var_hreg@var = var0
            var_hreg@diag_script = (/diag_script/)
            var_hreg@ncdf = regr_path
            ncdf_outfile = ncdf_write(var_hreg, regr_path)

            ; clean up

            if (isdefined("f")) then
                delete(f)
            end if

            delete(srclons)
            delete(srclats)

        end if

        if (.not.isvar("mmm")) then
            mmm = new(dimsizes(var_hreg), float)
            mmm = 0.0
            copy_VarMeta(var_hreg, mmm)

            modelstr = models@name(0)
            do i = 1, dim_MOD - 1
                modelstr = modelstr + "," + models@name(i)
            end do
            mmm@models = tochar(modelstr)
            mmm@experiment = models@experiment(0)

            timeval = toint(cd_calendar(mmm&time, 0))

            t1 = sprinti("%0.4i", timeval(0, 0))
            t2 = sprinti("%0.2i", timeval(0, 1))
            t3 = sprinti("%0.4i", timeval(dimsizes(mmm&time) - 1, 0))
            t4 = sprinti("%0.2i", timeval(dimsizes(mmm&time) - 1, 1))
            timestring = t1 + t2 + "-" + t3 + t4

            mmm@period = tochar(timestring)

            pts = new(dimsizes(var_hreg), integer)
            pts = 0

            delete(timeval)
        end if

        mmm = where(ismissing(var_hreg), mmm, mmm + var_hreg)
        pts = where(ismissing(var_hreg), pts, pts + 1)
    end do  ; imod-loop

    pts = where(pts.eq.0, pts@_FillValue, pts)
    mmm = mmm / pts

    mmm@var = var0
    mmm@diag_script = (/diag_script/)
    mmm@ncdf = regr_path

    ; *************************************************************************
    ; optionally mask results (only for 2-d fields)
    ; *************************************************************************

    lsmask = 0  ; default = no masking

    if (isatt(diag_script_info, "lsmask")) then
        lsmask = diag_script_info@lsmask
    end if

    if (.not.isStrSubset(field_type0, "T2") .and. (lsmask.ne.0)) then
        info_output("Masking only supported for 2-dim fields, disabling masking.", \
                            verbosity, 0)
        lsmask = 0
    end if

    if (lsmask.eq.1 .or. lsmask.eq.2) then
        f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
        a = f->LSMASK
        sftlf = byte2flt(landsea_mask(a, mmm&lat, mmm&lon))
        if (lsmask.eq.1) then  ; mask land grid cells
            sftlf = where(sftlf.ge.1, sftlf@_FillValue, 0.0)
        else  ; mask ocean grid cells
            sftlf = where(sftlf.eq.0, sftlf@_FillValue, 0.0)
        end if

        ; apply mask
        mmm = mmm + conform(mmm, sftlf, (/1, 2/))
    end if

    ; *************************************************************************
    ; optionally calculate standared deviation
    ; *************************************************************************

    calc_stddev = False  ; default = do not calculate standard deviation

    if (isatt(diag_script_info, "stddev")) then
        calc_stddev = diag_script_info@stddev
    end if

    if (calc_stddev) then

        std = new(dimsizes(mmm), float)
        std = 0.0
        copy_VarMeta(mmm, std)
        std@var = var0 + "-std"
        std@long_name = "standard devation " + mmm@long_name

        do imod = 0, dim_MOD - 1
            regr_path = clim_diag_file_path(var0, field_type0, imod, \
                                            diag_script, opt)

            ; data have already been regridded

            info_output("Reading in file = " + regr_path, verbosity, 1)
            var = ncdf_read(regr_path, var0)

            std = where(ismissing(var), std, std + (var - mmm)^2)

            delete(var)
        end do ; imod-loop

        ; pts is reused from above and contains the number of valid models that have been averaged
        pts = pts - 1
        pts = where(pts.eq.0, pts@_FillValue, pts)

        std = sqrt(std / pts)

        ; optionally also apply mask to standard deviation

        if (lsmask.eq.1 .or. lsmask.eq.2) then
            std = std + conform(std, sftlf, (/1, 2/))
        end if

    end if  ; if calc_stddev

    ; *************************************************************************
    ; save results to netCDF
    ; *************************************************************************

    workpath = getenv("ESMValTool_wrk_dir")

    regr_path = workpath + "/" + basename(diag_script) + "/" + var0 + "_" + field_type0 + \
                "_" + models@experiment(0) + "_" + "mmm_" + selected_grid + "_" + timestring + ".nc"
    ncdf_mmm_outfile = ncdf_write(mmm, regr_path)

    if (calc_stddev) then
        regr_path = workpath + "/" + basename(diag_script) + "/" + var0 + "-std_" + field_type0 + \
                    "_" + models@experiment(0) + "_" + "mmm_" + selected_grid + "_" + timestring + ".nc"
        ncdf_mmm_outfile = ncdf_write(std, regr_path)
    end if

    ; ===================================================================================

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
