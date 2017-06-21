;;#############################################################################
;; CLOUDS_TAYLOR
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME EMBRACE
;;#############################################################################
;; Description
;;    Calculates the performance of models in reproducing 2-d fields of annual
;;    mean or seasonal (DFJ, MAM, JJA, SON) mean cloud properties
;;    (Taylor diagramms). The code is based on
;;       1) perfmetrics_main.ncl
;;       2) perfmetrics_taylor.ncl
;;       3) perfmetrics_taylor_collect.ncl
;;    originally written by Franziska Frank (DLR, Germany).
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding
;;                     - coarsest = coarsest grid of all input datasets
;;                     - finest   = highest resolution grid of all input datasets
;;                     - ref      = grid of reference datasets
;;                     - 1x1      = 1 deg x 1 deg grid
;;    valid_fraction:  minimum fraction (0-1) of valid data points per grid cell
;;                     (for creating a mask)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    calcmm:          - False (default) = do not calculate multi-model mean
;;                     - True = calculate multi-model mean
;;    models_to_skip:  name(s) of data set(s) to be skipped when calculating the
;;                     multi-model mean (usually all observations)
;;    embracelegend:   - False (default) = include legend in plot, max. 2 columns
;;                              with model names in legend
;;                     - True = write extra file with legend, max. 7 model names
;;                              per column in legend, alternative observational
;;                              dataset(s) will be plotted as a red star and labeled
;;                              "altern. ref. dataset" in legend (only if dataset is
;;                              of class "OBS")
;;    estimate_obs_uncertainty: - True = estimate observational uncertainties from
;;                     mean values (assuming fractions of obs. RMSE from documentation
;;                     of the obs data); only available for "CERES-EBAF", "MODIS",
;;                     "MODIS-L3"
;;                     - False = do not estimate obs. uncertainties from mean values
;;    mask_ts_sea_ice: - True = mask T < 272 K as sea ice (only for variable "ts")
;;                     - False = no additional grid cells masked for variable "ts"
;;    styleset:        "CMIP5", "DEFAULT"
;;                     (if not set, CLOUDS_TAYLOR will create a color table and
;;                     symbols for plotting)
;;    ref_model:       name of reference data set (if this attribute is not set
;;                     the variable attribute "ref_model" defined in the namelist
;;                     is used; see below)
;;                     note: if neither diag_script_info@ref_model nor the variable
;;                           attribute ref_model is set, an error message is issued
;;    timemean:        time averaging
;;                     - annualclim (default) = 1 plot annual mean
;;                     - seasonalclim = 4 plots (DJF, MAM, JJA, SON)
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Optional variable attributes (defined in namelist)
;;    ref_model:       name of reference data set
;;                     note: if diag_script_info@ref_model is defined,
;;                           diag_script_info@ref_model will be used and
;;                           the variable attribute ref_model will be ignored
;;
;; Caveats
;;    - support for irregular grids is not fully implemented yet
;;    - only 2-dim variables are currently supported
;;    - observational uncertainties are regridded like standard variables
;;    - for derived variables (e.g. SW_CRE), also the original variables have
;;      to be specified in the nameslist for the reference model (for the reference
;;      model only) in order to estimate the observational uncertainty
;;      (if estimate_obs_uncertainty = True)
;;    - the variable to plot must be the *first* variable in the namelist,
;;      the observational uncertainty (or the variables from which a derived
;;      variable has been calculated from) the *second* (third, ...)
;;
;; Modification history
;;    20170620-A_laue_ax: added tags for reporting
;;    20161104-A_laue_ax: changed calculation of correlation and standard deviation
;;                        ratio from unweighted to grid area weighted
;;    20160901-A_laue_ax: added regridding option 1 deg x 1 deg
;;    20160818-A_laue_ax: added support for variable attribute "ref_model"
;;    20160404-A_laue_ax: added optional drawing of observational
;;                        uncertainties
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150415-A-laue_ax: written.
;;
;;############################################################################

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"

load "./plot_scripts/ncl/taylor_diagram_less_hardcoded.ncl"
load "./plot_scripts/ncl/aux_plotting.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ; ========================================================================
    ; ========================= write references =============================
    ; ========================================================================

    write_references(diag_script,      \  ; script name
                     "A_laue_ax",      \  ; authors
                     "",               \  ; contributors
                     "D_lauer13jclim", \  ; diag_references
                     "",               \  ; obs_references
                     (/"P_embrace"/))     ; proj_references

    dim_MOD = dimsizes(models@name)
    dim_VAR = 1
    var0 = variables(0)
    field_type0 = field_types(0)

    ; time averaging: at the moment, only "annualclim" and "seasonalclim"
    ; are supported

    timemean = "annualclim"  ; default
    numseas = 1              ; default
    season = (/"annual"/)

    if (isatt(diag_script_info, "timemean")) then
        timemean = diag_script_info@timemean
        if (timemean.eq."seasonalclim") then
            numseas = 4
            delete(season)
            season = (/"DJF", "MAM", "JJA", "SON"/)
        end if
    end if

end

begin

    vardeffile = "interface_data/" + variables(0) + "_info.tmp"
    loadscript(vardeffile)

end

begin
    ; ========================================================================
    ; ========================== initialization ==============================
    ; ========================================================================

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    diag_script_base = basename(diag_script)
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if

    ; check for reference model definition
    ; first try diag_script_info@ref_model
    ; if diag_script_info@ref_model is not present, check variable attribute
    ; ref_model

    if (isatt(diag_script_info, "ref_model")) then
        ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
        if (ismissing(ref_ind)) then
            error_msg("f", diag_script, "", "reference model missing")
        end if
        ref_model = diag_script_info@ref_model(0)
    else
        ref_ind = get_ref_model_idx(models, var_attr_ref(0))
        if (ismissing(ref_ind)) then
            error_msg("f", diag_script, "", "reference model missing")
        end if
        ref_model = models@name(ref_ind)
    end if

    if (isatt(diag_script_info, "estimate_obs_uncertainty")) then
        estimate_obs_uncertainty = diag_script_info@estimate_obs_uncertainty
    else
        estimate_obs_uncertainty = False
    end if

    if (isatt(diag_script_info, "mask_ts_sea_ice")) then
        mask_ts_sea_ice = diag_script_info@mask_ts_sea_ice
    else
        mask_ts_sea_ice = False
    end if

    if (isatt(diag_script_info, "filename_add")) then
        filename_add = "_" + diag_script_info@filename_add
    else
        filename_add = ""
    end if

    ; basename of diag_script
    diag_script_base = basename(diag_script)

    ; output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ; set regridding method
    if (.not.(isatt(diag_script_info, "grid"))) then
        error_msg("w", diag_script, "", "Regridding method not defined. " + \
                  "Using coarsest")
        selected_grid = "coarsest"
    else
        selected_grid = diag_script_info@grid
    end if

    ;; Check for irregular grid
    l_irreg = new(dim_MOD, logical)
    l_irreg = False
    do imod = 0, dim_MOD - 1
        data = read_data(imod, var0, field_type0)
        cnt = interface_get_dictkeys(var0, imod) + "_dimNames"
        dims = data@$cnt$
        ndims = dimsizes(dims)
        if (ndims.lt.2) then
            continue
        end if
        if (ndims.eq.2 .and. dims(0).eq."time") then
            continue
        end if
        if (dims(ndims - 1).ne."lon" .and. dims(ndims - 2).ne."lat") then
            l_irreg(imod) = True
            continue
        end if
    end do
    delete([/data, cnt, dims, ndims/])

    ;; Define target grid for regridding
    if (any(l_irreg) .and. selected_grid.ne."1x1") then
        error_msg("f", diag_script, "", "models are on irregular grid, " + \
                  "diag_script_info@grid = '1x1' must be used")
    end if

    ; defined the target grid for regridding

    models_ind = ispan(0, dim_MOD - 1, 1)
    info_output("Setting target grid:", verbosity, 2)
    if (selected_grid.eq."ref") then
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

    if (isatt(diag_script_info, "embracelegend")) then
        embracelegend = diag_script_info@embracelegend
    else
        embracelegend = False
    end if

    ;; Selected region
    if (isatt(diag_script_info, "region")) then
        region = select_region(diag_script_info@region)
    end if

    ; ========================================================================
    ; ====================== regridding + creating mask ======================
    ; ========================================================================

    opt = True
    opt@plot_type = "regridded"
    opt@time_avg = selected_grid

    numvars = dimsizes(variables)

    do ivar = 0, numvars - 1
        ; process all models for first variable
        if (ivar.eq.0) then
            modelf = 0
            modell = dim_MOD - 1
        ; process only reference model for error estimates
        else
            modelf = ref_ind
            modell = ref_ind
        end if

        do imod = modelf, modell
            info_output("Regridding " + models@name(imod), verbosity, 2)
            regr_path = clim_diag_file_path(variables(ivar), field_types(ivar), \
                                     imod, diag_script, opt)

            ;; Do not process data if already processed
            if (isfilepresent(regr_path) .and. \
                getenv("ESMValTool_force_calc").eq."False") then
                pre_var = ncdf_read(regr_path, variables(ivar))
                dnames = getVarDimNames(pre_var)
            else
                ;; Extract model corresponding data
                data = read_data(imod, variables(ivar), field_types(ivar))
                var = extract_data_extended(imod, var0, data, -1, 0, 0)
                dnames = getVarDimNames(var)
    
            if (l_irreg(imod)) then  ; irregular grid

                ;; Regrid using ESMF_regrid

                ;; ESMF_regrid cannot handle integer coordinates
                do ii = 0, dimsizes(dnames) - 1
                    if (dnames(ii).eq."time") then
                        continue
                    end if
                    if (typeof(var&$dnames(ii)$).eq."integer") then
                        var&$dnames(ii)$ := tofloat(var&$dnames(ii)$)
                    end if
                end do

                ;; Pass information of irregular grid as attributes of var
                ;;  - lat/lon grid information comes from another file
                ;; this is NOT CMOR standard! Reformat below.
                var@lat2d := extract_area(imod, data, "lat", var0)
                var@lon2d := extract_area(imod, data, "lon", var0)

                ;; Naming of intermediate regrid files: will be re-used if
                ;; already present - ensure unambiguous naming!
                ;; Because of the mask of missing values, even different
                ;; time-steps / model levels need their own wgt file
                stem   = interface_get_inpaths(imod) + \
                    "/ESMF_regrid_" + models@project(imod) + "_" + \
                    models@name(imod) + "_" + models@experiment(imod) + "_" + \
                    models@ensemble(imod)
                ending = var0 + "_" + models@start_year(imod) + \
                    "-" + models@end_year(imod) + ".nc"

                ;; Set options
                Opt = True
                Opt@SrcMask2D  = where(ismissing(var(0, :, :)), 0, 1)
                Opt@DstGridLat = grid&lat
                Opt@DstGridLon = grid&lon
                ;Opt@DstGridType    = "rectilinear"  ; problem with that Opt
                Opt@ForceOverwrite = True
                Opt@PrintTimings   = True
                Opt@Debug          = True
                Opt@InterpMethod   = "patch"  ; "bilinear", "patch", "conserve"

                ;; Save regrid files for later use: careful with naming!
                Opt@WgtFileName = stem + "_Wgt_" + var0 + ".nc"
                Opt@SrcFileName = stem + "_Src_" + ending
                Opt@DstFileName = stem + "_Dst_" + ending

                ;; Reuse regrid files if present
                if (isfilepresent(Opt@SrcFileName)) then
                    Opt@SkipSrcGen = True
                end if
                if (isfilepresent(Opt@DstFileName)) then
                    Opt@SkipDstGen = True
                end if
                if (isfilepresent(Opt@WgtFileName)) then
                    Opt@SkipWgtGen = True
                    var_hreg = ESMF_regrid_with_weights(var, \
                                                        Opt@WgtFileName, Opt)
                else
                    var_hreg = ESMF_regrid(var, Opt)
                end if
                delete(Opt)

            else  ; regular grid

                ;; Determine dimensions to regrid
                zregrid = False  ; plev
                xregrid = False  ; lon
                yregrid = False  ; lat
                dnames = getVarDimNames(var)
                if (any(dnames.eq."lon")) then
                    xregrid = True
                    if (dimsizes(var&lon).eq.dimsizes(grid&lon)) then
                        if (max(abs(var&lon - grid&lon)).lt.1.0d-5) then
                            xregrid = False
                        end if
                    end if
                end if
                if (any(dnames.eq."lat")) then
                    yregrid = True
                    if (dimsizes(var&lat).eq.dimsizes(grid&lat)) then
                        if (max(abs(var&lat - grid&lat)).lt.1.0d-5) then
                            yregrid = False
                        end if
                    end if
                end if

                ;; Regrid in the horizontal (if necessary)
                if (xregrid.and.yregrid) then
                    var_hreg = linint2_Wrap(var&lon, var&lat, var, True, \
                                            grid&lon, grid&lat, 0)
                else if (xregrid.and..not.yregrid) then
                    dx = ind(dnames.eq."lon")
                    var_hreg = linint1_n_Wrap(var&lon, var, True, grid&lon, 0, dx)
                    info_output("Latitude regridding not required", verbosity, 4)
                else if (.not.xregrid.and.yregrid) then
                    dx = ind(dnames.eq."lat")
                    var_hreg = linint1_n_Wrap(var&lat, var, False, grid&lat, 0, dx)
                    info_output("Longitude regridding not required", verbosity, 4)
                else if (.not.xregrid .and. .not.yregrid) then
                    var_hreg = var
                    info_output("Horizontal regridding not required", verbosity, 4)
                end if
                end if
                end if
                end if

            end if  ; regular or irregular grid?

                delete(data)
                delete(var)

                pre_var = var_hreg
                delete(var_hreg)

                ;; Output to NetCDF
                pre_var@var = variables(ivar)
                pre_var@diag_script = (/diag_script/)
                pre_var@ncdf = regr_path
                ncdf_outfile = ncdf_write(pre_var, regr_path)
            end if

            ;; Create a missing value mask for each model
            ;; The dim_avg_wgt function with option limit is used on the time
            ;; coordinate: in this way each grid box containing more than
            ;; valid_fraction of missing values along the time series is
            ;; set to missing.

            if ((imod.eq.0).and.(ivar.eq.0)) then
                global_mask = dim_avg_n_Wrap(pre_var, 0)  ; to get rid of time dim
                global_mask = 0.
            end if

            ww = pre_var&time
            ww = 1.
            limit  = toint(diag_script_info@valid_fraction * \
                           dimsizes(pre_var&time))

            if (limit.lt.1) then
                limit = 1
            end if

            ; mask sea ice (variable ts)?
            if (mask_ts_sea_ice) then
                if ((var0.eq."ts").and.(imod.eq.ref_ind).and.(ivar.eq.0)) then
                    ; temperatures below 272 K are sea ice ---> mask
                    pre_var = where(pre_var.lt.272.0, pre_var@_FillValue, pre_var)
                end if
            end if

            ;; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1
            if (dimsizes(dimsizes(pre_var)).eq.1) then
                tmp = dim_avg_wgt_Wrap(pre_var, ww, limit)
            else
                tmp = dim_avg_wgt_n_Wrap(pre_var, ww, limit, 0)
            end if

            tmp = where(ismissing(tmp), tmp, 0.)
            if ((all(ismissing(tmp))).and.(ivar.eq.0)) then
                error_msg("w", diag_script, "", models@name(imod) + \
                          "contains only missing values, excluding from global " + \
                          "masking.")
            else
                global_mask = global_mask + tmp
            end if

            delete(tmp)
            delete(ww)
            delete(pre_var)
        end do  ; imod
    end do  ; ivar

    ; ========================================================================
    ; ============================ statistics ================================
    ; ========================================================================

    ;; Define taylor file directory
    work_dir = getenv("ESMValTool_wrk_dir")
    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    ncdf_dir = get_output_dir(work_dir, diag_script_base) + "taylor_" \
                              + var0 + ".nc"

    ; calculate multi-model mean (if requested in cfg_clouds_taylor.ncl)

    flag_mean = False
    if (isatt(diag_script_info, "calcmm")) then
        if (diag_script_info@calcmm) then
            flag_mean = True
            dim_MOD = dim_MOD + 1
            mm_ind = 0

            modelsatname = array_append_record((/"multi-model-mean"/), \
                                               models@name, 0)
            delete(models@name)
            models@name = modelsatname
            delete(modelsatname)

            projectsatname = array_append_record((/"mmm"/), models@project, 0)
            delete(models@project)
            models@project = projectsatname
            delete(projectsatname)

            ; update index of reference data set
            ref_ind = ind(models@name.eq.ref_model)

            model_ind_woref = ind(ref_model.ne.models@name)

            ; check for model(s) to skip (typically observation(s)) and multi-model mean entry

            if (isatt(diag_script_info, "models_to_skip")) then
                do imod = 0, dim_MOD - 1
                    if (any(models@name(imod) .eq. diag_script_info@models_to_skip)) then
                        idx = ind(model_ind_woref.eq.imod)
                        if (.not.ismissing(idx)) then
                            model_ind_woref(idx) = -999
                        end if
                        print("skipped model " + models@name(imod) + " in multi-model mean calculation")
                    end if
                end do
            end if

            ; also exclude multi-model mean entry from averaging
            idx = ind(model_ind_woref.eq.mm_ind)
            if (.not.ismissing(idx)) then
                model_ind_woref(idx) = -999
            end if

            tmpidx = ind(model_ind_woref.ge.0)
            tmp = model_ind_woref(tmpidx)
            delete(model_ind_woref)
            model_ind_woref = tmp
            delete(tmp)
        end if
    end if

    ; debug output
;    system ("rm debug.nc")
;    debugfile = addfile("debug.nc","c")
;    debugfile->mask = global_mask

    ; read data
    if (isdefined("data")) then
        delete(data)
    end if

    do ivar = 0, numvars - 1
        ; process all models for first variable
        if (ivar.eq.0) then
            modelf = 0
            modell = dim_MOD - 1
        ; process only reference model for error estimates
        else
            modelf = ref_ind
            modell = ref_ind
        end if

        do imod = modelf, modell
            if (flag_mean) then
                if (imod.eq.0) then
                    new_path = "none"
                else
                    new_path = clim_diag_file_path(variables(ivar), field_types(ivar), \
                                  imod - 1, diag_script, opt)
                end if
            else
                new_path = clim_diag_file_path(variables(ivar), field_types(ivar), \
                                  imod, diag_script, opt)
            end if

            if (isfilepresent_esmval(new_path)) then
                info_output("Read in: " + new_path, verbosity, 1)
                data_temp = ncdf_read(new_path, variables(ivar))

                ; ivar.eq.0  ---> variable data
                ;     .ne.0  ---> uncertainty estimates (observations)
                ; The uncertainty estimates are assumed to be given as
                ; 1-sigma of normally distributed error estimates
                ; Note: uncertainty estimates must have "err" in their
                ;       variable name, otherwise variables are assumed
                ;       to be fields used for calculating derived variables
                ;       (in this case, fields are not squared before averaging)
                ; ---> square sigma before averaging over time, then
                ;      calculate square-root to get average sigma


                if ((ivar.ne.0).and.isStrSubset(variables(ivar), "err")) then
                    data_temp = data_temp * data_temp
                end if

                data1 = time_operations(data_temp, -1, -1, "average", \
                                    timemean, True)

                if ((ivar.ne.0).and.isStrSubset(variables(ivar), "err")) then
                    data1 = sqrt(data1)
                end if

                if (numseas.eq.1) then
                    data1 = data1 + global_mask
                else
                    do is = 0, numseas - 1
                        data1(is, :, :) = data1(is, :, :) + global_mask
                    end do
                end if

                if (ivar.eq.0) then
                    if (.not.isdefined("data")) then
                        dim_data = array_append_record((/dim_MOD/), dimsizes(data1), 0)
                        data = new(dim_data, float)
                        data!0 = "model"
                        data&model = models@name
                        if (numseas.eq.1) then
                            copy_VarCoords(data1, data(imod, :, :))
                        else
                            copy_VarCoords(data1, data(imod, :, :, :))
                        end if
                    end if
                    dim_data = dimsizes(data)
                    rank = dimsizes(dim_data)
                    if (numseas.eq.1) then
                        data(imod, :, :) = data1
                    else
                        data(imod, :, :, :) = data1
                    end if
                    delete(data_temp)
                    delete(data1)
                else
                    if (.not.isdefined("err")) then
                        dim_data = array_append_record((/numvars - 1/), dimsizes(data1), 0)
                        err = new(dim_data, float)
                        err!0 = "var"
                        err&var = variables(1:numvars - 1)
                    end if
                    if (numseas.eq.1) then
                        err(ivar - 1, :, :) = data1
                    else
                        err(ivar - 1, :, :, :) = data1
                    end if
                end if

                ; debug output

;                debugfile->$models@name(imod)$ = data1
            else
                if ((flag_mean).and.(imod.eq.mm_ind)) then
                else
                    ; skip this model as no data are available

                    info_output("No data available for model:" + \
                                models@name(imod), verbosity, 1)
                    info_output("Skipping", verbosity, 1)
                end if
            end if
        end do  ; imod loop
    end do  ; ivar loop

    ; number of used reference data sets

    dim_REF = 1 ; dimsizes(diag_script_info@ref_model)

    ; define result variable

    val = new((/dim_MOD - dim_REF, 2, numseas/), float)
    val!0 = "models"
    val!1 = "statistic"
    val!2 = "time"
    val&models = models@name(ind((models@name.ne.ref_model)))
    val&statistic = (/"stddev_ratio", "correlation"/)
    val&time = ispan(0, numseas - 1, 1)
    val@ncdf_dir = ncdf_dir

    ; extract data

    do is = 0, numseas - 1

        if (numseas.eq.1) then
            ref = data(ref_ind, :, :)
        else
            ref = data(ref_ind, is, :, :)
        end if

        ; loop over models
        ; (second model index, for models only, skip reference)

        mod_i = 0
        do imod = 0, dim_MOD - 1

            ; skip reference

            if (imod.eq.ref_ind) then
                print(" Refmodel = " + data&model(imod))
                continue
            end if

            ; multi-model mean

            if (flag_mean.and.imod.eq.mm_ind) then
                ; special case: if there is only one model, skip averaging

                if (dimsizes(model_ind_woref).eq.1) then
                    if (numseas.gt.1) then
                        var = data(model_ind_woref, :, :, :)
                    else
                        var = data(model_ind_woref, :, :)
                    end if
                else
                    if (numseas.gt.1) then
                        var = dim_avg_n_Wrap(data(model_ind_woref, is, :, :), 0)
                    else
                        var = dim_avg_n_Wrap(data(model_ind_woref, :, :), 0)
                    end if
                end if

            else  ; Given model
                if (numseas.eq.1) then
                    var = data(imod, :, :)
                else
                    var = data(imod, is, :, :)
                end if
            end if

            ; calculate rate of standard deviations

;            val(mod_i, 0, is) = calculate_metric(var, ref, "stddev_ratio_nowgt")
            val(mod_i, 0, is) = calculate_metric(var, ref, "stddev_ratio")

            ; calculate pattern correlation

;            val(mod_i, 1, is) = calculate_metric(var, ref, "correlation_nowgt")
            val(mod_i, 1, is) = calculate_metric(var, ref, "correlation")

            delete(var)
            mod_i = mod_i + 1

        end do  ; loop over models
    end do  ; loop over seasons

    ; attach attributes to the results

    val@title = "taylor metrics"
    val@long_name = "Taylor diagram"
    val@diag_script = (/diag_script/)
    val@var = "var"

    ; write NetCDF output

    ncdf_outfile = ncdf_write(val, val@ncdf_dir)

    ; ---------------------- calculate/estimate RMS of observations -------------------

    rmsobs = new((/numseas/), float)

    ; optionally, if there are no variable(s) available for the error estimate,
    ; the reference data can be used to create one for selected data sets

    if (estimate_obs_uncertainty) then
        if (any(ref_model.eq.(/"CERES-EBAF", "MODIS", "MODIS-L3"/))) then
            if (numvars.eq.1) then  ; no variable(s) for error estimate
                numvars = 2
                variables0 = (/variables(0), variables(0)/)
                delete(variables)
                variables = variables0
                dim_data = dimsizes(data)
                dim_data(0) = numvars - 1
                err = new(dim_data, float)
                if (numseas.eq.1) then
                    err(0, :, :) = data(ref_ind, :, :)
                else
                    err(0, :, :, :) = data(ref_ind, :, :, :)
                end if
            end if
        end if
    end if  ; if estimate_obs_uncertainty

    do is = 0, numseas - 1
        if (numseas.eq.1) then
            ref = data(ref_ind, :, :)
        else
            ref = data(ref_ind, is, :, :)
        end if

        ; This code is equivalent to the function "calculate_metric" (statistics.ncl)
        weights = map_area(ref&lat, ref&lon)
        ; optional: no weights --> reset weights
;        wgt1d = 1.0

        ; convert to 1-D arrays
        ref1d = ndtooned(ref)
        wgt1d = ndtooned(weights)
        wgt1d = wgt1d / dim_avg_n_Wrap(wgt1d, 0)

        rmsobs(is) = 0.0

        do iv = 0, numvars - 2
            if (numseas.eq.1) then
                var = err(iv, :, :)
            else
                var = err(iv, is, :, :)
            end if

            ; If the variable name contains "err", the variable is assumed to contain
            ; an estimate of the observational uncertainty (e.g. UWisc, GPCP-SG, and
            ; SSMI-MERIS provide error estimates, MODIS-L3-C6 provides error estimates
            ; for lwp and iwp but not for clt). In this case, there is nothing to do
            ; here.
            ; If the variable name does not contain "err", the variable is assumed to
            ; be a field (e.g. rsut, rsutcs) used to calculate a derived variable
            ; (e.g. CRE_SW). In this case, the observational uncertainty can be estimated
            ; if the fractional average error is known (e.g. CERES-EBAF, MODIS).

            if (isStrSubset(variables(iv + 1), "err")) then
                ; variable contains error estimate provided by data set
                ; ---> nothing to do
                factor = 1.0
            else
                ; variable(s) contain(s) fields used to derive the variable plotted
                ; --> observational uncertainty might be estimated from these fields
                ;     if the average error(s) of these fields are known

                ; by default, uncertainty cannot be estimated

                factor = 0.0

                ; for CERES-EBAF, we estimate the errors from the original variables
                ; (rsut+rsutcs / rlut+rlutcs) assuming variable specific factors

                if (ref_model.eq."CERES-EBAF") then
                    if (any(variables(iv + 1).eq.(/"rsut", "rsutcs"/))) then
                        ; From the "Technical Note on CERES EBAF Ed2.7":
                        ; rsutcs: 1% ~ 0.5 W/m2
                        ;         estimate total error = 2.6 W/m2 ---> ~5.2%
                        ; rsut: 4-5 W/m2 ---> ~4-5%
                        ;
                        ; ---> use 5% for rsut and rsutcs
                        factor = 0.05
                    end if
                    if (any(variables(iv + 1).eq.(/"rlut", "rlutcs"/))) then
                        ; From the Technical Note on CERES EBAF Ed2.7:
                        ; "For CERES, calibration uncertainty is 0.5% (1 sigma)
                        ; which for a typical global mean clear-sky LW flux corresponds
                        ; to ~1 W/m2. [...] The total error in TOA outgoing clear-sky LW radiation
                        ; in a region is estimated as sqrt(1^2 + 1.74^2 + 0.7^2 + 1^2 + 2.75^2)
                        ; or approximately 3.6 W/m2."
                        ;
                        ; --> here, we estimate the total error as 0.5% * 3.6 W/m2 / 1.0 W/m2 = 1.8%
                        factor = 0.018
                    end if
                end if

                if (ref_model.eq."MODIS") then
                    if (any(variables(iv + 1).eq.(/"clt"/))) then
                        ; From the technical note on obs4MIPS MODIS satellite data
                        ; "Moderate Resolution Imaging Spectroradiometer (MODIS)
                        ; Cloud Fraction Technical Document":
                        ;
                        ; Global day/night fractional agreement of cloud detection between
                        ; MODIS and CALIPSO lidar (CALIOP):
                        ; Aug 2006 (clear):  0.84
                        ; Aug 2006 (cloudy): 0.88
                        ; Feb 2006 (clear):  0.85
                        ; Feb 2006 (cloudy): 0.88
                        ;
                        ; ---> using fractional agreement as error estimate (1 sigma)
                        ; ---> average error estimate = sqrt((0.84^2 + 0.88^2 + 0.85^2 + 0.88^2) / 4) = 0.86
                        ; ---> factor = 1 - 0.86 = 0.14
                        factor = 0.14
                    end if
                end if

                if (ref_model.eq."MODIS-L3") then
                    if (any(variables(iv + 1).eq.(/"clt"/))) then
                        factor = 0.0  ; no error estimates for clt from MODIS_C6_L3
                    end if
                end if
            end if

            ; convert to 1-D arrays
            var1d = ndtooned(var * factor)

            ; --- RMS ---
            ; in case of derived variables (e.g., SW_CRE), the squared error
            ; estimates of all input variables are added before calculating
            ; the square root
 
            rmsobs(is) = rmsobs(is) + dim_avg_wgt_Wrap(var1d^2, wgt1d, 1)
        end do  ; iv-loop

        ; finish calculation of RMS
        rmsobs(is) = sqrt(rmsobs(is))

        ; normalize RMS by dividing by standard deviation of obs (ref)
        avg_ref = dim_avg_wgt_Wrap(ref1d, wgt1d, 1)
        rmsobs(is) =  rmsobs(is) / sqrt(dim_avg_wgt_Wrap((ref1d - avg_ref)^2, wgt1d, 1))

    end do  ; loop over seasons

    if (isvar("var")) then
        delete(var)
    end if
    if (isvar("ref")) then
        delete(ref)
    end if
    if (isvar("weights")) then
        delete(weights)
    end if
    if (isvar("var1d")) then
        delete(var1d)
    end if
    if (isvar("ref1d")) then
        delete(ref1d)
    end if
    if (isvar("wgt1d")) then
        delete(wgt1d)
    end if

    ; ========================================================================
    ; ============================= plotting =================================
    ; ========================================================================

    if (getenv("ESMValTool_write_plots").eq.True) then

        if (flag_mean) then
            nm = dim_MOD - 1
        else
            nm = dim_MOD
        end if

        climofiles = new(nm, string)

        do ii = 0, nm - 1
            climofiles(ii) = \
                interface_get_inpaths(ii) + "/" \
                + interface_get_infile(variables(0), field_types(0), ii)
        end do

        if (numvars.gt.1) then
            climofilesx = array_append_record(climofiles, \
                interface_get_inpaths(0) + "/" + \
                interface_get_infile(variables(1), field_types(1), 0), 0)
            delete(climofiles)
            climofiles = climofilesx
            delete(climofilesx)
        end if

        alltags = array_append_record(tags, (/"PT_taylor", "DM_global"/), 0)

        ; create outfile directory
        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "clouds_taylor")

        colors  = new(mod_i, string)
        markers = new(mod_i, integer)

        ratio = new((/mod_i, 1, numseas/), float)  ; standard deviation ratios
        ratio = val(:, 0, :)
        cc = new((/mod_i, 1, numseas/), float)    ; correlations
        cc = val(:, 1, :)

        if (isatt(diag_script_info, "styleset")) then
            colortab = project_style(diag_script_info, "colors")
            markertab = project_style(diag_script_info, "markers")
        else
            colortab = (/"(/0.00, 0.00, 0.59/)", "(/0.00, 0.39, 1.00/)", \
                         "(/0.20, 1.00, 1.00/)", "(/0.20, 0.88, 0.00/)", \
                        "(/1.00, 0.88, 0.00/)", "(/1.00, 0.59, 0.00/)", \
                         "(/1.00, 0.20, 0.00/)", "(/0.59, 0.00, 0.00/)", \
                         "(/0.78, 0.00, 0.78/)", "(/0.59, 0.00, 0.59/)", \
                         "(/0.90, 0.90, 0.90/)", "(/0.70, 0.70, 0.70/)", \
                         "(/0.50, 0.50, 0.50/)", "(/0.30, 0.30, 0.30/)"/)
            markertab = (/16, 4, 5, 0/)
        end if

        do is = 0, numseas - 1
            if (isvar("wks")) then
                delete(wks)
            end if

            if (isvar("legendwks")) then
                delete(legendwks)
            end if

            if (isvar("plot")) then
                delete(plot)
            end if

            wks = get_wks("dummy_for_wks", diag_script, var0 + "_" + season(is) + filename_add)

            ; create new marker: filled star

            mstring = "z"
            fontnum = 35
            size    = 1.5
            angle   = 0.0

            new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, \
                                     1.0, size, angle)

            if (embracelegend) then
                legendwks = get_wks("legend_wks", diag_script, "legend" + filename_add)
                new_index = NhlNewMarker(legendwks, mstring, fontnum, 0.0, 0.0, \
                            1.0, size, angle)
            else
                legendwks = wks
            end if

            i = 0
            idx = 0

            if (isatt(diag_script_info, "styleset")) then
                do while (i.lt.mod_i)

                    colors(i)  = colortab(i)
                    markers(i) = markertab(i)

                    ; override multi-model mean
                    if ((flag_mean).and.(i.eq.0)) then
                        colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
                        markers(i) = new_index
                    end if
                    ; override alternative obs (if "EMBRACE" legend)
                    if (embracelegend) then
                        if ((isStrSubset(str_lower(models@project(i)), "obs"))) then
                            if (flag_mean) then
                                colors(i)  = "(/1.00, 0.00, 0.00/)"  ; red
                            else
                                colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
                            end if
                            markers(i) = new_index
                            models@name(i) = "altern. ref. dataset"
                            val&models(i) = models@name(i)
                        end if
                    end if
                    i = i + 1
                end do
            else
                do while (i.lt.mod_i)
                    if ((flag_mean).and.(i.eq.0)) then
                        colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
                        markers(i) = new_index
                        i = i + 1
                    else
                        if (embracelegend) then
                            colors(i)  = colortab(idx)
                            if ((isStrSubset(str_lower(models@project(i)), "obs"))) then
                                if (flag_mean) then
                                    colors(i)  = "(/1.00, 0.00, 0.00/)"  ; red
                                else
                                    colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
                                end if
                                markers(i) = new_index
                                models@name(i) = "altern. ref. dataset"
                                val&models(i) = models@name(i)
                            else
                                markers(i) = markertab(mod(i,2))
                            end if
                            i = i + 1
                            if (i.ge.mod_i) then
                                break
                            end if
                            idx = idx + 1
                            if (idx.ge.dimsizes(colortab)) then
                                idx = dimsizes(colortab) - 1
                            end if
                        else
                            do n = 0, mod_i / dimsizes(colortab)
                                colors(i)  = colortab(idx)
                                markers(i) = markertab(n)
                                i = i + 1
                                if (i.ge.mod_i) then
                                    break
                                end if
                            end do
                            idx = idx + 1
                            if (idx.ge.dimsizes(colortab)) then
                                idx = dimsizes(colortab) - 1
                            end if
                        end if
                    end if
                end do
            end if

            ropts = True
            ropts@Colors          = colors
            ropts@Markers         = markers  ; marker styles
            ropts@gsMarkerSizeF   = 0.0125

            modelnames            = val&models
            ropts@caseLabels      = modelnames

            ropts@stnRad          = (/ 0.5, 1.5, 2.0, 2.5 /)
            ropts@OneX            = "1.00"
            ropts@ccRays          = (/ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, \
                                       0.9, 0.95, 0.99 /)
            ropts@ccRays_color    = "Black"
            ropts@centerDiffRMS   = True
            ropts@centerDiffRMS_color = "LightGray"  ; "PaleGreen1"

            ropts@printDataLabels = False  ; print numeric label for each data pt.
            ropts@reverseOrder    = True  ; reverse order of items in legend
            ropts@legendBox       = True  ; draw box around legend

            ; legend parameters

            ropts@legendXpos            = 0.575
            ropts@legendYpos            = -0.5
            ropts@legendWidth           = 0.2

            ropts@caseLabelsFontHeightF = 0.1
            ropts@plotSize              = 0.6

            if (embracelegend) then
                ropts@legendExtraFile = True  ; save legend to extra file
                ropts@modelsperrow = 7  ; 7 models per column in the legend
            else
                ropts@legendExtraFile = False
                if (dim_MOD.ge.20) then
                    ropts@caseLabelsFontHeightF = 0.06
                    ropts@plotSize              = 0.45
                end if
            end if

            ; plot

            ropts@tiMainString    = var0 + " (" + season(is) + ")"  ; title
            ; observational uncertainties as normalized RMS error
            ropts@rmsobs = rmsobs(is)

            plot = taylor_diagram(wks, legendwks, ratio(:, :, is), cc(:, :, is), ropts)

            info_output(" Wrote " + wks@fullname, verbosity, 1)

            ; add meta data to plot (for reporting)

            caption = "Taylor diagram for variable " + variables(0) + " (" + season(is) \
                      + "), reference = " + models@name(ref_ind) + "."
            id = diag_script + "_" + season(is)

            ESMValMD(wks@fullname, alltags, caption, id, variables(0), models@name, climofiles)

            if (embracelegend.and.(is.eq.(numseas-1))) then
                frame(legendwks)
            end if
        end do  ; is-loop (seasons)
    end if

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
