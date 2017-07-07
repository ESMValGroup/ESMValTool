;; CLOUDS_IPCC
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME EMBRACE
;;#############################################################################
;; Description
;;    Calculates the multi-model mean bias of annual mean 2-d cloud variables
;;    compared with a reference data set (observations). In addition,
;;    zonal averages of the individual models, the multi-model mean and the
;;    reference data set (observations) are calculated.
;;    If more than one variable is specified in the namelist, variables 2 to n
;;    are assumed to be observational error estimates. These error estimates are
;;    added to the zonal mean plot as light red shading.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref, 1x1)
;;    ref_model:       name of reference data set
;;    valid_fraction:  minimum fraction of valid data points
;;                     (for creating a mask)
;;    projection:      map projection, e.g., Mollweide, Mercator
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    mask_ts_sea_ice: - True = mask T < 272 K as sea ice (only for variable "ts")
;;                     - False = no additional grid cells masked for variable "ts"
;;    models_to_skip:  name(s) of data set(s) to be skipped when calculating the
;;                     multi-model mean (usually all observations)
;;    ncdf:            enable to output to netCDF; either use "default" or give
;;                     a full file name
;;    styleset:        style set for zonal mean plot ("CMIP5", "DEFAULT")
;;    timemean:        time averaging, i.e. "seasonalclim" (DJF, MAM, JJA, SON),
;;                     "annualclim" (annual mean)
;;
;; Required variable_info attributes (variable specific)
;;    long_name:       description of variable
;;    map_diff_Levels: contour levels for difference plot
;;    units:           variable units
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;    Specifying more than one data set for the observational uncertainties
;;    my lead to unexpected or undefined results.
;;
;; Modification history
;;    20170620-A_laue_ax: added tags for reporting
;;    20160920-A_laue_ax: added optional shading of observational uncertainties
;;                        to the zonal mean plot
;;    20160901-A_laue_ax: added regridding option 1 deg x 1 deg
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150428A-laue_ax: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"

load "./plot_scripts/ncl/contour_maps.ncl"

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
                     "",               \  ; diag_references
                     "",               \  ; obs_references
                     (/"P_embrace"/))     ; proj_references

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
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

    if (.not.(isatt(diag_script_info, "ref_model"))) then
        error_msg("f", diag_script, "", "no reference model is selected")
    end if

    ; set reference model

    ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
    ref_model = diag_script_info@ref_model(0)
    if (ismissing(ref_ind)) then
        error_msg("f", diag_script, "", "no adequate reference model provided")
    end if

    if (dimsizes(diag_script_info@ref_model).gt.1) then
        ref_ind2 = ind(models@name.eq.diag_script_info@ref_model(1))
        ref_model2 = diag_script_info@ref_model(1)
    else
        ref_ind2 = -1
    end if

    if (isatt(diag_script_info, "mask_ts_sea_ice")) then
        mask_ts_sea_ice = diag_script_info@mask_ts_sea_ice
    else
        mask_ts_sea_ice = False
    end if

    ; basename of diag_script

    diag_script_base = basename(diag_script)

    ; output plot directory

    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ; set regridding method

    if (.not.(isatt(diag_script_info, "grid"))) then
        error_msg("w", diag_script, "", "Regridding method not defined. " + \
                  "Using 1x1")
        selected_grid = "1x1"
    else
        selected_grid = diag_script_info@grid
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

    ; ========================================================================
    ; ============================ regridding ================================
    ; ========================================================================

    opt = True
    opt@plot_type = "regridded"
    opt@time_avg = selected_grid
    do imod = 0, dim_MOD - 1

        info_output("Regridding " + models@name(imod), verbosity, 2)
        regr_path = clim_diag_file_path(var0, field_type0, imod, \
                                        diag_script, opt)

        ; do not process data if already processed

        force = getenv("ESMValTool_force_calc")
        if (ismissing(force)) then
            force = "False"
        end if

        if (isfilepresent(regr_path) .and. force.eq."False") then

;            info_output("Reading in file = " + regr_path, verbosity, 1)
            pre_var = ncdf_read(regr_path, var0)
            dnames = getVarDimNames(pre_var)

        else

            ; extract model corresponding data

            data = read_data(imod, var0, field_type0)
            var = extract_data(imod, data, -1, 0, 0)
            delete(data)

            ; determine dimensions to regrid

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

            ; regrid in the horizontal (if necessary)

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

            delete(var)

            pre_var = var_hreg
            delete(var_hreg)

            ; output to NetCDF
            pre_var@var = var0
            pre_var@diag_script = (/diag_script/)
            pre_var@ncdf = regr_path
            ncdf_outfile = ncdf_write(pre_var, regr_path)

        end if

        ; create a missing value mask for each model

        ; The dim_avg_wgt function with option limit is used on the time
        ; coordinate: this way each grid box containing more than
        ; 'valid_fraction' of missing values along the time series is
        ; set to missing.

        if (imod.eq.0) then
            global_mask = dim_avg_n_Wrap(pre_var, 0)  ; to get rid of time dim
            global_mask = 0.
        end if

        ww = pre_var&time
        ww = 1.
        limit = toint(diag_script_info@valid_fraction * dimsizes(pre_var&time))

        if (limit.lt.1) then
            limit = 1
        end if

        ; mask sea ice (variable ts)?
        if (mask_ts_sea_ice) then
            if ((var0.eq."ts").and.(imod.eq.ref_ind)) then
                ; temperatures below 272 K are sea ice ---> mask
                pre_var = where(pre_var.lt.272.0, pre_var@_FillValue, pre_var)
            end if
        end if

        ; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1

        if (dimsizes(dimsizes(pre_var)).eq.1) then
            tmp = dim_avg_wgt_Wrap(pre_var, ww, limit)
        else
            tmp = dim_avg_wgt_n_Wrap(pre_var, ww, limit, 0)
        end if

        tmp = where(ismissing(tmp), tmp, 0.)
        if (all(ismissing(tmp))) then
            error_msg("w", diag_script, "", models@name(imod) + \
                      "contains only missing values, Excluding from global " + \
                      "masking.")
        else
            global_mask = global_mask + tmp
        end if
        delete(tmp)

        delete(ww)
        delete(pre_var)

    end do  ; imod

    ; additional variables (other than var0) are assumed to be uncertainty estimates
    ; of observation(s)

    do ivar = 1, dim_VAR - 1
        do imod = 0, dim_MOD - 1
            in_file = interface_get_fullpath(variables(ivar), field_types(ivar), imod)
            if (isfilepresent(in_file)) then
                info_output("Regridding " + models@name(imod), verbosity, 2)
                regr_path = clim_diag_file_path(variables(ivar), field_types(ivar), imod, \
                                                diag_script, opt)

                ; do not process data if already processed

                if (isfilepresent(regr_path) .and. getenv("ESMValTool_force_calc").eq."False") then
                    pre_var = ncdf_read(regr_path, variables(ivar))
                    pre_var = pre_var^2
                    dnames = getVarDimNames(pre_var)
                else
                    ; extract model corresponding data
                    data = read_data(imod, variables(ivar), field_types(ivar))
                    ; square error values before averaging
                    var = extract_data(imod, data, -1, 0, 0)
                    var = var^2
                    delete(data)
                    ; determine dimensions to regrid
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

                    ; regrid in the horizontal (if necessary)

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

                    delete(var)
                    pre_var = var_hreg
                    delete(var_hreg)
                    outvar = sqrt(pre_var)
                    copy_VarMeta(pre_var, outvar)

                    ; output to NetCDF
                    outvar@var = variables(ivar)
                    outvar@diag_script = (/diag_script/)
                    outvar@ncdf = regr_path
                    ncdf_outfile = ncdf_write(outvar, regr_path)

                    delete(outvar)
                end if

                if (isdefined("tmp")) then
                    delete(tmp)
                end if

                tmp = time_operations(pre_var, -1, -1, "average", timemean, True)
                err = dim_avg_n_Wrap(tmp, dimsizes(dimsizes(tmp)) - 1)
                delete(tmp)
                err = sqrt(err)  ; pre_var is still the squared error (squared for averaging)!

                if (.not.isdefined("err_zm")) then
                    tmp = dimsizes(err)
                    dims = array_append_record(dim_MOD+1, tmp, 0)
                    delete(tmp)
                    err_zm = new(dims, float)
                    err_zm!0 = "model"
                    err_zm&model = new(dim_MOD+1, string)
                    err_zm&model = ""
                    if (numseas.eq.1) then
                        err_zm!1 = "lat"
                        err_zm&lat = err&lat
                    else
                        err_zm!1 = "season"
                        err_zm&season = err&season
                        err_zm!2 = "lat"
                        err_zm&lat = err&lat
                    end if
                    delete(dims)
                end if

                err_zm&model(imod) = models@name(imod)
                if (numseas.gt.1) then
                    err_zm(imod, :, :) = err(:, :)
                else
                    err_zm(imod, :) = err(:)
                end if

                delete(err)

            end if  ; if filepresent
        end do ; imod-loop
    end do ; ivar loop

    opt@time_avg = selected_grid

    ; ========================================================================
    ; =========================== calculations ===============================
    ; ========================================================================

    ; **********************************
    ; *** calculate multi-model mean ***
    ; **********************************

    mm_ind = 0  ; index of multi-model mean
    dim_MOD = dim_MOD + 1
    modelsatname = array_append_record((/"multi-model-mean"/), \
                                       models@name, 0)
    delete(models@name)
    models@name = modelsatname
    delete(modelsatname)

    ; update index of reference data set(s)

    ref_ind  = ind(models@name.eq.diag_script_info@ref_model(0))
    if (ref_ind2.ge.0) then
        ref_ind2 = ind(models@name.eq.diag_script_info@ref_model(1))
    end if

    ; debug output

;    system ("rm debug.nc")
;    debugfile = addfile("debug.nc","c")
;    debugfile->mask = global_mask

    ; read data

    do imod = 0, dim_MOD - 1

        if (imod.eq.0) then
            new_path = "none"
        else
            new_path = clim_diag_file_path(var0, field_type0, imod - 1, \
                                           diag_script, opt)
        end if

        if (isfilepresent_esmval(new_path)) then
            info_output("Read in: " + new_path, verbosity, 1)
            data_temp = ncdf_read(new_path, var0)
            data1 = time_operations(data_temp, -1, -1, "average", \
                                    timemean, True)

            if (numseas.gt.1) then
                do is = 0, numseas - 1
                    data1(is, :, :) = data1(is, :, :) + global_mask
                end do
            else
                data1 = data1 + global_mask
            end if

            if (.not.isdefined("data")) then
                dim_data = array_append_record((/dim_MOD/), dimsizes(data1), 0)
                data = new(dim_data, float)
                data!0 = "model"
                data&model = models@name
            end if
            rank = dimsizes(dim_data)

            ; debug output

;            debugfile->$models@name(imod)$ = data1

            if (numseas.gt.1) then
                data(imod, :, :, :) = data1
            else
                data(imod, :, :) = data1
            end if

            delete(data_temp)
            delete(data1)

        else
            if (imod.eq.mm_ind) then
            else
                ; skip this model as no data are available

                info_output("No data available for model:" + \
                            models@name(imod), verbosity, 1)
                info_output("Skipping", verbosity, 1)
            end if
        end if
    end do  ; imod

    if (ref_ind2.ge.0) then
        model_ind_woref = ind((ref_model.ne.models@name).and. \
                              (ref_model2.ne.models@name))
    else
        model_ind_woref = ind(ref_model.ne.models@name)
    end if

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

    ; special case: if there is only one model, skip averaging

    if (dimsizes(model_ind_woref).eq.1) then
        if (numseas.gt.1) then
            data(mm_ind, :, :, :) = data(model_ind_woref, :, :, :)
        else
            data(mm_ind, :, :) = data(model_ind_woref, :, :)
        end if
    else
        if (numseas.gt.1) then
            data(mm_ind, :, :, :) = dim_avg_n_Wrap(data(model_ind_woref, :, :, :), 0)
        else
            data(mm_ind, :, :) = dim_avg_n_Wrap(data(model_ind_woref, :, :), 0)
        end if
    end if

    data&model(mm_ind) = models@name(mm_ind)

    ; differences between multi-model mean and reference data set

    if (numseas.gt.1) then
        diff = data(mm_ind, :, :, :) - data(ref_ind, :, :, :)
        copy_VarMeta(data(ref_ind, :, :, :), diff)
    else
        diff = data(mm_ind, :, :) - data(ref_ind, :, :)
        copy_VarMeta(data(ref_ind, :, :), diff)
    end if

;    debugfile->diff = diff

    ; calculate zonal means

    tmp = dim_avg_n_Wrap(data, dimsizes(dimsizes(data)) - 1)

    ; we order the zonal mean array in the following way so that
    ; the lines for the multi-model mean and reference model will
    ; be drawn on top of the lines for the individual models:
    ;   (1) individual model(s)
    ;   (2) reference model(s) (observations)
    ;   (3) multi-model mean

    dims = dimsizes(tmp)
    zonalmean = new(dims, float)
    copy_VarMeta(tmp, zonalmean)

    i = 0
    n = dimsizes(model_ind_woref) - 1

    ; first entries in "zonalmean" = individual models

    do i = 0, n
        if (numseas.gt.1) then
            zonalmean(i, :, :) = tmp(model_ind_woref(i), :, :)
        else
            zonalmean(i, :) = tmp(model_ind_woref(i), :)
        end if
        zonalmean&model(i) = tmp&model(model_ind_woref(i))
    end do

    model_ind_rest = new((/dim_MOD - n - 1/), integer)   ; observations only

    ; get indices of observations

    k = 0

    do i = 0, dim_MOD - 1
        if (any(model_ind_woref.eq.i)) then
        else
            if (i.ne.mm_ind) then
                model_ind_rest(k) = i
                k = k + 1
            end if
        end if
    end do

    ; now add observations to "zonal mean"

    n = n + 1

    do i = 0, k - 1
        if (numseas.gt.1) then
            zonalmean(n, :, :) = tmp(model_ind_rest(i), :, :)
        else
            zonalmean(n, :) = tmp(model_ind_rest(i), :)
        end if
        zonalmean&model(n) = tmp&model(model_ind_rest(i))
        n = n + 1
    end do

    ; last entry in "zonalmean" = multi-model mean

    if (numseas.gt.1) then
        zonalmean(n, :, :) = tmp(mm_ind, :, :)   ; multi-model mean
    else
        zonalmean(n, :) = tmp(mm_ind, :)   ; multi-model mean
    end if

    zonalmean&model(n) = tmp&model(mm_ind)

    ; update indices for array "zonalmean"

    mm_ind = n

    delete(tmp)

    ; ========================================================================
    ; ============================= plotting =================================
    ; ========================================================================

    if (getenv("ESMValTool_write_plots").eq.True) then

        climofiles = new(dim_MOD - 1, string)  ; dim_MOD = number of models + 1 (i.e. models + multi-model-mean)

        do ii = 0, dim_MOD - 2  ; loop only over individual models (skip multi-model-mean)
            climofiles(ii) = \
                interface_get_inpaths(ii) + "/" \
                + interface_get_infile(variables(0), field_types(0), ii)
        end do

        if (dim_VAR.gt.1) then  ; is an error estimate present?
            climofilesx = array_append_record(climofiles, \
                interface_get_inpaths(0) + "/" + \
                interface_get_infile(variables(1), field_types(1), 0), 0)
            delete(climofiles)
            climofiles = climofilesx
            delete(climofilesx)
        end if

        alltags = array_append_record(tags, (/"PT_geo", "PT_zonal", "ST_clim", "DM_global"/), 0)

        ; create outfile directory

        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "clouds_ipcc")

        diff@res_gsnMaximize       = True         ; use full page for the plot
        diff@res_cnFillOn          = True         ; color plot desired
        diff@res_cnLineLabelsOn    = False        ; contour lines
        diff@res_cnLinesOn         = False
        diff@res_tiMainOn          = False
        diff@res_mpPerimOn         = False
        diff@res_mpGridLineColor   = -1
        diff@res_mpGridAndLimbOn   = True

;        diff@res_gsnLeftString           = annots(imod)
        diff@res_gsnLeftStringFontHeightF = 0.015
;        diff@res_gsnLeftString           = models@name(imod)
        diff@res_cnLevelSelectionMode     = "ExplicitLevels"
        diff@res_mpOutlineOn              = True
        diff@res_cnLevels                 = variable_info@map_diff_Levels
        diff@res_mpFillOn                 = False
;        diff@res_mpLandFillColor          = "Black"
;        pal = read_colormap_file("diag_scripts/lib/ncl/rgb/qcm3.rgb")
;        diff@res_cnFillColors             = pal
        diff@res_lbLabelBarOn             = True
        diff@res_gsnRightString           = ""
        diff@res_mpFillDrawOrder          = "PostDraw"     ; draw map fill last
        diff@res_cnMissingValFillColor    = "Gray"
        diff@res_tmYLLabelsOn             = False
        diff@res_tmYLOn                   = False
        diff@res_tmYRLabelsOn             = False
        diff@res_tmYROn                   = False
        diff@res_tmXBLabelsOn             = False
        diff@res_tmXBOn                   = False
        diff@res_tmXTLabelsOn             = False
        diff@res_tmXTOn                   = False
;        diff@res_gsnRightString           = variable_info@long_name + " in " \
;                                            + variable_info@units
        diff@res_cnInfoLabelOn            = False      ; turn off cn info label
        diff@res_mpProjection             = diag_script_info@projection
        if (isatt(diag_script_info, "explicit_cn_levels")) then
            diff@res_cnLevelSelectionMode = "ExplicitLevels"
            diff@res_cnLevels = diag_script_info@explicit_cn_levels
        end if

        diff@var = var0  ; Overwrite existing entry
        if (isatt(variable_info, "long_name")) then
            diff@var_long_name = variable_info@long_name
        end if
        if (isatt(variable_info, "units")) then
            diff@var_units     = variable_info@units
        end if

        if (var0.eq."clt") then
            if (isatt(diff, "res_cnLevels")) then
                delete(diff@res_cnLevels)
            end if
            diff@res_cnLevels = fspan(-25, 25, 11)
        end if

        if ((var0.eq."lwp").or.(var0.eq."clivi")) then
            if (isatt(diff, "res_cnLevels")) then
                delete(diff@res_cnLevels)
            end if
            diff@res_cnLevels = fspan(-100, 100, 11)
            diff = diff * 1000
            zonalmean = zonalmean * 1000
            diff@var_units = "g m-2"
            variable_info@units = diff@var_units
            zonalmean@units = diff@var_units
        end if

        if ((var0.eq."pr").or.(var0.eq."pr-mmday")) then
            diff@res_mpLandFillColor          = "Black"
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-precip-delta.rgb")
            diff@res_cnFillColors             = pal
        end if

        ; ###########################################
        ; # Optional output to netCDF               #
        ; ###########################################
        ; optional output (controlled by diag_script_info)

        if (isatt(diag_script_info, "ncdf")) then
            diff@ncdf = diag_script_info@ncdf

            ; function in ~/interface_scripts/auxiliary.ncl
            ncdf_outfile = ncdf_write(diff, diag_script_info@ncdf)
        end if

        plots = new((/2, numseas/), graphic)

        do is = 0, numseas - 1

            ; --------------------------------------------------------------------
            ; create workspace

            if (isvar("wks")) then
                delete(wks)
            end if

            wks = get_wks("dummy_for_wks", diag_script, var0 + "_" + season(is))
;            drawNDCGrid(wks) ; debugging option

            ; --------------------------------------------------------------------
            ; plot contour map

            diff@res_gsnDraw         = False  ; Do not draw yet
            diff@res_gsnFrame        = False  ; Don't advance frame.

            diff@res_lbTitleString   = "~F8~D~F21~" + variable_info@long_name + \
                                       " (" + variable_info@units + ")"
            diff@res_lbTitlePosition      = "Bottom"
            diff@res_lbLabelFontHeightF   = 0.015
            diff@res_lbAutoManage         = False
;            diff@res_lbTopMarginF         = 0.1
            diff@res_lbTitleFontHeightF   = 0.015

            ; plot contour map

            if (numseas.gt.1) then
                 plots(0, is) = contour_map(wks, diff(is, :, :), var0)
            else
                plots(0, 0) = contour_map(wks, diff, var0)
            end if

            ; --------------------------------------------------------------------
            ; plot zonal means

            lat = (/zonalmean&lat/)

            ; the attribute "long_name" will be used as title string for the x-axis

            lat@long_name = "Latitude"
            lat@units = "degrees_north"

            ; the attribute "long_name" will be used as title string for the y-axis

            zonalmean@long_name = zonalmean@long_name + " (" + \
                                  zonalmean@units + ")"

            linethickness          = new(dim_MOD, float)
            linecolor              = new(dim_MOD, string)
            linedash               = new(dim_MOD, integer)

            n = dimsizes(model_ind_woref) - 1

            ; settings for all models that have been used to calculate the
            ; multi-model mean (= first entries in "zonalmean")

            linethickness(0:n)      = 1.0     
            linecolor(0:n)          = "(/0.5, 0.5, 0.5/)"
            linedash                = 0
;            linedash(n+1:dim_MOD-1) = 0

            ; for all other models, styleset defined in config file is
            ; used (if defined), otherwise use default values

            if (isatt(diag_script_info, "styleset")) then
                colortab     = project_style(diag_script_info, "colors")
                dashtab      = project_style(diag_script_info, "dashes")
                thicknesstab = project_style(diag_script_info, "thicks")
                do i = n + 1, dim_MOD - 1
                    if (isdefined("idx")) then
                        delete(idx)
                    end if
                    idx = ind(models@name.eq.zonalmean&model(i))
                    if (all(ismissing(idx))) then
                    else
                        linecolor(i) = colortab(idx(0))
                        linethickness(i) = thicknesstab(idx(0))
                        linedash(i) = dashtab(idx(0))
                    end if
                end do
            else
                linethickness(n+1:dim_MOD-1) = 4.0     ; reference data set(s) (observations)
                linethickness(mm_ind)  = 4.0           ; multi-model mean
                linecolor(n+1:dim_MOD-1) = "Black"     ; reference data set
                linecolor(mm_ind)          = "Red"     ; multi-model mean
            end if

            res                    = True
            res@xyDashPattern      = linedash
            res@xyMonoLineColor    = False
            res@xyLineThicknesses  = linethickness
            res@xyLineColors       = linecolor
;            res@tiMainString       = "...title..."      ; add title
            res@tmXBLabelFontHeightF = 0.023
            res@tmYLLabelFontHeightF = 0.023
            res@gsnDraw            = False  ; Do not draw yet
            res@gsnFrame           = False  ; Don't advance frame.

            ; this controls the size and location of the second plot

            res@vpWidthF           = 0.6
            res@vpHeightF          = 0.475

            if (numseas.gt.1) then
                plots(1, is) = gsn_csm_xy(wks, lat, zonalmean(:, is, :), res)
                if (isdefined("err_zm")) then
                    y = new((/2, dimsizes(err_zm&lat)/), float)
                    res@gsnXYFillColors = "(/1.0, 0.8, 0.8/)"
                    res@xyLineColors(:) = "transparent"
                    do imod = 0, dim_MOD - 1
                        erridx = ind(err_zm&model.eq.zonalmean&model(imod))
                        if (.not.all(ismissing(erridx))) then
                            y(0, :) = zonalmean(imod, is, :) - err_zm(erridx, is, :)
                            y(1, :) = zonalmean(imod, is, :) + err_zm(erridx, is, :)
                            errorshade = gsn_csm_xy(wks, err_zm&lat, y, res)
                            overlay(plots(1, is), errorshade)
                        end if
                    end do
                end if
            else
                plots(1, 0) = gsn_csm_xy(wks, lat, zonalmean, res)

                if (isdefined("err_zm")) then
                    y = new((/2, dimsizes(err_zm&lat)/), float)
                    res@gsnXYFillColors = "(/1.0, 0.8, 0.8/)"
                    res@xyLineColors(:) = "transparent"
                    do imod = 0, dim_MOD - 1
                        erridx = ind(err_zm&model.eq.zonalmean&model(imod))
                        if (.not.all(ismissing(erridx))) then
                            y(0, :) = zonalmean(imod, :) - err_zm(erridx, :)
                            y(1, :) = zonalmean(imod, :) + err_zm(erridx, :)
                            errorshade = gsn_csm_xy(wks, err_zm&lat, y, res)
                            overlay(plots(1, 0), errorshade)
                        end if
                    end do
                end if
            end if

            ; --------------------------------------------------------------------

            pres                    = True
            pres@gsnPanelCenter     = False
            pres@gsnPanelXF         = (/0.075, 0.625/)  ; hor. pos. of sub-plots
            pres@txString = season(is)

            outfile = panelling(wks, plots(:, is), 1, 2, pres)
            info_output(" Wrote " + wks@fullname, verbosity, 1)

            ; add meta data to plot (for reporting)

            caption = "Multi model mean bias (left) and zonal averages (right) for variable " + variables(0) \
                      + " (" + season(is) + "), reference = " + models@name(ref_ind) + "."
            id = diag_script + "_" + variables(0) + "_" + season(is)

            ESMValMD(wks@fullname, alltags, caption, id, variables(0), models@name, climofiles)

        end do ; is-loop (seasons)
    end if ; if create plots

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
