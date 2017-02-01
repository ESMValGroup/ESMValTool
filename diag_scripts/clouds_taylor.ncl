;;#############################################################################
;; CLOUDS_TAYLOR
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME EMBRACE
;;#############################################################################
;; Description
;;    Calculates the performance of models in reproducing 2-d fields of annual
;;    mean cloud properties (Taylor diagramms). The code is based on
;;       1) perfmetrics_main.ncl
;;       2) perfmetrics_taylor.ncl
;;       3) perfmetrics_taylor_collect.ncl
;;    originally written by Franziska Frank (DLR, Germany).
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref)
;;    ref_model:       name of reference data set
;;    valid_fraction:  minimum fraction of valid data points
;;                     (for creating a mask)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    calcmm:          include calculating of the multi-model mean
;;    ncdf:            enable to output to netCDF; either use "default" or give
;;                     a full file name
;;    styleset:        "CMIP5", "DEFAULT"
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;    none
;;
;; Modification history
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
    if (.not.(isatt(diag_script_info, "ref_model"))) then
        error_msg("f", diag_script, "", "no reference model is selected")
    end if

    ; set reference model
    ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
    ref_model = diag_script_info@ref_model
    if (ismissing(ref_ind)) then
        error_msg("f", diag_script, "", "no adequate reference model provided")
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

    ; defined the target grid for regridding

    models_ind = ispan(0, dim_MOD - 1, 1)
    info_output("Setting target grid:", verbosity, 2)
    if (selected_grid.eq."ref") then
        grid = find_destination_grid(ref_ind, \
                                     var0, field_type0, "coarsest")
    else
        grid = find_destination_grid(models_ind, \
                                     var0, field_type0, selected_grid)
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

        ;; Do not process data if already processed
        if (isfilepresent(regr_path) .and. \
            getenv("ESMValTool_force_calc").eq."False") then

;            info_output("Reading in file = " + regr_path, verbosity, 1)
            pre_var = ncdf_read(regr_path, var0)
            dnames = getVarDimNames(pre_var)

        else

            ;; Extract model corresponding data
            data = read_data(imod, var0, field_type0)
            var = extract_data(imod, data, -1, 0, 0)
            delete(data)

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
            delete(var)

            pre_var = var_hreg
            delete(var_hreg)

            ;; Output to NetCDF
            pre_var@var = var0
            pre_var@diag_script = (/diag_script/)
            pre_var@ncdf = regr_path
            ncdf_outfile = ncdf_write(pre_var, regr_path)

        end if

        ;; Create a missing value mask for each model
        ;; The dim_avg_wgt function with option limit is used on the time
        ;; coordinate: in this way each grid box containing more than
        ;; valid_fraction of missing values along the time series is
        ;; set to missing.

        if (imod.eq.0) then
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

        ;; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1
        if (dimsizes(dimsizes(pre_var)).eq.1) then
            tmp = dim_avg_wgt_Wrap(pre_var, ww, limit)
        else
            tmp = dim_avg_wgt_n_Wrap(pre_var, ww, limit, 0)
        end if

        tmp = where(ismissing(tmp), tmp, 0.)
        if (all(ismissing(tmp))) then
            error_msg("w", diag_script, "", models@name(imod) + \
                      "contains only values, Excluding from global " + \
                      "masking.")
        else
            global_mask = global_mask + tmp
        end if

        delete(tmp)
        delete(ww)
        delete(pre_var)

    end do  ; imod

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
            modelsatname = array_append_record((/"multi-model mean"/), \
                                               models@name, 0)
            delete(models@name)
            models@name = modelsatname
            delete(modelsatname)
            ; update index of reference data set
            ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
        end if
    end if

    ; debug output
;    system ("rm debug.nc")
;    debugfile = addfile("debug.nc","c")
;    debugfile->mask = global_mask

    ; read data

    do imod = 0, dim_MOD - 1
        if (flag_mean) then
            if (imod.eq.0) then
                new_path = "none"
            else
                new_path = clim_diag_file_path(var0, field_type0, imod - 1, \
                                               diag_script, opt)
            end if
        else
            new_path = clim_diag_file_path(var0, field_type0, imod, \
                                           diag_script, opt)

        end if

        if (isfilepresent_esmval(new_path)) then
            info_output("Read in: " + new_path, verbosity, 1)
            data_temp = ncdf_read(new_path, var0)
            data1 = time_operations(data_temp, -1, -1, "average", \
                                    timemean, True)
            if (numseas.eq.1) then
                data1 = data1 + global_mask
            else
                do is = 0, numseas - 1
                    data1(is, :, :) = data1(is, :, :) + global_mask
                end do
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

            if (numseas.eq.1) then
                data(imod, :, :) = data1
            else
                data(imod, :, :, :) = data1
            end if

            delete(data_temp)
            delete(data1)

        else
            if ((flag_mean).and.(imod.eq.0)) then
            else
                ; skip this model is no data are available

                info_output("No data available for model:" + \
                            models@name(imod), verbosity, 1)
                info_output("Skipping", verbosity, 1)
            end if
        end if
    end do  ; imod loop

    ; number of used reference data sets

    dim_REF = dimsizes(diag_script_info@ref_model)

    ; define result variable

    val = new((/dim_MOD - dim_REF, 2, numseas/), float)
;    val!0 = "diagnostics"
    val!0 = "models"
    val!1 = "statistic"
    val!2 = "time"
;    val&diagnostics = var0_new
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

            if (flag_mean.and.imod.eq.0) then
                model_ind_woref = ind(ref_model.ne.data&model)
                if (numseas.eq.1) then
                    var = dim_avg_n_Wrap(data(model_ind_woref, :, :), 0)
                else
                    var = dim_avg_n_Wrap(data(model_ind_woref, is, :, :), 0)
                end if
            else  ; Given model
                if (numseas.eq.1) then
                    var = data(imod, :, :)
                else
                    var = data(imod, is, :, :)
                end if
            end if

            ; calculate rate of standard deviations

            val(mod_i, 0, is) = calculate_metric(var, ref, "stddev_ratio_nowgt")

            ; calculate pattern correlation

            val(mod_i, 1, is) = calculate_metric(var, ref, "correlation_nowgt")

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

    ; ========================================================================
    ; ============================= plotting =================================
    ; ========================================================================

    if (getenv("ESMValTool_write_plots").eq.True) then
        ; create outfile directory
        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "clouds_taylor")

        ;; create workspace

        wks = get_wks("dummy_for_wks", diag_script, var0)

        ratio = new((/mod_i, 1, numseas/), float)  ; standard deviation ratios
        ratio = val(:, 0, :)
        cc = new((/mod_i, 1, numseas/), float)    ; correlations
        cc = val(:, 1, :)

        colortab = (/"(/0.00, 0.00, 0.59/)", "(/0.00, 0.39, 1.00/)", \
                     "(/0.20, 1.00, 1.00/)", "(/0.20, 0.88, 0.00/)", \
                     "(/1.00, 0.88, 0.00/)", "(/1.00, 0.59, 0.00/)", \
                     "(/1.00, 0.20, 0.00/)", "(/0.59, 0.00, 0.00/)", \
                     "(/0.78, 0.00, 0.78/)", "(/0.59, 0.00, 0.59/)", \
                     "(/0.90, 0.90, 0.90/)", "(/0.70, 0.70, 0.70/)", \
                     "(/0.50, 0.50, 0.50/)", "(/0.30, 0.30, 0.30/)"/)

        markertab = (/16, 4, 5, 0/)

        colors  = new(mod_i, string)
        markers = new(mod_i, integer)

        i = 0
        idx = 0

        do while (i.lt.mod_i)
            if ((flag_mean).and.(i.eq.0)) then

                ; create new marker: filled star

                mstring = "z"
                fontnum = 35
                size    = 1.5
                angle   = 0.0

                new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, \
                                         1.0, size, angle)

                colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
                markers(i) = new_index
                i = i + 1
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
        end do

        ropts = True
        ropts@gsnDraw         = False
        ropts@gsnFrame        = False
        ropts@Colors          = colors
        ropts@Markers         = markers  ; marker styles
        ropts@gsMarkerSizeF   = 0.0125

        modelnames            = val&models
        ropts@caseLabels      = modelnames

;        ropts@stnRad          = (/ 0.5, 1.5 /)
        ropts@stnRad          = (/ 0.5, 1.5, 2.0, 2.5 /)
        ropts@OneX            = "1.00"
;        ropts@ccRays          = (/ 0.6, 0.9 /)
        ropts@ccRays          = (/ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, \
                                   0.9, 0.95, 0.99 /)
;        ropts@ccRays_color    = "LightGray"
;        ropts@ccRays_color    = "PaleGreen1"
        ropts@ccRays_color    = "Black"
        ropts@centerDiffRMS   = True
        ropts@centerDiffRMS_color = "LightGray"  ; "PaleGreen1"

        ropts@printDataLabels = False  ; print numeric label for each data pt.
        ropts@reverseOrder    = True  ; reverse order of items in legend
        ropts@legendBox       = True  ; draw box around legend

        ; legend parameters

;        if (dim_MOD.ge.20) then
;            ropts@caseLabelsFontHeightF = 0.06
;            ropts@plotSize              = 0.45
;        else
            ropts@caseLabelsFontHeightF = 0.1
            ropts@plotSize              = 0.6
;        end if

        ropts@legendXpos            = 0.575
        ropts@legendYpos            = -0.5
        ropts@legendWidth           = 0.2

        ; plot

        plots = new(numseas, graphic)  ; collect individual maps in a graphic array

        do is = 0, numseas - 1
            ropts@tiMainString    = var0 + " (" + season(is) + ")"  ; title

            plots(is) = taylor_diagram(wks, ratio(:, :, is), cc(:, :, is), ropts)

            pres = True      ; needed to override panelling defaults
            pres@gsnPanelLabelBar = False  ; no general label bar desired here
        end do

;;        outfile = panelling(wks, plots, 2, 2, pres)

        info_output(" Wrote " + wks@fullname, verbosity, 1)
    end if

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
