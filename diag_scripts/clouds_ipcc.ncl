;;#############################################################################
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
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref)
;;    ref_model:       name of reference data set
;;    valid_fraction:  minimum fraction of valid data points
;;                     (for creating a mask)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    colormap:        e.g., WhiteBlueGreenYellowRed, rainbow
;;    ncdf:            enable to output to netCDF; either use "default" or give
;;                     a full file name
;;    projection:      map projection, e.g., Mollweide, Mercator
;;    styleset:        "CMIP5", "DEFAULT"
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
;;    none
;;
;; Modification history
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
    dim_VAR = 1
    var0 = variables(0)
    field_type0 = field_types(0)

    ; time averaging: at the moment, only "annualclim" and "seasonalclim"
    ; are supported

    timemean = "annualclim"  ; default
    numseas = 1              ; default
    season = (/"DJF", "MAM", "JJA", "SON"/)

    if (isatt(diag_script_info, "timemean")) then
        timemean = diag_script_info@timemean
        if (timemean.eq."seasonalclim") then
            numseas = 4
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

        ; do not process data if already processed

        if (isfilepresent(regr_path) .and. \
            getenv("ESMValTool_force_calc").eq."False") then

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

        ; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1

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
    ; =========================== calculations ===============================
    ; ========================================================================

    ; **********************************
    ; *** calculate multi-model mean ***
    ; **********************************

    mm_ind = 0  ; index of multi-model mean
    dim_MOD = dim_MOD + 1
    modelsatname = array_append_record((/"multi-model mean"/), \
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
        model_ind_woref = ind((ref_model.ne.data&model).and. \
                              (ref_model2.ne.data&model))
    else
        model_ind_woref = ind(ref_model.ne.data&model)
    end if

    if (numseas.gt.1) then
        data(mm_ind, :, :, :) = dim_avg_n_Wrap(data(model_ind_woref, :, :, :), 0)
    else
        data(mm_ind, :, :) = dim_avg_n_Wrap(data(model_ind_woref, :, :), 0)
    end if

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
    ;   (1) individual models
    ;   (2) reference model (observations)
    ;   (3) multi-model mean

    dims = dimsizes(tmp)
    zonalmean = new(dims, float)
    copy_VarMeta(tmp, zonalmean)

    i = 0

    do n = 0, dims(0) - 1  ; loop over models
        if ((n.ne.ref_ind).and.(n.ne.ref_ind2).and.(n.ne.mm_ind)) then
            if (numseas.gt.1) then
                zonalmean(i, :, :) = tmp(n, :, :)
            else
                zonalmean(i, :) = tmp(n, :)
            end if
            i = i + 1
        end if
    end do

    if (numseas.gt.1) then
        if (ref_ind2.ge.0) then
            zonalmean(i, :, :) = tmp(ref_ind2, :, :)  ; alternative reference data set
            i = i + 1
        end if
        zonalmean(i, :, :)     = tmp(ref_ind, :, :)   ; observations
        zonalmean(i + 1, :, :) = tmp(mm_ind, :, :)    ; multi-model mean
    else
        if (ref_ind2.ge.0) then
            zonalmean(i, :)     = tmp(ref_ind2, :)    ; alternative reference data set
            i = i + 1
        end if
        zonalmean(i, :)     = tmp(ref_ind, :)  ; observations
        zonalmean(i + 1, :) = tmp(mm_ind, :)   ; multi-model mean
    end if

    ; update indices for array "zonalmean"

    if (ref_ind2.ge.0) then
        ref_ind2 = i - 1
    end if
    ref_ind  = i
    mm_ind   = i + 1

    delete(tmp)

    ; ========================================================================
    ; ============================= plotting =================================
    ; ========================================================================

    if (getenv("ESMValTool_write_plots").eq.True) then
        ; create outfile directory

        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "clouds_ipcc")

        diff@res_gsnMaximize       = True         ; use full page for the plot
        diff@res_cnFillOn          = True         ; color plot desired
        diff@res_cnLineLabelsOn    = False        ; contour lines
        diff@res_cnLinesOn         = False
        diff@res_tiMainOn          = False
        diff@res_mpPerimOn         = False
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

        ; --------------------------------------------------------------------
        ; create workspace

        wks = get_wks("dummy_for_wks", diag_script, var0)
;        drawNDCGrid(wks) ; debugging option

        plots = new((/2, numseas/), graphic)

        ; --------------------------------------------------------------------
        ; plot contour map

        diff@res_gsnDraw         = False  ; Do not draw yet
        diff@res_gsnFrame        = False  ; Don't advance frame.

        diff@res_lbTitleString   = "~F8~D~F21~" + variable_info@long_name + \
                                   " (" + variable_info@units + ")"
        diff@res_lbTitlePosition      = "Bottom"
;        diff@res_lbLabelFontHeightF   = 0.05
        diff@res_lbAutoManage         = True
        diff@res_lbTopMarginF         = 0.1
;        diff@res_lbTitleFontHeightF   = 0.05
;        diff@res_pmLabelBarWidthF     = 0.8

        ; plot contour map

        if (numseas.gt.1) then
            do is = 0, numseas - 1
                plots(0, is) = contour_map(wks, diff(is, :, :), var0)
            end do
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

        linethickness          = new(dims(0), float)
        linethickness          = 1.0     ; default
        linethickness(ref_ind) = 4.0     ; reference data set (observations)
        if (ref_ind2.ge.0) then
            linethickness(ref_ind2)  = 4.0  ; alternative reference data set
        end if
        linethickness(mm_ind)  = 4.0     ; multi-model mean

        linecolor              = new(dim_MOD, string)
        linecolor              = "LightGray"  ; default
        if (ref_ind2.ge.0) then
            linecolor(ref_ind2)    = "RoyalBlue"  ; alternative reference data set
        end if
        linecolor(ref_ind)     = "Black"      ; reference data set
        linecolor(mm_ind)      = "Red"        ; multi-model mean

        res                    = True
        res@xyDashPattern      = 0               ; Make curves all solid
        res@xyMonoLineColor    = False
        res@xyLineThicknesses  = linethickness
        res@xyLineColors       = linecolor
;        res@tiMainString       = "...title..."      ; add title
        res@tmXBLabelFontHeightF = 0.023
        res@tmYLLabelFontHeightF = 0.023
        res@gsnDraw            = False  ; Do not draw yet
        res@gsnFrame           = False  ; Don't advance frame.

        ; this controls the size and location of the second plot

        res@vpWidthF           = 0.6
        res@vpHeightF          = 0.475

        if (numseas.gt.1) then
            do is = 0, numseas - 1
                plots(1, is) = gsn_csm_xy(wks, lat, zonalmean(:, is, :), res)
            end do
        else
            plots(1, 0) = gsn_csm_xy(wks, lat, zonalmean, res)
        end if

        ; --------------------------------------------------------------------

        pres                   = True
        pres@gsnPanelCenter    = False
        pres@gsnPanelXF        = (/0.075, 0.625/)  ; hor. pos. of sub-plots

        if (numseas.gt.1) then
            do is = 0, numseas - 1
                pres@txString = season(is)
                outfile = panelling(wks, plots(:, is), 1, 2, pres)
            end do
        else
            outfile = panelling(wks, plots(:, 0), 1, 2, pres)
        end if

        info_output(" Wrote " + wks@fullname, verbosity, 1)
    end if

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
