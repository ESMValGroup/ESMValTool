;; SEASONALITY_MM
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME CRESCENDO
;;#############################################################################
;; Description
;;    Calculates the multi-model mean seasonality (DJF - JJA), multi-model mean
;;    of absolute seasonality, multi-model mean bias in seasonality, and multi-
;;    model mean bias in absolute seasonality (similar to IPCC ch. 9 fig. 9.3).
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref, 1x1)
;;    ref_model:       name of reference data set
;;    valid_fraction:  minimum fraction of valid data points
;;                     (for creating a mask)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    ncdf:            enable to output to netCDF; either use "default" or give
;;                     a full file name
;;    projection:      map projection, e.g., Mollweide, Mercator
;;
;; Required variable_info attributes (variable specific)
;;    map_diff_Levels: contour levels for difference plot
;;    units:           variable units
;;
;; Optional variable_info attributes (variable specific)
;;    long_name:       description of variable
;;
;; Caveats
;;    none
;;
;; Modification history
;;    20170622-A_laue_ax: added tags for reporting
;;    20170320-A_laue_ax: written.
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
                     (/"P_crescendo"/))   ; proj_references

    dim_MOD = dimsizes(models@name)
    dim_VAR = 1
    var0 = variables(0)
    field_type0 = field_types(0)

    season = (/"DJF", "JJA"/)
    numseas = 2

    climofiles = new(dim_MOD, string)

    do ii = 0, dim_MOD - 1
        climofiles(ii) = \
            interface_get_inpaths(ii) + "/" \
            + interface_get_infile(variables(0), field_types(0), ii)
    end do

    alltags = array_append_record(tags, (/"PT_geo", "ST_clim", "ST_mean", "ST_diff", "DM_global"/), 0)
    allmodelnames = array_append_record(models@name, (/"multi-model-mean"/), 0)
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

    ; check for required diag_script_info and variable_info attributes

    exit_if_missing_atts(diag_script_info, (/"grid", "ref_model", "valid_fraction"/))

    reqatts = (/"map_diff_Levels", "units"/)
    if (.not.all(isatt(variable_info, reqatts))) then
        error_msg("f", diag_script, "", "missing required " + \
                  "variable_info attribute: " + \
                  reqatts(ind(.not.isatt(variable_info, reqatts))))
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

    ; update index of reference data set

    ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))

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
            if (.not.isdefined("data")) then
                dim_data = new((/4/), integer)
                dim_data(0) = dim_MOD
                dim_data(1) = numseas
                dim_data(2:3) = dimsizes(global_mask)
                data = new(dim_data, float)
                copy_VarMeta(global_mask, data(0, 0, :, :))
                data!0 = "model"
                data&model = models@name
                data!1 = "season"
                data&season = season
                delete(dim_data)
            end if

            data_temp = ncdf_read(new_path, var0)
            do is = 0, numseas - 1
                data1 = time_operations(data_temp, -1, -1, "average", \
                                    season(is), True)
                data(imod, is, :, :) = data1 + global_mask
            end do

            delete(data_temp)
            delete(data1)

        else
            if (imod.eq.mm_ind) then
            else
                ; skip this model is no data are available

                info_output("No data available for model:" + \
                            models@name(imod), verbosity, 1)
                info_output("Skipping", verbosity, 1)
            end if
        end if
    end do  ; imod

    model_ind_woref = ind(ref_model.ne.data&model)
    data(mm_ind, :, :, :) = dim_avg_n_Wrap(data(model_ind_woref, :, :, :), 0)

    ; seasonality = DJF - JJA

    seas_ref = data(ref_ind, 0, :, :) - data(ref_ind, 1, :, :)
    seas_mm = data(mm_ind, 0, :, :) - data(mm_ind, 1, :, :)
    seas_mm@diag_script = (/diag_script/)
    copy_VarMeta(data(ref_ind, 0, :, :), seas_ref)
    copy_VarMeta(seas_ref, seas_mm)
    delete(data)

    abs_seas = abs(seas_mm)
    bias_seas = seas_mm - seas_ref
    bias_abs_seas = abs_seas - abs(seas_ref)

    ; ========================================================================
    ; ============================= plotting =================================
    ; ========================================================================

    if (getenv("ESMValTool_write_plots").eq.True) then
        ; create outfile directory

        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "clouds_ipcc")

        seas_mm@res_gsnMaximize       = True         ; use full page for the plot
        seas_mm@res_cnFillOn          = True         ; color plot desired
        seas_mm@res_cnLineLabelsOn    = False        ; contour lines
        seas_mm@res_cnLinesOn         = False
        seas_mm@res_tiMainOn                 = True
        seas_mm@res_gsnLeftStringFontHeightF = 0.015
        seas_mm@res_cnLevelSelectionMode     = "ExplicitLevels"
        seas_mm@res_mpOutlineOn              = True
        seas_mm@res_mpFillOn                 = False
        seas_mm@res_lbLabelBarOn             = True
        seas_mm@res_gsnRightString           = ""
        seas_mm@res_mpFillDrawOrder          = "PostDraw"     ; draw map fill last
        seas_mm@res_cnMissingValFillColor    = "Gray"
        seas_mm@res_tmYLLabelsOn             = False
        seas_mm@res_tmYLOn                   = False
        seas_mm@res_tmYRLabelsOn             = False
        seas_mm@res_tmYROn                   = False
        seas_mm@res_tmXBLabelsOn             = False
        seas_mm@res_tmXBOn                   = False
        seas_mm@res_tmXTLabelsOn             = False
        seas_mm@res_tmXTOn                   = False
        seas_mm@res_cnInfoLabelOn            = False      ; turn off cn info label
        seas_mm@res_mpProjection             = diag_script_info@projection

        seas_mm@var = var0  ; Overwrite existing entry
        if (isatt(variable_info, "long_name")) then
            seas_mm@var_long_name = variable_info@long_name
        end if
        if (isatt(variable_info, "units")) then
            seas_mm@var_units     = variable_info@units
        end if

        ; ###########################################
        ; # Optional output to netCDF               #
        ; ###########################################
        ; optional output (controlled by diag_script_info)

        if (isatt(diag_script_info, "ncdf")) then
            seas_mm@ncdf = diag_script_info@ncdf
            ; function in ~/interface_scripts/auxiliary.ncl
            ncdf_outfile = ncdf_write(seas_mm, diag_script_info@ncdf)
        end if

        ; --------------------------------------------------------------------
        ; create workspace

        wks = get_wks("dummy_for_wks", diag_script, var0)
;        drawNDCGrid(wks) ; debugging option

        plots = new((/4/), graphic)

        ; --------------------------------------------------------------------
        ; plot contour map

        seas_mm@res_gsnDraw         = False  ; Do not draw yet
        seas_mm@res_gsnFrame        = False  ; Don't advance frame.
        seas_mm@res_mpPerimOn       = False

        seas_mm@res_lbTitleString   = "(" + variable_info@units + ")"
        seas_mm@res_lbTitlePosition      = "Bottom"

        seas_mm@res_lbLabelFontHeightF   = 0.014
        seas_mm@res_lbTopMarginF         = 0.1
        seas_mm@res_lbTitleFontHeightF   = 0.014

        seas_mm@res_tiMainFontHeightF    = 0.016

        seas_mm@res_tiMainString = "Multi Model Mean Bias"

        if (isatt(seas_mm, "res_cnLevels")) then
            delete(seas_mm@res_cnLevels)
        end if
        seas_mm@res_cnLevels = variable_info@map_diff_Levels

        copy_VarMeta(seas_mm, abs_seas)
        copy_VarMeta(seas_mm, bias_seas)
        copy_VarMeta(seas_mm, bias_abs_seas)

        if ((var0.eq."tas").or.(var0.eq."tas-degC")) then
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-tas.rgb")
            seas_mm@res_cnFillColors = pal
            pal2 = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-tas-absdelta.rgb")
            abs_seas@res_cnFillColors = pal2
            pal3 = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-tas-seasdelta.rgb")
            bias_seas@res_cnFillColors = pal3
            bias_abs_seas@res_cnFillColors = pal3

            delete(seas_mm@res_cnLevels)
            delete(abs_seas@res_cnLevels)
            delete(bias_seas@res_cnLevels)
            delete(bias_abs_seas@res_cnLevels)

            ; IPCC ch.9 fig. 9.3 scales for tas
            seas_mm@res_cnLevels = fspan(-50, 50, 21)
            abs_seas@res_cnLevels = fspan(5, 50, 10)
            bias_seas@res_cnLevels = fspan(-8, 8, 17)
            bias_abs_seas@res_cnLevels = fspan(-8, 8, 17)
        end if

        seas_mm@res_tiMainString = "Multi Model Mean"
        abs_seas@res_tiMainString = "Multi Model Mean of Absolute Seasonality"
        bias_seas@res_tiMainString = "Multi Model Mean Bias in Seasonality"
        bias_abs_seas@res_tiMainString = "Multi Model Mean Bias in Absolute Seasonality"

        ; panelling resources

        pres                   = True
        pres@gsnPanelCenter    = False
        pres@gsnPanelRowSpec = True                   ; tell panel what order to plt
        pres@gsnPanelYWhiteSpacePercent = 5
        pres@gsnPanelXWhiteSpacePercent = 5

        plots(0) = contour_map(wks, seas_mm, var0)
        plots(1) = contour_map(wks, abs_seas, var0)
        plots(2) = contour_map(wks, bias_seas, var0)
        plots(3) = contour_map(wks, bias_abs_seas, var0)

        gsn_panel(wks, plots, (/2, 2/), pres)

        ; add meta data to plot (for reporting)

        caption = "Multi model values for variable " + variables(0) \
                  + ", from top left to bottom right: mean, mean of absolute seasonality, " \
                  + "mean bias in seasonality, mean bias in absolute seasonality. Reference dataset: " \
                  + models@name(ref_ind) + ". Similar to IPCC AR5, fig. 9.3."

        id = diag_script + "_" + variables(0)

        contrib_authors = "A_laue_ax"

        ESMValMD(wks@fullname, alltags, caption, id, variables(0), allmodelnames, climofiles, diag_script, contrib_authors)

        info_output(" Wrote " + wks@fullname, verbosity, 1)
    end if

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
