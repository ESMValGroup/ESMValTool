;##############################################################################
;# CONTOUR PLOTS                                                              #
;##############################################################################
;# Author: Klaus-Dirk Gottschaldt (DLR, Germany)                              #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    *                                                                       #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    *                                                                       #
;# Required variable_info attributes (variable specific):                     #
;#    *                                                                       #
;# Caveats:                                                                   #
;# Modification history:                                                      #
;#    * 20151203-A_laue_ax added option to calculate/plot multi-model-mean    #
;#                         (diag_script_info@calcmm) and differences to       #
;#                         reference data set (diag_script_info@calcdiff)     #
;#    * 20151027-A_laue_ax moved call to 'write_references' to the beginning  #
;#                         of the code                                        #
;#    * 20150325-A_laue_ax modified reference tags used for acknowledgements  #
;#                         (projects, observations, etc)                      #
;#    * 20140528 A_gott_kl written based on ancyc_lat.ncl                     #
;#                                                                            #
;##############################################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"

;; fetch general parameters
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; References Information
    write_references(diag_script,      \  ; script name
                     "A_gott_kl",      \  ; authors
                     "",               \  ; contributors
                     "D_0000",         \  ; diag_references
                     "",               \  ; obs_references
                     (/"P_esmval"/))          ; proj_references

; A-laue_ax+
    flag_diff = False
    flag_mean = False

    if (isatt(diag_script_info, "calcmm")) then
        if (diag_script_info@calcmm) then
            flag_mean = True
        end if
    end if

    if (isatt(diag_script_info, "calcdiff")) then
        if (diag_script_info@calcdiff) then
            info_output("diag_script_info@calcdiff = True requires a reference model", verbosity, 1)
            ref_model = var_attr_ref(0)  ; produces an error if ref_model is not set in the namelist
            info_output("reference model = " + ref_model, verbosity, 1)
            idx_ref = get_ref_model_idx(models, var_attr_ref(0))
            flag_diff = True
        end if
    end if
; A-laue_ax-

    do idx_mod = 0, dim_MOD - 1
        data1 = read_data(idx_mod, var0, field_type0)
    end do
end

;; Load variable_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end


;; main code
begin
    ;; Check valid plot fields
    plot_fields = (/"T2Ms", "T2Ds", "C2Ms", "C2Ds"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_type0 + \
              " for this diag_script")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; determine destination grid (regridding.ncl)
    idx = ispan(0, dim_MOD - 1, 1)
    dest_grd = find_destination_grid(idx, var0, field_type0, "coarsest")

    ;; default time range & annots
    y1 = toint(models@start_year)
    y2 = toint(models@end_year)
    annots = project_style(diag_script_info, "annots")

    ;; Loop over models to regrid to a reference grid and create
    ;;     mask for missing values
    do idx_mod = 0, dim_MOD - 1

        ;; Extract data
        A0 = extract_data_extended(idx_mod, var0, data1, -1, 0, 0)
        ;; Use overlapping time period or as given in nml?
        if(isatt(diag_script_info, "range_option")) then
            y1(idx_mod) = get_start_year(idx_mod)
            y2(idx_mod) = get_end_year(idx_mod)
            temp = time_operations(A0, y1(idx_mod), y2(idx_mod), \
                                   "extract", "", 0)
            delete(A0)
            A0 = temp
            delete(temp)
        end if  ; Default: use as in nml
        ;; annual cycle (Note: Returns time dimension "month")
        A1 = time_operations(A0, -1, -1, "average", "annualclim", True)

        ;; initialize target fields ...
        if(idx_mod.eq.0) then
            ;; ... for regridded data
            if(.not.isatt(A1, "_FillValue")) then
                Fill = default_fillvalue(typeof(A1))
            else
                Fill = A1@_FillValue
            end if
            var_grd = new((/dim_MOD, dimsizes(dest_grd&lat),\
                          dimsizes(dest_grd&lon)/), typeof(A1), Fill)
            var_grd!0 = "models"
            var_grd!1 = "lat"
            var_grd!2 = "lon"
            var_grd&models = annots
            var_grd&lat = dest_grd&lat
            var_grd&lon = dest_grd&lon

            ;; ... for the mask
            mask0 = var_grd(0, :, :)
            mask0(:, :) = 1
        end if

        ;; regrid to ref
        A1 = where(ismissing(A1), Fill, A1)
        A1@_FillValue = Fill  ; ensure consistent _FillValue
        var_grd(idx_mod, :, :) = area_hi2lores_Wrap(A1&lon, A1&lat, A1, \
            True, 1, dest_grd&lon, dest_grd&lat, False)

        ;; Build the mask: missing in at least one model = missing in the mask
        mask0 = where(ismissing(var_grd(idx_mod, :, :)), mask0@_FillValue, \
                      mask0)

        ;; clean up
        delete([/A0, A1/])
    end do  ; idx_MOD, for regridding and mask

    info_output(diag_script + ": Mean value calculated after weighted " + \
        "averaging over time, regridding (lon, lat) to coarsest grid, " + \
        "masking all points that are missing in at least one model, " + \
        "considering area weighting", verbosity, 1)

; A-laue_ax+
    if (flag_mean) then
        ; add entries for "multi-model mean"

        dim_MOD = dim_MOD + 1
        modelsatname = array_append_record(annots, (/"multi-model mean"/), 0)
        delete(annots)
        annots = modelsatname
        delete(modelsatname)

        ; calculate multi-model-mean
        ; if a reference model is given (ref_model), this model will not be
        ; included in the multi-model-mean

        avgidx = ind(ispan(0, dim_MOD-2, 1).ne.idx_ref)
        mmm = dim_avg_n_Wrap(var_grd(avgidx, :, :), 0)
        mmm = mmm * mask0

        tmp = array_append_record(var_grd, conform(var_grd,mmm, (/1,2/)), 0)
        delete(var_grd)
        delete(mmm)
        var_grd = tmp
        delete(tmp)
    end if
; A-laue_ax-

    ;; Output dir
    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; General settings for this plot
    res = True
    res@cnFillOn = True
    res@cnLinesOn = False
    ;res@cnFillPalette =
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = variable_info@acp_cnlvl
    res@mpProjection = diag_script_info@mpProjection
    res@txFont = diag_script_info@font
    res@tmXBLabelFont = diag_script_info@font
    res@tmYLLabelFont = diag_script_info@font
    ;res@tmXBMode = "Explicit"
    ;res@tmXBValues = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    ;res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", \
                      ;"N", "D"/)
    if(isatt(diag_script_info, "lbLabelBarOn")) then
        res@lbLabelBarOn = diag_script_info@lbLabelBarOn
    else
        res@lbLabelBarOn = True
    end if
    res@lbLabelFont = diag_script_info@font
    res@lbTitleFont = diag_script_info@font
    res@lbTitleFontHeightF = 0.01  ; default: 0.025
    res@lbTitleOffsetF = -0.25  ; default: 0.03 (negative values move it down)
    res@lbTopMarginF = 0.1  ; default: 0.05 (negative values move it up)
    res@lbTitleJust = "BottomCenter"  ; default: CenterCenter
    ;res@lbTitleString = variable_info@long_name
    if(isatt(diag_script_info, "gsnZonalMean")) then
        res@gsnZonalMean = diag_script_info@gsnZonalMean
    end if

; A-laue_ax+
    if (flag_diff) then
        plot_types = 2
    else
        plot_types = 1
    end if

    do itype = 1, plot_types
    ; calculate differences for difference plots
    ; (mean plots already finished --> simply overwrite array var_grd)
    if (itype.eq.2) then
        obs = var_grd(idx_ref, :, :)
        do idx_mod = 0, dim_MOD - 1
            var_grd(idx_mod, :, :) = var_grd(idx_mod, :, :) - obs(:, :)
        end do
    end if
; A-laue_ax-

    ;; Loop over regridded models for plotting
    do idx_mod = 0, dim_MOD - 1

        ;; calculate mean
        A2 = where(ismissing(mask0), mask0, var_grd(idx_mod, :, :))
        copy_VarMeta(mask0, A2)
        time_mean = area_operations(A2, -90, 90, 0, 360, "average", True)
        if(time_mean.ge.1000) then
            fmt = tostring(toint(time_mean))
        else if(time_mean.lt.1000 .and. time_mean.ge.100) then
            fmt = sprintf("%5.1f", time_mean)
        else if(time_mean.lt.100 .and. time_mean.ge.10) then
            fmt = sprintf("%5.2f", time_mean)
        else if(time_mean.lt.10 .and. time_mean.ge.1) then
            fmt = sprintf("%5.3f", time_mean)
        else
            fmt = tostring(time_mean)
        end if
        end if
        end if
        end if

        ;; Define output workstation
        file_type = getenv("ESMValTool_output_file_type")
        if(ismissing(file_type)) then
            file_type = "PS"
        end if
        aux_title_info = ""
; A-laue_ax+
        if ((flag_mean).and.(idx_mod.eq.dim_MOD-1)) then
            imod = -1
            aux_title_info = "multi-model-mean"
        else
            imod = idx_mod
        end if
        if (itype.eq.2) then  ; difference plots
            aux_title_info = aux_title_info + "_diff"
        end if
;        output_filename = interface_get_figure_filename(diag_script_base,\
;            var0, field_type0, aux_title_info, idx_mod)
        output_filename = interface_get_figure_filename(diag_script_base,\
            var0, field_type0, aux_title_info, imod)
; A-laue_ax-
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(file_type, output_file_path)

        ;; Data to plot
        A3 = var_grd(idx_mod, :, :)

        ;; model specific ressources
        res@tiMainString = annots(idx_mod)
; A-laue_ax+
        if (itype.eq.2) then  ; difference plots
            res@tiMainString = res@tiMainString + " - " + annots(idx_ref)
        end if
; A-laue_ax-
        res@gsnLeftString = variable_info@short_name + \
            " [" + variable_info@units + "]"
        res@gsnRightString = fmt + " " + variable_info@units

; A-laue_ax+
        ; change color bar annotations for difference plots
        if (itype.eq.2) then  ; difference plots
            delete(res@cnLevels)
            res@cnLevels = variable_info@acp_dilvl
        end if
; A-laue_ax-

        ;; Plotting
        map = gsn_csm_contour_map(wks, A3, res)
        print("gv " + output_file_path + "." + file_type)

        ;; clean up
        delete([/wks, A2, A3, time_mean/])

    end do  ; idx_MOD, for plotting
; A-laue_ax+
    end do  ; itype (plot types)
; A-laue_ax-

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
