;##############################################################################
;# eyring13jgr_fig04.ncl                                                      #
;##############################################################################
;# Author: Irene Cionni (ENEA, Italy)                                         #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    * Calculates climatological annual mean tropospheric ozone columns      #
;#      and plots the geographical distribution.                              #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#     styleset: style set (“DEFAULT, "CMIP5")                                #
;#     font: overrides default font (e.g., 21, see                            #
;#         www.ncl.ucar.edu/Document/Graphics/Resources/tx.shtml#txFont)      #
;#     range_option: 0 = as in nml, 1 = overlapping time period               #
;#     lbLabelBarOn: plot a label bar (True, False)                           #
;#     e13fig04_ = "True"                                                     #
;#     e13fig04_list_chem_mod: list of models in the group “chem” (array of   #
;#         strings, default = (/"All"/))                                      #
;#     e13fig04_list_chem_mod_string: plotting label for group “chem”,        #
;#         e.g., "CMIP5"                                                      #
;#     e13fig04_list_nochem_mod: list of models in the group “nochem”         #
;#         (array of strings, default = (/""/))                               #
;#     e13fig01_list_nochem_mod_string: plotting label for group “nochem”,    #
;#         e.g., "NoChem"                                                     #
;#     e13fig04_diff_ref: name of reference model for difference plots,       #
;#         e.g., "AURA-MLS-OMI"                                               #
;#                                                                            #
;# Optional diag_script_info attributes:                                      #
;#     mpProjection: map projection                                           #
;#                                                                            #
;# Required variable_info attributes (variable specific):                     #
;#     variable_info@acp_cnlvl                                                #
;#     variable_info@map_ref_Levels                                           #
;#     variable_info@map_diff_Levels                                          #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    20151208-A-laue_ax: code clean-up, moved call to write_references to    #
;#                        beginning of code                                   #
;#    20140528-A_GoKl: written based on ancyc_lat.ncl                         #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
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
    dim_OBS = dimsizes(ind(models@project.eq."OBS"))
    OBS=models@name(ind(models@project.eq."OBS")) 
    

    var0 = variables(0)
    field_type0 = field_types(0)

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
    ;; References Information
    write_references(diag_script,      \  ; script name
                     "A_cion_ir",      \  ; authors
                     "",               \  ; contributors
                     "D_eyring13jgr",  \  ; diag_references
                     "E_auramls",      \  ; obs_references
                     "P_esmval")          ; proj_references

    ;********************diag_script_info*****************

    if ((isatt(diag_script_info,"e13fig04_list_chem_mod"))) then
        list_chem_mod = diag_script_info@e13fig04_list_chem_mod
        attrib_chem = "_CHEM"
        if (list_chem_mod.eq."All") then
            delete(list_chem_mod)
            list_chem_mod = models@name
            attrib_chem = " "
        end if
        ind_chem_mod = ind(models@name.eq.list_chem_mod)
    else
        list_chem_mod = models@name(ind(models@project.ne."OBS"))
        ind_chem_mod = ind(models@project.ne."OBS")
        attrib_chem = " "
    end if

    if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
        list_nochem_mod = diag_script_info@e13fig04_list_nochem_mod
        ind_nochem_mod = ind(models@name.eq.list_nochem_mod)
    end if

    if ((isatt(diag_script_info, "e13fig04_diff_ref"))) then
        OBS_ref = diag_script_info@e13fig04_diff_ref
        ind_OBS_ref = ind(models@name.eq.OBS_ref) 
    end if

    ;*************************************************

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
        if (isatt(diag_script_info, "range_option")) then
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
    if (isatt(diag_script, "mpProjection")) then
        res@mpProjection = diag_script_info@mpProjection
    else
        res@mpProjection = "CylindricalEquidistant"
    end if
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
    res@lbLabelFontHeightF = 0.025
    ;res@lbTitleString = variable_info@long_name

    if (isatt(diag_script_info, "gsnZonalMean")) then
        res@gsnZonalMean = diag_script_info@gsnZonalMean
    end if

    ;; Loop over regridded models for plotting
    do idx_mod = 0, dim_MOD - 1
        ;; calculate mean
        A2 = where(ismissing(mask0), mask0, var_grd(idx_mod, :, :))
        copy_VarMeta(mask0, A2)
        time_mean = area_operations(A2, -90, 90, 0, 360, "average", True)
        if (time_mean.ge.1000) then
            fmt = tostring(toint(time_mean))
        else if (time_mean.lt.1000 .and. time_mean.ge.100) then
            fmt = sprintf("%5.1f", time_mean)
        else if (time_mean.lt.100 .and. time_mean.ge.10) then
            fmt = sprintf("%5.2f", time_mean)
        else if (time_mean.lt.10 .and. time_mean.ge.1) then
            fmt = sprintf("%5.3f", time_mean)
        else
            fmt = tostring(time_mean)
        end if
        end if
        end if
        end if

        ;; Define output workstation
        file_type = getenv("ESMValTool_output_file_type")
        if (ismissing(file_type)) then
            file_type = "PS"
        end if
        aux_title_info = ""
        output_filename = interface_get_figure_filename(diag_script_base, \
            var0, field_type0, aux_title_info, idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(file_type, output_file_path)

        ;; Data to plot
        A3 = var_grd(idx_mod, :, :)

        ;; model specific ressources
        res@tiMainFontHeightF = 0.03
        res@tiMainString = annots(idx_mod)
        res@gsnLeftString = variable_info@short_name + \
            " [" + variable_info@units + "]"
        res@gsnRightString = fmt + " " + variable_info@units
        res@gsnRightStringFontHeightF = 0.023
        res@gsnLeftStringFontHeightF = 0.023

        ;; Plotting
        map = gsn_csm_contour_map(wks, A3, res)
        print("gv " + output_file_path + "." + file_type)

        ;; clean up
        delete([/wks, A2, A3, time_mean/])
    end do  ; idx_MOD, for plotting

    ;*********plot MMM plot lista_chem plot lista_nochem***********

    var_grd_mod=var_grd(ind_chem_mod, :, :)
    MMM_chem = dim_avg_n_Wrap(var_grd_mod, 0)

    if ((isatt(diag_script_info, "e13fig04_diff_ref"))) then
        var_grd_ref = var_grd(ind_OBS_ref, :, :)
    end if

    if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
        if(dimsizes(list_nochem_mod).gt.1) then
            var_grd_nochem_mod = var_grd(ind_nochem_mod, :, :)
            MMM_nochem = dim_avg_n_Wrap(var_grd_nochem_mod, 0)
            if ((isatt(diag_script_info, "e13fig04_diff_ref"))) then
                diff_var_grd_nochem = MMM_chem-var_grd_ref
            end if
        end if
    end if

    diff_var_grd = MMM_chem-var_grd_ref

    ;**************plot add**********

    ;; calculate mean MMM
    MMM_chem = where(ismissing(mask0), mask0, MMM_chem)
    copy_VarMeta(mask0, MMM_chem)
    time_mean = area_operations(MMM_chem, -90, 90, 0, 360, "average", True)

    if (time_mean.ge.1000) then
        fmt = tostring(toint(time_mean))
    else if (time_mean.lt.1000 .and. time_mean.ge.100) then
        fmt = sprintf("%5.1f", time_mean)
    else if (time_mean.lt.100 .and. time_mean.ge.10) then
        fmt = sprintf("%5.2f", time_mean)
    else if (time_mean.lt.10 .and. time_mean.ge.1) then
        fmt = sprintf("%5.3f", time_mean)
    else
        fmt = tostring(time_mean)
    end if
    end if
    end if
    end if

    ;; Define output workstation
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if
    aux_title_info = "MMM_chem"        
       
    output_filename = interface_get_figure_filename(diag_script_base, \
         var0, field_type0, aux_title_info, -1)
    output_file_path = output_dir + output_filename
    aux_title_info = "MMM_chem"

    aux_title_info1 = "panel"
    output_filename1 = interface_get_figure_filename(diag_script_base,\
        var0, field_type0, aux_title_info1, -1)
    output_file_path1 = output_dir + output_filename1
        
    wks1 = gsn_open_wks(file_type, output_file_path1)    
    ip = 0
    map_p = new(5, "graphic")

    if ((isatt(diag_script_info, "e13fig04_list_chem_mod_string"))) then
        res@tiMainString = diag_script_info@e13fig04_list_chem_mod_string
    else
        res@tiMainString = "MMM_chem"
    end if

    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = variable_info@short_name + \
        " [" + variable_info@units + "]"
    res@gsnRightString = fmt + " " + variable_info@units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023

    ;; Plotting
    map_p(ip) = gsn_csm_contour_map(wks1, MMM_chem, res)
    ip = ip + 1
;    print("gv " + output_file_path1 + "." + file_type)

    ;; clean up
    delete([/ time_mean/])

    ;*****************NO CHEM ***************************

    if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
        MMM_nochem = where(ismissing(mask0), mask0, MMM_nochem)
        copy_VarMeta(mask0, MMM_nochem)
        time_mean = area_operations(MMM_nochem, -90, 90, 0, 360, "average", True)
        if (time_mean.ge.1000) then
            fmt = tostring(toint(time_mean))
        else if (time_mean.lt.1000 .and. time_mean.ge.100) then
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
        if (ismissing(file_type)) then
            file_type = "PS"
        end if
        aux_title_info = "MMM_nochem"

        output_filename = interface_get_figure_filename(diag_script_base, \
            var0, field_type0, aux_title_info, -1)
        output_file_path = output_dir + output_filename

        if ((isatt(diag_script_info, "e13fig04_list_nochem_mod_string"))) then
            res@tiMainString = e13fig04_list_nochem_mod_string
        else
            res@tiMainString = "MMM_nochem"
        end if

        res@tiMainFontHeightF = 0.03
        res@gsnLeftString = variable_info@short_name + \
            " [" + variable_info@units + "]"
        res@gsnRightString = fmt + " " + variable_info@units
        res@gsnLeftStringFontHeightF = 0.023
        res@gsnRightStringFontHeightF = 0.023

        ;; Plotting
        map_p(ip) = gsn_csm_contour_map(wks1, MMM_nochem, res)
        ip = ip + 1
;        print("gv " + output_file_path1 + "." + file_type)

        ;; clean up
        delete([/time_mean/])
    end if

    ;****************************************************
    ;******************REF*******************************

    if ((isatt(diag_script_info, "e13fig04_diff_ref"))) then        
        var_grd_ref = where(ismissing(mask0), mask0, var_grd_ref)
        copy_VarMeta(mask0, var_grd_ref)
        time_mean = area_operations(var_grd_ref, -90, 90, 0, 360, "average", True)
        if (time_mean.ge.1000) then
            fmt = tostring(toint(time_mean))
        else if (time_mean.lt.1000 .and. time_mean.ge.100) then
            fmt = sprintf("%5.1f", time_mean)
        else if (time_mean.lt.100 .and. time_mean.ge.10) then
            fmt = sprintf("%5.2f", time_mean)
        else if (time_mean.lt.10 .and. time_mean.ge.1) then
            fmt = sprintf("%5.3f", time_mean)
        else
            fmt = tostring(time_mean)
        end if
        end if
        end if
        end if

        ;; Define output workstation
        file_type = getenv("ESMValTool_output_file_type")
        if (ismissing(file_type)) then
            file_type = "PS"
        end if
        aux_title_info = diag_script_info@e13fig04_diff_ref
        output_filename = interface_get_figure_filename(diag_script_base, \
            var0, field_type0, aux_title_info, -1)
        output_file_path = output_dir + output_filename
        res@tiMainString = diag_script_info@e13fig04_diff_ref 
        res@tiMainFontHeightF = 0.03
        res@gsnLeftString = variable_info@short_name + \
            " [" + variable_info@units + "]"
        res@gsnRightString = fmt + " " + variable_info@units
        res@gsnRightStringFontHeightF = 0.023
        res@gsnLeftStringFontHeightF = 0.023

        ;; Plotting
        map_p(ip) = gsn_csm_contour_map(wks1, var_grd_ref, res)
        ip = ip + 1
;        print("gv " + output_file_path + "." + file_type)

        ;; clean up
        delete([/time_mean/])

        ;*****************diff**************************************

        diff_var_grd=where(ismissing(mask0), mask0, diff_var_grd)
        copy_VarMeta(mask0, diff_var_grd)
        time_mean = area_operations(diff_var_grd, -90, 90, 0, 360, "average", True)

        if (time_mean.ge.1000) then
            fmt = tostring(toint(time_mean))
        else if (time_mean.lt.1000 .and. time_mean.ge.100) then
            fmt = sprintf("%5.1f", time_mean)
        else if (time_mean.lt.100 .and. time_mean.ge.10) then
            fmt = sprintf("%5.2f", time_mean)
        else if (time_mean.lt.10 .and. time_mean.ge.1) then
            fmt = sprintf("%5.3f", time_mean)
        else
            fmt = tostring(time_mean)
        end if
        end if
        end if
        end if

        ;; Define output workstation
        file_type = getenv("ESMValTool_output_file_type")
        if (ismissing(file_type)) then
            file_type = "PS"
        end if
        aux_title_info = "diff_MMM-"+diag_script_info@e13fig04_diff_ref
        output_filename = interface_get_figure_filename(diag_script_base, \
            var0, field_type0, aux_title_info, -1)
        output_file_path = output_dir + output_filename

        if ((isatt(diag_script_info, "e13fig04_list_chem_mod_string"))) then
            res@tiMainString = diag_script_info@e13fig04_list_chem_mod_string \
                + " - " + diag_script_info@e13fig04_diff_ref
        else
            res@tiMainString = "MMM_chem"
        end if

        res@tiMainFontHeightF = 0.03
        res@gsnLeftString = variable_info@short_name + \
            " [" + variable_info@units + "]"
        res@gsnRightString = fmt + " " + variable_info@units
        res@gsnRightStringFontHeightF = 0.023
        res@gsnLeftStringFontHeightF = 0.023

        ;; Plotting
        delete(res@cnLevelSelectionMode)
        delete(res@cnLevels)
        map_p(ip) = gsn_csm_contour_map(wks1, diff_var_grd, res)
        ip = ip + 1
;        print("gv " + output_file_path1 + "." + file_type)

        ;; clean up
        delete([/ time_mean/])

        if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
            diff_var_grd_nochem = where(ismissing(mask0), mask0, diff_var_grd_nochem)
            copy_VarMeta(mask0, diff_var_grd_nochem)
            time_mean = area_operations(diff_var_grd_nochem, -90, 90, 0, 360, "average", True)
            if (time_mean.ge.1000) then
                fmt = tostring(toint(time_mean))
            else if (time_mean.lt.1000 .and. time_mean.ge.100) then
                fmt = sprintf("%5.1f", time_mean)
            else if (time_mean.lt.100 .and. time_mean.ge.10) then
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
            if (ismissing(file_type)) then
                file_type = "PS"
            end if
            aux_title_info = "diff_MMM_nochem-"+diag_script_info@e13fig04_diff_ref
            output_filename = interface_get_figure_filename(diag_script_base, \
                var0, field_type0, aux_title_info,-1)
            output_file_path = output_dir + output_filename

            if ((isatt(diag_script_info, "e13fig04_list_nochem_mod_string"))) then
                res@tiMainString = e13fig04_list_nochem_mod_string + " - " \
                    + diag_script_info@e13fig04_diff_ref
            else
                res@tiMainString = "MMM_nochem"
            end if

            res@tiMainFontHeightF = 0.03
            res@gsnLeftString = variable_info@short_name + \
                " [" + variable_info@units + "]"
            res@gsnRightString = fmt + " " + variable_info@units
            res@gsnRightStringFontHeightF = 0.023
            res@gsnLeftStringFontHeightF = 0.023

            ;; Plotting
            delete(res@cnLevelSelectionMode)
            delete(res@cnLevels)
            map_p(ip) = gsn_csm_contour_map(wks1, diff_var_grd_nochem, res)
            ip = ip + 1
;            print("gv " + output_file_path1 + "." + file_type)

            ;; clean up
            delete([/ time_mean/])
        end if
    end if

    resP = True
    gsn_panel(wks1, map_p(0:ip - 1), (/1, ip/), resP) 
    print("gv " + output_file_path1 + "." + file_type)

    ;***********************************************************************
    ;********plot diff *****************************************************

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
