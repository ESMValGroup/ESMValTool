;;#############################################################################
;; GENERAL ROUTINES FOR GO PLOTTING
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;    function  get_title_suffix
;;    procedure remove_attrs
;;    procedure plot_two_by_one
;;    procedure plot_three_by_one_diff
;;    procedure two_by_one
;;    procedure three_by_one_diff
;;    procedure plot_three_by_one_vector
;;    procedure three_by_one_vector
;;    procedure plot_multipanel
;;    procedure multipanel
;;    procedure plot_multipanel_vector
;;    procedure multipanel_vector
;;    procedure seasonal_plot
;;    procedure xy_plot_wrapper
;;    procedure ts_line_wrapper
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/logging.ncl"


begin
    loadscript(exact_panel_position_file)
end

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon_panel_positions.ncl"
load "./diag_scripts/aux/SAMonsoon/exact_panel_positions_generic.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/xy_line.ncl"

undef("get_title_suffix")
function get_title_suffix(plot_settings[1]: logical)
local main_title_suffix
begin

    if (isatt(plot_settings, "main_title_suffix")) then
        main_title_suffix = plot_settings@main_title_suffix
    else
        main_title_suffix = ""
    end if
    return main_title_suffix

end

undef("remove_attrs")
procedure remove_attrs(p_settings[1]:logical,
                       attributes[*]:string)
local idx, settings
begin
    do idx = 0, dimsizes(attributes) - 1
        if (isatt(p_settings, attributes(idx))) then
            delete(p_settings@$attributes(idx)$)
        end if
    end do
end

undef("plot_two_by_one")
procedure plot_two_by_one(wks[1]:graphic,\
                          res[1]:logical,\
                           di[1]:logical,\
                plot_settings[1]:logical,\
             valid_statistics[*]:string,
                storage_name1[1]:string,\
                storage_name2[1]:string,\
                    debuginfo[1]:logical,\
                storage_vault[1]:logical,\
                      idx_mod[1]:integer)
begin
    res = panel_three_by_one(res, 0)  ; Update resource for first plot
    res@cnLevelSelectionMode = "ExplicitLevels"
    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if

    dummy_array = (/1e+20/)
    dummy_array@_FillValue = 1e+20

    remove_attrs(res, (/"cnLevels"/))
    cn_levels_string = str_sub_str("cn_levels_" + storage_name1 + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    remove_attrs(plot_settings, (/"idx"/))
    plot_settings@idx = plot_settings@idx_ref

    ;; -------------------------------------
    ;; Plot model entry (first plot)
    ;; -------------------------------------
    main_title_prefix = ""  ; model_info@model(idx_mod)
    main_title_suffix = ""
    res = panel_two_by_one(res, 0)  ; Update resource for second plot

    storage_record = str_join((/storage_name1, sprinti("%i", idx_mod)/), str_vault_sep)
    upper = retrieve_from_vault(storage_vault, storage_record)
    plot_settings@type = "mean"

    remove_attrs(plot_settings, (/"idx", "gsn_right_string"/))
    plot_settings@idx = idx_mod
    plot_settings@gsn_right_string = upper@gsn_right_string

    statistics = True
    statistics = compute_stat((/"yrs", "mean", "bob", "eio", "mean-corr"/), valid_statistics, upper, dummy_array)

    plot_settings@skip_default_title = False
    plot1 = single_contour(wks,\
                           upper,\
                           main_title_prefix,\
                           main_title_suffix,\
                           plot_settings,\
                           debuginfo,\
                           statistics,\
                           res)
    delete(statistics)

    ;; ---------------------------------------
    ;; Plot model entry (second plot)
    ;; ---------------------------------------
    main_title_prefix = ""
    main_title_suffix = ""
    res = panel_two_by_one(res, 1)  ; Update resource for third plot
    delete(res@cnLevels)
    cn_levels_string = str_sub_str("cn_levels_" + storage_name2 + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    storage_record = str_join((/storage_name2, sprinti("%i", idx_mod)/), str_vault_sep)
    lower = retrieve_from_vault(storage_vault, storage_record)
    plot_settings@type = "mean"

    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = idx_mod

    statistics = True
    statistics = compute_stat((/"yrs", "mean", "bob", "eio", "mean-corr"/), valid_statistics, lower, dummy_array)
    plot_settings@gsn_right_string = lower@gsn_right_string
    plot_settings@skip_default_title = True
    plot2 = single_contour(wks,\
                           lower,\
                           main_title_prefix,\
                           main_title_suffix,\
                           plot_settings,\
                           debuginfo,\
                           statistics,\
                           res)
    delete(statistics)

    txres = True
    if (any(idx_mod .eq. plot_settings@idx_ref)) then
        inset_top_text(wks, plot1, "Reference", txres)
        inset_top_text(wks, plot2, "Reference", txres)
    end if
    lbtitle = plot_settings@lb_units
    two_by_one_labelbar(wks, plot2, lbtitle)

    if (debuginfo) then
        txres@txFuncCode = "~"
        debugbox = write_info(debuginfo)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot1)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot2)
        drawNDCGrid(wks)
    end if
    draw(plot1)
    draw(plot2)
    if (debuginfo) then
        place_description(wks,\
                          debuginfo@description,\
                          debuginfo@description_ycoord)
    end if
end

undef("plot_three_by_one_diff")
procedure plot_three_by_one_diff(wks[1]:graphic,\
                                 res[1]:logical,\
                                  di[1]:logical,\
                       plot_settings[1]:logical,\
                    valid_statistics[*]:string,\
                        storage_name[1]:string,\
                           debuginfo[1]:logical,\
                       storage_vault[1]:logical,\
                              ref[*][*]:numeric,\
                             idx_mod[1]:integer)
begin
    res = panel_three_by_one(res, 0)  ; Update resource for first plot
    res@cnLevelSelectionMode = "ExplicitLevels"
    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if

    dummy_array = (/1e+20/)
    dummy_array@_FillValue = 1e+20

    ;; -------------------------------------
    ;; Plot reference plot (top plot)
    ;; -------------------------------------
    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    main_title_prefix = "(1) " + di@season + "-" + storage_name + " of "
    main_title_suffix = get_title_suffix(plot_settings)
    plot_settings@type = "ref"

    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = plot_settings@idx_ref

    statistics = True
    statistics = compute_stat((/"yrs", "mean", "bob", "eio"/), valid_statistics, ref, dummy_array)
    plot_ref = single_contour(wks,\
                              ref,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)
    delete(statistics)
    txres = True
    inset_top_text(wks, plot_ref, "Reference", txres)

    three_by_one_labelbar(wks, plot_ref, plot_settings)

    ;; -------------------------------------
    ;; Plot model entry (middle plot)
    ;; -------------------------------------
    main_title_prefix = "(2) " + di@season + "-" + storage_name + " of "
    res = panel_three_by_one(res, 1)  ; Update resource for second plot

    storage_record = str_join((/storage_name, sprinti("%i", idx_mod)/), str_vault_sep)
    curr = retrieve_from_vault(storage_vault, storage_record)
    plot_settings@type = "mean"

    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = idx_mod

    statistics = True
    statistics = compute_stat((/"yrs", "mean", "bob", "eio", "mean-corr"/), valid_statistics, curr, dummy_array)
    statistics = compute_stat((/"corr"/), valid_statistics, curr, ref)
    plot = single_contour(wks,\
                          curr,\
                          main_title_prefix,\
                          main_title_suffix,\
                          plot_settings,\
                          debuginfo,\
                          statistics,\
                          res)
    delete(statistics)

    ;; ---------------------------------------
    ;; Plot reference model diff (bottom plot)
    ;; ---------------------------------------
    main_title = "(2) - (1)"
    main_title_suffix = ""
    res = panel_three_by_one(res, 2)  ; Update resource for third plot
    delete(res@cnLevels)
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_diff_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$
    res@cnFillPalette = "BlWhRe"
    ;; Plot mean of differences of first and second data set, first
    ;; interpolate to the reference data set grid (conservative interpolation)
    diff_model_ref = get_model_minus_ref(curr, ref)
    plot_settings@type = "diff"
    plot_settings@skip_default_title = True
    statistics = True
    statistics = compute_stat((/"mean"/), valid_statistics, curr, dummy_array)
    statistics = compute_stat((/"rmse"/), valid_statistics, curr, ref)
    plot_diff = single_contour(wks,\
                               diff_model_ref,\
                               main_title,\
                               main_title_suffix,\
                               plot_settings,\
                               debuginfo,\
                               statistics,\
                               res)
    delete(statistics)
    delete(diff_model_ref)
    plot_settings@skip_default_title = False
    ;;;; remove the light_centered colour palette for climatological plots (OA)
    delete(res@cnFillPalette)
    if (debuginfo) then
        txres@txFuncCode = "~"
        debugbox = write_info(debuginfo)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_ref)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_diff)
        drawNDCGrid(wks)
    end if
    draw(plot_ref)
    draw(plot_diff)
    draw(plot)
    if (debuginfo) then
        place_description(wks,\
                          debuginfo@description,\
                          debuginfo@description_ycoord)
    end if
end

undef("two_by_one")
procedure two_by_one(storage_vault[1]: logical,
                     di[1]: logical,
                     plot_settings[1]: logical,
                     storage_name1[1]: string,
                     storage_name2[1]: string,
                     debuginfo[1]: logical,
                     valid_statistics[*]: string,
                     res[1]: logical)
; Description:
;   For each 2D-model field stored in the 'storage_vault', plot two separate
;   figures. The various logcials carry switches for the plot routine behaviour.
local aux_title_info,      cn_levels_string,   curr,\
    diag_description,    diag_script_base,   dim_MOD,\
    dim_VAR,             idx_mod,            dummy_array,\
    field_type0,         lbtitle,            main_title,\
    main_title_prefix,   output_dir,         output_filename,\
    output_file_path,    plot,               plot_diff,\
    plot_ref,            ref,                res,\
    storage_record,      textres,            var0,\
    wks,                 diff_model_ref,     txres
begin

    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    ;; Handle the case when using two variables
    if (dimsizes(variables) .gt. 1) then
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(plot_settings@diag_script)
    output_dir = config_user_info@plot_dir

    ;; --------------------------------
    ;; Static resources for these plots
    ;; --------------------------------
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False
    res@gsnFrame = False
    res@gsnDraw = False
    res@gsnAddCyclic = False
    do idx_mod = 0, dim_MOD - 1
        ;; ----------------------------------------
        ;; Define output workstation
        ;; ----------------------------------------
        aux_title_info = di@season + "-" + storage_name1 + "-" + storage_name2
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0 + var1,\
                                                        field_type0 + field_type1,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(output_file_type, output_file_path)
        plot_two_by_one(wks,\
                        res,\
                        di,\
                        plot_settings,\
                        valid_statistics,\
                        storage_name1,\
                        storage_name2,\
                        False,\
                        storage_vault,\
                        idx_mod)
        frame(wks)  ; Workstation output
        delete(wks)

        ;; Debug version
        if (debuginfo) then
            aux_title_info = di@season + "-" + storage_name1 + "-" + storage_name2 + "-debug"
            output_filename = interface_get_figure_filename(diag_script_base,\
                                                            var0 + var1,\
                                                            field_type0 + field_type1,\
                                                            aux_title_info,\
                                                            idx_mod)
            output_file_path = output_dir + output_filename
            wks_debug = gsn_open_wks(output_file_type, output_file_path)
            plot_two_by_one(wks_debug,\
                            res,\
                            di,\
                            plot_settings,\
                            valid_statistics,\
                            storage_name1,\
                            storage_name2,\
                            debuginfo,\
                            storage_vault,\
                            idx_mod)
            frame(wks_debug)  ; Workstation output
            delete(wks_debug)
        end if

    end do  ; idx_mod = 0, dim_MOD - 1

end
undef("three_by_one_diff")
procedure three_by_one_diff(storage_vault[1]: logical,
                            di[1]: logical,
                            plot_settings[1]: logical,
                            storage_name[1]: string,
                            debuginfo[1]: logical,
                            valid_statistics[*]: string,
                            res[1]: logical)
; Description:
;   For each 2D-model field stored in the 'storage_vault', plot a separate figure
;   with a reference plot at the top, the model in the middle, and a diff at the
;   bottom. The various logcials carry switches for the plot routine behaviour.
local aux_title_info,      cn_levels_string,   curr,\
    diag_description,    diag_script_base,   dim_MOD,\
    dim_VAR,             idx_mod,            dummy_array,\
    field_type0,         lbtitle,            main_title,\
    main_title_prefix,   output_dir,         output_filename,\
    output_file_path,    plot,               plot_diff,\
    plot_ref,            ref,                res,\
    storage_record,      textres,            var0,\
    wks,                 diff_model_ref,     txres
begin

    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    ;; Handle the case when using two variables
    if (dimsizes(variables) .gt. 1) then
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(plot_settings@diag_script)
    output_dir = config_user_info@plot_dir

    ;; --------------------------------
    ;; Static resources for these plots
    ;; --------------------------------
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False
    res@gsnFrame = False
    res@gsnDraw = False
    res@gsnAddCyclic = False

    storage_record = str_join((/storage_name, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    ref = retrieve_from_vault(storage_vault, storage_record)

    do idx_mod = 0, dim_MOD - 1
        if (any(idx_mod .eq. plot_settings@idx_ref)) then
            continue
        end if

        ;; ----------------------------------------
        ;; Define output workstation
        ;; ----------------------------------------
        aux_title_info = di@season + "-" + storage_name
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0 + var1,\
                                                        field_type0 + field_type1,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(output_file_type, output_file_path)
        plot_three_by_one_diff(wks,\
                               res,\
                               di,\
                               plot_settings,\
                               valid_statistics,\
                               storage_name,\
                               False,\
                               storage_vault,\
                               ref,\
                               idx_mod)
        frame(wks)  ; Workstation output
        delete(wks)

        ;; Debug version
        if (debuginfo) then
            aux_title_info = di@season + "-" + storage_name + "-debug"
            output_filename = interface_get_figure_filename(diag_script_base,\
                                                            var0 + var1,\
                                                            field_type0 + field_type1,\
                                                            aux_title_info,\
                                                            idx_mod)
            output_file_path = output_dir + output_filename
            wks_debug = gsn_open_wks(output_file_type, output_file_path)
            plot_three_by_one_diff(wks_debug,\
                                   res,\
                                   di,\
                                   plot_settings,\
                                   valid_statistics,\
                                   storage_name,\
                                   debuginfo,\
                                   storage_vault,\
                                   ref,\
                                   idx_mod)
            frame(wks_debug)  ; Workstation output
            delete(wks_debug)
        end if

    end do  ; idx_mod = 0, dim_MOD - 1

end

undef("plot_three_by_one_vector")
procedure plot_three_by_one_vector(wks[1]:graphic,
                                   res[1]:logical,
                                    di[1]:logical,
                         plot_settings[1]:logical,
                      valid_statistics[*]:string,
                          storage_name[1]:string,
                             debuginfo[1]:logical,
                         storage_vault[1]:logical,
                             ua_ref[*][*]:numeric,
                             va_ref[*][*]:numeric,
                          speed_ref[*][*]:numeric,
                               idx_mod[1]:integer)
begin
    res = panel_three_by_one(res, 0)  ; Update resource for first plot
    res@cnLevelSelectionMode = "ExplicitLevels"
    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if

    ;; ------------------------------
    ;; Plot reference plot (top plot)
    ;; ------------------------------
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    main_title_prefix = "(1) " + di@season + "-" + storage_name + " of "
    main_title_suffix = ""
    plot_settings@type = "ref"

    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = plot_settings@idx_ref

    dummy_array = (/1e+20/)
    dummy_array@_FillValue = 1e+20

    statistics = True
    statistics = compute_stat((/"yrs", "mean"/), valid_statistics, speed_ref, dummy_array)
    plot_ref = single_contour(wks,\
                              speed_ref,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)
    delete(statistics)
    txres = True
    inset_top_text(wks, plot_ref, "Reference", txres)

    max_speed_ref = max(speed_ref)
    if (storage_name .ne. "stddev") then
        plot_ref_v = single_vector(wks,\
                                   max_speed_ref,\
                                   speed_ref,\
                                   ua_ref,\
                                   va_ref,\
                                   main_title_prefix,\
                                   main_title_suffix,\
                                   plot_settings,\
                                   debuginfo,\
                                   res)
        overlay(plot_ref, plot_ref_v)
    end if
    three_by_one_labelbar(wks, plot_ref, plot_settings@lb_units)

    ;; -------------------------------------
    ;; Plot model entry plot (middle plot)
    ;; -------------------------------------
    main_title_prefix = "(2) " + di@season + "-" + storage_name + " of "
    res = panel_three_by_one(res, 1)  ; Update resource for second plot

    storage_record = str_join((/storage_name, "ua", sprinti("%i", idx_mod)/), str_vault_sep)
    ua = retrieve_from_vault(storage_vault, storage_record)

    storage_record = str_join((/storage_name, "va", sprinti("%i", idx_mod)/), str_vault_sep)
    va = retrieve_from_vault(storage_vault, storage_record)

    storage_record = str_join((/storage_name, sprinti("%i", idx_mod)/), str_vault_sep)
    speed = retrieve_from_vault(storage_vault, storage_record)
    plot_settings@type = "mean"

    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = idx_mod

    ;; Interpolate to reference grid for pattern correlation calculation
    diff_model_ref = get_model_minus_ref(speed, speed_ref)

    statistics = True
    statistics = compute_stat((/"yrs", "mean"/), valid_statistics, speed, dummy_array)
    statistics = compute_stat((/"corr"/), valid_statistics, speed, speed_ref)
    plot = single_contour(wks,\
                          speed,\
                          main_title_prefix,\
                          main_title_suffix,\
                          plot_settings,\
                          debuginfo,\
                          statistics,\
                          res)
    delete(statistics)
    if (storage_name .ne. "stddev") then
        plot_v = single_vector(wks,\
                               max_speed_ref,\
                               speed,\
                               ua,\
                               va,\
                               main_title_prefix,\
                               main_title_suffix,\
                               plot_settings,\
                               debuginfo,\
                               res)
        overlay(plot, plot_v)
    end if

    ;; --------------------------------------------
    ;; Plot model reference diff plot (bottom plot)
    ;; --------------------------------------------
    ;; Plot mean of differnces of first and second data set, first
    ;; interpolate to the reference data set grid (bilinear interpolation)

    diff_model_ref = get_model_minus_ref(speed, speed_ref)
    diff_model_ua_ref = get_model_minus_ref(ua, ua_ref)
    diff_model_va_ref = get_model_minus_ref(va, va_ref)

    main_title = "(1) - (2)"
    res = panel_three_by_one(res, 2)  ; Update resource for third plot
    delete(res@cnLevels)
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_diff_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    plot_settings@type = "diff"
    plot_settings@skip_default_title = True

    statistics = True
    statistics = compute_stat((/"mean"/), valid_statistics, diff_model_ref, dummy_array)
    statistics = compute_stat((/"rmse", "corr"/), valid_statistics, speed, speed_ref)
    plot_diff = single_contour(wks,\
                               diff_model_ref,\
                               main_title,\
                               main_title_suffix,\
                               plot_settings,\
                               debuginfo,\
                               statistics,\
                               res)
    delete(statistics)
    if (storage_name .ne. "stddev") then
        plot_diff_v = single_vector(wks,\
                                    max_speed_ref,\
                                    diff_model_ref,\
                                    diff_model_ua_ref,\
                                    diff_model_va_ref,\
                                    main_title,\
                                    main_title_suffix,\
                                    plot_settings,\
                                    debuginfo,\
                                    res)
        overlay(plot_diff, plot_diff_v)
    end if
    plot_settings@skip_default_title = False

    if (debuginfo) then
        debugbox = write_info(debuginfo)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_ref)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_diff)
        drawNDCGrid(wks)
    end if
    draw(plot_ref)
    draw(plot_diff)
    draw(plot)
    if (debuginfo) then
        place_description(wks,\
                          debuginfo@description,\
                          debuginfo@description_ycoord)
    end if
end

undef("three_by_one_vector")
procedure three_by_one_vector(storage_vault[1]: logical,
                              di[1]: logical,
                              plot_settings[1]: logical,
                              storage_name[1]: string,
                              debuginfo[1]: logical,
                              valid_statistics[*]: string,
                              res[1]: logical)
; Description:
;   For each set of three 2D-model fields (vector_x, vector_y and abs) stored
;   in the 'storage_vault', plot a separate figure with a reference plot at
;   the top (abs + vector f, the model in the middle, and a diff at the
;   bottom.  The various logcials carry switches for the plot routine behaviour.
local aux_title_info,      cn_levels_string,   curr,\
    diag_description,    diag_script_base,   dim_MOD,\
    dim_VAR,             idx_mod,                        \
    field_type0,         lbtitle,            main_title,\
    main_title_prefix,   output_dir,         output_filename,\
    output_file_path,    plot,               plot_diff,\
    plot_ref,            ref,                res,\
    storage_record,      textres,            var0,\
    wks,                 wks_debug,          txres
begin
    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(plot_settings@diag_script)
    output_dir = config_user_info@plot_dir

    ;; --------------------------------
    ;; Static resources for these plots
    ;; --------------------------------
    res@mpFillOn = False
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False
    res@gsnFrame = False
    res@gsnDraw = False
    res@gsnAddCyclic = False

    txres = True
    txres@txFuncCode = "~"

    ;; Fetch reference plot
    storage_record = str_join((/storage_name, "ua", sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    ua_ref = retrieve_from_vault(storage_vault, storage_record)

    storage_record = str_join((/storage_name, "va", sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    va_ref = retrieve_from_vault(storage_vault, storage_record)

    storage_record = str_join((/storage_name, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    speed_ref = retrieve_from_vault(storage_vault, storage_record)

    do idx_mod = 0, dim_MOD - 1
        if (idx_mod .eq. plot_settings@idx_ref) then
            continue
        end if

        ;; -------------------------
        ;; Define output workstation
        ;; -------------------------
        aux_title_info = di@season + "-" + storage_name
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0 + var1,\
                                                        field_type0 + field_type1,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(output_file_type, output_file_path)
        plot_three_by_one_vector(wks,\
                                 res,\
                                 di,\
                                 plot_settings,\
                                 valid_statistics,\
                                 storage_name,\
                                 False,\
                                 storage_vault,\
                                 ua_ref,\
                                 va_ref,\
                                 speed_ref,\
                                 idx_mod)
        frame(wks)
        delete(wks)

        ;; Debug version
        if (debuginfo) then
            aux_title_info = di@season + "-" + storage_name + "-debug"
            output_filename = interface_get_figure_filename(diag_script_base,\
                                                            var0 + var1,\
                                                            field_type0 + field_type1,\
                                                            aux_title_info,\
                                                            idx_mod)
            output_file_path = output_dir + output_filename
            wks_debug = gsn_open_wks(output_file_type, output_file_path)
            plot_three_by_one_vector(wks_debug,\
                                     res,\
                                     di,\
                                     plot_settings,\
                                     valid_statistics,\
                                     storage_name,\
                                     debuginfo,\
                                     storage_vault,\
                                     ua_ref,\
                                     va_ref,\
                                     speed_ref,\
                                     idx_mod)
            frame(wks_debug)
            delete(wks_debug)
        end if
    end do  ; idx_mod = start_idx, dim_MOD - 1
end

undef("plot_multipanel")
procedure plot_multipanel(cols[*]:integer,
                          rows[*]:float,
                          curr_idx[1]:integer,
                          curr_page[1]:integer,
                          res[1]:logical,
                          storage_name[1]:string,
                          storage_vault[1]:logical,
                          wks[1]:graphic,
                          di[1]:logical,
                          plot_settings[1]:logical,
                          valid_statistics[*]:string,
                          debuginfo[1]:logical,
                          figures_per_page[*]:integer,
                          model_panel_placement[*]:integer,
                          figure_panel_placement[*]:integer,
                          plot_array[*]:graphic,
                          type_specifier[1]:string,
                          no_figures_on_this_page[1]:integer)
begin
    ;; Update position, labelbar and title
    curr_figure_pos = figure_panel_placement(curr_idx)
    res = panel_n_by_cols(res, curr_figure_pos, rows, cols, figures_per_page(curr_page))

    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    main_title_prefix = ""
    main_title_suffix = ""

    ;; Fetch reference plot
    storage_record = str_join((/storage_name, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    ref = retrieve_from_vault(storage_vault, storage_record)

    ;; Reference plot
    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = plot_settings@idx_ref

    dummy_array = (/1e+20/)
    dummy_array@_FillValue = 1e+20

    statistics = True
    statistics = compute_stat((/"yrs", "mean", "bob", "eio", "mean-corr"/), valid_statistics, ref, dummy_array)
    plot_ref = single_contour(wks,\
                              ref,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)
    delete(statistics)
    txres = True
    txres@txFuncCode = "~"
    if (debuginfo) then
        debugbox = write_info(debuginfo)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_ref)
    end if

    idx_fig = figure_panel_placement(curr_idx)
    plot_array(idx_fig) = plot_ref

    ;; Skip past the reference plot
    curr_idx = curr_idx + 1

    lbres = True

    if (plot_settings@type .eq. "diff") then
        inset_top_text(wks, plot_ref, "REF", txres)
        inset_labelbar(wks, plot_ref, res, "REF", lbres)
        main_title_suffix = " - REF"
    else
        inset_top_text(wks, plot_ref, "Reference", txres)
        main_title_suffix = ""
    end if

    delete(res@cnLevels)
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + type_specifier + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$
    ;; -----------------------------
    ;; Create the non-reference plots
    ;; -----------------------------
    do curr_fig = 1, figures_per_page(curr_page) - 1

        main_title_prefix = ""
        idx_mod = model_panel_placement(curr_idx)
        idx_fig = figure_panel_placement(curr_idx)

        ;; Skip reference models
        if any((idx_mod .eq. plot_settings@idx_ref)) then
            continue
        end if

        ;; Update placement and labelbar colors
        res = panel_n_by_cols(res,\
                              figure_panel_placement(curr_idx),\
                              rows,\
                              cols,\
                              figures_per_page(curr_page))

        storage_record = str_join((/storage_name, sprinti("%i", idx_mod)/), str_vault_sep)
        curr = retrieve_from_vault(storage_vault, storage_record)

        statistics = True
        statistics = compute_stat((/"yrs"/), valid_statistics, curr, dummy_array)
        if (plot_settings@type .eq. "diff") then
            statistics = compute_stat((/"rmse"/), valid_statistics, curr, ref)
            diff_model_ref = get_model_minus_ref(curr, ref)
            delete(curr)
            curr = diff_model_ref
            delete(diff_model_ref)
        else
            statistics = compute_stat((/"corr"/), valid_statistics, curr, ref)
        end if
        statistics = compute_stat((/"mean", "bob", "eio", "mean-corr"/), valid_statistics, curr, dummy_array)

        if (isatt(plot_settings, "idx")) then
            delete(plot_settings@idx)
        end if
        plot_settings@idx = idx_mod
        plot = single_contour(wks,\
                              curr,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)

        if (debuginfo) then
            debugbox = write_info(debuginfo)
            am_infobox_id = place_debuginfo(wks, debugbox, txres, plot)
        end if

        plot_array(idx_fig) = plot

        ;; Update index to point to next field
        curr_idx = curr_idx + 1
        delete(curr)
        delete(statistics)

    end do  ; curr_fig=1, figures_per_page(curr_page) - 1

    plottype_lbres = False

    ;; ---------------------------------------------------------
    ;; Create an blank plot for shared labelbar placement (mean)
    ;; ---------------------------------------------------------
    blank_plot = add_blank_plot_title(wks,\
                                      di@season + "-" + plot_settings@type + plot_settings@part_of_header,\
                                      rows,\
                                      cols)
    ;; Create shared labelbar
    n_by_cols_labelbar(wks,\
                       blank_plot,\
                       plot_array(no_figures_on_this_page - 1),\
                       rows,\
                       cols,\
                       plot_settings@lb_units,\
                       plottype_lbres)

    ;;---------------------
    ;; Draw mean value plot
    ;;---------------------
    if (debuginfo) then
        drawNDCGrid(wks)
    end if
    draw(plot_array)
    draw(blank_plot)
    if (debuginfo) then
        place_description(wks,\
                          debuginfo@description,\
                          debuginfo@description_ycoord)
    end if
end

undef("multipanel")
procedure multipanel(storage_vault[1]: logical,
                     di[1]: logical,
                     plot_settings[1]: logical,
                     storage_name[1]: string,
                     debuginfo[1]: logical,
                     valid_statistics[*]: string,
                     res[1]: logical)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local aux_title_info,      blank_plot,              cn_levels_string,\
    cols,                curr,                    curr_fig,\
    curr_figure_pos,     curr_idx,                curr_idx_debug,            curr_page,\
    diag_script_base,    dim_MOD,                 dim_VAR,\
    dummy_array,         field_type0,             figure_panel_placement,    figures_per_page,\
    idx_fig,             idx_mod,                 lbres,\
    main_title_prefix,                            \
                         model_panel_placement,   no_figures_on_this_page,\
    output_dir,          output_filename,         output_file_path,\
    page_no,             plot,                    plot_array,\
    plot_mean_diff,      plot_ref,                plot_stddev,\
    plot_stddev_diff,    plottype_lbres,          ref,\
    res,                 rows,                    storage_record,\
    total_no_of_pages,   txres,                   type_specifier,\
    var0,                wks
begin
    ; No of models may differ from dimsizes(dim_MOD) if variables
    ; from different models are combined (e.g., 'ts' from HadISST
    ; and 'pr' from TRMM)
    ref_no = 0
    stored_ref_name = "model"
    dim_MOD = no_unique_vault_entries(storage_vault, stored_ref_name, ref_no)

    dim_VAR = dimsizes(variables)

    var0 = variables(0)
    field_type0 = field_types(0)
    if (dimsizes(variables) .gt. 1) then
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(plot_settings@diag_script)
    output_dir = config_user_info@plot_dir

    ;; --------------------------------
    ;; Static resources for these plots
    ;; --------------------------------
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False
    res@gsnFrame = False
    res@gsnDraw = False
    res@lbLabelBarOn = False
    res@gsnAddCyclic = False

    ;; --------------------------------------
    ;; Compute the layout of paneled figures
    ;; --------------------------------------
    figures_per_page = get_figures_per_page(dim_MOD,\
                                            max_figures_pp,\
                                            min_figures_pp)

    ;; Which model goes where across all pages
    model_panel_placement = new((/sum(figures_per_page)/), integer)

    ;; Which model goes where on each page?
    figure_panel_placement = new((/sum(figures_per_page)/), integer)
    place_models_on_pages(model_info,\
                          plot_settings@idx_ref,\
                          figures_per_page,\
                          model_panel_placement,\
                          figure_panel_placement)

    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir

    if (plot_settings@type .eq. "diff") then
        type_specifier = "_diff"
    else
        type_specifier = ""
    end if

    ;; ---------------------------
    ;; Loop over all output pages
    ;; ---------------------------
    curr_idx = 0
    curr_idx_debug = 0
    total_no_of_pages = dimsizes(figures_per_page)

    do curr_page = 0, total_no_of_pages - 1
        ;;---------------------------
        ;; Plot arrays for gsn_panels
        ;;---------------------------
        plot_array = new((/max_figures_pp/), graphic)

        no_figures_on_this_page = figures_per_page(curr_page)

        ;; Create a string to add to the figure output\
            ;; filename for mulitple pages
            if (total_no_of_pages .gt. 1) then
                page_no = "-page" + sprinti("%i", curr_page)
            else
                page_no = ""
            end if

        ;; -----------------------------------
        ;; Define output workstation for plots
        ;; -----------------------------------
        idx_mod = -1  ; No specific model defined
        aux_title_info = di@season + "-" + storage_name + type_specifier + page_no
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0 + var1,\
                                                        field_type0 + field_type1,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(output_file_type, output_file_path)

        cols = multipanel_get_no_cols(no_figures_on_this_page, max_cols)
        rows = multipanel_get_no_rows(no_figures_on_this_page, max_cols)

        plot_multipanel(cols,\
                        rows,\
                        curr_idx,\
                        curr_page,\
                        res,\
                        storage_name,\
                        storage_vault,\
                        wks,\
                        di,\
                        plot_settings,\
                        valid_statistics,\
                        False,\  ; No debuginfo
                        figures_per_page,\
                        model_panel_placement,\
                        figure_panel_placement,\
                        plot_array,\
                        type_specifier,\
                        no_figures_on_this_page)
        frame(wks)

        if (debuginfo) then
            aux_title_info = di@season + "-" + storage_name + type_specifier + page_no + "-debug"
            output_filename = interface_get_figure_filename(diag_script_base,\
                                                            var0 + var1,\
                                                            field_type0 + field_type1,\
                                                            aux_title_info,\
                                                            idx_mod)
            output_file_path = output_dir + output_filename
            wks_debug = gsn_open_wks(output_file_type, output_file_path)

            plot_multipanel(cols,\
                            rows,\
                            curr_idx_debug,\
                            curr_page,\
                            res,\
                            storage_name,\
                            storage_vault,\
                            wks_debug,\
                            di,\
                            plot_settings,\
                            valid_statistics,\
                            debuginfo,\
                            figures_per_page,\
                            model_panel_placement,\
                            figure_panel_placement,\
                            plot_array,\
                            type_specifier,\
                            no_figures_on_this_page)
            frame(wks_debug)
        end if
    end do  ; curr_fig = 1, figures_per_page(curr_page) - 1
end

undef("plot_multipanel_vector")
procedure plot_multipanel_vector(cols[*]:integer,
                                 rows[*]:float,
                                 curr_idx[1]:integer,
                                 curr_page[1]:integer,
                                 res[1]:logical,
                                 storage_name[1]:string,
                                 storage_vault[1]:logical,
                                 wks[1]:graphic,
                                 di[1]:logical,
                                 plot_settings[1]:logical,
                                 valid_statistics[*]:string,
                                 debuginfo[1]:logical,
                                 figures_per_page[*]:integer,
                                 model_panel_placement[*]:integer,
                                 figure_panel_placement[*]:integer,
                                 plot_array[*]:graphic,
                                 type_specifier[1]:string,
                                 no_figures_on_this_page[1]:integer)
begin
    ;; Update position, labelbar and title
    curr_figure_pos = figure_panel_placement(curr_idx)
    res = panel_n_by_cols(res, curr_figure_pos, rows, cols, figures_per_page(curr_page))

    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    main_title_prefix = ""
    main_title_suffix = ""

    ;; Fetch reference plot
    storage_record = str_join((/"mean", "ua", sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    ua_ref = retrieve_from_vault(storage_vault, storage_record)

    storage_record = str_join((/"mean", "va", sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    va_ref = retrieve_from_vault(storage_vault, storage_record)

    storage_record = str_join((/storage_name, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
    speed_ref = retrieve_from_vault(storage_vault, storage_record)
    max_speed_ref = max(speed_ref)

    if (isatt(plot_settings, "idx")) then
        delete(plot_settings@idx)
    end if
    plot_settings@idx = plot_settings@idx_ref

    dummy_array = (/1e+20/)
    dummy_array@_FillValue = 1e+20

    statistics = True
    statistics = compute_stat((/"yrs", "mean"/), valid_statistics, speed_ref, dummy_array)
    plot_ref = single_contour(wks,\
                              speed_ref,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)
    delete(statistics)
    if (storage_name .ne. "stddev") then
        plot_ref_v = single_vector(wks,\
                                   max_speed_ref,\
                                   speed_ref,\
                                   ua_ref,\
                                   va_ref,\
                                   main_title_prefix,\
                                   main_title_suffix,\
                                   plot_settings,\
                                   debuginfo,\
                                   res)
        overlay(plot_ref, plot_ref_v)
    end if

    txres = True
    txres@txFuncCode = "~"
    if (debuginfo) then
        debugbox = write_info(debuginfo)
        am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_ref)
    end if

    delete(res@cnLevels)
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + type_specifier + "_basic", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    idx_fig = figure_panel_placement(curr_idx)
    plot_array(idx_fig) = plot_ref

    ;; Skip past the reference plot
    curr_idx = curr_idx + 1

    lbres = True
    txres = True
    txres@txFuncCode = "~"

    if (plot_settings@type .eq. "diff") then
        inset_top_text(wks, plot_ref, "REF", txres)
        inset_labelbar(wks, plot_ref, res, "REF", lbres)
        main_title_suffix = " - REF"
    else
        inset_top_text(wks, plot_ref, "Reference", txres)
        main_title_suffix = ""
    end if

    ;; ------------------------------
    ;; Create the non-reference plots
    ;; ------------------------------
    do curr_fig = 1, figures_per_page(curr_page) - 1

        main_title_prefix = ""
        idx_mod = model_panel_placement(curr_idx)
        idx_fig = figure_panel_placement(curr_idx)

        ;; Update placement and labelbar colors
        res = panel_n_by_cols(res, figure_panel_placement(curr_idx), rows, cols, figures_per_page(curr_page))

        storage_record = str_join((/"mean", "ua", sprinti("%i", idx_mod)/), str_vault_sep)
        ua = retrieve_from_vault(storage_vault, storage_record)

        storage_record = str_join((/"mean", "va", sprinti("%i", idx_mod)/), str_vault_sep)
        va = retrieve_from_vault(storage_vault, storage_record)

        storage_record = str_join((/storage_name, sprinti("%i", idx_mod)/), str_vault_sep)
        speed = retrieve_from_vault(storage_vault, storage_record)

        statistics = True
        statistics = compute_stat((/"yrs"/), valid_statistics, speed, dummy_array)
        if (plot_settings@type .eq. "diff") then
            ;; Plot mean of differences of first and second data set, first
            ;; interpolate to the reference data set grid (bilinear interpolation)

            statistics = compute_stat((/"rmse"/), valid_statistics, speed, speed_ref)

            ;; ua/va/speed field interpolation
            diff_model_ref = get_model_minus_ref(speed, speed_ref)
            diff_model_ua_ref = get_model_minus_ref(ua, ua_ref)
            diff_model_va_ref = get_model_minus_ref(va, va_ref)

            delete(speed)
            speed = diff_model_ref

            delete(ua)
            ua = diff_model_ua_ref

            delete(va)
            va = diff_model_va_ref

            delete(diff_model_ref)
            delete(diff_model_ua_ref)
            delete(diff_model_va_ref)
        end if

        if (isatt(plot_settings, "idx")) then
            delete(plot_settings@idx)
        end if
        plot_settings@idx = idx_mod

        statistics = compute_stat((/"mean"/), valid_statistics, speed, dummy_array)
        statistics = compute_stat((/"corr"/), valid_statistics, speed, speed_ref)

        plot = single_contour(wks,\
                              speed,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)
        delete(statistics)
        if (storage_name .ne. "stddev") then
            plot_v = single_vector(wks,\
                                   max_speed_ref,\
                                   speed,\
                                   ua,\
                                   va,\
                                   main_title_prefix,\
                                   main_title_suffix,\
                                   plot_settings,\
                                   debuginfo,\
                                   res)
            overlay(plot, plot_v)
        end if

        if (debuginfo) then
            debugbox = write_info(debuginfo)
            am_infobox_id = place_debuginfo(wks, debugbox, txres, plot)
        end if

        plot_array(idx_fig) = plot

        ;; Update index to point to next field
        curr_idx = curr_idx + 1
        delete(ua)
        delete(va)
        delete(speed)

    end do  ; curr_fig=1, figures_per_page(curr_page) - 1

    plottype_lbres = False

    ;; ---------------------------------------------------------
    ;; Create an blank plot for shared labelbar placement (mean)
    ;; ---------------------------------------------------------
    blank_plot = add_blank_plot_title(wks,\
                                      di@season + "-" + plot_settings@type + plot_settings@part_of_header,\
                                      rows,\
                                      cols)
    ;; Create shared labelbar
    n_by_cols_labelbar(wks,\
                       blank_plot,\
                       plot_array(no_figures_on_this_page - 1),\
                       rows,\
                       cols,\
                       plot_settings@lb_units,\
                       plottype_lbres)

    ;;---------------------
    ;; Draw mean value plot
    ;;---------------------
    if (debuginfo) then
        drawNDCGrid(wks)
    end if
    draw(plot_array)
    draw(blank_plot)
    if (debuginfo) then
        place_description(wks,\
                          debuginfo@description,\
                          debuginfo@description_ycoord)
    end if
end

undef("multipanel_vector")
procedure multipanel_vector(storage_vault[1]: logical,
                            di[1]: logical,
                            plot_settings[1]: logical,
                            storage_name[1]: string,
                            debuginfo[1]: logical,
                            valid_statistics[*]: string,
                            res[1]: logical)
; Description:
;   Place the sets of three 2D-model fields (vector_x, vector_y and abs)
;   stored in the 'storage_vault' on a set of paneled figures.
;   First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local aux_title_info,          blank_plot,        cn_levels_string,\
    cols,                    curr_fig,            curr_figure_pos,\
    curr_idx,                curr_page,           diag_script_base,\
    dim_MOD,                 dim_VAR,             dummy_array,\
    field_type0,             field_type1,         figure_panel_placement,\
    figures_per_page,        idx_fig,             idx_mod,\
    lbres,                   main_title_prefix,\
                                                  model_panel_placement,\
    no_figures_on_this_page, output_dir,          output_filename,\
    output_file_path,                             page_no,\
    plot,                    plot_array,          plot_ref,\
    plottype_lbres,          res,                 rows,\
    speed,                   speed_ref,           storage_record,\
    total_no_of_pages,       txres,               type_specifier,\
    ua,                      ua_ref,              va,\
    var0,                    var1,                va_ref,\
    wks,                     plot_ref_v
begin
    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    if (dimsizes(variables) .gt. 1) then
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(plot_settings@diag_script)
    output_dir = config_user_info@plot_dir

    ;; --------------------------------
    ;; Static resources for these plots
    ;; --------------------------------
    res@mpFillOn = False
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False
    res@gsnFrame = False
    res@gsnDraw = False
    res@lbLabelBarOn = False
    res@gsnAddCyclic = False

    ;; --------------------------------------
    ;; Compute the layout of paneled figures
    ;; --------------------------------------
    figures_per_page = get_figures_per_page(dim_MOD,\
                                            max_figures_pp,\
                                            min_figures_pp)

    ;; Which model goes where across all pages
    model_panel_placement = new((/sum(figures_per_page)/), integer)

    ;; Which model goes where on each page?
    figure_panel_placement = new((/sum(figures_per_page)/), integer)
    place_models_on_pages(model_info,\
                          plot_settings@idx_ref,\
                          figures_per_page,\
                          model_panel_placement,\
                          figure_panel_placement)

    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir

    if (plot_settings@type .eq. "diff") then
        type_specifier = "_diff"
    else
        type_specifier = ""
    end if

    ;; ---------------------------
    ;; Loop over all output pages
    ;; ---------------------------
    curr_idx = 0
    curr_idx_debug = 0
    total_no_of_pages = dimsizes(figures_per_page)

    do curr_page = 0, total_no_of_pages - 1
        ;;---------------------------
        ;; Plot arrays for gsn_panels
        ;;---------------------------
        plot_array = new((/max_figures_pp/), graphic)

        no_figures_on_this_page = figures_per_page(curr_page)

        ;; Create a string to add to the figure output\
            ;; filename for mulitple pages
            if (total_no_of_pages .gt. 1) then
                page_no = "-page" + sprinti("%i", curr_page)
            else
                page_no = ""
            end if

        ;; -----------------------------------
        ;; Define output workstation for plots
        ;; -----------------------------------
        idx_mod = -1  ; No specific model defined
        aux_title_info = di@season + "-" + storage_name + type_specifier + page_no
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0 + var1,\
                                                        field_type0 + field_type1,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(output_file_type, output_file_path)

        cols = multipanel_get_no_cols(no_figures_on_this_page, max_cols)
        rows = multipanel_get_no_rows(no_figures_on_this_page, max_cols)

        plot_multipanel_vector(cols,\
                               rows,\
                               curr_idx,\
                               curr_page,\
                               res,\
                               storage_name,\
                               storage_vault,\
                               wks,\
                               di,\
                               plot_settings,\
                               valid_statistics,\
                               False,\  ; No debuginfo
                               figures_per_page,\
                               model_panel_placement,\
                               figure_panel_placement,\
                               plot_array,\
                               type_specifier,\
                               no_figures_on_this_page)
        frame(wks)

        if (debuginfo) then
            aux_title_info = di@season + "-" + storage_name + type_specifier + page_no + "-debug"
            output_filename = interface_get_figure_filename(diag_script_base,\
                                                            var0 + var1,\
                                                            field_type0 + field_type1,\
                                                            aux_title_info,\
                                                            idx_mod)
            output_file_path = output_dir + output_filename
            wks_debug = gsn_open_wks(output_file_type, output_file_path)

            plot_multipanel_vector(cols,\
                                   rows,\
                                   curr_idx_debug,\
                                   curr_page,\
                                   res,\
                                   storage_name,\
                                   storage_vault,\
                                   wks_debug,\
                                   di,\
                                   plot_settings,\
                                   valid_statistics,\
                                   debuginfo,\
                                   figures_per_page,\
                                   model_panel_placement,\
                                   figure_panel_placement,\
                                   plot_array,\
                                   type_specifier,\
                                   no_figures_on_this_page)
            frame(wks_debug)
        end if
    end do
end

undef("seasonal_plot")
procedure seasonal_plot(storage_vault[1]: logical,
                        di[1]: logical,
                        plot_settings[1]: logical,
                        storage_name[1]: string,
                        debuginfo[1]: logical)
; Description:
;   Produces a bar plot for seasonal climatology or interannual
;   variability. The various logcials carry switches for the
;   plot routine behaviour.
local am_infobox_id,             am_labelbar_id,         amres,\
    debugbox,                  diag_script_base,       dim_VAR,\
    field_type0,               idx_mod,                labelbar,\
    labels,                    lbres,                  output_dir,\
    output_filename,           output_file_path,       plot_array,\
    precip,                    precip_seasonal_maxY,   precip_seasonal_minY,\
    precip_seasonal_spreadY,   res,                    storage_record,\
    txres,                     var0,                   wks,\
    xaxis,                     var1,                   field_type1
begin
    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    if (dim_VAR .gt. 1) then  ; Assume two variables?
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    storage_record = storage_name
    precip = retrieve_from_vault(storage_vault, storage_record)
    units = precip@units

    precip_seasonal_maxY = max(precip)
    precip_seasonal_minY = min(precip)
    precip_seasonal_spreadY = precip_seasonal_maxY - precip_seasonal_minY

    ;; Plot arrays for gsn_plots
    plot_array  = new((/dim_MOD/), graphic)

    ;;--------------------------------
    ;; General resources for all plot
    ;;--------------------------------
    res                       = True
    res@gsnDraw               = False
    res@gsnFrame              = False
    res@gsnXYBarChart         = True
    res@gsnXYBarChartBarWidth = 0.15  ; change bar widths
    res@gsnYRefLine           = 0.    ; reference line

    ;; Specific plot resource for plot iii)
    if (plot_settings@type .eq. "iav") then
        res@tmXBFormat = "4f"
        xaxis = di@years
    else
        ;; Set months on x-axis for plots i) and ii)
        res@tmXBMode      = "Explicit"       ; Define your own tick mark labels.
        res@tmXBMinorOn   = False            ; No minor tick marks.
        res@tmXBValues    = ispan(0, 11, 1)  ; Location to put tick mark labels
        res@tmXBLabels    = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun",\
                              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/)
        res@tmXBLabelFont = 25            ; Change font of labels.
        res@tmXBLabelFontHeightF = 0.015  ; Change font height of labels.

        ;; Add some margins
        res@trXMinF = -1
        res@trXMaxF = 12

        xaxis = fspan(0, 11, 12)
    end if

    labels = new(dim_MOD, string)

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir

    ;;-----------------------------------------------
    ;; Define output workstation for seasonal anomaly
    ;;-----------------------------------------------
    idx_mod = -1  ; No specific model defined
    aux_filename_info = plot_settings@aux_filename_info
    output_filename = interface_get_figure_filename(diag_script_base,\
                                                    var0 + var1,\
                                                    field_type0 + field_type1,\
                                                    aux_filename_info,\
                                                    idx_mod)
    output_file_path = output_dir + output_filename
    wks = gsn_open_wks(output_file_type, output_file_path)

    ;;----------
    ;; Plot loop
    ;;----------
    do idx_mod = 0, dim_MOD - 1
        labels(idx_mod) = model_info@model(idx_mod)
        res@gsnXYBarChartColors = di@colors_seasonal(idx_mod)

        ;; Shift plot to the right and draw/redraw
        res@vpXF = 0.2 + idx_mod * 0.005

        res@trYMinF = precip_seasonal_minY - precip_seasonal_spreadY * 0.05
        res@trYMaxF = precip_seasonal_maxY + precip_seasonal_spreadY * 0.05

        ;; Set strings for titles and axis for plot i) and ii)
        res@tiMainString = set_string_once(plot_settings@title_string, idx_mod)
        res@tiYAxisString = set_string_once(plot_settings@yaxis_string, idx_mod)
        res@gsnLeftString = set_string_once(plot_settings@left_plot_subheader, idx_mod)
        res@gsnRightString = set_string_once(units, idx_mod)

        if (isatt(plot_settings, "xaxis_string")) then
            res@tiXAxisString = set_string_once(plot_settings@xaxis_string, idx_mod)
        end if

        plot_array(idx_mod) = gsn_csm_xy(wks,\
                                         xaxis,\
                                         precip(idx_mod, :),\
                                         res)

        ;; Remove borders such that next plot only draws the actual "bars"
        res = remove_figure_borders(res)
    end do  ; idx_mod = 0, dim_MOD - 1

    ;; Label bar resources
    lbres                   = True             ; labelbar only resources
    lbres@vpWidthF          = 0.30             ; labelbar width
    lbres@vpHeightF         = 0.024 * dim_MOD  ; labelbar height
    lbres@lbBoxMajorExtentF = 0.36             ; puts space between color boxes
    lbres@lbFillColors      = di@colors_seasonal
    lbres@lbMonoFillPattern = True             ; Solid fill pattern
    lbres@lbLabelJust       = "CenterLeft"     ; left justify labels
    lbres@lbPerimOn         = True
    lbres@lbPerimFill       = 0
    lbres@lbPerimFillColor  = "white"
    lbres@lbPerimColor      = "black"
    labelbar = gsn_create_labelbar(wks, dim_MOD, labels, lbres)

    ;; Place annotations
    amres                = True
    amres@amZone         = 3
    amres@amSide         = "Right"
    amres@amParallelPosF = 0.88
    amres@amOrthogonalPosF = -0.25
    am_labelbar_id = gsn_add_annotation(plot_array(dim_MOD - 1), labelbar, amres)

    ;; If requested, add a text info box and/or debug info box to each plot
    txres            = True
    txres@txFuncCode = "~"

    if (debuginfo) then
        debuginfo@years  = di@years_string
        debugbox = write_info(debuginfo)
        am_infobox_id = place_debuginfo(wks,\
                                        debugbox,\
                                        txres,\
                                        plot_array(dim_MOD - 1))
        drawNDCGrid(wks)
    end if

    draw(wks)
    if (debuginfo) then
        place_description(wks,\
                          debuginfo@description,\
                          debuginfo@description_ycoord)
    end if
    frame(wks)
end

undef("xy_plot_wrapper")
procedure xy_plot_wrapper(storage_vault[1]: logical,
                                     di[1]: logical,
                          plot_settings[1]: logical,
                           storage_name[1]: string,
                              debuginfo[1]: logical)
local annots,           annots_mmm,        avgstd,\
    avgstd_mmm,       colors,            colors_mmm,\
    dashes,           dashes_mmm,        diag_script_base,\
    dim_VAR,          field_type0,       field_type1,       i,  idx_mod,\
    imon,                                output_dir,\
    output_filename,  output_file_path,  res,\
    storage_record,   temp,              thicks,\
    thicks_mmm,       values,            values_stddev,\
    values_stddev,    var0,              var1,\
    wks,              wks_debug
begin
    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    if (dim_VAR .gt. 1) then  ; Assume two variables?
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    storage_record = storage_name
    values = retrieve_from_vault(storage_vault, storage_record)
    units = values@units

    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir
    ;;-----------------------------------------------
    ;; Define output workstation for xy plot
    ;;-----------------------------------------------
    idx_mod = -1  ; No specific model defined
    aux_filename_info = plot_settings@aux_filename_info
    output_filename = interface_get_figure_filename(diag_script_base,\
                                                    var0 + var1,\
                                                    field_type0 + field_type1,\
                                                    aux_filename_info,\
                                                    idx_mod)
    output_file_path = output_dir + output_filename
    wks = gsn_open_wks(output_file_type, output_file_path)

    aux_filename_info = plot_settings@aux_filename_info + "-debug"
    output_filename = interface_get_figure_filename(diag_script_base,\
                                                    var0 + var1,\
                                                    field_type0 + field_type1,\
                                                    aux_filename_info,\
                                                    idx_mod)
    output_file_path = output_dir + output_filename
    if (debuginfo) then
        wks_debug = gsn_open_wks(output_file_type, output_file_path)
    end if

    values&month@units = "Multi-year monthly mean"

    ;; Select colors and other plotting attributes
    ;; See ./diag_scripts/shared/plot/style.ncl
    colors = project_style(diag_script_info, "colors")
    dashes = project_style(diag_script_info, "dashes")
    thicks = project_style(diag_script_info, "thicks")
    annots = project_style(diag_script_info, "annots")
    avgstd = project_style(diag_script_info, "avgstd")
    ;; Select colors and other plotting attributes for multi-model mean
    if (diag_script_info@multi_model_mean .eq. "y") then
        ;; project_style evaluates metadata of variable "models"
        temp = models

        ;; -> keep original "models" in "temp" and restore later
        copy_VarMeta(models, temp)
        delete(models)

        ;; Use "models" to pass on attribute names
        models = getvaratts(temp)
        do i = 0, dimsizes(models) - 1
            ;; Define all original attributes again, but empty
            model_info@$models(i)$ = ""
        end do
        model_info@model = "model_mean"

        ;; See ./diag_scripts/shared/plot/style.ncl
        colors_mmm = project_style(diag_script_info, "colors")
        dashes_mmm = project_style(diag_script_info, "dashes")
        thicks_mmm = project_style(diag_script_info, "thicks")
        annots_mmm = project_style(diag_script_info, "annots")
        avgstd_mmm = project_style(diag_script_info, "avgstd")
        delete(models)
        models = temp   ; Restore original "models"
        copy_VarMeta(temp, models)
        delete(temp)
    end if

    ;; Calculate standard deviation of models
    if (di@multi_model_mean .ne. "y") then
        ;; define anyway, because fields are needed as parameters for xy_line
        values_stddev  = 0
    else
        ;; Mean, stddev, -1 * stddev, +1 * stddev
        values_stddev = new((/4, 12/), float)
        ;; See ./diag_scripts/shared/plot/style.ncl for which data shall be
        ;; included in the statistics
        temp = ind(avgstd .eq. 0)
        do imon = 0, 11
            values_stddev(0, imon) = \
                dim_avg_n_Wrap(values(temp, imon), (/0/))
            values_stddev(1, imon) = \
                dim_stddev_n_Wrap(values(temp, imon), (/0/))
            values_stddev(2, imon) = \
                values_stddev(0, imon) - values_stddev(1, imon)
            values_stddev(3, imon) = \
                values_stddev(0, imon) + values_stddev(1, imon)
        end do
        delete(temp)
    end if

    ;; Min and max values on x-axis
    res = True
    res@tmXTOn = False
    res@tmYROn = True
    if (di@supporting_gridlines .eq. "y") then
        res@tmXMajorGrid = True
        res@gsnYRefLine = 0
    end if
    res@trXMinF = min(values&month) - 0.05 * \
        (max(values&month) - min(values&month))
    res@trXMaxF = max(values&month) + 0.25 * \
        (max(values&month) - min(values&month))

    res@tmXBMode        = "Explicit"
    res@tmXBValues      = (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11/)
    res@tmXBLabels      = (/"J", "F", "M", "A", "M", "J",\
                           "J", "A", "S", "O", "N", "D"/)

    res@tiMainString = plot_settings@title_string
    res@tiYAxisString = plot_settings@yaxis_string
    res@gsnRightString = values@units
    res@gsnRightStringFontHeightF = 16
    res@gsnLeftString = plot_settings@left_plot_subheader
    res@gsnLeftStringFontHeightF = 16

    xy_line(wks, values, values&month, values_stddev, res, False)

    if (debuginfo) then
        xy_line(wks_debug, values, values&month, values_stddev, res, debuginfo)
    end if
end

undef("ts_line_wrapper")
procedure ts_line_wrapper(storage_vault[1]: logical,
                                     di[1]: logical,
                          plot_settings[1]: logical,
                           storage_name[1]: string,
                              debuginfo[1]: logical)
local    annots,           avgstd,            diag_script_base,\
    dim_VAR,          field_type0,       field_type1,\
    idx_mod,          imonth,            output_dir,\
    output_filename,  output_file_path,  res,\
    storage_record,   temp,              val_area_stddev,\
    values,           values_stddev,     var0,\
    var1,             wks,               wks_debug
begin
    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    if (dim_VAR .gt. 1) then  ; Assume two variables?
        var1 = variables(1)
        field_type1 = field_types(1)
    else
        var1 = ""
        field_type1 = ""
    end if

    storage_record = storage_name
print(storage_vault)
    values = retrieve_from_vault(storage_vault, storage_record)
    units = values@units

    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir
    ;;-----------------------------------------------
    ;; Define output workstation for xy plot
    ;;-----------------------------------------------
    idx_mod = -1  ; No specific model defined
    aux_filename_info = plot_settings@aux_filename_info
    output_filename = interface_get_figure_filename(diag_script_base,\
                                                    var0 + var1,\
                                                    field_type0 + field_type1,\
                                                    aux_filename_info,\
                                                    idx_mod)
    output_file_path = output_dir + output_filename
    wks = gsn_open_wks(output_file_type, output_file_path)

    aux_filename_info = plot_settings@aux_filename_info + "-debug"
    output_filename = interface_get_figure_filename(diag_script_base,\
                                                    var0 + var1,\
                                                    field_type0 + field_type1,\
                                                    aux_filename_info,\
                                                    idx_mod)
    output_file_path = output_dir + output_filename
    if (debuginfo) then
        wks_debug = gsn_open_wks(output_file_type, output_file_path)
    end if

    ;; Determine time range
    start_year = min(stringtoint(model_info@start_year))
    end_year = max(stringtoint(model_info@end_year))

    ;; Select attributes
    annots = project_style(diag_script_info, "annots")
    avgstd = project_style(diag_script_info, "avgstd")

    ;; Prepare calculation of standard deviation of models
    if (diag_script_info@multi_model_mean .eq. "y") then
        values_stddev = new((/4, end_year - start_year + 1/), float)
        val_area_stddev = new((/4, end_year - start_year + 1/), float)
    else  ; initialize anyway, because both are needed as parameters in call of xy_line
        values_stddev = 0
        val_area_stddev = 0
    end if

    ;; Calculate standard deviation of models
    if (diag_script_info@multi_model_mean .eq. "y") then
        ;; See style_$project$.ncl for which data shall
        ;; be included in the statistics
        temp = ind(avgstd .eq. 0)
        do imonth = 0, dimsizes(values&years) - 1
           ;; For extent
            values_stddev(0, imonth) = \
                dim_avg_n_Wrap(values(temp, imonth), (/0/))
            values_stddev(1, imonth) = \
                dim_stddev_n_Wrap(values(temp, imonth), (/0/))
            values_stddev(2, imonth) = \
                values_stddev(0, imonth) - values_stddev(1, imonth)
            values_stddev(3, imonth) = \
                values_stddev(0, imonth) + values_stddev(1, imonth)
        end do
        delete(temp)
    end if

    res = True
    res@tmXTOn = False
    res@tmYROn = True
    if (di@supporting_gridlines .eq. "y") then
        res@tmXMajorGrid = True
        res@gsnYRefLine = 0
    end if
    res@tmXBFormat = "4f"
    res@tiMainString = plot_settings@title_string
    res@tiYAxisString = plot_settings@yaxis_string
    res@tiXAxisString = plot_settings@xaxis_string
    res@gsnRightString = values@units
    res@gsnRightStringFontHeightF = 16
    res@gsnLeftString = plot_settings@left_plot_subheader
    res@gsnLeftStringFontHeightF = 16
    res@trXMinF = min(di@years) - 0.05 * (max(di@years) - \
                                          min(di@years))
    res@trXMaxF = max(di@years) + 0.25 * (max(di@years) - \
                                          min(di@years))

    xy_line(wks, values, values&years, values_stddev, res, False)

    if (debuginfo) then
        xy_line(wks_debug, values, di@years, values_stddev, res, debuginfo)
    end if
end
