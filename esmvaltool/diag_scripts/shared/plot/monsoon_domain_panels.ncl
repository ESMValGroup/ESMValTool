;;#############################################################################
;; GENERAL ROUTINES FOR MONSOON PANELLING
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;    plot_precip_domain
;;    precip_domain
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/logging.ncl"

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip_world.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon_panel_positions.ncl"
load "./diag_scripts/aux/SAMonsoon/exact_panel_positions_generic.ncl"
load "./diag_scripts/shared/plot/style.ncl"

begin
    loadscript(exact_panel_position_file)
end

;;#############################################################################
undef("plot_precip_domain")
procedure plot_precip_domain(cols[*]:integer,
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
;;
;; Arguments
;;                       cols: number of columns for this panel plot
;;                       rows: number of rows for this panel plot
;;                   curr_idx: current index
;;                  curr_page: current page (may be more than one)
;;                        res: valid ncl resources to apply to plot
;;               storage_name: name of first processed variable in the vault
;;              storage_vault: handle/pointer to vault with all models
;;                        wks: workstation
;;                         di: logical with diagnostic script info
;;              plot_settings: logical with general plot_settings
;;           valid_statistics: valid statistics (to be computed if defined)
;;                  debuginfo: info to put onto plot in debug mode
;;           figures_per_page: array with number of figures on each page
;;      model_panel_placement: where to place respective model
;;     figure_panel_placement: where to place respective figure on the page
;;                 plot_array: plot handles/pointers
;;             type_specifier: kind of plot, 'mean' or 'stddev'
;;    no_figures_on_this_page: no of figures on this page
;;
;; Return value
;;
;; Description
;;    Multipanel plot, plots all models on the current page. Top left entry
;;    is always the reference model. Used for the global domain plots only.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150702-A_eval_ma: written.
;;
local am_infobox_id,      blank_plot,         cn_levels_string,\
    curr,               curr_idx,           debugbox,\
    diff_model_ref,     dummy_array,        idx_fig,\
    idx_mod,            lb_prefix,          lbres,\
    main_title_prefix,  main_title_suffix,  plot,\
    plot_ref,           plottype_lbres,     ref,\
    statistics,         storage_record,     txres
begin

    funcname = "plot_precip_domain"
    scriptname = "plot_script/ncl/monsoon_domain_panels.ncl"
    enter_msg(scriptname, funcname)

    if (isatt(res, "cnLevels")) then
        delete(res@cnLevels)
    end if
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_mean", "-", "_")
    res@cnLevels = di@$cn_levels_string$

    if (isatt(plot_settings, "default_palette")) then
        if (isatt(res, "cnFillPalette")) then
            delete(res@cnFillPalette)
        end if
        res@cnFillPalette = plot_settings@default_palette
    end if

    main_title_prefix = ""
    main_title_suffix = ""

    storage_record = str_join((/storage_name, sprinti("%i", plot_settings@idx_ref)/), str_vault_sep)
    ref = retrieve_from_vault(storage_vault, storage_record)

    ;; Reference plot
    plot_settings@idx = plot_settings@idx_ref

    ;; Dummy array required for some plot routines
    dummy_array = (/(/1e20, 1e20/), (/1e20, 1e20/)/)
    dummy_array@_FillValue = 1e20

    statistics = True
    statistics = compute_stat((/"yrs"/), valid_statistics, ref, dummy_array)

    res = panel_n_by_cols(res, figure_panel_placement(curr_idx), rows, cols, figures_per_page(curr_page))
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
        if (storage_name .eq. "intensity") then
            inset_labelbar(wks, plot_ref, res, "REF", lbres)
        end if
    end if

    delete(res@cnLevels)
    cn_levels_string = str_sub_str("cn_levels_" + storage_name + type_specifier, "-", "_")
    res@cnLevels = di@$cn_levels_string$

    if (isatt(plot_settings, "diff_palette")) then
        if (isatt(res, "cnFillPalette")) then
            delete(res@cnFillPalette)
        end if
        res@cnFillPalette = plot_settings@diff_palette
    end if

    ;; -----------------------------
    ;; Create the non-reference plots
    ;; -----------------------------
    do curr_fig = 1, figures_per_page(curr_page) - 1

        main_title_prefix = ""
        idx_mod = model_panel_placement(curr_idx)
        idx_fig = figure_panel_placement(curr_idx)

        ;; Update placement and labelbar colors
        res = panel_n_by_cols(res, figure_panel_placement(curr_idx), rows, cols, figures_per_page(curr_page))

        storage_record = str_join((/storage_name, sprinti("%i", idx_mod)/), str_vault_sep)
        curr = retrieve_from_vault(storage_vault, storage_record)
        if (plot_settings@type .eq. "diff") then
            ;; Plot mean of differences of first and second data set, first
            ;; interpolate to the reference data set grid (bilinear interpolation)
            diff_model_ref = get_model_minus_ref(curr, ref)

            delete(curr)
            curr = diff_model_ref
            delete(diff_model_ref)

            main_title_suffix = " - REF"
        else
            main_title_suffix = ""
        end if

        statistics = True
        statistics = compute_stat((/"yrs"/), valid_statistics, curr, dummy_array)

        plot_settings@idx = idx_mod
        plot = single_contour(wks,\
                              curr,\
                              main_title_prefix,\
                              main_title_suffix,\
                              plot_settings,\
                              debuginfo,\
                              statistics,\
                              res)
        delete(statistics)
        if (debuginfo) then
            debugbox = write_info(debuginfo)
            am_infobox_id = place_debuginfo(wks, debugbox, txres, plot)
        end if

        plot_array(idx_fig) = plot
        delete(curr)

        ;; Update index to point to next field
        curr_idx = curr_idx + 1
    end do  ; curr_fig=1, figures_per_page(curr_page) - 1

    plottype_lbres = True

    ;; Copy over any attributes with this prefix
    lb_prefix = "lbres_"
    plottype_lbres = filter_attrs(plot_settings, plottype_lbres, lb_prefix)

    ;; Create an blank plot for shared labelbar placement
    blank_plot = add_blank_plot_title(wks,\
                                      "Monsoon " + plot_settings@part_of_header,\
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

    leave_msg(scriptname, funcname)

end

;;#############################################################################
undef("precip_domain")
procedure precip_domain(storage_vault [1] : logical,
                     di [1] : logical,
                     plot_settings [1] : logical,
                     storage_name [1] : string,
                     debuginfo [1] : logical,
                     valid_statistics [*] : string,
                     res [1] : logical)
;;
;; Arguments
;;                  storage_vault: handle/pointer to vault with all models
;;                             di: logical with diagnostic script info
;;                  plot_settings: logical with general plot_settings
;;                   storage_name: name of first processed variable in the vault
;;                      debuginfo: info to put onto plot in debug mode
;;               valid_statistics: valid statistics (to be computed if defined)
;;                            res: valid ncl resources to apply to plot
;;
;; Return value
;;
;; Description
;;    Determines how to place a number of contour plots in a grid across
;;    multiple pages. Loop over pages and call plot_multipanel(...) for
;;    each page to plot entries. Used for the domain plots only.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20150702-A_eval_ma: written.
;;
local aux_filename_info,         blank_plot,          cn,\
    cols,                      curr,                curr_fig,\
    curr_figure_pos,           curr_idx,            curr_page,\
    diag_script_base,          dim_MOD,             dim_VAR,\
                                                    field_type0,\
    figure_panel_placement,    figures_per_page,    fill_colors_string,\
    idx_fig,                   idx_mod,             idx_ref,\
    lbres,                     main_title_prefix,   max_cols,\
                                                    model_panel_placement,\
    no_figures_on_this_page,   output_dir,          output_filename,\
    output_file_path,                               page_no,\
    plot,                      plot_array,          plot_ref,                plottype_lbres,\
    ref,                       res,                 rows,\
    storage_record,            total_no_of_pages,   txres,\
    type_specifier,            var0,                wks
begin

    funcname = "precip_domain"
    scriptname = "plot_script/ncl/monsoon_domain_panels.ncl"
    enter_msg(scriptname, funcname)

    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir

    res@tiMainFontHeightF = 0.011

    ;;-------------------------------
    ;; Static resources for this plot
    ;;-------------------------------
    res@mpFillOn = False
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False
    res@gsnFrame = False
    res@gsnDraw = False
    res@lbLabelBarOn = False
    res@gsnAddCyclic = True
    ;res@cnSpanFillPalette = False

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

    type_specifier = "_" + plot_settings@type

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

        ;; ----------------------------------------
        ;; Define output workstation for mean plots
        ;; ----------------------------------------
        idx_mod = -1  ; No specific model defined
        aux_filename_info = di@summer_season + "-"\
                          + storage_name\
                          + type_specifier\
                          + page_no
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0,\
                                                        field_type0,\
                                                        aux_filename_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(output_file_type, output_file_path)

        max_cols = 2
        cols = multipanel_get_no_cols(no_figures_on_this_page, max_cols)
        rows = multipanel_get_no_rows(no_figures_on_this_page, max_cols)

        ;; Update position, labelbar and title
        curr_figure_pos = figure_panel_placement(curr_idx)
        res = panel_n_by_cols(res, curr_figure_pos, rows, cols, figures_per_page(curr_page))

        plot_precip_domain(cols,\
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
            aux_filename_info = di@summer_season + "-"\
                              + storage_name\
                              + type_specifier\
                              + page_no\
                              + "-debug"
            output_filename = interface_get_figure_filename(diag_script_base,\
                                                            var0,\
                                                            field_type0,\
                                                            aux_filename_info,\
                                                            idx_mod)
            output_file_path = output_dir + output_filename
            wks_debug = gsn_open_wks(output_file_type, output_file_path)

            plot_precip_domain(cols,\
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

    end do  ; curr_page = 0, total_no_of_pages - 1

    leave_msg(scriptname, funcname)

end
