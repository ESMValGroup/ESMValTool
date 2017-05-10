:mod:`monsoon_panels`
=====================
.. function:: get_title_suffix(plot_settings [1] : logical)

   :param  logical plot_settings: Logical with plot settings as attributes

   Return value
      main_title_suffix: Main title suffix
  
   Description
      Checks and returns the main title suffix attribute from
      plot_settings, empty string if missing
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  remove_attrs(p_settings [1]:logical, attributes [*]:string)

   :param logical p_settings: Logical with plot settings as attributes
   :param string attributes: array with attribute names

   Return value
  
   Description
      Removes supplied attributes if they are defined
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  plot_two_by_one(wks [1]:graphic, res [1]:logical, di [1]:logical, plot_settings [1]:logical, valid_statistics [*]:string, storage_name1 [1]:string, storage_name2 [1]:string, debuginfo [1]:logical, storage_vault [1]:logical, idx_mod [1]:integer)

   :param graphic wks: workstation
   :param logical res: resources for plot
   :param logical di: logical with diagnostic script info
   :param logical plot_settings: logical with general plot_settings
   :param string valid_statistics: valid statistics (to be computed if defined)
   :param string storage_name1: name of first processed variable in the vault
   :param string storage_name2: name of second processed variable in the vault
   :param logical debuginfo: info to put onto plot in debug mode
   :param logical storage_vault: handle/pointer to vault with all models
   :param integer idx_mod: index of current model

   Return value
  
   Description
      Plots two contour plots on top of each other
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  plot_three_by_one_diff(wks [1]:graphic, res [1]:logical, di [1]:logical, plot_settings [1]:logical, valid_statistics [*]:string, storage_name [1]:string, debuginfo [1]:logical, storage_vault [1]:logical, ref [*][*]:numeric, idx_mod [1]:integer)

   :param graphic wks: workstation
   :param logical res: resources for plot
   :param logical di: logical with diagnostic script info
   :param logical plot_settings: logical with general plot_settings
   :param string valid_statistics: valid statistics (to be computed if defined)
   :param string storage_name: name of first processed variable in the vault
   :param logical debuginfo: info to put onto plot in debug mode
   :param logical storage_vault: handle/pointer to vault with all models/variables idx_mod: index of current model  Return value  Description
   :param numeric ref:erence, middle model, bottom their difference  Caveats  References  Modification history 20150702-A_eval_ma: written.  local am_infobox_id,   cn_levels_string,   curr,\ debugbox,        diff_model_ref,     dummy_array,\ header_prefix,   main_title,         main_title_suffix,\ plot,            plot_diff,          plot_ref,\ res,             statistics,         storage_record,\ title_name,      txres  begin verbosity  = stringtointeger(getenv("ESMValTool_verbosity")) info_output("<<<<<<<< Entering plot_three_by_one_diff", verbosity, 8) res = panel_three_by_one(res, 0)  ; Update resource for first plot res@cnLevelSelectionMode = "ExplicitLevels" if (isatt(res, "cnLevels")) then delete(res@cnLevels) end if  dummy_array = (/1e+20/) dummy_array@_FillValue = 1e+20  ------------------------------------- Plot reference plot (top plot) ------------------------------------- if (isatt(res, "cnLevels")) then delete(res@cnLevels) end if cn_levels_string = str_sub_str("cn_levels_" + storage_name + "_basic", "-", "_") res@cnLevels = di@$cn_levels_string$  Title string might be a substring of the variable name used for storing the data in the vault title_name = keep_only_n_cols(plot_settings, storage_name) header_prefix = empty_string_if_missing_attribute(plot_settings, "prefix_to_header")  main_title_prefix = "(1) " + header_prefix + di@season + "-" + title_name + " of " main_title_suffix = get_title_suffix(plot_settings) plot_settings@type = "ref"  if (isatt(plot_settings, "idx")) then delete(plot_settings@idx) end if plot_settings@idx = plot_settings@idx_ref  Possibly switch from default colormap for non-diff plot if (isatt(plot_settings, "default_palette")) then if (isatt(res, "cnFillPalette")) then delete(res@cnFillPalette) end if res@cnFillPalette = plot_settings@default_palette end if  statistics = True statistics = compute_stat((/"yrs", "mean", "bob", "eio", "sahel", "sa34", "en34"/), valid_statistics, ref, dummy_array) plot_ref = single_contour(wks,\ ref,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ statistics,\ res) delete(statistics) txres = True inset_top_text(wks, plot_ref, "Reference", txres)  ;lbtitle = plot_settings@lb_units ;three_by_one_labelbar(wks, plot_ref, lbtitle) three_by_one_labelbar(wks, plot_ref, plot_settings)  ------------------------------------- Plot model entry (middle plot) ------------------------------------- header_prefix = empty_string_if_missing_attribute(plot_settings, "prefix_to_header")  title_name = keep_only_n_cols(plot_settings, storage_name) main_title_prefix = "(2) " + header_prefix + di@season + "-" + title_name + " of " res = panel_three_by_one(res, 1)  ; Update resource for second plot 
   :param integer idx_mod:)/), str_vault_sep) curr = retrieve_from_vault(storage_vault, storage_record) plot_settings@type = "mean"

.. function::  two_by_one(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name1 [1] : string, storage_name2 [1] : string, debuginfo [1] : logical, valid_statistics [*] : string, res [1] : logical)

   :param  logical storage_vault: handle/pointer to vault with all models
   :param  logical di: logical with diagnostic script info
   :param  logical plot_settings: logical with general plot_settings res: resources for plot
   :param  string storage_name1: name of first processed variable in the vault
   :param  string storage_name2: name of second processed variable in the vault
   :param  logical debuginfo: info to put onto plot in debug mode
   :param  string valid_statistics: valid statistics (to be computed if defined)
   :param  logical res: valid ncl resources to apply to plot

   Return value
  
   Description
      Loops models in storage_vault and calls plot_two_by_one(..)
      for each model
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  three_by_one_diff(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name [1] : string, debuginfo [1] : logical, valid_statistics [*] : string, res [1] : logical)

   :param  logical storage_vault: handle/pointer to vault with all models
   :param  logical di: logical with diagnostic script info
   :param  logical plot_settings: logical with general plot_settings res: resources for plot
   :param  string storage_name: name of first processed variable in the vault
   :param  logical debuginfo: info to put onto plot in debug mode
   :param  string valid_statistics: valid statistics (to be computed if defined)
   :param  logical res: valid ncl resources to apply to plot

   Return value
  
   Description
      Loops models in storage_vault and calls plot_three_by_one_diff(..)
      for each model
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  plot_three_by_one_vector(wks [1]:graphic, res [1]:logical, di [1]:logical, plot_settings [1]:logical, valid_statistics [*]:string, storage_names [1]:logical, debuginfo [1]:logical, storage_vault [1]:logical, ua_ref [*][*]:numeric, va_ref [*][*]:numeric, speed_ref [*][*]:numeric, idx_mod [1]:integer)

   :param graphic wks: workstation
   :param logical res: resources for plot
   :param logical di: logical with diagnostic script info
   :param logical plot_settings: logical with general plot_settings
   :param string valid_statistics: valid statistics (to be computed if defined) storage_name: name of first processed variable in the vault debuginfo: info to put onto plot in debug mode storage_vault: handle/pointer to vault with all models/variables ua_ref: ua wind reference model/obs va_ref: ua wind reference model/obs speed_ref: wind speed reference for model/obs idx_mod: index of current model  Return value  Description Plots three contour plots, top: reference, middle model, bottom their difference. On top of the contour the vector field is plotted  Caveats  References  Modification history 20150702-A_eval_ma: written.   local am_infobox_id,      cn_levels_string,   debugbox,\ diff_model_ref,     diff_model_ua_ref,  diff_model_va_ref,\ dummy_array,        main_title,         main_title_prefix,\ main_title_suffix,  max_speed_ref,      plot,\ plot_diff,          plot_diff_v,        plot_ref,\ plot_ref_v,         plot_v,             speed,\ statistics,         storage_record,     txres,\ ua,                 va begin verbosity  = stringtointeger(getenv("ESMValTool_verbosity")) info_output("<<<<<<<< Entering plot_three_by_one_vector", verbosity, 8) res = panel_three_by_one(res, 0)  ; Update resource for first plot res@cnLevelSelectionMode = "ExplicitLevels" if (isatt(res, "cnLevels")) then delete(res@cnLevels) end if  ------------------------------ Plot reference plot (top plot) ------------------------------
   :param logical storage_names:@base_name + "_basic", "-", "_") res@cnLevels = di@$cn_levels_string$  main_title_prefix = "(1) " + di@season + "-" + storage_names@base_name + " of " main_title_suffix = "" plot_settings@type = "ref"  if (isatt(plot_settings, "idx")) then delete(plot_settings@idx) end if plot_settings@idx = plot_settings@idx_ref  dummy_array = (/1e+20/) dummy_array@_FillValue = 1e+20  statistics = True statistics = compute_stat((/"yrs", "mean"/), valid_statistics, speed_ref, dummy_array) plot_ref = single_contour(wks,\ speed_ref,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\
   :param logical debuginfo:,\ statistics,\ res) delete(statistics) txres = True inset_top_text(wks, plot_ref, "Reference", txres)  if (isatt(di, "max_speed_ref")) then max_speed_ref = di@max_speed_ref else max_speed_ref = max(speed_ref) end if if (storage_names@base_name .ne. "stddev") then plot_ref_v = single_vector(wks,\ max_speed_ref,\ speed_ref,\ ua_ref,\ va_ref,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ res) overlay(plot_ref, plot_ref_v) end if three_by_one_labelbar(wks, plot_ref, plot_settings)  ------------------------------------- Plot model entry plot (middle plot) ------------------------------------- main_title_prefix = "(2) " + di@season + "-" + storage_names@base_name + " of " res = panel_three_by_one(res, 1)  ; Update resource for second plot  storage_record = str_join((/storage_names@storage_x_component, sprinti("%i", idx_mod)/), str_vault_sep)
   :param logical storage_vault:, storage_record)  storage_record = str_join((/storage_names@storage_y_component, sprinti("%i", idx_mod)/), str_vault_sep) va = retrieve_from_vault(storage_vault, storage_record)  storage_record = str_join((/storage_names@storage_xy_component, sprinti("%i", idx_mod)/), str_vault_sep) speed = retrieve_from_vault(storage_vault, storage_record) plot_settings@type = "mean"  if (isatt(plot_settings, "idx")) then delete(plot_settings@idx) end if plot_settings@idx = idx_mod  Interpolate to reference grid for pattern correlation calculation diff_model_ref = get_model_minus_ref(speed, speed_ref)  statistics = True statistics = compute_stat((/"yrs", "mean"/), valid_statistics, speed, dummy_array) statistics = compute_stat((/"corr"/), valid_statistics, speed, speed_ref) plot = single_contour(wks,\ speed,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ statistics,\ res) delete(statistics) if (storage_names@base_name .ne. "stddev") then plot_v = single_vector(wks,\ max_speed_ref,\ speed,\ ua,\ va,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ res) overlay(plot, plot_v) end if  -------------------------------------------- Plot model reference diff plot (bottom plot) -------------------------------------------- Plot mean of differnces of first and second data set, first interpolate to the reference data set grid (bilinear interpolation)  diff_model_ref = get_model_minus_ref(speed, speed_ref)
   :param numeric ua_ref: = get_model_minus_ref(ua, ua_ref)
   :param numeric va_ref: = get_model_minus_ref(va, va_ref)  main_title = "(2) - (1)" res = panel_three_by_one(res, 2)  ; Update resource for third plot delete(res@cnLevels) cn_levels_string = str_sub_str("cn_levels_" + storage_names@base_name+ "_diff_basic", "-", "_") res@cnLevels = di@$cn_levels_string$  plot_settings@type = "diff"  statistics = True statistics = compute_stat((/"mean"/), valid_statistics, diff_model_ref, dummy_array)
   :param numeric speed_ref:) plot_diff = single_contour(wks,\ diff_model_ref,\ main_title,\ main_title_suffix,\ plot_settings,\ debuginfo,\ statistics,\ res) delete(statistics) if (isatt(di, "max_speed_ref_diff")) then max_speed_ref_diff = di@max_speed_ref_diff else max_speed_ref_diff = max(diff_model_ref) end if if (storage_names@base_name .ne. "stddev") then plot_diff_v = single_vector(wks,\ max_speed_ref,\ diff_model_ref,\ diff_model_ua_ref,\ diff_model_va_ref,\ main_title,\ main_title_suffix,\ plot_settings,\ debuginfo,\ res) overlay(plot_diff, plot_diff_v) end if  if (debuginfo) then debugbox = write_info(debuginfo) am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_ref) am_infobox_id = place_debuginfo(wks, debugbox, txres, plot) am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_diff) drawNDCGrid(wks) end if draw(plot_ref) draw(plot_diff) draw(plot) if (debuginfo) then place_description(wks,\ debuginfo@description,\ debuginfo@description_ycoord) end if info_output(">>>>>>>> Leaving plot_three_by_one_vector", verbosity, 8) end  undef("three_by_one_vector") procedure three_by_one_vector(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_names [1] : logical, debuginfo [1] : logical, valid_statistics [*] : string, res [1] : logical)  Arguments storage_vault: handle/pointer to vault with all models di: logical with diagnostic script info plot_settings: logical with general plot_settings storage_name: name of first processed variable in the vault debuginfo: info to put onto plot in debug mode valid_statistics: valid statistics (to be computed if defined) res: valid ncl resources to apply to plot  Return value  Description Loops models in storage_vault and calls plot_three_by_one_vector(..) for each model producing a reference plot at the top (abs + vector) ditto for the model in the middle, and a diff at the bottom  Caveats  References  Modification history 20150702-A_eval_ma: written.  local aux_title_info,      cn_levels_string,   curr,\ diag_description,    diag_script_base,   dim_MOD,\
   :param integer idx_mod:,                        \ field_type0,         lbtitle,            main_title,\ main_title_prefix,   output_dir,         output_filename,\ output_file_path,    plot,               plot_diff,\ plot_ref,            ref,                res,\ storage_record,      textres,            var0,\ wks,                 wks_debug,          txres begin verbosity  = stringtointeger(getenv("ESMValTool_verbosity")) info_output("<<<<<<<< Entering three_by_one_vector", verbosity, 8) dim_MOD = dimsizes(models@name) dim_VAR = dimsizes(variables) var0 = variables(0) var1 = variables(1) field_type0 = field_types(0) field_type1 = field_types(1)

.. function::  plot_multipanel(cols[*]:integer, rows[*]:float, curr_idx[1]:integer, curr_page[1]:integer, res[1]:logical, storage_name[1]:string, storage_vault[1]:logical, wks[1]:graphic, di[1]:logical, plot_settings[1]:logical, valid_statistics[*]:string, debuginfo[1]:logical, figures_per_page[*]:integer, model_panel_placement[*]:integer, figure_panel_placement[*]:integer, plot_array[*]:graphic, type_specifier[1]:string, no_figures_on_this_page[1]:integer)

   :param integer cols: number of columns for this panel plot
   :param float rows: number of rows for this panel plot
   :param integer curr_idx: current index
   :param integer curr_page: current page (may be more than one)
   :param logical res: valid ncl resources to apply to plot
   :param string storage_name: name of first processed variable in the vault
   :param logical storage_vault: handle/pointer to vault with all models
   :param graphic wks: workstation
   :param logical di: logical with diagnostic script info
   :param logical plot_settings: logical with general plot_settings
   :param string valid_statistics: valid statistics (to be computed if defined)
   :param logical debuginfo: info to put onto plot in debug mode
   :param integer figures_per_page: array with number of figures on each page
   :param integer model_panel_placement: where to place respective model
   :param integer figure_panel_placement: where to place respective figure on the page
   :param graphic plot_array: plot handles/pointers
   :param string type_specifier: kind of plot, 'mean' or 'stddev'
   :param integer no_figures_on_this_page: no of figures on this page

   Return value
  
   Description
      Multipanel plot, plots all models on the current page. Top left entry
      is always the reference model.
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  multipanel(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name [1] : string, debuginfo [1] : logical, valid_statistics [*] : string, res [1] : logical)

   :param  logical storage_vault: handle/pointer to vault with all models
   :param  logical di: logical with diagnostic script info
   :param  logical plot_settings: logical with general plot_settings
   :param  string storage_name: name of first processed variable in the vault
   :param  logical debuginfo: info to put onto plot in debug mode
   :param  string valid_statistics: valid statistics (to be computed if defined)
   :param  logical res: valid ncl resources to apply to plot

   Return value
  
   Description
      Determines how to place a number of contour plots in a grid across
      multiple pages. Loop over pages and call plot_multipanel(...) for
      each page to plot entries.
  
   Caveats
  
   References
  
   Modification history
      20150702-A_eval_ma: written.
  
.. function::  plot_multipanel_vector(cols[*]:integer, rows[*]:float, curr_idx[1]:integer, curr_page[1]:integer, res[1]:logical, storage_names[1]:logical, storage_vault[1]:logical, wks[1]:graphic, di[1]:logical, plot_settings[1]:logical, valid_statistics[*]:string, debuginfo[1]:logical, figures_per_page[*]:integer, model_panel_placement[*]:integer, figure_panel_placement[*]:integer, plot_array[*]:graphic, type_specifier[1]:string, no_figures_on_this_page[1]:integer)

   :param integer cols: number of columns for this panel plot
   :param float rows: number of rows for this panel plot
   :param integer curr_idx: current index
   :param integer curr_page: current page (may be more than one)
   :param logical res: valid ncl resources to apply to plot storage_name: name of first processed variable in the vault storage_vault: handle/pointer to vault with all models wks: workstation di: logical with diagnostic script info plot_settings: logical with general plot_settings valid_statistics: valid statistics (to be computed if defined) debuginfo: info to put onto plot in debug mode figures_per_page: array with number of figures on each page model_panel_placement: where to place respective model figure_panel_placement: where to place respective figure on the page plot_array: plot handles/pointers type_specifier: kind of plot, 'mean' or 'stddev' no_figures_on_this_page: no of figures on this page  Return value  Description Multipanel plot for contour with vector overlay, plots all models on the current page. Top left entry is always the reference model.  Caveats  References  Modification history 20150702-A_eval_ma: written.  local am_infobox_id,       blank_plot,          cn_levels_string,\ curr_figure_pos,     curr_idx,            debugbox,\ diff_model_ref,      diff_model_ua_ref,   diff_model_va_ref,\ dummy_array,         header_prefix,       idx_fig,\ idx_mod,             lbres,               main_title_prefix,\ main_title_suffix,   max_speed_ref,       plot,\ plot_ref,            plot_ref_v,          plot_v,\ plottype_lbres,      res,                 speed,\ speed_ref,           statistics,          storage_record,\ txres,               ua,                  ua_ref,\ va,                  va_ref begin verbosity  = stringtointeger(getenv("ESMValTool_verbosity")) info_output("<<<<<<<< Entering plot_multipanel_vector", verbosity, 8) Update position, labelbar and title curr_figure_pos = figure_panel_placement(curr_idx) res = panel_n_by_cols(res, curr_figure_pos, rows, cols, figures_per_page(curr_page))  if (isatt(res, "cnLevels")) then delete(res@cnLevels) end if
   :param logical storage_names:@base_name + "_basic", "-", "_") res@cnLevels = di@$cn_levels_string$  main_title_prefix = "" main_title_suffix = ""  Fetch reference plot storage_record = str_join((/storage_names@storage_x_component, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep)
   :param logical storage_vault:, storage_record)  storage_record = str_join((/storage_names@storage_y_component, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep) va_ref = retrieve_from_vault(storage_vault, storage_record)  storage_record = str_join((/storage_names@storage_xy_component, sprinti("%i", plot_settings@idx_ref(0))/), str_vault_sep) speed_ref = retrieve_from_vault(storage_vault, storage_record) if (isatt(di, "max_speed_ref")) then max_speed_ref = di@max_speed_ref else max_speed_ref = max(speed_ref) end if  if (isatt(plot_settings, "idx")) then delete(plot_settings@idx) end if plot_settings@idx = plot_settings@idx_ref  dummy_array = (/1e+20/) dummy_array@_FillValue = 1e+20  statistics = True statistics = compute_stat((/"yrs", "mean"/), valid_statistics, speed_ref, dummy_array)
   :param graphic wks:,\ speed_ref,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ statistics,\ res)  delete(statistics) if (storage_names@base_name .ne. "stddev") then plot_ref_v = single_vector(wks,\ max_speed_ref,\ speed_ref,\ ua_ref,\ va_ref,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ res) overlay(plot_ref, plot_ref_v) end if  txres = True txres@txFuncCode = "~" if (debuginfo) then debugbox = write_info(debuginfo) am_infobox_id = place_debuginfo(wks, debugbox, txres, plot_ref) end if  delete(res@cnLevels) cn_levels_string = str_sub_str("cn_levels_" + storage_names@base_name + type_specifier + "_basic", "-", "_")
   :param logical di:@$cn_levels_string$  idx_fig = figure_panel_placement(curr_idx) plot_array(idx_fig) = plot_ref  Skip past the reference plot curr_idx = curr_idx + 1  lbres = True txres = True txres@txFuncCode = "~" 
   :param logical plot_settings:@type .eq. "diff") then inset_top_text(wks, plot_ref, "REF", txres) inset_labelbar(wks, plot_ref, res, "REF", lbres) main_title_suffix = " - REF" else inset_top_text(wks, plot_ref, "Reference", txres) main_title_suffix = "" end if  ------------------------------ Create the non-reference plots ------------------------------ do curr_fig = 1, figures_per_page(curr_page) - 1  main_title_prefix = "" idx_mod = model_panel_placement(curr_idx) idx_fig = figure_panel_placement(curr_idx)  Update placement and labelbar colors res = panel_n_by_cols(res, figure_panel_placement(curr_idx), rows, cols, figures_per_page(curr_page))  storage_record = str_join((/storage_names@storage_x_component, sprinti("%i", idx_mod)/), str_vault_sep) ua = retrieve_from_vault(storage_vault, storage_record)  storage_record = str_join((/storage_names@storage_y_component, sprinti("%i", idx_mod)/), str_vault_sep) va = retrieve_from_vault(storage_vault, storage_record)  storage_record = str_join((/storage_names@storage_xy_component, sprinti("%i", idx_mod)/), str_vault_sep) speed = retrieve_from_vault(storage_vault, storage_record)  statistics = True
   :param string valid_statistics:, speed, dummy_array) if (plot_settings@type .eq. "diff") then Plot mean of differences of first and second data set, first interpolate to the reference data set grid (bilinear interpolation)  statistics = compute_stat((/"rmse"/), valid_statistics, speed, speed_ref)  ua/va/speed field interpolation diff_model_ref = get_model_minus_ref(speed, speed_ref) diff_model_ua_ref = get_model_minus_ref(ua, ua_ref) diff_model_va_ref = get_model_minus_ref(va, va_ref)  delete(speed) speed = diff_model_ref  delete(ua) ua = diff_model_ua_ref  delete(va) va = diff_model_va_ref  if (isatt(di, "max_speed_ref_diff")) then max_speed_ref = di@max_speed_ref_diff else max_speed_ref = max(diff_model_ref) end if  delete(diff_model_ref) delete(diff_model_ua_ref) delete(diff_model_va_ref)  else statistics = compute_stat((/"corr"/), valid_statistics, speed, speed_ref) end if statistics = compute_stat((/"mean"/), valid_statistics, speed, dummy_array)  if (isatt(plot_settings, "idx")) then delete(plot_settings@idx) end if plot_settings@idx = idx_mod   plot = single_contour(wks,\ speed,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\
   :param logical debuginfo:,\ statistics,\ res) delete(statistics) if (storage_names@base_name .ne. "stddev") then plot_v = single_vector(wks,\ max_speed_ref,\ speed,\ ua,\ va,\ main_title_prefix,\ main_title_suffix,\ plot_settings,\ debuginfo,\ res) overlay(plot, plot_v) end if  if (debuginfo) then debugbox = write_info(debuginfo) am_infobox_id = place_debuginfo(wks, debugbox, txres, plot) end if  plot_array(idx_fig) = plot  Update index to point to next field curr_idx = curr_idx + 1 delete(ua) delete(va) delete(speed) 
   :param integer figures_per_page:(curr_page) - 1  plottype_lbres = False  --------------------------------------------------------- Create an blank plot for shared labelbar placement (mean) --------------------------------------------------------- header_prefix = empty_string_if_missing_attribute(plot_settings, "prefix_to_header") blank_plot = add_blank_plot_title(wks,\ header_prefix + di@season + "-" + plot_settings@type + plot_settings@part_of_header,\ rows,\ cols) Create shared labelbar n_by_cols_labelbar(wks,\ blank_plot,\ plot_array(no_figures_on_this_page - 1),\ rows,\ cols,\ plot_settings@lb_units,\ plottype_lbres)  --------------------- Draw mean value plot --------------------- if (debuginfo) then drawNDCGrid(wks) end if draw(plot_array) draw(blank_plot) if (debuginfo) then place_description(wks,\ debuginfo@description,\ debuginfo@description_ycoord) end if info_output(">>>>>>>> Leaving plot_multipanel_vector", verbosity, 8) end  undef("multipanel_vector") procedure multipanel_vector(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_names [1] : logical, debuginfo [1] : logical, valid_statistics [*] : string, res [1] : logical)  Arguments storage_vault: handle/pointer to vault with all models di: logical with diagnostic script info plot_settings: logical with general plot_settings storage_name: name of first processed variable in the vault debuginfo: info to put onto plot in debug mode valid_statistics: valid statistics (to be computed if defined) res: valid ncl resources to apply to plot  Return value  Description Determines how to place a number of contour plots in a grid across multiple pages. Loop over pages and call plot_multipanel(...) for each page to plot entries.  Caveats  References  Modification history 20150702-A_eval_ma: written.  local aux_title_info,          blank_plot,        cn_levels_string,\ cols,                    curr_fig,            curr_figure_pos,\ curr_idx,                curr_page,           diag_script_base,\ dim_MOD,                 dim_VAR,             dummy_array,\ field_type0,             field_type1,         figure_panel_placement,\ figures_per_page,        idx_fig,             idx_mod,\ lbres,                   main_title_prefix,\
   :param integer model_panel_placement:,\ no_figures_on_this_page, output_dir,          output_filename,\ output_file_path,                             page_no,\ plot,                    plot_array,          plot_ref,\ plottype_lbres,          res,                 rows,\ speed,                   speed_ref,           storage_record,\ total_no_of_pages,       txres,               type_specifier,\ ua,                      ua_ref,              va,\ var0,                    var1,                va_ref,\ wks,                     plot_ref_v begin dim_MOD = dimsizes(models@name) dim_VAR = dimsizes(variables) var0 = variables(0) field_type0 = field_types(0) if (dimsizes(variables) .gt. 1) then var1 = variables(1) field_type1 = field_types(1) else var1 = "" field_type1 = "" end if  'output_file_type' if fetched from ncl.interface if (.not. isdefined("output_file_type")) then output_file_type = "ps" end if  Output dir 'plot_dir' if fetched from ncl.interface diag_script_base = basename(plot_settings@diag_script) output_dir = get_output_dir(plot_dir, diag_script_base)  -------------------------------- Static resources for these plots -------------------------------- res@mpFillOn = False res@cnFillOn = True res@cnLinesOn = False res@cnLevelSelectionMode = "ExplicitLevels" res@cnMissingValFillColor = "Background" res@cnLineLabelsOn = False res@gsnFrame = False res@gsnDraw = False res@lbLabelBarOn = False res@gsnAddCyclic = False  -------------------------------------- Compute the layout of paneled figures -------------------------------------- figures_per_page = get_figures_per_page(dim_MOD,\ max_figures_pp,\ min_figures_pp)  Which model goes where across all pages model_panel_placement = new((/sum(figures_per_page)/), integer)  Which model goes where on each page?
   :param integer figure_panel_placement: = new((/sum(figures_per_page)/), integer) place_models_on_pages(models,\ plot_settings@idx_ref,\ figures_per_page,\ model_panel_placement,\ figure_panel_placement)  Output dir 'plot_dir' is fetched from ncl.interface diag_script_base = basename(diag_script) output_dir = get_output_dir(plot_dir, diag_script_base)  if (plot_settings@type .eq. "diff") then type_specifier = "_diff" else type_specifier = "" end if  --------------------------- Loop over all output pages --------------------------- curr_idx = 0 curr_idx_debug = 0 total_no_of_pages = dimsizes(figures_per_page)  do curr_page = 0, total_no_of_pages - 1 --------------------------- Plot arrays for gsn_panels ---------------------------
   :param graphic plot_array: = new((/max_figures_pp/), graphic)  no_figures_on_this_page = figures_per_page(curr_page)  Create a string to add to the figure output\ filename for mulitple pages if (total_no_of_pages .gt. 1) then page_no = "-page" + sprinti("%i", curr_page) else page_no = "" end if  ----------------------------------- Define output workstation for plots ----------------------------------- idx_mod = -1  ; No specific model defined if (isatt(di,"filter_name")) then
   :param string type_specifier: + page_no else aux_title_info = di@season + "-" + storage_names@base_name + type_specifier + page_no end if output_filename = interface_get_figure_filename(diag_script_base,\ var0 + var1,\ field_type0 + field_type1,\ aux_title_info,\ idx_mod) output_file_path = output_dir + output_filename wks = gsn_open_wks(output_file_type, output_file_path) 
   :param integer no_figures_on_this_page:, max_cols) rows = multipanel_get_no_rows(no_figures_on_this_page, max_cols)

.. function::  seasonal_plot(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name [1] : string, debuginfo [1] : logical)

   :param  logical storage_vault: handle/pointer to vault with all models
   :param  logical di: logical with diagnostic script info
   :param  logical plot_settings: logical with general plot_settings
   :param  string storage_name: name of first processed variable in the vault
   :param  logical debuginfo: info to put onto plot in debug mode

   Return value
  
   Description
.. function::  xy_plot_wrapper(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name [1] : string, debuginfo [1] : logical)

   :param  logical storage_vault: handle/pointer to vault with all models
   :param  logical di: logical with diagnostic script info
   :param  logical plot_settings: logical with general plot_settings
   :param  string storage_name: name of first processed variable in the vault
   :param  logical debuginfo: info to put onto plot in debug mode

   Return value
  
   Description
.. function::  ts_line_wrapper(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name [1] : string, debuginfo [1] : logical)

   :param  logical storage_vault: handle/pointer to vault with all models
   :param  logical di: logical with diagnostic script info
   :param  logical plot_settings: logical with general plot_settings
   :param  string storage_name: name of first processed variable in the vault
   :param  logical debuginfo: info to put onto plot in debug mode

   Return value
  
   Description
      Wrapper script for the plot script 'xy_line(..)' with a
      time series on the x-axis.
  
   Caveats
  
   References
  
   Modification history
      20150703-A_eval_ma: written.
  
.. function::  xy_line_overlap(storage_vault [1] : logical, di [1] : logical, plot_settings [1] : logical, storage_name [1] : string, debuginfo [1] : logical)
