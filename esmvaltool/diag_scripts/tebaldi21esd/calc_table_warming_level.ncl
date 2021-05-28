; #############################################################################
;  Produces entries for the table of warming level crossings
;  Author: Debeire Kevin (DLR, Germany)
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"

load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  ; 'datasets', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)
  dataset_names = metadata_att_as_array(info_items, "dataset")
  dataset_exps = metadata_att_as_array(info_items, "exp")
  dataset_ensembles = metadata_att_as_array(info_items, "ensemble")
  n_MOD = count_unique_values(dataset_names)
  MOD_list = get_unique_values(dataset_names)
  labels = new((/n_MOD/), string)
  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears"/))
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  offset = diag_script_info@offset
end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
end

begin
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if
  if (isatt(diag_script_info, "reference_run")) then
    reference_run_name = diag_script_info@reference_run
  else
    log_info("No reference run name given, set to default 'historical'")
    reference_run_name = "historical"
  end if
  attsh = True
  attsh@exp = reference_run_name
  item_hist = select_metadata_by_atts(info_items, attsh)
  dim_MOD_hist = ListCount(item_hist)
  warming_levels = diag_script_info@warming_levels
  scenario_list = diag_script_info@scenarios
  n_SCEN = dimsizes(diag_script_info@scenarios)
  n_WAR = dimsizes(diag_script_info@warming_levels)
  t_len = toint(eyears(1))-toint(syears(1)) + 1
  mean_table = new((/n_WAR, n_SCEN/), integer)
  five_perc_table = new((/n_WAR, n_SCEN/), integer)
  ninetyfive_perc_table = new((/n_WAR, n_SCEN/), integer)
  number_models_crossing_table = new((/n_WAR, n_SCEN/), integer)
  number_models_table =  new((/n_WAR, n_SCEN/), integer)
  do issp = 0, n_SCEN-1
    idx_scen = ind((dataset_exps .eq. diag_script_info@scenarios(issp)))
    n_MOD = dimsizes(idx_scen)
    do iwar = 0, n_WAR-1
      years_for_all_models = new((/n_MOD/), integer)
      all_models_smoothed_ts = new((/n_MOD, t_len/), float)
      do imod = 0, n_MOD - 1
        ; LOAD HISTORICAL EXP and GET historical baseline
        idx_histo = ind((dataset_exps .eq. reference_run_name) .and. \
                        (dataset_names .eq. dataset_names(idx_scen(imod))))
        Ahisto = read_data(info_items[idx_histo(0)])
        datahisto = time_operations(Ahisto, -1, -1, "average", "yearly", True)
        delete(Ahisto)
        datahisto_glob = coswgt_areaave(datahisto)
        ind_start = ind(datahisto&year .eq. \
                        max((/diag_script_info@begin_ref_year, \
                            datahisto&year(0)/)))
        ind_end = ind(datahisto&year .eq. diag_script_info@end_ref_year)
        histo_ref = dim_avg(datahisto_glob(ind_start : ind_end))
        delete([/datahisto, datahisto_glob/])

        A0 = read_data(info_items[idx_scen(imod)])
        dims = getvardims(A0)
        if (dimsizes(dims) .lt. 2) then
          error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + " dimensions, " + \
                    "need 2 or 3")
        end if
        idx = ind(dims .eq. "lat")
        if (ismissing(idx)) then
          error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
        end if
        idx = ind(dims .eq. "lon")
        if (ismissing(idx)) then
          error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
        end if

        ; Average over time (function in ./diag_scripts/shared/statistics.ncl)
        data0 = time_operations(A0, -1, -1, "average", "yearly", True)

        dim_file = dimsizes(data0)  ; file size of scenario file
        ; Average over globe/region
        ; (function in ./diag_scripts/shared/latlon.ncl)
        if (isatt(diag_script_info, "ts_minlat")) then
          data0_glob = area_operations(data0, diag_script_info@ts_minlat, \
                                       diag_script_info@ts_maxlat, \
                                       diag_script_info@ts_minlon, \
                                       diag_script_info@ts_maxlon, "average", \
                                       True)
        else
          data0_glob = coswgt_areaave(data0)
        end if
        running_avg = runave_n(data0_glob, 11, 0, 0) - histo_ref + offset
        t_len_mod = dimsizes(running_avg)
        all_models_smoothed_ts(imod, 0:t_len_mod-1) = running_avg
        idx_warming_level = ind(running_avg .ge. warming_levels(iwar))
        year_warming_level = toint(min(idx_warming_level) + \
                                   diag_script_info@syears(1))
        years_for_all_models(imod) = year_warming_level
        delete(idx_warming_level)
        delete(A0)
        delete(data0)
        delete(data0_glob)
        delete(running_avg)
        delete(histo_ref)
        delete(idx_histo)
      end do
      mean_ts_scenario = dim_avg_n(all_models_smoothed_ts, 0)
      std_ts_scenario = dim_stddev_n(all_models_smoothed_ts, 0)
      five_perc_ts_scenario = mean_ts_scenario - 1.645 * std_ts_scenario
      nin_fiv_perc_ts_scenario = mean_ts_scenario + 1.645 * std_ts_scenario
      year_crossing_mean = \
        toint(min(ind(mean_ts_scenario .ge. warming_levels(iwar)))) + \
        diag_script_info@syears(1)
      year_crossing_five_perc = \
        toint(min(ind(five_perc_ts_scenario .ge. warming_levels(iwar)))) + \
        diag_script_info@syears(1)
      year_crossing_ninety_five_perc = \
        toint(min(ind(nin_fiv_perc_ts_scenario .ge. warming_levels(iwar)))) + \
        diag_script_info@syears(1)
      mean_table(iwar, issp) = year_crossing_mean
      five_perc_table(iwar, issp) = year_crossing_five_perc
      ninetyfive_perc_table(iwar, issp) = year_crossing_ninety_five_perc
      number_models_crossing_table(iwar, issp) = num(.not.ismissing( \
                                                     years_for_all_models))
      number_models_table(iwar, issp) = n_MOD
      delete(all_models_smoothed_ts)
      delete(years_for_all_models)
    end do
    delete(idx_scen)
  end do
  to_plot = new((/5, n_WAR, n_SCEN/), integer)
  to_plot(0, :, :) = mean_table
  to_plot(1, :, :) = five_perc_table
  to_plot(2, :, :) = ninetyfive_perc_table
  to_plot(3, :, :) = number_models_crossing_table
  to_plot(4, :, :) = number_models_table
  to_plot!0 = "panel"
  to_plot&panel = (/"mean", "five_perc", "ninetyfive_perc", \
                    "num_model_crossing", "num_model"/)
  to_plot!1 = "warming"
  to_plot&warming = warming_levels
  to_plot!2 = "ssp"
  to_plot&ssp = scenario_list
  to_plot@units     = "None"

  ; Add to to_plot, as attributes without prefix
  to_plot@diag_script = (/DIAG_SCRIPT/)
  to_plot@var = var0

  ; Check units and adjust for plotting if necessary
  ; difference in K is the same as difference in degC, only change label

  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################
  to_plot@ncdf = variable_info[0]@diagnostic + ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)

  leave_msg(DIAG_SCRIPT, "")
end
