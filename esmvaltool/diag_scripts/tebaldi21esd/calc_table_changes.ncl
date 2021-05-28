; #############################################################################
;  Computes the entries for the table of changes given scenario and period
;  Author: Debeire Kevin (DLR, Germany)
; #############################################################################
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
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears", \
                       "begin_ref_year", "end_ref_year"/))
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
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
  n_PER = dimsizes(syears)
  period_list = (/tostring(syears(0)) + "-" + tostring(eyears(0)), \
                  tostring(syears(1)) + "-" + tostring(eyears(1))/)
  scenario_list = diag_script_info@scenarios
  n_SCEN= dimsizes(diag_script_info@scenarios)
  spread = diag_script_info@spread
  mean_table = new((/n_PER, n_SCEN/), float)
  five_perc_table = new((/n_PER, n_SCEN/), float)
  ninetyfive_perc_table = new((/n_PER, n_SCEN/), float)
  number_models_table =  new((/n_PER, n_SCEN/), float)
  do issp = 0, n_SCEN-1
    idx_scen = ind((dataset_exps .eq. diag_script_info@scenarios(issp)))
    n_MOD = dimsizes(idx_scen)
    do iper = 0, n_PER-1
      changes_for_all_models = new((/n_MOD/), float)
      do imod = 0, n_MOD - 1
        ; LOAD HISTORICAL EXP and GET historical baseline
        idx_histo = ind((dataset_exps .eq. reference_run_name) .and. \
                       (dataset_names .eq. dataset_names(idx_scen(imod))) \
                        .and. (dataset_ensembles .eq. \
                        dataset_ensembles(idx_scen(imod))))
        Ahisto = read_data(info_items[idx_histo(0)])
        datahisto = time_operations(Ahisto, -1, -1, "average", "yearly", True)
        delete(Ahisto)
        datahisto_glob = coswgt_areaave(datahisto)
        ind_start = ind(datahisto&year .eq. diag_script_info@begin_ref_year)
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
        ind_start_per = ind(data0_glob&year .eq. syears(iper))
        ind_end_per = ind(data0_glob&year .eq. min((/eyears(iper), \
                          data0_glob&year(dimsizes(data0_glob&year)-1)/)))
        scen_period_val = dim_avg(data0_glob(ind_start_per : ind_end_per))
        if var0 .eq. "pr" then
            changes_current_mod = 100 * (scen_period_val - histo_ref) / \
                                  histo_ref
        else
            changes_current_mod = scen_period_val - histo_ref
        end if
        changes_for_all_models(imod) = changes_current_mod
        delete(idx)
        delete(A0)
        delete(data0)
        delete(data0_glob)
        delete(histo_ref)
        delete(idx_histo)
      end do
      mean_table(iper, issp) = avg(changes_for_all_models)
      five_perc_table(iper, issp) = avg(changes_for_all_models) - \
                                    spread * stddev(changes_for_all_models)
      ninetyfive_perc_table(iper, issp) = avg(changes_for_all_models) + \
                                          spread * \
                                          stddev(changes_for_all_models)
      number_models_table(iper, issp) = n_MOD
      delete(changes_for_all_models)
    end do
    delete(idx_scen)
  end do
  to_plot = new((/4, n_PER, n_SCEN/), float)
  to_plot(0, :, :) = mean_table
  to_plot(1, :, :) = five_perc_table
  to_plot(2, :, :) = ninetyfive_perc_table
  to_plot(3, :, :) = tofloat(number_models_table)
  to_plot!0 = "panel"
  to_plot&panel = (/"mean", "five_perc", "ninetyfive_perc", "num_model"/)
  to_plot!1 = "periods"
  to_plot&periods = period_list
  to_plot!2 = "ssp"
  to_plot&ssp = scenario_list
  to_plot@units     = variable_info[0]@units
  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if (isatt(to_plot, "diag_script")) then  ; Add to existing entries
    tmp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(tmp, (/DIAG_SCRIPT/), 0)
    delete(tmp)
  else  ; Add as new attribute
    to_plot@diag_script = (/DIAG_SCRIPT/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if (isatt(variable_info, "long_name")) then
    to_plot@var_long_name = variable_info@long_name
  end if
  ; Check units and adjust
  if (to_plot@units.eq."K") then
    to_plot@units = "degC"
  end if
  if var0 .eq. "pr" then
    to_plot@units = "%"
  end if
  if (isatt(diag_script_info, "plot_units")) then
    to_plot = convert_units(to_plot, diag_script_info@plot_units)
  end if
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
