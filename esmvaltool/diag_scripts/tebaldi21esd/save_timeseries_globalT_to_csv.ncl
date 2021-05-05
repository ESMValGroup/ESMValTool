; #############################################################################
;  Computes changes between two periods and store results in a CSV file
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
  MOD_list= get_unique_values(dataset_names)
  print(MOD_list)
  labels = new((/n_MOD/),string)
  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears"/))
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
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
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
  scenario_list = diag_script_info@scenarios
  n_SCEN= dimsizes(diag_script_info@scenarios)
  DATA = True
  list_name = NewList("lifo")
  list_timeseries = NewList("lifo")
  ListAppend(list_name,"YEAR")
  years = ispan(syears(1),eyears(1),1)
  n_years= dimsizes(years)
  ListAppend(list_timeseries,years)
  all_format = "%4i"
  header_format = "%s"
  do i_scen= 0, n_SCEN-1
    current_scen = scenario_list(i_scen)
    idx_scen = ind((dataset_exps .eq. current_scen))
    print(idx_scen)
    print(current_scen) 
    n_MOD = dimsizes(idx_scen)
    do imod=0, n_MOD-1
      ; LOAD HISTORICAL EXP and GET historical baseline
      ;LOAD Histo file
      print(dataset_names(idx_scen(imod)))
      idx_histo= ind((dataset_exps .eq. reference_run_name) .and. (dataset_names .eq. dataset_names(idx_scen(imod))) .and. (dataset_ensembles .eq. dataset_ensembles(idx_scen(imod))) )
      Ahisto = read_data(info_items[idx_histo(0)])
      datahisto = time_operations(Ahisto, -1, -1, "average", "yearly", True)
      delete(Ahisto)
      datahisto_glob = coswgt_areaave(datahisto)
      ind_start = ind(datahisto&year .eq. max((/diag_script_info@begin_ref_year,datahisto&year(0)/)))
      ind_end = ind(datahisto&year .eq. diag_script_info@end_ref_year)
      histo_ref = dim_avg(datahisto_glob(ind_start : ind_end))
      delete([/datahisto, datahisto_glob/])

      ; LOAD scenario file
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
      if current_scen .eq. "ssp534-over" then
        data0rel = new((/n_years/),float)
        ssp585_name= "ssp585"+"_"+dataset_names(idx_scen(imod))+"_"+dataset_ensembles(idx_scen(imod))
        data_ssp585= DATA@$ssp585_name$
        data0rel(0:24)=data_ssp585(0:24)
        data0rel(25:)= data0_glob(0:60) - histo_ref
      else
        data0rel = new((/n_years/),float)
        data0rel(0:dimsizes(data0_glob)-1)= data0_glob - histo_ref
      end if
      ;print(data0rel) 
      delete(A0)
      delete(data0)
      delete(data0_glob)
      delete(histo_ref)
      delete(idx_histo)
      tmp_header = header_format
      tmp_format = all_format
      current_name = current_scen+"_"+dataset_names(idx_scen(imod))+"_"+dataset_ensembles(idx_scen(imod))
      ListAppend(list_name, current_scen+"_"+dataset_names(idx_scen(imod))+"_"+dataset_ensembles(idx_scen(imod)))
      DATA@$current_name$=data0rel  
      ListAppend(list_timeseries,DATA@$current_name$)
      all_format = tmp_format +",%3.3f"
      header_format = tmp_header +",%s"
      delete(tmp_format)
      delete(tmp_header)
      delete(data0rel)
    end do
    delete(idx_scen)
  end do
  write_table(work_dir+"/all_timeseries.csv", "w",list_name, header_format)
  write_table(work_dir+"/all_timeseries.csv", "a",list_timeseries, all_format)
  delete(DATA)
end
