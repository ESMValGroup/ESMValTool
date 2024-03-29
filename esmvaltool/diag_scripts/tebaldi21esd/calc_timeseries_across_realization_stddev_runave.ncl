; #############################################################################
;  Compute across realization Std. dev. for large ensemble experiment(SSP3-7.0)
;  Author: Kevin Debeire (DLR, Germany)
; #############################################################################
;
;  Description
;     Computes time series of ensemble spreads (inter-member stddev).
;     One dataset is used for resampling subsets of 10 members.
;
;  Required diag_script_info attributes (diagnostics specific)
;     scenarios: list with scenarios to be included in the
;                figure, e.g (/"ssp370"/)
;     syears: list with start years in time periods (e.g. start of historical
;             period and SSPs)
;     eyears: list with end years in time periods (end year of historical
;             period and SSPs)
;     begin_ref_year: start year of reference period (e.g. 1995)
;     end_ref_year: end year of reference period (e.g. 2014)
;     sampled_model: name of dataset on which to sample (ex. CanESM5)
;     n_samples: number of samples of size 10 to draw among all the
;                ensembles of sampled_model
;     runave_window: size window used for the centered running average
;
;  Modification history
;     20220803-A_debe_kevin: written for ESMValTool
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
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
  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears"/))
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
end

begin
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

  t_len = toint(eyears(1))-toint(syears(0))+1

  if (isatt(diag_script_info, "n_samples")) .and. \
     (isatt(diag_script_info, "sampled_model")) then
    n_samples = toint(diag_script_info@n_samples)
    to_plot = new((/n_MOD+n_samples, t_len/), float)
    to_plot!0 = "metric"
    to_plot&metric = new(n_MOD+n_samples, string, "std")
  else
    to_plot = new((/n_MOD, t_len/), float)
    to_plot!0 = "metric"
    to_plot&metric = new(n_MOD, string, "std")
  end if
  ; Get runave_window
  if (isatt(diag_script_info, "runave_window")) then
    runave_window = diag_script_info@runave_window
  else
    runave_window = 1
  end if
  to_plot!1 = "year"
  to_plot&year = ispan(toint(syears(0)), toint(eyears(1)), 1)

  do imod = 0, n_MOD - 1
    current_mod = MOD_list(imod)
    log_debug("processing " + info_items[imod]@dataset)
    ; find ssp runs from same model
    idx_scen = ind((dataset_names .eq. current_mod) .and. \
                   (dataset_exps .eq. diag_script_info@scenarios(0)))
    dim_ENS =  dimsizes(idx_scen)
    n_ENS = 10
    collect_n_ENS = new((/n_ENS, t_len/), double)
    ; find ssp runs from same model
    do iens = 0, n_ENS - 1   ; loop over ensembles
      ; LOAD AND AVERAGE SSP data
      A0 = read_data(info_items[idx_scen(iens)])
      ; Check dimensions
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

      ; LOAD AND AVERAGE HISTORICAL data
      idx_hist = ind((dataset_names .eq. current_mod) .and. \
                     (dataset_exps .eq. reference_run_name) .and. \
                     (dataset_ensembles .eq. \
                      dataset_ensembles(idx_scen(iens))))
      A1 = read_data(info_items[idx_hist])
      dims1 = getvardims(A1)
      if (dimsizes(dims1) .lt. 2) then
        error_msg("f", DIAG_SCRIPT, "", dimsizes(dims1) + " dimensions, " + \
                  "need 2 or 3")
      end if
      idx1 = ind(dims1 .eq. "lat")
      if (ismissing(idx1)) then
        error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
      end if
      idx = ind(dims1 .eq. "lon")
      if (ismissing(idx1)) then
        error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
      end if
      ; Average over time (function in ./diag_scripts/shared/statistics.ncl)
      data1 = time_operations(A1, -1, -1, "average", "yearly", True)
      dim_file1 = dimsizes(data1)  ; file size of scenario file
      ; Average over globe/region
      ; (function in ./diag_scripts/shared/latlon.ncl)
      if (isatt(diag_script_info, "ts_minlat")) then
        data1_glob = area_operations(data1, diag_script_info@ts_minlat, \
                                     diag_script_info@ts_maxlat, \
                                     diag_script_info@ts_minlon, \
                                     diag_script_info@ts_maxlon, "average", \
                                     True)
      else
        data1_glob = coswgt_areaave(data1)
      end if
      ensemble = dataset_ensembles(idx_scen(iens))
      model = current_mod
      data_glob = array_append_record(data1_glob, data0_glob, 0)
      asciiwrite(work_dir + "/" + var0 + "_timeseries_" + model + "_" + \
                 ensemble + ".txt", data_glob)
      ; store for all ens
      if var0 .eq. "pr" then
        collect_n_ENS(iens, :) = runave_n(data_glob, runave_window, 0, 0) \
                                 * 86400.0
      else
        collect_n_ENS(iens, :) = runave_n(data_glob, runave_window, 0, 0)
      end if
      delete(data_glob)
      delete(A0)
      delete(A1)
      delete(data0)
      delete(data0_glob)
      delete(data1)
      delete(data1_glob)
    end do
    to_plot(imod, :) = tofloat(dim_stddev_n(collect_n_ENS, 0))
    delete(idx_scen)
    delete(collect_n_ENS)
  end do
  ; SAMPLING across ensemble of SAMPLED MODEL
  if (isatt(diag_script_info, "sampled_model") .and. \
      isatt(diag_script_info, "n_samples")) then
    sampled_model = diag_script_info@sampled_model
    n_samples = diag_script_info@n_samples
    idx_samp_mod = ind((dataset_names .eq. sampled_model) .and. \
                       (dataset_exps .eq. diag_script_info@scenarios(0)))

    dim_ENS =  dimsizes(idx_samp_mod)
    n_ENS = 10
    random_setallseed(1234567890, 236484749)
    do i_samp = 0, n_samples - 1
      print("Drawing sample #" + i_samp)
      collect_n_ENS = new((/n_ENS, t_len/), double)
      ; randomly sample across ensemble for sampled_model
      rand_idx_50 = generate_sample_indices(dim_ENS, 0)
      rand_idx = idx_samp_mod(rand_idx_50(0:9))
      print("current sample indices " + rand_idx)
      do iens = 0, n_ENS - 1
        A0 = read_data(info_items[rand_idx(iens)])
        ; Check dimensions
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

        ; LOAD AND AVERAGE HISTORICAL data
        idx_hist = ind((dataset_names .eq. sampled_model) .and. \
                       (dataset_exps .eq. reference_run_name) .and. \
                       (dataset_ensembles .eq. \
                        dataset_ensembles(rand_idx(iens))))
        A1 = read_data(info_items[idx_hist])
        dims1 = getvardims(A1)
        if (dimsizes(dims1) .lt. 2) then
          error_msg("f", DIAG_SCRIPT, "", dimsizes(dims1) + " dimensions, " + \
                    "need 2 or 3")
        end if
        idx1 = ind(dims1 .eq. "lat")
        if (ismissing(idx1)) then
          error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
        end if
        idx = ind(dims1 .eq. "lon")
        if (ismissing(idx1)) then
          error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
        end if
        ; Average over time (function in ./diag_scripts/shared/statistics.ncl)
        data1 = time_operations(A1, -1, -1, "average", "yearly", True)
        dim_file1 = dimsizes(data1)
        ; Average over globe/region
        ; (function in ./diag_scripts/shared/latlon.ncl)
        if (isatt(diag_script_info, "ts_minlat")) then
          data1_glob = area_operations(data1, diag_script_info@ts_minlat, \
                                       diag_script_info@ts_maxlat, \
                                       diag_script_info@ts_minlon, \
                                       diag_script_info@ts_maxlon, "average", \
                                       True)
        else
          data1_glob = coswgt_areaave(data1)
        end if
        ensemble = dataset_ensembles(rand_idx(iens))
        model = current_mod
        data_glob = array_append_record(data1_glob, data0_glob, 0)

        ; store for all ens
        if var0 .eq. "pr" then
          collect_n_ENS(iens, :) = runave_n(data_glob, runave_window, 0, 0) \
                                   * 86400.0
        else
          collect_n_ENS(iens, :) = runave_n(data_glob, runave_window, 0, 0)
        end if
        delete(data_glob)
        delete(A0)
        delete(A1)
        delete(data0)
        delete(data0_glob)
        delete(data1)
        delete(data1_glob)
      end do
      to_plot(n_MOD+i_samp, :) = tofloat(dim_stddev_n(collect_n_ENS, 0))
      delete(collect_n_ENS)
    end do
  end if

  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if(isatt(to_plot, "diag_script")) then  ; Add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/DIAG_SCRIPT/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if(isatt(variable_info[0], "long_name")) then
    to_plot@var_long_name = variable_info[0]@long_name
  end if
  if(isatt(variable_info[0], "units")) then
    to_plot@units     = variable_info[0]@units
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
