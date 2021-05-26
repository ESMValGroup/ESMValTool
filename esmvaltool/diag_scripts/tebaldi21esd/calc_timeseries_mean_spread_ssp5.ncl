; #############################################################################
;  Computes timeseries of two SSP5 scenarios between 2015-2100
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

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears", \
                                           "begin_ref_year", "end_ref_year"/))
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
  ; ###########################################
  ; # Get data and average annual and globally#
  ; ###########################################
  ; get data from first dataset for historical period
  imod = 0  ; NCL array indices start from zero
  log_debug("processing " + info_items[imod]@dataset + "_" \
            + info_items[imod]@exp + "_" \
            + info_items[imod]@ensemble)

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(info_items[imod])

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
  data1 = time_operations(A0, -1, -1, "average", "yearly", True)

  dim_file = dimsizes(data1)  ; file size of historical file
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

  ; How many historical model runs?
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

  ref_data2 = new((/dim_MOD_hist/), float, data1_glob@_FillValue)
  data2_rel = new((/dim_MOD_hist, dim_file(0)/), float, \
                  data1_glob@_FillValue)

  ; How many scenarios?
  dim_scen = dimsizes(diag_script_info@scenarios)
  ; start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)
  nyears = toint(eyears(nperiods - 1)) - toint(syears(0)) + 1
  nyearsp1 = toint(eyears(nperiods - 2)) - toint(syears(nperiods - 2)) + 1
  nyearsp2 = toint(eyears(nperiods - 1)) - toint(syears(nperiods - 1)) + 1

  atts1 = True
  atts1@exp = diag_script_info@scenarios(0)
  item_ssp534 = select_metadata_by_atts(info_items, atts1)
  dim_MOD_534 = ListCount(item_ssp534)
  ; indicate array for scenarios
  data534_rel = new((/dim_MOD_534, nyearsp1/), float, \
                   data1_glob@_FillValue)
  i534 = 0

  if (dim_scen .gt. 1) then
    atts2 = True
    atts2@exp = diag_script_info@scenarios(1)
    item_ssp585 = select_metadata_by_atts(info_items, atts2)
    dim_MOD_585 = ListCount(item_ssp585)
    data585_rel = new((/dim_MOD_585, nyearsp2/), float, \
                     data1_glob@_FillValue)
    i585 = 0
  end if

  ; set runave_window parameter
  if (isatt(diag_script_info, "runave_window")) then
    runave_window = diag_script_info@runave_window
    half_window= toint(floor(runave_window/2.))
  else
    runave_window = 11
    half_window= toint(floor(runave_window/2.))
  end if

  ; set percent parameter
  if (isatt(diag_script_info, "percent")) then
    percent = diag_script_info@percent
  else
    percent = 0
  end if

  do imod = 0, dim_MOD_hist - 1
    if (imod .eq. 0) then
      data2 = data1
      delete(A0)
    else
      A0 = read_data(info_items[imod])
      data2 = time_operations(A0, -1, -1, "average", "yearly", True)
      delete(A0)
    end if
    ; Average over globe (function in
    ; ./diag_scripts/shared/statistics.ncl)
    data2_glob = coswgt_areaave(data2)
    ind_start = ind(data2&year .eq. diag_script_info@begin_ref_year)
    ind_end = ind(data2&year .eq. diag_script_info@end_ref_year)
    ref_data2(imod) = dim_avg(data2_glob(ind_start : ind_end))
    dim_data2g = dimsizes(data2_glob)
    if (dim_data2g .ne. dim_file(0)) then
      tdiff = dim_file(0) - dim_data2g
      data2_rel(imod, 0 : tdiff - 1) = data2_glob@_FillValue
      data2_rel(imod, tdiff :) =  data2_glob - ref_data2(imod)
    else
      data2_rel(imod, :) = data2_glob - ref_data2(imod)
    end if
    if  percent .eq. 1 then
      data2_rel(imod, :) = 100*(data2_rel(imod, :)/ref_data2(imod))
    end if
    ; clean up
    delete([/data2, data2_glob/])
    ; find all other runs from this model and calculate relative ts
    atts = True
    atts@dataset = info_items[imod]@dataset
    atts@exp = diag_script_info@scenarios(0)
    atts@ensemble = info_items[imod]@ensemble
    item_534 = select_metadata_by_atts(info_items, atts)
    if (ListCount(item_534) .ne. 0) then
      A1 = read_data(item_534[0])
      data3 = time_operations(A1, -1, -1, "average", "yearly", True)
      data3_glob = coswgt_areaave(data3)
      dim_data3g = dimsizes(data3_glob)
      if (dim_data3g .ne. nyearsp1) then
        if (dim_data3g .gt. nyearsp1) then
          error_msg("w", DIAG_SCRIPT, "", "Length of dataset " + \
                    info_items[imod]@dataset + " in recipe is " + \
                    "longer than period given in recipe (" + nyearsp1 + \
                    " years)")
        else
          data534_rel(i534, 0 : dim_data3g - 1) = (data3_glob - ref_data2(imod))
          data534_rel(i534, dim_data3g :) =  data3_glob@_FillValue
        end if
      else
        data534_rel(i534, :) = (data3_glob - ref_data2(imod))
      end if
      if  percent .eq. 1 then
        data534_rel(i534, :) = 100*(data534_rel(i534, :)/ref_data2(imod))
      end if

      i534 = i534 + 1
      ; clean up
      delete([/A1, data3, data3_glob/])
    else
      log_info("No corresponding model for historial run " + \
               info_items[imod]@dataset + " found in scenario " + \
               diag_script_info@scenarios(0) + \
               ", cannot calculate relative change.")
    end if
    delete(item_534)

    if (dim_scen .gt. 1) then
      atts@exp = diag_script_info@scenarios(1)
      item_585 = select_metadata_by_atts(info_items, atts)

      if (ListCount(item_585) .ne. 0) then
        A1 = read_data(item_585[0])
        data3 = time_operations(A1, -1, -1, "average", "yearly", True)
        data3_glob = coswgt_areaave(data3)
        dim_data3g = dimsizes(data3_glob)
        if (dim_data3g .ne. nyearsp2) then
          if (dim_data3g .gt. nyearsp2) then
            error_msg("w", diag_script, "", "Length of dataset " + \
                      info_items[imod]@dataset + " in recipe is " + \
                      "longer than period given in recipe (" + nyearsp2 + \
                      " years)")
          else
            data585_rel(i585, 0:dim_data3g - 1) = (data3_glob - ref_data2(imod))
            data585_rel(i585, dim_data3g:) =  data3_glob@_FillValue
          end if
        else
          data585_rel(i585, :) = (data3_glob - ref_data2(imod))
        end if
        if  percent .eq. 1 then
          data585_rel(i585, :) = 100*(data585_rel(i585, :)/ref_data2(imod))
        end if
        ; perform running average and extend ssp534 with ssp585 before 2040
        app_arr585 = array_append_record(data2_rel(imod,dim_data2g-half_window:dim_data2g-1),data585_rel(i585,:),0)
        tmp_arr585 = runave_n(app_arr585,runave_window,0,0)
        id_junction= nyearsp2-nyearsp1
        app_arr534 = array_append_record(data585_rel(i585,id_junction-half_window:id_junction-1),data534_rel(i534-1,:),0)
        tmp_arr534 = runave_n(app_arr534,runave_window,0,0)
        dim_tmp534= dimsizes(tmp_arr534)
        data534_rel(i534-1,:) = tmp_arr534(half_window:dim_tmp534-1)
        data585_rel(i585, :) = tmp_arr585(half_window:dim_data3g+4)
        delete(tmp_arr534)
        delete(tmp_arr585)
        delete(app_arr534)
        delete(app_arr585)
        i585 = i585 + 1
        ; clean up
        delete([/A1, data3, data3_glob/])
      else
        log_info("No corresponding model ensemble for historial run " + \
                 info_items[imod]@dataset + " found in scenario " + \
                 diag_script_info@scenarios(1) + \
                 ", cannot calculate relative change.")
      end if
      delete(item_585)
    end if

  end do  ; imod

  ; average and standard deviation over all models
  datahist_avg = dim_avg_n_Wrap(data2_rel, 0)
  datahist_std = dim_stddev_n_Wrap(data2_rel, 0)
  ; put all data from 1850-2100 into mean_val and std_val
  mean_val = new((/dim_scen+1, nyearsp2/), typeof(data534_rel), \
                 data534_rel@_FillValue)
  std_val = new((/dim_scen+1, nyearsp2/), typeof(data534_rel), \
                data534_rel@_FillValue)
  p0 = toint(eyears(0)) - toint(syears(0))
  p1 = toint(syears(1)) - toint(syears(2))
  p2 = toint(syears(2)) - toint(syears(2))
  mean_val(0, 0 : p0) = mean_val@_FillValue
  std_val(0, 0 : p0) = std_val@_FillValue

  data534_avg = dim_avg_n_Wrap(data534_rel, 0)
  data534_std = dim_stddev_n_Wrap(data534_rel, 0)
  mean_val(1,p1:) = data534_avg
  std_val(1,p1:) = data534_std
  if (dim_scen .gt. 1) then
    data585_avg = dim_avg_n_Wrap(data585_rel, 0)
    data585_std = dim_stddev_n_Wrap(data585_rel, 0)
    mean_val(2, p2 :) = data585_avg
    std_val(2, p2 :) = data585_std
  end if

  ; check number of models for timeperiods in scenarios
  ; if there are less than 5 models set to missing
  runs = new((/1 + dim_scen, nperiods/), integer)
  atts_runs = True
  atts_runs@exp = reference_run_name
  item_runs = select_metadata_by_atts(info_items, atts_runs)
  runs(0, 0) = ListCount(item_runs)
  do scen = 0, dim_scen - 1
    do p = 1, nperiods - 1
      atts_runs@exp = diag_script_info@scenarios(scen)
      item_runs = select_metadata_by_atts(info_items, atts_runs)
      nrp1 = NewList("fifo")
      nrp2 = NewList("fifo")
      do r = 0, ListCount(item_runs) - 1
        if (item_runs[r]@end_year .gt. toint(diag_script_info@eyears(p))) then
          ListAppend(nrp1, item_runs[r])
        elseif (item_runs[r]@end_year .eq. \
                toint(diag_script_info@eyears(p))) then
          ListAppend(nrp2, item_runs[r])
        end if
      end do
      tmp = ListCount(nrp1)
      ; check if only found missing, in that case set tmp to zero
      if tmp .eq. 0 then
        tmp = 0
      end if
      tmp2 = ListCount(nrp2)
      ; check if only found missing, in that case set tmp2 to zero
      if tmp2 .eq. 0 then
        tmp2 = 0
      end if
      runs(scen + 1, p) = tmp2 + tmp
      delete([/tmp, tmp2, nrp1, nrp2/])
      if (runs(scen + 1, p) .le. 4) .and. (p .ne. 0) then
        p2 = toint(eyears(p - 1)) - toint(syears(0))
        p3 = toint(eyears(p)) - toint(syears(0))
        mean_val(scen + 1, p2 : p3) = mean_val@_FillValue
        std_val(scen + 1, p2 : p3) = std_val@_FillValue
        log_info("Scenario " + diag_script_info@scenarios(scen) + \
                 " in period ending " + diag_script_info@eyears(p) + \
                 " has less than 5 models, set to missing.")
      end if
    end do
  end do
  ;FIX MODEL_NR
  dim_1=dimsizes(data534_rel)
  dim_2=dimsizes(data585_rel)
  runs(1,1) = dim_1(0)
  runs(2,1) = dim_2(0)

  dims_data = dimsizes(mean_val)
  to_plot = new((/3 * dims_data(0), dims_data(1)/), float)
  if (isatt(diag_script_info, "spread")) then
    spread = diag_script_info@spread
  else
    spread = 1.0
  end if
  do j = 0, dims_data(0) - 1
    ;  if variable sw or lw TOA flux: change direction
    if (var0 .eq. "rlut" .or. var0 .eq. "rsut") then
      to_plot(j * 3, :) = - mean_val(j, :)
    ; if variable is pr change unit to mm/day
    ;else if var0 .eq. "pr" then
    ;  to_plot(j * 3, :) = mean_val(j, :)*86400
    else
      to_plot(j * 3, :) = mean_val(j, :)
    end if
    ;end if
    to_plot((j * 3) + 1, :) = to_plot(j * 3, :) + spread * std_val(j, :)
    to_plot((j * 3) + 2, :) = to_plot(j * 3, :) - spread * std_val(j, :)
  end do
  copy_VarMeta(mean_val, to_plot)
  to_plot!0 = "metric"
  str_spread = sprintf("%5.2f", spread)
  to_plot&metric = (/"mean", "+" + str_spread + "std", "-" + \
                     str_spread + "std", "mean", "+" + \
                     str_spread + "std", "-" + str_spread + \
                     "std", "mean", "+" + str_spread + "std","-" + str_spread + "std"/)
  to_plot!1 = "year"
  to_plot&year = ispan(toint(syears(2)), toint(eyears(nperiods - 1)), 1)

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
  if  percent .eq. 1 then
    to_plot@units     = "%"
  end if

  runs!0 = "scenario"
  runs&scenario = diag_script_info@label
  runs!1 = "period"
  runs&period = diag_script_info@syears
  if(isatt(runs, "diag_script")) then  ; Add to existing entries
    temp = runs@diag_script
    delete(runs@diag_script)
    runs@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    runs@diag_script = (/DIAG_SCRIPT/)
  end if
  runs@var = "nr"
  runs@var_long_name = "number of model runs per scenario and period"
  runs@units     = 1
  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################
  to_plot@ncdf = variable_info[0]@diagnostic + ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  if(isatt(diag_script_info, "model_nr") .and. \
     diag_script_info@model_nr .eq. True) then
    ncdf_file2 = work_dir + "/nr_runs_" + variable_info[0]@diagnostic + ".nc"
    ncdf_outfile2 = ncdf_write(runs, ncdf_file2)
  end if

  leave_msg(DIAG_SCRIPT, "")
end
