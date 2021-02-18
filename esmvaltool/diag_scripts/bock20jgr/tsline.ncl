; #############################################################################
; PLOTS TIME SERIES
; Authors: Lisa Bock (DLR, Germany)
; #############################################################################
;
; Description
;    This script is for plotting a time series of the global mean (anomaly).
;
; Required diag_script_info attributes (diagnostic specific)
;    styleset: as in diag_scripts/shared/plot/style.ncl functions
;
; Optional diag_script_info attributes (diagnostic specific)
;    time_avg: type of time average (currently only "yearly" and "monthly" are
;              available).
;    ts_anomaly: calculates anomalies with respect to the defined period;
;                for each gird point by removing the mean for the given
;                calendar month (requiring at least 50% of the data to be
;                non-missing)
;    ref_start: start year of reference period for anomalies
;    ref_end: end year of reference period for anomalies
;    ref_value: if true, right panel with mean values is attached
;    ref_mask: if true, model fields will be masked by reference fields
;    region: name of domain
;    plot_units: variable unit for plotting
;    y_min: set min of y-axis
;    y_max: set max of y-axis
;    mean_nh_sh: if true, calculate first NH and SH mean
;    volcanoes: if true, lines of main volcanic eruptions will be added
;    header: if true, region name as header
;    write_stat: if true, write multi model statistics in nc-file
;
; Caveats
;
; Modification history
;    20190911-bock_lisa: included method of Jones et al., 2013 and
;                        added provenance
;    20181112-bock_lisa: code rewritten for ESMValTool v2.0
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  project0 = input_file_info[0]@project
  if (project0 .eq. "OBS") then
    project0 = input_file_info[1]@project
    exp0 = input_file_info[1]@exp
    obs_first = True
  else
    exp0 = input_file_info[0]@exp
    obs_first = False
  end if
  info_items = select_metadata_by_name(input_file_info, var0)
  datasetnames = metadata_att_as_array(info_items, "dataset")
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Time averaging
  if (isatt(diag_script_info, "time_avg")) then
    time_avg = diag_script_info@time_avg
    if (all(time_avg.ne.(/"yearly", "monthly"/))) then
        error_msg("f", DIAG_SCRIPT, "", \
                  "time averaging option " + time_avg + " not yet " + \
                  "implemented")
    end if
  else
    time_avg = "monthly"
  end if

  ; Check for required settings
  exit_if_missing_atts(diag_script_info, "styleset")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "time_avg", "monthly")
  set_default_att(diag_script_info, "ts_anomaly", "noanom")
  set_default_att(diag_script_info, "ref_value", True)
  set_default_att(diag_script_info, "ref_mask", False)
  set_default_att(diag_script_info, "region", "Global")
  set_default_att(diag_script_info, "mean_nh_sh", False)
  set_default_att(diag_script_info, "header", False)
  set_default_att(diag_script_info, "volcanoes", False)
  set_default_att(diag_script_info, "write_stat", False)
  set_default_att(diag_script_info, "stat", "MinMax")

  ; Determine time range
  start_year = min(metadata_att_as_array(info_items, "start_year"))
  end_year = max(metadata_att_as_array(info_items, "end_year"))
  all_years = ispan(start_year, end_year, 1)

  ; Create time coordinate
  if (time_avg.eq."monthly") then
    ntime = 12 * (end_year - start_year + 1)
    time = new(ntime, integer)
    do yy = start_year, end_year
        do mm = 1, 12
            time(12 * (yy - start_year) + mm - 1) = 100 * yy + mm
        end do
    end do
  elseif (time_avg.eq."yearly") then
    ntime = end_year - start_year + 1
    time = new(ntime, integer)
    time = ispan(start_year, end_year, 1)
  end if
  ntime_old = ntime

  ; get multi-model mean index
  mmm_ind = ind(datasetnames .eq. "MultiModelMean")

  ; Set index of the reference dataset and read it
  if (isatt(variable_info[0], "reference_dataset")) then
    ref_ind = ind(datasetnames.eq.variable_info[0]@reference_dataset)
    A_ref = read_data(info_items[ref_ind])
  end if

  ; Anomaly
  if (isatt(diag_script_info, "ts_anomaly")) then
    anom = diag_script_info@ts_anomaly
    if (anom .eq. "anom") then
      if (isatt(diag_script_info, "ref_start") .and. \
          isatt(diag_script_info, "ref_end")) then
        ref_start = diag_script_info@ref_start
        ref_end   = diag_script_info@ref_end
        if (ref_start.lt.start_year) then
          ref_start = start_year
        end if
        if ((ref_start.lt.start_year) .or. (ref_end.gt.end_year) .or. \
            (ref_end.lt.ref_start)) then
          error_msg("f", DIAG_SCRIPT, "", \
                    "period for reference years is not properly defined")
        end if
      else
        error_msg("f", DIAG_SCRIPT, "", \
                  "period for reference years is not defined " + \
                  "(needed for anomaly)")
      end if

      anom_ref = new((/dim_MOD/), double)
      anom_ref!0 = "model"
      anom_ref&model = datasetnames
    end if
  else
    anom = "noanom"
  end if

  ; Create model array
  model_arr = new((/dim_MOD, ntime/), double)
  model_arr!0 = "model"
  model_arr!1 = "time"
  model_arr&model = datasetnames
  model_arr&time = time
  model_arr@_FillValue = 1e+20
  model_arr = model_arr@_FillValue

  ; Loop over models
  do imod = 0, dim_MOD - 1

    log_info("Process dataset: " + datasetnames(imod))

    ; Read data
    A0 = read_data(info_items[imod])
    dnames = getVarDimNames(A0)

    ; Convert units for plotting (if required)
    if (isatt(diag_script_info, "plot_units")) then
      A0 = convert_units(A0, diag_script_info@plot_units)
    end if

    if (isatt(variable_info[0], "reference_dataset")) then

      ; masking with reference dataset
      if (diag_script_info@ref_mask .and. imod .ne. ref_ind) then
        ; Determine start/end year
        start_year = info_items[imod]@start_year
        end_year = info_items[imod]@end_year
        if (start_year .gt. diag_script_info@ref_start) then
          ref_start = start_year
        else
          ref_start = diag_script_info@ref_start
        end if
        A_ref_mask = time_operations(A_ref, start_year, \
                                     end_year, "extract", "", 0)
        A0_nomask_ref = time_operations(A0, ref_start, \
                                        ref_end, "extract", "", 0)
        A0 = where(A_ref_mask.eq.A_ref_mask@_FillValue, \
                   A_ref_mask@_FillValue, A0)
        delete(A_ref_mask)
      end if

    end if

    ; Anomaly
    if (anom .eq. "anom") then
      ; calculate monthly mean of ref period if 0.5 data points are available
      tmp = time_operations(A0, ref_start, ref_end, "extract", "", 0)
      limit  = toint(0.5 * dimsizes(tmp&time))
      do i = 0, dimsizes(tmp!0)-1
        tmp(i, :, :) = where(dim_num_n(.not.ismissing(tmp), 0).ge.limit, \
                             tmp(i, :, :), tmp@_FillValue)
      end do
      A0_monavg_ref = time_operations(tmp, ref_start, ref_end, "average", \
                                      "monthlyclim", True)

      ; calculate anomaly for each grid point
      do i = 0, dimsizes(A0&time) - 1
        A0(i, :, :) = A0(i, :, :) - A0_monavg_ref(mod(i, 12), :, :)
      end do

      if (diag_script_info@ref_value) then
        ; Calculate time average of ref period
        ; annual mean if at least 2 months of data is available
        ; if masking then take original unmasked dataset for ref mean
        if (diag_script_info@ref_mask .and. imod .ne. ref_ind) then
          tmp = A0_nomask_ref
          delete(A0_nomask_ref)
        end if
        date := cd_calendar(tmp&time, 0)
        year := date(:, 0)
        month := date(:, 1)
        weights = days_in_month(toint(year), toint(month))
        A0_timavg_ref = dim_avg_wgt_n_Wrap(tmp, weights, 2, 0)
        delete(year)
        delete(month)
        delete(date)
        delete(weights)
        delete(A0_monavg_ref)

        ; calculate global mean of reference period
        ; first for each hemisphere - if choosen
        if (diag_script_info@mean_nh_sh) then
          tmp1 = area_operations(A0_timavg_ref, -90., 0., 0., 360., \
                                 "average", True)
          tmp2 = area_operations(A0_timavg_ref, 0., 90., 0., 360., \
                                 "average", True)
          anom_ref_tmp = (tmp1 + tmp2) / 2.
          delete(tmp1)
          delete(tmp2)
        else
          anom_ref_tmp = area_operations(A0_timavg_ref, -90., 90., \
                                         0., 360., "average", True)
        end if
        delete(A0_timavg_ref)

        anom_ref(imod) = anom_ref_tmp
        ; delete(anom_ref_tmp)
      end if

      delete(tmp)

    end if

    ; Calculate time average of dataset
    ; annual mean if at least 2 months of data is available
    if (time_avg.eq."yearly") then
      do i = 0, dimsizes(A0!0) - 1
        A0(i, :, :) = where(dim_num_n(.not.ismissing(A0), 0).ge.2, \
                            A0(i, :, :), A0@_FillValue)
      end do
      A0_timavg = time_operations(A0, -1, -1, "average", \
                                  "yearly", True)
    else
      A0_timavg = A0
    end if
    delete(A0)

    ; calculate global mean of complete time series
    ; first for each hemisphere - if choosen
    if (diag_script_info@mean_nh_sh) then
      tmp1 = area_operations(A0_timavg, -90., 0., 0., 360., \
                             "average", True)
      tmp2 = area_operations(A0_timavg, 0., 90., 0., 360., \
                             "average", True)
      procmod = tmp1
      procmod = (tmp1 + tmp2) / 2.
      delete(tmp1)
      delete(tmp2)
    else
      procmod = area_operations(A0_timavg, -90., 90., 0., 360., \
                                "average", True)
    end if
    delete(A0_timavg)

    if (.not.isdefined("procmod")) then
      error_msg("f", diag_script_info, "", "cannot process this " + \
                "field type " + field_type0)
    end if

    ; Match time coordinate
    if (time_avg.eq."monthly") then
      date = cd_calendar(procmod&time, -1)
    elseif (time_avg.eq."yearly") then
      date = procmod&year
    end if
    idx1 = ind(date(0).eq.model_arr&time)
    idx2 = ind(date(dimsizes(date) - 1).eq.model_arr&time)
    model_arr(imod, idx1:idx2) = (/procmod/)
    if (imod.eq.0) then
      copy_VarAtts(procmod, model_arr)
    end if
    delete(procmod)
    delete(date)

  end do

  ; rm OBS from last entry
  if (obs_first) then
    tmp = model_arr(0, :)
    tmp1 =  info_items[0]
    do i = 0, dim_MOD - 3
      model_arr(i, :) = model_arr(i + 1, :)
      model_arr&model(i) = model_arr&model(i + 1)
      info_items[i] = info_items[i + 1]
    end do
    model_arr(dim_MOD - 2, :) = tmp
    model_arr&model(dim_MOD - 2) = variable_info[0]@reference_dataset
    info_items[dim_MOD - 2] = tmp1
    ref_ind = dim_MOD - 2
    delete(tmp)
    delete(tmp1)
  end if

  ; Convert time coordinate to years (required by the plot routine)
  if (time_avg.eq."monthly") then
    year = model_arr&time / 100
    xmin = toint(floor(min(year)))
    xmax = toint(ceil(max(year)))
    month = model_arr&time - 100 * year
    time = todouble(year + month/12.)
    delete(model_arr&time)
    model_arr&time = time
    delete(time)
    delete(year)
    delete(month)
  elseif (time_avg.eq."yearly") then
    xmin = toint(floor(min(model_arr&time)))
    xmax = toint(ceil(max(model_arr&time)))
    tmp = todouble(model_arr&time)
    delete(model_arr&time)
    model_arr&time = tmp
    delete(tmp)
  end if

  ; individual case for HadCRUT4 observations
  ; mean value comes from climatology file (absolute.nc).
  ; There are no missing values as in the anomaly data.
  if (anom .eq. "anom" .and. isatt(variable_info[0], "reference_dataset") \
      .and. variable_info[0]@reference_dataset .eq. "HadCRUT4") then
    if(ref_start .eq. 1961 .and. ref_end .eq. 1990) then
      anom_ref(ref_ind) = 14.0
    else
      anom_ref(ref_ind) = anom_ref@_FillValue
    end if
  end if

  ; calculate multi model statistics
  if (diag_script_info@write_stat) then
    if (isatt(variable_info[0], "reference_dataset")) then
      ind_mod = ind(datasetnames.ne."MultiModelMean" .and. \
                    datasetnames.ne.variable_info[0]@reference_dataset)
    else
      ind_mod = ind(datasetnames.ne."MultiModelMean")
    end if

    model_arr_stat = new((/ntime, 10/), double)
    model_arr_stat!1 = "stat"
    model_arr_stat!0 = "time"
    model_arr_stat&time = model_arr&time
    model_arr_stat&stat = (/"mean", "stddev", "spread", "min", "5%", "10%", \
                            "90%", "95%", "max", "ref"/)

    model_arr_stat(:, 0) = (/model_arr(mmm_ind, :)/)
    do i = 0, ntime-1
      statb = stat_dispersion(model_arr(ind_mod, i), False)
      ; standard deviation
      model_arr_stat(i, 1) = (/statb(1)/)
      ; spread
      model_arr_stat(i, 2) = (/statb(14)-statb(2)/)
      ; min
      model_arr_stat(i, 3) = (/statb(2)/)
      ; max
      model_arr_stat(i, 8) = (/statb(14)/)
      ; 5% quantile
      model_arr_stat(i, 4) = (/statb(24)/)
      ; 95% quantile
      model_arr_stat(i, 7) = (/statb(25)/)
      ; 10% quantile
      model_arr_stat(i, 5) = (/statb(3)/)
      ; 90% quantile
      model_arr_stat(i, 6) = (/statb(13)/)
    end do
    if (isatt(variable_info[0], "reference_dataset")) then
      model_arr_stat(:, 9) = (/model_arr(ref_ind, :)/)
    end if

  end if

  ; *************************************
  ; output to NetCDF
  ; *************************************
  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)
  out_path1 = out_path + "tsline_" + var0 + "_" + anom + ".nc"
  model_arr@ncdf = out_path
  model_arr@experiment = project0 + "_" + exp0
  model_arr@diag_script = DIAG_SCRIPT
  model_arr@var = var0
  ncdf_outfile = ncdf_write(model_arr, out_path1)
  ; statistics
  if (diag_script_info@write_stat) then
    out_path2 = out_path + "tsline_mm_" + var0 + "_" + anom + "_stat.nc"
    model_arr_stat@ncdf = out_path
    model_arr_stat@experiment = project0 + "_" + exp0
    model_arr_stat@diag_script = DIAG_SCRIPT
    model_arr_stat@var = var0
    model_arr_stat@units = model_arr@units
    ncdf_outfile = ncdf_write(model_arr_stat, out_path2)
  end if

  ; Define workstation
  outfile = config_user_info@plot_dir + var0 + "_" \
    + str_sub_str(diag_script_info@region, " ", "_") \
    + "_" + project0 + "_" + exp0 + "_" + anom + "_" \
    + start_year + "-" + end_year
  wks = gsn_open_wks(file_type, outfile)

  ; Set resources
  res = True
  res@trXMinF = xmin
  res@trXMaxF = xmax
  ; res@trXMaxF = xmax + 0.25 * (xmax - xmin)
  res@tmXBMode = "Explicit"
  if (xmax - xmin.gt.20) then
    res@tmXBValues = ispan(xmin, xmax, 10)
    res@tmXBLabels = ispan(xmin, xmax, 10)
    res@tmXBMinorValues =  ispan(xmin, xmax, 5)
  else
    res@tmXBValues = ispan(xmin, xmax, 5)
    res@tmXBLabels = ispan(xmin, xmax, 5)
    res@tmXBMinorValues =  ispan(xmin, xmax, 1)
  end if
  res@tmXBLabelAngleF = 45
  res@tmXBLabelJust     = "CenterRight"
  if (isatt(diag_script_info, "y_min")) then
    res@trYMinF = diag_script_info@y_min
  end if
  if (isatt(diag_script_info, "y_max")) then
    res@trYMaxF = diag_script_info@y_max
  end if

  res@tmXBMode         = "Manual"
  res@tmXBTickSpacingF = 20

  if (diag_script_info@header) then
    res@tiMainString = diag_script_info@region
  else
    res@tiMainString = ""
  end if

  if (isatt(variable_info[0], "long_name")) then
    if (var0 .eq. "tas") then
      varname = "Temperature"
    else
      varname = variable_info[0]@long_name
    end if
  else
    varname = var0
  end if

  if (model_arr@units .eq. "degC") then
    units = "~F34~0~F~ C"
  else
    units = model_arr@units
  end if

  if (anom .eq. "anom") then
    res@tiYAxisString = varname + " Anomaly" + " (" + units + ")"
  else
    res@tiYAxisString = varname + " (" + units + ")"
  end if

  res0 = True
  ref_start = diag_script_info@ref_start
  ref_end = diag_script_info@ref_end
  res0@tiYAxisString = tostring(ref_start) + "-" + tostring(ref_end) \
    + " Mean " + varname + " (" + units + ")"

  if (anom .eq. "anom") then
    xy_line_anom(wks, anom_ref, model_arr, model_arr&time, \
                 0, ref_start, ref_end, res, res0, info_items)
  else
    xy_line(wks, model_arr, model_arr&time, model_arr_stddev, res, info_items)
  end if

  log_info(" wrote " + outfile + "." + file_type)

  ; ***************************************
  ; add meta data to plot (for reporting)
  ; ***************************************

  if (diag_script_info@region .eq. "Global") then
    domain = "global"
  else
    domain = "reg"
  end if

  if (anom .eq. "anom") then
    statistics = "anomaly"
  else
    statistics = "mean"
  end if

  caption = "Time series of the " + statistics + " for variable " \
    + varname + ", similar to IPCC AR5, fig. 9.8."
  contrib_authors = (/"cionni_irene", "righi_mattia", \
                      "wenzel_sabrina", "bock_lisa"/)

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 outfile + "." + file_type, \
                 caption, \
                 statistics, \
                 domain, \
                 "times", \
                 contrib_authors, \
                 (/"flato13ipcc", "jones13jgr"/), \
                 metadata_att_as_array(input_file_info, "filename"))

  leave_msg(DIAG_SCRIPT, "")

end
