; #############################################################################
; WRAPPER SCRIPT FOR COLLECTING AND PLOTTING TIMESERIES
; Author: Lisa Bock (DLR, Germany)
; #############################################################################
; Description
;    Collects timeseries and their statistics previously calculated by
;    tsline.ncl and passes them to a new tsline plot
;
; Required diag_script_info attributes
;    styleset: as in diag_scripts/shared/plot/style.ncl functions
;
; Optional diag_script_info attributes
;    time_avg: type of time average (currently only "yearly" and "monthly" are
;              available).
;    ts_anomaly: calculates anomalies with respect to the defined period
;    ref_start: start year of reference period for anomalies
;    ref_end: end year of reference period for anomalies
;    region: name of domain
;    plot_units: variable unit for plotting
;    y_min: set min of y-axis
;    y_max: set max of y-axis
;    order: order in which experiments should be plotted
;    header: if true, region name as header
;    stat_shading: if true: shading of statistic range
;    ref_shading: if true: shading of reference period
;
; Caveats
;
; Modification history
;    20190204-A_bock_lisa: written.
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
  diag_script_base = basename(DIAG_SCRIPT)

  var_names = metadata_att_as_array(variable_info, "variable_group")
  var0 = "tas"
  ind_tas = ind(var_names .eq. "tas")
  info_items = select_metadata_by_name(input_file_info, var0)
  units0 = info_items[0]@units
  if(any(var_names .eq. "tasUnc1") .and. any(var_names .eq. "tasUnc2")) then
    obs_unc = True
    var1 = "tasUnc1"
    var2 = "tasUnc2"
    ind1 = ind(var_names .eq. var1)
    ind2 = ind(var_names .eq. var2)
    info_items1 = select_metadata_by_name(input_file_info, \
                                          variable_info[ind1]@short_name)
    info_items2 = select_metadata_by_name(input_file_info, \
                                          variable_info[ind2]@short_name)
  else
    obs_unc = False
  end if

  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)

  ; Check file type
  file_type = config_user_info@output_file_type
  if(ismissing(file_type)) then
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

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "time_avg", "monthly")
  set_default_att(diag_script_info, "ts_anomaly", "noanom")
  set_default_att(diag_script_info, "region", "Global")
  set_default_att(diag_script_info, "header", True)
  set_default_att(diag_script_info, "stat_shading", False)
  set_default_att(diag_script_info, "ref_shading", False)

  ; Determine time range
  start_year = diag_script_info@start_year
  end_year = diag_script_info@end_year
  all_years = ispan(start_year, end_year, 1)
  nyears = end_year - start_year + 1

  ref_start = diag_script_info@ref_start
  ref_end = diag_script_info@ref_end

  ; Anomaly
  if (isatt(diag_script_info, "ts_anomaly")) then
    anom = diag_script_info@ts_anomaly
    if (anom .eq. "anom") then
      if (isatt(diag_script_info, "ref_start") .and. \
          isatt(diag_script_info, "ref_end")) then
        ref_start = diag_script_info@ref_start
        ref_end   = diag_script_info@ref_end
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
      anom_ref = 0.

    end if
  else
    anom = "noanom"
  end if

  ; -----------------------------------------------
  ; Read pre-calculated multi model statistics
  ; -----------------------------------------------

  input_files = diag_script_info@input_files + "/" + \
    "tsline_mm_" + var0 + "_" + anom + "_stat.nc"
  input_files := tostring(input_files)

  ; Filter out non-existing files
  input_files := input_files(ind(isfilepresent(input_files)))

  nfiles = dimsizes(input_files)

  ; Loop over files in the list, read and append data
  do ii = 0, nfiles - 1

    log_info("Read in " + input_files(ii))

    data_temp = ncdf_read(input_files(ii), var0)

    log_info("Experiment: " + data_temp@experiment)

    if(isatt(diag_script_info, "order")) then
      ii_n = ind(data_temp@experiment .eq. diag_script_info@order)
      if (all(ismissing(ii_n))) then
        error_msg("f", DIAG_SCRIPT, "", "Collected experiment " + \
                  data_temp@experiment + " does not appear in " + \
                  "diag_script_info@order")
      elseif (ii_n .gt. nfiles-1)
        error_msg("f", DIAG_SCRIPT, "", "diag_script_info@order not " \
                  + "well defined (less experiments available)")
      end if
    else
      ii_n = ii
    end if

    if (ii .eq. 0) then
      time = data_temp&time
      ntime = dimsizes(time)

      dim_data   = (/nfiles, ntime/)
      data_mean  = new(dim_data, double)
      data_mean!0 = "experiment"
      data_mean!1 = "time"
      data_mean&experiment = new(nfiles, string, "exp")
      data_mean&time = time
      data_mean@var = var0
      if (isatt(diag_script_info, "ref")) then
        data_ref = new((/1, ntime/), double)
        data_ref!0 = "dataset"
        data_ref!1 = "time"
        data_ref&time = time
      end if
      data_stat = new((/2, nfiles, ntime/), double)
    end if

    ; Match time coordinate
    if (time_avg.eq."monthly") then
      date = cd_calendar(data_temp&time, -1)
    elseif (time_avg.eq."yearly") then
      date = data_temp&time
    end if
    idx1 = ind(date(0).eq.data_mean&time)
    idx2 = ind(date(dimsizes(date) - 1).eq.data_mean&time)
    delete(date)

    data_mean(ii_n, idx1:idx2)  = (/data_temp(:, 0)/)
    data_mean&experiment(ii_n) = data_temp@experiment

    if (ii .eq. 0) then
      data_mean@units = data_temp@units
      if (isatt(diag_script_info, "ref")) then
        data_ref(0, idx1:idx2)  = (/data_temp(:, 9)/)
        data_ref&dataset(0) = diag_script_info@ref
      end if
    end if

    ; add statistics to multi-model mean
    ; stddev
    data_stat(0, ii_n, idx1:idx2) = data_temp(:, 0) - data_temp(:, 1)
    data_stat(1, ii_n, idx1:idx2) = data_temp(:, 0) + data_temp(:, 1)

    delete(data_temp)

  end do

  if (obs_unc) then

    A_ref1 = read_data(info_items1[0])
    A_ref2 = read_data(info_items2[0])

    if (time_avg .eq. "yearly") then
      A_ref1 := time_operations(A_ref1, -1, -1, "average", "yearly", True)
      A_ref2 := time_operations(A_ref2, -1, -1, "average", "yearly", True)
    end if

    data_ref_stderr = new((/2, ntime/), double)
    data_ref_stderr!0 = "statistics"
    data_ref_stderr!1 = "time"
    data_ref_stderr&time = time

    data_ref_stderr(0, :) = (/data_ref(0, :) - A_ref1/)
    data_ref_stderr(1, :) = (/data_ref(0, :) + A_ref2/)

    delete(A_ref1)
    delete(A_ref2)

  else

    data_ref_stderr = 0

  end if

  ; -------------------------------------------
  ; NetCDF Output
  ; -------------------------------------------

  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)
  out_path1 = out_path + "tsline_collect_" + var0 + ".nc"
  data_mean@ncdf = out_path
  data_mean@diag_script = DIAG_SCRIPT
  data_mean@var = var0
  ncdf_outfile = ncdf_write(data_mean, out_path1)

  if (isatt(diag_script_info, "ref")) then
    out_path1 = out_path + "tsline_collect_" + var0 + "_ref.nc"
    data_ref@ncdf = out_path
    data_ref@diag_script = DIAG_SCRIPT
    data_ref@var = var0
    ncdf_outfile_ref = ncdf_write(data_ref, out_path1)
  end if

  ; -------------------------------------------
  ; Plotting
  ; -------------------------------------------

  ; Define workstation
  outfile = config_user_info@plot_dir + data_mean@var + "_" + \
    str_sub_str(diag_script_info@region, " ", "_") + \
    "_multimodel_" + anom + "_" + start_year + "-" + end_year
  wks = gsn_open_wks(file_type, outfile)

  ; Set resources
  res = True
  xmin = start_year
  xmax = end_year
  res@trXMinF = xmin
  res@trXMaxF = xmax
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
  if (isatt(diag_script_info, "y_min")) then
    res@trYMinF = diag_script_info@y_min
  end if
  if (isatt(diag_script_info, "y_max")) then
    res@trYMaxF = diag_script_info@y_max
  end if

  if (xmin .eq. 1850) then
    res@tmXBMode         = "Manual"
    res@tmXBTickSpacingF = 50
  elseif (xmin .eq. 1950) then
    res@tmXBMode         = "Manual"
    res@tmXBTickSpacingF = 20
  end if

  if (diag_script_info@header) then
    res@tiMainString = diag_script_info@region
  else
    res@tiMainString = ""
  end if

  res@tiMainString = diag_script_info@region

  if (isatt(variable_info[0], "long_name")) then
    if (var0 .eq. "tas") then
      varname = "Temperature"
    elseif (var0 .eq. "siconc")
      varname = "Sea Ice"
    else
      varname = variable_info[0]@long_name
    end if
  else
    varname = var0
  end if

  if (data_mean@units .eq. "degC") then
    units = "~F34~0~F~ C"
  else
    units = data_mean@units
  end if

  if (anom .eq. "anom") then
    res@tiYAxisString = varname + " Anomaly" + " (" + units + ")"
  else
    res@tiYAxisString = varname + " (" + units + ")"
  end if

  res0 = True

  if (isatt(diag_script_info, "ref")) then
    xy_line_collect(wks, data_mean, data_ref, data_ref_stderr, data_ref&time, \
                    data_stat, ref_start, ref_end, res, res0, info_items)
  else
    xy_line_collect(wks, data_mean, 0, 0, data_mean&time, \
                    data_stat, ref_start, ref_end, res, res0, info_items)
  end if

  log_info(" wrote " + outfile + "." + file_type)

  ; ***************************************
  ; add meta data to plot (for reporting)
  ; ***************************************

  if (diag_script_info@region .eq. "Global") then
    domain = (/"global"/)
  else
    domain = (/"reg"/)
  end if

  if (anom .eq. "anom") then
    statistics = (/"anomaly"/)
  else
    statistics = (/"mean"/)
  end if

  caption = "Multi model mean time series of the " + statistics + \
            " for variable " + varname + ", similar to IPCC AR5, fig. 10.21."

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 outfile + "." + file_type, \
                 caption, \
                 statistics, \
                 domain, \
                 "times", \
                 (/"bock_lisa"/), \
                 (/"flato13ipcc", "jones13jgr"/), \
                 metadata_att_as_array(input_file_info, "filename"))

  leave_msg(DIAG_SCRIPT, "")

end
