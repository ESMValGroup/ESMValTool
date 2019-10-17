; #############################################################################
; INTERANNUAL VARIABILITY MULTI-MODEL MEAN FOR STIPPLING SIGNIFICANCE IPCCCH12
; Author: Ruth Lorenz (ETH Zurich, Switzerland)
; CRESCENDO project
; #############################################################################
;
;  Description
;     Calculate Interannual variability from piControl runs for plotting
;     significance with stippling and hatching, save for each model
;     (iavmode: "each") or multi-model mean only.
;     Can either be calculated over full time period of piControl run
;     (no periodlength given)
;     or first averaged annually/seasonally over periodlength and then standard
;     deviation calculated based on averaged periods.
;
;  Required diag_script_info attributes (diagnostics specific)
;     time_avg:  time period to calculate IAV, e.g. annualclim, seasonalclim
;
;  Optional diag_script_info attributes (diagnostic specific)
;     iavmode: calculate multi-model mean of IAV over all models or calculate
;              and save IAV for each individual model? (mmm, each)
;     periodlength: length of periods to calculate IAV accross, depends on
;                   period lengths used in
;                   ch12_calc_map_diff_mmm_stippandhatch.ncl
;                   if not given whole time period calculated at once
;
;  Caveats
;     Needs lots of memory for 3D ocean variables
;
;  Modification history
;     20161219-A_lorenz_ruth: remove seasonal cycle before std if seasonal
;     20161024-A_lorenz_ruth: adapted to ESMValTool
;     20130501-A_sedlacek_jan: written for IPCC AR5 as get_natvar.ncl.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/statistics.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(DIAG_SCRIPT, "")

  ;  'datasets', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  ; Check required diag_script_info attributes
  req_atts = (/"time_avg"/)
  exit_if_missing_atts(diag_script_info, req_atts)
end

begin
  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if (isvar("MyParam")) then
    delete(MyParam)
  end if
  if (isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if
end

begin
  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; get data from first dataset
  imod = 0  ; NCL array indicies start from zero
  log_debug("processing " + info_items[imod]@dataset + "_" \
            + info_items[imod]@exp + "_" \
            + info_items[imod]@ensemble)

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(info_items[imod])

  ; Check dimensions
  dims = getvardims(A0)
  ndim = dimsizes(dims)

  if (ndim .lt. 3) then
    error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + " dimensions, " + \
              "need 3 or 4")
  end if
  idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
  end if
  lat_name = dims(idx)

  idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
  end if
  lon_name = dims(idx)

  if (ndim .gt. 3) then
    idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
    if (ismissing(idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no level dimension")
    end if
    lev_name = dims(idx)
  end if

  index = ispan(0, dim_MOD - 1, 1)
  if (diag_script_info@time_avg .eq. "seasonalclim") then
    dim_seas = 4
    diag_script_info@seasons = (/0, 1, 2, 3/)
  else
    dim_seas = 1
    diag_script_info@seasons = (/0/)
  end if

  Fill = default_fillvalue(typeof(A0))
  if (ndim .eq. 3) then
    data1 = new((/dim_MOD, dim_seas, dimsizes(A0&lat), \
                 dimsizes(A0&lon)/), typeof(A0), Fill)
    if ((.not. isatt(diag_script_info, "iavmode")) .or. \
        (diag_script_info@iavmode .ne. "each")) then
      IAV_mmm = new((/dim_seas, dimsizes(A0&lat), \
                     dimsizes(A0&lon)/), typeof(A0), Fill)
    end if
  elseif (ndim .eq. 4) then
    data1 = new((/dim_MOD, dim_seas, dimsizes(A0&$lev_name$), \
                 dimsizes(A0&lat), dimsizes(A0&lon)/), \
                typeof(A0), Fill)
    if ((.not. isatt(diag_script_info, "iavmode")) .or. \
        (diag_script_info@iavmode .ne. "each")) then
      IAV_mmm = new((/dim_seas, dimsizes(A0&$lev_name$), \
                     dimsizes(A0&lat)/), \
                    typeof(A0), Fill)
    end if
  end if

  do imod = 0, dim_MOD - 1
    if (imod .ne. 0) then
      delete(A0)
      A0 = read_data(info_items[imod])
      dims = getvardims(A0)
      idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
      end if
      lat_name = dims(idx)

      idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
      end if
      lon_name = dims(idx)

      if (ndim .eq. 4) then
        idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
        if (ismissing(idx)) then
          error_msg("f", DIAG_SCRIPT, "", "no level dimension")
        end if
        lev_name = dims(idx)
      end if
    end if

    ; Calculate IAV (function in ./diag_scripts/lib/ncl/statistics.ncl)
    ; cut the first 100 yr
    time = A0&time
    start_yr = cd_calendar(time(0), 0)
    ind_end = dimsizes(time) - 1
    end_yr = cd_calendar(time(ind_end), 0)
    if ((end_yr(0, 0) - start_yr(0, 0) + 1) .lt. 500) then
      log_info("warning: Control run is less than 500 years.")
    end if
    new_start_yr = toint((start_yr(0, 0) + 100))
    end_yr_int = toint(end_yr(0, 0))
    if (isatt(diag_script_info, "periodlength")) then
      length_of_period = toint(diag_script_info@periodlength)
      nr_periods = toint(floor((dimsizes(time) - 12 * 100) / \
                               (length_of_period * 12.)))
      if (nr_periods .lt. 1) then
        error_msg("w", DIAG_SCRIPT, "", "time range too short, " + \
                  "less than 1 period covered, continue with " + \
                  "next model")
        continue
      end if
      yr_possible = (dimsizes(time) - 12 * 100) / (length_of_period * 12.)
      rest = yr_possible - nr_periods
      start_yrs = ispan(new_start_yr, end_yr_int, length_of_period)
      if (ndim .eq. 3) then
        data_tmp = new((/nr_periods, dim_seas, dimsizes(A0&$lat_name$), \
                        dimsizes(A0&$lon_name$)/), typeof(A0), Fill)
      elseif (ndim .eq. 4) then
        data_tmp = new((/nr_periods, dim_seas, dimsizes(A0&$lev_name$), \
                        dimsizes(A0&$lat_name$), dimsizes(A0&$lon_name$)/), \
                       typeof(A0), Fill)
      end if

      do per = 0, nr_periods - 1
        if ((rest .gt. 0) .and. (per .eq. nr_periods - 1)) then
          new_end_yr = end_yr_int
        else
          new_end_yr = toint(start_yrs(per) + length_of_period - 1.0)
        end if
        ;  calculate seasonal/annual averages over periods
        if ((dim_seas .eq. 1) .and. \
            (diag_script_info@time_avg .eq. "annualclim")) then
          if (ndim .eq. 3) then
            data_tmp(per, 0, :, :) = \
              time_operations(A0, start_yrs(per), new_end_yr, "average", \
                              diag_script_info@time_avg, True)
          elseif (ndim .eq. 4) then
            data_tmp(per, 0, :, :, :) = \
              time_operations(A0, start_yrs(per), new_end_yr, "average", \
                              diag_script_info@time_avg, True)
          end if
        else
          if (ndim .eq. 3) then
            data_tmp(per, :, :, :) = \
              time_operations(A0, start_yrs(per), new_end_yr, "average", \
                              diag_script_info@time_avg, True)
          elseif (ndim .eq. 4) then
            data_tmp(per, :, :, :, :) = \
              time_operations(A0, start_yrs(per), new_end_yr, "average", \
                              diag_script_info@time_avg, True)
          end if
        end if
      end do
      data_dtr = dtrend_quadratic_msg_n(data_tmp, False, False, 0)
      if (typeof(data_dtr) .ne. typeof(data1)) then
        if ((typeof(data_dtr) .eq. "double") .and. \
            (typeof(data1) .eq. "float")) then
          tmp = data_dtr
          delete(data_dtr)
          data_dtr = doubletofloat(tmp)
        elseif ((typeof(data_dtr) .eq. "float") .and. \
                (typeof(data1) .eq. "double")) then
          tmp = data_dtr
          delete(data_dtr)
          data_dtr = floattodouble(tmp)
        else
          error_msg("f", DIAG_SCRIPT, "", "Type conversion issue, " + \
                    "data_dtr has a different type than data1 which is " + \
                    "neither float or double.")
        end if
      end if
      if (ndim .eq. 3) then
        data1(imod, :, :, :) = dim_stddev_n(data_dtr, 0)
      elseif (ndim .eq. 4) then
        data1(imod, :, :, :, :) = dim_stddev_n(data_dtr, 0)
      end if
      delete(start_yrs)
      delete(data_tmp)
      delete(data_dtr)
    else
      data1_tmp = interannual_variability(A0, new_start_yr, end_yr_int, \
                                          diag_script_info@time_avg, \
                                          "quadratic")
      if (ndim .eq. 3) then
        data1(imod, :, :, :) = data1_tmp
      elseif (ndim .eq. 4) then
        data1(imod, :, :, :, :) = data1_tmp
      end if
    end if
    delete(start_yr)
    delete(end_yr)
    delete(time)

    if ((isatt(diag_script_info, "iavmode")) .and. \
        (diag_script_info@iavmode .eq. "each")) then
      outfile = "IAV_" + info_items[imod]@exp + "_" + \
        info_items[imod]@dataset + "_" + info_items[imod]@ensemble + "_" + \
        var0 + "_" + diag_script_info@time_avg + ".nc"
      file_exist = isfilepresent(outfile)
      if (file_exist .and. diag_script_info@overwrite .eq. False) then
        continue
      end if
      if (imod .eq.0) then
        if (ndim .eq. 3) then
          data1!0 = "models"
          data1!1 = "season"
          data1!2 = "lat"
          data1!3 = "lon"
          data1&lon = A0&lon
        elseif (ndim .eq. 4) then
          data1!0 = "models"
          data1!1 = "season"
          data1!2 = lev_name
          data1&$lev_name$ = A0&$lev_name$
          data1!3 = "lat"
          data1!4 = "lon"
          data1&lon = A0&lon
        end if
        if (diag_script_info@time_avg .eq. "seasonalclim") then
          data1&season = (/ispan(0, dim_seas - 1, 1)/)
        else
          data1&season = (/0/)
        end if
        data1&lat = A0&lat
        data1@diag_script = (/DIAG_SCRIPT/)
        data1@var = "iav"
        if (isatt(variable_info[0], "long_name")) then
          data1@var_long_name = "iav of " + variable_info[0]@long_name
        end if
        if (isatt(variable_info[0], "units")) then
          data1@units     = variable_info[0]@units
        end if
      end if
      ; ###########################################
      ; # Output to netCDF                        #
      ; ###########################################
      if (dim_seas .eq. 1) then
        if (ndim .eq. 3) then
          write_data = data1(imod, 0, :, :)
        elseif (ndim .eq. 4) then
          write_data = data1(imod, 0, :, :, :)
        end if
      else
        if (ndim .eq. 3) then
          write_data = data1(imod, :, :, :)
        elseif (ndim .eq. 4) then
          write_data = data1(imod, :, :, :, :)
        end if
      end if
      ; Function in ~/interface_scripts/auxiliary.ncl
      ncdf_outfile = ncdf_write(write_data, work_dir + outfile)
    end if
  end do  ; imod

  if ((.not. isatt(diag_script_info, "iavmode")) .or. \
      (diag_script_info@iavmode .ne. "each")) then
    ; ###########################################
    ; # Calculate multi-model mean of IAV       #
    ; ###########################################
    do s = 0, dim_seas - 1
      if (ndim .eq. 3) then
        IAV_mmm(s, :, :) = \
          rm_single_dims(dim_avg_n(data1(:, s, :, :), 0) * sqrt(2.))
      elseif (ndim .eq. 4) then
        ; Calculate the zonal average
        tmp_zon = dim_avg_n_Wrap(data1, 4)  ; data1(mon, seas, lev, lat, lon)
        ; Calculate multi-model mean
        IAV_mmm(s, :, :) = \
          rm_single_dims(dim_avg_n(tmp_zon(:, s, :, :), 0) * sqrt(2.))
      end if
    end do

    if (ndim .eq. 3) then
      IAV_mmm!0 = "season"
      IAV_mmm!1 = "lat"
      IAV_mmm!2 = "lon"
      IAV_mmm&lon = A0&lon
    elseif (ndim .eq. 4) then
      IAV_mmm!0 = "season"
      IAV_mmm!1 = lev_name
      IAV_mmm&$lev_name$ = A0&$lev_name$
      IAV_mmm!2 = "lat"
    end if
    if (isatt(diag_script_info, "seasons")) then
      IAV_mmm&season = (/ispan(0, dim_seas - 1, 1)/)
    else
      IAV_mmm&season = (/0/)
    end if
    IAV_mmm&lat = A0&lat

    ; ###########################################
    ; # Output to netCDF                        #
    ; ###########################################
    outfile = "IAV_mmm_piControl_" + \
      var0 + "_" + diag_script_info@time_avg + ".nc"

    IAV_mmm@diag_script = (/DIAG_SCRIPT/)
    IAV_mmm@var = "iav"
    if (isatt(variable_info[0], "long_name")) then
      IAV_mmm@var_long_name = "iav of " + variable_info[0]@long_name
    end if
    if (isatt(variable_info[0], "units")) then
      IAV_mmm@units     = variable_info[0]@units
    end if
    IAV_mmm@comment = metadata_att_as_array(info_items, "dataset")
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(IAV_mmm, work_dir + outfile)
  end if

  ; collect meta-data
  nc_file = ncdf_outfile
  plot_file = "None"
  caption = "Inter-annual variability based on piControl runs."
  statistics = ("var")
  domains = ("global")
  plot_types = ("other")
  authors = (/"lorenz_ruth"/)
  references = (/"collins13ipcc"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, plot_file, caption, statistics, domains, \
                 plot_types, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")

end
