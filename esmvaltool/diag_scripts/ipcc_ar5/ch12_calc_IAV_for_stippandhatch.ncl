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
;     List possible caveats or limitations of this diagnostic
;     Features to-be-implemented shall also be mentioned here
;
;  Modification history
;     20161219-A_lore_ru: remove seasonal cycle before std if seasonal
;     20161024-A_lore_ru: adapted to ESMValTool
;     20130501-A_sedl_ja: written for IPCC AR5 as get_natvar.ncl.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"
load "./diag_scripts/shared/statistics.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ;  'diag_script' and 'variables' are fetched from the
  ;  above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ;  'datasets', 'variables' and 'field_types' are fetched from the
  ;  above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
  write_references(diag_script, \  ; script name
                   "A_lore_ru", \  ; authors
                   "", \           ; contributors
                   "D_collins13ipcc",    \  ; diag_references
                   "",          \  ; obs_references
                   "P_crescendo")   ; proj_references
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
  if (isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if
end

begin
  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; get data from first dataset
  imod = 0  ; NCL array indicies start from zero
  log_debug("processing " + dataset_info@dataset(imod) + "_"\
            + dataset_info@exp(imod) + "_" \
            + dataset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)

  ; Check dimensions
  dims = getvardims(A0)

  if (dimsizes(dims) .lt. 3) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 3 or 4")
  end if
  idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  lat_name = dims(idx)

  if (field_type0 .eq. "T2Ms") then
    idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
    if (ismissing(idx)) then
      error_msg("f", diag_script, "", "no lon dimension")
    end if
    lon_name = dims(idx)
  end if

  index = ispan(0, dim_MOD - 1, 1)
  if diag_script_info@time_avg .eq. "seasonalclim" then
    dim_seas = 4
    diag_script_info@seasons = (/0, 1, 2, 3/)
  else
    dim_seas = 1
    diag_script_info@seasons = (/0/)
  end if

  if (field_type0 .eq. "T3Ms") then
    idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
    if (ismissing(idx)) then
      error_msg("f", diag_script, "", "no level dimension")
    end if
    lev_name = dims(idx)
  end if

  Fill = default_fillvalue(typeof(A0))
  if (field_type0 .eq. "T2Ms") then
    data1 = new((/dim_MOD, dim_seas, dimsizes(A0&lat), \
                 dimsizes(A0&lon)/), typeof(A0), Fill)
    if ((.not. isatt(diag_script_info, "iavmode")) .or. \
        (diag_script_info@iavmode .ne. "each")) then
      IAV_mmm = new((/dim_seas, dimsizes(A0&lat), \
                     dimsizes(A0&lon)/), typeof(A0), Fill)
    end if
  else if (field_type0 .eq. "T3Ms") then
    data1 = new((/dim_MOD, dim_seas, dimsizes(A0&$lev_name$), \
                 dimsizes(A0&lat)/), \
                typeof(A0), Fill)
    if ((.not. isatt(diag_script_info, "iavmode")) .or. \
        (diag_script_info@iavmode .ne. "each")) then
      IAV_mmm = new((/dim_seas, dimsizes(A0&$lev_name$), \
                     dimsizes(A0&lat)/), \
                    typeof(A0), Fill)
    end if
  end if
  end if

  do imod = 0, dim_MOD - 1
    if (imod .ne. 0) then
      delete(A0)
      A0 = read_data(imod, var0, field_type0)
      dims = getvardims(A0)
      idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
      if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lat dimension")
      end if
      lat_name = dims(idx)
      if (field_type0 .eq. "T2Ms") then
        idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
        if (ismissing(idx)) then
          error_msg("f", diag_script, "", "no lon dimension")
        end if
        lon_name = dims(idx)
      end if
      if (field_type0 .eq. "T3Ms") then
        idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
        if (ismissing(idx)) then
          error_msg("f", diag_script, "", "no level dimension")
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
        error_msg("w", diag_script, "", "time range too short, " + \
                  "less than 1 period covered, continue with " + \
                  "next model")
        continue
      end if
      yr_possible = (dimsizes(time) - 12 * 100) / (length_of_period * 12.)
      rest = yr_possible - nr_periods
      start_yrs = ispan(new_start_yr, end_yr_int, length_of_period)
      if (field_type0 .eq. "T2Ms") then
        data_tmp = new((/nr_periods, dim_seas, dimsizes(A0&$lat_name$), \
                        dimsizes(A0&$lon_name$)/), typeof(A0), Fill)
      else if (field_type0 .eq. "T3Ms") then
        data_tmp = new((/nr_periods, dim_seas, dimsizes(A0&$lev_name$), \
                        dimsizes(A0&$lat_name$)/), \
                       typeof(A0), Fill)
      end if
      end if

      if (diag_script_info@time_avg .eq. "seasonalclim") then
        ; remove seasonal cycle
        deseas_data = rmMonAnnCycTLL(A0)
        delete(A0)
        A0 = deseas_data
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
          data_tmp(per, 0, :, :) = \
            time_operations(A0, start_yrs(per), new_end_yr, "average", \
                            diag_script_info@time_avg, True)
        else
          data_tmp(per, :, :, :) = \
            time_operations(A0, start_yrs(per), new_end_yr, "average", \
                            diag_script_info@time_avg, True)
        end if
      end do
      data_dtr = dtrend_quadratic_msg_n(data_tmp, False, False, 0)

      data1(imod, :, :, :) = dim_stddev_n(data_dtr, 0)
      delete(start_yrs)
      delete(data_tmp)
      delete(data_dtr)
      if (diag_script_info@time_avg .eq. "seasonalclim") then
        delete(deseas_data)
      end if
    else
      data1_tmp = interannual_variability(A0, new_start_yr, end_yr_int, \
                                          diag_script_info@time_avg, \
                                          "quadratic")
      data1(imod, :, :, :) = data1_tmp
    end if
    delete(start_yr)
    delete(end_yr)
    delete(time)

    if ((isatt(diag_script_info, "iavmode")) .and. \
        (diag_script_info@iavmode .eq. "each")) then
      outfile = "IAV_" + dataset_info@exp(imod) + "_" + \
        dataset_info@dataset(imod) +\
        "_" + dataset_info@ensemble(imod) + "_" + field_type0 + "_" + \
        var0 + "_" + diag_script_info@time_avg + ".nc"
      file_exist = isfilepresent(outfile)
      if (file_exist .and. diag_script_info@overwrite .eq. False) then
        continue
      end if
      if (imod .eq.0) then
        if (field_type0 .eq. "T2Ms") then
          data1!0 = "models"
          data1!1 = "season"
          data1!2 = "lat"
          data1!3 = "lon"
          data1&lon = A0&lon
        else if (field_type0 .eq. "T3Ms") then
          data1!0 = "models"
          data1!1 = "season"
          data1!2 = lev_name
          data1&$lev_name$ = A0&$lev_name$
          data1!3 = "lat"
        end if
        end if
        if (diag_script_info@time_avg .eq. "seasonalclim") then
          data1&season = (/ispan(0, dim_seas - 1, 1)/)
        else
          data1&season = (/0/)
        end if
        data1&lat = A0&lat
        data1@diag_script = (/diag_script/)
        data1@var = "iav"
        if (isatt(variable_info, "long_name")) then
          data1@var_long_name = "iav of " + variable_info@long_name
        end if
        if (isatt(variable_info, "units")) then
          data1@var_units     = variable_info@units
        end if
      end if
      ; ###########################################
      ; # Output to netCDF                        #
      ; ###########################################
      if (dim_seas .eq. 1) then
        write_data = data1(imod, 0, :, :)
      else
        write_data = data1(imod, :, :, :)
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
      IAV_mmm(s, :, :) = \
        rm_single_dims(dim_avg_n(data1(:, s, :, :), 0) * sqrt(2.))
    end do

    if (field_type0 .eq. "T2Ms") then
      IAV_mmm!0 = "season"
      IAV_mmm!1 = "lat"
      IAV_mmm!2 = "lon"
      IAV_mmm&lon = A0&lon
    else if (field_type0 .eq. "T3Ms") then
      IAV_mmm!0 = "season"
      IAV_mmm!1 = lev_name
      IAV_mmm&$lev_name$ = A0&$lev_name$
      IAV_mmm!2 = "lat"
    end if
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
    outfile = "IAV_mmm_piControl_" + field_type0 + "_" \
      + var0 + "_" + diag_script_info@time_avg + ".nc"

    IAV_mmm@diag_script = (/diag_script/)
    IAV_mmm@var = "iav"
    if (isatt(variable_info, "long_name")) then
      IAV_mmm@var_long_name = "iav of " + variable_info@long_name
    end if
    if (isatt(variable_info, "units")) then
      IAV_mmm@var_units     = variable_info@units
    end if
    IAV_mmm@comment = dataset_info@dataset
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(IAV_mmm, work_dir + outfile)
  end if
end
