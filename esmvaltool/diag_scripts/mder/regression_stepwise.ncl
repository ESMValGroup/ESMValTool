; #############################################################################
; mder/regression_stepwise.ncl
; #############################################################################
; MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; Author: Sabrina Wenzel (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     This routine computes weights for a set of models listed in the
;     recipe recipes/recipe_wenzel16jclim.xml following the multiple diagnostic
;     ensemble regression (MDER) described by Karpechko et al. 2013.
;     The script calls a variable of interest from different models (which
;     shall be weighted) and additionally several diagnostics which are
;     specified in the recipe for computig the weights.
;     Calls three plotscripts:
;     + Plots the absolute correlation coefficient between the variable of
;       interest and present-day diagnostic.
;     + Scatter plot of the climatological mean of the variable of interest
;       and the result of the MDER.
;     + Timeseries of the variable of interest including the multi-model
;       mean, the weighted and unweighted model ensemble.
;
; Required diag_script_info attributes:
;     p_time: time range for projections.
;     p_step: averaging time range for target period.
;     scal_time: a year to scale the timeseries with.
;     time_oper: select for operation in time_operation function.
;     time_opt: select for option in time_operation function.
;     calc_type: trend/pos.
;     wregr: exit criterion for stepwise regression (True/False).
;     wexit: exit criterion for stepwise regression (True/False).
;
; Optional diag_script_info attributes:
;     units: units of target variable.
;     loglin: perform log-linear regression (default: False).
;     smooth: smoothes time period with 1-2-1 filter.
;     iter: iteration for smoothing.
;     cross_validation_mode: perform cross validation.
;     p_value: p value for F-test in stepwise regression (default: p=0.05).
;     plot_mder_time_series: plot whole time series of MDER results.
;
; Caveats:
;     Adds input_file_info@units.
;     Modifies input_file_info@dataset and input_file_info@ensemble.
;     Modifies variable_info@units.
;
; Modification history
;     20180806_A-schl_ma: ported to v2.0.
;     201808??_A-senf_da: added option to plot log-linear regression.
;     201303??_A-wenz_sa: written.
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/latlon.ncl"

load "./diag_scripts/shared/plot/mder.ncl"
load "./diag_scripts/shared/plot/aux_plotting.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

load "./diag_scripts/shared/mder.ncl"
; load "./diag_scripts/aux/SeaIce/seaice.ncl"


begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Variable
  var0 = variables(0)
  target_var = var0
  target_var@units = variable_info@units
  field_type0 = field_types(0)
  if (var0 .eq. "sic") then
    target_var = "SIE"  ; sea ice extent
    target_var@units = diag_script_info@units
  end if

  ; Check if variable is multi-dimensional
  is_multi_dim = .not. isStrSubset(field_type0, "0")

  ; Extract observations from dataset list
  IDX_OBS = ind(input_file_info@project .eq. "OBS")
  DIM_OBS = dimsizes(IDX_OBS)
  OBS_NAME = input_file_info@dataset(IDX_OBS)

  ; Extract historical simulations from datasets
  IDX_HDAT = ind(input_file_info@exp .eq. "historical" .and. \
                 input_file_info@project .ne. "OBS")
  DIM_HDAT = dimsizes(IDX_HDAT)

  ; Extract individual future projections from datasets (= unique models)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens) then
    idx_datasets = UNIQ(input_file_info@dataset(DIM_OBS + DIM_HDAT:)) + \
                   DIM_OBS + DIM_HDAT
  else
    idx_datasets = ispan(DIM_OBS + DIM_HDAT, \
                         dimsizes(input_file_info@dataset) - 1, 1)
    idx_hist = UNIQ(input_file_info@dataset(DIM_OBS:DIM_OBS + DIM_HDAT))
    if (dimsizes(idx_hist) .gt. 1) then
      all_datasets = unique_labels_min((/"dataset", "ensemble"/))
      input_file_info@dataset := all_datasets
      delete(all_datasets)
    end if
    delete(idx_hist)
  end if

  ; Get dataset variables
  DIM_DAT = dimsizes(input_file_info@dataset(idx_datasets))
  input_file_info@ensemble := input_file_info@ensemble(idx_datasets)

end

begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Read years
  ny_max = max(toint(input_file_info@end_year))
  ny_min = min(toint(input_file_info@start_year(DIM_OBS:)))
  n_max = ny_max - ny_min + 1

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Required attributes
  req_atts = (/"p_time", "p_step", "scal_time", "time_oper", "time_opt", \
              "calc_type"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  p_time = diag_script_info@p_time
  p_step = diag_script_info@p_step
  scal_time = diag_script_info@scal_time
  time_oper = diag_script_info@time_oper
  time_opt = diag_script_info@time_opt
  calc_type = diag_script_info@calc_type

  ; For multi-dimensional variables, 'region' and 'area_oper' are required
  if (is_multi_dim) then
    exit_if_missing_atts(diag_script_info, (/"region", "area_oper"/))
    area_oper = diag_script_info@area_oper
    if (isstring(diag_script_info@region)) then
      region = select_region(diag_script_info@region)
    else if (dimsizes(diag_script_info@region .eq. 2)) then
      region = (/diag_script_info@region(0), \
                diag_script_info@region(1), \
                0.0, \
                360.0/)
      region@name = flt2string(region(0)) + "_-_" + flt2string(region(1))
      region@DM_tag = "DM_reg"
    else
      error_msg("f", diag_script, "", diag_script_info@region + \
                " is an invalid region type")
    end if
    end if
  end if

  ; Optional attributes
  if (.not. isatt(diag_script_info, "units")) then
    diag_script_info@units = variable_info@units
  end if
  loglin = False
  if (isatt(diag_script_info, "loglin")) then
    loglin = diag_script_info@loglin
  end if
  smooth = False
  if (isatt(diag_script_info, "smooth")) then
    smooth = diag_script_info@smooth
    if (isatt(diag_script_info, "iter")) then
      iter = diag_script_info@iter
    else
      error_msg("f", diag_script, "", "requires attribute iter in recipe " + \
                "when smooth is set to True")
    end if
  end if
  cross_validation_mode = False
  if (isatt(diag_script_info, "cross_validation_mode")) then
    cross_validation_mode = diag_script_info@cross_validation_mode
  end if

  ; Output dirs
  plot_dir = config_user_info@plot_dir
  work_dir = config_user_info@work_dir
  system("mkdir -p " + plot_dir)
  system("mkdir -p " + work_dir)
  plot_num = integertochar(ispan(97, 122, 1))
  filename = input_file_info@project(DIM_OBS) + "_" + target_var + \
    "-" + calc_type + "_" + input_file_info@exp(DIM_OBS + DIM_HDAT) + \
    "_" + (p_step + 1) + "ystep"

  ; Output arrays
  xx = new((/DIM_DAT, dimsizes(p_time)/), float)
  xx_val = new((/DIM_DAT, dimsizes(p_time), 2/), float)
  YY = new((/DIM_DAT + 3, n_max/), float)         ; incl. uMMM and wMMM
  t_mder_errors = new((/2, n_max/), float)
  mder_errors = new((/2, n_max/), float)
  yy_obs = new((/DIM_OBS, n_max/), float)
  weights_t = new((/1, DIM_DAT, dimsizes(p_time)/), float)
  cfint_t = new((/DIM_DAT, dimsizes(p_time)/), float)
  coeff_y_est = new((/DIM_DAT, dimsizes(p_time)/), double)  ; Type of coef
  Y_WEIGHT = new((/dimsizes(p_time)/), double)              ; Type of coef
  MCORR = new((/dimsizes(p_time)/), double)                 ; Type of coef

  ; ---------------------------------------------------------------------------
  ; Read and preprocess model data
  ; ---------------------------------------------------------------------------

  ; Call the variable of interest and save in one array YY
  do i = 0, DIM_DAT - 1
    idat = idx_datasets(i)
    log_debug(input_file_info@dataset(idat))

    ; Call projection simulations for regression analysis and time series,
    ; scale time series and average over 10-year period
    tmp = read_data(idat, var0, field_type0)
    idx_hist = ind(input_file_info@dataset(idat) .eq. \
                   input_file_info@dataset(:DIM_OBS + DIM_HDAT - 1))
    idx_fut = ind(input_file_info@dataset(idat) .eq. \
                  input_file_info@dataset(DIM_OBS + DIM_HDAT:)) + \
      DIM_OBS + DIM_HDAT

    ; Process multi-dimensional variables
    if (is_multi_dim) then
      tmp := area_operations(tmp, region(0), region(1), region(2), \
                             region(3), area_oper, True)
    end if

    if (dimsizes(idx_hist) .ge. 1) then
      ; Call historical simulations for time series
      tmp_hist = read_data(idx_hist(0), var0, field_type0)

      ; Process multi-dimensional variables
      if (is_multi_dim) then
        tmp_hist := area_operations(tmp_hist, region(0), region(1), \
                                    region(2), region(3), area_oper, True)
      end if

      ; Average over ensemble members if more than one exists (see recipe)
      if (avgens) then
        if (dimsizes(idx_fut) .gt. 1) then
          if (dimsizes(idx_fut) .eq. dimsizes(idx_hist)) then
            tmp_avg_fut = new((/dimsizes(idx_fut), dimsizes(tmp&time)/), \
                              typeof(tmp))
            tmp_avg_hist = new((/dimsizes(idx_hist), \
                                dimsizes(tmp_hist&time)/), typeof(tmp_hist))
            copy_VarCoords(tmp, tmp_avg_fut(0, :))
            copy_VarCoords(tmp_hist, tmp_avg_hist(0, :))
            tmp_avg_fut(0, :) = tmp
            tmp_avg_hist(0, :) = tmp_hist
            delete([/tmp, tmp_hist/])

            ; Loop over ensemble members, read the data and save in array
            do iens = 1, dimsizes(idx_fut) - 1
              tmp_fut = read_data(idx_fut(iens), var0, field_type0)
              tmp_hist = read_data(idx_hist(iens), var0, field_type0)

              ; Process multi-dimensional variables
              if (is_multi_dim) then
                tmp_fut := area_operations(tmp_fut, region(0), region(1), \
                                           region(2), region(3), area_oper, \
                                           True)
                tmp_hist := area_operations(tmp_hist, region(0), region(1), \
                                            region(2), region(3), area_oper, \
                                            True)
              end if

              tmp_avg_fut(iens, :) = (/tmp_fut/)
              tmp_avg_hist(iens, :) = (/tmp_hist/)
              delete([/tmp_fut, tmp_hist/])
            end do

            ; Average over ensemble members
            tmp = tofloat(dim_avg_n_Wrap(tmp_avg_fut, 0))
            tmp_hist = tofloat(dim_avg_n_Wrap(tmp_avg_hist, 0))
            copy_VarCoords(tmp_avg_fut(0, :), tmp)
            copy_VarCoords(tmp_avg_hist(0, :), tmp_hist)
            delete([/tmp_avg_fut, tmp_avg_hist/])
          else
            error_msg("f", diag_script, "", "inconsistent ensemble member " + \
                      "size for dataset " + input_file_info@dataset(idat))
          end if
        end if
      end if

      ; Adapt future values for cumulative quantities
      if (var0 .eq. "cumnbp" .or. var0 .eq. "cumfgco2") then
        do tt = 0, dimsizes(tmp&time) - 1
          tmp(tt) = tmp(tt) + tmp_hist(dimsizes(tmp_hist&time) - 1)
        end do
      end if

      ; Append future array to historical array
      tmpy = array_append_record(tmp_hist, tmp, 0)
      delete([/tmp, tmp_hist, idx_hist, idx_fut/])
    end if

    ; Convert units if desired
    if (isatt(variable_info, "plot_units")) then
      tmpy := convert_units(tmpy, variable_info@plot_units)
      variable_info@units = variable_info@plot_units
    end if

    ; Adapt scal_time(1) (necessary of larger than end_year)
    if (scal_time(1) .gt. \
        toint(input_file_info@end_year(idat - DIM_HDAT))) then
      scal_time(1) = toint(input_file_info@end_year(idat - DIM_HDAT))
    end if

    ; Calculate average over starting period for trend calculation
    trend_base = time_operations(tmpy, scal_time(0), scal_time(1), \
                                 time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      trend_base := time_operations(trend_base, scal_time(0), scal_time(1), \
                                    "average", "yearly", True)
    end if
    trend_base = tofloat(avg(trend_base))

    ; Loop over different target periods
    do pt = 0, dimsizes(p_time) - 1
      ; Average over var time series
      temp_mean = time_operations(tmpy, p_time(pt), p_time(pt) + p_step, \
                                  time_oper, time_opt, True)
      if (time_oper .eq. "extract") then
        temp_mean := time_operations(temp_mean, p_time(pt), \
                                     p_time(pt) + p_step, "average", \
                                     "yearly", True)
      end if

      ; Perform calculation type
      if (calc_type .eq. "trend") then
        xx(i, pt) = tofloat(avg(temp_mean) - trend_base)
        xx_val(i, pt, 0) = tofloat(avg(temp_mean) - trend_base)

        ; FIXME: Why substracting base value from error?
        xx_val(i, pt, 1) = tofloat(stddev(temp_mean) / \
                                   sqrt(dimsizes(temp_mean)) - trend_base)
        variable_info@units = variable_info@units + "/dec"
      else if (calc_type .eq. "pos") then
        xx(i, pt) = tofloat(avg(temp_mean))
        xx_val(i, pt, 0) = tofloat(avg(temp_mean))
        xx_val(i, pt, 1) = tofloat(stddev(temp_mean) / \
                                   sqrt(dimsizes(temp_mean)))
      else if (calc_type .eq. "int") then
        integrated_vals = temp_mean
        do tt = 0, dimsizes(temp_mean) - 2
          integrated_vals(tt + 1) = integrated_vals(tt) + temp_mean(tt + 1)
        end do
        xx(i, pt) = tofloat(avg(integrated_vals))
        xx_val(i, pt, 0) = tofloat(avg(integrated_vals))
        xx_val(i, pt, 1) = tofloat(stddev(integrated_vals) / \
                                   sqrt(dimsizes(integrated_vals)))
        delete([/integrated_vals, tt/])

      else
        error_msg("f", diag_script, "", "calc_type " + calc_type + \
                  " not supported")
      end if
      end if
      end if
      delete(temp_mean)
    end do

    ; Extract season for var time series
    temp_mean = time_operations(tmpy, -1, -1, time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      temp_mean := time_operations(temp_mean, -1, -1, "average", "yearly", \
                                   True)
    end if

    ; Smoothing var time series with 1-2-1 filter
    if (smooth) then
      temp_mean := filter121(temp_mean, iter)
    end if
    if (calc_type .eq. "trend") then
      YY(i + 3, 0:dimsizes(temp_mean) - 1) = tofloat(temp_mean - trend_base)
    else if (calc_type .eq. "pos") then
      YY(i + 3, 0:dimsizes(temp_mean) - 1) = tofloat(temp_mean)
    else if (calc_type .eq. "int") then
      YY(i + 3, 0) = temp_mean(0)
      do tt = 0, dimsizes(temp_mean) - 2
        YY(i + 3, tt + 1) = YY(i + 3, tt) + temp_mean(tt + 1)
      end do
    end if
    end if
    end if
    delete([/trend_base, tmpy, temp_mean/])
  end do

  ; Add attributes to models time series
  YY!1 = "time"
  YY&time = ispan(ny_min, ny_max, 1)      ; temp_mean&year
  YY&time@calendar = "standard"

  ; Save unique dataset names (modifies input_file_info@dataset!)
  input_file_info@dataset := input_file_info@dataset(idx_datasets)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess observational data
  ; ---------------------------------------------------------------------------

  ; Call OBS of the variable of interest and save in one array yy_obs
  do iobs = 0, DIM_OBS - 1
    log_debug(OBS_NAME(iobs))

    ; Call projection simulations for regression analysis and time series,
    ; scale time series and average over 10-year period
    tmpy = read_data(iobs, var0, field_type0)

    ; Process multi-dimensional variables
    if (is_multi_dim) then
      tmpy := area_operations(tmpy, region(0), region(1), region(2), \
                              region(3), area_oper, True)
    end if

    ; Convert units if desired
    if (isatt(variable_info, "plot_units")) then
      tmpy := convert_units(tmpy, variable_info@plot_units)
    end if

    ; Perform calculation type
    ymin = min(toint(input_file_info@start_year(iobs))) - ny_min
    obs_mean = time_operations(tmpy, -1, -1, time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      obs_mean := time_operations(obs_mean, -1, -1, "average", "yearly", True)
    end if
    if (calc_type .eq. "trend") then
      obs_base = time_operations(tmpy, scal_time(0), -1, "average", "yearly", \
                                 True)
      obs_base = tofloat(avg(obs_base))
      yy_obs(iobs, ymin:ymin + dimsizes(obs_mean) - 1) = \
        tofloat(obs_mean - obs_base)
      delete(obs_base)
    else if (calc_type .eq. "pos") then
      yy_obs(iobs, ymin:ymin + dimsizes(obs_mean) - 1) = tofloat(obs_mean)
    else if (calc_type .eq. "int") then
      yy_obs(iobs, ymin) = obs_mean(0)
      do tt = 0, dimsizes(obs_mean) - 2
        yy_obs(iobs, ymin + tt + 1) = yy_obs(iobs, ymin + tt) + \
          obs_mean(tt + 1)
      end do
      delete(tt)
    end if
    end if
    end if
    delete([/tmpy, ymin, obs_mean/])
  end do

  ; Add attributes to observations time series
  yy_obs!1 = "time"
  yy_obs&time = YY&time
  yy_obs&time@calendar = "standard"
  yy_obs@name = OBS_NAME

  ; ---------------------------------------------------------------------------
  ; Call diagnostics
  ; ---------------------------------------------------------------------------

  ; Get input directory (using ancestor diagnostics)
  input_files = diag_script_info@input_files
  input_dirs = input_files(ind(input_files .ne. \
                               str_match(input_files, "_info.ncl")))

  ; Abort when no diagnostics are given
  if (dimsizes(input_dirs) .lt. 1) then
    error_msg("f", diag_script, "", "no diagnositcs given (use ancestors " + \
              "key in recipe")
  end if

  ; Get diagnostics
  file_list = new(1, string)
  do idx_dir = 0, dimsizes(input_dirs) - 1
    file_list := array_append_record(\
      file_list, \
      systemfunc("ls " + input_dirs(idx_dir) + "/*.nc"), \
      0)
  end do
  file_list := file_list(1:)

  ; Create diagnostic arrays
  DIM_DIAG = dimsizes(file_list)
  diag_a = new((/DIM_DIAG, DIM_DAT/), float)
  diag_a_std = new((/DIM_DIAG, DIM_DAT/), float)
  obs0 = new((/DIM_DIAG/), float)
  obs0_std = new((/DIM_DIAG/), float)
  wdiag = new((/DIM_DIAG/), string)

  log_info("")
  ; Loop over files in the list and append data
  do idx_diag = 0, DIM_DIAG - 1
    data_temp = ncdf_read(file_list(idx_diag), "diag")

    ; Loop over datasets and save diagnostic data
    do i = 0, DIM_DAT - 1
      do idat = 0, dimsizes(data_temp&datasets) - 1

        ; Check if data_temp dataset labels are different
        if (isStrSubset(input_file_info@dataset(i), "_") .and. .not. \
            isStrSubset(data_temp&datasets(idat), "_")) then
          dataset_name = str_split(input_file_info@dataset(i), "_")
          input_file_info@dataset(i) = dataset_name(0)
          delete(dataset_name)
        end if

        ; Add data if datasets match
        if (input_file_info@dataset(i) .eq. data_temp&datasets(idat)) then
          diag_a(idx_diag, i) = data_temp(0, idat, 0)
          diag_a_std(idx_diag, i) = data_temp(0, idat, 1)
          if (isdim(data_temp, "diagnostics")) then
            wdiag(idx_diag) = tostring(data_temp&diagnostics(0))
          else
            error_msg("f", diag_script, "", "Cannot read diagnostic name " + \
                      "from data " + data_temp + ", need dimension " + \
                      "'diagnostics'")
          end if
        end if
      end do
    end do

    ; Save reference data
    if (isatt(data_temp, "ref_dataset"))then
      if (any(data_temp&datasets .eq. data_temp@ref_dataset)) then
        idx_ref = ind(data_temp&datasets .eq. data_temp@ref_dataset)
        obs0(idx_diag) = data_temp(0, idx_ref, 0)
        obs0_std(idx_diag) = data_temp(0, idx_ref, 1)
        delete(idx_ref)
      end if
    end if
    delete(data_temp)

    ; FIXME: degrees of freedom are hardcoded (26)
    log_info(wdiag(idx_diag) + ": OBS: +/- 95% conv.: " + obs0(idx_diag) + \
             " +/- " + (cdft_t(0.05 / 2.0, 26) * obs0_std(idx_diag)) + \
             " MODELS: " + (avg(diag_a(idx_diag, :))) + " +/- " + \
             (stddev(diag_a(idx_diag, :))))

    ; Only account for good diagnostics where no model is missing
    if (.not. any(ismissing(diag_a(idx_diag, :)))) then
      if (.not. isvar("good_diag_idx")) then
        good_diag_idx = idx_diag
      else
        tmp_gd = good_diag_idx
        good_diag_idx := new((/dimsizes(tmp_gd) + 1/), integer)
        good_diag_idx(:dimsizes(tmp_gd) - 1) = tmp_gd
        good_diag_idx(dimsizes(tmp_gd)) = idx_diag
        delete(tmp_gd)
      end if
    else
      missing_datasets := ind(ismissing(diag_a(idx_diag, :)))
      ; FIXME: use "f"
      error_msg("w", diag_script, "", \
                "MDER was not succesfull. In diagnostic " + wdiag(idx_diag) + \
                ", some datasets are missing: " + \
                diag_a&datasets(missing_datasets))
    end if
  end do
  log_info("")

  ; Save "good" diagnostics without missing datasets
  wdiag := wdiag(good_diag_idx)
  diag_all = diag_a(good_diag_idx, :)
  diag_all_std = diag_a_std(good_diag_idx, :)
  diag_all&diagnostics = wdiag
  diag_all_std&diagnostics = wdiag
  DIM_DIAG = dimsizes(good_diag_idx)
  DIM_DIAG_ALL = dimsizes(diag_a(:, 0))
  delete([/diag_a, diag_a_std, good_diag_idx/])

  ; ---------------------------------------------------------------------------
  ; Regression calculations
  ; ---------------------------------------------------------------------------

  ; Loop over different target periods
  do pt = 0, dimsizes(p_time) - 1
    log_info("Starting calculations for target period " + p_time(pt) + "-" + \
             (p_time(pt) + p_step))

    ; Skip calculation if values are missing
    if (any(ismissing(xx(:, pt)))) then
      continue
    end if

    ; Regression model: cross-validation cycle
    if (cross_validation_mode) then
      do i = 0, DIM_DAT - 1
        log_debug("Starting cross-validation cycle for dataset " + \
                  input_file_info@dataset(i))
        good = ind(input_file_info@dataset .ne. input_file_info@dataset(i))
        coef = regress_stepwise(diag_all(:, good), xx(good, pt))
        const = coef@const
        used_diags = toint(coef@idx_used_diags)
        n_terms = coef@n_terms
        YFIT = coef@YFIT
        mcorr = coef@MCORRELATION
        corr = coef@correlation
        log_debug("Cross-validation mode. Dataset: " + \
                  input_file_info@dataset(i) + ", constant: " + const)
        if (n_terms .gt. 0) then
          do it = 0, dimsizes(coef) - 1
            log_debug(wdiag(used_diags(it)) + ", " + coef(it) + ", " + \
                      corr(it))
          end do
        end if
        log_debug("Multiple correlation: " + mcorr)
        log_debug("")
        delete([/const, used_diags, n_terms, YFIT, mcorr, corr, coef, \
                good/])
      end do
    end if

    ; Building the final regression
    coef := regress_stepwise(diag_all, xx(:, pt))
    const = coef@const
    used_diags := toint(coef@idx_used_diags)
    n_terms := coef@n_terms
    YFIT := coef@YFIT
    mcorr := coef@MCORRELATION
    corr := coef@correlation

    ; Print final regression model
    log_info("The final model (" + scal_time(0) + "-" + \
             (p_time(pt) + p_step) + "):")
    log_info("    Constant: " + const)
    if (n_terms .gt. 0) then
      do it = 0, dimsizes(coef) - 1
        log_info("    " + wdiag(used_diags(it)) + ": " + coef(it) + \
                 ", (corr = " + corr(it) + ")")
      end do
    end if
    log_info("    " + "Multiple correlation: " + mcorr)
    log_info("")

    ; -------------------------------------------------------------------------
    ; Model weighting
    ; -------------------------------------------------------------------------

    x_diag := transpose(diag_all(used_diags, :))
    good_diag := ind(.not. ismissing(x_diag(:, 0)))
    x := x_diag(good_diag, :)

    ; tmp_dat = DIM_DAT
    tmp_dat := dimsizes(x(:, 0))

    ; Create Design Matrix D
    D := new((/tmp_dat, dimsizes(used_diags) + 1/), typeof(diag_all))
    D = 1
    do id = 0, dimsizes(x(:, 0)) - 1
      D(id, 1:) = x(id, :)
    end do

    ; Observations
    xobs := new((/1, dimsizes(used_diags) + 1/), typeof(diag_all), \
                getFillValue(diag_all))
    xobs = 1.0
    xobs(0, 1:) = obs0(used_diags)
    df = tmp_dat - dimsizes(xobs(0, :))

    ; Identity matrix
    Id := new((/tmp_dat, tmp_dat/), typeof(diag_all))
    Id = 0.0
    do idat = 0, tmp_dat - 1
      Id(idat, idat) = 1.0
    end do

    ; Auxiliary vectors
    m := new((/tmp_dat, 1/), float)    ; column vector
    m = 1.0
    n := new((/tmp_dat, 1/), float)    ; column vector
    n = (/1.0 / tmp_dat/)

    ; Future projections
    y := xx(good_diag, pt)

    ; Calculate model weights (see Karpechko et al. 2013, eq. (4))
    weights := (transpose(n) + (xobs(:, 1:) - transpose(n) # x) # \
                inverse_matrix(transpose(x) # x - transpose(x) # m # \
                               transpose(n) # x) # \
                (transpose(x) - transpose(x) # m # transpose(n)))

    ; Confidence interval
    cfint = cdft_t(0.05 / 2.0, df) ^ 2 * \
      (1.0 + xobs # inverse_matrix(transpose(D) # D) # transpose(xobs))
    SSE1 := y # (Id - D # inverse_matrix(transpose(D) # D) # transpose(D))
    sigmae := (1.0 / df) * SSE1
    tmp_cf := cfint(0, 0) * sigmae
    delete(cfint)

    ; Weighted mean and std
    y_weight := (weights # y)
    cf := sqrt(tmp_cf # y)
    uM = avg(y)
    uM_std = stddev(y) / sqrt(DIM_DAT)
    log_info("uMMM = " + uM + " +/- " + uM_std)
    log_info("y_weight = " + y_weight)
    log_info("cf = " + cf)
    log_info("")

    ; Save weights and confidence intervals in one array
    do idat = 0, DIM_DAT - 1
      do gdiad = 0, dimsizes(good_diag) - 1
        if (idat .eq. good_diag(gdiad))
          weights_t(:, idat, pt) = weights(:, gdiad)
          cfint_t(idat, pt) = tmp_cf(gdiad)
        end if
      end do
    end do

    ; Output
    log_info("Model weights:")
    do idat = 0, DIM_DAT - 1
      log_info(input_file_info@dataset(idat) + ": " + weights_t(:, idat, pt))
    end do
    log_info("")

    ; Save Y-estimates and multiple correllation coeficients
    coeff_y_est(:, pt) = coef@Yest
    MCORR(pt) = mcorr
    Y_WEIGHT(pt) = y_weight

    ; Save coefficients and terms of regression model in array
    ; Array sizes are first known at runtime (-> build at runtime)
    if (isvar("mder_coefs")) then
      max_dim = max((/dimsizes(coef), dimsizes(mder_coefs(pt - 1, 1:))/))
      tmp_mder_coefs = mder_coefs
      tmp_names = names_used_diags
      tmpt  = idx_used_diags
      delete([/mder_coefs, names_used_diags, idx_used_diags/])

      mder_coefs = new((/pt + 1, max_dim + 1/), typeof(tmp_mder_coefs))
      names_used_diags = new((/pt + 1, max_dim + 1/), typeof(tmp_names))
      idx_used_diags = new((/pt + 1, max_dim + 1/), typeof(tmpt))

      mder_coefs(:pt - 1, :dimsizes(tmp_mder_coefs(0, :)) - 1) = tmp_mder_coefs
      names_used_diags(:pt - 1, :dimsizes(tmp_names(0, :)) - 1) = tmp_names
      idx_used_diags(:pt-1, :dimsizes(tmpt(0, :)) - 1) = tmpt
      mder_coefs(pt, 0) = const
      mder_coefs(pt, 1:dimsizes(coef)) = coef
      names_used_diags(pt, 1:dimsizes(wdiag(used_diags))) = wdiag(used_diags)
      idx_used_diags(pt, 1:dimsizes(used_diags)) = used_diags

      delete([/tmp_mder_coefs, tmp_names, tmpt, max_dim/])
    else
      mder_coefs = new((/1, dimsizes(coef) + 1/), typeof(coef))
      names_used_diags = new((/1, dimsizes(coef) + 1/), string)
      idx_used_diags = new((/1, dimsizes(coef) + 1/), typeof(used_diags))

      mder_coefs(pt, 0) = const
      mder_coefs(pt, 1:) = coef
      names_used_diags(pt, 1:) = wdiag(used_diags)
      idx_used_diags(pt, 1:) = used_diags
    end if

    delete([/y, corr, mcorr, x, x_diag, good_diag, xobs, cf, tmp_cf, \
            weights, uM, uM_std/])
  end do    ; Loop over different target periods

  ; Iterate over target periods and save MDER results
  do pt = 0, dimsizes(p_time) - 1
    tmp_yy = new(dimsizes(YY), typeof(YY))
    tmp_mder_errors = new(dimsizes(mder_errors), typeof(mder_errors))
    do tt = 0, n_max - 2
      YY(0, tt) = avg(YY(3:, tt))

      ; MDER errors for time periods
      if (tt .ge. p_time(pt) - ny_min .and. \
          tt .lt. (p_time(pt) + p_step + 1 - ny_min)) then
        tmp_yy(1, tt) = sum(YY(3:, tt) * weights_t(0, :, pt))
        tmp_mder_errors(0, tt) = tmp_yy(1, tt) + \
          sqrt(abs(sum(YY(3:, tt) * cfint_t(:, pt))))
        tmp_mder_errors(1, tt) = tmp_yy(1, tt) - \
          sqrt(abs(sum(YY(3:, tt) * cfint_t(:, pt))))
      end if
      if (tt .eq. (p_time(pt) + p_step - ny_min)) then
        YY(1, tt - (p_step + 1) / 2)  = avg(tmp_yy(1, tt-p_step:tt))
        mder_errors(0, tt - (p_step + 1) / 2) = \
          avg(tmp_mder_errors(0, tt - p_step:tt))
        mder_errors(1, tt - (p_step + 1) / 2) = \
          avg(tmp_mder_errors(1, tt - p_step:tt))
      end if

      ; MDER errors for whole time span
      YY(2, tt) = sum(YY(3:, tt) * weights_t(0, :, pt))
      t_mder_errors(0, tt) = YY(2, tt) + \
        sqrt(abs(sum(YY(3:, tt) * cfint_t(:, pt))))
      t_mder_errors(1, tt) = YY(2, tt) - \
        sqrt(abs(sum(YY(3:, tt) * cfint_t(:, pt))))
    end do
  end do

  ; ---------------------------------------------------------------------------
  ; Plot preparations
  ; ---------------------------------------------------------------------------

  ; Plot regression line
  ;   y_mder_results(:, 0): const
  ;   y_mder_results(:, 1): fitted observations (with MDER coefficients)
  ;   y_mder_results(:, 2:): fitted models (with model weights)
  y_mder_results = new((/dimsizes(p_time), DIM_DAT + 2/), typeof(mder_coefs))
  y_mder_results_std = new((/dimsizes(p_time), DIM_DAT + 2/), \
                           typeof(mder_coefs))

  nfit = fspan(0, dimsizes(y_mder_results(0, 2:)) - 1, \
               dimsizes(y_mder_results(0, 2:)))
  ; nfit = fspan(0, 49, 50)
  wy = new((/dimsizes(p_time), DIM_DAT + 2/), typeof(mder_coefs))
  wy_std = new((/dimsizes(p_time), DIM_DAT + 2/), float)
  yfit_reg = new((/dimsizes(p_time), dimsizes(nfit)/), typeof(mder_coefs))
  yfit_std = new((/dimsizes(p_time), 2, dimsizes(nfit)/), typeof(mder_coefs))

  ; Iterate over target periods
  do ppt = 0, dimsizes(p_time) - 1
    xgood = ind(.not. ismissing(idx_used_diags(ppt, :)))
    hlp_tmp = new((/dimsizes(mder_coefs(ppt, xgood)), dimsizes(nfit)/), float)
    hlp = new((/1, dimsizes(mder_coefs(ppt, xgood)) + 1/), float)
    y_mder_results(ppt, :) = mder_coefs(ppt, 0)          ; const
    y_mder_results_std(ppt, :) = 0.0
    yfit_reg(ppt, :) = mder_coefs(ppt, 0)                 ; const
    hlp = 1.0

    wy(ppt, 0) = avg(xx(:, ppt))
    wy(ppt, 1) = Y_WEIGHT(ppt)
    wy(ppt, 2:) = xx(:, ppt)
    wy_std(ppt, 2:) = xx_val(:, ppt, 1)       ; cfint_t(:, ppt)

    do ii = 1, dimsizes(mder_coefs(ppt, :)) - 1
      if (.not. ismissing(mder_coefs(ppt, ii))) then
        y_mder_results(ppt, 1) = y_mder_results(ppt, 1) + \
          mder_coefs(ppt, ii) * obs0(idx_used_diags(ppt, ii))
        y_mder_results_std(ppt, 1) = y_mder_results_std(ppt, 1) + \
          (mder_coefs(ppt, ii) * (obs0_std(idx_used_diags(ppt, ii)))) ^ 2

        y_mder_results(ppt, 2:) = coeff_y_est(:, ppt)
        ; = y_mder_results(2:) + coef(ii) * diag_all(terms(ii), :)

        y_mder_results_std(ppt, 2:) = y_mder_results_std(ppt, 2:) + \
          (mder_coefs(ppt, ii) * \
           (diag_all_std(idx_used_diags(ppt, ii), :))) ^ 2

        minx = min(diag_all(idx_used_diags(ppt, ii), :)) - 1.1 * \
          (max(diag_all(idx_used_diags(ppt, ii), :)) - \
           min(diag_all(idx_used_diags(ppt, ii), :)))
        maxx = max(diag_all(idx_used_diags(ppt, ii), :)) + 1.1 * \
          (max(diag_all(idx_used_diags(ppt, ii), :)) - \
           min(diag_all(idx_used_diags(ppt, ii), :)))
        diag_var = minx + (maxx - minx) * nfit / dimsizes(nfit)
        hlp_tmp(ii - 1, :) = diag_var
        delete([/diag_var, minx, maxx/])
      end if
    end do

    y_mder_results_std = sqrt(y_mder_results_std)
    x_diag = transpose(diag_all(idx_used_diags(ppt, xgood), :))
    good_diag = ind(.not. ismissing(x_diag(:, 0)))
    x = x_diag(good_diag, :)

    ; Create Design matrix D and Identity matrix Id for pt = 1
    D := new((/dimsizes(x(:, 0)), dimsizes(idx_used_diags(ppt, xgood)) + 1/), \
             typeof(diag_all))
    D = 1.0
    do id = 0, dimsizes(x(:, 0)) - 1
      D(id, 1:) = x(id, :)
    end do
    df = tmp_dat - dimsizes(obs0(idx_used_diags(ppt, xgood))) + 1
    Id = new((/dimsizes(x(:, 0)), dimsizes(x(:, 0))/), typeof(diag_all))
    Id = 0.0
    do idat = 0, tmp_dat - 1
      Id(idat, idat) = 1.0
    end do

    ; Error
    SSE = xx(:, ppt) # (Id - D # inverse_matrix(transpose(D) # D) # \
                        transpose(D)) # transpose(xx(:, ppt))
    sigmae2 = sqrt((1.0 / df) * SSE)

    ; Calculate 95% confidence interval bounds
    if (.not. loglin) then
      miny = min(xx(:, ppt)) - 1.1 * (max(xx(:, ppt)) - min(xx(:, ppt)))
      maxy = max(xx(:, ppt)) + 1.1 * (max(xx(:, ppt)) - min(xx(:, ppt)))
      yfit_reg(ppt, :) = miny + (maxy - miny) * nfit / dimsizes(nfit)
      delete([/miny, maxy/])
      do ti = 0, dimsizes(nfit) - 1
        hlp(0, 1:) = (/hlp_tmp(:, ti)/)
        yfit_std(ppt, 0, ti) = yfit_reg(ppt, ti) + cdft_t(0.05 / 2.0, df) * \
          sigmae2 * sqrt(1.0 + hlp # inverse_matrix(transpose(D) # D) # \
                         transpose(hlp))
        yfit_std(ppt, 1, ti) = yfit_reg(ppt, ti) - cdft_t(0.05 /2.0, df) * \
          sigmae2 * sqrt(1.0 + hlp # inverse_matrix(transpose(D) # D) # \
                         transpose(hlp))
      end do

    ; Log-linear regression
    else
      isort_tmp = dim_pqsort_n(y_mder_results(ppt, 2:), 1, 0)
      isort = new(dimsizes(isort_tmp) + 2, integer)
      isort(0) = 0
      isort(1) = 1
      isort(2:) = isort_tmp + 2
      y_mder_results(ppt, :) = y_mder_results(ppt, isort)
      y_mder_results_std(ppt, :) = y_mder_results_std(ppt, isort)
      wy(ppt, :) = wy(ppt, isort)
      wy_std(ppt, :) = wy_std(ppt, isort)
      model_numbers = ispan(0, dimsizes(y_mder_results(ppt, :)) - 1, 1)

      ; Save order of datasets
      y_mder_results@model_numbers = model_numbers(isort)
      rc = regline_stats(y_mder_results(ppt, 2:), log(wy(ppt, 2:)))
      yfit_reg(ppt, :) = exp(rc@b(0) + rc@b(1) * y_mder_results(ppt, 2:))
      yfit_std(ppt, :, :) = \
        (/exp(rc@y95(0) + rc@b95(0) * y_mder_results(ppt, 2:)), \
         exp(rc@y95(1) + rc@b95(1) * y_mder_results(ppt, 2:))/)
      wy(ppt, 1) = exp(rc@b(0) + rc@b(1) * y_mder_results(ppt, 1))  ; yfit(obs)
      delete([/isort_tmp, isort, model_numbers, rc/])
    end if

    ; Clean up
    delete([/xgood, hlp_tmp, hlp, x_diag, good_diag, x, D, df, SSE, sigmae2/])
  end do

  ; Plot squared differences between ensemble mean and the total change in
  ; pseudoreality
  ; 1. #grades, 2. #pseudorealities
  grade_all = new((/DIM_DIAG_ALL, DIM_DAT - 1, DIM_DAT/), float)
  xx_pre = new((/DIM_DAT, dimsizes(p_time)/), float)
  grade_trans = dim_avg_n_Wrap(grade_all(:, :, :), 0)
  do im = 0, DIM_DAT - 1
    good = ind(input_file_info@dataset(im) .ne. input_file_info@dataset)
    do pt = 0, dimsizes(p_time) - 1
      xx_pre(im, pt) = avg(xx(good, pt))
    end do
    delete(good)
  end do

  ; Calculating statistics
  x25 = round(0.25 * DIM_DAT, 3) - 1
  x75 = round(0.75 * DIM_DAT, 3) - 1

  presq = new((/dimsizes(p_time), DIM_DAT/), typeof(mder_coefs))
  premdersq = new((/dimsizes(p_time), DIM_DAT/), typeof(mder_coefs))
  BSS = new((/dimsizes(p_time)/), float)

  ; Iterate over target periods
  do pt = 0, dimsizes(p_time) - 1
    presq(pt, :) = (xx(:, pt) - xx_pre(:, pt)) ^ 2
    premdersq(pt, :) = (xx(:, pt) - coeff_y_est(:, pt)) ^ 2
    qsort(presq(pt, :))
    qsort(premdersq(pt, :))

    tmpRDS = (/(/presq(pt, x25), presq(pt, x25), dim_median(presq(pt, :)), \
                presq(pt, x75), presq(pt, x75)/), \
              (/premdersq(pt, x25), premdersq(pt, x25), \
               dim_median(premdersq(pt, :)), premdersq(pt, x75), \
               premdersq(pt, x75)/)/)
    tmpAVG = (/avg(presq(pt, :)), avg(premdersq(pt, :))/)
    tmpPRE = (/presq(pt, :), premdersq(pt, :)/)

    if (all(isvar((/"preRDS", "preAVG", "prePRE"/)))) then
      preRDS := array_append_record(preRDS, sqrt(tmpRDS), 0)
      preAVG := array_append_record(preAVG, sqrt(tmpAVG), 0)
      prePRE := array_append_record(prePRE, sqrt(tmpPRE), 0)
    else
      preRDS = sqrt(tmpRDS)
      preAVG = sqrt(tmpAVG)
      prePRE = sqrt(tmpPRE)
    end if
    delete([/tmpRDS, tmpAVG, tmpPRE/])
    log_info("RMSE (" + p_time(pt) + "-" + (p_time(pt) + p_step) + \
             "), uMMM: " + sqrt(avg(presq(pt, :))) + ", MDER: " + \
             sqrt(avg(premdersq(pt, :))))
  end do

  ; Add uMMM and MDER to input_file_info@dataset
  tmpm = input_file_info@dataset
  input_file_info@dataset := new(DIM_DAT + 3, string)
  input_file_info@dataset(0) = "uMMM"
  input_file_info@dataset(1) = "MDER time dependent"
  input_file_info@dataset(2) = "MDER"
  input_file_info@dataset(3:) = tmpm
  delete(tmpm)
  DIM_DAT = DIM_DAT + 2

  ; ---------------------------------------------------------------------------
  ; Write regression model equations to ASCII table
  ; ---------------------------------------------------------------------------

  log_info("")

  ; Create arrays
  times_str = new((/dimsizes(p_time) + 1/), string)
  terms_str = new((/dimsizes(p_time) + 1/), string)
  mcorr_str = new((/dimsizes(p_time) + 1/), string)

  ; Header
  times_str(0) = "Target period"
  terms_str(0) = "Model equation"
  mcorr_str(0) = "Correlation R^2"

  ; Contents
  do pt = 0, dimsizes(p_time) - 1
    terms_str(pt + 1) = sprintf("%4.2f", mder_coefs(pt, 0))
    do ii = 1, dimsizes(mder_coefs(pt, :)) - 1
      if (.not. ismissing(mder_coefs(pt, ii))) then
        if (mder_coefs(pt, ii) .gt. 0.0) then
          terms_str(pt + 1) = terms_str(pt + 1) + " +" + \
            sprintf("%4.2f", mder_coefs(pt, ii)) + " x " + \
            names_used_diags(pt, ii)
        else
          terms_str(pt + 1) = terms_str(pt + 1) + " " + \
            sprintf("%4.2f", mder_coefs(pt, ii)) + " x " + \
            names_used_diags(pt, ii)
        end if
      end if
    end do
    p1_time = p_time(pt) + (p_step - 2000)
    times_str(pt + 1) = p_time(pt) + "-" + sprinti("%0.2hi", p1_time)
    mcorr_str(pt + 1) = sprintf("%.3f", MCORR(pt))
  end do

  ; Write files
  table_path = work_dir + filename + ".txt"
  ascii_path = work_dir + var0 + "_" + \
    input_file_info@exp(0 + DIM_OBS + DIM_HDAT) + "_" + (p_step + 1) + \
    "ystep_diagnostics.txt"
  write_table(table_path, "w", [/times_str, terms_str, mcorr_str/], \
              "%s    %s    %s")
  asciiwrite(ascii_path, names_used_diags)
  log_info("Wrote " + ascii_path)
  log_info("Wrote " + table_path)
  delete([/table_path, ascii_path/])

  ; ---------------------------------------------------------------------------
  ; Write NetCDF files
  ; ---------------------------------------------------------------------------

  dataset_list = input_file_info@dataset(3:)
  do pt = 0, dimsizes(p_time) - 1
    weights_out = weights_t(0, :, pt)
    weights_path = work_dir + filename + "_weights_" + \
      p_time(pt) + "-" + (p_time(pt) + p_step) + ".nc"
    weights_out!0 = "dataset"
    weights_out&dataset = dataset_list
    weights_out@diag_script = diag_script
    weights_out@var = "weight"
    weights_out@description = "Weights from MDER"
    if (target_var .eq. "SIE") then
      weights_out@description = "Weights for September Arctic sea ice extent"
    end if
    weights_file = ncdf_write(weights_out, weights_path)
    delete([/weights_out, weights_path, weights_file/])
  end do

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Line plot showing the time series including MDER weights
  plot_name_1 = filename + "_FIG1"
  wks = get_wks(0, diag_script, plot_name_1)
  MainStg = ""
  YStg = target_var + " [" + variable_info@units + "]"
  plot = weight_lin_plot(YY, mder_errors, t_mder_errors, yy_obs, MainStg, \
                         YStg, wks, "draw")
  log_info("Wrote " + plot_dir + plot_name_1 + "." + file_type)

  ; TODO: Add meta data
  ; climofiles = new(DIM_DAT, string)
  ; do ii = 0, DIM_DAT - 1
  ;   climofiles(ii) = interface_get_inpaths(ii) + "/" + \
  ;                    interface_get_infile(var0, field_type0, ii)
  ; end do
  ; alldatasetnames = input_file_info@dataset
  ; zz = diag_script + "_" + target_var
  ; contrib_authors = (/"A_wenz_sa, A_schl_ma"/)
  ; alltags = array_append_record(tags, (/"ST_mder", "ST_corr", "ST_range", \
  ;                                      "ST_stddev", "PT_times", "DM_sh"/), 0)
  ; caption = "Time series of " + target_var + " for " + \
  ;   input_file_info@ensemble(0) + \
  ;   " scenario, similar to Wenzel et al., JClim, 2016, fig. 5"

  ; plot_fname = wks@fullname
  ; ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
  ;          alldatasetnames, climofiles, diag_script, contrib_authors)
  ; delete([/alltags, caption, wks, plot_fname/])

  ; MDER scatter plots for all time periods
  do ppt = 0, dimsizes(p_time) - 1
    plot_name_2 = filename + "_FIG2" + plot_num(ppt)
    wks = get_wks(0, diag_script, plot_name_2)
    MainStg = plot_num(ppt) + ") " + p_time(ppt) + "-" + \
      (p_time(ppt) + p_step)
    XStg = terms_str(ppt + 1) + " [" + variable_info@units + "]"
    YStg = target_var + " [" + variable_info@units + "]"
    plot2 = mder_scatter_plot(y_mder_results(ppt, :), \
                              y_mder_results_std(ppt, :), \
                              wy(ppt, :), wy_std(ppt, :), \
                              yfit_reg(ppt, :), \
                              yfit_std(ppt, :, :), \
                              MainStg, XStg, YStg, wks, "draw")
    log_info("Wrote " + plot_dir + plot_name_2 + "." + file_type)

    ; TODO: Add meta data
    ; alltags = array_append_record(tags, (/"ST_mder", "ST_corr", \
    ;                                      "ST_range", "ST_stddev", \
    ;                                      "ST_corr", "PT_scatter", \
    ;                                      "T_EC", "DM_sh"/), 0)
    ; caption = "Scatterplot showing the correlation between future " + \
    ;           target_var + " and the quantity resulting from MDER " + \
    ;           "analysis on the x-axes, similar to Wenzel et al., " + \
    ;           "JClim, 2016, fig. 3"
    ; plot_fname = wks@fullname
    ; ESMValMD(plot_fname, alltags, caption, zz, var0, \
    ;          alldatasetnames, climofiles, diag_script, contrib_authors)
    ; delete([/alltags, caption, wks, plot_fname/])
  end do

  ; Boxplot of RMSE of different weigting methods inkl pseudoreality check
  plot_name_3 = filename + "_FIG3"
  YStg = "RMSE projection error [" + variable_info@units + "]"
  XBname = new((dimsizes(times_str) - 1) * 2, string)
  if (dimsizes(XBname) .eq. 2) then
    XBname = (/"uMMM", "MDER"/)
  else if (dimsizes(XBname) .eq. 4) then
    XBname(0::2) = (/"uMMM~C~" + times_str(1:)/)
    XBname(1::2) = (/"MDER~C~" + times_str(1:)/)
  else
    XBname(0::2) = (/" "/)
    XBname(1::2) = times_str(1:)
  end if
  end if

  wks = get_wks(0, diag_script, plot_name_3)
  MainStg = ""
  plot3 = squared_error_plot(preRDS, preAVG, prePRE, MainStg, XBname, YStg, \
                             wks, "draw")
  log_info("Wrote " + plot_dir + plot_name_3 + "." + file_type)

  ; TODO: Add meta data
  ; alltags = array_append_record(tags, (/"ST_mder", "ST_rmsd", \
  ;                                      "ST_stddev", "ST_mean", "PT_bar", \
  ;                                      "PT_errorbar", "DM_sh"/), 0)
  ; caption = "RMSE difference between the ensemble mean future " + \
  ;   "climatological mean " + target_var + "and its equivalent in " + \
  ;   "pseudoreality under " + input_file_info@ensemble(0) + " scenario, " + \
  ;   "similar to Wenzel et al., JClim, 2016, fig. 4"
  ; plot_fname = wks@fullname
  ; ESMValMD(plot_fname, alltags, caption, zz, var0, \
  ;          alldatasetnames, climofiles, diag_script, contrib_authors)
  ; delete([/alltags, caption, zz, climofiles, alldatasetnames, plot_fname/])

  leave_msg(diag_script, "")

end
