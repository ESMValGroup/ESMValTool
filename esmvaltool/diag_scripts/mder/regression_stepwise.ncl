; #############################################################################
; mder/regression_stepwise.ncl
; #############################################################################
;
; Description:
;     This routine computes weights for a set of climate models following the
;     multiple diagnostic ensemble regression (MDER) described by Karpechko et
;     al. 2013. The script calls a variable of interest from different models
;     (which shall be weighted) and additionally several diagnostics which are
;     specified in the recipe for computig the weights.
;     Calls three plotscripts:
;     + Timeseries of the variable of interest including the multi-model
;       mean, the weighted and unweighted model ensemble.
;     + Scatter plot of the climatological mean of the variable of interest
;       and the result of the MDER.
;     + Boxplot of RMSE in a pseudo-reality out-of-sample testing.
;
; Required diag_script_info attributes:
;     p_time: start years for future projections.
;     p_step: time range for future projections.
;     scal_time: Time range for base period used in anomaly calculations when
;                "calc_type = 'trend'".
;     time_oper: select for operation in time_operation function.
;     time_opt: select for option in time_operation function.
;     calc_type: trend/pos.
;     domain: domain of the diagnostic (used for provenance tracking).
;
; Optional diag_script_info attributes:
;     average_ens: average over ensemble members of the same model (default:
;                  False).
;     smooth: smoothes time period with 1-2-1 filter.
;     iter: number of iteration for smoothing (required when "smooth" is set).
;     cross_validation_mode: perform cross validation.
;     region: the region to be averaged (required for multidimensional input).
;     area_oper: operation in area_operation function (required for
;                multidimensional input).

; Optional variable_info attributes:
;     plot_units: convert variable to these units.
;
; Caveats:
;     "warning: in unique_labels_min (diag_scripts/shared/plot/style.ncl), Add
;     more attributes to prio to make labels unique! Continuing with non-unique
;     labels" is normal when using "average_ens = true".
;
; Modification history
;     20191121-schlund_manuel: added new provenance tracking.
;     20180806-schlund_manuel: ported to v2.0.
;     201303??-wenzel_sabrina: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/latlon.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/mder.ncl"
load "$diag_scripts/shared/plot/style.ncl"

load "$diag_scripts/shared/mder.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  ; Variable
  VAR0 = variable_info[0]
  var0 = VAR0@short_name
  DIM_VAR = ListCount(variable_info)
  if (DIM_VAR .gt. 1) then
    error_msg("w", DIAG_SCRIPT, "", "this diagnostic supports only one " + \
              "variable, processing " + VAR0@short_name)
  end if

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  ALL_DATASETS = metadata_att_as_array(INFO0, "dataset")

  ; Rank of input data
  dummy_data = read_data(INFO0[0])
  rank = dimsizes(dimsizes(dummy_data))
  if (rank .gt. 1) then
    is_multi_dim = True
  else
    is_multi_dim = False
  end if
  delete([/dummy_data, rank/])

  ; OBS
  OBS = get_obs_list(INFO0)
  if (ListCount(OBS) .lt. 1) then
    error_msg("f", DIAG_SCRIPT, "", "this diagnostic needs at least one " + \
              "OBS dataset")
  end if
  OBS_DATASETS = metadata_att_as_array(OBS, "dataset")
  DIM_OBS = dimsizes(OBS_DATASETS)

  ; CMIP5 models
  atts := True
  atts@project = "CMIP5"
  CMIP5 = select_metadata_by_atts(INFO0, atts)
  exps = metadata_att_as_array(CMIP5, "exp")
  future_exp = get_unique_values(exps(ind(exps .ne. "historical")))
  if (dimsizes(future_exp) .ne. 1) then
    error_msg("f", DIAG_SCRIPT, "", "this diagnostic supports only one " + \
              "projection experiment, got " + future_exp)
  end if

  ; Historical experiment
  atts@exp = "historical"
  HIST = select_metadata_by_atts(CMIP5, atts)
  hist_datasets = metadata_att_as_array(HIST, "dataset")

  ; Future experiment
  atts@exp = future_exp
  FUT = select_metadata_by_atts(CMIP5, atts)
  fut_datasets = metadata_att_as_array(FUT, "dataset")
  fail = False
  if (dimsizes(hist_datasets) .ne. dimsizes(fut_datasets)) then
    fail = True
  else
    if (all(hist_datasets .ne. fut_datasets)) then
      fail = True
    end if
  end if
  if (fail) then
    error_msg("f", DIAG_SCRIPT, "", "historical and future experiment do " + \
              "not cover the same models, got " + dimsizes(hist_datasets) + \
              " historical models and " + dimsizes(fut_datasets) + \
              " projection models")
  else
    MODELS = fut_datasets
  end if

  ; Ensemble averaging (if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens) then
    MODELS := get_unique_values(MODELS)
  else
    MODELS := unique_labels_min(FUT, (/"dataset", "ensemble"/))
  end if
  DIM_MOD = dimsizes(MODELS)

  ; Directories
  work_dir = config_user_info@work_dir
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

end

begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Read years
  min_year = 9999
  max_year = -9999
  do idat = 0, ListCount(CMIP5) - 1
    if (CMIP5[idat]@start_year .lt. min_year) then
      min_year = INFO0[idat]@start_year
    end if
    if (CMIP5[idat]@end_year .gt. max_year) then
      max_year = INFO0[idat]@end_year
    end if
  end do
  time_span = max_year - min_year + 1

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Required attributes
  req_atts = (/"p_time", "p_step", "scal_time", "time_oper", "time_opt", \
              "calc_type", "domain"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  p_time = diag_script_info@p_time
  p_step = diag_script_info@p_step
  scal_time = diag_script_info@scal_time
  time_oper = diag_script_info@time_oper
  time_opt = diag_script_info@time_opt
  calc_type = diag_script_info@calc_type
  dim_p_time = dimsizes(p_time)

  ; Provenance
  AUTHORS = (/"wenzel_sabrina", "schlund_manuel"/)
  DOMAIN = diag_script_info@domain
  STATISTICS = (/"mder"/)
  REFERENCES = (/"wenzel16jclim"/)

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
  smooth = False
  if (isatt(diag_script_info, "smooth")) then
    smooth = diag_script_info@smooth
    if (isatt(diag_script_info, "iter")) then
      iter = diag_script_info@iter
    else
      error_msg("f", diag_script, "", "this diagnotic requires attribute " + \
                "'iter' in recipe when 'smooth' is set to True")
    end if
  end if
  cross_validation_mode = False
  if (isatt(diag_script_info, "cross_validation_mode")) then
    cross_validation_mode = diag_script_info@cross_validation_mode
  end if

  ; Plot numbering
  plot_num = integertochar(ispan(97, 122, 1))
  filename = "CMPI5_" + var0 + \
    "-" + calc_type + "_" + future_exp + "_" + (p_step + 1) + "ystep"

  ; Output arrays
  xx = new((/DIM_MOD, dimsizes(p_time)/), float)
  xx_val = new((/DIM_MOD, dimsizes(p_time), 2/), float)
  YY = new((/DIM_MOD + 3, time_span/), float)         ; incl. uMMM and wMMM
  t_mder_errors = new((/2, time_span/), float)
  mder_errors = new((/2, time_span/), float)
  yy_obs = new((/DIM_OBS, time_span/), float)
  weights_t = new((/1, DIM_MOD, dimsizes(p_time)/), float)
  cfint_t = new((/DIM_MOD, dimsizes(p_time)/), float)
  coeff_y_est = new((/DIM_MOD, dimsizes(p_time)/), double)  ; Type of coef
  Y_WEIGHT = new((/dimsizes(p_time)/), double)              ; Type of coef
  MCORR = new((/dimsizes(p_time)/), double)                 ; Type of coef

  ; ---------------------------------------------------------------------------
  ; Read and preprocess model data
  ; ---------------------------------------------------------------------------

  ; Call the variable of interest and save in one array YY
  ANCESTORS = new(DIM_MOD * 2, string)
  do imod = 0, DIM_MOD - 1
    log_info(MODELS(imod))

    ; Average over ensemble members if desired
    if (avgens) then
      atts := True
      atts@dataset = MODELS(imod)
      info_hist := select_metadata_by_atts(HIST, atts)
      info_fut := select_metadata_by_atts(FUT, atts)
    else
      info_hist := NewList("lifo")
      info_fut := NewList("lifo")
      ListAppend(info_hist, HIST[imod])
      ListAppend(info_fut, FUT[imod])
    end if

    ; Get data
    data_hist = get_average(info_hist)
    data_fut = get_average(info_fut)
    info_hist := info_hist[0]
    info_fut := info_fut[0]

    ; Get filenames
    ANCESTORS(2 * imod) = info_hist@filename
    ANCESTORS(2 * imod + 1) = info_fut@filename

    ; Process multi-dimensional variables
    if (is_multi_dim) then
      data_hist := area_operations(data_hist, region(0), region(1), \
                                   region(2), region(3), area_oper, True)
      data_fut := area_operations(data_fut, region(0), region(1), \
                                  region(2), region(3), area_oper, True)
    end if

    ; Append future array to historical array
    data_y = array_append_record(data_hist, data_fut, 0)
    delete([/data_hist, data_fut/])

    ; Convert units if desired
    if (isatt(VAR0, "plot_units")) then
      data_y := convert_units(data_y, VAR0@plot_units)
    end if

    ; Adapt scal_time(1) (necessary of larger than end_year)
    if (scal_time(1) .gt. info_hist@end_year) then
      scal_time(1) = info_hist@end_year
    end if

    ; Calculate average over starting period for trend calculation
    trend_base = time_operations(data_y, scal_time(0), scal_time(1), \
                                 time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      trend_base := time_operations(trend_base, scal_time(0), scal_time(1), \
                                    "average", "yearly", True)
    end if
    trend_base = tofloat(avg(trend_base))

    ; Loop over different target periods
    do pt = 0, dim_p_time - 1
      ; Average over var time series
      y_mean = time_operations(data_y, p_time(pt), p_time(pt) + p_step, \
                               time_oper, time_opt, True)
      if (time_oper .eq. "extract") then
        y_mean := time_operations(y_mean, p_time(pt), \
                                  p_time(pt) + p_step, "average", \
                                  "yearly", True)
      end if

      ; Perform calculation type
      if (calc_type .eq. "trend") then
        xx(imod, pt) = tofloat(avg(y_mean) - trend_base)
        xx_val(imod, pt, 0) = tofloat(avg(y_mean) - trend_base)

        ; FIXME: Why substracting base value from error?
        xx_val(imod, pt, 1) = tofloat(stddev(y_mean) / \
                                      sqrt(dimsizes(y_mean)) - trend_base)
      else if (calc_type .eq. "pos") then
        xx(imod, pt) = tofloat(avg(y_mean))
        xx_val(imod, pt, 0) = tofloat(avg(y_mean))
        xx_val(imod, pt, 1) = tofloat(stddev(y_mean) / \
                                      sqrt(dimsizes(y_mean)))
      else if (calc_type .eq. "int") then
        integrated_vals = y_mean
        do tt = 0, dimsizes(y_mean) - 2
          integrated_vals(tt + 1) = integrated_vals(tt) + y_mean(tt + 1)
        end do
        xx(imod, pt) = tofloat(avg(integrated_vals))
        xx_val(imod, pt, 0) = tofloat(avg(integrated_vals))
        xx_val(imod, pt, 1) = tofloat(stddev(integrated_vals) / \
                                      sqrt(dimsizes(integrated_vals)))
        delete([/integrated_vals, tt/])
      else
        error_msg("f", diag_script, "", "calc_type " + calc_type + \
                  " not supported")
      end if
      end if
      end if
      delete(y_mean)
    end do

    ; Extract season for var time series
    y_mean = time_operations(data_y, -1, -1, time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      y_mean := time_operations(y_mean, -1, -1, "average", "yearly", \
                                True)
    end if

    ; Smoothing var time series with 1-2-1 filter
    if (smooth) then
      y_mean := filter121(y_mean, iter)
    end if
    if (calc_type .eq. "trend") then
      YY(imod + 3, 0:dimsizes(y_mean) - 1) = tofloat(y_mean - trend_base)
    else if (calc_type .eq. "pos") then
      YY(imod + 3, 0:dimsizes(y_mean) - 1) = tofloat(y_mean)
    else if (calc_type .eq. "int") then
      YY(imod + 3, 0) = y_mean(0)
      do tt = 0, dimsizes(y_mean) - 2
        YY(imod + 3, tt + 1) = YY(i + 3, tt) + y_mean(tt + 1)
      end do
    end if
    end if
    end if
    delete([/trend_base, data_y, y_mean/])
  end do

  ; Consider ensemble averaging
  if (avgens) then
    INFO0 := get_unique_items(INFO0, "dataset")
    CMIP5 := get_unique_items(CMIP5, "dataset")
    HIST := get_unique_items(HIST, "dataset")
    FUT := get_unique_items(FUT, "dataset")
  end if

  ; Add attributes to models time series
  YY!1 = "time"
  YY&time = ispan(min_year, max_year, 1)      ; y_mean&year
  YY&time@calendar = "standard"

  ; ---------------------------------------------------------------------------
  ; Read and preprocess observational data
  ; ---------------------------------------------------------------------------

  ; Call OBS of the variable of interest and save in one array yy_obs
  do iobs = 0, DIM_OBS - 1
    log_info(OBS_DATASETS(iobs))
    info = OBS[iobs]
    ANCESTORS := array_append_record(ANCESTORS, info@filename, 0)

    ; Read observational data
    data_obs = read_data(info)

    ; Process multi-dimensional variables
    if (is_multi_dim) then
      data_obs := area_operations(data_obs, region(0), region(1), region(2), \
                                  region(3), area_oper, True)
    end if

    ; Convert units if desired
    if (isatt(VAR0, "plot_units")) then
      data_obs := convert_units(data_obs, VAR0@plot_units)
      VAR0@units = VAR0@plot_units
    end if

    ; Perform calculation type
    ymin = min(info@start_year) - min_year
    obs_mean = time_operations(data_obs, -1, -1, time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      obs_mean := time_operations(obs_mean, -1, -1, "average", "yearly", True)
    end if
    if (calc_type .eq. "trend") then
      VAR0@units = VAR0@units + "/dec"
      obs_base = time_operations(data_obs, scal_time(0), -1, "average", \
                                 "yearly", True)
      obs_base = tofloat(avg(obs_base))
      yy_obs(iobs, ymin:ymin + dimsizes(obs_mean) - 1) = tofloat(obs_mean - \
                                                                 obs_base)
      delete(obs_base)
    else if (calc_type .eq. "pos") then
      yy_obs(iobs, ymin:ymin + dimsizes(obs_mean) - 1) = tofloat(obs_mean)
    else if (calc_type .eq. "int") then
      yy_obs(iobs, ymin) = obs_mean(0)
      do tt = 0, dimsizes(obs_mean) - 2
        yy_obs(iobs, ymin + tt + 1) = yy_obs(iobs, ymin + tt) + \
          obs_mean(tt + 1)
      end do
      VAR0@units = integrate_time_units(VAR0@units, 0)
      delete(tt)
    end if
    end if
    end if
    delete([/data_obs, ymin, obs_mean/])
  end do

  ; Add attributes to observations time series
  yy_obs!0 = "name"
  yy_obs!1 = "time"
  yy_obs&name := OBS_DATASETS
  yy_obs&time = YY&time
  yy_obs&time@calendar = "standard"

  ; ---------------------------------------------------------------------------
  ; Call diagnostics
  ; ---------------------------------------------------------------------------

  ; Get input directory (using ancestor diagnostics)
  input_files := diag_script_info@input_files
  var_files = str_match(input_files, "_info.ncl")
  do ivar = 0, dimsizes(var_files) - 1
    input_files := input_files(ind(input_files .ne. var_files(ivar)))
  end do
  input_dirs = input_files
  if (dimsizes(input_dirs) .lt. 1) then
    error_msg("f", DIAG_SCRIPT, "", "no diagnositcs given (use ancestors " + \
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
  ANCESTORS := array_append_record(ANCESTORS, file_list, 0)

  ; Create diagnostic arrays
  DIM_DIAG = dimsizes(file_list)
  diag_a = new((/DIM_DIAG, DIM_MOD/), float)
  diag_a_std = new((/DIM_DIAG, DIM_MOD/), float)
  obs0 = new((/DIM_DIAG/), float)
  obs0_std = new((/DIM_DIAG/), float)
  wdiag = new((/DIM_DIAG/), string)

  log_info("")
  ; Loop over files in the list and append data
  do idiag = 0, DIM_DIAG - 1
    log_info("Reading " + file_list(idiag))
    data_temp = ncdf_read(file_list(idiag), "diag")

    ; Loop over datasets and save diagnostic data
    do imod = 0, DIM_MOD - 1
      do idat = 0, dimsizes(data_temp&datasets) - 1

        ; Add data if datasets match
        if (MODELS(imod) .eq. data_temp&datasets(idat)) then
          diag_a(idiag, imod) = data_temp(0, idat, 0)
          diag_a_std(idiag, imod) = data_temp(0, idat, 1)
          if (isdim(data_temp, "diagnostics")) then
            wdiag(idiag) = tostring(data_temp&diagnostics(0))
          else
            error_msg("f", DIAG_SCRIPT, "", "Cannot read diagnostic name " + \
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
        obs0(idiag) = data_temp(0, idx_ref, 0)
        obs0_std(idiag) = data_temp(0, idx_ref, 1)
        delete(idx_ref)
      end if
    end if
    delete(data_temp)

    ; FIXME: degrees of freedom are hardcoded (26)?
    log_info(wdiag(idiag) + ": OBS: +/- 95% conv.: " + obs0(idiag) + \
             " +/- " + (cdft_t(0.05 / 2.0, 26) * obs0_std(idiag)) + \
             " MODELS: " + (avg(diag_a(idiag, :))) + " +/- " + \
             (stddev(diag_a(idiag, :))))

    ; Only account for good diagnostics where no model is missing
    if (.not. any(ismissing(diag_a(idiag, :)))) then
      if (.not. isvar("good_diag_idx")) then
        good_diag_idx = idiag
      else
        tmp_gd = good_diag_idx
        good_diag_idx := new((/dimsizes(tmp_gd) + 1/), integer)
        good_diag_idx(:dimsizes(tmp_gd) - 1) = tmp_gd
        good_diag_idx(dimsizes(tmp_gd)) = idiag
        delete(tmp_gd)
      end if
    else
      missing_datasets := ind(ismissing(diag_a(idiag, :)))
      error_msg("w", DIAG_SCRIPT, "", "in diagnostic " + wdiag(idiag) + \
                ", some datasets are missing: " + \
                diag_a&datasets(missing_datasets))
    end if
  end do
  if (.not. isvar("good_diag_idx")) then
    error_msg("f", DIAG_SCRIPT, "", "no diagnostic without missing " + \
              "datasets given")
  end if

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
      do imod = 0, DIM_MOD - 1
        model = MODELS(imod)
        log_debug("Starting cross-validation cycle for dataset " + model)
        good = ind(MODELS .ne. model)
        coef = regress_stepwise(diag_all(:, good), xx(good, pt))
        const = coef@const
        used_diags = toint(coef@idx_used_diags)
        n_terms = coef@n_terms
        YFIT = coef@YFIT
        mcorr = coef@MCORRELATION
        corr = coef@correlation
        log_debug("Cross-validation mode. Dataset: " + model + \
                  ", constant: " + const)
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

    ; tmp_dat = DIM_MOD
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
    cfint1 = cdft_t(0.05 / 2.0, df) ^ 2 * \
      (1.0 + xobs # inverse_matrix(transpose(D) # D) # transpose(xobs))
    SSE1 := y # (Id - D # inverse_matrix(transpose(D) # D) # transpose(D))
    sigmae1 := (1.0 / df) * SSE1
    tmp_cf := cfint1(0, 0) * sigmae1
    delete([/SSE1, sigmae1, cfint1/])

    ; Weighted mean and std
    y_weight := (weights # y)
    cf := sqrt(tmp_cf # y)
    uM = avg(y)
    uM_std = stddev(y) / sqrt(DIM_MOD)
    log_info("uMMM = " + uM + " +/- " + uM_std)
    log_info("y_weight = " + y_weight)
    log_info("cf = " + cf)
    log_info("")

    ; Save weights and confidence intervals in one array
    do imod = 0, DIM_MOD - 1
      do gdiad = 0, dimsizes(good_diag) - 1
        if (imod .eq. good_diag(gdiad))
          weights_t(:, imod, pt) = weights(:, gdiad)
          cfint_t(imod, pt) = tmp_cf(gdiad)
        end if
      end do
    end do

    ; Output
    log_info("Model weights:")
    do imod = 0, DIM_MOD - 1
      log_info(MODELS(imod) + ": " + weights_t(:, imod, pt))
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
      tmpt = idx_used_diags
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
  end do        ; Loop over different target periods

  ; Iterate over target periods and save MDER results
  do pt = 0, dimsizes(p_time) - 1
    tmp_yy = new(dimsizes(YY), typeof(YY))
    tmp_mder_errors = new(dimsizes(mder_errors), typeof(mder_errors))
    do tt = 0, time_span - 2
      YY(0, tt) = avg(YY(3:, tt))

      ; MDER errors for time periods
      if (tt .ge. p_time(pt) - min_year .and. \
          tt .lt. (p_time(pt) + p_step + 1 - min_year)) then
        tmp_yy(1, tt) = sum(YY(3:, tt) * weights_t(0, :, pt))
        tmp_mder_errors(0, tt) = tmp_yy(1, tt) + \
          sqrt(abs(sum(YY(3:, tt) * cfint_t(:, pt))))
        tmp_mder_errors(1, tt) = tmp_yy(1, tt) - \
          sqrt(abs(sum(YY(3:, tt) * cfint_t(:, pt))))
      end if
      if (tt .eq. (p_time(pt) + p_step - min_year)) then
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
  y_mder_results = new((/dimsizes(p_time), DIM_MOD + 2/), typeof(mder_coefs))
  y_mder_results_std = new((/dimsizes(p_time), DIM_MOD + 2/), \
                           typeof(mder_coefs))

  nfit = fspan(0, dimsizes(y_mder_results(0, 2:)) - 1, \
               dimsizes(y_mder_results(0, 2:)))
  ; nfit = fspan(0, 49, 50)
  wy = new((/dimsizes(p_time), DIM_MOD + 2/), typeof(mder_coefs))
  wy_std = new((/dimsizes(p_time), DIM_MOD + 2/), float)
  yfit_reg = new((/dimsizes(p_time), dimsizes(nfit)/), typeof(mder_coefs))
  yfit_std = new((/dimsizes(p_time), 2, dimsizes(nfit)/), typeof(mder_coefs))

  ; Iterate over target periods
  do ppt = 0, dimsizes(p_time) - 1
    xgood = ind(.not. ismissing(idx_used_diags(ppt, :)))
    hlp_tmp = new((/dimsizes(mder_coefs(ppt, xgood)), dimsizes(nfit)/), float)
    hlp = new((/1, dimsizes(mder_coefs(ppt, xgood)) + 1/), float)
    y_mder_results(ppt, :) = mder_coefs(ppt, 0)           ; const
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

    ; Create Design matrix D and Identity matrix Id
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
    SSE2 = xx(:, ppt) # (Id - D # inverse_matrix(transpose(D) # D) # \
                         transpose(D))
    sigmae2 = (1.0 / df) * SSE2
    delete(SSE2)

    ; Calculate 95% confidence interval bounds
    miny = min(xx(:, ppt)) - 1.1 * (max(xx(:, ppt)) - min(xx(:, ppt)))
    maxy = max(xx(:, ppt)) + 1.1 * (max(xx(:, ppt)) - min(xx(:, ppt)))
    yfit_reg(ppt, :) = miny + (maxy - miny) * nfit / dimsizes(nfit)
    delete([/miny, maxy/])
    do ti = 0, dimsizes(nfit) - 1
      hlp(0, 1:) = (/hlp_tmp(:, ti)/)
      cfint2 = cdft_t(0.05 / 2.0, df) ^ 2 * ( \
        1.0 + hlp # inverse_matrix(transpose(D) # D) # transpose(hlp))
      cfint2 := cfint2(0, 0) * sigmae2
      cfint2 := sqrt(cfint2 # xx(:, ppt))

      yfit_std(ppt, 0, ti) = yfit_reg(ppt, ti) + cfint2
      yfit_std(ppt, 1, ti) = yfit_reg(ppt, ti) - cfint2

      delete(cfint2)
    end do

    ; Clean up
    delete([/xgood, hlp_tmp, hlp, x_diag, good_diag, x, D, df, sigmae2/])
  end do

  ; Plot squared differences between ensemble mean and the total change in
  ; pseudoreality
  ; 1. #grades, 2. #pseudorealities
  grade_all = new((/DIM_DIAG_ALL, DIM_MOD - 1, DIM_MOD/), float)
  xx_pre = new((/DIM_MOD, dimsizes(p_time)/), float)
  grade_trans = dim_avg_n_Wrap(grade_all(:, :, :), 0)
  do imod = 0, DIM_MOD - 1
    good = ind(MODELS(imod) .ne. MODELS)
    do pt = 0, dimsizes(p_time) - 1
      xx_pre(imod, pt) = avg(xx(good, pt))
    end do
    delete(good)
  end do

  ; Calculating statistics
  x25 = round(0.25 * DIM_MOD, 3) - 1
  x75 = round(0.75 * DIM_MOD, 3) - 1

  presq = new((/dimsizes(p_time), DIM_MOD/), typeof(mder_coefs))
  premdersq = new((/dimsizes(p_time), DIM_MOD/), typeof(mder_coefs))
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

  ; ---------------------------------------------------------------------------
  ; Write files
  ; ---------------------------------------------------------------------------

  if (diag_script_info@write_netcdf) then

    ; Ascii
    table_path = work_dir + filename + ".txt"
    ascii_path = work_dir + var0 + "_" + future_exp + "_" + \
      (p_step + 1) + "ystep_diagnostics.txt"
    write_table(table_path, "w", [/times_str, terms_str, mcorr_str/], \
                "%s    %s    %s")
    asciiwrite(ascii_path, names_used_diags)
    log_info("Wrote " + ascii_path)
    log_info("Wrote " + table_path)
    delete([/table_path, ascii_path/])

    ; NetCDF
    do pt = 0, dimsizes(p_time) - 1
      weights_out = weights_t(0, :, pt)
      weights_path = work_dir + filename + "_weights_" + \
        p_time(pt) + "-" + (p_time(pt) + p_step) + ".nc"
      weights_out!0 = "dataset"
      weights_out&dataset = MODELS
      weights_out@diag_script = DIAG_SCRIPT
      weights_out@var = "weight"
      weights_out@description = "Weights from MDER"
      weights_file = ncdf_write(weights_out, weights_path)

      ; Provenance tracking
      plot_path = "missing"
      plot_type = ""
      caption = "MDER model weights for target variable " + var0
      log_provenance(weights_file, plot_path, caption, STATISTICS, \
                     DOMAIN, plot_type, AUTHORS, REFERENCES, ANCESTORS)
      delete([/weights_out, weights_path, weights_file/])
    end do
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  if (diag_script_info@write_plots) then

    ; Add uMMM and MDER to array
    datasets := new(DIM_MOD + 3, string)
    datasets(0) = "uMMM"
    datasets(1) = "MDER time dependent"
    datasets(2) = "MDER"
    datasets(3:) = MODELS
    YY!0 = "dataset"
    YY&dataset := datasets
    NEW_INFO = NewList("lifo")
    do i = 0, 2
      ListAppend(NEW_INFO, True)
      NEW_INFO[i]@project = "CMIP5"
      NEW_INFO[i]@dataset = datasets(i)
    end do
    do imod = 0, DIM_MOD - 1
      ListAppend(NEW_INFO, FUT[imod])
    end do
    do iobs = 0, DIM_OBS - 1
      ListAppend(NEW_INFO, OBS[iobs])
    end do
    delete(datasets)

    ; Line plot showing the time series including MDER weights
    plot_name_1 = filename + "_FIG1"
    path_1 = plot_dir + plot_name_1 + "." + file_type
    wks = gsn_open_wks(file_type, plot_dir + plot_name_1)
    MainStg = future_exp
    YStg = var0 + " [" + VAR0@units + "]"
    plot = weight_lin_plot(YY, mder_errors, t_mder_errors, yy_obs, MainStg, \
                           YStg, wks, "draw", NEW_INFO)
    log_info("Wrote " + path_1)

    ; Write ncdf file
    new_path = work_dir + plot_name_1 + ".nc"
    YY@var = var0
    YY@diag_script = DIAG_SCRIPT
    YY@ncdf = new_path
    outfile = ncdf_write(YY, new_path)

    ; Provenance tracking
    caption = "Time series of " + var0 + " for " + future_exp + "."
    plot_type = (/"times"/)
    log_provenance(outfile, path_1, caption, STATISTICS, DOMAIN, plot_type, \
                   AUTHORS, REFERENCES, ANCESTORS)

    ; MDER scatter plots for all time periods
    y_mder_results!0 = "time"
    y_mder_results!1 = "dataset"
    y_mder_results&time = p_time
    data_dim = (/"constant", "OBS"/)
    data_dim := array_append_record(data_dim, MODELS, 0)
    y_mder_results&dataset = data_dim
    delete(data_dim)
    do ppt = 0, dimsizes(p_time) - 1
      plot_name_2 = filename + "_FIG2" + plot_num(ppt)
      path_2 = plot_dir + plot_name_2 + "." + file_type
      wks = gsn_open_wks(file_type, plot_dir + plot_name_2)
      MainStg = plot_num(ppt) + ") " + p_time(ppt) + "-" + \
        (p_time(ppt) + p_step) + " (" + future_exp + ")"
      XStg = terms_str(ppt + 1) + " [" + VAR0@units + "]"
      YStg = var0 + " [" + VAR0@units + "]"
      plot2 = mder_scatter_plot(y_mder_results(ppt, :), \
                                y_mder_results_std(ppt, :), \
                                wy(ppt, :), wy_std(ppt, :), \
                                yfit_reg(ppt, :), \
                                yfit_std(ppt, :, :), \
                                MainStg, XStg, YStg, wks, "draw", INFO0)
      log_info("Wrote " + path_2)

      ; Write ncdf file
      new_path = work_dir + plot_name_2 + ".nc"
      y_mder_results@var = var0
      y_mder_results@diag_script = DIAG_SCRIPT
      y_mder_results@ncdf = new_path
      outfile = ncdf_write(y_mder_results, new_path)

      ; Provenance tracking
      caption = "Scatterplot showing the correlation between " + \
                "future " + var0 + " and the quantity resulting from MDER " + \
                "analysis on the X axis for " + future_exp + "."
      plot_type = (/"scatter"/)
      log_provenance(outfile, path_2, caption, STATISTICS, DOMAIN, plot_type, \
                     AUTHORS, REFERENCES, ANCESTORS)
    end do

    ; Boxplot of RMSE of different weigting methods incl pseudoreality check
    plot_name_3 = filename + "_FIG3"
    path_3 = plot_dir + plot_name_3 + "." + file_type
    YStg = "RMSE projection error [" + VAR0@units + "]"
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

    wks = gsn_open_wks(file_type, plot_dir + plot_name_3)
    MainStg = future_exp
    plot3 = squared_error_plot(preRDS, preAVG, prePRE, MainStg, XBname, YStg, \
                               wks, "draw")
    log_info("Wrote " + path_3)

    ; Write ncdf file
    new_path = work_dir + plot_name_3 + ".nc"
    prePRE!0 = "type"
    prePRE!1 = "pseudo_reality_obs"
    do pt = 0, dimsizes(p_time) - 1
      if (isvar("dim_type")) then
        dim_type := array_append_record(dim_type, (/"uMMM " + p_time(pt), \
                                                   "wMMM " + p_time(pt)/), 0)
      else
        dim_type := (/"uMMM " + p_time(pt), "wMMM " + p_time(pt)/)
      end if
    end do
    prePRE&type = dim_type
    dim_PRE = dimsizes(prePRE)
    prePRE&pseudo_reality_obs = ispan(0, dim_PRE(1) - 1, 1)
    prePRE@var = var0
    prePRE@diag_script = DIAG_SCRIPT
    prePRE@ncdf = new_path
    outfile = ncdf_write(prePRE, new_path)

    ; Provenance tracking
    caption = "RMSE difference between the ensemble mean future " + \
      "climatological mean " + var0 + " and its equivalent in " + \
      "pseudo-reality under " + future_exp + "."
    plot_type = (/"box"/)
    log_provenance(outfile, path_3, caption, STATISTICS, DOMAIN, plot_type, \
                   AUTHORS, REFERENCES, ANCESTORS)
  end if

  leave_msg(DIAG_SCRIPT, "")

end
