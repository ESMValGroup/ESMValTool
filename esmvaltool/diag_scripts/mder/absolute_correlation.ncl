; #############################################################################
; mder/absolute_correlation.ncl
; #############################################################################
;
; Description:
;     This routine calculates absolute correlation coefficients between all
;     diagnostics, calculated for a present-day period.
;
;     Calls one plotscript:
;     + Plots the absolute correlation coefficient between the variable of
;       interest and present-day diagnostic.
;
; Required diag_script_info attributes:
;     p_time: start years for future projections.
;     p_step: time range for future projections.
;     scal_time: Time range for base period used in anomaly calculations when
;                "calc_type = 'trend'".
;     time_oper: operation in time_operation function.
;     time_opt: option in time_operation function.
;     calc_type: trend/pos/int.
;     domain: domain of the diagnostic (used for provenance tracking).
;
; Optional diag_script_info attributes:
;     average_ens: average over ensemble members of the same model (default:
;                  False).
;     region: the region to be averaged (required for multidimensional input).
;     area_oper: operation in area_operation function (required for
;                multidimensional input).
;
; Optional variable_info attributes:
;     plot_units: convert variable to these units.
;
; Caveats:
;     "warning: in unique_labels_min (diag_scripts/shared/plot/style.ncl), Add
;     more attributes to prio to make labels unique! Continuing with non-unique
;     labels" is normal when using "average_ens = true".
;
; Modification history:
;     20191121-schlund_manuel: added new provenance tracking.
;     20180726-schlund_manuel: ported to v2.0.
;     201307??-wenzel_sabrina: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/latlon.ncl"

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
  PLOT_TYPE = "errorbar"
  STATISTICS = (/"corr"/)
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
      error_msg("f", DIAG_SCRIPT, "", diag_script_info@region + \
                " is an invalid region type")
    end if
    end if
  end if

  ; Plot numbering
  plot_num = integertochar(ispan(97, 122, 1))

  ; Output arrays
  xx = new((/dim_p_time, DIM_MOD/), float)
  future_diag = new((/DIM_MOD/), float)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Call the variable of interest and merge to one timeline
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
      VAR0@units = VAR0@plot_units
    end if

    ; Calculate future diagnostic for diagnostic scatter plot
    future_temp = time_operations(data_y, p_time(0), \
                                  (/p_time(dim_p_time - 1) + p_step/), \
                                  time_oper, time_opt, True)
    if (time_oper .eq. "extract") then
      future_diag(imod) = tofloat(avg(time_operations( \
        future_temp, p_time(0), (/p_time(dim_p_time - 1) + p_step/), \
        "average", "yearly", True)))
    else
      future_diag(imod) = tofloat(avg(future_temp))
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
        xx(pt, imod) = tofloat(avg(y_mean)) - trend_base
      else if (calc_type .eq. "pos") then
        xx(pt, imod) = tofloat(avg(y_mean))
      else if (calc_type .eq. "int") then
        integrated_vals = y_mean
        do tt = 0, dimsizes(y_mean) - 2
          integrated_vals(tt + 1) = integrated_vals(tt) + y_mean(tt + 1)
        end do
        xx(pt, imod) = tofloat(avg(integrated_vals))
        delete([/integrated_vals, tt/])
      else
        error_msg("f", DIAG_SCRIPT, "", "calc_type " + calc_type + \
                  " not supported")
      end if
      end if
      end if
      delete(y_mean)
    end do
    delete([/data_y, trend_base/])
  end do

  ; Consider ensemble averaging
  if (avgens) then
    INFO0 := get_unique_items(INFO0, "dataset")
    HIST := get_unique_items(HIST, "dataset")
    FUT := get_unique_items(FUT, "dataset")
  end if

  ; Adapt units for "int" calculation
  if (calc_type .eq. "int") then
    VAR0@units = integrate_time_units(VAR0@units, 0)
  end if

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
  diag_all = new((/DIM_DIAG, DIM_MOD/), float)
  wdiag = new((/DIM_DIAG/), string)

  ; Loop over files in the list and append data
  do idiag = 0, DIM_DIAG - 1
    log_info("Reading " + file_list(idiag))
    data_temp = ncdf_read(file_list(idiag), "diag")

    ; Loop over datasets and save diagnostic data
    do imod = 0, DIM_MOD - 1
      do idat = 0, dimsizes(data_temp&datasets) - 1

        ; Add data if datasets match
        if (MODELS(imod) .eq. data_temp&datasets(idat)) then
          diag_all(idiag, imod) = data_temp(0, idat, 0)
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
    delete(data_temp)
  end do
  diag_all!0 = "diag"
  diag_all!1 = "dataset"
  diag_all&diag = wdiag
  diag_all&dataset = MODELS

  ; Create scatterplots for each diagnostic
  res_scat = True
  future_diag@future_exp = future_exp
  var0@units = VAR0@units
  scatterplots = diag_scatter(res_scat, diag_all, future_diag, var0, FUT)

  ; ---------------------------------------------------------------------------
  ; Correlation calculations
  ; ---------------------------------------------------------------------------

  ; New arrays
  corr = new((/dim_p_time, DIM_DIAG/), float)
  diag_sign = new((/dim_p_time, DIM_DIAG/), string)

  ; Loop over different target periods
  do pt = 0, dim_p_time - 1
    corr_sig = new((/2, DIM_DIAG/), float)
    r_sig0  = new((/DIM_DIAG/), float)
    r_sig1  = new((/DIM_DIAG/), float)
    r_wdiag  = new((/DIM_DIAG/), string)

    ; Iterate over diagnostics
    do idiag = 0, DIM_DIAG - 1
      hlp = diag_all(idiag, :)

      ; Only use models with data (filter missing values)
      good = ind(hlp .ne. 999.0)
      ngood = dimsizes(good)

      ; Calculate correlations
      corr(pt, idiag) = escorc(hlp(good), xx(pt, good))
      rc = regline(hlp(good), xx(pt, good))
      df = rc@nptxy - 2
      z = 0.5 * log((1.0 + abs(corr(pt, idiag))) / \
                    (1.0 - abs(corr(pt, idiag))))
      corr_sig(0, idiag) = tanh(z - 1.96 / (sqrt(ngood - 3.0)))
      corr_sig(1, idiag) = tanh(z + 1.96 / (sqrt(ngood - 3.0)))

      t0 = abs(corr(pt, idiag)) / sqrt((1 - corr(pt, idiag) ^ 2) / \
                                       (ngood-2))
      sign0 = 1.0 - betainc(df / (df + rc@tval ^ 2), df / 2.0, 0.5)

      log_info(wdiag(idiag) + " corr: " + corr(pt, idiag) + \
               "  datasets: " + ngood + "  stddv: " + sign0(0))
      delete([/hlp, good, ngood, rc, df/])
    end do

    ; Sort abs(coerr. coef.) to decreasing order in array
    acorr = abs(corr(pt, :))
    qsort(acorr)
    acorr = acorr(::-1)
    do idiag = 0, DIM_DIAG - 1
      idx_diag = ind(abs(corr(pt, :)) .eq. acorr(idiag))
      r_wdiag(idiag) = wdiag(idx_diag)
      r_sig0(idiag) = corr_sig(0, idx_diag)
      r_sig1(idiag) = corr_sig(1, idx_diag)
      diag_sign(pt, idiag) = where(corr(pt, idx_diag) .gt. 0.0, "red", "blue")
      delete(idx_diag)
    end do
    backup_wdiag = wdiag
    delete([/wdiag, corr_sig/])
    wdiag = r_wdiag
    corr_sig = new((/2, DIM_DIAG/), float)
    corr_sig(0, :) = r_sig0
    corr_sig(1, :) = r_sig1
    delete([/r_wdiag, r_sig0, r_sig1/])
    log_info(wdiag + ": " + acorr)
    acorr!0 = "diag"
    acorr&diag = wdiag

    ; -------------------------------------------------------------------------
    ; Plot
    ; -------------------------------------------------------------------------

    if (config_user_info@write_plots) then
      filename  = var0 + "-" + calc_type + "_" + future_exp + "_" + \
        p_time(pt) + "-" + (/p_time(pt) + p_step/)
      wks = gsn_open_wks(file_type, plot_dir + filename)
      MainStg = plot_num(pt) + ") " + p_time(pt) + "-" + \
        (/p_time(pt) + p_step/) + " (" + future_exp + ")"
      YStg = "Absolute Correlation Coefficient"
      plot = error_bar_plot(acorr, corr_sig, wdiag, diag_sign(pt, :), \
                            MainStg, YStg, wks, "draw")

      ; Write ncdf file
      new_path = work_dir + filename + ".nc"
      acorr@var = var0
      acorr@diag_script = DIAG_SCRIPT
      acorr@ncdf = new_path
      outfile_acorr = ncdf_write(acorr, new_path)

      ; Provenance tracking
      plot_path = plot_dir + filename + "." + file_type
      caption = "Absolute values of the correlation coefficient between " + \
        var0 + " and multiple present-day diagnostics."
      log_provenance(outfile_acorr, plot_path, caption, STATISTICS, \
                     DOMAIN, PLOT_TYPE, AUTHORS, REFERENCES, ANCESTORS)
    end if

    ; Restore old order
    wdiag = backup_wdiag
    delete([/corr_sig, acorr, backup_wdiag/])
  end do

  leave_msg(DIAG_SCRIPT, "")

end
