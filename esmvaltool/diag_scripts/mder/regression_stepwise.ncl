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
; Required recipe attributes:
;     p_time : time range for projections.
;     p_step : averaging time range for target period.
;     scal_time : a year to scale the timeseries with.
;     time_oper : select for operation in time_operation function.
;     time_opt : select for option in time_operation function.
;     calc_type : trend/pos
;     wregr : exit criterion for stepwise regression (True/False)
;     wexit : exit criterion for stepwise regression (True/False)
;
; Optional recipe attributes:
;     units : units of target variable.
;     loglin : perform log-linear regression (default: False).
;     smooth : smoothes time period with 1-2-1 filter.
;     iter : iteration for smoothing.
;     cross_validation_mode : perform cross validation.
;     p_value : p value for F-test in stepwise regression (default: p=0.05)
;
; Caveats:
;     Adds dataset_info@units.
;     Modifies dataset_info@dataset and dataset_info@ensemble.
;     Modifies variable_info@units.
;
; Modification history
;   * 20180806_A-schl_ma: ported to v2.0
;   * 201808??_A-senf_da: added option to plot log-linear regression
;   * 201303??_A-wenz_sa: written
;
; #############################################################################

load "interface_scripts/interface.ncl"

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

  ; Check field type: need 1-d field
  ; supported_vars = (/"sic"/)  ; List of non-1-d variables that are supported
  supported_vars = (/""/)
  if (.not. isStrSubset(field_type0, "0")) then
    if (.not. any(var0 .eq. supported_vars)) then
    error_msg("f", diag_script, "", "Variable " + var0 + " not supported yet.")
    end if
  end if

  ; Extract observations from model-list
  IDX_OBS = str_match_ind_ic(dataset_info@project, "OBS")   ; Case insensitive
  DIM_OBS = dimsizes(IDX_OBS)
  OBS_NAME = dataset_info@dataset(IDX_OBS)

  ; Extract historical simulations from datasets
  IDX_HDAT_tmp = ind(input_file_info@start_year .ne. \
                     max(tofloat(input_file_info@start_year)))
  IDX_HDAT = IDX_HDAT_tmp(ind(IDX_HDAT_tmp .ne. IDX_OBS))
  DIM_HDAT = dimsizes(IDX_HDAT)
  delete(IDX_HDAT_tmp)

  ; Extract individual future projections from datasets (= unique models)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens) then
    idx_datasets = UNIQ(dataset_info@dataset(DIM_OBS + DIM_HDAT:)) + \
                   DIM_OBS + DIM_HDAT
  else
    idx_datasets = ispan(DIM_OBS + DIM_HDAT, \
                         dimsizes(dataset_info@dataset) - 1, 1)
    idx_hist = UNIQ(dataset_info@dataset(DIM_OBS:DIM_OBS + DIM_HDAT))
    if (dimsizes(idx_hist) .gt. 1) then
      all_datasets = unique_labels_min((/"dataset", "ensemble"/))
      dataset_info@dataset := all_datasets
      delete(all_datasets)
    end if
    delete(idx_hist)
  end if

  ; Get dataset variables
  DIM_DAT = dimsizes(dataset_info@dataset(idx_datasets))
  dataset_info@ensemble := dataset_info@ensemble(idx_datasets)

end

begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Read latitude range and level from input files
  f = addfile(input_file_info@filename(0), "r")
  if (isfilevar(f, "plev") .and. isfilevar(f, "lat")) then
    plev = read_data(0, "plev", field_type0)
    lat = read_data(0, "lat", field_type0)
    if (dimsizes(plev) .eq. 1) then
      lev = plev(0)
    else
      error_msg("f", diag_script, "", "ambiguous pressure levels given, " + \
                "only one supported")
    end if
    if (dimsizes(lat) .gt. 1) then
      lat_range = (/min(lat), max(lat)/)
    else
      lat_range = read_data(0, "lat_bnds", field_type0)
    end if
  else
    if (isatt(f, "plev") .and. \
        isatt(f, "lat_range_0") .and. \
        isatt(f, "lat_range_1")) then
      lev = f@plev
      lat_range = (/f@lat_range_0, f@lat_range_1/)
    else
      error_msg("f", diag_script, "", "cannot read plev and latrange")
    end if
  end if
  log_info("Retrieved pressure level " + lev + " and latitude range " + \
           lat_range(0) + " to " + lat_range(1) + " from preprocessed data")

  ; Read years
  ny_max = max(toint(input_file_info@end_year))
  ny_min = min(toint(input_file_info@start_year(DIM_OBS:)))
  n_max = ny_max - ny_min + 1

  ; Required attributes
  if (isatt(diag_script_info, "p_time")) then
    p_time = diag_script_info@p_time
  else
    error_msg("f", diag_script, "", "requires attribute p_time in recipe")
  end if
  if (isatt(diag_script_info, "p_step"))then
    p_step = diag_script_info@p_step
  else
    error_msg("f", diag_script, "", "requires attribute p_step in recipe")
  end if
  if (isatt(diag_script_info, "scal_time")) then
    scal_time = diag_script_info@scal_time
  else
    error_msg("f", diag_script, "", "requires attribute scal_time in recipe")
  end if
  if (isatt(diag_script_info, "time_oper"))then
    time_oper = diag_script_info@time_oper
  else
    error_msg("f", diag_script, "", "requires attribute time_oper in recipe")
  end if
  if (isatt(diag_script_info, "time_opt"))then
    time_opt = diag_script_info@time_opt
  else
    error_msg("f", diag_script, "", "requires attribute time_opt in recipe")
  end if
  if (isatt(diag_script_info, "calc_type"))then
    calc_type = diag_script_info@calc_type
  else
    error_msg("f", diag_script, "", "requires attribute calc_type in recipe")
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

  ; Output dir
  diag_script_base = basename(diag_script)
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  plot_num = integertochar(ispan(97, 122, 1))
  plot_file = dataset_info@project(DIM_OBS) + "_" + target_var + \
    "-" + calc_type + "_" + dataset_info@exp(DIM_OBS + DIM_HDAT) + \
    "_" + (p_step + 1) + "ystep"

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Output arrays
  xx = new((/DIM_DAT, dimsizes(p_time)/), float)
  xx_val = new((/DIM_DAT, dimsizes(p_time), 2/), float)
  yy = new((/DIM_DAT + 3, n_max/), float)         ; incl. uMMM and wMMM
  yyw = new((/2, n_max/), float)
  yywt = new((/2, n_max/), float)
  yy_obs = new((/DIM_OBS, n_max/), float)
  weights_t = new((/1, DIM_DAT, dimsizes(p_time)/), float)
  cfint_t = new((/DIM_DAT, dimsizes(p_time)/), float)
  coeff_y_est = new((/DIM_DAT, dimsizes(p_time)/), double)
  Y_WEIGHT = new((/dimsizes(p_time)/), double)
  MCORR = new((/dimsizes(p_time)/), double)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess model data
  ; ---------------------------------------------------------------------------

  ; Call the variable of interest and save in one array yy
  do i = 0, DIM_DAT - 1
    idat = idx_datasets(i)
    log_debug(dataset_info@dataset(idat))

    ; Call projection simulations for regression analysis and time series,
    ; scale time series and average over 10-year period
    tmp = read_data(idat, var0, field_type0)
    idx_hist = ind(dataset_info@dataset(idat) .eq. \
                   dataset_info@dataset(:DIM_OBS + DIM_HDAT - 1))
    idx_fut = ind(dataset_info@dataset(idat) .eq. \
                  dataset_info@dataset(DIM_OBS + DIM_HDAT:)) + \
      DIM_OBS + DIM_HDAT

    ; Calculate derived var (i.e. make it 1-d)
    ; TODO: Sea-ice extent
    ; if (var0 .eq. "sic") then  ; Calculate sea ice extent
    ;   if (isatt(diag_script_info, "region")) then
    ;     region = diag_script_info@region
    ;   else
    ;     error_msg("w", diag_script, "", "No region for sea ice extent " + \
    ;               "selected in cfg file, selecting Arctic.")
    ;     region = "Arctic"
    ;   end if
    ; tmp := sea_ice_area(tmp, i, data, region, "all", "extent", True)
    ; end if

    if (dimsizes(idx_hist) .ge. 1) then
      ; Call historical simulations for time series
      tmp_hist = read_data(idx_hist(0), var0, field_type0)

      ; TODO: Sea-ice extent
      ; if (var0 .eq. "sic") then
      ;   tmp_hist := sea_ice_area(tmp_hist, idx_hist(0), data, region, \
      ;                            "all", "extent", True)
      ; end if

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

              ; TODO: calculate sea ice extent
              ; if (var0 .eq. "sic") then
              ;   tmp_fut := sea_ice_area(tmp_fut, idx_fut(iens), data, \
              ;                           region, "all", "extent", True)
              ;   tmp_hist := sea_ice_area(tmp_hist, idx_hist(iens), data, \
              ;                            region, "all", "extent", True)
              ; end if

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
                      "size for dataset " + dataset_info@dataset(idat))
          end if
        end if
      end if

      ; Append future array to historical array
      tmpy = array_append_record(tmp_hist, tmp, 0)
      delete([/tmp, tmp_hist, idx_hist, idx_fut/])
    end if

    ; Calculate average over starting period for trend calculation
    tmpxx = time_operations(tmpy, scal_time(0), scal_time(1), "extract", \
                            time_opt, True)
    tmpyy = time_operations(tmpxx, scal_time(0), scal_time(1), "average", \
                            "yearly", True)
    tmppd = tofloat(avg(tmpyy))

    ; Loop over different target periods
    do pt = 0, dimsizes(p_time) - 1
      ; Average over var time series
      tmps = time_operations(tmpy, p_time(pt), p_time(pt) + p_step, \
                             time_oper, time_opt, True)
      tmpj = time_operations(tmps, p_time(pt), p_time(pt) + p_step, \
                             "average", "yearly", True)

      ; Perform calculation type
      if (calc_type .eq. "trend") then
        xx(i, pt) = tofloat(avg(tmpj) - tmppd)
        xx_val(i, pt, 0) = tofloat(avg(tmpj) - tmppd)
        xx_val(i, pt, 1) = tofloat(stddev(tmpj) / sqrt(dimsizes(tmpj)) - tmppd)
        variable_info@units = diag_script_info@units + "/dec"
      else if (calc_type .eq. "pos") then
        xx(i, pt) = tofloat(avg(tmpj))
        xx_val(i, pt, 0) = tofloat(avg(tmpj))
        xx_val(i, pt, 1) = tofloat(stddev(tmpj) / sqrt(dimsizes(tmpj)))
        variable_info@units = diag_script_info@units
      else
        error_msg("f", diag_script, "", "calc_type " + calc_type + \
                  " not supported")
      end if
      end if
      delete([/tmps, tmpj/])
    end do

    ; Extract season for var time series
    tmpst = time_operations(tmpy, -1, -1, time_oper, time_opt, True)
    tmpsc = time_operations(tmpst, -1, -1, "average", "yearly", True)

    ; Smoothing var time series with 1-2-1 filter
    if (smooth) then
      tmpsm = filter121(tmpsc, iter)
    else
      tmpsm = tmpsc
    end if
    if (calc_type .eq. "trend") then
      yy(i + 3, 0:dimsizes(tmpsm) - 1) = tofloat(tmpsm - tmppd)
    else if (calc_type .eq. "pos") then
      yy(i + 3, 0:dimsizes(tmpsm) - 1) = tofloat(tmpsm)
    end if
    end if
    delete([/tmppd, tmpyy, tmpxx, tmpy, tmpsm, tmpsc, tmpst/])
  end do

  ; Add attributes to models time series
  yy!1 = "time"
  yy&time = ispan(ny_min, ny_max, 1)      ; tmpsc&year
  yy&time@calendar = "standard"

  ; Save unique dataset names (modifies dataset_info@dataset!)
  tmpm = dataset_info@dataset
  dataset_info@dataset := tmpm(idx_datasets)
  delete(tmpm)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess observational data
  ; ---------------------------------------------------------------------------

  ; Call OBS of the variable of interest and save in one array yy_obs
  do iobs = 0, DIM_OBS - 1
    log_debug(OBS_NAME(iobs))

    ; Call projection simulations for regression analysis and time series,
    ; scale time series and average over 10-year period
    tmpy = read_data(iobs, var0, field_type0)

    ; Calculate derived var (i.e. make it 1-d)
    ; TODO: Sea-ice extent
    ; if (var0 .eq. "sic") then  ; Calculate sea ice extent
    ;   tmpy = sea_ice_area(tmpy, iobs, data, region, "all", "extent", True)
    ; end if

    ; Scale var time series
    tmpoy = time_operations(tmpy, scal_time(0), -1, "average", "yearly", True)
    tmppd = tofloat(avg(tmpoy))
    tmpyy = time_operations(tmpy, -1, -1, time_oper, time_opt, True)
    tmpsm = time_operations(tmpyy, -1, -1, "average", "yearly", True)
    ymin = min(toint(input_file_info@start_year(iobs))) - ny_min

    ; Perform calculation type
    if (calc_type .eq. "trend") then
      yy_obs(iobs, ymin:ymin + dimsizes(tmpsm) - 1) = tofloat(tmpsm - tmppd)
    else if (calc_type .eq. "pos") then
      yy_obs(iobs, ymin:ymin + dimsizes(tmpsm) - 1) = tofloat(tmpsm)
    else
      error_msg("f", diag_script, "", "calc_type " + calc_type + \
                " not supported")
    end if
    end if
    ; yy_obs(iobs, :dimsizes(tmpsm) - 1) = tofloat(tmpsm)
    delete([/tmpy, tmpyy, tmpsm, ymin, tmppd, tmpoy/])
  end do

  ; Add attributes to observations time series
  yy_obs!1 = "time"
  yy_obs&time = yy&time
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

  ; Loop over files in the list and append data
  do idx_diag = 0, DIM_DIAG - 1
    data_temp = ncdf_read(file_list(idx_diag), "diag")

    ; Loop over datasets and save diagnostic data
    do i = 0, DIM_DAT - 1
      do idat = 0, dimsizes(data_temp&datasets) - 1
        if (dataset_info@dataset(i) .eq. data_temp&datasets(idat)) then
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
             " DATASETS: " + (avg(diag_a(idx_diag, :))) + " +/- " + \
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
    if (cross_validation_mode .eq. True) then
      do i = 0, DIM_DAT - 1
        log_debug("Starting cross-validation cycle for dataset " + \
                  dataset_info@dataset(i))
        good = ind(dataset_info@dataset .ne. dataset_info@dataset(i))
        coef = regress_stepwise(diag_all(:, good), xx(good, pt))
        const = coef@const
        terms = toint(coef@terms)
        nterms = coef@nterms
        YFIT = coef@YFIT
        mcorr = coef@MCORRELATION
        corr = coef@correlation
        log_debug("Cross-validation mode. Dataset: " + \
                  dataset_info@dataset(i) + ", constant: " + const)
        if (nterms .gt. 0) then
          do it = 0, dimsizes(coef) - 1
            log_debug(wdiag(terms(it)) + ", " + coef(it) + ", " + corr(it))
          end do
        end if
        log_debug("Multiple correlation: " + mcorr)
        log_debug("")
        delete([/const, terms, nterms, YFIT, mcorr, corr, coef, good/])
      end do
    end if

    ; Building the final regression
    coef := regress_stepwise(diag_all, xx(:, pt))
    const = coef@const
    terms := toint(coef@terms)
    nterms := coef@nterms
    YFIT := coef@YFIT
    mcorr := coef@MCORRELATION
    corr := coef@correlation

    ; Print final regression model
    log_info("The final model (" + scal_time(0) + "-" + \
             (p_time(pt) + p_step) + "):")
    log_info("Constant: " + const)
    if (nterms .gt. 0) then
      do it = 0, dimsizes(coef) - 1
        log_debug(wdiag(terms(it)) + ", " + coef(it) + ", " + corr(it))
      end do
    end if
    log_info("Multiple correlation: " + mcorr)
    log_info("")

    ; -------------------------------------------------------------------------
    ; Model weighting
    ; -------------------------------------------------------------------------

    x_diag := transpose(diag_all(terms, :))
    good_diag := ind(.not. ismissing(x_diag(:, 0)))
    x := x_diag(good_diag, :)

    ; tmp_dat = DIM_DAT
    tmp_dat := dimsizes(x(:, 0))

    ; Create Design Matrix D
    D := new((/tmp_dat, dimsizes(terms) + 1/), typeof(diag_all))
    D = 1
    do id = 0, dimsizes(x(:, 0)) - 1
      D(id, 1:) = x(id, :)
    end do

    ; Observations
    xobs := new((/1, dimsizes(terms) + 1/), typeof(diag_all), \
                getFillValue(diag_all))
    xobs = 1.0
    xobs(0, 1:) = obs0(terms)
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
    log_info("Model weights: ")
    do idat = 0, DIM_DAT - 1
      log_info(dataset_info@dataset(idat) + ": " + weights_t(:, idat, pt))
    end do
    log_info("")

    ; Save Y-estimates and multiple correllation coeficients
    coeff_y_est(:, pt) = coef@Yest
    MCORR(pt) = mcorr
    Y_WEIGHT(pt) = y_weight

    ; Save coefficients and terms of regression model in array
    if (isvar("mderrm")) then
      idx = max((/dimsizes(coef), dimsizes(mderrm(0, 1:))/))
      tmprm = mderrm
      tmptr = termsrm
      tmpt  = termst
      delete([/mderrm, termsrm, termst/])

      mderrm = new((/pt + 1, idx + 1/), typeof(tmprm))
      termsrm = new((/pt + 1, idx + 1/), typeof(tmptr))
      termst = new((/pt + 1, idx + 1/), typeof(tmpt))

      mderrm(:pt - 1, :dimsizes(tmprm(0, :)) - 1) = tmprm
      termsrm(:pt - 1, :dimsizes(tmptr(0, :)) - 1) = tmptr
      termst(:pt-1, :dimsizes(tmpt(0, :)) - 1) = tmpt
      mderrm(pt, 0) = const
      mderrm(pt, 1:dimsizes(coef)) = coef
      termsrm(pt, 1:dimsizes(wdiag(terms))) = wdiag(terms)
      termst(pt, 1:dimsizes(terms)) = terms
      delete([/tmprm, tmptr, tmpt, idx/])
    else
      mderrm = new((/pt + 1, dimsizes(coef) + 1/), double)
      termsrm = new((/pt + 1, dimsizes(coef) + 1/), string)
      termst = new((/pt + 1, dimsizes(coef) + 1/), typeof(terms))

      mderrm(pt, 0) = const
      mderrm(pt, 1:) = coef
      termsrm(pt, 1:) = wdiag(terms)
      termst(pt, 1:) = terms
    end if

    delete([/y, corr, mcorr, x, x_diag, good_diag, xobs, cf, tmp_cf, \
            weights, uM, uM_std/])
  end do    ; Loop over different target periods

  ; Iterate over target periods and save MDER results
  do pt = 0, dimsizes(p_time) - 1
    tmp_yy = new(dimsizes(yy), typeof(yy))
    tmp_yywt = new(dimsizes(yywt), typeof(yywt))
    do tt = 0, n_max - 2
      yy(0, tt) = avg(yy(3:, tt))
      if (tt .ge. p_time(pt) - ny_min .and. \
          tt .lt. (p_time(pt) + p_step + 1 - ny_min)) then
        tmp_yy(1, tt) = sum(yy(3:, tt) * weights_t(0, :, pt))
        tmp_yywt(0, tt) = tmp_yy(1, tt) + sqrt(abs(sum(yy(3:, tt) * \
                                               cfint_t(:, pt))))
        tmp_yywt(1, tt) = tmp_yy(1, tt) - sqrt(abs(sum(yy(3:, tt) * \
                                               cfint_t(:, pt))))
      end if
      if (tt .eq. (p_time(pt) + p_step - ny_min)) then
        yy(1, tt - (p_step + 1) / 2)  = avg(tmp_yy(1, tt-p_step:tt))
        yywt(0, tt - (p_step + 1) / 2) = avg(tmp_yywt(0, tt - p_step:tt))
        yywt(1, tt - (p_step + 1) / 2) = avg(tmp_yywt(1, tt - p_step:tt))
      end if
      yy(2, tt) = sum(yy(3:, tt) * weights_t(0, :, pt))
      yyw(0, tt) = yy(2, tt) + sqrt(abs(sum(yy(3:, tt) * cfint_t(:, pt))))
      yyw(1, tt) = yy(2, tt) - sqrt(abs(sum(yy(3:, tt) * cfint_t(:, pt))))
    end do
  end do

  ; ---------------------------------------------------------------------------
  ; Plot preparations
  ; ---------------------------------------------------------------------------

  ; Plot regression line
  diag_reg = new((/dimsizes(p_time), DIM_DAT + 2/), double)
  nfit = fspan(0, dimsizes(diag_reg(0, 2:)) - 1, dimsizes(diag_reg(0, 2:)))
  ; nfit = fspan(0, 49, 50)
  diag_std = new((/dimsizes(p_time), DIM_DAT + 2/), double)
  wy = new((/dimsizes(p_time), DIM_DAT + 2/), double)
  wy_std = new((/dimsizes(p_time), DIM_DAT + 2/), float)
  yfit_reg = new((/dimsizes(p_time), dimsizes(nfit)/), double)
  yfit_std = new((/dimsizes(p_time), 2, dimsizes(nfit)/), double)

  ; Iterate over target periods
  do ppt = 0, dimsizes(p_time) - 1
    xgood = ind(.not. ismissing(termst(ppt, :)))
    hlp_tmp = new((/dimsizes(mderrm(ppt, xgood)), dimsizes(nfit)/), float)
    hlp = new((/1, dimsizes(mderrm(ppt, xgood)) + 1/), float)
    diag_reg(ppt, :) = mderrm(ppt, 0)         ; const
    diag_std(ppt, :) = 0.0
    yfit_reg(ppt, :) = mderrm(ppt, 0)         ; const
    hlp = 1.0

    wy(ppt, 0) = avg(xx(:, ppt))
    wy(ppt, 1) = Y_WEIGHT(ppt)
    wy(ppt, 2:) = xx(:, ppt)
    wy_std(ppt, 2:) = xx_val(:, ppt, 1)       ; cfint_t(:, ppt)

    do ii = 1, dimsizes(mderrm(ppt, :)) - 1
      if (.not. ismissing(mderrm(ppt, ii))) then
        diag_reg(ppt, 1) = diag_reg(ppt, 1) + mderrm(ppt, ii) * \
                                              obs0(termst(ppt, ii))
        diag_std(ppt, 1) = diag_std(ppt, 1) + (mderrm(ppt, ii) * \
                                               (obs0_std(termst(ppt, ii)))) ^ 2
        diag_reg(ppt, 2:) = coeff_y_est(:, ppt)

        ; diag_reg(2:) + coef(ii) * diag_all(terms(ii), :)
        diag_std(ppt, 2:) = diag_std(ppt, 2:) + \
          (mderrm(ppt, ii) * (diag_all_std(termst(ppt, ii), :))) ^ 2

        minx = min(diag_all(termst(ppt, ii), :)) - 1.1 * \
          (max(diag_all(termst(ppt, ii), :)) - \
           min(diag_all(termst(ppt, ii), :)))
        maxx = max(diag_all(termst(ppt, ii), :)) + 1.1 * \
          (max(diag_all(termst(ppt, ii), :)) - \
           min(diag_all(termst(ppt, ii), :)))
        diag_var = minx + (maxx - minx) * nfit / dimsizes(nfit)
        hlp_tmp(ii - 1, :) = diag_var
        delete([/diag_var, minx, maxx/])
      end if
    end do

    diag_std = sqrt(diag_std)
    x_diag = transpose(diag_all(termst(ppt, xgood), :))
    good_diag = ind(.not. ismissing(x_diag(:, 0)))
    x = x_diag(good_diag, :)

    ; Create Design matrix D and Identity matrix Id for pt = 1
    D := new((/dimsizes(x(:, 0)), dimsizes(termst(ppt, xgood)) + 1/), \
             typeof(diag_all))
    D = 1.0
    do id = 0, dimsizes(x(:, 0)) - 1
      D(id, 1:) = x(id, :)
    end do
    df = tmp_dat - dimsizes(obs0(termst(ppt, xgood))) + 1
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
      isort_tmp = dim_pqsort_n(diag_reg(ppt, 2:), 1, 0)
      isort = new(dimsizes(isort_tmp) + 2, integer)
      isort(0) = 0
      isort(1) = 1
      isort(2:) = isort_tmp + 2
      diag_reg(ppt, :) = diag_reg(ppt, isort)
      diag_std(ppt, :) = diag_std(ppt, isort)
      wy(ppt, :) = wy(ppt, isort)
      wy_std(ppt, :) = wy_std(ppt, isort)
      model_numbers = ispan(0, dimsizes(diag_reg(ppt, :)) - 1, 1)

      ; Save order of datasets
      diag_reg@model_numbers = model_numbers(isort)
      rc = regline_stats(diag_reg(ppt, 2:), log(wy(ppt, 2:)))
      yfit_reg(ppt, :) = exp(rc@b(0) + rc@b(1) * diag_reg(ppt, 2:))
      yfit_std(ppt, :, :) = (/exp(rc@y95(0) + rc@b95(0) * diag_reg(ppt, 2:)), \
                             exp(rc@y95(1) + rc@b95(1) * diag_reg(ppt, 2:))/)
      wy(ppt, 1) = exp(rc@b(0) + rc@b(1) * diag_reg(ppt, 1))   ; yfit(obs)
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
    good = ind(dataset_info@dataset(im) .ne. dataset_info@dataset)
    do pt = 0, dimsizes(p_time) - 1
      xx_pre(im, pt) = avg(xx(good, pt))
    end do
    delete(good)
  end do

  ; Calculating statistics
  x25 = round(0.25 * DIM_DAT, 3) - 1
  x75 = round(0.75 * DIM_DAT, 3) - 1

  presq = new((/dimsizes(p_time), DIM_DAT/), double)
  premdersq = new((/dimsizes(p_time), DIM_DAT/), double)
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

  ; Add uMMM and MDER to dataset_info@dataset
  tmpm = dataset_info@dataset
  dataset_info@dataset := new(DIM_DAT + 3, string)
  dataset_info@dataset(0) = "uMMM"
  dataset_info@dataset(1) = "MDER time dependent"
  dataset_info@dataset(2) = "MDER"
  dataset_info@dataset(3:) = tmpm
  delete(tmpm)
  DIM_DAT = DIM_DAT + 2

  ; ---------------------------------------------------------------------------
  ; Write regression model equations to ASCII table
  ; ---------------------------------------------------------------------------

  RegStg = new((/dimsizes(p_time)/), string)
  PStg = new((/dimsizes(p_time)/), string)
  do pt = 0, dimsizes(p_time) - 1
    RegStg(pt) = sprintf("%4.2f", mderrm(pt, 0))
    do ii = 1, dimsizes(mderrm(pt, :)) - 1
      if (.not. ismissing(mderrm(pt, ii))) then
        if (mderrm(pt, ii) .gt. 0.0) then
          RegStg(pt) = RegStg(pt) + " +" + sprintf("%4.2f", mderrm(pt, ii)) + \
            " x " + termsrm(pt, ii)
        else
          RegStg(pt) = RegStg(pt) + " " + sprintf("%4.2f", mderrm(pt, ii)) + \
            " x " + termsrm(pt, ii)
        end if
      end if
    end do
    p1_time = p_time(pt) + p_step - 2000      ; p_time(0)
    PStg(pt) = p_time(pt) + "-" + sprinti("%0.2hi", p1_time)
  end do

;   header = (/"Target Period             Model Equation  R^2",\
;              "                          "+models@experiment(0+DIM_OBS+dim_hMOD)+"         "/)
;   hlist = [/header/]
;   alist = [/PStg,RegStg,MCORR/]
;
;   write_table(plot_dir+plot_file+".txt", "w", hlist, "%s")
;   write_table(plot_dir+plot_file+".txt", "a", alist, "%s %s %f")
;   delete([/header,hlist,alist/])
;
;   asciiwrite(climo_dir + "/" + target_var + "_" + \
;              models@experiment(0+DIM_OBS+dim_hMOD) + "_" +(p_step+1) + \
;              "ystep_diagnostics.txt", termsrm)
;   ;;--------------------------------------------------------------------------
;   ;;Write weights to NetCDF
;   model_list = dataset_info@dataset(3:)
;   do pt = 0, dimsizes(p_time)-1
;     weights_out = weights_t(0,:,pt)
;     weights_filename = plot_dir + plot_file + "_weights_" + \
;                        p_time(pt) + "-" + (p_time(pt)+p_step)
;     weights_out!0           = "model"
;     weights_out&model       = model_list
;     weights_out@diag_script = diag_script
;     weights_out@var         = "weight"
;     weights_out@description = "Weights from MDER"
;     if (target_var.eq."SIE") then
;       weights_out@description = "Weights for September Arctic sea ice extent"
;     end if
;     weights_file            = ncdf_write(weights_out, weights_filename+".nc")
;     delete([/weights_out, weights_filename, weights_file/])
;   end do
;   ;;--------------------------------------------------------------------------
;   ;;PLOT
;
;   ;; add meta data to plot (for reporting)
;   climofiles = new(DIM_DAT, string)
;   do ii = 0, DIM_DAT - 1
;       climofiles(ii) = interface_get_inpaths(ii) + "/" \
;                        + interface_get_infile(variables(0), field_types(0), ii)
;   end do
;   allmodelnames = dataset_info@dataset
;   zz = diag_script + "_" + variables(0)
;   contrib_authors = (/"A_wenz_sa"/)
;
;   ;; Line plot showing the time series including MDER weights
;   wks     = get_wks(0, diag_script, plot_file + "_FIG1")
;   MainStg = ""
;   YStg    = target_var+" ["+variable_info@units+"]"
;   plot    = weight_lin_plot(yy,yywt,yyw,yy_obs,MainStg,YStg,wks,"draw")
;
;   alltags = array_append_record(tags, (/"ST_mder", "ST_corr", "ST_range", \
;                                         "ST_stddev", "PT_times", "DM_sh"/), 0)
;   caption = "Time series of " + variables(0) + " for " + models@ensemble(0) + \
;             " scenario, similar to Wenzel et al., JClim, 2016, fig. 5"
;
;   plot_fname = wks@fullname
;   ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
;            allmodelnames, climofiles, diag_script, contrib_authors)
;   delete([/alltags, caption, wks, plot_fname/])
;
;   do ppt=0, dimsizes(p_time)-1
;     wks     = get_wks(0, diag_script, plot_file + "_FIG2" + plot_num(ppt))
;     MainStg = plot_num(ppt)+") "+p_time(ppt)+"-"+(p_time(ppt)+p_step)
;     XStg    = RegStg(ppt)+" ["+variable_info@units+"]"
;     YStg    = target_var+" ["+variable_info@units+"]"
; print(diag_reg(ppt,:))
; print(diag_std(ppt,:))
; print(wy(ppt,:))
; print(wy_std(ppt,:))
; print(yfit_reg(ppt,:))
; print(yfit_std(ppt,:,:))
;     plot2   = mder_scatter_plot(diag_reg(ppt,:), diag_std(ppt,:), \
;                                 wy(ppt,:), wy_std(ppt,:), \
;                                 yfit_reg(ppt,:), yfit_std(ppt,:,:), \
;                                 MainStg, XStg, YStg, wks, "draw")
;     alltags = array_append_record(tags, (/"ST_mder", "ST_corr", "ST_range", \
;                                           "ST_stddev", "ST_corr", "PT_scatter",\
;                                           "T_EC", "DM_sh"/), 0)
;     caption = "Scatterplot showing the correlation between future " + \
;               variables(0) + " and the quantity resulting from MDER analysis" +\
;               " on the x-axes, similar to Wenzel et al., JClim, 2016, fig. 3"
;     plot_fname = wks@fullname
;     ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
;              allmodelnames, climofiles, diag_script, contrib_authors)
;     delete([/alltags, caption, wks, plot_fname/])
;   end do  ; ppt
; exit
;
;   YStg   = "RMSE projection error ["+variable_info@units+"]"
;   XBname = new(dimsizes(PStg)*2,string)
;
;   if (dimsizes(XBname).eq.2) then
;       XBname = (/"uMMM","MDER"/)
;   else if (dimsizes(XBname).eq.4) then
;       XBname(0::2) = (/"uMMM~C~"+PStg/)
;       XBname(1::2) = (/"MDER~C~"+PStg/)
;   else
;       XBname(0::2) = (/" "/)
;       XBname(1::2) = PStg
;   end if
;   end if
;
;   ;; Boxplot of RMSE of different weigting methods inkl pseudoreality check
;   wks     = get_wks(0, diag_script, plot_file + "_FIG3")
;   MainStg = ""
;   Xtitle  = PStg(0)
;   do pt = 1, dimsizes(p_time)-1
;     Xtitle = Xtitle+"  "+PStg(pt)
;   end do
;   plot3   = sqrt_error_plot(preRDS,preAVG,prePRE,\
;                            MainStg,XBname,Xtitle,YStg,wks,"draw")
;
;   alltags = array_append_record(tags, (/"ST_mder", "ST_rmsd", "ST_stddev", \
;                                 "ST_mean", "PT_bar", "PT_errorbar", "DM_sh"/),0)
;   caption = "RMSE difference between the ensemble mean future climatological" +\
;             " mean " + variables(0) + "and its equivalent in pseudoreality " +\
;             "under " + models@ensemble(0) + " scenario, similar to Wenzel " +\
;             "et al., JClim, 2016, fig. 4"
;   plot_fname = wks@fullname
;   ESMValMD(plot_fname, alltags, caption, zz, variables(0), \
;            allmodelnames, climofiles, diag_script, contrib_authors)
;   delete([/alltags, caption, zz, climofiles, allmodelnames, plot_fname/])
;
;   ;; Remove temporary diagnostic list
;   system("rm -f " + temp_dir)

  leave_msg(diag_script, "")

end
