; #############################################################################
; mder/absolute_correlation.ncl
; #############################################################################
; MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; Author: Sabrina Wenzel (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     This routine calculates absolute correlation coefficients s between all
;     diagnostics, calculated for a present-day period, called in the namelist
;     and a future variable of interest.
;
;     Calls one plotscript:
;     Plots the absolute correlation coefficient between the variable of
;     interest and present-day diagnostic.
;
; Required recipe attributes:
;     p_time : time range for projections.
;     p_step : averaging time range for target period.
;     scal_time : a year to scale the timeseries with.
;     time_oper : select for opperation in time_opperation function.
;     time_opt : select for option in time_opperation function.
;     calc_type : trend/pos
;
; Optional recipe attributes:
;
; Caveats:
;     Modifies dataset_info@dataset and dataset_info@ensemble.
;
; Modification history:
;   * 20180726_A-schl_ma: ported to v2.0
;   * 201307??_A-wenz_sa: written
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/latlon.ncl"

load "./diag_scripts/shared/plot/mder.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

load "./diag_scripts/shared/mder.ncl"


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
  if (avgens .eq. "True") then
    idx_datasets = UNIQ(dataset_info@dataset(DIM_OBS + DIM_HDAT:)) + \
                   DIM_OBS + DIM_HDAT
    avgens = True
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
    avgens = False
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

  ; Global attributes
  if (isatt(diag_script_info, "p_time"))then
    p_time = diag_script_info@p_time
  else
    error_msg("f", diag_script, "", "requires attribute p_time in recipe")
  end if
  if (isatt(diag_script_info, "p_step"))then
    p_step = diag_script_info@p_step
  else
    error_msg("f", diag_script, "", "requires attribute p_step in recipe")
  end if
  if (isatt(diag_script_info, "scal_time"))then
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

  ; Output dir
  diag_script_base = basename(diag_script)
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  plot_num = integertochar(ispan(97, 122, 1))

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Output arrays
  xx = new((/dimsizes(p_time), DIM_DAT/), float)
  future_diag = new((/DIM_DAT/), float)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Call the variable of interest and merge to one timeline
  do i = 0, DIM_DAT - 1
    idat = idx_datasets(i)

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

    ; Calculate future diagnostic for diagnostic scatter plot
    future_tmps = time_operations(tmpy, p_time(0), \
                                  (/p_time(dimsizes(p_time) - 1) + p_step/), \
                                  time_oper, time_opt, True)
    future_diag(i) = tofloat(avg(time_operations( \
        future_tmps, p_time(0), (/p_time(dimsizes(p_time)-1) + p_step/), \
        "average", "yearly", True)))
    delete(future_tmps)

    ; Loop over different target periods
    do pt = 0, dimsizes(p_time) - 1
      ; Average over var time series
      tmps = time_operations(tmpy, p_time(pt), p_time(pt) + p_step, \
                             time_oper, time_opt, True)
      tmpj = time_operations(tmps, p_time(pt), p_time(pt) + p_step, \
                             "average", "yearly", True)

      ; Calculate the scaling factor between scal_time(0) and scal_time(1)
      if (scal_time(1) .gt. \
          toint(input_file_info@end_year(idat - DIM_HDAT))) then
        scal_time(1) = toint(input_file_info@end_year(idat - DIM_HDAT))
      end if

      ; Calculate average over starting period for trend calculation
      tmpxx = time_operations(tmpy, scal_time(0), scal_time(1), "extract", \
                              time_opt, True)
      tmpyy = time_operations(tmpxx, scal_time(0), scal_time(1), "average", \
                              "yearly", True)
      tmppd = tofloat(avg(tmpyy))

      ; Perform calculation type
      if (calc_type .eq. "trend") then
        xx(pt, i) = tofloat(avg(tmpj)) - tmppd
      else if (calc_type .eq. "pos") then
        xx(pt, i) = tofloat(avg(tmpj))
      else
        error_msg("f", diag_script, "", "calc_type " + calc_type + \
                  " not supported")
      end if
      end if
      delete([/tmppd, tmps, tmpj, tmpyy, tmpxx/])
    end do
    delete(tmpy)
  end do

  ; Save unique dataset names (modifies dataset_info@dataset!)
  tmpm = dataset_info@dataset
  dataset_info@dataset := tmpm(idx_datasets)
  delete(tmpm)

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
  diag_all = new((/DIM_DIAG, DIM_DAT/), float)
  wdiag = new((/DIM_DIAG/), string)           ; Array with diagnostic names

  ; Loop over files in the list and append data
  do idx_diag = 0, DIM_DIAG - 1
    data_temp = ncdf_read(file_list(idx_diag), "diag")

    ; Loop over datasets and save diagnostic data
    do i = 0, DIM_DAT - 1
      do idat = 0, dimsizes(data_temp&datasets) - 1
        if (dataset_info@dataset(i) .eq. data_temp&datasets(idat)) then
          diag_all(idx_diag, i) = data_temp(0, idat, 0)
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
    delete(data_temp)
  end do
  diag_all&diagnostics = wdiag

  ; Create scatterplots for each diagnostic
  res_scat = True
  scatterplots = diag_scatter(res_scat, diag_all, future_diag, target_var)

  ; ---------------------------------------------------------------------------
  ; Correlation calculations
  ; ---------------------------------------------------------------------------

  ; New arrays
  corr = new((/dimsizes(p_time), DIM_DIAG/), float)
  diag_sign = new((/dimsizes(p_time), DIM_DIAG/), string)

  ; Loop over different target periods
  do pt = 0, dimsizes(p_time) - 1
    corr_sig = new((/2, DIM_DIAG/), float)
    r_sig0  = new((/DIM_DIAG/), float)
    r_sig1  = new((/DIM_DIAG/), float)
    r_wdiag  = new((/DIM_DIAG/), string)

    ; Iterate over diagnostics
    do idx_diag = 0, DIM_DIAG - 1
      hlp = diag_all(idx_diag, :)

      ; Only use models with data (filter missing values)
      good = ind(hlp .ne. 999.0)
      ngood = dimsizes(good)

      ; Calculate correlations
      corr(pt, idx_diag) = escorc(hlp(good), xx(pt, good))
      rc = regline(hlp(good), xx(pt, good))
      df = rc@nptxy - 2
      z = 0.5 * log((1.0 + abs(corr(pt, idx_diag))) / \
                    (1.0 - abs(corr(pt, idx_diag))))
      corr_sig(0, idx_diag) = tanh(z - 1.96 / (sqrt(ngood - 3.0)))
      corr_sig(1, idx_diag) = tanh(z + 1.96 / (sqrt(ngood - 3.0)))

      t0 = abs(corr(pt, idx_diag)) / sqrt((1 - corr(pt, idx_diag) ^ 2) / \
                                          (ngood-2))
      sign0 = 1.0 - betainc(df / (df + rc@tval ^ 2), df / 2.0, 0.5)

      log_info(wdiag(idx_diag) + " corr: " + corr(pt, idx_diag) + \
               "  datasets: " + ngood + "  stddv: " + sign0(0))
      delete([/hlp, good, ngood, rc, df/])
    end do

    ; Sort abs(coerr. coef.) to decreasing order in array
    acorr = abs(corr(pt, :))
    qsort(acorr)
    acorr = acorr(::-1)
    do idx_diag = 0, DIM_DIAG - 1
      i_diag = ind(abs(corr(pt, :)) .eq. acorr(idx_diag))
      r_wdiag(idx_diag) = wdiag(i_diag)
      r_sig0(idx_diag) = corr_sig(0, i_diag)
      r_sig1(idx_diag) = corr_sig(1, i_diag)
      diag_sign(pt, idx_diag) = where(corr(pt, i_diag) .gt. 0.0, "red", "blue")
      delete(i_diag)
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

    filename  = target_var + "-" + calc_type + "_" + \
      dataset_info@exp(DIM_OBS + DIM_HDAT) + "_" + \
      p_time(pt) + "-" + (/p_time(pt) + p_step/)
    wks = gsn_open_wks(file_type, plot_dir + filename)
    MainStg = plot_num(pt) + ") " + p_time(pt) + "-" + (/p_time(pt) + p_step/)
    YStg = "Absolute Correlation Coefficient"
    plot = error_bar_plot(acorr, corr_sig, wdiag, diag_sign(pt, :), MainStg, \
                          YStg, wks, "draw")

    ; Add meta data (TODO)
    ; climofiles = new(DIM_DAT, string)
    ; do ii = 0, DIM_DAT - 1
    ;   climofiles(ii) = interface_get_inpaths(ii) + "/" + \
    ;     interface_get_infile(var0, field_type0, ii)
    ; end do
    ; alltags = array_append_record(tags, \
    ;                               (/"ST_mder","ST_spectrum", "ST_corr", \
    ;                                "ST_range", "PT_errorbar", "PT_scatter", \
    ;                                "DM_et", "DM_global", "DM_polar", \
    ;                                "DM_sh", "DM_trop"/), \
    ;                               0)
    ; alldatasetnames = dataset_info@dataset
    ; caption = "Absolute values of the correlation coefficient between " + \
    ;   var0  + " and pressent-day diagnostics, similar to Wenzel et al., " + \
    ;   "JClim, 2016, fig. 2"
    ; id = diag_script + "_" + var0
    ; contrib_authors = (/"A_wenz_sa", "A_schl_ma"/)
    ;
    ; ESMValMD(plot_dir + filename + "." + file_type, alltags, caption, id, \
    ;          var0, alldatasetnames, climofiles, diag_script, contrib_authors)
    ;
    ; delete([/alltags, caption, id, climofiles, alldatasetnames, wks, plot, \
    ;         MainStg, acorr/])

    ; Restore old order
    wdiag = backup_wdiag
    delete([/corr_sig, acorr, backup_wdiag/])
  end do

  leave_msg(diag_script, "")

end
