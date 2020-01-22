; #############################################################################
; shared/mder.ncl
; #############################################################################
; AUXILIARY FUNCTIONS FOR MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; #############################################################################
; Please consider using of extending existing routines before adding new
; ones. Check the header of each routine for documentation.
;
; Contents:
;     function get_obs_list: get list of OBS datasets
;     function diag_file_name: filename of calculated diagnostic.
;     function integrate_time_units: temporally integrates units.
;     function mder_output_desired: check if MDER output is desired.
;     function reg_multlin_stats_rap: costumized version of NCL's
;                                     reg_multlin_stats.
;     function regress_stepwise: does stepwise regression on diagnostics.
;
; #############################################################################


; #############################################################################
undef("get_obs_list")
function get_obs_list(input_list: list)
;
; Arguments:
;     input_list: List of input metadata
;
; Return value:
;     list: All OBS datasets.
;
; Description:
;     Returns OBS and OBS6 datasets.
;
; Modification history:
;     20191120-schlund_manuel: written.
;
local atts, obs, obs6

begin

  funcname = "get_obs_list"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  atts := True
  atts@project = "OBS"
  obs = select_metadata_by_atts(input_list, atts)
  atts@project = "OBS6"
  obs6 = select_metadata_by_atts(input_list, atts)
  do iobs = 0, ListCount(obs6) - 1
    ListAppend(obs, obs6[iobs])
  end do

  leave_msg(scriptname, funcname)
  return(obs)

end


; #############################################################################
undef("diag_file_name")
function diag_file_name(filename[1]: string, \
                        opt[1]: logical)
;
; Arguments:
;     filename: Name of the current file.
;     opt: options with the optional attributes
;          opt@wdiag (diagnostic name).
;          opt@period_avg (year-, month- or season name).
;          opt@level (specify plev range).
;          opt@region (specify lat and/or lon range).
;
; Return value:
;     file_name: returns the file name of the calculated climatological
;                diagnostic.
; Description:
;     Provides a file name for the calculation results of diagsnostics.
;
; Modification history:
;     20180717-schlund_manuel: ported to v2.0.
;     20140723-wenzel_sabrina: written.
;
local outname

begin

  funcname = "diag_file_name"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; File name
  outname = basename(systemfunc("basename " + filename))

  ; Add extra info to the file
  if (isatt(opt, "wdiag")) then
    outname = outname + "_" + opt@wdiag
  end if
  if (isatt(opt, "period_avg")) then
    if (dimsizes(opt@period_avg) .eq. 2) then
      outname = outname + "_" + str_sub_str(opt@period_avg(0) + \
                                            opt@period_avg(1), " ", "")
    else
      outname = outname + "_" + opt@period_avg
    end if
  end if
  if (isatt(opt, "level")) then
    if (dimsizes(opt@level) .eq. 2) then
      outname = outname + "_" + str_sub_str(opt@level(0) + \
                                            opt@level(1), " ", "")
    else
      outname = outname + "_" + str_sub_str(opt@level, " ", "")
    end if
  end if
  if (isatt(opt, "region")) then
    if (dimsizes(opt@region) .eq. 4) then
      outname = outname + "_" + str_sub_str(opt@region(0) + opt@region(1) + \
                                            opt@region(2) + opt@region(3), \
                                            " ", "")
    else if (dimsizes(opt@region) .eq. 2) then
      outname = outname + "_" + str_sub_str(opt@region(0) + opt@region(1), \
                                            " ", "")
    else
      outname = outname + "_" + str_sub_str(opt@region, " ", "")
    end if
    end if
  end if
  outname = outname + ".nc"

  ; Return file name
  leave_msg(scriptname, funcname)
  return(outname)

end


; #############################################################################
undef("integrate_time_units")
function integrate_time_units(units[1] : string, \
                              opt[1] : integer)
;
; Arguments:
;     opt: integer determining how to handle multiple time units.
;          opt < 0: remove last units.
;          opt = 0: throw error if there are multiple time units.
;          opt = n > 1: remove n-th time units.
;
; Return value:
;     New units.
;
; Description:
;     Temporally integrates units, i.e. removes reciprocal time units.
;
; Caveat:
;     Units with higher exponents than -1 cannot be processed, e.g. "m s-2" or
;     "m/s2".
;
; Modification history:
;     20180927-schlund_manuel: written.
;
local allowed_time_units, allowed_prefixes, allowed_suffixes, iu, ip, \
      matches, dim_time, all_matches, all_lengths

begin

  funcname = "integrate_time_units"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Allowed time units
  allowed_time_units = (/"s", "min", "h", "hr", "d", "day", "mon", "y", "yr", \
                        "dec"/)
  allowed_prefixes = (/"/"/)
  allowed_suffixes = (/"-1", "^-1"/)

  ; Find units by prefix
  do iu = 0, dimsizes(allowed_time_units) - 1
    do ip = 0, dimsizes(allowed_prefixes) - 1
      time = allowed_prefixes(ip) + allowed_time_units(iu)
      matches := str_index_of_substr(units, time, 0)
      if (.not. all(ismissing(matches))) then
        dim_time := new(dimsizes(matches), integer)
        dim_time = strlen(time)
        if (.not. isvar("all_matches")) then
          all_matches = matches
          all_lengths = dim_time
        else
          all_matches := array_append_record(all_matches, matches, 0)
          all_lengths := array_append_record(all_lengths, dim_time, 0)
        end if
      end if
    end do
  end do

  ; Find units by suffix
  do iu = 0, dimsizes(allowed_time_units) - 1
    do is = 0, dimsizes(allowed_suffixes) - 1
      time = allowed_time_units(iu) + allowed_suffixes(is)
      matches := str_index_of_substr(units, time, 0)
      if (.not. all(ismissing(matches))) then
        dim_time := new(dimsizes(matches), integer)
        dim_time = strlen(time)
        if (.not. isvar("all_matches")) then
          all_matches = matches
          all_lengths = dim_time
        else
          all_matches := array_append_record(all_matches, matches, 0)
          all_lengths := array_append_record(all_lengths, dim_time, 0)
        end if
      end if
    end do
  end do

  ; Process different options
  if (isvar("all_matches")) then
    if (opt .lt. 0) then
      max_ind = max(all_matches)
    else if (opt .eq. 0) then
      max_ind = all_matches(0)
      if (dimsizes(all_matches) .gt. 1) then
        if (.not. all(all_matches .eq. max_ind)) then
          error_msg("f", scriptname, funcname, "temporal unit integration " + \
                    "failed for units " + units + ", ambiguous time " + \
                    "units given. Use different 'opt' parameter to enable " + \
                    "this conversion.")
        end if
      end if
    else
      unique_matches := get_unique_values(all_matches)
      qsort(unique_matches)
      if (opt .gt. dimsizes(unique_matches)) then
        error_msg("w", scriptname, funcname, "cannot remove time units " + \
                  "at position " + opt + ", removing last appearance")
        max_ind = max(all_matches)
      else
        max_ind = unique_matches(opt - 1)
      end if
    end if
    end if

    ; Remove time units
    max_len = max(where(all_matches .eq. max_ind, all_lengths, 0))
    if (max_len .eq. strlen(units)) then
      new_units = "1"
    else
      new_units = stringtochar(units)
      new_units := chartostring(new_units(:max_ind - 1)) + \
        chartostring(new_units(max_ind + max_len:))
    end if
  else
    error_msg("w", scriptname, funcname, "no reciprocal time units " + \
              "found in " + units)
  end if

  ; Return file name
  leave_msg(scriptname, funcname)
  return(str_squeeze(new_units))

end


; #############################################################################
undef("mder_output_desired")
function mder_output_desired()
;
; Return value:
;     True if mder output is desired, False if not.
;
; Description:
;     Checks if MDER output of diagnostic is desired.
;
; Modification history:
;     20180906-schlund_manuel: written.
;
local out

begin

  funcname = "mder_output_desired"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Check diag_script_info@wdiag
  if (isatt(diag_script_info, "wdiag") .and. \
      .not. any(diag_script_info@wdiag .eq. "")) then
    if (config_user_info@write_netcdf) then
      if (.not. isatt(diag_script_info, "ref_dataset")) then
        error_msg("w", scriptname, funcname, "Skipping MDER output: " + \
                  "attribute 'ref_dataset' not given in recipe")
        out = False
      else
        log_info("Writing MDER output with the following diagnostics: " + \
                 diag_script_info@wdiag)
        out = True
        out@work_dir = config_user_info@work_dir
        system("mkdir -p " + out@work_dir)
      end if
    else
      error_msg("w", scriptname, funcname, "Skipping MDER output: " + \
                "'write_netcdf' is set to False in user configuration")
      out = False
    end if
  else
    log_info("Skipping MDER output: 'wdiag' is not given in recipe")
    out = False
  end if

  ; Return boolean
  leave_msg(scriptname, funcname)
  return(out)

end


; #############################################################################
undef("reg_multlin_stats_rap")
function reg_multlin_stats_rap(Y[*]: numeric, \
                               XP:numeric, \
                               opt)
;
; Arguments:
;     Y: dependent variable (size NY). Missing values (_FillValue) are not
;        allowed.
;     XP: one [*] or more [*][*] independent variables. Missing values
;         (_FillValue) are not allowed. The size of the leftmost dimension of
;         XP must be the same size as Y[*]. The rightmost dimension for [*][*]
;         contains the independent variables.
;     opt: options [type logical].
;          opt = False: no options.
;          opt = True.
;          opr@print_data = True: print input Y and XP in table form.
;          opr@print_anova= True: print ANOVA information.
;
; Return value:
;     See reg_mullin_stats.
;
; Description:
;     This function is based on the NCL-function "reg_multlin_stats". The
;     original can be found in:
;         "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl".
;     See:
;         Wilks, D.S. (2006): Statistical Methods in the Atmospheric Sciences
;                             Academic Press (2nd Edition).
;         Davis, J.C. (2002): Statistics and Data Analysis in Geology, Wiley
;                             (3rd Edition), 462-470.
;     Very readable discussion:
;         http://reliawiki.org/index.php/Multiple_Linear_Regression_Analysis.
;
; Caveats:
;     In this version the following attributes are not available:
;         b@stderr = stderr.
;         b@tval = tval.
;         b@pval = pval.
;
; Modification history:
;     20180808-schlund_manuel: ported to v2.0.
;     2013????-wenzel_sabrina: written.
;
local dimXP, N, NP, NY, M, T, Y, X, b, B, rankXP, Yavg, Yvar, Ystd, Xavg, \
      Xvar, Xstd, Yest, n, m, Yres, SST, SSE, SSR, MST, MSE, MSR, SE, Frat, \
      r2, r2a, r, fuv, XXt, XXTi, varcovx, df, stderr, tval, pval, opt_wm, \
      fmt_XXt, fmt_XXti, fmt_covx, opt_pd, fmt_data

begin

  funcname = "reg_multlin_stats_rap"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Arrays
  dimXP = dimsizes(XP)
  NX = dimXP(0)              ; number of rows ('observations')
  N = dimsizes(Y)
  rankXP = dimsizes(dimXP)

  ; Check input
  if (N .ne. NX) then
    error_msg("f", scriptname, funcname, "number of rows of Y and XP must " + \
              "match, NY = " + N + ", NXP = " + NX)
  end if
  if (rankXP.gt.2) then
    error_msg("f", scriptname, funcname, "rank of XP > 2: rankXP = " + \
              rankXP)
  end if
  ; if (isatt(Y, "_FillValue") .and. any(ismissing(Y))) then
  ;   error_msg("f", scriptname, funcname, "Y has missing values, not allowed")
  ; end if
  ; if (isatt(XP, "_FillValue") .and. any(ismissing(XP))) then
  ;   error_msg("f", scriptname, funcname, "XP has missing values, " + \
  ;             "not allowed")
  ; end if

  ; Get number of predictor variables
  if (rankXP .eq. 2) then
    NP = dimXP(1)
  else
    NP = 1
  end if

  ; Print data
  if (opt .and. isatt(opt, "print_data") .and. opt@print_data) then
    log_info("")
    opt_pd = True
    opt_pd@title  = "----- reg_multlin_stats: Y, XP -----"
    data = new((/N, NP + 1/), typeof(Y))
    data(:, 0) = (/Y/)

    if (NP.eq.1) then
      data(:, 1) = (/XP/)
    else
      data(:, 1:) = (/XP/)
    end if

    fmt_data   = "f11.2"
    if (isatt(opt, "fmt_data")) then
      fmt_data  = opt@fmt_data          ; e.g.: "f8.3", "e13.5", etc.
    end if
    write_matrix(data, (NP + 1) + fmt_data, opt_pd)
    log_info("")
    delete(data)
  end if

  ; Create the necessary arrays. X is a "design matrix":
  ;     http://en.wikipedia.org/wiki/Design_matrix
  ; Most multiple regression models include a constant term. This ensures that
  ; the model will be "unbiased", i.e. the mean of the residuals will be
  ; exactly zero.

  ; Extra column for design matrix (required ordering for reg_multlin)
  M = NP + 1
  X = new((/M, N/), typeof(XP), getFillValue(XP))
  X(0, :) = 1.0

  if (NP .eq. 1) then
    X(1, :) = XP                ; One predictor
  else
    do m = 0, NP - 1
      X(m + 1, :) = XP(:, m)    ; m refers to predictor variables
    end do
  end if

  ; Calculate the partial regression coefficients: b -> beta_hat
  ;     b = (X'X)^-1 X'y (theoretically)
  ; http://www.ncl.ucar.edu/Document/Functions/Built-in/reg_multlin.shtml

  ; Partial regression coef (unstandardized...depends on units)
  b = reg_multlin(Y, X, False)
  delete(b@constant)

  ; Assorted calculations
  Yavg = avg(Y)
  Yvar = variance(Y)            ; Unbiased estimate the variance
  Ystd = sqrt(Yvar)

  Xavg = new(NP, typeof(b), "No_FillValue")
  Xvar = new(NP, typeof(b), "No_FillValue")
  Xstd = new(NP, typeof(b), "No_FillValue")

  ; Standardized regression coef.: describes what the partial regression coef.
  ; would equal if all variables had the same standard deviation.
  bstd = new(M, typeof(b), "No_FillValue")

  do n = 0, NP - 1
    Xavg(n) = avg(X(n + 1, :))
    Xvar(n) = variance(X(n + 1, :))
    Xstd(n) = sqrt(Xvar(n))
    bstd(n + 1) = b(n + 1) * Xstd(n) / Ystd         ; don't do b(0)
  end do

  ; Calculate model estimates: Yest = b(0) + b(1)* + b(2)* + ...
  Yest = new(N, typeof(Y), getFillValue(Y))
  do n = 0, N - 1
    Yest(n) = b(0) + sum(b(1:) * X(1:, n))          ; array syntax
  end do

  ; ANOVA table info. (Eg, Wilks pg 185 & pg 197)
  ; Sums of squares , degrees of freedom, mean squares, F statistic
  ; Other statistics are also calculated.
  ;     http://reliawiki.org/index.php/
  ;     ANOVA_Calculations_in_Multiple_Linear_Regression

  ; Mean squares are obtained by dividing the sum of squares by their
  ; associated degrees of freedom.
  dof  = N - NP - 1
  Yres = Yest-Y                 ; residuals (array operation)

                                ; SS ==> Sum of Squares
  SST = sum((Y - Yavg) ^ 2)     ; total SS [dof = N - 1]
  SSR = sum((Yest - Yavg) ^ 2)  ; regression SS [dof = NP]
  SSE = sum(Yres ^ 2)           ; residual SS [dof = N - NP - 1]

  MST = SST / (N - 1)           ; Mean SST [dof = N - 1]
  MSR = SSR / NP                ; Mean SSR [dof = NP]
  MSE = SSE / dof               ; Mean SSE [dof = (N - NP - 1)]

  RSE = sqrt(MSE)               ; residual standard error

  ; F-Test is an overall (global) test of the models fit. If at least one
  ; variable has a significant coefficient, then the model fit should be
  ; significant. The F-statistic might be interpreted as the variation
  ; explained by the regression relative to the variation not explained by the
  ; regression.
  ; R dof=(M,N-NP-1)
  Frat = MSR / MSE

  ; The square of the Pearson correlation coefficient between the observed
  ; and modeled (predicted) data values of the dependent variable. It is
  ; another measure of 'goodness of fit.' biased high, particularly for small N
  ; r2 = coefficient of determination. It is
  r2 = SSR / SST

  r2a = r2 - (1 - r2) * NP / tofloat(dof)   ; adjusted r2... better for small N
  r = sqrt(r2)                              ; multiple (overall) correlation
  fuv = 1.0 - r2                            ; fraction of variance unexplained

  ; Determine (a) standard error of coefficients; (b) t-values; (c) p-values

  XXt     = X#transpose(X)    ; (M,M); opposite of text books
  XXti    = inverse_matrix(XXt)
  varcovx = MSE*XXti          ; variance-covariance matrix
                              ; of estimated regression coefficients

  ; Debug log
  if (config_user_info@log_level .eq. "debug") then
    log_debug("")
    opt_wm = True
    opt_wm@title  = "---------   XXt  ---------"
    fmt_XXt = "f11.2"
    if (isatt(opt, "fmt_XXt")) then
      fmt_XXt  = opt@fmt_XXt                ; "f8.2", "f12.4", "e13.5", ...
    end if
    write_matrix(XXt, M + fmt_XXt, opt_wm)
    opt_wm@title  = "---------   XXti ---------"
    fmt_XXti = "f11.2"
    if (isatt(opt, "fmt_XXti")) then
      fmt_XXti = opt@fmt_XXti
    end if
    write_matrix(XXti, M + fmt_XXti, opt_wm)
    opt_wm@title  = "--------- varcovx -----k----"
    fmt_covx = "f11.2"
    if (isatt(opt, "fmt_covx")) then
      fmt_covx = opt@fmt_covx
    end if
    write_matrix(varcovx, M + fmt_covx, opt_wm)
    log_debug("")
  end if

  stderr = new(M, typeof(Y))
  tval = new(M, typeof(Y))
  pval = new(M, typeof(Y))
  df  = N - NP

  ; do m = 0, M - 1
  ;   stderr(m) = sqrt(varcovx(m, m))
  ;   tval(m) = b(m) / stderr(m)
  ;   pval(m) = betainc(df / (df + tval(m) ^ 2), df / 2.0, 0.5)
  ; end do

  ; Print ANOVA information
  if (opt .and. isatt(opt, "print_anova") .and. opt@print_anova) then
    log_info("")
    log_info("--------- ANOVA information ---------")
    log_info("")
    log_info("SST = " + SST + ", SSR = " + SSR + ", SSE = " + SSE)
    log_info("MST = " + MST + ", MSR = " + MSR + ", MSE = " + MSE + \
             ", RSE = " + RSE)
    log_info("F-statistic = " + Frat + ", dof = (" + NP + "," + \
             (N - NP - 1) + ")")
    log_info("---------                   ---------")
    log_info("r2 = " + r2 + ", r = " + r + ", r2a = " + r2a + ", fuv = " + fuv)
    log_info("---------                   ---------")
    log_info("stderr, tval, pval: " + stderr + ", " + tval + ", " + pval)
    log_info("")
  end if

  ; Associate assorted information as attributes of 'b'
  b@long_name = "multiple regression coefficients"
  b@model  = "Yest = b(0) + b(1)*X1 + b(2)*X2 + ... + b(M)*XM"

  b@N = N                   ; # of 'observations'
  b@NP = NP                 ; # of predictors
  b@M = M                   ; design matrix size
  b@bstd = bstd             ; standardized coefficients

  b@SST = SST               ; [1]
  b@SSE = SSE
  b@SSR = SSR

  b@MST = MST
  b@MSE = MSE
  b@MSE_dof = dof
  b@MSR = MSR

  b@RSE = RSE               ; [1]
  b@RSE_dof = N - M - 1     ; [1]

  b@F = Frat                ; [1]
  b@F_dof = (/NP, dof/)     ; [2]

  df1 = b@F_dof(0)
  df2 = b@F_dof(1)
  b@F_pval = ftest(MSR, df1 + 1, MSE, df2 + 1, 0) * 0.5   ; [1]

  b@r2 = r2                 ; [1]
  b@r = r                   ; [1]
  b@r2a = r2a               ; [1]

  b@fuv = fuv

  b@Yest = Yest             ; [NY]
  b@Yavg = Yavg             ; [1]
  b@Ystd = Ystd             ; [1]

  b@Xavg = Xavg             ; [1]
  b@Xstd = Xstd             ; [1]

  ; b@stderr = stderr         ; [M]
  ; b@tval = tval             ; [M]
  ; b@pval = pval             ; [M]

  leave_msg(scriptname, funcname)
  return(b)

end


; #############################################################################
undef("regress_stepwise")
function regress_stepwise(xarr_f[*][*]: numeric, \
                          yarr_f[*]: numeric)
;
; Arguments:
;     xarr_f: array containing all diagnostics for each model.
;     yarr_f: array with variable to which should be regressed.
;
; Return value:
;     coef: coefficients of the regression (beta vector).
;           coef@const: constant of regression (beta_0).
;           coef@idx_used_diags: array with indices of diagnostics used in the
;                                regression model.
;           coef@n_terms: number of terms in the regression model (= number of
;                         used diagnostics).
;           coef@YFIT: array with y values from linear regression.
;           coef@MCORRELATION: correlation coefficient of multicorrelation.
;           coef@correlation: correlation coefficients of individual
;                             correlations.
;
; Description:
;     Performs a multi diagnostic ensemble regression following Karpechko et
;     al., 2013.
;
; Caveats:
;     Needs attributes diag_script_info@wregr, diag_script_info@wexit and
;     optionally diag_scripts@p_value.
;
; Modification history:
;     20180807-schlund_manuel: ported to v2.0.
;     201303??-wenzel_sabrina: written.
;
local p_value, im, wregr, wexit, dim_diag, dim_dat, f_value_b, ssr0, sst, \
      used_diags, iiter, igd, jgd, nii, cvf, xarr, yarr

begin

  funcname = "regress_stepwise"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Required attributes
  req_atts = (/"wregr", "wexit"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  wregr = diag_script_info@wregr
  wexit = diag_script_info@wexit

  ; Optional attributes
  if (isatt(diag_script_info, "p_value")) then
    p_value = diag_script_info@p_value
  else
    p_value = 0.05
  end if

  ; Convert float arrays to double arrays
  xarr = todouble(xarr_f)
  yarr = todouble(yarr_f)
  copy_VarMeta(xarr_f, xarr)
  copy_VarMeta(yarr_f, yarr)

  ; Define some arrays
  dim_diag = dimsizes(xarr(:, 0))     ; Number of diagnostics
  dim_dat = dimsizes(xarr(0, :))      ; Number of datasets (i.e. models)
  diag_names = xarr&diagnostics       ; Diagnostic names

  used_diags = -1                     ; Array for regression variables
  ssr0 = 0.0d
  sst = sum((yarr - avg(yarr)) ^ 2)
  opt = False                         ; Options for multi linear regression

  ; Look-up table for critical value (= 1-p percentiles for a F-distribution
  ; with 1 degree of freedom in the nominator and varying (0-99 = index of the
  ; array) degrees of freedom in the denominator)
  if (p_value .eq. 0.1) then
    cvf = (/39.8635, 8.52632, 5.53832, 4.54477, 4.06042, 3.77595, 3.58943, \
           3.45792, 3.36030, 3.28502, 3.22520, 3.17655, 3.13621, 3.10221, \
           3.07319, 3.04811, 3.02623, 3.00698, 2.98990, 2.97466, 2.96096, \
           2.94858, 2.93736, 2.92712, 2.91774, 2.90913, 2.90119, 2.89385, \
           2.88703, 2.88069, 2.87478, 2.86926, 2.86409, 2.85922, 2.85465, \
           2.85035, 2.84628, 2.84244, 2.83880, 2.83535, 2.83208, 2.82897, \
           2.82600, 2.82317, 2.82047, 2.81790, 2.81544, 2.81308, 2.81082, \
           2.80866, 2.80658, 2.80458, 2.80266, 2.80082, 2.79904, 2.79733, \
           2.79568, 2.79409, 2.79255, 2.79107, 2.78963, 2.78825, 2.78691, \
           2.78560, 2.78434, 2.78312, 2.78194, 2.78080, 2.77968, 2.77860, \
           2.77755, 2.77654, 2.77555, 2.77458, 2.77364, 2.77272, 2.77184, \
           2.77098, 2.77013, 2.76931, 2.76851, 2.76773, 2.76696, 2.76622, \
           2.76550, 2.76479, 2.76410, 2.76342, 2.76276, 2.76211, 2.76149, \
           2.76087, 2.76026, 2.75968, 2.75909, 2.75853, 2.75797, 2.75743, \
           2.75689, 2.75637/)
  else if (p_value .eq. 0.05) then
    cvf = (/161.448, 18.5128, 10.1280, 7.70865, 6.60789, 5.98738, 5.59145, \
           5.31766, 5.11735, 4.96460, 4.84434, 4.74723, 4.66719, 4.60011, \
           4.54308, 4.49400, 4.45132, 4.41388, 4.38075, 4.35124, 4.32479, \
           4.30095, 4.27935, 4.25968, 4.24170, 4.22520, 4.21001, 4.19597, \
           4.18296, 4.17088, 4.15962, 4.14910, 4.13926, 4.13002, 4.12134, \
           4.11316, 4.10546, 4.09817, 4.09128, 4.08475, 4.07855, 4.07265, \
           4.06705, 4.06171, 4.05661, 4.05175, 4.04710, 4.04265, 4.03839, \
           4.03431, 4.03039, 4.02663, 4.02302, 4.01954, 4.01620, 4.01298, \
           4.00987, 4.00687, 4.00398, 4.00120, 3.99849, 3.99588, 3.99337, \
           3.99092, 3.98856, 3.98627, 3.98405, 3.98190, 3.97981, 3.97778, \
           3.97581, 3.97390, 3.97203, 3.97023, 3.96847, 3.96676, 3.96509, \
           3.96348, 3.96189, 3.96035, 3.95885, 3.95739, 3.95597, 3.95457, \
           3.95320, 3.95189, 3.95059, 3.94932, 3.94809, 3.94687, 3.94570, \
           3.94454, 3.94341, 3.94230, 3.94123, 3.94017, 3.93913, 3.93811, \
           3.93711, 3.93615/)
  else
    error_msg("f", scriptname, funcname, "p_value " + p_value + " not " + \
              "supported yet, function is hardcoded for p=0.05 and p=0.1")
  end if
  end if

  ; Start stepwise regression
  iiter = 1
  do while (wregr .or. wexit)
    ; .or.(f_prob .ge. f_value .and. f_prob .lt. f_value_b)
    log_debug("stepwise regression: iteration " + iiter)

    ; 1. Diagnostics to be tried in forward selection
    n_good_diag = dim_diag + 1 - dimsizes(used_diags)
    if (n_good_diag .le. 0) then
      log_debug("No diagnostics left to test: stopping iteration")
      break
    end if

    ; Array of untested "good" diagnostics
    good_diags = new(n_good_diag, integer)

    ; Loop over good diagnostics
    jgd  = 0
    do igd = 0, dim_diag - 1
      nii = num(used_diags .eq. igd)

      ; Only include the diagnostic if it is not used already
      if (nii .ne. 1) then
        good_diags(jgd) = igd
        jgd = jgd + 1
      end if
    end do
    log_debug("The following diagnostics are tried at forward step " + \
              iiter + ": " + diag_names(good_diags))

    ; 2. Trying the diagnostics in forward selection
    ssr1 = new(n_good_diag, double)
    sse1 = new(n_good_diag, double)
    F_val = new(n_good_diag, double)
    yfit := new((/n_good_diag, dim_dat/), double)

    ; Iterate over good diagnostics
    do igd = 0, n_good_diag - 1
      if (dimsizes(used_diags) .eq. 1) then
        test_diag = xarr(igd, :)
        rc = regline_stats(test_diag(:), yarr(:))
      else
        test_diag = new((/dimsizes(used_diags), dim_dat/), double)
        if (dimsizes(used_diags) .eq. 2) then
          test_diag(:dimsizes(used_diags) - 2, :) = \
            xarr(used_diags(1:):used_diags(1:), :)
        else
          test_diag(:dimsizes(used_diags) - 2, :) = \
            xarr(used_diags(1:), :)
        end if
        test_diag(dimsizes(used_diags) - 1, :) = xarr(good_diags(igd), :)
        rc  = reg_multlin_stats(yarr, transpose(test_diag), opt)
        MCORRELATION = rc@r
      end if
      sse1(igd) = rc@SSE
      ssr1(igd) = rc@SSR
      F_val(igd) = rc@F - ssr0
      yfit(igd, 0:dimsizes(rc@Yest) - 1) = rc@Yest
      delete([/test_diag, rc/])
    end do
    log_debug("ssr1 - ssr0 = " + (ssr1 - ssr0))

    ; 3. Defining which diagnostic is the best
    im = maxind(ssr1 - ssr0)
    log_debug("The following diagnostic has the largest increase in SSR: " + \
              diag_names(good_diags(im)) + ": " + (ssr1(im) - ssr0))

    ; 4. Testing whether it should be included
    dfe = (dim_dat - 1 - dimsizes(used_diags))
    f_value = (ssr1(im) - ssr0) / (sse1(im) / dfe)
    f_prob  = cvf(dfe - 1)
    delete(F_val)

    ; F-test
    if (f_value .le. f_prob) then
      log_debug("F-test shows that F-value is smaller than theoretical " + \
                "F-value: " + f_value + ", " + f_prob + ", " + ssr1(im) + \
                ", " + ssr0 + ", " + sse1(im) + ", " + dfe)
      log_debug("The diagnostic " + diag_names(good_diags(im)) + \
                " is not included into the regression model, skipping " + \
                "backward elimination and stopping iteration")

      ; If the diagnostic is not included then backward elimination step is
      ; not needed.
      wexit = False
      break   ; continue
    end if
    if (f_value .gt. f_prob) then
      log_debug("F-test shows that F-value is larger than theoretical " + \
                "F-value: " + f_value + ", " + f_prob + ", " + ssr1(im) + \
                ", " + ssr0 + ", " + sse1(im) + ", " + dfe)
      log_debug("The diagnostic " + diag_names(good_diags(im)) + \
                " is included into the regression model.")

      ; Include tested diagnostic (good_diags) to used diagnostics (used_diags)
      tmp = new((dimsizes(used_diags) + dimsizes(good_diags(im))), integer)
      tmp(0:dimsizes(used_diags) - 1) = used_diags
      tmp(dimsizes(used_diags):dimsizes(tmp) - 1) = good_diags(im)
      used_diags := tmp
      delete(tmp)
      ssr0 = ssr1(im)
      sse0 = sse1(im)
    end if

    ; Only one diagnostic used: break after more than two iterations
    if (dimsizes(used_diags) .eq. 2) then
      delete([/good_diags, ssr1, sse1/])
      if (iiter .ge. 2) then
        log_debug("Only one diagnotic used and more than 2 iterations: " + \
                  "stopping iteration")
        break
      else
        log_debug("Only one diagnostic used, skipping backward elimination")
        wregr = False
        iiter = iiter + 1
        continue
      end if
    end if

    ; 5. Trying the diagnostics in backward elimination
    wregr = False
    log_debug("The following diagnostics are tried at backward step " + \
              iiter + ": " + diag_names(used_diags(1:)))
    ssr1 := new(dimsizes(used_diags) - 1, double)
    F_val = new(dimsizes(used_diags) - 1, double)
    yfit := new((/dimsizes(used_diags) - 1, dim_dat/), double)

    ; Iterate over used diagnostics
    do igd = 1, dimsizes(used_diags) - 1
      ngood = ispan(1, dimsizes(used_diags) - 1, 1)
      good = ind(ngood .ne. igd)

      if (.not. all(ismissing(good))) then
        good = good + 1
        test_diag = xarr(used_diags(good), :)
        if (dimsizes(ngood) .eq. 1) then
          rc = regline_stats(test_diag, yarr)
        else
          rc = reg_multlin_stats_rap(yarr, transpose(test_diag), opt)
        end if
        sse2 = rc@SSE
        ssr1(igd - 1) = rc@SSR
        F_val(igd - 1) = ssr0 - rc@F
        yfit(igd - 1, 0:dimsizes(rc@Yest) - 1) = rc@Yest
        delete([/rc, test_diag/])
      else
        ssr1(igd - 1) = 0
        F_val(igd - 1) = (ssr0 - ssr1(im)) / (sse0 / dfe)
      end if
      delete([/good, ngood/])
    end do
    log_debug("ssr1 = " + ssr1)

    ; 6. Defining which diagnostic is the worst
    im := minind(ssr0 - ssr1)
    log_debug("The following diagnostic has the smallest decrease in SSR: " + \
              diag_names(used_diags(im+1)) + ": " + (ssr0 - ssr1(im)))

    ; 7. Testing whether it should be excluded
    dfe = (dim_dat - dimsizes(used_diags))
    f_value_b = (ssr0 - ssr1(im)) / (sse0 / dfe)
    f_prob = cvf(dfe - 1)
    delete(F_val)

    ; F-test
    if (f_value_b .le. f_prob) then
      log_debug("F-test shows that F-value is smaller than theoretical " + \
                "F-value: " + f_value_b + ", " + f_prob + ", " + ssr0 + \
                ", " + ssr1(im) + ", " + sse0 + ", " + dfe)
      log_debug("The diagnostic " + diag_names(good_diags(im + 1)) + \
                " is excluded from the regression model.")
      good = ind(used_diags .ne. used_diags(im + 1))
      tmp = used_diags(good)
      used_diags := tmp
      delete([/tmp, good/])
      ssr0 = ssr1(im)
    end if
    if (f_value_b .gt. f_prob) then
      log_debug("F-test shows that F-value is larger than theoretical " + \
                "F-value: " + f_value_b + ", " + f_prob + ", " + ssr0 + \
                ", " + ssr1(im) + ", " + sse0 + ", " + dfe)
      log_debug("The diagnostic " + diag_names(good_diags(im)) + \
                " is left in the regression model.")
    end if
    delete([/good_diags, ssr1, sse1/])

    ; if iiter equals 2 then stop
    ; if (iiter .eq. dim_diag - 2) then
    ;   ...
    ; end if

    ; if the regression model includes two terms then stop
    if (dimsizes(used_diags) .eq. 3)
      log_debug("Regression model includes two terms, stopping regression")
      break
    end if
    iiter = iiter + 1
  end do      ; Stepwise regression

  ; Final model
  if (.not. wregr) then
    log_debug("regress_stepwise correctly returned regression model")
    idx_used_diags = used_diags(1:)
    n_terms = dimsizes(used_diags) - 1
    finald = xarr(used_diags(1:), :)
    rc = reg_multlin_stats(yarr, transpose(finald), opt)
    MCORRELATION = rc@r2
    correlation = rc@pval
    coef = rc(1:)
    YFIT = rc@Yest
    const = rc(0)
  else
    error_msg("w", scriptname, funcname, "regress_stepwise returned " + \
              "empty regression model")
    const = avg(yarr)
    idx_used_diags := 0
    YFIT = new(dim_dat, double)
    YFIT = const
    rc = reg_multlin_stats(yarr, transpose(yfit), opt)
    MCORRELATION = rc@r
    correlation = escorc(yfit, yarr)
    coef = 0.0
    n_terms = 0
    coef@Yest = new(dim_dat, double)
    coef@Yest(:) = const
  end if

  ; Return elements
  coef@const = const
  coef@idx_used_diags = idx_used_diags
  coef@n_terms = n_terms
  coef@YFIT = YFIT
  coef@MCORRELATION = MCORRELATION
  coef@correlation = correlation

  ; Return data
  leave_msg(scriptname, funcname)
  return(coef)

end
