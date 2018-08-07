; #############################################################################
; mder.ncl
; #############################################################################
; AUXILIARY FUNCTIONS FOR MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; #############################################################################
; Please consider using of extending existing routines before adding new
; ones. Check the header of each routine for documentation.
;
; Contents:
;
;   * function diag_file_name: filename of calculated climatological diagnostic
;   * function f_cvf: computes cutoff value
;   * function regress_stepwise: does stewise regression on
;                                diagnostics
;   * function r_metric: calculates performance grade for each model
;
; #############################################################################

load "interface_scripts/interface.ncl"


; #############################################################################
undef("diag_file_name")
function diag_file_name(dataset_idx[1]: integer, \
                        opt[1]: logical)
;
; Arguments:
;     dataset_idx : index of the current model
;     opt : options, needs the attributes
;         @wdiag (diagnostic name)
;         @period_avg (year-, month- or season name)
;         @level (specify plev range)
;         @region (specify lat and/or lon range)
;
; Return value:
;     file_name : returns the file name of the calculated climatological
;                 diagnostic
; Description:
;     Provides a file name for the calculation results of diagsnostics.
;
; Modification history:
;   * 20180717_A-schl_ma: ported to v2.0
;   * 20140723_A-wenz_sa: written
;
local outname

begin

  funcname = "diag_file_name"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; File name
  outname = basename(systemfunc("basename " + \
                                input_file_info@filename(dataset_idx)))

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
undef("f_cvf")
function f_cvf(i[1]: numeric, \
               p[1]: numeric)
;
; Arguments:
;     i: Index.
;     p: A non-negative scalar, in the interval [0.0, 1.0], of
;         type float or double that specifies the probability of
;         occurance or success.
;
; Return value:
;     Value of specified f-distribution.
;
;     DFN: A positive scalar of type integer, float or double that
;         specifies the degrees of freedom of the F distribution
;         numerator.
;     DFD: A positive scalar of type integer, float or double that
;         specifies the degrees of freedom of the F distribution
;         denominator.
;
; Description:
;     Quick and dirty version:
;     This function computes the cutoff value (v) such that:
;         Probability(X > v) = p
;     where X is a random variable from the F distribution with
;     (dfn) and (dfd) degrees of freedom.
;
; Caveats:
;     FIXME (needs to be more general).
;
; Modification history:
;   * 20180807_A-schl_ma: ported to v2.0
;   * 201502??_A-wenz_sa: written
;
local tmp

begin

  funcname = "f_cvf"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  if (p .eq. 0.05) then
    tmp = (/161.448, 18.5128, 10.1280, 7.70865, 6.60789, 5.98738, 5.59145, \
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
  end if

  cvf = tmp(i)
  leave_msg(scriptname, funcname)
  return(cvf)

end


; #############################################################################
undef("regress_stepwise")
function regress_stepwise(xarr[*][*]:numeric, \
                          yarr[*]:numeric)
;
; Arguments:
;     xarr: array containing all diagnostics for each model
;     yarr: array with variable to which should be regressed
;
; Return value:
;     const = constant of regression (beta_0)
;     terms = array with diagn used for regression model (beta vector)
;     nterms = dimensionsize of terms
;     YFIT = array with y values from linear regression
;     MCORRELATION = correlation coefficient of multicorrelation
;     correlation  = correlation coefficient of individual correlations
;
; Description:
;     Performs a multi diagnostic ensemble regression following Karpechko et
;     al., 2013.
;
; Modification history:
;   * 20180807_A-schl_ma: ported to v2.0
;   * 201303??-A_schl_ma: written
;
local pvalue, im, wregr, wexit, nterms_in, npoints, f_value_b, ssr0,\
      sst, iiud, iiter, igd, jgd, nii, cvf

begin

  funcname = "regress_stepwise"
  scriptname = "diag_scripts/shared/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Local attributes for function called from diag_script_info
;   if (isatt(diag_script_info, "wregr") .and. \
;       isatt(diag_script_info, "wexit")) then
;     wregr = diag_script_info@wregr      ; 'no' for exit criteria
;     wexit = diag_script_info@wexit      ; 'false'for exit criteria
;   else
;     error_msg("f", scriptname, funcname, "Exit criteria need to be " + \
;               "defined in 'wregr' and 'wexit' in recipe")
;   end if

;   if (isatt(diag_script_info,"pvalue"))then
;     pvalue = diag_script_info@pvalue
;   else
;     pvalue = 0.05
;   end if

;   ;;Define some arrays
;   nterms_in = dimsizes(xarr(:,0))     ;;num of diagnostics
;   npoints   = dimsizes(xarr(0,:))     ;;num of models
;   wterm = xarr&diagnostics            ;;array with diagnostic names
;   nterms = 0
;   f_value_b = 999.;new(1,float);
;   ;f_prob = 999.;new(1,float);
;   ;f_value = 999.;new(1,float);

;   iiud = -1                           ;;Array for regression variables (=num of 'used diagnostics')
;   ssr0 = 0.
;   sst = sum((yarr-avg(yarr))^2)
;   opt = True
;   opt = False                         ;;Options for multi linear regression


;   ; Look-up table for critical value
;   ; 0.1 prob level
; ;  cvf = (/39.8635,8.52632,5.53832,4.54477,4.06042,3.77595,3.58943,3.45792,\
; ;             3.36030,3.28502,3.22520,3.17655,3.13621,3.10221,3.07319,3.04811,\
; ;             3.02623,3.00698,2.98990,2.97466,2.96096,2.94858,2.93736,2.92712,\
; ;             2.91774,2.90913,2.90119,2.89385,2.88703,2.88069,2.87478,2.86926,\
; ;             2.86409,2.85922,2.85465,2.85035,2.84628,2.84244,2.83880,2.83535,\
; ;             2.83208,2.82897,2.82600,2.82317,2.82047,2.81790,2.81544,2.81308,\
; ;             2.81082,2.80866,2.80658,2.80458,2.80266,2.80082,2.79904,2.79733,\
; ;             2.79568,2.79409,2.79255,2.79107,2.78963,2.78825,2.78691,2.78560,\
; ;             2.78434,2.78312,2.78194,2.78080,2.77968,2.77860,2.77755,2.77654,\
; ;             2.77555,2.77458,2.77364,2.77272,2.77184,2.77098,2.77013,2.76931,\
; ;             2.76851,2.76773,2.76696,2.76622,2.76550,2.76479,2.76410,2.76342,\
; ;             2.76276,2.76211,2.76149,2.76087,2.76026,2.75968,2.75909,2.75853,\
; ;             2.75797,2.75743,2.75689,2.75637/)

;   ; 0.05 prob level
;   cvf = (/161.448,18.5128,10.1280,7.70865,6.60789,5.98738,5.59145,5.31766,\
;           5.11735,4.96460,4.84434,4.74723,4.66719,4.60011,4.54308,4.49400,\
;           4.45132,4.41388,4.38075,4.35124,4.32479,4.30095,4.27935,4.25968,\
;           4.24170,4.22520,4.21001,4.19597,4.18296,4.17088,4.15962,4.14910,\
;           4.13926,4.13002,4.12134,4.11316,4.10546,4.09817,4.09128,4.08475,\
;           4.07855,4.07265,4.06705,4.06171,4.05661,4.05175,4.04710,4.04265,\
;           4.03839,4.03431,4.03039,4.02663,4.02302,4.01954,4.01620,4.01298,\
;           4.00987,4.00687,4.00398,4.00120,3.99849,3.99588,3.99337,3.99092,\
;           3.98856,3.98627,3.98405,3.98190,3.97981,3.97778,3.97581,3.97390,\
;           3.97203,3.97023,3.96847,3.96676,3.96509,3.96348,3.96189,3.96035,\
;           3.95885,3.95739,3.95597,3.95457,3.95320,3.95189,3.95059,3.94932,\
;           3.94809,3.94687,3.94570,3.94454,3.94341,3.94230,3.94123,3.94017,\
;           3.93913,3.93811,3.93711,3.93615/)

;   iiter = 1

;   ;;Start stepwise regression::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
;   do while (wregr.eq."no".or.wexit.eq."false");.or.(f_prob.ge.f_value.AND.f_prob.lt.f_value_b)

;     ;;1. Diagnostics to be tried in forward selection:
;     ngd  = nterms_in+1-dimsizes(iiud) ;;num of untested 'good diagnostics'
;     if (ngd.le.0) then
;       break
;     end if
;     iigd = new(ngd,integer)           ;;array of untested 'good diagnostics'
;     jgd  = 0

;     do igd = 0,nterms_in-1            ;;Loop over good diagnostics
;       nii = num(iiud .eq. igd)
;       ;;Only include the diagnostic if it is not used already
;       if (nii.ne.1) then
;         iigd(jgd) = igd
;         jgd = jgd+1
;       end if
;     end do

;     info_output("The following diagnostics are tried at "+iiter+" forward step:",verbosity, 2)
;     if (verbosity.ge.2) then
;       print("info:  "+wterm(iigd))
;     end if

;     ;;2. Trying the diagnostics in forward selection:
;     ssr1  = new(ngd,float)
;     sse1  = new(ngd,float)
;     F_val = new(ngd,float)

;     if (isvar("yfit")) then
;       delete(yfit)
;     end if
;     yfit  = new((/ngd,npoints/),float)

;     do igd = 0,ngd-1
;       if (dimsizes(iiud).eq.1) then
;         testd = xarr(igd,:)
;         rc  = regline_stats(testd(:),yarr(:))
;       else
;         testd = new((/dimsizes(iiud),npoints/),float)
;         if (dimsizes(iiud).eq.2) then
;           testd(:dimsizes(iiud)-2,:) = xarr(iiud(1:):iiud(1:),:)
;         else
;           testd(:dimsizes(iiud)-2,:) = xarr(iiud(1:),:)
;         end if
;         testd(dimsizes(iiud)-1,:) = xarr(iigd(igd),:)
;         rc  = reg_multlin_stats(yarr,transpose(testd),opt)
;         MCORRELATION = rc@r
;       end if
;       sse1(igd)   = rc@SSE
;       ssr1(igd)   = rc@SSR
;       F_val(igd)  = rc@F-ssr0
;       yfit(igd,0:dimsizes(rc@Yest)-1) = rc@Yest
;       delete(testd)
;       delete(rc)
;     end do

;     if (verbosity.ge.2) then
;       print("info:  ssr1-ssr0 = "+(ssr1 - ssr0))
;     end if

;     ;;3. Defining which diagnostic is the best:
;     im = maxind(ssr1-ssr0)

;     info_output("The following diagnostic has the largest increase in SSR:", verbosity, 2)
;     if (verbosity.ge.2) then
;       print("info:  "+wterm(iigd(im))+": "+(ssr1(im)-ssr0))
;     end if

;     ;;4. Testing whether it should be included
;     ;;F-test
;     dfe = (npoints-1-(1+dimsizes(iiud)-1))
;     f_value = (ssr1(im)-ssr0)/(sse1(im)/dfe)
;     f_prob  = cvf(dfe-1)
;     delete(F_val)

;     if (f_value .gt. f_prob) then
;       info_output("F-test shows that F-value is larger than theoretical F-value:",verbosity, 2)
;       info_output(f_value+" , "+f_prob+" , "+ssr1(im)+" , "+ssr0+" , "+sse1(im)+" , "+dfe,verbosity, 2)
;       info_output("The diagnostic "+wterm(iigd(im))+" is included into the model.",verbosity, 2)

;       ;;include tested diagnostic (iigd) to 'used diagnostics' (iiud)
;       tmp = new((dimsizes(iiud)+dimsizes(iigd(im))),integer)
;       tmp(0:dimsizes(iiud)-1) = iiud
;       tmp(dimsizes(iiud):dimsizes(tmp)-1) = iigd(im)
;       delete(iiud)
;       iiud = tmp
;       delete(tmp)

;       ssr0 = ssr1(im)
;       sse0 = sse1(im)
;     end if

;     if (f_value .lt. f_prob) then
;       info_output("F-test shows that F-value is smaller than theoretical F-value:",verbosity, 2)
;       info_output(f_value+" , "+f_prob+" , "+ssr1(im)+" , "+ssr0+" , "+sse1(im)+" , "+dfe,verbosity, 2)
;       info_output("The diagnostic "+wterm(iigd(im))+" is not included into the model.",verbosity, 2)

;       ;;IF the diagnostic is not included then backward elimination step is not needed.
;       wexit = "true"
;       break;continue                         ;GOTO ,toexit
;     end if

;     ;;break if only one diagnostic is used
;     if (dimsizes(iiud).eq.2) then
;       delete(iigd)
;       delete(ssr1)
;       delete(sse1)
;       if (iiter.ge.2) then
;         break
;       else
;         info_output("CONTINUE",verbosity, 2)
;         wregr = "yes"
;         iiter = iiter+1
;         continue
;       end if
;     end if

;     ;;5. Trying the diagnostics in backward elimination:
;     wregr = "yes"
;     info_output("The following diagnostics are tried at "+iiter+" backward step:",verbosity, 2)
;     if (verbosity.ge.2) then
;       print("info:  "+wterm(iiud(1:)))
;     end if

;     delete(ssr1)
;     delete(yfit)
;     ssr1  = new(dimsizes(iiud)-1,float)
;     F_val = new(dimsizes(iiud)-1,float)
;     yfit  = new((/dimsizes(iiud)-1,npoints/),float)

;     do igd = 1,dimsizes(iiud)-1
;       ngood = ispan(1,dimsizes(iiud)-1,1)       ;;Num of 'good diagnostics'
;       good = ind(ngood.ne.igd)

;       if (.not.all(ismissing(good))) then
;         good  = good+1
;         testd = xarr(iiud(good),:)

;         if (dimsizes(ngood).eq.1) then
;           rc  = regline_stats(testd,yarr)
;         else
;           rc = reg_multlin_stats_rap(yarr,transpose(testd),opt)
;         end if

;         sse2         = rc@SSE
;         ssr1(igd-1)  = rc@SSR
;         F_val(igd-1) = ssr0-rc@F
;         yfit(igd-1,0:dimsizes(rc@Yest)-1)  = rc@Yest
;         delete(rc)
;         delete(testd)
;       else
;         ssr1(igd-1)  = 0
;         F_val(igd-1) = (ssr0-ssr1(im))/(sse0/dfe)
;       end if

;       delete(good)
;       delete(ngood)
;     end do

;     if (verbosity.ge.2) then
;       print("info:  ssr1 = "+ssr1)
;     end if

;     ;;6. Defining which diagnostic is the worst:
;     delete(im)
;     im = minind(ssr0-ssr1)

;     info_output("The following diagnostic has the smallest decrease in SSR:",verbosity, 2)
;     if (verbosity.ge.2) then
;       print("info:  "+wterm(iiud(im+1))+": "+(ssr0-ssr1(im)))
;     end if

;     ;;7. Testing whether it should be excluded
;     ;;F-test
;     dfe       = (npoints-1-(dimsizes(iiud)-1))
;     f_value_b = (ssr0-ssr1(im))/(sse0/dfe)
;     f_prob    = cvf(dfe-1)
;     delete(F_val)

;     if (f_value_b .le. f_prob) then
;       info_output("F-test shows that F-value is smaller than theoretical F-value:",verbosity, 2)
;       info_output(f_value_b+" , "+f_prob+" , "+ssr0+" , "+ssr1(im)+" , "+sse0+" , "+dfe,verbosity, 2)
;       info_output("The diagnostic "+wterm(iiud(im+1))+" is excluded from the model.",verbosity, 2)

;       good = ind(iiud.ne.iiud(im+1))
;       tmp = iiud(good)
;       delete(iiud)
;       iiud = tmp
;       delete(tmp)
;       delete(good)
;       ssr0 = ssr1(im)
;     end if

;     if (f_value_b .gt. f_prob) then
;       info_output("F-test shows that F-value is larger than theoretical F-value:",verbosity, 2)
;       info_output(f_value_b+" , "+f_prob+" , "+ssr0+" , "+ssr1(im)+" , "+sse0+" , "+dfe,verbosity, 2)
;       info_output("The diagnostic "+wterm(iiud(im+1))+" is left in the model.",verbosity, 2)
;     end if

;     delete(iigd)
;     delete(ssr1)
;     delete(sse1)

;     ;;if iiter eq 2 then stop
;     ;if (iiter.eq.nterms_in-2) then

;     ;;if the regression model includes two terms then stop
;     if (dimsizes(iiud).eq.3)
;       ;info_output("No diagnostic matches the criteria",verbosity, 1)
;       break
;     end if

;     iiter = iiter+1

;   end do ;;End stepwise regression:::::

;   ;;Final model
;   if (wregr .eq. "yes") then
;     terms   = iiud(1:)
;     nterms  = dimsizes(iiud)-1
;     finald  = xarr(iiud(1:),:)
;     rc      = reg_multlin_stats(yarr,transpose(finald),opt)
;     MCORRELATION = rc@r2
;     correlation  = rc@pval
;     coef  = rc(1:)
;     YFIT  = rc@Yest
;     const = rc(0)
;   else
;     const = avg(yarr)
;     terms = 0.
;     YFIT  = new(npoints,float)
;     YFIT  = const
;     rc    = reg_multlin_stats(yarr,transpose(yfit),opt)
;     MCORRELATION = rc@r
;     correlation  = escorc(yfit,yarr)
;     coef   = 0.
;     nterms = 0
;     coef@Yest = new(npoints,float)
;     coef@Yest(:) = const
;   end if

;   ;;Return elements
;   coef@const        = const
;   coef@terms        = terms
;   coef@nterms       = nterms
;   coef@YFIT         = YFIT
;   coef@MCORRELATION = MCORRELATION
;   coef@correlation  = correlation

;   ; Return data
;   leave_msg(scriptname, funcname)
;   return(coef)

end


; ; #############################################################################
; undef("r_metric")
; function r_metric(mu_m[*]:numeric,\
;                   mu_r[*]:numeric,\
;                  sigma[*]:numeric)
; ;;
; ;; Author: Sabrina Wenzel (DLR, Germany)
; ;;
; ;; Description:
; ;;    * calculates a grade for each models according to its performance #
; ;;      compared with observations, following Waugh and Eyring (2008)
; ;;
; ;; Return value:
; ;;     metric = array with gradees for each model
; ;;
; ;; Arguments:
; ;;    mu_m  = array with climatological mean diagnostics for all models
; ;;    mu_r  = array with climatological mean diagnostics from OBS
; ;;    sigma = observational uncertainty (stddev(OBS))
; ;;
; ;; Caveats:
; ;;
; ;; Modification history:
; ;;    * 201411-A_WeSa: written by Sabrina Wenzel
; local neg,mu_m,mu_r,sigma
; begin

;   funcname = "r_metric"
;   scriptname = "diag_scripts/shared/mder.ncl"
;   enter_msg(scriptname, funcname)

;   if (.not.ismissing(sigma)) then
;     if (sigma.ne.0.) then
;       metric = 1.-abs(mu_m-mu_r)/(3.*sigma)
;     else
;       metric = new(dimsizes(mu_m),typeof(mu_m))
;     end if
;   else
;     metric = new(dimsizes(mu_m),typeof(mu_m))
;   end if

;   neg = ind(metric.lt.0.)
;   if (.not.all(ismissing(neg))) then
;     metric(neg) = 0.
;   end if

;   leave_msg(scriptname, funcname)
;   return(metric)

; end
