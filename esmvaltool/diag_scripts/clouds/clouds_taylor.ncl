; ############################################################################
; CLOUDS_TAYLOR
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates the performance of models in reproducing 2-d fields of annual
;   mean or seasonal (DFJ, MAM, JJA, SON) mean cloud properties
;   (Taylor diagramms). The code is based on
;     1) perfmetrics_main.ncl
;     2) perfmetrics_taylor.ncl
;     3) perfmetrics_taylor_collect.ncl
;   originally written by Franziska Frank (DLR, Germany).
;   Note: This code requires that all data are provided on the same grid.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   embracelegend:   - False (default) = include legend in plot, max.
;                      2 columns with dataset names in legend
;                    - True = write extra file with legend, max. 7 dataset
;                      names per column in legend, alternative observational
;                      dataset(s) will be plotted as a red star and labeled
;                      "altern. ref. dataset" in legend (only if dataset is
;                      of class "OBS")
;   estimate_obs_uncertainty: - True = estimate observational uncertainties
;                      from mean values (assuming fractions of obs. RMSE from
;                      documentation of the obs data); only available for
;                      "CERES-EBAF", "MODIS", "MODIS-L3"
;                    - False = do not estimate obs. uncertainties from mean
;                      values
;   filename_add:    legacy feature: arbitrary string to be added to all
;                    filenames of plots and netcdf output produced
;                    (default = "")
;   mask_ts_sea_ice: - True = mask T < 272 K as sea ice (only for
;                      variable "ts")
;                    - False = no additional grid cells masked for
;                      variable "ts"
;   styleset:        "CMIP5", "DEFAULT"
;                    (if not set, CLOUDS_TAYLOR will create a color table and
;                    symbols for plotting)
;   timemean:        time averaging
;                    - annualclim (default) = 1 plot annual mean
;                    - seasonalclim = 4 plots (DJF, MAM, JJA, SON)
;   valid_fraction:  used for creating sea ice mask (mask_ts_sea_ice = true):
;                    fraction of valid time steps required to mask grid cell
;                    as valid data
;
; Required variable_info attributes (variable specific)
;   reference_dataset: name of reference data set
;
; Optional variable attributes (variable specific)
;   none
;
; Caveats
;   KNOWN LIMITATIONS
;     1) only 2-dim variables are currently supported
;     2) observational uncertainties are regridded like standard variables
;     3) for derived variables (e.g. SW_CRE), also the original variables have
;        to be specified in the nameslist for the reference model (for the
;         reference dataset only) in order to estimate the observational
;         uncertainty (if estimate_obs_uncertainty = True)
;     4) the variable to plot must be the *first* variable in the namelist,
;        the observational uncertainty (or the variables from which a derived
;        variable has been calculated from) the *second* (third, ...)
;
; Modification history
;   20190221-lauer_axel: added provenance to output (v2.0)
;   20181120-lauer_axel: adapted code to multi-variable capable framework
;   20180923-lauer_axel: added writing of results to netcdf
;   20180611-lauer_axel: code rewritten for ESMValTool v2.0
;   20170620-lauer_axel: added tags for reporting
;   20161104-lauer_axel: changed calculation of correlation and standard
;                        deviation ratio from unweighted to grid area weighted
;   20160901-lauer_axel: added regridding option 1 deg x 1 deg
;   20160818-lauer_axel: added support for variable attribute "ref_model"
;   20160404-lauer_axel: added optional drawing of observational
;                        uncertainties
;   20151027-lauer_axel: moved call to 'write_references' to the beginning
;                        of the code
;   20150415-lauer_axel: written.
;
; ###########################################################################


load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/taylor_plot.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  set_default_att(diag_script_info, "embracelegend", False)
  set_default_att(diag_script_info, "estimate_obs_uncertainty", False)
  set_default_att(diag_script_info, "mask_ts_sea_ice", False)
  set_default_att(diag_script_info, "timemean", "annualclim")
  set_default_att(diag_script_info, "valid_fraction", 0.5)

  estimate_obs_uncertainty = diag_script_info@estimate_obs_uncertainty

  variables = metadata_att_as_array(variable_info, "short_name")
  numvars = dimsizes(variables)

  ; find "main" variable and if present uncertainty estimates or auxiliary
  ; variables that can optionally be used to estimate the uncertainties

  mainvarind = 0  ; default = main variable is first (and only) variable
  errvarind = -1  ; default = no error estimate available

  if (numvars.gt.1) then
    ; uncertainty estimates are identified by "err" in their short_name
    do i = 0, numvars - 1
      if (isStrSubset(variables(i), "err")) then
        errvarind = i
        break
      end if
    end do
    ; now find the main variable:
    ; it is assumed that the main variable is either the first non-error
    ; variable or the first "derived" variable
    idx = ind(ispan(0, numvars - 1, 1) .ne. errvarind)
    mainvarind = idx(0)  ; first non-error variable
    ; now check for possibly derived variables
    do ii = 1, dimsizes(idx) - 1
      i = idx(ii)
      ; the "derived" attribute is stored in input_file_info only
      do j = 0, ListCount(input_file_info) - 1
        if (input_file_info[j]@short_name .eq. variables(i)) then
          if isatt(input_file_info[j], "derive") then
            if (input_file_info[j]@derive) then
              mainvarind = i
              break
            end if
          end if
        end if
      end do
    end do
  else
  end if

  var0 = variable_info[mainvarind]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)
  if (isatt(variable_info[mainvarind], "reference_dataset")) then
    refname = variable_info[mainvarind]@reference_dataset
  end if
  if (isatt(variable_info[mainvarind], "alternative_dataset")) then
    refname2 = variable_info[mainvarind]@alternative_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  infiles = metadata_att_as_array(input_file_info, "filename")

  ; check for reference model definition

  if (.not.isvar("refname")) then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
  end if

  ; print info on variable(s) to log file

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT)
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info("main variable: " + variables(mainvarind))
  if (errvarind .ge. 0) then
    log_info("error variable: " + variables(errvarind))
  else
    log_info("error variable: none")
  end if
  tmp = ispan(0, numvars - 1, 1)
  auxind = ind(tmp .ne. mainvarind .and. tmp .ne. errvarind)
  if(.not. all(ismissing(auxind))) then
    if (estimate_obs_uncertainty) then
      log_info("aux variable: " + variables(auxind))
    end if
  else
    if (estimate_obs_uncertainty) then
      log_info("no aux variable(s): cannot estimate obs uncertainty")
    end if
  end if
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  delete(tmp)
  delete(auxind)

  ; optional uncertainty estimates are only supported for the reference
  ; dataset, i.e. max 1 dataset

  if (estimate_obs_uncertainty) then
    auxind = ind(ispan(0, numvars - 1, 1) .ne. mainvarind)
    if (.not. all(ismissing(auxind))) then
      do i = 0, dimsizes(auxind) - 1
        tmpinfo = select_metadata_by_name(input_file_info, variables(i))
        items = ListCount(tmpinfo)
        if (items .gt. 1) then
          error_msg("f", DIAG_SCRIPT, "", "uncertainty estimate only " + \
                    "supported for 1 dataset (i.e. reference dataset)")
        end if
        tmpname = metadata_att_as_array(tmpinfo, "dataset")
        if (tmpname .ne. refname) then
          error_msg("f", DIAG_SCRIPT, "", "uncertainty estimate only " + \
                    "supported for the reference dataset (" + refname + ")")
        end if
        delete(tmpinfo)
      end do
    end if
    delete(auxind)
  end if

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  timemean = diag_script_info@timemean
  numseas  = 1          ; default
  season   = (/"annual"/)

  if (timemean.eq."seasonalclim") then
    numseas = 4
    delete(season)
    season = (/"DJF", "MAM", "JJA", "SON"/)
  end if

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  if (config_user_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if

end

begin
  ; ========================================================================
  ; ========================== initialization ==============================
  ; ========================================================================

  ; get reference model

  ref_ind = ind(names .eq. refname)
  if (ismissing(ref_ind)) then
    error_msg("f", DIAG_SCRIPT, "", "reference dataset (" \
              + refname + ") is missing")
  end if

  ; get multi-model mean index (if present)

  mm_ind = ind(names .eq. "MultiModelMean")

  if (ismissing(mm_ind)) then
    mm_ind = -1
  end if

  mask_ts_sea_ice = diag_script_info@mask_ts_sea_ice

  if (isatt(diag_script_info, "filename_add")) then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  embracelegend = diag_script_info@embracelegend

  ; ========================================================================
  ; ============================ statistics ================================
  ; ========================================================================

  if ((var0 .eq. "ts") .and. (mask_ts_sea_ice))
    ; Create a missing value mask (from reference data)
    ; The dim_avg_wgt function with option limit is used on the time
    ; coordinate: in this way each grid box containing more than
    ; valid_fraction of missing values along the time series is
    ; set to missing.

    data_temp = read_data(info0[ref_ind])

    ww = data_temp&time
    ww = 1.
    limit  = toint(diag_script_info@valid_fraction * \
                   dimsizes(data_temp&time))

    if (limit.lt.1) then
      limit = 1
    end if

    ; temperatures below 272 K are sea ice ---> mask
    data_temp = where(data_temp.lt.272.0, data_temp@_FillValue, data_temp)

    ; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1
    if (dimsizes(dimsizes(data_temp)).eq.1) then
      tmp = dim_avg_wgt_Wrap(data_temp, ww, limit)
    else
      tmp = dim_avg_wgt_n_Wrap(data_temp, ww, limit, 0)
    end if

    global_mask = where(ismissing(tmp), tmp, 0.)

    delete(tmp)
    delete(ww)
    delete(data_temp)
  end if

  ;  ; debug output
  ;  system ("rm debug.nc")
  ;  debugfile = addfile("debug.nc","c")
  ;  debugfile->mask = global_mask

  ; read data and calculate time average

  ierridx = 0

  do ivar = 0, numvars - 1

    ; ignore auxiliary variables if estimate_obs_uncertainty = False

    if (.not.estimate_obs_uncertainty .and. ivar .ne. mainvarind .and. \
        ivar .ne. errvarind) then
      continue
    end if

    ; process all models for first variable
    if (ivar.eq.mainvarind) then
      modelf = 0
      modell = dim_MOD - 1
    ; process error estimates / auxiliary variables for error estimate
    ; (reference model only)
    else
      modelf = 0
      modell = 0
    end if

    info = select_metadata_by_name(input_file_info, variables(ivar))

    do imod = modelf, modell
      data_temp = read_data(info[imod])

      ; The uncertainty estimates are assumed to be given as
      ; 1-sigma of normally distributed error estimates
      ; Note: uncertainty estimates must have "err" in their
      ;       variable name, otherwise variables are assumed
      ;       to be fields used for calculating derived variables
      ;       (in this case, fields are not squared before averaging)
      ; ---> square sigma before averaging over time, then
      ;      calculate square-root to get average sigma

      if (ivar.eq.errvarind) then
        data_temp = data_temp * data_temp
      end if

      data1 = time_operations(data_temp, -1, -1, "average", \
                              timemean, True)

      if (ivar.eq.errvarind) then
        data1 = sqrt(data1)
      end if

      if (isdefined("global_mask")) then
        if (numseas.eq.1) then
          data1 = data1 + global_mask
        else
          do is = 0, numseas - 1
            data1(is, :, :) = data1(is, :, :) + global_mask
          end do
        end if
      end if

      if (ivar.eq.mainvarind) then
        if (.not.isdefined("data")) then
          dim_data = array_append_record((/dim_MOD/), dimsizes(data1), 0)
          data = new(dim_data, float)
          data!0 = "model"
          data&model = names
          if (numseas.eq.1) then
            copy_VarCoords(data1, data(imod, :, :))
          else
            copy_VarCoords(data1, data(imod, :, :, :))
          end if
        end if
        dim_data = dimsizes(data)
        rank = dimsizes(dim_data)
        if (numseas.eq.1) then
          data(imod, :, :) = data1
        else
          data(imod, :, :, :) = data1
        end if
        delete(data_temp)
        delete(data1)
      else
        if (.not.isdefined("err")) then
          dim_data = array_append_record((/numvars - 1/), dimsizes(data1), 0)
          err = new(dim_data, float)
          err!0 = "var"
          err&var = variables(ind(variables .ne. var0))
        end if
        if (numseas.eq.1) then
          err(ierridx, :, :) = data1
        else
          err(ierridx, :, :, :) = data1
        end if
      end if

      ; debug output

    ; debugfile->$input_file_info@dataset(imod)$ = data1
    end do  ; imod loop
    if (ivar .ne. mainvarind) then
      ierridx = ierridx + 1
    end if
  end do  ; ivar loop

  ; number of used reference data sets

  dim_REF = 1
  if (isvar("refname2")) then
    dim_REF = dim_REF + 1
  end if

  ; define result variable

  val = new((/dim_MOD - dim_REF, 2, numseas/), float)
  val!0 = "models"
  val!1 = "statistic"
  val!2 = "time"
  val&models = names(ind((names .ne. refname)))
  val&statistic = (/"stddev_ratio", "correlation"/)
  val&time = ispan(0, numseas - 1, 1)

  ; extract data

  do is = 0, numseas - 1

    if (numseas.eq.1) then
      ref = data(ref_ind, :, :)
    else
      ref = data(ref_ind, is, :, :)
    end if

    ; loop over models
    ; (second model index, for models only, skip reference)

    mod_i = 0
    do imod = 0, dim_MOD - 1

      ; skip reference

      if (imod.eq.ref_ind) then
        log_info("Ref_dataset = " + data&model(imod))
        continue
      end if

      if (numseas.eq.1) then
        var = data(imod, :, :)
      else
        var = data(imod, is, :, :)
      end if

      ; calculate rate of standard deviations

      ; val(mod_i, 0, is) = calculate_metric(var, ref, "stddev_ratio_nowgt")
      val(mod_i, 0, is) = calculate_metric(var, ref, "stddev_ratio")

      ; calculate pattern correlation

      ; val(mod_i, 1, is) = calculate_metric(var, ref, "correlation_nowgt")
      val(mod_i, 1, is) = calculate_metric(var, ref, "correlation")

      delete(var)
      mod_i = mod_i + 1

    end do  ; loop over datasets
  end do  ; loop over seasons

  ; attach attributes to the results

  val@title = "taylor metrics"
  val@long_name = "Taylor diagram"
  val@diag_script = (/DIAG_SCRIPT/)
  val@var = "var"

  ; ------------------- calculate/estimate RMS of observations ---------------

  rmsobs = new((/numseas/), float)

  ; optionally, if there are no variable(s) available for the error estimate,
  ; the reference data can be used to create one for selected data sets

  if (estimate_obs_uncertainty) then
    if (any(refname.eq.(/"CERES-EBAF", "MODIS", "MODIS-L3"/))) then
      if (numvars.eq.1) then  ; no variable(s) for error estimate
        numvars = 2
        variables0 = (/variables(0), variables(0)/)
        delete(variables)
        variables = variables0
        dim_data = dimsizes(data)
        dim_data(0) = numvars - 1
        err = new(dim_data, float)
        if (numseas.eq.1) then
          err(0, :, :) = data(ref_ind, :, :)
        else
          err(0, :, :, :) = data(ref_ind, :, :, :)
        end if
      end if
    end if
  end if  ; if estimate_obs_uncertainty

  do is = 0, numseas - 1
    if (numseas.eq.1) then
      ref = data(ref_ind, :, :)
    else
      ref = data(ref_ind, is, :, :)
    end if

    ; This code is equivalent to the function "calculate_metric"
    ; (statistics.ncl)
    weights = map_area(ref&lat, ref&lon)
    ; optional: no weights --> reset weights
    ; wgt1d = 1.0

    ; convert to 1-D arrays
    ref1d = ndtooned(ref)
    wgt1d = ndtooned(weights)
    wgt1d = wgt1d / dim_avg_n_Wrap(wgt1d, 0)

    rmsobs(is) = 0.0

    if (errvarind .ge. 0 .or. estimate_obs_uncertainty) then
      ierridx = 0
      do iv = 0, numvars - 1
        if (iv .eq. mainvarind) then
          continue
        end if

        if (numseas.eq.1) then
          var = err(ierridx, :, :)
        else
          var = err(ierridx, is, :, :)
        end if

        ; If the variable name contains "err", the variable is assumed to
        ; contain an estimate of the observational uncertainty (e.g. UWisc,
        ; GPCP-SG, and SSMI-MERIS provide error estimates, MODIS-L3-C6 provides
        ; error estimates for lwp and iwp but not for clt). In this case, there
        ; is nothing to do here.
        ; If the variable name does not contain "err", the variable is assumed
        ; to be a field (e.g. rsut, rsutcs) used to calculate a derived
        ; variable (e.g. swcre). In this case, the observational uncertainty
        ; can be estimated if the fractional average error is known
        ; (e.g. CERES-EBAF, MODIS).

        if (errvarind .ge. 0) then
          ; variable contains error estimate provided by data set
          ; ---> nothing to do
          factor = 1.0
        else
          ; variable(s) contain(s) fields used to derive the variable plotted
          ; --> observational uncertainty might be estimated from these fields
          ;    if the average error(s) of these fields are known

          ; by default, uncertainty cannot be estimated

          factor = 0.0

          ; for CERES-EBAF, we estimate the errors from the original variables
          ; (rsut+rsutcs / rlut+rlutcs) assuming variable specific factors

          if (refname.eq."CERES-EBAF") then
            if (any(variables(iv).eq.(/"rsut", "rsutcs"/))) then
              ; From the "Technical Note on CERES EBAF Ed2.7":
              ; rsutcs: 1% ~ 0.5 W/m2
              ;      estimate total error = 2.6 W/m2 ---> ~5.2%
              ; rsut: 4-5 W/m2 ---> ~4-5%
              ;
              ; ---> use 5% for rsut and rsutcs
              factor = 0.05
            end if
            if (any(variables(iv).eq.(/"rlut", "rlutcs"/))) then
              ; From the Technical Note on CERES EBAF Ed2.7:
              ; "For CERES, calibration uncertainty is 0.5% (1 sigma)
              ; which for a typical global mean clear-sky LW flux corresponds
              ; to ~1 W/m2. [...] The total error in TOA outgoing clear-sky LW
              ; radiation in a region is estimated as
              ; sqrt(1^2 + 1.74^2 + 0.7^2 + 1^2 + 2.75^2) or approximately
              ; 3.6 W/m2."
              ;
              ; --> here, we estimate the total error as
              ;     0.5% * 3.6 W/m2 / 1.0 W/m2 = 1.8%
              factor = 0.018
            end if
          end if

          if (refname.eq."MODIS") then
            if (any(variables(iv).eq.(/"clt"/))) then
              ; From the technical note on obs4MIPS MODIS satellite data
              ; "Moderate Resolution Imaging Spectroradiometer (MODIS)
              ; Cloud Fraction Technical Document":
              ;
              ; Global day/night fractional agreement of cloud detection
              ; between MODIS and CALIPSO lidar (CALIOP):
              ; Aug 2006 (clear):  0.84
              ; Aug 2006 (cloudy): 0.88
              ; Feb 2006 (clear):  0.85
              ; Feb 2006 (cloudy): 0.88
              ;
              ; ---> using fractional agreement as error estimate (1 sigma)
              ; ---> average error estimate =
              ;      sqrt((0.84^2 + 0.88^2 + 0.85^2 + 0.88^2) / 4) = 0.86
              ; ---> factor = 1 - 0.86 = 0.14
              factor = 0.14
            end if
          end if

          if (refname.eq."MODIS-L3") then
            if (any(variables(iv).eq.(/"clt"/))) then
              factor = 0.0  ; no error estimates for clt from MODIS_C6_L3
            end if
          end if
        end if

        ; convert to 1-D arrays
        var1d = ndtooned(var * factor)

        ; --- RMS ---
        ; in case of derived variables (e.g., SW_CRE), the squared error
        ; estimates of all input variables are added before calculating
        ; the square root

        rmsobs(is) = rmsobs(is) + dim_avg_wgt_Wrap(var1d ^ 2, wgt1d, 1)

        ierridx = ierridx + 1

      end do  ; iv-loop
    end if

    ; finish calculation of RMS
    rmsobs(is) = sqrt(rmsobs(is))

    ; normalize RMS by dividing by standard deviation of obs (ref)
    avg_ref = dim_avg_wgt_Wrap(ref1d, wgt1d, 1)
    rmsobs(is) =  rmsobs(is) / sqrt(dim_avg_wgt_Wrap((ref1d - avg_ref) ^ 2, \
                                                     wgt1d, 1))

  end do  ; loop over seasons

  if (isvar("var")) then
    delete(var)
  end if
  if (isvar("ref")) then
    delete(ref)
  end if
  if (isvar("weights")) then
    delete(weights)
  end if
  if (isvar("var1d")) then
    delete(var1d)
  end if
  if (isvar("ref1d")) then
    delete(ref1d)
  end if
  if (isvar("wgt1d")) then
    delete(wgt1d)
  end if

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  nummods = dim_MOD - dim_REF

  colors  = new(nummods, string)
  markers = new(nummods, integer)

  ratio = new((/nummods, 1, numseas/), float)  ; standard deviation ratios
  ratio = val(:, 0, :)
  cc = new((/nummods, 1, numseas/), float)   ; correlations
  cc = val(:, 1, :)

  if (isatt(diag_script_info, "styleset")) then
    colortab = project_style(info0, diag_script_info, "colors")
    markertab = project_style(info0, diag_script_info, "markers")
  else
    colortab = (/"(/0.00, 0.00, 0.59/)", "(/0.00, 0.39, 1.00/)", \
                 "(/0.20, 1.00, 1.00/)", "(/0.20, 0.88, 0.00/)", \
                 "(/1.00, 0.88, 0.00/)", "(/1.00, 0.59, 0.00/)", \
                 "(/1.00, 0.20, 0.00/)", "(/0.59, 0.00, 0.00/)", \
                 "(/0.78, 0.00, 0.78/)", "(/0.59, 0.00, 0.59/)", \
                 "(/0.90, 0.90, 0.90/)", "(/0.70, 0.70, 0.70/)", \
                 "(/0.50, 0.50, 0.50/)", "(/0.30, 0.30, 0.30/)"/)
    markertab = (/16, 4, 5, 0/)
  end if

  plotfile = new((/numseas/), string)

  do is = 0, numseas - 1
    if (isvar("wks")) then
      delete(wks)
    end if

    if (isvar("legendwks")) then
      delete(legendwks)
    end if

    if (isvar("plot")) then
      delete(plot)
    end if

    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_taylor_" + var0 \
                  + "_" + season(is) + filename_add)

    ; create new marker: filled star

    mstring = "z"
    fontnum = 35
    size   = 1.5
    angle  = 0.0

    new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, \
                             1.0, size, angle)

    if (embracelegend) then
      legendwks = get_wks("legend_wks", DIAG_SCRIPT, "clouds_taylor_" \
                          + "legend" + filename_add)
      new_index = NhlNewMarker(legendwks, mstring, fontnum, 0.0, 0.0, \
                               1.0, size, angle)
    else
      legendwks = wks
    end if

    i = 0
    idx = 0

    if (isatt(diag_script_info, "styleset")) then
      do ii = 0, dim_MOD - 1

        if (ii.eq.ref_ind) then
          continue
        end if

        colors(i)  = colortab(i)
        markers(i) = markertab(i)

        ; override multi-model mean (if present)
        if (ii .eq. mm_ind) then
          colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
          markers(i) = new_index
        end if
        ; override alternative obs (if "EMBRACE" legend)
        if (embracelegend) then
          if ((isStrSubset(str_lower(input_file_info@project(i)), "obs"))) \
            then
            if (mm_ind .ge. 0) then
              colors(i)  = "(/1.00, 0.00, 0.00/)"  ; red
            else
              colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
            end if
            markers(i) = new_index
            val&models(i) = "altern. ref. dataset"
          end if
        end if
        i = i + 1
      end do
    else
      do ii = 0, dim_MOD - 1

        if (ii.eq.ref_ind) then
          continue
        end if

        if (ii .eq. mm_ind) then
          colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
          markers(i) = new_index
          i = i + 1
        else
          if (embracelegend) then
            colors(i)  = colortab(idx)
            if ((isStrSubset(str_lower(input_file_info@project(i)), "obs"))) \
              then
              if (mm_ind .ge. 0) then
                colors(i)  = "(/1.00, 0.00, 0.00/)"  ; red
              else
                colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
              end if
              markers(i) = new_index
              val&models(i) = "altern. ref. dataset"
            else
              markers(i) = markertab(mod(i, 2))
            end if
            i = i + 1
            if (i.ge.dim_MOD) then
              break
            end if
            idx = idx + 1
            if (idx.ge.dimsizes(colortab)) then
              idx = dimsizes(colortab) - 1
            end if
          else
            do n = 0, dim_MOD / dimsizes(colortab)
              colors(i)  = colortab(idx)
              markers(i) = markertab(n)
              i = i + 1
              if (i.ge.dim_MOD) then
                break
              end if
            end do
            idx = idx + 1
            if (idx.ge.dimsizes(colortab)) then
              idx = dimsizes(colortab) - 1
            end if
          end if
        end if
      end do
    end if

    ropts = True
    ropts@Colors        = colors
    ropts@Markers       = markers  ; marker styles
    ropts@gsMarkerSizeF = 0.0125

    modelnames          = val&models
    ropts@caseLabels    = modelnames

    ropts@stnRad        = (/ 0.5, 1.5, 2.0, 2.5 /)
    ropts@OneX          = "1.00"
    ropts@ccRays        = (/ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, \
                           0.9, 0.95, 0.99 /)
    ropts@ccRays_color  = "Black"
    ropts@centerDiffRMS = True
    ropts@centerDiffRMS_color = "LightGray"  ; "PaleGreen1"

    ropts@printDataLabels = False  ; print numeric label for each data pt.
    ropts@reverseOrder  = True  ; reverse order of items in legend
    ropts@legendBox     = True  ; draw box around legend

    ; legend parameters

    ropts@legendXpos        = 0.575
    ropts@legendYpos        = -0.5
    ropts@legendWidth       = 0.2

    ropts@caseLabelsFontHeightF = 0.1
    ropts@plotSize              = 0.6

    if (embracelegend) then
      ropts@legendExtraFile = True  ; save legend to extra file
      ropts@modelsperrow    = 7     ; 7 models per column in the legend
    else
      ropts@legendExtraFile = False
      if (dim_MOD.ge.20) then
        ropts@caseLabelsFontHeightF = 0.06
        ropts@plotSize              = 0.45
      end if
    end if

    ; plot
    if (write_plots) then
      ropts@tiMainString   = var0 + " (" + season(is) + ")"  ; title
      ; observational uncertainties as normalized RMS error
      ropts@rmsobs = rmsobs(is)

      plot = taylor_diagram(wks, legendwks, ratio(:, :, is), cc(:, :, is), \
                            ropts)

      plotfile(is) = wks@fullname
      log_info("Wrote " + plotfile)

      if (embracelegend.and.(is.eq.(numseas-1))) then
        frame(legendwks)
      end if
    else
      plotfile(is) = ""
    end if  ; if write_plots
  end do  ; is-loop (seasons)

  ; write netCDF output

  if (any(rmsobs .gt. 0.0)) then
    val@RMSE_observations = rmsobs
  end if

  nc_filename = work_dir + "clouds_taylor_" + var0 + filename_add + ".nc"
  nc_outfile = ncdf_write(val, nc_filename)

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s)
  ; ------------------------------------------------------------------------

  statistics = "clim"
  domain = "global"
  plottype = "taylor"

  do is = 0, numseas - 1
    ; note: because function log_provenance does not yet support to attach
    ;       different captions to netcdf (contains all seasons) and plots
    ;       (contain one season each), the caption cannot specifiy the
    ;       season plotted; using "annual" or "DJF/MAM/JJA/SON" instead.

    caption = "Taylor diagram for variable " + var0 + " (" + allseas \
              + "), reference = " + refname + "."

    log_provenance(nc_outfile, plotfile, caption, statistics, domain, \
                   plottype, "", "", infiles)
  end do

  leave_msg(DIAG_SCRIPT, "")

end
