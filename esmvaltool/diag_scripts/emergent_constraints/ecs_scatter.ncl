; ECS_SCATTER
; #############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME CRESCENDO
; #############################################################################
; Description
;   Calculates equilibrium climate sensitivity (ECS) versus:
;     - southern ITCZ index (Tian, 2015)
;     - tropical mid-tropospheric humidity asymmetry index (Tian, 2015)
;     - lower tropospheric mixing index (LTMI) (Sherwood et al., 2014)
;     - covariance of shortwave cloud reflection (Brient and Schneider, 2016)
;     - climatological Hadley cell extent (Lipat et al., 2017)
;
;    References:
;        - Brient, F., and T. Schneider, J. Climate, 29, 5821-5835,
;          doi: 10.1175/JCLIM-D-15-0897.1 (2016)
;        - Lipat, B. R., et al., Geophys. Res. Lett., 44, 5739-5748,
;          doi: 10.1002/2017GL073151 (2017)
;        - Sherwood et al., nature, 505, 37-42, doi: 10.1038/nature12829 (2014)
;        - Tian, B., Geophys. Res. Lett., 42, 4133-4141,
;          doi: 10.1002/2015GL064119 (2015)
;
; Required diag_script_info attributes (diagnostics specific)
;   diag: emergent constraint to calculate ("itczidx", "humidx", "ltmi",
;         "covrefl", "shhc", "sherwood_s", "sherwood_d")
;   ecs_file:  path + file name of netCDF containing precalculated
;              ECS values (only necessary if 'output_diag_only' is not True)
;
; Optional diag_script_info attributes (diagnostic specific)
;   calcmm:         - False (default) = do not calculate multi-model mean
;                   - True = calculate multi-model mean
;   legend_outside: True, False
;   output_diag_only: - False (default) = output full emergent constraint
;                                         (X and Y axis, plots)
;                     - True = output only X axis (no plots)
;   output_attributes: Additional attributes for the output files
;   output_models_only: Only write models (no reference datasets) to netcdf
;                       files.
;   predef_minmax:  True, False
;   styleset:       "CMIP5" (if not set, ECS_SCATTER will create a color table
;                   and symbols for plotting)
;   suffix: string to add to output filenames (e.g."cmip3")
;
; Required variable attributes
;   reference_dataset: name of reference data set
;
; Optional variable_info attributes
;   none
;
; Caveats
;   none
;
; Modification history
;    20200102-schlund_manuel: added Sherwood S and D index.
;    20190226-lauer_axel: added output of provenance (v2.0)
;    20181022-lauer_axel: code rewritten for ESMValTool v2.0
;    20170818-lauer_axel: preparation of code for reporting (tags)
;    20170810-lauer_axel: added emergent constraint "covariance of shortwave
;                         cloud reflection" (Brient and Schneider, 2016)
;    20161007-lauer_axel: written.
;
; #############################################################################
; List of functions
;   function covrefl: covariance of shortwave cloud reflectivity
;   function humidx: tropical mid-tropopheric humidity asymmetry index
;   function itczidx: southern ITCZ index
;   function ltmi: lower tropospheric mixing index (LTMI)
;   function sherwood_d: Sherwood D index
;   function sherwood_s: Sherwood S index
;   function shhc: Southern hemisphere Hadley cell extent
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"

; =============================================================================
; Function covrefl
; Calculate covariance of shortwave cloud reflectivity.
;
; Arguments
; - dataset: names of datasets (must be of same size as variables)
; - variables: names of available variables (must be of same size as dataset)
; - ldebug: output debug data as netCDF (true, false)
;
; Reference
; - Brient and Schneider (2016)
; =============================================================================

function covrefl(dataset[*]:string, variables[*]:string, ldebug:logical)
begin
  ; get indices associated with specific variables
  ; sst is special since it can be "ts" or alternatively "tos"
  sst_idx = ind(variables.eq."ts")
  if (ismissing(sst_idx)) then
    sst_idx = ind(variables.eq."tos")
  end if
  hur_idx = ind(variables.eq."hur")
  rsut_idx = ind(variables.eq."rsut")
  rsutcs_idx = ind(variables.eq."rsutcs")
  rsdt_idx = ind(variables.eq."rsdt")

  all_idx = (/sst_idx, hur_idx, rsut_idx, rsutcs_idx, rsdt_idx/)

  ; -------------------------------------------------------------------
  ; if needed, regrid all variables to common grid
  ; (e.g. needed if observations are taken from different sources)
  ; for definition of the common grid, the grid used for relative
  ; humidity is used
  ; -------------------------------------------------------------------

  ; get grid of variable 'relative humidity' (hur)

  atts = True
  atts@short_name = variables(hur_idx)
  atts@dataset = dataset(hur_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  rh500 = read_data(info[0])
  delete(info)

  if (.not.isatt(rh500, "_FillValue")) then
    rh500@_FillValue = -999.
  end if

  tsteps0 = dimsizes(rh500&time)

  do iv = 0, dimsizes(variables) - 1

    ; relative humidity has already been read and does not need to be
    ; processed because rh500 defines the reference grid

    if (iv.eq.hur_idx) then
      continue
    end if

    ; do not process any (possibly present) variables that we do not need
    ; for the calculations

    if (.not.any(all_idx .eq. iv)) then
      continue
    end if

    ; extract data of variable #'iv'
    atts = True
    atts@short_name = variables(iv)
    atts@dataset = dataset(iv)
    info = select_metadata_by_atts(input_file_info, atts)
    delete(atts)

    var = read_data(info[0])
    delete(info)

    ; determine dimensions to regrid
    xregrid = False  ; lon
    yregrid = False  ; lat
    dnames = getVarDimNames(var)
    if (any(dnames.eq."lon")) then
      xregrid = True
      if (dimsizes(var&lon).eq.dimsizes(rh500&lon)) then
        if (max(abs(var&lon - rh500&lon)).lt.1.0e-3) then
          xregrid = False
        end if
      end if
    end if
    if (any(dnames.eq."lat")) then
      yregrid = True
      if (dimsizes(var&lat).eq.dimsizes(rh500&lat)) then
        if (max(abs(var&lat - rh500&lat)).lt.1.0e-3) then
          yregrid = False
        end if
      end if
    end if
    ; regrid in the horizontal (if necessary)
    if (xregrid.and.yregrid) then
      var_hreg = linint2_Wrap(var&lon, var&lat, var, True, \
                              rh500&lon, rh500&lat, 0)
    else if (xregrid.and..not.yregrid) then
      dx = ind(dnames.eq."lon")
      var_hreg = linint1_n_Wrap(var&lon, var, True, rh500&lon, 0, dx)
    else if (.not.xregrid.and.yregrid) then
      dx = ind(dnames.eq."lat")
      var_hreg = linint1_n_Wrap(var&lat, var, False, rh500&lat, 0, dx)
    else if (.not.xregrid .and. .not.yregrid) then
      var_hreg = var
    end if
    end if
    end if
    end if

    delete(var)
    delete(dnames)

    ; check if number of time steps are identical for all variables

    tsteps = dimsizes(rh500&time)
    if (tsteps.ne.tsteps0) then
      infostr = variables(hur_idx) + ": " + sprinti("%i", tsteps0) + ", " \
                + variables(iv) + ": " + sprinti("%i", tsteps)
      model1 = names(imod(hur_idx))
      model2 = names(imod(iv))
      error_msg("f", DIAG_SCRIPT, "", "number of time steps for " \
                + "different variables not consistent (model = " \
                + model1 + ", " + model2 + "). " + infostr)
    end if

    ; rename (regridded) field to something meaningful

    if (iv.eq.rsut_idx)
      rsut = var_hreg
    end if
    if (iv.eq.rsutcs_idx)
      rsutcs = var_hreg
    end if
    if (iv.eq.rsdt_idx)
      rsdt = var_hreg
    end if
    if (iv.eq.sst_idx)
      sst = var_hreg
    end if

    delete(var_hreg)

  end do  ; iv-loop (regridding if needed)

  ; -------------------------------------------------------------------

  ; 500 hPa relative humidity

  ; calculate shortwave cloud forcing

  sw_cre = rsutcs
  sw_cre = sw_cre - rsut

  delete(rsut)
  delete(rsutcs)

  ; generate land-/sea-mask
  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
  lsdata = a->LSMASK
  lsm0  = landsea_mask(lsdata, rh500&lat, rh500&lon)
  lsm = where(lsm0.eq.0, 0, lsm0@_FillValue)
  copy_VarMeta(lsm0, lsm)
  delete(lsdata)
  delete(lsm0)

  ; lower quartile of the mid-tropospheric (500 hPa) relative humidity
  ; over the tropical ocean, i.e. 30S-30N
  rh500 = rh500 + conform(rh500, lsm, (/1, 2/))
  delete(lsm)
  rh500(:, {-90.:-30.}, :) = rh500@_FillValue
  rh500(:, {30.:90.}, :) = rh500@_FillValue

  if (ldebug) then
    tmp = get_unique_values(dataset)
    if (dimsizes(tmp) .eq. 1) then
      refstr = tmp
    else
      refstr = tmp(0)
      do i = 1, dimsizes(tmp) - 1
        refstr = refstr + "_" + tmp(i)
      end do
    end if
    fname = config_user_info@work_dir + "/debug_" + refstr + ".nc"
    system("rm " + fname)
    debugfile = addfile(fname, "c")
    debugout = time_operations(rh500, -1, -1, "average", "annualclim", True)
    statb = stat_dispersion(debugout, False)
    q25 = statb(6)
    debugout = where(debugout.le.q25, 0.0, debugout@_FillValue)
    debugfile->mask = debugout
    debugfile->rh500mask = rh500
  end if

  dims = dimsizes(rh500)
  do t = 0, dims(0) - 1
    statb = stat_dispersion(rh500(t, :, :), False)
    q25 = statb(6)
    ; mask all grid cells with relative humidity values larger than the
    ; lower quartile
    rh500(t, :, :) = where(rh500(t, :, :).le.q25, 0.0, rh500@_FillValue)
  end do
  delete(dims)

  if (ldebug) then
    debugfile->rh500mask_final = rh500
  end if

  ; apply rh500 mask to shortwave cloud forcing, solar insolation and sst
  sw_cre = sw_cre + rh500(:, :, :)
  sst = sst + rh500(:, :, :)
  rsdt = rsdt + rh500(:, :, :)
  delete(rh500)

  if (ldebug) then
    debugfile->sst = sst
    debugfile->sw_cre = sw_cre
    debugfile->rsdt = rsdt
  end if

  ; calculate averages (missing values will be ignored and grid boxes
  ; area-weighted) over TLC areas (determined by rh500 mask)
  t1 = area_operations(sw_cre, -90., 90., 0., 360., "average", True)
  t2 = area_operations(rsdt, -90., 90., 0., 360., "average", True)
  t3 = area_operations(sst, -90., 90., 0., 360., "average", True)

  delete(sw_cre)
  delete(rsdt)
  delete(sst)

  ; shortwave cloud reflectivity (in %)
  alpha = t1
  alpha = -100. * alpha / t2

  ; deseasonalize time series
  ; calculate climatological seasonal cycles
  seasclim0 = time_operations(alpha, -1, -1, "average", "monthlyclim", True)
  seasclim1 = time_operations(t3, -1, -1, "average", "monthlyclim", True)

  ; subtract climatological seasonal cycles

  alpha_dseas = alpha
  sst_dseas = t3

  do i = 0, 11
    alpha_dseas(i::12) = alpha_dseas(i::12) - seasclim0(i)
    sst_dseas(i::12) = sst_dseas(i::12) - seasclim1(i)
  end do

  delete(seasclim0)
  delete(seasclim1)
  delete(alpha)

  ; calculate regression coefficient alpha_dseas vs. sst_dseas
  rc = regline(sst_dseas, alpha_dseas)

  if (ldebug) then
    debugfile->alpha_dseas = alpha_dseas
    debugfile->sst_dseas = sst_dseas
    t1 = -t1 / t2
    debugfile->alpha = t1
  end if

  delete(t1)
  delete(t2)
  delete(t3)

  delete(alpha_dseas)
  delete(sst_dseas)

  return(rc)
end

; =============================================================================
; Function humidx
; Calculate tropical mid-tropopheric humidity asymmetry index.
;
; Arguments
; - dataset: name of dataset
; - variables: names of available variables
; - ldebug: output debug data as netCDF (true, false)
; - nh_ref: average NH tropical Pacific, reference dataset
; - sh_ref: average SH tropical Pacific, reference dataset
;
; Reference
; - Tian (2015)
; =============================================================================

function humidx(dataset:string, variables[*]:string, ldebug:logical, \
                nh_ref:numeric, sh_ref:numeric)
begin
  hus_idx = ind(variables.eq."hus")

  atts = True
  atts@short_name = variables(hus_idx)
  atts@dataset = dataset
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  A0 = read_data(info[0])

  ; calculate climatological annual mean
  q500 = time_operations(A0, -1, -1, "average", "annualclim", True)

  ; average SH tropical Pacific (120E-80W, 30S-0)
  sh_avg = area_operations(q500, -30., 0., 120., 280., "average", True)
  ; average NH tropical Pacific (120E-80W, 0-20N)
  nh_avg = area_operations(q500, 0., 20., 120., 280., "average", True)

  delete(q500)

  result = (sh_avg - sh_ref) / sh_ref * 100.0 \
    - (nh_avg - nh_ref) / nh_ref * 100.0
  result@standard_name = "Tropical midtropospheric humidity asymmetry index"
  result@units = "%"

  ; process error estimate if available (for observations)

  minmax = new((/2/), float)
  minmax@standard_name = result@standard_name + " error estimate"
  minmax@units = result@units

  ierr = ind(variables.eq."husStderr")

  if (.not.ismissing(ierr)) then
    atts = True
    atts@short_name = variables(ierr)
    atts@dataset = dataset
    info = select_metadata_by_atts(input_file_info, atts)
    delete(atts)
    n = ListCount(info)

    if (n .eq. 1) then
      A0err = read_data(info[0])
      delete(info)

      ; *********************
      ; *** mean + stderr ***
      ; *********************

      A0 = A0 + A0err

      ; calculate climatological annual mean
      q500 = time_operations(A0, -1, -1, "average", "annualclim", True)

      ; average SH tropical Pacific (120E-80W, 30S-0)
      sh1 = area_operations(q500, -30., 0., 120., 280., "average", True)
      ; average NH tropical Pacific (120E-80W, 0-20N)
      nh1 = area_operations(q500, 0., 20., 120., 280., "average", True)
      delete(q500)

      minmax(1) = (sh1 - sh_ref) / sh_ref * 100.0 \
        - (nh1 - nh_ref) / nh_ref * 100.0

      ; *********************
      ; *** mean - stderr ***
      ; *********************

      A0 = A0 - 2.0 * A0err
      delete(A0err)
      A0 = where(A0.ge.0.0, A0, 0.0)

      ; calculate climatological annual mean
      q500 = time_operations(A0, -1, -1, "average", "annualclim", True)

      ; average SH tropical Pacific (120E-80W, 30S-0)
      sh0 = area_operations(q500, -30., 0., 120., 280., "average", True)
      ; average NH tropical Pacific (120E-80W, 0-20N)
      nh0 = area_operations(q500, 0., 20., 120., 280., "average", True)
      delete(q500)

      minmax(0) = (sh0 - sh_ref) / sh_ref * 100.0 \
        - (nh0 - nh_ref) / nh_ref * 100.0
    end if  ; if dataset for error estimate available
  end if  ; if error variable present

  delete(A0)

  return([/ result, nh_avg, sh_avg, minmax /])
end

; =============================================================================
; Function itczidx
; Calculate southern ITCZ index.
;
; Arguments
; - dataset: name of dataset
; - variables: names of available variables
; - ldebug: output debug data as netCDF (true, false)
; - SEPavg_ref: average over southeastern Pacific for reference dataset
;   in mm/day
;
; Reference
; - Tian (2015)
; =============================================================================

function itczidx(dataset:string, variables[*]:string, ldebug:logical, \
                 SEPavg_ref:numeric)
begin
  pr_idx = ind(variables.eq."pr")

  atts = True
  atts@short_name = variables(pr_idx)
  atts@dataset = dataset
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  pr = read_data(info[0])
  delete(info)

  ; average over time (annual mean)
  pr_mean = time_operations(pr, -1, -1, "average", "annualclim", True)

  ; average over southeastern Pacific (150W-100W, 30S-0)
  SEPavg = area_operations(pr_mean, -30., 0., 210., 260., "average", True)
  delete(pr_mean)
  SEPavg = SEPavg * 3600.0 * 24.0  ; convert from kg/m2/s to mm/day
  result = SEPavg - SEPavg_ref

  result@standard_name = "Southern ITCZ index"
  result@units = "mm day-1"

  ; process error estimate if available (for observations)

  minmax = new((/2/), float)
  minmax@standard_name = result@standard_name + " error estimate"
  minmax@units = result@units

  ierr = ind(variables.eq."prStderr")

  if (.not.ismissing(ierr)) then
    atts = True
    atts@short_name = variables(ierr)
    atts@dataset = dataset
    info = select_metadata_by_atts(input_file_info, atts)
    delete(atts)
    n = ListCount(info)

    if (n .eq. 1) then
      pr_err = read_data(info[0])
      delete(info)

      ; mean + stderr
      pr = pr + pr_err
      tmp1 = time_operations(pr, -1, -1, "average", "annualclim", True)
      minmax(1) = area_operations(tmp1, -30., 0., 210., 260., "average", \
                                  True)
      delete(tmp1)

      ; max(0, mean - stderr)
      pr = pr - 2.0 * pr_err
      delete(pr_err)
      pr = where(pr.ge.0.0, pr, 0.0)
      tmp0 = time_operations(pr, -1, -1, "average", "annualclim", True)
      minmax(0) = area_operations(tmp0, -30., 0., 210., 260., "average", \
                                  True)
      delete(tmp0)

      minmax = minmax * 3600.0 * 24.0  ; convert from kg/m2/s to mm/day
      minmax = minmax - SEPavg_ref

    end if  ; if dataset for error estimate available
  end if  ; if error variable present

  delete(pr)

  return([/ result, minmax /])

end

; =============================================================================
; Function ltmi
; Calculate lower tropospheric mixing index (LTMI) for one dataset.
;
; Arguments
; - dataset: names of datasets
; - variables: names of available variables
; - ldebug: output debug data as netCDF (true, false)
;
; Reference
; - Sherwood et al. (2014)
; =============================================================================

function ltmi(dataset[*]:string, variables[*]:string, ldebug:logical)
begin
  hur_idx = ind(variables.eq."hur")
  ta_idx = ind(variables.eq."ta")
  wap_idx = ind(variables.eq."wap")

  ; vertical velocity (at 850, 700, 600, 500, 400 hPa)

  atts = True
  atts@short_name = variables(wap_idx)
  atts@dataset = dataset(wap_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  omega = read_data(info[0])
  delete(info)

  omega_mean = time_operations(omega, -1, -1, "average", "annualclim", True)

  ; grid boxes with upward motion at 500 hPa (annual mean)
  ascend = where(omega_mean({50000}, :, :).lt.0, \
                 -1.0 * omega_mean({50000}, :, :), omega_mean@_FillValue)
  ascend@_FillValue = omega_mean@_FillValue
  ascend!0 = "lat"
  ascend!1 = "lon"
  ascend&lat = omega_mean&lat
  ascend&lon = omega_mean&lon

  ; generate land-/sea-mask
  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
  lsdata = a->LSMASK
  lsm0  = landsea_mask(lsdata, omega_mean&lat, omega_mean&lon)
  lsm = where(lsm0.eq.0, 0, lsm0@_FillValue)
  copy_VarMeta(lsm0, lsm)
  delete(lsdata)
  delete(lsm0)
  delete(omega_mean)

  ; upper quartile of the annual mean mid-tropospheric (500 hPa) ascent rate
  ; in ascending regions (over the tropical ocean 30S-30N only)

  ascend = ascend + lsm
  ascend({-90.:-30.}, :) = ascend@_FillValue
  ascend({30.:90.}, :) = ascend@_FillValue

  statb = stat_dispersion(ascend, False)
  q25 = statb(10)

  ; mask all grid cells with ascend rates lower than the upper quartile
  omgmask = where(ascend.ge.q25, 0.0, ascend@_FillValue)
  omgmask@_FillValue = ascend@_FillValue

  ; temperature
  atts = True
  atts@short_name = variables(ta_idx)
  atts@dataset = dataset(ta_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  ta = read_data(info[0])
  delete(info)

  ta_mean = time_operations(ta, -1, -1, "average", "annualclim", True)
  delete(ta)

  ; relative humidity
  atts = True
  atts@short_name = variables(hur_idx)
  atts@dataset = dataset(hur_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  rh = read_data(info[0])
  delete(info)

  rh_mean = time_operations(rh, -1, -1, "average", "annualclim", True)
  delete(rh)

  ; -------------------------------------------------------------------------
  ; if needed, regrid (horizontally) all fields to common grid
  ; (needed since some models provide variables on inconsistent grids)
  ; for definition of the common grid, the grid used for omega (mask) is used
  ; -------------------------------------------------------------------------

  do i = 0, 1
    if (i.eq.0) then
      var = ta_mean
      delete(ta_mean)
    else
      var = rh_mean
      delete(rh_mean)
    end if

    xregrid = False  ; lon
    yregrid = False  ; lat
    dnames = getVarDimNames(var)
    if (any(dnames.eq."lon")) then
      xregrid = True
      if (dimsizes(var&lon).eq.dimsizes(omega&lon)) then
        if (max(abs(var&lon - omega&lon)).lt.1.0d-5) then
          xregrid = False
        end if
      end if
    end if
    if (any(dnames.eq."lat")) then
      yregrid = True
      if (dimsizes(var&lat).eq.dimsizes(omega&lat)) then
        if (max(abs(var&lat - omega&lat)).lt.1.0d-5) then
          yregrid = False
        end if
      end if
    end if
    ; regrid horizontally (if necessary)
    if (xregrid.and.yregrid) then
      var_hreg = linint2_Wrap(var&lon, var&lat, var, True, \
                              omega&lon, omega&lat, 0)
    else if (xregrid.and..not.yregrid) then
      dx = ind(dnames.eq."lon")
      var_hreg = linint1_n_Wrap(var&lon, var, True, omega&lon, 0, dx)
    else if (.not.xregrid.and.yregrid) then
      dx = ind(dnames.eq."lat")
      var_hreg = linint1_n_Wrap(var&lat, var, False, omega&lat, 0, dx)
    else if (.not.xregrid .and. .not.yregrid) then
      var_hreg = var
    end if
    end if
    end if
    end if

    delete(var)
    delete(dnames)

    if (i.eq.0) then
      ta_mean = var_hreg
    else
      rh_mean = var_hreg
    end if

    delete(var_hreg)
  end do

  ; -------------------------------------------------------------------------

  ; calculate "S" = small-scale lower-tropospheric mixing
  ; (CMOR standard unit for relative humidity (hur) is percent (%))

  ; ((rh(700)-rh(850))/ 100% - (T(700)-T(850)) / 9) / 2

  S = ((rh_mean({70000}, :, :) - rh_mean({85000}, :, :)) / 100.0 \
       - (ta_mean({70000}, :, :) - ta_mean({85000}, :, :)) / 9.0) / 2.0

  ; apply udraft mask + land-/sea-mask

  S = S + omgmask + lsm
  S!0 = "lat"
  S!1 = "lon"
  S&lat = ta_mean&lat
  S&lon = ta_mean&lon

  S_avg = area_operations(S, -90., 90., 0., 360., "average", True)

  ; average omega at 850 hPa, 700 hPa
  omega1 = (omega(:, {85000}, :, :) + omega(:, {70000}, :, :)) / 2.0
  ; average omega at 600 hPa, 500 hPa, 400 hPa
  omega2 = (omega(:, {60000}, :, :) + omega(:, {50000}, :, :) \
            + omega(:, {40000}, :, :)) / 3.0
  delta = omega2 - omega1

  D1 = where((delta.gt.0.).and.(omega1.lt.0.), delta, 0.)
  D2 = where(omega2.lt.0., -1. * omega2, 0.)

  D1!0 = "time"
  D1!1 = "lat"
  D1!2 = "lon"
  D1&time = omega&time
  D1&lat = omega&lat
  D1&lon = omega&lon

  D2!0 = "time"
  D2!1 = "lat"
  D2!2 = "lon"
  D2&time = omega&time
  D2&lat = omega&lat
  D2&lon = omega&lon

  ; average over tropical ocean regions: 160W-30E, 30S-30N
  ; (ocean grid cells only)

  lsmC = conform_dims(dimsizes(D1), lsm, (/1, 2/))
  D1 = D1 + lsmC
  D2 = D2 + lsmC
  delete(lsmC)

  D1_avg = avg(area_operations(D1, -30., 30., 200., 30., "average", True))
  D2_avg = avg(area_operations(D2, -30., 30., 200., 30., "average", True))

  if (D2_avg.ne.0.) then
    D_avg = D1_avg / D2_avg
  else
    D_avg = D@_FillValue
  end if

  ; optionally write debug output

  if (ldebug) then
    tmp = get_unique_values(dataset)
    if (dimsizes(tmp) .eq. 1) then
      refstr = tmp
    else
      refstr = tmp(0)
      do i = 1, dimsizes(tmp) - 1
        refstr = refstr + "_" + tmp(i)
      end do
    end if
    fname = config_user_info@work_dir + "/debug_" + refstr + ".nc"
    system("rm " + fname)
    debugfile = addfile(fname, "c")
    debugfile->S = S
    debugfile->omgmask = omgmask
    debugfile->ascend = ascend
    debugfile->lsm = lsm
    debugfile->delta = delta
    debugfile->D1 = D1
    debugfile->D2 = D2
  end if

  result = D_avg + S_avg  ; lower tropospheric mixing index

  result@standard_name = "LTMI"
  result@long_name = "lower_tropospheric_mixing_index"
  result@units = "1"
  result@S = S_avg
  result@D = D_avg

  return(result)
end

; =============================================================================
; Function sherwood_d
; Calculate Sherwood D index.
;
; Arguments
; - dataset: names of datasets
; - variables: names of available variables
; - ldebug: output debug data as netCDF (true, false)
;
; Reference
; - Sherwood et al. (2014)
; =============================================================================

function sherwood_d(dataset[*]:string, variables[*]:string, ldebug:logical)
begin
  wap_idx = ind(variables.eq."wap")

  ; vertical velocity (at 850, 700, 600, 500, 400 hPa)

  atts = True
  atts@short_name = variables(wap_idx)
  atts@dataset = dataset(wap_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  omega = read_data(info[0])
  delete(info)

  omega_mean = time_operations(omega, -1, -1, "average", "annualclim", True)

  ; grid boxes with upward motion at 500 hPa (annual mean)
  ascend = where(omega_mean({50000}, :, :).lt.0, \
                 -1.0 * omega_mean({50000}, :, :), omega_mean@_FillValue)
  ascend@_FillValue = omega_mean@_FillValue
  ascend!0 = "lat"
  ascend!1 = "lon"
  ascend&lat = omega_mean&lat
  ascend&lon = omega_mean&lon

  ; generate land-/sea-mask
  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
  lsdata = a->LSMASK
  lsm0  = landsea_mask(lsdata, omega_mean&lat, omega_mean&lon)
  lsm = where(lsm0.eq.0, 0, lsm0@_FillValue)
  copy_VarMeta(lsm0, lsm)
  delete(lsdata)
  delete(lsm0)
  delete(omega_mean)

  ; upper quartile of the annual mean mid-tropospheric (500 hPa) ascent rate
  ; in ascending regions (over the tropical ocean 30S-30N only)

  ascend = ascend + lsm
  ascend({-90.:-30.}, :) = ascend@_FillValue
  ascend({30.:90.}, :) = ascend@_FillValue

  statb = stat_dispersion(ascend, False)
  q25 = statb(10)

  ; mask all grid cells with ascend rates lower than the upper quartile
  omgmask = where(ascend.ge.q25, 0.0, ascend@_FillValue)
  omgmask@_FillValue = ascend@_FillValue

  ; average omega at 850 hPa, 700 hPa
  omega1 = (omega(:, {85000}, :, :) + omega(:, {70000}, :, :)) / 2.0
  ; average omega at 600 hPa, 500 hPa, 400 hPa
  omega2 = (omega(:, {60000}, :, :) + omega(:, {50000}, :, :) \
            + omega(:, {40000}, :, :)) / 3.0
  delta = omega2 - omega1

  D1 = where((delta.gt.0.).and.(omega1.lt.0.), delta, 0.)
  D2 = where(omega2.lt.0., -1. * omega2, 0.)

  D1!0 = "time"
  D1!1 = "lat"
  D1!2 = "lon"
  D1&time = omega&time
  D1&lat = omega&lat
  D1&lon = omega&lon

  D2!0 = "time"
  D2!1 = "lat"
  D2!2 = "lon"
  D2&time = omega&time
  D2&lat = omega&lat
  D2&lon = omega&lon

  ; average over tropical ocean regions: 160W-30E, 30S-30N
  ; (ocean grid cells only)

  lsmC = conform_dims(dimsizes(D1), lsm, (/1, 2/))
  D1 = D1 + lsmC
  D2 = D2 + lsmC
  delete(lsmC)

  D1_avg = avg(area_operations(D1, -30., 30., 200., 30., "average", True))
  D2_avg = avg(area_operations(D2, -30., 30., 200., 30., "average", True))

  if (D2_avg.ne.0.) then
    result = D1_avg / D2_avg
  else
    result = D@_FillValue
  end if

  ; optionally write debug output

  result@standard_name = "Sherwood_D"
  result@long_name = "Sherwood D index"
  result@units = "1"

  return(result)
end

; =============================================================================
; Function sherwood_s
; Calculate Sherwood S index
;
; Arguments
; - dataset: names of datasets
; - variables: names of available variables
; - ldebug: output debug data as netCDF (true, false)
;
; Reference
; - Sherwood et al. (2014)
; =============================================================================

function sherwood_s(dataset[*]:string, variables[*]:string, ldebug:logical)
begin
  hur_idx = ind(variables.eq."hur")
  ta_idx = ind(variables.eq."ta")
  wap_idx = ind(variables.eq."wap")

  ; vertical velocity (at 850, 700, 600, 500, 400 hPa)

  atts = True
  atts@short_name = variables(wap_idx)
  atts@dataset = dataset(wap_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  omega = read_data(info[0])
  delete(info)

  omega_mean = time_operations(omega, -1, -1, "average", "annualclim", True)

  ; grid boxes with upward motion at 500 hPa (annual mean)
  ascend = where(omega_mean({50000}, :, :).lt.0, \
                 -1.0 * omega_mean({50000}, :, :), omega_mean@_FillValue)
  ascend@_FillValue = omega_mean@_FillValue
  ascend!0 = "lat"
  ascend!1 = "lon"
  ascend&lat = omega_mean&lat
  ascend&lon = omega_mean&lon

  ; generate land-/sea-mask
  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
  lsdata = a->LSMASK
  lsm0  = landsea_mask(lsdata, omega_mean&lat, omega_mean&lon)
  lsm = where(lsm0.eq.0, 0, lsm0@_FillValue)
  copy_VarMeta(lsm0, lsm)
  delete(lsdata)
  delete(lsm0)
  delete(omega_mean)

  ; upper quartile of the annual mean mid-tropospheric (500 hPa) ascent rate
  ; in ascending regions (over the tropical ocean 30S-30N only)

  ascend = ascend + lsm
  ascend({-90.:-30.}, :) = ascend@_FillValue
  ascend({30.:90.}, :) = ascend@_FillValue

  statb = stat_dispersion(ascend, False)
  q25 = statb(10)

  ; mask all grid cells with ascend rates lower than the upper quartile
  omgmask = where(ascend.ge.q25, 0.0, ascend@_FillValue)
  omgmask@_FillValue = ascend@_FillValue

  ; temperature
  atts = True
  atts@short_name = variables(ta_idx)
  atts@dataset = dataset(ta_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  ta = read_data(info[0])
  delete(info)

  ta_mean = time_operations(ta, -1, -1, "average", "annualclim", True)
  delete(ta)

  ; relative humidity
  atts = True
  atts@short_name = variables(hur_idx)
  atts@dataset = dataset(hur_idx)
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  rh = read_data(info[0])
  delete(info)

  rh_mean = time_operations(rh, -1, -1, "average", "annualclim", True)
  delete(rh)

  ; -------------------------------------------------------------------------
  ; if needed, regrid (horizontally) all fields to common grid
  ; (needed since some models provide variables on inconsistent grids)
  ; for definition of the common grid, the grid used for omega (mask) is used
  ; -------------------------------------------------------------------------

  do i = 0, 1
    if (i.eq.0) then
      var = ta_mean
      delete(ta_mean)
    else
      var = rh_mean
      delete(rh_mean)
    end if

    xregrid = False  ; lon
    yregrid = False  ; lat
    dnames = getVarDimNames(var)
    if (any(dnames.eq."lon")) then
      xregrid = True
      if (dimsizes(var&lon).eq.dimsizes(omega&lon)) then
        if (max(abs(var&lon - omega&lon)).lt.1.0d-5) then
          xregrid = False
        end if
      end if
    end if
    if (any(dnames.eq."lat")) then
      yregrid = True
      if (dimsizes(var&lat).eq.dimsizes(omega&lat)) then
        if (max(abs(var&lat - omega&lat)).lt.1.0d-5) then
          yregrid = False
        end if
      end if
    end if
    ; regrid horizontally (if necessary)
    if (xregrid.and.yregrid) then
      var_hreg = linint2_Wrap(var&lon, var&lat, var, True, \
                              omega&lon, omega&lat, 0)
    else if (xregrid.and..not.yregrid) then
      dx = ind(dnames.eq."lon")
      var_hreg = linint1_n_Wrap(var&lon, var, True, omega&lon, 0, dx)
    else if (.not.xregrid.and.yregrid) then
      dx = ind(dnames.eq."lat")
      var_hreg = linint1_n_Wrap(var&lat, var, False, omega&lat, 0, dx)
    else if (.not.xregrid .and. .not.yregrid) then
      var_hreg = var
    end if
    end if
    end if
    end if

    delete(var)
    delete(dnames)

    if (i.eq.0) then
      ta_mean = var_hreg
    else
      rh_mean = var_hreg
    end if

    delete(var_hreg)
  end do

  ; -------------------------------------------------------------------------

  ; calculate "S" = small-scale lower-tropospheric mixing
  ; (CMOR standard unit for relative humidity (hur) is percent (%))

  ; ((rh(700)-rh(850))/ 100% - (T(700)-T(850)) / 9) / 2

  S = ((rh_mean({70000}, :, :) - rh_mean({85000}, :, :)) / 100.0 \
       - (ta_mean({70000}, :, :) - ta_mean({85000}, :, :)) / 9.0) / 2.0

  ; apply udraft mask + land-/sea-mask

  S = S + omgmask + lsm
  S!0 = "lat"
  S!1 = "lon"
  S&lat = ta_mean&lat
  S&lon = ta_mean&lon

  result = area_operations(S, -90., 90., 0., 360., "average", True)

  result@standard_name = "Sherwood_S"
  result@long_name = "Sherwood S index"
  result@units = "1"

  return(result)
end

; =============================================================================
; Function shhc
; Calculate Southern hemisphere Hadley cell extent for one dataset.
;
; Arguments
; - imod: dataset index
; - ldebug: output debug data as netCDF (true, false)
;
; Reference
; - Lipat et al. (2017)
; =============================================================================

function shhc(dataset:string, variables[*]:string, ldebug:logical)
begin
  va_idx = ind(variables.eq."va")

  atts = True
  atts@short_name = variables(va_idx)
  atts@dataset = dataset
  info = select_metadata_by_atts(input_file_info, atts)
  delete(atts)

  va = read_data(info[0])
  delete(info)

  ; check units of meridional wind
  if (va@units.ne."m s-1") then
    error_msg("f", DIAG_SCRIPT, "", "meridional wind 'va' not in expected " \
              + "units (m s-1)")
  end if

  ; average over DJF
  va_seas = time_operations(va, -1, -1, "average", "DJF", True)

  delete(va)

  ; if needed, invert levels
  n = dimsizes(va_seas&plev)
  if (va_seas&plev(0) .gt. va_seas&plev(n-1)) then
    tmp = va_seas(::-1, :, :)
    tmp&plev = va_seas&plev(::-1)
    delete(va_seas)
    va_seas = tmp
    delete(tmp)
  end if
  delete(n)

  ; create dummy surface pressure field (meridional wind is on pressure levels)
  n = dimsizes(va_seas)
  ps = new((/n(1), n(2)/), float)
  ps = 101325.  ; Pa
  delete(n)

  ; calculate zonal average stream function
  ; note: pressure must be between 5 and 1005 hPa (exclusive) and
  ;       monotonically increasing
  i0 = ind((va_seas&plev .gt. 500.) .and. (va_seas&plev .lt. 100500.))
  psi = zonal_mpsi_Wrap(va_seas(i0, :, :), va_seas&lat, va_seas&plev(i0), ps)
  delete(i0)
  delete(va_seas)
  delete(ps)
  ; 500 hPa stream function
  psi500 = int2p_n_Wrap(psi&plev, psi, (/50000./), 1, 0)

  ; optionally write debug output

  if (ldebug) then
    fname = config_user_info@work_dir + "/debug_" + dataset + ".nc"
    system("rm " + fname)
    debugfile = addfile(fname, "c")
    debugfile->psi = psi
  end if

  delete(psi)

  ; calculate Hadley cell (HC) edge:
  ;   - identify first two grid cells from the equator where psi500 changes
  ;     sign (note: sign change from minus to plus as we are interested in the
  ;     downward branch)
  ;   - interpolate linearly between these two grid points to find the latitude
  ;     of the zero crossing

  n = dimsizes(psi500&lat)
  if (psi500&lat(0).le.psi500&lat(n-1)) then
    istart = max(ind(psi500&lat .lt. 0.))
    iend = 0
    idelta = -1
  else
    istart = min(ind(psi500&lat .lt. 0.))
    iend = n - 1
    idelta = 1
  end if
  delete(n)

  if (psi500(0, istart) .le.0) then
    sgn0 = -1
  else
    sgn0 = 1
  end if

  i = istart
  do j = min((/istart, iend/)), max((/istart, iend/)) - 1
    if (psi500(0, i) .le.0) then
      sgn1 = -1
    else
      sgn1 = 1
    end if

    if ((sgn0 .eq. -1) .and. (sgn1 .eq. 1)) then
      break
    end if

    sgn0 = sgn1

    i = i + idelta
  end do

  delta_y = tofloat(psi500(0, i) - psi500(0, i - idelta))
  delta_x = tofloat(psi500&lat(i) - psi500&lat(i - idelta))

  result = -1.0 * psi500(0, i - idelta) * delta_x / delta_y \
    + tofloat(psi500&lat(i - idelta))

  if (ldebug) then
    debugfile->psi500 = psi500
  end if

  return(result)
end

; =============================================================================
; =============================================================================
; =============================================================================
;                                main program
; =============================================================================
; =============================================================================
; =============================================================================

begin

  enter_msg(DIAG_SCRIPT, "")

  variables = metadata_att_as_array(input_file_info, "short_name")
  variables := get_unique_values(variables)
  dim_VAR = dimsizes(variables)

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + variables + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  tmp = metadata_att_as_array(input_file_info, "dataset")

  names = get_unique_values(tmp)  ; note: get_unique_values sorts array!

  delete(tmp)
  dim_MOD = dimsizes(names)

  infiles = metadata_att_as_array(input_file_info, "filename")
  var0 = variables(0)

  ; write references

  if (isatt(diag_script_info, "diag")) then
    diag = diag_script_info@diag
  else
    error_msg("f", DIAG_SCRIPT, "", "no diagnostic selected in recipe " + \
              "('diag' not set).")
  end if

  caption = "Scatterplot of equilibrium climate sensitivity (ECS) vs. "

  ; check that required variables are there

  ; southern ITCZ index (Tian, Geophys. Res. Lett., 2015)

  if (diag .eq. "itczidx") then
    varlist = (/"pr"/)
    caption = caption + "southern ITCZ index (Tian et al., 2015)."
  end if

  ; tropical mid-tropospheric humidity asymmetry
  ; (Tian, Geophys. Res. Lett., 2015)
  if (diag .eq. "humidx") then
    varlist = (/"hus"/)
    caption = caption + "tropical mid-tropospheric humidity asymetry " \
      + "(Tian et al., 2015)."
  end if

  ; lower tropospheric mixing index (LTMI) (Sherwood et al., Nature, 2014)
  if (diag .eq. "ltmi") then
    varlist = (/"ta", "hur", "wap"/)
    caption = caption + "lower tropospheric mixing index (LTMI) " \
      + "(Sherwood et al., 2014)."
  end if

  ; Sherwood D index (Sherwood et al., Nature, 2014)
  if (diag .eq. "sherwood_d") then
    varlist = (/"wap"/)
    caption = caption + "Sherwood D index (Sherwood et al., 2014)."
  end if

  ; Sherwood S index (Sherwood et al., Nature, 2014)
  if (diag .eq. "sherwood_s") then
    varlist = (/"ta", "hur", "wap"/)
    caption = caption + "Sherwood S index (Sherwood et al., 2014)."
  end if

  ; SH Hadley cell extent (SHHC) (Lipat et al., 2017, Geophys. Res. Lett.)
  if (diag .eq. "shhc") then
    varlist = (/"va"/)
    caption = caption + "Southern Hemisphere Hadley cell extent " \
      + "(Lipat et al., 2017)."
  end if

  ; covariance of shortwave cloud reflection
  ; (Brient and Schneider, J. Climate, 2016)
  if (diag .eq. "covrefl") then
    varlist = (/"ts", "hur", "rsut", "rsutcs", "rsdt"/)
    sst_idx = ind(variables.eq.varlist(0))
    if (ismissing(sst_idx)) then
      varlist(0) = "tos"
    end if
    caption = caption + "tropical low clouds (TLC) reflection feedback " \
      + "(Brient and Schneider, 2016)."
  end if

  if (.not. isvar("varlist")) then
    error_msg("f", DIAG_SCRIPT, "", "selected diagnostic not implemented")
  end if

  idx = new(dimsizes(varlist), integer)

  do i = 0, dimsizes(varlist) - 1
    idx(i) = ind(variables.eq.varlist(i))
  end do

  if (any(ismissing(idx))) then
    errstr = "diagnostic " + diag + " requires the following variable(s): " \
             + str_join(varlist, ", ")
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  delete(idx)

  ; make sure path for (optional) netcdf output exists

  if (config_user_info@write_netcdf.eq."True") then
    write_nc = True
    work_dir = config_user_info@work_dir + "/"
    ; Create work dir
    system("mkdir -p " + work_dir)
  else
    write_nc = False
  end if

  flag_mean = False
  mmm_ind = -1
  if (isatt(diag_script_info, "calcmm")) then
    if (diag_script_info@calcmm) then
      flag_mean = True
      modelsatname = array_append_record(names, \
                                         (/"MultiModelMean"/), 0)
      mmm_ind = dim_MOD
      dim_MOD = dim_MOD + 1
      delete(names)
      names = modelsatname
      delete(modelsatname)
    end if
  end if

  ; reference dataset(s)

  refname = new(dim_VAR, string)
  do i = 0, dim_VAR - 1
    short_name = variables(i)
    short_name_data = select_metadata_by_name(input_file_info, short_name)
    do idat = 0, ListCount(short_name_data) - 1
      if (.not. isatt(short_name_data[idat], "reference_dataset")) then
        error_msg("f", DIAG_SCRIPT, "", "reference dataset not defined " + \
                  "for dataset " + short_name_data[idat]@filename)
      end if
      if (isdefined("short_name_ref")) then
        if (short_name_ref .ne. short_name_data[idat]@reference_dataset) then
          error_msg("f", DIAG_SCRIPT, "", "Differing reference datasets " + \
                    "variable " + short_name + ", got " + short_name_ref + \
                    " and " + short_name_data[idat]@reference_dataset)
        end if
      else
        short_name_ref = short_name_data[idat]@reference_dataset
        log_info("Using reference dataset " + short_name_ref + " for " + \
                 "variable " + short_name)
      end if
    end do
    refname(i) = short_name_ref
    delete(short_name_ref)
  end do
  dim_REF = dimsizes(refname)

  ; dataset indices with no reference dataset(s) and no multi-model mean

  model_ind = ispan(0, dim_MOD - 1, 1)
  n = dimsizes(refname)
  ref_ind = ispan(0, n - 1, 1)

  do i = 0, n - 1
    ri = ind(names .eq. refname(i))
    model_ind(ri) = -1
    ref_ind(i) = ri
  end do

  model_ind_woref = ind(model_ind.ge.0)
  if (flag_mean) then
    model_ind_worefmmm = model_ind_woref(ind(model_ind_woref.ne.mmm_ind))
  else
    model_ind_worefmmm = model_ind_woref
  end if

  delete(model_ind)

  ; ========================== initialization ==============================

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "calcmm", False)
  set_default_att(diag_script_info, "ecs_file", "ECS.nc")
  set_default_att(diag_script_info, "legend_outside", False)
  set_default_att(diag_script_info, "predef_minmax", False)
  set_default_att(diag_script_info, "output_diag_only", False)
  set_default_att(diag_script_info, "output_models_only", False)

  legend_outside = diag_script_info@legend_outside
  ecsfile = diag_script_info@ecs_file

  ; ============================ calculations ==============================

  diagdata = new(dim_MOD, float)
  diagminmax = new((/2/), float)
  debug_out = False

  ; southern ITCZ index (Tian, 2015)

  if (diag.eq."itczidx") then
    pr_idx = ind(variables.eq."pr")
    ref0 = ref_ind(pr_idx)

    result = itczidx(names(ref0), variables, debug_out, 0.0)

    SEPavg_ref = result[0]
    diagminmax = result[1]
    delete(result)

    do imod = 0, dim_MOD - 1
      ; skip multi-model mean

      if (imod.eq.mmm_ind) then
        continue
      end if

      result = itczidx(names(imod), variables, debug_out, SEPavg_ref)

      diagdata(imod) = result[0]
      if (.not.all(ismissing(result[1]))) then
        diagminmax = result[1]
      end if
      delete(result)
    end do

  end if  ; diag.eq."itczidx"

  ; tropical mid-tropopheric humidity asymmetry index (Tian, 2015)

  if (diag.eq."humidx") then
    sh_avg = new(dim_MOD, float)
    nh_avg = new(dim_MOD, float)

    hus_idx = ind(variables.eq."hus")
    ref0 = ref_ind(hus_idx)
    result = humidx(names(ref0), variables, debug_out, 1.0, 1.0)
    nh_avg(ref0) = result[1]
    sh_avg(ref0) = result[2]
    delete(result)

    do imod = 0, dim_MOD - 1
      ; skip multi-model mean entry

      if (imod.eq.mmm_ind) then
        continue
      end if

      result = humidx(names(imod), variables, debug_out, nh_avg(ref0), \
                      sh_avg(ref0))

      diagdata(imod) = result[0]

      if (.not.all(ismissing(result[3]))) then
        diagminmax = result[3]
      end if
      delete(result)

    end do  ; imod-loop
  end if  ; diag.eq."humidx"

  ; LTMI (Sherwood et al., 2014)

  if (diag.eq."ltmi") then
    obs_done = False
    datasets = new(dim_VAR, string)
    do imod = 0, dim_MOD - 1
      ; skip multi-model mean entry
      if (imod.eq.mmm_ind) then
        continue
      end if

      if (any(ref_ind.eq.imod)) then
        do_obs = True
        datasets = names(ref_ind)
        if (.not.isvar("ref0")) then
          ref0 = imod
        end if
      else
        do_obs = False
        datasets(:) = names(imod)
      end if

      ; all observations (possibly from different sources = different datasets)
      ; are processed at the same time --> skip all observations if
      ; observations have been processed already

      if (do_obs.and.obs_done) then
        continue
      end if

      diagdata(imod) = ltmi(datasets, variables, debug_out)

    end do  ; imod-loop

    diagdata@standard_name = "Lower tropospheric mixing index"
    diagdata@units = "1"
  end if

  ; Sherwood D index (Sherwood et al., 2014)

  if (diag.eq."sherwood_d") then
    obs_done = False
    datasets = new(dim_VAR, string)
    do imod = 0, dim_MOD - 1
      ; skip multi-model mean entry
      if (imod.eq.mmm_ind) then
        continue
      end if

      if (any(ref_ind.eq.imod)) then
        do_obs = True
        datasets = names(ref_ind)
        if (.not.isvar("ref0")) then
          ref0 = imod
        end if
      else
        do_obs = False
        datasets(:) = names(imod)
      end if

      ; all observations (possibly from different sources = different datasets)
      ; are processed at the same time --> skip all observations if
      ; observations have been processed already

      if (do_obs.and.obs_done) then
        continue
      end if

      diagdata(imod) = sherwood_d(datasets, variables, debug_out)

    end do  ; imod-loop

    diagdata@standard_name = "Sherwood_D_index"
    diagdata@units = "1"
  end if

  ; Sherwood S index (Sherwood et al., 2014)

  if (diag.eq."sherwood_s") then
    obs_done = False
    datasets = new(dim_VAR, string)
    do imod = 0, dim_MOD - 1
      ; skip multi-model mean entry
      if (imod.eq.mmm_ind) then
        continue
      end if

      if (any(ref_ind.eq.imod)) then
        do_obs = True
        datasets = names(ref_ind)
        if (.not.isvar("ref0")) then
          ref0 = imod
        end if
      else
        do_obs = False
        datasets(:) = names(imod)
      end if

      ; all observations (possibly from different sources = different datasets)
      ; are processed at the same time --> skip all observations if
      ; observations have been processed already

      if (do_obs.and.obs_done) then
        continue
      end if

      diagdata(imod) = sherwood_s(datasets, variables, debug_out)

    end do  ; imod-loop

    diagdata@standard_name = "Sherwood_S_index"
    diagdata@units = "1"
  end if

  ; covariance of shortwave cloud reflectivity (Brient and Schneider, 2016)

  if (diag.eq."covrefl") then
    obs_done = False
    datasets = new(dim_VAR, string)

    do imod = 0, dim_MOD - 1

      ; skip multi-model mean entry

      if (imod.eq.mmm_ind) then
        continue
      end if

      if (any(ref_ind.eq.imod)) then
        do_obs = True
        datasets = names(ref_ind)
        if (.not.isvar("ref0")) then
          ref0 = imod
        end if
      else
        do_obs = False
        datasets(:) = names(imod)
      end if

      ; all observations (possibly from different sources = different datasets)
      ; are processed at the same time --> skip all observations if
      ; observations have been processed already

      if (do_obs.and.obs_done) then
        continue
      end if

      diagdata(imod) = covrefl(datasets, variables, debug_out)

    end do

    diagdata@standard_name = "Regression sw cloud relectivity - sst"
    diagdata@units = "%/K"
  end if

  ; Southern Hemisphere Hadley cell extent (Lipat et al., 2017)

  if (diag.eq."shhc") then
    va_idx = ind(variables.eq."va")
    ref0 = ref_ind(va_idx)
    do imod = 0, dim_MOD - 1
      ; skip multi-model mean entry
      if (imod.eq.mmm_ind) then
        continue
      end if
      diagdata(imod) = shhc(names(imod), variables, debug_out)
    end do  ; imod-loop

    ; this is the "Southern Hemisphere Hadley cell extent"

    diagdata@standard_name = "SH Hadley cell extent"
    diagdata@units = "degrees"
  end if

  ; ==================== calculate MMM if desired ==================

  if (flag_mean) then
    diagdata(mmm_ind) = avg(diagdata(model_ind_worefmmm))
  end if

  ; ------------------------------------------------------------------------
  ; Skip further calculation/output if desired
  ; ------------------------------------------------------------------------

  if (diag_script_info@output_diag_only) then
    log_info("Returning only X axis of emergent constraint because " + \
             "'output_diag_only' is set to True")

    outfile = str_sub_str(diagdata@standard_name, " ", "_") + "_ref"
    n = dimsizes(ref_ind)
    do i = 0, n - 1
      outfile = outfile + "_" + names(ref_ind(i))
    end do
    if (isatt(diag_script_info, "suffix")) then
      outfile = outfile + "_" + diag_script_info@suffix
    end if

    diagdata!0 = "model"
    diagdata&model = names
    diagdata@var = diag
    diagdata@diag_script = DIAG_SCRIPT

    ; Additional attributes
    if (isdefined("output_attributes")) then
      atts = getvaratts(output_attributes)
      do iatt = 0, dimsizes(atts) - 1
        diagdata@$atts(iatt)$ = output_attributes@$atts(iatt)$
      end do
    end if

    ; Reference datasets
    do iref = 0,  dim_REF - 1
      if (.not. isvar("ref_str")) then
        ref_str = refname(iref)
      else
        ref_str = ref_str + "|" + refname(iref)
      end if
    end do
    if (isvar("ref_str")) then
      diagdata@reference_dataset = ref_str
    end if

    nc_filename = work_dir + outfile + ".nc"
    if (diag_script_info@output_models_only) then
      nc_outfile = ncdf_write(diagdata(model_ind_woref), nc_filename)
    else
      nc_outfile = ncdf_write(diagdata, nc_filename)
    end if

    plotfile = new(1, string)
    log_provenance(nc_outfile, plotfile, caption, "other", "reg", "", "", "", \
                   infiles)

    leave_msg(DIAG_SCRIPT, "")
    exit

  end if

  ; ==================== read precalculated ECS from file ==================

  if (isfilepresent(ecsfile)) then
    log_info("Read in: " + ecsfile)
    ecs_tmp = ncdf_read(ecsfile, "ecs")
  else
    error_msg("f", DIAG_SCRIPT, "", "file with precalculated ECS data (" + \
              ecsfile + ") not found.")
  end if

  ; now copy the precalculated ECS values read from the external file to the
  ; array ECS

  ECS = new(dim_MOD, float)
  ECS = ecs_tmp@_FillValue
  ECS@_FillValue = ecs_tmp@_FillValue
  ECS@units = ecs_tmp@units
  ECS@standard_name = ecs_tmp@standard_name

  exit_flag = False

  do i = 0, dim_MOD - 1
    if (any(i.eq.ref_ind).or.(i.eq.mmm_ind)) then
      continue
    end if
    name = names(i)
    idx = ind(name.eq.ecs_tmp&model)
    if (ismissing(idx)) then
      log_info("no pre-calcuated ECS value for model " + name + " available")
      exit_flag = True
    else
      ECS(i) = ecs_tmp(idx)
      log_info("pre-calcuated ECS value for model " + name + " = " + ECS(i))
    end if
    delete(idx)
  end do

  if (exit_flag) then
    error_msg("f", DIAG_SCRIPT, "", "precalculated ECS data (" + \
              ecsfile + ") not available for all models.")
  end if

  f = addfile(ecsfile, "r")
  climofiles = array_append_record(infiles, getfilepath(f), 0)

  ; ==================== calculate MMM if desired ==================

  if (flag_mean) then
    ECS(mmm_ind) = avg(ECS(model_ind_worefmmm))
  end if

  ; ============================ plot results ==============================

  if (config_user_info@write_plots.eq."True") then
    outfile = str_sub_str(diagdata@standard_name, " ", "_") + \
      "-" + str_sub_str(ECS@standard_name, " ", "_") + "_ref"
    n = dimsizes(ref_ind)
    do i = 0, n - 1
      outfile = outfile + "_" + names(ref_ind(i))
    end do

    if (isatt(diag_script_info, "suffix")) then
      outfile = outfile + "_" + diag_script_info@suffix
    end if

    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, outfile)
    wks@legendfile = outfile + "_legend"

    if (isatt(diag_script_info, "styleset")) then
      my_info = NewList("lifo")
      do i = 0, dim_MOD - 1
        ListAppend(my_info, new(1, logical))
        my_info[i]@dataset = names(i)
      end do
      colortab = project_style(my_info, diag_script_info, "colors")
      markertab = project_style(my_info, diag_script_info, "markers")
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

    nmod = dimsizes(model_ind_woref)

    colors  = new(nmod, string)
    markers = new(nmod, integer)

    ; create new marker: filled star

    mstring = "z"
    fontnum = 35
    size    = 1.5
    angle   = 0.0

    new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, 1.0, size, angle)

    if (isatt(diag_script_info, "styleset")) then
      colors = colortab(model_ind_woref)
      markers = markertab(model_ind_woref)
      i = ind(model_ind_woref.eq.mmm_ind)
      if (.not. all(ismissing(i))) then
        colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
        markers(i) = new_index
      end if
    else
      i = 0
      idx = 0
      do while (i.lt.nmod)
        imod = model_ind_woref(i)
        if (imod.eq.mmm_ind) then
          colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
          markers(i) = new_index
          i = i + 1
        else
          do n = 0, (nmod - 1) / dimsizes(colortab)
            colors(i)  = colortab(idx)
            markers(i) = markertab(n)
            i = i + 1
            if (i.ge.nmod) then
              break
            end if
          end do
          idx = idx + 1
          if (idx.ge.dimsizes(colortab)) then
            idx = dimsizes(colortab) - 1
          end if
        end if
      end do
    end if

    ; Collect data

    data_arr = new((/2, nmod/), float)
    data_arr!0 = "statistic"
    data_arr!1 = "model"
    data_arr&statistic = (/diagdata@standard_name, ECS@standard_name/)
    data_arr&model = names(model_ind_woref)
    data_arr(0, :) = (/diagdata(model_ind_woref)/)
    data_arr(1, :) = (/ECS(model_ind_woref)/)
    data_arr@units = (/diagdata@units, ECS@units/)
    data_arr@legend_outside = legend_outside
    data_arr@colors = colors
    data_arr@markers = markers
    data_arr@annots = names(model_ind_woref)

    yregmin = 999.9
    yregmax = -999.9

    xmin = min(data_arr(0, :))
    xmax = max(data_arr(0, :))
    xdelta = xmax - xmin

    if (diagdata(ref0) .lt. xmin) then
      data_arr@trXMinF = diagdata(ref0) - 0.02 * xdelta
    end if

    if (diagdata(ref0) .gt. xmax) then
      data_arr@trXMaxF = diagdata(ref0) + 0.02 * xdelta
    end if

    nx = dimsizes(diagdata(model_ind_worefmmm))
    ymin = min(ECS)
    ymax = max(ECS)

    xmin = min(diagdata)
    xmax = max(diagdata)
    xdelta = xmax - xmin
    x0 = xmin - 0.5 * xdelta
    x1 = xmax + 0.5 * xdelta

    if (nx.ge.3) then
      x = diagdata(model_ind_worefmmm)
      y = ECS(model_ind_worefmmm)
      ii = dim_pqsort_n(x, 1, 0)
      xx = x(ii)                ; temporary 'work' arrays
      yy = y(ii)

      rc = regline_stats(xx, yy)  ; requires NCL 6.2.0 or higher

      ; calculate confidence intervals (25%, 75%) of regression

      xxx = fspan(x0, x1, 50)
      yyy = rc@b(0) + xxx * rc@b(1)
      ; t-value given the one-sided probability and the degrees of freedom
      tval = cdft_t(0.25, nx - 2)

      mean_se_fit = 1.0 / nx + (xxx - rc@xave) ^ 2 / sum((x - rc@xave) ^ 2)
      mean_conf_upper = yyy + tval * sqrt(rc@MSE * mean_se_fit)
      mean_conf_lower = yyy - tval * sqrt(rc@MSE * mean_se_fit)

      ; calculate prediction intervals (25%, 75%)

      prediction_upper = yyy + tval * sqrt(rc@MSE * (1.0 + mean_se_fit))
      prediction_lower = yyy - tval * sqrt(rc@MSE * (1.0 + mean_se_fit))

      yregmin = min(mean_conf_lower)
      yregmax = max(mean_conf_upper)

      data_arr@trYMinF = min((/ymin, yregmin/))
      data_arr@trYMaxF = max((/ymax, yregmax/))

      delete(yyy)
      delete(mean_se_fit)
      delete(x)
      delete(y)
      delete(ii)
    end if

    ; If requested by user, use predefined min/max for x- and y-axis
    ; depending on diagnostic

    if (diag_script_info@predef_minmax) then
      data_arr@trYMinF = 2.0     ; CMIP3/CMIP5
      data_arr@trYMaxF = 6.5     ; CMIP3/CMIP5: 5.0, CMIP6: 6.5
      if (diag.eq."itczidx") then
        data_arr@trXMinF = -1.0    ; southern ITCZ index
        data_arr@trXMaxF = 3.0     ; southern ITCZ index
      end if
      if (diag.eq."humidx") then
        data_arr@trXMinF = -15.0   ; tropical mid-tropos. humidity asym. index
        data_arr@trXMaxF = 40.0    ; tropical mid-tropos. humidity asym. index
      end if
      if (diag.eq."sherwood_d") then
        data_arr@trXMinF = 0.0     ; Sherwood D index
        data_arr@trXMaxF = 0.6     ; Sherwood D index
      end if
      if (diag.eq."sherwood_s") then
        data_arr@trXMinF = 0.1     ; Sherwood S index
        data_arr@trXMaxF = 0.8     ; Sherwood S index
      end if
      if (diag.eq."ltmi") then
        data_arr@trXMinF = 0.4     ; lower tropospheric mixing index (LTMI)
        data_arr@trXMaxF = 1.2     ; lower tropospheric mixing index (LTMI)
      end if
      if (diag.eq."shhc") then
        data_arr@trXMinF = -40.0   ; SH HC extent
        data_arr@trXMaxF = -30.0   ; SH HC extent
      end if
      if (diag.eq."covrefl") then
        data_arr@trXMinF = -2.0    ; covariance of TLC sw reflection with SST
        data_arr@trXMaxF = 1.0     ; covariance of TLC sw reflection with SST
      end if

      ymin = min((/ymin, data_arr@trYMinF/))
      ymax = max((/ymax, data_arr@trYMaxF/))
    end if

    ; draw scatter plot
    data_arr@res_tiMainString = ""
    plot = scatterplot_markers(wks, data_arr, var0, input_file_info)

    ydelta = ymax - ymin
    y0 = min((/ymin - 0.5 * ydelta, yregmin/))
    y1 = max((/yregmax, ymax + 0.5 * ydelta/))

    ; plot observational uncertainty (if available)

    if (isdefined("diagminmax")) then
      x = (/diagminmax(0), diagminmax(0), diagminmax(1), diagminmax(1), \
           diagminmax(0)/)
      y = (/y0, y1, y1, y0, y0/)
      res = True
      res@tfPolyDrawOrder = "Draw"
      res@gsFillColor = (/0.9, 0.9, 0.9/)
      referr = gsn_add_polygon(wks, plot, x, y, res)
      delete(x)
      delete(y)
      delete(res)
    end if

    ; draw line for reference data

    x = (/diagdata(ref0), diagdata(ref0)/)
    y = (/y0, y1/)
    res = True
    res@gsLineColor = (/0.75, 0.75, 0.75/)
    res@gsLineThicknessF = 4.0
    res@tfPolyDrawOrder = "Draw"

    ref = gsn_add_polyline(wks, plot, x, y, res)

    delete(x)
    delete(y)

    ; add regression line (drawn in the background)
    ; (can only be calculated if there are at least 3 models)

    if (nx.ge.3) then
      x = (/x0, x1/)
      y = (/x0 * rc + rc@yintercept, x1 * rc + rc@yintercept/)
      res = True
      res@gsLineColor = (/1.00, 0.00, 0.00/)
      res@gsLineThicknessF = 4.0
      res@tfPolyDrawOrder = "Draw"
      regl = gsn_add_polyline(wks, plot, x, y, res)

      res@gsLineThicknessF = 2.0
      res@gsLineDashPattern = 1
      reglc1 = gsn_add_polyline(wks, plot, xxx, mean_conf_upper, res)
      reglc2 = gsn_add_polyline(wks, plot, xxx, mean_conf_lower, res)
      res@gsLineDashPattern = 2
      reglp1 = gsn_add_polyline(wks, plot, xxx, prediction_upper, res)
      reglp2 = gsn_add_polyline(wks, plot, xxx, prediction_lower, res)

      delete(xxx)
      delete(mean_conf_upper)
      delete(mean_conf_lower)
      delete(prediction_upper)
      delete(prediction_lower)
    end if

    draw(plot)
    frame(wks)

    plotfile = wks@fullname

    log_info("Wrote " + plotfile)
  else
    plotfile = ""
  end if  ; if config_user_info@write_plots.eq."True"

  ; ======================== netcdf output ===========================

  outfile = str_sub_str(diagdata@standard_name, " ", "_") + "_ref"
  n = dimsizes(ref_ind)
  do i = 0, n - 1
    outfile = outfile + "_" + names(ref_ind(i))
  end do
  if (isatt(diag_script_info, "suffix")) then
    outfile = outfile + "_" + diag_script_info@suffix
  end if

  diagdata!0 = "model"
  diagdata&model = names
  diagdata@var = diag
  diagdata@diag_script = DIAG_SCRIPT

  ; Additional attributes
  if (isdefined("output_attributes")) then
    atts = getvaratts(output_attributes)
    do iatt = 0, dimsizes(atts) - 1
      diagdata@$atts(iatt)$ = output_attributes@$atts(iatt)$
    end do
  end if

  ; Reference datasets
  do iref = 0,  dim_REF - 1
    if (.not. isvar("ref_str")) then
      ref_str = refname(iref)
    else
      ref_str = ref_str + "|" + refname(iref)
    end if
  end do
  if (isvar("ref_str")) then
    diagdata@reference_dataset = ref_str
  end if

  nc_filename = work_dir + outfile + ".nc"
  if (diag_script_info@output_models_only) then
    nc_outfile = ncdf_write(diagdata(model_ind_woref), nc_filename)
  else
    nc_outfile = ncdf_write(diagdata, nc_filename)
  end if

  nc_filename@existing = "append"
  ECS!0 = "model"
  ECS&model = names
  ECS@var = "ECS"
  ECS@diag_script = DIAG_SCRIPT

  if (diag_script_info@output_models_only) then
    nc_outfile = ncdf_write(ECS(model_ind_woref), nc_filename)
  else
    nc_outfile = ncdf_write(ECS, nc_filename)
  end if

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file
  ; ------------------------------------------------------------------------

  statistics = ("other")
  domain = ("reg")
  plottype = ("scatter")

  log_provenance(nc_outfile, plotfile, caption, statistics, domain, \
                 plottype, "", "", climofiles)

  leave_msg(DIAG_SCRIPT, "")

end
