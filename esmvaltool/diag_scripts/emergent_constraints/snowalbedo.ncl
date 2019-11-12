; #############################################################################
; snowalbedo (ipcc ar5 fig. 9.45a)
; #############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME CRESCENDO
; #############################################################################
; Description
;   Calculates springtime snow-albedo feedback values in climate change
;   versus springtime values in the seasonal cycle in transient climate
;   change experiments.
;   Figure resembles IPCC AR5 Fig. 9.45a (Flato et al., 2013).
;
;   References:
;     - Flato et al., Evaluation of climate models, in: Climate Change 2013:
;       the Physical Science Basis, 2013.
;     - Hall, A., and X. Qu, Using the current seasonal cycle to constrain
;       snow albedo feedback in future climate change, Geophys. Res. Lett.,
;       33, L03502, doi:10.1029/2005GL025127, 2006.
;
; Required settings (scripts)
;   exp_presentday:  name of present-day experiment (e.g. "historical")
;   exp_future:      name of climate change experiment (e.g. "rcp45")
;
; Optional settings (scripts)
;   diagminmax:      observational uncertainty (min and max)
;   legend_outside:  True, False
;   styleset:        "CMIP5" (if not set, this diagnostic will create its own
;                    color table and symbols for plotting)
;   suffix:          string to be added to output filenames
;   xmax:            upper limit of x-axis (default = automatic)
;   xmin:            lower limit of x-axis (default = automatic)
;   ymax:            upper limit of y-axis (default = automatic)
;   ymin:            lower limit of y-axis (default = automatic)
;
; Required settings (variables)
;   ref_model:       name of reference data set
;
; Optional settings (variables)
;   none
;
; Caveats
;   - requires NCL version 6.4 or higher to plot confidence bands of
;     regression
;   - no regridding: rsdscs and albsfc fields are required to be on the same
;     grid
;   - it remains unclear which experiments and time periods have to be used
;     to exactly reproduce fig. 9.45a in IPCC AR5
;
; Modification history
;   20191104-lauer_axel: rewritten for v2.0
;   20171220-lauer_axel: written
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"

; =============================================================================
; Function snowalbedo_precalc
; Calculate variables needed to obtain springtime snow-albedo feedback values:
;   - seasonal cycle of surface albedo weighted by April incoming solar
;     radiation and averaged over the region 30N-90N, 0-360 (land grid
;     cells only)
;   - seasonal cycle of 2m temperature averaged over the region 30N-90N, 0-360
;     (land grid cells only)
;
; Arguments
; - tas: 2-m temperature
; - alb: surface albedo
; - rad: downwelling shortwave radiation at surface
; - lsdata: NCL build-in land-sea mask
;
; Reference
; - Hall et al. (2006)
; =============================================================================

function snowalbedo_precalc(tas, alb, rad, lsdata)
begin
  ; check that dimensions of "rsdt" and "alb" are identical

  dims_r = dimsizes(rad)
  dims_a = dimsizes(alb)

  dimensions_fine = False

  if (dimsizes(dims_r).eq.dimsizes(dims_a)) then
    if (all(dims_r.eq.dims_a)) then
      dimensions_fine = True
    end if
  end if

  delete(dims_r)
  delete(dims_a)

  if (.not.dimensions_fine) then
    error_msg("f", DIAG_SCRIPT, "", "dimensions of variables alb and rsdt " \
              + "do not match.")
  end if

  rad_seas = time_operations(rad, -1, -1, "average", "monthlyclim", True)
  rad_apr = rad_seas(3, :, :)
  delete(rad_seas)

  ; -------------------------------------------------------------------

  ; generate and apply land-mask

  lsm0 = landsea_mask(lsdata, alb&lat, alb&lon)
  lsm  = where(lsm0.eq.1, 0, lsm0@_FillValue)
  copy_VarMeta(lsm0, lsm)
  delete(lsm0)
  alb = alb + conform(alb, lsm, (/1, 2/))

  ; ---------------------------------------------------------------------
  ; weight surface albedo by April incoming solar radiation

  ; apply land-mask to April incoming solar radiation
  rad_apr = rad_apr + lsm

  ; extract region 30N-90N, 0-360
  rad_reg = area_operations(rad_apr, 30., 90., 0., 360., "extract", False)
  delete(rad_apr)
  alb_reg = area_operations(alb, 30., 90., 0., 360., "extract", False)

  ; calculate grid box areas
  area_reg = map_area(rad_reg&lat, rad_reg&lon)

  ; calculate weights
  weights = rad_reg * area_reg / sum(rad_reg * area_reg)
  delete(area_reg)
  delete(rad_reg)

  ; calculate weighted mean albedo

  albmean = dim_sum_n_Wrap(alb_reg * conform(alb_reg, weights, (/1, 2/)), \
                           (/1, 2/))
  albmean!0 = "time"
  albmean&time = alb_reg&time
  delete(weights)
  delete(alb_reg)

  ; ---------------------------------------------------------------------

  delete(lsm)

  lsm0  = landsea_mask(lsdata, tas&lat, tas&lon)
  lsm = where(lsm0.eq.1, 0, lsm0@_FillValue)
  copy_VarMeta(lsm0, lsm)
  delete(lsm0)
  tas = tas + conform(tas, lsm, (/1, 2/))
  delete(lsm)

  ; calculate average temperature (missing values will be ignored and grid
  ; boxes area-weighted) over NH land grid cells (polewards of 30 deg N)
  t2  = area_operations(tas, 30., 90., 0., 360., "average", True)

  alb_seas = time_operations(albmean, -1, -1, "average", "monthlyclim", True)
  t2_seas  = time_operations(t2, -1, -1, "average", "monthlyclim", True)

  delete(t2)

  return([/alb_seas, t2_seas/])

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

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT)
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  variables = metadata_att_as_array(variable_info, "short_name")
  dim_VAR = dimsizes(variables)

  tmp_dim_MOD = new(dim_VAR, integer)
  do i = 0, dim_VAR - 1
    var0 = variable_info[i]@short_name
    info0 = select_metadata_by_name(input_file_info, var0)
    tmp_dim_MOD(i) = ListCount(info0)
    delete(info0)
  end do

  dim_MOD = tmp_dim_MOD(0)

  if (any(tmp_dim_MOD.ne.dim_MOD)) then
    error_msg("f", DIAG_SCRIPT, "", "number of datasets has to be " \
              + "identical for all variables.")
  end if

  names = new((/dim_VAR, dim_MOD/), string)
  exps = new((/dim_VAR, dim_MOD/), string)
  ensembles = new((/dim_VAR, dim_MOD/), string)

  do i = 0, dim_VAR - 1
    var0 = variable_info[i]@short_name
    info0 = select_metadata_by_name(input_file_info, var0)
    names(i, :) = metadata_att_as_array(info0, "dataset")
    exps(i, :) = metadata_att_as_array(info0, "exp")
    ensembles(i, :) = metadata_att_as_array(info0, "ensemble")
    delete(info0)
  end do

  ; check if required variables are available

  tas_idx = ind(variables.eq."tas")
  rad_idx = ind(variables.eq."rsdt")
  alb_idx = ind(variables.eq."alb")

  if (ismissing(tas_idx)) then
    error_msg("f", DIAG_SCRIPT, "", "variable tas required but not found.")
  end if
  if (ismissing(alb_idx)) then
    error_msg("f", DIAG_SCRIPT, "", "variable albsfc required but not found.")
  end if
  if (ismissing(rad_idx)) then
    error_msg("f", DIAG_SCRIPT, "", "variable rsdt required but not found.")
  end if

  ; check if present-day and future experiment names are defined

  if (isatt(diag_script_info, "exp_presentday")) then
    exp_presentday = diag_script_info@exp_presentday
  else
    error_msg("f", DIAG_SCRIPT, "", "'exp_presentday' not defined in recipe" \
              + " (scripts).")
  end if

  if (isatt(diag_script_info, "exp_future")) then
    exp_future = diag_script_info@exp_future
  else
    error_msg("f", DIAG_SCRIPT, "", "'exp_future' not defined in recipe" \
              + " (scripts).")
  end if

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "legend_outside", False)

  ; reference datasets

  refname = new(dim_VAR, string)
  ref_ind = new(dim_VAR, integer)
  do i = 0, dim_VAR - 1
    if (isatt(variable_info[i], "reference_dataset")) then
      refname(i) = variable_info[i]@reference_dataset
      ref_ind(i) = ind(names(i, :).eq.refname(i))
    else
      error_msg("f", DIAG_SCRIPT, "", "reference dataset not defined for " \
                + "variable " + variables(i))
    end if
  end do

  ; check that each model provides two datasets: exp_presentday and exp_future

  ; model indices with no reference dataset

  do i = 0, dim_VAR - 1
    ind_tmp = ind(exps(i, :).eq.exp_presentday)
    ind_tmp = where(ind_tmp.eq.ref_ind(i), -1, ind_tmp)
    idx = ind(ind_tmp.ge.0)
    ind_present = ind_tmp(idx)
    delete(ind_tmp)
    delete(idx)
    ind_tmp = ind(exps(i, :).eq.exp_future)
    ind_tmp = where(ind_tmp.eq.ref_ind(i), -1, ind_tmp)
    idx = ind(ind_tmp.ge.0)
    ind_future = ind_tmp(idx)
    delete(ind_tmp)
    delete(idx)

    if ((dimsizes(ind_present).ne.dim_MOD / 2) .or. \
        (dimsizes(ind_future).ne.dim_MOD / 2)) then
      error_msg("f", DIAG_SCRIPT, "", "variable " + variables(i) \
                + ": number of model datasets for exp " + exp_presentday \
                + " and exp " + exp_future + " is not identical.")
    end if

    do j = 0, dimsizes(ind_present) - 1
      itest = ind(names(i, ind_future).eq.names(i, ind_present(j)))
      if (ismissing(itest)) then
        error_msg("f", DIAG_SCRIPT, "", "variable " + variables(i) \
                  + ", dataset " + names(i, ind_present(j)) \
                  + ": no data for experiment " + exp_future + " found.")
      end if
    end do
  end do

  num_models_present = dimsizes(ind_present)
  delete(ind_future)

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  climofiles = metadata_att_as_array(input_file_info, "filename")
end

begin

  ; ========================================================================
  ; ============================ calculations ==============================
  ; ========================================================================

  ydata = new(num_models_present + 1, float)
  xdata = new(num_models_present + 1, float)
  ds_names = new(num_models_present + 1, string)

  obs_done = False

  a = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
  lsdata = a->LSMASK

  do ii = 0, num_models_present

    ; --------------------------------------------------
    ; calculate seasonal cycle delta_albedo_s / delta_T_s
    ; (differences between April and May)
    ; unit = %/K
    ; --------------------------------------------------

    do jj = 0, dim_VAR - 1
      atts = True
      atts@short_name = variables(jj)

      if (ii .eq. num_models_present) then
        ; all observations (possibly from different sources)
        ; are processed at the same time (at ii .eq. num_models_present)

        atts@dataset = refname(jj)
        ds_names(ii) = "ref"
      else  ; model dataset
        imod = ind_present(ii)
        atts@dataset = names(0, imod)
        atts@exp = exp_presentday
        atts@ensemble = ensembles(0, imod)
        ds_names(ii) = names(0, imod) + "_" + ensembles(0, imod)
      end if

      info = select_metadata_by_atts(input_file_info, atts)
      delete(atts)

      if (variables(jj) .eq. "tas") then
        tas = read_data(info[0])
      end if
      if (variables(jj) .eq. "alb") then
        alb = read_data(info[0])
      end if
      if (variables(jj) .eq. "rsdt") then
        rad = read_data(info[0])
      end if

      delete(info)
    end do

    result = snowalbedo_precalc(tas, alb, rad, lsdata)
    alb_seas = result[0]
    t2_seas = result[1]
    delete(result)

    xdata(ii) = (alb_seas(3) - alb_seas(4)) * 100.0 / (t2_seas(3) - t2_seas(4))

    delete(tas)
    delete(alb)
    delete(rad)

    ; --------------------------------------------------
    ; calculate climate change delta_alpha_s / delta_T_s
    ; note: models only, no observations
    ; --------------------------------------------------

    if (ii .lt. num_models_present) then  ; i.e. model dataset
      do jj = 0, dim_VAR - 1
        atts = True
        atts@short_name = variables(jj)
        imod = ind_present(ii)
        atts@dataset = names(0, imod)
        atts@exp = exp_future
        atts@ensemble = ensembles(0, imod)

        info = select_metadata_by_atts(input_file_info, atts)
        delete(atts)

        if (variables(jj) .eq. "tas") then
          tas = read_data(info[0])
        end if
        if (variables(jj) .eq. "alb") then
          alb = read_data(info[0])
        end if
        if (variables(jj) .eq. "rsdt") then
          rad = read_data(info[0])
        end if

        delete(info)
      end do

      result = snowalbedo_precalc(tas, alb, rad, lsdata)

      delete(tas)
      delete(alb)
      delete(rad)

      alb_seas_future = result[0]
      t2_seas_future = result[1]
      delete(result)

      delta_t = t2_seas_future(3) - t2_seas(3)
      if (delta_t.ne.0.0) then
        ydata(ii) = (alb_seas_future(3) - alb_seas(3)) * 100.0 / delta_t
      else
        ydata(ii) = ydata@_FillValue
      end if

      delete(alb_seas_future)
      delete(t2_seas_future)
    end if

    delete(alb_seas)
    delete(t2_seas)

    print("-------------------------------------------")
    print(ds_names(ii) + ": " + xdata(ii) + "%/K, " + ydata(ii) + "%/K")

  end do  ; loop over datasets

  delete(lsdata)

  xdata@standard_name = "seasonal cycle"
  ydata@standard_name = "climate change"
  xdata@units = "%/K"
  ydata@units = "%/K"

  ; ========================================================================
  ; ============================ plot results ==============================
  ; ========================================================================

  outfile = "snow_albedo_feedback_ref"
  do i = 0, dim_VAR - 1
    outfile = outfile + "_" + names(i, ref_ind(i))
  end do

  if (isatt(diag_script_info, "suffix")) then
    outfile = outfile + "_" + diag_script_info@suffix
  end if

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, outfile)
  plot_filename = wks@fullname
  wks@legendfile = outfile + "_legend"

  if (isatt(diag_script_info, "styleset")) then
    info0 = select_metadata_by_name(input_file_info, var0)
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

  nx = dimsizes(xdata)

  colors  = new(nx, string)
  markers = new(nx, integer)

  ; create new marker: filled star

  mstring = "z"
  fontnum = 35
  size    = 1.5
  angle   = 0.0

  new_index = NhlNewMarker(wks, mstring, fontnum, 0.0, 0.0, 1.0, size, angle)

  if (isatt(diag_script_info, "styleset")) then
    colors(0:nx - 2) = colortab(ind_present)
    markers(0:nx - 2) = markertab(ind_present)
    i = ind(ds_names.eq."MultiModelMean")
    if (.not. ismissing(i)) then
      colors(i)  = "(/0.00, 0.00, 0.00/)"  ; black
      markers(i) = new_index
    end if
  else
    i = 0
    j = 0
    do n = 0, nx - 2
      if (ds_names(n).eq."MultiModelMean") then
        colors(n)  = "(/0.00, 0.00, 0.00/)"  ; black
        markers(n) = new_index
      else
        colors(n) = colortab(i)
        markers(n) = markertab(j)
        i = i + 1
        if (i.ge.dimsizes(colortab)) then
          i = 0
          j = j + 1
          if (j.ge.dimsizes(markertab)) then
            log_info(DIAG_SCRIPT + " (info): not enough colors/symbols " \
                     + "for all datasets, reusing symbols/colors")
            j = 0
          end if
        end if
      end if
    end do
  end if

  ; reference data: red star
  colors(nx - 1) = "(/1.00, 0.00, 0.00/)"  ; red
  markers(nx - 1) = new_index

  ; Collect data (annotations, colors and markers are defined for
  ; model datasets only)

  data_arr = new((/2, nx/), float)
  data_arr!0 = "statistic"
  data_arr!1 = "model"
  data_arr&statistic = (/xdata@standard_name, ydata@standard_name/)
  data_arr&model = ds_names

  data_arr(0, :) = (/xdata/)
  data_arr(1, :) = (/ydata/)
  data_arr@units = (/xdata@units, ydata@units/)
  data_arr@legend_outside = diag_script_info@legend_outside
  data_arr@colors = colors(0: nx - 2)
  data_arr@markers = markers(0: nx - 2)
  data_arr@annots = ds_names(0: nx - 2)

  yregmin = 999.9
  yregmax = -999.9

  ymin = min(ydata)
  ymax = max(ydata)

  xmin = min(xdata)
  xmax = max(xdata)
  xdelta = xmax - xmin
  x0 = xmin - 0.5 * xdelta
  x1 = xmax + 0.5 * xdelta

  delete(ii)

  if (nx.ge.3) then
    x = xdata(0:nx - 2)  ; w/o reference dataset
    y = ydata(0:nx - 2)  ; w/o reference dataset
    ii = dim_pqsort_n(x, 1, 0)
    xx = x(ii)           ; temporary 'work' arrays
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

  ; if requested by user, use predefined min/max for x- and y-axis

  if (isatt(diag_script_info, "xmax")) then
    data_arr@trXMaxF = -0.3
  end if
  if (isatt(diag_script_info, "xmin")) then
    data_arr@trXMinF = -1.7
  end if
  if (isatt(diag_script_info, "ymax")) then
    data_arr@trYMaxF = -0.3
  end if
  if (isatt(diag_script_info, "ymin")) then
    data_arr@trYMinF = -1.7
  end if

  data_arr@trXReverse = True  ; flip x axis
  data_arr@trYReverse = True  ; flip y axis

  ; draw scatter plot
  data_arr@res_tiMainString = ""
  plot = scatterplot_markers(wks, data_arr(:, 0:nx - 2), var0, input_file_info)

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

  x = (/xdata(nx - 1), xdata(nx - 1)/)
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

  ; mandatory netCDF output

  nc_filename = work_dir + outfile + ".nc"

  data_arr@diag_script = (/DIAG_SCRIPT/)
  data_arr@var = "snowalbedofeedback"
  nc_outfile = ncdf_write(data_arr, nc_filename)

  ; ------------------------------------------------------------------------
  ; write provenance to netCDF output (and plot file)
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "other"/)
  domain = "nhext"
  plottype = "scatter"
  caption = "Scatterplot of simulated springtime snow-albedo feedback " \
            + "values versus springtime values in the seasonal cycle. " \
            + "Resembles IPCC AR5 fig. 9.45a."

  log_provenance(nc_outfile, plot_filename, caption, statistics, \
                 domain, plottype, "", "", climofiles)

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
