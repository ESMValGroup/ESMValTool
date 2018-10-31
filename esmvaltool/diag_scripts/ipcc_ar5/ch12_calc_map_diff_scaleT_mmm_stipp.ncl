; #############################################################################
; Multi-model mean change scaled by global T change map with significance
; Author: Ruth Lorenz (ETH, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Description
;    Pattern scaling. Annual mean or seasonal change scaled by global T
;    change per model multi-model mean with significance in different
;    periods (e.g. 2081-2100 and 2181-2200 with respect to 1986-2005)
;
; Required diag_script_info attributes (diagnostics specific)
;    scenarios: list with scenarios to be included in the
;               figure, e.g (/"rcp26","rcp45","rcp60","rcp85"/)
;    periods: list with start years of periods to be included
;    	       e.g. (/"2046","2081","2180"/)
;    time_avg: list with seasons or annual to be included
;              e.g (/"annualclim"/)
;
; Optional diag_script_info attributes (diagnostic specific)
;    seasons: list with seasons index to be included if
;             averagetime is "seasonalclim" (required for
;             seasonalclim), DJF:0, MAM:1, JJA:2, SON:3
;    percent: 0 or 1, difference in percent = 1, default = 0
;    explicit_cn_levels: explicit levels for all contour plots
;    max_vert: maximum number of plots in vertical
;              default determined by number of scenarios
;    max_hori: maximum number of plots in horizontal
;              default determined by number of periods
;    title: Figure title
;    label: labels on top of each panel, e.g. scenario + rcp,
;           loop to plot is seasons, scenarios, periods
;    colormap: specify colormap for contour plots
;    units: unit string in ncl formatting for legend title
;
; Caveats
;    Only handles one variable at a time
;
; Modification history
;    20180618-A_lore_ru: refactored for v2, fixed calculation of stippling
;    20161027-A_lore_ru: written for ESMValTool
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/ensemble.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ; 'diag_script' and 'variables' are fetched from the
  ; above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; 'models', 'variables' and 'field_types' are fetched from the
  ; above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  if (dim_VAR .eq. 2) then
    var0 = variables(0)
    var1 = variables(1)
  else
    var0 = variables(0)
    var1 = variables(0)
  end if
  field_type0 = field_types(0)
  if (isatt(diag_script_info, "seasons")) then
    dim_seas = dimsizes(diag_script_info@seasons)
  else
    dim_seas = 1
  end if
  ; ###########################################
  ; # References Information                  #
  ; ###########################################
  write_references(diag_script, \  ; script name
                   "A_lore_ru", \  ; authors
                   "", \  ; contributors
                   "D_collins13ipcc",    \  ; diag_references
                   "",          \  ; obs_references
                   "P_crescendo")   ; proj_references
  ; Check required diag_script_info attributes
  req_atts = (/"scenarios", "periods", "time_avg"/)
  exit_if_missing_atts(diag_script_info, req_atts)

end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if(isvar("MyParam")) then
      delete(MyParam)
  end if
  if(isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; get data from first model
  imod = 0  ; NCL array indices start from zero
  log_info("processing " + dataset_info@dataset(imod) + "_"\
           + dataset_info@exp(imod) + "_"\
           + dataset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)
  Fill = default_fillvalue(typeof(A0))  ; 1e20

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  nlat = dimsizes(A0&lat)
  idx = ind(dims .eq. "lon")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lon dimension")
  end if
  nlon = dimsizes(A0&lon)

  temp_avg = new((/dim_MOD/), typeof(A0), Fill)
  temp_avg!0 = "models"
  data1 = new((/dim_MOD,  dim_seas, nlat, nlon/), \
              typeof(A0), Fill)
  ; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
  if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
      (dim_seas .ne. 4)) then
    tmp_seas = time_operations(A0, -1, -1, "average", \
                               diag_script_info@time_avg, True)
    if (dim_seas .eq. 1) then
      data1(imod, 0, :, :) = tmp_seas(diag_script_info@seasons, :, :)
    else
      data1(imod, :, :, :) = tmp_seas(diag_script_info@seasons, :, :)
    end if
  else if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
           (dim_seas .eq. 4)) then
    data1(imod, :, :, :) = time_operations(A0, -1, -1, "average", \
                                           diag_script_info@time_avg, True)
  else
    data1(imod, 0, :, :) = time_operations(A0, -1, -1, "average", \
                                           diag_script_info@time_avg, True)
  end if
  end if
  rad = (4.0 * atan(1.0) / 180.0)
  do imod = 0, dim_MOD - 1
    A0_temp = read_data(imod, var1, field_type0)
    temp_tavg = dim_avg_n_Wrap(A0_temp, 0)
    latw = cos(data1&lat * rad)
    temp_avg(imod) = wgt_areaave_Wrap(temp_tavg, latw, 1.0, 0)
    delete(temp_tavg)
    delete(latw)
    if (imod .ne. 0) then
      A0_temp = read_data(imod, var1, field_type0)
      temp_tavg = dim_avg_n_Wrap(A0_temp, 0)
      latw = cos(temp_tavg&lat * rad)
      temp_avg(imod) = wgt_areaave_Wrap(temp_tavg, latw, 1.0, 0)
      A0 = read_data(imod, var0, field_type0)
      ; Average over time,
      ; function in ./diag_scripts/lib/ncl/statistics.ncl)
      if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
          (dim_seas .ne. 4)) then
        tmp_seas = time_operations(A0, -1, -1, "average", \
                                   diag_script_info@time_avg, True)
        if (dim_seas .eq. 1) then
          data1(imod, 0, :, :) = tmp_seas(diag_script_info@seasons, :, :)
        else
          data1(imod, :, :, :) = tmp_seas(diag_script_info@seasons, :, :)
        end if
        delete(tmp_seas)
      else if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
               (dim_seas .eq. 4)) then
        data1(imod, :, :, :) = time_operations(A0, -1, -1, "average", \
                                               diag_script_info@time_avg, True)
      else
        data1(imod, 0, :, :) = time_operations(A0, -1, -1, "average", \
                                               diag_script_info@time_avg, True)
      end if
      end if
    end if
    delete(A0_temp)
    delete(A0)
  end do
  delete(imod)

  ; ##############################################
  ; # Calculate change from hist to periods      #
  ; ##############################################
  ; loop over rcps and periods
  dim_scen = dimsizes(diag_script_info@scenarios)
  dim_per = dimsizes(diag_script_info@periods)

  to_plot = new((/dim_per * dim_seas, nlat, nlon/), typeof(data1), Fill)
  if (isatt(diag_script_info, "percent")) then
    percent = diag_script_info@percent
  else
    percent = 0
  end if
  to_plot_signif = new((/dim_per * dim_seas, nlat, nlon/), typeof(data1), Fill)
  perc95 = new((/dim_per * dim_seas,  nlat, nlon/), typeof(data1), Fill)
  i = 0
  do seas = 0, dim_seas - 1  ; loop over seasons
    do per = 0, dim_per - 1  ; loop over periods
      do rcp = 0, dim_scen - 1  ; loop over rcps
        idx_rcp = ind((dataset_info@exp .eq. diag_script_info@scenarios(rcp)) \
                      .and. (dataset_info@start_year .eq. \
                      diag_script_info@periods(per)))
        proj = data1(idx_rcp, seas, :, :)
        temp_proj = temp_avg(idx_rcp)
        ; **cut low values ***
        proj = where(abs(proj) .gt. 1.e-14, proj, proj@_FillValue)
        ; find historical runs from same models as in rcp
        idx_hist = find_model_index(dataset_info@dataset(idx_rcp), \
                                    dataset_info@dataset, "historical", \
                                    dataset_info@exp, \
                                    dataset_info@ensemble(idx_rcp), \
                                    dataset_info@ensemble, 0)
        base = data1(idx_hist, seas, :, :)
        temp_base = temp_avg(idx_hist)
        ; **cut low values ***
        base = where(abs(base) .gt. 1.e-14, base, base@_FillValue)
        ; scale each model by global T change
        log_debug("Scale each model by global T change")
        dim_mod = dimsizes(idx_rcp)
        if (dim_mod .lt. 2) then
          error_msg("f", diag_script, "", \
                    "Only one model found, multi-model needs at least two.")
        end if
        var_diff_scal = new((/dim_mod, nlat, nlon/), typeof(data1), Fill)
        do imod = 0, dim_mod - 1
          if (percent .eq. 1) then
            var_diff_scal(imod, :, :) = \
              ((100 * (proj(imod, :, :) - base(imod, :, :))) / \
               where(base(imod, :, :) .ne. 0.,\
               base(imod, :, :), base@_FillValue)) / \
              (temp_proj(imod) - temp_base(imod))
          else
            var_diff_scal(imod, :, :) = (proj(imod, :, :) - \
                                         base(imod, :, :)) / \
                                        (temp_proj(imod) - temp_base(imod))
          end if
        end do
        delete(base)
        delete(proj)
        delete(temp_proj)
        delete(temp_base)
        if rcp .eq. 0 then
          var_diff_scal_all_mod = var_diff_scal
        else
          tmp_scal_all_mod = var_diff_scal_all_mod
          delete(var_diff_scal_all_mod)
          var_diff_scal_all_mod = \
            array_append_record(tmp_scal_all_mod, var_diff_scal, 0)
          delete(tmp_scal_all_mod)
        end if
        delete(idx_hist)
        delete(idx_rcp)
        delete(var_diff_scal)
      end do  ; rcp
      var_diff_scal_all_mod!0 = "models"
      var_diff_scal_all_mod!1 = "lat"
      var_diff_scal_all_mod!2 = "lon"
      var_diff_scal_all_mod&lat = data1&lat
      var_diff_scal_all_mod&lon = data1&lon
      ; average over rcps
      log_debug("Average over models and rcps")
      to_plot(i, :, :) = dim_avg_n_Wrap(var_diff_scal_all_mod, 0)
      log_debug("determine significance")
      ; calculate standard deviation and 95-percentile (assuming gaussian
      ; distribution) over all models and scenarios
      nvalid = dim_num_n(.not.ismissing(var_diff_scal_all_mod(:, nlat / 2, \
                                                              nlon / 2)), 0)
      printVarSummary(nvalid)
      print(nvalid)
      perc95(i, :, :) = dim_stddev_n_Wrap(var_diff_scal_all_mod, 0) * 1.96 / \
        sqrt(nvalid)  ; standard error
      ; compare change to perc95
      to_plot_signif(i, :, :) = where(abs(to_plot(i, :, :)) .gt. \
                                      abs(perc95(i, :, :)), 1., 0)
      i = i + 1
      delete([/var_diff_scal_all_mod/])
    end do  ; per
  end do  ; seas
  to_plot!0 = "panel"
  to_plot&panel = diag_script_info@label
  to_plot!1 = "lat"
  to_plot&lat = data1&lat
  to_plot!2 = "lon"
  to_plot&lon = data1&lon

  if (percent .eq. 1) then
    to_plot@units     = "% per ~F35~J~F~C"
  end if

  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if (isatt(to_plot, "diag_script")) then  ; Add to existing entries
    tmp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(tmp, (/diag_script/), 0)
    delete(tmp)
  else  ; Add as new attribute
    to_plot@diag_script = (/diag_script/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if (isatt(variable_info, "long_name")) then
    to_plot@var_long_name = variable_info@long_name
  end if
  ; check unit and adjust for plotting if necessary
  if (isatt(variable_info, "units")) then
    if ((var0 .eq. "mrro" .or. var0 .eq. "evspsbl") \
        .and. variable_info@units .eq. "kg m-2 s-1") then
      ; convert [kg m-2 s-1] to [mm day-1]
      to_plot = to_plot * 60 * 60 * 24
      variable_info@units = "mm day-1"
      to_plot@var_units   = variable_info@units
    else if (var0 .eq. "mrsos" .and. variable_info@units .eq. "kg m-2") then
      ; convert [kg m-2] to [mm]
      variable_info@units = "mm"
      to_plot@var_units   = variable_info@units
    else if (var0 .eq. "hfls" .and. variable_info@units .eq. "W m-2") then
      ; convert [W m-2] to [mm day-1]
      ; ET[kg m**-2 s**-1] = LH / lambda
      ;                    = LE / (2.5 * 10 ^ 6 [J kg**-1])
      ; [kg m**-2 s**-1] * 24 * 60 * 60 = [mm / day]
      to_plot = to_plot / (2.5 * 10 ^ 6) * 24 * 60 * 60
      variable_info@units = "mm day-1"
      to_plot@var_units   = variable_info@units
    else if (variable_info@units .eq. "Pa") then
      ; convert to hPa
      to_plot = to_plot / 100
      to_plot@var_units     = "hPa"
    else
      to_plot@var_units     = variable_info@units
    end if
    end if
    end if
    end if
  end if
  ; swap downward/upward radiation fluxes according to IPCC Chapter12
  if ((var0 .eq. "rsut") .or. (var0 .eq. "rlut")) then
    to_plot = - to_plot
  end if
  copy_VarCoords(to_plot, to_plot_signif)

  to_plot_signif@diag_script = (/diag_script/)
  to_plot_signif@var = "signif"
  to_plot_signif@var_long_name = "significance: where average change over" + \
                                 "all realizations is larger than " + \
                                 "95-percentile of distribution of models."
  to_plot_signif@var_units     = 1

  copy_VarCoords(to_plot, perc95)

  perc95@diag_script = (/diag_script/)
  perc95@var = "95perc"
  perc95@var_long_name = "95-percentile of distribution of models"
  perc95@var_units     = 1

  ; ###########################################
  ; # Output to netCDF               #
  ; ###########################################
  to_plot@ncdf = variable_info@diagnostic + ".nc"
  ncdf_file = work_dir + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  ncdf_file@existing = "append"
  ncdf_outfile = ncdf_write(to_plot_signif, ncdf_file)
  ncdf_outfile = ncdf_write(perc95, ncdf_file)

  ; copy preproc info file to work dir
  system("cp " + diag_script_info@input_files + " " + work_dir)
end
