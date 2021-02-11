; #############################################################################
;  xco2_analysis/main.ncl
; #############################################################################
; Description
;    Plotting timeseries and histograms for growth rate and seasonal cycle.
;
; Required diag_script_info attributes:
;    styleset: styleset to use for plotting colors, linestyles...
;    region: latitude range for averaging
;    masking: different masking options are available to use on dataset:
;              "none" - no masking
;              "obs" - observational masking
;    ensemble_mean: if true calculates multi-model mean only accounting for
;                   the ensemble member named in "ensemble_refs"
;
; Optional diag_script_info attributes:
;    output_file_type: output file type for plots. Default: png
;    ensemble_refs: list of model-ensemble pairs to denote which ensemble
;                   member to use for calculating multi-model mean. required if
;                   ensemble_mean = true
;    var_plotname: String formatting how variable should be named in plots
;                  defaults to short_name if not assigned
;
; Caveats
;
; Modification history
;    20201116-gier_bettina: Added provenance, clean up
;    20200226-gier_bettina: Adapted to version 2
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"

load "$diag_scripts/xco2_analysis/stat.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  AUTHORS = (/"gier_bettina"/)
  REFERENCES = (/"gier20bg"/)

  ; Variable
  var0 = variable_info[0]@short_name

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  DATASETS := array_append_record(DATASETS, "multi-model mean", 0)
  experiments = metadata_att_as_array(INFO0, "exp")
  ensembles = metadata_att_as_array(INFO0, "ensemble")
  dim_MOD = dimsizes(DATASETS)
  ALL_FILES = metadata_att_as_array(INFO0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

end
begin
  ; Maximum amount of missing values per year
  min_nmonth = 7

  ; Prepare region
  lat_min = diag_script_info@region(0)
  lat_max = diag_script_info@region(1)

  if lat_min.eq.(-90) .and. lat_max.eq.90 then
    region = "global"
  else if lat_min.eq.(30) .and. lat_max.eq.60 then
    region = "nhmidlat"
  else if lat_min.eq.(-60) .and. lat_max.eq.(-30) then
    region = "shmidlat"
  else if lat_min.eq.(-30) .and. lat_max.eq.(30) then
    region = "trop"
  else if lat_min.eq.(0) .and. lat_max.eq.(90) then
    region = "nh"
  else if lat_min.eq.(-90) .and. lat_max.eq.0 then
    region = "sh"
  else
    region = "lat_" + tostring(lat_min) + "_" + tostring(lat_max)
  end if
  end if
  end if
  end if
  end if
  end if
  DOMAIN = (/region/)

  ; Plot file type
  if (isatt(diag_script_info, "output_file_type")) then
    file_type = diag_script_info@output_file_type
  elseif (isatt(config_user_info, "output_file_type")) then
    file_type = config_user_info@output_file_type
  else
    file_type = "png"
  end if

  ; Output directories
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ; Determine start + end year
  start_years = metadata_att_as_array(INFO0, "start_year")
  start_year = min(start_years)
  end_years = metadata_att_as_array(INFO0, "end_year")
  end_year = max(end_years)
  nyear = end_year - start_year + 1
  ntime = nyear*12
  time = new(ntime, float)
  do yy = start_year, end_year
      do mm = 0, 11
          time(12 * (yy - start_year) + mm) = yy + (mm + 0.5)/12.
      end do
  end do

  if (isatt(INFO0[0], "reference_dataset")) then
    ref_ind = ind(DATASETS.eq.INFO0[0]@reference_dataset)
  else
    ref_ind = -999
  end if

  mod_inds = ind(DATASETS.ne.INFO0[0]@reference_dataset)

  ; Need array that maps obs to index 0, the other models to following indices
  mapping_array = new(dimsizes(DATASETS), integer)
  if ref_ind.ne.-999 then
    mapping_array(ref_ind) = 0
  end if
  mapping_array(mod_inds) = ispan(1, dimsizes(DATASETS)-1, 1)

  ; For MMM calculation
  if ref_ind.eq."-999" then
    subtract_mmm = 1
  else
    subtract_mmm = 2
  end if

  ; Check if ensembles members should be handled as ensemble or separately
  calc_ensemble_mean = diag_script_info@ensemble_mean
  if calc_ensemble_mean then
    ensemble_refs = diag_script_info@ensemble_refs
  end if
  ; find unique models
  unique_models = get_unique_values(DATASETS)
  dim_unique_MOD = dimsizes(unique_models)

  model_arr = new((/dim_MOD, ntime/), float)
  model_arr!0 = "model"
  model_arr!1 = "time"
  model_arr&model = array_append_record("OBS", DATASETS(mod_inds), 0)
  model_arr&time = time
  model_arr@var = var0
  model_arr@diag_script = (/DIAG_SCRIPT/)

  gr_arr = model_arr
  gr_arr@var = "gr"
  sca_arr = model_arr
  sca_arr@var = "sca"

  growth_series = new((/dim_MOD, nyear/), float)
  growth_series!0 = "model"
  growth_series!1 = "time"
  growth_series&model = array_append_record("OBS", DATASETS(mod_inds), 0)
  growth_series&time = ispan(start_year, end_year, 1)

  amp_series = growth_series
  amp_series@long_name = "Seasonal Cycle Amplitude"
  amp_series@var_yaxis = True

  opt_mask = "_" + str_join(diag_script_info@masking, "_")

  ; Formatted varname for plots
  if (isatt(diag_script_info, "var_plotname")) then
      var0_plotname = diag_script_info@var_plotname
  else
      var0_plotname = var0
  end if

end
begin
  ; First read obs to have it ready for masking
  obs_data = read_data(INFO0[ref_ind])

  model_arr@units = obs_data@units
  gr_arr@units = obs_data@units + " yr-1"
  growth_series@units = gr_arr@units
  amp_series@units = obs_data@units

  do imod = 0, dim_MOD - 1
    if DATASETS(imod).eq."multi-model mean" then
      if calc_ensemble_mean then
        var0_mod = mmm_array / tofloat(dim_unique_MOD - subtract_mmm)
      else
        var0_mod = mmm_array / tofloat(dimsizes(DATASETS) - subtract_mmm)
      end if
      copy_VarMeta(mmm_array, var0_mod)
    else
      ; Check if it's an ensemble member
      if calc_ensemble_mean then
        if imod.eq.ref_ind then
          mmm_member = False
        elseif .not. calc_ensemble_mean then
          mmm_member = True
        elseif all(ensemble_refs(:, 0).ne.DATASETS(imod))
          mmm_member = True
        else
          ensemble_ref_ind = ind(ensemble_refs(:, 0).eq.DATASETS(imod))
          if ensemble_refs(ensemble_ref_ind, 1).eq. ensembles(imod) then
            mmm_member = True
          else
            mmm_member = False
          end if
        end if
      else
        mmm_member = True
      end if
      if (imod.ne.ref_ind)then
        tmp = read_data(INFO0[imod])
        var0_mod = tofloat(tmp)
        copy_VarMeta(tmp, var0_mod)
        delete(tmp)

        if mmm_member
          if .not. isdefined("mmm_array") then
            mmm_array =  area_hi2lores_Wrap(var0_mod&lon, var0_mod&lat, \
                                            var0_mod, True, 1, obs_data&lon, \
                                            obs_data&lat, False)
          else
            var_reg = area_hi2lores(var0_mod&lon, var0_mod&lat, var0_mod, \
                                    True, 1, mmm_array&lon, mmm_array&lat, \
                                    False)
            mmm_array = mmm_array + (/var_reg/)
            delete(var_reg)
          end if
        end if
      end if
    end if

    if imod.eq.ref_ind then
      var0_reg = obs_data
    else
      if any(diag_script_info@masking .eq. "obs") then
        var0_reg = area_hi2lores_Wrap(var0_mod&lon, var0_mod&lat, var0_mod, \
                                      True, 1, obs_data&lon, obs_data&lat, \
                                      False)
        var0_reg = where(ismissing(obs_data), obs_data, var0_reg)
        delete(var0_mod)
      else
        var0_reg = var0_mod
        delete(var0_mod)
      end if
    end if

    var0_avg = area_operations(var0_reg, lat_min, lat_max, 0, 360, \
                               "average", True)
    model_arr(mapping_array(imod), :) = (/var0_avg/)
    gr_arr(mapping_array(imod), :) = (/calc_gr(var0_avg, "monthly", \
                                               min_nmonth)/)
    mnth_filled = calc_gr(var0_avg, "monthlyfilled", min_nmonth)
    sca_arr(mapping_array(imod), :) = (/(var0_avg - \
                                         dim_cumsum(mnth_filled/12., 1)) - \
                                       avg(var0_avg - \
                                           dim_cumsum(mnth_filled/12., 1))/)
    growth_series(mapping_array(imod), :) = (/calc_gr(var0_avg, "yearly", \
                                                      min_nmonth)/)
    amp_series(mapping_array(imod), :) = (/calc_sca(var0_avg, min_nmonth)/)
    delete([/var0_avg, var0_reg, mnth_filled/])
  end do  ; End loop over models

  ; ---------------------------------------------------------------------------
  ; -------------------------- PLOTTING ---------------------------------------
  ; ---------------------------------------------------------------------------

  ; ---------------------------------------------------------------------------
  ; ---------------------- TIMESERIES PANELS-----------------------------------
  ; ---------------------------------------------------------------------------

  filename_ts_panels = var0 + "_" + experiments(0) + "_" + region \
    + "_time_series_panels_" + start_year + "-" + end_year + opt_mask

  outfile_netcdf_ts_panels = work_dir + "/" + filename_ts_panels + ".nc"
  model_arr@ncdf = outfile_netcdf_ts_panels
  gr_arr@ncdf = outfile_netcdf_ts_panels
  sca_arr@ncdf = outfile_netcdf_ts_panels

  outfile_netcdf_ts_panels@existing = "overwrite"
  ncdf_outfile = ncdf_write(model_arr, outfile_netcdf_ts_panels)
  outfile_netcdf_ts_panels@existing = "append"
  ncdf_outfile = ncdf_write(gr_arr, outfile_netcdf_ts_panels)
  ncdf_outfile = ncdf_write(sca_arr, outfile_netcdf_ts_panels)

  outfile_ts_panels = plot_dir + "/" + filename_ts_panels
  wks = gsn_open_wks(file_type, outfile_ts_panels)

  colors = project_style(INFO0, diag_script_info, "colors")
  colors := array_append_record(colors, "red", 0)
  ; Thicker lines for OBS (first) and MMM (last)
  thicks = new(dim_MOD, float)
  thicks(0) = 6.
  thicks(1:dim_MOD-2) = 3.
  thicks(dim_MOD-1) = 6.
  dashes = project_style(INFO0, diag_script_info, "dashes")
  dashes := array_append_record(dashes, 0, 0)

  labels_ini = model_arr&model
  if dim_unique_MOD .ne. dim_MOD then
    inds = ispan(0, dim_MOD - 1, 1)
    do imod = 0, dim_MOD - 1
      other_sets = DATASETS(ind(inds.ne.imod))
      if any(other_sets.eq.DATASETS(imod)) then
        labels_ini(imod) = labels_ini(imod) + " " + ensembles(imod)
        get_r = str_split(ensembles(imod), "ri")
        dashes(imod) = toint(get_r(0)) - 1
      end if
      delete(other_sets)
    end do
  end if

  ; Plot series
  res = True
  res@tiYAxisString = ""
  res@tiYAxisFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02
  res@tiYAxisString = var0_plotname + " [" + model_arr@units + "]"
  res@tiXAxisFontHeightF = 0.02
  res@tmXBLabelFontHeightF = 0.02
  res@pmTickMarkDisplayMode = "Conditional"
  res@tiMainFontHeightF = 0.02
  res@tiMainString = ""
  res@gsnMaximize = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@trXMinF = start_year
  res@trXMaxF = end_year + 1
  res@vpHeightF = 0.3
  res@vpWidthF = 0.8
  res@tmLabelAutoStride = True
  res@xyLineColors      = colors  ; change line color
  res@xyLineThicknesses = thicks  ; make 2nd lines thicker
  res@xyMonoDashPattern = False
  res@xyDashPatterns    = dashes
  res@tmXBFormat = "f"

  plots = new(3, graphic)

  res@trYMinF = min(model_arr) - 0.05 * (max(model_arr) - min(model_arr))
  res@trYMaxF = max(model_arr) + 0.05 * (max(model_arr) - min(model_arr))

  plots(0) = gsn_csm_xy(wks, model_arr&time, model_arr, res)
  ; Redraw obs to make it more visible on top
  respoly = True
  respoly@gsLineThicknessF = thicks(0)
  plots@poly0 = gsn_add_polyline(wks, plots(0), model_arr&time, \
                                 model_arr(0, :), respoly)

  res@tiYAxisString = "Growth Rate [" + model_arr@units + " yr~S~-1~N~]"
  res@trYMinF = min(gr_arr) - 0.05 * (max(gr_arr) - min(gr_arr))
  res@trYMaxF = max(gr_arr) + 0.05 * (max(gr_arr) - min(gr_arr))

  plots(1) = gsn_csm_xy(wks, gr_arr&time, gr_arr, res)
  plots@poly1 = gsn_add_polyline(wks, plots(1), gr_arr&time, \
                                 gr_arr(0, :), respoly)

  res@trYMinF = min(sca_arr) - 0.05 * (max(sca_arr) - min(sca_arr))
  res@trYMaxF = max(sca_arr) + 0.05 * (max(sca_arr) - min(sca_arr))
  res@tiYAxisString = "Seasonal Cycle [" + model_arr@units + "]"
  plots(2) = gsn_csm_xy(wks, sca_arr&time, sca_arr, res)
  plots@poly2 = gsn_add_polyline(wks, plots(2), sca_arr&time, \
                                 sca_arr(0, :), respoly)

  ; Calc annual mean GR correlations
  gr_corrs = new(dim_MOD-1, float)
  gr_corrs!0 = "models"
  gr_corrs&models = labels_ini(1:)

  ; Calc SCA correlations
  sca_corrs = gr_corrs

  do imod = 1, dim_MOD - 1
    gr_corrs(imod-1) = escorc(growth_series(0, :), growth_series(imod, :))
    corr = escorc(sca_arr(0, :), sca_arr(imod, :))
    sca_corrs(imod-1) = corr
  end do
  print("Annual Mean Growth Rate correlation to observations")
  print(gr_corrs&models + ": " + gr_corrs)

  print("SCA Correlation to observations")
  print(sca_corrs&models + ": " + sca_corrs)

  lgres = True
  lgres@vpWidthF = 0.25
  lgres@vpHeightF = 0.5
  lgres@lgLineColors = colors(::-1)
  lgres@lgLabelFontHeightF = 0.1
  lgres@lgMonoDashIndex = False
  lgres@lgDashIndexes = dashes(::-1)
  lgres@lgMonoLineThickness = False
  lgres@lgLineThicknesses = thicks(::-1)
  lgres@lgPerimOn = False
  lgres@lgItemType = "Lines"

  lgres@lgAutoManage = False
  lgres@lgLabelFontHeightF = 0.06

  labels = labels_ini(::-1)
  nlabels = dim_MOD
  legend = gsn_create_legend(wks, nlabels, labels, lgres)

  anres = True
  anres@amJust = "TopLeft"
  ; Parallel -> x-direction
  anres@amParallelPosF = 0.53
  anres@amOrthogonalPosF = -0.5
  annoid1 = gsn_add_annotation(plots(0), legend, anres)

  res1                     = True
  res2                     = True
  res1@gsnMaximize         = True
  res2@gsnAttachPlotsXAxis = True

  newplot = gsn_attach_plots(plots(0), plots(1:), res1, res2)

  draw(plots(0))

  frame(wks)
  delete(res)

  ; Provenance
  log_provenance(ncdf_outfile, \
                 outfile_ts_panels, \
                 "Timeseries, growth rate and seasonal cycle of " \
                 + var0 + " for " + region + ", " + start_year + "-" \
                 + end_year + ". Using masking: " + opt_mask \
                 + "Similar to Gier et al. 2020 Fig 3", \
                 (/"mean", "diff", "detrend"/), \
                 DOMAIN, \
                 "times", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

  ; ---------------------------------------------------------------------------
  ; ---------------------- BARPLOTS -------------------------------------------
  ; ---------------------------------------------------------------------------

  ; Set colors so first = OBS is black and last = MMM is red
  bp_colors = new(dimsizes(growth_series&model), string)
  bp_colors = "skyblue"
  bp_colors(0) = "grey40"
  bp_colors(dimsizes(growth_series&model) - 1) = "red"

  ; ---- Growth Rate
  filename_grbp = var0 + "_" + experiments(0) \
    + "_" + region + "_" + (start_year) + "-" + (end_year) \
    + "_barplot_grRates" + opt_mask
  outfile_grbp = plot_dir + "/" + filename_grbp
  outfile_netcdf_grbp = work_dir + "/" + filename_grbp + ".nc"

  wks = gsn_open_wks(file_type, outfile_grbp)

  gr_Rates = new((/dimsizes(growth_series&model), 2/), float)
  gr_Rates(:, 0) = dim_avg_n(growth_series, 1)
  gr_Rates(:, 1) = dim_stddev_n(growth_series, 1)
  gr_Rates@units = model_arr@units + " yr~S~-1~N~"
  gr_Rates!0 = "models"
  gr_Rates&models = growth_series&model
  gr_Rates!1 = "statistics"
  gr_Rates&statistics = (/"mean", "stddev"/)
  gr_Rates@var = "gr"
  gr_Rates@ncdf = outfile_netcdf_grbp
  gr_Rates@diag_script = (/DIAG_SCRIPT/)
  printVarSummary(gr_Rates)
  ncdf_outfile = ncdf_write(gr_Rates, outfile_netcdf_grbp)

  ; Sort by size
  sort_by = gr_Rates(:, 0)
  qsort(sort_by)
  perm_idx = new(dimsizes(sort_by), integer)
  do idim = 0, dimsizes(sort_by) - 1
    perm_idx(idim) = ind(sort_by(idim).eq.gr_Rates(:, 0))
  end do
  gr_Rates_sort = gr_Rates(perm_idx, :)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize = True
  res@gsnXYBarChart = True
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, dimsizes(growth_series&model) - 1, 1)
  res@tmXBLabels = labels_ini(perm_idx)
  res@tmXBLabelAngleF = 60.
  res@tmXBLabelJust = "CenterRight"
  res@tmXBLabelFontHeightF = 0.014
  res@gsnYRefLine = 0.
  res@trYMinF = 0.
  res@trYMaxF = max(gr_Rates(:, 0) + gr_Rates(:, 1)) + 0.15
  res@trXMinF = -1.
  res@trXMaxF = dimsizes(growth_series&model)
  res@tiYAxisString = "Average growth rate [" + gr_Rates@units + "]"
  res@tiYAxisFontHeightF = 0.02
  res@gsnXYBarChartBarWidth = 0.75
  res@gsnXYBarChartColors = bp_colors(perm_idx)

  x = ispan(0, dimsizes(growth_series&model)-1, 1)

  plot_grbp = gsn_csm_xy(wks, x, gr_Rates_sort(:, 0), res)

  ; Add errorbars

  err_wdth = 0.33 * res@gsnXYBarChartBarWidth

  lres = True
  lres@gsLineThicknessF = 2.

  do imod = 0, dimsizes(growth_series&model) - 1
    plot_grbp@$unique_string("errbar")$ = \
      gsn_add_polyline(wks, plot_grbp, (/imod, imod/), \
                       (/gr_Rates_sort(imod, 0) + gr_Rates_sort(imod, 1),  \
                         gr_Rates_sort(imod, 0) - gr_Rates_sort(imod, 1)/), \
                       lres)
    plot_grbp@$unique_string("errbar")$ = \
      gsn_add_polyline(wks, plot_grbp, (/imod + err_wdth, imod - err_wdth/), \
                       (/gr_Rates_sort(imod, 0) + gr_Rates_sort(imod, 1), \
                        gr_Rates_sort(imod, 0) + gr_Rates_sort(imod, 1)/), \
                       lres)
    plot_grbp@$unique_string("errbar")$ = \
      gsn_add_polyline(wks, plot_grbp, (/imod + err_wdth, imod - err_wdth/), \
                       (/gr_Rates_sort(imod, 0) - gr_Rates_sort(imod, 1), \
                        gr_Rates_sort(imod, 0) - gr_Rates_sort(imod, 1)/), \
                       lres)
  end do

  draw(plot_grbp)
  frame(wks)

  delete(res)
  delete(lres)

  ; Provenance
  log_provenance(ncdf_outfile, \
                 outfile_grbp, \
                 "Growth Rate histogram of " \
                 + var0 + " for " + region + ", " + start_year + "-" \
                 + end_year + ". Using masking: " + opt_mask \
                 + "Similar to Gier et al. 2020, Fig 4.", \
                 (/"mean", "stddev"/), \
                 DOMAIN, \
                 "histogram", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

  ; ---------------------------------------------------------------------------
  ; ---- Seasonal Cycle Amplitude ---------------------------------------------
  filename_sabp = var0 + "_" + experiments(0) + "_" + region \
    + "_" + (start_year) + "-" + (end_year) + "_barplot_SA" + opt_mask
  outfile_sabp = plot_dir + "/" + filename_sabp
  outfile_netcdf_sabp = work_dir + "/" + filename_sabp + ".nc"

  wks = gsn_open_wks(file_type, outfile_sabp)

  sa_Rates = new((/dimsizes(amp_series&model), 2/), float)
  sa_Rates(:, 0) = dim_avg_n(amp_series, 1)
  sa_Rates(:, 1) = dim_stddev_n(amp_series, 1)
  sa_Rates!0 = "models"
  sa_Rates&models = amp_series&model
  sa_Rates!1 = "statistics"
  sa_Rates&statistics = (/"mean", "stddev"/)
  sa_Rates@units = model_arr@units
  sa_Rates@var = "sca"
  sa_Rates@ncdf = outfile_netcdf_sabp
  sa_Rates@diag_script = (/DIAG_SCRIPT/)
  ncdf_outfile = ncdf_write(sa_Rates, outfile_netcdf_grbp)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize = True
  res@gsnXYBarChart = True
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, dimsizes(amp_series&model) - 1, 1)
  res@tmXBLabels = amp_series&model
  res@tmXBLabelAngleF = 60.
  res@tmXBLabelJust = "CenterRight"
  res@tmXBLabelFontHeightF = 0.014
  res@gsnYRefLine = 0.
  res@trYMinF = 0.
  res@trYMaxF = max(sa_Rates(:, 0) + sa_Rates(:, 1)) + 0.15
  res@trXMinF = -1.
  res@trXMaxF = dimsizes(amp_series&model)
  res@tiMainString = "Mean Seasonal Cycle Amplitude from " + \
                     tostring(start_year) + " - " + tostring(end_year) + \
                     ", " + region
  res@tiYAxisString = "Seasonal Cycle Amplitude [" + sa_Rates@units + "]"
  res@tiYAxisFontHeightF = 0.02
  res@gsnXYBarChartBarWidth = 0.75
  res@gsnXYBarChartColors = bp_colors

  x = ispan(0, dimsizes(amp_series&model)-1, 1)

  plot_sabp = gsn_csm_xy(wks, x, sa_Rates(:, 0), res)

  ; Add errorbars

  err_wdth = 0.33 * res@gsnXYBarChartBarWidth

  lres = True
  lres@gsLineThicknessF = 2.

  do imod = 0, dimsizes(amp_series&model) - 1
    plot_sabp@$unique_string("errbar")$ = \
      gsn_add_polyline(wks, plot_sabp, (/imod, imod/), \
                       (/sa_Rates(imod, 0) + sa_Rates(imod, 1),  \
                        sa_Rates(imod, 0) - sa_Rates(imod, 1)/), lres)
    plot_sabp@$unique_string("errbar")$ = \
      gsn_add_polyline(wks, plot_sabp, (/imod + err_wdth, imod - err_wdth/), \
                       (/sa_Rates(imod, 0) + sa_Rates(imod, 1), \
                        sa_Rates(imod, 0) + sa_Rates(imod, 1)/), lres)
    plot_sabp@$unique_string("errbar")$ = \
      gsn_add_polyline(wks, plot_sabp, (/imod + err_wdth, imod - err_wdth/), \
                       (/sa_Rates(imod, 0) - sa_Rates(imod, 1), \
                        sa_Rates(imod, 0) - sa_Rates(imod, 1)/), lres)
  end do

  draw(plot_sabp)
  frame(wks)

  delete(plot_sabp)
  delete(res)
  delete(lres)
  delete(err_wdth)

  ; Provenance
  log_provenance(ncdf_outfile, \
                 outfile_sabp, \
                 "Seasonal Cycle Amplitude histogram of " \
                 + var0 + " for " + region + ", " + start_year + "-" \
                 + end_year + ". Using masking: " + opt_mask, \
                 (/"mean", "stddev"/), \
                 DOMAIN, \
                 "histogram", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

end
