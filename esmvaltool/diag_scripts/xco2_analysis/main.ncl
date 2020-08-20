; #############################################################################
;  xco2_analysis/global_maps.ncl
; #############################################################################
; Description
;    Plotting global maps for XCO2 data and SCA
;
; Required diag_script_info attributes:
;
; Optional diag_script_info attributes:
;
; Caveats
;
; Modification history
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
  REFERENCES = (/"acknow_author"/)

  ; Variable
  var0 = variable_info[0]@short_name

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  DATASETS := array_append_record(DATASETS, "multi-model mean", 0)
  experiments = metadata_att_as_array(INFO0, "exp") ; Careful! OBS also has exp: esm-hist
  ensembles = metadata_att_as_array(INFO0, "ensemble")
  print(ensembles)
  dim_MOD = dimsizes(DATASETS)
  ALL_FILES = metadata_att_as_array(INFO0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0  + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

end
begin
  ; Minimum amount of available months allowed for calculating SCA
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
    region = "tropics"
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
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ;; Determine start + end year
  start_years = metadata_att_as_array(INFO0, "start_year")
  start_year = min(start_years)
  end_years = metadata_att_as_array(INFO0, "end_year")
  end_year = max(end_years)
  nyear = end_year - start_year + 1
  ntime = nyear*12
  time = new(ntime, float)
  do yy = start_year, end_year
      do mm = 0, 11
          time(12* (yy - start_year) + mm) = yy + (mm+0.5)/12.
      end do
  end do

  if (isatt(INFO0[0], "reference_dataset")) then
    ref_ind = ind(DATASETS.eq.INFO0[0]@reference_dataset)
  else
    ref_ind = -999
  end if

  mod_inds = ind(DATASETS.ne.INFO0[0]@reference_dataset)

  ; Need array that maps obs to index 0, the other models to following indizes
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

  ; Check if ensembles members should be handled as ensemble or seperately
  ensemble_mean = True ;diag_script_info@ensemble_mean
  if ensemble_mean then
    ; find unique models
    unique_models = get_unique_values(DATASETS)
    dim_unique_MOD = dimsizes(unique_models)
  end if


  model_arr = new((/dim_MOD, ntime/), float)
  model_arr!0 = "model"
  model_arr!1 = "time"
  model_arr&model = array_append_record("OBS", DATASETS(mod_inds), 0)
  model_arr&time = time

  gr_arr = model_arr
  sca_arr = model_arr

  growth_series = new((/dim_MOD, nyear/), float)
  growth_series!0 = "model"
  growth_series!1 = "time"
  growth_series&model = array_append_record("OBS", DATASETS(mod_inds), 0)
  growth_series&time = ispan(start_year, end_year, 1)

  amp_series = growth_series
  amp_series@units = "ppmv"
  amp_series@long_name = "Seasonal Cycle Amplitude"
  amp_series@var_yaxis = True

  opt_mask = "_" + str_join(diag_script_info@masking, "_")

    ;; Formatted varname for plots
  if var0.eq."xco2" then
      var0_plotname = "XCO~B~2~N~"
  else
      var0_plotname = var0
  end if

  print("Ensemble refs")
  test = diag_script_info@ensemble_refs
  print(test(:, 0))

end
begin
  ; First read obs to have it ready for masking
  obs_data = read_data(INFO0[ref_ind])
  obs_data = convert_units(obs_data, "ppmv")

  do imod = 0, dim_MOD - 1
    if DATASETS(imod).eq."multi-model mean" then
      var0_mod = mmm_array / tofloat(dimsizes(DATASETS) - subtract_mmm)
      copy_VarMeta(mmm_array, var0_mod)
    else
      if imod.ne.ref_ind then
        tmp = read_data(INFO0[imod])
        var0_mod = tofloat(tmp)
        copy_VarMeta(tmp, var0_mod)
        delete(tmp)
        var0_mod = convert_units(var0_mod, "ppmv")

        if .not. isdefined("mmm_array") then
          mmm_array =  area_hi2lores_Wrap(var0_mod&lon, var0_mod&lat, var0_mod, \
                                          True,1, obs_data&lon, obs_data&lat, False)
        else
          var_reg = area_hi2lores(var0_mod&lon, var0_mod&lat, var0_mod, \
                                  True,1, mmm_array&lon, mmm_array&lat, False)
          mmm_array = mmm_array + (/var_reg/)
          delete(var_reg)
        end if
      end if
    end if

    if imod.eq.ref_ind then
      var0_reg = obs_data
    else
      if any(diag_script_info@masking .eq. "obs") then
        var0_reg = area_hi2lores_Wrap(var0_mod&lon, var0_mod&lat, var0_mod, True, \
                                  1, obs_data&lon, obs_data&lat, False)
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
    gr_arr(mapping_array(imod), :) = (/calc_gr(var0_avg, "monthlyfilled", min_nmonth)/)
    sca_arr(mapping_array(imod), :) = (/(var0_avg - dim_cumsum(gr_arr(mapping_array(imod), :)/12., 1)) - \
                                      avg(var0_avg - dim_cumsum(gr_arr(mapping_array(imod), :)/12., 1))/)
    growth_series(mapping_array(imod), :) = (/calc_gr(var0_avg, "yearly", min_nmonth)/)
    amp_series(mapping_array(imod), :) = (/calc_sca(var0_avg, min_nmonth)/)
    delete([/var0_avg, var0_reg/])
  end do ; End loop over models

  ; -------------------------------------------------------------------------------
  ; -------------------------- PLOTTING -------------------------------------------
  ; -------------------------------------------------------------------------------
  ; --------------------- TIMESERIES ----------------------------------------------
  ; -------------------------------------------------------------------------------

  outfile_ts = plot_dir + "/" + var0 + "_" + experiments(0) \
               + "_" + region + "_time_series_" \
               + start_year + "-" + end_year + opt_mask
  wks = gsn_open_wks(file_type, outfile_ts)

  res = True
  res@tiYAxisString = var0_plotname + " [ppmv]";"  [" + model_arr@units + "]"
  res@tiMainString = "" ;"Time series of " + var0 + " for "\
                      ;+ experiments(0) + " " + region

  ts_plot = xy_line_gen(wks, model_arr, model_arr&time, 0, res, False, INFO0)
  frame(wks)
  delete(res)

  ; -------------------------------------------------------------------------------
  ; ---------------------- SERIES PANELS-------------------------------------------
  ; -------------------------------------------------------------------------------

  outfile_ts_panels = plot_dir + "/" + var0 + "_" + experiments(0) \
               + "_" + region + "_time_series_panels" \
               + start_year + "-" + end_year + opt_mask
  wks = gsn_open_wks(file_type, outfile_ts_panels)

  colors = project_style(INFO0, diag_script_info, "colors")
  colors := array_append_record(colors, "red", 0)
  thicks = project_style(INFO0, diag_script_info, "thicks")
  thicks := array_append_record(thicks, 6., 0)
  dashes = project_style(INFO0, diag_script_info, "dashes")
  dashes := array_append_record(dashes, 0, 0)

  ;; Plot series
  res = True
  res@tiYAxisString = ""
  res@tiYAxisFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02
  res@tiYAxisString = var0_plotname + " [ppmv]";"  [" + model_arr@units + "]"
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

  plots = new(3, graphic)



  res@trYMinF = min(model_arr) - 0.05 * (max(model_arr) - min(model_arr))
  res@trYMaxF = max(model_arr) + 0.05 * (max(model_arr) - min(model_arr)) ;0.22

  plots(0) = gsn_csm_xy(wks, model_arr&time, model_arr, res)

  res@tiYAxisString = "Growth Rate [ppmv yr~S~-1~N~]"
  res@trYMinF = min(gr_arr) - 0.05 * (max(gr_arr) - min(gr_arr))
  res@trYMaxF = max(gr_arr) + 0.05 * (max(gr_arr) - min(gr_arr)) ;0.22

  plots(1) = gsn_csm_xy(wks, gr_arr&time, gr_arr, res)

  res@trYMinF = min(sca_arr) - 0.05 * (max(sca_arr) - min(sca_arr))
  res@trYMaxF = max(sca_arr) + 0.05 * (max(sca_arr) - min(sca_arr)) ;0.22
  res@tiYAxisString = "SCA [ppmv]"
  plots(2) = gsn_csm_xy(wks, sca_arr&time, sca_arr, res)

  ; Add correlation to plot 2!
  do imod = 1, dim_MOD - 1
    corr = escorc(sca_arr(:, 0), sca_arr(:, imod))
    ;print(corr)
  end do


  lgres = True
  lgres@vpWidthF = 0.15
  lgres@vpHeightF = 0.5
  lgres@lgLineColors = colors(::-1)
  lgres@lgLabelFontHeightF = 0.1
  lgres@lgMonoDashIndex     = False
  lgres@lgDashIndexes     = dashes(::-1)
  lgres@lgMonoLineThickness = False
  lgres@lgLineThicknesses = thicks(::-1);*2
  lgres@lgPerimOn = False
  lgres@lgItemType         = "Lines"

  labels = model_arr&model(::-1)
  nlabels = dim_MOD
  legend = gsn_create_legend(wks, nlabels, labels, lgres)

  anres = True
  anres@amJust = "TopLeft"
  ; Parallel -> x-direction
  anres@amParallelPosF = 0.53 ;0.5
  anres@amOrthogonalPosF = -0.5 ;0.6
  annoid1 = gsn_add_annotation(plots(0), legend, anres)

  res1                     = True
  res2                     = True
  res1@gsnMaximize         = True
  res2@gsnAttachPlotsXAxis = True

  newplot = gsn_attach_plots(plots(0), plots(1:), res1, res2)

  draw(plots(0))

  frame(wks)
  delete(res)

  ; -------------------------------------------------------------------------------
  ; ---------------------- BARPLOTS -----------------------------------------------
  ; -------------------------------------------------------------------------------


  ; Set colors so first = OBS is black and last = MMM is red
  bp_colors = new(dimsizes(growth_series&model), string)
  bp_colors = "skyblue"
  bp_colors(0) = "grey40"
  bp_colors(dimsizes(growth_series&model) - 1) = "red"

  ; ---- Growth Rate
  outfile_grbp = plot_dir + "/" + var0 + "_" + experiments(0) \
               + "_" + region + "_" + (start_year) + "-" + (end_year)\
               + "_barplot_grRates" + opt_mask

  wks = gsn_open_wks(file_type, outfile_grbp)

  gr_Rates = new((/dimsizes(growth_series&model), 2/), float)
  gr_Rates(:, 0) = dim_avg_n(growth_series, 1)
  gr_Rates(:, 1) = dim_stddev_n(growth_series, 1)
  asciiwrite(work_dir + "/GR_wstd.txt", gr_Rates)

  ; Sort by size
  sort_by = gr_Rates(:, 0)
  qsort(sort_by)
  perm_idx = new(dimsizes(sort_by), integer)
  do idim = 0, dimsizes(sort_by) -1
    perm_idx(idim) = ind(sort_by(idim).eq.gr_Rates(:,0))
  end do
  gr_Rates_sort = gr_Rates(perm_idx, :)

  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@gsnMaximize = True
  res@gsnXYBarChart = True
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, dimsizes(growth_series&model) - 1, 1)
  res@tmXBLabels = growth_series&model(perm_idx)
  res@tmXBLabelAngleF = 60.
  res@tmXBLabelJust = "CenterRight"
  res@tmXBLabelFontHeightF = 0.014
  res@gsnYRefLine = 0.
  res@trYMinF = 0.
  res@trYMaxF = max(gr_Rates(:,0) + gr_Rates(:,1)) + 0.15
  res@trXMinF = -1.
  res@trXMaxF = dimsizes(growth_series&model)
  res@tiMainString = "" ;"Mean Growth Rate from " + tostring(start_year) + \
                     ;" - " + tostring(end_year) + ", " + region
  res@tiYAxisString = "Average growth rate [ppmv yr~S~-1~N~]"
  res@tiYAxisFontHeightF = 0.02
  res@gsnXYBarChartBarWidth = 0.75
  res@gsnXYBarChartColors = bp_colors(perm_idx)
  ; res@gsnAboveYRefLineColor = "skyblue"

  x = ispan(0, dimsizes(growth_series&model)-1, 1)

  plot_grbp = gsn_csm_xy(wks, x, gr_Rates_sort(:, 0), res)

  ;; Add errorbars

  err_wdth = 0.33 * res@gsnXYBarChartBarWidth

  lres = True
  lres@gsLineThicknessF = 2.

  do imod = 0, dimsizes(growth_series&model) - 1
    plot_grbp@$unique_string("errbar")$ = gsn_add_polyline(wks, plot_grbp, \
          (/imod, imod/), (/gr_Rates_sort(imod, 0) + gr_Rates_sort(imod, 1),  \
          gr_Rates_sort(imod, 0) - gr_Rates_sort(imod, 1)/), lres)
    plot_grbp@$unique_string("errbar")$ = gsn_add_polyline(wks, plot_grbp, \
          (/imod + err_wdth, imod - err_wdth/), (/gr_Rates_sort(imod, 0) \
          + gr_Rates_sort(imod, 1),gr_Rates_sort(imod, 0) \
          + gr_Rates_sort(imod, 1)/), lres)
    plot_grbp@$unique_string("errbar")$ = gsn_add_polyline(wks, plot_grbp, \
          (/imod + err_wdth, imod - err_wdth/), (/gr_Rates_sort(imod, 0) \
          - gr_Rates_sort(imod, 1),gr_Rates_sort(imod, 0) \
          - gr_Rates_sort(imod, 1)/), lres)
  end do

  draw(plot_grbp)
  frame(wks)

  delete(res)
  delete(lres)

  ; ---- Seasonal Cycle Amplitude
  outfile_sabp = plot_dir + "/" + var0 + "_" + experiments(0) \
               + "_" + region + "_" + (start_year) + "-" + (end_year)\
               + "_barplot_SA" + opt_mask

  wks = gsn_open_wks(file_type, outfile_sabp)

  sa_Rates = new((/dimsizes(amp_series&model), 2/), float)
  sa_Rates(:, 0) = dim_avg_n(amp_series, 1)
  sa_Rates(:, 1) = dim_stddev_n(amp_series, 1)
  asciiwrite(work_dir + "/SCA_wstd.txt", sa_Rates)

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
  res@trYMaxF = max(sa_Rates(:,0) + sa_Rates(:,1)) + 0.15
  res@trXMinF = -1.
  res@trXMaxF = dimsizes(amp_series&model)
  res@tiMainString = "Mean Seasonal Cycle Amplitude from " + tostring(start_year) + \
                     " - " + tostring(end_year) + ", " + region
  res@tiYAxisString = "Seasonal Cycle Amplitude [ppmv]"
  res@tiYAxisFontHeightF = 0.02
  res@gsnXYBarChartBarWidth = 0.75
  res@gsnXYBarChartColors = bp_colors
  ; res@gsnAboveYRefLineColor = "skyblue"

  x = ispan(0, dimsizes(amp_series&model)-1, 1)

  plot_sabp = gsn_csm_xy(wks, x, sa_Rates(:, 0), res)

  ;; Add errorbars

  err_wdth = 0.33 * res@gsnXYBarChartBarWidth

  lres = True
  lres@gsLineThicknessF = 2.

  do imod = 0, dimsizes(amp_series&model) - 1
    plot_sabp@$unique_string("errbar")$ = gsn_add_polyline(wks, plot_sabp, \
          (/imod, imod/), (/sa_Rates(imod, 0) + sa_Rates(imod, 1),  \
          sa_Rates(imod, 0) - sa_Rates(imod, 1)/), lres)
    plot_sabp@$unique_string("errbar")$ = gsn_add_polyline(wks, plot_sabp, \
          (/imod + err_wdth, imod - err_wdth/), (/sa_Rates(imod, 0) \
          + sa_Rates(imod, 1),sa_Rates(imod, 0) \
          + sa_Rates(imod, 1)/), lres)
    plot_sabp@$unique_string("errbar")$ = gsn_add_polyline(wks, plot_sabp, \
          (/imod + err_wdth, imod - err_wdth/), (/sa_Rates(imod, 0) \
          - sa_Rates(imod, 1),sa_Rates(imod, 0) \
          - sa_Rates(imod, 1)/), lres)
  end do

  draw(plot_sabp)
  frame(wks)

  delete(plot_sabp)
  delete(res)
  delete(lres)
  delete(err_wdth)

end
