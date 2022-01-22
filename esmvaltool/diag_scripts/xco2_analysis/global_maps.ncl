; #############################################################################
;  xco2_analysis/global_maps.ncl
; #############################################################################
; Description
;    Plotting global maps for XCO2 data and Seasonal Cycle Amplitude
;
; Required diag_script_info attributes:
;   contour_max_level: maximum seasonal cycle value displayed for contour plot
;
; Optional diag_script_info attributes:
;    output_file_type: output file type for plots. Default: png
;
; Caveats
;
; Modification history
;    20201118-gier_bettina: Added provenance, clean up
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
  DOMAIN = (/"global"/)

  ; Variable
  var0 = variable_info[0]@short_name

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  DATASETS := array_append_record(DATASETS, "multi-model mean", 0)
  experiments = metadata_att_as_array(INFO0, "exp")
  dim_MOD = dimsizes(DATASETS)
  ALL_FILES = metadata_att_as_array(INFO0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

end

begin
  ; Maximum amount of missing values per year
  min_nmonth = 7

  ; Plot file type
  if (isatt(diag_script_info, "output_file_type")) then
    file_type = diag_script_info@output_file_type
  elseif (isatt(config_user_info, "output_file_type")) then
    file_type = config_user_info@output_file_type
  else
    file_type = "png"
  end if

  ; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

end

begin
  ; Prepare arrays

  ; Make array to collect timeseries data
  start_years = metadata_att_as_array(INFO0, "start_year")
  start_year = min(start_years)
  end_years = metadata_att_as_array(INFO0, "end_year")
  end_year = max(end_years)
  nyear = (end_year - start_year + 1)
  time = new(nyear*12, float)
  do yy = start_year, end_year
    do mm = 0, 11
      time(12 * (yy - start_year) + mm) = yy + (mm + 0.5)/12.
      end do
  end do

end

begin
  ; First read observations for observational masking
  if (isatt(INFO0[0], "reference_dataset")) then
    ref_ind = ind(DATASETS.eq.INFO0[0]@reference_dataset)
  else
    ref_ind = -999
  end if

  mod_inds = ind(DATASETS.ne.INFO0[0]@reference_dataset)
  mod_names = DATASETS(mod_inds)

  ; Need array that maps obs to index 0, the other models to following indices
  mapping_array = new(dim_MOD, integer)
  if ref_ind.ne.-999 then
    mapping_array(ref_ind) = 0
  end if
  mapping_array(mod_inds) = ispan(1, dim_MOD - 1, 1)

  ; For MMM calculation
  if ref_ind.eq."-999" then
    subtract_mmm = 1
  else
    subtract_mmm = 2
  end if

  if ref_ind.ne."-999" then
    obs_data =  read_data(INFO0[ref_ind])
  end if

  do imod = 0, dim_MOD - 1
    if imod.eq.ref_ind then
      continue
    end if

    if imod.eq.(dim_MOD-1) then
      ; Compute MMM
      mod_var = mmm_array / (dimsizes(DATASETS) - subtract_mmm)
      copy_VarMeta(mmm_array, mod_var)

    else
      tmp = read_data(INFO0[imod])
      mod_var = tofloat(tmp)
      copy_VarMeta(tmp, mod_var)
      delete(tmp)
      var0_units = mod_var@units

      if .not. isdefined("mmm_array") then
        mmm_array =  area_hi2lores_Wrap(mod_var&lon, mod_var&lat, mod_var, \
                                        True, 1, obs_data&lon, obs_data&lat, \
                                        False)
      else
        var_reg = area_hi2lores(mod_var&lon, mod_var&lat, mod_var, \
                                True, 1, mmm_array&lon, mmm_array&lat, False)
        mmm_array = mmm_array + (/var_reg/)
      end if
    end if  ; End MMM check and computation

    mod_nomv = area_hi2lores_Wrap(mod_var&lon, mod_var&lat, mod_var, True, \
                                  1, obs_data&lon, obs_data&lat, False)
    mod_mv = where(ismissing(obs_data), obs_data, mod_nomv)

    ; Calculate growth, seasonal amplitude on a grid level
    ; Generate new arrays to hold values
    nyr = end_year - start_year + 1
    mod_nomv_growth = new((/nyr, dimsizes(mod_nomv&lat), \
                           dimsizes(mod_nomv&lon)/), float)
    mod_nomv_growth!0 = "year"
    mod_nomv_growth&year = ispan(start_year, \
                                 end_year, 1)
    mod_nomv_growth!1 = "lat"
    mod_nomv_growth&lat = mod_nomv&lat
    mod_nomv_growth!2 = "lon"
    mod_nomv_growth&lon = mod_nomv&lon
    mod_nomv_growth@diag_script = DIAG_SCRIPT
    mod_nomv_seasamp = mod_nomv_growth
    mod_mv_growth = mod_nomv_growth
    mod_mv_seasamp = mod_nomv_growth

    mod_nomv_growth@var = var0 + "_growth_rate_no_missing_values"
    mod_nomv_seasamp@var = var0 + "_seasonal_amplitude_no_missing_values"
    mod_mv_growth@var = var0 + "_growth_rate_with_missing_values"
    mod_mv_seasamp@var = var0 + "_seasonal_amplitude_with_missing_values"

    ; Only calculate obs version once
    if .not. isdefined("obs_growth") then
      obs_growth = mod_nomv_growth
      obs_seasamp = mod_nomv_growth
      obs_growth@var = var0 + "_growth_rate_obs"
      obs_seasamp@var = var0 + "_seasonal_amplitude_obs"

      ; Gridded GR and SCA
      obs_monthly = new(dimsizes(obs_data), float)
      copy_VarAtts(obs_data, obs_monthly)

      do itim = 6, dimsizes(obs_data&time)-7
          obs_monthly(itim, :, :) = obs_data(itim + 6, :, :) \
                                    - obs_data(itim - 6, :, :)
      end do

      do itim = 0, nyear - 1
        ; GR
        obs_cum_mv = new(dimsizes(obs_monthly(0, :, :)), float)
        do ilat = 0, dimsizes(obs_monthly(0, :, 0)) - 1
          do ilon = 0, dimsizes(obs_monthly(0, 0, :)) - 1
            obs_cum_mv(ilat, ilon) = num(ismissing(\
                    obs_monthly(itim*12:itim*12+11, ilat, ilon)))
          end do
        end do
        obs_growth(itim, :, :) = \
          where(obs_cum_mv.le.min_nmonth, \
                dim_avg_n(obs_monthly(itim*12:itim*12+11, :, :), 0), \
                obs_growth(itim, :, :))

        ; SCA
        obs_detrended = new(dimsizes(obs_monthly), float)
        obs_cum_mv_sca = new(dimsizes(obs_monthly(0, :, :)), float)

        do ilat = 0, dimsizes(obs_monthly(0, :, 0)) - 1
          do ilon = 0, dimsizes(obs_monthly(0, 0, :)) - 1
            gr_yr = where(ismissing(obs_monthly(itim*12:itim*12+11, ilat, \
                                                ilon)), \
                          obs_growth(itim, ilat, ilon), \
                          obs_monthly(itim*12:itim*12+11, ilat, ilon))
            obs_detrended(itim*12:itim*12+11, ilat, ilon) = \
              obs_data(itim*12:itim*12+11, ilat, ilon) - \
              dim_cumsum_n(gr_yr/12., 1, 0)
            delete(gr_yr)

            obs_cum_mv_sca(ilat, ilon) = num(ismissing(\
                    obs_detrended(itim*12:itim*12+11, ilat, ilon)))
          end do
        end do
        obs_seasamp(itim, :, :) = \
          where(obs_cum_mv_sca.le.min_nmonth, \
                dim_max_n(obs_detrended(itim*12:itim*12+11, :, :), 0) - \
                dim_min_n(obs_detrended(itim*12:itim*12+11, :, :), 0), \
                obs_seasamp(itim, :, :))
      end do

      obs_map_lon = dim_avg_n(obs_seasamp, 2)

      coslat = cos(obs_seasamp&lat/180*3.14159)
      obs_map_latavg = dim_avg_wgt_n(obs_map_lon, coslat, 1, 1)

    end if

    ; Gridded GR and SCA
    mod_mv_monthly = new(dimsizes(mod_mv), float)
    copy_VarAtts(mod_mv, mod_mv_monthly)
    mod_nomv_monthly = new(dimsizes(mod_nomv), float)
    copy_VarAtts(mod_nomv, mod_nomv_monthly)
    do itim = 6, dimsizes(obs_data&time) - 7
      mod_mv_monthly(itim, :, :) = mod_mv(itim+6, :, :) - mod_mv(itim-6, :, :)
      mod_nomv_monthly(itim, :, :) = \
        mod_nomv(itim+6, :, :) - mod_nomv(itim-6, :, :)
    end do

    do itim = 0, nyear - 1
      ; GR
      mod_cum_mv = new(dimsizes(mod_mv(0, :, :)), float)
      mod_cum_nomv = new(dimsizes(mod_nomv(0, :, :)), float)
      do ilat = 0, dimsizes(obs_monthly(0, :, 0)) - 1
        do ilon = 0, dimsizes(obs_monthly(0, 0, :)) - 1
          mod_cum_mv(ilat, ilon) = num(ismissing(\
                     mod_mv_monthly(itim*12:itim*12+11, ilat, ilon)))
          mod_cum_nomv(ilat, ilon) = num(ismissing(\
                     mod_nomv_monthly(itim*12:itim*12+11, ilat, ilon)))
        end do
      end do

      mod_mv_growth(itim, :, :) = \
        where(mod_cum_mv.le.min_nmonth, \
              dim_avg_n(mod_mv_monthly(itim*12:itim*12+11, :, :), 0), \
              mod_mv_growth(itim, :, :))
      mod_nomv_growth(itim, :, :) = \
        where(mod_cum_nomv.le.min_nmonth, \
              dim_avg_n(mod_nomv_monthly(itim*12:itim*12+11, :, :), 0), \
              mod_nomv_growth(itim, :, :))

      ; SCA
      mod_mv_detrended = new(dimsizes(mod_mv), float)
      mod_nomv_detrended = new(dimsizes(mod_nomv), float)
      mod_cum_mv_sca = new(dimsizes(mod_mv_monthly(0, :, :)), float)
      mod_cum_nomv_sca = new(dimsizes(mod_nomv_monthly(0, :, :)), float)

      do ilat = 0, dimsizes(obs_monthly(0, :, 0)) - 1
        do ilon = 0, dimsizes(obs_monthly(0, 0, :)) - 1

          gr_yr = where(ismissing(mod_mv_monthly(itim*12:itim*12+11, \
                                                 ilat, ilon)), \
                        mod_mv_growth(itim, ilat, ilon), \
                        mod_mv_monthly(itim*12:itim*12+11, ilat, ilon))
          mod_mv_detrended(itim*12:itim*12+11, ilat, ilon) = \
            mod_mv(itim*12:itim*12+11, ilat, ilon) \
            - dim_cumsum_n(gr_yr/12., 1, 0)
          delete(gr_yr)

          gr_yr = where(ismissing(mod_nomv_monthly(itim*12:itim*12+11, \
                                                   ilat, ilon)), \
                        mod_nomv_growth(itim, ilat, ilon),\
                        mod_nomv_monthly(itim*12:itim*12+11, ilat, ilon))
          mod_nomv_detrended(itim*12:itim*12+11, ilat, ilon) = \
            mod_nomv(itim*12:itim*12+11, ilat, ilon) \
            - dim_cumsum_n(gr_yr/12., 1, 0)
          delete(gr_yr)

          mod_cum_mv_sca(ilat, ilon) = num(ismissing(\
                     mod_mv_detrended(itim*12:itim*12+11, ilat, ilon)))
          mod_cum_nomv_sca(ilat, ilon) = num(ismissing(\
                     mod_nomv_detrended(itim*12:itim*12+11, ilat, ilon)))
        end do
      end do
      mod_mv_seasamp(itim, :, :) = \
        where(mod_cum_mv_sca.le.min_nmonth, \
              dim_max_n(mod_mv_detrended(itim*12:itim*12+11, :, :), 0) \
              - dim_min_n(mod_mv_detrended(itim*12:itim*12+11, :, :), 0), \
              mod_mv_seasamp(itim, :, :))
      mod_nomv_seasamp(itim, :, :) = \
        where(mod_cum_nomv_sca.le.min_nmonth, \
              dim_max_n(mod_nomv_detrended(itim*12:itim*12+11, :, :), 0) \
              - dim_min_n(mod_nomv_detrended(itim*12:itim*12+11, :, :), 0), \
              mod_nomv_seasamp(itim, :, :))
    end do

    ; #####################################################################
    ; #####################  PLOTTING  ####################################
    ; #####################################################################
    ; Set up the different colormaps to be used
    ; Non-divergent
    non_div_cm = read_colormap_file("MPL_viridis")
    non_div_cm = non_div_cm(::-1, :)
    ; Divergent
    div_cm = read_colormap_file("ncl_default")

    ; Plot seasonal cycle amplitude maps with and without obs sampling
    if DATASETS(imod) .eq. "multi-model mean" then
      fname_dataset = "multi-model_mean"
    else
      fname_dataset = DATASETS(imod)
    end if
    filename = experiments(0) + "_" + fname_dataset + "_xco2_mean_seas_cyc" \
      + "_" + start_year + "-" + end_year
    outfile = plot_dir + "/" + filename
    outfile_netcdf = work_dir + "/" + filename + ".nc"
    wks = gsn_open_wks(file_type, outfile)
    res = True
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnRightString = ""
    res@lbLabelBarOn = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnStringFontHeightF = 0.02
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    res@cnConstFLabelOn = False
    res@cnLineLabelsOn      = False
    res@tmLabelAutoStride = True
    res@lbTitleOn = True
    res@lbTitlePosition = "Bottom"
    res@lbLabelFontHeightF = 0.02
    res@lbTitleFontHeightF = 0.025
    res@lbBottomMarginF = -0.2
    res@lbTopMarginF = 0.4
    res@lbTitleOffsetF = 0.2
    res@lbLabelBarOn = False

    seasamp_plots = new(4, graphic)

    res@gsnLeftString = "OBS to sampled " + DATASETS(imod) + " difference"
    res@lbTitleString = "Seasonal cycle amplitude difference [" \
      + var0_units + "]"
    seasamp_obsdiff = dim_avg_n_Wrap(mod_mv_seasamp, 0)
    seasamp_obsdiff = dim_avg_n_Wrap(obs_seasamp, 0) \
      - dim_avg_n_Wrap(mod_mv_seasamp, 0)
    seasamp_obsdiff@var = "obs_SCA_diff"
    seasamp_diff = dim_avg_n_Wrap(mod_mv_seasamp, 0) \
      - dim_avg_n_Wrap(mod_nomv_seasamp, 0)
    seasamp_diff@var = "SCA_diff_missing_values"
    seasamp_diff!0 = "lat"
    seasamp_diff&lat = mod_mv_growth&lat
    seasamp_diff!1 = "lon"
    seasamp_diff&lon = mod_mv_growth&lon
    maxlev_obs = decimalPlaces(max(abs(seasamp_obsdiff)), 0, False)
    maxlev_diff = decimalPlaces(max(abs(seasamp_diff)), 0, False)
    maxlev = max((/maxlev_obs, maxlev_diff/))
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -maxlev
    res@cnMaxLevelValF = maxlev
    res@cnLevelSpacingF = 1.
    mean_plot = avg(seasamp_obsdiff)
    res@gsnRightString = sprintf("%1.2f", mean_plot) + " " \
      + var0_units
    res@cnFillPalette = div_cm
    plot_seasamp_diff = gsn_csm_contour_map(wks, seasamp_obsdiff, res)

    res@gsnLeftString = "Unsampled to sampled " + DATASETS(imod) \
                        + " difference"
    res@lbTitleString = "Seasonal cycle amplitude difference [" \
      + var0_units + "]"
    mean_plot = avg(seasamp_diff)
    res@gsnRightString = sprintf("%1.2f", mean_plot) + " " \
      + var0_units

    seasamp_plots(3) = gsn_csm_contour_map(wks, seasamp_diff, res)

    ; to have a common label bar,plots should be set to the same interval
    ; b/c the label bar is drawn from the interval of the first plot.
    res@cnLevelSelectionMode =  "ManualLevels"
    res@cnMinLevelValF       = 0.
    res@cnMaxLevelValF       = 12.
    ; toint(max((/max(mod_mv_seasamp), \
    ; max(mod_nomv_seasamp), max(obs_seasamp)/))) + 1
    res@cnLevelSpacingF      = 1.

    res@gsnLeftString = "Sampled " + DATASETS(imod)
    delete(res@cnFillPalette)
    res@cnFillPalette = non_div_cm   ; set color map

    mean_plot = avg(mod_mv_seasamp)
    res@gsnRightString = sprintf("%1.2f", mean_plot) + " " \
      + var0_units
    seasamp_plots(0) = \
      gsn_csm_contour_map(wks, dim_avg_n_Wrap(mod_mv_seasamp, 0), res)

    mean_plot = avg(mod_nomv_seasamp)
    res@gsnRightString = sprintf("%1.2f", mean_plot) + " " \
      + var0_units
    res@gsnLeftString = "Unsampled " + DATASETS(imod)
    seasamp_plots(1) = \
      gsn_csm_contour_map(wks, dim_avg_n_Wrap(mod_nomv_seasamp, 0), res)

    res@gsnLeftString = "OBS"

    mean_plot = avg(obs_seasamp)
    res@gsnRightString = sprintf("%1.2f", mean_plot) + " " \
      + var0_units
    seasamp_plots(2) = \
      gsn_csm_contour_map(wks, dim_avg_n_Wrap(obs_seasamp, 0), res)

    pres = True
    pres@gsnMaximize = True
    pres@gsnFrame = False
    pres@gsnPanelRowSpec = True
    pres@gsnPanelCenter = False
    pres@gsnPanelLabelBar = True
    pres@gsnPanelBottom = 0.4
    pres@lbTitleOn = True
    pres@lbTitlePosition = "Bottom"
    pres@lbTitleString = "Mean seasonal cycle amplitude [" \
      + var0_units + "]"
    pres@lbLabelFontHeightF = 0.012
    pres@lbTitleFontHeightF = 0.015
    pres@lbBottomMarginF = -0.2
    pres@lbTopMarginF = 0.3
    pres@lbTitleOffsetF = 0.2
    pres@gsnPanelYWhiteSpacePercent = 2.

    gsn_panel(wks, seasamp_plots(:2), (/2, 1/), pres)

    pres@gsnPanelTop = 0.4
    pres@gsnPanelBottom = 0.01

    gsn_panel(wks, (/plot_seasamp_diff, seasamp_plots(3)/), (/2/), pres)
    frame(wks)
    delete(pres)

  ; Write data to netcdf
  outfile_netcdf@existing = "overwrite"
  ncdf_outfile = ncdf_write(obs_seasamp, outfile_netcdf)
  outfile_netcdf@existing = "append"
  ncdf_outfile = ncdf_write(mod_nomv_seasamp, outfile_netcdf)
  ncdf_outfile = ncdf_write(mod_mv_seasamp, outfile_netcdf)
  ncdf_outfile = ncdf_write(seasamp_diff, outfile_netcdf)
  ncdf_outfile = ncdf_write(seasamp_obsdiff, outfile_netcdf)

  ; Provenance
  log_provenance(ncdf_outfile, \
                 outfile, \
                 "Maps of mean Seasonal Cycle Amplitude of " + var0 + " " \
                 + DATASETS(imod) + " for " + start_year + "-" \
                 + end_year + ". Top: SCA with observational sampling (left)" \
                 + ", without sampling (right). Middle: SCA of observations." \
                 + "Bottom: Difference between observations and models " \
                 + "(left) and sampled and unsampled model (right)." \
                 + " Similar to Gier et al. 2020, Figure 6.", \
                 (/"mean", "diff"/), \
                 DOMAIN, \
                 "geo", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

    ; ------------------------- Collect for panel plots
    seasamp_avg = dim_avg_n_Wrap(mod_nomv_seasamp, 0)
    if .not. isdefined("sca_p_array") then
        sca_p_dims = array_append_record(dimsizes(mod_inds), \
                                         dimsizes(seasamp_avg), 0)
        dim = 0
        sca_p_array = new(sca_p_dims, float)
        sca_p_array!0 = "models"
        sca_p_array&models = DATASETS(mod_inds)
        sca_p_array!1 = "lat"
        sca_p_array&lat = seasamp_avg&lat
        sca_p_array!2 = "lon"
        sca_p_array&lon = seasamp_avg&lon
        sca_p_array@var = "mean_seasonal_cycle_amplitude"
        sca_p_array@diag_script = DIAG_SCRIPT
        delete(sca_p_dims)
    end if
    sca_p_array(dim, :, :) = (/seasamp_avg/)
    dim = dim + 1
    delete(seasamp_avg)

    delete(mod_var)
    delete(mod_mv)
    delete(mod_nomv)
    delete(res@cnFillPalette)
  end do  ; Loop over models

  ; ---------------- Make panel plot

  do imod = 0, dimsizes(mod_inds) - 1

    if .not. isdefined("sca_panel_plots") then
      filename_scap = var0 + "_" + experiments(0) \
                     + "_SCApanels_" \
                     + start_year + "-" + end_year
      outfile_scap = plot_dir + "/" + filename_scap
      outfile_netcdf_scap = work_dir + "/" + filename_scap + ".nc"
      wks_scap = gsn_open_wks(file_type, outfile_scap)
      sca_panel_plots = new(dimsizes(mod_inds), graphic)
      res_sca_panel_plots = True
      res_sca_panel_plots@cnFillOn = True
      res_sca_panel_plots@cnLinesOn = False
      res_sca_panel_plots@lbLabelBarOn = False
      res_sca_panel_plots@gsnDraw = False
      res_sca_panel_plots@gsnFrame = False
      res_sca_panel_plots@cnLevelSelectionMode = "ManualLevels"
      res_sca_panel_plots@cnMinLevelValF = 0.
      res_sca_panel_plots@cnMaxLevelValF = diag_script_info@contour_max_level
      res_sca_panel_plots@cnLevelSpacingF = 1.
      res_sca_panel_plots@gsnStringFontHeightF = 0.03
      res_sca_panel_plots@tmXBLabelFontHeightF = 0.025
      res_sca_panel_plots@tmYLLabelFontHeightF = 0.025
      res_sca_panel_plots@cnConstFLabelOn = False
      res_sca_panel_plots@cnLineLabelsOn      = False
      res_sca_panel_plots@tmLabelAutoStride = True
      res_sca_panel_plots@gsnZonalMean = True
      res_sca_panel_plots@gsnZonalMeanXMinF = 0
      res_sca_panel_plots@gsnZonalMeanXMaxF = \
        diag_script_info@contour_max_level
      res_sca_panel_plots@amOrthogonalPosF = .2
      res_sca_panel_plots@cnFillPalette = non_div_cm   ; set color map
      res_sca_panel_plots@gsnRightStringOrthogonalPosF = -0.005
      res_sca_panel_plots@gsnLeftStringOrthogonalPosF = -0.005
    end if
    res_sca_panel_plots@gsnLeftString = mod_names(imod)
    mean_sca_mod = avg(sca_p_array(imod, :, :))
    res_sca_panel_plots@gsnRightString = \
      sprintf("%1.2f", mean_sca_mod) + " " + var0_units
    sca_panel_plots(imod) = \
      gsn_csm_contour_map(wks_scap, sca_p_array(imod, :, :), \
                          res_sca_panel_plots)
  end do

  pres = True
  pres@gsnMaximize = True
  pres@gsnFrame = False
  pres@gsnPanelCenter = False
  pres@gsnPanelLabelBar = True
  pres@gsnPanelBottom = 0.1
  pres@lbTitleOn = True
  pres@lbTitlePosition = "Bottom"
  pres@lbTitleString = \
    "Mean seasonal cycle amplitude [" + var0_units + "]"
  pres@lbLabelFontHeightF = 0.012
  pres@lbTitleFontHeightF = 0.015
  pres@lbBottomMarginF = -1
  pres@lbTopMarginF = 0.8
  pres@lbTitleOffsetF = 0.2

  ncol = toint(floor(sqrt(dimsizes(mod_inds))))
  nrow = toint(ceil((dimsizes(mod_inds)/tofloat(ncol))))
  gsn_panel(wks_scap, sca_panel_plots, (/nrow, ncol/), pres)
  frame(wks_scap)
  delete(res_sca_panel_plots)
  delete(pres)

  ncdf_outfile = ncdf_write(sca_p_array, outfile_netcdf_scap)

  ; Provenance
  log_provenance(ncdf_outfile, \
                 outfile_scap, \
                 "Maps of mean Seasonal Cycle Amplitude of " \
                 + var0 + " for " + start_year + "-" \
                 + end_year + ". Similar to Gier et al 2020, Fig 5.", \
                 (/"mean", "detrend"/), \
                 DOMAIN, \
                 "geo", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

  leave_msg(DIAG_SCRIPT, "")
end
