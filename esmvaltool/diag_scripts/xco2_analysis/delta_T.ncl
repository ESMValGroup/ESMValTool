; #############################################################################
;  xco2_analysis/delta_T.ncl
; #############################################################################
; Description
;    Computes IAV of growth rate and plots against growing season temperature.
;
; Required diag_script_info attributes:
;    region: region to average over
;    masking: the kind of masking to apply prior to region average
;    var_order: First main variable, then temperature variable to compare
;
; Optional diag_script_info attributes:
;    styleset: styleset for color cording panels
;    output_file_type: output file type for plots
;    var_plotname: String formatting how variable should be named in plots
;                  defaults to short_name if not assigned
;
; Caveats
;
; Modification history
;    20201116-gier_bettina: Added provenance, clean up
;    20200227-gier_bettina: Adapted to version 2
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
load "$diag_scripts/xco2_analysis/carbon_plots.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  AUTHORS = (/"gier_bettina"/)
  REFERENCES = (/"gier20bg"/)

  ; Variable
  var0 = diag_script_info@var_order(0)
  var1 = diag_script_info@var_order(1)
  var2 = diag_script_info@var_order(2)

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  INFO1 = select_metadata_by_name(input_file_info, var1)
  INFO2 = select_metadata_by_name(input_file_info, var2)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  ; Rename MultiModelMean to multi-model mean for publication
  DATASETS = where(DATASETS.eq."MultiModelMean", "multi-model mean", DATASETS)
  DATASETS2 = metadata_att_as_array(INFO1, "dataset")
  DATASETS2 = where(DATASETS2.eq."MultiModelMean", "multi-model mean", \
                    DATASETS2)
  ALL_FILES = metadata_att_as_array(INFO0, "filename")
  ALL_FILES := array_append_record(ALL_FILES, \
                                   metadata_att_as_array(INFO1, "filename"), 0)
  experiments = metadata_att_as_array(INFO0, "exp")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + " and " + var1 + ")")
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

  ; Output plot directory
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
  mapping_array1 = new(dimsizes(DATASETS), integer)
  if ref_ind.ne.-999 then
    mapping_array1(ref_ind) = 0
  end if
  mapping_array1(mod_inds) = ispan(1, dimsizes(DATASETS)-1, 1)

  ; Need array that maps obs to index 0, the other models to following indizes
  mapping_array2 = new(dimsizes(DATASETS2) + 1, integer)
  mapping_array2(dimsizes(DATASETS2)) = 0
  mapping_array2(:dimsizes(DATASETS2)-1) = ispan(1, dimsizes(DATASETS2), 1)

  ; set up needed arrays
  growth_series = new((/dimsizes(DATASETS), nyear/), float)
  growth_series!0 = "model"
  growth_series!1 = "time"
  growth_series&model = DATASETS(mapping_array1)
  growth_series&time = ispan(start_year, end_year, 1)

  delta_t = new((/dimsizes(DATASETS2)+1, nyear/), float)
  delta_t!0 = "model"
  delta_t!1 = "time"
  delta_t&model = array_append_record("OBS", DATASETS2, 0)
  delta_t&time = growth_series&time

  ; Temp Anomalies with reference to 2000 - 2010
  if start_year.lt.2003 then
    int_yr = 2003
  else
    int_yr = start_year
  end if
  if end_year.gt.2010 then
    fin_yr = 2010
  else
    fin_yr = end_year
  end if
  if (fin_yr - int_yr).lt.3 then
    fin_yr = end_year
    int_yr = start_year
  end if

end

begin
  ; Read Sat data
  obs_data = read_data(INFO0[ref_ind])
  obs_avg = area_operations(obs_data, lat_min, lat_max, 0, 360, \
                            "average", True)
  growth_series(0, :) = (/calc_gr(obs_avg, "yearly", min_nmonth)/)
  var0_units = obs_data@units

  ; Prepare tas data
  MODIS_path = diag_script_info@auxiliary_data_dir \
    + "/Land_Cover_Class_1degree.nc4"
  MODIS_dum = addfile(MODIS_path, "r")
  MOD_map = byte2flt(MODIS_dum->land_cover_class)
  MOD_map = where(MOD_map.eq.0 .or. MOD_map.eq.15 .or. MOD_map.eq.16,\
                  MOD_map@_FillValue, 1.)
  delete(MODIS_dum)
  ; This map is lon -180, 180, we want 0, 360
  MOD_map = lonFlip(MOD_map)

  ; ---------------------------------------------
  ; --------------------- First prepare XCO2 data
  do imod = 0, dimsizes(DATASETS) - 1
    if imod.eq.ref_ind then
      continue
    end if

    tmp = read_data(INFO0[imod])
    var0_mod = tofloat(tmp)
    copy_VarMeta(tmp, var0_mod)
    delete(tmp)

    var0_reg = area_hi2lores_Wrap(var0_mod&lon, var0_mod&lat, var0_mod, True, \
                                  1, obs_data&lon, obs_data&lat, False)
    var0_reg = where(ismissing(obs_data), obs_data, var0_reg)
    var0_avg = area_operations(var0_reg, lat_min, lat_max, 0, 360, \
                               "average", True)
    growth_series(mapping_array1(imod), :) = (/calc_gr(var0_avg, "yearly", \
                                                       min_nmonth)/)
    delete([/var0_mod, var0_avg, var0_reg/])
  end do

  ; ---------------------- Secondly tas data
  do imod = 0, dimsizes(DATASETS2)
    if imod.eq.(dimsizes(DATASETS2)) then
      tmp = read_data(INFO2[0])
      var1_units = tmp@units
    else
      tmp = read_data(INFO1[imod])
    end if

    var1_mod = tofloat(tmp)
    copy_VarMeta(tmp, var1_mod)
    delete(tmp)

    time_coord = cd_calendar(var1_mod&time, 0)
    temp_anom = var1_mod

    do imnth = 1, 12
      ref_inds = ind((time_coord(:, 1).eq.imnth) .and. \
                     (time_coord(:, 0).lt.fin_yr) .and. \
                     (time_coord(:, 0).ge.int_yr))
      mean = dim_avg_n_Wrap(var1_mod(ref_inds, :, :), 0)
      ind_mnth = ind(time_coord(:, 1).eq.imnth)
      mnth_anom = tofloat(var1_mod(ind_mnth, :, :))
      mnth_anom = mnth_anom - conform(mnth_anom, mean, (/1, 2/))
      temp_anom(ind_mnth, :, :) = (/mnth_anom/)
      delete(ref_inds)
      delete(ind_mnth)
      delete(mnth_anom)
      delete(mean)
    end do

    delete(time_coord)

    ; Apply MODIS mask to get rid of ocean, snow/ice and
    ; sparsely vegetated regions
    temp_hreg = linint2_Wrap(MOD_map&longitude, MOD_map&latitude, MOD_map, \
                             True, var1_mod&lon, var1_mod&lat, 0)
    temp_var = temp_anom*conform(temp_anom, temp_hreg, (/1, 2/))
    copy_VarMeta(var1_mod, temp_var)
    delete(temp_hreg)

    temp_area = area_operations(temp_var, lat_min, lat_max, 0, \
                                360, "average", False)

    ; First determine which part of the year to average over
    if lat_min.ge.0 then
    ; Schneising14 growing season: SH Dec-May, NH April-sept
      mnth_string = "AMJJAS"
    elseif lat_max.le.0 then
      mnth_string = "DJFMAM"
    else
      mnth_string = "annualclim"
    end if

    del_T = new(nyear, "float")
    do ii = 0, nyear - 1
      del_T(ii) = time_operations(temp_area, start_year + ii, \
                                  start_year + ii, "average", \
                                  mnth_string, False)
    end do
    delta_t(mapping_array2(imod), :) = (/del_T/)

    delete(del_T)
    delete(temp_var)
    delete(temp_anom)
    delete(var1_mod)

  end do

  ; -------------------------------------- PLOT -------------------------------
  obs_pres = True
  delta_t@units = var1_units
  delta_t@var = var1
  delta_t@diag_script = DIAG_SCRIPT
  delta_t@long_name = "Growing Season Temperature"
  growth_series@units = var0_units + " yr-1"
  growth_series@var = "gr"
  growth_series@diag_script = DIAG_SCRIPT
  growth_series@long_name = "Growth Rate"
  delta_t@region = region

  filename_dtgr_detr = var0 + "_" + var1 + "_" + experiments(0) \
    + "_" + region + "_" + (start_year) + "-" + (end_year) \
    + "_dt_GR_detr_panels"
  outfile_dtgr_detr = plot_dir + "/" + filename_dtgr_detr
  outfile_netcdf_dtgr_detr = work_dir + "/" + filename_dtgr_detr + ".nc"

  wks = gsn_open_wks(file_type, outfile_dtgr_detr)
  delta_t_detr = new(dimsizes(delta_t), float)
  copy_VarMeta(delta_t, delta_t_detr)
  growth_series_detr = new(dimsizes(growth_series), float)
  copy_VarMeta(growth_series, growth_series_detr)
  do imod = 0, dimsizes(growth_series&model) - 1
    delta_t_detr(imod, :) = var_detrend(delta_t&time, delta_t(imod, :))
    growth_series_detr(imod, :) = var_detrend(growth_series&time, \
                                              growth_series(imod, :))
  end do
  delta_t_detr@long_name = "IAV of growing season temperature"
  growth_series_detr@long_name = "IAV of growth rate"
  var_detr_GR = "Detrended GR"
  plt_dtgr_panels = panel_plots(wks, delta_t_detr, growth_series_detr, \
                                "Detrended ~F8~D~F~ T", var_detr_GR, \
                                obs_pres, INFO0)
  delete(wks)
  delete(plt_dtgr_panels)

  outfile_netcdf_dtgr_detr@existing = "overwrite"
  ncdf_outfile = ncdf_write(growth_series_detr, outfile_netcdf_dtgr_detr)
  outfile_netcdf_dtgr_detr@existing = "append"
  ncdf_outfile = ncdf_write(delta_t_detr, outfile_netcdf_dtgr_detr)

  log_provenance(ncdf_outfile, \
                 outfile_netcdf_dtgr_detr, \
                 "Sensitivity of interannual variability of " \
                 + var0 + " growth rate in the " + region + ", " + \
                 start_year + "-" + end_year + "to the interannual " \
                 + "variability of growing season temperature. " \
                 + "Similar to Gier et al. 2020, Fig C1.", \
                 (/"mean", "stddev"/), \
                 DOMAIN, \
                 "scatter", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

  leave_msg(DIAG_SCRIPT, "")

end
