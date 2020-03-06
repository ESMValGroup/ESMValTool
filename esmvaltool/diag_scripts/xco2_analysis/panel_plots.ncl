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
load "$diag_scripts/xco2_analysis/carbon_plots.ncl"

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

  if diag_script_info@masking[0] .ne. "none" then
    opt_mask = "_" + str_join(diag_script_info@masking, "_")
  else
    opt_mask = ""
  end if

    ;; Formatted varname for plots
  if var0.eq."xco2" then
      var0_plotname = "XCO~B~2~N~"
  else
      var0_plotname = var0
  end if

  ; Prepare arrays
  amp_series = new((/dim_MOD, nyear/), float)
  amp_series!0 = "model"
  amp_series&model = DATASETS
  amp_series!1 = "time"
  amp_series&time = ispan(start_year, end_year, 1)
  amp_series@units = "ppmv"
  amp_series@long_name = "Seasonal Cycle Amplitude"
  amp_series@var_yaxis = True
  amp_series@region = region

  var0_yr = amp_series
  var0_yr@long_name = var0_plotname

  varname2 = "SCA"
  varname2@mean = diag_script_info@plot_var2_mean

end
begin
  ; First read obs to have it ready for masking
  obs_data = read_data(INFO0[ref_ind])
  obs_data = convert_units(obs_data, "ppmv")
  time_mnth = cd_calendar(obs_data&time, 0)

  ; Prepare sat masks for usage if needed
  if any(diag_script_info@masking .eq. "sciamachy") then
    ; Sciamachy Masks! 2003-2008
    scia_ind = max(ind(time_mnth(:,0).eq.2008))
    scia_l_mask = obs_data(0, :, :)
    scia_l_mask = 0.

    do i_tim = 0, dimsizes(obs_data&time(:scia_ind)) - 1
        mnth_mv = where(ismissing(obs_data(i_tim, :, :)), 0, 1)
        scia_l_mask = scia_l_mask + mnth_mv
    end do
    scia_l_mask = scia_l_mask / tofloat(dimsizes(obs_data&time(:scia_ind)))
    scia_l_mask = where(scia_l_mask.ge.0.5, 1, scia_l_mask@_FillValue)
  end if

  if any(diag_script_info@masking .eq. "gosat") then
    gosat_ind = min(ind(time_mnth(:,0).eq.2013))
    gosat_l_mask = obs_data(0, :, :)
    gosat_l_mask = 0.

    do i_tim = gosat_ind, dimsizes(obs_data&time) - 1
        mnth_mv = where(ismissing(obs_data(i_tim, :, :)), 0, 1)
        gosat_l_mask = gosat_l_mask + mnth_mv
    end do
    gosat_l_mask = gosat_l_mask / tofloat(dimsizes(obs_data&time(gosat_ind:)))
    gosat_l_mask = where(gosat_l_mask.ge.0.5, 1, gosat_l_mask@_FillValue)
  end if

  ; Prepare Land-Sea Mask
  f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
  a = f->LSMASK
  sftlf = byte2flt(landsea_mask(a, obs_data&lat, obs_data&lon))
  sftlf = where(sftlf.gt.1., 1., sftlf)
  sftlf = where(sftlf.eq.0, sftlf@_FillValue, sftlf)

  ; Read data
  do imod = 0, dim_MOD -1
    ; Load data
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

    ; Go through all masks and apply -> check for which are applied to obs as well!
    if imod .ne. ref_ind then
      var0_reg = area_hi2lores_Wrap(var0_mod&lon, var0_mod&lat, var0_mod, True, \
                                  1, obs_data&lon, obs_data&lat, False)
      delete(var0_mod)
      if any(diag_script_info@masking .eq. "obs") then
        var0_reg = where(ismissing(obs_data), obs_data, var0_reg)
      end if
      if any(diag_script_info@masking .eq. "sciamachy") then
        var0_reg = var0_reg * conform(var0_reg, scia_l_mask, (/1, 2/))
      end if
      if any(diag_script_info@masking .eq. "gosat") then
        var0_reg = var0_reg * conform(var0_reg, gosat_l_mask, (/1, 2/))
      end if
    else
      var0_reg = obs_data
    end if

    if any(diag_script_info@masking .eq. "land") then
      var0_reg = var0_reg * conform(var0_reg, sftlf, (/1, 2/))
    end if

    ; Different averages
    if diag_script_info@area_avg .eq. "lat-first" then
      var0_lonavg = dim_avg_n_Wrap(var0_reg, 2)
      printVarSummary(var0_lonavg)
      var0_sca_lonavg = new((/nyear, dimsizes(var0_lonavg&lat)/), float)
      delete(var0_reg)

      do ilat = 0, dimsizes(var0_lonavg&lat) -1
        gr = calc_gr(var0_lonavg(:, ilat), "monthlyfilled", min_nmonth)
        do itim = 0, nyear - 1
          slope = dim_cumsum(gr(itim*12:itim*12+11)/12., 1)
          detrended = var0_lonavg(itim*12:itim*12+11, ilat) - slope
          mvnum = num(ismissing(detrended))
          if mvnum .le. min_nmonth then
              var0_sca_lonavg(itim, ilat) = max(detrended) - min(detrended)
          end if
          delete(detrended)
          delete(slope)
        end do

        ;var0_sca_lonavg(:, ilat) = (/calc_sca(var0_lonavg(:, ilat), min_nmonth)/)
      end do

      coords = ind_nearest_coord((/lat_min, lat_max/), var0_lonavg&lat, 0)
      cos_lats = tofloat(cos(var0_lonavg&lat(coords(0):coords(1))/180*3.14159))
      var0_lonavg_xco2 = dim_avg_wgt_n(var0_lonavg(:,coords(0): coords(1)), cos_lats, 1, 1)
      copy_VarMeta(var0_lonavg, var0_lonavg_xco2)

      var0_yr(mapping_array(imod), :) = (/time_operations(var0_lonavg_xco2, -1, -1, "average", "yearly", True)/)
      amp_series(mapping_array(imod), :) = dim_avg_wgt(var0_sca_lonavg(:,coords(0): coords(1)) , \
                cos_lats, 1)
    else
      var0_avg = area_operations(var0_reg, lat_min, lat_max, 0, 360, "average", True)
      var0_yr(mapping_array(imod), :) = (/time_operations(var0_avg, -1, -1, "average", "yearly", True)/)
      amp_series(mapping_array(imod), :) = (/calc_sca(var0_avg, min_nmonth)/)
    end if

  end do



  ; -------------------------------------------------------------------------------
  ; ---------------------- PANEL PLOTS---------------------------------------------
  ; -------------------------------------------------------------------------------
  obs_pres = diag_script_info@obs_in_panel
  if .not. obs_pres then
    var0_yr = var0_yr(1:, :)
    amp_series = amp_series(1:, :)
  end if
  outfile_dtgr_detr = plot_dir + "/" + var0 + "_" + experiments(0) \
               + "_" + region + "_" + (start_year) + "-" + (end_year) + "_" + var0 + "_SCA_panels" \
               + opt_mask
  wks = gsn_open_wks(file_type, outfile_dtgr_detr)
  plt_xco2_SCA_panels = panel_plots(wks, var0_yr, amp_series, var0, \
                                    varname2, obs_pres, INFO0)
  delete(wks)
  delete(plt_xco2_SCA_panels)

end
