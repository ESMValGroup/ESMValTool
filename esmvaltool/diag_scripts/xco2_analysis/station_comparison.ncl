; #############################################################################
;  xco2_analysis/station_comparison.ncl
; #############################################################################
; Description
;    Sampling Station CO2 data like XCO2 datasets and comparing the two
;
; Required diag_script_info attributes:
;    var_order: in this case xco2, co2, co2s - column averaged with obs dataset
;               first, then 2D variable, followed by surface stations
;
; Optional diag_script_info attributes:
;    output_file_type: output file type for plots. Default: png
;    var_plotnames: String formatting how variables should be named in plots
;                   defaults to short_name if not assigned
;    overwrite_altitudes: Give other altitude values than the ones attached in
;                         the station data. Valid if altitude changes and
;                         timeseries spans range with different sample
;                         altitude. Caveat: If used, need to give altitude
;                         values for all stations.
;    output_map: boolean if stations to be displayed on map. As this requires
;                finetuning, currently only implemented for station set of
;                (ASK, CGO, HUN, LEF, WIS) following the paper. Change code
;                if others are desired.
;
; Caveats
;    Currently assumes same time for all input variables (mmm computation)
;
; Modification history
;    20201210-gier_bettina: Added provenance and clean up
;    20200220-gier_bettina: Adapted to version 2
;    20180809-gier_bettina: written.
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
  var0 = diag_script_info@var_order(0)
  var1 = diag_script_info@var_order(1)
  var2 = diag_script_info@var_order(2)
  dim_VAR = ListCount(variable_info)

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  INFO1 = select_metadata_by_name(input_file_info, var1)
  INFO2 = select_metadata_by_name(input_file_info, var2)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  DATASETS1 = metadata_att_as_array(INFO1, "dataset")
  ; ESRL station name is in version
  DATASETS2 = metadata_att_as_array(INFO2, "version")
  IN_FILES0 = metadata_att_as_array(INFO0, "filename")
  MOD_FILES = array_append_record(IN_FILES0, \
                                  metadata_att_as_array(INFO1, "filename"), 0)
  STAT_FILES = metadata_att_as_array(INFO2, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + " and " + var1 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

end

begin
  ; Plot file type
  if (isatt(diag_script_info, "output_file_type")) then
    file_type = diag_script_info@output_file_type
  elseif (isatt(config_user_info, "output_file_type")) then
    file_type = config_user_info@output_file_type
  else
    file_type = "png"
  end if

  ; Formatted varname for plots
  if (isatt(diag_script_info, "var_plotnames")) then
      var_plotnames = diag_script_info@var_plotnames
  else
      var_plotnames = (/var0, var1/)
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

end

begin
  ; Read Sat data
  if (isatt(INFO0[0], "reference_dataset")) then
    ref_ind = ind(DATASETS.eq.INFO0[0]@reference_dataset)
  else
    error_msg("w", DIAG_SCRIPT, "", "this diagnostic requires reference " + \
              "datasets.")
  end if
  sat_data = read_data(INFO0[ref_ind])
  var_unit = sat_data@units

  ; read multi-model means
  xco2_mmm = ind(DATASETS.eq."MultiModelMean")
  mod_xco2 = read_data(INFO0[xco2_mmm])

  co2_mmm = ind(DATASETS1.eq."MultiModelMean")
  mod_co2 = read_data(INFO1[co2_mmm])

  ; Make array to hold stationdata
  station_array = new((/dimsizes(DATASETS2), ntime/), float)
  station_array!0 = "station"
  station_array&station = DATASETS2
  station_array@lat = new(dimsizes(DATASETS2), float)
  station_array@lon = new(dimsizes(DATASETS2), float)
  station_array@alt = new(dimsizes(DATASETS2), float)
  station_array!1 = "time"
  station_array&time = time
  station_array@var = var2
  station_array@units = var_unit
  station_array@diag_script = DIAG_SCRIPT

  ; Make arrays to hold masked station data and corresponding sat measurement
  station_array_masked = station_array
  station_array_masked = station_array
  station_array_masked@var = station_array@var + "_masked"
  sat_array = station_array
  sat_array@var = var0
  co2_array = station_array
  xco2_array = station_array
  sat_array@var = var0 + "_model"
  offsets = new(dimsizes(DATASETS2), float)

  ; Reading in station data
  do istat = 0, dimsizes(DATASETS2) - 1
    stat_data = read_data(INFO2[istat])
    station_array(istat, :) = (/stat_data/)
    station_array@lat(istat) = tofloat(stat_data&lat(0))
    station_array@lon(istat) = tofloat(stat_data&lon(0))
    if (isatt(diag_script_info, "overwrite_altitudes")) then
      station_array@alt(istat) = \
        tofloat(diag_script_info@overwrite_altitudes(istat))
    else
      station_array@alt(istat) = tofloat(stat_data@altitude)
    end if

    ; Sample Station data according to Sat Obs
    add_lat = 2.5
    add_lon = 2.5
    latmin = station_array@lat(istat) - add_lat
    latmax = station_array@lat(istat) + add_lat
    lonmin = station_array@lon(istat) - add_lon
    lonmax = station_array@lon(istat) + add_lon
    if lonmin.lt.0 then
        lonmin = 360.0 - abs(lonmin)
    end if
    if lonmax .gt.360 then
        lonmax = lonmax - 360.0
    end if
    obs_area = area_operations(sat_data, latmin, latmax, \
                               lonmin, lonmax, "average", True)
    sat_array(istat, :) = (/obs_area/)
    xco2_area = area_operations(mod_xco2, latmin, latmax, \
                                lonmin, lonmax, "average", True)

    ; Calculate model offset through regression
    sat_rc = regline(obs_area&time, obs_area)
    mod_rc = regline(xco2_area&time, xco2_area)
    offset_mod = (obs_area&time(0) + \
                  (obs_area&time(dimsizes(obs_area&time)-1) \
                  - obs_area&time(0))/2) * sat_rc + sat_rc@yintercept \
      - ((xco2_area&time(0) + (xco2_area&time(dimsizes(xco2_area&time) - 1)\
          - xco2_area&time(0))/2)*mod_rc + mod_rc@yintercept)
    offsets(istat) = tofloat(offset_mod)
    xco2_array(istat, :) = (/xco2_area/) + tofloat(offset_mod)

    atmos = stdatmus_z2tdp(station_array@alt(istat))
    p_level = atmos(2)*100
    if p_level.gt.100000 then
        p_level = 100000
    end if

    mod_co2h_data = mod_co2(:, {p_level}, :, :)
    co2_area = area_operations(mod_co2h_data, latmin, latmax, \
                               lonmin, lonmax, "average", True)
    co2_array(istat, :) = (/co2_area/) + tofloat(offset_mod)
    station_array_masked(istat, :) = \
      where(ismissing(obs_area), \
            station_array_masked(istat, :), \
            station_array(istat, :))
    delete(obs_area)
  end do  ; End reading station data

  stobs_ts = new((/dimsizes(DATASETS2), 4, dimsizes(station_array&time)/), \
                 float)
  stobs_ts!0 = "station"
  stobs_ts&station = station_array&station
  stobs_ts!1 = "dataset"
  stobs_ts&dataset = (/"Station " + var_plotnames(1), \
                       "Satellite " + var_plotnames(0), \
                       "Multi-model Mean " + var_plotnames(1), \
                       "Multi-Model Mean " + var_plotnames(0)/)
  stobs_ts!2 = "time"
  stobs_ts&time = station_array&time

  ; Loop over stations and plot co2 timeseries
  do istat = 0, dimsizes(DATASETS2) - 1

    filename_stations_masked = var0 + "_surface_flasks_" + \
      station_array&station(istat) + "_" + (start_year) + "-" + \
      (end_year) + "_masked"
    outfile_stations_masked = plot_dir + "/" + filename_stations_masked
    outfile_stations_masked_netcdf = work_dir + "/" \
      + filename_stations_masked + ".nc"

    outfile_stations_masked_netcdf@existing = "overwrite"
    ncdf_outfile = ncdf_write(station_array(istat, :), \
                              outfile_stations_masked_netcdf)
    outfile_stations_masked_netcdf@existing = "append"
    ncdf_outfile = ncdf_write(station_array_masked(istat, :), \
                              outfile_stations_masked_netcdf)

    wks = gsn_open_wks(file_type, outfile_stations_masked)
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tiYAxisString =  var_plotnames(1) + " [" + var_unit + "]"
    res@tiMainString = "Time series for surface flask measurement"
    res@xyLineThicknessF       = 2.0      ; default is 1
    res@xyMonoDashPattern      = True     ; force all solid lines
    res@pmLegendOrthogonalPosF = -1.05
    res@pmLegendParallelPosF   = 0.25
    res@pmLegendDisplayMode    = "Always"
    res@pmLegendWidthF         = 0.05
    res@pmLegendHeightF        = 0.2
    res@lgLabelFontHeightF     = .015
    res@xyDashPattern = 0
    res@tmXBFormat = "f"  ; Don't put unnecessary 0's
    res@xyExplicitLegendLabels = station_array&station(istat) + " (" + \
      station_array@lat(istat) + ", " + station_array@lon(istat) + ")"
    res@xyLineColors = (/"blue", "black", "red", "orange"/)
    res@tiMainString = "Time series masked " + station_array&station(istat) \
      + " (" + station_array@lat(istat) + ", " + station_array@lon(istat) + ")"
    res@pmLegendDisplayMode = "Never"
    statm_plot = gsn_csm_xy(wks, station_array_masked&time, \
                            station_array_masked(istat, :), res)
    res@xyDashPattern = 2
    res@pmLegendDisplayMode = "Never"
    stat_plot = gsn_csm_xy(wks, station_array&time, \
                           station_array(istat, :), res)
    NhlAddOverlay(statm_plot, stat_plot, statm_plot)
    draw(statm_plot)
    frame(wks)
    delete(stat_plot)
    delete(statm_plot)
    delete(res@xyExplicitLegendLabels)

    print(INFO2[istat]@filename)

    ; Provenance
    log_provenance(ncdf_outfile, \
                   outfile_stations_masked, \
                   "Time series of station " + station_array&station(istat) \
                   + ", dotted denoting areas masked like observations.", \
                   (/""/), \
                   "reg", \
                   "times", \
                   AUTHORS, \
                   REFERENCES, \
                   (/INFO2[istat]@filename, INFO0[ref_ind]@filename/))

    filename_stations_obs = var0 + "_surface_flasks_" + \
      station_array&station(istat) + "_" + (start_year) + "-" + \
      (end_year) + "_maskedvsobs"
    outfile_stations_obs_netcdf = work_dir + "/" \
      + filename_stations_obs + ".nc"

    outfile_stations_obs = plot_dir + "/" + filename_stations_obs
    wks = gsn_open_wks(file_type, outfile_stations_obs)
    res@tiMainString = "Time series " + station_array&station(istat) + \
      " (" + station_array@lat(istat) + ", " + station_array@lon(istat) + \
      ") Mod Offset: " + round(offsets(istat), 3)
    res@xyDashPattern = 0
    res@pmLegendDisplayMode = "Always"
    stobs = new((/4, dimsizes(station_array&time)/), float)
    stobs!0 = "dataset"
    stobs&dataset = (/station_array&station(istat), \
                     "Satellite " + var_plotnames(0), \
                     "Multi-model Mean " + var_plotnames(1), \
                     "Multi-Model Mean " + var_plotnames(0)/)
    stobs!1 = "time"
    stobs&time = station_array&time
    stobs(0, :) = (/station_array(istat, :)/)
    stobs(1, :) = (/sat_array(istat, :)/)
    stobs(2, :) = (/co2_array(istat, :)/)
    stobs(3, :) = (/xco2_array(istat, :)/)
    stobs@var = var2
    stobs@diag_script = DIAG_SCRIPT
    stobs@units = var_unit
    stobs_ts(istat, :, :) = (/stobs/)
    ncdf_outfile = ncdf_write(stobs, outfile_stations_obs_netcdf)

    res@xyExplicitLegendLabels = stobs&dataset
    stat_plot = gsn_csm_xy(wks, stobs&time, stobs, res)
    draw(stat_plot)
    frame(wks)
    delete(res@xyExplicitLegendLabels)

    ; Provenance
    log_provenance(ncdf_outfile, \
                   outfile_stations_obs, \
                   "Time series of station " + station_array&station(istat) \
                   + ", for " + var0 + " and " + var1 + ".", \
                   (/"mean"/), \
                   "reg", \
                   "times", \
                   AUTHORS, \
                   REFERENCES, \
                   array_append_record(MOD_FILES, INFO2[istat]@filename, 0))

  end do
  delete(res)

  if ((isatt(diag_script_info, "output_map")) \
      .and. (diag_script_info@output_map)) then
    ; Make Mapplot and put timeseries on it

    ; Empty dummy map plot:

    filename_map = var0 + "_surface_flasks_" + (start_year) \
      + "-" + (end_year) + "_map"
    outfile_map_netcdf = work_dir + "/" + filename_map + ".nc"
    outfile_map = plot_dir + "/" + filename_map
    stobs_ts@var = var2
    stobs_ts@diag_script = DIAG_SCRIPT
    stobs_ts@units = var_unit
    ncdf_outfile = ncdf_write(stobs_ts, outfile_map_netcdf)

    wks = gsn_open_wks(file_type, outfile_map)

    res               = True                        ; plot mods desired
    res@gsnDraw       = False                       ; don't draw yet
    res@gsnFrame      = False                       ; don't advance frame

    resm = res
    resm@mpFillOn = True
    resm@mpLandFillColor = "grey92"
    plot = gsn_csm_map_ce(wks, resm)

    ; Make Station timeseries plots
    sta_ts_plots = new(dimsizes(DATASETS2), graphic)

    do istat = 0, dimsizes(DATASETS2) - 1
      res@tiMainString = station_array&station(istat) + " (" + \
           decimalPlaces(station_array@lat(istat), 1, True) + \
           ", " + decimalPlaces(station_array@lon(istat), 1, True) + \
           ") Offset: " + decimalPlaces(offsets(istat), 1, True)
      res@vpWidthF          = 0.155       ; Make plots quite a bit smaller.
      res@vpHeightF         = 0.15

      res@trXMinF = start_year
      res@trXMaxF = end_year + 1
      res@trYMinF = min(stobs_ts(istat, :, :)) - 2
      res@trYMaxF = max(stobs_ts(istat, :, :)) + 2
      res@xyDashPattern = 0
      res@xyLineColors = (/"blue", "black", "red", "orange"/)
      res@tmXBFormat = "f"  ; Don't put unnecessary 0's
      res@pmLegendDisplayMode = "Never"

      res@tiYAxisString = var_plotnames(1) + " [" + var_unit + "]"
      res@xyLineThicknessF       = 2.0      ; default is 1
      res@xyMonoDashPattern      = True     ; force all solid lines
      res@tiMainFontHeightF = 0.01
      res@tiYAxisFontHeightF = 0.01
      res@tmXBLabelFontHeightF = 0.007
      res@tmYLLabelFontHeightF = 0.007
      res@tmYLLabelDeltaF = -0.5
      res@tmXBLabelDeltaF = -0.5
      res@tiDeltaF = 0.9

      sta_ts_plots(istat) = gsn_csm_xy(wks, stobs_ts&time, \
                                       stobs_ts(istat, :, :), res)

      mres = True
      mres@gsMarkerIndex = 5
      mres@gsMarkerSizeF    = 0.01
      mres@gsMarkerColor    = "black"
      mres@gsMarkerThicknessF = 5.
      plot@$unique_string("marker")$ = \
        gsn_add_polymarker(wks, plot, station_array@lon(istat), \
                           station_array@lat(istat), mres)
    end do

    txres = True
    txres@txFontHeightF = 0.01

    ; Manual text placement! If using different station please change.
    plot@$unique_string("test")$ = \
      gsn_add_text(wks, plot, stobs_ts&station(0) + " (" \
                   + station_array@alt(0) \
                   + " m)", station_array@lon(0) + 20, \
                   station_array@lat(0) - 7, txres)

    plot@$unique_string("test")$ = \
      gsn_add_text(wks, plot, stobs_ts&station(1) + " (" \
                   + station_array@alt(1) \
                   + " m)", station_array@lon(1) + 10, \
                   station_array@lat(1) - 7, txres)

    plot@$unique_string("test")$ = \
      gsn_add_text(wks, plot, stobs_ts&station(2) + " (" \
                   + station_array@alt(2) \
                   + " m)", station_array@lon(2) + 5, \
                   station_array@lat(2) + 7, txres)

    plot@$unique_string("test")$ = \
      gsn_add_text(wks, plot, stobs_ts&station(3) + " (" \
                   + station_array@alt(3) \
                   + " m)", station_array@lon(3) - 20, \
                   station_array@lat(3) - 7, txres)

    plot@$unique_string("test")$ = \
      gsn_add_text(wks, plot, stobs_ts&station(4) + " (" \
                   + station_array@alt(4) \
                   + " m)", station_array@lon(4), \
                   station_array@lat(4) - 7, txres)

    amres                  = True
    amres@amJust           = "BottomLeft"

    ; ASK
    amres@amOrthogonalPosF =  0.43   ; 0.5 is the bottom edge of the plot.
    amres@amParallelPosF   = -0.10   ; -0.5 is the left edge of the plot.
    anno1 = gsn_add_annotation(plot, sta_ts_plots(0), amres)

    ; CGO
    amres@amOrthogonalPosF =  0.45
    amres@amParallelPosF  = 0.18
    anno2 = gsn_add_annotation(plot, sta_ts_plots(1), amres)

    ; HUN
    amres@amOrthogonalPosF =  -0.06
    amres@amParallelPosF  =  -0.19
    anno3 = gsn_add_annotation(plot, sta_ts_plots(2), amres)

    ; LEF
    amres@amOrthogonalPosF =  0.30
    amres@amParallelPosF  =  -0.445
    anno3 = gsn_add_annotation(plot, sta_ts_plots(3), amres)

    ; WIS
    amres@amOrthogonalPosF =  -0.05
    amres@amParallelPosF  =  0.21
    anno3 = gsn_add_annotation(plot, sta_ts_plots(4), amres)

    ; Make shared legend

    gres = True
    gres@YPosPercent = 97.
    gres@XPosPercent = 1.
    gres@LineLabelWhiteSpacePercent = 1.

    lineres = True
    lineres@lgLineColors = (/"blue", "black", "red", "orange"/)
    lineres@lgLineThicknesses = 2.5
    lineres@LineLengthPercent = 4.

    textres = True
    textres@lgLabelFontHeights = 0.01
    textres@lgLabels = stobs_ts&dataset  ; legend labels (required)

    plot = simple_legend(wks, plot, gres, lineres, textres)
    draw(plot)
    frame(wks)

        ; Provenance
    log_provenance(ncdf_outfile, \
                   outfile_map, \
                   "Time series for satellite, multi-model mean and station " \
                   + var0 + " and " + var1 \
                   + ". Similar to Gier et al, 2020 Fig. 2", \
                   (/"mean"/), \
                   "global", \
                   "times", \
                   AUTHORS, \
                   REFERENCES, \
                   array_append_record(MOD_FILES, STAT_FILES, 0))
  end if

  leave_msg(DIAG_SCRIPT, "")
end
