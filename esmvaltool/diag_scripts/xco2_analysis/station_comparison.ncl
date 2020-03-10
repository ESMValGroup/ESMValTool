; #############################################################################
;  xco2_analysis/station_comparison.ncl
; #############################################################################
; Description
;    Sampling Station CO2 data like XCO2 datasets and comparing the two
;
; Required diag_script_info attributes:
;    path: path to station data - move this to recipe!
;    stationlist: ascii data for stations -> make this into recipe cfg
;    station_varname: short_name of variable for station data
;    var_order: in this case xco2, co2 - column averaged with obs dataset from recipe first
;
; Optional diag_script_info attributes:
;
; Caveats
;    Currently assumes same time for all input variables (mmm computation)
;
; Modification history
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
  REFERENCES = (/"acknow_author"/)
  DOMAIN = (/"global"/) ;diag_script_info@region

  ; Variable
  var0 = diag_script_info@var_order(0)
  var1 = diag_script_info@var_order(1)
  dim_VAR = ListCount(variable_info)

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  INFO1 = select_metadata_by_name(input_file_info, var1)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  ALL_FILES = metadata_att_as_array(INFO0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + " and " + var1 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

end

begin
  ; Maximum amount of missing months allowed for calculating SCA
  min_nmonth = 7 ;5

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
    IYEAR = min(start_years)
    end_years = metadata_att_as_array(INFO0, "end_year")
    FYEAR = max(end_years)
    nyear = FYEAR - IYEAR + 1
    ntime = nyear*12
    time = new(ntime, float)
    do yy = IYEAR, FYEAR
        do mm = 0, 11
            time(12* (yy - IYEAR) + mm) = yy + (mm+0.5)/12.
        end do
    end do

end

begin
  varp = "CO2"
  ; Read Sat data
  if (isatt(INFO0[0], "reference_dataset")) then
    ref_ind = ind(DATASETS.eq.INFO0[0]@reference_dataset)
  end if
  sat_data = read_data(INFO0[ref_ind])

    ; Make mmm array
    mod_xco2_hreg = sat_data
    mod_xco2_hreg = 0.

    ; Compute multi-model mean (regrid and avg)
    do ivar = 0, dim_VAR - 1
      if ivar.eq.0 then
        INFO = INFO0
      else
        INFO = INFO1
        ref_ind = -999
      end if

      dim_MOD = dimsizes(metadata_att_as_array(INFO, "dataset"))

        ; Loop over xco2 + co2, regrid
        do imod = 0, dim_MOD - 1
        ; Sat data not inside MMM
        if ivar.eq.0 .and. imod.eq.ref_ind then
          continue
        end if

        mod_var= read_data(INFO[imod])
        var_reg = area_hi2lores_Wrap(mod_var&lon, mod_var&lat, mod_var, \
                                             True,1, sat_data&lon, sat_data&lat, False)

        if ivar.eq.1 then
          if .not. isdefined("mod_co2") then
            print("defining mod_co2")
            mod_co2 = var_reg
          else
            if dimsizes(var_reg&plev) .ne. dimsizes(mod_co2&plev) then
              var_zreg = int2p_n_Wrap(var_reg&plev, var_reg, mod_co2&plev, -1, 1)
              print("adding to modco2 with diff plev")
              mod_co2 = mod_co2 + (/tofloat(var_zreg)/)
              delete(var_zreg)
            else
              print("adding to modco2 with same plev")
              mod_co2 = mod_co2 +  (/tofloat(var_reg)/)
            end if
          end if
        else if ivar.eq.0 then
          mod_xco2_hreg = mod_xco2_hreg + (/tofloat(var_reg)/)
        end if
        end if
        delete([/mod_var, var_reg/])
      end do ; Loop over models

    end do ; Loop over variables

    dim_MOD = dimsizes(metadata_att_as_array(INFO0, "dataset"))

    mod_xco2_hreg = mod_xco2_hreg / tofloat(dim_MOD -1 )
    mod_co2 = mod_co2 / tofloat(dim_MOD -1)

    mod_xco2 = where(ismissing(sat_data), sat_data, mod_xco2_hreg)
    copy_VarMeta(sat_data, mod_xco2)

    ; Convert units to ppmv
    mod_co2 = convert_units(mod_co2, "ppmv")
    mod_xco2 = convert_units(mod_xco2, "ppmv")
    sat_data = convert_units(sat_data, "ppmv")

  ; Print stuffs!
  print("Mod CO2 Avg " + avg(mod_co2))
  print("Sat XCO2 Avg " + avg(sat_data))
  print("Mod XCO2 Avg " + avg(mod_xco2))



    ; Read station info
    log_info("Reading station information")
    stationlist = diag_script_info@path + diag_script_info@stationlist
    stat = readAsciiTable(stationlist, 1, "string", 0)
    st_code = str_squeeze(str_get_field(stat(:, 0), 1, ","))
    st_name = str_squeeze(str_get_field(stat(:, 0), 2, ","))
    st_lat  = tofloat(str_get_field(stat(:, 0), 3, ","))
    st_lon  = tofloat(str_get_field(stat(:, 0), 4, ","))
    st_alt  = toint(str_get_field(stat(:, 0), 5, ","))
    YEAR1   = toint(str_get_field(stat(:, 0), 6, ","))
    YEAR2   = toint(str_get_field(stat(:, 0), 7, ","))
    nstat   = dimsizes(st_code)
    delete(stat)

    ; Make array to hold stationdata
    station_array = new((/nstat, ntime/), float)
    station_array!0 = "station"
    station_array&station = st_code
    station_array&station@name = st_name
    station_array&station@lat = decimalPlaces(st_lat, 1, True)
    station_array&station@lon = decimalPlaces(st_lon, 1, True)
    station_array&station@alt = st_alt
    station_array!1 = "time"
    station_array&time = time

    ; Make arrays to hold masked station data and corresponding sat measurement
    station_array_masked = station_array
    sat_array = station_array
    co2_array = station_array
    xco2_array = station_array
    offsets = new(nstat, float)


    obsfiles = diag_script_info@path + "OBS_ESRL_surface_flask_" + st_code + \
               "_T0M_" + diag_script_info@station_varname + "_" + YEAR1 + "01-" + \
               YEAR2 + "12.nc"
    ; Reading in station data
    do istat=0, nstat-1
        var = addfile(obsfiles(istat), "r")
        statdata = var->$diag_script_info@station_varname$

        ; Extract correct timestep
        timestat = cd_calendar(statdata&time, 0)
        ind_i = ind(timestat(:,0).eq.IYEAR)
        ind_f = ind(timestat(:,0).eq.FYEAR)
        station_array(istat, :) = (/statdata(min(ind_i):max(ind_f))/)
        delete([/timestat, ind_i, ind_f, var, statdata/])

        ; Sample Station data according to Sat Obs
        add_lat = 2.5
        add_lon = 2.5
        ;lat = ind_nearest_coord(st_lat(istat), sat_data&lat, 0)
        ;lon = ind_nearest_coord(st_lon(istat), sat_data&lon, 0)
        latmin = st_lat(istat) - add_lat
        latmax = st_lat(istat) + add_lat
        lonmin = st_lon(istat) - add_lon
        lonmax = st_lon(istat) + add_lon
        if lonmin.lt.0 then
            lonmin = 360.0 - abs(lonmin)
        end if
        if lonmax .gt.360 then
            lonmax = lonmax - 360.0
        end if
        obs_area = area_operations(sat_data, latmin, latmax, \
                                   lonmin, lonmax, "average", True)
        ;obs_area = dim_avg_n_Wrap(sat_data(:, {latmin:latmax}, :), (/1,2/))
        sat_array(istat, :) = (/obs_area/)
        xco2_area = area_operations(mod_xco2, latmin, latmax, \
                                   lonmin, lonmax, "average", True)

        ; Calculate model offset through rgression
        sat_rc = regline(obs_area&time, obs_area)
        mod_rc = regline(xco2_area&time, xco2_area)
        offset_mod = (obs_area&time(0) + \
                     (obs_area&time(dimsizes(obs_area&time)-1) \
                      - obs_area&time(0))/2) * sat_rc + sat_rc@yintercept \
                     - ((xco2_area&time(0) \
                         + (xco2_area&time(dimsizes(xco2_area&time) - 1)\
                            - xco2_area&time(0))/2)*mod_rc + mod_rc@yintercept)
        offsets(istat) = tofloat(offset_mod)
        xco2_array(istat,:) = (/xco2_area/) + tofloat(offset_mod)

        atmos = stdatmus_z2tdp(st_alt(istat))
        p_level = atmos(2)*100
        if p_level.gt.100000 then
            p_level = 100000
        end if

        mod_co2h_data = mod_co2(:, {p_level}, :, :)
        ;mod_co2h = where(ismissing(sat_data), sat_data, mod_co2h_data)
        ;copy_VarMeta(mod_co2h_data, mod_co2h)
        co2_area = area_operations(mod_co2h_data, latmin, latmax, \
                                   lonmin, lonmax, "average", True)
        co2_array(istat, :) = (/co2_area/) + tofloat(offset_mod)
        station_array_masked(istat, :) = \
            where(ismissing(obs_area), \
            station_array_masked(istat, :), \
            station_array(istat, :))
         delete(obs_area)

     end do


     stobs_ts = new((/nstat, 4, dimsizes(station_array&time)/), float)
     stobs_ts!0 = "station"
     stobs_ts&station = st_code
     stobs_ts!1 = "dataset"
     stobs_ts&dataset = (/"Station CO~B~2~N~", "Satellite XCO~B~2~N~", \
                           "Multi-model Mean CO~B~2~N~", "Multi-Model Mean XCO~B~2~N~"/)
     stobs_ts!2 = "time"
     stobs_ts&time = station_array&time

     ; Loop over stations and plot co2 timeseries
     do istat = 0, nstat-1
         outfile_stations = plot_dir + "/" + var0 + "_surface_flasks_" + st_code(istat) + "_" + (IYEAR) + "-" + (FYEAR)
         wks = gsn_open_wks(file_type, outfile_stations)
         res = True
         res@gsnDraw = False
         res@gsnFrame = False
         res@tiYAxisString =  varp + " [ ppmv ]"
         res@tiMainString = "Time series for surface flask measurement"
         res@xyLineThicknessF       = 2.0      ; default is 1
         res@xyMonoDashPattern      = True     ; force all solid lines
                                ;---Make legend
         res@pmLegendOrthogonalPosF = -1.05
         res@pmLegendParallelPosF   = 0.25
         res@pmLegendDisplayMode    = "Always"
         res@pmLegendWidthF         = 0.05
         res@pmLegendHeightF        = 0.2
         res@lgLabelFontHeightF     = .015
         res@xyDashPattern = 0

                                ;---Change the legend line labels
         res@xyExplicitLegendLabels = station_array&station(istat) + " (" + \
             station_array&station@lat(istat) + \
             ", " + station_array&station@lon(istat) + ")"

         ;res@pmLegendDisplayMode    = "Always"
         stat_plot = gsn_csm_xy(wks,station_array&time, station_array(istat,:), res)
         res = True
         draw(stat_plot)
         frame(wks)
         delete(stat_plot)

         outfile_stations_masked = plot_dir + "/" + var0 + "_surface_flasks_" + \
             st_code(istat) + "_" + (IYEAR) + "-" + \
             (FYEAR) + "_masked2"
         wks = gsn_open_wks(file_type, outfile_stations_masked)
         res@tiMainString = "Time series masked " + st_name(istat) + " (" + \
             station_array&station@lat(istat) + \
             ", " + station_array&station@lon(istat) + ")"
         res@pmLegendDisplayMode = "Never"
         statm_plot = gsn_csm_xy(wks,station_array_masked&time, station_array_masked(istat,:), res)
         res@xyDashPattern = 2
         res@pmLegendDisplayMode = "Never"
         stat_plot = gsn_csm_xy(wks,station_array&time, station_array(istat,:), res)
         NhlAddOverlay(statm_plot, stat_plot, statm_plot)
         draw(statm_plot)
         frame(wks)
         delete(stat_plot)
         delete(statm_plot)
         delete(res@xyExplicitLegendLabels)

         outfile_stations_obs = plot_dir + "/" + var0 + "_surface_flasks_" + \
             st_code(istat) + "_" + (IYEAR) + "-" + \
             (FYEAR) + "_maskedvsobs"
         wks = gsn_open_wks(file_type, outfile_stations_obs)
         res@tiMainString = "Time series " + st_name(istat) + " (" + \
             station_array&station@lat(istat) + \
             ", " + station_array&station@lon(istat) + \
             ") Mod Offset: " + round(offsets(istat), 3)
         ;res@xyMonoDashPattern = False
         res@xyDashPattern = 0;(/0, 2, 0, 0, 0/)                         ; Make curves all solid
         res@xyLineColors = (/"blue", "black", "red","orange"/) ; 3 different colors
         res@tmXBFormat = "f" ; Don't put unnecessary 0's
         res@pmLegendDisplayMode = "Always"
         stobs = new((/4, dimsizes(station_array&time)/), float)
         stobs!0 = "dataset"
         stobs&dataset = (/st_name(istat), "Satellite XCO~B~2~N~", \
                           "Multi-model Mean CO~B~2~N~", "Multi-Model Mean XCO~B~2~N~"/)
         stobs!1 = "time"
         stobs&time = station_array&time
         stobs(0, :) = (/station_array(istat, :)/)
         stobs(1, :) = (/sat_array(istat, :)/)
         stobs(2, :) = (/co2_array(istat, :)/)
         stobs(3, :) = (/xco2_array(istat, :)/)
         stobs_ts(istat, :, :) = (/stobs/)
         res@xyExplicitLegendLabels = stobs&dataset
         stat_plot = gsn_csm_xy(wks, stobs&time, stobs, res)
         draw(stat_plot)
         frame(wks)
         delete(res@xyExplicitLegendLabels)
     end do

     ; Compute SCA

     stat_sca = new((/nstat, nyear/), float)
     stat_sca!0 = "station"
     stat_sca&station = station_array&station
     stat_sca!1 = "time"
     stat_sca&time = ispan(IYEAR, FYEAR, 1)
     statm_sca = stat_sca
     sat_sca = stat_sca
     do istat = 0, nstat - 1
         stat_sca(istat, :) = (/calc_sca(station_array(istat, :), 5)/)
         statm_sca(istat, :) = (/calc_sca(station_array_masked(istat, :), 5)/)
         sat_sca(istat, :) = (/calc_sca(sat_array(istat, :), 5)/)
     end do

     ; Loop over stations and plot sca timeseries
     do istat = 0, nstat-1
         outfile_stats = plot_dir + "/" + var0 + "_surface_flasks_" + \
             st_code(istat) + "_" + (IYEAR) + "-" + \
             (FYEAR) + "_SCA_comparison"
         wks = gsn_open_wks(file_type, outfile_stats)
         res@trXMinF = IYEAR
         res@trXMaxF = FYEAR
         res@tiYAxisString = "SCA [ppmv]"
         res@trYMinF = min((/sat_sca(istat, :), stat_sca(istat, :), statm_sca(istat, :)/))
         res@trYMaxF = max((/sat_sca(istat, :), stat_sca(istat, :), statm_sca(istat, :)/))
         res@xyDashPattern = 0
         res@xyLineColor = "black"
         res@tiMainString = "SCA comparison " + st_name(istat) + " (" + \
             station_array&station@lat(istat) + \
             ", " + station_array&station@lon(istat) + ")"
         res@pmLegendDisplayMode = "Never"
         stat_plot = gsn_csm_xy(wks,stat_sca&time, stat_sca(istat,:), res)
         resL = True
         resL@gsLineColor = "black"
         resL@gsLineThicknessF = 2
         resL@gsLineDashPattern = 2
         stat_plot@$unique_string("lines")$ = gsn_add_polyline(wks, \
              stat_plot, statm_sca&time, statm_sca(istat,:), resL)
         resL@gsLineDashPattern = 0
         resL@gsLineColor = "blue"
         stat_plot@$unique_string("lines")$ = gsn_add_polyline(wks, \
              stat_plot, sat_sca&time, sat_sca(istat,:), resL)
         draw(stat_plot)
         frame(wks)
         delete(stat_plot)
     end do

     delete(res)

     ; Make Mapplot and put timeseries on it

     ; Empty dummy map plot:

     outfile_map = plot_dir + "/" + var0 + "_surface_flasks_" + (IYEAR) + "-" + \
             (FYEAR) + "_map"
     wks = gsn_open_wks(file_type, outfile_map)               ; send graphics to PNG file

     res               = True                        ; plot mods desired
     res@gsnDraw       = False                       ; don't draw yet
     res@gsnFrame      = False                       ; don't advance frame
     plot = gsn_csm_map_ce(wks,res)                 ; draw global map
     ;drawNDCGrid(wks)

     ; Make Station timeseries plots
     sta_ts_plots = new(nstat, graphic)

     do istat = 0, nstat - 1
         res@tiMainString = st_code(istat) + " (" + \
             station_array&station@lat(istat) + \
             ", " + station_array&station@lon(istat) + \
             ") Offset: " + decimalPlaces(offsets(istat), 1, True)
         ;res@xyMonoDashPattern = False
             res@vpWidthF          = 0.155       ; Make plots quite a bit smaller.
             res@vpHeightF         = 0.15

         res@trXMinF = IYEAR
         res@trXMaxF = FYEAR+1
         res@trYMinF = min(stobs_ts(istat, :, :)) - 2
         res@trYMaxF = max(stobs_ts(istat, :, :)) + 2
         res@xyDashPattern = 0;(/0, 2, 0, 0, 0/)                         ; Make curves all solid
         res@xyLineColors = (/"blue", "black", "red","orange"/) ; 3 different colors
         res@tmXBFormat = "f" ; Don't put unnecessary 0's
         res@pmLegendDisplayMode = "Never" ;"Always"
         ;stobs_ts(istat, :, :) = (/stobs/)
         ;res@xyExplicitLegendLabels = stobs_ts&dataset

         res@tiYAxisString =  "XCO~B~2~N~ [ppmv]"
         res@xyLineThicknessF       = 2.0      ; default is 1
         res@xyMonoDashPattern      = True     ; force all solid lines
         res@tiMainFontHeightF = 0.01
         res@tiYAxisFontHeightF = 0.01
         res@tmXBLabelFontHeightF = 0.007;0.025
         res@tmYLLabelFontHeightF = 0.007;0.025
         res@tmYLLabelDeltaF = -0.5 ;offset to regular position, default 0
         res@tmXBLabelDeltaF = -0.5
         res@tiDeltaF = 0.9 ;Default 1.5
                                ;---Make legend
         ;res@pmLegendOrthogonalPosF = -1.05
         ;res@pmLegendParallelPosF   = 0.25
         ;res@pmLegendDisplayMode    = "Always"
         ;res@pmLegendWidthF         = 0.05/5.
         ;res@pmLegendHeightF        = 0.2/5.
         ;res@lgLabelFontHeightF     = .015/5.
         sta_ts_plots(istat) = gsn_csm_xy(wks, stobs_ts&time, stobs_ts(istat, :, :), res)

         mres = True
         mres@gsMarkerIndex = 5
         mres@gsMarkerSizeF    = 0.01
         mres@gsMarkerColor    = "black"
         mres@gsMarkerThicknessF = 5.
         plot@$unique_string("marker")$ = gsn_add_polymarker(wks, plot,\
               st_lon(istat), st_lat(istat), mres)
     end do

     txres = True
     txres@txFontHeightF = 0.01

     ; Manual text placement!
     plot@$unique_string("test")$ = gsn_add_text(wks, plot, \
               stobs_ts&station(0) + " (" + st_alt(0) + " m)", \
               st_lon(0) + 20, st_lat(0) - 7, txres)

     plot@$unique_string("test")$ = gsn_add_text(wks, plot, \
               stobs_ts&station(1) + " (" + st_alt(1) + " m)", \
               st_lon(1) + 10, st_lat(1) - 7, txres)

     plot@$unique_string("test")$ = gsn_add_text(wks, plot, \
               stobs_ts&station(2) + " (" + st_alt(2) + " m)", \
               st_lon(2) - 20, st_lat(2)-7, txres)

     plot@$unique_string("test")$ = gsn_add_text(wks, plot, \
               stobs_ts&station(3) + " (" + st_alt(3) + " m)", \
               st_lon(3), st_lat(3) - 7, txres)

     plot@$unique_string("test")$ = gsn_add_text(wks, plot, \
               stobs_ts&station(4) + " (" + st_alt(4) + " m)", \
               st_lon(4) + 5, st_lat(4) + 7, txres)

     amres                  = True
     amres@amJust           = "BottomLeft"  ; Corner of plot for positioning

     ; ASK
     amres@amOrthogonalPosF =  0.43   ;  0.5 is the bottom edge of the plot.
                                   ; Need to leave room for tickmarks.
     amres@amParallelPosF   = -0.10   ; -0.5 is the left edge of the plot.
     anno1 = gsn_add_annotation(plot, sta_ts_plots(0), amres)  ; Attach first plot

     ; CGO
     amres@amOrthogonalPosF =  0.45   ;  0.5 is the bottom edge of the plot.
     amres@amParallelPosF  = 0.18
     anno2 = gsn_add_annotation(plot, sta_ts_plots(1), amres)  ; Attach first plot

     ; LEF
     amres@amOrthogonalPosF =  0.30   ;  0.5 is the bottom edge of the plot.
     amres@amParallelPosF  =  -0.445
     anno3 = gsn_add_annotation(plot, sta_ts_plots(2), amres)  ; Attach first plot

     ; WIS
     amres@amOrthogonalPosF =  -0.05   ;  0.5 is the bottom edge of the plot.
     amres@amParallelPosF  =  0.21
     anno3 = gsn_add_annotation(plot, sta_ts_plots(3), amres)  ; Attach first plot

     ; HUN
     amres@amOrthogonalPosF =  -0.06   ;  0.5 is the bottom edge of the plot.
     amres@amParallelPosF  =  -0.19
     anno3 = gsn_add_annotation(plot, sta_ts_plots(4), amres)  ; Attach first plot


     ;; Make shared legend

  gres = True
  gres@YPosPercent = 97.    ; expressed as %, 0->100, sets position of top border of legend
                            ;  when gres@Position is set to its default setting of "Top" (Default = 95.)
  gres@XPosPercent = 1.      ; expressed as %, 0->100, sets position of left border of legend(Default = 5.)
  gres@LineLabelWhiteSpacePercent = 1.

  lineres = True
  lineres@lgLineColors = (/"blue", "black", "red","orange"/) ; line colors
  lineres@lgLineThicknesses = 2.5                        ; line thicknesses
  lineres@LineLengthPercent = 4.                         ; expressed as %, 0->100, length of line

  textres = True
  textres@lgLabelFontHeights = 0.01
  textres@lgLabels = stobs_ts&dataset  ; legend labels (required)

  plot = simple_legend(wks,plot,gres,lineres,textres)


;---Drawing the "base" plot will draw all plots.
     draw(plot)
     frame(wks)
     ;print(stobs_ts)


  leave_msg(DIAG_SCRIPT, "")
end
