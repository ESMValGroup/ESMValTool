; #############################################################################
; austral_jet_diags.ncl
; #############################################################################
; Author: Sabrina Zechlau (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     This script is based on the diagnostic eyring13jclim.ncl but
;     calculates trends and annual mean for various present day diagnostiocs
;     such as near-global (82.5°S-82.5°N) ozone at 50 hPa and temperature at
;     100 hPa, September-October-November-December (SOND) ozone at 50 hPa
;     over Antarctica (60°S-90°S), October-November-December-January (ONDJ)
;     temperature at 100 hPa over Antarctica (60°S-90°S), DJF SH jet
;     position at 850 hPa, and DJF upper tropospheric tropical (30°S-30°N)
;     temperatures at 250 hPa. The trends are calculated over 1979-2005 for
;     the past and over 2006-2050 for the future.
;     The output will optional be saved in NetCDF-files for further
;     processing of the MDER diagnostics.
;
; Required recipe attributes:
;     styleset: style set ("DEFAULT", "CMIP5").
;     e13fig10_season: seasons (n-element array of strings, "ANN", "JAN",
;         "FEB", "MAR", "DJF", "SON", etc.).
;
;  Optional recipe attributes:
;      rsondes: a set of additional observations can be added to the plot
;          but will not be saved for the MDER diagnostics.
;      wdiag: if the diagnostics will be used for the MDER a diagnostic name
;          name is needed.
;      wdiag_title: this name will apper in the figure title.
;
; Caveats:
;
; Modification history:
;   * 20180716-A_schl_ma: ported to v2.0
;   * 20180329-A-wenz_sa: modified original script (eyring13jclim.ncl) to
;                         austral-jet_diags.ncl
;   * 20170816-A_bock_ls: added tags for reporting
;   * 20151209-A-laue_ax: code clean-up, added call to write_references,
;                         several bugfixes for processing zonal wind (ua)
;                         and handling missing variable attribute "units"
;   * 201510??-A-cion_ir: written
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"

load "./diag_scripts/shared/plot/style.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

load "./diag_scripts/shared/mder.ncl"


begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Write references (TODO)
  ; write_references(diag_script,       \
  ;                  "A_zech_sa",       \
  ;                  "A_schl_ma",       \
  ;                  "D_wenzel16jclim", \
  ;                  "",                \
  ;                  "P_embrace")

end

begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Season
  if(isatt(diag_script_info, "e13fig10_season")) then
    season = diag_script_info@e13fig10_season
  else
    if (isStrSubset(field_type0, "z") .eq. False) then
      error_msg("f", diag_script, "", "plot requires attribute " + \
                "'e13fig10_season' in recipe")
    end if
  end if

  ; Additional observations
  if (isatt(diag_script_info, "rsondes")) then
    if(isatt(diag_script_info, "e13fig10_yr_min")) then
      year_min = diag_script_info@e13fig10_yr_min
    else
      error_msg("f", diag_script, "", "plot requires attribute " + \
                " 'e13fig10_yr_min' in recipe")
    end if
    if(isatt(diag_script_info, "e13fig10_yr_max")) then
      year_max = diag_script_info@e13fig10_yr_max
    else
      error_msg("f", diag_script, "", "plot requires attribute " + \
                " 'e13fig10_yr_max' in recipe")
    end if
  end if

  ; Attributes for MDER analysis
  if (.not. isatt(diag_script_info, "wdiag")) then
    diag_script_info@wdiag = ""
  end if
  if (.not. isatt(diag_script_info, "wdiag_title")) then
    diag_script_info@wdiag_title = (/"", ""/)
  end if

  ; Get datasets (consider ensemble averaging if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens .eq. "True") then
    idx_datasets = UNIQ(dataset_info@dataset)
  else
    idx_datasets = ispan(0, dimsizes(dataset_info@dataset) - 1, 1)
  end if

  ; Variable
  var0 = variables(0)
  field_type0 = field_types(0)

  ; Array sizes
  dim_DAT = dimsizes(idx_datasets)

  ; Trend array
  val_trend = new((/2, dim_DAT/), "float")
  val_trend!0 = "case"
  val_trend&case = (/"val", "stddev"/)
  val_trend!1 = "dataset"
  val_trend&dataset = dataset_info@dataset(idx_datasets)

  ; Mean array
  val_mean = new((/2, dim_DAT/), "float")
  val_mean!0 = "case"
  val_mean&case = (/"val", "stddev"/)
  val_mean!1 = "dataset"
  val_mean&dataset = dataset_info@dataset(idx_datasets)

  ; Get plev and lat range from preprocessed data
  plev = read_data(0, "plev", field_type0)
  lat = read_data(0, "lat", field_type0)
  if (dimsizes(plev) .eq. 1) then
    lev = plev(0)
  else
    error_msg("f", diag_script, "", "ambiguous pressure levels given, " + \
              "only one supported")
  end if
  if (dimsizes(lat) .gt. 1) then
    lat_range = (/min(lat), max(lat)/)
  else
    lat_range = read_data(0, "lat_bnds", field_type0)
  end if

  ; ---------------------------------------------------------------------------
  ; Read and precprocess data
  ; ---------------------------------------------------------------------------

  ; Iterate over datasets
  do i = 0, dim_DAT - 1
    idat = idx_datasets(i)
    log_info(dataset_info@dataset(idat))

    ; Read data
    A0 = read_data(idat, var0, field_type0)

    ; Average over ensemble members if desired
    if (avgens .eq. True) then
      ens = ind(dataset_info@dataset(idat) .eq. dataset_info@dataset)
      rank = dimsizes(dimsizes(A0))

      ; Check if there is more than one ensemble member
      if (dimsizes(ens) .gt. 1) then
        dims = getVarDimNames(A0)
        if (rank .eq. 4) then
          tmpA0 = new((/dimsizes(ens), dimsizes(A0&$dims(0)$), \
                      dimsizes(A0&$dims(1)$), dimsizes(A0&$dims(2)$), \
                      dimsizes(A0&$dims(3)$)/), typeof(A0))
          copy_VarCoords(A0, tmpA0(0, :, :, :, :))
          tmpA0(0, :, :, :, :) = (/A0/)
        else if (rank .eq. 3) then
          tmpA0 = new((/dimsizes(ens), dimsizes(A0&$dims(0)$), \
                      dimsizes(A0&$dims(1)$), dimsizes(A0&$dims(2)$)/), \
                      typeof(A0))
          copy_VarCoords(A0, tmpA0(0, :, :, :))
          tmpA0(0, :, :, :) = (/A0/)
        else if (rank .eq. 2) then
          tmpA0 = new((/dimsizes(ens), dimsizes(A0&$dims(0)$), \
                      dimsizes(A0&$dims(1)$)/), typeof(A0))
          copy_VarCoords(A0, tmpA0(0, :, :))
          tmpA0(0, :, :) = (/A0/)
        else if (rank .eq. 1) then
          tmpA0 = new((/dimsizes(ens), dimsizes(A0&$dims(0)$)/), typeof(A0))
          copy_VarCoords(A0, tmpA0(0, :))
          tmpA0(0, :) = (/A0/)
        else
          error_msg("f", diag_script, "", "unsupported rank (" + rank + \
                    ") of variable " + var0)
        end if
        end if
        end if
        end if
        delete(A0)

        ; Average over ensemble members
        do iens = 1, dimsizes(ens) - 1
          tmp_data = read_data(idat + iens, var0, field_type0)
          if (rank .eq. 4) then
            tmpA0(iens, :, :, :, :) = (/tmp_data/)
          else if (rank .eq. 3) then
            tmpA0(iens, :, :, :) = (/tmp_data/)
          else if (rank .eq. 2) then
            tmpA0(iens, :, :) = (/tmp_data/)
          else if (rank .eq. 1) then
            tmpA0(iens, :) = (/tmp_data/)
          end if
          end if
          end if
          end if
          delete(tmp_data)
        end do
        A0 = dim_avg_n_Wrap(tmpA0, 0)
        delete(tmpA0)
      end if
      delete(ens)
    end if

    ; -------------------------------------------------------------------------
    ; Process spatial data
    ; -------------------------------------------------------------------------

    ; Process ua (jet position)
    if (var0 .eq. "ua") then
      tmpx   = dim_avg_Wrap(A0)
      dummy1 = new(dimsizes(tmpx&time), typeof(A0&lat))

      ; 2nd degree polynomial fit at maximum for all times to get jet position
      do tt = 0, dimsizes(tmpx&time) - 1
        idx = maxind(tmpx(tt, :))
        scal_x = tmpx(tt, idx-1:idx+1)
        scal_y = tmpx&lat(idx-1:idx+1)
        pf = lspoly(scal_x, scal_y, 1.0, 3)
        yfit = pf(0) + pf(1) * scal_x + pf(2) * scal_x ^ 2
        dummy1(tt) = yfit(1)
        delete([/pf, idx, yfit, scal_x, scal_y/])
      end do
      dummy1!0 = "time"
      dummy1&time = A0&time
      variable_info@units = "deg"
      delete(tmpx)

    ; Process mmstf (SH Hadley cell boundary)
    else if (var0 .eq. "mmstf") then
      tmp = A0
      dummy1 = new(dimsizes(A0&time), typeof(A0&lat))

      ; 1st degree polynomial fit for all times
      do tt = 0, dimsizes(A0&time) - 1
        idxmin = minind(tmp(tt, :))
        idxmax = maxind(tmp(tt, :))
        tmpp = tmp(tt, idxmin:idxmax)
        idx1 = min(ind(tmpp .ge. 0.0))
        idx2 = max(ind(tmpp .le. 0.0))
        if (idx2 .gt. (idx1 - 1)) then
          idx2 = idx1 - 1
        end if
        scal_x = tmpp(idx2:idx1)
        scal_y = tmpp&lat(idx2:idx1)
        pf = lspoly(scal_x, scal_y, 1, 2)
        yfit = pf(0) + pf(1) * scal_x
        dummy1(tt) = yfit(1)
        delete([/pf, idx1, idx2, idxmin, idxmax, yfit, tmpp, scal_x, scal_y/])
      end do
      dummy1!0 = "time"
      dummy1&time = A0&time
      variable_info@units = "deg"
      delete(tmp)

    ; Other variables (no operations needed, done in preprocessing)
    else
        dummy1 = A0
    end if
    end if
    delete(A0)

    ; Convert units of tro3
    if (var0 .eq. "tro3") then
      dummy1 = dummy1 * 10 ^ -3
      variable_info@units = "ppmv"
    end if

    ; -------------------------------------------------------------------------
    ; Process temporal data
    ; -------------------------------------------------------------------------

    if (season .ne. "yearly") then
      tmp = time_operations(dummy1, \
                            toint(input_file_info@start_year(idat)), \
                            toint(input_file_info@end_year(idat)), \
                            "extract", season, True)
      dummy2 = time_operations(tmp, \
                               toint(input_file_info@start_year(idat)), \
                               toint(input_file_info@end_year(idat)), \
                               "average", "yearly", True)
      delete(tmp)
    else
      dummy2 = time_operations(dummy1, \
                               toint(input_file_info@start_year(idat)), \
                               toint(input_file_info@end_year(idat)), \
                               "average", season, True)
    end if

    ; Trend
    ttmp = ispan(toint(input_file_info@start_year(idat)), \
                 toint(input_file_info@end_year(idat)), 1)
    dummy3 = dummy2 - dim_avg_n(dummy2, 0)
    rc = regline(ttmp, dummy3)
    if (var0 .ne. "tpp") then
      val_trend(0, i) = tofloat(rc) * 10
      val_trend(1, i) = tofloat(rc@rstd) * 10
    else
      val_trend(0, i) = tofloat(rc)
      val_trend(1, i) = tofloat(rc@rstd)
    end if

    ; Mean
    val_mean(0, i) = tofloat(avg(dummy2))
    val_mean(1, i) = tofloat(stddev(dummy2) / sqrt(dimsizes(dummy2)))
    delete([/dummy1, dummy2, dummy3/])
    delete([/ttmp, rc/])
  end do

  ; ---------------------------------------------------------------------------
  ; Process additional observations
  ; ---------------------------------------------------------------------------

  if (isatt(diag_script_info, "rsondes")) then
    rsondes = diag_script_info@rsondes
    rsondes_f = diag_script_info@rsondes_file
    dim_rsondes = dimsizes(rsondes)
    tmp_vtrend = val_trend
    tmp_vt = val_mean
    delete([/val_trend, val_mean/])
    dim_DAT = dim_DAT + dim_rsondes
    dataset_names = new(dim_DAT, string)
    dataset_names(:dim_rsondes-1) = rsondes
    dataset_names(dim_rsondes:) = dataset_info@dataset(idx_dataset)

    ; Trend
    val_trend   = new((/2, dim_DAT/), "float")
    val_trend(:, dim_rsondes:) = (/tmp_vtrend/)
    val_trend!0 = "case"
    val_trend&case = (/"val", "stddev"/)
    val_trend!1 = "dataset"
    val_trend&dataset = dataset_names

    ; Mean
    val_mean = new((/2, dim_DAT/), "float")
    val_mean(:, dim_rsondes:) = (/tmp_vt/)
    val_mean!0 = "case"
    val_mean&case = (/"val", "stddev"/)
    val_mean!1 = "dataset"
    val_mean&dataset = dataset_names
    delete([/tmp_vtrend, tmp_vt/])

    ; Read files
    do irs = 0, dim_rsondes - 1
      rs_fils = systemfunc("ls " + rsondes_f(irs))
      tmp = addfile(rs_fils, "r")
      A0  = tmp->$var0$
      delete([/rs_fils, tmp/])
      tt = cd_calendar(A0&time, 0)
      rst_max = max(tt(:, 0))
      rst_min = min(tt(:, 0))
      if (year_min .ge. rst_min) then
        t_min = year_min
      else
        t_min = rst_min
      end if
      if (year_max .ge. rst_max) then
        t_max = year_max
      else
        t_max = rst_max
      end if

      ; Read data and process data
      ttmp = ispan(toint(t_min), toint(t_max), 1)
      A0&time@calendar = "standard"
      A0!1 = "plev"
      A0!2 = "station"
      dummy1 = tofloat(dim_avg_Wrap(A0(time|:, {plev|lev}, station|:)))
      dummy1!0 = "time"
      dummy1&time = A0&time
      if (season.ne."yearly") then
          tmp_t = time_operations(dummy1, toint(t_min), toint(t_max), \
                                  "extract", season, True)
          dummy2 = time_operations(tmp_t, toint(t_min), toint(t_max), \
                                   "average", "yearly", True)
          delete(tmp_t)
      else
          dummy2 = time_operations(dummy1, toint(t_min), toint(t_max), \
                                   "average", season, True)
      end if

      ; Trend
      dummy3 = dummy2 - dim_avg_n(dummy2, 0)
      rc = regline(ttmp, dummy3)
      val_trend(0, irs) = tofloat(rc) * 10
      val_trend(1, irs) = tofloat(rc@rstd) * 10

      ; Mean
      val_mean(0, irs) = tofloat(avg(dummy2))
      val_mean(1, irs) = tofloat(stddev(dummy2))
      delete([/dummy1, dummy2, dummy3/])
      delete([/ttmp, rc, A0/])
    end do
  else
    dataset_names = dataset_info@dataset(idx_datasets)
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Basic plotting settings
  if (config_user_info@write_plots) then
    delete(dataset_info@dataset)
    dataset_info@dataset = dataset_names
    colors = project_style(diag_script_info, "colors")
    thicks = project_style(diag_script_info, "thicks")
    lg_labels = unique_labels_min((/"dataset", "ensemble"/))
    plot_file1 = config_user_info@plot_dir + diag_script + "_" + var0 + \
      "_" + diag_script_info@e13fig10_season(0) + "_trend"
    plot_file2 = config_user_info@plot_dir + diag_script + "_" + var0 + \
      "_" + diag_script_info@e13fig10_season(0) + "_climmean"

    ; Create directory
    system("mkdir -p " + config_user_info@plot_dir)

    ; Plotting instances
    wks1 = gsn_open_wks(file_type, plot_file1)
    wks2 = gsn_open_wks(file_type, plot_file2)
    dum0 = new(dim_DAT, graphic)
    dum1 = new(dim_DAT, graphic)
    error_bar0 = new(dim_DAT, graphic)
    error_bar1 = new(dim_DAT, graphic)
    ti = ispan(1, dimsizes(val_trend(0, :)), 1)

    ; Get plot descriptions
    if (lat_range(0).eq.90.0 .and. lat_range(1).eq.-90.0 .or. \
        lat_range(0).eq.82.5 .and. lat_range(1).eq.-82.5) then
      lat_string = "NG"
    else
      if (lat_range(0) .ge. 0.0) then
        lat_string = "(" + round(lat_range(0), 3) + "N"
      else if (lat_range(0) .lt. 0.0) then
        lat_string = "(" + -1 * round(lat_range(0), 3) + "S"
      end if
      end if
      if (lat_range(1) .ge. 0.0) then
        lat_string = lat_string + round(lat_range(1), 3) + "N)"
      else if (lat_range(1) .lt. 0.0) then
        lat_string = lat_string + round(lat_range(1), 3) + "S)"
      end if
      end if
    end if
    if (season .eq. "yearly") then
      season = "ANN"
    end if

    ; -------------------------------------------------------------------------
    ; Plot 1: Trend
    ; -------------------------------------------------------------------------

    res = True
    if ((isatt(variable_info, "units"))) then
        res@tiYAxisString = "(" + variable_info@units + "/dec)"
    else
        res@tiYAxisString = "(per dec)"
    end if
    res@tiMainString = lev/100 + " hPa " + season + " " + var0 + " " + \
      lat_string + " trend (" + diag_script_info@wdiag_title(0) + ")"
    res@gsnFrame = False                    ; don't draw yet
    res@gsnDraw = False                     ; don't advance frame
    res@gsnMaximize = True
    res@xyMarker = thicks(0)                ; choose type of marker (16)
    res@xyMonoMarkerColor = False
    res@xyMarkerColors = colors(0)          ; Marker color
    res@xyMarkerSizeF = 0.01                ; Marker size (default 0.01)
    res@tiMainFontHeightF = 0.022
    res@trYMinF = min(val_trend(0, :)) - 0.5 * abs(min(val_trend(0, :)))
    res@trYMaxF = max(val_trend(0, :)) + 0.5 * abs(max(val_trend(0, :)))
    res@trXMinF = 0
    res@trXMaxF = dim_DAT + 1
    res@vpWidthF = 0.9                      ; Change the aspect ratio, but
    res@vpHeightF = 0.4                     ; make plot as large as possible.
    res@gsnYRefLine = 0.0
    res@tmXBLabelAngleF = 90.0
    res@tmXBMode = "Explicit"               ; explicitly set Y-axis labels
    res@tmXBLabelJust = "CenterRight"
    res@tmXBValues = ti
    res@tmXBLabels = lg_labels
    res@tmXBLabelFontHeightF = 0.015

    ; Plot trend
    plot1 = gsn_csm_xy(wks1, \
                       (/ti(0), ti(0)/), \
                       (/val_trend(0, 0), val_trend(0, 0)/), \
                       res)
    mkres = True
    mkres@gsMarkerSizeF = 0.01
    mkres@gsMarkerThicknessF = 1.0
    do t = 0, dim_DAT - 1
      mkres@gsMarkerColor = colors(t)
      mkres@gsMarkerIndex = 16
      dum0(t) = gsn_add_polymarker(wks1, plot1, ti(t), val_trend(0, t), mkres)
    end do

    ; Plot error bars
    polyres = True
    do t = 0, dim_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar0(t) = gsn_add_polyline(wks1, plot1, (/ti(t), ti(t)/), \
                                       (/val_trend(0, t) - val_trend(1, t), \
                                       val_trend(0, t) + val_trend(1, t)/), \
                                       polyres)
    end do
    draw(plot1)
    frame(wks1)

    ; ; Add meta data (TODO)
    ; climofiles = new(dim_DAT, string)
    ; do ii = 0, dim_DAT - 1
    ;     climofiles(ii) = interface_get_inpaths(ii) + "/" + \
    ;                      interface_get_infile(var0, field_type0, ii)
    ; end do
    ; id = diag_script + "_" + var0
    ; contrib_authors = (/"A_zech_sa", "A_schl_ma"/)
    ; alldatasetnames = dataset_info@dataset
    ; alltags = array_append_record(tags, \
    ;                               (/"DM_sh", "ST_mean", "ST_clim", \
    ;                               "ST_range", "PT_errorbar", \
    ;                               "PT_scatter"/), 0)
    ; caption = "Trends in " + season + " " + var0 + " anomalies at " \
    ;           + lev/100 + " hPa for " + lat_string + \
    ;           ", similar to Wenzel et al., JClim, 2016, SI"
    ; ESMValMD(plot_file1 + "." + file_type, alltags, caption, id, var0, \
    ;          alldatasetnames, climofiles, diag_script, contrib_authors)
    ; delete([/alltags, caption/])

    ; -------------------------------------------------------------------------
    ; Plot 2: Climatological mean
    ; -------------------------------------------------------------------------

    res@tiMainString = lev/100 + " hPa " + season + " " + var0 + " " + \
      lat_string + " mean (" + diag_script_info@wdiag_title(1) + ")"
    res@tiYAxisString = "(" + variable_info@units + ")"
    res@trYMinF = min(val_mean(0, :)) - 0.1 * abs(min(val_mean(0, :) - \
                                                      val_mean(1, :)))
    res@trYMaxF = max(val_mean(0, :)) + 0.1 * abs(max(val_mean(0, :) + \
                                                      val_mean(1, :)))

    ; Plot mean
    plot2 = gsn_csm_xy(wks2, \
                       (/ti(0), ti(0)/), \
                       (/val_mean(0, 0), val_mean(0, 0)/), \
                       res)
    mkres = True
    mkres@gsMarkerSizeF = 0.01
    mkres@gsMarkerThicknessF = 1.0
    do t = 0, dim_DAT - 1
      mkres@gsMarkerColor = colors(t)
      mkres@gsMarkerIndex = 16    ; thicks(t)
      dum1(t) = gsn_add_polymarker(wks2, plot2, ti(t), val_mean(0, t), mkres)
    end do

    ; Plot error bars
    polyres = True
    do t = 0, dim_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar1(t) = gsn_add_polyline(wks2, plot2, (/ti(t), ti(t)/), \
                                       (/val_mean(0, t) - val_mean(1, t), \
                                       val_mean(0, t) + val_mean(1, t)/), \
                                       polyres)
    end do
    draw(plot2)
    frame(wks2)

    ; ; Add meta data (TODO)
    ; alltags = array_append_record(tags, \
    ;                               (/"DM_sh", "ST_mean", "ST_clim", \
    ;                               "ST_range", "PT_errorbar", \
    ;                               "PT_scatter"/), 0)
    ; caption = "Climatological mean in " + season + " " + var0 + \
    ;           " anomalies at " + lev/100 + " hPa for " + lat_string + \
    ;           ", similar to Wenzel et al., JClim, 2016, SI"
    ; ESMValMD(plot_file2 + "." + file_type, alltags, caption, id, var0, \
    ;          alldatasetnames, climofiles, diag_script, contrib_authors)
    ; delete([/alltags, caption, id, climofiles, alldatasetnames/])
  end if

  ; ---------------------------------------------------------------------------
  ; Write NETCDF file
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_netcdf) then
    tmp_d = diag_script_info@wdiag
    do dd = 0, dimsizes(tmp_d) - 1
      delete(diag_script_info@wdiag)
      diag_script_info@wdiag = tmp_d(dd)

      ; Iterate over all datasets
      do i = 0, dim_DAT - 1
        idat = idx_datasets(i)
        log_debug("Processing " + dataset_info@dataset(i))
        new_path = diag_file_path(idat, diag_script, diag_script_info)

        ; Write files
        if (isStrSubset(diag_script_info@wdiag, "_t")) then
          val = val_trend(case|:, dataset|i)
          val@var = var0
          val@units = "(" + variable_info@units + "/dec)"
          val@diag_script = (/diag_script/)
          copy_VarAtts(diag_script_info, val)
          val@dataset = dataset_info@dataset(i)
          val@start_year = input_file_info@start_year(idat)
          val@end_year = input_file_info@end_year(idat)
          val@ncdf = new_path
          ncdf_outfile = ncdf_write(val, new_path)
          delete(val)
        else if(isStrSubset(diag_script_info@wdiag, "_c")) then
          val = val_mean(case|:, dataset|i)
          val@var = var0
          val@units = "(" + variable_info@units + "/dec)"
          val@diag_script = (/diag_script/)
          copy_VarAtts(diag_script_info, val)
          val@dataset = dataset_info@dataset(i)
          val@start_year = input_file_info@start_year(idat)
          val@end_year = input_file_info@end_year(idat)
          val@ncdf = new_path
          ncdf_outfile = ncdf_write(val, new_path)
          delete(val)
        end if
        end if
        delete(new_path)
      end do
    end do
  end if

  leave_msg(diag_script, "")

end
