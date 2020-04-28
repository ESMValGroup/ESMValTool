; #############################################################################
; austral_jet/main.ncl
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
; Required diag_script_info attributes:
;     styleset: style set ("DEFAULT", "CMIP5").
;     season: seasons (n-element array of strings, "ANN", "JAN", "FEB", "MAR",
;                      "DJF", "SON", etc.).
;
; Optional diag_script_info attributes:
;     average_ens: average over ensemble members of the same model (default:
;                  False).
;     rsondes: a set of additional observations can be added to the plot
;              but will not be saved for the MDER diagnostics.
;     rsondes_file: paths to additional observations, required when rsondes is
;                   set
;     rsondes_yr_min: minimum year for plot, required when rsondes is set.
;     rsondes_yr_max: maximum year for plot, required when rsondes is set.
;     wdiag: if the diagnostics will be used for the MDER a diagnostic name
;            name is needed.
;     wdiag_title: this name will apper in the figure title.
;     derive_var: derive variable (one of "tpp", "mmstf")
;     derive_latrange: latitude range for derived variable, required when
;                      derive_var is set.
;     derive_lev: level for derived variable, required when derive_var is set.
;
; Caveats:
;     "warning: in unique_labels_min (diag_scripts/shared/plot/style.ncl), Add
;     more attributes to prio to make labels unique! Continuing with non-unique
;     labels" is normal when using "average_ens = true".
;
; Modification history:
;     20191121-schlund_manuel: added new provenance tracking.
;     20180716-schlund_manuel: ported to v2.0.
;     20180329-wenzel_sabrina: modified original script (eyring13jclim.ncl) to
;                              austral-jet_diags.ncl.
;     20170816-bock_lisa: added tags for reporting.
;     20151209-lauer_axel: code clean-up, added call to write_references,
;                          several bugfixes for processing zonal wind (ua)
;                          and handling missing variable attribute "units".
;     201510??-cionni_irene: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/set_operators.ncl"

load "$diag_scripts/shared/plot/style.ncl"

load "$diag_scripts/shared/mder.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  ; Variable
  VAR0 = variable_info[0]
  var0 = VAR0@short_name
  DIM_VAR = ListCount(variable_info)
  if (DIM_VAR .gt. 2) then
    error_msg("w", DIAG_SCRIPT, "", "this diagnostic supports at most two " + \
              "variables, got " + DIM_VAR)
  end if

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  DATASETS = metadata_att_as_array(INFO0, "dataset")

  ; Ensemble averaging (if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens) then
    DATASETS := get_unique_values(DATASETS)
  else
    DATASETS := unique_labels_min(INFO0, (/"dataset", "ensemble"/))
  end if
  DIM_DAT = dimsizes(DATASETS)

  ; Directories
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  ; Provenance
  AUTHORS = (/"cionni_irene", "lauer_axel", "bock_lisa", "wenzel_sabrina", \
             "schlund_manuel"/)
  DOMAIN = "sh"
  PLOT_TYPE = "errorbar"
  REFERENCES = (/"wenzel16jclim"/)
  ANCESTORS = metadata_att_as_array(INFO0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

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

  ; Required attributes
  req_atts = (/"season"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  season = diag_script_info@season

  ; Additional observations
  if (isatt(diag_script_info, "rsondes")) then
    if(isatt(diag_script_info, "rsondes_yr_min")) then
      year_min = diag_script_info@rsondes_yr_min
    else
      error_msg("f", DIAG_SCRIPT, "", "plot requires attribute " + \
                "'rsondes_yr_min' in recipe if 'rsondes' is given")
    end if
    if(isatt(diag_script_info, "rsondes_yr_max")) then
      year_max = diag_script_info@rsondes_yr_max
    else
      error_msg("f", DIAG_SCRIPT, "", "plot requires attribute " + \
                "'rsondes_yr_max' in recipe if 'rsondes' is given")
    end if
  end if

  ; Attributes for MDER analysis
  if (.not. isatt(diag_script_info, "wdiag")) then
    diag_script_info@wdiag = (/"", ""/)
  end if
  if (.not. isatt(diag_script_info, "wdiag_title")) then
    diag_script_info@wdiag_title = (/"", ""/)
  end if

  ; Output arrays
  val_trend = new((/2, DIM_DAT/), "float")
  val_trend!0 = "case"
  val_trend&case = (/"val", "stddev"/)
  val_trend!1 = "dataset"
  val_trend&dataset = DATASETS
  val_mean = new((/2, DIM_DAT/), "float")
  val_mean!0 = "case"
  val_mean&case = (/"val", "stddev"/)
  val_mean!1 = "dataset"
  val_mean&dataset = DATASETS

  ; ---------------------------------------------------------------------------
  ; Derive variables (using NCL)
  ; ---------------------------------------------------------------------------

  if (isatt(diag_script_info, "derive_var")) then
    new_var = diag_script_info@derive_var
    var0 = new_var
    log_info("Deriving variable " + new_var)

    ; Check attributes
    if (isatt(diag_script_info, "derive_latrange")) then
      lat_range = diag_script_info@derive_latrange
    else
      error_msg("f", DIAG_SCRIPT, "", "when deriving a variable the " + \
                "attribute derive_latrange is required")
    end if
    if (isatt(diag_script_info, "derive_lev")) then
      lev = diag_script_info@derive_lev
    else
      error_msg("f", DIAG_SCRIPT, "", "when deriving a variable the " + \
                "attribute derive_lev is required")
    end if
    log_info("Retrieved pressure level " + lev + " and latitude range " + \
             lat_range(0) + " to " + lat_range(1) + " from recipe")

    ; Tropopause pressure tpp
    if (new_var .eq. "tpp") then
      INFO0 := select_metadata_by_name(input_file_info, "ta")
      VAR0@short_name = new_var
      VAR0@standard_name = "tropopause_air_pressure"
      VAR0@long_name = "Tropopause Air Pressure"
      VAR0@units = "hPa"

      ; Iterate over datasets
      dim_info = ListCount(INFO0)
      do idat = 0, dim_info - 1
        info = INFO0[idat]
        f = addfile(info@filename, "r")
        if isfilevar(f, new_var) then
          log_info("Variable " + new_var + " for " + info@dataset + \
                   " already derived, skipping")
          copy_VarAtts(VAR0, INFO0[idat])
          delete(f)
          continue
        end if

        ; Read data
        A0 = read_data(info)
        copy_VarAtts(VAR0, INFO0[idat])

        ; Order plev: lowest pressure first (in hPa)
        pr = A0&plev
        pr = pr / 100.0
        if (pr(0) .eq. max(pr)) then
          pin = pr(::-1)
          tin = A0({time|:}, {lat|:}, {lon|:}, {plev|::-1})
        else
          pin = pr
          tin = A0({time|:}, {lat|:}, {lon|:}, {plev|:})
        end if

        ; Calculate tropopause pressure and copy metedata
        outvar = trop_wmo(pin, tin, 0, False)
        outvar!0 = A0!0
        outvar!1 = A0!2
        outvar!2 = A0!3
        outvar&time = A0&time
        outvar&lat = A0&lat
        outvar&lon = A0&lon
        outvar@var = new_var
        outvar := area_operations(outvar, lat_range(0), lat_range(1), 0.0, \
                                  360.0, "average", True)
        outvar@diag_script = DIAG_SCRIPT
        ncdf_outfile = ncdf_write(outvar, info@filename)
        delete([/A0, pr, pin, tin, outvar, ncdf_outfile/])
      end do

    ; Zonal mean meridional stream function mmstf
    else if (new_var .eq. "mmstf") then
      INFO0 := select_metadata_by_name(input_file_info, "ps")
      INFO1 := select_metadata_by_name(input_file_info, "va")
      VAR0@short_name = new_var
      VAR0@standard_name = "zonal_meridional_stream_function"
      VAR0@long_name = "Zonal Meridional Stream Function"
      VAR0@units = "kg s-1"

      ; Iterate over datasets
      dim_info = ListCount(INFO0)
      do idat = 0, dim_info - 1
        info0 = INFO0[idat]
        info1 = INFO1[idat]

        ; Check datasets
        if (isatt(info0, "ensemble") .and. isatt(info1, "ensemble")) then
          fail = ((info0@dataset .ne. info1@dataset) .or. \
                  (info0@ensemble .ne. info1@ensemble))
        else
          fail = (info0@dataset .ne. info1@dataset)
        end if
        if (fail) then
          error_msg("f", DIAG_SCRIPT, "", "derivation of 'mmstf' not " + \
                    "possible, 'va' and 'ps' do not cover the same datasets")
        end if

        ; Check if variable was already derived
        f = addfile(info0@filename, "r")
        if isfilevar(f, new_var) then
          log_info("Variable " + new_var + " for " + info0@dataset + \
                   " already derived, skipping")
          copy_VarAtts(VAR0, INFO0[idat])
          delete(f)
          continue
        end if

        ; Read data
        ps = read_data(info0)
        va = read_data(info1)
        copy_VarAtts(VAR0, INFO0[idat])
        if (dimsizes(va&lat) .ne. dimsizes(ps&lat)) then
          tmp = ps({time|:}, {lat|va&lat}, {lon|:})
          delete(ps)
          ps = tmp
          delete(tmp)  ; ps := ps(...) causes memory leaks
        end if

        ; Preprocess plev (see docu of zonal_mpsi_Wrap)
        pr = va&plev
        pidx = ind(pr .gt. 500.0 .and. pr .lt. 100500.0)
        pr := pr(pidx)
        vt2b = va(:, pidx, :, :)
        if (pr(0) .eq. max(pr)) then
          pr := pr(::-1)
          tmp = vt2b({time|:}, {plev|::-1}, {lat|:}, {lon|:})
          delete(vt2b)
          vt2b = tmp
          delete(tmp)  ; vt2b := vt2b(...) causes memory leaks
        end if

        ; Calculate zonal meridional stream function
        outvar = zonal_mpsi_Wrap(vt2b, va&lat, pr, ps)
        outvar!0 = va!0
        outvar!1 = va!1
        outvar!2 = va!2
        outvar&time = va&time
        outvar&plev = pr
        outvar&lat = va&lat
        outvar@var = new_var
        outvar@diag_script = DIAG_SCRIPT
        ncdf_outfile = ncdf_write(outvar, info0@filename)
        delete([/va, ps, pr, pidx, vt2b, outvar, ncdf_outfile/])
      end do
    else
      error_msg("f", DIAG_SCRIPT, "", "don't know how to derive variable " + \
                new_var)
    end if
    end if

  ; Get plev and lat range from preprocessed data if variable is not derived
  else
    info = INFO0[0]
    f = addfile(info@filename, "r")
    if (isfilevar(f, "plev") .and. isfilevar(f, "lat")) then
      plev = f->plev
      lat = f->lat
      if (dimsizes(plev) .eq. 1) then
        lev = plev(0)
      else
        error_msg("f", DIAG_SCRIPT, "", "ambiguous pressure levels given, " + \
                  "only one supported")
      end if
      if (dimsizes(lat) .gt. 1) then
        lat_range = (/min(lat), max(lat)/)
      else
        lat_range = f->lat_bnds
      end if
    else
      if (isatt(f, "plev") .and. \
          isatt(f, "lat_range_0") .and. \
          isatt(f, "lat_range_1")) then
        lev = f@plev
        lat_range = (/f@lat_range_0, f@lat_range_1/)
      else
        error_msg("f", DIAG_SCRIPT, "", "cannot read plev and latrange")
      end if
    end if
    log_info("Retrieved pressure level " + lev + " and latitude range " + \
             lat_range(0) + " to " + lat_range(1) + " from preprocessed data")
    delete(f)
  end if

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Iterate over datasets
  do idat = 0, DIM_DAT - 1
    log_info(DATASETS(idat))

    ; Average over ensemble members if desired
    if (avgens) then
      atts = True
      atts@dataset = DATASETS(idat)
      info := select_metadata_by_atts(INFO0, atts)
    else
      info := NewList("lifo")
      ListAppend(info, INFO0[idat])
    end if

    ; Get data
    A0 = get_average(info)
    info := info[0]

    ; -------------------------------------------------------------------------
    ; Process spatial data
    ; -------------------------------------------------------------------------

    ; Process mmstf (SH Hadley cell boundary)
    if (var0 .eq. "mmstf") then
      tmp = A0({time|:}, {plev|lev}, {lat|lat_range(0):lat_range(1)})
      dummy1 = new(dimsizes(A0&time), typeof(A0&lat))

      ; 1st degree polynomial fit for all times
      do tt = 0, dimsizes(A0&time) - 1
        idxmin = minind(tmp(tt, :))
        idxmax = maxind(tmp(tt, :))
        tmpp = tmp(tt, idxmin:idxmax)
        idx1 = min(ind(tmpp .ge. 0.0))
        idx2 = max(ind(tmpp .le. 0.0))
        if (ismissing(idx1)) then
          idx1 = 1
        end if
        if (ismissing(idx2)) then
          idx2 = dimsizes(tmpp) - 1
        end if
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
      VAR0@units = "degrees"
      delete(tmp)

    ; Other variables (no operations needed, done in preprocessing)
    else
        dummy1 = A0
    end if
    delete(A0)

    ; Convert units of tro3
    if (var0 .eq. "tro3") then
      dummy1 = dummy1 * 10e-3
      VAR0@units = "ppmv"
    end if

    ; -------------------------------------------------------------------------
    ; Process temporal data
    ; -------------------------------------------------------------------------

    if (season .ne. "yearly") then
      tmp = time_operations(dummy1, info@start_year, info@end_year, \
                            "extract", season, True)
      dummy2 = time_operations(tmp, info@start_year, info@end_year, \
                               "average", "yearly", True)
      delete(tmp)
    else
      dummy2 = time_operations(dummy1, info@start_year, info@end_year, \
                               "average", season, True)
    end if

    ; Trend
    ttmp = ispan(info@start_year, info@end_year, 1)
    dummy3 = dummy2 - dim_avg_n(dummy2, 0)
    rc = regline(ttmp, dummy3)
    if (var0 .ne. "tpp") then
      val_trend(0, idat) = tofloat(rc) * 10
      val_trend(1, idat) = tofloat(rc@rstd) * 10
    else
      val_trend(0, idat) = tofloat(rc)
      val_trend(1, idat) = tofloat(rc@rstd)
    end if

    ; Mean
    val_mean(0, idat) = tofloat(avg(dummy2))
    val_mean(1, idat) = tofloat(stddev(dummy2) / sqrt(dimsizes(dummy2)))
    delete([/dummy1, dummy2, dummy3, ttmp, rc/])
  end do

  ; Consider ensemble averaging
  if (avgens) then
    INFO0 := get_unique_items(INFO0, "dataset")
  end if

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
    DIM_DAT = DIM_DAT + dim_rsondes
    dataset_names = new(DIM_DAT, string)
    dataset_names(:dim_rsondes-1) = rsondes
    dataset_names(dim_rsondes:) = DATASETS
    DATASETS := dataset_names
    delete(dataset_names)

    ; Trend
    val_trend = new((/2, DIM_DAT/), "float")
    val_trend(:, dim_rsondes:) = (/tmp_vtrend/)
    val_trend!0 = "case"
    val_trend&case = (/"val", "stddev"/)
    val_trend!1 = "dataset"
    val_trend&dataset = DATASETS

    ; Mean
    val_mean = new((/2, DIM_DAT/), "float")
    val_mean(:, dim_rsondes:) = (/tmp_vt/)
    val_mean!0 = "case"
    val_mean&case = (/"val", "stddev"/)
    val_mean!1 = "dataset"
    val_mean&dataset = DATASETS
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
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Basic plotting settings
  if (config_user_info@write_plots) then
    colors = project_style(INFO0, diag_script_info, "colors")
    thicks = project_style(INFO0, diag_script_info, "thicks")
    lg_labels = unique_labels_min(INFO0, (/"dataset", "ensemble"/))
    plot_file1 = config_user_info@plot_dir + DIAG_SCRIPT + "_" + var0 + \
      "_" + season + "_trend"
    plot_file2 = config_user_info@plot_dir + DIAG_SCRIPT + "_" + var0 + \
      "_" + season + "_climmean"

    ; Plotting instances
    wks1 = gsn_open_wks(file_type, plot_file1)
    wks2 = gsn_open_wks(file_type, plot_file2)
    dum0 = new(DIM_DAT, graphic)
    dum1 = new(DIM_DAT, graphic)
    error_bar0 = new(DIM_DAT, graphic)
    error_bar1 = new(DIM_DAT, graphic)
    ti = ispan(1, dimsizes(val_trend(0, :)), 1)

    ; Get plot descriptions
    if (lat_range(0).le.-82.5 .and. lat_range(1).ge.82.5) then
      lat_string = "NG"
    else
      if (lat_range(0) .ge. 0.0) then
        lat_string = "(" + round(lat_range(0), 3) + "N"
      else if (lat_range(0) .lt. 0.0) then
        lat_string = "(" + -1 * round(lat_range(0), 3) + "S"
      end if
      end if
      lat_string = lat_string + "-"
      if (lat_range(1) .ge. 0.0) then
        lat_string = lat_string + round(lat_range(1), 3) + "N)"
      else if (lat_range(1) .lt. 0.0) then
        lat_string = lat_string + -1 * round(lat_range(1), 3) + "S)"
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
    if ((isatt(VAR0, "units"))) then
        res@tiYAxisString = "(" + VAR0@units + "/dec)"
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
    res@trXMaxF = DIM_DAT + 1
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
    do t = 0, DIM_DAT - 1
      mkres@gsMarkerColor = colors(t)
      mkres@gsMarkerIndex = 16
      dum0(t) = gsn_add_polymarker(wks1, plot1, ti(t), val_trend(0, t), mkres)
    end do

    ; Plot error bars
    polyres = True
    do t = 0, DIM_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar0(t) = gsn_add_polyline(wks1, plot1, (/ti(t), ti(t)/), \
                                       (/val_trend(0, t) - val_trend(1, t), \
                                       val_trend(0, t) + val_trend(1, t)/), \
                                       polyres)
    end do
    draw(plot1)
    frame(wks1)

    ; Write ncdf file
    new_path = work_dir + "trend_" + season + "_" + var0 + ".nc"
    val_trend@var = var0
    val_trend@diag_script = DIAG_SCRIPT
    val_trend@ncdf = new_path
    outfile_trend = ncdf_write(val_trend, new_path)

    ; Provenance tracking
    plot_path = plot_file1 + "." + file_type
    caption = "Trends in " + season + " " + var0 + "."
    statistics := (/"trend", "stddev"/)
    log_provenance(outfile_trend, plot_path, caption, statistics, \
                   DOMAIN, PLOT_TYPE, AUTHORS, REFERENCES, ANCESTORS)

    ; -------------------------------------------------------------------------
    ; Plot 2: Climatological mean
    ; -------------------------------------------------------------------------

    res@tiMainString = lev/100 + " hPa " + season + " " + var0 + " " + \
      lat_string + " mean (" + diag_script_info@wdiag_title(1) + ")"
    res@tiYAxisString = "(" + VAR0@units + ")"
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
    do t = 0, DIM_DAT - 1
      mkres@gsMarkerColor = colors(t)
      mkres@gsMarkerIndex = 16    ; thicks(t)
      dum1(t) = gsn_add_polymarker(wks2, plot2, ti(t), val_mean(0, t), mkres)
    end do

    ; Plot error bars
    polyres = True
    do t = 0, DIM_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar1(t) = gsn_add_polyline(wks2, plot2, (/ti(t), ti(t)/), \
                                       (/val_mean(0, t) - val_mean(1, t), \
                                       val_mean(0, t) + val_mean(1, t)/), \
                                       polyres)
    end do
    draw(plot2)
    frame(wks2)

    ; Write ncdf file
    new_path = work_dir + "mean_" + season + "_" + var0 + ".nc"
    val_mean@var = var0
    val_mean@diag_script = DIAG_SCRIPT
    val_mean@ncdf = new_path
    outfile_mean = ncdf_write(val_mean, new_path)

    ; Provenance tracking
    plot_path = plot_file2 + "." + file_type
    caption = "Climatological mean in " + season + " " + var0 + "."
    statistics := (/"mean", "clim", "stddev"/)
    log_provenance(outfile_mean, plot_path, caption, statistics, \
                   DOMAIN, PLOT_TYPE, AUTHORS, REFERENCES, ANCESTORS)
  end if

  ; ---------------------------------------------------------------------------
  ; Write NETCDF file for MDER
  ; ---------------------------------------------------------------------------

  output = mder_output_desired()
  if (output) then
    work_dir := output@work_dir
    opt = diag_script_info
    plot_path = "missing"
    plot_type = ""

    ; Iterate over all desired diagnostics
    do dd = 0, dimsizes(diag_script_info@wdiag) - 1
      opt@wdiag := diag_script_info@wdiag(dd)
      log_info(opt@wdiag)

      ; Iterate over all datasets
      do idat = 0, DIM_DAT - 1
        log_debug("Processing " + DATASETS(idat))
        info = INFO0[idat]
        new_path = work_dir + diag_file_name(info@filename, opt)

        ; Write files
        if (isStrSubset(opt@wdiag, "_t")) then
          val = val_trend(case|:, dataset|idat)
          val@units = "(" + VAR0@units + "/dec)"
          statistics := (/"trend", "stddev"/)
          caption = "Trend in " + season + " " + var0 + " of " + \
                    DATASETS(idat) + "."
        else if(isStrSubset(opt@wdiag, "_c")) then
          val = val_mean(case|:, dataset|idat)
          val@units = "(" + VAR0@units + ")"
          statistics := (/"mean", "clim", "stddev"/)
          caption = "Climatological mean in " + season + " " + var0 + \
                    " of " + DATASETS(idat) + "."
        else
          error_msg("f", DIAG_SCRIPT, "", "Diagnostic " + opt@wdiag + \
                    "not supported for MDER output")
        end if
        end if
        val@var = var0
        val@diag_script = (/DIAG_SCRIPT/)
        copy_VarAtts(diag_script_info, val)
        val@wdiag := opt@wdiag
        val@dataset = DATASETS(idat)
        val@start_year = info@start_year
        val@end_year = info@end_year
        val@ncdf = new_path
        ncdf_outfile = ncdf_write(val, new_path)
        delete([/val, new_path/])

        ; Provenance tracking
        log_provenance(ncdf_outfile, plot_path, caption, statistics, \
                       DOMAIN, plot_type, AUTHORS, REFERENCES, \
                       info@filename)
      end do
    end do
  end if

  leave_msg(DIAG_SCRIPT, "")

end
