; #############################################################################
; austral_jet/asr.ncl
; #############################################################################
; Author: Sabrina Zechlau (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     Calculates the absorbed shortwave correlation gradient and saves output
;     for MDER analysis.
;
; Required diag_script_info attributes:
;     season: specify season of year.
;
; Optional diag_script_info attributes:
;     wdiag: if the diagnostics will be used for the MDER a diagnostic name
;            name is needed.
;     wdiag_title: this name will apper in the figure title.
;
; Caveats:
;     Modifies input_file_info@dataset.
;
; Modification history
;     20180725_A-schl_ma: ported to v2.0
;     201303??_A-wenz_sa: written.
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
  if (DIM_VAR .gt. 1) then
    error_msg("w", DIAG_SCRIPT, "", "this diagnostic supports only one " + \
              "variable, processing " + VAR0@short_name)
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

  ; Output arrays
  val_grad = new((/2, DIM_DAT/), "float")
  val_grad!0 = "case"
  val_grad&case = (/"val", "stddev"/)
  val_grad!1 = "dataset"
  val_grad&dataset = DATASETS

  ; Get lat range from preprocessed data
  f = addfile(INFO0[0]@filename, "r")
  lat = f->lat
  if (dimsizes(lat) .gt. 1) then
    lat_range = (/min(lat), max(lat)/)
  else
    lat_range = f->lat_bnds
  end if
  log_info("Retrieved latitude range " + lat_range(0) + " to " + \
           lat_range(1) + " from preprocessed data")
  delete(f)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Iterate over datasets
  do idat = 0, DIM_DAT-1
    log_info(DATASETS(idat))

    ; Average over ensemble members if desired
    if (avgens) then
      atts = True
      atts@dataset = DATASETS(idat)
      info := select_metadata_by_atts(INFO0, atts)
    else
      info = NewList("lifo")
      ListAppend(info, INFO0[idat])
    end if

    ; Get data
    A0 = get_average(info)
    info := info[0]

    ; -------------------------------------------------------------------------
    ; Process temporal data
    ; -------------------------------------------------------------------------

    dummy0 = time_operations(A0, info@start_year, info@end_year, \
                             "average", season, True)

    ; -------------------------------------------------------------------------
    ; Process spatial data
    ; -------------------------------------------------------------------------

    grad = new(dimsizes(dummy0&month), float)
    do tt = 0, dimsizes(dummy0&month) - 1
      dummy2 = area_operations(dummy0(tt, :, :), -20.0, -50.0, 0.0, 360.0, \
                               "average", True)
      dummy3 = area_operations(dummy0(tt, :, :), -50.0, -90.0, 0.0, 360.0, \
                               "average", True)
      grad(tt) = dummy2 - dummy3
      delete([/dummy2, dummy3/])
    end do

    val_grad(0, idat) = tofloat(avg(grad))
    val_grad(1, idat) = tofloat(stddev(grad))

    delete([/dummy0, grad, A0/])
  end do

  ; Consider ensemble averaging
  if (avgens) then
    INFO0 := get_unique_items(INFO0, "dataset")
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Basic plotting settings
  if (config_user_info@write_plots) then
    colors = project_style(INFO0, diag_script_info, "colors")
    thicks = project_style(INFO0, diag_script_info, "thicks")
    lg_labels = unique_labels_min(INFO0, (/"dataset", "ensemble"/))
    plot_file = config_user_info@plot_dir + var0 + "_" + season

    ; Create directory
    system("mkdir -p " + config_user_info@plot_dir)

    ; Plotting instances
    wks = gsn_open_wks(file_type, plot_file)
    dum0 = new(DIM_DAT, graphic)
    error_bar0 = new(DIM_DAT, graphic)
    ti = ispan(1, dimsizes(val_grad(0, :)), 1)

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

    ; Plot appearance
    res = True
    res@tiYAxisString = "W m-2"               ; add units title
    res@tiMainString = season + " ASR " + lat_string + " gradient (" + \
      diag_script_info@wdiag_title(0) + ")"
    res@gsnFrame = False                      ; don't draw yet
    res@gsnDraw  = False                      ; don't advance frame
    res@gsnMaximize = True
    res@xyMarker = thicks(0)                  ; choose type of marker
    res@xyMonoMarkerColor = False
    res@xyMarkerColors = colors(0)            ; Marker color
    res@xyMarkerSizeF = 0.01                  ; Marker size (default 0.01)
    res@tiMainFontHeightF = 0.022
    res@trYMinF = min(val_grad(0, :)) - 0.5 * abs(min(val_grad(0, :)))
    res@trYMaxF = max(val_grad(0, :)) + 0.5 * abs(max(val_grad(0, :)))
    res@trXMinF = 0
    res@trXMaxF = DIM_DAT + 1
    res@vpWidthF  = 0.9                       ; Change the aspect ratio, but
    res@vpHeightF = 0.4                       ; make plot as large as possible.
    res@gsnYRefLine = 0.0
    res@tmXBLabelAngleF = 90.0
    res@tmXBMode = "Explicit"                 ; explicitly set Y-axis labels
    res@tmXBLabelJust = "CenterRight"
    res@tmXBValues = ti
    res@tmXBLabels = lg_labels
    res@tmXBLabelFontHeightF = 0.015

    ; Plot data
    plot = gsn_csm_xy(wks, \
                      (/ti(0), ti(0)/), \
                      (/val_grad(0, 0), val_grad(0, 0)/), \
                      res)
    mkres = True
    mkres@gsMarkerSizeF = 0.01
    mkres@gsMarkerThicknessF = 1.0
    do t = 0, DIM_DAT - 1
      mkres@gsMarkerColor = colors(t)
      mkres@gsMarkerIndex = 16
      dum0(t) = gsn_add_polymarker(wks, plot(0), ti(t), val_grad(0, t), mkres)
    end do

    ; Plot error bars
    polyres = True
    do t = 0, DIM_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar0(t) = gsn_add_polyline(wks, plot, (/ti(t), ti(t)/), \
                                       (/val_grad(0, t) - val_grad(1, t), \
                                       val_grad(0, t) + val_grad(1, t)/), \
                                       polyres)
    end do
    draw(plot)
    frame(wks)

    ; ; Add meta data (TODO)
    ; climofiles = new(DIM_DAT, string)
    ; do ii = 0, DIM_DAT - 1
    ;   climofiles(ii) = interface_get_inpaths(ii) + "/" + \
    ;                    interface_get_infile(var0, field_type0, ii)
    ; end do
    ; id = diag_script + "_" + var0
    ; contrib_authors = (/"A_zech_sa", "A_schl_ma"/)
    ; alldatasetnames = DATASETS
    ; alltags = array_append_record(tags,
    ;                               (/"DM_sh", "ST_mean", "ST_clim", \
    ;                               "ST_range", "PT_errorbar", "PT_scatter", \
    ;                               "T_clouds"/), 0)
    ; caption = season + " ASR gradient for " + lat_string + \
    ;           ", similar to Wenzel et al., JClim, 2016, SI"
    ; ESMValMD(plot_file + "." + file_type, alltags, caption, id, var0, \
    ;          alldatasetnames, climofiles, diag_script, contrib_authors)
    ; delete([/alltags, caption, id, climofiles, alldatasetnames/])
  end if

  ; ---------------------------------------------------------------------------
  ; Write NETCDF file for MDER
  ; ---------------------------------------------------------------------------

  output = mder_output_desired()
  if (output) then
    work_dir = output@work_dir
    opt = diag_script_info
    log_info(opt@wdiag)

    ; Iterate over all datasets
    do idat = 0, DIM_DAT - 1
      log_debug("Processing " + DATASETS(idat))
      info = INFO0[idat]
      new_path = work_dir + diag_file_name(info@filename, opt)

      ; Write files
      val = val_grad(case|:, dataset|idat)
      val@var = var0
      val@units = "W m-2"
      val@diag_script = (/DIAG_SCRIPT/)
      copy_VarAtts(diag_script_info, val)
      val@dataset = DATASETS(idat)
      val@start_year = info@start_year
      val@end_year = info@end_year
      val@ncdf = new_path
      ncdf_outfile = ncdf_write(val, new_path)
      delete([/val, new_path/])
    end do
  end if

  leave_msg(DIAG_SCRIPT, "")

end
