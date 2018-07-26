; #############################################################################
; mder/asr.ncl
; #############################################################################
; Author: Sabrina Zechlau (DLR, Germany)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;     Calculates the absorbed shortwave correlation gradient and saves output
;     for MDER analysis.
;
; Required recipe attributes:
;     season : specify season of year.
;
; Optional recipe attributes:
;     wdiag: if the diagnostics will be used for the MDER a diagnostic name
;         name is needed.
;     wdiag_title: this name will apper in the figure title.
;
; Caveats:
;     Modifies dataset_info@dataset.
;
; Modification history
;   * 20180725_A-schl_ma: ported to v2.0
;   * 201303??_A-wenz_sa: written.
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"

load "./diag_scripts/shared/plot/style.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

load "./diag_scripts/shared/mder.ncl"


begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Variable
  var0 = variables(0)
  field_type0 = field_types(0)

end

begin

  ; Plot file type
  plot_dir = config_user_info@plot_dir
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Season
  if (isatt(diag_script_info, "season")) then
    season = diag_script_info@season
  else
    if (isStrSubset(field_type0, "z") .eq. False) then
      error_msg("f", diag_script, "", "plot requires attribute " + \
                "'season' in recipe")
    end if
  end if

  ; Get datasets (consider ensemble averaging if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens .eq. "True") then
    avgens = True
    idx_datasets = UNIQ(dataset_info@dataset)
  else
    avgens = False
    idx_datasets = ispan(0, dimsizes(dataset_info@dataset) - 1, 1)
  end if

  ; Array sizes
  DIM_DAT = dimsizes(idx_datasets)

  ; Output arrays
  val_grad = new((/2, DIM_DAT/), "float")
  val_grad!0 = "case"
  val_grad&case = (/"val", "stddev"/)
  val_grad!1 = "dataset"
  val_grad&dataset = dataset_info@dataset(idx_datasets)

  ; Get lat range from preprocessed data
  lat = read_data(0, "lat", field_type0)
  if (dimsizes(lat) .gt. 1) then
    lat_range = (/min(lat), max(lat)/)
  else
    lat_range = read_data(0, "lat_bnds", field_type0)
  end if
  log_info("Retrieved latitude range " + lat_range(0) + " to " + \
           lat_range(1) + " from preprocessed data")

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Iterate over datasets
  do i = 0, DIM_DAT-1
    idat = idx_datasets(i)
    log_info(dataset_info@dataset(idat))

    ; Read data
    A0 = read_data(idat, var0, field_type0)

    ; Average over ensemble members if desired
    if (avgens) then
      ens = ind(dataset_info@dataset(idat) .eq. dataset_info@dataset)
      rank = dimsizes(dimsizes(A0))

      ; Check if there is more than one ensemble member
      if (dimsizes(ens) .gt. 1) then
        dims = getVarDimNames(A0)
        if (rank .eq. 3) then
          tmpA0 = new((/dimsizes(ens), dimsizes(A0&$dims(0)$), \
                      dimsizes(A0&$dims(1)$), dimsizes(A0&$dims(2)$)/), \
                      typeof(A0))
          copy_VarCoords(A0, tmpA0(0, :, :, :))
          tmpA0(0, :, :, :) = (/A0/)
        else
          error_msg("f", diag_script, "", "unsupported rank (" + rank + \
                    ") of variable " + var0)
        end if
        delete(A0)

        ; Average over ensemble members
        do iens = 1, dimsizes(ens) - 1
          tmp_data = read_data(idat + iens, var0, field_type0)
            if (rank .eq. 3) then
              tmpA0(iens, :, :, :) = (/tmp_data/)
            end if
          delete(tmp_data)
        end do
        A0 = dim_avg_n_Wrap(tmpA0, 0)
        delete(tmpA0)
      end if
      delete(ens)
    end if

    ; -------------------------------------------------------------------------
    ; Process temporal data
    ; -------------------------------------------------------------------------

    dummy0 = time_operations(A0, \
                             toint(input_file_info@start_year(idat)), \
                             toint(input_file_info@end_year(idat)), \
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

    val_grad(0, i) = tofloat(avg(grad))
    val_grad(1, i) = tofloat(stddev(grad))

    delete([/dummy0, grad, A0/])
  end do
  dataset_names = dataset_info@dataset(idx_datasets)

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Basic plotting settings
  if (config_user_info@write_plots) then
    dataset_info@dataset := dataset_names   ; Used by unique_labels_min
    colors = project_style(diag_script_info, "colors")
    thicks = project_style(diag_script_info, "thicks")
    lg_labels = unique_labels_min((/"dataset", "ensemble"/))
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
    res@tiYAxisString = "(index)"             ; add units title
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
    ; alldatasetnames = dataset_names
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

  if (config_user_info@write_netcdf) then
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)
    opt = True
    opt@wdiag = diag_script_info@wdiag
    log_info(opt@wdiag)

    ; Iterate over all datasets
    do i = 0, DIM_DAT - 1
      idat = idx_datasets(i)
      log_debug("Processing " + dataset_names(i))
      new_path = work_dir + diag_file_name(idat, opt)

      ; Write files
      val = val_grad(case|:, dataset|i)
      val@var = var0
      val@units = "(index)"
      val@diag_script = (/diag_script/)
      copy_VarAtts(diag_script_info, val)
      val@dataset = dataset_names(i)
      val@start_year = input_file_info@start_year(idat)
      val@end_year = input_file_info@end_year(idat)
      val@ncdf = new_path
      ncdf_outfile = ncdf_write(val, new_path)
      delete([/val, new_path/])
    end do
  end if

  leave_msg(diag_script, "")

end
