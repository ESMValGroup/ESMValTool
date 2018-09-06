; #############################################################################
; seaice/ts_line.ncl
; #############################################################################
; TIMESERIES PLOTS FOR SEA ICE AREA AND EXTENT
; Author: Daniel Senftleben (DLR, Germany)
; ESMVal project
; #############################################################################
;
; Description:
;     Creates line plot for sea ice area and extent versus time.
;
; Required diag_script_info attributes:
;     multi_model_mean: additionally plot MMM and standard deviation.
;     month: selected month or annual mean.
;     region: region to be plotted ( "Arctic" or "Antarctic").
;
; Optional diag_script_info attributes:
;     fill_pole_hole: fill observational hole at North pole, Default: False.
;     styleset: specify line and color style.
;     xy_line_legend: False: do not draw a legend inside the plot.
;     legend_outside: True: draw legend in an extra plot.
;     wdiag: if the diagnostics will be used for the MDER a diagnostic name
;            name is needed.
;     wdiag_title: this name will apper in the figure title.
;
; Caveats:
;     Modifies input_file_info@dataset.
;
; Modification history:
;     20180725-A_schl_ma: ported to v2.0.
;     20180321-A_laue_ax: bugfix: handle case of wdiag not defined;
;                         reintroduced call to write_references
;                         (bugs introduced with MDER changes).
;     20170512-A_senf_da: Added reporting.
;     20160906-A_senf_da: Added option to fill polar observation hole.
;     20151027-A_laue_ax: Moved call to 'write_references' to the beginning
;                         of the code.
;     20150325-A_laue_ax: Modified reference tags used for acknowledgements
;                         (projects, observations, etc.).
;     20140312-A_righ_ma: Adapted to new time_operations structure.
;     20140129-A_senf_da: Modularized to fit to new structure.
;     20130603-A_gott_kl: Implemented into ./diag_scripts.
;     20130405-A_gott_kl: Written based on previous code by Melanie Braeu.
;
; #############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/aux_plotting.ncl"
load "./diag_scripts/shared/plot/xy_line.ncl"

; load "./diag_scripts/shared/meta_data.ncl"

; load "./diag_scripts/calc_seasonal.ncl"

load "./diag_scripts/shared/mder.ncl"


begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Variable
  var0 = variables(0)
  field_type0 = field_types(0)

  ; Write references (TODO)
  ; write_references(diag_script, \
  ;                  "A_gott_kl", \
  ;                  (/"A_hueb_mi", "A_eval_ma", "A_brae_me", "A_senf_da"/), \
  ;                  "D_0001", \
  ;                  "", \
  ;                  (/"P_embrace", "P_climval"/))

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

  ; Check field type
  plot_fields = (/"T2Ms"/)
  ind_f = any(field_type0 .eq. plot_fields)
  if (all(ismissing(ind_f))) then
    error_msg("f", diag_script, "", "cannot use field " + field_type0)
  end if
  delete([/plot_fields, ind_f/])

  ; Required attributes
  req_atts = (/"region"/)
  exit_if_missing_atts(diag_script_info, req_atts)

  ; Check region selection
  if (diag_script_info@region .ne. "Arctic" .and. \
      diag_script_info@region .ne. "Antarctic") then
    error_msg("f", diag_script, "", "no (valid) region selected")
  end if

  ; Check fill_pole_hole option
  if (.not. isatt(diag_script_info, "fill_pole_hole")) then
    diag_script_info@fill_pole_hole = False
  else
    if (diag_script_info@fill_pole_hole .and. \
        diag_script_info@region .eq. "Antarctic") then
      diag_script_info@fill_pole_hole = False
      error_msg("w", diag_script, "", "No pole hole at Antarctic " + \
                "- setting diag_script_info@fill_pole_hole to False")
    end if
  end if

  ; Attributes for MDER analysis
  if (.not.isatt(diag_script_info, "wdiag")) then
    diag_script_info@wdiag = ""
  end if
  if (.not.isatt(diag_script_info, "wdiag_title")) then
    diag_script_info@wdiag_title = (/"", ""/)
  end if

  ; Get datasets (consider ensemble averaging if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens) then
    idx_datasets = UNIQ(input_file_info@dataset)
  else
    idx_datasets = ispan(0, dimsizes(input_file_info@dataset) - 1, 1)
  end if
  dataset_names = input_file_info@dataset(idx_datasets)

  ; Determine time range
  start_year = min(stringtoint(input_file_info@start_year))
  end_year = max(stringtoint(input_file_info@end_year))
  all_years = ispan(start_year, end_year, 1)

  ; For annotation
  smonths = (/"January", "February", "March", "April", \
              "May", "June", "July", "August", \
              "September", "October", "November", "December"/)
  i2months = (/"01", "02", "03", "04", "05", "06",\
               "07", "08", "09", "10", "11", "12"/)

  ; Output file name specifications
  ; (passed to get_wks for final output file name)
  if (diag_script_info@month .eq. "A") then
      i2month = "annual"
  else
      i2month = smonths(toint(diag_script_info@month) - 1)
  end if
  outfile_ext = "extent_" + var0 + "_" + \
    diag_script_info@region + "_" + i2month + "_" + \
    sprinti("%0.4i", start_year) + "-" + \
    sprinti("%0.4i", end_year)
  outfile_area = "area_" + var0 + "_" + \
    diag_script_info@region + "_" + i2month + "_" + \
    sprinti("%0.4i", start_year) + "-" + \
    sprinti("%0.4i", end_year)

  ; Get unique list of all models
  annots = project_style(diag_script_info, "annots")

  ; Array sizes
  DIM_DAT = dimsizes(idx_datasets)
  DIM_ALL_DAT = dimsizes(annots)

  ; Array for collecting sea-ice extent for all datasets
  val_ext = new((/DIM_ALL_DAT, end_year - start_year + 1/), float)
  val_ext!1 = "time"
  val_ext&time = int2dble(all_years)

  ; Array for collecting sea-ice area for all datasets
  val_area = new((/DIM_ALL_DAT, end_year - start_year + 1/), float)
  val_area!1 = "time"
  val_area&time = int2dble(all_years)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Prepare calculation of standard deviation of datasets
  if (diag_script_info@multi_model_mean .eq. "y") then
    val_ext_stddev = new((/4, end_year - start_year + 1/), float)
    val_area_stddev = new((/4, end_year - start_year + 1/), float)

  ; Initialize anyway: both are needed as parameters for xy_line
  else
    val_ext_stddev = 0
    val_area_stddev = 0
  end if

  ; Iterate over all datasets
  do i = 0, DIM_ALL_DAT - 1
    log_debug("Processing " + annots(i))

    ; Read data
    sic = read_data(i, var0, field_type0)

    ; Check units (Subsequent code expects %/100)
    if (sic@units.eq."%") then
      sic = sic / 100.0
      sic@units = "%/100"
    end if

    ; Further processing depends on whether sic is on a lat-lon grid
    if(all(iscoord(sic, (/"lat", "lon"/)))) then
      log_debug("sic is on a lat - lon grid")
      area_temp = map_area(sic&lat, sic&lon)

      ; Make 2D to be compatible with irregular grid handling below
      area_lat = conform_dims(dimsizes(area_temp), sic&lat, 0)
      area_lon = conform_dims(dimsizes(area_temp), sic&lon, 1)
      if (isatt(sic, "_FillValue")) then
        if (isatt(area_temp, "_FillValue")) then
          area_temp = where(area_temp .eq. area_temp@_FillValue, \
                            sic@_FillValue, area_temp)
        end if
        area_temp@_FillValue = sic@_FillValue
      end if
    else  ; irregular grid
      log_debug("sic is not on a lat-lon grid")

      ; Get area data & lat-lon from another file
      aux_vars = read_fx_data(i, "areacello")
      area_temp = aux_vars[0]
      if (ListCount(aux_vars) .gt. 1) then
        area_lat = aux_vars[1]
        area_lon = aux_vars[2]
      else
        area_lat = conform_dims(dimsizes(area_temp), area_temp&lat, 0)
        area_lon = conform_dims(dimsizes(area_temp), area_temp&lon, 1)
      end if
      delete(aux_vars)
    end if

    ; Time axis
    sic_temp = time_operations(sic, \
                               toint(input_file_info@start_year(i)),\
                               toint(input_file_info@end_year(i)), \
                               "extract", "", 0)
    sic_temp!0 = "time"
    delete(sic)
    sic = sic_temp
    copy_VarMeta(sic_temp, sic)
    delete(sic_temp)

    ; Fill pole hole
    if (diag_script_info@fill_pole_hole) then
      latstep = 0.1
      max_lat = 90.0
      lat1D = ndtooned(area_lat)
      do itime = 0, dimsizes(sic&time) - 1  ; Hole depends on time
        actual_lat = max_lat
        data1D = ndtooned(sic(itime, :, :))
        if (all(ismissing(data1D))) then    ; Time step missing
          continue
        end if
        do while(True)                      ; Will break when hole is filled
          i1D = ind(lat1D .ge. actual_lat .and. \
                    lat1D .lt. actual_lat + latstep)
          if (.not. all(ismissing(i1D))) then
            i2D = ind_resolve(ind(lat1D.ge.actual_lat), \
                              dimsizes(area_lat))
            dim_i2D = dimsizes(i2D)
            npts = dim_i2D(0)               ; Number of elements
            ndim = dim_i2D(1)               ; Rank of "lat"
            if (any(ismissing(data1D(i1D)))) then
              do n = 0, npts - 1
                ring = sic(itime, i2D(n, 0), i2D(n, 1))
                sic(itime, i2D(n, 0), i2D(n, 1)) = where(ismissing(ring), \
                                                         1.0, ring)
                delete(ring)
              end do
            else                          ; Out of the hole!
              delete([/i1D, i2D, dim_i2D, npts, ndim/])
              break
            end if
            delete([/i2D, dim_i2D, npts, ndim/])
          end if
          delete([/i1D/])
          actual_lat = actual_lat - latstep
        end do                              ; while
        delete(data1D)
      end do                                ; itime
      delete(lat1D)
    end if

    ; Select region -> name corresponds to the pole region, but entire
    ; hemisphere will be evaluated and converted to million sqkm
    if (diag_script_info@region .eq. "Antarctic") then
      area_temp = where(area_lat .lt. 0.0, area_temp / 1.0e12, 0.0)
    else if (diag_script_info@region .eq. "Arctic")
      area_temp = where(area_lat .gt. 0.0, area_temp / 1.0e12, 0.0)
    end if
    end if

    ; Calculate sea ice area
    sic_area = sic * conform(sic, area_temp, (/1, 2/))
    copy_VarMeta(sic, sic_area)

    ; Calculate sea ice extent
    extent_temp = where(sic .ge. 0.15, 1.0, 0.0)
    sic_ext = conform(sic, area_temp, (/1, 2/)) * extent_temp
    copy_VarMeta(sic, sic_ext)
    delete([/extent_temp, sic/])

    ; Calculate weighted annual mean or select month
    if (diag_script_info@month .eq. "A") then     ; weighted annual mean
      sic_ann_ext  = time_operations( \
        sic_ext, \
        toint(input_file_info@start_year(i)), \
        toint(input_file_info@end_year(i)), \
        "average", "yearly", True)
      sic_ann_area = time_operations( \
        sic_area, \
        toint(input_file_info@start_year(i)), \
        toint(input_file_info@end_year(i)), \
        "average", "yearly", True)
    else                                          ; select month
      month = toint(diag_script_info@month)
      date_ext = cd_calendar(sic_ext&time, 0)
      sic_ann_ext = sic_ext(ind(date_ext(:, 1).eq.month), :, :)
      date_area = cd_calendar(sic_area&time, 0)
      sic_ann_area = sic_area(ind(date_area(:, 1).eq.month), :, :)
      delete([/date_ext, date_area, month/])
    end if
    delete([/sic_ext, sic_area/])

    ; Sum up
    temp = dim_sum_Wrap(sic_ann_ext)
    delete(sic_ann_ext)
    sic_sum_ext = dim_sum_Wrap(temp)
    delete(temp)
    temp = dim_sum_Wrap(sic_ann_area)
    delete(sic_ann_area)
    sic_sum_area = dim_sum_Wrap(temp)
    delete(temp)

    ; Put into global array
    ia = toint(input_file_info@start_year(i)) - start_year
    ie = toint(input_file_info@end_year(i)) - start_year
    val_ext(i, ia:ie) = tofloat(sic_sum_ext)
    val_area(i, ia:ie) = tofloat(sic_sum_area)

    ; Clean up
    delete([/area_lon, area_lat, area_temp, sic_sum_ext, sic_sum_area, ia, \
           ie/])
  end do

  ; Average over ensemble members if desired
  if (avgens) then
    sic_ext = new((/DIM_DAT, end_year - start_year + 1/), typeof(val_ext))
    sic_area = new((/DIM_DAT, end_year - start_year + 1/), typeof(val_area))

    dims_ext = getVarDimNames(val_ext)
    dims_area = getVarDimNames(val_area)

    ; Iterate over datasets
    do i = 0, DIM_DAT - 1
      idat = idx_datasets(i)

      ; Check if there is more than one ensemble member
      ens = ind(input_file_info@dataset(idat) .eq. input_file_info@dataset)
      if (dimsizes(ens) .gt. 1) then
        tmpA0_ext = val_ext(ens, :)
        tmpA0_area = val_area(ens, :)
        sic_ext(i, :) = dim_avg_n_Wrap(tmpA0_ext, 0)
        sic_area(i, :) = dim_avg_n_Wrap(tmpA0_area, 0)
        delete([/tmpA0_ext, tmpA0_area/])
      else
        sic_ext(i, :) = val_ext(idat, :)
        sic_area(i, :) = val_area(idat, :)
      end if
      delete(ens)
    end do
    delete([/val_ext, val_area/])

    ; Array for collecting sea-ice extent for all datasets
    val_ext = sic_ext
    val_ext!1 = "time"
    val_ext&time = int2dble(all_years)

    ; Array for collecting sea-ice area for all datasets
    val_area = sic_area
    val_area!1 = "time"
    val_area&time = int2dble(all_years)
    delete([/sic_ext, sic_area, dims_ext, dims_area/])
  else
    if (dimsizes(UNIQ(input_file_info@dataset)) .eq. 1) then
      dataset_names := unique_labels_min((/"dataset", "ensemble"/))
    end if
  end if

  ; Calculate standard deviation of datasets
  if (diag_script_info@multi_model_mean .eq. "y") then
    ; See style_$project$.ncl for which data shall
    ; be included in the statistics
    avgstd = project_style(diag_script_info, "avgstd")
    temp = new((/DIM_DAT/), integer)
    do i = 0, DIM_DAT - 1
      idat = idx_datasets(i)
      temp(i) = avgstd(idat)
    end do
    do imonth = 0, dimsizes(val_ext&time) - 1
      ; For extent
      val_ext_stddev(0, imonth) = \
        dim_avg_n_Wrap(val_ext(temp, imonth), (/0/))
      val_ext_stddev(1, imonth) = \
        dim_stddev_n_Wrap(val_ext(temp, imonth), (/0/))
      val_ext_stddev(2, imonth) = \
        val_ext_stddev(0, imonth) - val_ext_stddev(1, imonth)
      val_ext_stddev(3, imonth) = \
        val_ext_stddev(0, imonth) + val_ext_stddev(1, imonth)

      ; For area
      val_area_stddev(0, imonth) = \
        dim_avg_n_Wrap(val_area(temp, imonth), (/0/))
      val_area_stddev(1, imonth) = \
        dim_stddev_n_Wrap(val_area(temp, imonth), (/0/))
      val_area_stddev(2, imonth) = \
        val_area_stddev(0, imonth) - val_area_stddev(1, imonth)
      val_area_stddev(3, imonth) = \
        val_area_stddev(0, imonth) + val_area_stddev(1, imonth)
    end do
    delete(temp)
  end if

  ; Some attributes for extent plot
  val_ext!1 = "time"
  val_ext&time = int2dble(all_years)
  if (diag_script_info@month .eq. "A") then   ; weighted annual mean
    val_ext&time@units = "Annual"
  else                                      ; select a month
    val_ext&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_ext@units = "10~S~6~N~ km~S~2~N~"
  val_ext@long_name = "Sea Ice Extent"

  ; Some attributes for area plot
  val_area!1 = "time"
  val_area&time = int2dble(all_years)
  if diag_script_info@month .eq. "A" then   ; weighted annual mean
    val_area&time@units = "Annual"
  else                                      ; select a month
    val_area&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_area@units = "10~S~6~N~ km~S~2~N~"
  val_area@long_name = "Sea Ice Area"

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------
  if (config_user_info@write_plots) then
    input_file_info@dataset := dataset_names   ; Used by unique_labels_min

    if (dimsizes(val_ext&time) .eq. 1) then
      log_debug("Skipping tsline: plot needs at least two time values, has 1")
    else
      ; Create plot variables
      wks_ext = get_wks(0, diag_script, outfile_ext)
      wks_area = get_wks(0, diag_script, outfile_area)

      ; Define plot resources
      res = True
      res@tmXBLabelAngleF = 45
      res@tmXBLabelJust = "CenterRight"

      ; Min and max values on x-axis, leave space on right for legend
      res@trXMaxF = max(val_ext&time) + 0.25 * (max(val_ext&time) - \
                                                min(val_ext&time))
      res@trXMinF = min(val_ext&time) - 0.05 * (max(val_ext&time) - \
                                                min(val_ext&time))
      res@tmXBMode = "Automatic"
      res@tmXBFormat = "f"
      if (isatt(diag_script_info, "xy_line_legend")) then
        if (.not.diag_script_info@xy_line_legend) then
          res@trXMaxF = max(val_ext&time)
          res@trXMinF = min(val_ext&time)
        end if
      end if

      ; Resources specific for extent plot
      res_ext = True
      copy_VarMeta(res, res_ext)
      res_ext@tiMainString = val_ext&time@units + " " + \
        diag_script_info@region + " " + val_ext@long_name
      res_ext@tiYAxisString = val_ext@units

      ; Resources specific for area plot
      res_area = True
      copy_VarMeta(res, res_area)
      res_area@tiMainString  = val_area&time@units + " " + \
        diag_script_info@region + " " + val_area@long_name
      res_area@tiYAxisString = val_area@units

      ; Draw the plots
      debuginfo = False
      xy_line(wks_ext, val_ext, val_ext&time, val_ext_stddev, res_ext, \
              debuginfo)
      if (isatt(diag_script_info, "legend_outside")) then
        ; Legend can only once be drawn in an extra file (hlu error),
        ; but is same for both area and extent anyway
        if (diag_script_info@legend_outside) then
          val_area@legend_outside = True
          wks_area@legendfile = get_outfile_name(outfile_area + "_legend")
        end if
      else
        diag_script_info@legend_outside = False
      end if
      xy_line(wks_area, val_area, val_area&time, val_area_stddev, res_area, \
              debuginfo)

      ; ; Add meta data (TODO)
      ; climofiles = new(DIM_DAT, string)
      ; do ii = 0, DIM_DAT - 1
      ;     climofiles(ii) = interface_get_inpaths(ii) + "/" + \
      ;                      interface_get_infile(var0, field_type0, ii)
      ; end do
      ; plot_fname = wks_ext@fullname
      ; alltags = array_append_record(tags, \
      ;                               (/"DM_polar", "PT_seas", "ST_mean", \
      ;                               "ST_stddev"/), 0)
      ; caption = "Timeseries of sea ice extent."
      ; if (diag_script_info@fill_pole_hole) then
      ;   caption = caption + " Polar hole in observations was filled " + \
      ;               "assuming 100% sea ice concentration."
      ; end if
      ; id = "#idSeaIceTslineExtent"
      ; contrib_authors = (/"A_gott_kl", "A_hueb_mi", "A_eval_ma", \
      ;                     "A_brae_me", "A_senf_da", "A_wenz_sa"/)
      ; ESMValMD(plot_fname, alltags, caption, id, var0, dataset_names, \
      ;          climofiles, diag_script, contrib_authors)
      ; delete([/plot_fname, alltags, caption, id/])
      ;
      ; plot_fname = wks_area@fullname
      ; alltags = array_append_record(tags, \
      ;                               (/"DM_polar", "PT_seas", "ST_mean",\
      ;                               "ST_stddev"/), 0)
      ; caption = "Timeseries of sea ice area."
      ; if (diag_script_info@fill_pole_hole) then
      ;   caption = caption + " Polar hole in observations was filled " + \
      ;     "assuming 100% sea ice concentration."
      ; end if
      ; id = "#idSeaIceTslineArea"
      ; ESMValMD(plot_fname, alltags, caption, id, var0, dataset_names, \
      ;          climofiles, diag_script, contrib_authors)
      ; delete([/plot_fname, alltags, caption, id/])
    end if  ; More than 1 time step
  end if    ; If write_plots

  ; ---------------------------------------------------------------------------
  ; MDER output
  ; ---------------------------------------------------------------------------

  ; Setup arrays
  val_trend = new((/2, DIM_DAT/), "float")
  val_trend!0 = "case"
  val_trend&case = (/"val", "stddev"/)
  val_trend!1 = "dataset"
  val_trend&dataset = dataset_names
  val_mean = new((/2, DIM_DAT/), "float")
  val_mean!0 = "case"
  val_mean&case = (/"val", "stddev"/)
  val_mean!1 = "dataset"
  val_mean&dataset = dataset_names

  ; Calculate trend and mean
  do idat = 0, DIM_DAT - 1
    rc = regline(val_ext&time, val_ext(idat, :))
    val_trend(0, idat) = tofloat(rc)
    val_trend(1, idat) = tofloat(rc@rstd)
    val_mean(0, idat) = tofloat(avg(val_ext(idat, :)))
    val_mean(1, idat) = tofloat(stddev(val_ext(idat, :)) / \
                                sqrt(dimsizes(val_ext(idat, :))))
  end do

  ; Basic plot settings
  if (config_user_info@write_plots) then
    delete(res)
    dum0 = new(DIM_DAT, graphic)
    dum1 = new(DIM_DAT, graphic)
    error_bar0 = new(DIM_DAT, graphic)
    error_bar1 = new(DIM_DAT, graphic)
    plot_num = integertochar(ispan(97, 122, 1))
    colors = project_style(diag_script_info, "colors")
    thicks = project_style(diag_script_info, "thicks")
    lg_labels = unique_labels_min((/"dataset", "ensemble"/))
    wks_ext1 = get_wks(0, diag_script, outfile_ext + "_trend")
    wks_ext2 = get_wks(0, diag_script, outfile_ext + "_climmean")
    ti = ispan(1, dimsizes(val_trend(0, :)), 1)

    ; Appearance
    res = True
    res@tiYAxisString = " ( " + val_ext@units + "/dec)"
    res@tiMainString = diag_script_info@wdiag_title(0)
    res@gsnFrame = False                        ; don't draw yet
    res@gsnDraw = False                         ; don't advance frame
    res@gsnMaximize = True
    res@xyMarker = thicks(0)
    res@xyMonoMarkerColor = False
    res@xyMarkerColors = colors(0)              ; Marker color
    res@xyMarkerSizeF = 0.01                    ; Marker size (default 0.01)
    res@tiMainFontHeightF = 0.022
    res@trYMinF = min(val_trend(0, :)) - 0.5 * abs(min(val_trend(0, :)))
    res@trYMaxF = max(val_trend(0, :)) + 0.5 * abs(max(val_trend(0, :)))
    res@trXMinF = 0                             ; set X-axis min
    res@trXMaxF = DIM_DAT+1
    res@vpWidthF = 0.9
    res@vpHeightF = 0.4
    res@gsnYRefLine = 0.0
    res@tmXBLabelAngleF = 90.0
    res@tmXBMode = "Explicit"                   ; explicitly set Y-axis labels
    res@tmXBLabelJust = "CenterRight"
    res@tmXBValues = ti
    res@tmXBLabels = lg_labels
    res@tmXBLabelFontHeightF = 0.015

    ; Plot trend
    plot1 = gsn_csm_xy(wks_ext1, \
                       (/ti(0), ti(0)/), \
                       (/val_trend(0, 0), val_trend(0, 0)/), \
                       res)
    mkres = True
    mkres@gsMarkerSizeF = 0.01
    mkres@gsMarkerThicknessF = 1.0
    do t = 0, DIM_DAT - 1
      mkres@gsMarkerColor = colors(t)
      mkres@gsMarkerIndex = 16
      dum0(t) = gsn_add_polymarker(wks_ext1, plot1, ti(t), val_trend(0, t), \
                                   mkres)
    end do

    ; Plot error bars
    polyres   = True
    do t = 0, DIM_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar0(t) = gsn_add_polyline(wks_ext1, \
                                       plot1, \
                                       (/ti(t), ti(t)/), \
                                       (/val_trend(0, t) - val_trend(1, t), \
                                         val_trend(0, t) + val_trend(1, t)/), \
                                       polyres)
    end do
    draw(plot1)
    frame(wks_ext1)

    ; ; Add meta data (TODO)
    ; climofiles = new(DIM_DAT, string)
    ; do ii = 0, DIM_DAT - 1
    ;     climofiles(ii) = interface_get_inpaths(ii) + "/" + \
    ;                      interface_get_infile(var0, field_type0, ii)
    ; end do
    ; id = diag_script + "_" + var0
    ; contrib_authors = (/"A_zech_sa", "A_schl_ma"/)
    ; alldatasetnames = dataset_names
    ; alltags = array_append_record(tags, \
    ;                               (/"DM_sh", "ST_mean", "ST_clim", \
    ;                               "ST_range", "PT_errorbar", \
    ;                               "PT_scatter"/), 0)
    ; caption = "Trends in " + " " + var0 + \
    ;           ", similar to Wenzel et al., JClim, 2016, SI"
    ; plot_fname = wks_ext1@fullname
    ; ESMValMD(plot_fname, alltags, caption, id, var0, \
    ;          alldatasetnames, climofiles, diag_script, contrib_authors)
    ; delete([/alltags, caption, plot_fname/])

    ; Plot mean
    res@tiMainString = diag_script_info@wdiag_title(1)
    res@tiYAxisString = " (" + val_ext@units + ")"
    res@trYMinF = min(val_mean(0, :)) - 0.2 * abs(min(val_mean(0, :)))
    res@trYMaxF = max(val_mean(0, :)) + 0.2 * abs(max(val_mean(0, :)))
    plot2 = gsn_csm_xy(wks_ext2, \
                       (/ti(0), ti(0)/), \
                       (/val_mean(0, 0), val_mean(0, 0)/), \
                       res)
    mkres = True
    mkres@gsMarkerSizeF = 0.01
    mkres@gsMarkerThicknessF = 1.0
    do t = 0, DIM_DAT - 1
      mkres@gsMarkerColor = colors(t)
    mkres@gsMarkerIndex = 16
      dum1(t) = gsn_add_polymarker(wks_ext2, plot2, ti(t), val_mean(0, t), \
                                   mkres)
    end do

    ; Plot error bars
    polyres = True
    do t = 0, DIM_DAT - 1
      polyres@gsLineColor = colors(t)
      error_bar1(t) = gsn_add_polyline(wks_ext2, \
                                       plot2, \
                                       (/ti(t), ti(t)/), \
                                       (/val_mean(0, t) - val_mean(1, t), \
                                         val_mean(0, t) + val_mean(1, t)/), \
                                       polyres)
    end do
    draw(plot2)
    frame(wks_ext2)

    ; ; Add meta data (TODO)
    ; alltags = array_append_record(tags, \
    ;                               (/"DM_sh", "ST_mean", "ST_clim", \
    ;                               "ST_range", "PT_errorbar", \
    ;                               "PT_scatter"/), 0)
    ; caption = "Climatological mean in " + var0 + " anomalies" + \
    ;           ", similar to Wenzel et al., JClim, 2016, SI"
    ; plot_fname = wks_ext2@fullname
    ; ESMValMD(plot_fname, alltags, caption, id, var0, \
    ;          alldatasetnames, climofiles, diag_script, contrib_authors)
    ; delete([/alltags, caption, id, climofiles, plot_fname/])
  end if

  ; ---------------------------------------------------------------------------
  ; Write NETCDF file for MDER
  ; ---------------------------------------------------------------------------

  output = mder_output_desired()
  if (output) then
    work_dir = output@work_dir
    opt = diag_script_info

    ; Iterate over all desired diagnostics
    do dd = 0, dimsizes(diag_script_info@wdiag) - 1
      opt@wdiag := diag_script_info@wdiag(dd)
      log_info(opt@wdiag)

      ; Iterate over all datasets
      do i = 0, DIM_DAT - 1
        idat = idx_datasets(i)
        log_debug("Processing " + dataset_names(i))
        new_path = work_dir + diag_file_name(idat, opt)

        ; Write files
        if (isStrSubset(opt@wdiag, "_t")) then
          val = val_trend(case|:, dataset|i)
          val@units = "(" + variable_info@units + ")"
        else if(isStrSubset(opt@wdiag, "_c")) then
          val = val_mean(case|:, dataset|i)
          val@units = "(" + variable_info@units + "/yr)"
        else
          error_msg("f", diag_script, "", "Diagnostic " + opt@wdiag + \
                    "not supported for MDER output")
        end if
        end if

        val@var = var0
        val@diag_script = (/diag_script/)
        copy_VarAtts(diag_script_info, val)
        val@wdiag := opt@wdiag
        val@dataset = dataset_names(i)
        val@start_year = input_file_info@start_year(idat)
        val@end_year = input_file_info@end_year(idat)
        val@ncdf = new_path
        ncdf_outfile = ncdf_write(val, new_path)
        delete([/val, new_path/])
      end do
    end do
  end if

  leave_msg(diag_script, "")

end
