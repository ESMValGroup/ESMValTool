; #############################################################################
; seaice_tsline.ncl
; #############################################################################
; LINE PLOTS OF SEA ICE AREA AND EXTENT vs TIME
; Author: Daniel Senftleben (DLR, Germany)
; ESMVal project
; #############################################################################
;
; Description
;     Creates line plot for sea ice area and extent versus time.
;
; Required recipe attributes:
;     month: selected month or annual mean.
;     region: region to be plotted ("Arctic" or "Antarctic").
;
; Optional recipe attributes:
;     fill_pole_hole: fill observational hole at North pole, Default: False
;     styleset: specify line and color style
;     xy_line_legend: False: do not draw a legend inside the plot
;     legend_outside: True: draw legend in an extra plot
;
; Caveats:
;     CMOR units expected: areacello in m2, sic in %
;
; Modification history
;   * 20180723-A_schl_ma: Extended diagnostic to use it in MDER.
;   * 20180620-A_senf_da: Ported to ESMValTool v2.
;   * 20160906-A_senf_da: Added option to fill polar observation hole.
;   * 20151027-A_laue_ax: Moved call to 'write_references' to the beginning
;                         of the code.
;   * 20150325-A_laue_ax: Modified reference tags used for acknowledgements
;                         (projects, observations, etc).
;   * 20140312-A_righ_ma: Adapted to new time_operations structure.
;   * 20140129-A_senf_da: Modularized to fit to new structure.
;   * 20130419-A_gott_kl: Implemented into ./diag_scripts.
;   * 20130405-A_gott_kl: Implemented.
;
; #############################################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/aux_plotting.ncl"
load "./diag_scripts/shared/plot/xy_line.ncl"

load "./diag_scripts/shared/mder.ncl"

begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  DIM_DAT = dimsizes(dataset_info@dataset)
  DIM_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ; Write references (TODO)
  ; write_references(diag_script, \
  ;                  "A_gott_kl", \
  ;                  (/"A_hueb_mi", "A_eval_ma", "A_brae_me", "A_senf_da"/), \
  ;                  "D_0001", \
  ;                  "", \
  ;                  (/"P_embrace", "P_climval"/))

end

begin

  ; Plot file type
  file_type = config_user_info@output_file_type
  if(ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Check field type
  plot_fields = (/"T2Ms"/)
  ind_f = any(field_type0 .eq. plot_fields)
  if (all(ismissing(ind_f))) then
    error_msg("f", diag_script, "", "can't use field " + field_type0)
  end if
  delete([/plot_fields, ind_f/])

  ; Check region selection
  if (diag_script_info@region.ne."Arctic" .and. \
      diag_script_info@region.ne."Antarctic") then
    error_msg("f", diag_script, "", "no (valid) region selection")
  end if

  ; Check fill_pole_hole option
  if (.not.isatt(diag_script_info, "fill_pole_hole")) then
    diag_script_info@fill_pole_hole = False
  else
    if (diag_script_info@fill_pole_hole .and. \
        diag_script_info@region.eq."Antarctic") then
      diag_script_info@fill_pole_hole = False
      error_msg("w", diag_script, "", "No pole hole at Antarctic " + \
                "- setting diag_script_info@fill_pole_hole to False")
    end if
  end if

  ; Determine time range
  start_year = min(stringtoint(dataset_info@start_year))
  end_year = max(stringtoint(dataset_info@end_year))
  all_years = ispan(start_year, end_year, 1)

  ; Select attributes
  annots = project_style(diag_script_info, "annots")
  avgstd = project_style(diag_script_info, "avgstd")

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

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Array for collecting sea-ice extent for all datasets
  val_ext = new((/dimsizes(annots), end_year - start_year + 1/), float)
  val_ext!1 = "time"
  val_ext&time = int2dble(all_years)

  ; Array for collecting sea-ice area for all datasets
  val_area = new((/dimsizes(annots), end_year - start_year + 1/), float)
  val_area!1 = "time"
  val_area&time = int2dble(all_years)

  ; Prepare calculation of standard deviation of datasets
  if (diag_script_info@multi_model_mean .eq. "y") then
    val_ext_stddev = new((/4, end_year - start_year + 1/), float)
    val_area_stddev = new((/4, end_year - start_year + 1/), float)
  else  ; initialize anyway: both are needed as parameters for xy_line
    val_ext_stddev = 0
    val_area_stddev = 0
  end if

  ; Loop over datasets
  do idat = 0, dimsizes(annots) - 1

    log_info("processing " + annots(idat) + \
             "  ############################################")
    ; See ./interface_scripts/data_handling.ncl
    sic = read_data(idat, var0, field_type0)
    if (sic@units.eq."%") then
      sic = sic / 100.0  ; Subsequent code expects %/100
      sic@units = "%/100"
    end if

    ; Further processing depends on whether sic is on a lat-lon grid
    if(all(iscoord(sic, (/"lat", "lon"/)))) then
      log_info("SeaIce_tsline.ncl: " + \
               "sic is on a lat - lon grid")
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
      log_info("SeaIce_tsline.ncl: " + \
               "sic is not on a lat-lon grid")

      ; Get area data & lat-lon from another file
      ; (dataset specific: ~/reformat_scripts folder)
      aux_vars  = read_fx_data(idat, "areacello")
      area_lon  = aux_vars[2]  ; lon
      area_lat  = aux_vars[1]  ; lat
      area_temp = aux_vars[0]  ; areacello
      delete(aux_vars)
    end if

    ; Time axis
    sic_temp = time_operations(sic, toint(dataset_info@start_year(idat)),\
                               toint(dataset_info@end_year(idat)), \
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
        if (all(ismissing(data1D))) then  ; Time step missing
          continue
        end if
        do while(1.lt.2)  ; Will break when hole is filled
          i1D = ind(lat1D.ge.actual_lat .and. \
                    lat1D.lt.actual_lat+latstep)
          if (.not. all(ismissing(i1D))) then
            i2D     = ind_resolve(ind(lat1D.ge.actual_lat), \
                                  dimsizes(area_lat))
            dim_i2D = dimsizes(i2D)
            npts    = dim_i2D(0)  ; Number of elements
            ndim    = dim_i2D(1)  ; Rank of "lat"
            if (any(ismissing(data1D(i1D)))) then
              do n = 0, npts - 1
                ring = sic(itime, i2D(n, 0), i2D(n, 1))
                sic(itime, i2D(n, 0), i2D(n, 1)) = \
                  where(ismissing(ring), 1., ring)
                delete(ring)
              end do
            else  ; Out of the hole!
              delete([/i1D, i2D, dim_i2D, npts, ndim/])
              break
            end if
            delete([/i2D, dim_i2D, npts, ndim/])
          end if
          delete([/i1D/])
          actual_lat = actual_lat - latstep
        end do  ; while
        delete(data1D)
      end do  ; itime
      delete(lat1D)
    end if

    ; Select region -> name corresponds to the pole region, but entire
    ; hemisphere will be evaluated and converted to million sqkm
    if (diag_script_info@region .eq. "Antarctic") then
                                ; Set area_temp to 0 in NH
      area_temp = where(area_lat .lt. 0.0, area_temp / 1.0e12, 0.0)
    else if (diag_script_info@region .eq. "Arctic")
                                ; Set area_temp to 0 in SH
      area_temp = where(area_lat .gt. 0.0, area_temp / 1.0e12, 0.0)
    end if
    end if

    ; Calculate sea ice area
    sic_area = sic * conform(sic, area_temp, (/1, 2/))
    copy_VarMeta(sic, sic_area)

    ; Calculate sea ice extent
    extent_temp = where(sic .ge. 0.15, 1., 0.)
    sic_ext = conform(sic, area_temp, (/1, 2/)) * extent_temp
    delete(extent_temp)
    copy_VarMeta(sic, sic_ext)
    delete(sic)

    ; Calculate weighted annual mean or select month
    if (diag_script_info@month .eq. "A") then    ; weighted annual mean
      sic_ann_ext  = time_operations(sic_ext, \
                                     toint(dataset_info@start_year(idat)), \
                                     toint(dataset_info@end_year(idat)), \
                                     "average", "yearly", True)
      sic_ann_area = time_operations(sic_area, \
                                     toint(dataset_info@start_year(idat)),\
                                     toint(dataset_info@end_year(idat)), \
                                     "average", "yearly", True)
    else  ; Select month
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
    ia = toint(dataset_info@start_year(idat)) - start_year
    ie = toint(dataset_info@end_year(idat)) - start_year
    val_ext(idat, ia:ie) = tofloat(sic_sum_ext)
    val_area(idat, ia:ie) = tofloat(sic_sum_area)

    ; Clean up
    delete([/area_lon, area_lat, area_temp, sic_sum_ext, sic_sum_area, \
           ia, ie/])

  end do  ; idat

  ; Calculate standard deviation of datasets
  if (diag_script_info@multi_model_mean .eq. "y") then
    ; See style_$project$.ncl for which data shall
    ; be included in the statistics
    temp = ind(avgstd .eq. 0)
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
  if diag_script_info@month .eq. "A" then   ; weighted annual mean
    val_ext&time@units = "Annual"
  else  ; select a month
    val_ext&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_ext@units = "10~S~6~N~ km~S~2~N~"
  val_ext@long_name = "Sea Ice Extent"

  ; Some attributes for area plot
  val_area!1 = "time"
  val_area&time = int2dble(all_years)
  if diag_script_info@month .eq. "A" then   ; weighted annual mean
    val_area&time@units = "Annual"
  else  ; select a month
    val_area&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_area@units = "10~S~6~N~ km~S~2~N~"
  val_area@long_name = "Sea Ice Area"

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  if (dimsizes(val_ext&time) .eq. 1) then
    log_info(" Skipping tsline: " + \
             "plot needs at least two time values, has 1")
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

    ; Draw the plots (procedure in ./plot_scripts/ncl/xy_line.ncl)
    debuginfo = False
    xy_line(wks_ext, val_ext, val_ext&time, val_ext_stddev, res_ext, \
            debuginfo)
    if (isatt(diag_script_info, "legend_outside")) then
      ; Legend can only once be drawn in an extra file (hlu error),
      ; but is same for both area and extent anyway
      if (diag_script_info@legend_outside) then
        val_area@legend_outside = True
      end if
    else
      diag_script_info@legend_outside = False
    end if
    xy_line(wks_area, val_area, val_area&time, val_area_stddev, res_area, \
            debuginfo)

  end if  ; More than 1 time step

  ; ---------------------------------------------------------------------------
  ; Write NETCDF file for MDER
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_netcdf) then
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)
    opt = True

    ; Iterate over all desired diagnostics
    do dd = 0, dimsizes(diag_script_info@wdiag) - 1
      opt@wdiag = diag_script_info@wdiag(dd)
      log_info(opt@wdiag)

      ; Iterate over all datasets
      do i = 0, DIM_DAT - 1
        idat = idx_datasets(i)
        log_debug("Processing " + dataset_info@dataset(idat))
        new_path = work_dir + diag_file_name(idat, opt)

        ; Write files
        if (isStrSubset(opt@wdiag, "_t")) then
          val = val_trend(case|:, dataset|i)
          val@var = var0
          val@units = "(" + variable_info@units + "/dec)"
          val@diag_script = (/diag_script/)
          copy_VarAtts(diag_script_info, val)
          val@dataset = dataset_info@dataset(idat)
          val@start_year = input_file_info@start_year(idat)
          val@end_year = input_file_info@end_year(idat)
          val@ncdf = new_path
          delete(val)
        else if(isStrSubset(opt@wdiag, "_c")) then
          val = val_mean(case|:, dataset|i)
          val@var = var0
          val@units = "(" + variable_info@units + "/dec)"
          val@diag_script = (/diag_script/)
          copy_VarAtts(diag_script_info, val)
          val@dataset = dataset_info@dataset(idat)
          val@start_year = input_file_info@start_year(idat)
          val@end_year = input_file_info@end_year(idat)
          val@ncdf = new_path
          ncdf_outfile = ncdf_write(val, new_path)
        end if
        end if
        delete([/val, new_path/])
      end do
    end do
  end if

  leave_msg(diag_script, "")

end
