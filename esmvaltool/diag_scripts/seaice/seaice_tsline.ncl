; #############################################################################
; LINE PLOTS OF SEA ICE AREA AND EXTENT vs TIME
; Author: Daniel Senftleben (DLR, Germany)
; ESMVal project
; #############################################################################
;
; Description
;   Creates line plot for sea ice area and extent versus time.
;
; Required diag_script_info attributes (diagnostics specific)
;   month: selected month or annual mean.
;   region: region to be plotted ( "Arctic" or "Antarctic").
;
; Optional diag_script_info attributes (diagnostics specific)
;   fill_pole_hole: fill observational hole at North pole, Default: False
;   styleset: specify line and color style
;   EMs_in_lg: create legend label for each individual ensemble member
;              (Default: "True")
;   xy_line_legend: False: do not draw a legend inside the plot
;   legend_outside: True: draw legend in an extra plot
;
; Required variable_info attributes (variable specific)
;
; Optional variable_info attributes (variable specific)
;
; Caveats
;   CMOR units expected: areacello in m2, sic in %
;
; Modification history
;   20190227-senftleben_daniel: added provenance.
;   20190218-senftleben_daniel: ported to ESMValTool v2.0a1
;   20180620-senftleben_daniel: ported to ESMValTool v2
;   20160906-senftleben_daniel: added option to fill polar observation hole
;   20151027-lauer_axel: moved call to 'write_references' to the beginning
;   20150325-lauer_axel: modified reference tags used for acknowledgements
;   20140312-lauer_axel: adapted to new time_operations structure
;   20140129-senftleben_daniel: modularized to fit to new structure
;   20130419-gottschaldt_klaus-dirk: implemented into ./diag_scripts
;   20130405-gottschaldt_klaus-dirk: written
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/seaice/seaice_aux.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  ; Get metadata items
  att = True
  att@mip = "OImon"
  info = select_metadata_by_atts(input_file_info, att)  ; variable
  var0 = info[0]@short_name
  datasets = metadata_att_as_array(info, "dataset")
  n_datasets = dimsizes(datasets)
  att@mip = "fx"
  areainfo = select_metadata_by_atts(input_file_info, att)  ; area
  areadatasets = metadata_att_as_array(areainfo, "dataset")
  delete(att)

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"region", "month"/))

  ; Check region selection
  if (all(diag_script_info@region.ne.(/"Arctic", "Antarctic"/))) then
    error_msg("f", DIAG_SCRIPT, "", "no (valid) region selection")
  end if

  ; Check fill_pole_hole option
  if (.not.isatt(diag_script_info, "fill_pole_hole")) then
    diag_script_info@fill_pole_hole = False
  else
    if (diag_script_info@fill_pole_hole) then
      if (diag_script_info@region.eq."Antarctic") then
        diag_script_info@fill_pole_hole = False
        error_msg("w", DIAG_SCRIPT, "", "No pole hole at Antarctic " + \
                  "- setting diag_script_info@fill_pole_hole to False")
      end if
    end if
  end if

  ; Check multi model mean option
  if (isatt(diag_script_info, "multi_model_mean")) then
    multi_model_mean = diag_script_info@multi_model_mean
  else
    multi_model_mean = False
  end if

  ; Determine time range
  start_years = metadata_att_as_array(info, "start_year")
  start_year  = min(start_years)
  end_years   = metadata_att_as_array(info, "end_year")
  end_year    = max(end_years)
  all_years   = ispan(start_year, end_year, 1)

  ; Select attributes
  avgstd = project_style(info, diag_script_info, "avgstd")

  ; Create output plot directory
  if (diag_script_info@write_plots) then
    plot_dir = config_user_info@plot_dir
    system("mkdir -p " + plot_dir)
  end if

  ; Plot file type
  file_type = config_user_info@output_file_type

  ; For annotation
  smonths = (/"January", "February", "March", "April", \
              "May", "June", "July", "August", \
              "September", "October", "November", "December"/)
  i2months = sprinti("%0.2i", ispan(1, 12, 1))

  ; Output file name specifications (passed to get_wks for output file name)
  system("mkdir -p " + config_user_info@work_dir)
  if (diag_script_info@month .eq. "A") then
    i2month = "annual"
  else
    i2month = smonths(toint(diag_script_info@month) - 1)
  end if
  out_ext = "extent_" + var0 + "_" \
    + diag_script_info@region + "_" + i2month + "_" \
    + sprinti("%0.4i", start_year) + "-" \
    + sprinti("%0.4i", end_year)
  out_area = "area_" + var0 + "_" \
    + diag_script_info@region + "_" + i2month + "_" \
    + sprinti("%0.4i", start_year) + "-" \
    + sprinti("%0.4i", end_year)

  ; Array for collecting sea-ice extent for all datasets
  val_ext = new((/n_datasets, end_year - start_year + 1/), float)
  val_ext!1 = "time"
  val_ext&time = int2dble(all_years)

  ; Array for collecting sea-ice area for all datasets
  val_area = new((/n_datasets, end_year - start_year + 1/), float)
  val_area!1 = "time"
  val_area&time = int2dble(all_years)

  ; Prepare calculation of standard deviation of datasets
  if (multi_model_mean) then
    val_ext_stddev = new((/4, end_year - start_year + 1/), float)
    val_area_stddev = new((/4, end_year - start_year + 1/), float)
  else  ; initialize anyway: both are needed as parameters for xy_line
    val_ext_stddev = 0
    val_area_stddev = 0
  end if

  ; Loop over datasets
  do imod = 0, n_datasets - 1

    log_info("Processing dataset " + datasets(imod))

    sic = read_data(info[imod])

    ; Check grid type
    l_irreg = True
    if (sic!1.eq."lat" .and. sic!2.eq."lon") then
      if (dimsizes(dimsizes(sic&lat)).eq.1 .and. \
          dimsizes(dimsizes(sic&lon)).eq.1) then
        l_irreg = False
      end if
    end if

    ; Look for area variable
    if (any(datasets(imod).eq.areadatasets)) then

      ; Read areacello
      log_info("Using variable areacello for the area")
      area = read_data(areainfo[imod])
      locinfo = areainfo[imod]
      locinfo@short_name = "lat"  ; read auxiliary latitude coordinate
      arealat = read_data(locinfo)
      delete(locinfo)

      ; Check grid consistency
      if (.not.isconform(sic(0, :, :), area)) then
        if (l_irreg) then
          error_msg("f", DIAG_SCRIPT, "", "areacello not conformal")
        else
          delete(area)
          delete(arealat)
          log_info("areacello has inconsistent dimensions, " + \
                   "calculating area for a regular grid instead")
          area = map_area(sic&lat, sic&lon)
          arealat = area&lat
        end if
      end if

    else

      if (.not.l_irreg) then
        log_info("Calculating area for a regular grid")
        area = map_area(sic&lat, sic&lon)
        arealat = area&lat
      else
        error_msg("f", DIAG_SCRIPT, "", "cannot calculate area")
      end if

    end if

    ; Calcluate sea ice area
    sic_sum_ext  = sea_ice_area(sic, area, arealat, l_irreg, \
                                diag_script_info@region, \
                                start_years(imod), end_years(imod), \
                                diag_script_info@month, "extent", \
                                diag_script_info@fill_pole_hole)

    ; Calculate sea ice extent
    sic_sum_area = sea_ice_area(sic, area, arealat, l_irreg, \
                                diag_script_info@region, \
                                start_years(imod), end_years(imod), \
                                diag_script_info@month, "area", \
                                diag_script_info@fill_pole_hole)

    ; Put into global array
    ia = toint(start_years(imod)) - start_year
    ie = toint(end_years(imod)) - start_year
    val_ext(imod, ia:ie) = tofloat(sic_sum_ext)
    val_area(imod, ia:ie) = tofloat(sic_sum_area)

    ; Clean up
    delete([/sic, area, arealat, sic_sum_ext, sic_sum_area, ia, ie/])

  end do  ; imod

  ; Calculate standard deviation of datasets
  if (multi_model_mean) then
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
  val_ext!0 = "model"
  val_ext&model = metadata_att_as_array(info, "dataset")
  val_ext!1 = "time"
  val_ext&time = int2dble(all_years)
  if (diag_script_info@month .eq. "A") then   ; weighted annual mean
    val_ext&time@units = "Annual"
  else  ; select a month
    val_ext&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_ext@units = "10~S~6~N~ km~S~2~N~"
  val_ext@long_name = "Sea Ice Extent"

  ; Some attributes for area plot
  val_area!0 = "model"
  val_area&model = metadata_att_as_array(info, "dataset")
  val_area!1 = "time"
  val_area&time = int2dble(all_years)
  if (diag_script_info@month .eq. "A") then   ; weighted annual mean
    val_area&time@units = "Annual"
  else  ; select a month
    val_area&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_area@units = "10~S~6~N~ km~S~2~N~"
  val_area@long_name = "Sea Ice Area"

  ; Plotting
  if (dimsizes(val_ext&time) .eq. 1) then
    log_info(" Skipping tsline: " + \
             "plot needs at least two time values, has 1")
  else

    ; Create plot variables
    wks_ext = gsn_open_wks(file_type, config_user_info@plot_dir + out_ext)
    wks_area = gsn_open_wks(file_type, config_user_info@plot_dir + out_area)

    ; Define plot resources
    res                 = True
    res@tmXBLabelAngleF = 45
    res@tmXBLabelJust   = "CenterRight"

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
    xy_line(wks_ext, val_ext, val_ext&time, val_ext_stddev, res_ext, info)
    if (isatt(diag_script_info, "legend_outside")) then
      ; Legend can only once be drawn in an extra file (hlu error),
      ; but is same for both area and extent anyway
      if (diag_script_info@legend_outside) then
        val_area@legend_outside = True
        wks_area@legendfile = config_user_info@plot_dir + out_ext + "_legend"
      end if
    else
      diag_script_info@legend_outside = False
    end if
    xy_line(wks_area, val_area, val_area&time, val_area_stddev, res_area, info)

    ; Provenance
    if (diag_script_info@region.eq."Arctic") then
      domain = "nh"
    elseif (diag_script_info@region.eq."Antarctic") then
      domain = "sh"
    end if
    infiles = metadata_att_as_array(info, "filename")

    val_ext@var = var0
    val_ext@diag_script = DIAG_SCRIPT
    caption_ext = "Timeseries of " + \
      smonths(toint(diag_script_info@month) - 1) + " " + \
      diag_script_info@region + " sea ice extent in " + val_ext@units
    nc_ext = \
      ncdf_write(val_ext, config_user_info@work_dir + out_ext + ".nc")

    log_provenance(nc_ext, out_ext, caption_ext, \
                   (/"mean", "stddev", "clim"/), domain, "times", \
                   (/"senftleben_daniel", "gottschaldt_klaus-dirk"/), \
                   "stroeve07grl", infiles)

    val_area@var = var0
    val_area@diag_script = DIAG_SCRIPT
    caption_area = "Timeseries of " + \
      smonths(toint(diag_script_info@month) - 1) + " " + \
      diag_script_info@region + " sea ice area in " + val_area@units
    nc_area = \
      ncdf_write(val_area, config_user_info@work_dir + out_area + ".nc")

    log_provenance(nc_area, out_area, caption_area, \
                   (/"mean", "stddev", "clim"/), domain, "times", \
                   (/"senftleben_daniel", "gottschaldt_klaus-dirk"/), \
                   "stroeve07grl", infiles)

  end if  ; More than 1 time step

  leave_msg(DIAG_SCRIPT, "")

end
