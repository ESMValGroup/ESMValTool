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
;   20190227-A_senf_da: Added provenance.
;   20190218-A_senf_da: Ported to ESMValTool v2.0a1
;   20180620-A_senf_da: Ported to ESMValTool v2
;   20160906-A_senf_da: Added option to fill polar observation hole
;   20151027-A_laue_ax: Moved call to 'write_references' to the beginning
;                       of the code
;   20150325-A_laue_ax: Modified reference tags used for acknowledgements
;                       (projects, observations, etc)
;   20140312-A_righ_ma: Adapted to new time_operations structure.
;   20140129-A_senf_da: Modularized to fit to new structure.
;   20130419-A_gott_kl: Implemented into ./diag_scripts.
;   20130405-A_gott_kl: Implemented.
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

  n_variables = ListCount(variable_info)
  var0 = variable_info[0]@short_name

end


begin

  ; Get metadata items
  info = select_metadata_by_name(input_file_info, var0)
  n_datasets = ListCount(info)

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
  start_years = metadata_att_as_array(input_file_info, "start_year")
  start_year  = min(start_years)
  end_years   = metadata_att_as_array(input_file_info, "end_year")
  end_year    = max(end_years)
  all_years   = ispan(start_year, end_year, 1)

  ; Select attributes
  avgstd = project_style(info, diag_script_info, "avgstd")

  ; For annotation
  smonths = (/"January", "February", "March", "April", \
              "May", "June", "July", "August", \
              "September", "October", "November", "December"/)
  i2months = sprinti("%0.2i", ispan(1, 12, 1))

  ; Output file name specifications
  ; (passed to get_wks for final output file name)
  if (diag_script_info@month .eq. "A") then
    i2month = "annual"
  else
    i2month = smonths(toint(diag_script_info@month) - 1)
  end if
  outfile_ext = "extent_" + var0 + "_"\
    + diag_script_info@region + "_" + i2month + "_"\
    + sprinti("%0.4i", start_year) + "-"\
    + sprinti("%0.4i", end_year)
  outfile_area = "area_" + var0 + "_"\
    + diag_script_info@region + "_" + i2month + "_"\
    + sprinti("%0.4i", start_year) + "-"\
    + sprinti("%0.4i", end_year)

  ; Start processing ##########################################
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

    sic = read_data(info[imod])

    ; Calcluate sea ice area / extent
    sic_sum_ext  = sea_ice_area(sic, imod, info, diag_script_info@region, \
                                start_years(imod), end_years(imod), \
                                diag_script_info@month, "extent", \
                                diag_script_info@fill_pole_hole)
    sic_sum_area = sea_ice_area(sic, imod, info, diag_script_info@region, \
                                start_years(imod), end_years(imod), \
                                diag_script_info@month, "area", \
                                diag_script_info@fill_pole_hole)

    ; Put into global array
    ia = toint(start_years(imod)) - start_year
    ie = toint(end_years(imod)) - start_year
    val_ext(imod, ia:ie) = tofloat(sic_sum_ext)
    val_area(imod, ia:ie) = tofloat(sic_sum_area)

    ; Clean up
    delete([/sic, sic_sum_ext, sic_sum_area, ia, ie/])

  end do  ; imod

  ; Calculate standard deviation of datasets
  if (multi_model_mean) then
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
  val_ext!0 = "model"
  val_ext&model = metadata_att_as_array(info, "dataset")
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
  val_area!0 = "model"
  val_area&model = metadata_att_as_array(info, "dataset")
  val_area!1 = "time"
  val_area&time = int2dble(all_years)
  if diag_script_info@month .eq. "A" then   ; weighted annual mean
    val_area&time@units = "Annual"
  else  ; select a month
    val_area&time@units = smonths(stringtoint(diag_script_info@month) - 1)
  end if
  val_area@units = "10~S~6~N~ km~S~2~N~"
  val_area@long_name = "Sea Ice Area"

  ; *********** PLOTTING ************
  if (dimsizes(val_ext&time) .eq. 1) then
    log_info(" Skipping tsline: " + \
             "plot needs at least two time values, has 1")
  else
    ; Create plot variables
    wks_ext = get_wks(0, DIAG_SCRIPT, outfile_ext)
    wks_area = get_wks(0, DIAG_SCRIPT, outfile_area)

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
      end if
    else
      diag_script_info@legend_outside = False
    end if
    xy_line(wks_area, val_area, val_area&time, val_area_stddev, res_area, info)

    ; Provenance
    statistics = (/"mean", "stddev", "clim"/)
    domain = diag_script_info@region
    plottype = "tsline"
    authors = (/"senf_da", "gott_kl"/)
    references = "stroeve07grl"
    infiles = metadata_att_as_array(info, "filename")
    val_ext@var = var0
    val_ext@diag_script = DIAG_SCRIPT
    if (isatt(wks_ext, "fullname")) then
      outfile_ext = wks_ext@fullname
    else
      outfile_ext = outfile_ext + "." + output_type()
    end if
    caption_ext = "Timeseries of " + \
      smonths(toint(diag_script_info@month) - 1) + " " + \
      diag_script_info@region + " sea ice extent in " + val_ext@units
    nc_file_ext = ncdf_write(val_ext, outfile_ext + ".nc")
    log_provenance(nc_file_ext, outfile_ext, caption_ext, statistics, \
                   domain, plottype, authors, references, infiles)
    val_area@var = var0
    val_area@diag_script = DIAG_SCRIPT
    if (isatt(wks_area, "fullname")) then
      outfile_area = wks_area@fullname
    else
      outfile_area = outfile_area + "." + output_type()
    end if
    caption_area = "Timeseries of " + \
      smonths(toint(diag_script_info@month) - 1) + " " + \
      diag_script_info@region + " sea ice area in " + val_area@units
    nc_file_area = ncdf_write(val_area, outfile_area + ".nc")
    log_provenance(nc_file_area, outfile_area, caption_area, statistics, \
                   domain, plottype, authors, references, infiles)

  end if  ; More than 1 time step

  leave_msg(DIAG_SCRIPT, "")

end
