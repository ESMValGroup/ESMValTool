; #############################################################################
; russell_figure3b.ncl (Subantarctic Front)
;
; Russell, J.L.,et al., 2018, J. Geophysical Research – Oceans, 123, 3120-3143
; https://doi.org/10.1002/2017JC013461 (figure 3b - Subantarctic Front)
;
; Author:  Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;
; #############################################################################
; Description
;   - Uses original grid (no regridding).
;   - Changes i-j and x-y coordinates to lat - lon coordinates
;   - Takes the time average of temperature
;   - Makes sure the temperature values are in Kelvin
;   - extracts the temperature closest to and less than 400m depth
;   - aranges the temperature array in ascending order of lon (dataset)
;   - creates as a contour map of ascending sorted temperature
;     (no actual plotting of contour map)
;   - extracts the position of 4 degree celsius isotherm from the above plot
;   - increase the lon range from 0 - 360 to -360 to 360
;   - plots the isoline after removing stray points from isoline
;   - overlays all the plots on the first plot
;   - creates legend
;
; Required diag_script_info attributes
;
;   - styleset     : CMIP5         - default
;    - ncdf        : default
;
; Required Preprocessor attributes ( no_preprocessor )
;   - none (no preprocessing required)
;     in the recipe do not keep a preprocessor in variable section
;
; Caveats
;
;   - IPSL models do not work as they have 360+ range in longitude.
;   - CNRM models give a plotting warning, but the plot looks fine.
;
;
; Modification history
;
;   20190510- russell_joellen, pandde_amarjiit - written and
;              implemented for ESMValTool v2
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  datasetnames = metadata_att_as_array(info_items, "dataset")
  start_years_data = metadata_att_as_array(info_items, "start_year")
  end_years_data = metadata_att_as_array(info_items, "end_year")
  inputfile_paths = metadata_att_as_array(info_items, "filename")
  nDatasets = ListCount(info_items)
end

begin

  ; Set annotations
  annots = project_style(input_file_info, diag_script_info, "annots")
  colors = project_style(input_file_info, diag_script_info, "colors")
  dashes = project_style(input_file_info, diag_script_info, "dashes")
  thicks = project_style(input_file_info, diag_script_info, "thicks")

  plotpath = config_user_info@plot_dir + "/Russell18jgr_fig3_Subantarctic" \
    + "-Fronts_" + sprinti("%0.4i", min(toint(start_years_data))) \
    + "-" + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  plot = new(nDatasets, graphic)
  plot1 = new(nDatasets, graphic)
  res                     = True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@trYMaxF             = -40.
  res@trYMinF             = -70.
  res@trXMaxF             = 300.
  res@trXMinF             = -60.
  res@vpXF                = 0.1                   ; orientation on page
  res@vpYF                = 0.875
  res@vpWidthF            = 0.625                     ; width
  res@vpHeightF           = 0.815                    ; height
  res@gsnLeftString       = "Russell et al -2018 - Figure 3 b"
  res@tmXBMinorPerMajor   = 1
  res@tmXBMaxTicks        = 15
  res@tmXBTickStartF      = -60.
  res@tmXBTickSpacingF    = 30.
  res@tmXBTickEndF        = 300.
  res@tmYLTickStartF      = -70.
  res@tmYLTickSpacingF    = 2.
  res@tmYLTickEndF        = -40.
  res@tmYLMinorPerMajor   = 1
  res@tmYLMaxTicks        = 20
  res@tmYLMode            = "Manual"
  res@tmXBLabelFontHeightF = 0.0125
  res@tmYLLabelFontHeightF = 0.0145
  res@gsnRightStringFontHeightF = 15.
  res@gsnLeftStringFontHeightF  = 15.
  res@tiMainString   = " Subantarctic Fronts"

  res1 = True
  res1@gsnDraw = False
  res1@gsnFrame = False

  do iii = 0, nDatasets - 1
    dataset_time = read_data(info_items[iii])
    if (max(dataset_time) .gt. 273.0) then
      dataset_time@units = "K"
      if (min(dataset_time) .eq. 0) then  ; some models have 0 as fillvalues
        dataset_time = where(dataset_time .eq.0, dataset_time@_FillValue, \
                             dataset_time)
      end if
    else
      dataset_time = dataset_time + 273.0
      dataset_time@units = "converted to K"
    end if

    dataset_lev = dim_avg_n_Wrap(dataset_time, 0)
    delete(dataset_time)
    lev_data = dataset_lev&lev
    lev_a = closest_val(400.0, lev_data)

    if (lev_data(lev_a) .lt. 400.0) then
      lev_new = lev_data(lev_a)
      dataset_1 = dataset_lev(lev_a, :, :)
    else
      lev_new = lev_data(lev_a - 1)
      dataset_1 = dataset_lev(lev_a - 1, :, :)
    end if
    delete(dataset_lev)
    delete(lev_data)
    delete(lev_new)

    res@xyDashPatterns              = dashes(iii)
    res@xyLineColors                = colors(iii)
    res@xyExplicitLegendLabels      = annots(iii)
    res@xyLineThicknessF            = thicks(iii)

    if (iscoord(dataset_1, "i")) then  ; changes i-j coordinates to lat-lon
      delete(dataset_1&i)
      delete(dataset_1&j)
      infile_path = inputfile_paths(iii)
      infile_iii  = addfile(infile_path, "r")
      area_lon    = infile_iii->lon
      area_lat    = infile_iii->lat
      dataset_1!0   = "lat"
      dataset_1!1   = "lon"
      dataset_1&lat = tofloat(area_lat(:, 10))
      dataset_1&lon = tofloat(area_lon(10, :))
      delete(area_lat)
      delete(area_lon)
    end if

    if(.not.(iscoord(dataset_1, "lat"))) then
      infile_path = inputfile_paths(iii)
      infile_iii  = addfile(infile_path, "r")
      area_lon = infile_iii->lon
      area_lat = infile_iii->lat
      a1 = closest_val(-35.0, area_lat(:, 0))
      a2 = closest_val(-80.0, area_lat(:, 0))
      ; sort the temperature variable in ascending order of lon
      ip0 = dim_pqsort(area_lon(0, :), 1)
      dataset = dataset_1(a1:a2, ip0)
      dataset@lon2d = area_lon(a1:a2, ip0)
      dataset@lat2d = area_lat(a1:a2, ip0)
      delete(area_lat)
      delete(area_lon)
    else
      a1 = closest_val(-35.0, dataset_1&lat)
      a2 = closest_val(-80.0, dataset_1&lat)
      ; sort the temperature variable in ascending order of lon
      ip0 = dim_pqsort(dataset_1&lon, 1)
      dataset = dataset_1(a1:a2, ip0)
    end if
    delete(ip0)
    delete(dataset_1)

    plot1(iii) = gsn_csm_contour(wks, dataset, res1)    ; create the plot
    delete(dataset)

    nc_filename = config_user_info@work_dir + "russell18jgr_fig3b_" \
      + "subantarctic-front-position_" + annots(iii) + "_" \
      + (start_years_data(iii)) + "-" + (end_years_data(iii)) + ".nc"

    isoline = get_isolines(plot1(iii), 277.15)
    b = isoline@start_point(0)
    e = b + isoline@n_points(0) - 1

    ; NCL picks a few isolated points on the isoline in some of the models
    ; we are discarding them if the difference in lon is more than 20 degrees

    if (isoline@segment_count .gt.1) then
        do segmod = 1, isoline@segment_count - 1
        start_point1 = isoline@start_point(segmod)
        end_point1 = isoline@start_point(segmod-1) \
          + isoline@n_points(segmod-1) - 1

        if (abs(isoline(1, start_point1) - isoline(1, e)) .lt. 20) then
            e = e + isoline@n_points(segmod)
        end if
        end do
    end if
    last_index = e-b

    y = new((3*(e + 1 - b)), typeof(isoline))
    x = new((3*(e + 1 - b)), typeof(isoline))
    x@units = "Degrees_east"
    y@units = "Degrees_north"

    ; making the isoline repeat 3 times to make lon range from -360 to 360
    if (isoline(1, 11) .gt. isoline(1, 1)) then
        y(0 : last_index)                    = isoline(0, b:e)
        y(last_index+1 : 2*last_index + 1)   = isoline(0, b:e)
        y(2*last_index+2 : 3*last_index + 2) = isoline(0, b:e)

        x(0 : last_index)                    = isoline(1, b:e) - 360
        x(last_index+1 : 2*last_index + 1)   = isoline(1, b:e)
        x(2*last_index+2 : 3*last_index + 2) = isoline(1, b:e) + 360

    else
        y(0 : last_index)                    = isoline(0, b:e:-1)
        y(last_index+1 : 2*last_index + 1)   = isoline(0, b:e:-1)
        y(2*last_index+2 : 3*last_index + 2) = isoline(0, b:e:-1)

        x(0 : last_index)                    = isoline(1, b:e:-1) - 360
        x(last_index+1 : 2*last_index + 1)   = isoline(1, b:e:-1)
        x(2*last_index+2 : 3*last_index + 2) = isoline(1, b:e:-1) + 360

    end if

    plot(iii) = gsn_csm_xy(wks, x, y, res)

    if (iii .ne. 0) then
      overlay(plot(0), plot(iii))
    end if

    out_var = new((/2, dimsizes(x) /), float)
    out_var(0, :) = x
    out_var(1, :) = y
    out_var!0 = "points"
    out_var!1 = "xy"
    out_var&points = (/0, 1/)
    out_var&xy = ispan(0, (dimsizes(x) - 1), 1)
    out_var@var = "position_of_subantarctic_front"
    out_var@diag_script = DIAG_SCRIPT
    out_var@description = "out_var(0,:) are lon positions of isoline " \
      + " and out_var(1,:) is the lat position of isoline - (277.15K)"
    ncdf_outfile = ncdf_write(out_var, nc_filename)
    delete(isoline)
    delete(x)
    delete(y)
    delete(b)
    delete(e)
    delete(last_index)
    delete(out_var)
    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 3b", \
                   "mean", \
                   "sh", \
                   "geo", \
                   "russell_joellen", \
                   "russell18jgr", \
                   inputfile_paths(iii))
  end do
  draw(plot(0))

  legend = create "Legend" legendClass wks
  "vpXF"                     : 0.625                   ; orientation on page
  "vpYF"                     : 0.925
  "vpWidthF"                 : 0.5                     ; width
  "vpHeightF"                : 0.725                   ; height
  "lgPerimOn"                : False                   ; no perimeter
  "lgItemCount"              : dimsizes(annots)        ; how many
  "lgLineLabelStrings"       : annots                  ; labels
  "lgLabelsOn"               : False                   ; no default lables
  "lgLineLabelFontHeightF"   : 0.0085                  ; font height
  "lgDashIndexes"            : dashes                  ; line paterns
  "lgLineColors"             : colors
  "lgMonoLineLabelFontColor" : True                    ; one label color
  end create
  draw(legend)

end
