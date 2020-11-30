; #############################################################################
;  russell18jgr_fig9b.ncl
;
; Russell, J.L.,et al., 2018, J. Geophysical Research – Oceans, 123, 3120-3143
;  https://doi.org/10.1002/2017JC013461 (figure 9b)
;
; Author:  Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;          Pandde Amarjiit (University of Arizona, USA)
; ESMVal project
; #############################################################################
; Description
;  - Uses original grid (no regridding).
;    (for Total carbon flux part of script)
;  - Maskes the fgco2 data points on the land using ESMValTool preprocessor.
;  - if var coordinates are lat-lon then the script calculates the areacello
;      otherwise, the script reads the areacello variable from input file.
;  - Multiplies time averaged fgco2 with areacello to get flux per cell.
;  - flux per lat = sum of all cells across a lat.
;  - total carbon flux in Southern ocean = sum of flux per lat from 30S - 90S.
;
;    (for Southern westerly band part of the script)
;  - Uses ESMValTool land masking (preprocessor) to mask land data points.
;  - Calculates longitudinal averaged tauu/tauuo values.
;  - Calculates the first lat north and south of 50S where tauu is zero.
;  - Subtracts the north lat value from south lat value to get
;      Latitudinal width of Southern Hemisphere Westerly Band.
;  - Calculates line of linear regression between lat width and
;      total carbon flux.
;  - Plots the line of linear regression.
;  - Plots each model's Carbon flux and lat width as markers.
;
; Required diag_script_info attributes (configured for russell figure 9a)
;      None
;
; Required variable_info_attributtes:
;   - preprocessor
;   - mip
;   - project
;   - exp
;   - ensemble
;   - start_year
;   - end_year
;   - additional datasets
;
; Required preprocessor attributtes:
;   preprocessor_time - for tauuo and hfds
;    climate_statistics:
;      operator: mean
;      period: full
;
; Required diag_script_info attributes (diagnostics specific)
;
;    - styleset : CMIP5 - default
;    - ncdf     : default
;
; Caveats
;  - the lat width calculation does not match exactly with one in paper
;
; Modification history
;
;   20190510 - russell_joellen, pandde_amarjiit - written
;             and implemented for ESMValTool v2.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function
load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  fgco2_items = select_metadata_by_name(input_file_info, "fgco2")

  var_names = metadata_att_as_array(variable_info, "short_name")
  var_ind = str_match_ind(var_names, "tauu")
  tauu_var = var_names(var_ind)
  ; to run either tauu or tauuo
  tauu_items = select_metadata_by_name(input_file_info, tauu_var)
  tauu_inputfile_paths = metadata_att_as_array(tauu_items, "filename")

  areacello_items = select_metadata_by_name(input_file_info, "areacello")

  fgco2_inputfile_paths = metadata_att_as_array(fgco2_items, "filename")
  start_years_data = metadata_att_as_array(fgco2_items, "start_year")
  end_years_data   = metadata_att_as_array(fgco2_items, "end_year")
  nDatasets = ListCount(tauu_items)
  datasetnames = metadata_att_as_array(fgco2_items, "dataset")

end

begin

  plotpath = config_user_info@plot_dir + "/russell18jgr-fig9b_" \
            + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
            + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  wks@fullname = plotpath
  plots = new(nDatasets, graphic)
  lat_width = new(nDatasets, double)
  carbon_flux_iii = new(nDatasets, double)

  data1 = new((/2, dimsizes(datasetnames)/), double)
  aaa = new(dimsizes(datasetnames), double)

  ; Set annotations
  annots = project_style(tauu_items, diag_script_info, "annots")
  markers = project_style(tauu_items, diag_script_info, "markers")
  colors = project_style(tauu_items, diag_script_info, "colors")

  res                            = True             ; plot mods desired
  res@gsnMaximize                = False            ; maximize plot in frame
  res@gsnFrame                   = False
  res@gsnDraw                    = True
  res@tiYAxisString              = "Southern ocean carbon uptake (Pg/yr)"
  res@tiMainString               = "Russell et at 2018 - fig 9b "
  res@vpHeightF                  = 0.45
  res@vpWidthF                   = 0.6
  res@tmYLMode                   = "Automatic"
  res@tmYLMaxTicks               = 8
  res@tmXBMode                   = "Manual"
  res@tmXBMinorPerMajor          = 0
  res@tiYAxisFontHeightF         = 0.0125
  res@tiXAxisFontHeightF         = 0.0125
  res@tmXBLabelFontHeightF       = 0.01
  res@tmYLLabelFontHeightF       = 0.01
  res@tmYROn                     = False    ; Turn off top tickmarks
  res@tmXTOn                     = False    ; Turn off top tickmarks
  res@vpYF                       = 0.9
  res@vpXF                       = 0.08
  res@tmXMajorGrid               = True     ; Add vertical grid lines
  res@tmXMajorGridLineColor      = "grey"
  res@tmYMajorGrid               = True     ; Add horizontal grid lines
  res@tmYMajorGridLineColor      = "grey"
  res@tiXAxisString = "Latitudinal width of Southern Hemisphere Westerly Band"

  do iii = 0, nDatasets - 1

    fgco2_dataset = read_data(fgco2_items[iii])
    tauu_dataset = read_data(tauu_items[iii])
    infile_path = fgco2_inputfile_paths(iii)
    infile_iii  = addfile(infile_path, "r")

    if (iscoord(fgco2_dataset, "lat")) then
      ; some models have fgco2 on different grid than areacello
      var_lat = infile_iii->lat
      var_lon = infile_iii->lon
      radius_earth = 6.37e06
      deg2rad_convF = 0.0174533
      dlat = abs(var_lat(20) - var_lat(19))
      ; some models have closer lat points near poles.
      ; hence taking difference of 20th and 19th lat points
      dlon = abs(var_lon(20) - var_lon(19))
      dist_x_deg_earth = radius_earth * deg2rad_convF * dlat
      clat = cos(var_lat*deg2rad_convF)
      dx = dlon * radius_earth * deg2rad_convF * clat
      ; dx = radius of earth * cos(lat of this data point)
      ;     * (lon1 - lon2 (in radians))
      dy = dlat * radius_earth * deg2rad_convF
      ; dy = radius of earth *(lat1 - lat2 (in radians))
      dxdy = tofloat(dx*dy)     ; area of cell = dx * dy
      areacello_2d = new(dimsizes(fgco2_dataset), float)
      areacello_2d = conform(areacello_2d, dxdy, 0)
      delete(var_lon)
      delete(dx)
      delete(dy)
      delete(dxdy)
      delete(clat)
      delete(dlon)
      delete(dlat)
      fgco2_dataset = fgco2_dataset * 0.000031536
      ; unit conversion from kg/s to Pg/yr
      carbon_flux = fgco2_dataset * areacello_2d
      ; flux per cell   = (flux per area) * (area per cell)
      carbon_flux_per_lat = dim_sum_n_Wrap(carbon_flux, 1)
      ; flux per lat    = sum of flux per cell on lon dimension
      a = closest_val(-30.0, var_lat)
      carbon_flux_iii(iii) = sum(carbon_flux_per_lat(0:a))
      delete(var_lat)

    else
      fx_var = read_data(areacello_items[iii])
      if (all(ismissing(fx_var))) then
        error_msg("f", "russell18jgr-fig9b.ncl", " ", "areacello file of " + \
                  datasetnames(iii) + "not found in the recipe. If " + \
                  "areacello file available, please copy the dataset name " + \
                  "to additional dataset section of areacello.")
      end if
      areacello_2d = fx_var
      area_lat = infile_iii->lat
      fgco2_dataset = fgco2_dataset * 0.000031536
      ; unit conversion from kg/s to Pg/yr
      carbon_flux = fgco2_dataset * areacello_2d
      ; flux per cell   = (flux per area) * (area per cell)
      carbon_flux_per_lat = dim_sum_n_Wrap(carbon_flux, 1)
      ; flux per lat    = sum of flux per cell on lon dimension
      a = closest_val(-30.0, area_lat(:, 0))
      carbon_flux_iii(iii) = sum(carbon_flux_per_lat(0:a))
      delete(area_lat)
      delete(fx_var)
    end if

    tauu_lon = dim_avg_n_Wrap(tauu_dataset, 1)
    delete(tauu_dataset)

    if (iscoord(tauu_lon, "lat")) then   ; extracting lat values
      var_lat = tauu_lon&lat
    else
      infile_path = tauu_inputfile_paths(iii)
      infile_iii = addfile(infile_path, "r")
      area_lat = infile_iii->lat
      var_lat = dim_avg_n_Wrap(area_lat, 1)
      delete(area_lat)
    end if

    a1 = closest_val(-50.0, var_lat)
    a2 = closest_val(-75.0, var_lat)
    counter_lat = True
    do lat_ind = a1, (dimsizes(var_lat) / 2)
      ; loops over each lat from 50S equatorward till it becomes negative
      if (counter_lat) then
        if(tauu_lon(lat_ind).ge. 0 .and. tauu_lon(lat_ind+1).lt.0) then
          lat1 = var_lat(lat_ind)
          lat2 = var_lat(lat_ind+1)
          final_lat = lat1 - ((lat1 - lat2) * tauu_lon(lat_ind)) \
            / (tauu_lon(lat_ind)-tauu_lon(lat_ind+1))
          counter_lat = False
          delete(lat1)
          delete(lat2)
        end if
      end if
    end do
    do lat_ind1 = a2, a1
      ; loops over each lat from 50S poleward till tauu becomes negative
      if ((tauu_lon(lat_ind1) .lt. 0).and.(tauu_lon(lat_ind1+1).ge.0)) then
        lat1 = var_lat(lat_ind1)
        lat2 = var_lat(lat_ind1+1)
        lower_lat = lat1 - ((lat1 - lat2) * tauu_lon(lat_ind1) \
                            / (tauu_lon(lat_ind1) - tauu_lon(lat_ind1+1)))
        delete(lat1)
        delete(lat2)
      end if
    end do
    delete(var_lat)
    lat_width(iii) = final_lat - lower_lat
    delete(lower_lat)
    delete(final_lat)
    delete(infile_iii)
    delete(carbon_flux)
    delete(carbon_flux_per_lat)
    delete(areacello_2d)
    delete(fgco2_dataset)
    delete(tauu_lon)
  end do

  res@trYMaxF           =  decimalPlaces(max(carbon_flux_iii + 0.1), 1, True)
  res@trYMinF           =  decimalPlaces(min(carbon_flux_iii - 0.1), 1, True)
  xval_start_point      = (round(min(2*lat_width), 3)/2.0) - 0.5
  xval_end_point        = (round(max(2*lat_width), 3)/2.0) + 0.5
  res@trXMaxF           = xval_end_point
  res@trXMinF           = xval_start_point
  res@tmXBTickStartF    = xval_start_point
  res@tmXBTickSpacingF  = 0.5
  res@tmXBTickEndF      = xval_end_point
  res@xyDashPatterns    = 1        ; solid line
  res@xyLineThicknessF  = 1        ; set second line to 2

  aaa = fspan(xval_start_point, xval_end_point, nDatasets)
  ; array of x coordinates for line of best fit
  rc = regline(lat_width, carbon_flux_iii)
  ; calculates the slope and y-intercept for line of best fit
  linereg = rc@yintercept + (rc * aaa)
  ; calculates the y coordinates of line of best fit
  linereg = rc@yintercept + rc * aaa
  data1(0, :) = carbon_flux_iii
  data1(1, :) = linereg
  data1!0 = "i"
  data1&i = (/1, 2/)

  plot  = gsn_csm_xy(wks, aaa, linereg, res)       ; create plot
  mres                    = True
  mres@gsMarkerSizeF      = 0.0125         ; choose marker size

  do i = 0, dimsizes(lat_width) - 1
    mres@gsMarkerColor   = colors(i)       ; set marker colors
    mres@gsMarkerIndex   = markers(i)      ; choose marker types
    id = unique_string("mark")             ; create unique id
    ; add marker to plot
    plot@$id$ = gsn_add_polymarker(wks, plot, lat_width(i), \
                                   carbon_flux_iii(i), mres)
  end do

  draw(plot)    ; draw the plot

  legend = create "Legend" legendClass wks
  "vpXF"                     : 0.675                ; orientation on page
  "vpYF"                     : 0.925
  "vpWidthF"                 : 0.31                 ; width
  "vpHeightF"                : 0.725                ; height
  "lgPerimOn"                : False                ; no perimeter
  "lgItemCount"              : dimsizes(annots)     ; how many
  "lgLabelStrings"           : annots               ; labels
  "lgLabelsOn"               : True                 ; no default lables
  "lgLabelFontHeightF"       : 0.001                ; font height
  "lgItemType"               : "markers"            ; line paterns
  "lgMarkerColors"           : colors
  "lgMarkerIndexes"          : markers              ; one label color
  end create
  draw(legend)

  do idd = 0, dimsizes(lat_width) - 1
    nc_filename = config_user_info@work_dir + "russell18jgr_fig9b_" \
                  + annots(idd) + "_" + (start_years_data(idd)) + "-" \
                  + (end_years_data(idd)) + ".nc"
    outvar = data1(:, idd)
    outvar@var   = "carbon-flux_lat-width"
    outvar@diag_script = "russell18jgr_fig9b.ncl"
    outvar@model_name = annots(idd)
    outvar@regline_y_coord = linereg
    outvar@regline_x_coord = aaa
    outvar@description = "Total carbon flux and lat width of southern " \
      + "westerly band for dataset : " + annots(idd)  \
      + " for years " + (start_years_data(idd)) + "-" \
      + (end_years_data(idd)) + " are in 0th and 1st "\
      + "dimension respectively. Line of best fit's " \
      + "coordinates are added as attributes regline_coords"

    ncdf_outfile = ncdf_write(outvar, nc_filename)

    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 9b", \
                   "mean", \
                   "sh", \
                   "scatter", \
                   "russell_joellen", \
                   "russell18jgr", \
                   (/tauu_inputfile_paths(idd), fgco2_inputfile_paths(idd)/))

  end do
end
