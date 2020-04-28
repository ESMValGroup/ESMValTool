; ############################################################################
;  russell18jgr_fig9a.ncl
;
; Russell, J.L.,et al., 2018, J. Geophysical Research – Oceans, 123, 3120-3143.
;  https://doi.org/10.1002/2017JC013461 (figure 9a)
;
; Author:   Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;          Pandde Amarjiit (University of Arizona, USA)
; ESMVal project
; #############################################################################
; Description
;
;  - Uses original grid (no regridding).
;     (for total heat flux part of script)
;  - Maskes the hfds data points on the land using ESMValTool preprocessor.
;  - if var coordinates are lat-lon then the script calculates the areacello
;      otherwise, the script reads the areacello variable from input file.
;  - Multiplies time averaged hfds with areacello to get heat flux per cell.
;  - flux per lat = sum of all cells across a lat.
;  - total heat flux in Southern ocean = sum of flux per lat from 30S - 90S.
;
;    (for Southern westerly band part of the script)
;  - Uses ESMValTool land masking (preprocessor) to mask land data points.
;  - Calculates longitudinal averaged tauu/tauuo values.
;  - Calculates the first lat north and south of 50S where tauu is zero.
;  - Subtracts the north lat value from south lat value to get
;      Latitudinal width of Southern Hemisphere Westerly Band.
;  - Calculates line of linear regression between lat width and
;      total heat flux.
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
;    - time_average
;
; Required diag_script_info attributes (diagnostics specific)
;
;    - styleset : CMIP5 - default
;    - ncdf     : default
;
; Caveats
;  - Does not work for models without hfds, these models will be made
;       compatible in future.
;
; Modification history
;
; 20190510 - russell_joellen, pandde_amarjiit - written and
;      implemented for ESMValTool v2.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function
load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  nVariables = ListCount(variable_info)
  hfds_items = select_metadata_by_name(input_file_info, "hfds")

  var_names = metadata_att_as_array(variable_info, "short_name")
  var_ind = str_match_ind(var_names, "tauu")
  tauu_var = var_names(var_ind)
  tauu_items = select_metadata_by_name(input_file_info, tauu_var)
  tauu_inputfile_paths = metadata_att_as_array(tauu_items, "filename")
  hfds_inputfile_paths = metadata_att_as_array(hfds_items, "filename")
  start_years_data = metadata_att_as_array(tauu_items, "start_year")
  end_years_data   = metadata_att_as_array(tauu_items, "end_year")
  nDatasets = ListCount(tauu_items)
  areacello_items = select_metadata_by_name(input_file_info, "areacello")
  areadatasets = metadata_att_as_array(areacello_items, "dataset")

  if (nVariables .ne. 3) then
    error_msg("f", "russell_fig9a.ncl", " ", "Number of variables for" + \
              " this diag script must be 3 (hfds, areacello & tauu/tauuo)" + \
              ". Please make sure one of tauu or tauuo is commented out")
  end if
  if (nDatasets .lt. 3) then
    error_msg("f", "russell_fig9a.ncl", " ", "Minimum number of " + \
              "datasets required for this diag script must be 3. " + \
              "Please add more datasets.")
  end if
  datasetnames = metadata_att_as_array(tauu_items, "dataset")
end

begin

  plotpath = config_user_info@plot_dir + "/russell18jgr-fig9a_" \
            + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
            + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  wks@fullname = plotpath
  plots = new(nDatasets, graphic)
  lat_width = new(nDatasets, double)
  heat_flux_iii = new(nDatasets, double)

  data1 = new((/2, dimsizes(datasetnames)/), double)
  aaa = new(dimsizes(datasetnames), double)

  ; Set annotations
  annots = project_style(tauu_items, diag_script_info, "annots")
  markers = project_style(tauu_items, diag_script_info, "markers")
  colors = project_style(tauu_items, diag_script_info, "colors")

  res                       = True
  res@gsnMaximize           = False
  res@gsnFrame              = False
  res@gsnDraw               = True
  res@tiYAxisString         = "Southern ocean heat uptake (PW)"
  res@tiMainString          = " Russell et al 2018 - Figure 9a "
  res@vpHeightF             = 0.45
  res@vpWidthF              = 0.65
  res@tmYLMinorPerMajor     = 0
  res@tmXBMinorPerMajor     = 0
  res@tmYLMode              = "Automatic"
  res@tmXBMode              = "Manual"
  res@tiYAxisFontHeightF    = 0.0125
  res@tiXAxisFontHeightF    = 0.0125
  res@tmXBLabelFontHeightF  = 0.01     ; fontsize of tickmark labels (x-axis)
  res@tmYLLabelFontHeightF  = 0.01     ; fontsize of tickmark labels (x-axis)
  res@tmYROn                = False    ; Turn off top tickmarks
  res@tmXTOn                = False    ; Turn off top tickmarks
  res@vpYF                  = 0.9
  res@vpXF                  = 0.1
  res@tmXMajorGrid          = True     ; Add vertical grid lines
  res@tmXMajorGridLineColor = "grey"
  res@tmYMajorGrid          = True     ; Add horizontal grid lines
  res@tmYMajorGridLineColor = "grey"
  res@tiXAxisString = "Latitudinal width of Southern Hemisphere Westerly Band"

  do iii = 0, nDatasets - 1

    hfds_dataset =  read_data(hfds_items[iii])
    infile_path = hfds_inputfile_paths(iii)
    infile_iii  = addfile(infile_path, "r")

    tauu_dataset = read_data(tauu_items[iii])
    tauu_lon = dim_avg_n_Wrap(tauu_dataset, 1)
    delete(tauu_dataset)

    fx_var = read_data(areacello_items[iii])
    ; checking for areacello
    if (all(ismissing(fx_var))) then
      error_msg("f", "russell18jgr-fig9a.ncl", " ", "areacello file of " + \
                datasetnames(iii) + "not found in the recipe. If areacello" + \
                " file available, please copy the dataset name to " + \
                "additional dataset section of areacello.")
    end if
    areacello_2d = fx_var

    if (iscoord(hfds_dataset, "lat")) then
      var_lat = hfds_dataset&lat
    else
      area_lat = infile_iii->lat
      var_lat  = area_lat(:, 0)
      delete(area_lat)
    end if
    ; unit conversion from W to PW
    heat_flux = hfds_dataset * areacello_2d / (10.0 ^ 15)
    ; flux per cell   = (flux per area) * (area per cell)
    heat_flux_per_lat = dim_sum_n_Wrap(heat_flux, 1)
    ; flux per lat    = sum of flux per cell on lon dimension
    a = closest_val(-30.0, var_lat)
    ; finds the closest lat grid point to 30S
    heat_flux_iii(iii) = sum(heat_flux_per_lat(0:a))
    delete(var_lat)
    delete(fx_var)

    if (iscoord(tauu_lon, "lat")) then  ; extracting lat values
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
    do lat_ind = a1, (dimsizes(var_lat)/2)
      ; loops over each lat from 50S equatorward till it becomes negative
      if (counter_lat) then
        if(tauu_lon(lat_ind).ge. 0 .and. tauu_lon(lat_ind+1).lt.0) then
          lat1 = var_lat(lat_ind)
          lat2 = var_lat(lat_ind+1)
          final_lat = lat1 - ((lat1 - lat2) * tauu_lon(lat_ind)) \
            / (tauu_lon(lat_ind)-tauu_lon(lat_ind+1))
          counter_lat = False
        end if
      end if
    end do
    delete(lat1)
    delete(lat2)
    do lat_ind1 = a2, a1
      ; loops over each lat from 50S poleward till tauu becomes negative
      if ((tauu_lon(lat_ind1) .lt. 0).and.(tauu_lon(lat_ind1+1).ge.0)) then
        lat1 = var_lat(lat_ind1)
        lat2 = var_lat(lat_ind1+1)
        lower_lat = lat1 - ((lat1 - lat2) * tauu_lon(lat_ind1) \
                            / (tauu_lon(lat_ind1) - tauu_lon(lat_ind1+1)))
      end if
    end do
    delete(lat1)
    delete(lat2)
    delete(var_lat)
    lat_width(iii) = final_lat - lower_lat
    delete(lower_lat)
    delete(final_lat)
    delete(infile_iii)
    delete(heat_flux)
    delete(heat_flux_per_lat)
    delete(areacello_2d)
    delete(hfds_dataset)
    delete(tauu_lon)

  end do

  res@trYMaxF      = decimalPlaces(max(heat_flux_iii+0.1), 1, True)
  res@trYMinF      = decimalPlaces(min(heat_flux_iii-0.1), 1, True)
  xval_start_point = (round(min(lat_width - 0.5), 3))
  xval_end_point   = (round(max(lat_width + 0.5), 3))

  res@tmYLMaxTicks = 10
  res@trXMaxF      = xval_end_point
  res@trXMinF      = xval_start_point
  res@tmXBTickStartF    = xval_start_point
  res@tmXBTickSpacingF  = 1.0
  res@tmXBTickEndF      = xval_end_point
  res@xyMarkLineMode    = "Lines"
  res@xyDashPatterns    = 1            ; solid line
  res@xyLineThicknessF  = 1

  aaa = fspan(xval_start_point, xval_end_point, nDatasets)
  ; array of x coordinates for line of best fit
  rc    = regline(lat_width, heat_flux_iii)
  ; calculates the slope and y-intercept for line of best fit
  linereg = rc@yintercept + (rc * aaa)
  ; calculates the y coordinates of line of best fit
  data1(0, :) = heat_flux_iii
  data1(1, :) = lat_width
  data1!0 = "i"
  data1&i = (/1, 2/)
  plot  = gsn_csm_xy(wks, aaa, linereg, res)   ; create plot

  mres                    = True
  mres@gsMarkerSizeF      = 0.0125             ; choose marker size

  do i = 0, dimsizes(lat_width) - 1
    ; overlaying the markers on the plot
    mres@gsMarkerColor   = colors(i)            ; set marker colors
    mres@gsMarkerIndex   = markers(i)           ; choose marker types
    id = unique_string("mark")                  ; create unique id
    ; add markers to plot
    plot@$id$ = gsn_add_polymarker(wks, plot, lat_width(i), \
                                   heat_flux_iii(i), mres)
  end do

  draw(plot)    ; draw the plot

  legend = create "Legend" legendClass wks
  "vpXF"                     : 0.675                  ; orientation on page
  "vpYF"                     : 0.925
  "vpWidthF"                 : 0.31                   ; width
  "vpHeightF"                : 0.725                  ; height
  "lgPerimOn"                : False                  ; no perimeter
  "lgItemCount"              : dimsizes(annots)       ; how many
  "lgLabelStrings"           : annots                 ; labels
  "lgLabelsOn"               : True                   ; no default lables
  "lgLabelFontHeightF"       : 0.001                  ; font height
  "lgItemType"               : "markers"              ; line paterns
  "lgMarkerColors"           : colors
  "lgMarkerIndexes"          : markers                ; one label color
  end create
  draw(legend)

  do idd = 0, dimsizes(lat_width) - 1
    nc_filename = config_user_info@work_dir + "russell18jgr_fig9a_" \
                  + annots(idd) + "_" + (start_years_data(idd)) + "-" \
                  + (end_years_data(idd)) + ".nc"
    outvar = data1(:, idd)
    outvar@var   = "heat-flux_lat-width"
    outvar@diag_script = "russell18jgr_fig9a.ncl"
    outvar@model_name = annots(idd)
    outvar@regline_y_coord = linereg
    outvar@regline_x_coord = aaa

    outvar@description = "Total heat flux and lat width of southern westerly" \
      + " band for dataset : " + annots(idd) + " for years " \
      + (start_years_data(idd)) + "-" + (end_years_data(idd)) \
      + " are in 0th and 1st dimension respectively. Line of best fit's " \
      + "coordinates are added as attributes regline_coords"

    ncdf_outfile = ncdf_write(outvar, nc_filename)

    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 9a", \
                   "mean", \
                   "sh", \
                   "scatter", \
                   "russell_joellen", \
                   "russell18jgr", \
                   (/tauu_inputfile_paths(idd), hfds_inputfile_paths(idd)/))

  end do
end
