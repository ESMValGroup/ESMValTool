; #############################################################################
;  russell_figure7i.ncl
;
;  Based on Figure 7i - Russell, J.L.,et al., 2018, J. Geophysical Research –
;    Oceans, 123, 3120-3143. https://doi.org/10.1002/2017JC013461 (figure 7i)
;
; Author:  Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
; ESMVal project
; #############################################################################
; Description
;
;   - Uses original grid (no regridding).
;   - Maskes the fgco2 data points on the land using ESMValTool preprocessor.
;   - if var coordinates are lat-lon then the script calculates the areacello
;       otherwise, the script reads the areacello variable from input file.
;   - Multiplies time averaged fgco2 with areacello to get flux per cell
;   - flux per lat = sum of all cells across a lat
;   - integrated flux at lat = cumulative sum of flux per lat from pole
;       till that lat.
;   - Plots integrated flux of fgco2 vs lat as xy line plot
;
; Required diag_script_info attributes (diagnostics specific)
;
;   - styleset : CMIP5 - default
;
;  Caveats:
;
; Modification history
;
;  20190510 - russell_joellen, pandde_amarjiit - written and
;                    implemented for ESMValTool v2.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  fgco2_items = select_metadata_by_name(input_file_info, "fgco2")
  datasetnames = metadata_att_as_array(fgco2_items, "dataset")
  start_years_data = metadata_att_as_array(fgco2_items, "start_year")
  end_years_data = metadata_att_as_array(fgco2_items, "end_year")
  inputfile_paths = metadata_att_as_array(fgco2_items, "filename")
  areacello_items = select_metadata_by_name(input_file_info, "areacello")
  areadatasets = metadata_att_as_array(areacello_items, "dataset")
  nDatasets = ListCount(fgco2_items)

end

begin

  ; Set annotations
  annots = project_style(input_file_info, diag_script_info, "annots")
  colors = project_style(input_file_info, diag_script_info, "colors")
  dashes = project_style(input_file_info, diag_script_info, "dashes")

  var0 = variable_info[0]@short_name

  plotpath = config_user_info@plot_dir + "/Russell_figure7i_" + var0 + "_" \
    + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
    + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  plot = new(nDatasets, graphic)

    res = True
    res@tmXBLabelFontHeightF       = 0.008
    res@tmYLLabelFontHeightF       = 0.008
    res@gsnFrame                   = False
    res@gsnDraw                    = False
    res@trXMinF                    = -80.
    res@trXMaxF                    = -30.
    res@trYMaxF                    = 1.
    res@trYMinF                    = -1.2
    res@vpHeightF                  = 0.65
    res@vpWidthF                   = 0.65
    res@tiMainString               = "Integrated Flux"
    res@gsnYRefLine                = 0
    res@gsnYRefLineDashPattern     = 2
    res@gsnYRefLineColor           = "grey"
    res@vpYF                       = 0.9
    res@vpXF                       = 0.08
    res@gsnRightString             = "Units - ( PgC/yr )"
    res@gsnRightStringFontHeightF  = 13.
    res@gsnLeftStringFontHeightF   = 13.
    res@tmXBMinorPerMajor          = 0
    res@tmXBTickStartF             = -80.
    res@tmXBTickSpacingF           = 5.
    res@tmXBTickEndF               = -30.
    res@tmXBMode                   = "Manual"
    res@tmYLMinorPerMajor          = 1
    res@tmYLTickStartF             = -1.2
    res@tmYLTickSpacingF           = 0.2
    res@tmYLTickEndF               = 1.0
    res@tmYLMode                   = "Manual"
    res@tiYAxisFontHeightF         = 0.0175
    res@tiYAxisOffsetXF            = 0.01

  do iii = 0, nDatasets - 1
    dataset = read_data(fgco2_items[iii])
    infile_path = inputfile_paths(iii)
    infile_iii  = addfile(infile_path, "r")

    if (iscoord(dataset, "lat")) then
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
      ; dx = radius of earth * cos(lat of this data point) *
      ;   (lon1 - lon2 (in radians))
      dy = dlat * radius_earth * deg2rad_convF
      ; dy = radius of earth *(lat1 - lat2 (in radians))
      dxdy = tofloat(dx*dy)   ; area of cell = dx * dy
      areacello_2d = new(dimsizes(dataset), float)
      areacello_2d = conform(areacello_2d, dxdy, 0)
      delete(var_lon)
      delete(dx)
      delete(dy)
      delete(dxdy)
      delete(clat)
      delete(dlon)
      delete(dlat)
      dataset = dataset * -0.000031536
      ; unit conversion from kg/s to Pg/yr
      carbon_flux = dataset * areacello_2d
      ; flux per cell   = (flux per area) * (area per cell)
      carbon_flux_per_lat = dim_sum_n_Wrap(carbon_flux, 1)
      ; flux per lat    = sum of flux per cell on lon dimension
      int_carbon_flux = cumsum(carbon_flux_per_lat, 2)
      ; integrated flux = cumulative sum of flux per lat
      int_carbon_flux!0 = "lat"
      int_carbon_flux&lat = var_lat
      delete(var_lat)
    else
      dataset_compare = datasetnames(iii) + "$"
      fx_ind = str_match_ind_regex(areadatasets, dataset_compare)
      fx_var = read_data(areacello_items[fx_ind])

      if (all(ismissing(fx_var))) then
        fx_variable = "areacello"
        error_msg("f", "russell_fig-7i.ncl", " ", "areacello file for " + \
                  datasetnames(iii) + " not found in the metadata file," + \
                  " please add 'fx_files: [areacello]' to the variable " + \
                  "dictionary in the recipe or add the location of " + \
                  " file to config-user.yml")
      end if
      areacello_2d = fx_var
      delete(fx_var)
      dataset = dataset * -0.000031536
      ; unit conversion from kg/s to Pg/yr
      carbon_flux = dataset * areacello_2d
      ; flux per cell   = (flux per area) * (area per cell)
      carbon_flux_per_lat = dim_sum_n_Wrap(carbon_flux, 1)
      ; flux per lat    = sum of flux per cell on lon dimension
      int_carbon_flux = cumsum(carbon_flux_per_lat, 2)
      ; integrated flux = cumulative sum of flux per lat
      int_carbon_flux!0 = "lat"
      area_lat = infile_iii->lat
      int_carbon_flux&lat = tofloat(area_lat(:, 0))
      delete(area_lat)

    end if
    delete(infile_iii)
    delete(carbon_flux)
    delete(carbon_flux_per_lat)
    delete(areacello_2d)
    delete(dataset)

    int_carbon_flux@var = var0
    int_carbon_flux@diag_script = "russell18jgr-fig7i.ncl"
    res@xyDashPatterns              = dashes(iii)
    res@xyLineColors                = colors(iii)
    res@xyExplicitLegendLabels      = annots(iii)
    res@gsnLeftString              = "Russell et al -2018 - Figure 7i "

    plot(iii) = gsn_csm_xy(wks, int_carbon_flux&lat, int_carbon_flux, res)

    nc_filename = config_user_info@work_dir + "russell_figure-7i_" + var0 + \
      "_" + annots(iii) + "_" + (start_years_data(iii)) + "-" + \
      (end_years_data(iii)) + ".nc"

    ncdf_outfile = ncdf_write(int_carbon_flux, nc_filename)
    delete(int_carbon_flux)
    if (iii .ne. 0) then
      overlay(plot(0), plot(iii))
    end if

    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 7i", \
                   "mean", \
                   "sh", \
                   "zonal", \
                   "russell_joellen", \
                   "russell18jgr", \
                   infile_path)

  end do
  draw(plot(0))

  legend = create "Legend" legendClass wks
  "vpXF"                     : 0.625                  ; orientation on page
  "vpYF"                     : 0.925
  "vpWidthF"                 : 0.5                    ; width
  "vpHeightF"                : 0.725                  ; height
  "lgPerimOn"                : False                  ; no perimeter
  "lgItemCount"              : dimsizes(annots)       ; how many
  "lgLineLabelStrings"       : annots                 ; labels
  "lgLabelsOn"               : False                  ; no default lables
  "lgLineLabelFontHeightF"   : 0.0085                 ; font height
  "lgDashIndexes"            : dashes                 ; line paterns
  "lgLineColors"             : colors
  "lgMonoLineLabelFontColor" : True                   ; one label color
  end create
  draw(legend)

end
