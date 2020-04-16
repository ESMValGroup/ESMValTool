; #############################################################################
; russell18jgr-fig5g.ncl
;
;  Russell, J.L.,et al., 2018, J. Geophysical Research – Oceans, 123,
;   3120-3143.  https://doi.org/10.1002/2017JC013461 (figure 5g)
;
; Author:  Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;
; #############################################################################
; Description
;   - Uses original grid (no regridding).
;   - Calculates the monthly climatology values of sic
;   - if var coordinates are lat-lon then the script calculates the areacello
;       otherwise, the script reads the areacello variable from input file.
;   - Multiplies areacello and sic/100 to find the total area of each
;       cell covered by ice
;   - Adds all cells below equator to get total sea ice area for each month
;   - Plots total area in 10^12 m^2 vs months
;   - overlays all the plots on the first one
;   - Draws the legend
;
;  Areacello is calculated for lat-lon model as many lat-lon models
;          have atmo grid for sic
;
; Required diag_script_info attributes (diagnostics specific)
;    (Do not change)
;
;    styleset : CMIP5 - Default
;
; Optional diag_script_info attributes (diagnostic specific)
;
;
; Caveats
;
;   - no caveats known on May 10, 2019
;
; Modification history
;
;  20190510 - russell_joellen, pandde_amarjiit - written and
;              implemented for ESMValTool v2.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function
load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")
  att = True
  att@mip = "OImon"
  info = select_metadata_by_atts(input_file_info, att)  ; variable
  var0 = info[0]@short_name

  info_items = select_metadata_by_name(input_file_info, var0)
  sic_datasets = metadata_att_as_array(info_items, "dataset")
  inputfile_paths = metadata_att_as_array(info_items, "filename")
  start_years_data = metadata_att_as_array(info_items, "start_year")
  end_years_data = metadata_att_as_array(info_items, "end_year")
  nDatasets = ListCount(info_items)

  att@mip = "fx"
  volinfo = select_metadata_by_atts(input_file_info, att)  ; area
  voldatasets = metadata_att_as_array(volinfo, "dataset")
  volfile_paths = metadata_att_as_array(volinfo, "filename")
  delete(att)
  if (dimsizes(sic_datasets) .ne. dimsizes(voldatasets)) then
    error_msg("f", DIAG_SCRIPT, " ", "areacello files for " + \
              "russell18jgr-fig5g.ncl do not match with sic datasets. " + \
              "Please check to make sure all the sic files have a " + \
              "areacello file in the recipe.")
  end if
end

begin

  ; Set annotations
  annots = project_style(input_file_info, diag_script_info, "annots")
  colors = project_style(input_file_info, diag_script_info, "colors")
  dashes = project_style(input_file_info, diag_script_info, "dashes")
  thicks = project_style(input_file_info, diag_script_info, "thicks")

  plotpath = config_user_info@plot_dir + "/russell18jgr-fig5g_" + var0 + "_" \
    + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
    + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  wks@fullname = plotpath
  plots = new(nDatasets, graphic)

  monthss = ispan(0, 11, 1)
  areacello_lat = False      ; counter for warning

  res                           = True
  res@tmXBLabelFontHeightF      = 0.008
  res@tmYLLabelFontHeightF      = 0.008
  res@gsnFrame                  = False
  res@gsnDraw                   = False
  res@trYMaxF                   = 24
  res@trYMinF                   = 0.
  res@tmXBMode                  = "Explicit"
  res@tmXBValues                = (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11/)
  res@vpHeightF                 = 0.65
  res@vpWidthF                  = 0.65
  res@gsnYRefLine               = 0
  res@gsnYRefLineDashPattern    = 2
  res@gsnYRefLineColor          = "grey"
  res@xyMarkLineMode            = "MarkLines"
  res@xyMarkers                 = 4
  res@vpYF                      = 0.9
  res@vpXF                      = 0.08
  res@gsnRightStringFontHeightF = 15.
  res@tiXAxisString             = "months"
  res@gsnLeftString             = "Russell et al -2018 - Figure 5 g"
  res@gsnLeftStringFontHeightF  = 17.
  res@tiXAxisFontHeightF        = 0.0175
  res@tiYAxisFontHeightF        = 0.0175
  res@tiYAxisOffsetXF           = 0.01
  res@tiYAxisString        = "Area under sea ice ( 10~S~12 ~N~ m~S~2~N~ )"
  res@tmXBLabels  = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", \
                     "Aug", "Sep", "Oct", "Nov", "Dec"/)

  do iii = 0, nDatasets - 1
    dataset = read_data(info_items[iii])  ; reading data
    sic12m = clmMonTLL(dataset)           ; create monthly climatology
    delete(dataset)
    if (max(sic12m) .lt. 5.0) then
        sic12m = sic12m * 100.0
        ; if variable is not in percentage
    end if
    dim_sic_12m = dimsizes(sic12m)
    a1 = round((dim_sic_12m(1) /2), 3)
    ; finding the index of grid point half of dimension sizes
    delete(dim_sic_12m)
    infile_path = inputfile_paths(iii)
    infile_iii  = addfile(infile_path, "r")

    if (iscoord(sic12m, "lat")) then
      areacello_lat = True
      var_lon = infile_iii->lon
      var_lat = infile_iii->lat
      radius_earth = 6.37e06
      deg2rad_convF = 0.0174533
      dlat = abs(var_lat(20) - var_lat(19))
      ; some models have closer lat points near poles.
      dlon = abs(var_lon(20) - var_lon(19))
      dist_x_deg_earth = radius_earth * deg2rad_convF * dlat
      clat = cos(var_lat*deg2rad_convF)

      dx = dlon * radius_earth * deg2rad_convF * clat
      ; dx = radius of earth * cos(lat of this data point) * (lon1 - lon2)
      dy = dlat * radius_earth * deg2rad_convF
      ; dy = radius of earth *(lat1 - lat2 (in radians))
      dxdy = tofloat(dx*dy)
      ; area of cell = dx * dy
      areacello_2d = new(dimsizes(sic12m(0, :, :)), float)
      areacello_2d = conform(areacello_2d, dxdy, 0)
      delete(var_lon)
      delete(dx)
      delete(dy)
      delete(dxdy)
      delete(clat)
      delete(dlon)
      delete(dlat)
      delete(var_lat)

    else
      dataset_compare = sic_datasets(iii) + "$"
      fx_ind = str_match_ind_regex(voldatasets, dataset_compare)
      fx_var = read_data(volinfo[fx_ind])
      if (all(ismissing(fx_var))) then
        error_msg("f", DIAG_SCRIPT, " ", "areacello file for " + \
                  annots(iii) + " not found, please add dataset " \
                  + "name in the additional dataset section of areacello" \
                  + " or remove the dataset from sic section.")
      end if

      areacello_2d = fx_var
      delete(fx_var)
    end if

    areacello_3d = new(dimsizes(sic12m), float)
    ; making areacello same dimensions as var12m
    areacello_3d = conform(areacello_3d, areacello_2d, (/1, 2/))
    delete(infile_iii)
    delete(areacello_2d)
    copy_VarCoords(sic12m, areacello_3d)
    sic12m = sic12m * areacello_3d
    delete(areacello_3d)
    sic_12m_lon_avg =  dim_sum_Wrap(sic12m(:, 0:a1, :))
    ; take total sum of southern hemisphere
    delete(sic12m)
    var_final = dim_sum_Wrap(sic_12m_lon_avg)  ; unit conversion
    delete(sic_12m_lon_avg)
    var_final = var_final / (10.0 ^ 14)
    ; 10^14 = 10^12 for unit conversion &
    ; 100 for converting sic from percentage to decimal

    res@xyDashPatterns            = dashes(iii)
    res@xyLineColors              = colors(iii)
    res@xyMarkerColors            = colors(iii)
    res@xyExplicitLegendLabels    = annots(iii)
    res@xyLineThicknessF          = thicks(iii)
    var_final@var = var0
    var_final@diag_script = DIAG_SCRIPT

    plots(iii) = gsn_csm_xy(wks, monthss, var_final, res)

    nc_filename = config_user_info@work_dir + "russell18jgr-fig5g_" + var0 \
      + "_" + annots(iii) + "_" + (start_years_data(iii)) + "-"  \
      + (end_years_data(iii)) + ".nc"

    ncdf_outfile = ncdf_write(var_final, nc_filename)
    delete(var_final)

    if (iii .ne. 0) then
        overlay(plots(0), plots(iii))
    end if

    ; Call provenance logger
    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 5g", \
                   "mean", \
                   "sh",   \
                   "times",  \
                   "russell_joellen", \
                   "russell18jgr", \
                   infile_path)

  end do

  draw(plots(0))

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

  if (areacello_lat) then  ; warning for lat lon areacello manual colection
    error_msg("w", DIAG_SCRIPT, " ", "All the models having lat-lon " + \
              "coordinates had a manual calculation of areacello this is " + \
              "not an error, no responce is needed by user for plotting ")
  end if

end
