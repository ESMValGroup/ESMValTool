; #############################################################################
;  russell18jgr_fig7h.ncl
;
;  Based on Figure 7h - Russell, J.L.,et al., 2018, J. Geophysical Research –
;   Oceans, 123, 3120-3143. https://doi.org/10.1002/2017JC013461 (figure 7h)
;
; Author:  Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;
; #############################################################################
; Description
;
;  - Uses original grid (no regridding).
;  - Maskes the fgco2 data points on the land using ESMValTool preprocessor.
;  - Calculates time and longitudinal averaged fgco2 values.
;  - Plots averaged fgco2 values vs lat as xy line plot.
;
; Required diag_script_info attributes (diagnostics specific)
;
;    - styleset : CMIP5 - default
;
; Required Preprocessor attributes ( no_preprocessor )
;     - None (no preprocessing required)
;       but for esmvaltool just add (time_average: false)
;
;  Caveats:
;
;    1) i-j models are assigned lat - lon arrays assuming the grid has
;         constant lat across j and constant lon across i.
;    2) CNRM-CM5 model doesnt work with ncdf_write function as it has
;         x,y coordinate variable instead of lat-lon or i-j or rlat-rlon
;
;
; Modification history
;
; 20190510 -russell_joellen, pandde_amarjiit - written and
;    implemented for ESMValTool v2.
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

  plotpath = config_user_info@plot_dir + "/Russell_figure7h_" + var0 + "_" \
    + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
    + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  plot = new(nDatasets, graphic)

  res                           = True
  res@tmXBLabelFontHeightF      = 0.009
  res@tmYLLabelFontHeightF      = 0.009
  res@gsnFrame                  = False
  res@gsnDraw                   = False
  res@trXMinF                   = -80.
  res@trXMaxF                   = -30.
  res@trYMaxF                   =  40.
  res@trYMinF                   = -60.
  res@vpHeightF                 = 0.65
  res@vpWidthF                  = 0.65
  res@tiMainString              = "Zonal-mean Flux"
  res@gsnYRefLine               = 0
  res@gsnYRefLineDashPattern    = 2
  res@gsnYRefLineColor          = "grey"
  res@vpYF                      = 0.9
  res@vpXF                      = 0.08
  res@gsnRightString            = "Units - gC/ (m~S~2~N~ * yr)"
  res@gsnRightStringFontHeightF = 15.
  res@gsnLeftString             = "Russell et al -2018 - Figure 7 h"
  res@gsnLeftStringFontHeightF  = 17.
  res@tmYLMinorPerMajor         = 1
  res@tmYLTickStartF            = -60.
  res@tmYLTickSpacingF          = 10.
  res@tmYLTickEndF              = 40.
  res@tmYLMode                  = "Manual"
  res@tmXBMinorPerMajor         = 0
  res@tmXBTickStartF            = -80.
  res@tmXBTickSpacingF          = 5.
  res@tmXBTickEndF              = -30.
  res@tmXBMode                  = "Manual"
  res@tiYAxisFontHeightF        = 0.0175
  res@tiYAxisOffsetXF           = 0.01
  res@tiXAxisString             = "Latitude"

  do iii = 0, nDatasets - 1
    dataset = read_data(info_items[iii])
    infile_path = inputfile_paths(iii)

    if (iscoord(dataset, "i")) then
      ; to convert i-j coordinate variable to lat -lon
      delete(dataset&i)
      delete(dataset&j)
      infile_iii  = addfile(infile_path, "r")
      area_lon = infile_iii->lon
      area_lat = infile_iii->lat

      dataset!0   = "lat"
      dataset&lat = tofloat(area_lat(:, 0))
      dataset!1   = "lon"
      dataset&lon = tofloat(area_lon(0, :))
      delete(area_lat)
      delete(area_lon)
    end if
    dataset = dataset * -31536000000.
    ;  unit conversion factor from kg /(m^2 * sec) to g / (m^2 * yr)
    if(.not.(iscoord(dataset, "lat"))) then

      infile_iii  = addfile(infile_path, "r")
      area_lon = infile_iii->lon
      area_lat = infile_iii->lat
      var_lat = area_lat(:, 10)
      delete(area_lat)
      delete(area_lon)
    end if
    var_lon_avg = dim_avg_n_Wrap(dataset, 1)          ; taking zonal average
    var_lon_avg@var = var0
    var_lon_avg@diag_script = "russell18jgr_fig7h.ncl"
    res@xyDashPatterns              = dashes(iii)
    res@xyLineColors                = colors(iii)
    res@xyExplicitLegendLabels      = annots(iii)

    if (iscoord(var_lon_avg, "lat")) then
      plot(iii) = gsn_csm_xy(wks, var_lon_avg&lat, var_lon_avg, res)
    else
      plot(iii) = gsn_csm_xy(wks, var_lat, var_lon_avg, res)
      delete(var_lat)
    end if
    if (iii .ne. 0) then
      overlay(plot(0), plot(iii))
    end if

    nc_filename = config_user_info@work_dir + "russell18jgr_fig-7h_" + var0 \
      + "_" + datasetnames(iii) + "_" + (start_years_data(iii)) + "-" + \
      (end_years_data(iii)) + ".nc"

    ncdf_outfile = ncdf_write(var_lon_avg, nc_filename)
    delete(var_lon_avg)
    delete(dataset)
    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 7h", \
                   "mean", \
                   "sh", \
                   "zonal", \
                   "russell_joellen", \
                   "russell18jgr", \
                   infile_path)

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
