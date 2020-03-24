; ############################################################################
;  russell18jgr_fig2.ncl
;  Based on Figure 2 - Russell, J.L.,et al., 2018, J. Geophysical Research –
;    Oceans, 123, 3120-3143. https://doi.org/10.1002/2017JC013461 (figure 2)
;
;  Author: Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;
; ############################################################################
; Description
;    - Uses original grid (no regridding).
;    - Uses ESMValTool land masking (preprocessor) to mask land data points.
;    - Calculates longitudinal averaged tauu/tauuo values
;    - Plots averaged tauu values vs lat as xy line plot
;
;  Note - The figure in russell18jgr paper was made using tauuo but here
;         we can use both tauu and tauuo, it can be changed in the recipe.
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
;   preprocessor_time_land: (for tauu)
;    climate_statistics:
;      operator: mean
;      period: full
;    mask_landsea:
;      mask_out: "land"
;   preprocessor_time - for tauuo
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
;    - uses the same color, dash and thichkness for datasets of same models
;       (for ex. GFDL-ESM2M tauu and tauuo datasets use the same line)
;    - not tested for any observation models.
;
; Modification history
;  20190510 - A_russell_joellen, A_pandde_amarjiit - written and implemented
;             for ESMValToolv2
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  start_years_data = metadata_att_as_array(input_file_info, "start_year")
  end_years_data = metadata_att_as_array(input_file_info, "end_year")
  inputfile_paths = metadata_att_as_array(input_file_info, "filename")
  var_name = metadata_att_as_array(input_file_info, "short_name")
  nDatasets = ListCount(input_file_info)

end

begin

  ; Set annotations
  annots = project_style(input_file_info, diag_script_info, "annots")
  colors = project_style(input_file_info, diag_script_info, "colors")
  dashes = project_style(input_file_info, diag_script_info, "dashes")
  thicks = project_style(input_file_info, diag_script_info, "thicks")
  plotpath = config_user_info@plot_dir + "/russell18jgr_fig2_"  \
    + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
    + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  plot = new(nDatasets, graphic)
  res                             = True
  res@tmXBLabelFontHeightF        = 0.008
  res@tmYLLabelFontHeightF        = 0.008
  res@gsnFrame                    = False
  res@gsnDraw                     = False
  res@trXMinF                     = -80.
  res@trXMaxF                     = -30.
  res@trYMaxF                     = 0.25
  res@trYMinF                     = -0.1
  res@vpHeightF                   = 0.6
  res@vpWidthF                    = 0.6
  res@tiMainString                = " "
  res@gsnYRefLine                 = 0
  res@gsnYRefLineDashPattern      = 2
  res@gsnYRefLineColor            = "grey"
  res@vpYF                        = 0.85
  res@vpXF                        = 0.08
  res@gsnRightString              = "Units - (Pa)"
  res@gsnRightStringFontHeightF   = 15.
  res@gsnLeftString               = "Russell et al -2018 - Figure 2"
  res@gsnLeftStringFontHeightF    = 17.
  res@tmYLMinorPerMajor           = 4
  res@tmYLTickStartF              = -0.1
  res@tmYLTickSpacingF            = 0.1
  res@tmYLTickEndF                = 0.25
  res@tmYLMode                    = "Manual"
  res@tmXBMinorPerMajor           = 1
  res@tmXBTickStartF              = -80.
  res@tmXBTickSpacingF            = 5.
  res@tmXBTickEndF                = -30.
  res@tmXBMode                    = "Manual"
  res@tiYAxisFontHeightF          = 0.013
  res@tiXAxisFontHeightF          = 0.013
  res@tiYAxisOffsetXF             = 0.0
  res@tiYAxisString              = "Surface eastward wind stress"
  res@tiXAxisString              = "Latitude"

  do iii = 0, nDatasets - 1

    infile_path = inputfile_paths(iii)
    dataset = read_data(input_file_info[iii])        ; reading data
    var_lon_avg = dim_avg_n_Wrap(dataset, 1)         ; taking zonal average
    var_lon_avg@var = var_name(iii)
    var_lon_avg@diag_script = "russell18jgr-fig2.ncl"
    res@xyDashPatterns              = dashes(iii)
    res@xyLineColors                = colors(iii)
    res@xyExplicitLegendLabels      = annots(iii)
    res@xyLineThicknessF            = thicks(iii)

    if(.not.(iscoord(dataset, "lat"))) then
      ; extracting lat array from file for plotting
      infile_iii  = addfile(infile_path, "r")
      area_lat = infile_iii->lat
      var_lat = area_lat(:, 10)
      delete(area_lat)
    end if

    if (iscoord(var_lon_avg, "lat")) then
      plot(iii) = gsn_csm_xy(wks, var_lon_avg&lat, var_lon_avg, res)
    else
      plot(iii) = gsn_csm_xy(wks, var_lat, var_lon_avg, res)
      delete(var_lat)
    end if
    nc_filename = config_user_info@work_dir + "russell18jgr_fig2_" \
      + var_name(iii) + "_" + annots(iii) + "_" \
      + (start_years_data(iii)) + "-" + (end_years_data(iii))+".nc"

    ncdf_outfile = ncdf_write(var_lon_avg, nc_filename)
    delete(var_lon_avg)
    delete(dataset)

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 plotpath, \
                 "Russell et al 2018 figure 2", \
                 "mean", \
                 "sh", \
                 "geo", \
                 "russell_joellen", \
                 "russell18jgr", \
                 infile_path)

    if (iii .ne. 0) then
        overlay(plot(0), plot(iii))
    end if

  end do

  draw(plot(0))

  legend = create "Legend" legendClass wks
  "vpXF"                     : 0.62                 ; orientation on page
  "vpYF"                     : 0.925
  "vpWidthF"                 : 0.5                 ; width
  "vpHeightF"                : 0.725                ; height
  "lgPerimOn"                : False                ; no perimeter
  "lgItemCount"              : dimsizes(annots)     ; how many
  "lgLineLabelStrings"       : annots               ; labels
  "lgLabelsOn"               : False                ; no default lables
  "lgLineLabelFontHeightF"   : 0.00825              ; font height
  "lgDashIndexes"            : dashes               ; line paterns
  "lgLineColors"             : colors
  "lgLeftMarginF"            : -0.09
  "lgMonoLineLabelFontColor" : True                 ; one label color
  "lgLineDashSegLenF"        : 0.11
  end create
  draw(legend)
end
