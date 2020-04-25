; #############################################################################
;  russell18jgr_fig5.ncl
;
;  Russell, J.L.,et al., 2018, J. Geophysical Research – Oceans, 123,
;     3120-3143. https://doi.org/10.1002/2017JC013461 (figure 5)
;
; Author:   Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;          Pandde Amarjiit (University of Arizona, USA)
; ESMVal project
; #############################################################################
; Description
;  - Uses original grid (no regridding).
;  - Creates monthly climatology of sic data
;  - replaces grid cells with more than 15% sic
;         with 1 for September & 2 for March
;  - Plots March & September extent 15% or higher of sic
;         as red and blue respectively
;  - Panels multiple plots: starts a new page after max_vert*max_hori
;
; Required Preprocessor attributes ( no_preprocessor )
;  - None (no preprocessing required)
;     in the recipe do not keep a preprocessor in variable section
;
; Required diag_script_info attributes (configured for russell figure 5)
;         (Do not change)
;  - max_lat : plot ranges for SH
;
; Optional diag_script_info attributes (diagnostic specific)
;  - max_vert, max_hori: allowed dimensions of a panel page
;          (number of plots per column/row)
;
; Caveats
;  - BNU-ESM appears to have a rotated grid
;
; Modification history
;  20190510- russell_joellen, pandde_amarjiit - written and
;                    implemented for ESMValTool v2
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function
load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  var_atts = variable_info[0]
  info_items = select_metadata_by_name(input_file_info, var0)
  inputfile_paths = metadata_att_as_array(info_items, "filename")
  start_years_data = metadata_att_as_array(info_items, "start_year")
  end_years_data = metadata_att_as_array(info_items, "end_year")
  nDatasets = ListCount(info_items)
end

begin

  plotpath = config_user_info@plot_dir + "/russell18jgr-fig5_" + var0 + "_" \
            + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
            + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  wks@fullname = plotpath
  plots = new(nDatasets, graphic)

  ; Set annotations
  annots = project_style(input_file_info, diag_script_info, "annots")
  ; panel parameters
  nvert = 1  ; default
  nhori = 1  ; default
  if (isatt(diag_script_info, "max_vert")) then
      nvert = toint(diag_script_info@max_vert)
      ; Maximum allowed number of plots per page (vertical)
  end if

  if (isatt(diag_script_info, "max_hori")) then
      nhori = toint(diag_script_info@max_hori)
      ; Maximum allowed number of plots per page (horizontal)
  end if

  cnLevels = (/ 1.01, 2.001/)
  colorshere = (/ "blue4", "red"/)

  res                              = True
  res@lbLabelBarOn                 = False    ; individual bars
  res@cnLevelSelectionMode         = "ExplicitLevels"
  res@cnLevels                     = cnLevels
  res@cnFillPalette                = colorshere
  res@mpGridLatSpacingF            = 10
  res@mpFillDrawOrder              = "PostDraw"  ; covers antartica with grey
  res@mpLandFillColor              =  (/128, 128, 128/) / 256.0
  res@cnLinesOn                    = False
  res@cnLineLabelsOn               = False
  res@cnFillOn                     = True
  res@gsnCenterString              = "(Blue - September & Red - March)"
  res@gsnDraw                      = False       ; Draw and Frame is
  res@gsnFrame                     = False     ; called in panelling
  res@gsnLeftString                = "Southern Ocean Max Min Sea ice extent"
  res@gsnPolar                     = "SH"
  res@mpMaxLatF                    = diag_script_info@max_lat  ; Default is 0
  res@gsnAddCyclic                 = True
  res@gsnRightStringOrthogonalPosF = 0.14
  res@gsnLeftStringOrthogonalPosF  = 0.14
  res@gsnRightStringFontHeightF    = 13.
  res@gsnLeftStringFontHeightF     = 15.

  do iii = 0, nDatasets - 1

    dataset = read_data(info_items[iii])
    sic12m = clmMonTLL(dataset)  ; create monthly climatology
    delete(dataset)
    if (max(sic12m) .lt. 5.0) then
        sic12m = sic12m * 100.0
        ; if variable is not in percentage
    end if

    sic_march = sic12m(2, :, :)      ; average data for march month
    sic_september = sic12m(8, :, :)  ; average data for september month
    sic_september = where((sic_september .gt. 15), 1, sic12m@_FillValue)
    ; masking all the grid points which have sic less than 15%
    dataset = where((sic_march .gt. 15), 2, sic_september)
    delete(sic_september)
    copy_VarCoords(sic_march, dataset)
    delete(sic_march)
    infile_path = inputfile_paths(iii)
    infile_iii  = addfile(infile_path, "r")

    if(.not.(iscoord(sic12m, "lat"))) then

      area_lon = infile_iii->lon
      area_lat = infile_iii->lat
      dataset@lat2d = area_lat
      dataset@lon2d = area_lon
      delete(area_lat)
      delete(area_lon)
    else
      ; if lat lon arrays lost in where command
      var_lon = infile_iii->lon
      var_lat = infile_iii->lat
      dataset!0 = "lat"
      dataset!1 = "lon"
      dataset&lat = var_lat
      dataset&lon = var_lon
      delete(var_lon)
      delete(var_lat)
    end if

    delete(sic12m)
    range = start_years_data(iii) + " - " + end_years_data(iii)
    res@tiMainString                 = annots(iii)
    res@gsnRightString               = "annual mean " + range
    dataset@var = var0
    dataset@diag_script = "russell18jgr-fig5"

    plots(iii) = gsn_csm_contour_map_polar(wks, dataset, res)

    nc_filename = config_user_info@work_dir + "russell18jgr_fig5_" + var0  \
      + "_" + annots(iii) + "_" + (start_years_data(iii)) + "-" \
      + (end_years_data(iii)) + ".nc"

    ncdf_outfile = ncdf_write(dataset, nc_filename)
    delete(dataset)

    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 5 -polar", \
                   "mean", \
                   "sh", \
                   "geo", \
                   "russell_joellen", \
                   "russell18jgr", \
                   infile_path)

  end do

  ; Draw the panel
  pres = True
  pres@gsnPanelLabelBar = False
  outfile = panelling(wks, plots, nvert, nhori, pres)

end
