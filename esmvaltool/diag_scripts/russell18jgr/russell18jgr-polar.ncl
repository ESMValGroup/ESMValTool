; ############################################################################
; russell_polcon_landmasked.ncl (polor contour plot with land values masked)
;
; Russell, J.L.,et al., 2018, J. Geophysical Research – Oceans, 123, 3120-3143
;  https://doi.org/10.1002/2017JC013461 (figures 1, 7a, & 8a)
;
; Author:  Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;          Pandde Amarjiit (University of Arizona, USA)
;
; adapted from sea ice polcon made by Daniel Senftleben (DLR, Germany)
;
; ###########################################################################
; Description
;   - Uses original grid (no regridding).
;   - Panels multiple plots: starts a new page after max_vert*max_hori.
;   - Uses the preprocessor land masked values of tauu, ph and fgco2.
;   - The script takes a time average of the masked data.
;   - Plots a polar contour plot.
;   - Adapted from the SouthernOcean_polcon.ncl & SeaIce Polar contor plot.
;
; Required diag_script_info attributes (configured to recreate figure 1, 7 & 8)
;   (Does not need changes)
;
;   - max_lat : plot ranges for SH.
;   - labelBar_end_type: determines the type of label bars.
;
; Optional diag_script_info attributes (configured to recreate figure 1, 7 & 8)
;   (Does not need changes)
;
;   - max_vert, max_hori: allowed dimensions of a panel page.
;                        (number of plots per column/row)
;   - grid_min, grid_max, grid_step: min, max and step size for color scale.
;   - colormap: NCL colour map name.
;   - colors:  array of rgb value as colormap.
;
; Required Preprocessor attributes ((configured to recreate figure 1, 7 & 8))
;   preprocessor_time_land: (for tauu)
;     climate_statistics:
;       operator: mean
;       period: full
;     mask_landsea:
;       mask_out: "land"
;
; Caveats
;  - fgco2 of CNRM-CM5 model doesnt work with ncdf_write function as it has
;      x,y coordinate variable instead of lat-lon or i-j or rlat-rlon
;
;
; Modification history
;
;  20190610 - russell_joellen, pandde_amarjiit - written and
;              implemented for ESMValTool v2
;
; #########################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function
load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  var_atts = variable_info[0]
  ; info_items = select_metadata_by_name(input_file_info, var0)
  var_name = metadata_att_as_array(input_file_info, "short_name")

  inputfile_paths = metadata_att_as_array(input_file_info, "filename")
  start_years_data = metadata_att_as_array(input_file_info, "start_year")
  end_years_data = metadata_att_as_array(input_file_info, "end_year")
  nDatasets = ListCount(input_file_info)

end

begin

  plotpath = config_user_info@plot_dir + "/Russell_polar-contour_" + var0 + \
     "_" + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
     + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  wks@fullname = plotpath
  plots = new(nDatasets, graphic)

  ; Set annotations
  annots = project_style(input_file_info, diag_script_info, "annots")
  ; Set contour parameters
  grid_min    =  0.   ; default
  grid_max    =  1.   ; default
  grid_step   =  0.1  ; default
  if (isatt(diag_script_info, "grid_min")) then
      grid_min = tofloat(diag_script_info@grid_min)
  end if
  if (isatt(diag_script_info, "grid_max")) then
      grid_max = tofloat(diag_script_info@grid_max)
  end if
  if (isatt(diag_script_info, "grid_step")) then
      grid_step = tofloat(diag_script_info@grid_step)
  end if
  nsteps = round((grid_max - grid_min) / grid_step, 3) + 1
  cnLevels = fspan(grid_min, grid_max, nsteps)
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

  ; Define color table
  if (isatt(diag_script_info, "colors")) then
    colors = diag_script_info@colors
    colors = colors / 256.0
  else
    if (isatt(diag_script_info, "colormap")) then
      colormap = diag_script_info@colormap
      gsn_define_colormap(wks, colormap)
      ColorStart = 2
      ColorEnd = -1
      if (colormap .eq. "BlWhRe") then
        ColorEnd = 95
        ColorStart = 7
      end if
    else
      colormap = "nrl_sirkes_nowhite"
      gsn_define_colormap(wks, colormap)
      ColorStart = 10
      ColorEnd = 240
    end if
  end if

  if (diag_script_info@labelBar_end_type .eq. "ExcludeOuterBoxes") then
    type_of_outer_boxes = "ExcludeOuterBoxes"
    outerBoxes_present = False
  else
    outerBoxes_present = True
    type_of_outer_boxes = "TriangleBothEnds"
  end if

  strUnits = " (" + var_atts@units + ")"  ; to show units in title
  if (isatt(diag_script_info, "new_units")) then
    strUnits = " (" + diag_script_info@new_units + ") "
  end if

  res                             = True
  res@lbLabelBarOn                = True         ; individual bars
  res@lbOrientation               = "Vertical"

  if (outerBoxes_present) then
    res@lbBoxEndCapStyle          = type_of_outer_boxes
  else
    res@cnLabelBarEndStyle        = type_of_outer_boxes
  end if

  if (isatt(diag_script_info, "grid_color")) then
    res@mpGridLineColor           = diag_script_info@grid_color
  else
    res@mpGridLineColor           = "green"
  end if

  res@cnLevelSelectionMode        = "ExplicitLevels"
  res@cnLevels                    = cnLevels
  res@cnFillOn                    = True       ; turn on color fill

  if (isatt(diag_script_info, "colors")) then
    res@cnFillPalette            = colors
  else
    res@gsnSpreadColors          = True         ; use full colormap
    res@gsnSpreadColorStart      = ColorStart   ; color to start
    res@gsnSpreadColorEnd        = ColorEnd     ; color to end
  end if

  res@mpGridLineDashPattern        = 0
  res@cnLinesOn                    = True       ; turn on color lines
  res@cnLineLabelsOn               = True       ; turn on contour line labels
  res@gsnAddCyclic                 = True       ; to make the plot complete
  res@mpFillOn                     = True
  res@mpLandFillColor              = (/128, 128, 128/) / 256.0
  res@gsnDraw                      = False
  res@gsnFrame                     = False
  res@gsnPolar                     = "SH"
  res@mpMaxLatF                    = diag_script_info@max_lat  ; Default is 0
  res@gsnContourZeroLineThicknessF = 2.
  res@gsnContourNegLineDashPattern = 2
  res@gsnRightStringOrthogonalPosF = 0.11
  res@mpFillDrawOrder          = "PostDraw"  ; to cover antartica with grey

  do iii = 0, nDatasets - 1
    dataset = read_data(input_file_info[iii])  ; reading data
    dataset@var = var_name(iii)
    dataset@diag_script = DIAG_SCRIPT
    if (isatt(diag_script_info, "unitCorrectionalFactor")) then
      ; to convert the units in datasets
      correction_factor = diag_script_info@unitCorrectionalFactor
      dataset = dataset * correction_factor
      delete(correction_factor)
    end if

    if(.not.(iscoord(dataset, "lat"))) then

      infile_path = inputfile_paths(iii)
      infile_iii  = addfile(infile_path, "r")
      area_lon = infile_iii->lon
      area_lat = infile_iii->lat
      dataset@lat2d = area_lat
      dataset@lon2d = area_lon
      delete(area_lat)
      delete(area_lon)
    end if

    range = start_years_data(iii) + " - " + end_years_data(iii)
    res@tiMainString                 = annots(iii)
    res@gsnRightString               = "annual mean " + range
    res@gsnLeftString                = var_name(iii) + strUnits

    plots(iii) = gsn_csm_contour_map_polar(wks, dataset, res)
    nc_filename = config_user_info@work_dir + "russell18jgr_polar_" \
      + var_name(iii) + "_" + annots(iii) + "_" + (start_years_data(iii)) \
      + "-" + (end_years_data(iii)) + ".nc"

    ncdf_outfile = ncdf_write(dataset, nc_filename)
    delete(dataset)
    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 polar plot " + var0, \
                   "mean", \
                   "sh", \
                   "geo", \
                   "russell_joellen", \
                   "russell18jgr", \
                   inputfile_paths(iii))
  end do

    ; Draw the panel
    pres = True
    pres@gsnPanelLabelBar = False
    outfile = panelling(wks, plots, nvert, nhori, pres)
end
