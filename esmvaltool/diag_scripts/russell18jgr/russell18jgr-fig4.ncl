; ############################################################################
; russell18jgr_fig4.ncl
;
;  Russell, J.L.,et al., 2018, J. Geophysical Research - Oceans, 123,
;  3120-3143. https://doi.org/10.1002/2017JC013461 (figure 4)
;
; Author:  Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;
; ############################################################################
; Description
;   - Uses original grid (no regridding).
;   - Changes i-j and x-y coordinates to lat-lon coordinates
;   - Uses the time averaged data from the closes lon to 69W (291E)
;   - Converts units if unitCorrectionalFactor attribute is present
;   - Calculates total transport from volcello file
;   - If volcello file not found then skips directly to plotting
;   - Volcello(m^3) is divided by the horizontal distance between two data
;     points along the same latitude, which results in east-west
;     cross-section area of each grid cell
;   - total transport is sum of transport of each cell between 75S and 49S
;   - Panels all the plots together
;
; Required diag_script_info attributes
;
;   - styleset    : CMIP5         - default
;   - ncdf        : default
;
; Optional diag_script_info attributes
;
;   - max_vert, max_hori: allowed dimensions of a panel page
;                        (number of plots per column/row)
;
;   - unitCorrectionalFactor : 100 (to convert from m/s to cm/s)
;   - new_units              : "cm/s"
;
; Required Preprocessor attributes ((configured to recreate figure 4))
;   - None (no preprocessing required)
;     in the recipe do not keep a preprocessor in variable section
;
;
; Caveats
;
;   - MIROC-ESM and BNU-ESM does not work as depth variable is not called lev
;   - total transport is very small in case of NorESM1-M and ME as volcello
;     values look incorrect(very small)
;   - MRI models have 0 in land points instead of missing values.
;     So, no grey region produced on land.
;
; Modification history
;
;  20190610-russell_joellen, pandde_amarjiit: written and implemented for
;           ESMValTool v2
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/contour_maps.ncl"  ; load plot function
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  att = True
  att@mip = "Omon"
  info = select_metadata_by_atts(input_file_info, att)  ; variable
  var0 = info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  inputfile_paths = metadata_att_as_array(info_items, "filename")
  start_years_data = metadata_att_as_array(info_items, "start_year")
  end_years_data = metadata_att_as_array(info_items, "end_year")
  nDatasets = ListCount(info_items)
  annots = project_style(input_file_info, diag_script_info, "annots")
  uo_dataset_name = metadata_att_as_array(info_items, "dataset")

  att@mip = "fx"
  volinfo = select_metadata_by_atts(input_file_info, att)  ; area
  voldatasets = metadata_att_as_array(volinfo, "dataset")
  volfile_paths = metadata_att_as_array(volinfo, "filename")
  delete(att)

end

begin

  colors_new = (/ (/15.0, 69, 168/), (/36, 118, 205/), (/57, 162, 245/),   \
                  (/96, 190, 250/), (/131, 212, 253/), (/146, 230, 253/),  \
                  (/161, 241, 255/), (/188, 246, 255/), (/205, 226, 229/), \
                  (/234, 231, 211/), (/251, 246, 190/), (/255, 232, 154/), \
                  (/252, 224, 97/), (/254, 173, 26/), (/251, 136, 10/),    \
                  (/238, 91, 12/), (/209, 49, 7/), (/178, 0, 0/) /)
  colors_new = colors_new / 256.0

  plotpath = config_user_info@plot_dir + "/Russell18jgr-fig4_" \
    + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
    + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  file_type = output_type()
  wks = gsn_open_wks(file_type, plotpath)
  wks@fullname = plotpath
  plots = new(nDatasets, graphic)

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

  warn_ij = True

  res                           = True           ; plot mods desired
  res@gsnDraw                   = False          ; draw plot
  res@gsnFrame                  = False          ; advance frame
  res@cnFillOn                  = True           ; turn on color fill
  res@cnLinesOn                 = True           ; turn on contour lines
  res@cnLineLabelsOn            = True           ; turn on contour labels
  res@cnLineLabelsOn            = True           ; draw contour labels
  res@cnLineLabelFontHeightF    = 0.006
  res@cnMissingValFillColor     = "grey44"
  res@cnFillPalette             = colors_new
  res@cnLevelSelectionMode      = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF            = -40.0        ; set min contour level
  res@cnMaxLevelValF            = 40.0         ; set max contour level
  res@cnLevelSpacingF           = 5.0          ; set contour spacing
  res@lbLabelBarOn              = True         ; color bar
  res@lbOrientation             = "vertical"   ; vertical label bars
  res@pmLabelBarOrthogonalPosF  = -0.04        ; shift left-right
  res@lbBoxMinorExtentF         = 0.2          ; slim width
  res@tiXAxisString             = "Latitude"   ; title
  res@tiYAxisString             = "Depth (m)"
  res@gsnCenterString           = ""
  res@gsnLeftStringFontHeightF  = 13.
  res@gsnRightStringFontHeightF = 13.
  res@vpWidthF                  = 0.60         ; resize
  res@vpHeightF                 = 0.45
  res@trXMinF                   = -74.
  res@trXMaxF                   = -50.
  res@tmXBTickSpacingF          = 2.0          ; tick mark interval
  res@tmXBMinorPerMajor         = 1
  res@tmXBTickStartF            = -74.
  res@tmXBTickEndF              = -50.
  res@tmXBMode                  = "Manual"     ; adjust X axis
  res@gsnXAxisIrregular2Linear  = True         ; adjust Y axis
  res@gsnYAxisIrregular2Linear  = True
  res@trYReverse                = True         ; reverse Y axis
  res@tiMainFontHeightF         = 0.022        ; fontsize of the main title
  res@txFontHeightF             = 0.02         ; fontsize of the subtitles
  res@tmXBLabelFontHeightF      = 0.01         ; fontsize of tickmark labels
  res@tmYLLabelFontHeightF      = 0.01         ; fontsize of tickmark labels
  res@tiXAxisFontHeightF        = 0.015
  res@tiYAxisFontHeightF        = 0.017
  res@lbBoxEndCapStyle          = "TriangleBothEnds"
  res@gsnContourZeroLineThicknessF = 2      ; doubles thickness of zero contour
  res@gsnContourNegLineDashPattern = 1      ; dash pattern 1 for neg. contours

  do iii = 0, nDatasets - 1

    dataset_all_time = read_data(info_items[iii])
    dataset = dim_avg_n_Wrap(dataset_all_time, 0)
    delete(dataset_all_time)
    dataset@var = var0
    dataset@diag_script = DIAG_SCRIPT
    dataset@model_name = annots(iii)

    strUnits = " (" + dataset@units + ")"  ; to show units in title
    if (isatt(diag_script_info, "new_units")) then
      strUnits = " (" + diag_script_info@new_units + ") "
    end if

    dataset_compare = uo_dataset_name(iii) + "$"
    fx_ind = str_match_ind_regex(voldatasets, dataset_compare)
    if (ismissing(fx_ind)) then
      fx_var = new(1, float)
      fx_var = fx_var@_FillValue
    else
      fx_var = read_data(volinfo[fx_ind])
    end if

    if (iscoord(dataset, "i")) then  ; changes i-j coordinates to lat-lon

      delete(dataset&i)
      delete(dataset&j)
      infile_path = inputfile_paths(iii)
      infile_iii  = addfile(infile_path, "r")
      area_lon = infile_iii->lon
      area_lat = infile_iii->lat
      dataset!1  = "lat"
      dataset!2  = "lon"
      dataset&lat = tofloat(area_lat(:, 10))
      dataset&lon = tofloat(area_lon(10, :))
      delete(area_lat)
      delete(area_lon)
      warn_ij = True

    end if

    if(.not.(iscoord(dataset, "lat"))) then
      infile_path = inputfile_paths(iii)
      infile_iii  = addfile(infile_path, "r")
      area_lon = infile_iii->lon
      area_lat = infile_iii->lat
      var_lat = area_lat(:, 10)
      var_lon = area_lon(10, :)

      delete(area_lat)
      delete(area_lon)

      if(max(var_lon) .gt. 291) then  ; finds the closest lon grid point to 69W
        a = closest_val(291.0, var_lon)
      else
        a = closest_val(-69.0, var_lon)
      end if

      exact_lon = var_lon(a)
      var_final = dataset(:, :, a)  ; taking the lon value closest to 69W
      b1 = closest_val(-76, var_lat)
      b2 = closest_val(-48, var_lat)

      if (.not. all((ismissing(fx_var)))) then

        drlon = tofloat(abs(var_lon(a) - var_lon(a + 1))) * 6.37e06 * 0.0174533
        ; drlon = horizontal distance between 2 data poins along the same lon
        volcello_3d = fx_var
        volcello_3d = volcello_3d/drlon
        if(abs(volcello_3d&lev(0)) .gt. abs(volcello_3d&lev(1))) then
          volcello_2d = volcello_3d(::-1, :, a)
          ; if volcello has decreasing lev
        else
          volcello_2d = volcello_3d(:, :, a)
        end if
        if (any(dimsizes(volcello_2d) .ne. dimsizes(var_final))) then
          ; if dimensions of volcello are not same as uo
          volcello_2d_new = linint1_n_Wrap(volcello_2d&rlat, volcello_2d, \
                                           False, var_final&rlat, 0, 1)
          delete(volcello_2d)
          volcello_2d = volcello_2d_new
          delete(volcello_2d_new)
        end if

        transportpercell = var_final * volcello_2d / (10 ^ 6)
        ; to convert from m^3/s to sverdrups
        copy_VarCoords(var_final, transportpercell)
        transportperlat = dim_sum_n_Wrap(transportpercell, 0)
        ; sum of transport on each lat
        transportperlat = transportperlat / tofloat(cos(var_lat * 0.0174533))
        ; to compensate for changing horizontal distance with lat
        totaltransport = sum(transportperlat(b1:b2))
        delete(drlon)
        delete(volcello_2d)
        delete(transportpercell)
        delete(transportperlat)
        delete(volcello_3d)

      else
        totaltransport = "missing volume file"
        error_msg("w", "russell18jgr-fig4.ncl", " ", " volcello file for " \
                  + annots(iii) + " not found in the recipe. If volcello " \
                  + "file available, please copy the dataset name from " \
                  + "additional dataset section of uo to additional " \
                  + "dataset section of volcello. This is a warning, not " \
                  + "an error, it just skips the transport calcultions")
      end if
      delete(var_lat)
      delete(var_lon)

    else
      ; finds the closest lon grid point to 69 W
      if(max(dataset&lon) .gt. 291) then
        a = closest_val(291.0, dataset&lon)
      else
        a = closest_val(-69.0, dataset&lon)
      end if
      var_final = dataset(:, :, a)
      b1 = closest_val(-76, var_final&lat)
      b2 = closest_val(-48, var_final&lat)
      exact_lon = dataset&lon(a)  ; taking the lon value closest to 69W

      if (.not. all(ismissing(fx_var))) then
        dlon = tofloat(abs(dataset&lon(a) - dataset&lon(a+1))) * 6.37e06 \
              * 0.0174533
        ; dlon = horizontal distance between 2 data poins along the same lon
        volcello_3d = fx_var
        volcello_3d = volcello_3d/dlon

        if(abs(volcello_3d&lev(0)) .gt. abs(volcello_3d&lev(1))) then
          volcello_2d = volcello_3d(::-1, :, a)
          ; if volcello has decreasing lev
        else
          volcello_2d = volcello_3d(:, :, a)
        end if
        if (any(dimsizes(volcello_2d) .ne. dimsizes(var_final))) then
          ; if dimensions of volcello are not same as uo
          volcello_2d_new = linint1_n_Wrap(volcello_2d&lat, volcello_2d, \
                                           False, var_final&lat, 0, 1)
          delete(volcello_2d)
          volcello_2d = volcello_2d_new
          delete(volcello_2d_new)
        end if
        transportpercell = var_final * volcello_2d / (10 ^ 6)
        ; to convert from m^3/s to sverdrups
        transportperlat = dim_sum_n_Wrap(transportpercell, 0)
        transportperlat = transportperlat \
          / tofloat(cos(dataset&lat * 0.0174533))
        ; to compensate for changing horizontal distance with lat
        totaltransport = sum(transportperlat(b1:b2))

        delete(dlon)
        delete(volcello_2d)
        delete(transportpercell)
        delete(transportperlat)
        delete(volcello_3d)
      else
        totaltransport = "missing volume file"
        error_msg("w", "russell18jgr-fig4.ncl", " ", " volcello file for " \
                  + annots(iii) + " not found in the recipe. If volcello " \
                  + "file available, please copy the dataset name from " \
                  + "additional dataset section of uo to additional " \
                  + "dataset section of volcello. This is a warning, not " \
                  + "an error, it just skips the transport calcultions")
      end if
    end if
    delete(fx_var)

    if (exact_lon .gt. 200) then  ; editing the lon value for printing as W
      exact_lon = 360.0 - exact_lon
    else
      exact_lon = exact_lon * -1
    end if

    if(typeof(totaltransport) .eq. "float") then
      res@gsnLeftString = "Net transport : " \
                          + sprintf("%4.1f", totaltransport) + "Sv"
    else
      res@gsnLeftString = " no volcello file "
    end if

    if (isatt(diag_script_info, "unitCorrectionalFactor")) then
      ; to convert the units in datasets
      var_final = var_final * diag_script_info@unitCorrectionalFactor
    end if

    res@gsnRightString  = "Drake passage (" + sprintf("%4.2f", exact_lon) \
                          + " W)"
    res@tiMainString    = "Section velocity of " + annots(iii) + strUnits

    plots(iii) = gsn_csm_contour(wks, var_final(:, b1:b2), res)

    nc_filename = config_user_info@work_dir + "russell18jgr_fig4_" + var0 \
      + "_" + annots(iii) + "_" + (start_years_data(iii)) + "-" \
      + (end_years_data(iii)) + ".nc"

    ncdf_outfile = ncdf_write(dataset, nc_filename)
    delete(dataset)
    delete(var_final)
    delete(exact_lon)
    delete(totaltransport)
    ; Call provenance logger
    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 4", \
                   "mean", \
                   "sh", \
                   "geo", \
                   "russell_joellen", \
                   "russell18jgr", \
                   (/inputfile_paths(iii), volfile_paths(iii)/))

  end do

  ; Draw the panel
  pres = True
  pres@gsnPanelLabelBar = False
  outfile = panelling(wks, plots, nvert, nhori, pres)
  if (warn_ij) then
    error_msg("w", "russell18jgr-fig4.ncl", " ", " All models having " \
              + "i-j coordinates have been assigned lat-lon coordinates. " \
              + "This is not an error, no response is needed by user for " \
              + "plotting. ")
  end if
end
