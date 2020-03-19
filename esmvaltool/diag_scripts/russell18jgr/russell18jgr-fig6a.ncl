; #############################################################################
;
; russell18jgr_fig6a.ncl
;
;  Based on Figure 6a - Russell, J.L.,et al., 2018, J. Geophysical Research –
;   Oceans, 123, 3120-3143. https://doi.org/10.1002/2017JC013461 (figure 6a)
;
;  Author: Pandde Amarjiit (University of Arizona, USA)
;          Russell Joellen (University of Arizona, USA)
;          Goodman Paul    (University of Arizona, USA)
;
; #############################################################################
; Description
;
;    - Calculates the time average of thetao, so and vo
;    - Regrids temperature and salinity onto the vo grid
;    - Calculate the potential density of all cells relative to
;        0, 2000, 4000 decibars (or 1977m and 3948m)
;    - extracts the volcello of the lat closest to 30S
;    - divide volcello by latitudinal distance between two data points on
;       same lon, this gives us north-south cross-sectional area of each grid
;    - volume transport(m^3/s) in each cell is calculated by vo(m/s)
;        * cross sectional area(m^2)
;    - make lxx variables that have 1 for specific density layer and
;        missing values in rest
;    - multiply lxx with volume transport per cell to get volume transport
;        of cells in that layer (lvxx)
;    - total sum of lvxx gives the net volume transported in that density layer
;    - plots the volume transport per layer as a bar chart
;
;  Density layers defined as per : (Talley, L.D., 2003. Shallow, intermediate
;           and deep overturning components of the global heat budget.
;            Journal of Physical Oceanography 33, 530–560)
;
; Required Preprocessor attributes ( no_preprocessor )
;     - None (no preprocessing required)
;     in the recipe do not keep a preprocessor in variable section
;
; Required diag_script_info attributes
;
;     styleset = "CMIP5"  - default
;     ncdf     = "default"
;
; Optional diag_script_info attributes (diagnostic specific)
;
; Caveats
;
;  - MIROC-ESM and BNU-ESM do not work as depth variable is not called lev
;  - MRI models does not work as the data has 0 as fillvalue instead of 1e+20
;  - CCSM4 and CESM1-CAM5 dont work as the units for so is 1,
;      which is not accepted by ESMValTool
;  - Transport is very small in case of NorESM1-M and ME as volcello
;      values look incorrect(very small)
;
;
; Modification history
;  20190510 - russell_joellen, pandde_amarjiit - written and
;            implemented for ESMValTool v2.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"  ; load metadata
load "$diag_scripts/shared/plot/style.ncl"  ; load plot style functions
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  vo_items = select_metadata_by_name(input_file_info, "vo")
  so_items = select_metadata_by_name(input_file_info, "so")
  thetao_items = select_metadata_by_name(input_file_info, "thetao")
  volcello_items = select_metadata_by_name(input_file_info, "volcello")

  vo_datasets = metadata_att_as_array(vo_items, "dataset")
  start_years_data = metadata_att_as_array(vo_items, "start_year")
  end_years_data = metadata_att_as_array(vo_items, "end_year")
  vo_inputfile_paths = metadata_att_as_array(vo_items, "filename")
  thetao_inputfile_paths = metadata_att_as_array(thetao_items, "filename")
  so_inputfile_paths = metadata_att_as_array(so_items, "filename")
  nDatasets = ListCount(vo_items)
  nVolcello = ListCount(volcello_items)
  nVariables = ListCount(variable_info)
  dim_models = dimsizes(vo_datasets)
  if (nVolcello .ne. nDatasets) then
    error_msg("f", "russell18jgr-fig6a.ncl", " ", "volcello files " + \
              "for russell18jgr-fig6a.ncl do not match with vo datasets." + \
              " Please do not add additional variables in variable groups ")
  end if
end

begin

  plotpath = config_user_info@plot_dir + "/Russell_figure-6a_" \
            + sprinti("%0.4i", min(toint(start_years_data))) + "-" \
            + sprinti("%0.4i", max(toint(end_years_data)))

  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  wks = gsn_open_wks(output_type(), plotpath)
  plots = new(nDatasets, graphic)
  plot_1 = new(nDatasets, graphic)
  plot_talley = new(nDatasets, graphic)

  y_val  = (/-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5/)
  y1_val = fspan(0.125, 10.875, 44)
  yaxis_labels = (/-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11/)
  talley_SO_Zonal_all_Levels = (/0.0, -7.34, -1.9, 9.71, 2.37, -5.86, -10.02, \
                                -11.11, -2.84, 9.96, 16.49, 0.51 /)

  res                   = True           ; resources for plotting
  res@gsnXYBarChart     = True
  res@gsnDraw           = False
  res@gsnFrame          = False
  res@vpXF              = 0.1
  res@vpYF              = 0.75
  res@vpHeightF         = 0.5
  res@vpWidthF          = 0.4
  res@gsnMaximize       = True
  res@trYReverse        = True
  res@tmXBTickStartF    = -20.
  res@tmXBTickSpacingF  = 4.
  res@tmXBTickEndF      = 20.
  res@tmXBMode          = "Manual"
  res@tmXBMinorPerMajor = 1
  res@trXMinF           = -20.
  res@trXMaxF           = 20.
  res@trYMinF           = -1.
  res@trYMaxF           = 11.
  res@gsnXRefLine       = 0.0
  res@gsnXYBarChartColors2 = ("red")
  res@xyLineColors      = (/"black", "black"/)
  res@tmYLMode          = "Explicit"
  res@tmYLLabelsOn      = True
  res@tmYLValues        = yaxis_labels
  res@tmYLLabelFont     = 5
  res@tmYLLabelFontHeightF = 0.008
  res@tmYLLabels = (/"~F0~Net~F~", "~F0~Surface~F~", "26.10s~B~0",  \
                    "26.40s~B~0", "26.90s~B~0", "27.10s~B~0", "27.40s~B~0", \
                    "36.8s~B~2", "45.8s~B~4", "45.86s~B~4", "45.92s~B~4", \
                    "46.0s~B~4", "~F0~Bottom~F~" /)
  res@gsnRightStringFontHeightF = 0.009
  res@gsnLeftStringFontHeightF = 0.009

  do iii = 0, dim_models - 1

    fx_var = read_data(volcello_items[iii])

    if (all(ismissing(fx_var))) then
    ; to give fatal error if volcello is missing
      fx_variable = "volcello"
      error_msg("f", "russell18jgr-fig6.ncl", " ", "volcello file for "  \
                + vo_datasets(iii) \
                + " not found in the metadata file, please add "\
                + "'fx_files: [volcello]' to the variable dictionary in the " \
                + "recipe or add the location of file to input directory " \
                + "in config-user.yml ")
    end if

    dataset_so_time = read_data(so_items[iii])
    dataset_so = dim_avg_n_Wrap(dataset_so_time, 0)
    delete(dataset_so_time)
    dataset_thetao_time = read_data(thetao_items[iii])
    dataset_thetao = dim_avg_n_Wrap(dataset_thetao_time, 0)
    delete(dataset_thetao_time)

    if (max(dataset_thetao) .gt. 250) then  ; making sure temperature is in K
      dataset_thetao = dataset_thetao - 273.15
    end if

    rho0 = new(2, double)
    assignFillValue(dataset_thetao, rho0)

    vo_file = addfile(vo_inputfile_paths(iii), "r")
    thetao_file = addfile(thetao_inputfile_paths(iii), "r")
    var_test_lat = vo_file->lat                 ; extracting lat array of vo
    var_test_lon = vo_file->lon                 ; extracting lon array of vo

    if((iscoord(dataset_so, "lat"))) then
      a = closest_val(-30.0, var_test_lat)
      ; getting index for lat closest to 30S
      exact_lat   =  var_test_lat(a)            ; lat value closest to -30
      ; interpolate thetao and so onto vo grid
      theta_inter = linint1_n_Wrap(dataset_thetao&lat, dataset_thetao, \
                                   False, exact_lat, 0, 1)
      theta_new   = linint1_n_Wrap(dataset_thetao&lon, theta_inter, \
                                   True, var_test_lon, 0, 2)
      delete(theta_inter)
      delete(dataset_thetao)
      so_inter = linint1_n_Wrap(dataset_so&lat, dataset_so, False, \
                                exact_lat, 0, 1)
      so_new   = linint1_n_Wrap(dataset_so&lon, so_inter, True, \
                                var_test_lon, 0, 2)
      delete(dataset_so)
      delete(so_inter)
      ; potential density calculation
      rho0_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 0.0)
      rho2_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 1977.0)
      rho4_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 3948.0)

      delete(theta_new)
      delete(so_new)

    elseif (iscoord(dataset_so, "rlat")) then

      theta_lat    = thetao_file->lat       ; extracting lat array of thetao
      theta_lon    = thetao_file->lon       ; extracting lon arrays of thetao
      lat1 = closest_val(-15.0, var_test_lat(:, 0))
      ; index value of vo-lat closest to -10
      lat2 = closest_val(-45.0, var_test_lat(:, 0))
      ; index value of vo-lat closest to -50
      a = closest_val(-30.0, var_test_lat(:, 0))
      ; index value of vo-lat closest to -30
      exact_lat   =  var_test_lat(a, 0)          ; lat value closest to -30
      ; interpolate thetao and so onto vo grid
      theta_inter = linint1_n_Wrap(theta_lat(lat2:lat1, 1), \
                                   dataset_thetao(:, lat2:lat1, :), False, \
                                   exact_lat, 0, 1)
      delete(dataset_thetao)

      so_inter   = linint1_n_Wrap(theta_lat(lat2:lat1, 1), \
                                  dataset_so(:, lat2:lat1, :), False, \
                                  exact_lat, 0, 1)
      delete(dataset_so)

      ; checking for monotonic nature of lon arrays, as NCL
      ;   needs monotonically increasing lon for interpolations.
      if ((isMonotonic(theta_lon(a, :)) .eq. 1) .and. \
          (isMonotonic(var_test_lon(a, :)) .eq. 1)) then
        theta_new = linint1_n_Wrap(theta_lon(a, :), theta_inter, \
                                   True, var_test_lon(a, :), 0, 2)
      else
        theta_new = theta_inter
      end if
      delete(theta_inter)

      if ((isMonotonic(theta_lon(a, :)) .eq. 1) .and. \
          (isMonotonic(var_test_lon(a, :)) .eq. 1)) then
        so_new = linint1_n_Wrap(theta_lon(a, :), so_inter, True, \
                                var_test_lon(a, :), 0, 2)
      else
        so_new = so_inter
      end if
      delete(theta_lat)
      delete(theta_lon)
      delete(so_inter)

      ; potential density calculation
      rho0_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 0.)
      rho2_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 1977.0)
      rho4_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 3948.0)
      delete(theta_new)
      delete(so_new)
    else
      theta_lat    = thetao_file->lat      ; extracting lat array of thetao
      theta_lon    = thetao_file->lon      ; extracting lon arrays of thetao

      lat1 = closest_val(-15.0, var_test_lat(:, 0))
      ; index value of vo-lat closest to -10
      lat2 = closest_val(-45.0, var_test_lat(:, 0))
      ; index value of vo-lat closest to -50
      a = closest_val(-30.0, var_test_lat(:, 0))
      ; index value of vo-lat closest to -30
      exact_lat   =  var_test_lat(a, 0)   ; lat value closest to -30
      ; interpolate thetao and so onto vo grid

      theta_inter = linint1_n_Wrap(theta_lat(lat2:lat1, 1), \
                                   dataset_thetao(:, lat2:lat1, :), False, \
                                   exact_lat, 0, 1)
      delete(dataset_thetao)

      so_inter = linint1_n_Wrap(theta_lat(lat2:lat1, 1), \
                                dataset_so(:, lat2:lat1, :), False, \
                                exact_lat, 0, 1)
      delete(dataset_so)

      ; checking for monotonic nature of lon arrays, as NCL needs
      ;      monotonically increasing lon for interpolations.
      if ((isMonotonic(theta_lon(a, :)) .eq. 1) .and. \
          (isMonotonic(var_test_lon(a, :)) .eq. 1)) then
        theta_new = linint1_n_Wrap(theta_lon(a, :), theta_inter, True, \
                                   var_test_lon(a, :), 0, 2)
      else
        theta_new = theta_inter
      end if

      if ((isMonotonic(theta_lon(a, :)) .eq. 1) .and. \
          (isMonotonic(var_test_lon(a, :)) .eq. 1)) then
        so_new = linint1_n_Wrap(theta_lon(a, :), so_inter, True, \
                                var_test_lon(a, :), 0, 2)
      else
        so_new = so_inter
      end if

      delete(theta_inter)
      delete(so_inter)
      delete(theta_lat)
      delete(theta_lon)

      ; potential density calculation
      rho0_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 0.)
      rho2_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 1977.0)
      rho4_30s = rho_mwjf(theta_new(:, 0, :), so_new(:, 0, :), 3948.0)
      delete(theta_new)
      delete(so_new)

    end if

    dataset_vo_time = read_data(vo_items[iii])
    dataset_vo = dim_avg_n_Wrap(dataset_vo_time, 0)
    delete(dataset_vo_time)
    delete(var_test_lon)
    volcello = fx_var
    delete(fx_var)

    ; transport calculation
    if(iscoord(dataset_vo, "lat")) then

      var_tmp2 = dataset_vo(:, a, :)  ; y velocity of lat closest to -30

      dlat = tofloat(abs(dataset_vo&lat(a) - dataset_vo&lat(a + 1)) \
                     * 6.37e06 * 0.0174533)
      ; dlat is the north-south distance between 2 consecutive data latitudes
      volumecello_3d = volcello/dlat
      ; north-south crossectional area is volcello / dlat
      if(abs(volcello&lev(0)) .gt. abs(volcello&lev(1))) then
        ; extracting volume of cells closest to 30S
        volcello_2d = volumecello_3d(::-1, a, :)
        ; reversing lev if it is in descending order.
      else
        volcello_2d = volumecello_3d(:, a, :)
      end if
      transportpercell = var_tmp2 * volcello_2d / (10 ^ 6)
      ; unit conversion from m^3/s to sverdrups
      copy_VarCoords(var_tmp2, transportpercell)

    elseif (iscoord(dataset_vo, "rlat")) then

      var_tmp2 = dataset_vo(:, a, :)   ; y transport of lat closest to -30
      drlat = tofloat(abs(var_test_lat(a, 0) - var_test_lat(a+1, 0)) \
                      * 6.37e06 * 0.0174533)
      ; drlat is north-south distance between 2 consecutive data latitudes
      volumecello_3d = volcello/drlat
      ; north-south crossectional area is volcello / drlat
      if(abs(volcello&lev(0)) .gt. abs(volcello&lev(1))) then
        ; extracting volume of cells closest to 30S
        volcello_2d = volumecello_3d(::-1, a, :)
        ; reversing lev if it is in descending order.
      else
        volcello_2d = volumecello_3d(:, a, :)
      end if
      transportpercell = var_tmp2 * volcello_2d / (10 ^ 6)
      ; unit conversion from m^3/s to sverdrups
      copy_VarCoords(var_tmp2, transportpercell)

    else

      var_tmp2 = dataset_vo(:, a, :)   ; y transport of lat closest to 30S
      dlat = tofloat(abs(var_test_lat(a, 0) - var_test_lat(a + 1, 0)) \
                     * 6.37e06 * 0.0174533)
      ; dlat is the north-south distance between 2 consecutive data latitudes
      volumecello_3d = volcello / dlat
      ; north-south crossectional area is volcello / dlat

      if(abs(volcello&lev(1)) .gt. abs(volcello&lev(2))) then
        ; extracting volume of cells closest to 30S
        volcello_2d = volumecello_3d(::-1, a, :)
        ; reversing lev if it is in descending order.
      else
        volcello_2d = volumecello_3d(:, a, :)
      end if
      transportpercell = var_tmp2 * volcello_2d / (10 ^ 6)
      ; unit conversion from m^3/s to sverdrups
      copy_VarCoords(var_tmp2, transportpercell)

    end if
    delete(volcello)
    delete(var_test_lat)
    delete(dataset_vo)
    delete(volcello_2d)
    delete(volumecello_3d)
    delete(var_tmp2)

    rho0_30s = 1000. * (rho0_30s - 1.0)
    rho2_30s = 1000. * (rho2_30s - 1.0)
    rho4_30s = 1000. * (rho4_30s - 1.0)

    ; making masks for density levels based on : Talley, L.D., 2003.
    l10 = where((rho0_30s .lt. 26.10), 1, rho0@_FillValue)
    l11 = where((rho0_30s .lt. 24.900), 1, rho0@_FillValue)
    l12 = where((rho0_30s .ge. 24.900) .and. (rho0_30s .lt. 25.300), 1, \
                rho0@_FillValue)
    l13 = where((rho0_30s .ge. 25.300) .and. (rho0_30s .lt. 25.700), 1, \
                rho0@_FillValue)
    l14 = where((rho0_30s .ge. 25.700) .and. (rho0_30s .lt. 26.100), 1, \
                rho0@_FillValue)

    l20 = where((rho0_30s .ge. 26.100) .and. (rho0_30s .lt. 26.400), 1, \
                rho0@_FillValue)
    l21 = where((rho0_30s .ge. 26.100) .and. (rho0_30s .lt. 26.175), 1, \
                rho0@_FillValue)
    l22 = where((rho0_30s .ge. 26.175) .and. (rho0_30s .lt. 26.250), 1, \
                rho0@_FillValue)
    l23 = where((rho0_30s .ge. 26.250) .and. (rho0_30s .lt. 26.325), 1, \
                rho0@_FillValue)
    l24 = where((rho0_30s .ge. 26.325) .and. (rho0_30s .lt. 26.400), 1, \
                rho0@_FillValue)

    l30 = where((rho0_30s .ge. 26.400) .and. (rho0_30s .lt. 26.900), 1, \
                rho0@_FillValue)
    l31 = where((rho0_30s .ge. 26.400) .and. (rho0_30s .lt. 26.525), 1, \
                rho0@_FillValue)
    l32 = where((rho0_30s .ge. 26.525) .and. (rho0_30s .lt. 26.650), 1, \
                rho0@_FillValue)
    l33 = where((rho0_30s .ge. 26.650) .and. (rho0_30s .lt. 26.775), 1, \
                rho0@_FillValue)
    l34 = where((rho0_30s .ge. 26.775) .and. (rho0_30s .lt. 26.900), 1, \
                rho0@_FillValue)

    l40 = where((rho0_30s .ge. 26.900) .and. (rho0_30s .lt. 27.100), 1, \
                rho0@_FillValue)
    l41 = where((rho0_30s .ge. 26.900) .and. (rho0_30s .lt. 26.950), 1, \
                rho0@_FillValue)
    l42 = where((rho0_30s .ge. 26.950) .and. (rho0_30s .lt. 27.000), 1, \
                rho0@_FillValue)
    l43 = where((rho0_30s .ge. 27.000) .and. (rho0_30s .lt. 27.050), 1, \
                rho0@_FillValue)
    l44 = where((rho0_30s .ge. 27.050) .and. (rho0_30s .lt. 27.100), 1, \
                rho0@_FillValue)

    l50 = where((rho0_30s .ge. 27.100) .and. (rho0_30s .lt. 27.400), 1, \
                rho0@_FillValue)
    l51 = where((rho0_30s .ge. 27.100) .and. (rho0_30s .lt. 27.175), 1, \
                rho0@_FillValue)
    l52 = where((rho0_30s .ge. 27.175) .and. (rho0_30s .lt. 27.250), 1, \
                rho0@_FillValue)
    l53 = where((rho0_30s .ge. 27.250) .and. (rho0_30s .lt. 27.325), 1, \
                rho0@_FillValue)
    l54 = where((rho0_30s .ge. 27.325) .and. (rho0_30s .lt. 27.400), 1, \
                rho0@_FillValue)

    l60 = where((rho0_30s .ge. 27.400) .and. (rho2_30s .lt. 36.800), 1, \
                rho0@_FillValue)
    l61 = where((rho0_30s .ge. 27.400) .and. (rho0_30s .lt. 27.500), 1, \
                rho0@_FillValue)
    l62 = where((rho0_30s .ge. 27.500) .and. (rho2_30s .lt. 36.700), 1, \
                rho0@_FillValue)
    l63 = where((rho2_30s .ge. 36.700) .and. (rho2_30s .lt. 36.750), 1, \
                rho0@_FillValue)
    l64 = where((rho2_30s .ge. 36.750) .and. (rho2_30s .lt. 36.800), 1, \
                rho0@_FillValue)

    l70 = where((rho2_30s .ge. 36.800) .and. (rho4_30s .lt. 45.800), 1, \
                rho0@_FillValue)
    l71 = where((rho2_30s .ge. 36.800) .and. (rho2_30s .lt. 36.850), 1, \
                rho0@_FillValue)
    l72 = where((rho2_30s .ge. 36.850) .and. (rho2_30s .lt. 36.900), 1, \
                rho0@_FillValue)
    l73 = where((rho2_30s .ge. 36.900) .and. (rho2_30s .lt. 36.950), 1, \
                rho0@_FillValue)
    l74 = where((rho2_30s .ge. 36.950) .and. (rho4_30s .lt. 45.800), 1, \
                rho0@_FillValue)

    l80 = where((rho4_30s .ge. 45.800) .and. (rho4_30s .lt. 45.860), 1, \
                rho0@_FillValue)
    l81 = where((rho4_30s .ge. 45.800) .and. (rho4_30s .lt. 45.815), 1, \
                rho0@_FillValue)
    l82 = where((rho4_30s .ge. 45.815) .and. (rho4_30s .lt. 45.830), 1, \
                rho0@_FillValue)
    l83 = where((rho4_30s .ge. 45.830) .and. (rho4_30s .lt. 45.845), 1, \
                rho0@_FillValue)
    l84 = where((rho4_30s .ge. 45.845) .and. (rho4_30s .lt. 45.860), 1, \
                rho0@_FillValue)

    l90 = where((rho4_30s .ge. 45.860) .and. (rho4_30s .lt. 45.920), 1, \
                rho0@_FillValue)
    l91 = where((rho4_30s .ge. 45.860) .and. (rho4_30s .lt. 45.875), 1, \
                rho0@_FillValue)
    l92 = where((rho4_30s .ge. 45.875) .and. (rho4_30s .lt. 45.890), 1, \
                rho0@_FillValue)
    l93 = where((rho4_30s .ge. 45.890) .and. (rho4_30s .lt. 45.905), 1, \
                rho0@_FillValue)
    l94 = where((rho4_30s .ge. 45.905) .and. (rho4_30s .lt. 45.920), 1, \
                rho0@_FillValue)

    l100 = where((rho4_30s .ge. 45.920) .and. (rho4_30s .lt. 46.000), 1, \
                 rho0@_FillValue)
    l101 = where((rho4_30s .ge. 45.920) .and. (rho4_30s .lt. 45.940), 1, \
                 rho0@_FillValue)
    l102 = where((rho4_30s .ge. 45.940) .and. (rho4_30s .lt. 45.960), 1, \
                 rho0@_FillValue)
    l103 = where((rho4_30s .ge. 45.960) .and. (rho4_30s .lt. 45.980), 1, \
                 rho0@_FillValue)
    l104 = where((rho4_30s .ge. 45.980) .and. (rho4_30s .lt. 46.000), 1, \
                 rho0@_FillValue)

    l110 = where((rho4_30s .ge. 46.000), 1, rho0@_FillValue)
    l111 = where((rho4_30s .ge. 46.000) .and. (rho4_30s .lt. 46.050), 1, \
                 rho0@_FillValue)
    l112 = where((rho4_30s .ge. 46.050) .and. (rho4_30s .lt. 46.100), 1, \
                 rho0@_FillValue)
    l113 = where((rho4_30s .ge. 46.100) .and. (rho4_30s .lt. 46.150), 1, \
                 rho0@_FillValue)
    l114 = where((rho4_30s .ge. 46.150), 1, rho0@_FillValue)

    delete(rho0_30s)
    delete(rho2_30s)
    delete(rho4_30s)

    ; assignning filling values from rho0 to masked layers
    assignFillValue(rho0, l10)
    assignFillValue(rho0, l20)
    assignFillValue(rho0, l30)
    assignFillValue(rho0, l40)
    assignFillValue(rho0, l50)
    assignFillValue(rho0, l60)
    assignFillValue(rho0, l70)
    assignFillValue(rho0, l80)
    assignFillValue(rho0, l90)
    assignFillValue(rho0, l100)
    assignFillValue(rho0, l110)
    assignFillValue(rho0, l11)
    assignFillValue(rho0, l21)
    assignFillValue(rho0, l31)
    assignFillValue(rho0, l41)
    assignFillValue(rho0, l51)
    assignFillValue(rho0, l61)
    assignFillValue(rho0, l71)
    assignFillValue(rho0, l81)
    assignFillValue(rho0, l91)
    assignFillValue(rho0, l101)
    assignFillValue(rho0, l111)
    assignFillValue(rho0, l12)
    assignFillValue(rho0, l22)
    assignFillValue(rho0, l32)
    assignFillValue(rho0, l42)
    assignFillValue(rho0, l52)
    assignFillValue(rho0, l62)
    assignFillValue(rho0, l72)
    assignFillValue(rho0, l82)
    assignFillValue(rho0, l92)
    assignFillValue(rho0, l102)
    assignFillValue(rho0, l112)
    assignFillValue(rho0, l13)
    assignFillValue(rho0, l23)
    assignFillValue(rho0, l33)
    assignFillValue(rho0, l43)
    assignFillValue(rho0, l53)
    assignFillValue(rho0, l63)
    assignFillValue(rho0, l73)
    assignFillValue(rho0, l83)
    assignFillValue(rho0, l93)
    assignFillValue(rho0, l103)
    assignFillValue(rho0, l113)
    assignFillValue(rho0, l14)
    assignFillValue(rho0, l24)
    assignFillValue(rho0, l34)
    assignFillValue(rho0, l44)
    assignFillValue(rho0, l54)
    assignFillValue(rho0, l64)
    assignFillValue(rho0, l74)
    assignFillValue(rho0, l84)
    assignFillValue(rho0, l94)
    assignFillValue(rho0, l104)
    assignFillValue(rho0, l114)
    delete(rho0)

    lv = new((/12/), double)     ; lv is array of big blue bars
    assignFillValue(l10, lv)
    lv(1)  = sum(l10 * transportpercell)
    lv(2)  = sum(l20 * transportpercell)
    lv(3)  = sum(l30 * transportpercell)
    lv(4)  = sum(l40 * transportpercell)
    lv(5)  = sum(l50 * transportpercell)
    lv(6)  = sum(l60 * transportpercell)
    lv(7)  = sum(l70 * transportpercell)
    lv(8)  = sum(l80 * transportpercell)
    lv(9)  = sum(l90 * transportpercell)
    lv(10) = sum(l100 * transportpercell)
    lv(11) = sum(l110 * transportpercell)
    lv(0)  = 0
    lv(0)  = sum(lv)

    lv0 = new((/44/), double)     ; lv0 is array of small red bars
    assignFillValue(l104, lv0)
    lv0(0)   = sum(l11 * transportpercell)
    lv0(1)   = sum(l12 * transportpercell)
    lv0(2)   = sum(l13 * transportpercell)
    lv0(3)   = sum(l14 * transportpercell)
    lv0(4)   = sum(l21 * transportpercell)
    lv0(5)   = sum(l22 * transportpercell)
    lv0(6)   = sum(l23 * transportpercell)
    lv0(7)   = sum(l24 * transportpercell)
    lv0(8)   = sum(l31 * transportpercell)
    lv0(9)   = sum(l32 * transportpercell)
    lv0(10)  = sum(l33 * transportpercell)
    lv0(11)  = sum(l34 * transportpercell)
    lv0(12)  = sum(l41 * transportpercell)
    lv0(13)  = sum(l42 * transportpercell)
    lv0(14)  = sum(l43 * transportpercell)
    lv0(15)  = sum(l44 * transportpercell)
    lv0(16)  = sum(l51 * transportpercell)
    lv0(17)  = sum(l52 * transportpercell)
    lv0(18)  = sum(l53 * transportpercell)
    lv0(19)  = sum(l54 * transportpercell)
    lv0(20)  = sum(l61 * transportpercell)
    lv0(21)  = sum(l62 * transportpercell)
    lv0(22)  = sum(l63 * transportpercell)
    lv0(23)  = sum(l64 * transportpercell)
    lv0(24)  = sum(l71 * transportpercell)
    lv0(25)  = sum(l72 * transportpercell)
    lv0(26)  = sum(l73 * transportpercell)
    lv0(27)  = sum(l74 * transportpercell)
    lv0(28)  = sum(l81 * transportpercell)
    lv0(29)  = sum(l82 * transportpercell)
    lv0(30)  = sum(l83 * transportpercell)
    lv0(31)  = sum(l84 * transportpercell)
    lv0(32)  = sum(l91 * transportpercell)
    lv0(33)  = sum(l92 * transportpercell)
    lv0(34)  = sum(l93 * transportpercell)
    lv0(35)  = sum(l94 * transportpercell)
    lv0(36)  = sum(l101 * transportpercell)
    lv0(37)  = sum(l102 * transportpercell)
    lv0(38)  = sum(l103 * transportpercell)
    lv0(39)  = sum(l104 * transportpercell)
    lv0(40)  = sum(l111 * transportpercell)
    lv0(41)  = sum(l112 * transportpercell)
    lv0(42)  = sum(l113 * transportpercell)
    lv0(43)  = sum(l114 * transportpercell)
    if(any(ismissing(lv))) then  ; if all cells are missing then print 0
      aa = ind(ismissing(lv))
      lv(aa) = 0
      delete(aa)
    end if
    delete([/ l10, l11, l12, l13, l14, l20, l21, l22, l23, l24, l30, l31 /])
    delete([/ l32, l33, l34, l40, l41, l42, l43, l44, l50, l51, l52, l53 /])
    delete([/ l54, l60, l61, l62, l63, l64, l70, l71, l72, l73, l74, l80 /])
    delete([/ l81, l82, l83, l84, l90, l91, l92, l93, l94, l100, l101 /])
    delete([/ l102, l103, l104, l110, l111, l112, l113, l114 /])
    delete(transportpercell)
    exact_lat = exact_lat * -1
    strUnits = ""

    res@gsnLeftString = "(" + start_years_data(iii) + " - " + \
                        end_years_data(iii) + ") at (" \
                        + sprintf("%4.2f", exact_lat) + "S)"
    res@gsnRightString  = "Net Transport out of Southern ocean = " + \
                          sprintf("%4.2f", lv(0)) + "Sv"
    res@tiMainString    = vo_datasets(iii)
    res@xyLineColors    = (/"black", "black"/)
    res@gsnXYBarChartColors2 = ("red")

    plot_1(iii) = gsn_csm_xy(wks, lv0, y1_val, res)  ; plotting of red bars

    res@gsnXYBarChartColors2 = ("blue4")
    res@xyLineColors         = (/"blue4", "blue4"/)
    plots(iii) = gsn_csm_xy(wks, lv, y_val, res)    ; plotting of blue bars
    txres                 = True   ; to print volume transport of blue bars
    txres@gsnFrame        = False
    txres@txFontHeightF   = 0.009
    gsn_text_ndc(wks, sprintf("%4.3f", lv(0)), 0.27, 0.845, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(1)), 0.27, 0.78, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(2)), 0.27, 0.715, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(3)), 0.27, 0.655, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(4)), 0.27, 0.58, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(5)), 0.27, 0.51, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(6)), 0.27, 0.435, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(7)), 0.27, 0.375, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(8)), 0.27, 0.31, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(9)), 0.27, 0.2475, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(10)), 0.27, 0.175, txres)
    gsn_text_ndc(wks, sprintf("%4.3f", lv(11)), 0.27, 0.11, txres)

    res@gsnXYBarChartColors2 = ("White")
    res@xyLineColors      = (/"magenta2", "magenta2"/)
    plot_talley(iii) = gsn_csm_xy(wks, talley_SO_Zonal_all_Levels, y_val, res)
    ; to plot the magenta for talley values of transport
    overlay(plots(iii), plot_1(iii))
    overlay(plot_talley(iii), plots(iii))
    draw(plot_talley(iii))
    frame(wks)
    out_var = new((/56/), double)
    out_var(0:11) = (/lv/)
    out_var(12:55) = (/lv0/)
    out_var!0 = "i"
    out_var&i = ispan(0, 55, 1)
    out_var@var = "transport_per_layer"
    out_var@diag_script = "russell18jgr_fig6a.ncl"
    out_var@description = "Transport in main layers(blue bars) in i(0-11)" + \
      " and transport in sub layers (red bars) in i(12-55) at " + \
      exact_lat + "S of model " + vo_datasets(iii)
    delete(lv)
    delete(lv0)
    delete(exact_lat)

    nc_filename = config_user_info@work_dir + "russell18jgr-figure6a_" \
      + vo_datasets(iii) + "_" + (start_years_data(iii)) + "-" + \
      (end_years_data(iii)) + ".nc"
    ncdf_outfile = ncdf_write(out_var, nc_filename)

    log_provenance(ncdf_outfile, \
                   plotpath, \
                   "Russell et al 2018 figure 6 part a", \
                   "mean", \
                   "sh",   \
                   (/"bar", "vert"/),  \
                   "russell_joellen", \
                   "russell18jgr", \
                   (/ vo_inputfile_paths(iii), thetao_inputfile_paths(iii), \
                    so_inputfile_paths(iii)/))
  end do

end
