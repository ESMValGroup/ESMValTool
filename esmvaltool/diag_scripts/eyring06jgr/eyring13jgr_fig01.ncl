; #############################################################################
; SCRIPT TO REPRODUCE FIG. 1 OF EYRING ET AL. (JGR, 2013)
; Authors: Irene Cionni (ENEA, Italy) and Veronika Eyring (DLR, Germany)
; CCMI Project
; #############################################################################
; #############################################################################
; # Author: Irene Cionni (ENEA, Italy)
; # ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
; #
; # Description:
; #     Calculates seasonal cycles of zonally averaged total ozone columns.
; #
; # Required diag_script_info attributes (diagnostics specific):
; #     rgb_file: path + filename of color table
; #         (e.g., "diag_scripts/lib/ncl/rgb/eyring_toz.rgb")
; #     styleset: style set (“DEFAULT, "CMIP5")
; #     multimean = "True"
; #     list_chem_mod: list of models in the group “chem” (array of
; #         strings, default = (/"All"/))
; #     list_chem_mod_string: plotting label for group “chem”,
; #         e.g., "CMIP5"
; #     list_nochem_mod: list of models in the group “nochem”
; #         (array of strings, default = (/""/))
; #     list_nochem_mod_string: plotting label for gro “nochem”,
; #         e.g., "NoChem"
; #     diff_ref: name of reference model for difference plots,
; #         e.g., "NIWA"
; ############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  ; Get variables and datasets
  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD1 = ListCount(info0)
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")
  y1 = (/ metadata_att_as_array(info0, "start_year")/)
  y2 = (/ metadata_att_as_array(info0, "end_year") /)

  ; Check for reference dataset definition
  if (variable_info[0]@reference_dataset.eq."None") then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset is specified")
  end if
  iref = ind(names.eq.variable_info[0]@reference_dataset)

  ; Check variable attributes
  if(isatt(variable_info, "units")) then
    units = variable_info@units
  else
    units = "DU"
  end if

  if(isatt(variable_info, "long_name")) then
    long_name = variable_info@long_name
  else
    long_name = var0
  end if

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "Min", 185.)
  set_default_att(diag_script_info, "Max", 500.)
  set_default_att(diag_script_info, "Spacing", 25.)
  set_default_att(diag_script_info, "Mindiff", -50.)
  set_default_att(diag_script_info, "Maxdiff", 50.)
  set_default_att(diag_script_info, "Spacingdiff", 5.)
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "start_year", max(y1))
  set_default_att(diag_script_info, "end_year", min(y2))
  set_default_att(diag_script_info, "list_chem_mod", "All")

  set_default_att(diag_script_info, "list_chem_mod_string", "CMIP5")
  set_default_att(diag_script_info, "list_nochem_mod_string", "")
  set_default_att(diag_script_info, "diff_ref", "NIWA_BS")

  flag_mod = where(projects.ne."OBS6" .and. projects.ne."OBS" \
                   .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  Max = diag_script_info@Max
  Min = diag_script_info@Min
  Spacing = diag_script_info@Spacing

  Maxdiff = diag_script_info@Maxdiff
  Mindiff = diag_script_info@Mindiff
  Spacingdiff = diag_script_info@Spacingdiff
  multimean = diag_script_info@multimean

  if (multimean)
    dimTot = dim_MOD + dim_OBS + 1
  else
    dimTot = dim_MOD + dim_OBS
  end if
  list_chem_mod = diag_script_info@list_chem_mod
  attrib_chem = "_CHEM"
  if (list_chem_mod.eq."All") then
    delete(list_chem_mod)
    list_chem_mod = names
    attrib_chem = " "
  end if

  ind_chem_mod = ind((names.eq.list_chem_mod).and.(flag_mod.gt.0))
  if ((isatt(diag_script_info, "list_nochem_mod_string"))) then
    list_nochem_mod = diag_script_info@list_nochem_mod
    ind_nochem_mod = ind((names.eq.list_nochem_mod).and.(flag_mod.gt.0))
  end if
  OBS_ref = names(iref)

  ; make sure path for (mandatory) netcdf output exists
  year_start = toint(diag_script_info@start_year)
  year_end = toint(diag_script_info@end_year)

  ; Create output directories
  system("mkdir -p " + config_user_info@work_dir)
  system("mkdir -p " + config_user_info@plot_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

end

begin
  plotpath2 = config_user_info@plot_dir + "eyring13f01_panelmod" + var0
  wks2 = gsn_open_wks(file_type, plotpath2)
  plotpath1 = config_user_info@plot_dir + "eyring13f01_multi" + var0
  wks1 = gsn_open_wks(file_type, plotpath1)

  ; output_filename = wks@fullname
  A0 = read_data(info0[iref])
  lat0 =  A0&lat
  delete(A0)
  tmpData = new((/dim_MOD + dim_OBS, 12, dimsizes(lat0), 2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "month"
  tmpData!2 = "lat"
  tmpData!3 = "stat"

  masktmpData = new((/dim_MOD + dim_OBS, 12, dimsizes(lat0), 2/), "double")
  masktmpData!0 = "line"
  masktmpData!1 = "month"
  masktmpData!2 = "lat"
  masktmpData!3 = "stat"
  masktmpData = (/ 1 /)

  XXdata = new((/dimTot, 12, dimsizes(lat0), 2/), "double")
  XXdata!0 = "line"
  XXdata!1 = "month"
  XXdata!2 = "lat"
  XXdata!3 = "stat"

  maskXXdata = new((/dimTot, 12, dimsizes(lat0), 2/), "double")
  maskXXdata!0 = "line"
  maskXXdata!1 = "month"
  maskXXdata!2 = "lat"
  maskXXdata!3 = "stat"
  maskXXdata = (/ 1 /)

  XXref = new((/ 12, dimsizes(lat0), 2/), "double")
  XXref!0 = "month"
  XXref!1 = "lat"

  years1 = new(dimTot, "integer")
  years2 = new(dimTot, "integer")

  plot = new(dimTot, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  plotMod = new(dimTot, typeof(names))

  ik = 0

  do imod = 0, dim_MOD1 - 1
    years1(imod) = y1(imod)
    years2(imod) = y2(imod)
    plotMod(imod) = names(imod)
    A0 = read_data(info0[imod])
    dummy1 = time_operations(A0, year_start, year_end, "average", \
                             "monthlyclim", True)
    dummy2_b = time_operations(A0, year_start, year_end, "extract", \
                               "", True)
    date = yyyymm_time(year_start, year_end, "integer")
    yyyy = date / 100
    mm = date - yyyy * 100
    weights = days_in_month(yyyy, mm)
    dummy2 = new((/12, dimsizes(lat0)/), "double")
    do ilat = 0, dimsizes(lat0) - 1
      do im = 0, 11
        if (.not.all(ismissing(dummy2_b(im :: 12, ilat)))) then
          dummy2(im, ilat) = \
            dim_stddev_wgt_Wrap(dummy2_b(im :: 12, ilat), weights(im :: 12), 1)
        end if
      end do
    end do

    dummy2!0 = "month"
    dummy2!1 = "lat"

    delete(A0)
    delete(weights)
    delete(date)
    delete(yyyy)
    delete(mm)
    delete(dummy2_b)
    dummy2!0 = dummy1!0

    tmpData(imod, :, :, 0) = (/dummy1/)
    tmpData(imod, :, :, 1) = (/dummy2/)

    masktmpData(imod, :, :, 0) = where(ismissing(dummy1), \
                                       masktmpData@_FillValue, \
                                       masktmpData(imod, :, :, 0))
    masktmpData(imod, :, :, 1) = where(ismissing(dummy2), \
                                       masktmpData@_FillValue, \
                                       masktmpData(imod, :, :, 1))
    delete(dummy1)
    delete(dummy2)
  end do
  XXdata(:dim_MOD1 - 1, :, :, :) = (/tmpData/)
  maskXXdata(:dim_MOD1 - 1, :, :, :) = (/masktmpData/)
  XXref = (/tmpData(iref, :, :, :)/)

  XXdata!0 = "line"
  XXdata!1 = "month"
  XXdata!2 = "lat"
  XXdata!3 = "stat"
  XXdata&lat = lat0
  XXdata&month = ispan(1, 12, 1)

  maskXXdata!0 = "line"
  maskXXdata!1 = "month"
  maskXXdata!2 = "lat"
  maskXXdata!3 = "stat"
  maskXXdata&lat = lat0
  maskXXdata&month = ispan(1, 12, 1)

  var_grd_mod = XXdata(ind_chem_mod, :, :, 0)

  MMM_chem = dim_avg_n_Wrap(var_grd_mod, 0)
  MMM_chem&month = ispan(1, 12, 1)
  maskMMM_chem = MMM_chem
  maskMMM_chem = (/ 1. /)
  maskMMM_chem = where(ismissing(MMM_chem), maskMMM_chem@_FillValue, \
                       maskMMM_chem)
  maskMMM_chem = maskMMM_chem * MMM_chem
  rad = 4.0 * atan(1.0) / 180.0
  coslat   = cos(lat0 * rad)
  coslat!0 = "lat"
  coslat&lat = lat0

  area_MMM_chem = dim_avg_wgt_Wrap(maskMMM_chem, coslat, 1)
  time_MMM_chem = avg(area_MMM_chem)
  delete(area_MMM_chem)
  if (time_MMM_chem .ge. 1000) then
    fmt_MMM_chem = tostring(toint(time_MMM_chem))
  else if (time_MMM_chem .lt. 1000 .and. time_MMM_chem .ge. 100) then
    fmt_MMM_chem = sprintf("%5.1f", time_MMM_chem)
  else if (time_MMM_chem .lt. 100 .and. time_MMM_chem .ge. 10) then
    fmt_MMM_chem = sprintf("%5.2f", time_MMM_chem)
  else if (time_MMM_chem .lt. 10 .and. time_MMM_chem .ge. 1) then
    fmt_MMM_chem = sprintf("%5.3f", time_MMM_chem)
  else
    fmt_MMM_chem = tostring(time_MMM_chem)
  end if
  end if
  end if
  end if

  if ((isatt(diag_script_info, "diff_ref"))) then
    var_grd_ref = XXdata(iref, :, :, 0)
    var_grd_ref&month = ispan(1, 12, 1)
    diff_var_grd = MMM_chem(:, :)
    diff_var_grd = (/MMM_chem - var_grd_ref/)
    diff_maskMMM_chem = diff_var_grd
    diff_maskMMM_chem = (/1./)
    diff_maskMMM_chem = where(ismissing(diff_var_grd), \
                              diff_maskMMM_chem@_FillValue, \
                              diff_maskMMM_chem)
    diff_maskMMM_chem = diff_maskMMM_chem * diff_var_grd
    diff_area_MMM_chem = dim_avg_wgt_Wrap(diff_maskMMM_chem, coslat, 1)
    diff_time_MMM_chem = avg(diff_area_MMM_chem)
    delete(diff_area_MMM_chem)
    if (diff_time_MMM_chem .ge. 1000) then
      diff_fmt_MMM_chem = tostring(toint(diff_time_MMM_chem))
    else if (diff_time_MMM_chem .lt. 1000 .and. \
             diff_time_MMM_chem .ge. 100) then
      diff_fmt_MMM_chem = sprintf("%5.1f", diff_time_MMM_chem)
    else if(diff_time_MMM_chem.lt.100 .and. diff_time_MMM_chem.ge.10) then
      diff_fmt_MMM_chem = sprintf("%5.2f", diff_time_MMM_chem)
    else if (diff_time_MMM_chem.lt.10 .and. diff_time_MMM_chem.ge.1) then
      diff_fmt_MMM_chem = sprintf("%5.3f", diff_time_MMM_chem)
    else
      diff_fmt_MMM_chem = tostring(diff_time_MMM_chem)
    end if
    end if
    end if
    end if
  end if
  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    if (dimsizes(list_nochem_mod) .gt. 1) then
      var_grd_nochem_mod = XXdata(ind_nochem_mod, :, :, 0)
      MMM_nochem = dim_avg_n_Wrap(var_grd_nochem_mod, 0)
      MMM_nochem&month = ispan(1, 12, 1)
      if ((isatt(diag_script_info, "diff_ref"))) then
        diff_var_grd_nochem = MMM_nochem
        diff_var_grd_nochem = (/MMM_nochem - var_grd_ref/)
        diff_var_grd_nochem&month = ispan(1, 12, 1)
        diff_maskMMM_nochem = diff_var_grd_nochem
        diff_maskMMM_nochem = (/1./)
        diff_maskMMM_nochem = where(ismissing(diff_var_grd_nochem), \
                                    diff_maskMMM_nochem@_FillValue, \
                                    diff_maskMMM_nochem)
        diff_maskMMM_nochem = diff_maskMMM_nochem * diff_var_grd_nochem

        diff_area_MMM_nochem = dim_avg_wgt_Wrap(diff_maskMMM_nochem, \
                                                coslat, 1)

        diff_time_MMM_nochem = avg(diff_area_MMM_nochem)
        delete(diff_area_MMM_nochem)
        if (diff_time_MMM_nochem .ge. 1000) then
          diff_fmt_MMM_nochem = tostring(toint(diff_time_MMM_nochem))
        else if (diff_time_MMM_nochem .lt. 1000 .and. \
                 diff_time_MMM_nochem .ge. 100) then
          diff_fmt_MMM_nochem = sprintf("%5.1f", diff_time_MMM_nochem)
        else if(diff_time_MMM_nochem .lt. 100 .and. \
                diff_time_MMM_nochem .ge. 10) then
          diff_fmt_MMM_nochem = sprintf("%5.2f", diff_time_MMM_nochem)
        else if (diff_time_MMM_nochem .lt. 10 .and. \
                 diff_time_MMM_nochem .ge. 1) then
          diff_fmt_MMM_nochem = sprintf("%5.3f", diff_time_MMM_nochem)
        else
          diff_fmt_MMM_nochem = tostring(diff_time_MMM_nochem)
        end if
        end if
        end if
        end if
        end if
        maskMMM_nochem = MMM_nochem
        maskMMM_nochem = (/1./)
        maskMMM_nochem = where(ismissing(MMM_nochem), \
                               maskMMM_nochem@_FillValue, \
                               maskMMM_nochem)
        maskMMM_nochem = maskMMM_nochem * MMM_nochem
        area_MMM_nochem = dim_avg_wgt_Wrap(maskMMM_nochem, coslat, 1)
        time_MMM_nochem = avg(area_MMM_nochem)
        delete(area_MMM_nochem)
        if (time_MMM_nochem .ge. 1000) then
          fmt_MMM_nochem = tostring(toint(time_MMM_nochem))
        else if (time_MMM_nochem .lt. 1000 .and. time_MMM_nochem .ge. 100) then
          fmt_MMM_nochem = sprintf("%5.1f", time_MMM_nochem)
        else if(time_MMM_nochem.lt.100 .and. time_MMM_nochem.ge.10) then
          fmt_MMM_nochem = sprintf("%5.2f", time_MMM_nochem)
        else if (time_MMM_nochem.lt.10 .and. time_MMM_nochem.ge.1) then
          fmt_MMM_nochem = sprintf("%5.3f", time_MMM_nochem)
        else
          fmt_MMM_nochem = tostring(time_MMM_nochem)
        end if
        end if
        end if
        end if
        end if
  end if

  if (multimean) then
    XXdata(dimTot-1, :, :, 0) = (/dim_avg_n(XXdata(index_mod, :, :, 0), 0)/)
    XXdata(dimTot-1, :, :, 1) = (/dim_stddev_n(XXdata(index_mod, :, :, 0), 0)/)
    maskXXdata(dimTot - 1, :, :, 0) = where(ismissing(XXdata(dimTot - 1,
                                                             :, :, 0)), \
                                            maskXXdata@_FillValue, \
                                            maskXXdata(dimTot - 1, :, :, 0))
    plotMod(dimTot - 1) = "MMM"
    years1(dimTot - 1) = max(y1)
    years2(dimTot - 1) = min(y2)
  end if
  X_MMM_std = (/dim_stddev_n(XXdata(index_mod, :, :, :), 0)/)
  maskXXdata = maskXXdata * XXdata
  area_mean = dim_avg_wgt_n_Wrap(maskXXdata, coslat, 1, 2)
  time_mean = dim_avg_n_Wrap(area_mean, 1)
  delete(area_mean)
  if (time_mean(0, 0).ge.1000) then
    fmt = tostring(toint(time_mean))
  else if (time_mean(0, 0).lt.1000 .and. time_mean(0, 0).ge.100) then
    fmt = sprintf("%5.1f", time_mean)
  else if (time_mean(0, 0).lt.100 .and. time_mean(0, 0).ge.10) then
    fmt = sprintf("%5.2f", time_mean)
  else if (time_mean(0, 0).lt.10 .and. time_mean(0, 0).ge.1) then
    fmt = sprintf("%5.3f", time_mean)
  else
    fmt = tostring(time_mean)
  end if
  end if
  end if
  end if

  plot_d = new(dimTot, "graphic")
  plot_d1 = new(dimTot, "graphic")
  map_p = new(7, "graphic")
  ; Plotting

  cmap = (/"White", "Gray1", "Gray75", "maroon4", "DarkOrchid4", \
           "RoyalBlue4", "RoyalBlue", "CadetBlue3", "CadetBlue1",\
           "DarkOliveGreen3", "Chartreuse3", "GreenYellow", "Yellow",\
           "Orange", "OrangeRed", "OrangeRed4", "Black"/)
  gsn_define_colormap(wks1, cmap)
  gsn_define_colormap(wks2, cmap)
  do ipanel = 0, dimTot - 1
    res = True                  ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False
    res@lbLabelBarOn = True
    res@tiMainString = plotMod(ipanel) + " " + years1(ipanel) + "-" \
      + years2(ipanel)
    res@gsnLeftString = var0 + " [" + units + "]"
    res@gsnRightString = fmt(ipanel, 0) + " " + units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023

    res@tiYAxisFont = "helvetica"
    res@tiMainFont = "helvetica"
    res@tmXBLabelFont = "helvetica"
    res@tmYLLabelFont = "helvetica"
    res@lbLabelFont = "helvetica"
    res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
    res@cnLevelSpacingF = Spacing                    ; contour interval
    res@cnMinLevelValF = Min                   ; min level
    res@cnMaxLevelValF = Max                 ; max level
    res@cnLineLabelsOn = False                 ; turn off line labels
    res@cnFillOn = True                  ; turn on color fill
    res@cnInfoLabelOn = False                 ; turn off contour label
    res@gsnSpreadColors = True                  ; use full range of colors
    res@lbLabelAutoStride = True                  ; optimal labels
    res@gsnLeftStringFontHeightF = 0.015
    res@tmXBMode = "Explicit"            ; explicit labels
    res@gsnAddCyclic = True
    res@tmXBValues = (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", \
                       "A", "S", "O", "N", "D", "J"/)    ; labels themselves
    res@tmYLValues = (/ -90., -60., -30., 0., 30., 60., 90./)
    res@tmYLLabels = (/"-90", "-60", "-30", "0", "30", "60", "90"/)
    res@vpHeightF = 0.4                   ; change aspect ratio of plot
    res@vpWidthF = 0.75
    res@lbLabelFontHeightF = 0.011
    res@tiYAxisString = "Latitude"
    plot_d1(ipanel) = gsn_csm_lat_time(wks1, XXdata(line|ipanel, lat|:,\
                                       month|:, stat|0), res)
    delete(res)
  end do
  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  Pres@amJust = "TopLeft"
  if (dimTot.gt.3) then
    gsn_panel(wks1, plot_d1, (/1 + dimTot / 3, 3/), Pres)
  else
    gsn_panel(wks1, plot_d1, (/1, 3/), Pres)
  end if
  frame(wks1)
  ; plot panel multi
  res = True                  ; plot mods desired
  res@gsnDraw = False
  res@gsnFrame = False
  res@lbLabelBarOn = True

  res@gsnRightStringFontHeightF = 0.023
  res@gsnLeftStringFontHeightF = 0.023

  res@tiYAxisFont = "helvetica"
  res@tiMainFont = "helvetica"
  res@tmXBLabelFont = "helvetica"
  res@tmYLLabelFont = "helvetica"
  res@lbLabelFont = "helvetica"
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnLevelSpacingF = Spacing
  res@cnMinLevelValF = Min
  res@cnMaxLevelValF = Max
  res@cnLineLabelsOn = False
  res@cnFillOn = True
  res@gsnSpreadColors = True
  res@lbLabelAutoStride  = True
  res@gsnLeftStringFontHeightF = 0.015
  res@tmXBMode = "Explicit"
  res@gsnAddCyclic = True
  res@tmXBValues = (/0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
  res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S",\
                     "O", "N", "D", "J"/)
  res@tmYLValues = (/-90., -60., -30., 0., 30., 60., 90./)
  res@tmYLLabels = (/"-90", "-60", "-30", "0", "30", "60", "90"/)
  res@vpHeightF = 0.4
  res@vpWidthF = 0.75
  res@lbLabelFontHeightF = 0.011
  res@tiYAxisString = "Latitude"
  ip = 0

  if ((isatt(diag_script_info, "list_chem_mod_string"))) then
    res@tiMainString = diag_script_info@list_chem_mod_string
  else
    res@tiMainString = "MMM_chem"
  end if

  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame

  res@tiMainFontHeightF = 0.03
  res@gsnLeftString = var0 + \
    " [" + units + "]"
  res@gsnRightString = fmt_MMM_chem + " " + units
  res@gsnRightStringFontHeightF = 0.023
  res@gsnLeftStringFontHeightF = 0.023
  MMM_chem&lat = lat0
  ; Plotting
  res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
  res@cnLevelSpacingF      = Spacing                    ; contour interval
  res@cnMinLevelValF       = Min                   ; min level
  res@cnMaxLevelValF       = Max                 ; max level
  map_p(ip) = gsn_csm_lat_time(wks2, \
                               MMM_chem(lat|:, month|:), res)
  ip = ip + 1
  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    if ((isatt(diag_script_info, "list_nochem_mod_string"))) then
      res@tiMainString = nochem_mod_string
    else
      res@tiMainString = "MMM_nochem"
    end if

    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + \
      " [" + units + "]"
    res@gsnRightString = fmt_MMM_nochem + " " + units
    res@gsnLeftStringFontHeightF = 0.023
    res@gsnRightStringFontHeightF = 0.023
    res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
    res@cnLevelSpacingF = Spacing                    ; contour interval
    res@cnMinLevelValF = Min                   ; min level
    res@cnMaxLevelValF = Max                 ; max level

    ; Plotting
    MMM_nochem&lat = lat0
    map_p(ip) = gsn_csm_lat_time(wks2, \
                                 MMM_nochem(lat|:, month|:), res)
    ip = ip + 1
  end if
  if ((isatt(diag_script_info, "diff_ref"))) then
    res@tiMainString = OBS_ref
    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + " [" + units + "]"
    res@gsnRightString = diff_fmt_MMM_chem + " " + units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023
    res@cnLevelSelectionMode = "ExplicitLevels"
    ; Plotting
    var_grd_ref&lat = lat0
    map_p(ip) = gsn_csm_lat_time(wks2, \
                                 var_grd_ref(lat|:, month|:), res)
    ip = ip + 1
  end if
  if ((isatt(diag_script_info, "diff_ref"))) then
    res@tiMainString = OBS_ref
    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + " [" + units + "]"
    res@gsnRightString = fmt(iref, 0) + " " + units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnLevelSpacingF = Spacing
    res@cnMinLevelValF = Min
    res@cnMaxLevelValF = Max

    ; Plotting
    map_p(ip) = gsn_csm_lat_time(wks1, XXdata(line|iref, lat|:, month|:,\
                                 stat|0), res)
    ip = ip + 1
    if ((isatt(diag_script_info, "list_chem_mod_string"))) then
      res@tiMainString = diag_script_info@list_chem_mod_string + " - " \
                         + OBS_ref
    else
      res@tiMainString =  "MMM_chem - " + OBS_ref
    end if
    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + " [" + units + "]"
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnLevelSpacingF = Spacingdiff
    res@cnMinLevelValF = Mindiff
    res@cnMaxLevelValF = Maxdiff
    res@gsnRightString = diff_fmt_MMM_chem + " " + units
    diff_var_grd&lat = lat0
    map_p(ip) = gsn_csm_lat_time(wks2, \
                                 diff_var_grd(lat|:, month|:), res)
    ip = ip + 1
    if ((isatt(diag_script_info, "list_nochem_mod"))) then
      if ((isatt(diag_script_info, "list_nochem_mod_string"))) then
        res@tiMainString = list_nochem_mod_string + " - " \
                           + OBS_ref
      else
        res@tiMainString = "MMM_nochem - " + OBS_ref
      end if

      res@tiMainFontHeightF = 0.03
      res@gsnLeftString = var0 + " [" + units + "]"
      res@gsnRightString = diff_fmt_MMM_chem + " " + units
      res@gsnRightStringFontHeightF = 0.023
      res@gsnLeftStringFontHeightF = 0.023
      diff_var_grd_nochem&lat = lat0
      map_p(ip) = gsn_csm_lat_time(wks2, \
                                   diff_var_grd_nochem(lat|:, month|:), res)

    end if
  end if

  gsn_panel(wks2, map_p, (/4, 2/), Pres)

  frame(wks2)

  month_names = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",\
                  "AUG", "SEP", "OCT", "NOV", "DEC"/)
  nc_filename = config_user_info@work_dir + "eyring13f01_" + var0 + ".nc"
  Xdata = XXdata(line|:, month|:, lat|:, stat|:)
  Xdata@var = var0
  Xdata@diag_script = "eyring13f01"
  Xdata&stat = (/"mean", "stdev"/)
  Xdata&line = plotMod
  delete(Xdata&month)
  Xdata&month = month_names
  Xdata&lat = lat0
  nc_outfile = ncdf_write(Xdata, nc_filename)
  delete(Xdata)
; OUTPUT FILE
  Xdata2 = new((/7, 12, dimsizes(lat0)/), "double")
  mmod = new((/7/), "string")
  nc_filename1 = config_user_info@work_dir + "eyring13f01_MMM" \
    + var0 + ".nc"
  it = 0
  Xdata2(it, :, :) = MMM_chem(month|:, lat|:)
  mmod(it) = "MMM_chem"
  it = it + 1
  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    Xdata2(it, :, :) = MMM_nochem(month|:, lat|:)
    mmod(it) = "MMM_nochem"
    it = it + 1
  end if
  Xdata2(it, :, :) = XXdata(line|iref, month|:, lat|:, stat|0)
  mmod(it) = "OBS_ref"
  it = it + 1
  Xdata2(it, :, :) = diff_var_grd(month|:, lat|:)
  mmod(it) = "diff_MMM_chem"
  it = it + 1
  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    Xdata2(it, :, :) = diff_var_grd_nochem(month|:, lat|:)
    mmod(it) = "diff_MMM_nochem"
    it = it + 1
  end if

  Xdata = new((/it, 12, dimsizes(lat0)/), "double")
  Xdata!0 = "line"
  Xdata!1 = "month"
  Xdata!2 = "lat"
  mmod1 = new(it, "string")
  Xdata = (/Xdata2(:it - 1, :, :)/)
  Xdata@var = var0
  Xdata@diag_script = "eyring13f01"
  mmod1 = mmod(:it - 1)
  Xdata&line = mmod1
  Xdata&month = month_names
  Xdata&lat = lat0
  nc_outfile1 = ncdf_write(Xdata, nc_filename1)
    ; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Annual cycle of " + var0 + " bias", \
                 (/"mean"/), \
                 (/"global"/),\
                 "zonal", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
