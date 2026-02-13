; ###########  ###############################################################
; # ANNUAL CYCLE CONTOUR PLOTS OF REGIONAL MEANS                              #
; #############################################################################
; # Author: Irene Cionni(ENEA, Italy)                                         #
; # CRESCENDO project
; #############################################################################
; #
; # Description
; #    * Creates annual cycle regional plots for individual models and
; #       multi model over the selected project&experiment&mip
; #       - Read regions
; #       - select the regions
; #       - calculate annual cycle for each model
; #       - evaluate the MMM over selected project&experiment&mip
; #       - plot the individual models annual cycle
; #       - plot the MMMs annual cycle and standard deviations
; #
; # Optional diag_script_info attributes:
; #
; #     *  styleset, default "CMIP5"
; #     *  fig938_region_label:
; #       (/"WNA","ENA","CAM","TSA","SSA","EUM","NAF","CAF","SAF",
; #         "NAS","CAS","EAS","SAS","SEA","AUS"/)
; #        these reagions are defined in function select_region
; #        (./diag_scripts/regional_downscaling/regional_function.ncl)
; #        according to a personal comunication of Sin Chan Chou, default "WNA"
; #     *  fig938_project_MMM: projects to average, default "CMIP5"
; #     *  fig938_experiment_MMM: experiments to average, default "historical"
; #     *  fig938_mip_MMM: mip to average, default "Amon"
; #     *  fig938_names_MMM: names in legend, default "CMIP5"
; #     *  fig938_colors_MMM: colors, default "red"
; #     *  fig938_YMin        : minimum Y Axis
; #     *  fig938_YMax        : maximum Y Axis
; #     *  fig938_diff: difference to reference data (true) or
; #                     absolute annual cycle (false), default: true
; #
; #        fig938_mip_MMM, fig938_experiment_MMM, fig938_project_MMM,
; #        fig938_names_MMM, and fig938_colors_MMM
; #        must have the same number of elements i.e.
; #        fig938_project_MMM=(/"CMIP5", "CMIP3"/)
; #        fig938_experiment_MMM=(/"historical", "historical"/)
; #        fig938_mip_MMM=(/"Amon", "Amon"/)
; #        fig938_names_MMM=(/"CMIP5","CMIP3"/)
; #        fig938_colors_MMM:=(/"red","blue"/)
; #
; #
; #
; # Caveats:
; #
; # Comments:
; #      Regions borders are polygons defined in function select_region_938
; #      in regional_function.ncl
; #
; # Modification history:
; # 20230127-weigel_katja: Update for absolute annual cycle including
; #                        new fig938_diff attribute
; #                        and bug fixes in index of data sets/models and
; #                        color for regions
; # 20221026-weigel_katja: header updated
; # 20220314-weigel_katja: header updated
; # original file anncycplot.ncl
; #                                                             #
; ###############################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/regional_downscaling/regional_function.ncl"


begin
  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  units = variable_info[0]@units
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD1 = ListCount(info0)
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")
  mips = metadata_att_as_array(info0, "mip")
  y1 = metadata_att_as_array(info0, "start_year")
  y2 = metadata_att_as_array(info0, "end_year")
  long_name = variable_info[0]@long_name
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
    do imod = 0, dimsizes(names) - 1
      if (names(imod) .eq. refname) then
        iref = imod
      end if
    end do
  else
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset  " + \
              "(variable_info[0]@reference_dataset) needs to be defined.")
  end if
  if(isatt(variable_info, "units")) then
    units = variable_info@units
  else
    units = "K"
  end if
  if(isatt(variable_info, "long_name")) then
    LONG_NAME = variable_info@long_name
  else
    LONG_NAME = var0
  end if

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "fig938_region_label", "WNA")
  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "fig938_MMM", "True")
  set_default_att(diag_script_info, "fig938_project_MMM", "CMIP5")
  set_default_att(diag_script_info, "fig938_experiment_MMM", "historical")
  set_default_att(diag_script_info, "fig938_mip_MMM", "Amon")
  set_default_att(diag_script_info, "fig938_names_MMM", \
                  diag_script_info@fig938_project_MMM)
  set_default_att(diag_script_info, "fig938_colors_MMM", "red")
  set_default_att(diag_script_info, "refModel", "ERA-Interim")
  set_default_att(diag_script_info, "fig938_diff", "True")
  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")
  lin_mar_mod = new((/dimsizes(colors)/), "string")
  lin_mar_mod = (/"Lines"/)

  if (isatt(diag_script_info, "refModel")) then
    idx_ref = ind(names .eq. diag_script_info@refModel)
    if (all(ismissing(idx_ref))) then
      error_msg("f", diag_script, "", "refModel is not included in " + \
                "model list (namelist)")
    end if
    p_ref = "ref_"
  else
    idx_ref = -1
    p_ref = ""
  end if
  region_label = tostring(diag_script_info@fig938_region_label)
  dim_reg = dimsizes(region_label)
  latrange = new((/dimsizes(region_label), 2/), "double")
  lonrange = new((/dimsizes(region_label), 2/), "double")
  do ir = 0, dimsizes(region_label) - 1
    region = select_region1(region_label(ir))
    latrange(ir, :) = region(0:1)
    lonrange(ir, :) = region(2:3)
  end do
  flag_mod = where(projects.ne."OBS" .and.projects.ne."OBS6"\
                   .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if
  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end
begin
  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_38_" + var0)
  wks_r = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_38_regions_" + var0)
  wks_l = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_38_legend_" + var0)
  if (diag_script_info@fig938_MMM.eq."True")
    wks_MMM = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_38_MMM_" + var0)
  end if
  ; ******************************************************************
  ; common resources for  plots
  pan = True
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@vpHeightF = 0.6
  res@vpWidthF = 0.8
  ; res@txFontHeightF = 0.018
  months               = ispan(0, 12, 1)

  res@tiMainFontHeightF      = 0.022
  res@tmXBMode               = "Explicit"
  res@tmXBValues             = ispan(0, 12, 1)
  res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                     "O", "N", "D", "J"/)
  res@tmXBLabelFontHeightF = 0.035
  res@tmYLLabelFontHeightF = 0.035
  res@tiXAxisString = " "
  res@tiYAxisString = " "
  res@xyLineDashSegLenF = 0.5
; ------------res1--------------------------------------------------

  pan1 = True
  res1 = True
  res1@gsnDraw = False
  res1@gsnFrame = False
  ; res1@txFontHeightF          = 0.02
  res1@tiMainFontHeightF = 0.02
  res1@tmXBMode = "Explicit"
  res1@tmXBValues = ispan(0, 12, 1)
  res1@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                      "O", "N", "D", "J"/)
  res1@vpHeightF = 0.6
  res1@vpWidthF = 0.8
  res1@tmXBLabelFontHeightF = 0.035
  res1@tmYLLabelFontHeightF = 0.035
  res1@tiXAxisString = " "
  res1@tiYAxisString = " "
  res1@xyLineDashSegLenF = 0.5
  ; ************************select MMM*********************
  dim_MMM = 1
  if (diag_script_info@fig938_MMM.eq."True")
    if (dimsizes(diag_script_info@fig938_project_MMM).eq.\
        dimsizes(diag_script_info@fig938_experiment_MMM).eq.\
        dimsizes(diag_script_info@fig938_mip_MMM)) then
      dim_MMM = dimsizes(diag_script_info@fig938_project_MMM)
      indic_s = new((/dimsizes(diag_script_info@fig938_project_MMM),\
                    dimsizes(names)/), "integer")
      indic_not_s = new((/dimsizes(diag_script_info@fig938_project_MMM),\
                         dimsizes(names)/), "integer")
      do ik = 0, dimsizes(diag_script_info@fig938_project_MMM) - 1
        ii_s = ind((projects.eq.diag_script_info@fig938_project_MMM(ik)).and.\
                   (exps.eq.diag_script_info@fig938_experiment_MMM(ik)).and.\
                   (mips.eq.diag_script_info@fig938_mip_MMM(ik)))
        indic_s(ik, :dimsizes(ii_s) - 1) = ii_s
        delete(ii_s)

        ii_not_s = \
          ind((projects.ne.diag_script_info@fig938_project_MMM(ik))\
              .or.(exps.ne.diag_script_info@fig938_experiment_MMM(ik))\
              .or.(mips.ne.diag_script_info@fig938_mip_MMM(ik)))
        indic_not_s(ik, :dimsizes(ii_not_s)-1) = ii_not_s
        delete(ii_not_s)
      end do
    else
      error_msg("f", diag_script, "", "diag_script_info@fig938_project_MMM" + \
                "diag_script_info@fig938_experiment_MMM and" + \
                " diag_script_info@fig938_experiment_MMM" + \
                " must have the same dimension")
    end if
  end if
  if (dimsizes(diag_script_info@fig938_project_MMM).gt.1) then
    indic_not_sela = get_unique_difference(indic_s, indic_not_s)
  else
    ndim = ndtooned(indic_not_s)
    ise = ind(.not.ismissing(ndim))
    indic_not_sela = ndim(ise)
    delete(ise)
    delete(ndim)
  end if
  plot = new(dim_reg, "graphic")
  plot_MMM = new(dim_reg, "graphic")
  plot_STD = new((/dim_reg, dim_MMM/), "graphic")

  ; map = new(dimsizes(region_label), "graphic")
  map_r = new(1, "graphic")
  mres_r = True
  mres_r@gsnDraw = False
  mres_r@gsnFrame = False
  mres_r@lbLabelBarOn = False
  mres_r@cnLevelSelectionMode = "ManualLevels"
  mres_r@cnMinLevelValF = 0.
  mres_r@cnMaxLevelValF = 200.0
  mres_r@cnMaxLevelCount = 42.0
  mres_r@cnLevelSpacingF = max((/5.0, 200.0/tofloat(dim_reg + 2.0)/))
  mres_r@cnInfoLabelOn = False
  mres_r@cnLineLabelsOn = False
  mres_r@cnFillOn = True
  mres_r@cnFillOpacityF = 0.5
  cmap_r = read_colormap_file("GMT_wysiwygcont")
  mres_r@cnFillPalette = cmap_r(::-1, :)
  mres_r@mpFillOn = True
  txres               = True
  txres@txFontHeightF = 0.015
  txres@txFont = 22
  txres@txJust        = "CenterLeft"

  random_setallseed(36484749, 9494848)
  unf   = random_uniform(0, 1, (/dim_reg + 2/))
  ip = dim_pqsort(unf, 1)

  do ir = 0, dim_reg - 1
    case1 = region_label(ir)
    if(idx_ref.ne.-1) then
      if (diag_script_info@fig938_MMM.eq."True")
        if ((diag_script_info@fig938_diff.eq."True").and.(idx_ref.ne.-1)) then
          indic_s = where(indic_s.eq.idx_ref, indic_s@_FillValue, indic_s)
          indic_not_sel = where(indic_not_sela.eq.idx_ref,\
                                indic_not_sela@_FillValue, indic_not_sela)
          if(.not.(any(ismissing(indic_not_sel)))) then
            indic_not_sel1 = indic_not_sel(ind(.not.ismissing(indic_not_sel)))
            delete(indic_not_sel)
            indic_not_sel = indic_not_sel1
            delete(indic_not_sel1)
          end if
        else
          indic_not_sel = indic_not_sela
        end if
      end if

      if (diag_script_info@fig938_diff.eq."True") then
        pan@txString = LONG_NAME + " bias vs. " + names(idx_ref)
      else
        pan@txString = LONG_NAME
      end if
      D = new((/dim_MOD1, 13/), "double")
      models = new((/dim_MOD1/), "string")
      C = new((/dim_MOD1 - 1, 13/), "double")
      models1 = new((/dim_MOD1 - 1/), "string")
      colors1 = new((/dim_MOD1 - 1/), typeof(colors))
      dashes1 = new((/dim_MOD1 - 1/), typeof(dashes))
      thicks1 = new((/dim_MOD1 - 1/), typeof(thicks))
      markers1 = new((/dim_MOD1 - 1/), typeof(markers))
      lin_mar_mod1 = new((/dim_MOD1 - 1/), "string")
      ; (b) difference of each model from 'reference'
      ; CASE 1
      A0 = read_data(info0[idx_ref])
      ; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ; Get infos for map plot
      lsdata = A0(0, :, :)
      cell_fill = (/select_region_938(region_label(ir), lsdata)/)
      cell_fill!0 = "lat"
      cell_fill!1 = "lon"
      cell_fill&lat = lsdata&lat
      cell_fill&lon = lsdata&lon
      irc = (ip(ir + 1)) * 200.0/tofloat(dim_reg + 2)

      if (ir.eq.0) then
        cell_fill2 = where(ismissing(cell_fill), -1., irc)
        cell_fill2!0 = "lat"
        cell_fill2!1 = "lon"
        cell_fill2&lat = lsdata&lat
        cell_fill2&lon = lsdata&lon
      else
        cell_fill2 = where(ismissing(cell_fill), cell_fill2, irc)
      end if
      ; *****************************************************************
      Atmp1 = A0
      Atmp = area_operations(Atmp1, latrange(ir, 0), latrange(ir, 1), \
                             lonrange(ir, 0), lonrange(ir, 1), "average",\
                             True)
      delete(A0)
      if (Atmp@units.eq."K") then
        A = convert_units(Atmp, "degC")  ; get Celsius
        UNITS = "~S~o~N~C"
      else
        if (Atmp@units.eq."kg m-2 s-1") then
          A = convert_units(Atmp, "mm/day")
          UNITS = "mm/day"
        else
          A = (/Atmp/)
          UNITS = units
        end if
      end if
      ; end if
      delete(Atmp)
      delete(Atmp1)

      jmod = 0
      do imod = 0, dim_MOD1 - 1
        models(imod) = names(imod)
        if (imod .ne. idx_ref) then
          B0 = read_data(info0[imod])
          Btmp1 = B0
          Btmp = area_operations(Btmp1, latrange(ir, 0), latrange(ir, 1), \
                                 lonrange(ir, 0), lonrange(ir, 1), "average",\
                                 True)
          delete(B0)
          if (Btmp@units.eq."K") then
            B = convert_units(Btmp, "degC")  ; get Celsius
            UNITS = "~S~o~N~C"
          else
            if (Btmp@units.eq."kg m-2 s-1") then
              B = convert_units(Btmp, "mm/day")
              UNITS = "mm/day"
            else
              B = (/Btmp/)
              UNITS = units
            end if
          end if
          delete(Btmp)
          delete(Btmp1)

          ; ---------------------------------------------------------
          ; difference plot of time(x) vs. latitude(y)

          if (diag_script_info@fig938_diff.eq."True") then
            C(jmod, 0:11) = (/A - B/)  ; get diff values
            C(jmod, 12) = (/A(0) - B(0)/)
            D(imod, 0:11) = (/A - B/)
            D(imod, 12) = (/A(0) - B(0)/)
          else
            C(jmod, 0:11) = (/B/)
            C(jmod, 12) = (/B(0)/)
            D(imod, 0:11) = (/B/)
            D(imod, 12) = (/B(0)/)
          end if

          delete(B)
          colors1(jmod) = colors(imod)
          dashes1(jmod) = dashes(imod)
          thicks1(jmod) = thicks(imod)
          markers1(jmod) = markers(imod)
          models1(jmod) = names(imod)
          lin_mar_mod1(jmod) = "Lines"
          jmod = jmod + 1
        end if
      end do  ; imod
    else
      if (diag_script_info@fig938_MMM.eq."True")
        if ((diag_script_info@fig938_diff.eq."True").and.(idx_ref.ne.-1)) then
          indic_s = where(indic_s.eq.idx_ref, indic_s@_FillValue, indic_s)
          indic_not_sel = where(indic_not_sela.eq.idx_ref,\
                                indic_not_sela@_FillValue, indic_not_sela)
          if(.not.(any(ismissing(indic_not_sel)))) then
            indic_not_sel1 = indic_not_sel(ind(.not.ismissing(indic_not_sel)))
            delete(indic_not_sel)
            indic_not_sel = indic_not_sel1
            delete(indic_not_sel1)
          end if
        else
          indic_not_sel = indic_not_sela
        end if
      end if
      lin_mar_mod1 = lin_mar_mod
      C = new((/dim_MOD1, 13/), "double")
      do imod = 0, dim_MOD1 - 1
        ; CASE 1
        pan@txString = LONG_NAME
        A0 = read_data(info0[imod])
        Atmp1 = A0
        Atmp = area_operations(Atmp1, latrange(ir, 0), latrange(ir, 1), \
                               lonrange(ir, 0), lonrange(ir, 1), "average",\
                               True)
        delete(A0)
        if (Atmp@units.eq."K") then
          A = convert_units(Atmp, "degC")  ; get Celsius
          UNITS = "~S~o~N~C"
        else
          if (Atmp@units.eq."kg m-2 s-1") then
            A = convert_units(Atmp, "mm/day")
            UNITS = "mm/day"
          else
            A = (/Atmp/)
            UNITS = units
          end if
        end if
        delete(Atmp)
        delete(Atmp1)
        C(imod, 0:11) = A
        C(imod, 12) = A(0)
        lin_mar_mod1(imod) = "Lines"
      end do

      colors1 = colors
      dashes1 = dashes
      thicks1 = thicks
      markers1 = markers
      models1 = names
    end if

    ; case contour plots of time(x) vs. latitude(y)

    ; res@tiMainString =variable_info@long_name
    if (diag_script_info@fig938_diff.eq."True") then
      E = C
    else
      E = D
      E(idx_ref, 0:11) = (/A/)
      E(idx_ref, 12) = (/A(0)/)
    end if
    if (diag_script_info@fig938_diff.eq."True") then
      Ehlp = reshape(E, (/13*(dim_MOD1 - 1)/))
    else
      Ehlp = reshape(E, (/13*dim_MOD1/))
    end if
    if (any(.not.ismissing(Ehlp))) then
      min1 = min(Ehlp(ind(.not.ismissing(Ehlp))))
      max1 = max(Ehlp(ind(.not.ismissing(Ehlp))))
    else
      min1 = 0.0
      max1 = 1.0
    end if
    if (isatt(diag_script_info, "fig938_YMin")) then
      res@trYMinF = diag_script_info@Figure_938_YMin
    else
      res@trYMinF = min1
    end if
    if (isatt(diag_script_info, "fig938_YMax")) then
      res@trYMaxF = diag_script_info@Figure_938_YMax
    else
      res@trYMaxF = max1
    end if
    if (diag_script_info@fig938_diff.eq."True") then
      res@xyLineColors = colors1  ; line colors
      res@xyLineThicknesses = thicks1        ; line thicknesses
      res@xyDashPatterns = dashes1       ; line patterns
    else
      res@xyLineColors = colors  ; line colors
      res@xyLineThicknesses = thicks        ; line thicknesses
      res@xyDashPatterns = dashes       ; line patterns
    end if

    res@gsnLeftString = case1
    res@gsnLeftStringFontHeightF = 0.04
    res@gsnRightString = UNITS
    res@gsnRightStringFontHeightF = 0.04
    res@tfPolyDrawOrder = "Draw"

    plot(ir) = gsn_csm_xy(wks, months, E, res)

    ; ******************single region  output file**************************

    nc_filename = work_dir + "fig938_" + var0 + ".nc"
    E@var = var0 + "_" + case1
    E@diag_script = "fig938"
    E!0 = "line"
    E!1 = "month"
    month_n = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
               "O", "N", "D", "J"/)

    if (diag_script_info@fig938_diff.eq."True") then
      E&line = models1
    else
      E&line = models
    end if
    E&month = month_n

    if (ir .eq. 0) then
      E@existing = "overwrite"
    else
      E@existing = "append"
    end if

    nc_outfile = ncdf_write(E, nc_filename)

    ; ******************single region  output file**************************

    if (diag_script_info@fig938_MMM.eq."True")
      if ((diag_script_info@fig938_diff.ne."True").and.(idx_ref .ne. -1)) then
        if(.not.all(ismissing(indic_not_sel))) then
          if (any(indic_not_sel.eq.idx_ref)) then
            dim_extr = dimsizes(indic_not_sel)
          else
            dim_extr = dimsizes(indic_not_sel) + 1
          end if
        else
          dim_extr = 1
        end if
      else
        if(.not.all(ismissing(indic_not_sel))) then
          dim_extr = dimsizes(indic_not_sel)
        else
          dim_extr = 0
        end if
      end if

      MMM = new((/dimsizes(diag_script_info@fig938_project_MMM) +\
                dim_extr, 13/), typeof(C))
      STD = new((/dimsizes(diag_script_info@fig938_project_MMM) + \
                dim_extr, 13/), typeof(C))
      names_MMM = new((/dimsizes(diag_script_info@fig938_project_MMM) + \
                      dim_extr/), "string")
      colors_MMM = new((/dimsizes(diag_script_info@fig938_project_MMM) + \
                       dim_extr/), typeof(colors))
      thicks_MMM = new((/dimsizes(diag_script_info@fig938_project_MMM) + \
                       dim_extr/), typeof(thicks))
      dashes_MMM = new((/dimsizes(diag_script_info@fig938_project_MMM) + \
                       dim_extr/), typeof(dashes))
      lin_mar_mod_MMM = new((/dimsizes(diag_script_info@fig938_project_MMM) + \
                            dim_extr/), typeof(lin_mar_mod1))
      colori_def = (/"red", "blue", "green", "magenta", "yellow"/)
      do iin = 0, dimsizes(diag_script_info@fig938_project_MMM) - 1
        aa = indic_s(iin, ind(.not.ismissing(indic_s(iin, :))))
        if (idx_ref .ne. -1)
          D!0 = "case"
          MMM(iin, :) = dim_avg_n_Wrap(D(aa, :), 0)
          STD(iin, :) = dim_stddev_n_Wrap(D(aa, :), 0)
        else
          C!0 = "case"
          MMM(iin, :) = dim_avg_n_Wrap(C(aa, :), 0)
          STD(iin, :) = dim_stddev_n_Wrap(C(aa, :), 0)
        end if
        delete(aa)
        if(isatt(diag_script_info, "fig938_colors_MMM")) then
          colors_MMM(iin) = diag_script_info@fig938_colors_MMM(iin)
        else
          colors_MMM(iin) = colori_def(iin)
        end if
        thicks_MMM(iin) = 3.
        dashes_MMM(iin) = 0
        lin_mar_mod_MMM(iin) = "Lines"
        if(isatt(diag_script_info, "fig938_names_MMM")) then
          names_MMM(iin) = diag_script_info@fig938_names_MMM(iin)
        else
          names_MMM(iin) = diag_script_info@fig938_project_MMM(iin) + " " + \
                           diag_script_info@fig938_experiment_MMM(iin) + \
                           " " + diag_script_info@fig938_mip_MMM(iin)
        end if
      end do
      if(.not.all(ismissing(indic_not_sel))) then
        if (idx_ref .ne. -1)
          MMM(dimsizes(diag_script_info@fig938_project_MMM):, :) = \
            (/E(indic_not_sel(ind(.not.ismissing(indic_not_sel))), :)/)
          colors_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            colors(indic_not_sel(ind(.not.ismissing(indic_not_sel))))
          thicks_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = 3.
          dashes_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            dashes(indic_not_sel(ind(.not.ismissing(indic_not_sel))))
          lin_mar_mod_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            lin_mar_mod(indic_not_sel(ind(.not.ismissing(indic_not_sel))))
          names_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            models(indic_not_sel(ind(.not.ismissing(indic_not_sel))))
        else
          MMM(dimsizes(diag_script_info@fig938_project_MMM):, :) = \
            (/C(indic_not_sel, :)/)
          colors_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            colors1(indic_not_sel)
          thicks_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = 3.
          dashes_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            dashes1(indic_not_sel)
          lin_mar_mod_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            lin_mar_mod1(indic_not_sel)
          names_MMM(dimsizes(diag_script_info@fig938_project_MMM):) = \
            models1(indic_not_sel)
        end if
      end if
    end if
    if(isatt(diag_script_info, "fig938_YMin")) then
      res1@trYMinF = diag_script_info@fig938_YMin
    else
      res1@trYMinF = min1
    end if

    if(isatt(diag_script_info, "fig938_YMax")) then
      res1@trYMaxF = diag_script_info@fig938_YMax
    else
      res1@trYMaxF = max1
    end if

    res1@xyLineColors = colors_MMM  ; line colors
    res1@xyLineThicknesses = thicks_MMM        ; line thicknesses
    res1@xyDashPatterns = dashes_MMM       ; line patterns
    res1@gsnLeftString = case1
    res1@gsnRightString = UNITS
    res1@gsnLeftStringFontHeightF = 0.04
    res1@gsnRightStringFontHeightF = 0.04
    res1@tfPolyDrawOrder = "Draw"
    plot_MMM(ir) = gsn_csm_xy(wks_MMM, months, MMM, res1)
    xp = new((/2 * 13/), typeof(MMM))
    yp = new((/2 * 13, dimsizes(diag_script_info@fig938_project_MMM)/), \
             typeof(MMM))
    do k = 0, 12
      yp(k, :) = MMM(:dimsizes(diag_script_info@fig938_project_MMM) - 1, k) + \
        STD(:dimsizes(diag_script_info@fig938_project_MMM) - 1, k)
      xp(k) = months(k)
      xp(2 * 13 - 1 - k) = months(k)
      yp(2 * 13 - 1 - k, :) = \
        MMM(:dimsizes(diag_script_info@fig938_project_MMM) - 1, k) \
        - STD(:dimsizes(diag_script_info@fig938_project_MMM) - 1, k)
    end do
    do iin = 0, dimsizes(diag_script_info@fig938_project_MMM) - 1
      gsres = True
      gsres@gsFillColor = colors_MMM(iin)
      gsres@gsFillOpacityF = 0.3
      plot_STD(ir, iin) = gsn_add_polygon(wks_MMM, plot_MMM(ir), xp,\
                                          yp(:, iin), gsres)
    end do

    ; ******************single region  output file**************************

    nc_filename_MMM = work_dir + "fig938_" + var0 + "_MMM.nc"
    MMM@var = var0 + "_" + case1
    MMM@diag_script = "fig938"
    month_n = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                     "O", "N", "D", "J"/)
    MMM!0 = "line"
    MMM!1 = "month"
    MMM&line = names_MMM
    MMM&month = month_n

    if (ir .eq. 0) then
      MMM@existing = "overwrite"
    else
      MMM@existing = "append"
    end if

    nc_outfile = ncdf_write(MMM, nc_filename_MMM)

    ; ******************single region  output file**************************

    nc_filename_STD = work_dir + "fig938_" + var0 + "_STD.nc"
    STD@var = var0 + "_" + case1
    STD@diag_script = "fig938"
    month_n = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                     "O", "N", "D", "J"/)
    STD!0 = "line"
    STD!1 = "month"
    STD&line = names_MMM
    STD&month = month_n

    if (ir .eq. 0) then
      STD@existing = "overwrite"
    else
      STD@existing = "append"
    end if

    nc_outfile = ncdf_write(STD, nc_filename_STD)

    ; ******************single region  output file**************************

    delete(min1)
    delete(max1)
    delete(C)
    delete(MMM)
    delete(STD)
    delete(xp)
    delete(yp)
    delete(indic_not_sel)
    delete(A)
    delete(E)
    delete(Ehlp)
  end do  ; region

  ; Panel plots
  pan = True
  pan@gsnMaximize = True
  pan@gsnFrame = False
  pan@gsnPaperOrientation = "portrait"
  n_p = dimsizes(plot) / 5
  if (mod(dimsizes(plot), 5) .ne. 0) then
    n_p = n_p + 1
  end if
  gsn_panel(wks, plot, (/n_p, 5/), pan)
  frame(wks)

  ; ***********add legend****************************
  if (diag_script_info@fig938_diff.eq."True") then
    add_legenda_page(wks_l, models1, colors1, dashes1, lin_mar_mod1,\
                     thicks1, "0.9")
  else
    add_legenda_page(wks_l, models, colors, dashes, lin_mar_mod,\
                     thicks, "0.9")
  end if
  frame(wks_l)

  if (diag_script_info@fig938_MMM.eq."True")
    pan1 = True
    pan1@gsnMaximize = True
    pan1@gsnFrame = False
    pan1@gsnPaperOrientation = "portrait"
    n_p = dimsizes(plot_MMM) / 5
    if (mod(dimsizes(plot_MMM), 5) .ne. 0) then
      n_p = n_p + 1
    end if
    gsn_panel(wks_MMM, plot_MMM, (/n_p, 5/), pan1)
    add_legend(wks_MMM, names_MMM, colors_MMM,\
               dashes_MMM, lin_mar_mod_MMM, thicks_MMM, "0.2")
    ; end if
  end if
  frame(wks_MMM)
  delete(plot_MMM)
  delete(res1)
  delete(pan1)

  map_r = gsn_csm_contour_map(wks_r, cell_fill2, mres_r)

  do ir = 0, dim_reg - 1
    txres@txFontColor   = "black"
    region_array_hlp = select_region1(region_label(ir))
    dum = gsn_add_text(wks_r, map_r, region_label(ir),\
                       region_array_hlp(2) + 5.0, region_array_hlp(0) + 15.0,\
                       txres)
    delete(region_array_hlp)
  end do
  draw(map_r)
  frame(wks_r)

  ; -------------------------------------------------------------
  ; Add provenance
  ; -------------------------------------------------------------

  caption = "Mean seasonal cycle for " + var0 + " over " + \
            "land in different regions"

  log_provenance(nc_filename, \
                 wks@fullname, \
                 caption, \
                 (/"clim"/), \
                 (/"reg"/), \
                 "seas", \
                 (/"cionni_irene"/), \
                 (/"flato13ipcc"/), \
                 metadata_att_as_array(info0, "filename"))

  if (diag_script_info@fig938_MMM .eq. "True")

    log_provenance(nc_filename_MMM, \
                   wks_MMM@fullname, \
                   caption, \
                   (/"clim"/), \
                   (/"reg"/), \
                   "seas", \
                   (/"cionni_irene"/), \
                   (/"flato13ipcc"/), \
                   metadata_att_as_array(info0, "filename"))

  end if


end
