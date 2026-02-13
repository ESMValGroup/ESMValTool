; #############################################################################
; # SEASONAL OR ANNUAL BIAS                                                   #
; #############################################################################
; # Author: Irene Cionni(ENEA, Italy)                                         #
; # CRESCENDO project
; #############################################################################
;
; # Description:                                                              #
; #    * Create seasonal or annual bias box plots at selected regions         #
; #       - Read regions
; #        -create mask map using region poligons
; #       - for all the models calculate the bias reference_dataset
; #       - sort biases for each selected project&experiment&mip
; #       - evaluate percentiles (5th 25th 50th 75th 95th)
; #       - plot the first selected project&experiment&mip as box-and-whisker
; #       - plot the others selected projects&experiment&mip as markers
; #       - compare the root-mean-square error of first project&experiment&mip
; #         with the second if exist labeled with red for the regions where
; #         the first error is larger than the second.
; #
; # Required variable_info attributes (variable specific):
; #
; #    *  reference_dataset: reference dataset name
; #
; # Optional diag_script_info attributes:
; #
; #    *  styleset, default "CMIP5"
; #    *  fig939_season : seasons i.e. (/"DJF","JJA","ANN"/)
; #    *  fig939_region_label :
; #      (/"ALA","CGI","WNA","CNA","ENA","CAM","AMZ","NEB","WSA","SSA","NEU",
; #       "CEU","MED","SAH","WAF","EAF","SAF","NAS","WAS","CAS","TIB","EAS",
; #       "SAS","SEA","NAU","SAU"/)
; #        these region are defined in Seneviratne et al., 2012: Appendix 3.A
; #        Notes and technical details on Chapter 3 figures.
; #       In: Managing the Risks of Extreme Events and Disasters to Advance
; #          Climate Change Adaptation(IPCC). <www.ipcc.ch>.
; #      * fig939_MMM, default: True
; #      * fig939_project_MMM: projects to average
; #      * fig939_experiment_MMM: experiments to average
; #      * fig939_mip_MMM: mip to average
; #      * Dimensions of fig939_project_MMM, fig939_experiment_MMM and
; #         fig939_mip_MMM must have the number of elements i.e.
; #         fig939_project_MMM=(/"CMIP5", "CMIP3"/)
; #         fig939_experiment_MMM=(/"historical", "historical"/)
; #         fig939_mip_MMM=(/"Amon", "Amon"/)
; #
; #      * Names for legend
; #          fig939_names_MMM   :names in legend  i.e. (/"CMIP5","CMIP3"/
; #      * Vertical lines divide groups of region
; #         fig939_vert_line_pos
; #          i.e. (/6,10,13,17,24,26/)
; #      * labels of vertical lines
; #        fig939_vert_line_label
; #        e.g. (/"North America","South America","Europe","Africa",
; #                "Asia","Australia"/)
; #        fig939_mode    : True= cumulative mode
; #        fig939_YMin        : minimum Y Axis
; #        fig939_YMax        : maximum Y Axis
; #        fig939_percentage  : default: False
; # Caveats:
; #      This script requires NCL version 6.40
; #
; # Comments:
; #      Regions borders are polygons defined in function select_region_srex
; #      in regional_function.ncl
; #
; # Modification history:
; # 20230110-weigel_katja: fixed calendar, calendar starting point needs to
; #                        including time: "months since 1-1-1 00:00:00"
; #                        and converting month number to time Axis
; #                        with -1.0
; # 20230111-weigel_katja: added option for percentage difference
; # 20230110-weigel_katja: fixed calendar, calendar starting point needs to
; #                        including time: "months since 1-1-1 00:00:00"
; #                        and converting month number to time Axis
; #                        with -1.0
; # 20221026-weigel_katja: header updated
; # 20221026-weigel_katja: fixed calendar, time_operations expects
; #                        "gregorian" instead of "Gregorian"
; # 20221026-weigel_katja: removed unused variable fig939_colors_MMM
; # 20221026-weigel_katja: removed unused variable refModel
; #                        (reference_dataset instead)
; # 20220314-weigel_katja: header updated
; ###############################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/regional_downscaling/regional_function.ncl"


begin
  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  UNITS = variable_info[0]@units
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
        idx_ref = imod
      end if
    end do
  else
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset  " + \
              "(variable_info[0]@reference_dataset) needs to be defined.")
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

  set_default_att(diag_script_info, "fig939_region_label", "ALAs")
  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "fig939_MMM", "True")
  set_default_att(diag_script_info, "fig939_project_MMM", "CMIP5")
  set_default_att(diag_script_info, "fig939_experiment_MMM", "historical")
  set_default_att(diag_script_info, "fig939_mip_MMM", "Amon")
  ; set_default_att(diag_script_info, "fig939_colors_MMM", "red")
  set_default_att(diag_script_info, "fig939_names_MMM", \
                  diag_script_info@fig939_project_MMM)
  set_default_att(diag_script_info, "fig939_season", "DJF")
  set_default_att(diag_script_info, "fig939_percentage", "False")

  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")
  lin_mar_mod = new((/dimsizes(colors)/), "string")
  lin_mar_mod = (/"Lines"/)
  season = diag_script_info@fig939_season
  mode = new(1, logical)
  if(isatt(diag_script_info, "fig939_mode")) then
    mode = True
  else
    mode = False
  end if
  region_label = tostring(diag_script_info@fig939_region_label)
  dim_reg = dimsizes(region_label)
  dim_seas = dimsizes(season)

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
  dim_MMM = 1

  project_MMM = diag_script_info@fig939_project_MMM
  exp_MMM = diag_script_info@fig939_experiment_MMM
  mip_MMM = diag_script_info@fig939_mip_MMM
  if (diag_script_info@fig939_MMM.eq."True")
    if (dimsizes(project_MMM).eq.dimsizes(exp_MMM).eq.dimsizes(mip_MMM)) then
      dim_MMM = dimsizes(project_MMM)
      indic_s = new((/dimsizes(project_MMM),\
                    dimsizes(names)/), "integer")
      indic_not_s = new((/dimsizes(project_MMM),\
                        dimsizes(names)/), "integer")
      do ik = 0, dimsizes(diag_script_info@fig939_project_MMM) - 1
        ii_s = ind((projects.eq.project_MMM(ik))\
                   .and.(exps.eq.exp_MMM(ik))\
                   .and.(mips.eq.mip_MMM(ik)))
        indic_s(ik, :dimsizes(ii_s) - 1) = ii_s
        delete(ii_s)
        ii_not_s = ind((projects.ne.project_MMM(ik)).or.(exps.ne.exp_MMM(ik))\
                       .or.(mips.ne.mip_MMM(ik)))
        indic_not_s(ik, :dimsizes(ii_not_s) - 1) = ii_not_s
        delete(ii_not_s)
      end do
    else
      error_msg("f", diag_script, "", "diag_script_info@fig939_project_MMM" + \
                "diag_script_info@fig939_experiment_MMM and" + \
                " diag_script_info@fig939_experiment_MMM" + \
                " must have the same dimension")
    end if
  end if

  if (dimsizes(project_MMM).gt.1) then
    indic_not_sela = get_unique_difference(indic_s, indic_not_s)
  else
    ndim = ndtooned(indic_not_s)
    ise = ind(.not.ismissing(ndim))
    indic_not_sela = ndim(ise)
    delete(ise)
    delete(ndim)
  end if
  MMM = new((/dim_seas, dim_reg, dimsizes(project_MMM) +\
            dimsizes(indic_not_sela), 6/), "double")
  names_MMM = new((/dimsizes(project_MMM)/), "string")
  MMM_rmse = new((/dim_seas, dim_reg, dimsizes(project_MMM)/), "double")
  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_39_" + var0)
  wks_r = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_39_regions_" + var0)
  wks_l = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_39_legend_" + var0)
  ; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Amask = read_data(info0[idx_ref])
  lsdata = Amask(0, :, :)
  lsdata = where(ismissing(lsdata), 1., lsdata)
  ; *****************************************************************
  map_r = new(1, "graphic")
  plot = new(dim_seas, "graphic")
  dum_ref = new(dim_seas, "graphic")
  dum1 = new((/dim_seas, dimsizes(project_MMM)/),\
             "graphic")
  dum2 = new((/dim_seas, dimsizes(project_MMM)/),\
             "graphic")
  dum3 = new((/dim_seas, dimsizes(project_MMM)/),\
             "graphic")
  dum4 = new((/dim_seas, dimsizes(project_MMM)/),\
             "graphic")
  dum5 = new((/dim_seas, dimsizes(project_MMM)/),\
             "graphic")
  dum_lines = new((/dim_seas, 30/), "graphic")
  dum_txt = new((/dim_seas, 30/), "graphic")
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
  gsres                = True
  gsres@gsFillOpacityF = 0.7
  txres               = True
  txres@txFontHeightF = 0.012
  txres@txFont = 22
  txres@txJust        = "CenterLeft"

  random_setallseed(36484749, 9494848)
  unf   = random_uniform(0, 1, (/dim_reg + 2/))
  ip = dim_pqsort(unf, 1)

  do iseas = 0, dim_seas - 1
    do ir = 0, dim_reg - 1
      case1 = region_label(ir)
      cell_fill = (/select_region_srex(region_label(ir), lsdata)/)
      cell_fill!0 = "lat"
      cell_fill!1 = "lon"
      cell_fill&lat = lsdata&lat
      cell_fill&lon = lsdata&lon

      if (iseas.eq.0) then
        irc = (ip(ir + 1)) * 200.0/tofloat(dim_reg + 2.0)

        if (ir.eq.0) then
          cell_fill2 = where(ismissing(cell_fill), -1., irc)
          cell_fill2!0 = "lat"
          cell_fill2!1 = "lon"
          cell_fill2&lat = lsdata&lat
          cell_fill2&lon = lsdata&lon
        else
          cell_fill2 = where(ismissing(cell_fill), cell_fill2, irc)
        end if
      end if
      cell_fill = where(ismissing(cell_fill), -1., 1.)
      if (diag_script_info@fig939_MMM.eq."True")
        indic_s = where(indic_s.eq.idx_ref, indic_s@_FillValue, indic_s)
        indic_not_sel = where(indic_not_sela.eq.idx_ref,\
                              indic_not_sela@_FillValue, indic_not_sela)
        if(.not.ismissing(all(indic_not_sel))) then
          indic_not_sel1 = indic_not_sel(ind(.not.ismissing(indic_not_sel)))
          delete(indic_not_sel)
          indic_not_sel = indic_not_sel1
          delete(indic_not_sel1)
        end if
      end if
      D = new((/dim_MOD+1/), "double")
      C = new((/dim_MOD/), "double")
      models1 = new((/dim_MOD/), "string")
      colors1 = new((/dim_MOD/), typeof(colors))
      dashes1 = new((/dim_MOD/), typeof(dashes))
      thicks1 = new((/dim_MOD/), typeof(thicks))
      markers1 = new((/dim_MOD/), typeof(markers))
      lin_mar_mod1 = new((/dim_MOD/), "string")
      ; (b) difference of each model from 'reference'
      ; CASE 1 In this plot only case1 exsist!!!!!
      ; case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
      A0 = read_data(info0[idx_ref])
      if (season(iseas).ne."ANN") then
        if (mode) then
          A0!0 = "time"
          time1 = todouble(A0&time - 1.0)
          time1!0 = "time"
          ; ys=input_file_info[idx_ref]@start_year
          time1@units = "months since 1-1-1 00:00:00"
          time1@calendar = "gregorian"
          delete(A0&time)
          A0&time = time1
          delete(time1)
          Atmp3 = time_operations(A0, -1, -1, "extract", season(iseas), True)
          date3 = cd_calendar(Atmp3&time, 0)
          year3 = date3(:, 0)
          month3 = date3(:, 1)
          weights3 = days_in_month(toint(year3), toint(month3))
          Atmp4 = Atmp3
          Atmp4 = (/Atmp3 * conform(Atmp3, weights3, 0)/)
          Atmp1 = dim_sum_n_Wrap(Atmp4, 0)
          delete(Atmp3)
          delete(Atmp4)
          delete(date3)
          delete(year3)
          delete(month3)
          delete(weights3)
        else
          A0!0 = "time"
          time1 = todouble(A0&time - 1.0)
          time1!0 = "time"
          ; ys=input_file_info[idx_ref]@start_yea
          time1@units = "months since 1-1-1 00:00:00"
          time1@calendar = "gregorian"
          delete(A0&time)
          A0&time = time1
          Atmp1 = time_operations(A0, -1, -1, "average", season(iseas), True)
        end if
      else
        if (mode) then
          A0!0 = "time"
          time1 = todouble(A0&time - 1.0)
          time1!0 = "time"
          ; ys=input_file_info[idx_ref]@start_year
          time1@units = "months since 1-1-1 00:00:00"
          time1@calendar = "gregorian"
          delete(A0&time)
          A0&time = time1
          Atmp3 = A0
          date3 = cd_calendar(Atmp3&time, 0)
          year3 = date3(:, 0)
          month3 = date3(:, 1)
          weights3 = days_in_month(toint(year3), toint(month3))
          Atmp4 = Atmp3
          Atmp4 = (/Atmp3 * conform(Atmp3, weights3, 0)/)
          Atmp1 = dim_sum_n_Wrap(Atmp4, 0)
          delete(Atmp3)
          delete(Atmp4)
          delete(date3)
          delete(year3)
          delete(month3)
          delete(weights3)
        else
          A0!0 = "time"
          time1 = todouble(A0&time - 1.0)
          time1!0 = "time"
          ; ys=input_file_info[idx_ref]@start_year
          time1@units = "months since 1-1-1 00:00:00"
          time1@calendar = "gregorian"
          delete(A0&time)
          A0&time = time1
          Atmp1 = time_operations(A0, -1, -1, "average", "annualclim", True)
        end if
      end if
      ; ******************mask region*********************
      Atmp1 = mask(Atmp1, cell_fill.eq.-1, False)
      ; **************************************************
      Atmp = area_operations(Atmp1, -90., 90., 0., 360., "average", True)
      delete(A0)
      A = Atmp
      delete(Atmp)
      delete(Atmp1)
      jmod = 0
      do imod = 0, dim_MOD - 1
        if (imod .ne. idx_ref) then
          B0 = read_data(info0[imod])
          if (season(iseas).ne."ANN") then
            if (mode) then
              B0!0 = "time"
              time1 = todouble(B0&time - 1.0)
              time1!0 = "time"
              ; ys = input_file_info[idx_ref]@start_year
              time1@units = "months since 1-1-1 00:00:00"
              time1@calendar = "gregorian"
              delete(B0&time)
              B0&time = time1
              Btmp3 = time_operations(B0, -1, -1, "extract",\
                                      season(iseas), True)

              date3 = cd_calendar(Btmp3&time, 0)
              year3 = date3(:, 0)
              month3 = date3(:, 1)
              weights3 = days_in_month(toint(year3), toint(month3))
              Btmp4 = Btmp3
              Btmp4 = (/Btmp3 * conform(Btmp3, weights3, 0)/)
              Btmp1 = dim_sum_n_Wrap(Btmp4, 0)
              delete(Btmp3)
              delete(Btmp4)
              delete(date3)
              delete(year3)
              delete(month3)
              delete(weights3)
            else
              B0!0 = "time"
              time1 = todouble(B0&time - 1.0)
              time1!0 = "time"
              ; ys=input_file_info[idx_ref]@start_year
              time1@units = "months since 1-1-1 00:00:00"
              time1@calendar = "gregorian"
              delete(B0&time)
              B0&time = time1
              Btmp1 = time_operations(B0, -1, -1, "average",\
                                      season(iseas), True)
            end if
          else
            if (mode) then
              Btmp3 = B0
              Btmp3!0 = "time"
              time1 = todouble(Btmp3&time - 1.0)
              time1!0 = "time"
              time1@units = "months since 1-1-1 00:00:00"
              time1@calendar = "gregorian"
              delete(Btmp3&time)
              Btmp3&time = time1
              date3 = cd_calendar(Btmp3&time, 0)
              year3 = date3(:, 0)
              month3 = date3(:, 1)
              weights3 = days_in_month(toint(year3), toint(month3))
              Btmp4 = Btmp3
              Btmp4 = (/Btmp3 * conform(Btmp3, weights3, 0)/)
              Btmp1 = dim_sum_n_Wrap(Btmp4, 0)
              delete(Btmp3)
              delete(Btmp4)
              delete(date3)
              delete(year3)
              delete(month3)
              delete(weights3)
            else
              B0!0 = "time"
              time1 = todouble(B0&time - 1.0)
              time1!0 = "time"
              ; ys = input_file_info[idx_ref]@start_year
              time1@units = "months since 1-1-1 00:00:00"
              time1@calendar = "gregorian"
              delete(B0&time)
              B0&time = time1
              Btmp1 = time_operations(B0, -1, -1, "average", "annualclim",\
                                      True)
            end if
          end if
          ; ******************mask region**********************
          Btmp1 = mask(Btmp1, cell_fill.eq.-1, False)
          ; **************************************************
          Btmp = area_operations(Btmp1, -90., 90., 0., 360., "average", True)
          delete(B0)
          B = Btmp
          delete(Btmp)
          delete(Btmp1)
          ; --------------------------------------------------------
          ; Bias
          if (mode) then
            print("cumulative mode")
            C(jmod) = (/ 100 * (avg(B) - avg(A)) / avg(A)/)
            D(imod) = (/100 * (avg(B) - avg(A)) / avg(A)/)
          else
            C(jmod) = (/avg(B) - avg(A)/)  ; get diff values
            D(imod) = (/avg(B) - avg(A)/)
          end if

          if (UNITS.eq."K") then
            UNITS = "~S~o~N~C"
          else
            if (mode) then
              UNITS = "%"
            else
              UNITS = variable_info[0]@units
            end if
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
      delete(A)
      if (diag_script_info@fig939_MMM.eq."True")
        colori_def = (/"red", "orange", "green", "blue", "purple"/)
        scatters_def = (/0, 5, 16, 4, 7, 8, 12/)
        do iin = 0, dimsizes(diag_script_info@fig939_project_MMM) - 1
          if (.not.all(ismissing(indic_s(iin, :)))) then
            aa = indic_s(iin, ind(.not.ismissing(indic_s(iin, :))))
          end if
          F = D(aa)
          qsort(F)
          dimt = dimsizes(F)
          x5 = round(.05 * dimt, 3) - 1
          x25 = round(.25 * dimt, 3) - 1
          x50 = round(.50 * dimt, 3) - 1
          x75 = round(.75 * dimt, 3) - 1
          x95 = round(.95 * dimt, 3) - 1
          x5 = where(x5.lt.0, 0, x5)
          x25 = where(x25.lt.0, 0, x25)
          x50 = where(x50.lt.0, 0, x50)
          x75 = where(x75.lt.0, 0, x75)
          x95 = where(x95.lt.0, 0, x95)
          MMM_rmse(iseas, ir, iin) = (sum(F ^ 2) / dimt) ^ 0.5
          MMM(iseas, ir, iin, :) = (/F(x5), F(x25), F(x50),\
                                    F(x75), F(x95), stddev(F)/)
          delete(aa)
          delete(F)
          delete(dimt)

          ; if(isatt(diag_script_info, "fig939_markers_MMM")) then
          ;   markers_MMM(iin) = diag_script_info@fig939_markers_MMM(iin)
          ; else
          ;   markers_MMM(iin) = scatters_def(iin)
          ; end if

          ; markers_MMM(iin) = 0

          if(isatt(diag_script_info, "fig939_names_MMM")) then
            names_MMM(iin) = diag_script_info@fig939_names_MMM(iin)
          else
            names_MMM(iin) = diag_script_info@fig939_project_MMM(iin) + " " +\
                             diag_script_info@fig939_experiment_MMM(iin) +\
                             " " + diag_script_info@fig939_mip_MMM(iin)
          end if
        end do
        if (.not.ismissing(all(indic_not_sela))) then
          MMM(iseas, ir, dimsizes(diag_script_info@fig939_project_MMM):, 0) =\
           (/D(indic_not_sela)/)
        end if
      end if
      delete(C)
      delete(D)
    end do
    res = True                         ; plot mods desired
    ; add * if first project root mean square error is larger than the second
    region_label1 = region_label
    colors_reg = new(dimsizes(region_label), "string")
    do irg = 0, dimsizes(region_label) - 1
      if (dim_MMM.gt.1) then
        if (MMM_rmse(iseas, irg, 0).gt.MMM_rmse(iseas, irg, 1)) then
          region_label1(irg) = region_label(irg) + "*"
          colors_reg(irg) = "red"
        else
          region_label1(irg) = region_label(irg)
          colors_reg(irg) = "blue"
        end if
      else
        colors_reg(irg) = "blue"
      end if
    end do

    ind_blue = ind(colors_reg.eq."blue")
    ind_red = ind(colors_reg.eq."red")
    res@ind_blue = ind_blue
    res@ind_red = ind_red
    res@tmXBLabels = region_label1
    res@tmXBLabelFontHeightF = 0.017
    res@tiMainString  = "bias, " + long_name + \
                        " ("+UNITS+"), " + season(iseas)
    if(isatt(diag_script_info, "fig939_YMin")) then
      res@trYMinF = diag_script_info@fig939_YMin
    else
      res@trYMinF = min(MMM) - min(MMM) / 5
    end if

    if (isatt(diag_script_info, "fig939_YMax")) then
      res@trYMaxF = diag_script_info@fig939_YMax
    else
      res@trYMaxF = max(MMM) + max(MMM) / 5
    end if
    res@tmXBLabelAngleF = 90.
    res@vpHeightF = 0.5
    res@vpWidthF = 0.9
    ; res@gsnYRefLine=0.0
    ; res@gsnYRefLineColor="gray11"
    ; res@gsnYRefLineDashPattern=1.
    plot(iseas) = box_plot_cr(wks, ispan(1, dim_reg, 1),\
                              MMM(iseas, :, 0, :), False, res, False)
    if (isatt(diag_script_info, "fig939_vert_line_pos")) then
      lines_x = diag_script_info@fig939_vert_line_pos
      do ilin = 0, dimsizes(lines_x) - 1
        dres = True
        dres@gsLineColor = "gray11"
        dres@gsLineDashPattern = 2.
        dum_lines(iseas, ilin) = gsn_add_polyline(wks, plot(iseas),\
                                                  (/lines_x(ilin) + 0.5,\
                                                    lines_x(ilin) + 0.5/),\
                                                  (/min(MMM) - avg(MMM) / 10,\
                                                    max(MMM) + \
                                                    avg(MMM) / 10/),\
                                                  dres)
      end do
      if (isatt(diag_script_info, "fig939_vert_line_label")) then
        lab_x = diag_script_info@fig939_vert_line_label

        do itxt = 0, dimsizes(lab_x) - 1
          txres = True
          txres@txFontHeightF = 0.02
          txres@txAngleF = 90.
          txres@txJust = "BottomRight"
          dum_txt(iseas, itxt) = gsn_add_text(wks, plot(iseas), lab_x(itxt),\
                                              lines_x(itxt) + 0.4,\
                                              min(MMM) - min(MMM) / 4,\
                                              txres)
        end do
      end if
    end if
    rres = True
    rres@gsLineColor = "gray11"
    rres@gsLineDashPattern = 1.
    dum_ref(iseas) = gsn_add_polyline(wks, plot(iseas), (/0., 100./),\
                                      (/0., 0./), rres)
    do iin = 1, dimsizes(diag_script_info@fig939_project_MMM) - 1
      mres = True                     ; marker mods desired
      mres@gsMarkerIndex = scatters_def(iin)         ; polymarker style
      mres@gsMarkerSizeF = 0.008 * 26 / dim_reg
      mres@gsMarkerThicknessF = 2.
      mres@gsMarkerColor = "purple"                    ; polymarker color
      dum1(iseas, iin) = gsn_add_polymarker(wks, plot(iseas),\
                                            ispan(1, dim_reg, 1),\
                                            MMM(iseas, :, iin, 0), mres)
      mres@gsMarkerColor = "blue"                    ; polymarker color
      dum2(iseas, iin) = gsn_add_polymarker(wks, plot(iseas),\
                                            ispan(1, dim_reg, 1),\
                                            MMM(iseas, :, iin, 1),\
                                            mres)
      mres@gsMarkerColor = "green"                  ; polymarker color
      dum3(iseas, iin) = gsn_add_polymarker(wks, plot(iseas),\
                                            ispan(1, dim_reg, 1), \
                                            MMM(iseas, :, iin, 2), \
                                            mres)
      mres@gsMarkerColor = "orange"                ; polymarker color
      dum4(iseas, iin) = gsn_add_polymarker(wks, plot(iseas),\
                                            ispan(1, dim_reg, 1),\
                                            MMM(iseas, :, iin, 3),\
                                            mres)
      mres@gsMarkerColor = "red"                    ; polymarker color
      dum5(iseas, iin) = gsn_add_polymarker(wks, plot(iseas), \
                                            ispan(1, dim_reg, 1),\
                                            MMM(iseas, :, iin, 4), mres)
    end do
    delete(ind_blue)
    delete(ind_red)
    delete(res@ind_blue)
    delete(res@ind_red)
  end do

  ; ******************single region  output file**************************

  print(work_dir)
  nc_filename_MMM = work_dir + "fig939_" + var0 + "_MMM.nc"
  MMM@var = var0
  MMM@diag_script = "fig939"
  MMM!0 = "season"
  MMM!1 = "region"
  MMM!2 = "line"
  MMM!3 = "stat"
  dsizes_x = dimsizes(MMM)
  nline     = dsizes_x(2)
  MMM&line = new((/nline/), "string")  ; (/"None", "None", "None"/)
  MMM&season = season
  MMM&region = region_label
  MMM&stat = (/"x5", "x25", "x50", "x75", "x95", "stddev"/)
  nc_outfile = ncdf_write(MMM, nc_filename_MMM)

  pan = True
  pan@gsnMaximize = True
  pan@gsnFrame = False
  pan@gsnPaperOrientation = "portrait"
  n_p = dimsizes(plot) / 3
  if (mod(dimsizes(plot), 3) .ne. 0) then
    n_p = n_p + 1
  end if
  gsn_panel(wks, plot, (/3, n_p/), pan)
  frame(wks)
  res1 = True                         ; plot mods desired
  res1@tmXTLabels = names_MMM
  res1@tmXTValues = ispan(1, dimsizes(diag_script_info@fig939_project_MMM)\
                          + 1, 1)
  res1@vpHeightF = 0.5
  res1@vpWidthF = 0.5
  lg_y = new((/1, 6/), "float")
  lg_y(0, 0:4) = (/5., 25., 50., 75., 95./)
  plot_lg = box_plot_lg(wks_l, 1, lg_y, False, res1, False)
  duml1 = new(dimsizes(project_MMM), "graphic")
  duml2 = new(dimsizes(project_MMM), "graphic")
  duml3 = new(dimsizes(project_MMM), "graphic")
  duml4 = new(dimsizes(project_MMM), "graphic")
  duml5 = new(dimsizes(project_MMM), "graphic")
  do iin = 1, dimsizes(project_MMM) - 1
    mresl = True                     ; marker mods desired
    mresl@gsMarkerIndex = scatters_def(iin)             ; polymarker style
    ; mresl@gsMarkerSizeF = 60./dim_reg                   ; polymarker size
    mresl@gsMarkerSizeF = 0.025 * 26 / dim_reg
    mresl@gsMarkerThicknessF = 3.
    mresl@gsMarkerColor = "purple"                    ; polymarker color
    duml1(iin) = gsn_add_polymarker(wks_l, plot_lg, iin + 1, lg_y(0, 0), mresl)
    mresl@gsMarkerColor = "blue"                    ; polymarker color
    duml2(iin) = gsn_add_polymarker(wks_l, plot_lg, iin + 1, lg_y(0, 1), mresl)
    mresl@gsMarkerColor = "green"                    ; polymarker color
    duml3(iin) = gsn_add_polymarker(wks_l, plot_lg, iin + 1, lg_y(0, 2), mresl)
    mresl@gsMarkerColor = "orange"                    ; polymarker color
    duml4(iin) = gsn_add_polymarker(wks_l, plot_lg, iin + 1, lg_y(0, 3), mresl)
    mresl@gsMarkerColor = "red"                    ; polymarker color
    duml5(iin) = gsn_add_polymarker(wks_l, plot_lg, iin + 1, lg_y(0, 4), mresl)
  end do
  draw(plot_lg)
  delete(pan)
  delete(n_p)
  frame(wks_l)
  pan = True
  pan@gsnMaximize = True
  pan@gsnFrame = False
  pan@gsnPaperOrientation = "portrait"
  map_r = gsn_csm_contour_map(wks_r, cell_fill2, mres_r)

  do ir = 0, dim_reg - 1
    txres@txFontColor   = "black"
    cregion_array = (/select_region_srex_poly(region_label(ir))/)
    dum = gsn_add_text(wks_r, map_r, region_label(ir),\
                       cregion_array(0), cregion_array(1), txres)
    delete(cregion_array)
  end do
  draw(map_r)
  frame(wks_r)
;  info_output("wrote " + output_dir + output_filename, verbosity, 1)
;  info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

  ; -------------------------------------------------------------
  ; Add provenance
  ; -------------------------------------------------------------

  caption = "Seasonal- and annual mean biases of " + var0 + " over " + \
            "land in different regions"

  log_provenance(nc_filename_MMM, \
                 wks@fullname, \
                 caption, \
                 (/"diff", "mean", "perc"/), \
                 (/"reg"/), \
                 "box", \
                 (/"cionni_irene"/), \
                 (/"flato13ipcc", "seneviratne12ipcc"/), \
                 metadata_att_as_array(info0, "filename"))

end
