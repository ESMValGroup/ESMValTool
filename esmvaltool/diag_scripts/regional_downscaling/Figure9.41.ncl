; #############################################################################
; #############################################################################
; # Author: Irene Cionni (ENEA, Italy)                                        #
; # ESMVal project
; #############################################################################
; #    * Creates plot ranked monthly data of observations vs models.          #
; #  3 pages are created:
; #  1st page Centred (data are adjusted to get a zero mean in model data with
; #         respect to the diagonal)
; #  2nd page  Uncentred (data not adjusted)
; #  3th page  Full models legend
; #
; # Required variable_info attributes (variable specific):
; #     *Reference Dataset
; #      variable_info[0]@reference_dataset
; # Optional diag_script_info attributes (diagnostics specific):
; #     *Region Label
; #      diag_script_info@fig941_region_label, default "MEDs"
; #     *styleset, default "CMIP5"
; #
; #Optional diag_script_info attributes
; #                                                                           #
; #                                                                           #
; # Caveats:                                                                  #
; #                                                                           #
; # Modification history:
; # 20221026-weigel_katja: header updated
; # 20220314-weigel_katja: header updated
; #############################################################################
; #############################################################################
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

  set_default_att(diag_script_info, "fig941_region_label", "MEDs")
  set_default_att(diag_script_info, "styleset", "CMIP5")
  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")
  lin_mar_mod = new((/dimsizes(colors)/), "string")
  lin_mar_mod = (/"Lines"/)
  region_label = tostring(diag_script_info@fig941_region_label)
  dim_reg = dimsizes(region_label)
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
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_41_" + var0)
  wks_r = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_41_raw_" + var0)
  wks_l = get_wks("dummy_for_wks", DIAG_SCRIPT, "Figure9_41_legend_" + var0)
  ; +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Amask = read_data(info0[idx_ref])
  lsdata = Amask(0, :, :)
  ; *****************************************************************
  plot = new(dimsizes(region_label), "graphic")
  plot1 = new(dimsizes(region_label), "graphic")
  dum_diag = new(dimsizes(region_label), "graphic")
  dum_diag1 = new(dimsizes(region_label), "graphic")
  do ir = 0, dim_reg - 1
    case1 = region_label(ir)
    cell_fill = (/select_region_srex(region_label(ir), lsdata)/)
    cell_fill!0 = "lat"
    cell_fill!1 = "lon"
    cell_fill&lat = lsdata&lat
    cell_fill&lon = lsdata&lon
    cell_fill = where(ismissing(cell_fill), -1., 1.)
    models1 = new((/dim_MOD/), "string")
    colors1 = new((/dim_MOD/), typeof(colors))
    dashes1 = new((/dim_MOD/), typeof(dashes))
    thicks1 = new((/dim_MOD/), typeof(thicks))
    markers1 = new((/dim_MOD/), typeof(markers))
    lin_mar_mod1 = new((/dim_MOD/), "string")
    A0 = read_data(info0[idx_ref])
    time = A0&time
    Atmp1 = A0
    ; ******************mask region*********************
    Atmp1 = (/ mask(A0, cell_fill.eq.-1, False)/)
    ; **************************************************

    Atmp = area_operations(Atmp1, -90., 90., 0., 360., "average", True)
    delete(A0)
    A = Atmp
    delete(Atmp)
    qsort(A)  ; rank the observation in ascending order

    ; -----------------------------------------------------------
    ; netcdf output file

    nc_filename = work_dir + "fig941_" + var0 + "_" + \
                             region_label(ir) + "_obs.nc"
    A@var = var0
    A@diag_script = "fig941"
    nc_outfile = ncdf_write(A, nc_filename)

    D = new((/dim_MOD, dimsizes(time)/), "double")
    C = new((/dim_MOD, dimsizes(time)/), "double")

    jmod = 0
    do imod = 0, dim_MOD - 1
      if (imod .ne. idx_ref) then
        B0 = read_data(info0[imod])
        B0&time = time
        ; ******************mask region**********************
        Btmp1 = B0
        Btmp1 = (/mask(B0, cell_fill.eq.-1, False)/)
        ; **************************************************
        Btmp = area_operations(Btmp1, -90., 90., 0., 360., "average", True)
        delete(B0)
        B = Btmp
        qsort(B)
        delete(Btmp)
        delete(Btmp1)
        ; --------------------------------------------------------
        ; Bias
        C(jmod, :) = B
        D(jmod, :) = B - avg(B) + avg(A)
        delete(B)
        colors1(jmod) = colors(imod)
        dashes1(jmod) = dashes(imod)
        thicks1(jmod) = thicks(imod)+1
        markers1(jmod) = markers(imod)
        models1(jmod) = names(imod)
        lin_mar_mod1(jmod) = "Lines"
        jmod = jmod + 1
      end if
    end do  ; imod

    UNITS = variable_info[0]@units
    if (UNITS.eq."K") then
      UNITS = "~S~o~N~C"
      C = C - 273.15
      D = D - 273.15
      A = A - 273.15
    end if

    ; -----------------------------------------------------------
    ; netcdf output file
    nc_filename = work_dir + "fig941_" + var0 + "_" \
                           + region_label(ir) + "_models.nc"
    D!0 = "models"
    D!1 = "time"
    D&models =  models1
    D&time = A&time
    D@var = var0
    D@diag_script = "fig941"
    D@existing = "overwrite"
    nc_outfile = ncdf_write(D, nc_filename)

    ; -----------------------------------------------------------
    ; plots of time(x) vs. latitude(y)
    min1 = min(C)
    max1 = max(C)
    res = True
    res@gsnDraw = False                          ; don't draw
    res@gsnFrame = False                 ; don't advance frame
    ; res@tiMainString =variable_info@long_name
    res@tiYAxisString = "RAW MODELS"
    res@tiXAxisString = "OBS"
    res@trYMinF = min1
    res@trYMaxF = max1
    res@trXMinF = min1
    res@trXMaxF = max1
    res@xyLineColors = colors1  ; line colors
    res@xyLineThicknesses = thicks1        ; line thicknesses
    res@xyDashPatterns = dashes1       ; line patterns
    res@gsnLeftString = case1
    res@gsnRightString = UNITS
    res@vpHeightF = 0.4
    res@vpWidthF = 0.8
    res@tiMainFontHeightF      = 0.02
    res@txFontHeightF = 0.02
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    ; res@tfPolyDrawOrder = "Draw"
    plot(ir) = gsn_csm_xy(wks_r, A, C, res)
    lnres = True
    dum_diag(ir) = gsn_add_polyline(wks_r, plot(ir), (/min1, max1/), \
                                    (/min1, max1/), lnres)
    res@tiYAxisString = "MODELS"
    res@tiXAxisString = "OBS"
    plot1(ir) = gsn_csm_xy(wks, A, D, res)
    dum_diag1(ir) = gsn_add_polyline(wks, plot1(ir), (/min1, max1/), \
                                     (/min1, max1/), lnres)
    delete(C)
    delete(A)
    delete(D)
  end do

  ; Panel plots
  pan = True
  pan@gsnMaximize = True
  pan@gsnFrame = False
  pan@gsnPaperOrientation = "portrait"
  n_p = dimsizes(plot) / 3
  if (mod(dimsizes(plot), 3) .ne. 0) then
    n_p = n_p + 1
  end if
  pan@txString = "Centered " + LONG_NAME + " bias vs. " + \
                 names(idx_ref)
  gsn_panel(wks, plot1, (/3, n_p/), pan)

  ; ***********add legend****************************
  frame(wks)
  pan@txString = "Uncentered " + LONG_NAME + " bias vs. " + \
                 names(idx_ref)
  gsn_panel(wks_r, plot, (/3, n_p/), pan)
  frame(wks_r)
  add_legenda_page(wks_l, models1, colors1, dashes1, lin_mar_mod1,\
                   thicks1, "0.99")
  frame(wks_l)

  ; -------------------------------------------------------------
  ; Add provenance
  ; -------------------------------------------------------------

  caption = "Modelled versus observed monthly mean temperature"

  log_provenance(nc_filename, \
                 wks@fullname, \
                 caption, \
                 (/"corr"/), \
                 (/"reg"/), \
                 "line", \
                 (/"cionni_irene"/), \
                 (/"flato13ipcc"/), \
                 metadata_att_as_array(info0, "filename"))

end
