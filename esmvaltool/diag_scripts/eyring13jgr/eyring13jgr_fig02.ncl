; #############################################################################
; eyring13jgr_fig02.ncl
; #############################################################################
; Author: Irene Cionni (ENEA, ITALY)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;   Calculates time series of area-weighted total ozone averaged over given
;   domains and seasons.
;
; Required diag_script_info attributes:
;   e13fig02_latrange: extremes of latitude band(s)
;
; Optional diag_script_info attributes:
;   e13fig02_legend: switch on or off the full legend (True, False)
;   e13fig02_legend_MMM: switch on or off the legend of MMM groups
;   e13fig02_season: season(s) - ANN = annual, DJF = December, January,
;     February, etc
;   list_chem_mod: list of models in the group “chem” (array of strings,
;     default = (/"All"/))
;   nochem_mod: list of models in the group “nochem” (array of strings,
;     default = (/""/))
;   styleset: style scheme for line types and colors (e.g. "CMIP5")
;
; Required variable_info attributes:
;   none
;
; Optional variable_info attributes:
;   none
;
; Caveats:
;   none
;
; Modification history:
;   20190515-A_laue_ax: code rewritten for v2, optional masking/regridding/
;                       level selection (for 3-dim variables) now done by
;                       backend
;   20170807-A_bock_ls: added tags for reporting
;   20160929-A_laue_ax: added optional masking of missing values before
;                       calculating averages
;   20160826-A_laue_ax: added ensemble member to model name in legend if
;                       ensemble is other than r1i1p1
;   20151203-A_laue_ax: code clean-up, moved call to write_references to
;                       beginning of code
;   201510??-A-cion_ir: written
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"

load "$diag_scripts/shared/plot/misc_function.ncl"

; ************************************************
; ************************************************

begin
  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD1 = ListCount(info0)
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "e13fig02_season", "ANN")
  set_default_att(diag_script_info, "e13fig02_legend", True)
  set_default_att(diag_script_info, "e13fig02_legend_MMM", True)
  set_default_att(diag_script_info, "styleset", "CMIP5")

  if (.not.isatt(diag_script_info, "e13fig02_latrange")) then
    error_msg("f", DIAG_SCRIPT, "", "no latitude bands " + \
              "(e13fig03_latrage) defined.")
  end if

  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

  flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
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
  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f02_" + var0)
  output_filename = wks@fullname

  plot = new(dimsizes(diag_script_info@e13fig02_season), "graphic")
  plot_num = integertochar(ispan(97, 122, 1))

  ;-----------------------------------------------------------------

  ; get unique model experiments (excluding observations)
  dim_EXP = 0  ; dimsizes of different experiments
  do n = 0, dimsizes(exps) - 1
    if (projects(n) .ne. "OBS" .and. projects(n) .ne. "obs4mips") then
      if (isvar("label_exp")) then  ; n always ge 1 in that case
        if (exps(n) .ne. exps(n - 1)) then
          dim_EXP = dim_EXP + 1
          ltmp = array_append_record(label_exp, exps(n), 0)
          delete(label_exp)
          label_exp = ltmp
          delete(ltmp)
        end if
      else
        label_exp = exps(n)
        dim_EXP = 1
      end if
    end if
  end do

  y1 = metadata_att_as_array(info0, "start_year")
  y2 = metadata_att_as_array(info0, "end_year")

  enssize = dimsizes(ind(exps(0) .eq. exps))
  yearsize = max(toint(y2) - min(toint(y1))) + 1
  dum = new((/dimsizes(diag_script_info@e13fig02_season), (dim_EXP) * 3/), \
    "graphic")

  if ((isatt(diag_script_info, "list_chem_mod"))) then
    list_chem_mod = diag_script_info@list_chem_mod
    attrib_chem = "_CHEM"
    if (list_chem_mod.eq."All") then
      delete(list_chem_mod)
      list_chem_mod = names(index_mod)
      attrib_chem = " "
    end if
  else
    list_chem_mod = names
    attrib_chem = " "
  end if

  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    list_nochem_mod = diag_script_info@list_nochem_mod
  end if
  latrange = diag_script_info@e13fig02_latrange
  season = diag_script_info@e13fig02_season
  if ((isatt(diag_script_info, "e13fig02_XMin"))) then
    XMin = diag_script_info@e13fig02_XMin
  else
    XMin= new(dimsizes(season), "float")
  end if
  if ((isatt(diag_script_info, "e13fig02_XMax"))) then
    XMax = diag_script_info@e13fig02_XMax
  else
    XMax=new(dimsizes(season), "float")
  end if
  if ((isatt(diag_script_info, "e13fig02_YMin"))) then
    YMin = diag_script_info@e13fig02_YMin
  else
    YMin=new(dimsizes(season), "float")
  end if
  if ((isatt(diag_script_info, "e13fig02_YMax"))) then
    YMax = diag_script_info@e13fig02_YMax
  else
    YMax=new(dimsizes(season), "float")
  end if

  color_fill_exp=(/(/0.90, 0.90, 0.90/), \
                   (/0.105, 0.67, 0.8/), \
                   (/0.88, 0.88, 0.88/), \
                   (/0.66, 0.66, 0.66/), \
                   (/0.55, 0.55, 0.55/), \
                   (/0.33, 0.33, 0.33/), \
                   (/0.22, 0.22, 0.22/)/)

  index_mark=(/4, 6, 7, 5, 12/)
  plot_poly = new(dimsizes(diag_script_info@e13fig02_season) * dim_EXP, \
    "graphic")
  if (dim_OBS.gt.0) then
    plot_obs = new(dimsizes(diag_script_info@e13fig02_season) * dim_OBS, \
      "graphic")
  end if
  plot_line = new(dimsizes(diag_script_info@e13fig02_season) * dim_EXP \
                  * (dimsizes(list_nochem_mod) + 1), "graphic")

  ;----------------------------------------------------------------------------

  ip = 0
  iobs1 = 0
  in = 0
  let = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", \
          "n", "o", "p", "q", "o", "s", "t", "u", "v", "w", "x", "y", "z"/)
  iobs2 = 0

  ;; domain for provenance
  DM_tag = new(dimsizes(season) , string)

  ; filename(s) for netCDF output
  ; (one season/region per file as writing all results to one file
  ; does not work)

  nc_outfile = new(dimsizes(season), "string")

  do ipan = 0, dimsizes(season) - 1  ; panel loop
    month=(/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", \
            "SEP", "OCT", "NOV", "DEC"/)
    ind_month = season(ipan)
    do im = 0, 11
      if (season(ipan).eq.month(im)) then
        delete(ind_month)
        ind_month = im + 1
      end if
    end do

    ;; tagging
    DM_tag(ipan) = "reg"
    if (latrange(ipan, 1).eq.90) then
      if (latrange(ipan, 0).eq.-90) then
        DM_tag(ipan) = "glob"
      else if (latrange(ipan, 0).eq.60) then
        DM_tag(ipan) = "polar"
      else if (latrange(ipan, 0).eq.0) then
        DM_tag(ipan) = "nh"
      end if
      end if
      end if
    else if (latrange(ipan, 1).eq.30) then
      if (latrange(ipan, 0).eq.30.) then
        DM_tag(ipan) = "trop"
      end if
    else if (latrange(ipan, 1).eq.0) then
      if (latrange(ipan, 0).eq.-90) then
        DM_tag(ipan) = "sh"
      end if
    else if (latrange(ipan, 1).eq.-60) then
      if (latrange(ipan, 0).eq.-90.) then
        DM_tag(ipan) = "polar"
      end if
    end if
    end if
    end if
    end if

    if (dimsizes(list_chem_mod) .gt. 0) then
      chem_t = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), "float")
      chem_t_error = new((/dim_EXP, yearsize/), "float")
      colors_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
        typeof(colors))
      dashes_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
        typeof(dashes))
      thicks_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
        typeof(thicks))
      lines_chem  = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
        typeof(colors))

      time_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), \
        "float")
      time_c = ispan(min(toint(y1)), max(toint(y2)), 1)
      time_chem!0 = "exp"
      time_chem!2 = "time"
      time_chem&time = time_c

      chem_t!0 = "exp"
      chem_t!1 = "model"
      chem_t!2 = "time"
      chem_t&exp = label_exp
      chem_t&model = array_append_record(list_chem_mod, "MMM_chem", 0)
      chem_t&time = time_c

      chem_t@var = var0 + "_chem_" + season(ipan) + "_" + DM_tag(ipan)
      chem_t@diag_script = (/DIAG_SCRIPT/)
    end if

    if (dimsizes(list_nochem_mod) .gt. 0) then
      nochem_t = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), \
        "float")
      nochem_t_error = new((/dim_EXP, yearsize/), "float")
      colors_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
        typeof(colors))
      dashes_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
        typeof(dashes))
      thicks_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
        typeof(thicks))
      lines_nochem  = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
        typeof(colors))

      time_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), \
        "float")
      time_nc = ispan(min(toint(y1)), max(toint(y2)), 1)
      time_nochem!0 = "exp"
      time_nochem!2 = "time"
      time_nochem&time = time_nc

      nochem_t!0 = "exp"
      nochem_t!1 = "model"
      nochem_t!2 = "time"
      nochem_t&exp = label_exp
      nochem_t&model = array_append_record(list_nochem_mod, "MMM_nochem", 0)
      nochem_t&time = time_nc

      nochem_t@var = var0 + "_nochem_" + season(ipan) + "_" + DM_tag(ipan)
      nochem_t@diag_script = (/DIAG_SCRIPT/)
    end if

    if (dim_OBS.gt.0) then
      obs_t       = new((/dim_OBS, yearsize/), "float")
      colors_obs  = new((/dim_OBS/), typeof(colors))
      markers_obs = new((/dim_OBS/), typeof(dashes))
      thicks_obs  = new((/dim_OBS/), typeof(thicks))
      lines_obs   = new((/dim_OBS/), typeof(colors))

      time_obs = new((/dim_OBS, yearsize/), "float")
      time_o = ispan(min(toint(y1)), max(toint(y2)), 1)
      time_obs!1 = "time"
      time_obs&time = time_o
      obs_t!0 = "obs"
      obs_t!1 = "time"
      obs_t&obs = list_obs_mod
      obs_t&time = time_o

      obs_t@var = var0 + "_obs_" + season(ipan) + "_" + DM_tag(ipan)
      obs_t@diag_script = (/DIAG_SCRIPT/)
    end if

    ; ----------------------------chem------------------------------------

    nochem = 0
    chem = 0
    iobs = 0

    do imod = 0, dim_MOD1 - 1
      ; Extract data
      A0 = read_data(info0[imod])
      dummy1 = area_operations(A0(time|:, lat|:, lon|:), latrange(ipan, 0), \
                               latrange(ipan, 1), 0., 360., "average", True)
      if (var0.eq."tro3") then
         dummy1 = dummy1 * 10^-3
      end if

      delete(A0)

      if (season(ipan).ne."ANN") then
        tmp = time_operations(dummy1, toint(y1(imod)), toint(y2(imod)), \
                              "extract", tostring(ind_month), True)
        dummy2 = time_operations(tmp, toint(y1(imod)), toint(y2(imod)), \
                                 "average", "yearly", True)
        delete(tmp)
      else
        dummy2 = time_operations(dummy1, toint(y1(imod)), toint(y2(imod)), \
                                 "average", "yearly", True)
      end if

      iexp = ind(label_exp.eq.exps(imod))

      if (any(names(imod).eq.list_chem_mod)) then
        chem_t(iexp, chem, {toint(y1(imod)):toint(y2(imod))}) = tofloat(dummy2)
        time_chem(iexp, chem, {toint(y1(imod)):toint(y2(imod))}) = \
          ispan(toint(y1(imod)), toint(y2(imod)), 1)
        colors_chem(iexp, chem) = colors(imod)
        dashes_chem(iexp, chem) = dashes(imod)
        thicks_chem(iexp, chem) = 2.
        lines_chem(iexp, chem) = "Lines"
        chem = chem + 1
      end if

      if (any(names(imod).eq.list_nochem_mod)) then
        nochem_t(iexp, nochem, {toint(y1(imod)):toint(y2(imod))}) = \
          tofloat(dummy2)
        time_nochem(iexp, nochem, {toint(y1(imod)):toint(y2(imod))}) = \
          ispan(toint(y1(imod)), toint(y2(imod)), 1)
        colors_nochem(iexp, nochem) = colors(imod)
        dashes_nochem(iexp, nochem) = dashes(imod)
        thicks_nochem(iexp, nochem) = 1.
        lines_nochem(iexp, nochem) = "Lines"
        nochem = nochem + 1
      end if

      if (any(names(imod).eq.list_obs_mod)) then
        obs_t(iobs, {toint(y1(imod)):toint(y2(imod))}) = tofloat(dummy2)
        time_obs(iobs, {toint(y1(imod)):toint(y2(imod))}) = \
          ispan(toint(y1(imod)), toint(y2(imod)), 1)
        colors_obs(iobs) = colors(imod)
        markers_obs(iobs) = markers(imod)
        thicks_obs(iobs) = 3.
        lines_obs(iobs) = "Markers"
        iobs = iobs + 1
      end if

      delete(dummy1)
      delete(dummy2)
    end do  ; (do imod)

    do iexp = 0, dim_EXP - 1
      if (chem.gt.0) then
        nn = dimsizes(list_chem_mod) - 1
        chem_t(iexp, nn + 1, :) = \
          tofloat(dim_avg_n(chem_t(iexp, 0:nn, :), 0))
        chem_t_error(iexp, :) = \
          tofloat(dim_stddev_n(chem_t(iexp, 0:nn, :), 0))
        time_chem(iexp, nn + 1, \
          {toint(min(time_chem(iexp, 0:nn, :))):toint(max( \
          time_chem(iexp, 0:nn, :)))}) = \
          ispan(toint(min(time_chem(iexp, 0:nn, :))), \
          toint(max(time_chem(iexp, 0:nn, :))) ,1)

        ; styles for MMM (chem)
        dashes_chem(iexp, nn + 1) = 0
        thicks_chem(iexp, nn + 1) = 3.
        lines_chem(iexp, nn + 1) = "Lines"
        colors_chem(iexp, nn + 1) = "black"  ; default

        if (((label_exp).eq."historical").or.((label_exp).eq."refC1")) then
          colors_chem(iexp, nn + 1) = "black"
        end if
        if ((label_exp).eq."rcp26") then
          colors_chem(iexp, nn + 1) = "(/0., 0., 1./)"
        end if
        if ((label_exp).eq."rcp45") then
          colors_chem(iexp, nn + 1) = "(/0, 4757, 0.7373, 1./)"
        end if
        if ((label_exp).eq."rcp60") then
          colors_chem(iexp, nn + 1) = "(/1., 0.5098, 0.1765/)"
        end if
        if ((label_exp).eq."rcp85") then
          colors_chem(iexp, nn + 1) = "(/1., 0., 0./)"
        end if
      end if

      ; ****** no chem **********************************

      if (nochem.gt.0) then
        nn = dimsizes(list_nochem_mod) - 1
        nochem_t(iexp, nn + 1, :) = \
          tofloat(dim_avg_n(nochem_t(iexp, 0:nn, :), 0))
        nochem_t_error(iexp, :) = \
          tofloat(dim_stddev_n(nochem_t(iexp, 0:nn, :), 0))
        time_nochem(iexp, nn + 1, \
          {toint(min(time_nochem(iexp, 0:nn, :))):toint(max( \
          time_nochem(iexp, 0:nn, :)))}) = \
          ispan(toint(min(time_nochem(iexp, 0:nn, :))), \
          toint(max(time_nochem(iexp, 0:nn, :))), 1)

        ; styles for MMM (nochem)
        dashes_nochem(iexp, nn + 1) = 1
        thicks_nochem(iexp, nn + 1) = 3
        lines_nochem(iexp, nn + 1) = "Lines"
        colors_nochem(iexp, nn + 1) = "blue"

        if (((label_exp).eq."historical").or.((label_exp).eq."refC1")) then
          colors_nochem(iexp, nn + 1) = "black"
        end if
        if ((label_exp).eq."rcp26") then
          colors_nochem(iexp, nn + 1) = "(/0., 0., 1./)"
        end if
        if ((label_exp).eq."rcp45") then
          colors_nochem(iexp, nn + 1) = "(/0, 4757, 0.7373, 1./)"
        end if
        if ((label_exp).eq."rcp60") then
          colors_nochem(iexp, nn + 1) = "(/1., 0.5098, 0.1765/)"
        end if
        if ((label_exp).eq."rcp85") then
          colors_nochem(iexp, nn + 1) = "(/1., 0., 0./)"
        end if
      end if
    end do

    ; --------PLOT---------------------------------------

    ti = ispan(1, (dim_EXP) * 3, 1)

    if (latrange(ipan, 0).ge.0.) then
      latstring1 = latrange(ipan,0) + "N"
    else
      latstring1 = -latrange(ipan,0) + "S"
    end if

    if (latrange(ipan, 1).ge.0.) then
      latstring2 = latrange(ipan, 1) + "N"
    else
      latstring2 = -latrange(ipan,1) + "S"
    end if

    res = True  ; plot mods desired
    res@tiYAxisString = variable_info[0]@units  ; add units title
    res@tiXAxisString = "Year"               ; y-axis label
    res@tiMainString  = let(ipan) + ") " + variable_info[0]@long_name \
      + " " + season(ipan) + " " + latstring1 + "-" + latstring2
    res@gsnFrame      = False                ; don't draw yet
    res@gsnDraw       = False

    if (.not.ismissing(YMax(ipan))) then
      res@trYMaxF = YMax(ipan)
    else
      vv = (/max(chem_t), max(nochem_t)/)
      res@trYMaxF = max((/max(chem_t), max(nochem_t)/))
    end if

    if (.not.ismissing(YMin(ipan))) then
      res@trYMinF = YMin(ipan)
    else
      res@trYMinF = min((/min(chem_t), min(nochem_t)/))
    end if

    if (.not.ismissing(XMax(ipan))) then
      res@trXMaxF = XMax(ipan)
    end if

    if (.not.ismissing(XMin(ipan))) then
      res@trXMinF = XMin(ipan)
    end if

    res@vpHeightF         = 0.4                  ; change aspect ratio of plot
    res@vpWidthF          = 0.5
    res@xyLineThicknesses = ndtooned(thicks_chem)
    res@xyLineColors      = ndtooned(colors_chem)  ; line color
    res@xyDashPatterns    = ndtooned(dashes_chem)
    res@xyMarkLineModes   = ndtooned(lines_chem)

    plot(ipan) = gsn_csm_xy(wks, onedtond(ndtooned(time_chem(:, :, :)),  \
      (/dim_EXP * (dimsizes(list_chem_mod) + 1), yearsize/)), \
      onedtond(ndtooned(chem_t),                              \
      (/dim_EXP * (dimsizes(list_chem_mod) + 1), yearsize/)), res)

    ; *************** add error *****************

    dim_dummy = dimsizes(time_chem(0, dimsizes(list_chem_mod), :))
    xp = new((/dim_EXP, 2 * dim_dummy/), typeof(time_chem))
    yp = new((/dim_EXP, 2 * dim_dummy/), typeof(chem_t))

    do k = 0, dim_dummy - 1
      xp(:, k) = time_chem(:, dimsizes(list_chem_mod), k)
      xp(:, dim_dummy + k) = \
        time_chem(:, dimsizes(list_chem_mod), dim_dummy - 1 - k)
      yp(:,k) = chem_t(:, dimsizes(list_chem_mod), k) - chem_t_error(:, k)
      yp(:, dim_dummy + k) = \
        chem_t(:, dimsizes(list_chem_mod), dim_dummy-1-k) \
        + chem_t_error(:, dim_dummy - 1 - k)
    end do

    do iexp = 0, dim_EXP - 1
      gsres = True
      gsres@gsFillColor = color_fill_exp(iexp, :)
      gsres@tfPolyDrawOrder = "Draw"
      plot_poly(ip) = gsn_add_polygon(wks, plot(ipan), xp(iexp,:), \
        yp(iexp,:), gsres)
      ip = ip + 1
    end do

    ; ************ add nochem **************************

    if (nochem.gt.0) then
      do iexp = 0, dim_EXP - 1
        do imod = 0, dimsizes(list_nochem_mod)
          plres = True
          plres@gsLineColor = colors_nochem(iexp, imod)
          plres@gsLineDashPattern = dashes_nochem(iexp, imod)
          plres@gsLineThicknessF = thicks_nochem(iexp, imod)
          plot_line(in) = gsn_add_polyline(wks, plot(ipan), \
            time_nochem(iexp, imod, :), nochem_t(iexp, imod, :), plres)
          in = in + 1
        end do
      end do
    end if

    ; *************** add obs *******************

    do iobs1 = 0, dim_OBS - 1
      pmres = True
      pmres@gsMarkerSizeF = 0.01
      pmres@gsMarkerThicknessF = thicks_obs(iobs1)
      pmres@gsMarkerColor = colors_obs(iobs1)
      pmres@gsMarkerIndex = markers_obs(iobs1)
      plot_obs(iobs2) = gsn_add_polymarker(wks,plot(ipan), \
        time_obs(iobs1, :), obs_t(iobs1, :), pmres)
      iobs2 = iobs2 + 1
    end do

    ; ************** netCDF output **************

    ; put all results in one array for output

    ; make space for MMM (if chem_t / nochem_t are present)
    if (chem.gt.0) then
      chem = chem + 1
    end if
    if (nochem.gt.0) then
      nochem = nochem + 1
    end if

    outfield = new((/dim_EXP, chem + nochem + dim_OBS, yearsize/), "float")
    outfield@var = var0
    outfield@diag_script = (/DIAG_SCRIPT/)
    outfield!0 = "exp"
    outfield!1 = "dataset"
    outfield!2 = "time"

    outfield&exp = ""
    n = 0

    datasets = new((/chem + nochem + dim_OBS/), "string")

    if (chem.gt.0) then
      outfield(:, n:n + chem - 1, :) = (/chem_t/)
      outfield&exp = (/chem_t&exp/)
      datasets(n:n + chem - 1) = (/chem_t&model/)
      outfield&time = tofloat(chem_t&time)
      n = n + chem
    end if
    if (nochem.gt.0) then
      outfield(:, n:n + nochem - 1, :) = (/nochem_t/)
      outfield&exp = (/nochem_t&exp/)
      datasets(n:n + nochem - 1) = (/nochem_t&model/)
      outfield&time = tofloat(nochem_t&time)
      n = n + nochem
    end if
    if (dim_OBS.gt.0) then
      outfield(0, n:n + dim_OBS - 1, :) = (/obs_t/)
      datasets(n:n + dim_OBS - 1) = (/obs_t&obs/)
      outfield&time = tofloat(obs_t&time)
    end if

    ; filename for netCDF output (panel plot)
    nc_filename = work_dir + "eyring13f01_panel_" + var0 + "_" \
      + season(ipan) + "_" + DM_tag(ipan) + ".nc"

    outfield&dataset = datasets
    nc_outfile(ipan) = ncdf_write(outfield, nc_filename)

    delete(outfield)
    delete(datasets)

    ; ******************************************

    delete(nochem_t)
    delete(chem_t)
    delete(xp)
    delete(yp)
    delete(ind_month)
  end do  ; loop over all "seasons"

  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2

  gsn_panel(wks, plot, (/tointeger(dimsizes(season) / 2), 2/), Pres)

  if (dim_OBS.gt.0) then
    add_legenda(wks, list_obs_mod, colors_obs, markers_obs, lines_obs, \
      thicks_obs, "0.35")
  end if

  if (diag_script_info@e13fig02_legend_MMM) then
    chem_att = new(dim_EXP, "string")
    nochem_att = new(dim_EXP, "string")
    do iexp = 0, dim_EXP - 1
      chem_att(iexp) = attrib_chem
      nochem_att(iexp) = "_NOCHEM"
    end do

    if (nochem.gt.0) then
      add_legenda(wks, label_exp + nochem_att, \
        colors_nochem(:, dimsizes(list_nochem_mod)), \
        dashes_nochem(:, dimsizes(list_nochem_mod)), \
        lines_nochem(:, dimsizes(list_nochem_mod)),  \
        thicks_nochem(:, dimsizes(list_nochem_mod)), "0.08")
    end if
  end if

  ; ****************** a second page ********************

  if (diag_script_info@e13fig02_legend) then
    ; add ensemble member to model name in legend

    modelnames = list_chem_mod

    do i = 0, dim_MOD1 - 1
      if (isvar("idx")) then
        delete(idx)
      end if
      idx = ind(modelnames.eq.names(i))
      if (all(ismissing(idx))) then
        continue
      end if
      do j = 0, dimsizes(idx) - 1
        modelnames(idx(j)) = modelnames(idx(j)) + " ("  \
          + ensembles(idx(j)) + ")"
      end do
    end do

    add_legenda(wks, modelnames, \
      colors_chem(0, :dimsizes(list_chem_mod) - 1), \
      dashes_chem(0, :dimsizes(list_chem_mod) - 1), \
      lines_chem(0, :dimsizes(list_chem_mod) - 1),  \
      thicks_chem(0, :dimsizes(list_chem_mod) - 1), "0.3")

    if (nochem.gt.0) then
      add_legenda(wks, list_nochem_mod, \
        colors_nochem(0, :dimsizes(list_nochem_mod) - 1), \
        dashes_nochem(0, :dimsizes(list_nochem_mod) - 1), \
        lines_nochem(0, :dimsizes(list_nochem_mod) - 1),  \
        thicks_nochem(0, :dimsizes(list_nochem_mod) - 1), "0.05")
    end if
  end if

  if (config_user_info@write_plots.eq."True") then
    frame(wks)
    log_info(" Wrote " + output_filename)
  end if

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "mean"/)
  domain = DM_tag
  plottype = "times"
  contrib_authors = (/"cion_ir", "laue_ax"/)
  diag_ref = "eyring13jgr"

  caption = "Time series of area-weighted "  + var0 + \
            " over given domains and seasons" + \
            ", similar to IPCC AR5, fig. 9.10 and Eyring et al., " + \
            "JGR, 2013, fig. 2."

  climofiles = metadata_att_as_array(info0, "filename")

  do i = 0, dimsizes(season) - 1
    log_provenance(nc_outfile(i), output_filename, caption, statistics, \
                   domain, plottype, contrib_authors, diag_ref, climofiles)
  end do

  ; ------------------------------------------------------------------------

  leave_msg(DIAG_SCRIPT, "")

end
