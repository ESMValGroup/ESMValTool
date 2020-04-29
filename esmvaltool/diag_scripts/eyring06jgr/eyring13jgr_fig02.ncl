; ############################################################################
; eyring13jgr_fig02.ncl                                                      #
; ############################################################################
; Author: Irene Cionni (ENEA, Italy)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; ############################################################################
;
;  Description:
;     Calculates time series of area-weighted total ozone from 1960-2005
;     averaged over given domains and seasons.
;
;  Required diag_script_info attributes (diagnostics specific):
;     * min: extremes of latitude band
;     * max: extremes of latitude band
;     * season: ANN=Annual, DJF=December, January, February, etc
;     * lev: selection of level if original fiels is T3*, is
;       not required if variables already in T2*s fields
;
;  Optional diag_script_info attributes (diagnostics specific):
;    * list_chem_mod: names of models to apply a MMM, "All" if not need
;    * list_nochem_mod: names of models to apply a MMM, "None" if not need
;    * legend: switch on or off the full legend.Opt "True","False"
;    * legend_MMM: switch on or off the legend of MMM groups
;    * apply_masking: mask missing values from all data sets before
;                     calculating averages(True/False) (T2Ms fields only)
;    * grid: if masking is selected, all data must be regridded to a common
;            grid (coarsest, finest, ref, 1x1; default = 1x1)
;    * refmodel: if masking is selected and grid = "ref", a reference model
;            has to be specified
; ############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")
  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD1 = ListCount(info0)
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")
  print(ensembles)
  y1 = (/ metadata_att_as_array(info0, "start_year")/)
  y2 = (/ metadata_att_as_array(info0, "end_year") /)

  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
    do imod = 0, dimsizes(names) - 1
      if (names(imod) .eq. refname) then
        iref = imod
      end if
    end do
  end if
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

  set_default_att(diag_script_info, "season", "ANN")
  set_default_att(diag_script_info, "styleset", "CMIP5")
  set_default_att(diag_script_info, "latmin", -90.)
  set_default_att(diag_script_info, "latmax", 90.)
  set_default_att(diag_script_info, "multimean", False)
  set_default_att(diag_script_info, "show_trend", "False")
  set_default_att(diag_script_info, "xMin", 1960)
  set_default_att(diag_script_info, "xMax", 2001)
  set_default_att(diag_script_info, "yMin", -20)
  set_default_att(diag_script_info, "yMax", 20)
  set_default_att(diag_script_info, "legend", True)
  set_default_att(diag_script_info, "legend_MMM", True)
  set_default_att(diag_script_info, "list_chem_mod", "All")
  set_default_att(diag_script_info, "list_chem_mod_string", "CMIP5")
  set_default_att(diag_script_info, "apply_masking", False)

  seasons = diag_script_info@season
  latrange = new((/dimsizes(seasons), 2/), "float")
  latrange(:, 0) = diag_script_info@latmin
  latrange(:, 1) = diag_script_info@latmax
  seasons = diag_script_info@season
  trend = diag_script_info@show_trend
  XMin = diag_script_info@xMin
  XMax = diag_script_info@xMax
  YMin = diag_script_info@yMin
  YMax = diag_script_info@yMax

  colors = project_style(info0, diag_script_info, "colors")
  dashes = project_style(info0, diag_script_info, "dashes")
  thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
  markers = project_style(info0, diag_script_info, "markers")

  color_fill_exp = (/(/0.90, 0.90, 0.90/), \
                    (/0.105, 0.67, 0.8/), \
                    (/0.88, 0.88, 0.88/), \
                    (/0.66, 0.66, 0.66/), \
                    (/0.55, 0.55, 0.55/), \
                    (/0.33, 0.33, 0.33/), \
                    (/0.22, 0.22, 0.22/)/)

  index_mark = (/4, 6, 7, 5, 12/)

  flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  multimean = diag_script_info@multimean

  dim_EXP = 1  ; dimsizes of different experiments
  label_exp = exps(0)
  do n = 1, dimsizes(exps) - 1
    if (projects(n).ne."OBS") then
      if (exps(n).ne.exps(n - 1)) then
        dim_EXP = dim_EXP + 1
        ltmp = array_append_record(label_exp, exps(n), 0)
        delete(label_exp)
        label_exp = ltmp
        delete(ltmp)
      end if
    end if
  end do

  enssize = dimsizes(ind(exps(0).eq.exps))
  yearsize = max(toint(y2)) - min(toint(y1)) + 1
  dum = new((/dimsizes(seasons), (dim_EXP) * 3/), graphic)
  list_chem_mod = diag_script_info@list_chem_mod
  attrib_chem = "_CHEM"
  print(list_chem_mod)
  if (any(list_chem_mod.eq."All")) then
    delete(list_chem_mod)
    list_chem_mod = names
    attrib_chem = " "
  end if
  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    list_nochem_mod = diag_script_info@list_nochem_mod
  else
    list_nochem_mod = "0"
  end if

  if (multimean)
    dimTot = dim_MOD + dim_OBS + 1
  else
    dimTot = dim_MOD + dim_OBS
  end if
  ind_chem_mod = new(dimsizes(list_chem_mod), "integer")
  do ic = 0, dimsizes(list_chem_mod) - 1
    a = ind(names.eq.list_chem_mod(ic))
    if (flag_mod(a).gt.0) then
      ind_chem_mod(ic) = a
      delete(a)
    end if
  end do
  if ((isatt(diag_script_info, "list_nochem_mod"))) then
    list_nochem_mod = diag_script_info@list_nochem_mod
    ind_nochem_mod = ind((names.eq.list_nochem_mod).and.(flag_mod.gt.0))
  end if
  apply_masking = diag_script_info@apply_masking

  if (isatt(variable_info[0], "reference_dataset")) then
    OBS_ref = variable_info[0]@reference_dataset
    do iobs = 0, dimsizes(OBS_ref) - 1
      ind_OBS_ref = ind(names.eq.OBS_ref(iobs))
      if (.not.ismissing(ind_OBS_ref)) then
        break
      end if
    end do
    delete(OBS_ref)
    OBS_ref = names(ind_OBS_ref)
  end if

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end
begin
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f02_" + var0)
  wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f02_legend_" + var0)
  plot = new(dimsizes(seasons), "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  plot_poly = new(dimsizes(seasons) * dim_EXP, "graphic")
  if (dim_OBS.gt.0) then
    plot_obs = new(dimsizes(seasons) * dim_OBS, "graphic")
  end if
  plot_line = new(dimsizes(seasons) * dim_EXP \
                  * (dimsizes(list_nochem_mod) + 1), "graphic")

  ip = 0
  iobs1 = 0
  in = 0
  let = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",\
          "n", "o", "p", "q", "o", "s", "t", "u", "v", "w", "x", "y", "z"/)
  iobs2 = 0

  do  ipan = 0, dimsizes(seasons) - 1  ; panel loop
    month = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", \
              "SEP", "OCT", "NOV", "DEC"/)
    ind_month = seasons(ipan)
    do im = 0, 11
      if (seasons(ipan).eq.month(im)) then
        delete(ind_month)
        ind_month = im + 1
      end if
    end do

    chem_t   = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), "float")
    nochem_t = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/),\
                   "float")
    chem_t_error   = new((/dim_EXP, yearsize/), "float")
    nochem_t_error = new((/dim_EXP, yearsize/), "float")

    if (dim_OBS.gt.0) then
      obs_t = new((/dim_OBS, yearsize/), "float")
      colors_obs = new((/dim_OBS/), typeof(colors))
      markers_obs = new((/dim_OBS/), typeof(dashes))
      thicks_obs = new((/dim_OBS/), typeof(thicks))
      lines_obs = new((/dim_OBS/), typeof(colors))
    end if

    colors_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
                        typeof(colors))
    colors_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
                        typeof(colors))

    dashes_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
                        typeof(dashes))
    dashes_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
                        typeof(dashes))
    thicks_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
                        typeof(thicks))
    thicks_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
                        typeof(thicks))
    lines_chem    = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
                        typeof(colors))
    lines_nochem  = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
                        typeof(colors))

    iexp = 0
    nochem = 0
    chem = 0
    iobs = 0

    time_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), float)
    time_c = ispan(min(y1), max(y2), 1)
    time_chem!0 = "exp"
    time_chem!2 = "time"
    time_chem&time = time_c

    time_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), \
                      float)

    time_c = ispan(min(y1), max(y2), 1)
    time_nochem!0 = "exp"
    time_nochem!2 = "time"
    time_nochem&time = time_c
    if (dim_OBS.gt.0) then
      time_obs = new((/dim_OBS, yearsize/), float)
      time_o = ispan(min(y1), max(y2), 1)
      time_obs!1 = "time"
      time_obs&time = time_o
      obs_t!1 = "time"
      obs_t&time = time_o
    end if

    chem_t!2 = "time"
    nochem_t!2 = "time"
    chem_t&time = time_c
    nochem_t&time = time_c
    years1 = y1
    years2 = y2
    plotMod = names
    do imod = 0, dim_MOD1 - 1
    ; Extract data
      rad = 4.0 * atan(1.0) / 180.0
      A0 = read_data(info0[imod])
      if(all(isdim(A0, (/"time", "lat", "lon"/)))) then
        dummy1 = area_operations(A0(time|:, lat|:, lon|:), latrange(ipan, 0), \
                                 latrange(ipan, 1), 0., 360., "average", True)
        coslat = cos(A0&lat * rad)
        coslat!0 = "lat"
        coslat&lat = A0&lat
        delete(coslat)
        delete(A0)
      else
        if(all(isdim(A0, (/"time", "lat"/)))) then
          coslat = cos(A0&lat * rad)
          coslat!0 = "lat"
          coslat&lat = A0&lat
          dummy1 = dim_avg_wgt_Wrap(A0(:, {latrange(ipan, 0): \
                                           latrange(ipan, 1)}), \
                                    coslat({latrange(ipan, 0): \
                                            latrange(ipan, 1)}), 1)
          delete(coslat)
          delete(A0)
        end if
      end if
      if (seasons(ipan).ne."ANN") then
        dummy2 = time_operations(dummy1, y1(imod),  \
                                 y2(imod), "extract", \
                                 tostring(ind_month), True)
        delete(dummy1)
      else
        dummy2 = month_to_annual(dummy1, 1)
        delete(dummy1)
      end if
      if (imod.ge.1) then
        if (projects(imod).ne."OBS") then
          if (exps(imod).ne.exps(imod - 1)) then
            iexp = iexp + 1
            nochem = 0
            chem = 0
          end if
        end if
      end if

      if (any(names(imod).eq.list_chem_mod)) then
        chem_t(iexp, chem, {y1(imod):y2(imod)}) = \
          tofloat(dummy2)
        time_chem(iexp, chem, {y1(imod):y2(imod)}) = \
          ispan(y1(imod), y2(imod), 1)
        colors_chem(iexp, chem) = colors(imod)
        dashes_chem(iexp, chem) = dashes(imod)
        thicks_chem(iexp, chem) = 2.
        lines_chem(iexp, chem) = "Lines"
        chem = chem + 1
      end if

      if (any(names(imod).eq.list_nochem_mod)) then
        nochem_t(iexp, nochem, {y1(imod):y2(imod)}) = \
          tofloat(dummy2)
        time_nochem(iexp, nochem, {y1(imod):y2(imod)}) = \
          ispan(y1(imod), y2(imod), 1)
        colors_nochem(iexp, nochem) = colors(imod)
        dashes_nochem(iexp, nochem) = dashes(imod)
        thicks_nochem(iexp, nochem) = 1.
        lines_nochem(iexp, nochem) = "Lines"
        nochem = nochem + 1
      end if

      if (dim_OBS.gt.0) then
        if (any(names(imod).eq.list_obs_mod)) then
          obs_t(iobs, {y1(imod):y2(imod)}) = \
            tofloat(dummy2)
          time_obs(iobs, {y1(imod):y2(imod)}) = \
            ispan(y1(imod), y2(imod), 1)
          colors_obs(iobs) = colors(imod)
          markers_obs(iobs) = markers(imod)
          thicks_obs(iobs) = 3.
          lines_obs(iobs) = "Markers"
          iobs = iobs + 1
        end if
      end if

      delete(dummy2)
    end do  ; (do imod)
    dim_chem_mod = dimsizes(list_chem_mod)
    do iexp = 0, dim_EXP - 1
      chem_t(iexp, dim_chem_mod, :) = \
        tofloat(dim_avg_n(chem_t(iexp, :dim_chem_mod - 1, :), 0))
      chem_t_error(iexp, :) = \
        tofloat(dim_stddev_n(chem_t(iexp, :dim_chem_mod - 1, :), 0))
      if (chem.gt.0) then
        time_chem(iexp, dim_chem_mod, {toint(min(time_chem(iexp, \
                                                           :dim_chem_mod - 1,\
                                                           :))): \
                                       toint(max(time_chem(iexp, \
                                                           :dim_chem_mod - 1,\
                                                           :)))}) \
         = ispan(toint(min(time_chem(iexp, :dim_chem_mod - 1, :))),\
                 toint(max(time_chem(iexp, :dim_chem_mod - 1, :))), 1)
      end if

      dashes_chem(iexp, dim_chem_mod) = 0
      thicks_chem(iexp, dim_chem_mod) = 3.
      lines_chem(iexp, dim_chem_mod) = "Lines"

      if (((label_exp).eq."historical").or.((label_exp).eq."refC1")) then
        colors_chem(iexp, dim_chem_mod) = "black"
      end if
      if ((label_exp).eq."rcp26") then
        colors_chem(iexp, dim_chem_mod) = "(/0., 0., 1./)"
      end if
      if ((label_exp).eq."rcp45") then
        colors_chem(iexp, dim_chem_mod) = "(/0, 4757, 0.7373, 1./)"
      end if
      if ((label_exp).eq."rcp60") then
        colors_chem(iexp, dim_chem_mod) = "(/1., 0.5098, 0.1765/)"
      end if
      if ((label_exp).eq."rcp85") then
        colors_chem(iexp, dim_chem_mod) = "(/1., 0., 0./)"
      end if

      dim_nochem_mod = dimsizes(list_nochem_mod)
      nochem_t(iexp, dim_nochem_mod, :) = \
        tofloat(dim_avg_n(nochem_t(iexp, :dim_nochem_mod - 1, :), 0))
      nochem_t_error(iexp, :) = \
        tofloat(dim_stddev_n(nochem_t(iexp, :dim_nochem_mod - 1, :), 0))
      if (nochem.gt.0) then
        time_nochem(iexp, dim_nochem_mod, \
                    {toint(min(time_nochem(iexp, :dim_nochem_mod - 1, :))):\
                     toint(max(time_nochem(iexp, :dim_nochem_mod - 1, :)))})\
          = ispan(toint(min(time_nochem(iexp, :dim_nochem_mod - 1, :))), \
                  toint(max(time_nochem(iexp, :dim_nochem_mod - 1, :))), 1)
      end if

      dashes_nochem(iexp, dim_nochem_mod) = 1
      thicks_nochem(iexp, dim_nochem_mod) = 3
      lines_nochem(iexp, dim_nochem_mod) = "Lines"

      if (((label_exp).eq."historical").or.((label_exp).eq."refC1")) then
        colors_nochem(iexp, dimsizes(list_nochem_mod)) = "black"
      end if
      if ((label_exp).eq."rcp26") then
        colors_nochem(iexp, dim_nochem_mod) = "(/0., 0., 1./)"
      end if
      if ((label_exp).eq."rcp45") then
        colors_nochem(iexp, dim_nochem_mod) = "(/0, 4757, 0.7373, 1./)"
      end if
      if ((label_exp).eq."rcp60") then
        colors_nochem(iexp, dim_nochem_mod) = "(/1., 0.5098, 0.1765/)"
      end if
      if ((label_exp).eq."rcp85") then
        colors_nochem(iexp, dim_nochem_mod) =  "(/1., 0., 0./)"
      end if
    end do

    ti = ispan(1, (dim_EXP) * 3, 1)

    if (latrange(ipan, 0).ge.0.) then
      latstring1 = latrange(ipan, 0) + "N"
    else
      latstring1 = -latrange(ipan, 0) + "S"
    end if

    if (latrange(ipan, 1).ge.0.) then
      latstring2 = latrange(ipan, 1) + "N"
    else
      latstring2 = -latrange(ipan, 1) + "S"
    end if

    res = True  ; plot mods desired
    res@tiYAxisString = units  ; add units title
    res@tiXAxisString = "Year"               ; y-axis label
    res@tiMainString = let(ipan) + ") " + var0 + " " + seasons(ipan) + \
      " " + latstring1 + "-" + latstring2
    res@gsnFrame = False                ; don't draw yet
    res@gsnDraw = False

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
    res@vpHeightF = 0.4                    ; change aspect ratio of plot
    res@vpWidthF = 0.5
    res@xyLineThicknesses = ndtooned(thicks_chem)
    res@xyLineColors = ndtooned(colors_chem)  ; line color
    res@xyDashPatterns = ndtooned(dashes_chem)
    res@xyMarkLineModes = ndtooned(lines_chem)

    plot(ipan) = gsn_csm_xy(wks, onedtond(ndtooned(time_chem(:, :, :)), \
                            (/dim_EXP * (dimsizes(list_chem_mod) + 1), \
                             yearsize/)), onedtond(ndtooned(chem_t), \
                            (/dim_EXP * (dimsizes(list_chem_mod) + 1), \
                              yearsize/)), res)

    dim_dummy = dimsizes(time_chem(0, dim_chem_mod, :))
    xp = new((/dim_EXP, 2 * dim_dummy/), typeof(time_chem))
    yp = new((/dim_EXP, 2 * dim_dummy/), typeof(chem_t))

    do k = 0, dim_dummy - 1
      xp(:, k) = time_chem(:, dim_chem_mod, k)
      xp(:, dim_dummy + k) = time_chem(:, dim_chem_mod, dim_dummy - 1 - k)
      yp(:, k) = chem_t(:, dim_chem_mod, k) - chem_t_error(:, k)
      yp(:, dim_dummy + k) = chem_t(:, dim_chem_mod, dim_dummy-1-k) + \
        chem_t_error(:, dim_dummy - 1 - k)
    end do

    do iexp = 0, dim_EXP - 1
      if (label_exp(iexp).ne."OBS") then
        gsres = True
        gsres@gsFillColor = color_fill_exp(iexp, :)
        gsres@tfPolyDrawOrder = "Draw"
        plot_poly(ip) = gsn_add_polygon(wks, plot(ipan), xp(iexp, :), \
                                        yp(iexp, :), gsres)
        ip = ip + 1
      end if
    end do

    if (nochem.gt.0) then
      do iexp = 0, dim_EXP - 1
        do imod = 0, dim_nochem_mod
          plres = True
          plres@gsLineColor = colors_nochem(iexp, imod)
          plres@gsLineDashPattern = dashes_nochem(iexp, imod)
          plres@gsLineThicknessF = thicks_nochem(iexp, imod)
          plot_line(in) = gsn_add_polyline(wks, plot(ipan), \
                                           time_nochem(iexp, imod, :),\
                                           nochem_t(iexp, imod, :), plres)
          in = in + 1
        end do
      end do
    end if
    if (dim_OBS.gt.0) then
      do iobs1 = 0, dim_OBS - 1
        pmres = True
        pmres@gsMarkerSizeF = 0.01
        pmres@gsMarkerThicknessF = thicks_obs(iobs1)
        pmres@gsMarkerColor = colors_obs(iobs1)
        pmres@gsMarkerIndex = markers_obs(iobs1)
        plot_obs(iobs2) = gsn_add_polymarker(wks, plot(ipan), \
                                             time_obs(iobs1, :), \
                                             obs_t(iobs1, :),\
                                             pmres)
        iobs2 = iobs2 + 1
      end do
    end if

    delete(nochem_t)
    delete(chem_t)
    delete(xp)
    delete(yp)
    delete(ind_month)
  end do  ; loop over all "seasons"

  Pres = True
  Pres@gsnFrame = False
  Pres@gsnPanelBottom = 0.2
  if (dimsizes(seasons).gt.1) then
    gsn_panel(wks, plot, (/tointeger(dimsizes(seasons) / 2), 2/), Pres)
  else
    gsn_panel(wks, plot, (/tointeger(dimsizes(seasons)), 1 /), Pres)
  end if
  frame(wks)
  if (dim_OBS.gt.0) then
    add_legend(wks, list_obs_mod, colors_obs, markers_obs, lines_obs,\
               thicks_obs, "0.95")
  end if

  if (diag_script_info@legend_MMM.eq."True") then
    chem_att = new(dim_EXP, "string")
    nochem_att = new(dim_EXP, "string")
    do iexp = 0, dim_EXP - 1
      chem_att(iexp) = attrib_chem
      nochem_att(iexp) = "_NOCHEM"
    end do

    nochempos = 0.80 - (dim_chem_mod / 3) * 0.05
    if (nochem.gt.0) then
      add_legend(wks1, label_exp + nochem_att, \
                 colors_nochem(:, dim_nochem_mod), \
                 dashes_nochem(:, dim_nochem_mod), \
                 lines_nochem(:, dim_nochem_mod),  \
                 thicks_nochem(:, dim_nochem_mod), \
                 tostring(nochempos))
    end if

    if (chem.gt.0) then
      add_legend(wks1, label_exp + chem_att, \
                 colors_chem(:, dim_chem_mod), \
                 dashes_chem(:, dim_chem_mod), \
                 lines_chem(:, dim_chem_mod), \
                 thicks_chem(:, dim_chem_mod), \
                 "0.90")
    end if
  end if

  if (diag_script_info@legend.eq."True") then
    modelnames = names(index_mod)

    do i = 0, dimsizes(modelnames) - 1
      if (isvar("idx")) then
        delete(idx)
      end if
      idx = ind(list_chem_mod.eq.modelnames(i))
      if (all(ismissing(idx))) then
        continue
      end if
      do j = 0, dimsizes(idx) - 1
        if (ensembles(idx(j)).ne."r1i1p1") then
          if (isStrSubset(ensembles(idx(j)), "r") .and. \
              isStrSubset(ensembles(idx(j)), "i") .and. \
              isStrSubset(ensembles(idx(j)), "p")) then
            modelnames(idx(j)) = modelnames(idx(j)) + " ("  \
                                 + ensembles(idx(j)) + ")"
          end if
        end if
      end do
    end do

    add_legend(wks1, list_chem_mod, colors_chem(0, :dim_chem_mod - 1), \
               dashes_chem(0, :dim_chem_mod - 1), \
               lines_chem(0, :dim_chem_mod - 1),  \
               thicks_chem(0, :dim_chem_mod - 1), "0.85")
    if (nochem.gt.0) then
      add_legend(wks1, list_nochem_mod, \
                 colors_nochem(0, :dim_nochem_mod - 1), \
                 dashes_nochem(0, :dim_nochem_mod - 1), \
                 lines_nochem(0, :dim_nochem_mod - 1),  \
                 thicks_nochem(0, :dim_nochem_mod - 1), \
                 tostring(nochempos - 0.05))
    end if
  end if

  frame(wks1)

end
