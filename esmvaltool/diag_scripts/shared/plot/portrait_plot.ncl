; #############################################################################
; PLOTS SCRIPTS FOR PORTRAIT PLOTS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
;    function portrait_plot
;    function circle_plot
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"


; #############################################################################

undef("portrait_plot")
function portrait_plot(wks_in[1],
                       source,
                       varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source = (ndiag, nmod)
;    source!0 = diagnostics
;    source!1 = models
;
;    source = (ndiag, nmod, nref)
;    source(:,:,0) = reference dataset
;    source(:,:,1) = alternative dataset (optional)
;    source!0 = diagnostics
;    source!1 = models
;
; Return value
;    A graphic variable.
;
; Caveats
;
; Modification history
;    20151027-righi_mattia: added option for displaying rankings.
;    20140605-righi_mattia: modified with flexible plot shapes.
;    20140204-frank_franziska: extended.
;    20140114-righi_mattia: written.
;
local funcname, scriptname, wks, wks_in, data, source, res, atts, \
    base, varname
begin

  funcname = "portrait_plot"
  scriptname = "diag_scripts/shared/plot/portrait_plot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Use -999 as fill value (easier to display as text)
  data@_FillValue = -999.

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq.  defaults)) then
    grade = att2var(data, "var")
      DIAG_SCRIPT = att2var(data, "diag_script")
  else
    grade = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, grade)

  ; Check source dimensionality
  grade_dim = dimsizes(data)
  ndiag = grade_dim(0)
  nmod = grade_dim(1)

  if (dimsizes(grade_dim).eq.3) then
    if (grade_dim(2).le.2) then
      nref = grade_dim(2)
    else
      log_info("Maximum 2 references supported")
      log_info("Only the first 2 will be plotted")
      nref = 2
      var = data(:, :, 0:1)
      delete(data)
      data = var
      delete(var)
    end if
  elseif (dimsizes(grade_dim).eq.2) then
    nref = 1
    var = new((/ndiag, nmod, nref/), typeof(data))
    var(:, :, 0) = data
    delete(data)
    data = var
    delete(var)
  else
    error_msg("f", scriptname, funcname, \
              "input variable must be eithed 2D or 3D")
  end if

  ; Set levels
  if (isatt(data, "label_bounds")) then

    ; User provided
    chars = \
      tochar(tostring(abs(data@label_bounds(1) - \
             round(data@label_bounds(1), 1)) * 100))
    five = tochar("5")
    zero = tochar("0")
    if (chars(1).eq.five(0)).or.(chars(1).eq.zero(0)) then
      dec_pl = 100.
    else
      dec_pl = 10.
    end if
    upperbound = tofloat(round(data@label_bounds(1) * dec_pl, 1)) / dec_pl
    lowerbound = tofloat(round(data@label_bounds(0) * dec_pl, 1)) / dec_pl

  else

    ; Flexible in regard to the range of the data values
    bound  = max(abs(data(:, :, 0)))
    if (all(ismissing(data))) then
      delete(bound)
      bound = 0.1
      error_msg("w", scriptname, funcname, "all data values are missing!")
    end if
    upperbound = tofloat(round(bound * 10, 1)) / 10.
    upperbound = where(upperbound.lt.0.1, 0.1, upperbound)
    if (min(data).lt.0) then
      lowerbound = -1. * upperbound
    else
      lowerbound = 0
    end if

  end if

  ; Label scale
  if (isatt(data, "label_scale")) then
    scale = data@label_scale
  else
    if (((upperbound - lowerbound) / 0.5).le.6) then
      scale = 0.1
    else
      scale = 0.5
    end if
  end if

  ; Set levels
  if (isatt(data, "label_lo")) then
    label_lo = data@label_lo
  else
    label_lo = True
  end if
  if (isatt(data, "label_hi")) then
    label_hi = data@label_hi
  else
    label_hi = True
  end if
  nlevels = round((upperbound - lowerbound + scale) / scale, 3)
  levels = new(nlevels, "float")
  levels(0) = lowerbound
  do ii = 1, nlevels - 1
    levels(ii) = levels(ii - 1) + scale
  end do
  levels = round(100. * levels, 0) / 100.

  ; Display metrics value on the plot
  flag_values = False
  if (isatt(data, "disp_values")) then
    if (data@disp_values) then
      flag_values = True
      text_v = new((/ndiag, nmod, nref/), graphic)
    end if
  end if

  ; Diplay model rankings on the plot
  flag_ranks = False
  if (isatt(data, "disp_rankings")) then
    if (data@disp_rankings) then
      flag_ranks = True
      text_r = new((/ndiag, nmod, nref/), graphic)
      ranks = new((/ndiag, nmod, nref/), integer)
    end if
  end if

  ; Set a colormap
  if (isatt(data, "colormap")) then
    gsn_define_colormap(wks, data@colormap)
  else
    gsn_define_colormap(wks, "ncl_default")
    data@colormap = "ncl_default"
  end if
  tmp = read_colormap_file(data@colormap)
  if (isatt(data, "cm_interval")) then
    cm_start = data@cm_interval(0)
    cm_end = data@cm_interval(1)
  else
    cm_start = 2
    cm_size = dimsizes(tmp)
    cm_end = cm_size(0) - 1
  end if
  cmap = tmp(cm_start:cm_end, :)
  ncmap = dimsizes(cmap(:, 0))
  delete(tmp)

  ; Sample colormap and extract colors
  if (label_lo.and.label_hi) then  ; 2 more colors for the triangles
    pos = new((/nlevels + 1/), float)
  end if
  if (label_lo.xor.label_hi) then  ; 1 more color for upper/lower triangle
    pos = new((/nlevels/), float)
  end if
  if (.not.label_lo.and. .not.label_hi)  ; no more colors
    pos = new((/nlevels - 1/), float)
  end if
  pos(0) = cm_start
  npos = dimsizes(pos)
  do ii = 1, npos - 1
    pos(ii) = pos(ii - 1) + (ncmap - 1) / tofloat(npos - 1)
  end do
  delete([/cmap, ncmap/])

  ; Reverse order if required
  if (isatt(data, "cm_reverse")) then
    if (data@cm_reverse) then
      pos = pos(::-1)
    end if
  end if

  ; Default resources
  res = True
  res@gsnDraw  = False
  res@gsnFrame = False
  res@tmBorderThicknessF = False
  res@tmXBOn = False
  res@tmXTOn = False
  res@tmYLOn = False
  res@tmYROn = False
  res@trXMinF = 0.
  res@trXMaxF = 1.0
  res@trYMinF = 0.
  res@trYMaxF = 1.0

  ; Set plot aspect depending on ndiag and nmod: longst dimension on y-axis
  xbsize = 1. / nmod
  ybsize = 1. / ndiag
  if (ndiag.ge.nmod) then
    res@vpWidthF = 0.6 * nmod / tofloat(ndiag)
    res@vpHeightF = 0.6
  else
    res@vpWidthF = 0.6
    if (isatt(diag_script_info, "scale_aspect")) then
      ; code for ESA CCI paper Lauer et al. (2016)
      scale_aspect = diag_script_info@scale_aspect
      res@vpHeightF = 0.6 * ndiag / tofloat(nmod) * scale_aspect
    else
      ; standard code
      res@vpHeightF = 0.6 * ndiag / tofloat(nmod)
      if (2 * ndiag.le.nmod) then
        res@vpHeightF = 2 * res@vpHeightF
      end if
    end if
  end if
  res@vpYF = 0.2 + res@vpHeightF

  ; Title
  res@tiMainFontHeightF = 0.020
  if (isatt(data, "metric")) then
    res@tiMainString = data@metric
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Draw a blank plot, polygons will be added below
  plot = gsn_blank_plot(wks, res)

  ; Text resources
  resTxt = True
  resTxt@txFontHeightF = 0.010
  resTxt@txJust = "CenterCenter"

  ; Polygons resources
  resPol = True
  resPol@gsLineThicknessF = 0.5
  resPol@gsEdgesOn   = True
  resPol@gsLineColor = "black"
  poly = new((/ndiag, nmod, nref/), graphic)

  ; Calculate rankings
  if (flag_ranks) then
    if (.not.isatt(data, "rank_order")) then
      error_msg("f", scriptname, funcname, \
                "rank_order attribute must be specified if " + \
                "ranking display is desired")
    end if
    if (all(data@rank_order.ne.(/-1, 1/))) then
      error_msg("f", scriptname, funcname, \
                "rank_order attribute can be either 1 (increasing) " + \
                "or -1 (decreasing)")
    end if
    do idiag = 0, ndiag - 1
      do iref = 0, nref - 1
        pp = dim_pqsort(data(idiag, :, iref), data@rank_order)
        rk = 1
        do imod = 0, nmod - 1
          if (.not.ismissing(data(idiag, pp(imod), iref))) then
            ranks(idiag, pp(imod), iref) = rk
            rk = rk + 1
          end if
        end do
        delete(rk)
        delete(pp)
      end do
    end do
  end if

  ; Draw polygons
  do idiag = 0, ndiag - 1

    ; Set triangle or square
    ;  - draw a square if only 1 reference model is available or if the
    ;    alternative model contains only missing values
    ;  - draw a triangle if an alternative reference is available and
    ;    contains valid values
    l_square = True
    if (nref.gt.1) then
      if (.not.all(ismissing(data(idiag, :, 1)))) then
        l_square = False
      end if
    end if

    do imod = 0, nmod - 1

      ; Find corresponding color for primary reference
      curr = data(idiag, imod, 0)
      if (ismissing(curr)) then
        resPol@gsFillColor = 0
      else
        lidx = max(ind(curr.ge.levels))
        if (label_lo) then
          lidx = lidx + 1
          lidx = where(ismissing(lidx), 0, lidx)
        end if
        if (.not.label_hi .and. lidx.gt.npos - 1) then
          error_msg("f", scriptname, funcname, \
                    "metric value (" + curr + ") above the " + \
                    "label upper bound (" + levels(nlevels - 1) + \
                    "): change label_bounds or set " + \
                    "label_hi = True")
        end if
        if (ismissing(lidx)) then
          error_msg("f", scriptname, funcname, \
                    "metric value (" + curr + ") below the " + \
                    "label lower bound (" + levels(0) + \
                    "): change label_bounds or set " + \
                    "label_lo = True")
        end if
        resPol@gsFillColor = round(pos(lidx), 3)
      end if

      if (l_square) then
        xx = (/imod * xbsize, (imod + 1) * xbsize, \
              (imod + 1) * xbsize, imod * xbsize, imod * xbsize/)
        yy = (/idiag * ybsize, idiag * ybsize, (idiag + 1) * ybsize, \
              (idiag + 1) * ybsize, idiag * ybsize/)
        poly(idiag, imod, 0) = gsn_add_polygon(wks, plot, xx, yy, resPol)
        if (flag_values) then
          xt = (xx(0) + xx(1)) / 2.
          yt = (yy(1) + yy(2)) / 2
          text_v(idiag, imod, 0) = \
            gsn_add_text(wks, plot, sprintf("%4.2f", data(idiag, imod, 0)), \
                         xt, yt, resTxt)
          delete(xt)
          delete(yt)
        end if
        if (flag_ranks.and..not.ismissing(data(idiag, imod, 0))) then
          xt = (xx(0) + xx(1)) / 2.
          yt = (yy(1) + yy(2)) / 2.
          text_r(idiag, imod, 0) = \
            gsn_add_text(wks, plot, "#" + ranks(idiag, imod, 0), \
                         xt, yt, resTxt)
          delete(xt)
          delete(yt)
        end if
        delete(xx)
        delete(yy)

      else

        ; Primary reference => lower triangle
        xx = (/imod * xbsize, (imod + 1) * xbsize, \
              (imod + 1) * xbsize, imod * xbsize/)
        yy = (/idiag * ybsize, idiag * ybsize, (idiag + 1) * ybsize, \
              idiag * ybsize/)
        poly(idiag, imod, 0) = gsn_add_polygon(wks, plot, xx, yy, resPol)
        if (flag_values) then
          xt = sum(xx(0:2)) / 3.
          yt = sum(yy(0:2)) / 3.
          text_v(idiag, imod, 0) = \
            gsn_add_text(wks, plot, sprintf("%4.2f", data(idiag, imod, 0)), \
                         xt, yt, resTxt)
          delete(xt)
          delete(yt)
        end if
        if (flag_ranks.and..not.ismissing(data(idiag, imod, 0))) then
          xt = sum(xx(0:2)) / 3.
          yt = sum(yy(0:2)) / 3.
          text_r(idiag, imod, 0) = \
            gsn_add_text(wks, plot, "#" + ranks(idiag, imod, 0), \
                         xt, yt, resTxt)
          delete(xt)
          delete(yt)
        end if
        delete(xx)
        delete(yy)

        ; Find corresponding color for alternative reference
        curr = data(idiag, imod, 1)
        if (ismissing(curr)) then
          resPol@gsFillColor = 0
        else
          lidx = max(ind(curr.ge.levels))
          if (label_lo) then
            lidx = lidx + 1
            lidx = where(ismissing(lidx), 0, lidx)
          end if
          if (.not.label_hi .and. lidx.gt.npos - 1) then
            error_msg("f", scriptname, funcname, \
                      "metric value (" + curr + ") above the " + \
                      "label upper bound (" + \
                      levels(nlevels - 1) + \
                      "): change label_bounds or set " + \
                      "label_hi = True")
          end if
          if (ismissing(lidx)) then
            error_msg("f", scriptname, funcname, \
                      "metric value (" + curr + ") below the " + \
                      "label lower bound (" + levels(0) + \
                      "): change label_bounds or set " + \
                      "label_lo = True")
          end if
          resPol@gsFillColor = round(pos(lidx), 3)
        end if

        ; Alternative reference => upper triangle
        xx = (/imod * xbsize, (imod + 1) * xbsize, \
              imod * xbsize, imod * xbsize/)
        yy = (/idiag * ybsize, (idiag + 1) * ybsize, \
              (idiag + 1) * ybsize, idiag * ybsize/)
        poly(idiag, imod, 1) = gsn_add_polygon(wks, plot, xx, yy, resPol)
        if (flag_values) then
          xt = sum(xx(0:2)) / 3.
          yt = sum(yy(0:2)) / 3.
          text_v(idiag, imod, 1) = \
            gsn_add_text(wks, plot, sprintf("%4.2f", data(idiag, imod, 1)), \
                         xt, yt, resTxt)
          delete(xt)
          delete(yt)
        end if
        if (flag_ranks.and..not.ismissing(data(idiag, imod, 1))) then
          xt = sum(xx(0:2)) / 3.
          yt = sum(yy(0:2)) / 3.
          text_r(idiag, imod, 1) = \
            gsn_add_text(wks, plot, "#" + ranks(idiag, imod, 1), \
                         xt, yt, resTxt)
          delete(xt)
          delete(yt)
        end if
        delete(xx)
        delete(yy)

      end if
      delete(resPol@gsFillColor)
    end do
  end do

  res@poly = poly
  if (isdefined("text_v")) then
    res@text_v = text_v
  end if
  if (isdefined("text_r")) then
    res@text_v = text_r
  end if

  ; Add legend text for the diagnostics
  if(isatt(diag_script_info, "scale_font")) then
    scale_font = diag_script_info@scale_font
  else
    scale_font = 1.
  end if
  resTxt@txFontHeightF = min((/0.012 * scale_font, \
                              0.57 * scale_font / max((/ndiag, nmod/))/))
  resTxt@txJust = "CenterRight"
  step = res@vpHeightF / ndiag
  xpos = 0.195
  ypos = 0.2 + step / 2.
  do idiag = 0, ndiag - 1
    gsn_text_ndc(wks, data&diagnostics(idiag), xpos, ypos, resTxt)
    ypos = ypos + step
  end do

  ; Add text for the models
  resTxt@txAngleF = 90.
  step = res@vpWidthF / nmod
  resTxt@txJust = "CenterRight"
  xpos = 0.2 + step / 2.
  ypos = 0.19
  do imod = 0, nmod - 1
    gsn_text_ndc(wks, data&models(imod), xpos, ypos, resTxt)
    xpos = xpos + step
  end do

  ; Draw a labelbar
  resTxt@txFontHeightF = 0.015 * scale_font
  resTxt@txAngleF = 0.
  resTxt@txJust = "CenterLeft"
  height = res@vpHeightF / npos
  xp1 = 0.2 + res@vpWidthF + 0.01
  xp2 = 0.2 + res@vpWidthF + 0.04
  xpm = (xp1 + xp2) / 2.
  ypos = 0.2

  il = 0
  do ii = 0, npos - 1
    resPol@gsFillColor = pos(ii)
    if (ii.eq.0) then
      if (label_lo) then
        gsn_polygon_ndc(wks, (/xpm, xp1, xp2, xpm/),\
                        (/ypos, ypos + height, \
                         ypos + height, ypos/), resPol)
      else
        gsn_polygon_ndc(wks, (/xp1, xp1, xp2, xp2, xp1/),\
                        (/ypos, ypos + height, \
                         ypos + height, ypos, ypos/), resPol)
        gsn_text_ndc(wks, sprintf("%4.2f", levels(il)), \
                     xp2 + 0.01, ypos, resTxt)
        il = il + 1
      end if
    elseif (ii.eq.(npos - 1)) then
      if (label_hi) then
        gsn_polygon_ndc(wks, (/xp1, xpm, xp2, xp1/), \
                        (/ypos, ypos + height, ypos, ypos/), resPol)
        gsn_text_ndc(wks, sprintf("%4.2f", levels(il)), \
                     xp2 + 0.01, ypos, resTxt)
        il = il + 1
      else
        gsn_polygon_ndc(wks, (/xp1, xp1, xp2, xp2, xp1/),\
                        (/ypos, ypos + height, \
                         ypos + height, ypos, ypos/), resPol)
        gsn_text_ndc(wks, sprintf("%4.2f", levels(il)), \
                     xp2 + 0.01, ypos, resTxt)
        il = il + 1
        gsn_text_ndc(wks, sprintf("%4.2f", levels(il)), \
                     xp2 + 0.01, ypos + height, resTxt)
        il = il + 1
      end if
    else
      gsn_polygon_ndc(wks, (/xp1, xp1, xp2, xp2, xp1/), \
                      (/ypos, ypos + height, ypos + height, ypos, ypos/), \
                      resPol)
      gsn_text_ndc(wks, sprintf("%4.2f", levels(il)), \
                   xp2 + 0.01, ypos, resTxt)
      il = il + 1
    end if
    ypos = ypos + height
  end do
  draw(plot)
  frame(wks)

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################

undef("circle_plot")
function circle_plot(wks_in[1],
                     source,
                     varname[1]: string,
                     ensname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source = (nmod, 2)
;    source!0 = models
;    source!1 = statistic
;    source(:, 0) = SMPI
;    source(:, 1) = 95% conf radius
;
; Return value
;    A graphic variable.
;
; Caveats
;
; Modification history
;    20170424-gier_bettina: written.
;
local funcname, scriptname, verbosity, wks, wks_in, data, source, res, atts, \
  base, varname
begin

  funcname = "circle_plot"
  scriptname = "plot_scripts/ncl/portrait_plot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Use -999 as fill value (easier to display as text)
  data@_FillValue = -999.

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname.eq.defaults)) then
    performance_index = att2var(data, "var")
    DIAG_SCRIPT = att2var(data, "diag_script")
  else
    performance_index = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, performance_index)

  ; Make list of markers for different projects
  symbols_1 = array_append_record(inttochar(ispan(97, 122, 1)), \
                                  inttochar(ispan(65, 90, 1)), 0)

  ; Greek alphabet
  symbols_2 = new(dimsizes(symbols_1), string)
  do ii = 0, dimsizes(symbols_1) - 1
    symbols_2(ii) = "~F5~" + symbols_1(ii) + "~F~"
  end do

  n_projects = 1
  x_min = min((/0.38, min(data(:, 0) - data(:, 1))-0.05/))
  x_max = max(data(:, 0) + data(:, 1)) + 0.25
  y_max = n_projects + 2
  y_min = n_projects - 2

  ; Default resources
  res = True
  res@xyXStyle = "Log"
  res@gsnDraw  = False
  res@gsnFrame = False
  res@tmBorderThicknessF = False
  res@tmXBBorderOn = False
  res@tmXTBorderOn = False
  res@tmYLBorderOn = False
  res@tmYRBorderOn = False
  res@tmXBOn = False
  res@tmXTOn = False
  res@tmYLOn = False
  res@tmYROn = False
  res@trXMinF = x_min
  res@trXMaxF = x_max
  res@trYMinF = y_min
  res@trYMaxF = y_max
  res@vpXF = 0.1
  res@vpWidthF = 0.9
  res@vpYF = 1.0
  res@vpHeightF = 0.4
  res@tiXAxisString = "I~S~2"
  res@tiXAxisOffsetYF = 0.1

  res_lines = True
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineThicknessF = 2.
  res_lines@tfPolyDrawOrder = "PreDraw"

  res_circles = True
  res_circles@gsEdgesOn = True
  res_circles@gsEdgeThicknessF = 1.5

  res_text = True
  res_text@txFontHeightF = 0.02

  res_marker = True
  res_marker@gsMarkerSizeF = 0.02
  res_marker@gsMarkerIndex = 1
  res_marker@gsMarkerColor = "black"
  res_marker@tfPolyDrawOrder = "Draw"
  res_marker@txFontHeightF = 0.015
  res_marker@txFontThicknessF = 2

  plot = gsn_csm_xy(wks, (/x_min, x_max/), (/1., 1./), res)

  ; Draw baseline
  plot@$unique_string("dum")$ = \
    gsn_add_polyline(wks, plot, (/x_min, x_max/), (/1., 1./), res_lines)

  ; Draw ticks on baseline
  x_ticks = 0.1 * ispan(toint(x_min * 10), toint(x_max * 10), 1)
  x_big_ticks = 0.1 * ispan(0, toint(x_max * 10), 5)
  x_numbers = ispan(toint(x_min) + 1, toint(x_max), 1)

  do ticks = 0, dimsizes(x_ticks) - 1
    plot@$unique_string("dum")$ = \
      gsn_add_polyline(wks, plot, (/x_ticks(ticks), x_ticks(ticks)/), \
                       (/0.95, 1.05/), res_lines)
  end do

  do ticks = 0, dimsizes(x_big_ticks) - 1
    plot@$unique_string("dum")$ = \
      gsn_add_polyline(wks, plot, \
                       (/x_big_ticks(ticks), x_big_ticks(ticks)/), \
                       (/0.85, 1.15/), res_lines)
  end do

  do ticks = 0, dimsizes(x_numbers)-1
    plot@$unique_string("dum")$ = \
      gsn_add_text(wks, plot, tostring(x_numbers(ticks)), \
                   x_numbers(ticks), 0.2, res_text)
  end do

  ; Add models as circles, transform x-coords into ndc space
  ; to keep circles in log
  n = 50  ; nr of points to draw circle with
  do imod = 0, dimsizes(data&models) - 1
    x_in = (/data(imod, 0), data(imod, 0) + data(imod, 1), data(imod, 0)/)
    y_in = (/1., 1., 1.26/)
    x_ndc = new(3, float)
    y_ndc = new(3, float)
    datatondc(plot, x_in, y_in, x_ndc, y_ndc)
    mod_r = x_ndc(1) - x_ndc(0)
    x_1 = fspan(-mod_r, mod_r, n)
    y_1 = sqrt(mod_r ^ 2 - x_1 ^ 2)
    x = array_append_record(x_1, x_1(::-1), 0)
    y = array_append_record(y_1, -y_1(::-1), 0)
    x = x + x_ndc(0)
    y = y + y_ndc(0)
    if (data&models(imod).eq."multi-model-mean") then
      res_circles@gsFillColor = "black"
      res_circles@gsFillOpacityF = 0.8
    else
      res_circles@gsFillColor = "orange"
      res_circles@gsFillOpacityF = 0.3
    end if
    gsn_polygon_ndc(wks, x, y, res_circles)
    gsn_polymarker_ndc(wks, x_ndc(0), y_ndc(0), res_marker)
    delete([/mod_r, x_1, y_1, x, y, x_ndc, y_ndc, x_in, y_in/])
  end do

  ; Plot average I^2 and Project name (Multi Model Mean at the end)
  x_in = (/1., 1./)  ; arbitrary value
  y_in = (/1., 1./)
  x_ndc = new(2, float)
  y_ndc = new(2, float)
  datatondc(plot, x_in, y_in, x_ndc, y_ndc)
  res_text@txJust = "CenterLeft"
  gsn_text_ndc(wks, ensname, 0.01, y_ndc(0), res_text)
  delete([/x_in, y_in, x_ndc, y_ndc/])

  ; Text symbols after all circles are drawn, in case they fall inside
  x_in = data(:, 0)
  y_in = new(dimsizes(data(:, 0)), float)
  y_in(:) = max(data(:, 1)) + 1.5
  x_ndc = new(dimsizes(data(:, 0)), float)
  y_ndc = new(dimsizes(data(:, 0)), float)
  datatondc(plot, x_in, y_in, x_ndc, y_ndc)

  do imod = 0, dimsizes(data&models) - 1
    gsn_text_ndc(wks, symbols_2(imod), x_ndc(imod), y_ndc(imod), res_marker)
  end do
  delete([/x_in, y_in, x_ndc, y_ndc/])

  ; Draw legend
  txres = True
  txres@txFontHeightF = 0.015
  txres@txJust = "CenterLeft"
  ypos = fspan(0.55, 0.1, 20)
  half = toint(dimsizes(data&models) / 2)
  do imod = 0, half
    gsn_text_ndc(wks, symbols_2(imod) + ": " + \
                 data&models(imod), 0.01, ypos(imod), txres)
  end do
  do imod = half + 1, dimsizes(data&models) - 1
    gsn_text_ndc(wks, symbols_2(imod) + ": " + \
                 data&models(imod), 0.3, ypos(imod - half - 1), txres)
  end do

  draw(plot)
  frame(wks)

  leave_msg(scriptname, funcname)
  return(plot)

end
