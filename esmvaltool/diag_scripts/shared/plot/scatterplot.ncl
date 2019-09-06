; #############################################################################
; PLOT SCRIPTS FOR SCATTER PLOTS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
;    function scatterplot
;    function scatterplot3D
;    function scatterplot_markers
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

; #############################################################################
undef("scatterplot")
function scatterplot(wks_in[1],
                     source,
                     varname[1]: string,
                     reflines: logical,
                     stats: logical,
                     inlist:list)
;
; Arguments
;    wks_in: workstation ("graphic" object or default will be used).
;    source: data to be plotted (see source prototype above) or a NetCDF
;            filename with data.
;    varname: variable name in the file.
;    reflines: show/hide reference lines (1:1 and +/- factor 2).
;    stats: show/hide summary of statistical values on the plot.
;    inlist: list of input_file_info items
;
; Source prototype
;    source = (2, npoints)
;    source(0, :) = x-axis values
;    source(1, :) = y-axis values
;    source!0 = datasets
;
; Return value
;    A graphic object.
;
; Description
;    Creates a scatter plot and optionally adds lines for the factor-2 range
;    and statistical information.
;
; Caveats
;    Linear axes are used as default, log axis can be optionally provided by
;    the calling diag_script (option: diag_script_info@scatter_log).
;
; Modification history
;    20140228-righi_mattia: written.
;
local funcname, scriptname, scattermode, loc_res, resL, linec, lineh, linel, \
  resT, nstats, st, text, step, xpos, ypos, ii, tmp, ratio, fac2p, str,
  scatter_log
begin

  funcname = "scatterplot"
  scriptname = "diag_scripts/shared/plot/scatterplot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF
  if (typeof(source).eq."string") then
    data = ncdf_read(source, varname)
  else
    data = source
  end if

  ; Get scatter_log option
  scatter_log = False
  if (isatt(diag_script_info, "scatter_log")) then
    if (diag_script_info@scatter_log .eq. "True") then
      scatter_log = True
    end if
  end if

  ; Style options
  colors  = project_style(inlist, diag_script_info, "colors")

  ; Define x and y points
  xpts = data(0, :)
  ypts = data(1, :)

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname.eq.defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, var)
  gsn_define_colormap(wks, array_append_record((/"white", "black"/), \
                      colors, 0))

  ; Default resources
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@xyMarkLineMode = "Markers"
  res@xyMarker = 9
  res@tiMainFont = 22
  res@tiMainString  = var
  if (isatt(data, "long_name")) then
    res@gsnCenterString = data@long_name
  else
    res@gsnCenterString = var
  end if
  if (isatt(data, "units")) then
    res@gsnCenterString  = res@gsnCenterString + " [" + \
      format_units(data@units) + "]"
  end if
  res@tiXAxisString = data&datasets(0)
  res@tiYAxisString = data&datasets(1)
  if (scatter_log) then
    res@trXAxisType = "LogAxis"
    res@trYAxisType = "LogAxis"
    xpts = where(xpts.le.0, xpts@_FillValue, xpts)
    ypts = where(ypts.le.0, ypts@_FillValue, ypts)
    res@trXMinF = 10 ^ ( \
      min((/toint(log10(min(xpts))), toint(log10(min(ypts)))/)) - 1)
    res@trXMaxF = 10 ^ ( \
      max((/toint(log10(max(xpts))), toint(log10(max(ypts)))/)) + 1)
    res@tmXBMinorPerMajor = 8
    res@tmXTMinorPerMajor = 8
    res@tmYLMinorPerMajor = 8
    res@tmYRMinorPerMajor = 8
  else
    res@trXMinF = 0.
    res@trXMaxF = 1.10 * max((/max(xpts), max(ypts)/))
  end if
  res@trYMinF = res@trXMinF
  res@trYMaxF = res@trXMaxF

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  plot = gsn_csm_xy(wks, xpts, ypts, res)

  ; Add reference lines
  if (reflines) then

    resL = True
    resL@gsLineThicknessF = 1.5
    linec = gsn_add_polyline(wks_in, plot, \
                             (/res@trXMinF, res@trXMaxF/), \
                             (/res@trYMinF, res@trYMaxF/), resL)
    plot@linec = linec
    resL@gsLineDashPattern = 1
    if (scatter_log) then
      lineh = gsn_add_polyline(wks_in, plot, \
                               (/res@trXMinF, res@trXMaxF / 2./), \
                               (/2 * res@trYMinF, res@trYMaxF/), resL)
      linel = gsn_add_polyline(wks_in, plot, \
                               (/2 * res@trXMinF, res@trXMaxF/), \
                               (/res@trYMinF, res@trYMaxF / 2/), resL)
    else
      lineh = gsn_add_polyline(wks_in, plot, \
                               (/res@trXMinF, res@trXMaxF / 2./), \
                               (/res@trYMinF, res@trYMaxF/), resL)
      linel = gsn_add_polyline(wks_in, plot, \
                               (/res@trXMinF, res@trXMaxF/), \
                               (/res@trYMinF, res@trYMaxF / 2/), resL)
    end if
    plot@lineh = lineh
    plot@linel = linel

  end if

  ; Add legend with statistical values
  if (stats) then

    resT               = True
    resT@txFontHeightF = 0.016
    resT@txJust        = "CenterLeft"

    nstats = 7
    st = new(nstats, string)
    text = new(nstats, graphic)

    ; Define labels horizontal position
    if (scatter_log) then
      step = (log10(res@trYMaxF) - log10(res@trYMinF)) / 20.
      xpos = 10 ^ (log10(res@trYMinF) + step)
    else
      step = (res@trYMaxF - res@trYMinF) / 20.
      xpos = res@trYMinF + step
    end if

    ; Define labels vertical position
    ypos = new(nstats, float)
    do ii = 0, nstats - 1
      if (scatter_log) then
        tmp = log10(res@trYMaxF) - step * (ii + 1)
        ypos(ii) = 10 ^ tmp
        delete(tmp)
      else
        ypos(ii) = res@trYMaxF - step * (ii + 1)
      end if
    end do

    ; Define substrings for labeling
    x_substr = "~B" + strlen(res@tiXAxisString) + "~" + res@tiXAxisString + \
      empty_str(strlen(res@tiXAxisString) + 2)

    y_substr = "~B" + strlen(res@tiYAxisString) + "~" + res@tiYAxisString + \
      empty_str(strlen(res@tiYAxisString) + 2)

    ; Number of points
    st(0) = "N~B3~pts     = " + sprintf("%5.0f", num(.not. ismissing(xpts)))

    ; x-pts average +/- stddev
    st(1) = res@tiXAxisString + " = " + \
      sprintf("%4.2f", avg(xpts)) + " ~F18~S~F21~ " + \
      sprintf("%4.2f", stddev(xpts))

    ; y-pts average +/- stddev
    st(2) = res@tiYAxisString + " = " + \
      sprintf("%4.2f", avg(ypts)) + " ~F18~S~F21~ " + \
      sprintf("%4.2f", stddev(ypts))

    ; Ratio of the averages
    st(3) = res@tiYAxisString + " / " + res@tiXAxisString + " = " + \
      sprintf("%4.2f", avg(ypts) / avg(xpts))

    ; Fraction of model points within factor 2 of the observations
    ratio = xpts / ypts
    fac2p = tofloat(num(ratio .ge. 0.5 .and. ratio .le. 2.0)) \
      / tofloat(num(.not. ismissing(ratio))) * 100
    st(4) = "FAC2 = " + sprintf("%4.2f", fac2p) + "%"

    ; Root-mean-square error
    st(5) = "RMSE = " + sprintf("%5.2f", dim_rmsd(xpts, ypts))

    ; Correlation coefficient
    corr = escorc(xpts, ypts)
    st(6) = "R~S1~2  = " + sprintf("%3.2f", corr)

    ; Add labels
    do ii = 0, nstats - 1
      text(ii) = gsn_add_text(wks_in, plot, st(ii), xpos, ypos(ii), resT)
      str = "text" + ii
      plot@$str$ = text(ii)
    end do

  end if

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("scatterplot3D")
function scatterplot3D(wks_in[1],
                       source,
                       varname[1]:string,
                       levels[*]:numeric,
                       inlist:list)
;
; Arguments
;    wks_in: workstation ("graphic" object or default will be used).
;    source: data to be plotted (see source prototype above) or a NetCDF
;            filename with data.
;    varname: variable name in the file.
;    levels: levels for the color coding of the third dimension.
;
; Source prototype
;    source = (3, npoints)
;    source(0, :) = x-axis values
;    source(1, :) = y-axis values
;    source(2, :) = z-axis values
;    source!0 = statistic
;    source!1 = datasets
;
; Return value
;    A graphic object.
;
; Description
;    Creates a 3D scatter plot using color coded markers as third dimension
;    and statistical information.
;
; Caveats
;
; Modification history
;    20151103-righi_mattia: written.
;
local funcname, scriptname, data, xpts, ypts, zpts, defaults, var, wks, res, \
  res_new, colors, ncolors, resM, markers, resT, labels, imod, idx, xp, yp, \
  xb, yb, resB, lbid, resA, annoid
begin

  funcname = "scatterplot3D"
  scriptname = "diag_scripts/shared/plot/scatterplot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF
  if (typeof(source).eq."string") then
    data = ncdf_read(source, varname)
  else
    data = source
  end if

  ; Define x and y points
  xpts = data(0, :)
  ypts = data(1, :)
  zpts = data(2, :)

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname.eq.defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, var)

  ; Create a 2D scatter plot
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@vpHeightF = 0.5
  res@xyMarkLineMode = "Markers"
  res@xyMarkerSizeF = 1.e-6
  res@xyMarker = 16
  res@xyMarkerColors = 14
  res@tiMainFontHeightF = 0.025
  res@tiYAxisFontHeightF = 0.025
  res@tiXAxisFontHeightF = 0.025
  res@txFontHeightF = 0.01
  res@tmXBLabelFontHeightF = 0.02
  res@tmYLLabelFontHeightF = 0.02
  res@tiMainString = varname
  res@tiXAxisString = data&statistic(0)
  res@tiYAxisString = data&statistic(1)
  if (isatt(data, "units")) then
    res@tiXAxisString = res@tiXAxisString + \
      " [" + format_units(data@units(0)) + "]"
    res@tiYAxisString = res@tiYAxisString + \
      " [" + format_units(data@units(1)) + "]"
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  plot = gsn_csm_xy(wks, (/xpts, xpts/), (/ypts, ypts/), res)

  ; Get nice spacing through color map for marker colors
  getvalues wks
    "wkColorMapLen" : clen  ; # of colors in color map
  end getvalues
  colors = round(fspan(2, clen - 1, dimsizes(levels) - 1), 3)
  ncolors = dimsizes(colors)

  ; Add color-coded markers and text
  resM = True
  resM@gsMarkerSizeF = 0.015
  resM@gsMarkerIndex = 16
  markers = new(dimsizes(data(0, :)), graphic)

  resT = True
  resT@txFontHeightF = 0.012
  resT@txFontColor = "black"
  labels = new(dimsizes(data(0, :)), graphic)

  do imod = 0, dimsizes(data(0, :)) - 1

    ; Use a black marker for the reference dataset (i.e., MVI = 0)
    if (zpts(imod).ne.0) then
      idx = max(ind(zpts(imod).ge.levels))
      resM@gsMarkerColor = colors(idx)
    else
      resM@gsMarkerColor = "black"
    end if
    markers(imod) = \
      gsn_add_polymarker(wks, plot, xpts(imod), ypts(imod), resM)
    delete(resM@gsMarkerColor)
    str = "mark" + imod
    plot@$str$ = markers(imod)

    ; Draw text label (right or left justified depending on the position)
    xp = xpts(imod)
    yp = ypts(imod)
    xb = 0.  ; initialize
    yb = 0.  ; initialize
    datatondc(plot, xp, yp, xb, yb)
    if (xb.le.0.7) then
      resT@txJust = "CenterLeft"
      xb = xb + 0.01
    else
      resT@txJust = "CenterRight"
      xb = xb - 0.01
    end if
    ndctodata(plot, xb, yb, xp, yp)
    labels(imod) = gsn_add_text(wks, plot, data&datasets(imod), xp, yp, resT)
    str = "text" + imod
    plot@$str$ = labels(imod)

  end do

  ; Add legend
  resB                    = True
  resB@lbAutoManage       = False
  resB@lbFillColors       = colors
  resB@vpWidthF           = 0.15
  resB@vpHeightF          = res@vpHeightF
  resB@lbMonoFillPattern  = True
  resB@lbLabelFontHeightF = 0.018
  resB@lbOrientation      = "vertical"
  resB@lbPerimOn          = False
  resB@lbLabelAlignment   = "ExternalEdges"
  resB@lbTitleFontHeightF = 0.018
  resB@lbTitleString = data&statistic(2)
  if (isatt(data, "units")) then
    resB@lbTitleString = resB@lbTitleString + \
      " [" + format_units(data@units(2)) + "]"
  end if
  lbid = gsn_create_labelbar(wks, ncolors, sprintf("%4.2f", levels), resB)

  resA                    = True
  resA@amParallelPosF     =  0.65
  annoid = gsn_add_annotation(plot, lbid, resA)

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("scatterplot_markers")
function scatterplot_markers(wks_in[1],
                             source,
                             varname[1],
                             inlist:list)
;
; Arguments
;    wks_in: workstation ("graphic" object or default will be used).
;    source: data to be plotted (see source prototype above) or a NetCDF
;            filename with data.
;    varname: variable name in the file.
;    inlist: list of input_file_info items
;
; Source prototype
;    source = (2, npoints)
;    source(0, :) = x-axis values
;    source(1, :) = y-axis values
;    source!0 = statistic
;    source!1 = model
;
; Return value
;    A graphic object.
;
; Description
;    Creates a scatter plot using different markers/colors for each point,
;    according to the given styleset (default is used if not specified).
;
; Caveats
;
; Modification history
;    20151105-righi_mattia: written.
;
local funcname, scriptname, data, xpts, ypts, defaults, var, wks, annots, \
  colors, markers, res, res_new, lg_outfile, styles, lgres
begin

  funcname = "scatterplot_markers"
  scriptname = "diag_scripts/shared/plot/scatterplot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF
  if (typeof(source).eq."string") then
    data = ncdf_read(source, varname)
  else
    data = source
  end if

  ; Define x and y points (extra dimension for markers/colors specification)
  xpts = data(0, :)
  ypts = data(1, :)

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname.eq.defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, var)

  ; Define style
  ; new: if data array has attributes annots, colors, markers, use those
  ;      instead of function project_style
  if (isatt(data, "annots")) then
    annots = data@annots
  else
    annots = project_style(inlist, diag_script_info, "annots")
  end if
  if (isatt(data, "colors")) then
    colors = data@colors
  else
    colors  = project_style(inlist, diag_script_info, "colors")
  end if
  if (isatt(data, "markers")) then
    markers = data@markers
  else
    markers = project_style(inlist, diag_script_info, "markers")
  end if

  ; Create a 2D scatter plot
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@xyMarkLineMode = "Markers"
  res@xyMarkerSizeF = 0.02
  res@xyMarkerThicknessF = 3.
  res@xyMonoMarkerColor = False
  res@xyMarkerColors = colors
  res@xyMarkers = markers
  res@tiMainFontHeightF = 0.025
  res@tiYAxisFontHeightF = 0.02
  res@tiXAxisFontHeightF = 0.02
  res@tiXAxisString = data&statistic(0)
  res@tiYAxisString = data&statistic(1)
  if (isatt(data, "units")) then
    res@tiXAxisString = res@tiXAxisString + \
      " [" + format_units(data@units(0)) + "]"
    res@tiYAxisString = res@tiYAxisString + \
      " [" + format_units(data@units(1)) + "]"
  end if

  if (isatt(data, "trYMinF")) then
    res@trYMinF = data@trYMinF        ; min value on y-axis
  end if
  if (isatt(data, "trYMaxF")) then
    res@trYMaxF = data@trYMaxF        ; max value on y-axis
  end if
  if (isatt(data, "trXMinF")) then
    res@trXMinF = data@trXMinF        ; min value on x-axis
  end if
  if (isatt(data, "trXMaxF")) then
    res@trXMaxF = data@trXMaxF        ; max value on x-axis
  end if
  if (isatt(data, "trXReverse")) then
    res@trXReverse = data@trXReverse  ; flip the x axis
  end if
  if (isatt(data, "trYReverse")) then
    res@trYReverse = data@trYReverse  ; flip the y axis
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  plot = gsn_csm_xy(wks, transpose((/xpts, xpts/)), \
                    transpose((/ypts, ypts/)), res)

  if (data@legend_outside) then

    lg_outfile = wks@legendfile
    styles = True
    styles@colors = colors
    styles@markers = markers
    create_legend_lines(annots, styles, lg_outfile, "markers")

  else

    lgres                    = True
    lgres@lgAutoManage       = False
    lgres@lbTopMarginF       = 0.
    lgres@lbBottomMarginF    = 0.
    lgres@lbLeftMarginF      = 0.
    lgres@lbRightMarginF     = 0.
    lgres@lgItemType         = "Markers"
    lgres@lgMarkerColors     = colors(::-1)
    lgres@lgMarkerIndexes    = markers(::-1)
    lgres@lgPerimOn          = False
    lgres@lgLabelJust        = "CenterLeft"
    lgres@lgLabelFontHeightF = 0.05
    lgres@vpWidthF           = 0.18
    lgres@vpHeightF          = 0.02 * dimsizes(annots)
    gsn_legend_ndc(wks, dimsizes(annots), annots(::-1), 0.8, 0.78, lgres)

  end if

  leave_msg(scriptname, funcname)
  return(plot)

end
