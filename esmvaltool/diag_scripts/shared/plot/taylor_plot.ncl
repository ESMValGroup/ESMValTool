; #############################################################################
; TAYLOR DIAGRAM ROUTINES
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;
;    function taylor_plot
;    function taylor_diagram
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/set_operators.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"

; #############################################################################
undef("taylor_plot")
function taylor_plot(wks_in[1],
                     source,
                     varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype:
;    source = (ndiag, nmod, statistic)
;    source(:,:,0) = normalized standard deviation
;    source(:,:,1) = correlation
;    source!0 = diagnostics
;    source!1 = models
;    source!2 = statistic
;
; Return value:
;    A graphic variable.
;
; Caveats
;    The taylor plot is drawn in two different ways, depending on ndiag:
;      ndiag = 1: models are drawn with different colors and markers, a
;                 separate legend file is created.
;      ndiag > 1: variables are drawn with different colors, models are
;                 marked with a numbers. Two legends are drawn within the
;                 plot frame: one for the variables (markers) and one for the
;                 models (numbers).
;
; Modification history
;    20150505-righi_mattia: written based on the original NCL code.
;
local funcname, scriptname
begin

  funcname = "taylor_plot"
  scriptname = "diag_scripts/shared/plot/taylor_plot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if(typeof(source) .eq. "string") then
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
  wks = get_wks(wks_in, DIAG_SCRIPT, varname)

  ; Set main quantities
  RATIO      = data(:, :, 0)
  CC         = data(:, :, 1)
  dimR = dimsizes(RATIO)
  ndiags = dimR(0)  ; # of diagnostics
  nmodels = dimR(1)  ; # of models

  ; Select plot-type: single variable or multiple variable
  if (ndiags.eq.1) then
    singlediag = True
  else
    singlediag = False
  end if

  ; Set default colors
  default_colors = (/"red", "blue", "green", "darkorange", "cyan", \
                     "magenta", "black"/)
  default_markers = (/16, 4, 6, 7, 5, 2, 3, 12/)

  ; Default resources
  res                 = True
  res@gsnDraw         = False
  res@gsnFrame        = False
  res@txFontHeightF   = 0.018

  ; One diagnostics, single legend for models using markers and colors
  if (singlediag) then
    if (nmodels.gt.dimsizes(default_colors) * dimsizes(default_markers)) then
      error_msg("f", "taylor_plot.ncl", "taylor_plot", \
                "too many input models (" + nmodels + "), " + \
                "can't draw more than " + dimsizes(default_colors) * \
                dimsizes(default_markers))
    end if
    markers = new(nmodels, integer)
    colors = new(nmodels, string)
    ii = 0
    do idm = 0, dimsizes(default_markers) - 1
      do idc = 0, dimsizes(default_colors) - 1
        if (ii.lt.nmodels) then
          markers(ii) = default_markers(idm)
          colors(ii) = default_colors(idc)
        end if
        ii = ii + 1
      end do
    end do

  ; Multiple diagnostics, one legend for model (with numbers) and one
  ; legend for variable (with colors)
  else
    if (ndiags.gt.dimsizes(default_colors)) then
      error_msg("f", "taylor_plot.ncl", "taylor_plot", \
                "too many input variables (" + ndiags + "), " + \
                "can't draw more than " + dimsizes(default_colors))
    end if
    markers = new(ndiags, integer)
    markers = 16
    colors = default_colors(0:ndiags - 1)
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Set coordinates for plotting
  X = new((/ndiags, nmodels/), typeof(RATIO))
  Y = new((/ndiags, nmodels/), typeof(RATIO))

  do nc = 0, ndiags - 1
    angle = acos(CC(nc, :))
    X(nc, :) = RATIO(nc, :) * cos(angle)
    Y(nc, :) = RATIO(nc, :) * sin(angle)
  end do

  xyMin = 0.
  xyOne = 1.00
  xyMax = 1.65
  xyMax_Panel = xyMax + 0.10

  ; Create base plot
  resXY = True
  resXY@gsnDraw = False
  resXY@gsnFrame = False
  resXY@tmYLBorderOn = False
  resXY@tmXBBorderOn = False
  resXY@tiYAxisString = "Normalized standard deviation"
  resXY@tiYAxisFontHeightF = res@txFontHeightF

  ; Set X-axis labels
  resXY@tmXBMode = "Explicit"
  resXY@tmXBValues = (/0.0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.5/)
  resXY@tmXBLabels = (/"", "0.25", "0.50", "0.75", "REF", "1.25", "1.50"/)
  resXY@tmXBMajorLengthF = 0.015
  resXY@tmXBLabelFontHeightF = res@txFontHeightF
  resXY@tmXBMinorOn = False
  resXY@trXMaxF = xyMax_Panel

  ; Set Y-axis labels
  resXY@tmYLMode = "Manual"
  resXY@tmYLMinorOn = False
  resXY@tmYLMajorLengthF  = resXY@tmXBMajorLengthF
  resXY@tmYLLabelFontHeightF = res@txFontHeightF
  resXY@tmYLMode = "Explicit"
  resXY@tmYLValues = (/0.0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.5/)
  resXY@tmYLLabels = (/"0.00", "0.25", "0.50", "0.75", "1.00", \
                       "1.25", "1.50"/)
  resXY@trYMaxF = xyMax_Panel
  resXY@tmYRBorderOn = False
  resXY@tmYROn = False
  resXY@tmXTBorderOn = False
  resXY@tmXTOn = False

  ; Set correlation axis
  resXY@xyLineThicknessF = 2.
  npts = 200
  xx = fspan(xyMin, xyMax, npts)
  yy = sqrt(xyMax ^ 2 - xx ^ 2)
  sLabels = (/"0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", \
              "0.7", "0.8", "0.9", "0.95", "0.99", "1.0"/)
  cLabels = stringtofloat(sLabels)
  rad = 4. * atan(1.0) / 180.
  angC = acos(cLabels) / rad

  ; Set title (only if ndiag = 1)
  if (singlediag) then
    resXY@tiMainFont = 22
    resXY@tiMainString = data&diagnostics(0)
  end if

  ; Draw plot
  taylor  = gsn_xy(wks, xx, yy, resXY)

  ; Draw X and Y axes
  resL  = True
  resL@gsLineThicknessF  = resXY@xyLineThicknessF
  resL@gsLineDashPattern = 0
  dum0 = gsn_add_polyline(wks, taylor, (/0, 0/), (/0, xyMax/), resL)
  dum1 = gsn_add_polyline(wks, taylor, (/0, xyMax/), (/0, 0/), resL)
  taylor@$unique_string("dum")$ = dum0
  taylor@$unique_string("dum")$ = dum1

  ; Draw stddev ratio lines
  xx = fspan(xyMin, xyOne, npts)
  yy = sqrt(xyOne - xx ^ 2)
  resL@gsLineDashPattern = 1
  resL@gsLineColor = "LightGray"
  dum2 = gsn_add_polyline(wks, taylor, xx, yy, resL)  ; REF
  taylor@$unique_string("dum")$ = dum2

  resL@gsLineThicknessF  = 1
  StnRad = (/0.5, 1.5/)
  nStnRad = dimsizes(StnRad)
  dum3  = new(nStnRad, graphic)
  do ii = 0, nStnRad - 1
    rr = StnRad(ii)
    xx = fspan(xyMin, rr, npts)
    yy = sqrt(rr ^ 2 - xx ^ 2)
    dum3(ii) = gsn_add_polyline(wks, taylor, xx, yy, resL)
  end do
  taylor@$unique_string("dum")$ = dum3
  delete(xx)
  delete(yy)

  ; Get style info from taylor, use for correlation axis
  getvalues taylor
    "tmYLLabelFont"        : tmYLLabelFont
    "tmYLLabelFontHeightF" : tmYLLabelFontHeightF
  end getvalues

  ; Draw correlation labels and ticks
  radC = xyMax
  xC = radC * cos(angC * rad)
  yC = radC * sin(angC * rad)
  xC = xC + 0.020 * cos(rad * angC)
  yC = yC + 0.060 * sin(rad * angC)
  resT = True
  resT@txFontHeightF = res@txFontHeightF
  resT@tmYLLabelFont = tmYLLabelFont
  resT@txAngleF = -45.
  dum4 = gsn_add_text(wks, taylor, "Correlation", 1.30, 1.30, resT)
  taylor@$unique_string("dum")$ = dum4
  resT@txAngleF = 0.0
  resT@txFontHeightF = res@txFontHeightF * 0.50

  resP = True
  resP@gsLineThicknessF = 2.
  resT@txJust = "CenterLeft"
  resT@txFontHeightF = res@txFontHeightF
  tmEnd = 0.975
  radTM = xyMax * tmEnd
  xTM = new(2, "float")
  yTM = new(2, "float")

  dum5 = new(dimsizes(sLabels), graphic)
  dum6 = dum5

  do ii = 0, dimsizes(sLabels) - 1
    resT@txAngleF = angC(ii)
    dum5(ii) = gsn_add_text(wks, taylor, sLabels(ii), xC(ii), yC(ii), resT)
    xTM(0) = xyMax * cos(angC(ii) * rad)
    yTM(0) = xyMax * sin(angC(ii) * rad)
    xTM(1) = radTM * cos(angC(ii) * rad)
    yTM(1) = radTM * sin(angC(ii) * rad)
    dum6(ii) = gsn_add_polyline(wks, taylor, xTM, yTM, resP)
  end do
  taylor@$unique_string("dum")$ = dum5
  taylor@$unique_string("dum")$ = dum6

  mTM = (/0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, \
         0.75, 0.85, 0.91, 0.92, 0.93, 0.94, 0.96, 0.97, 0.98/)
  angmTM = acos(mTM) / rad
  radmTM = xyMax * (1. - (1. - tmEnd) * 0.5)
  dum7 = new(dimsizes(mTM), graphic)
  do ii = 0, dimsizes(mTM) - 1
    xTM(0) = xyMax * cos(angmTM(ii) * rad)
    yTM(0) = xyMax * sin(angmTM(ii) * rad)
    xTM(1) = radmTM * cos(angmTM(ii) * rad)
    yTM(1) = radmTM * sin(angmTM(ii) * rad)
    dum7(ii) = gsn_add_polyline(wks, taylor, xTM, yTM, resP)
  end do
  taylor@$unique_string("dum")$ = dum7

  ; Draw radial lines
  RL = (/0.6, 0.9/)
  angRL = acos(RL) / rad

  resRL = True
  resRL@gsLineDashPattern = 2
  resRL@gsLineThicknessF = 1
  resRL@gsLineColor = "LightGray"

  dum8 = new(dimsizes(angRL), graphic)
  do ii = 0, dimsizes(angRL) - 1
    xRL = xyMax * cos(angRL(ii) * rad)
    yRL = xyMax * sin(angRL(ii) * rad)
    dum8(ii) = gsn_add_polyline(wks, taylor, (/0, xRL/), (/0, yRL/), resRL)
  end do
  taylor@$unique_string("dum")$ = dum8

  ; Draw concentric about REF on X-axis
  resCC                    = True
  resCC@xyLineThicknessF   = 1.
  resCC@xyLineDashPattern  = 2
  resCC@gsLineColor        = "LightGray"
  dx   = 0.25
  ncon = 4
  npts = 100
  ang  = fspan(180, 360, npts) * rad
  dum9 = new(ncon, graphic)
  do ii = 1, ncon
    rr  = ii * dx
    xx  = 1. + rr * cos(ang)
    yy  = fabs(rr * sin(ang))
    if (ii.le.2) then
      dum9(ii - 1) = gsn_add_polyline(wks, taylor, xx, yy, resCC)
    end if
    if (ii.eq.3) then
      n3 = floattointeger(0.77 * npts)
      dum9(ii - 1) = gsn_add_polyline(wks, taylor, xx(0:n3), yy(0:n3), resCC)
    end if
    if (ii.eq.4) then
      n4 = floattointeger(0.61 * npts)
      dum9(ii - 1) = gsn_add_polyline(wks, taylor, xx(0:n4), yy(0:n4), resCC)
    end if
  end do
  delete(ang)
  delete(xx)
  delete(yy)
  taylor@$unique_string("dum")$ = dum9

  ; Add markers
  resM = True
  resM@gsMarkerThicknessF = 1.
  resM@gsMarkerSizeF = 0.015
  if (singlediag) then
    resM@gsMarkerThicknessF = 2.
    dum10 = new(nmodels, graphic)
    do idm = 0, nmodels - 1
      resM@gsMarkerIndex = markers(idm)
      resM@gsMarkerColor = colors(idm)
      dum10(idm) = gsn_add_polymarker(wks, taylor, X(0, idm), \
                                      Y(0, idm), resM)
    end do
    taylor@$unique_string("dum")$ = dum10
  else
    dum10 = new((ndiags * nmodels), graphic)
    resPT = True
    resPT@txJust = "BottomCenter"
    resPT@txFontThicknessF  = 1.2
    resPT@txFontHeightF = res@txFontHeightF
    markerTxYOffset = 0.04
    dum11 = dum10
    do idd = 0, ndiags - 1
      resM@gsMarkerIndex = markers(idd)
      resM@gsMarkerColor = colors(idd)
      resPT@txFontColor = resM@gsMarkerColor
      do idm = 0, nmodels - 1
        dum10(idd * nmodels + idm) = \
          gsn_add_polymarker(wks, taylor, X(idd, idm), Y(idd, idm), resM)
        dum11(idd * nmodels + idm) = \
          gsn_add_text(wks, taylor, idm + 1, X(idd, idm), \
                       Y(idd, idm) + markerTxYOffset, resPT)
      end do
    end do
    taylor@$unique_string("dum")$ = dum10
    taylor@$unique_string("dum")$ = dum11
  end if

  ; Add legend
  resLG  = True
  resLG@lgMarkerColors = colors
  resLG@lgMarkerIndexes = markers
  resLG@lgMarkerSizeF = resM@gsMarkerSizeF
  resLG@lgLabelFontHeightF = 0.15
  resLG@vpWidthF = 0.15
  resLG@lgPerimOn = False
  resLG@lgItemType = "Markers"
  resLG@lgJustification = "TopRight"
  if (singlediag) then
    ; Models (as separate file)
    styles = True
    styles@colors = colors
    styles@markers = markers
    thicks = new(nmodels, float)
    thicks = resM@gsMarkerThicknessF
    styles@thicks = thicks
    sizes = new(nmodels, float)
    sizes = resM@gsMarkerSizeF
    styles@sizes = sizes
    lg_outfile = wks@legendfile
    create_legend_lines(data&models, styles, lg_outfile, "markers")
  else
    ; Diagnostics
    resLG@vpHeightF = 0.030 * ndiags
    resLG@lgItemOrder = ispan(ndiags - 1, 0, 1)
    lbid = gsn_create_legend(wks, ndiags, data&diagnostics, resLG)
    amres = True
    amres@amParallelPosF =  0.40
    amres@amOrthogonalPosF = -0.45
    annoid1 = gsn_add_annotation(taylor, lbid, amres)

    ; Models
    resT@txJust = "CenterLeft"
    resT@txFontHeightF = 0.015
    dum12 = new(nmodels, graphic)
    ypos = 1.50
    do idm = 1, nmodels
      str = idm + " - " + data&models(idm - 1)
      dum12(idm - 1) = gsn_add_text(wks, taylor, str, 0.1, ypos, resT)
      ypos  = ypos - 0.06
    end do
    taylor@$unique_string("dum")$ = dum12
  end if

  draw(taylor)
  frame(wks)
  taylor@outfile = wks@name

  leave_msg(scriptname, funcname)
  return(taylor)

end

; #############################################################################
function taylor_diagram(wks:graphic,
                        legendwks:graphic,
                        RATIO[*][*]:numeric,
                        CC[*][*]:numeric,
                        rOpts:logical)
;
; Arguments
;    wks: workstation.
;    legendwks: workstation for the legend.
;    RATIO: ratio of the standard deviations.
;    CC: correlation.
;    rOpts: graphical resources.
;
; Return value:
;    A graphic variable.
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimR, nCase, nVar, X, Y, nc, angle, xyMin, xyOne, ceil_max_val, xyMax, \
  xyMax_Panel, FontHeightF, rxy, tempo_values, tempo_values1, temp_string, \
  temp1_string, temp_string_y, temp_string_y1, npts, xx, yy, sLabels, \
  cLabels, rad, angC, taylor, rsrRes, dum0, dum1, dum2, nStnRad, rr, radC, \
  xC, yC, txRes, dum4, plRes, tmEnd, radTM, xTM, yTM, dum5, dum6, i, mTM, \
  angmTM, radmTM, dum7, angRL, rlRes, dum8, xRL, yRL, respl, dx, ncon, npts, \
  ang, dum9, dist_n, Markers, Colors, gsRes, ptRes, markerTxYOffset, dum10, \
  dum11, n, lgres, nModel, lbid, lbid2, amres, annoid1, nVar, \
  varLabelsFontHeightF, txres, delta_y, ys, dum12
begin

  dimR = dimsizes(RATIO)
  nCase = dimR(0)  ; # of cases [models]
  nVar = dimR(1)   ; # of variables

  ; x/y coordinates for plotting
  X = new((/nCase, nVar/), typeof(RATIO))
  Y = new((/nCase, nVar/), typeof(RATIO))

  do nc = 0, nCase - 1
    angle = acos(CC(nc, :))
    X(nc, :) = RATIO(nc, :) * cos(angle)
    Y(nc, :) = RATIO(nc, :) * sin(angle)
  end do

  ; fix to nearest x.5
  xyMin = 0.
  xyOne = 1.00
  rmax = ceil(max(RATIO))
  if (rmax.lt.1.5) then
    rmax = 1.5
  end if
  if (rmax.gt.3.0) then
    rmax = 3.0
  end if
  xyMax = rmax + 0.1
  xyMax_Panel = xyMax + 0.10

  if (rOpts .and. isatt(rOpts, "txFontHeightF")) then
    FontHeightF = rOpts@txFontHeightF
  else
    FontHeightF = 0.0175
  end if

  ; ----------------------------------------------------------------
  ; Part 1:
  ; base plot: Based upon request of Mark Stevens
  ; basic x-y and draw the 1.0 observed and the outer curve at 1.65
  ; ----------------------------------------------------------------

  rxy                   = True
  rxy@gsnDraw           = False
  rxy@gsnFrame          = False
  rxy@vpHeightF         = 0.65
  rxy@vpWidthF          = 0.65
  if (rOpts .and. isatt(rOpts, "plotSize")) then
    rxy@vpHeightF         = rOpts@plotSize
    rxy@vpWidthF          = rOpts@plotSize
    rxy@vpXF              = 0.125
    rxy@vpYF              = 0.825
  end if
  rxy@tmYLBorderOn      = False
  rxy@tmXBBorderOn      = False

  rxy@tiYAxisString      = "Standardized Deviations (Normalized)"
  rxy@tiYAxisFontHeightF = FontHeightF

  rxy@tmXBMode          = "Explicit"
  delta = 0.5
  if (rmax.le.1.5) then
    delta = 0.25
  end if
  tempo_values = \
    decimalPlaces(fspan(0.0, rmax, toint(rmax / delta) + 1), 2, True)
  tempo_values1 = tempo_values(ind(tempo_values .lt. xyMax))
  delete(tempo_values)
  rxy@tmXBValues = tempo_values1
  delete(tempo_values1)
  temp_string = tostring(sprintf("%4.2f", rxy@tmXBValues))
  temp_string(0) = "    "
  temp1_string = where(temp_string .eq. "1.00", "REF", temp_string)
  rxy@tmXBLabels = temp1_string
  if (rOpts .and. isatt(rOpts, "OneX")) then
    rxy@tmXBLabels = temp_string
  end if
  delete(temp_string)
  delete(temp1_string)

  rxy@tmXBMajorLengthF  = 0.015
  rxy@tmXBLabelFontHeightF = FontHeightF
  rxy@tmXBMinorOn       = False
  rxy@trXMaxF           = xyMax_Panel
  rxy@tmYLMode          = "Manual"
  rxy@tmYLMinorOn       = False
  rxy@tmYLMajorLengthF  = rxy@tmXBMajorLengthF
  rxy@tmYLLabelFontHeightF = FontHeightF
  rxy@tmYLMode          = "Explicit"
  rxy@tmYLValues = rxy@tmXBValues
  temp_string_y = rxy@tmXBLabels
  temp_string_y(0) = "0.00"
  temp_string_y1 = where(temp_string_y .eq. "REF", "1.00", temp_string_y)
  delete(temp_string_y)
  rxy@tmYLLabels = temp_string_y1
  delete(temp_string_y1)
  rxy@trYMaxF           = xyMax_Panel

  rxy@tmYRBorderOn      = False
  rxy@tmYROn            = False
  rxy@tmXTBorderOn      = False
  rxy@tmXTOn            = False
  rxy@xyDashPatterns    = (/0 /)
  rxy@xyLineThicknesses = (/2./)
  rxy@gsnFrame          = False
  npts    = 100  ; arbitrary
  xx      = fspan(xyMin, xyMax, npts)
  yy      = sqrt(xyMax ^ 2 - xx ^ 2)  ; outer correlation line (xyMax)
  sLabels = (/"0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", \
              "0.7", "0.8", "0.9", "0.95", "0.99", "1.0" /)
  cLabels = stringtofloat(sLabels)
  rad     = 4. * atan(1.0) / 180.
  angC    = acos(cLabels) / rad  ; angles: correlation labels
  if (rOpts .and. isatt(rOpts, "tiMainString")) then
    rxy@tiMainString = rOpts@tiMainString
    if (isatt(rOpts, "tiMainFontHeightF")) then
      rxy@tiMainFontHeightF = rOpts@tiMainFontHeightF
    else
      rxy@tiMainFontHeightF = 0.0225
    end if
  end if

  taylor = gsn_xy(wks, xx, yy, rxy)

  ; --------------------------------------------------------------
  ; draw observational uncertainties given as normalized RMS error
  ; around the observational reference point at (1, 0)
  ; --------------------------------------------------------------

  if (rOpts .and. isatt(rOpts, "rmsobs")) then
    npts = 100                                     ; arbitrary
    rmsobj = new(2, graphic)

    respl0                    = True
    respl0@gsLineThicknessF   = 4.0
    respl0@gsLineDashPattern  = 0
    respl0@gsLineColor        = (/0.50, 0.95, 0.47/)
    respl0@gsFillColor        = (/0.85, 1.00, 0.80/)
    rr = rOpts@rmsobs

    ; uncertainty circle covers whole plotting area
    if (rr.ge.xyMax) then
      ang  = fspan(180, 360, npts) * rad
      xx  = 0.995 * xyMax * cos(ang)
      yy  = fabs(0.995 * xyMax * sin(ang))
      rmsobj(0) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
      rmsobj(1) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
    else  ; uncertainty circle covers only part of the plotting area
      dist = rr + 1.0
      ; uncertainty circle has to be cropped
      ; need to find intersection between uncertainty circule and outer arc
      if (dist .gt. xyMax) then
        xtmp = (xyMax ^ 2 - rr ^ 2 + 1.0) / 2.0
        ytmp = sqrt(xyMax ^ 2 - xtmp ^ 2)
        xtmp = xtmp - 1.0
        if (xtmp.eq.0.0) then
          phitmp = 3.1415927
        else if (xtmp.gt.0.0) then
          phitmp = 2.0 * 3.1415927 - atan(ytmp / xtmp)
        else
          phitmp = 3.1415927 + atan(ytmp / xtmp)
        end if
        end if
        ang  = fspan(3.1415927, phitmp, npts)
        xx1 = 1.0 + rr * cos(ang)
        yy1 = fabs(rr * sin(ang))

        xtmp = xx1(npts - 1)
        ytmp = yy1(npts - 1)
        if (xtmp.eq.0.0) then
          phitmp = 3.1415927
        else if (xtmp.gt.0.0) then
          phitmp = 2.0 * 3.1415927 - atan(ytmp / xtmp)
        else
          phitmp = 3.1415927 + atan(ytmp / xtmp)
        end if
        end if

        ang2 = fspan(phitmp, 0.0, npts)

        xx2 = 0.995 * xyMax * cos(ang2)
        yy2 = 0.995 * fabs(xyMax * sin(ang2))

        if (isvar("xx")) then
          delete(xx)
        end if
        if (isvar("yy")) then
          delete(yy)
        end if
        xx  = array_append_record(xx1, xx2, 0)
        yy  = array_append_record(yy1, yy2, 0)
        rmsobj(0) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
        rmsobj(1) = gsn_add_polyline(wks, taylor, xx1, yy1, respl0)
        delete(xx1)
        delete(xx2)
        delete(yy1)
        delete(yy2)
        delete(ang2)
      else  ; uncertainty circle does not need any adjustments
        ang  = fspan(180, 360, npts) * rad
        xx  = 1.0 + rr * cos(ang)
        yy  = fabs(rr * sin(ang))
        rmsobj(0) = gsn_add_polygon(wks, taylor, xx, yy, respl0)
        rmsobj(1) = gsn_add_polyline(wks, taylor, xx, yy, respl0)
      end if  ; if uncertainty circle has to be cropped

      delete(ang)
      delete(xx)
      delete(yy)

    end if  ; if uncertainty circle covers whole plotting area
  end if  ; if an uncertainty estimate is available

  rsrRes  = True
  rsrRes@gsLineThicknessF  = rxy@xyLineThicknesses(0)
  rsrRes@gsLineDashPattern = 0
  dum0 = gsn_add_polyline(wks, taylor, (/0.,   0. /), (/0., xyMax/), rsrRes)
  dum1 = gsn_add_polyline(wks, taylor, (/0., xyMax/), (/0.,   0. /), rsrRes)

  xx   = fspan(xyMin, xyOne, npts)
  yy   = sqrt(xyOne - xx ^ 2)
  rsrRes@gsLineDashPattern = 1
  rsrRes@gsLineThicknessF = rxy@xyLineThicknesses(0)
  dum2 = gsn_add_polyline(wks, taylor, xx, yy, rsrRes)
  delete(xx)
  delete(yy)

  if (rOpts .and. isatt(rOpts, "stnRad")) then
    rsrRes@gsLineThicknessF = 1
    nStnRad = dimsizes(rOpts@stnRad)

    dum3 = new(nStnRad, graphic)
    do n = 0, nStnRad - 1
      rr = rOpts@stnRad(n)
      if (rr.gt.rmax) then
        continue
      end if
      xx = fspan(xyMin, rr, npts)
      yy = sqrt(rr ^ 2 - xx ^ 2)
      dum3(n) = gsn_add_polyline(wks, taylor, xx, yy, rsrRes)
    end do
    taylor@$unique_string("dum")$ = dum3

    delete(xx)
    delete(yy)
  end if

  getvalues taylor
    "tmYLLabelFont"        : tmYLLabelFont
    "tmYLLabelFontHeightF" : tmYLLabelFontHeightF
  end getvalues

  ; ----------------------------------------------------------------
  ; Part 2:
  ; Correlation labels
  ; ----------------------------------------------------------------
  radC    = 1.02 * xyMax
  xC      = radC * cos(angC * rad)
  yC      = radC * sin(angC * rad)

  txRes               = True
  txRes@txFontHeightF = FontHeightF
  txRes@tmYLLabelFont = tmYLLabelFont
  txRes@txAngleF      = -45.
  if (.not.isatt(rOpts, "drawCorLabel") .or. rOpts@drawCorLabel) then
    xdum = 1.15 * xyMax * cos(45.0 * rad)
    ydum = 1.15 * xyMax * sin(45.0 * rad)
    dum4 = gsn_add_text(wks, taylor, "Correlation", xdum, ydum, txRes)
    taylor@$unique_string("dum")$ = dum4
  end if
  txRes@txAngleF      = 0.0
  txRes@txFontHeightF = FontHeightF * 0.50

  plRes               = True
  plRes@gsLineThicknessF = 2.
  txRes@txFontHeightF = FontHeightF

  tmEnd = 0.975
  radTM = xyMax * tmEnd
  xTM   = new(2, "float")
  yTM   = new(2, "float")

  dum5 = new(dimsizes(sLabels), graphic)
  dum6 = dum5

  do i = 0, dimsizes(sLabels) - 1
    txRes@txAngleF = angC(i)
    if (angC(i).le.45.0) then
      txRes@txJust = "CenterLeft"
    else
      txRes@txJust = "BottomCenter"
    end if
    dum5(i) = \
      gsn_add_text(wks, taylor, sLabels(i), xC(i), yC(i), txRes)
    xTM(0) = xyMax * cos(angC(i) * rad)
    yTM(0) = xyMax * sin(angC(i) * rad)
    xTM(1) = radTM * cos(angC(i) * rad)
    yTM(1) = radTM * sin(angC(i) * rad)
    dum6(i) = gsn_add_polyline(wks, taylor, xTM, yTM, plRes)
  end do
  mTM = (/0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, \
         0.75, 0.85, 0.91, 0.92, 0.93, 0.94, 0.96, 0.97, 0.98/)
  angmTM = acos(mTM) / rad
  radmTM = xyMax * (1. - (1. - tmEnd) * 0.5)

  dum7 = new(dimsizes(mTM), graphic)

  do i = 0, dimsizes(mTM) - 1
    xTM(0) = xyMax * cos(angmTM(i) * rad)
    yTM(0) = xyMax * sin(angmTM(i) * rad)
    xTM(1) = radmTM * cos(angmTM(i) * rad)
    yTM(1) = radmTM * sin(angmTM(i) * rad)
    dum7(i) = gsn_add_polyline(wks, taylor, xTM, yTM, plRes)
  end do

  if (rOpts .and. isatt(rOpts, "ccRays")) then
    angRL = acos(rOpts@ccRays) / rad
    rlRes = True
    rlRes@gsLineDashPattern = 2
    rlRes@gsLineThicknessF = 1
    if (isatt(rOpts, "ccRays_color")) then
      rlRes@gsLineColor    =  rOpts@ccRays_color
    end if

    dum8 = new(dimsizes(angRL), graphic)
    do i = 0, dimsizes(angRL) - 1
      xRL = xyMax * cos(angRL(i) * rad)
      yRL = xyMax * sin(angRL(i) * rad)
      dum8(i) = gsn_add_polyline(wks, taylor, (/0, xRL/), (/0, yRL/), rlRes)
    end do
    taylor@$unique_string("dum")$ = dum8
  end if

  ; ----------------------------------------------------------------
  ; Part 3:
  ; Concentric about 1.0 on XB axis
  ; ----------------------------------------------------------------
  if (rOpts .and. isatt(rOpts, "centerDiffRMS") .and. rOpts@centerDiffRMS) then
    respl                    = True
    respl@gsLineThicknessF   = 2.0
    respl@gsLineColor        = "Black"
    if (isatt(rOpts, "centerDiffRMS_color")) then
      respl@gsLineColor =  rOpts@centerDiffRMS_color
    end if

    dx   = delta
    ncon = 4
    npts = 100  ; arbitrary
    ang  = fspan(180, 360, npts) * rad
    dum9 = new(ncon, graphic)

    do n = 1, ncon
      rr = n * dx            ; radius from 1.0 [OBS] abscissa
      xx = 1. + rr * cos(ang)
      yy = fabs(rr * sin(ang))
      dist_n = dx * n + 1.0
      if (dist_n .gt. max(rxy@tmXBValues)) then
        xtmp = (xyMax ^ 2 - rr ^ 2 + 1.0) / 2.0
        ytmp = sqrt(xyMax ^ 2 - xtmp ^ 2)
        xtmp = xtmp - 1.0
        if (xtmp.eq.0.0) then
          phitmp = 3.1415927
        else if (xtmp.gt.0.0) then
          phitmp = 2.0 * 3.1415927 - atan(ytmp / xtmp)
        else
          phitmp = 3.1415927 + atan(ytmp / xtmp)
        end if
        end if
        ang2  = fspan(3.1415927, phitmp, npts)
        xx2  = 1.0 + rr * cos(ang2)
        yy2  = fabs(rr * sin(ang2))
        dum9(n - 1) = gsn_add_polyline(wks, taylor, xx2, yy2, respl)
        delete(ang2)
        delete(xx2)
        delete(yy2)
      else
        dum9(n - 1) = gsn_add_polyline(wks, taylor, xx, yy, respl)
      end if
      delete(dist_n)

    end do
    delete(ang)
    delete(xx)
    delete(yy)
    taylor@$unique_string("dum")$ = dum9

  end if

  ; draw filled black circle at (1, 0) to mark the reference point
  if (rOpts .and. isatt(rOpts, "gsMarkerSizeF")) then
    gsMarkerSizeF = rOpts@gsMarkerSizeF
  else
    gsMarkerSizeF = 0.0085
  end if
  gsRes = True
  gsRes@gsMarkerIndex = 16
  gsRes@gsMarkerColor = "Black"
  gsRes@gsMarkerSizeF = gsMarkerSizeF
  getvalues taylor
    "vpXF"      : vpx
    "vpYF"      : vpy
    "vpWidthF"  : vpw
    "vpHeightF" : vph
  end getvalues
  dotx = vpx + 1.0 / xyMax_Panel * vpw
  doty = vpy - vph
  gsn_polymarker_ndc(wks, dotx, doty, gsRes)
  dotobj = new(1, graphic)
  dotobj = gsn_add_polymarker(wks, taylor, 1.0, 0.0, gsRes)

  ; ---------------------------------------------------------------
  ; Part 4:
  ; generic resources that will be applied to all users data points
  ; of course, these can be changed
  ; http://www.ncl.ucar.edu/Document/Graphics/Resources/gs.shtml
  ; ---------------------------------------------------------------
  if (rOpts .and. isatt(rOpts, "Markers")) then
    Markers = rOpts@Markers
  else
    Markers = (/4, 6, 8,  0, 9, 12, 7, 2, 11, 16/)
  end if

  if (rOpts .and. isatt(rOpts, "Colors")) then
    Colors  = rOpts@Colors
  else
    Colors  = (/"red", "blue", "green", "cyan", "orange",  \
                "torquoise", "brown", "yellow", "purple", "black"/)
  end if

  if (rOpts .and. isatt(rOpts, "gsMarkerThicknessF")) then
    gsMarkerThicknessF = rOpts@gsMarkerThicknessF
  else
    gsMarkerThicknessF = 2.0
  end if

  if (rOpts .and. isatt(rOpts, "gsMarkerSizeF")) then
    gsMarkerSizeF = rOpts@gsMarkerSizeF
  else
    gsMarkerSizeF = 0.0085
  end if

  gsRes = True
  gsRes@gsMarkerThicknessF = gsMarkerThicknessF
  gsRes@gsMarkerSizeF      = gsMarkerSizeF

  ptRes = True
  ptRes@txJust = "BottomCenter"
  ptRes@txFontThicknessF   = 1.2
  ptRes@txFontHeightF      = 0.0125
  if (rOpts .and. isatt(rOpts, "txFontHeightF")) then
      ptRes@txFontHeightF = rOpts@txFontHeightF
  end if

  markerTxYOffset  = 0.0175
  if (rOpts .and. isatt(rOpts, "markerTxYOffset")) then
    markerTxYOffset = rOpts@markerTxYOffset
  end if

  dum10 = new((nCase*nVar), graphic)
  dum11 = dum10

  do n = 0, nCase - 1
    gsRes@gsMarkerIndex = Markers(n)
    gsRes@gsMarkerColor = Colors(n)
    ptRes@txFontColor   = gsRes@gsMarkerColor
    do i = 0, nVar - 1
      dum10(n * nVar + i) = \
        gsn_add_polymarker(wks, taylor, X(n, i), Y(n, i), gsRes)
      printdl = True
      if (rOpts .and. isatt(rOpts, "printDataLabels")) then
        if (.not.rOpts@printDataLabels) then
          printdl = False
        end if
      end if
      if (printdl) then
        dum11(n * nVar + i) = \
          gsn_add_text(wks, taylor, (i + 1), X(n, i), \
                       Y(n, i) + markerTxYOffset, ptRes)
      end if
    end do
  end do

  ; ---------------------------------------------------------------
  ; Part 5: ; add case legend and variable labels
  ; ---------------------------------------------------------------
  extrafile = False

  if (rOpts .and. isatt(rOpts, "legendExtraFile")) then
    if (rOpts@legendExtraFile) then
      extrafile = True
    end if
  end if

  if (rOpts .and. isatt(rOpts, "caseLabels")) then

    if (isatt(rOpts, "caseLabelsFontHeightF")) then
      caseLabelsFontHeightF = rOpts@caseLabelsFontHeightF
      else
        caseLabelsFontHeightF = 0.05
      end if

      lgres                    = True
      lgres@lgMarkerColors     = Colors
      lgres@lgMarkerIndexes    = Markers
      lgres@lgMarkerSizeF      = gsMarkerSizeF
      lgres@lgItemType         = "Markers"
      lgres@lgLabelFontHeightF = caseLabelsFontHeightF

      nModel                   = dimsizes(rOpts@caseLabels)
      if (isatt(rOpts, "legendWidth")) then
        lgres@vpWidthF       = rOpts@legendWidth
      else
        lgres@vpWidthF       = 0.15
      end if

      if (nModel.gt.20) then
        lgres@vpWidthF = lgres@vpWidthF * 0.5
        lgres@lgMarkerSizeF = lgres@lgMarkerSizeF * 0.75
      end if

      lgres@lgPerimOn          = False
      if (isatt(rOpts, "legendBox")) then
        if (rOpts@legendBox) then
          lgres@lgPerimOn = True
          lgres@lgRightMarginF = 0.1
        end if
      end if

      if (nModel.le.20) then
        if (isatt(rOpts, "reverseOrder")) then
          if (rOpts@reverseOrder) then
            lgres@lgItemOrder = ispan(nModel - 1, 0, 1)
          end if
        end if
        if (isatt(rOpts, "legendHeight")) then
          lgres@vpHeightF      = rOpts@legendHeight
        else
          lgres@vpHeightF      = 0.030 * nCase
        end if
        lbid = gsn_create_legend(legendwks, nModel, rOpts@caseLabels, lgres)
        amres = True
        if (isatt(rOpts, "legendXpos")) then
          amres@amParallelPosF = rOpts@legendXpos
          amres@amJust = "TopLeft"
        else
          amres@amParallelPosF = 0.35
        end if
        if (isatt(rOpts, "legendYpos")) then
          amres@amOrthogonalPosF = rOpts@legendYpos
          amres@amJust = "TopLeft"
        else
          amres@amOrthogonalPosF = -0.35
        end if
        if (extrafile) then
          gsn_legend_ndc(legendwks, nModel, rOpts@caseLabels(0:nModel - 1), \
                         0.78, 0.75, lgres)
        else
          annoid1 = gsn_add_annotation(taylor, lbid, amres)
        end if
      else
        if (isatt(rOpts, "modelsperrow")) then
          modelsperrow = rOpts@modelsperrow
        else
          modelsperrow = (nModel + 1) / 2
        end if
        do lr = 0, nModel / modelsperrow
          ix0 = lr * modelsperrow
          if (ix0.ge.nCase) then
            break
          end if
          ix1 = min((/nCase - 1, ix0 + modelsperrow - 1/))
          dumdelta = ix1 - ix0 + 1

          if (isatt(rOpts, "reverseOrder")) then
            if (rOpts@reverseOrder) then
              if (isatt(lgres, "lgItemOrder")) then
                delete(lgres@lgItemOrder)
              end if
              lgres@lgItemOrder = ispan(dumdelta - 1, 0, 1)
            end if
          end if

          lgres@vpHeightF      = 0.026 * dumdelta

          delete(lgres@lgMarkerColors)
          delete(lgres@lgMarkerIndexes)
          lgres@lgMarkerColors  = Colors(ix0:ix1)
          lgres@lgMarkerIndexes = Markers(ix0:ix1)
          lgres@lgPerimOn = False

          if (extrafile) then
            legend_x = 0.05 + 0.11 * lr
          else
            legend_x = 0.78 + 0.11 * lr
          end if

          gsn_legend_ndc(legendwks, dumdelta, rOpts@caseLabels(ix0:ix1), \
                         legend_x, 0.75, lgres)
        end do
    end if
  end if

  if (rOpts .and. isatt(rOpts, "varLabels")) then
    nVar = dimsizes(rOpts@varLabels)

    if (isatt(rOpts, "varLabelsFontHeightF")) then
      varLabelsFontHeightF = rOpts@varLabelsFontHeightF
    else
      varLabelsFontHeightF = 0.013
    end if

    txres = True
    txres@txFontHeightF = varLabelsFontHeightF
    txres@txJust = "CenterLeft"

    delta_y = 0.06
    if (rOpts .and. isatt(rOpts, "varLabelsYloc")) then
      ys  = rOpts@varLabelsYloc
    else
      ys = max((/nVar * delta_y, 0.30/))
    end if

    do i = 1, nVar
      if (i.eq.1) then
        dum12 = new(nVar, graphic)
      end if

      dum12(i - 1) = \
        gsn_add_text(wks, taylor, i + " - " + rOpts@varLabels(i - 1), \
                     0.125, ys, txres)
        ys = ys - delta_y
      end do

      taylor@$unique_string("dum")$ = dum12
  end if

  taylor@$unique_string("dum")$ = dum0
  taylor@$unique_string("dum")$ = dum1
  taylor@$unique_string("dum")$ = dum2
  taylor@$unique_string("dum")$ = dum5
  taylor@$unique_string("dum")$ = dum6
  taylor@$unique_string("dum")$ = dum7
  taylor@$unique_string("dum")$ = dum10
  taylor@$unique_string("dum")$ = dum11

  if (.not.isatt(rOpts, "taylorDraw") .or. \
      (isatt(rOpts, "taylorDraw") .and. rOpts@taylorDraw)) then
    draw(taylor)
  end if
  if (.not.isatt(rOpts, "taylorFrame") .or. \
     (isatt(rOpts, "taylorFrame") .and. rOpts@taylorFrame)) then
    frame(wks)
  end if

  return(taylor)

end
