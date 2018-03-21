;;#############################################################################
;; TAYLOR DIAGRAM ROUTINES
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;
;;    function taylor_plot
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/logging.ncl"

load "./diag_scripts/shared/set_operators.ncl"

load "./diag_scripts/shared/plot/aux_plotting.ncl"


;##############################################################################

undef("taylor_plot")
function taylor_plot(wks_in[1],
                     source,
                     varname[1]: string)
;;
;; Arguments
;;    wks_in: workstations (graphic object or default will be used).
;;    source: data to be plotted or a NetCDF filename with data.
;;    varname: variable name in the file.
;;
;; Source prototype:
;;    source = (ndiag, nmod, statistic)
;;    source(:,:,0) = normalized standard deviation
;;    source(:,:,1) = correlation
;;    source!0 = diagnostics
;;    source!1 = models
;;    source!2 = statistic
;;
;; Return value:
;;    A graphic variable.
;;
;; Caveats
;;    The taylor plot is drawn in two different ways, depending on ndiag:
;;      ndiag = 1: models are drawn with different colors and markers, a 
;;                 separate legend file is created.
;;      ndiag > 1: variables are drawn with different colors, models are
;;                 marked with a numbers. Two legends are drawn within the
;;                 plot frame: one for the variables (markers) and one for the
;;                 models (numbers).
;;
;; Modification history:
;;    20150505-A_righ_ma: written based on the original NCL code.
;;
local funcname, scriptname
begin

    funcname = "taylor_plot"
    scriptname = "diag_scripts/shared/plot/taylor_plot.ncl"
    enter_msg(scriptname, funcname)

    ;; Get data, either directly or via netCDF file
    if(typeof(source) .eq. "string") then
        data = ncdf_read(source, varname)
    else
        data = source
        copy_VarMeta(source, data)
    end if

    ;; Use -999 as fill value (easier to display as text)
    data@_FillValue = -999.

    ;; Retrieve basic metadata from data
    defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
    if (any(varname .eq.  defaults)) then
        grade = att2var(data, "var")
        diag_script = att2var(data, "diag_script")
    else
        grade = varname
    end if

    ;; Check if a valid wks has been provided, otherwise invoke default
    wks = get_wks(wks_in, diag_script, varname)

    ;; Set main quantities
    RATIO      = data(:, :, 0)
    CC         = data(:, :, 1)
    dimR = dimsizes(RATIO)
    ndiags = dimR(0)  ; # of diagnostics
    nmodels = dimR(1)  ; # of models

    ;; Select plot-type: single variable or multiple variable
    if (ndiags.eq.1) then
        singlediag = True
    else
        singlediag = False
    end if

    ;; Set default colors
    default_colors = (/"red", "blue", "green", "darkorange", "cyan", \
                       "magenta", "black"/)
    default_markers = (/16, 4, 6, 7, 5, 2, 3, 12/)

    ;; Default resources
    res                 = True
    res@gsnDraw         = False
    res@gsnFrame        = False
    res@txFontHeightF   = 0.018

    ;; One diagnostics, single legend for models using markers and colors
    if (singlediag) then
        if (nmodels.gt.dimsizes(default_colors) * \
            dimsizes(default_markers)) then
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
    ;; Multiple diagnostics, one legend for model (with numbers) and one
    ;; legend for variable (with colors)
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

    ;; Override defaults with "res_" attributes of "data"
    res_new = att2var(data, "res_")
    copy_VarMeta(res_new, res)

    ;; Set coordinates for plotting
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

    ;; Create base plot
    resXY = True
    resXY@gsnDraw = False
    resXY@gsnFrame = False
    resXY@tmYLBorderOn = False
    resXY@tmXBBorderOn = False
    resXY@tiYAxisString = "Normalized standard deviation"
    resXY@tiYAxisFontHeightF = res@txFontHeightF

    ;; Set X-axis labels
    resXY@tmXBMode = "Explicit"
    resXY@tmXBValues = (/0.0, 0.25, 0.50, 0.75, 1.00, 1.25, 1.5/)
    resXY@tmXBLabels = (/"", "0.25", "0.50", "0.75", "REF", "1.25", "1.50"/)
    resXY@tmXBMajorLengthF = 0.015
    resXY@tmXBLabelFontHeightF = res@txFontHeightF
    resXY@tmXBMinorOn = False
    resXY@trXMaxF = xyMax_Panel

    ;; Set Y-axis labels
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

    ;; Set correlation axis
    resXY@xyLineThicknessF = 2.
    npts = 200
    xx = fspan(xyMin, xyMax, npts)
    yy = sqrt(xyMax^2 - xx^2)
    sLabels = (/"0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", \
                "0.7", "0.8", "0.9", "0.95", "0.99", "1.0"/)
    cLabels = stringtofloat(sLabels)
    rad = 4. * atan(1.0) / 180.
    angC = acos(cLabels) / rad

    ;; Set title (only if ndiag = 1)
    if (singlediag) then
        resXY@tiMainFont = 22
        resXY@tiMainString = data&diagnostics(0)
    end if

    ;; Draw plot
    taylor  = gsn_xy(wks, xx, yy, resXY)

    ;; Draw X and Y axes
    resL  = True
    resL@gsLineThicknessF  = resXY@xyLineThicknessF
    resL@gsLineDashPattern = 0
    dum0 = gsn_add_polyline(wks, taylor, (/0, 0/), (/0, xyMax/), resL)
    dum1 = gsn_add_polyline(wks, taylor, (/0, xyMax/), (/0, 0/), resL)
    taylor@$unique_string("dum")$ = dum0
    taylor@$unique_string("dum")$ = dum1

    ;; Draw stddev ratio lines
    xx = fspan(xyMin, xyOne, npts)
    yy = sqrt(xyOne - xx^2)
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
        yy = sqrt(rr^2 - xx^2)
        dum3(ii) = gsn_add_polyline(wks, taylor, xx, yy, resL)
    end do
    taylor@$unique_string("dum")$ = dum3
    delete(xx)
    delete(yy)

    ;; Get style info from taylor, use for correlation axis
    getvalues taylor
        "tmYLLabelFont"        : tmYLLabelFont
        "tmYLLabelFontHeightF" : tmYLLabelFontHeightF
    end getvalues

    ;; Draw correlation labels and ticks
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

    ;; Draw radial lines
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

    ;; Draw concentric about REF on X-axis
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
            dum9(ii - 1) = \
                gsn_add_polyline(wks, taylor, xx(0:n3), yy(0:n3), resCC)
        end if
        if (ii.eq.4) then
            n4 = floattointeger(0.61 * npts)
            dum9(ii - 1) = \
                gsn_add_polyline(wks, taylor, xx(0:n4), yy(0:n4), resCC)
        end if
    end do
    delete(ang)
    delete(xx)
    delete(yy)
    taylor@$unique_string("dum")$ = dum9

    ;; Add markers
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
                    gsn_add_polymarker(wks, taylor, \
                                       X(idd, idm), Y(idd, idm), resM)
                dum11(idd * nmodels + idm) = \
                    gsn_add_text(wks, taylor, idm + 1, X(idd, idm), \
                                 Y(idd, idm) + markerTxYOffset, resPT)
            end do
        end do
        taylor@$unique_string("dum")$ = dum10
        taylor@$unique_string("dum")$ = dum11
    end if

    ;; Add legend
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
        ;; Models (as separate file)
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
        ;; Diagnostics
        resLG@vpHeightF = 0.030 * ndiags
        resLG@lgItemOrder = ispan(ndiags - 1, 0, 1)
        lbid = gsn_create_legend(wks, ndiags, data&diagnostics, resLG)
        amres = True
        amres@amParallelPosF =  0.40
        amres@amOrthogonalPosF = -0.45
        annoid1 = gsn_add_annotation(taylor, lbid, amres)

        ;; Models
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
