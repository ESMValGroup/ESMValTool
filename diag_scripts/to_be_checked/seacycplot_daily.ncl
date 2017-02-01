;*********************************************************************
; seacycplot.ncl
;*********************************************************************
; NCL script
; plot_seas_cycle.ncl (set 10)
; Mark Stevens April 2002
;*********************************************************************
; The original NCL script, plot_seas_cycle.ncl, was written to make...
; annual cycle line plots of global means (set10)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; This CCMVal script has been furthered developed as a part of the
; ESMValTool, an ESM Evaluation Tool within the EMBRACE project. Package
; leaders are Klaus Wyser (SMHI) and Veronika Eyring (DLR).
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/write_infobox.ncl"
load "./ncl_code/apply_mask.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;*********************************************************************

begin
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET10 - seacycplot_daily (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields  = (/"T2Ds"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;supported_info_transforms = (/"masks"/); , "latrange", "lonrange"/)
    supported_info_transforms = (/"masks", "scd_latrange", "scd_lonrange"/)

    wkdir      = plot_dir + "/seacycplot_daily/"
    file_type  = "PS"
    color_type = "COLOR"

    if(isvar("latrange")) then
        delete(latrange)
    end if
    if (isatt(info, "scd_latrange") \
           .and. inlist("scd_latrange", supported_info_transforms)) then
        latrange = info@scd_latrange
    else
        latrange = (/-90, 90/)
    end if

    if(isvar("lonrange")) then
        delete(lonrange)
    end if
    if (isatt(info, "scd_lonrange") \
           .and. inlist("scd_lonrange", supported_info_transforms)) then
        lonrange = info@scd_lonrange
    else
        lonrange = (/-360, 360/)
    end if

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;;******************************************************************
    ;;--------------------- Open files for plots -----------------------
    if (color_type .eq. "COLOR") then
        wks = gsn_open_wks(file_type, wkdir + "seacycplot_" + variable + "_c")
    else
        wks = gsn_open_wks(file_type, wkdir + "seacycplot_" + variable)
    end if

    ;;******************************************************************
    ;;------------------------------------------------------------------
    restick                = True
    restick@ttmFormat      = "%c"
    restick@ttmAxis        = "XB"
    restick@ttmMajorStride = 33
    time                   = ispan(1, 366, 1)
    time!0                 = "time"
    time&time              = time
    time@bounds            = "time_bnds"
    time@units             = "days since 1998-1-1"
    time@calendar          = "gregorian"
    time@axis              = "T"
    time@long_name         = "time"
    time@standard_name     = "time"

    res = True
    time_axis_labels(time, res, restick)

    gblmean   = new((/dim_MOD, 366/), float)
    gblmean!0 = "line"
    gblmean!1 = "month"
    lgLabels  = new((/dim_MOD/), string)

    do imod = 0, dim_MOD - 1
        ;; CASE 1  A(time, lat, lon)
        case1 = models@name(imod) + " " + models@case_name(imod)
        if ((stringtoint(models@end_date(imod)) \
                  - stringtoint(models@start_date(imod))) .gt. 0) then
            yrs_ave1 = models@start_date(imod) + "-" + models@end_date(imod)
            lgLabels(imod) = case1 + " (" + yrs_ave1 + ")"
        else
            lgLabels(imod) = case1
        end if

        A0 = extract_data(imod, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        if (isatt(info, "masks") .and. inlist("masks", supported_info_transforms)) then
            applyMask(A0, info@masks)
        end if

        ;; compute daily climatology
        if isStrSubset(field_number, "T") then
            tmp         = extract_time_range(A0, models@start_date(imod), \
                                                 models@end_date(imod))
            gregdate    = cd_calendar(tmp&time, -2)
            yyyyddd     = yyyymmdd_to_yyyyddd(gregdate)
            pr_runave   = runave_n_Wrap(tmp, 30, 0, 0)
            delete(tmp)
            climdayave  = clmDayTLL(pr_runave, yyyyddd)
            delete(pr_runave)
            delete(A0)
            A0 = climdayave
            delete(climdayave)
        end if

        A0&lat@units = "degrees_north"
        lat1 = A0&lat
        nlat1 = dimsizes(lat1)

        ;;****************************
        ;; compute weights for averging
        if (mod(nlat1, 2) .eq. 1) then
            ;; no gw in file
            gw1_tmp = latGauWgt(nlat1 - 1, "lat", "gaussian weights", "none")
            gw1     = new(nlat1, float)
            gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
            gw1(((nlat1/2) + 1):(nlat1 - 1)) = gw1_tmp((nlat1/2):(nlat1 - 2))
            delete(gw1_tmp)
        else
            ;; no gw in file
            gw1 = latGauWgt(nlat1, "lat", "gaussian weights", "none")
        end if
        gw1!0   = "lat"
        gw1&lat = lat1
        ;; ==================================

        ;;***************************************************************
        do m = 0, 365         ; average data
            A = A0(year_day|m, {lat|latrange(0):latrange(1)}, {lon|lonrange(0):lonrange(1)})
            gblmean(imod, m) = wgt_areaave(A, gw1({lat|latrange(0):latrange(1)}), 1., 0)

            delete(A)
        end do
        ;; ==================================

        gblmean(imod, 365) = gblmean(imod, 0)

        delete(A0)
        delete(lat1)
        delete(nlat1)
        delete(gw1)
    end do  ; imod

    ;;******************************************************************
    colors       = CCMVal2_ColorTable(models@name)
    dashs        = CCMVal2_DashTable(models@name)
    LineThick    = new((/dim_MOD/), integer)
    LineThick(:) = 2
    LineMode     = new((/dim_MOD/), string)
    LineMode(:)  = "MarkLines"

    ;;------------------------------------------------------------------
    ;; common plot resources
    res@gsnMaximize         = True
    res@gsnPaperOrientation = "portrait"
    res@tiMainString        = "Annual Cycle Global Mean Climatology"
    res@tiMainFontHeightF   = 0.022
    res@tiMainFont          = "helvetica-bold"
    res@tiYAxisFontHeightF  = 0.015
    if (color_type .eq. "COLOR") then
        res@xyLineColors = colors
    else
        res@xyMonoLineColor = True
        res@xyLineColor     = "black"
    end if
    res@xyLineThicknesses      = LineThick
    res@xyDashPatterns         = dashs
    res@xyMarkLineModes        = LineMode
    res@xyMarker               = 1
    res@xyMarkerSizeF          = 0.01
    res@pmLegendZone           = 7
    res@pmLegendDisplayMode    = "Always"
    res@pmLegendSide           = "Bottom"
    res@pmLegendWidthF         = 0.15
    res@pmLegendOrthogonalPosF = -.018
    res@pmLegendParallelPosF   = 0.5
    if (dim_MOD .le. 4) then
        res@pmLegendHeightF    = 0.08
        res@lgLabelFontHeightF = 0.016
    else
        if (dim_MOD .le. 8) then
            res@pmLegendHeightF    = 0.15
            res@lgLabelFontHeightF = 0.012
        else
            res@pmLegendHeightF    = 0.30
            res@lgLabelFontHeightF = 0.008
        end if
    end if
    res@lgPerimOn            = True
    res@txFontHeightF        = 0.020
    res@tmXBMode             = "Explicit"
    res@tmXBLabelFontHeightF = 0.015
    res@tmXBMinorOn          = False
;    res@tmXBLabels = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", \
;                       "Aug", "Sep", "Oct", "Nov", "Dec", "Jan"/)
;    res@tmXBValues = months
    res@xyExplicitLegendLabels = lgLabels

    ;;-------------------------------------------------------------------------
    res@tiYAxisString = info@units
    res@gsnLeftString = info@long_name

    ;; ********************************
    ;; create text info box and place it
    ;; ********************************
    txres               = True
    txres@txFontHeightF = 0.012
    txres@txFuncCode    = "~"
    txres@txPerimOn     = True
    infobox             = write_infobox(info, supported_info_transforms)
    infobox_id          = gsn_create_text(wks, infobox, txres)
    ;; =================================

    res@gsnDraw  = False                    ; Turn off draw and frame so
    res@gsnFrame = False                    ; we can attach some text.

    plot = gsn_csm_xy(wks, time, gblmean, res)

    amres                  = True
    amres@amZone           = 8
    amres@amJust           = 4     ; CenterCenter
    amres@amSide           = 1     ; Bottom
    amres@amOrthogonalPosF = -0.15  ; Move text up a little
    ann_id = gsn_add_annotation(plot, infobox_id, amres)
    draw(plot)
    frame(wks)

    delete(gblmean)
    delete(res)
    delete(plot)
    delete(lgLabels)
    delete(colors)
    delete(dashs)
    delete(LineThick)
    delete(LineMode)

end
