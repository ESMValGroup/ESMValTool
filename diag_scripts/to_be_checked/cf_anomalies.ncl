;*********************************************************************
; cf_anomalies.ncl
;*********************************************************************
; NCL script
; cf_anomalies.ncl
; Martin Evaldsson August 2012
;*********************************************************************
; Compare standard deviation of interannual anomalies for a variable
; by plotting annual bars for each model vs the reference
;*********************************************************************
load "./ncl_code/apply_mask.ncl"
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/write_infobox.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;*********************************************************************

begin
    plot_type = "cf_anomalies"
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET10 - " + plot_type + " (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    supported_info_transforms = (/"masks", "season", "cf_latrange", "cf_lonrange"/)

    ;;******************************************************************
    ;;check the input parameters
    if (isatt(info, "season")) then
        season = info@season
    else
        season = "ANN"
    end if

    if (isatt(info, "refModel")) then
        refModel = info@refModel
        ;; The 'models' array has been loaded in main.ncl
        idx_ref = ind(models@name .eq. info@refModel)
        if (all(ismissing(idx_ref))) then
            print(refModel)
            print("fatal:missing reference model - refModel isn't included in model list (namelist)")
            status_exit(1)
        end if
    else
        print("fatal:missing reference model - info@refModel must be specified for " + plot_type)
        status_exit(1)
    end if

    if(isvar("latrange")) then
        delete(latrange)
    end if
    if (isatt(info, "cf_latrange")) then
        latrange = info@cf_latrange
    else
        latrange = (/ -90, 90/)
    end if

    if(isvar("lonrange")) then
        delete(lonrange)
    end if
    if (isatt(info, "cf_lonrange")) then
        lonrange = info@cf_lonrange
    else
        lonrange = (/-360, 360/)
    end if

    wkdir = plot_dir + "/" + plot_type + "/"
    compare = "MODEL"
    file_type = "newps"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;;******************************************************************
    ;;------------------------------------------------------------------
    res = True

    ;; find size of array needed to store data and fetch time axis
    ;;------------------------------------------------------------------
    modeldata       = extract_data(0, data, -1, 0, 0)
    modeldata_time  = cd_calendar(modeldata&time(0::12), 0)
    modeldata_years = modeldata_time(:, 0)
    anom_plots      = new((/dim_MOD, dimsizes(modeldata_years)/), "float")
    delete(modeldata)
    delete(modeldata_time)
    ;;==================================================================

    do imod = 0, dim_MOD - 1
        A0 = extract_data(imod, data, -1, 0, 0)

        ;;cosine of latitude for averaging (i.e., weights)
        rad        = 4.0 * atan(1.0) / 180.0
        coslat     = cos(A0&lat * rad)
        coslat!0   = "lat"
        coslat&lat = A0&lat
        wgt_lat    = coslat({latrange(0):latrange(1)})

        if (isatt(info, "masks") .and. \
            inlist("masks", supported_info_transforms)) then
            applyMask(A0, info@masks)
        end if

        seasonal_anomaly = spatial_and_seasonal_average(A0({lat|latrange(0):latrange(1)}, \
                                  {lon|lonrange(0):lonrange(1)}, time|:), wgt_lat, season)
        std_anom = dim_standardize_Wrap(seasonal_anomaly, 0)
        anom_plots(imod, :) = std_anom

        delete(seasonal_anomaly)
        delete(A0)
        delete(coslat)
        delete(wgt_lat)
    end do  ; imod

    ;; prepare figure attributes
    ;; ---------------------------
    sres                       = True
    sres@vpWidthF              = 0.7
    sres@vpHeightF             = 0.5
    sres@gsnYRefLine           = 0.    ; reference line
    sres@vpXF                  = .15
    sres@trXMinF               = 0.4
    sres@trYMinF               = -1.8
    sres@trYMaxF               = 1.8
    sres@gsnDraw               = True
    sres@gsnFrame              = False
    sres@gsnXYBarChart         = True
    sres@gsnXYBarChartBarWidth = 0.15  ; change bar widths

    ;; compute the years on the x-axis (only include six years to avoid clutter)
    sres@tmXBMode   = "Explicit"        ; explicit labels
    no_years        = dimsizes(modeldata_years)
    sres@trXMaxF    = no_years + 0.6
    sres@tmXBValues = ispan(1, no_years - 1, 1)

    stride_year = toint(floor(no_years/6.))  ; round off and convert to integer
    firstYear   = toint(modeldata_years(0))

    ;; create x-axis labels (years). This has to be done explicitly.
    ;; -----------------------------------
    tmXBLabels_int = ispan(firstYear, firstYear + no_years - 1, 1)
    tmXBLabels_str = tostring(tmXBLabels_int)
    do tmind = 0, dimsizes(modeldata_years) - 1
        if (tmind % stride_year .ne. 0) then
            tmXBLabels_str(tmind) = ""
        end if
    end do
    sres@tmXBLabels = tmXBLabels_str
    ;; -----------------------------------

    sres@tmXBLabelFontHeightF = 0.0195
    sres@tmXTLabelFontHeightF = 0.0195
    sres@tmYLLabelFontHeightF = 0.0195

    sres@tiMainFontHeightF    = 0.025
    sres@tiMainFont           = "helvetica"
    if (isatt(info, "title")) then
        sres@tiMainString = info@title
    else
        sres@tiMainString = "interannual anomalies for " + info@long_name
    end if
    sres@gsnRightString      = ""
    sres@tiYAxisString       = "Std. deviation"
    sres@tiYAxisFontHeightF  = 0.0195
    colors                   = CCMVal2_ColorTable(models@name)
    useColors                = colors
    sres@gsnXYBarChartColors = useColors(0)
    ;;***************************************************************

    ;; prepare legend attributes
    ;; -----------------------------------
    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.3           ; labelbar width
    lbres@vpHeightF          = 0.1           ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
    lbres@lbFillColors       = useColors
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
    lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    lbres@lbPerimOn          = False
    lbres@lgPerimColor       = "white"
    ;;***************************************************************

    xaxis = fspan(1 - 0.075, no_years - 0.075, no_years)
    labels = new(2, string)

    sres@gsnDraw = False                     ; Turn off draw and frame so
    sres@gsnFrame = False                    ; we can attach some text.

    ;; loop although only one one model (+reference) are currently supported
    do imod = 0, dim_MOD - 1
        if (imod .eq. idx_ref) then
            continue
        end if
        path_filename = wkdir + plot_type + "_" + models@name(imod) + "_" + variable
        wks = gsn_open_wks(file_type, path_filename)

        ;; draw each bar separately (plot1 and plot2 respectively)
        bar_displacement         = 0.0
        sres@gsnXYBarChartColors = (/useColors(0)/)
        plot1                    = gsn_csm_xy(wks, xaxis + bar_displacement, \
                                              anom_plots(imod, :), sres)

        bar_displacement         = bar_displacement + 0.075
        sres@gsnXYBarChartColors = (/useColors(1)/)
        plot2                    = gsn_csm_xy(wks, xaxis + bar_displacement, \
                                              anom_plots(idx_ref, :), sres)

        labels(0)   = models@name(imod)
        labels(1)   = models@name(idx_ref)

        ;; draw right labelbar column
        labelbar_id = gsn_create_labelbar(wks, 2, labels, lbres)

        ;; create text info box
        txres               = True
        txres@txFontHeightF = 0.012
        txres@txFuncCode    = "~"
        txres@txPerimOn     = True
        info@years          = models@start_date(imod) + "-" + models@end_date(imod)
        infobox             = write_infobox(info, supported_info_transforms)

        ;; create the text info box
        infobox_id = gsn_create_text(wks, infobox, txres)
    end do   ; imod

    ;; place annotations
    amres                = True
    amres@amZone         = 7
    amres@amJust         = "BottomRight"
    amres@amParallelPosF = 1.0
    am_labelbar_id       = gsn_add_annotation(plot2, labelbar_id, amres)

    amres@amZone           = 8
    amres@amOrthogonalPosF = -0.2
    amres@amParallelPosF   = 0.0
    amres@amJust           = "BottomLeft"
    am_infobox_id          = gsn_add_annotation(plot2, infobox_id, amres)

    draw(wks)
    frame(wks)
end
