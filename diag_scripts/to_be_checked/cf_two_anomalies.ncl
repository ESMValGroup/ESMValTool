;*********************************************************************
; cf_two_anomalies.ncl
;*********************************************************************
; NCL script
; cf_two_anomalies.ncl
; Martin Evaldsson Augusti 2012
;*********************************************************************
; This plot_type is a try to visualise two separate variables in the same graph.
; In this case, the standard deviation of interannual anomalies for two variables
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/write_infobox.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./ncl_code/apply_mask.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

;*********************************************************************

begin
    plot_type = "cf_two_anomalies"
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET10 - " + plot_type + " (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields  = (/"T2Ms"/)
    ind_f        = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal: cf_two_anomalies.ncl: can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    supported_info_transforms = (/"masks_pr", "latrange", "lonrange", "season"/)
    if (isatt(info, "season")) then
        season = info@season
    else
        season = "ANN"
    end if

    if (isatt(info, "cf_latrange")) then
        latrange = info@cf_latrange
    else
        latrange = (/-90, 90/)
    end if

    if (isatt(info, "cf_lonrange")) then
        lonrange = info@cf_lonrange
    else
        lonrange = (/-360, 360/)
    end if
    wkdir     = plot_dir + "/" + plot_type + "/"
    file_type = "PS"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;; *************************************************************
    ;; find size of array needed to store data and fetch time axis
    ;;------------------------------------------------------------------
    modeldata       = extract_data(0, data, -1, 0, 0)
    modeldata_time  = cd_calendar(modeldata&time(0::12), 0)
    modeldata_years = modeldata_time(:, 0)
    delete(modeldata)
    delete(modeldata_time)

    anom_plots = new((/dim_MOD, dimsizes(modeldata_years), 2/), "float")
    ;; = =================================================================

    do imod = 0, dim_MOD - 1
;       if ((stringtoint(models@end_date(imod))-stringtoint(models@start_date(imod))) .gt. 0) then
;          yrs_ave1 = models@start_date(imod) + "-" + models@end_date(imod)
;          lgLabels(imod) = case1 + " ("+yrs_ave1+")"
;       else
;          lgLabels(imod) = case1
;       end if
        A0 = extract_data(imod, data, -1, 0, 0)

        ;;cosine of latitude for averaging (i.e., weights)
        rad        = 4.0 * atan(1.0) / 180.0
        coslat     = cos(A0&lat * rad)
        coslat!0   = "lat"
        coslat&lat = A0&lat
        wgt_lat    = coslat({latrange(0):latrange(1)})

        if (isatt(info, "masks_pr") \
            .and. inlist("masks_pr", supported_info_transforms)) then
            ;applyMask(A0(:,:,:, 0), info@masks) ; no mask for mhi
            applyMask(A0(:, :, :, 1), info@masks_pr)
        end if

        var1                      = A0(:, :, :, 0)
        seasonal_anomaly          = spatial_and_seasonal_average(var1({lat|latrange(0):latrange(1)}, \
                                          {lon|lonrange(0):lonrange(1)}, time|:), wgt_lat, season)
        std_anom                  = dim_standardize_Wrap(seasonal_anomaly, 0)
        anom_plots(imod, :, 0)     = std_anom
        anom_plots@standard_name1 = A0@standard_name1

        var2                      = A0(:, :, :, 1)
        seasonal_anomaly          = spatial_and_seasonal_average(var2({lat|latrange(0):latrange(1)}, \
                                          {lon|lonrange(0):lonrange(1)}, time|:), wgt_lat, season)
        std_anom                  = dim_standardize_Wrap(seasonal_anomaly, 0)
        anom_plots(imod, :, 1)     = std_anom
        anom_plots@standard_name2 = A0@standard_name2

        delete(seasonal_anomaly)
        delete(var1)
        delete(var2)
        delete(A0)
        delete(coslat)
        delete(wgt_lat)
    end do  ; imod

    sres                       = True
    sres@vpWidthF              = 0.7
    sres@vpHeightF             = 0.5
    sres@gsnYRefLine           = 0.         ; reference line
    sres@vpXF                  = .15
    sres@trXMinF               = 0.4
    sres@trYMinF               = -1.8
    sres@trYMaxF               = 1.8
    sres@gsnDraw               = True
    sres@gsnFrame              = False
    sres@gsnXYBarChart         = True
    sres@gsnXYBarChartBarWidth = 0.15       ; change bar widths

    sres@tmXBMode              = "Explicit"  ; explicit labels
    no_years                   = dimsizes(modeldata_years)
    sres@trXMaxF               = no_years + 0.6
    sres@tmXBValues            = ispan(1, no_years - 1, 1)

    ;; round off and convert to integer type
    stride_year = toint(floor(no_years/6.))
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

    sres@tiMainFontHeightF   = 0.025
    sres@tiMainFont          = "helvetica"
    sres@gsnRightString      = ""
    sres@tiYAxisString       = "Std. deviation"
    sres@tiYAxisFontHeightF  = 0.0195
    useColors                = (/"red", "lightblue"/)
    sres@gsnXYBarChartColors = useColors(0)

    ;; define plot labels
    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.3           ; labelbar width
    lbres@vpHeightF          = 0.1           ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.36          ; puts space between color boxes
    lbres@lbFillColors       = (/"red", "lightblue"/)
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
    lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    lbres@lbPerimOn          = False
    lbres@lgPerimColor       = "white"
    xaxis = fspan(1 - 0.075, no_years - 0.075, no_years)
    labels = new(2, string)

    do imod = 0, dim_MOD - 1
        path_filename = wkdir + plot_type + "_" + models@name(imod) + "_" + variable
        wks           = gsn_open_wks(file_type, path_filename)

        if (isatt(info, "title")) then
            sres@tiMainString = info@title
        else
            sres@tiMainString = "interannual anomalies for " + models@name(imod)
        end if

        ;; draw each bar separately (plot1 and plot2 respectively)
        bar_displacement         = 0.0
        sres@gsnXYBarChartColors = (/useColors(0)/)
        plot1                    = gsn_csm_xy(wks, xaxis + bar_displacement, anom_plots(imod, :, 0), sres)

        bar_displacement         = bar_displacement + 0.075
        sres@gsnXYBarChartColors = (/useColors(1)/)
        plot2                    = gsn_csm_xy(wks, xaxis + bar_displacement, anom_plots(imod, :, 1), sres)

        labels(0) = anom_plots@standard_name1
        labels(1) = anom_plots@standard_name2
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

    end do  ; imod

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
