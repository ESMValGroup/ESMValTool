;*********************************************************************
; pdf.ncl
;*********************************************************************
; NCL script
; plot_seas_cycle.ncl (set 10)
; Mark Stevens April 2002
;
; The CCMVal script has been furthered developed as a part of the
; ESMValTool, an ESM Evaluation Tool within the EMBRACE project. Package
; leaders are Klaus Wyser (SMHI) and Veronika Eyring (DLR).
;
; contact: Martin.Evaldsson@smhi.se
;*********************************************************************
; Plot a probability density function for daily data
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
    plot_type = "pdf"
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET10 - " + plot_type + " (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ds"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    supported_info_transforms = (/"masks", "years"/)

    wkdir = plot_dir + "/" + plot_type + "/"
    file_type = "PS"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;;******************************************************************
    ;;--------------------- Open files for plots -----------------------
    wks = gsn_open_wks(file_type, wkdir + plot_type + "_" + variable)

    ;;******************************************************************
    ;;------------------------------------------------------------------
    res = True

    global_mean = new((/dim_MOD, 366/), float)
    global_mean!0 = "line"
    global_mean!1 = "month"
    global_mean@_FillValue = 1e+20
    lgLabels = new((/dim_MOD/), string)

    do imod = 0, dim_MOD - 1
        case1 = models@name(imod) + " " + models@case_name(imod)
        if ((stringtoint(models@end_date(imod)) \
                          - stringtoint(models@start_date(imod))) .gt. 0) then
            yrs_ave1       = models@start_date(imod) + "-" + models@end_date(imod)
            lgLabels(imod) = case1 + " (" + yrs_ave1 + ")"
        else
            lgLabels(imod) = case1
        end if
        A0 = extract_data(imod, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        if (isatt(info, "masks") .and. \
                 inlist("masks", supported_info_transforms)) then
            applyMask(A0, info@masks)
        end if

        ;; compute daily climatology
        if isStrSubset(field_number, "T") then
            tmp         = extract_time_range(A0, models@start_date(imod), \
                                             models@end_date(imod))
            gregdate    = cd_calendar(tmp&time, -2)
            yyyyddd     = yyyymmdd_to_yyyyddd(gregdate)
            climdayave  = clmDayTLL(tmp, yyyyddd)
            delete(tmp)
            delete(A0)
            A0          = climdayave
            delete(climdayave)
            delete(gregdate)
            delete(yyyyddd)
        end if

        A0&lat@units = "degrees_north"
        lat1 = A0&lat
        nlat1 = dimsizes(lat1)

        ;;****************************
        ;; compute weights for averging
        if (mod(nlat1, 2).eq.1) then
            ;; no gw in file
            gw1_tmp = latGauWgt(nlat1 - 1, "lat", "gaussian weights", "none")
            gw1 = new(nlat1, float)
            gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
            gw1(((nlat1/2) + 1):(nlat1 - 1)) = gw1_tmp((nlat1/2):(nlat1 - 2))
            delete(gw1_tmp)
        else
            ;; no gw in file
            gw1 = latGauWgt(nlat1, "lat", "gaussian weights", "none")
        end if
        gw1!0 = "lat"
        gw1&lat = lat1
        ;; ==================================

        ;;***************************************************************
        do m = 0, 365         ; average data
            A = A0(m, :,:)
            global_mean(imod, m) = wgt_areaave(A, gw1, 1., 0)
            delete(A)
        end do
        ;; ==================================

        global_mean(imod, 365) = global_mean(imod, 0)

        delete(A0)
        delete(lat1)
        delete(nlat1)
        delete(gw1)

        ;; ********************************
        ;; create text info box
        ;; ********************************
        txres               = True
        txres@txFontHeightF = 0.012
        txres@txFuncCode    = "~"
        txres@txPerimOn     = True
        info@years          = models@start_date(imod) + "-" + models@end_date(imod)
        infobox             = write_infobox(info, supported_info_transforms)
        delete(info@years)
        ;; ==================================

        ;; put the text info box into the figure
        gsn_text_ndc(wks, infobox, 0.3, 0.10, txres)

    end do  ; imod

    ;; prepare and plot bins
    numberOfBins = info@number_of_bins
    ap      = pdfx(global_mean(0, :), numberOfBins, False)
    nBin    = ap@nbins          ; retrieve the number of bins
    xx      = new((/dim_MOD, nBin/), typeof(ap))
    yy      = new((/dim_MOD, nBin/), typeof(ap))

    do imod = 0, dim_MOD - 1
        ap = pdfx(global_mean(imod, :), numberOfBins, False)
        xx(imod, :) = ap@bin_center     ; assign appropriate "x" axis values
        yy(imod, :) = (/ ap /)
        delete(ap)
    end do

    ; set plot properties
    res                          = True
    res@xyLineThicknesses        = (/2.0, 2.0, 2.0/)
    res@xyMonoDashPattern        = True              ; all solid
    res@tiYAxisString            = "PDF (%)"
    colors                       = CCMVal2_ColorTable(models@name)
    res@xyLineColors             = colors

    res@gsnXYBarChart            = True              ; Create bar plot
    res@gsnXYBarChartOutlineOnly = True

    res@pmLegendDisplayMode      = "Never"  ; turn off legend
    res@pmLegendSide             = "Top"    ; Change location of
    res@pmLegendParallelPosF     = .85      ; move units right
    res@pmLegendOrthogonalPosF   = -0.15    ; move units down
    res@pmLegendWidthF           = 0.125    ; Change width and
    res@pmLegendHeightF          = 0.10     ; height of legend.
    res@pmLegendZone             = 4

    res@lgPerimOn                = True         ; turn off/on box around
    res@lgLabelFontHeightF       = .010         ; label font height
    res@lgTitleExtentF           = 0.05
    res@lgPerimFill              = "SolidFill"  ; fill the legend
    res@lgPerimFillColor         = "white"      ; with white
    res@tiXAxisString            = info@units
    res@xyExplicitLegendLabels   = lgLabels     ; create explicit labels

    res@tiMainDirection          = "Across"
    res@tiMainString             = "PDF for " + info@long_name  ; title
    res@gsnFrame                 = False       ; Don't advance frame

    plot = gsn_csm_xy(wks, xx, yy, res)      ; plot without legend

    res@tmXBBorderOn      = False            ; turn of axes and replot such that
    res@tmXBMinorPerMajor = 1                ; the legend is plotted on top of
    res@tmXBOn            = False            ; the axes
    res@tmXTBorderOn      = False
    res@tmXTOn            = False
    res@tmYLBorderOn      = False
    res@tmYLMinorPerMajor = 1
    res@tmYLOn            = False
    res@tmYRBorderOn      = False
    res@tmYROn            = False
    res@pmLegendDisplayMode    = "Always"  ; turn on legend
    res@tiMainOn          = False
    res@tiXAxisOn         = False
    res@tiYAxisOn         = False
    res@gsnFrame          = True          ; Advance frame

    plot = gsn_csm_xy(wks, xx, yy, res)  ; re-plot with legend (on top of axes)

    delete(xx)
    delete(yy)
    delete(res)
    delete(global_mean)
    delete(plot)
    delete(lgLabels)
end
