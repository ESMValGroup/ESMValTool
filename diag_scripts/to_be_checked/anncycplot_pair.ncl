;*********************************************************************
; anncycplot_pair.ncl
;*********************************************************************
; NCL script
; plot_ann_cycle.ncl
; Mark Stevens Sept 2001
;*********************************************************************
; The original NCL script, plot_ann_cycle.ncl, was written to make...
; annual cycle contour plots of zonal means (set8)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; This script uses only first two models in the name_list and displays
; individual models and their difference.
;
; This script takes input from the var_att.ncl files
; info attributes required:
; acp_dcnlvl            vector of contour levels covering range of
;                       values for difference field
; acp_cnlvl             vector of contour levels for original field
; rgb_file              color map filename
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/calc_cycle.ncl"
load "./plot_type/uniform_coordinates.ncl"
;*********************************************************************

begin
    print("+++++++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET8 - anncycplot_pair (var: " + variable + ")")
    print("+++++++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "T1Ms", "C2Ms", "C2Ds", "C1Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal: can't use field, " + field_number + " for this plot_type")
        exit
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/anncycplot_pair/"
    compare = "MODEL"
    file_type = "PS"

    ;; The 'models' array has been loaded in main.ncl
    case1 = models@name(0) + " " + models@case_name(0)
    case2 = models@name(1) + " " + models@case_name(1)

    i = 0
    if ((stringtoint(models@end_date(i)) \
         - stringtoint(models@start_date(i))) .gt. 0) then
        yrs_ave1 = models@start_date(i) + "-" + models@end_date(i)
    end if
    i = 1
    if ((stringtoint(models@end_date(i)) \
         - stringtoint(models@start_date(i))) .gt. 0) then
        yrs_ave2 = models@start_date(i) + "-" + models@end_date(i)
    end if

    ;;******************************************************************
    ;;check the input parameters
    if (isatt(info, "acp_cnlvl") .and. isatt(info, "acp_dcnlvl")) then
        cntrs = info@acp_cnlvl
        dcntrs = info@acp_dcnlvl
    else
        print("fatal: plot requires variable attributes 'acp_cnlvl' and 'acp_dcnlvl'")
        exit
    end if

    if (isatt(info, "rgb_file")) then
        color_type = "COLOR"
        rgb_file = "./rgb/" + info@rgb_file
    else
        color_type = "MONO"
        rgb_file = ""
    end if

    ;;******************************************************************
    ;; CASE 1    A(time, lat)
    imod = 0
    A0 = extract_data(imod, data, -1, 0, 0)
    uniform_coordinates(A0, field_number)

    ;; convert to T2Ms, C2Ms
    if ((field_number .eq. "T2Ds") .or. (field_number .eq. "C2Ds")) then
        tmp = day_to_month(A0, new_field_no)
        delete(A0)
        A0 = tmp
        delete(tmp)
    end if

    A = calc_cycle(A0)
    A&lat@units = "degrees_north"
    lat1 = A&lat
    nlat1 = dimsizes(lat1)
    delete(A0)

    ;; CASE 2
    imod = 1
    B0 = extract_data(imod, data, -1, 0, 0)
    uniform_coordinates(B0, field_number)

    ;; convert to T2Ms, C2Ms
    if ((field_number .eq. "T2Ds") .or. (field_number .eq. "C2Ds"))then
        tmp = day_to_month(B0, new_field_no)
        delete(B0)
        B0 = tmp
        delete(tmp)
    end if

    B = calc_cycle(B0)
    B&lat@units = "degrees_north"
    lat2 = B&lat
    nlat2 = dimsizes(lat2)
    delete(B0)

    if (all(ismissing(A)) .or. all(ismissing(B))) then
        print("fatal: " + variable + " not present in input files.")
        delete(A)
        delete(B)
        exit
    end if

    ;;******************************************************************
    ;;---------------------- Open files for plots ----------------------
    if (color_type .eq. "COLOR") then
        wks = gsn_open_wks(file_type, wkdir + "anncycplot_pair_" \
                           + variable + "_c")
        cmap = RGBtoCmap(rgb_file)
        gsn_define_colormap(wks, cmap)
    else
        wks = gsn_open_wks(file_type, wkdir + "anncycplot_pair_" + variable)
    end if

    ;;******************************************************************
    ;;------------------------------------------------------------------
    ;; common resources for contour plots
    plot = new(3, graphic)
    res                        = True
    res@gsnDraw                = False
    res@gsnFrame               = False
    res@txFontHeightF          = 0.025
    res@sfXArray               = ispan(0, 12, 1)
    res@tiMainFontHeightF      = 0.03
    res@tmXBMode               = "Explicit"
    res@tmXBValues             = ispan(0, 12, 1)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                       "O", "N", "D", "J"/)
    if (color_type .eq. "COLOR") then
        res@cnFillOn = True
        res@cnLinesOn = False
        res@lbTitleOn = True
        res@lbLabelFontHeightF = 0.018
        res@lbTitleFontHeightF = 0.02
        res@lbBoxMinorExtentF = 0.18
        res@cnLevelSelectionMode = "ExplicitLevels"
    else
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnFillOn = False
        res@cnLinesOn = True
        res@cnLineLabelsOn = True
        res@cnLineLabelInterval = 1
        res@cnLineLabelFontHeightF = 0.02
        res@cnInfoLabelFontHeightF = 0.02
        cnlvls = 15
        dcnlvls = 12
    end if

    ;;******************************************************************
    ;;------------------------------------------------------------------
    ;; case contour plots of time(x) vs. latitude(y)
    min1 = min(A)
    max1 = max(A)
    min2 = min(B)
    max2 = max(B)

    if (color_type .eq. "COLOR") then
        res@cnLevels = cntrs
        res@lbLabelStride = 2
        res@lbOrientation = "Horizontal"
        res@lbTitleString = "MIN = " + sprintf("%6.2f", min1) +  \
                            " MAX = " + sprintf("%6.2f", max1)
    else
        minab = min((/min1, min2/))
        maxab = max((/max1, max2/))
        mnmxint = nice_mnmxintvl(minab, maxab, cnlvls, False)
        res@cnMinLevelValF = mnmxint(0)
        res@cnMaxLevelValF = mnmxint(1)
        res@cnLevelSpacingF = mnmxint(2)
    end if

    res@gsnLeftString = info@long_name
    res@gsnRightString = info@units
    if (isvar("yrs_ave1")) then
        res@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
    else
        res@tiMainString = case1
    end if
    res@sfYArray = lat1
    plot(0) = gsn_csm_lat_time(wks, A(lat|:, time|:), res)
    if (color_type .eq. "MONO") then
        plot(0) = ZeroNegDashLineContour(plot(0))
    else
        delete(res@lbTitleString)
    end if
    delete(res@sfYArray)
    delete(res@tiMainString)

    if (isvar("yrs_ave2")) then
        res@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
    else
        res@tiMainString = case2
    end if
    res@sfYArray = lat2
    if (color_type .eq. "COLOR") then
        res@lbTitleString = "MIN = " + sprintf("%6.2f", min2) +  \
                            " MAX = " + sprintf("%6.2f", max2)
    end if
    plot(1) = gsn_csm_lat_time(wks, B(lat|:, time|:), res)
    if (color_type .eq. "MONO") then
        plot(1) = ZeroNegDashLineContour(plot(1))
    else
        delete(res@cnLevels)
        delete(res@lbTitleString)
    end if

    delete(res@sfYArray)
    delete(min1)
    delete(max1)
    delete(min2)
    delete(max2)

    ;;----------------------------------------------------------------------
    ;; difference plot of time(x) vs. latitude(y)
    ;; check for different number of latitudes and then
    ;; lineary interpolate to the smaller number of latitudes
    if (nlat1 .ne. nlat2) then                 ; lat grids different
        if (nlat1 .gt. nlat2) then
            tmp1 = linint1(lat1, A, False, lat2, 0)   ; a(time, lat)
            C = B                                  ; copy dims, coords
            C = (/tmp1 - B/)                         ; get diff values
            delete(tmp1)
            res@sfYArray = lat2
        else
            tmp2 = linint1(lat2, B, False, lat1, 0)   ; b(time, lat)
            C = A                                  ; copy dims, coords
            C = (/A - tmp2/)                         ; get diff values
            delete(tmp2)
            res@sfYArray = lat1
        end if
    else                                       ; same grid latitudes
        C = A                                    ; copy dims, coords
        C = (/A - B/)                              ; get diff values
        res@sfYArray = lat1
    end if

    mind = min(C)
    maxd = max(C)

    res@tiMainString = case1 + " - " + case2
    if (color_type .eq. "COLOR") then
        res@cnLevels = dcntrs
        res@lbLabelStride = 1
        res@lbOrientation = "Vertical"
        res@lbTitleString = "MIN = " + sprintf("%6.2f", mind) +  \
                            " MAX = " + sprintf("%6.2f", maxd)
    else
        mnmxint = nice_mnmxintvl(mind, maxd, dcnlvls, False)
        res@cnMinLevelValF = mnmxint(0)
        res@cnMaxLevelValF = mnmxint(1)
        res@cnLevelSpacingF = mnmxint(2)
    end if

    plot(2) = gsn_csm_lat_time(wks, C(lat|:, time|:), res)
    if (color_type .eq. "COLOR") then
        delete(res@lbTitleString)
        delete(res@cnLevels)
    else
        plot(2) = ZeroNegDashLineContour(plot(2))
    end if
    delete(mind)
    delete(maxd)

    pan = True
    pan@gsnMaximize = True
    pan@gsnFrame = False
    pan@gsnPaperOrientation = "portrait"
    gsn_panel(wks, plot, (/2, 2/), pan)
    frame(wks)

    delete(A)
    delete(B)
    delete(C)
    delete(lat1)
    delete(lat2)
    delete(nlat1)
    delete(nlat2)
    delete(plot)
    delete(res)
    delete(pan)
    delete(cntrs)
    delete(dcntrs)

end
