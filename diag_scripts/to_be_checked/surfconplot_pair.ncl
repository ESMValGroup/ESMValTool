;*********************************************************************
; surfconplot_pair.ncl
;*********************************************************************
; NCL script
; plot_surfaces_cons.ncl
; Mark Stevens Sept 2001
; set 5
; Rich Neale Jan 2008
; Jan 2002: updated with significance tests
; Mar 2003: isccp fields
; Jan 2008: ERA40 reanalysis
;*********************************************************************
; The original NCL script, plot_surfaces_cons.ncl, was written to make...
; longitude vs latitude 2D-field contour plots (set5)
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
; scp_dcnlvl  vector of contour levels covering range of
;             values for difference field
; scp_cnlvl   vector of contour levels for original field
; rgb_file    color map filename
; seasons     (optional) array of for averaging. (DJF, MAM,
;             JJA, SON or ANN for annual mean)
;             Note: if not set, use default (/"DJF", "JJA", "ANN"/)
; scp_tcnlvl  (optional) vector of contour levels for tropics field
; scp_tdcnlvl (optional) vector of contour levels covering
;             range of values for tropics difference field
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
;*********************************************************************

begin
    print("+++++++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET5 - surfconplot_pair (var: " + variable + ")")
    print("+++++++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "C2Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/surfconplot_pair/"
    compare = "MODEL"
    file_type = "PS"
    tick_marks = False

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

    ;; ******************************************************************
    ;; check the input parameters
    if (isatt(info, "scp_cnlvl") .and. isatt(info, "scp_dcnlvl")) then
        cntrs = info@scp_cnlvl
        dcntrs = info@scp_dcnlvl
    else
        print("fatal:plot requires variable attributes 'scp_cnlvl' and 'scp_dcnlvl'")
        status_exit(1)
    end if

    if (isatt(info, "scp_tcnlvl") .and. isatt(info, "scp_tdcnlvl")) then
        tropVar = True
        Tcntrs = info@scp_tcnlvl
        dTcntrs = info@scp_tdcnlvl
    else
        tropVar = False
    end if

    if (isatt(info, "seasons")) then
        seasons = info@seasons
    else
        seasons = (/"ANN", "DJF", "JJA"/)
    end if

    if (isatt(info, "rgb_file")) then
        color_type = "COLOR"
        rgb_file = "./rgb/" + info@rgb_file
    else
        color_type = "MONO"
        rgb_file = ""
    end if

    ;; ******************************************************************
    ;; CASE 1
    imod = 0
    A0 = extract_data(imod, data, -1, 0, 0)
    uniform_coordinates(A0, field_number)
    ;; convert to T2Ms
    if (field_number .eq. "T2Ds") then
        tmp = day_to_month(A0, "T2Ds")
        delete(A0)
        A0 = tmp
        delete(tmp)
    end if

    lat1 = A0&lat
    nlat1 = dimsizes(lat1)
    lon1 = A0&lon
    if (lon1(0) .lt. 0.0) then
        lon1 = lon1 + 180.0
        A0&lon = (/lon1/)
    end if

    ;; gw1 = inptr1->gw
    if (mod(nlat1, 2) .eq. 1) then
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

    ;; CASE 2
    imod = 1
    B0 = extract_data(imod, data, - 1, 0, 0)
    uniform_coordinates(B0, field_number)
    ;; convert to T2Ms
    if (field_number .eq. "T2Ds") then
        tmp = day_to_month(B0, "T2Ds")
        delete(B0)
        B0 = tmp
        delete(tmp)
    end if

    lat2 = B0&lat
    nlat2 = dimsizes(lat2)
    lon2 = B0&lon
    if (lon2(0) .lt. 0.0) then
        lon2 = lon2 + 180.0
        B0&lon = (/lon2/)
    end if

    ;; gw2 = inptr2->gw
    if (mod(nlat2, 2) .eq. 1) then
        ;; no gw in file
        gw2_tmp = latGauWgt(nlat2 - 1, "lat", "gaussian weights", "none")
        gw2 = new(nlat2, float)
        gw2(0:(nlat2/2)) = gw2_tmp(0:(nlat2/2))
        gw2(((nlat2/2) + 1):(nlat2 - 1)) = gw2_tmp((nlat2/2):(nlat2 - 2))
        delete(gw2_tmp)
    else
        ;; no gw in file
        gw2 = latGauWgt(nlat2, "lat", "gaussian weights", "none")
    end if
    gw2!0 = "lat"
    gw2&lat = lat2

    ;; lat/lon limis for global tropics plots
    minlat = -38
    maxlat = 38
    minlon = -70
    maxlon = 290

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;; loop through seasons
    do s = 0, dimsizes(seasons) - 1
        season = seasons(s)

        A = calc_SEASON(A0, season)
        A&lat@units = "degrees_north"
        A&lon@units = "degrees_east"

        B = calc_SEASON(B0, season)
        B&lat@units = "degrees_north"
        B&lon@units = "degrees_east"

        ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (all(ismissing(A)) .or. all(ismissing(B))) then
            print("surfconplot_pair.ncl: " + variable + " not present in input files.")
            delete(A)
            delete(B)
            continue
        end if
        ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ;; compute area weighted mean of the fields
        gblmean_A = wgt_areaave(A, gw1, 1.0, 0)
        gblmean_B = wgt_areaave(B, gw2, 1.0, 0)

        ;; compute RMSE of the two fields after interpolating to the
        ;; same grid if needed, and the mean of the difference of the
        ;; interpolated fields
        if (nlat1 .eq. nlat2) then            ; same grids
            gblrmse = wgt_arearmse(A, B, gw1, 1.0, 0)
            troprmse = wgt_arearmse(A({minlat:maxlat}, :), B({minlat:maxlat}, :), \
            gw1({minlat:maxlat}), 1.0, 0)
            C = A
            C = (/A - B/)
            gblmean_C = wgt_areaave(C, gw1, 1.0, 0)
            tropmean_C = wgt_areaave(C({minlat:maxlat}, :), gw1({minlat:maxlat}), 1.0, 0)
        else                                   ; different grids
            if (nlat1 .gt. nlat2) then           ; interp A to nlat2
                tmp = B
                tmp = linint2(lon1, lat1, A, True, lon2, lat2, 0)
                gblrmse = wgt_arearmse(tmp, B, gw2, 1.0, 0)
                troprmse = wgt_arearmse(tmp({minlat:maxlat}, :), B({minlat:maxlat}, :), \
                gw2({minlat:maxlat}), 1.0, 0)
                C = B
                C = (/tmp - B/)
                gblmean_C = wgt_areaave(C, gw2, 1.0, 0)
                tropmean_C = wgt_areaave(C({minlat:maxlat}, :), gw2({minlat:maxlat}), 1.0, 0)
            else                                ; interp B to nlat1
                tmp = A
                tmp = linint2(lon2, lat2, B, True, lon1, lat1, 0)
                gblrmse = wgt_arearmse(A, tmp, gw1, 1.0, 0)
                troprmse = wgt_arearmse(A({minlat:maxlat}, :), tmp({minlat:maxlat}, :), \
                gw1({minlat:maxlat}), 1.0, 0)
                C = A
                C = (/A - tmp/)
                gblmean_C = wgt_areaave(C, gw1, 1.0, 0)
                tropmean_C = wgt_areaave(C({minlat:maxlat}, :), gw1({minlat:maxlat}), 1.0, 0)
            end if
            delete(tmp)
        end if

        ;; ***************************************************************
        ;; -------------------- Open files for plots ---------------------
        if (color_type .eq. "COLOR") then
            cmap = RGBtoCmap(rgb_file)
            wks = gsn_open_wks(file_type, wkdir + "surfconplot_pair_" \
            + season + "_" + variable + "_c")
            gsn_define_colormap(wks, cmap)
            if (tropVar) then
                wksT = gsn_open_wks(file_type, wkdir + "surfconplot_pair_" \
                + season + "_" + variable + "_TROP_c")
                gsn_define_colormap(wksT, cmap)
            end if
        else     ; MONO
            wks = gsn_open_wks(file_type, wkdir + "surfconplot_pair_" \
            + season + "_" + variable)
            if (tropVar) then
                wksT = gsn_open_wks(file_type, wkdir + "surfconplot_pair_" \
                + season + "_" + variable + "_TROP")
            end if
        end if

        ;; ***************************************************************
        ;; -------------------- GLOBAL TROPICAL PLOTS --------------------
        if (tropVar) then
            plotT = new((/3/), "graphic")     ; tropics plots

            ;; ---------------------------------------------
            ;; common resources for global tropics plots
            trop           = True
            trop@gsnDraw   = False
            trop@gsnFrame  = False
            trop@mpProjection = "Mercator"
            trop@mpShapeMode = "FreeAspect"
            trop@mpLimitMode = "LatLon"
            trop@mpMaxLatF = maxlat
            trop@mpMinLatF = minlat
            trop@mpMaxLonF = maxlon
            trop@mpMinLonF = minlon
            trop@mpCenterLatF = (maxlat - minlat) / 2 + minlat
            trop@mpCenterLonF = (maxlon - minlon) / 2 + minlon
            trop@vpWidthF  = 1.00
            trop@vpHeightF = 0.27
            trop@vpXF      = 0.05
            trop@mpGridAndLimbOn = False
            trop@mpPerimOn = True
            trop@mpPerimDrawOrder = "PostDraw"
            trop@mpGeophysicalLineThicknessF = 1.0
            trop@tiMainFont = "Helvetica"
            trop@txFontHeightF = 0.018
            trop@tiMainFontHeightF = 0.022
            if (color_type .eq. "COLOR") then
                trop@lbLabelBarOn  = True
                trop@lbLabelStride = 1
                trop@lbTitleOn     = True
                trop@lbOrientation = "Horizontal"
                trop@pmLabelBarWidthF = 0.95
                trop@lbLabelFontHeightF = 0.022
                trop@lbTitleFontHeightF = 0.022
                trop@mpFillOn      = False
                trop@cnLevelSelectionMode = "ExplicitLevels"
                trop@cnFillOn      = True
                trop@cnLinesOn     = False
                trop@cnLineLabelsOn = False
                trop@cnMissingValFillColor = "Background"
            else
                trop@cnLevelSelectionMode = "ManualLevels"
                trop@mpFillOn      = True
                trop@cnLinesOn = True
                trop@cnLineLabelsOn = True
                trop@cnLineLabelInterval = 1
                trop@cnLineLabelFontHeightF = .008
                trop@cnInfoLabelFontHeightF = .010
                trop@cnInfoLabelOrthogonalPosF = 0.0
                tropcnlvls = 15
            end if

            ;; common resources for tropics panel plots
            panT = True
            panT@gsnMaximize = True
            panT@gsnFrame = False
            panT@txFontHeightF = 0.020
            panT@txFont = "helvetica-bold"
            panT@gsnPaperOrientation = "portrait"

            tropmean = wgt_areaave(A({minlat:maxlat}, :), gw1({minlat:maxlat}), 1.0, 0)
            min1 = min(A({minlat:maxlat}, :))
            max1 = max(A({minlat:maxlat}, :))
            min2 = min(B({minlat:maxlat}, :))
            max2 = max(B({minlat:maxlat}, :))
            if (color_type .eq. "COLOR") then
                trop@cnLevels = Tcntrs
                trop@lbTitleString = "Min = " + sprintf("%6.2f", min1) \
                                      + " Max = " + sprintf("%6.2f", max1)
            else
                minab = min((/min1, min2/))
                maxab = min((/max1, max2/))
                mnmxint = nice_mnmxintvl(minab, maxab, tropcnlvls, False)
                trop@cnMinLevelValF = mnmxint(0)
                trop@cnMaxLevelValF = mnmxint(1)
                trop@cnLevelSpacingF = mnmxint(2)
            end if
            if (isvar("yrs_ave1")) then
                trop@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
            else
                trop@tiMainString = case1
            end if
            trop@gsnCenterString = "mean= " + sprintf("%6.2f", tropmean)
            trop@gsnLeftString = info@long_name
            trop@gsnRightString = info@units

            plotT(0) = gsn_csm_contour_map_other(wksT, A, trop)
            if (color_type .eq. "MONO") then
                plotT(0) = ZeroNegDashLineContour(plotT(0))
            else
                delete(trop@lbTitleString)
            end if
            delete(trop@tiMainString)
            delete(trop@gsnCenterString)

            tropmean = wgt_areaave(B({minlat:maxlat}, :), gw2({minlat:maxlat}), 1.0, 0)
            if (isvar("yrs_ave2")) then
                trop@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
            else
                trop@tiMainString = case2
            end if
            trop@gsnCenterString = "mean= " + sprintf("%6.2f", tropmean)
            if (color_type .eq. "COLOR") then
                trop@lbTitleString = "Min = " + sprintf("%6.2f", min2) \
                                      + " Max = " + sprintf("%6.2f", max2)
            end if

            plotT(1) = gsn_csm_contour_map_other(wksT, B, trop)
            if (color_type .eq. "MONO") then
                plotT(1) = ZeroNegDashLineContour(plotT(1))
            else
                delete(trop@lbTitleString)
                delete(trop@cnLevels)
            end if
            delete(trop@tiMainString)
            delete(trop@gsnCenterString)
        end if

        ;; ***************************************************************
        ;; ---------------------- GLOBAL CE PLOTS ------------------------
        plot = new((/3/), "graphic")      ; global plots

        ;; -------------------------------------------------------------
        ;; common resources for global contour plots
        res = True
        res@mpCenterLonF = 180
        res@mpGeophysicalLineThicknessF = 1.0
        res@gsnDraw = False           ; Do not draw plot
        res@gsnFrame = False          ; Do not advance frame
        if (tick_marks .eq. "True") then
            res@gsnTickMarksOn = True
        else
            res@gsnTickMarksOn = False
        end if
        res@tiMainFont = "Helvetica"
        res@vpXF = 0.05
        res@vpWidthF = 1.00
        res@vpHeightF = 0.70
        res@mpPerimOn = True
        res@mpPerimDrawOrder = "PostDraw"
        if (color_type .eq. "COLOR") then
            res@mpFillOn = False
            res@lbLabelBarOn = True
            res@lbTitleOn = True
            res@lbBoxMinorExtentF = 0.15
            res@lbLabelStride = 1
            res@txFontHeightF = 0.020
            res@tiMainFontHeightF = 0.025
            res@lbOrientation = "Vertical"
            res@lbLabelFontHeightF = 0.015
            res@lbTitleFontHeightF = 0.018
            res@cnLevelSelectionMode = "ExplicitLevels"
            res@cnFillOn = True
            res@cnLinesOn = False
            res@cnMissingValFillColor = "Background"
            res@cnLineLabelsOn = False
        else
            res@cnLevelSelectionMode = "ManualLevels"
            res@cnLinesOn = True
            res@cnLineLabelsOn = True
            res@cnLineLabelInterval = 1
            res@cnLineLabelFontHeightF = .010
            res@cnInfoLabelOrthogonalPosF = -0.06
            res@cnInfoLabelFontHeightF = .010
            res@txFontHeightF = 0.020
            res@tiMainFontHeightF = 0.025
            gblcnlvls  = 15
            diffcnlvls = 12
        end if

        min1 = min(A)
        max1 = max(A)
        min2 = min(B)
        max2 = max(B)

        if (color_type .eq. "COLOR") then
            res@cnLevels = cntrs
            res@lbTitleString = "Min = " + sprintf("%6.2f", min1) \
                                 + " Max = " + sprintf("%6.2f", max1)
            if (variable .eq. "TREFHT_LAND" .or. variable .eq. "PRECIP_LAND" .or. \
                variable .eq. "PRECT_LAND") then
                res@cnRasterModeOn = True
            else
                res@cnRasterModeOn = False
            end if
        else
            minab = min((/min1, min2/))
            maxab = min((/max1, max2/))
            mnmxint = nice_mnmxintvl(minab, maxab, gblcnlvls, False)
            res@cnMinLevelValF = mnmxint(0)
            res@cnMaxLevelValF = mnmxint(1)
            res@cnLevelSpacingF = mnmxint(2)
        end if

        if (isvar("yrs_ave1")) then
            res@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
        else
            res@tiMainString = case1
        end if
        res@gsnLeftString = info@long_name
        res@gsnRightString = info@units
        res@gsnCenterString = "mean= " + sprintf("%6.2f", gblmean_A)

        plot(0) = gsn_csm_contour_map_ce(wks, A, res)   ; plot model data
        if (color_type .eq. "MONO") then
            plot(0) = ZeroNegDashLineContour(plot(0))
        else
            delete(res@lbTitleString)
        end if
        delete(res@tiMainString)
        delete(res@gsnCenterString)

        if (isvar("yrs_ave2")) then
            res@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
        else
            res@tiMainString = case2
        end if

        if (color_type .eq. "COLOR") then
            res@gsnCenterString = "mean= " + sprintf("%6.2f", gblmean_B)
            res@lbTitleString = "Min = " + sprintf("%6.2f", min2) \
                                 + " Max = " + sprintf("%6.2f", max2)
        else
            res@gsnCenterString = "mean= " + sprintf("%6.2f", gblmean_B)
        end if

        plot(1) = gsn_csm_contour_map_ce(wks, B, res)   ; plot model or obs
        if (color_type .eq. "MONO") then
            plot(1) = ZeroNegDashLineContour(plot(1))
        else
            delete(res@lbTitleString)
            delete(res@cnLevels)
        end if
        delete(res@tiMainString)
        delete(res@gsnCenterString)
        delete(res@gsnLeftString)
        delete(min1)
        delete(max1)
        delete(min2)
        delete(max2)

        ;; --------------------------------------------------------------
        ;; common resources for panel plots
        pan = True
        pan@gsnMaximize = True
        pan@gsnFrame = False
        pan@gsnPaperOrientation = "portrait"
        pan@gsnPanelYWhiteSpacePercent = 2

        ;; ------------------------ DIFFERENCE PLOTS --------------------------
        ;; MAKE DIFFERENCE PLOTS
        ;; tropical plots
        if (tropVar) then
            mind = min(C({minlat:maxlat}, :))
            maxd = max(C({minlat:maxlat}, :))
            trop@tiMainString = case1 + " - " + case2
            trop@gsnLeftString = "mean = " + sprintf("%6.2f", tropmean_C)
            trop@gsnCenterString = "rmse = " + sprintf("%6.2f", troprmse)

            if (color_type .eq. "COLOR") then
                trop@lbTitleString = "Min = " + sprintf("%6.2f", mind) \
                                      + " Max = " + sprintf("%6.2f", maxd)
                trop@cnLevels = dTcntrs
                plotT(2) = gsn_csm_contour_map_other(wksT, C, trop)
                delete(trop@cnLevels)
                delete(trop@lbTitleString)
            else
                mnmxint = nice_mnmxintvl(mind, maxd, diffcnlvls, False)
                trop@cnMinLevelValF = mnmxint(0)
                trop@cnMaxLevelValF = mnmxint(1)
                trop@cnLevelSpacingF = mnmxint(2)
                plotT(2) = gsn_csm_contour_map_other(wksT, C, trop)
                plotT(2) = ZeroNegDashLineContour(plotT(2))
            end if
            delete(trop@tiMainString)
            delete(trop@gsnLeftString)
            delete(trop@gsnCenterString)
            delete(trop@gsnRightString)

            delete(mind)
            delete(maxd)
        end if

        ;; Global CE plots
        mind = min(C)
        maxd = max(C)
        res@tiMainString = case1 + " - " + case2
        res@gsnLeftString = "mean = " + sprintf("%6.2f", gblmean_C)
        res@gsnCenterString = "rmse = " + sprintf("%6.2f", gblrmse)

        if (color_type .eq. "COLOR") then
            res@lbTitleString = "Min = " + sprintf("%6.2f", mind) \
                                 + " Max = " + sprintf("%6.2f", maxd)
            res@cnLevels = dcntrs
            plot(2) = gsn_csm_contour_map_ce(wks, C, res)   ; plot differences
            delete(res@cnLevels)
        else
            mnmxint = nice_mnmxintvl(mind, maxd, diffcnlvls, False)
            res@cnMinLevelValF = mnmxint(0)
            res@cnMaxLevelValF = mnmxint(1)
            res@cnLevelSpacingF = mnmxint(2)
            plot(2) = gsn_csm_contour_map_ce(wks, C, res)
            plot(2) = ZeroNegDashLineContour(plot(2))
        end if
        delete(res@tiMainString)
        delete(res@gsnLeftString)
        delete(res@gsnRightString)
        delete(res@gsnCenterString)

        delete(mind)
        delete(maxd)
        delete(C)
        delete(A)
        delete(B)

        ;; --------------------------------------------------------------
        ;; Panel plots
        if (tropVar) then   ; Tropics
            panT@txString = season
            gsn_panel(wksT, plotT, (/3, 1/), panT)
            frame(wksT)
        end if

        if (color_type .eq. "MONO") then
            pan@txFontHeightF = 0.015
            pan@txString = season
        end if
        gsn_panel(wks, plot, (/3, 1/), pan)

        if (color_type .eq. "COLOR") then
            txres = True
            txres@txFontHeightF = 0.022
            txres@txFont = "helvetica-bold"
            gsn_text_ndc(wks, season, 0.76, 0.97, txres)
        end if
        frame(wks)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do  ; s
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete(A0)
    delete(B0)
    delete(lat1)
    delete(lat2)
    delete(nlat1)
    delete(nlat2)
    delete(lon1)
    delete(lon2)
    delete(gw1)
    delete(gw2)
    delete(plot)
    delete(res)
    delete(pan)
    delete(txres)
    delete(cntrs)
    delete(dcntrs)
    if tropVar then
        delete(plotT)
        delete(trop)
        delete(panT)
        delete(Tcntrs)
        delete(dTcntrs)
    end if
    delete(tropVar)
    delete(seasons)
end
