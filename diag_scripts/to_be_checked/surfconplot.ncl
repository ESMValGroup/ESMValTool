;*********************************************************************
; surfconplot.ncl
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
; The switch 'refModel' stands for the reference model. If present,
; differences with the other models are displayed. If omitted,
; all models are displayed individually.
;
; This script takes input from the var_att.ncl files
; info attributes required:
; scp_dcnlvl            vector of contour levels covering range of
;                       values for difference field
; scp_cnlvl             vector of contour levels for original field
; refModel              reference model name (used with scp_dcnlvl)
; rgb_file              color map filename
; seasons               (optional) array of for averaging. (DJF, MAM,
;                       JJA, SON or ANN for annual mean)
;                       Note: if not set, use default (/"DJF", "JJA", "ANN"/)
; scp_tdcnlvl           (optional) vector of contour levels for tropics field
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
;*********************************************************************

begin
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET5 - surfconplot (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "C2Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal: surfconplot.ncl: can't use field, " \
                      + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/surfconplot/"
    compare = "MODEL"
    file_type = "PS"
    tick_marks = False

    ;; Models array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;; ******************************************************************
    ;; check the input parameters
    if (isatt(info, "refModel")) then
        idx_ref = ind(models@name .eq. info@refModel)
        if (all(ismissing(idx_ref))) then
            print("fatal: surfconplot.ncl: refModel isn't included in model list (namelist)")
            status_exit(1)
        end if
        if (isatt(info, "scp_dcnlvl")) then
            cntrs = info@scp_dcnlvl
        else
            print("fatal: surfconplot.ncl: plot requires variable attribute 'scp_dcnlvl'")
            status_exit(1)
        end if
        p_ref = "ref_"
    else
        idx_ref = -1
        if (isatt(info, "scp_cnlvl")) then
            cntrs = info@scp_cnlvl
        else
            print("fatal: surfconplot.ncl: plot requires variable attribute 'scp_cnlvl'")
            status_exit(1)
        end if
        p_ref = ""
    end if

    if (isatt(info, "scp_tdcnlvl")) then
        tropVar = True
        Tcntrs = info@scp_tdcnlvl
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

    ;; lat/lon limis for global tropics plots
    minlat = -38
    maxlat = 38
    minlon = -70
    maxlon = 290

    ;; ------------------------ GLOBAL TROPICAL PLOTS -------------------
    if (tropVar) then

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
            trop@cnLinesOn     = True
            trop@cnLineLabelsOn = True
            trop@cnLineLabelInterval = 1
            trop@cnLineLabelFontHeightF = .008
            trop@cnInfoLabelFontHeightF = .010
            trop@cnInfoLabelOrthogonalPosF = 0.0
            tropcnlvls = 15
        end if
    end if

    ;; ------------------------- GLOBAL CE PLOTS ------------------------
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
        gblcnlvls = 15
        diffcnlvls = 12
    end if

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;;loop through seasons
    do s = 0, dimsizes(seasons) - 1
        season = seasons(s)

        ;; ***************************************************************
        ;; ------------------ Open files for plots -----------------------
        if (color_type .eq. "COLOR") then
            cmap = RGBtoCmap(rgb_file)
            wks = gsn_open_wks(file_type, wkdir + "surfconplot_" \
                                 + p_ref + season + "_" + variable + "_c")
            gsn_define_colormap(wks, cmap)
            if (tropVar) then
                wksT = gsn_open_wks(file_type, wkdir + "surfconplot_" + p_ref \
                                 + season + "_" + variable + "_TROP_c")
                gsn_define_colormap(wksT, cmap)
            end if
        else     ; MONO
            wks = gsn_open_wks(file_type, wkdir + "surfconplot_" + p_ref \
                                 + season + "_" + variable)
            if (tropVar) then
                wksT = gsn_open_wks(file_type, wkdir + "surfconplot_" + p_ref \
                                 + season + "_" + variable + "_TROP")
            end if
        end if

        if (idx_ref .ne. -1) then
            ;; (b) difference of each model from 'reference'
            plot = new(dim_MOD - 1, "graphic")
            plotT = new(dim_MOD - 1, "graphic")

            ;; CASE 1
            case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
            A0 = extract_data(idx_ref, data, -1, 0, 0)
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

            ;; Handle odd number of latitudes explicitly
            if (mod(nlat1, 2) .eq. 1) then
                gw1_tmp = latGauWgt(nlat1 - 1, "lat", "gaussian weights", "none")
                gw1 = new(nlat1, float)
                gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
                gw1(((nlat1/2) + 1):(nlat1 - 1)) = gw1_tmp((nlat1/2):(nlat1 - 2))
                delete(gw1_tmp)
            else
                ;; An even number of latitudes is ok
                gw1 = latGauWgt(nlat1, "lat", "gaussian weights", "none")
            end if
            gw1!0 = "lat"
            gw1&lat = lat1

            A = calc_SEASON(A0, season)
            A&lat@units = "degrees_north"
            A&lon@units = "degrees_east"
            delete(A0)

            ;; compute area weighted mean of the fields
            gblmean_A = wgt_areaave(A, gw1, 1.0, 0)

            jmod = 0
            do imod = 0, dim_MOD - 1
                if (imod .ne. idx_ref) then
                    ;; CASE 2
                    case2 = models@name(imod) + " " + models@case_name(imod)
                    B0 = extract_data(imod, data, -1, 0, 0)
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

                    B = calc_SEASON(B0, season)
                    B&lat@units = "degrees_north"
                    B&lon@units = "degrees_east"
                    delete(B0)

                    ;; compute area weighted mean of the fields
                    gblmean_B = wgt_areaave(B, gw2, 1.0, 0)

                    ;; compute RMSE of the two fields after interpolating to the
                    ;; same grid if needed, and the mean of the difference of the
                    ;; interpolated fields
                    if (nlat1 .eq. nlat2) then            ; same grids
                        gblrmse = wgt_arearmse(A, B, gw1, 1.0, 0)
                        troprmse = wgt_arearmse(A({minlat:maxlat}, :), \
                                                B({minlat:maxlat}, :), \
                        gw1({minlat:maxlat}), 1.0, 0)
                        C = A
                        C = (/A - B/)
                        gblmean_C = wgt_areaave(C, gw1, 1.0, 0)
                        tropmean_C = wgt_areaave(C({minlat:maxlat}, :), \
                                                    gw1({minlat:maxlat}), 1.0, 0)
                    else                                   ; different grids
                        if (nlat1 .gt. nlat2) then           ; interp A to nlat2
                            tmp = B
                            tmp = linint2(lon1, lat1, A, True, lon2, lat2, 0)
                            gblrmse = wgt_arearmse(tmp, B, gw2, 1.0, 0)
                            troprmse = wgt_arearmse(tmp({minlat:maxlat}, :), \
                                                      B({minlat:maxlat}, :), \
                            gw2({minlat:maxlat}), 1.0, 0)
                            C = B
                            C = (/tmp - B/)
                            gblmean_C = wgt_areaave(C, gw2, 1.0, 0)
                            tropmean_C = wgt_areaave(C({minlat:maxlat}, :), \
                                                        gw2({minlat:maxlat}), 1.0, 0)
                        else                                ; interp B to nlat1
                            tmp = A
                            tmp = linint2(lon2, lat2, B, True, lon1, lat1, 0)
                            gblrmse = wgt_arearmse(A, tmp, gw1, 1.0, 0)
                            troprmse = wgt_arearmse(A({minlat:maxlat}, :), \
                                                    tmp({minlat:maxlat}, :), \
                            gw1({minlat:maxlat}), 1.0, 0)
                            C = A
                            C = (/A - tmp/)
                            gblmean_C = wgt_areaave(C, gw1, 1.0, 0)
                            tropmean_C = wgt_areaave(C({minlat:maxlat}, :), \
                                                        gw1({minlat:maxlat}), 1.0, 0)
                        end if
                        delete(tmp)
                    end if

                    ;; ---------------- GLOBAL TROPICAL PLOTS ---------------
                    if (tropVar) then
                        mind = min(C({minlat:maxlat}, :))
                        maxd = max(C({minlat:maxlat}, :))
                        trop@tiMainString = case1 + " - " + case2
                        trop@gsnLeftString = "mean = " + sprintf("%6.2f", tropmean_C)
                        trop@gsnCenterString = "rmse = " + sprintf("%6.2f", troprmse)
                        trop@gsnRightString = info@units

                        if (color_type .eq. "COLOR") then
                            trop@lbTitleString = "Min = " + sprintf("%6.2f", mind) \
                                                  + " Max = " + sprintf("%6.2f", maxd)
                            trop@cnLevels = Tcntrs
                            plotT(jmod) = gsn_csm_contour_map_other(wksT, C, trop)
                            delete(trop@cnLevels)
                            delete(trop@lbTitleString)
                        else
                            mnmxint = nice_mnmxintvl(mind, maxd, diffcnlvls, False)
                            trop@cnMinLevelValF = mnmxint(0)
                            trop@cnMaxLevelValF = mnmxint(1)
                            trop@cnLevelSpacingF = mnmxint(2)
                            plotT(jmod) = gsn_csm_contour_map_other(wksT, C, trop)
                            plotT(jmod) = ZeroNegDashLineContour(plotT(jmod))
                        end if
                        delete(trop@tiMainString)
                        delete(trop@gsnLeftString)
                        delete(trop@gsnCenterString)
                        delete(trop@gsnRightString)

                        delete(mind)
                        delete(maxd)
                    end if

                    ;; ---------------- GLOBAL CE PLOTS ---------------------
                    mind = min(C)
                    maxd = max(C)
                    res@tiMainString = case1 + " - " + case2
                    res@gsnLeftString = "mean = " + sprintf("%6.2f", gblmean_C)
                    res@gsnCenterString = "rmse = " + sprintf("%6.2f", gblrmse)
                    res@gsnRightString = info@units

                    if (color_type .eq. "COLOR") then
                        res@lbTitleString = "Min = " + sprintf("%6.2f", mind) \
                                             + " Max = " + sprintf("%6.2f", maxd)
                        res@cnLevels = cntrs
                        ;; plot differences
                        plot(jmod) = gsn_csm_contour_map_ce(wks, C, res)
                        delete(res@cnLevels)
                    else
                        mnmxint = nice_mnmxintvl(mind, maxd, diffcnlvls, False)
                        res@cnMinLevelValF = mnmxint(0)
                        res@cnMaxLevelValF = mnmxint(1)
                        res@cnLevelSpacingF = mnmxint(2)
                        plot(jmod) = gsn_csm_contour_map_ce(wks, C, res)
                        plot(jmod) = ZeroNegDashLineContour(plot(jmod))
                    end if
                    delete(res@tiMainString)
                    delete(res@gsnLeftString)
                    delete(res@gsnRightString)
                    delete(res@gsnCenterString)

                    delete(mind)
                    delete(maxd)

                    delete(lat2)
                    delete(nlat2)
                    delete(lon2)
                    delete(gw2)
                    delete(B)
                    delete(gblmean_B)

                    delete(C)
                    delete(gblrmse)
                    delete(troprmse)
                    delete(gblmean_C)
                    delete(tropmean_C)

                    jmod = jmod + 1
                end if
            end do  ; imod

            delete(lat1)
            delete(nlat1)
            delete(lon1)
            delete(gw1)
            delete(A)
            delete(gblmean_A)

        else
            ;; (a) full field
            plot = new(dim_MOD, "graphic")
            plotT = new(dim_MOD, "graphic")

            do imod = 0, dim_MOD - 1
                ;; CASE 1
                case1 = models@name(imod) + " " + models@case_name(imod)
                if ((stringtoint(models@end_date(imod)) \
                          - stringtoint(models@start_date(imod))) .gt. 0) then
                    yrs_ave1 = models@start_date(imod) + "-" + models@end_date(imod)
                end if
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

                A = calc_SEASON(A0, season)
                A&lat@units = "degrees_north"
                A&lon@units = "degrees_east"
                delete(A0)

                ;; compute area weighted mean of the fields
                gblmean_A = wgt_areaave(A, gw1, 1.0, 0)

                ;; ---------------- GLOBAL TROPICAL PLOTS ---------------
                if (tropVar) then
                    tropmean = wgt_areaave(A({minlat:maxlat}, :), \
                                              gw1({minlat:maxlat}), 1.0, 0)
                    min1 = min(A({minlat:maxlat}, :))
                    max1 = max(A({minlat:maxlat}, :))
                    if (color_type .eq. "COLOR") then
                        trop@cnLevels = Tcntrs
                        trop@lbTitleString = "Min = " + sprintf("%6.2f", min1) \
                                              + " Max = " + sprintf("%6.2f", max1)
                    else
                        mnmxint = nice_mnmxintvl(min1, max1, tropcnlvls, False)
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

                    plotT(imod) = gsn_csm_contour_map_other(wksT, A, trop)

                    if (color_type .eq. "MONO") then
                        plotT(imod) = ZeroNegDashLineContour(plotT(imod))
                    else
                        delete(trop@lbTitleString)
                    end if
                    delete(trop@tiMainString)
                    delete(trop@gsnCenterString)

                    delete(min1)
                    delete(max1)
                end if

                ;; ---------------- GLOBAL CE PLOTS ---------------------
                min1 = min(A)
                max1 = max(A)

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
                    mnmxint = nice_mnmxintvl(min1, max1, gblcnlvls, False)
                    res@cnMinLevelValF = mnmxint(0)
                    res@cnMaxLevelValF = mnmxint(1)
                    res@cnLevelSpacingF = mnmxint(2)
                end if

                if (isvar("yrs_ave1")) then
                    res@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
                    delete(yrs_ave1)
                else
                    res@tiMainString = case1
                end if
                res@gsnLeftString = info@long_name
                res@gsnRightString = info@units
                res@gsnCenterString = "mean= " + sprintf("%6.2f", gblmean_A)

                plot(imod) = gsn_csm_contour_map_ce(wks, A, res)   ; plot model data
                if (color_type .eq. "MONO") then
                    plot(imod) = ZeroNegDashLineContour(plot(imod))
                else
                    delete(res@lbTitleString)
                end if
                delete(res@tiMainString)
                delete(res@gsnCenterString)

                delete(min1)
                delete(max1)

                delete(lat1)
                delete(nlat1)
                delete(gw1)
                delete(A)
                delete(gblmean_A)
            end do  ; imod

        end if

        ;; --------------------------------------------------------------
        ;; Panel plots
        ;; common resources for tropics panel plots
        if (tropVar) then   ; Tropics
            panT = True
            panT@gsnMaximize = True
            panT@gsnFrame = False
            panT@txFontHeightF = 0.020
            panT@txFont = "helvetica-bold"
            panT@gsnPaperOrientation = "portrait"
            panT@txString = season

            n_p = dimsizes(plotT) / 2
            if (mod(dimsizes(plotT), 2) .ne. 0) then
                n_p = n_p + 1
            end if
            gsn_panel(wks, plotT, (/n_p, 2/), pan)    ; panel plot
            frame(wksT)
            delete(panT)
        end if

        ;; common resources for panel plots
        pan = True
        pan@gsnMaximize = True
        pan@gsnFrame = False
        pan@gsnPaperOrientation = "portrait"
        pan@gsnPanelYWhiteSpacePercent = 2
        pan@txFontHeightF = 0.018
        pan@txString = season

        n_p = dimsizes(plot) / 2
        if (mod(dimsizes(plot), 2) .ne. 0) then
            n_p = n_p + 1
        end if
        gsn_panel(wks, plot, (/n_p, 2/), pan)    ; panel plot
        frame(wks)

        delete(plot)
        delete(plotT)
        delete(pan)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do  ; s
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete(res)
    delete(cntrs)
    delete(seasons)
    if tropVar then
        delete(trop)
        delete(Tcntrs)
    end if
    delete(tropVar)
end
