; vertconplot_plev.ncl
;***************************************
; NCL script
; plot_vertical_cons.ncl
; Mark Stevens, Sept 2001
;***************************************
load "./plot_type/functions_contrib.ncl"
;load "./plot_type/calc_vertical.ncl"
load "./plot_type/calc_seasonal.ncl"
;***************************************

begin
    wkdir = getenv("ESMValTool_plot_dir")

    if(isatt(info, "rgb_file")) then
        color_type = "COLOR"
        rgb_file = "./rgb/" + info@rgb_file
    else
        color_type = "MONO"
        rgb_file = ""
    end if

    compare = "MODEL"
    plot_type = "PS"
    version = "beta"
    wkdir = wkdir + "/vertconplot_plev/"
    time_stamp = False

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

    if(isatt(info, "vcp_cnlvl") .and. isatt(info, "vcp_dcnlvl") \
                                .and. isatt(info, "vcp_ptop")) then
        con = info@vcp_cnlvl
        dcon = info@vcp_dcnlvl
        ptop = info@vcp_ptop
    else
        print("fatal:plot requires variable attributes 'vcp_cnlvl' and 'vcp_dcnlvl' and 'vcp_ptop'")
        status_exit(1)
    end if

    if(isatt(info, "vcp_LogScale")) then
        logScale = info@vcp_LogScale
    else
        logScale = False
    end if

    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET4 - vertconplot_plev (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    ;; ***********************************************************
    ;; CASE 1
    A0 = extract_data(0, data, -1, 0, 0)
    if isdim(A0, "lon") then
        tmp = dim_avg_Wrap(A0)
        delete(A0)
        A0 = tmp
        delete(tmp)
    end if
    A0&plev@units = "mb"

    lat1 = A0&lat
    nlat1 = dimsizes(lat1)
    if (lat1(0) .gt. lat1(1)) then
        A0 = A0(:, :, ::-1)
        lat1 = lat1(::-1)
    end if

    p1 = A0&plev
    nlev1 = dimsizes(p1)
    if (p1(0) .lt. p1(1)) then
        A0 = A0(:, ::-1, :)
        p1 = p1(::-1)
    end if
    p1&plev = p1
    ip1 = ind_nearest_coord(ptop, p1, 0)  ; nearest index of top level
    np1 = ip1 + 1                         ; number of pressures to plot

    ;; CASE 2
    B0 = extract_data(1, data, -1, 0, 0)
    if isdim(B0, "lon") then
        tmp = dim_avg_Wrap(B0)
        delete(B0)
        B0 = tmp
        delete(tmp)
    end if
    B0&plev@units = "mb"

    lat2 = B0&lat
    nlat2 = dimsizes(lat2)
    if (lat2(0) .gt. lat2(1)) then
        B0 = B0(:, :, ::-1)
        lat2 = lat2(::-1)
    end if

    p2 = B0&plev
    nlev2 = dimsizes(p2)
    if (p2(0) .lt. p2(1)) then
        B0 = B0(:, ::-1, :)
        p2 = p2(::-1)
    end if
    p2&plev = p2
    ip2 = ind_nearest_coord(ptop, p2, 0)  ; nearest index of top level
    np2 = ip2 + 1              ; number of pressures to plot

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    seasons = (/"ANN", "DJF", "JJA"/)

    do s = 0, 2
        season = seasons(s)

        A = calc_SEASON(A0, season)
        A&lat@units = "degrees_north"

        B = calc_SEASON(B0, season)
        B&lat@units = "degrees_north"

        if (all(ismissing(A)) .or. all(ismissing(B))) then
            print("vertconplot_plev.ncl: " + variable + " not present in input files.")
            delete(A)
            delete(B)
            continue
        end if

        ;; ******************************************************************
        ;; ------------------------- Open files for plots ----------------------
        if (color_type .eq. "COLOR") then
            wks = gsn_open_wks(plot_type, wkdir + "vertconplot_plev_" \
                                                + season + "_" + variable + "_c")
            cmap = RGBtoCmap(rgb_file)
            gsn_define_colormap(wks, cmap)
        else
            wks = gsn_open_wks(plot_type, wkdir + "vertconplot_plev_" \
                                                + season + "_" + variable)
        end if

        ;; *********************************************************************
        ;; common resources for contour plots
        cnplot = new(3, "graphic")
        cres                        = True
        cres@gsnDraw                = False
        cres@gsnFrame               = False
        cres@trXReverse             = True
        cres@tiYAxisString          = "Pressure (mb)"
        cres@txFontHeightF = 0.023
        if (color_type .eq. "COLOR") then
            cres@cnLevelSelectionMode = "ExplicitLevels"
            cres@cnFillOn             = True
            cres@cnLinesOn            = False
            cres@cnLineLabelsOn       = False
            cres@lbTitleOn = True
            cres@lbLabelFontHeightF = 0.018
            cres@lbTitleFontHeightF = 0.02
        else
            cres@cnLevelSelectionMode = "ManualLevels"
            cres@cnLinesOn              = True
            cres@cnLineLabelsOn         = True
            cres@cnLineLabelInterval    = 1
            cres@cnLineLabelFontHeightF = .02
            cres@cnInfoLabelFontHeightF = .02
            cnlvls = 15
            diffcnlvls = 12
        end if

        ;; -------------------------------------------
        ;; common resoureces for panel plot
        pan = True
        pan@gsnMaximize = True
        pan@gsnPaperOrientation = "portrait"
        pan@txString = season
        pan@txFontHeightF = 0.018
        pan@gsnFrame = False

        ;; ******************************************************************
        ;; ------------------------------------------------------------------
        ;;  for comparison with OBS data interpolate the model data to these
        ;;  16 constant pressure levels
        ;; plvlO = (/30., 50., 70., 100., 150., 200., 250., 300., 400., 500., 600., 700., \
        ;;          775., 850., 925., 1000./)
        ;; plvlO@units = "mb"
        ;;  for comparison with another model interpolate the model data to these
        ;;  25 constant pressure levels
        ;; plvlM = (/30., 50., 70., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., \
        ;;          600., 650., 700., 750., 800., 850., 875., 900., 925., 950., 975., 1000./)
        ;; plvlM@units = "mb"
        ;; if (compare .eq. "OBS") then
        ;;   plev = plvlO(::-1)               ; flip to bottom to top
        ;; else
        ;;     plev = plvlM(::-1)               ; flip to bottom to top
        ;; end if
        ;; itop = ind(plev .eq. ptop)       ; index of top level
        cres@trYMinF = max((/p1(ip1), p2(ip2), ptop/))

        ;; itop1 = itop
        ;; itop2 = itop
        ;; data1 = A(0:itop1, :)          ; model data(lev, lat)
        ;; data2 = B(0:itop2, :)          ; model or obs data(lev, lat)
        data1 = A(0:ip1, :)          ; model data(lev, lat)
        data2 = B(0:ip2, :)          ; model or obs data(lev, lat)

        ;; ----------------------------------------------------------
        ;; difference plot of latitude(x) vs. pressure(y)
        ;; check for different number of latitudes and then
        ;; lineary interpolate to the smaller number of latitudes
        if (nlat1 .ne. nlat2) then           ; lat grids different
            if (nlat1 .gt. nlat2) then
                data1_lat = linint1_Wrap(lat1, data1, False, lat2, 0)
                data2_lat = data2
            else
                data1_lat = data1
                data2_lat = linint1_Wrap(lat2, data2, False, lat1, 0)
            end if
        else                                 ; same grid latitudes
            data1_lat = data1
            data2_lat = data2
        end if

        ;; check for different number of latitudes and then
        ;; lineary interpolate to the smaller number of latitudes
        dumP1 = p1(0:ip1)
        dumP2 = p2(0:ip2)
        if (np1 .ne. np2) then               ; pressure grids different
            if (np1 .gt. np2) then
                tmp1 = data2_lat(lat|:, plev|:)
                tmp1 = int2p(dumP1, data1_lat(lat|:, plev|:), dumP2, 2)
                tmp1_rev = tmp1(plev|:, lat|:)
                C = data2_lat
                C = (/tmp1_rev - data2_lat/)               ; get diff values
                delete(tmp1)
                delete(tmp1_rev)
            else
                tmp2 = data1_lat(lat|:, plev|:)
                tmp2 = int2p(dumP2, data2_lat(lat|:, plev|:), dumP1, 2)
                tmp2_rev = tmp2(plev|:, lat|:)
                C = data1_lat                    ; copy dims, coords
                C = (/data1_lat - tmp2_rev/)           ; get diff values
                delete(tmp2)
                delete(tmp2_rev)
            end if
        else                                 ; same grid pressure
            dumPDiff = dumP1
            dumPDiff = (/dumP1 - dumP2/)
            idxP = ind(dumPDiff .ne. 0.0)
            if(all(ismissing(idxP))) then
                C = data1_lat                      ; copy dims, coords
                C = (/data1_lat - data2_lat/)        ; get diff values
            else
                tmp1 = data2_lat(lat|:, plev|:)
                tmp1 = int2p(dumP1, data1_lat(lat|:, plev|:), dumP2, 2)
                tmp1_rev = tmp1(plev|:, lat|:)
                C = data2_lat
                C = (/tmp1_rev - data2_lat/)               ; get diff values
                delete(tmp1)
                delete(tmp1_rev)
            end if
            delete(dumPDiff)
            delete(idxP)
        end if
        delete(data1_lat)
        delete(data2_lat)
        delete(dumP1)
        delete(dumP2)

        ;; ----------------------------------------------------------
        ;; case plots of latitude(x) vs. pressure/height(y)
        cres@trYLog = logScale

        min1 = min(data1)
        max1 = max(data1)
        min2 = min(data2)
        max2 = max(data2)

        if (color_type .eq. "COLOR") then
            cres@lbLabelStride = 2
            cres@lbOrientation = "Horizontal"
            cres@cnLevels = con
            cres@lbTitleString = "MIN = " + sprintf("%6.2f", min1) \
                                  + "  MAX = " + sprintf("%6.2f", max1)
        else
            minab = min((/min1, min2/))
            maxab = min((/max1, max2/))
            mnmxint = nice_mnmxintvl(minab, maxab, cnlvls, False)
            cres@cnMinLevelValF = mnmxint(0)
            cres@cnMaxLevelValF = mnmxint(1)
            cres@cnLevelSpacingF = mnmxint(2)
        end if

        if (isvar("yrs_ave1")) then
            cres@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
        else
            cres@tiMainString = case1
        end if
        cnplot(0)  = gsn_csm_pres_hgt(wks, data1, cres)
        if (color_type .eq. "MONO") then
            cnplot(0) = ZeroNegDashLineContour(cnplot(0))
        else
            delete(cres@lbTitleString)
        end if
        delete(cres@tiMainString)

        if(isvar("yrs_ave2")) then
            cres@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
        else
            cres@tiMainString = case2
        end if
        if(color_type .eq. "COLOR") then
            cres@lbTitleString = "MIN = " + sprintf("%6.2f", min2) \
                                  + "  MAX = " + sprintf("%6.2f", max2)
        end if

        cnplot(1)  = gsn_csm_pres_hgt(wks, data2, cres)
        if(color_type .eq. "MONO") then
            cnplot(1) = ZeroNegDashLineContour(cnplot(1))
        else
            delete(cres@cnLevels)
            delete(cres@lbLabelStride)
            delete(cres@lbOrientation)
            delete(cres@lbTitleString)
        end if
        delete(cres@tiMainString)
        delete(min1)
        delete(max1)
        delete(min2)
        delete(max2)

        mind = min(C)
        maxd = max(C)

        cres@tiMainString = case1 + " - " + case2
        if(color_type .eq. "COLOR") then
            cres@lbLabelStride = 1
            cres@lbOrientation = "Vertical"
            cres@lbBoxMinorExtentF = 0.18
            cres@cnLevels = dcon
            cres@lbTitleString = "MIN = " + sprintf("%6.2f", mind) \
                                  + "  MAX = " + sprintf("%6.2f", maxd)
        else
            mnmxint = nice_mnmxintvl(mind, maxd, diffcnlvls, False)
            cres@cnMinLevelValF = mnmxint(0)
            cres@cnMaxLevelValF = mnmxint(1)
            cres@cnLevelSpacingF = mnmxint(2)
        end if

        cnplot(2) = gsn_csm_pres_hgt(wks, C, cres)    ; contour plot

        if(color_type .eq. "COLOR") then
            delete(cres@lbTitleString)
            delete(cres@lbOrientation)
            delete(cres@lbLabelStride)
            delete(cres@lbBoxMinorExtentF)
            delete(cres@cnLevels)
        else
            cnplot(2) = ZeroNegDashLineContour(cnplot(2))
        end if

        if(time_stamp .eq. "True") then
            pan@gsnPanelBottom = 0.05
            gsn_panel(wks, cnplot, (/2, 2/), pan)         ; panel plot
            infoTimeStamp(wks, 0.011, "DIAG Version: " + version)
        else
            gsn_panel(wks, cnplot, (/2, 2/), pan)         ; panel plot
        end if
        frame(wks)

        delete(cres@trYMinF)
        delete(A)
        delete(B)
        delete(C)
        delete(data1)
        delete(data2)
        delete(mind)
        delete(maxd)
        delete(cnplot)
        delete(cres)
        delete(pan)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do  ; s
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete(A0)
    delete(B0)
    delete(lat1)
    delete(lat2)
    delete(nlat1)
    delete(nlat2)
    delete(con)
    delete(dcon)
    delete(ptop)
end
