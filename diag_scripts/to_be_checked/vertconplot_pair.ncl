;*********************************************************************
; vertconplot_pair.ncl
;*********************************************************************
; NCL script
; plot_vertical_cons.ncl
; Mark Stevens, Sept 2001
;*********************************************************************
; The original NCL script, plot_vertical_cons.ncl, was written to make...
; latitude vs pressure/height zonal mean contour plots (set4)
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
; vcp_ptop              top pressure limit for plot
; vcp_dcnlvl            vector of contour levels covering range of
;                       values for difference field
; vcp_cnlvl             vector of contour levels for original field
; rgb_file              color map filename
; seasons               (optional) array of for averaging. (DJF, MAM,
;                       JJA, SON or ANN for annual mean)
;                       Note: if not set, use default (/"DJF", "JJA", "ANN"/)
; vcp_LogScale          (optional) set True for y-axis logscale
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
;*********************************************************************

begin
    print("++++++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET4 - vertconplot_pair (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T3M", "T2Mz", "T2Dz", "C3M", "C2Mz"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/vertconplot_pair/"
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

    ;; ******************************************************************
    ;; check the input parameters
    if (isatt(info, "vcp_cnlvl") .and. isatt(info, "vcp_dcnlvl") \
                                 .and. isatt(info, "vcp_ptop")) then
        con = info@vcp_cnlvl
        dcon = info@vcp_dcnlvl
        ptop = info@vcp_ptop
    else
        print("fatal:plot requires variable attributes 'vcp_cnlvl' and 'vcp_dcnlvl' and 'vcp_ptop'")
        status_exit(1)
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

    if (isatt(info, "vcp_LogScale")) then
        logScale = info@vcp_LogScale
    else
        logScale = False
    end if

    ;; ******************************************************************
    ;; CASE 1
    imod = 0
    A0 = extract_data(imod, data, -1, 0, 0)
    uniform_coordinates(A0, field_number)
    if ((field_number .eq. "T3M") .or. (field_number .eq. "C3M")) then
        tmp = dim_avg_Wrap(A0)
        delete(A0)
        A0 = tmp
        delete(tmp)
    end if
    if (field_number .eq. "T2Dz") then
        tmp = day_to_month(A0, "T2Dz")
        delete(A0)
        A0 = tmp
        delete(tmp)
    end if
    A0&plev@units = "mb"
    plvl1 = A0&plev

    ;; nearest index of top level
    itop1 = ind_nearest_coord(ptop, plvl1, 0)
    np1 = itop1 + 1         ; number of pressures to plot

    lat1 = A0&lat
    nlat1 = dimsizes(lat1)

    ;; CASE 2
    imod = 1
    B0 = extract_data(imod, data, -1, 0, 0)
    uniform_coordinates(B0, field_number)
    if ((field_number .eq. "T3M") .or. (field_number .eq. "C3M")) then
        tmp = dim_avg_Wrap(B0)
        delete(B0)
        B0 = tmp
        delete(tmp)
    end if
    if (field_number .eq. "T2Dz") then
        tmp = day_to_month(B0, "T2Dz")
        delete(B0)
        B0 = tmp
        delete(tmp)
    end if
    B0&plev@units = "mb"
    plvl2 = B0&plev

    itop2 = ind_nearest_coord(ptop, plvl2, 0)  ; nearest index of top level
    np2 = itop2 + 1                            ; number of pressures to plot

    lat2 = B0&lat
    nlat2 = dimsizes(lat2)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;; loop through seasons
    do s = 0, dimsizes(seasons) - 1
        season = seasons(s)

        A = calc_SEASON(A0, season)
        A&lat@units = "degrees_north"

        B = calc_SEASON(B0, season)
        B&lat@units = "degrees_north"

        ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (all(ismissing(A)) .or. all(ismissing(B))) then
            print("vertconplot_pair.ncl: " + variable + " not present in input files.")
            delete(A)
            delete(B)
            continue
        end if

        ;; ***************************************************************
        ;; -------------------- Open files for plots ---------------------
        if (color_type .eq. "COLOR") then
            wks = gsn_open_wks(file_type, wkdir + "vertconplot_pair_" \
                                                + season + "_" + variable + "_c")
            cmap = RGBtoCmap(rgb_file)
            gsn_define_colormap(wks, cmap)
        else
            wks = gsn_open_wks(file_type, wkdir + "vertconplot_pair_" \
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
        pan               = True
        pan@gsnMaximize   = True
        pan@gsnPaperOrientation = "portrait"
        pan@txString      = season
        pan@txFontHeightF = 0.018
        pan@gsnFrame      = False

        cres@trYMinF = max((/ptop, plvl1(itop1), plvl2(itop2)/))

        data1 = A(0:itop1, :)          ; data(lev, lat)
        data2 = B(0:itop2, :)          ; data(lev, lat)

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
        dumP1 = plvl1(0:itop1)
        dumP2 = plvl2(0:itop2)
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
            if (all(ismissing(idxP))) then
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
            cres@lbLabelStride  = 2
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
        cnplot(0) = gsn_csm_pres_hgt(wks, data1, cres)
        if (color_type .eq. "MONO") then
            cnplot(0) = ZeroNegDashLineContour(cnplot(0))
        else
            delete(cres@lbTitleString)
        end if
        delete(cres@tiMainString)

        if (isvar("yrs_ave2")) then
            cres@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
        else
            cres@tiMainString = case2
        end if
        if (color_type .eq. "COLOR") then
            cres@lbTitleString = "MIN = " + sprintf("%6.2f", min2) \
                                  + "  MAX = " + sprintf("%6.2f", max2)
        end if
        cnplot(1) = gsn_csm_pres_hgt(wks, data2, cres)
        if (color_type .eq. "MONO") then
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
        if (color_type .eq. "COLOR") then
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
        if (color_type .eq. "COLOR") then
            delete(cres@lbTitleString)
            delete(cres@lbOrientation)
            delete(cres@lbLabelStride)
            delete(cres@lbBoxMinorExtentF)
            delete(cres@cnLevels)
        else
            cnplot(2) = ZeroNegDashLineContour(cnplot(2))
        end if

        gsn_panel(wks, cnplot, (/2, 2/), pan)         ; panel plot
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
    delete(seasons)

end
