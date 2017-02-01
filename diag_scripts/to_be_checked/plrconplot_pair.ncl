;*********************************************************************
; plrconplot_pair.ncl
;*********************************************************************
; NCL script
; plot_polar_cons.ncl
; Mark Stevens, Sept 2001
; Jan 2002 updated with significance plots
;*********************************************************************
; The original NCL script, plot_polar_cons.ncl, was written to make...
; longitude vs latitude polar 2D-field contour plots (set7)
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
; pcp_dcnlvl  vector of contour levels covering range of
;             values for difference field
; pcp_cnlvl   vector of contour levels for original field
; rgb_file    color map filename
; seasons     (optional) array of for averaging. (DJF, MAM,
;             JJA, SON or ANN for annual mean)
;             Note: if not set, use default (/"DJF", "JJA", "ANN"/)
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
;*********************************************************************

begin
    print("+++++++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET7 - plrconplot_pair (var: " + variable + ")")
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

    wkdir = plot_dir + "/plrconplot_pair/"
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
    if (isatt(info, "pcp_cnlvl") .and. isatt(info, "pcp_dcnlvl")) then
        Pcntrs = info@pcp_cnlvl
        dPcntrs = info@pcp_dcnlvl
    else
        print("fatal:plot requires variable attributes 'pcp_cnlvl' and 'pcp_dcnlvl'")
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
            print("plrconplot_pair.ncl: " + variable \
                         + " not present in input files.")
            delete(A)
            delete(B)
            continue
        end if
        ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ;; Interpolate the fields to the same grids if needed, and compute
        ;; the difference of the fields.
        if (nlat1 .eq. nlat2) then             ; same grids
            C = A
            C = (/A - B/)
        else                                   ; different grids
            if (nlat1 .gt. nlat2) then           ; interp a to nlat2
                tmp = linint2(lon1, lat1, A, True, lon2, lat2, 0)
                C = B
                C = (/tmp - B/)
            else                                ; interp b to nlat1
                tmp = linint2(lon2, lat2, B, True, lon1, lat1, 0)
                C = A
                C = (/A - tmp/)
            end if
            delete(tmp)
        end if

        ;; ***************************************************************
        ;; ------------------ Open file for plots ------------------------
        if (color_type .eq. "COLOR") then
            cmap = RGBtoCmap(rgb_file)
            wksNP = gsn_open_wks(file_type, wkdir + "plrconplot_pair_" \
                                           + season + "_" + variable + "_NP_c")
            wksSP = gsn_open_wks(file_type, wkdir + "plrconplot_pair_" \
                                           + season + "_" + variable + "_SP_c")
            gsn_define_colormap(wksNP, cmap)
            gsn_define_colormap(wksSP, cmap)
        else     ; Black and white
            wksNP = gsn_open_wks(file_type, wkdir + "plrconplot_pair_" \
                                           + season + "_" + variable + "_NP")
            wksSP = gsn_open_wks(file_type, wkdir + "plrconplot_pair_" \
                                           + season + "_" + variable + "_SP")
        end if

        ;; ***************************************************************
        ;; --------------------------------------------------------------
        plotNP = new((/3/), "graphic")    ; NH polar plots
        plotSP = new((/3/), "graphic")    ; SH polar plots

        ;; ---------------------------------------------------------------
        ;; common resources for north and south polar plots
        spmaxlat = -55
        npminlat =  50
        polar = True
        polar@gsnDraw  = False
        polar@gsnFrame = False
        polar@gsnAddCyclic = True
        if (tick_marks .eq. "True") then
            polar@gsnTickMarksOn = True
        else
            polar@gsnTickMarksOn = False
        end if
        polar@mpGridAndLimbOn   = False
        polar@tiMainFontHeightF = 0.025
        polar@txFontHeightF = 0.022
        polar@mpGeophysicalLineThicknessF = 1.5
        polar@mpPerimOn = True
        polar@mpPerimDrawOrder = "PostDraw"
        if (color_type .eq. "COLOR") then
            polar@mpFillOn = True           ; land fill off
            polar@mpLandFillColor      = "Background"
            polar@cnLevelSelectionMode = "ExplicitLevels"
            polar@cnFillOn           = True
            polar@cnLinesOn          = False
            polar@cnLineLabelsOn     = False
            polar@lbLabelFontHeightF = 0.017
            polar@lbTitleFontHeightF = 0.018
            polar@lbBoxMinorExtentF  = 0.18
        else
            polar@mpFillOn =  True           ; land fill on
            polar@cnLevelSelectionMode = "ManualLevels"
            polar@cnLinesOn              = True
            polar@cnLineLabelsOn         = True
            polar@cnLineLabelInterval    = 1
            polar@cnLineLabelFontHeightF = .018
            polar@cnInfoLabelFontHeightF = .018
            cnlvls                       = 15
            dcnlvls                      = 12
        end if

        ;; ----------------------- North Polar Plots ---------------------
        if (color_type .eq. "COLOR") then
            if (variable .eq. "SNOWHICE" .or. variable .eq. "SNOWHLND" .or. \
                variable .eq. "TREFHT_LAND") then
                polar@cnRasterModeOn = True
            else
                polar@cnRasterModeOn = False
            end if
            if (variable .eq. "ICEFRAC" .or. variable .eq. "SNOWHICE" .or. \
                variable .eq. "SNOWHLND") then
                if (variable .eq. "ICEFRAC") then
                    polar@mpFillDrawOrder = "PostDraw"
                    polar@lbLabelStride = 1
                else
                    polar@mpFillDrawOrder = "Draw"
                    polar@lbLabelStride = 1
                end if
            else
                polar@mpFillDrawOrder = "Draw"
                polar@lbLabelStride = 2
            end if
        end if

        polar@gsnLeftString = info@long_name
        polar@gsnRightString = info@units
        polar@gsnPolar = "NH"
        polar@mpMinLatF = npminlat

        ;; NP CASE 1
        if (isvar("yrs_ave1")) then
            polar@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
        else
            polar@tiMainString = case1
        end if

        meannp1 = wgt_areaave(A({npminlat:}, :), gw1({npminlat:}), 1., 0)
        minnp1 = min(A({npminlat:}, :))
        maxnp1 = max(A({npminlat:}, :))

        meannp2 = wgt_areaave(B({npminlat:}, :), gw2({npminlat:}), 1., 0)
        minnp2 = min(B({npminlat:}, :))
        maxnp2 = max(B({npminlat:}, :))

        if (color_type .eq. "COLOR") then
            polar@lbOrientation   = "Horizontal"
            polar@lbTitleOffsetF = -.25
            polar@cnLevels = Pcntrs
            polar@lbTitleString = "MEAN= " + sprintf("%6.2f", meannp1) \
                                   + "  Min= " + sprintf("%6.2f", minnp1) \
                                   + "  Max= " + sprintf("%6.2f", maxnp1)
        else
            polar@cnInfoLabelParallelPosF = 0.65
            minab = min((/minnp1, minnp2/))
            maxab = min((/maxnp1, maxnp2/))
            mnmxint = nice_mnmxintvl(minab, maxab, cnlvls, False)
            polar@cnMinLevelValF = mnmxint(0)
            polar@cnMaxLevelValF = mnmxint(1)
            polar@cnLevelSpacingF = mnmxint(2)
        end if

        plotNP(0) = gsn_csm_contour_map_polar(wksNP, A, polar)
        if (color_type .eq. "MONO") then
            plotNP(0) = ZeroNegDashLineContour(plotNP(0))
        end if

        ;; NP CASE 2
        if (isvar("yrs_ave2")) then
            polar@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
        else
            polar@tiMainString = case2
        end if

        if (color_type .eq. "COLOR") then
            polar@lbTitleString = "MEAN= " + sprintf("%6.2f", meannp2) \
                                   + "  Min= " + sprintf("%6.2f", minnp2) \
                                   + "  Max= " + sprintf("%6.2f", maxnp2)
        else
            polar@cnInfoLabelParallelPosF = 1.10
        end if

        plotNP(1) = gsn_csm_contour_map_polar(wksNP, B, polar)
        if (color_type .eq. "MONO") then
            plotNP(1) = ZeroNegDashLineContour(plotNP(1))
        end if
        delete(polar@mpMinLatF)

        ;; ------------------------- South Polar Plots -------------------
        polar@gsnPolar = "SH"
        polar@mpMaxLatF = spmaxlat

        meansp1 = wgt_areaave(A({:spmaxlat}, :), gw1({:spmaxlat}), 1., 0)
        minsp1 = min(A({:spmaxlat}, :))
        maxsp1 = max(A({:spmaxlat}, :))
        meansp2 = wgt_areaave(B({:spmaxlat}, :), gw2({:spmaxlat}), 1., 0)
        minsp2 = min(B({:spmaxlat}, :))
        maxsp2 = max(B({:spmaxlat}, :))

        ;; SP CASE 1
        if (isvar("yrs_ave1")) then
            polar@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
        else
            polar@tiMainString = case1
        end if

        if (color_type .eq. "COLOR") then
            if (variable .eq. "PSL") then
                polar@cnLevels = (/Pcntrs - 18/)
            end if
            polar@lbTitleString = "MEAN= " + sprintf("%6.2f", meansp1) \
                                   + "  Min= " + sprintf("%6.2f", minsp1) \
                                   + "  Max= " + sprintf("%6.2f", maxsp1)
        else
            polar@cnInfoLabelParallelPosF = 0.65
            minab = min((/minsp1, minsp2/))
            maxab = min((/maxsp1, maxsp2/))
            mnmxint = nice_mnmxintvl(minab, maxab, cnlvls, False)
            polar@cnMinLevelValF = mnmxint(0)
            polar@cnMaxLevelValF = mnmxint(1)
            polar@cnLevelSpacingF = mnmxint(2)
        end if

        plotSP(0) = gsn_csm_contour_map_polar(wksSP, A, polar)
        if (color_type .eq. "MONO") then
            plotSP(0) = ZeroNegDashLineContour(plotSP(0))
        end if

        ;; SP CASE 2
        if (isvar("yrs_ave2")) then
            polar@tiMainString = case2 + " (yrs " + yrs_ave2 + ")"
        else
            polar@tiMainString = case2
        end if

        if (color_type .eq. "COLOR") then
            polar@lbTitleString = "MEAN= " + sprintf("%6.2f", meansp2) \
                                   + "  Min= " + sprintf("%6.2f", minsp2) \
                                   + "  Max= " + sprintf("%6.2f", maxsp2)
        else
            polar@cnInfoLabelParallelPosF = 1.10
        end if

        plotSP(1) = gsn_csm_contour_map_polar(wksSP, B, polar)
        if (color_type .eq. "MONO") then
            plotSP(1) = ZeroNegDashLineContour(plotSP(1))
        else
            delete(polar@lbTitleString)
            delete(polar@lbOrientation)
            delete(polar@cnLevels)
            delete(polar@lbTitleOffsetF)
        end if
        delete(polar@mpMaxLatF)
        delete(polar@tiMainString)

        ;; MAKE DIFFERENCE PLOTS
        ;; NORTH POLE
        polar@gsnPolar = "NH"
        polar@mpMinLatF = npminlat
        polar@tiMainString = case1 + " - " + case2

        mind = min(C({npminlat:}, :))
        maxd = max(C({npminlat:}, :))

        if (color_type .eq. "COLOR") then
            polar@cnLevels = dPcntrs
            polar@lbLabelStride = 1
            polar@lbOrientation = "Vertical"
            polar@lbTitleString = "MIN = " + sprintf("%6.2f", mind) \
                                   + " MAX = " + sprintf("%6.2f", maxd)
            plotNP(2) = gsn_csm_contour_map_polar(wksNP, C, polar)
            delete(polar@cnLevels)
        else
            polar@cnInfoLabelParallelPosF = 0.80
            mnmxint = nice_mnmxintvl(mind, maxd, dcnlvls, False)
            polar@cnMinLevelValF = mnmxint(0)
            polar@cnMaxLevelValF = mnmxint(1)
            polar@cnLevelSpacingF = mnmxint(2)
            plotNP(2) = gsn_csm_contour_map_polar(wksNP, C, polar)
            plotNP(2) = ZeroNegDashLineContour(plotNP(2))
        end if
        delete(polar@tiMainString)
        delete(polar@mpMinLatF)

        ;; SOUTH POLE
        polar@gsnPolar = "SH"              ; south pole
        polar@mpMaxLatF = spmaxlat
        polar@tiMainString = case1 + " - " + case2
        polar@gsnLeftString = A@long_name
        polar@gsnRightString = A@units

        mind = min(C({:spmaxlat}, :))
        maxd = max(C({:spmaxlat}, :))
        if (color_type .eq. "COLOR") then
            polar@cnLevels = dPcntrs
            polar@lbTitleString = "MIN = " + sprintf("%6.2f", mind) \
                                   + " MAX = " + sprintf("%6.2f", maxd)
            plotSP(2) = gsn_csm_contour_map_polar(wksSP, C, polar)
            delete(polar@lbTitleString)
            delete(polar@lbOrientation)
            delete(polar@cnLevels)
        else
            polar@cnInfoLabelParallelPosF = 0.80
            mnmxint = nice_mnmxintvl(mind, maxd, dcnlvls, False)
            polar@cnMinLevelValF = mnmxint(0)
            polar@cnMaxLevelValF = mnmxint(1)
            polar@cnLevelSpacingF = mnmxint(2)
            plotSP(2) = gsn_csm_contour_map_polar(wksSP, C, polar)
            plotSP(2) = ZeroNegDashLineContour(plotSP(2))
        end if
        delete(polar@tiMainString)

        delete(polar@gsnLeftString)
        delete(polar@gsnRightString)
        delete(polar@mpMaxLatF)
        delete(mind)
        delete(maxd)
        delete(A)
        delete(B)
        delete(C)

        ;; ----------------------------------------------------
        ;; PANEL PLOTS
        ;; common resources for planel plots
        panP = True
        panP@gsnMaximize = True
        panP@gsnPaperOrientation = "portrait"
        panP@txFontHeightF = 0.020
        panP@txFont = "helvetica-bold"
        panP@gsnFrame = False
        panP@txString = season
        gsn_panel(wksNP, plotNP, (/2, 2/), panP)
        gsn_panel(wksSP, plotSP, (/2, 2/), panP)

        frame(wksNP)
        frame(wksSP)

        delete(panP)
        delete(polar)
        delete(plotNP)
        delete(plotSP)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do  ; s
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete(lat1)
    delete(lat2)
    delete(nlat1)
    delete(nlat2)
    delete(lon1)
    delete(lon2)
    delete(gw1)
    delete(gw2)
    delete(A0)
    delete(B0)
    delete(seasons)
    delete(Pcntrs)
    delete(dPcntrs)

end
