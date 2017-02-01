;*********************************************************************
; vertconplot.ncl
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
; The switch 'refModel' stands for the reference model. If present,
; differences with the other models are displayed. If omitted,
; all models are displayed individually.
;
; This script takes input from the var_att.ncl files
; info attributes required:
; vcp_ptop              top pressure limit for plot
; vcp_dcnlvl            vector of contour levels covering range of
;                       values for difference field
; vcp_cnlvl             vector of contour levels for original field
; refModel              reference model name (used with vcp_dcnlvl)
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
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET4 - vertconplot (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T3M", "T2Mz", "T2Dz", "C3M", "C2Mz"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/vertconplot/"
    compare = "MODEL"
    file_type = "PS"

   ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;; ******************************************************************
    ;; check the input parameters
    if (isatt(info, "refModel")) then
        idx_ref = ind(models@name .eq. info@refModel)
        if (all(ismissing(idx_ref))) then
            print("fatal:refModel isn't included in model list (namelist)")
            status_exit(1)
        end if
        if (isatt(info, "vcp_dcnlvl") .and. isatt(info, "vcp_ptop")) then
            con = info@vcp_dcnlvl
            ptop = info@vcp_ptop
        else
            print("fatal:plot requires variable attributes 'vcp_dcnlvl' and 'vcp_ptop'")
            status_exit(1)
        end if
        p_ref = "ref_"
    else
        idx_ref = -1
        if (isatt(info, "vcp_cnlvl") .and. isatt(info, "vcp_ptop")) then
            con = info@vcp_cnlvl
            ptop = info@vcp_ptop
        else
            print("fatal:plot requires variable attributes 'vcp_cnlvl' and 'vcp_ptop'")
            status_exit(1)
        end if
        p_ref = ""
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
    ;; common resources for contour plots
    cres                        = True
    cres@gsnDraw                = False
    cres@gsnFrame               = False
    cres@trXReverse             = True
    cres@trYLog                 = logScale
    cres@tiYAxisString          = "Pressure (mb)"
    cres@txFontHeightF = 0.023
    if (color_type .eq. "COLOR") then
        cres@cnLevelSelectionMode = "ExplicitLevels"
        cres@cnFillOn             = True
        cres@cnLinesOn            = False
        cres@cnLineLabelsOn       = False
        cres@lbTitleOn            = True
        cres@lbLabelFontHeightF   = 0.018
        cres@lbTitleFontHeightF   = 0.02
    else
        cres@cnLevelSelectionMode = "ManualLevels"
        cres@cnLinesOn            = True
        cres@cnLineLabelsOn       = True
        cres@cnLineLabelInterval  = 1
        cres@cnLineLabelFontHeightF = .02
        cres@cnInfoLabelFontHeightF = .02
        if (idx_ref .ne. -1) then
            cnlvls = 15
        else
            cnlvls = 12
        end if
    end if
    cres@gsnLeftString = info@long_name
    cres@gsnRightString = info@units

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;; loop through seasons
    do s = 0, dimsizes(seasons) - 1
        season = seasons(s)
        print("========== SEASON: " + season + " ==========")

        ;; ***************************************************************
        ;; --------------------- Open files for plots --------------------
        if (color_type .eq. "COLOR") then
            wks = gsn_open_wks(file_type, wkdir + "vertconplot_" + p_ref \
                                                + season + "_" + variable + "_c")
            cmap = RGBtoCmap(rgb_file)
            gsn_define_colormap(wks, cmap)
        else
            wks = gsn_open_wks(file_type, wkdir + "vertconplot_" + p_ref \
                                                + season + "_" + variable)
        end if

        if (idx_ref .ne. -1) then
            ;; (b) difference of each model from 'reference'
            cnplot = new(dim_MOD - 1, "graphic")

            ;; CASE 1
            case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
            A0 = extract_data(idx_ref, data, -1, 0, 0)
            uniform_coordinates(A0, field_number)
            ;; convert to T2Mz/C2Mz
            if ((field_number .eq. "T3M") .or. (field_number .eq. "C3M")) then
                tmp = dim_avg_Wrap(A0)
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if
            ;; convert to T2Mz
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
            np1 = itop1 + 1                      ; number of pressures to plot

            lat1 = A0&lat
            nlat1 = dimsizes(lat1)

            A = calc_SEASON(A0, season)
            A&lat@units = "degrees_north"
            data1 = A(0:itop1, :)         ; data(lev, lat)
            delete(A0)
            delete(A)

            jmod = 0
            do imod = 0, dim_MOD - 1
                if (imod .ne. idx_ref) then
                    ;;  CASE 2
                    case2 = models@name(imod) + " " + models@case_name(imod)
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

                    ;; nearest index of top level
                    itop2 = ind_nearest_coord(ptop, plvl2, 0)
                    np2 = itop2 + 1               ; number of pressures to plot

                    lat2 = B0&lat
                    nlat2 = dimsizes(lat2)

                    B = calc_SEASON(B0, season)
                    B&lat@units = "degrees_north"
                    data2 = B(0:itop2, :)          ; data(lev, lat)
                    delete(B0)
                    delete(B)

                    ;; ------------------------------------------------------
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

                    ;; check for different number of pressure_levels and then
                    ;; interpolate to the smaller number of pressure_levels (int2p)
                    dumP1 = plvl1(0:itop1)
                    dumP2 = plvl2(0:itop2)
                    ;; pressure grids different
                    if (np1 .ne. np2) then
                        if (np1 .gt. np2) then
                            tmp1 = data2_lat(lat|:, plev|:)
                            tmp1 = int2p(dumP1, data1_lat(lat|:, plev|:), dumP2, 2)
                            tmp1_rev = tmp1(plev|:, lat|:)
                            C = data2_lat
                            C = (/tmp1_rev - data2_lat/)     ; get diff values
                            delete(tmp1)
                            delete(tmp1_rev)
                        else
                            tmp2 = data1_lat(lat|:, plev|:)
                            tmp2 = int2p(dumP2, data2_lat(lat|:, plev|:), dumP1, 2)
                            tmp2_rev = tmp2(plev|:, lat|:)
                            C = data1_lat                    ; copy dims, coords
                            C = (/data1_lat - tmp2_rev/)     ; get diff values
                            delete(tmp2)
                            delete(tmp2_rev)
                        end if
                    else                                 ; same grid pressure
                        dumPDiff = dumP1
                        dumPDiff = (/dumP1 - dumP2/)
                        idxP = ind(dumPDiff .ne. 0.0)
                        if (all(ismissing(idxP))) then
                            C = data1_lat                    ; copy dims, coords
                            C = (/data1_lat - data2_lat/)    ; get diff values
                        else
                            tmp1 = data2_lat(lat|:, plev|:)
                            tmp1 = int2p(dumP1, data1_lat(lat|:, plev|:), dumP2, 2)
                            tmp1_rev = tmp1(plev|:, lat|:)
                            C = data2_lat
                            C = (/tmp1_rev - data2_lat/)      ; get diff values
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

                    mind = min(C)
                    maxd = max(C)

                    ;; ------------------------------------------------------
                    ;; case plots of latitude(x) vs. pressure/height(y)
                    cres@tiMainString = case1 + " - " + case2
                    if (color_type .eq. "COLOR") then
                        cres@lbLabelStride = 2
                        cres@lbOrientation = "Horizontal"
                        cres@lbBoxMinorExtentF = 0.18
                        cres@cnLevels = con
                        cres@lbTitleString = "MIN = " + sprintf("%6.2f", mind) \
                                              + " MAX = " + sprintf("%6.2f", maxd)
                    else
                        mnmxint = nice_mnmxintvl(mind, maxd, cnlvls, False)
                        cres@cnMinLevelValF = mnmxint(0)
                        cres@cnMaxLevelValF = mnmxint(1)
                        cres@cnLevelSpacingF = mnmxint(2)
                    end if

                    cnplot(jmod) = gsn_csm_pres_hgt(wks, C, cres)    ; contour plot

                    if (color_type .eq. "COLOR") then
                        delete(cres@lbTitleString)
                        delete(cres@lbOrientation)
                        delete(cres@lbLabelStride)
                        delete(cres@lbBoxMinorExtentF)
                        delete(cres@cnLevels)
                    else
                        cnplot(jmod) = ZeroNegDashLineContour(cnplot(jmod))
                    end if

                    delete(C)
                    delete(data2)
                    delete(lat2)
                    delete(nlat2)
                    delete(plvl2)
                    delete(np2)
                    delete(mind)
                    delete(maxd)

                    jmod = jmod + 1
                end if
            end do  ; imod

            delete(data1)
            delete(lat1)
            delete(nlat1)
            delete(plvl1)
            delete(np1)

        else
            ;; (a) full field
            cnplot = new(dim_MOD, "graphic")

            do imod = 0, dim_MOD - 1
                print(models@name(imod))
                ;; CASE 1
                case1 = models@name(imod) + " " + models@case_name(imod)
                if ((stringtoint(models@end_date(imod)) \
                          - stringtoint(models@start_date(imod))) .gt. 0) then
                    yrs_ave1 = models@start_date(imod) + "-" + models@end_date(imod)
                end if
                A0 = extract_data(imod, data, -1, 0, 0)
                uniform_coordinates(A0, field_number)
                ;; convert to T2Mz/C2Mz
                if ((field_number .eq. "T3M") .or. (field_number .eq. "C3M")) then
                    tmp = dim_avg_Wrap(A0)
                    delete(A0)
                    A0 = tmp
                    delete(tmp)
                end if
                ;; convert to T2Mz
                if (field_number .eq. "T2Dz") then
                    tmp = day_to_month(A0, "T2Dz")
                    delete(A0)
                    A0 = tmp
                    delete(tmp)
                end if

                A0&plev@units = "mb"
                plvl = A0&plev

                ;; nearest index of top level
                itop = ind_nearest_coord(ptop, plvl, 0)

                A = calc_SEASON(A0, season)
                A&lat@units = "degrees_north"
                data1 = A(0:itop, :)         ; data(lev, lat)
                delete(A0)
                delete(A)
                delete(plvl)
                delete(itop)

                ;; ------------------------------------------------------------
                ;; case plots of latitude(x) vs. pressure/height(y)
                min1 = min(data1)
                max1 = max(data1)

                if (color_type .eq. "COLOR") then
                    cres@lbLabelStride  = 2
                    cres@lbOrientation = "Horizontal"
                    cres@cnLevels = con
                    cres@lbTitleString = "MIN = " + sprintf("%6.2f", min1) \
                                          + " MAX = " + sprintf("%6.2f", max1)
                else
                    mnmxint = nice_mnmxintvl(min1, max1, cnlvls, False)
                    cres@cnMinLevelValF = mnmxint(0)
                    cres@cnMaxLevelValF = mnmxint(1)
                    cres@cnLevelSpacingF = mnmxint(2)
                end if

                if (isvar("yrs_ave1")) then
                    cres@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
                    delete(yrs_ave1)
                else
                    cres@tiMainString = case1
                end if

                cnplot(imod)  = gsn_csm_pres_hgt(wks, data1, cres)
                if (color_type .eq. "MONO") then
                    cnplot(imod) = ZeroNegDashLineContour(cnplot(imod))
                else
                    delete(cres@cnLevels)
                    delete(cres@lbLabelStride)
                    delete(cres@lbOrientation)
                    delete(cres@lbTitleString)
                end if
                delete(cres@tiMainString)
                delete(min1)
                delete(max1)
                delete(data1)
            end do  ; imod
        end if

        ;; -------------------------------------------
        ;; common resources for panel plot
        pan = True
        pan@gsnMaximize = True
        pan@gsnPaperOrientation = "portrait"
        pan@txString = season
        pan@txFontHeightF = 0.018
        pan@gsnFrame = False

        n_p = dimsizes(cnplot) / 3
        if (mod(dimsizes(cnplot), 3) .ne. 0) then
            n_p = n_p + 1
        end if
        gsn_panel(wks, cnplot, (/n_p, 3/), pan)    ; panel plot
        frame(wks)

        delete(cnplot)
        delete(pan)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do  ; s
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete(cres)
    delete(con)
    delete(seasons)
end
