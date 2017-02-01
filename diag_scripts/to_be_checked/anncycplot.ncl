;*********************************************************************
; anncycplot.ncl
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
; diagnostic software package
;
; The CCMVal script has been furthered developed as a part of the
; ESMValTool, an ESM Evaluation Tool within the EMBRACE project. Package
; leaders are Klaus Wyser (SMHI) and Veronika Eyring (DLR).
;
; The switch 'refModel' stands for the reference model. If present,
; differences with the other models are displayed. If omitted,
; all models are displayed individually.
;
; This script takes input from the var_att.ncl files
; info attributes required:
; acp_dcnlvl            vector of contour levels covering range of
;                       values for difference field
; acp_cnlvl             vector of contour levels for original field
; refModel              reference model name (used with acp_dcnlvl)
; rgb_file              color map filename
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/calc_cycle.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/write_infobox.ncl"
load "./ncl_code/info_output.ncl"
;*********************************************************************

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering anncycplot.ncl", verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot SET8 - anncycplot (var: " + variable + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "T1Ms", "C2Ms", "C2Ds", "C1Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    supported_info_transforms = (/"latrange", "lonrange"/)

    wkdir = plot_dir + "/anncycplot/"
    compare = "MODEL"
    file_type = "PS"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

;******************************************************************
;check the input parameters

    if (isatt(info, "refModel")) then
        idx_ref = ind(models@name .eq. info@refModel)
        if (all(ismissing(idx_ref))) then
            print("fatal:refModel isn't included in model list (namelist)")
            status_exit(1)
        end if
        if (isatt(info, "acp_dcnlvl")) then
            cntrs = info@acp_dcnlvl
        else
            print("fatal:plot requires variable attribute 'acp_dcnlvl'")
            status_exit(1)
        end if
        p_ref = "ref_"
    else
        idx_ref = -1
        if (isatt(info, "acp_cnlvl")) then
            cntrs = info@acp_cnlvl
        else
            print("fatal:plot requires variable attribute 'acp_cnlvl'")
            status_exit(1)
        end if
        p_ref = ""
    end if

    if(isvar("latrange")) then
        delete(latrange)
    end if
    if(isatt(info, "ac_latrange")) then
        latrange = info@ac_latrange
    end if

    if(isvar("lonrange")) then
        delete(lonrange)
    end if
    if(isatt(info, "ac_lonrange")) then
        lonrange = info@ac_lonrange
    end if

    if (isatt(info, "rgb_file")) then
        color_type = "COLOR"
        rgb_file = "./rgb/" + info@rgb_file
    else
        color_type =  "MONO"
        rgb_file = ""
    end if

;******************************************************************
; common resources for contour plots
    res                        = True
    res@gsnDraw                = False
    res@gsnFrame               = False
    res@txFontHeightF          = 0.018
    res@sfXArray               = ispan(0, 12, 1)
    res@tiMainFontHeightF      = 0.022
    res@tmXBMode               = "Explicit"
    res@tmXBValues             = ispan(0, 12, 1)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
    "O", "N", "D", "J"/)
    if (color_type .eq. "COLOR") then
        res@cnFillOn  = True
        res@cnLinesOn = False
        res@lbTitleOn = True
        res@lbLabelFontHeightF = 0.018
        res@lbTitleFontHeightF = 0.02
        res@lbBoxMinorExtentF = 0.18
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels = cntrs
        res@lbLabelStride = 2
        res@lbOrientation = "Horizontal"
    else
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnFillOn  = False
        res@cnLinesOn = True
        res@cnLineLabelsOn = True
        res@cnLineLabelInterval = 1
        res@cnLineLabelFontHeightF =  0.02
        res@cnInfoLabelFontHeightF =  0.02
        cnlvls = 15
        dcnlvls = 12
    end if
    res@gsnLeftString = info@long_name
    if(isvar("lonrange") .and. inlist("lonrange", supported_info_transforms)) then
        res@gsnLeftString = res@gsnLeftString + " (" + lonrange(0) \
                                              + "<lon<" + lonrange(1) + ")"
        res@gsnLeftStringFontHeightF = 14
    end if
    res@gsnRightString = info@units

;******************************************************************
;----------------------- Open files for plots ----------------------
    if (color_type .eq. "COLOR") then
        wks = gsn_open_wks(file_type, wkdir + "anncycplot_" \
                                            + p_ref + variable + "_c")
        cmap = RGBtoCmap(rgb_file)
        gsn_define_colormap(wks, cmap)
    else
        wks = gsn_open_wks(file_type, wkdir + "anncycplot_" + p_ref + variable)
    end if

    if (idx_ref .ne. -1) then
        ;; (b) difference of each model from 'reference'
        plot = new(dim_MOD - 1, "graphic")

        ;; CASE 1
        case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
        A0 = extract_data(idx_ref, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        ;; convert to T2Ms, C2Ms
        if ((field_number .eq. "T2Ds") .or. (field_number .eq. "C2Ds"))then
            tmp = day_to_month(A0, field_number)
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if

        ;; Select a longitude sub-set to plot
        if (iscoord(A0, "lon") .and. isvar("lonrange") \
                                .and. inlist("lonrange", supported_info_transforms)) then
            Atmp = A0(time|:, lat|:, {lon|lonrange(0):lonrange(1)})
            delete(A0)
            A0 = Atmp
            delete(Atmp)
        end if
        A = calc_cycle(A0)
        A&lat@units = "degrees_north"
        lat1 = A&lat
        nlat1 = dimsizes(lat1)
        delete(A0)

        jmod = 0
        do imod = 0, dim_MOD - 1
            if (imod .ne. idx_ref) then
                ;; CASE 2
                case2 = models@name(imod) + " " + models@case_name(imod)
                B0 = extract_data(imod, data, -1, 0, 0)
                uniform_coordinates(B0, field_number)

                ;; convert to T2Ms, C2Ms
                if ((field_number .eq. "T2Ds") .or. (field_number .eq. "C2Ds"))then
                    tmp = day_to_month(B0, field_number)
                    delete(B0)
                    B0 = tmp
                    delete(tmp)
                end if

                ;; Select a longitude sub-set to plot
                if (iscoord(B0, "lon") .and. isvar("lonrange") \
                                        .and. inlist("lonrange", supported_info_transforms)) then
                    Btmp = B0(time|:, lat|:, {lon|lonrange(0):lonrange(1)})
                    delete(B0)
                    B0 = Btmp
                    delete(Btmp)
                end if
                B = calc_cycle(B0)
                B&lat@units = "degrees_north"
                lat2 = B&lat
                nlat2 = dimsizes(lat2)
                delete(B0)

                ;;---------------------------------------------------------
                ;; difference plot of time(x) vs. latitude(y)
                ;; check for different number of latitudes and then
                ;; lineary interpolate to the smaller number of latitudes

                if (nlat1 .ne. nlat2) then                 ; lat grids different
                    if (nlat1 .gt. nlat2) then
                        tmp = linint1(lat1, A, False, lat2, 0)   ; a(time, lat)
                        C = B                               ; copy dims, coords
                        C = (/tmp - B/)                     ; get diff values
                        res@sfYArray = lat2
                    else
                        tmp = linint1(lat2, B, False, lat1, 0)   ; b(time, lat)
                        C = A                             ; copy dims, coords
                        C = (/A - tmp/)                   ; get diff values
                        res@sfYArray = lat1
                    end if
                    delete(tmp)
                else                                       ; same grid latitudes
                    C = A                                    ; copy dims, coords
                    C = (/A - B/)                              ; get diff values
                    res@sfYArray = lat1
                end if

                mind = min(C)
                maxd = max(C)

                res@tiMainString = case1 + " - " + case2
                if (color_type .eq. "COLOR") then
                    res@lbLabelStride = 1
                    res@lbTitleString = "MIN = " + sprintf("%6.2e", mind) \
                                        + " MAX = " + sprintf("%6.2e", maxd)
                else
                    mnmxint = nice_mnmxintvl(mind, maxd, dcnlvls, False)
                    res@cnMinLevelValF = mnmxint(0)
                    res@cnMaxLevelValF = mnmxint(1)
                    res@cnLevelSpacingF = mnmxint(2)
                end if

                if(isvar("latrange") .and. inlist("latrange", \
                                                  supported_info_transforms))then
                    lat1!0 = "lat"
                    lat1&lat = lat1
                    delete(res@sfYArray)
                    res@sfYArray = lat1({lat|latrange(0):latrange(1)})
                    plot(imod) = gsn_csm_lat_time(wks, C({lat|latrange(0):latrange(1)}, time|:), res)
                else
                    plot(jmod) = gsn_csm_lat_time(wks, C(lat|:, time|:), res)
                end if
                if (color_type .eq. "MONO") then
                    plot(jmod) = ZeroNegDashLineContour(plot(jmod))
                else
                    delete(res@lbTitleString)
                end if

                delete(mind)
                delete(maxd)
                delete(lat2)
                delete(nlat2)
                delete(B)
                delete(C)
                delete(res@sfYArray)

                jmod = jmod + 1
            end if
        end do  ; imod

        delete(lat1)
        delete(nlat1)
        delete(A)

    else
        ;; (a) full field
        plot = new(dim_MOD, "graphic")

        do imod = 0, dim_MOD - 1
            ;; CASE 1
            case1 = models@name(imod) + " " + models@case_name(imod)
            if ((stringtoint(models@end_date(imod)) \
                      - stringtoint(models@start_date(imod))) .gt. 0) then
                yrs_ave1 = models@start_date(imod) + "-" + models@end_date(imod)
            end if
            A0 = extract_data(imod, data, -1, 0, 0)
            uniform_coordinates(A0, field_number)

            ;; convert to T2Ms, C2Ms
            if ((field_number .eq. "T2Ds") .or. (field_number .eq. "C2Ds"))then
                tmp = day_to_month(A0, field_number)
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if

            ;; Select a longitude sub-set to plot
            if (iscoord(A0, "lon") .and. isvar("lonrange") \
                                    .and. inlist("lonrange", supported_info_transforms)) then
                Atmp = A0(time|:, lat|:, {lon|lonrange(0):lonrange(1)})
                delete(A0)
                A0 = Atmp
                delete(Atmp)
            end if
            A = calc_cycle(A0)
            A&lat@units = "degrees_north"
            lat1 = A&lat
            delete(A0)

            ;;------------------------------------------------------------
            ;; case contour plots of time(x) vs. latitude(y)
            min1 = min(A)
            max1 = max(A)

            if (color_type .eq. "COLOR") then
                res@lbLabelStride = 2
                res@lbTitleString = "MIN = " + sprintf("%6.2e", min1) \
                                    + " MAX = " + sprintf("%6.2e", max1)
            else
                mnmxint = nice_mnmxintvl(min1, max1, cnlvls, False)
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
            if(isvar("latrange") .and. inlist("latrange", \
                                               supported_info_transforms)) then
                lat1!0 = "lat"
                lat1&lat = lat1
                res@sfYArray = lat1({lat|latrange(0):latrange(1)})
                plot(imod) = gsn_csm_lat_time(wks, A({lat|latrange(0):latrange(1)}, time|:), res)
            else
                res@sfYArray = lat1
                plot(imod) = gsn_csm_lat_time(wks, A(lat|:, time|:), res)
            end if
            if (color_type .eq. "MONO") then
                plot(imod) = ZeroNegDashLineContour(plot(imod))
            else
                delete(res@lbTitleString)
            end if
            delete(res@sfYArray)
            delete(res@tiMainString)

            delete(min1)
            delete(max1)
            delete(A)
            delete(lat1)
        end do  ; imod

    end if

;----------------------------------------------------
; PANEL PLOTS
    pan = True
    pan@gsnMaximize = True
    pan@gsnFrame = False
    pan@gsnPaperOrientation = "portrait"
    n_p = dimsizes(plot) / 3
    if (mod(dimsizes(plot), 3) .ne. 0) then
        n_p = n_p + 1
    end if
    gsn_panel(wks, plot, (/n_p, 3/), pan)
    frame(wks)

    delete(plot)
    delete(res)
    delete(pan)
;   delete(cntrs)

    info_output(">>>>>>>> Leaving anncycplot.ncl", verbosity, 4)
end
