;*********************************************************************
; seadiffplot.ncl
;*********************************************************************
; NCL script
; plot_season_diff.ncl
; Mark Stevens, Sept 2001
;*********************************************************************
; The original NCL script, plot_season_diff.ncl, was written to make...
; longitude vs latitude DJF-JJA differences contour plots (set9)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; This script takes input from the var_att.ncl files
; info attributes required:
; sdp_cnlvl             vector of contour levels for covering range of
;                       values for difference field
; rgb_file              color map filename
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
;*********************************************************************

begin

    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET9 - seadiffplot (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "C2Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/seadiffplot/"
    compare = "MODEL"
    file_type = "PS"

    dim_MOD = dimsizes(models@name)

    ;; ******************************************************************
    ;;  check the input parameters
    if (isatt(info, "sdp_cnlvl")) then
        cntrs = info@sdp_cnlvl
    else
        print("fatal:plot requires variable attribute 'sdp_cnlvl'")
        status_exit(1)
    end if

    if (isatt(info, "rgb_file")) then
        color_type = "COLOR"
        rgb_file = "./rgb/" + info@rgb_file
    else
        color_type = "MONO"
        rgb_file = ""
    end if

    ;; -------------------------------------------------------------------
    ;;  common resources for plots
    res = True
    res@mpCenterLonF = 180
    res@mpGeophysicalLineThicknessF = 2.0
    res@gsnDraw = False
    res@gsnFrame = False
    res@txFontHeightF = 0.017
    res@tiMainString = ""
    if (color_type .eq. "COLOR") then
        res@lbTitleOn = True
        res@lbOrientation = "Vertical"
        res@lbBoxMinorExtentF = 0.15
        res@lbLabelFontHeightF = 0.013
        res@lbTitleFontHeightF = 0.013
        res@mpFillOn = False
        res@cnFillOn = True
        res@cnLinesOn = False
        res@cnLineLabelsOn = False
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels = cntrs
    else
        ;; resources for black and white plots
        res@mpFillOn = True
        res@cnFillOn = False
        res@cnLinesOn = True
        res@cnLineLabelInterval = 1
        res@cnLineLabelsOn = True
        res@cnLineLabelFontHeightF = .010
        res@cnInfoLabelOrthogonalPosF = -0.06     ; move info label
        res@cnInfoLabelFontHeightF = .010
        res@cnLevelSelectionMode = "ManualLevels"
        cnlvls = 15
    end if
    res@gsnLeftString = info@long_name

    ;; *****************************************************************
    ;; ----------------------- Open files for plots -------------------
    if (color_type .eq. "COLOR") then
        wks = gsn_open_wks(file_type, wkdir + "seadiffplot_" + variable + "_c")
        cmap = RGBtoCmap(rgb_file)
        gsn_define_colormap(wks, cmap)
    else
        wks = gsn_open_wks(file_type, wkdir + "seadiffplot_" + variable)
    end if

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

        ;; convert to T2Ms
        if (field_number .eq. "T2Ds") then
            tmp = day_to_month(A0, field_number)
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if

        A0&lat@units = "degrees_north"
        A0&lon@units = "degrees_east"
        data_DJF = calc_SEASON(A0, "DJF")
        data_JJA = calc_SEASON(A0, "JJA")

        A = data_DJF
        A = (/data_DJF - data_JJA/)
        delete(A0)
        delete(data_DJF)
        delete(data_JJA)

        min1 = min(A)
        max1 = max(A)

        ;; CASE 1 PLOT
        if (isvar("yrs_ave1")) then
            res@gsnCenterString = case1 + " (yrs " + yrs_ave1 + ")"
            delete(yrs_ave1)
        else
            res@gsnCenterString = case1
        end if

        if (color_type .eq. "COLOR") then
            res@lbTitleString = "MIN = " + sprintf("%6.2f", min1) \
                              + " MAX = " + sprintf("%6.2f", max1)
            plot(imod) = gsn_csm_contour_map_ce(wks, A, res)
            delete(res@lbTitleString)
        else
            mnmxint = nice_mnmxintvl(min1, max1, cnlvls, False)
            res@cnMinLevelValF = mnmxint(0)
            res@cnMaxLevelValF = mnmxint(1)
            res@cnLevelSpacingF = mnmxint(2)
            plot(imod) = gsn_csm_contour_map(wks, A, res)
            plot(imod) = ZeroNegDashLineContour(plot(imod))
        end if

        delete(res@gsnCenterString)
        delete(min1)
        delete(max1)
        delete(A)
    end do  ; imod

    ;; ----------------------------------------------------------------
    ;; common resources for panel plot
    pan = True
    pan@gsnMaximize = True
    pan@txString = "DJF-JJA"
    pan@txFont = "helvetica-bold"
    pan@txFontHeightF = 0.02
    pan@gsnFrame = False
    pan@gsnPaperOrientation = "portrait"

    n_p = dimsizes(plot) / 2
    if (mod(dimsizes(plot), 2) .ne. 0) then
        n_p = n_p + 1
    end if
    gsn_panel(wks, plot, (/n_p, 2/), pan)    ; panel plot
    frame(wks)

    delete(plot)
    delete(res)
    delete(pan)
    delete(cntrs)
end
