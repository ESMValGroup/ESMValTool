;*********************************************************************
; zonlnplot.ncl
;*********************************************************************
; NCL script
; plot_zonal_lines.ncl
; Mark Stevens, Sept 2001
;*********************************************************************
; The original NCL script, plot_zonal_lines.ncl, was written to make...
; line plots of zonal means (set3)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; This script takes input from the var_att.ncl files
; info attributes required:
; refModel   reference model name
; seasons    (optional) array of for averaging. (DJF, MAM,
;            JJA, SON or ANN for annual mean)
;            Note: if not set, use default (/"DJF", "JJA", "ANN"/)
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
;*********************************************************************

begin
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET3 - zonlnplot (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "T1Ms", "C2Ms", "C2Ds", "C1Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/zonlnplot/"
    file_type = "PS"
    color_type = "COLOR"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;; ******************************************************************
    ;; check the input parameters
    if (isatt(info, "refModel")) then
        idx_ref = ind(models@name .eq. info@refModel)
        if (all(ismissing(idx_ref))) then
            print("zonlnplot.ncl: refModel isn't included in model list (namelist)")
            exit
        end if
    else
        print("fatal:plot requires variable attribute 'refModel'")
        status_exit(1)
    end if

    if (isatt(info, "seasons")) then
        seasons = info@seasons
    else
        seasons = (/"ANN", "DJF", "JJA"/)
    end if

    ;; ******************************************************************
    colors = CCMVal2_ColorTable(models@name)
    dashs = CCMVal2_DashTable(models@name)
    lgLabels = models@name + " " + models@case_name
    LineThick = CCMVal_Thick(models@name)

    Ccolors = new((/dim_MOD - 1/), string)
    Cdashs = new((/dim_MOD - 1/), integer)
    CLineThick = new((/dim_MOD - 1/), integer)
    CLineThick(:) = 2

    ;; ==================================================================
    ;; setup data holders for plots
    dimXY = new((/dim_MOD/), "integer")
    do imod = 0, dim_MOD - 1
        A0 = extract_data(imod, data, -1, 0, 0)
        lat = A0&lat
        dimXY(imod) = dimsizes(lat)
        delete(A0)
        delete(lat)
    end do

    nMax = max(dimXY)
    Xdata = new((/dim_MOD, nMax/), "float")
    Ydata = new((/dim_MOD, nMax, dimsizes(seasons)/), "float")
    CXdata = new((/dim_MOD - 1, nMax/), "float")
    CYdata = new((/dim_MOD - 1, nMax, dimsizes(seasons)/), "float")
    Xdata!0 = "line"
    Xdata!1 = "lat"

    ;; get data into arrays (X) & (Y)
    do imod = 0, dim_MOD - 1

        ;; assign latitude data (X)
        A0DF = extract_data(imod, data, -1, 0, 0)
        A0 = convert_doubletofloat(A0DF)
        delete(A0DF)
        uniform_coordinates(A0, field_number)
        ;; convert to T2Ms
        if (field_number .eq. "T2Ds") then
            tmp = day_to_month(A0, "T2Ds")
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if
        if isStrSubset(field_number, "T") then
            tmp = extract_time_range(A0, models@start_date(imod), \
                                         models@end_date(imod))

            delete(A0)
            A0 = tmp
            delete(tmp)
        end if
        if (field_number .ne. "T1Ms") then
            tmp = dim_avg_Wrap(A0)
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if

        lat = A0&lat
        Xdata(imod, 0:dimXY(imod) - 1) = (/lat/)

        ;;  assign data (Y)
        do s = 0, dimsizes(seasons) - 1  ; loop through seasons
            season = seasons(s)
            if (models@name(imod) .eq. "CAM3.5" .or. models@name(imod) .eq. "WACCM" \
                                  .or. models@name(imod) .eq. "UMETRAC") then
                A = calc_SEASON(A0(0:227, :), season)
            else
                A = calc_SEASON(A0, season)
            end if
            Ydata(imod, 0:dimXY(imod) - 1, s) = (/A/)
            delete(A)
        end do  ; s
        delete(A0)
        delete(lat)
    end do  ; imod

    ;; assign data into array (C)
    ;; CASE 1
    Ydata1 = Ydata(idx_ref, 0:dimXY(idx_ref) - 1, :)
    lat1 = Xdata(idx_ref, 0:dimXY(idx_ref) - 1)
    nlat1 = dimXY(idx_ref)
    Ydata1!0 = "lat"
    Ydata1!1 = "season"
    Ydata1&lat = lat1
    Ydata1&season = ispan(1, dimsizes(seasons), 1)

    jmod = 0
    do imod = 0, dim_MOD - 1
        if (imod .ne. idx_ref) then
            ;; CASE 2
            Ydata2 = Ydata(imod, 0:dimXY(imod) - 1, :)
            lat2 = Xdata(imod, 0:dimXY(imod) - 1)
            nlat2 = dimXY(imod)
            Ydata2!0 = "lat"
            Ydata2!1 = "season"
            Ydata2&lat = lat2
            Ydata2&season = ispan(1, dimsizes(seasons), 1)

            ;; check for different number of latitudes and thenn
            ;; lineary interpolate to the smaller number of latitudes
            if (nlat1 .ne. nlat2) then           ; lat grids different
                if (nlat1 .gt. nlat2) then
                    C = (/Ydata2/)
                    tmp = linint1_Wrap(lat1, Ydata1(season|:, lat|:), False, lat2, 0)   ; case1
                    C = (/tmp(lat|:, season|:) - Ydata2/)
                    CXdata(jmod, 0:dimXY(imod) - 1) = (/lat2/)
                    CYdata(jmod, 0:dimXY(imod) - 1, :) = (/C/)
                else
                    C = (/Ydata1/)
                    tmp = linint1_Wrap(lat2, Ydata2(season|:, lat|:), False, lat1, 0)   ; case2
                    C = (/Ydata1 - tmp(lat|:, season|:)/)
                    CXdata(jmod, 0:dimXY(idx_ref) - 1) = (/lat1/)
                    CYdata(jmod, 0:dimXY(idx_ref) - 1, :) = (/C/)
                end if
                delete(tmp)
            else                                 ; same grid latitudes
                C = (/Ydata1 - Ydata2/)
                CXdata(jmod, 0:dimXY(idx_ref) - 1) = (/lat1/)
                CYdata(jmod, 0:dimXY(idx_ref) - 1, :) = (/C/)
            end if
            delete(C)
            delete(Ydata2)
            delete(lat2)

            Ccolors(jmod) = colors(imod)
            Cdashs(jmod) = dashs(imod)

            jmod = jmod + 1
        end if
    end do  ; imod

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;; loop through seasons
    do s = 0, dimsizes(seasons) - 1
        season = seasons(s)

        ;; ***************************************************************
        ;; --------------------- Open files for plots --------------------
        if (color_type .eq. "COLOR") then
            wks = gsn_open_wks(file_type, wkdir + "zonlnplot_" + season + "_" + variable + "_c")
        else
            wks = gsn_open_wks(file_type, wkdir + "zonlnplot_" + season + "_" + variable)
        end if

        ;; ***************************************************************
        ;;  common resources for xy plots
        plot = new(2, "graphic")
        res                        = True
        res@gsnDraw                = False               ; Do not draw plot
        res@gsnFrame               = False               ; Do not advance frame
        res@trXReverse             = False
        res@tiXAxisString          = "latitude"
        res@tiXAxisFontHeightF     = 0.020
        res@tiYAxisFontHeightF     = 0.020
        res@tiMainFontHeightF      = 0.025
        res@pmLegendSide           = "Right"            ; Change location of

        res@pmLegendOrthogonalPosF = 0.02               ; of the legend
        res@pmLegendWidthF         = 0.15               ; Change width and
        res@lgPerimOn              = True
        if (dim_MOD .le. 4) then
            res@pmLegendHeightF     = 0.10               ; height of legend.
            res@lgLabelFontHeightF  = .020
            res@pmLegendParallelPosF = 0.10               ; of the legend
        else
            if (dim_MOD .le. 8) then
                res@pmLegendHeightF  = 0.18
                res@lgLabelFontHeightF = 0.016
                res@pmLegendParallelPosF = 0.18
            else
                res@pmLegendHeightF  = 0.55
                res@lgLabelFontHeightF = 0.013
                res@pmLegendParallelPosF = 0.30
            end if
        end if
        res@txFontHeightF          = 0.015
        res@xyDashPatterns         = dashs
        res@xyLineThicknesses      = LineThick
        if (color_type .eq. "COLOR") then
            res@xyLineColors        = colors
        end if

        res@tiMainString           = season
        res@pmLegendDisplayMode    = "Always"
        res@tiYAxisString          = info@long_name + " (" + info@units + ")"
        res@xyExplicitLegendLabels = lgLabels
        res@xyMonoLineColor        = False

        if(info@units .eq. "hPa") then
            res@trYReverse = True
        end if

        sYdata = Ydata(:, :, s)
        plot(0) = gsn_csm_xy(wks, Xdata, sYdata, res)   ; cases line plot

        ;; overplot ccmVal levels...
        lls = (/300, 250, 200, 170, 150, 130, 115, 100, 90, 80, 70, 50./)
        nl = dimsizes(lls)
        do i = 0, nl - 1
            add_line(wks, plot(0), (/-90, 90./), (/lls(i), lls(i)/), "black", 2)
        end do

        delete(res@tiMainString)
        delete(res@xyExplicitLegendLabels)       ; no label
        delete(res@pmLegendDisplayMode)          ; default is off
        delete(res@xyDashPatterns)
        delete(res@xyLineThicknesses)
        if (color_type .eq. "COLOR") then
            delete(res@xyLineColors)
            res@xyLineColors = Ccolors
        end if

        sCYdata = CYdata(:, :, s)

        res@xyDashPatterns       = Cdashs
        res@xyLineThicknesses    = CLineThick
        res@tiMainString = models@name(idx_ref) + " - ModelX"
        plot(1) = gsn_csm_xy(wks, CXdata, sCYdata, res)      ; diff line plot
        delete(res@tiMainString)

        ;; ----------------------------------------------
        ;; common resources for panel plots
        pan = True
        pan@gsnMaximize = True
        pan@gsnPaperOrientation = "portrait"
        pan@gsnFrame = False
        gsn_panel(wks, plot, (/2, 1/), pan)         ; panel plot
        frame(wks)

        delete(sYdata)
        delete(sCYdata)
        delete(res)
        delete(plot)
        delete(pan)

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end do  ; s
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    delete(seasons)
    delete(lgLabels)
    delete(colors)
    delete(dashs)
    delete(LineThick)
    delete(Ccolors)
    delete(Cdashs)
    delete(CLineThick)
    delete(Xdata)
    delete(Ydata)
    delete(CXdata)
    delete(CYdata)
end
