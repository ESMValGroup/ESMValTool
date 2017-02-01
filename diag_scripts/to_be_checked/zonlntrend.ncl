;*********************************************************************
; zonlntrend.ncl
;*********************************************************************
; modified by andrew gettelman to take 2D data, calculate trends
; and plot zonal mean trend.
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
; refModel  reference model name
; seasons  (optional) array of for averaging. (DJF, MAM,
;   JJA, SON or ANN for annual mean)
;   Note: if not set, use default (/"DJF", "JJA", "ANN"/)
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
;*********************************************************************
begin
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot zonlntrend (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "T1Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/zonlntrend/"
    file_type = "PS"
    color_type = "COLOR"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;; ******************************************************************
    ;; check the input parameters
    if (isatt(info, "refModel")) then
        idx = ind(models@name .eq. info@refModel)
        idx_ref = idx(0)
        if (all(ismissing(idx_ref))) then
            print("fatal:refModel isn't included in model list (namelist)")
            status_exit(1)
        end if
    else
        print("fatal:plot requires variable attribute 'refModel'")
        status_exit(1)
    end if

    seasdt = (/ (/0, 12./), (/-1, 2./), (/2, 5./), (/5, 8/), (/8, 11./) /)
    seasons = (/"ANN", "DJF", "MAM", "JJA", "SON"/)
    nmos = floattoint((seasdt(:, 1) - seasdt(:, 0)))
    seasdt = seasdt / 12.

    ns = dimsizes(seasons)

    ;; same as above, but want y + s(0, s):y + s(1, s)
    ;; where  ANN           DJF        MAM     JJA    SON
    ;; 0      0/12.       -1./12.   2/12.   5/12.  8/12
    ;; 1      12./12.       2./12.     5/12.   8/12.  11/12
    ;; only use if there are nmos values in each.
    ;; ******************************************************************
    colors = CCMVal2_ColorTable(models@name)
    dashs = CCMVal2_DashTable(models@name)

    ;; hack to control dashs for namelist_3per_refb2
    ;; dashs = (/0, 0, 1, 1, 2, 2, 0, 0, 0, 1, 1, 1, 2, 2, 2/)
    lgLabels = models@name + " " + models@case_name
    LineThick = new((/dim_MOD/), integer)
    LineThick(:) = 2

    Ccolors = new((/dim_MOD - 1/), string)
    Cdashs = new((/dim_MOD - 1/), integer)
    CLineThick = new((/dim_MOD - 1/), integer)
    CLineThick(:) = 2

    ;; ==================================================================
    ;; setup data holders for plots
    dimXY = new((/dim_MOD/), "integer")
    dimT = new((/dim_MOD/), "integer")
    do imod = 0, dim_MOD - 1
        A0 = extract_data(imod, data, -1, 0, 0)
        lat = A0&lat
        dimXY(imod) = dimsizes(lat)

        delete(lat)
        time = A0&time
        dimT(imod) = dimsizes(time)
        delete(A0)
        delete(time)

    end do

    nMax = max(dimXY)
    Xdata = new((/dim_MOD, nMax/), "float")
    Ydata = new((/dim_MOD, nMax, dimsizes(seasons)/), "float")
    ;; significance of trends or sigma?
    Sdata = new((/dim_MOD, nMax, dimsizes(seasons)/), "float")
    CXdata = new((/dim_MOD - 1, nMax/), "float")
    CYdata = new((/ dim_MOD - 1, nMax, dimsizes(seasons)/), "float")
    Xdata!0 = "line"
    Xdata!1 = "lat"

    tMax = max(dimT)
    Tdata = new((/dim_MOD, tMax/), "float")
    Tdata!0 = "line"
    Tdata!1 = "time"

    aTdata = new((/dim_MOD, tMax/), "float")

    ;; get data into arrays (X) & (Y)
    do imod = 0, dim_MOD - 1
        print("zonlntrend.ncl: " + models@name(imod))
        ;; assign latitude data (X)
        A0 = extract_data(imod, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        ;; convert to T2Ms
        if (field_number .eq. "T2Ds") then
            tmp = day_to_month(A0, "T2Ds")
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if

        lat = A0&lat
        Xdata(imod, 0:dimXY(imod) - 1) = (/lat/)

        ;; convert to years, months, etc (0 = years, 1 = months)
        ;; note_cd_calendar requires NCL v6.0
        ti = A0&time
        ti2 = cd_calendar(ti, 0)
        ;; first: make array in years...
        tt = ti2(:, 0) + (ti2(:, 1) - 1) / 12.
        xx = dimsizes(tt)

        Tdata(imod, 0:xx(0) - 1) = (/tt/)

        ;; assign data (Y)
        ;; number of years for each model.
        tmp = ispan(floattoint(min(Tdata(imod, :))), \
                    floattoint(max(Tdata(imod, :))), 1) * 1.

        nyr = dimsizes(tmp)
        aTdata(imod, 0:nyr(0) - 1) = tmp
        delete(tmp)

        tsdim = dimsizes(A0)

        do s = 0, dimsizes(seasons) - 1  ; loop through seasons
            ;; calc trends each point.
            ;; calc seasonal time series (from TSLINE) (ann is just all months)
            ;; first need seasonal average timeseries (from tlsine)
            ;; tmpts(:, :, :) = tmpts@_FillValue
            if (dimsizes(tsdim) .eq. 3) then
                ;; conform would be a better way, this assumes y, x, t
                tmpts = new((/tsdim(1), tsdim(2), nyr/), "float")

                tmpts!0 = A0!1
                tmpts!1 = A0!2
                tmpts!2 = "year"
                tmpts&lat = A0&lat
                tmpts&lon = A0&lon
                tmpts&year = aTdata(imod, 0:nyr(0) - 1)
                tmpts@_FillValue = A0@_FillValue

                ;; array for trends
                tmptnd = tmpts(:, :, 0)

            end if

            if (dimsizes(tsdim) .eq. 2) then
                ;; conform would be a better way, this assumes y, x, t
                tmpts = new((/tsdim(1), nyr/), "float")

                tmpts!0 = A0!1
                tmpts!1 = "year"
                tmpts&lat = A0&lat
                tmpts&year = aTdata(imod, 0:nyr(0) - 1)
                tmpts@_FillValue = A0@_FillValue

                ;; array for trends
                tmptnd = tmpts(:, 0)

            end if

            ;; annual average
            ;; tmp = Ydata(imod, :)
            ;; tmp!0 = "time"
            ;; tmp&time = Xdata(imod, :)
            ;; take zonal mean of trends and assign to Y data 's'
            ;; Ydata is x, y
            ;; reform time array to years.
            tmpA0 = A0
            tmpA0&time = tt

            do y = 0, nyr - 1
                if (aTdata(imod, y) + seasdt(s, 0) .lt. max(aTdata(imod, :))) then

                    ;; average over season and puyt into x, y, year array for trends
                    if (dimsizes(tsdim) .eq. 3) then
                        tst = tmpA0({lat|:}, {lon|:}, {time | aTdata(imod, y) \
                               + seasdt(s, 0):aTdata(imod, y) + seasdt(s, 1) - 0.01})

                        n = dimsizes(tst)
                        if (n(2) .eq. nmos(s)) then
                            tmpts(:, :, y) = dim_avg_Wrap(tst)
                        end if
                    end if

                    if (dimsizes(tsdim) .eq. 2) then
                        tst = tmpA0({lat|:}, {time | aTdata(imod, y) \
                               + seasdt(s, 0):aTdata(imod, y) + seasdt(s, 1) - 0.01})

                        n = dimsizes(tst)
                        if (n(1) .eq. nmos(s)) then
                            tmpts(:, y) = dim_avg_Wrap(tst)
                        end if
                    end if

                    delete(tst)
                    delete(n)
                end if

            end do  ; end year loop

            delete(tmpA0)

            ;; now calculate trends to reduce to 2D x, y
            if (dimsizes(tsdim) .eq. 3) then
                do i = 0, tsdim(2) - 1
                    do j = 0, tsdim(1) - 1
                        rc = regline(aTdata(imod, 0:nyr(0) - 1), tmpts(j, i, :))
                        tmptnd(j, i) = rc
                    end do
                end do
                ;; now take zonal mean of trends and put into array
                A = dim_avg_Wrap(tmptnd({lat|:}, {lon|:}))

            end if

            ;; zonal mean data
            if (dimsizes(tsdim) .eq. 2) then
                do j = 0, tsdim(1) - 1
                    rc = regline(aTdata(imod, 0:nyr(0) - 1), tmpts(j, :))
                    tmptnd(j) = rc
                end do

                A = tmptnd

            end if

            Ydata(imod, 0:dimXY(imod) - 1, s) = (/A/)

            ;; add variance around latitude circle for Sdata?
            delete(A)
            delete(rc)
            delete(tmpts)
            delete(tmptnd)

        end do  ; seasonal loop
        delete(A0)
        delete(lat)
        delete(ti)
        delete(ti2)
        delete(tt)
        delete(xx)
    end do  ; imod loop

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
                    ;; case1
                    tmp = linint1_Wrap(lat1, Ydata1(season|:, lat|:), False, lat2, 0)
                    C = (/tmp(lat|:, season|:) - Ydata2/)
                    CXdata(jmod, 0:dimXY(imod) - 1) = (/lat2/)
                    CYdata(jmod, 0:dimXY(imod) - 1, :) = (/C/)
                else
                    C = (/Ydata1/)
                    ;; case2
                    tmp = linint1_Wrap(lat2, Ydata2(season|:, lat|:), False, lat1, 0)
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
            wks = gsn_open_wks(file_type, wkdir + "zonlntrend_" \
                                                + season + "_" + variable + "_c")
        else
            wks = gsn_open_wks(file_type, wkdir + "zonlntrend_" \
                                                + season + "_" + variable)
        end if

        ;; ***************************************************************
        ;; common resources for xy plots
        plot = new(2, "graphic")
        res                        = True
        res@gsnDraw                = False               ; Do not draw plot
        res@gsnFrame               = False               ; Do not advance frame
        res@trXReverse             = True
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
                res@pmLegendHeightF  = 0.35
                res@lgLabelFontHeightF = 0.016
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
        delete(res@tiMainString)
        delete(res@xyExplicitLegendLabels)       ; no label
        delete(res@pmLegendDisplayMode)          ; default is off
        delete(res@xyDashPatterns)
        delete(res@xyLineThicknesses)
        if (color_type .eq. "COLOR") then
            delete(res@xyLineColors)
            res@xyLineColors         = Ccolors
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
    delete(seasdt)
    delete(aTdata)
    delete(Tdata)
end
