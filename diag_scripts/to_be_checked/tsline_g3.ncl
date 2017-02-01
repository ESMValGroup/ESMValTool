;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# it is based on E06FIG12B.ncl by                                                #
;# Irene Cionni (irene.cionni@enea.it) & Veronika Eyring (veronika.eyring@dlr.de) #
;##################################################################################
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; Description:
; The goal is to take a variable of arbitrary dimensions, average it to 1D in
; space and 1D in time, and plot it as a line.
;
; info attributes required:
; either read data and sort by variables, or read type...

; ts_lonrange         array of pressure levels for time series plots if 2D or 3D
; ts_latrange         array of latitudes for time series plots if 2D or 3D
; ts_levrange         array of pressure levels for time series plots if 3D
; ts_ymax        (optional) y axis maximum for time series plots
; ts_ymin        (optional) y axis minimum for time series plots
; ts_xmax_year       (optional) start year for time series plots
; ts_xmin_year       (optional) end year for time series plots
; ts_anom(optional) plot anomalies.
; ts_smooth (optional) smooth...

;Note: cd_calendar requires NCL v6.0
;  A) Remove annual cycle option (anomalies) DONE: needs to be 'right' however
;     small adjustment for variable season/month length (minor)

;*******************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"

begin
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot tsline (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    ;; fields supported
    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T3M", "T2Ms", "T2Ds", "T1Ms", "T2Mz", "T0M"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    diri_out = plot_dir + "/tsline_g3/"
    file_type = "PS"

    color_type = "COLOR"

   ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    A = extract_data(0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(A))
    delete(A)

    ;; ******************************************************************
    ;; check the input parameters
    if(rank .gt. 1) then

        if(isvar("latrange")) then
            delete(latrange)
        end if
        if(isatt(info, "ts_latrange")) then
            latrange = info@ts_latrange
        else
            print("fatal:plot requires variable attribute 'ts_latrange'")
            status_exit(1)
        end if

        if(isvar("lonrange")) then
            delete(lonrange)
        end if
        if(isatt(info, "ts_lonrange")) then
            lonrange = info@ts_lonrange
        else
            if (isStrSubset(field_number, "z") .eq. False) then
                print("fatal:plot requires variable attribute 'ts_lonrange'")
                status_exit(1)
            end if
        end if

        if(field_number .eq. "T3M" .or. field_number .eq. "T2Mz") then
            if(isvar("levrange")) then
                delete(levrange)
            end if
            if(isatt(info, "ts_levrange")) then
                levrange = info@ts_levrange
            else
                print("fatal:plot requires variable attribute 'ts_levrange'")
                status_exit(1)
            end if
        end if

    end if  ; rank = 1

    multimean = False
    if(isatt(info, "ts_multimean")) then
        multimean = info@ts_multimean
    end if

    if(isatt(info, "units")) then
        UNITS = info@units
    end if

    ;; optional input parameters
    if(isatt(info, "ts_ymax")) then
        yMax_lev = info@ts_ymax
    end if
    if(isatt(info, "ts_ymin")) then
        yMin_lev = info@ts_ymin
    end if
    if(isatt(info, "ts_xmax_year")) then
        xMax_year = info@ts_xmax_year
    end if
    if(isatt(info, "ts_xmin_year")) then
        xMin_year = info@ts_xmin_year
    end if

    if(isvar("anom")) then
        delete(anom)
    end if
    if(isatt(info, "ts_anom")) then
        anom = info@ts_anom
    end if

    smooth = False
    if(isatt(info, "ts_smooth")) then
        smooth = info@ts_smooth
    end if

    ;; ******************************************************************
    ;; seasonal averaging
    ;; ******************************************************************
    seasdt = (/ (/-1, 2./), (/2, 5/), (/5, 8/), (/8, 11./) /) / 12.
    seasons = (/"DJF", "MAM", "JJA", "SON"/)
    ns = dimsizes(seasons)

    ;; same as above, but want y+s(0, s):y+s(1, s)
    ;; where  DJF        MAM     JJA    SON
    ;; 0      -11./12.   2/12.   5/12.  8/12
    ;; 1      2./12.     5/12.   8/12.  11/12
    ;; only use if there are 3 values in each.
    plot_num = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "l", "m", \
                 "n", "o", "p", "q", "r", "s", "t", "u", "v", "z"/)

    ;; find time range (if not set?)
    ;; radians for cosine of latitude weights
    rad    = 4.0 * atan(1.0) / 180.0

    ;; ==================================================================
    ;; setup data holders for plots
    ;; need to have dimensions (x & y): so need to set
    ;; A- time range  : xMin_year, xMax_year
    ;; B- value range : yMin_lev, yMax_lev
    ;; or can find them from data.
    ;; for now assume they are fixed.
    ;;not necessary? Produces a length of x axis.
    dimX = new((/dim_MOD/), "integer")
    do imod = 0, dim_MOD - 1
        A0 = extract_data(imod, data, -1, 0, 0)
        time = A0&time
        dimX(imod) = dimsizes(time)
        delete(A0)
        delete(time)
    end do
    nMax = max(dimX)
    Xdata = new((/dim_MOD, nMax/), "float")
    Ydata = new((/dim_MOD, nMax/), "float")

    Xdata!0 = "line"
    Xdata!1 = "time"

    ;; get data into arrays (X) & (Y)
    ;; loop for models
    do imod = 0, dim_MOD - 1

        ;; get data
        A0 = extract_data(imod, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        ;; convert to 1D if necessary
        if(rank .gt. 1) then

            ;; convert to T2Ms
            if (field_number .eq. "T2Ds") then
                tmp = day_to_month(A0, "T2Ds")
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if

            ;; cosine of latitude for averaging (=weights)
            coslat   = cos(A0&lat * rad)
            coslat!0 = "lat"
            coslat&lat = A0&lat

            ;; average the data over dimensions (may need more details here)
            if (field_number .ne. "T1Ms") then

                if (field_number .eq. "T3M") then
                    tmplt = dim_avg_Wrap(A0(time | :, {plev | levrange(0):levrange(1)}, \
                                            {lat | latrange(0):latrange(1)}, {lon | lonrange(0):lonrange(1)}))
                    tmpz = dim_avg_wgt_Wrap(tmplt(:, :, :), \
                                               coslat({latrange(0):latrange(1)}), 1)
                    tmp = dim_avg_Wrap(tmpz)
                    delete(tmplt)
                    delete(tmpz)
                end if
                if (field_number .eq. "T2Ms") then
                    tmplt = dim_avg_Wrap(A0(time | :, {lat | latrange(0):latrange(1)}, \
                                                     {lon | lonrange(0):lonrange(1)}))
                    tmp = dim_avg_wgt_Wrap(tmplt, coslat({latrange(0):latrange(1)}), 1)
                    delete(tmplt)
                end if
                if (field_number .eq. "T2Mz") then
                    tmpz = dim_avg_wgt_Wrap(A0(time | :, {plev | levrange(0):levrange(1)}, \
                                               {lat | latrange(0):latrange(1)}), \
                                               coslat({latrange(0):latrange(1)}), 1)
                    tmp = dim_avg_Wrap(tmpz)
                    delete(tmpz)
                end if

            end if

            if (field_number .eq. "T1Ms") then
                tmpz = dim_avg_wgt_Wrap(A0(time | :, {lat | latrange(0):latrange(1)}), \
                                             coslat({latrange(0):latrange(1)}), 1)
                tmp = tmpz
                delete(tmpz)
            end if

            delete(A0)
            A0 = tmp
            delete(tmp)

        end if  ; rank =1 loop

        ;; assign time data (X)
        ;; convert to years, months, etc (0=years, 1=months)
        ;; note_cd_calendar requires NCL v6.0
        ti = A0&time
        ti2 = cd_calendar(ti, 0)

        ;; first: make array in years...
        tt = ti2(:, 0) + (ti2(:, 1) - 1) / 12.

        Xdata(imod, 0:dimX(imod) - 1) = (/tt/)

        ;; assign data (Y)
        Ydata(imod, 0:dimX(imod) - 1) = (/A0/)

        delete(A0)
        delete(ti)
        delete(ti2)
        delete(tt)
        if rank .gt. 1 then
            delete(coslat)
        end if

    end do  ; imod

    ;; MAKE MULTIMODEL MEAN HERE.... Xdata, Ydata and dim_MOD = dim_MOD + 1
    ;; need to average appropropriate times!  Quick and dirty: take from model 0
    if (multimean .eq. True) then

        ;; shift names
        tmpn = models@name
        delete(models@name)
        models@name = new((/dim_MOD + 1/), "string")
        models@name(0) = "MEAN"
        models@name(1:dim_MOD) = tmpn

        delete(tmpn)

        tmpx = Xdata
        delete(Xdata)
        tmpy = Ydata
        delete(Ydata)

        ;; create bigger arrays and reload
        Xdata = new((/dim_MOD + 1, nMax/), "float")
        Ydata = new((/dim_MOD + 1, nMax/), "float")
        Xdata(1:dim_MOD, :) = (/tmpx(:, :)/)
        Ydata(1:dim_MOD, :) = (/tmpy(:, :)/)

        ;; copy X array from model 0
        Xdata(0, :) = tmpx(0, :)
        ;; average Y array for each time...
        do it = 0, dimX(0) - 1
            tmpx1d = ndtooned(tmpx)
            tmpy1d = ndtooned(tmpy)
            Ydata(0, it) = avg(tmpy1d(ind(tmpx1d .eq. tmpx(0, it))))
        end do

        dim_MOD = dim_MOD + 1
    end if

    ;; average seasonal and annual
    ;; sYdata, aYdata
    aYdata = new((/dim_MOD, nMax/), "float")
    aXdata = new((/dim_MOD, nMax/), "float")

    sYdata = new((/dim_MOD, nMax, ns/), "float")

    ;; trend lines
    aYreg = new((/dim_MOD, nMax/), "float")
    sYreg = new((/dim_MOD, nMax, ns/), "float")
    ;; slope
    arc = new((/dim_MOD/), "float")
    src = new((/dim_MOD, ns/), "float")
    ;; significance level (diff from zero w/ t - tst)
    asig = new((/dim_MOD/), "float")
    ssig = new((/dim_MOD, ns/), "float")

    ;; loop for models
    do imod = 0, dim_MOD - 1
        ;; x array (one value per year)
        tmp = ispan(floattoint(min(Xdata(imod, :))), \
                    floattoint(max(Xdata(imod, :))), 1) * 1.
        nyr = dimsizes(tmp)
        aXdata(imod, 0:nyr(0) - 1) = tmp
        delete(tmp)
        ;; annual average
        tmp = Ydata(imod, :)
        tmp!0 = "time"
        tmp&time = Xdata(imod, :)
        do y = 0, nyr - 1

            ;; should use 'month_to_annual_weighted' function...
            n = 0

            ttt = ind(floattoint(tmp&time) .eq. floattoint(aXdata(imod, y)))
            aa = tmp(ttt)

            delete(ttt)

            n = dimsizes(aa)

            if (n(0) .gt. 9.) then
                ;; aYdata(imod, y) = avg(tmp({aXdata(imod, y):aXdata(imod, y) + 0.99}))
                aYdata(imod, y) = avg(aa)

            end if
            delete(aa)
            delete(n)

            ;; loop for seasons
            ;; same as above, but want y+s(0, s):y+s(1, s)
            ;; where  DJF        MAM     JJA    SON
            ;; 0      -11./12.   2/12.   5/12.  8/12
            ;; 1      2./12.     5/12.   8/12.  11/12
            ;; only use if there are 3 values in each.
            do s = 0, ns - 1
                if (aXdata(imod, y) + seasdt(s, 0) .lt. max(aXdata(imod, :))) then
                    n = 0

                    ti = tmp&time
                    tr = aXdata(imod, y)

                    ttt = ind(ti .ge. tr + seasdt(s, 0) .and. ti .lt. tr + seasdt(s, 1))

                    ssz = dimsizes(ttt)
                    if (ssz(0) .gt. 1) then
                        aa = tmp(ttt)

                        n = dimsizes(aa)

                        if (n(0) .ge. 2) then
                            sYdata(imod, y, s) = avg(aa)
                        end if
                        delete(aa)
                        delete(n)
                    end if
                    delete(ti)
                    delete(tr)
                    delete(ttt)
                end if
            end do  ; seasonal loop

        end do  ; annual loop

        delete(tmp)

        ;; remove anomalies
        if (anom .eq. True) then

            ;; annual
            tmp = avg(aYdata(imod, :))
            aYdata(imod, :) = aYdata(imod, :) - tmp
            delete(tmp)

            ;; seasonal
            do s = 0, ns - 1
                tmp = avg(sYdata(imod, :, s))
                sYdata(imod, :, s) = sYdata(imod, :, s) - tmp
                delete(tmp)
            end do
        end if

        ;; calculate trends here using regline (seasonal and annual).
        ;; save lines for plotting (add line). may also want data.
        ;; annual
        rc = regline(aXdata(imod, :), aYdata(imod, :))
        df   = rc@nptxy - 2
        aYreg(imod, :) = rc * aXdata(imod, :) + rc@yintercept
        arc(imod) = rc
        asig(imod) = (1. - betainc(df/(df + rc@tval^2), df / 2.0, 0.5)) * 100.

        delete(rc)
        delete(df)

        ;; seasonal
        do s = 0, ns - 1
            rc = regline(aXdata(imod, :), sYdata(imod, :, s))
            df   = rc@nptxy - 2
            sYreg(imod, :, s) = rc * aXdata(imod, :) + rc@yintercept
            src(imod, s) = rc
            ssig(imod, s) = (1. - betainc(df/(df + rc@tval^2), df / 2.0, 0.5)) * 100.
            delete(rc)
            delete(df)
        end do  ; seasonal lines

        ;; write trends and significance
        print(models@name(imod) + "  ANN  Trend (" + info@units + "/yr)=" \
                                + arc(imod) + "  Significance (%)=" + asig(imod))

    end do     ; model loop

    ;; ================================================
    ;; METRICS
    ;; Make sure time period is the same!
    ;; ================================================
    ;; ******************************************************************
    ;; Determine which models are 'observations'
    obsi = new((/dim_MOD/), "integer")
    obsi(:) = 0
    do imod = 0, dim_MOD - 1
        if (models@name(imod) .eq. "ERAI") then
            obsi(imod) = 1
        end if
        if (models@name(imod) .eq. "ERA40") then
            obsi(imod) = 1
        end if
        if (models@name(imod) .eq. "NCEP") then
            obsi(imod) = 1
        end if
        if (models@name(imod) .eq. "NCEP2") then
            obsi(imod) = 1
        end if
        if (models@name(imod) .eq. "JRA25") then
            obsi(imod) = 1
        end if
        if (models@name(imod) .eq. "RICH-ERA40") then
            obsi(imod) = 1
        end if
        if (models@name(imod) .eq. "HALOE") then
            obsi(imod) = 1
        end if
    end do
    nobs = num(obsi .eq. 1)
    obsii = ind(obsi .eq. 1)
    ;******************************************************************

    if (isatt(info, "refModel")) then

        ;; logic:  test 3 metrics.
        ;; set up arrays
        gm = new((/dim_MOD/), "float")
        gc = new((/dim_MOD/), "float")
        gv = new((/dim_MOD/), "float")
        gsum = new((/dim_MOD/), "float")
        gradetx = new((/dim_MOD + 1/), "string")

        ng = 3.

        ;; select 'obs' (as a model)
        ii = -99
        do imod = 0, dim_MOD - 1
            if (models@name(imod) .eq. info@refModel) then
                ii = imod
            end if
        end do

        ;; loop for models (also calculate for 'obs' should be perfect)
        ref_avg = avg(aYdata(ii, :))
        ref_anom = aYdata(ii, :) - ref_avg
        ref_sd = stddev(ref_anom)

        ;; ******************************************************************
        ;; multi - obs averages (new for _g3)
        ref_avg3 = ref_avg
        ref_anom3 = ref_anom
        ref_sd3 = ref_sd

        ;; warning: requires the same years!
        if(nobs .gt. 1) then
            ref_avg3 = avg(aYdata(obsii, :))
            ref_sd3 = sqrt(sum((dim_avg(aYdata(obsii, :)) - ref_avg3)^2.) / (nobs - 1.))
            do y = 0, nMax - 1
                ref_anom3(y) = avg(aYdata(obsii, y) - avg(ref_avg3))
            end do
        end if

        ref_avg = ref_avg3
        ref_anom = ref_anom3
        ref_sd = ref_sd3
;******************************************************************

        do imod = 0, dim_MOD - 1

            ;; remove mean
            tmp = avg(aYdata(imod, :))
            ts = aYdata(imod, :) - tmp
            delete(tmp)

            ;; A: 'mean' following waugh & eyring
            ;; or select an acceptable 'delta'
            if(variable .eq. "ptp") then
                dd = 10.
            else
                dd = ng
            end if

            gm(imod) = 1. - abs(avg(aYdata(imod, :)) - ref_avg) / ref_sd / dd

            if (gm(imod) .lt. 0.) then
                gm(imod) = 0.
            end if

            ;;  B: "correlation"
            ;;  need to be able to match years!
            ;;  gc = (corr(anom_mod, anom_obs) + 1) / 2
            gg = ind(.not. ismissing(ref_anom))
            gc(imod) = (esccr(ts(gg), ref_anom(gg), 0) + 1.) / 2.

            ;;  C: "variance"
            ;;  gv = 1 - 1/ng * (sigma(mod)-sigma(obs))/sigma(obs)
            gv(imod) = 1. - abs(stddev(ts) - ref_sd) / ref_sd / ng
            if (gv(imod) .lt. 0.) then
                gv(imod) = 0.
            end if

;SUM
;       if(variable .eq. "ptp") then
;               gsum(imod) = gm(imod)/4. + gc(imod)/2. + gv(imod)/4.
;       else
            gsum(imod) = (gm(imod) + gc(imod) + gv(imod)) / 3.
;       end if

            ;; Text
            gradetx(imod) = sprintf("%5.2f", gm(imod)) + " " + sprintf("%5.2f", gc(imod)) \
                                     + " " + sprintf("%5.2f", gv(imod)) \
                                     + " " + sprintf("%5.2f", gsum(imod)) \
                                     + " " + models@name(imod)
        end do
        gradetx(dim_MOD) = " GM   " + " GC   " + " GV   " + " GSUM " + " Model"
    end if

    ;; end metrics
    ;; ================================================

    ;; -------------------------------------------
    ;; SMOOTHING
    ;; -------------------------------------------
    if (smooth .eq. True) then
        ;; smooth data if desired with a 5 year running mean....
        wgt = (/0.2, 0.2, 0.2, 0.2, 0.2/)
        tmp = wgt_runave_Wrap(aYdata, wgt, 0)
        aYdata = tmp
        delete(tmp)

        ;; smooth montly data
        do s = 0, ns - 1
            tmp = wgt_runave_Wrap(sYdata(:, :, s), wgt, 0)
            sYdata(:, :, s) = tmp
            delete(tmp)
        end do
    end if
;-------------------------------------------

    ;; ******************************************************************
    ;; The function CCMVal_Color and CCMVal_Dash  in misc_function.ncl
    ;; associate at any model a color and dash pattern
    ;; ******************************************************************
    ;; ******************************************************************
    colors = CCMVal2_ColorTable(models@name)
    dashs = CCMVal2_DashTable(models@name)
    lgLabels = models@name  ; + " " + models@case_name
    TLineThick = new((/dim_MOD/), integer)
    TLineThick(:) = 1
    LineThick = CCMVal_Thick(models@name) + 1

    Ccolors = new((/dim_MOD - 1/), string)
    Cdashs = new((/dim_MOD - 1/), integer)
    CLineThick = new((/dim_MOD - 1/), integer)
    CLineThick(:) = 4

    ;; make string for region
    rgstr = info@long_name
    if(isvar("lonrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", lonrange(0)) + "-" \
                      + sprintf("%3.0f", lonrange(1)) + "lon"
    end if
    if(isvar("latrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", latrange(0)) + "-" \
                      + sprintf("%3.0f", latrange(1)) + "lat"
    end if
    if(isvar("levrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", levrange(0)) + "-" \
                      + sprintf("%3.0f", levrange(1)) + "hPa"
    end if

    ;; ***************************************************************
    ;; --------------------- Open files for plots --------------------
    if(anom .eq. True) then
        wks = gsn_open_wks(file_type, diri_out + "tsline_g3_" + variable + "_ANN_anom")
    else
        wks = gsn_open_wks(file_type, diri_out + "tsline_g3_" + variable + "_ANN_c")
    end if

    ;; ***************************************************************
    ;;  common resources for xy plots
    plot = new(2, "graphic")
    res                        = True
    res@gsnDraw                = False               ; Do not draw plot
    res@gsnFrame               = False               ; Do not advance frame
    if(info@units .eq. "hPa") then
        res@trYReverse             = True
    else
        res@trYReverse             = False
    end if
    res@tiXAxisString          = "Year"
    res@tiXAxisFontHeightF     = 0.020
    res@tiYAxisFontHeightF     = 0.020
    res@tiMainFontHeightF      = 0.025

    res@vpHeightF       = 0.4                    ; change aspect ratio of plot
    res@vpWidthF        = 0.8

    res@pmLegendSide           = "Right"            ; Change location of

    res@pmLegendOrthogonalPosF = 0.01               ; of the legend
    res@pmLegendWidthF         = 0.2               ; Change width and
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
            res@pmLegendHeightF  = 0.5
            res@lgLabelFontHeightF = 0.016
            res@pmLegendParallelPosF = 0.50
        end if
    end if
    res@txFontHeightF          = 0.015
    res@xyDashPatterns         = dashs
    res@xyLineThicknesses      = LineThick
    if (color_type .eq. "COLOR") then
        res@xyLineColors        = colors
    end if

    res@tiMainString           = "ANN " + rgstr
    res@pmLegendDisplayMode    = "Always"
    res@tiYAxisString          = info@long_name + " (" + info@units + ")"
    res@xyExplicitLegendLabels = lgLabels
    res@xyMonoLineColor        = False

    ;; ************************************************************
    ;; XAxis
    if (isvar("xMin_year") .eq. True) then
        res@trXMinF       = xMin_year
    else
        res@trXMinF       = min(aXdata)
    end if
    if (isvar("xMax_year") .eq. True) then
        res@trXMaxF       = xMax_year
    else
        res@trXMaxF       = max(aXdata)
    end if

    ;; ************************************************************
    ;; YAxis
    res@trYMinF       = min(aYdata)
    res@trYMaxF       = max(aYdata)

    plot(0) = gsn_csm_xy(wks, aXdata, aYdata, res)   ; cases line plot

    ;; add trend lines
    res@xyLineThicknesses      = TLineThick
    do imod = 0, dim_MOD - 1
        add_line(wks, plot(0), aXdata(imod, :), aYreg(imod, :), colors(imod), dashs(imod))
    end do

    ;; ----------------------------------------------
    ;; common resources for panel plots
    pan = True
    pan@gsnMaximize = True
    pan@gsnPaperOrientation = "portrait"
    pan@gsnFrame = False
    gsn_panel(wks, plot, (/2, 1/), pan)         ; panel plot
    frame(wks)

    delete(plot)

    ;; seasonal loop here...
    plot = new(4, "graphic")

    if(anom .eq. True) then
        wks = gsn_open_wks(file_type, diri_out + "tsline_g3_" + variable + "_SEAS_anom")
    else
        wks = gsn_open_wks(file_type, diri_out + "tsline_g3_" + variable + "_SEAS_c")
    end if

    do s = 0, ns - 1
        ;; ***************************************************************
        ;; --------------------- Open files for plots --------------------
        res@tiMainString           = seasons(s) + rgstr
        res@xyLineThicknesses      = LineThick

        ;; ************************************************************
        ;; YAxis
        res@trYMinF       = min(sYdata(:, :, s))
        res@trYMaxF       = max(sYdata(:, :, s))

        plot(s) = gsn_csm_xy(wks, aXdata, sYdata(:, :, s), res)

        ;; add trend lines
        res@xyLineThicknesses      = TLineThick
        do imod = 0, dim_MOD - 1
            add_line(wks, plot(s), aXdata(imod, :), sYreg(imod, :, s), colors(imod), dashs(imod))
        end do

    end do  ; seasonal plots

    ;; ----------------------------------------------
    ;; common resources for panel plots
    pan = True
    pan@gsnMaximize = True
    pan@gsnPaperOrientation = "portrait"
    pan@gsnFrame = False
    gsn_panel(wks, plot, (/2, 2/), pan)         ; panel plot
    frame(wks)

    delete(Ydata)
    delete(Xdata)
    delete(res)
    delete(plot)
    delete(pan)
    delete(seasdt)
    delete(seasons)
end
