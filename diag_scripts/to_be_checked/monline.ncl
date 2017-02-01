;################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified#
;# in the CCMVal-2 data request available at                                    #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf   #
;# it is based on tsline.ncl by                                                 #
;# andrew gettelman                                                             #
;# Author: Andrew gettelman, august 2009                                        #
;################################################################################
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;
; Description:
; The goal is to take a variable of arbitrary dimensions, average it to 0D in
; space and 1D in time, then find monthly means and plot.
;
; info attributes required:
; either read data and sort by variables, or read type...

; ml_lonrange         array of longitudes for time series plots if 2D or 3D
; ml_latrange         array of latitudes for time series plots if 2D or 3D
; ml_levrange         array of pressure levels for time series plots if 3D
; ml_ymax        (optional) y axis maximum for time series plots
; ml_ymin        (optional) y axis minimum for time series plots
; ml_xmax_year       (optional) start year for time series plots
; ml_xmin_year       (optional) end year for time series plots
; ml_anom(optional) plot anomalies.

;Note: cd_calendar requires NCL v6.0
;  A) Remove annual cycle option (anomalies) DONE: needs to be 'right' however
;     small adjustment for variable season/month length (minor)


;*******************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"

begin

    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot monline (var: " + variable + ")")
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

    diri_out = plot_dir + "/monline/"
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
        if(isatt(info, "ml_latrange")) then
            latrange = info@ml_latrange
        else
            print("fatal:plot requires variable attribute 'ml_latrange'")
            status_exit(1)
        end if

        if(isvar("lonrange")) then
            delete(lonrange)
        end if
        if(isatt(info, "ml_lonrange")) then
            lonrange = info@ml_lonrange
        else
            if (isStrSubset(field_number, "z") .eq. False) then
                print("fatal:plot requires variable attribute 'ml_lonrange'")
                status_exit(1)
            end if
        end if

        if(field_number .eq. "T3M" .or. field_number .eq. "T2Mz") then
            if(isvar("levrange")) then
                delete(levrange)
            end if
            if(isatt(info, "ml_levrange")) then
                levrange = info@ml_levrange
            else
                print("fatal:plot requires variable attribute 'ml_levrange'")
                status_exit(1)
            end if
        end if

    end if  ; rank = 1

    if(isatt(info, "units")) then
        UNITS = info@units
    end if

    ;;  optional input parameters
    if(isatt(info, "ml_ymax")) then
        yMax_lev = info@ml_ymax
    end if
    if(isatt(info, "ml_ymin")) then
        yMin_lev = info@ml_ymin
    end if
    if(isatt(info, "ml_xmax_year")) then
        xMax_year = info@ml_xmax_year
    end if
    if(isatt(info, "ml_xmin_year")) then
        xMin_year = info@ml_xmin_year
    end if

    if(isvar("anom")) then
        delete(anom)
    end if
    if(isatt(info, "ml_anom")) then
        anom = info@ml_anom
    else
        anom = False
    end if

    ;;******************************************************************
    ;;monthly averaging
    ;;******************************************************************

    ;;monames = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",  \
    ;;            "Sep", "Oct", "Nov", "Dec"/)

    monames = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"/)
    nmos = dimsizes(monames)

    moarr = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)

    plot_num = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "l", "m", "n", \
                 "o", "p", "q", "r", "s", "t", "u", "v", "z"/)

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
    ;; not necessary? Produces a length of x axis.
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

    ;;  get data into arrays (X) & (Y)
    ;;loop for models
    do imod = 0, dim_MOD - 1

        ;;get data
        A0 = extract_data(imod, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        ;;convert to 1D if necessary
        if(rank .gt. 1) then

            ;; convert to T2Ms
            if (field_number .eq. "T2Ds") then
                tmp = day_to_month(A0, "T2Ds")
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if

            ;; cosine of latitude for averaging (=weights)
            coslat = cos(A0&lat * rad)
            coslat!0 = "lat"
            coslat&lat = A0&lat

            ;; average the data over dimensions (may need more details here)
            if (field_number .ne. "T1Ms") then

                if (field_number .eq. "T3M") then
                    tmplt = dim_avg_Wrap(A0(time | :, {plev | levrange(0):levrange(1)}, \
                                          {lat | latrange(0):latrange(1)}, {lon | lonrange(0):lonrange(1)}))
                    tmpz = dim_avg_wgt_Wrap(tmplt(:, :, :), coslat({latrange(0):latrange(1)}), 1)
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

                delete(A0)
                A0 = tmp
                delete(tmp)
            end if
        end if  ; rank =1 loop

        ;; assign time data (X)
        ;; convert to years, months, etc (0=years, 1=months)
        ;; note_cd_calendar requires NCL v6.0
        ti = A0&time
        ti2 = cd_calendar(ti, 0)

        ;; first: make X array in months
        tt = ti2(:, 1)

        ;; hack for waccm (old) and CAM3.5 data
        ;; CAREFUL: Corrected in read data!
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

    ;; average monthly
    mYdata = new((/dim_MOD, 12/), "float")

    mYdata!1 = "month"
    mYdata&month = moarr

    ;; standard deviation
    mSdata = new((/dim_MOD, 12/), "float")
    mSdata!1 = "month"
    mSdata&month = moarr

    mNdata = new((/dim_MOD, 12/), "float")
    mNdata!1 = "month"
    mNdata&month = moarr

    ;; loop for models
    do imod = 0, dim_MOD - 1
        do mi = 0, nmos - 1
            tmp = ind(Xdata(imod, :) .eq. mi + 1 .and. (.not. ismissing(Ydata(imod, :))))
            mYdata(imod, mi) = avg(Ydata(imod, tmp))
            mSdata(imod, mi) = stddev(Ydata(imod, tmp))
            mNdata(imod, mi) = num(.not. ismissing(Ydata(imod, tmp)))
            delete(tmp)
        end do
        ;; remove mean
        if (anom .eq. True) then
            mYdata(imod, :) = mYdata(imod, :) - avg(mYdata(imod, :))
        end if
    end do

    ;; ADD HALOE Data for H2O (from E06FIG07.ncl)
    if (isatt(info, "ml_obs_name") .and. isatt(info, "ml_obs_file") \
        .and. isatt(info, "ml_obs_stdv_file")) then
        f_obs = addfile(info@ml_obs_file, "r")
        ;; Observation by HALOE are already processed:mean
        v = f_obs->$variable$
        delete(f_obs)
        f_obs = addfile(info@ml_obs_stdv_file, "r")
        vName = variable + "@STDEV"

        ;; Observation by HALOE  are already processed:standard deviation
        v_STD = f_obs->$vName$
        delete(f_obs)
        v!0 = "month"
        v!1 = "plev"
        v!2 = "lat"
        v_STD!0 = "month"
        v_STD!1 = "plev"
        v_STD!2 = "lat"
        if ((v@units .eq. "m^3/m^3") .and. (variable .eq. "HCl")) then
            v = v * 10^9
            v_STD = v_STD * 10^9
            v@units = "ppbv"
            v_STD@units = "ppbv"
        end if
        if (v@units .eq. "m^3/m^3") then
            v = v * 10^6
            v_STD = v_STD * 10^6
            v@units = "ppm"
            v_STD@units = "ppm"
        end if

        ;; average to 1-D from T2Mz
        ;; cosine of latitude for averaging (=weights)
        coslat   = cos(v&lat * rad)
        coslat!0 = "lat"
        coslat&lat = v&lat

        tmpz = dim_avg_wgt_Wrap(v(month | :, {plev | levrange(0):levrange(1)}, \
        {lat | latrange(0):latrange(1)}), coslat({latrange(0):latrange(1)}), 1)
        tmp = dim_avg_Wrap(tmpz)
        delete(tmpz)
        delete(v)
        v = tmp
        delete(tmp)

        tmpz = dim_avg_wgt_Wrap(v_STD(month | :, {plev | levrange(0):levrange(1)}, \
        {lat | latrange(0):latrange(1)}), coslat({latrange(0):latrange(1)}), 1)
        tmp = dim_avg_Wrap(tmpz)
        delete(tmpz)
        delete(v_STD)
        v_STD = tmp
        delete(tmp)

        v_N = (/7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7./)

        delete(coslat)

        ;; remove anomalies
        if (anom .eq. True) then
            v(:) = v(:) - avg(v(:))
        end if

        ;; add to arrays
        mod_names = new((/dim_MOD + 1/), "string")
        mod_names(0) = info@ml_obs_name
        mod_names(1:dim_MOD) = models@name

        tmp = new((/dim_MOD + 1, 12/), "float")
        copy_VarMeta(mYdata, tmp)
        tmp(0, :) = v
        tmp(1:dim_MOD, :) = (/mYdata(:, :)/)
        delete(mYdata)
        mYdata = tmp
        delete(tmp)

        tmp = new((/dim_MOD + 1, 12/), "float")
        copy_VarMeta(mSdata, tmp)
        tmp(0, :) = v_STD
        tmp(1:dim_MOD, :) = (/mSdata(:, :)/)
        delete(mSdata)
        mSdata = tmp
        delete(tmp)

        tmp = new((/dim_MOD + 1, 12/), "float")
        copy_VarMeta(mNdata, tmp)
        tmp(0, :) = (/v_N/)
        tmp(1:dim_MOD, :) = (/mNdata(:, :)/)
        delete(mNdata)
        mNdata = tmp
        delete(tmp)

    else
        mod_names = models@name
    end if

    ;; make string for region
    rgstr = info@long_name
    if(isvar("lonrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", lonrange(0)) + "-" + sprintf("%3.0f", lonrange(1)) + "lon"
    end if
    if(isvar("latrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", latrange(0)) + "-" + sprintf("%3.0f", latrange(1)) + "lat"
    end if
    if(isvar("levrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", levrange(0)) + "-" + sprintf("%3.0f", levrange(1)) + "hPa"
    end if

    ;; write out array to netcdf file.
    mYdata!0 = "models"
    mYdata@info = rgstr

    mSdata!0 = "models"
    mSdata@info = rgstr

    mNdata!0 = "models"
    mNdata@info = rgstr

    mod_names!0 = "models"

    fout = diri_out + "monline_" + variable + ".nc"

    system("/bin/rm -f " + fout)   ; remove any pre-existing file

    ncdf = addfile(fout, "c")  ; open output netCDF file

    ;;  output variables directly
    ovarn = variable + "_ann_cycle"
    ncdf->$ovarn$ = mYdata
    osvarn = variable + "_ann_cycle_stdv"
    ncdf->$osvarn$ = mSdata
    onvarn = variable + "_nmonths"
    ncdf->$onvarn$ = mNdata

    nn = num(.not. ismissing(mod_names))

    ;; *****************************************************************
    ;;  create global attributes of the file
    ;; *****************************************************************
    setfileoption(ncdf, "DefineMode", True)
    fAtt              = True            ; assign file attributes
    fAtt@creation_date= systemfunc("date")
    do im = 0, nn - 1
        str = unique_string("model_")
        fAtt@$str$= mod_names(im)
    end do
                                ; copy attributes from the source file
    fileattdef(ncdf, fAtt)            ; copy file attributes
    setfileoption(ncdf, "DefineMode", False)

;      ncdf->model_names = mod_names

;================================================
;METRICS
;Make sure time period is the same!
;================================================

    if (isatt(info, "ml_stdv_plot")) then
        modn = num(.not. ismissing(mod_names))

        ;; logic:  test 3 metrics.
        ;; set up arrays
        gm = new((/modn/), "float")
        gc = new((/modn/), "float")
        gv = new((/modn/), "float")
        gsum = new((/modn/), "float")
        gradetx = new((/modn + 1/), "string")

        ng = 3.

        ;; select 'obs' (as a model)
        ii = -99
        do imod = 0, modn - 1
            if (mod_names(imod) .eq. info@ml_stdv_plot) then
                ii = imod
            end if
        end do

        print(models@name(ii))
        ;; loop for models (also calculate for 'obs' should be perfect)

        ref_avg = avg(mYdata(ii, :))
        ref_anom = mYdata(ii, :) - ref_avg

        do imod = 0, modn - 1

            ;; remove mean
            tmp = avg(mYdata(imod, :))
            ts = mYdata(imod, :) - tmp
            delete(tmp)

            ;; A: 'mean' following waugh & eyring
            ;; or select an acceptable 'delta'
            if(variable .eq. "ptp") then
                dd = 10.
            else
                dd = ng
            end if

            gm(imod) = 1. - abs(avg(mYdata(imod, :)) - ref_avg) / stddev(ref_anom) / dd

            if (gm(imod) .lt. 0.) then
                gm(imod) = 0.
            end if

            ;; B: "correlation"
            ;; need to be able to match years!
            gg = ind(.not. ismissing(ref_anom))
            gc(imod) = (esccr(ts(gg), ref_anom(gg), 0) + 1.) / 2.

            ;; C: "variance"
            ;;  gv = 1 - 1/ng * (sigma(mod)-sigma(obs))/sigma(obs)
            gv(imod)= 1. - abs((stddev(ts) - stddev(ref_anom))) / stddev(ref_anom) / ng
            if (gv(imod) .lt. 0.) then
                gv(imod) = 0.
            end if

            ;; SUM
            gsum(imod) = (gm(imod) + gc(imod) + gv(imod)) / 3.

            ;; Text
            gradetx(imod) = sprintf("%5.2f", gm(imod)) + " " + sprintf("%5.2f", gc(imod))+ \
            " " + sprintf("%5.2f", gv(imod)) + " " + sprintf("%5.2f", gsum(imod)) + " " + mod_names(imod)

        end do

        gradetx(modn) = " GM   " + " GC   " + " GV   " + " GSUM " + " Model"

        print(gradetx)
    end if
    ;; end metrics
    ;; ================================================

    ;; ******************************************************************
    ;; The function CCMVal_Color and CCMVal_Dash  in misc_function.ncl
    ;;  associate at any model a color and dash pattern
    ;; ******************************************************************
    ;; ******************************************************************
    colors = CCMVal2_ColorTable(mod_names)
    dashs = CCMVal2_DashTable(mod_names)
    lgLabels = mod_names
    nn = num(.not. ismissing(mod_names))
    LineThick= CCMVal_Thick(mod_names)

    Ccolors = new((/nn - 1/), string)
    Cdashs = new((/nn - 1/), integer)
    CLineThick = new((/nn - 1/), integer)
    CLineThick(:) = 2

    ;; ***************************************************************
    ;; --------------------- Open files for plots --------------------
    if(anom .eq. True) then
        wks = gsn_open_wks(file_type, diri_out + "monline_" + variable + "_ANN_anom")
    else
        wks = gsn_open_wks(file_type, diri_out + "monline_" + variable + "_ANN_c")
    end if

    ;; ***************************************************************
    ;;  common resources for xy plots
    plot = new(2, "graphic")
    res                        = True
    res@gsnDraw                = False               ; Do not draw plot
    res@gsnFrame               = False               ; Do not advance frame
    res@trXReverse             = False
    res@tiXAxisString          = "Month"
    res@tiXAxisFontHeightF     = 0.020
    res@tiYAxisFontHeightF     = 0.020
    res@tiMainFontHeightF      = 0.025

    res@pmLegendSide           = "Right"            ; Change location of

    res@pmLegendOrthogonalPosF = 0.02               ; of the legend
    res@pmLegendWidthF         = 0.15               ; Change width and
    res@lgPerimOn              = True
    if (dim_MOD .le. 4) then
        res@pmLegendHeightF      = 0.10               ; height of legend.
        res@lgLabelFontHeightF   = .020
        res@pmLegendParallelPosF = 0.10               ; of the legend
    else
        if (dim_MOD .le. 8) then
            res@pmLegendHeightF      = 0.18
            res@lgLabelFontHeightF   = 0.016
            res@pmLegendParallelPosF = 0.18
        else
            res@pmLegendHeightF      = 0.5
            res@lgLabelFontHeightF   = 0.012
            res@pmLegendParallelPosF = 0.50
        end if
    end if
    res@txFontHeightF          = 0.015
    res@xyDashPatterns         = dashs
    res@xyLineThicknesses      = LineThick
    if (color_type .eq. "COLOR") then
        res@xyLineColors        = colors
    end if

    res@tiMainString           = rgstr
    res@pmLegendDisplayMode    = "Always"
    res@tiYAxisString          = info@long_name + " (" + info@units + ")"
    res@xyExplicitLegendLabels = lgLabels
    res@xyMonoLineColor        = False

    ;; ************************************************************
    ;; XAxis
    Xarr = moarr

    res@vpHeightF = 0.4                    ; change aspect ratio of plot
    res@vpWidthF  = 0.8

    res@trXMaxF  = 13
    res@trXMinF  = 0

    res@tmXBMode  = "Explicit"
    res@tmXBValues = Xarr
    res@tmXBLabels = monames

    ;; ************************************************************
    ;; YAxis
    res@trYMinF       = min(mYdata)
    res@trYMaxF       = max(mYdata)

    ;; ************************************************************
    ;; PLOT
    res@tfPolyDrawOrder    = "Predraw"               ; put line on top

    plot(0) = gsn_csm_xy(wks, Xarr, mYdata, res)   ; cases line plot

    ;; ************************************************************
    ;; plot standard deviation
    if (isatt(info, "ml_stdv_plot"))

        ii = -99
        do imod = 0, dim_MOD - 1
            if (mod_names(imod) .eq. info@ml_stdv_plot) then
                ii = imod
            end if
        end do

        if(ii .lt. 0) then
            print("fatal:ml_stdv_plot model (" + info@ml_stdv_plot + ") not found, exiting")
            status_exit(1)
        end if

        ;;DEFINE A POLYGON BASED on selected 'model' (0=obs)
        xp = new((/2*nmos/), float)
        yp = new((/2*nmos/), float)
        do k = 0, nmos - 1
            yp(k) = mYdata(ii, k) + 2. * mSdata(ii, k)
            xp(k) = mYdata&month(k)
            xp(2 * nmos - 1-k) = mYdata&month(k)
            yp(2 * nmos - 1-k) = mYdata(ii, k) - 2. * mSdata(ii, k)
        end do

        gsres                   = True                        ; poly res
        gsres@gsFillColor       = "Gray"                 ; color chosen
        dummy = gsn_add_polygon(wks, plot(0), xp, yp, gsres)
    end if

    ;;----------------------------------------------
    ;; common resources for panel plots
    pan = True
    pan@gsnMaximize = True
    pan@gsnPaperOrientation = "portrait"
    pan@gsnFrame = False
    gsn_panel(wks, plot, (/2, 1/), pan)         ; panel plot
    frame(wks)

    delete(mYdata)
    delete(mSdata)
    delete(mNdata)
    delete(Ydata)
    delete(Xdata)
    delete(res)
    delete(plot)
    delete(pan)
end
;#################################################################################
