;##############################################################################
;# PLOTS MONTHLY MEAN TIME SERIES
;# Author: Andrew Gettelmann (NCAR, USA)
;# ESMVal project
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    * The goal is to take a variable of arbitrary dimensions, average it    #
;#      to 0D in space and 1D in time, then find monthly means and plot.      #
;#    * This script is based on tsline.ncl                                    #
;#    * The CCMVal script script has been furthered developed as a part of    #
;#      the ESMValTool                                                        #
;#                                                                            #
;# Required variable_info attributes (variable specific):                     #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    * ml_lonrange: array of longitudes for time series plots if 2D or 3D    #
;#    * ml_latrange: array of latitudes for time series plots if 2D or 3D     #
;#    * ml_levrange: array of pressure levels for time series plots if 3D     #
;#                                                                            #
;# Optional diag_script_info attributes (diagnostic specific)                 #
;#    * ml_ymax: y axis maximum for time series plots                         #
;#    * ml_ymin: y axis minimum for time series plots                         #
;#    * ml_xmax_year: start year for time series plots                        #
;#    * ml_xmin_year: end year for time series plots                          #
;#    * ml_anom: plot anomalies                                               #
;#                                                                            #
;# Caveats                                                                    #
;#    * cd_calendar requires NCL v6.0                                         #
;#    * remove annual cycle option (anomalies) DONE: needs to be 'right'      #
;#      however small adjustment for variable season/month length (minor)     #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/calc_seasonal.ncl"
;load "./diag_scripts/CCMVal_FUNCTIONS/misc_function.ncl"
load "./diag_scripts/lib/ncl/misc_function.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/uniform_coordinates.ncl"

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot SET8 - " + diag_script + " (var: " + \
                variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do
end

; Load var_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;; fields supported
    plot_fields = (/"T3M", "T2Ms", "T2Ds", "T1Ms", "T2Mz", "T0M"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field " + field_type0)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; Output dir
    ;; 'plot_dir' if fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; 'output_file_type' if fetched from ncl.interface
    if (isdefined("output_file_type")) then
        file_type = output_file_type
    else
        file_type = "PS"
    end if

    file_type = str_lower(file_type)

    color_type = "COLOR"

    A = extract_data(0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(A))
    delete(A)

    ;; ******************************************************************
    ;; check the input parameters
    if(rank .gt. 1) then

        if(isvar("latrange")) then
            delete(latrange)
        end if
        if(isatt(diag_script_info, "ml_latrange")) then
            latrange = diag_script_info@ml_latrange
        else
            error_msg("f", diag_script, "", "plot requires " + \
                      "diag_script_info attribute ml_latrange")
        end if

        if(isvar("lonrange")) then
            delete(lonrange)
        end if
        if(isatt(diag_script_info, "ml_lonrange")) then
            lonrange = diag_script_info@ml_lonrange
        else
            if (isStrSubset(field_type0, "z") .eq. False) then
                error_msg("f", diag_script, "", "plot requires " + \
                          "diag_script_info attribute ml_lonrange")
            end if
        end if

        if(field_type0 .eq. "T3M" .or. field_type0 .eq. "T2Mz") then
            if(isvar("levrange")) then
                delete(levrange)
            end if
            if(isatt(diag_script_info, "ml_levrange")) then
                levrange = diag_script_info@ml_levrange
            else
                error_msg("f", diag_script, "", "plot requires " + \
                          "diag_script_info attribute ml_levrange")
            end if
        end if

    end if  ; rank = 1

    if(isatt(variable_info, "units")) then
        UNITS = variable_info@units
    end if

    ;;  optional input parameters
    if(isatt(diag_script_info, "ml_ymax")) then
        yMax_lev = diag_script_info@ml_ymax
    end if
    if(isatt(diag_script_info, "ml_ymin")) then
        yMin_lev = diag_script_info@ml_ymin
    end if
    if(isatt(diag_script_info, "ml_xmax_year")) then
        xMax_year = diag_script_info@ml_xmax_year
    end if
    if(isatt(diag_script_info, "ml_xmin_year")) then
        xMin_year = diag_script_info@ml_xmin_year
    end if

    if(isvar("anom")) then
        delete(anom)
    end if
    if(isatt(diag_script_info, "ml_anom")) then
        anom = diag_script_info@ml_anom
    else
        anom = False
    end if

    ;;******************************************************************
    ;;monthly averaging
    ;;******************************************************************

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
        uniform_coordinates(A0, field_type0)

        ;;convert to 1D if necessary
        if(rank .gt. 1) then

            ;; convert to T2Ms
            if (field_type0 .eq. "T2Ds") then
                tmp = day_to_month(A0, "T2Ds")
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if

            ;; cosine of latitude for averaging (= weights)
            coslat = cos(A0&lat * rad)
            coslat!0 = "lat"
            coslat&lat = A0&lat

            ;; average the data over dimensions (may need more details here)
            if (field_type0 .ne. "T1Ms") then

                if (field_type0 .eq. "T3M") then
                    tmplt = dim_avg_Wrap( \
                        A0(time | :, {plev | levrange(0):levrange(1)}, \
                           {lat | latrange(0):latrange(1)}, \
                           {lon | lonrange(0):lonrange(1)}))
                    tmpz = dim_avg_wgt_Wrap( \
                        tmplt(:, :, :), coslat({latrange(0):latrange(1)}), 1)
                    tmp = dim_avg_Wrap(tmpz)
                    delete(tmplt)
                    delete(tmpz)
                end if
                if (field_type0 .eq. "T2Ms") then
                    tmplt = dim_avg_Wrap( \
                        A0(time | :, {lat | latrange(0):latrange(1)}, \
                           {lon | lonrange(0):lonrange(1)}))
                    tmp = dim_avg_wgt_Wrap( \
                        tmplt, coslat({latrange(0):latrange(1)}), 1)
                    delete(tmplt)
                end if
                if (field_type0 .eq. "T2Mz") then
                    tmpz = dim_avg_wgt_Wrap( \
                        A0(time | :, {plev | levrange(0):levrange(1)}, \
                           {lat | latrange(0):latrange(1)}), \
                        coslat({latrange(0):latrange(1)}), 1)
                    tmp = dim_avg_Wrap(tmpz)
                    delete(tmpz)
                end if

                delete(A0)
                A0 = tmp
                delete(tmp)
            end if
        end if  ; rank = 1 loop

        ;; assign time data (X)
        ;; convert to years, months, etc (0 = years, 1 = months)
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
            tmp = ind(Xdata(imod, :) .eq. mi + 1 .and. \
                      (.not. ismissing(Ydata(imod, :))))
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
    if (isatt(diag_script_info, "ml_obs_name") \
        .and. isatt(diag_script_info, "ml_obs_file") \
        .and. isatt(diag_script_info, "ml_obs_stdv_file")) then

        f_obs = addfile(diag_script_info@ml_obs_file, "r")
        ;; Observation by HALOE are already processed:mean
        v = f_obs->$var0$
        delete(f_obs)
        f_obs = addfile(diag_script_info@ml_obs_stdv_file, "r")
        vName = var0 + "@STDEV"

        ;; Observation by HALOE  are already processed:standard deviation
        v_STD = f_obs->$vName$
        delete(f_obs)
        v!0 = "month"
        v!1 = "plev"
        v!2 = "lat"
        v_STD!0 = "month"
        v_STD!1 = "plev"
        v_STD!2 = "lat"
        if ((v@units .eq. "m^3/m^3") .and. (var0 .eq. "HCl")) then
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
        ;; cosine of latitude for averaging (= weights)
        coslat   = cos(v&lat * rad)
        coslat!0 = "lat"
        coslat&lat = v&lat

        tmpz = dim_avg_wgt_Wrap(v(month | :, {plev | levrange(0):levrange(1)},\
                                             {lat | latrange(0):latrange(1)}),\
                                coslat({latrange(0):latrange(1)}), 1)
        tmp = dim_avg_Wrap(tmpz)
        delete(tmpz)
        delete(v)
        v = tmp
        delete(tmp)

        tmpz = dim_avg_wgt_Wrap( \
            v_STD(month | :, {plev | levrange(0):levrange(1)}, \
                  {lat | latrange(0):latrange(1)}), \
            coslat({latrange(0):latrange(1)}), 1)
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
        mod_names(0) = diag_script_info@ml_obs_name
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
    rgstr = variable_info@long_name
    if(isvar("lonrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", lonrange(0))\
                      + "-" + sprintf("%3.0f", lonrange(1)) + "lon"
    end if
    if(isvar("latrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", latrange(0)) + "-"\
                             + sprintf("%3.0f", latrange(1)) + "lat"
    end if
    if(isvar("levrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", levrange(0)) + "-" \
                             + sprintf("%3.0f", levrange(1)) + "hPa"
    end if

    ;; write out array to netcdf file.
    mYdata!0 = "models"
    mYdata@info = rgstr

    mSdata!0 = "models"
    mSdata@info = rgstr

    mNdata!0 = "models"
    mNdata@info = rgstr

    mod_names!0 = "models"

    diag_script_base = basename(diag_script)
    output_filename = interface_get_figure_filename(diag_script_base, \
                                                    var0, \
                                                    field_type0,\
                                                    "", \
                                                    -1)

    output_nc_filename = output_filename + ".nc"
    output_nc_path = output_dir + output_nc_filename

    print("output_filename = " + output_nc_path)
    system("/bin/rm -f " + output_nc_path)
    ncdf = addfile(output_nc_path, "c")

    ;;  output variables directly
    ovarn = var0 + "_ann_cycle"
    ncdf->$ovarn$ = mYdata
    osvarn = var0 + "_ann_cycle_stdv"
    ncdf->$osvarn$ = mSdata
    onvarn = var0 + "_nmonths"
    ncdf->$onvarn$ = mNdata

    nn = num(.not. ismissing(mod_names))

    ;; *****************************************************************
    ;;  create global attributes of the file
    ;; *****************************************************************
    setfileoption(ncdf, "DefineMode", True)
    fAtt = True            ; assign file attributes
    fAtt@creation_date = systemfunc("date")
    do im = 0, nn - 1
        str = unique_string("model_")
        fAtt@$str$= mod_names(im)
    end do
                                ; copy attributes from the source file
    fileattdef(ncdf, fAtt)      ; copy file attributes
    setfileoption(ncdf, "DefineMode", False)

    ;; ================================================
    ;; METRICS
    ;; Make sure time period is the same!
    ;; ================================================
    if (isatt(diag_script_info, "ml_stdv_plot")) then
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
            if (mod_names(imod) .eq. diag_script_info@ml_stdv_plot) then
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
            if(var0 .eq. "ptp") then
                dd = 10.
            else
                dd = ng
            end if

            gm(imod) = 1. - abs(avg(mYdata(imod, :)) - ref_avg) / \
                stddev(ref_anom) / dd

            if (gm(imod) .lt. 0.) then
                gm(imod) = 0.
            end if

            ;; B: "correlation"
            ;; need to be able to match years!
            gg = ind(.not. ismissing(ref_anom))
            gc(imod) = (esccr(ts(gg), ref_anom(gg), 0) + 1.) / 2.

            ;; C: "variance"
            ;;  gv = 1 - 1/ng * (sigma(mod)-sigma(obs))/sigma(obs)
            gv(imod) = 1. - abs((stddev(ts) - stddev(ref_anom))) / \
                stddev(ref_anom) / ng
            if (gv(imod) .lt. 0.) then
                gv(imod) = 0.
            end if

            ;; SUM
            gsum(imod) = (gm(imod) + gc(imod) + gv(imod)) / 3.

            ;; Text
            gradetx(imod) = sprintf("%5.2f", gm(imod)) + " " \
                + sprintf("%5.2f", gc(imod)) + " " \
                + sprintf("%5.2f", gv(imod)) + " " \
                + sprintf("%5.2f", gsum(imod)) + " " \
                + mod_names(imod)
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
    LineThick = CCMVal_Thick(mod_names)

    Ccolors = new((/nn - 1/), string)
    Cdashs = new((/nn - 1/), integer)
    CLineThick = new((/nn - 1/), integer)
    CLineThick(:) = 2

    if(anom .eq. True) then
        filename_aux = "ANN_anom"
    else
        filename_aux = "ANN_c"
    end if

    output_filename = interface_get_figure_filename(diag_script_base, \
                                                    var0, \
                                                    field_type0,\
                                                    filename_aux, -1)
    output_filename = output_filename + "." + file_type

    wks = gsn_open_wks(file_type, output_dir + output_filename)

    ;; Common resources for xy plots
    plot = new(2, "graphic")
    res                        = True
    res@gsnDraw                = False  ; Do not draw plot
    res@gsnFrame               = False  ; Do not advance frame
    res@trXReverse             = False
    res@tiXAxisString          = "Month"
    res@tiXAxisFontHeightF     = 0.020
    res@tiYAxisFontHeightF     = 0.020
    res@tiMainFontHeightF      = 0.025

    res@pmLegendSide           = "Right"  ; Change location of

    res@pmLegendOrthogonalPosF = 0.02  ; of the legend
    res@pmLegendWidthF         = 0.15  ; Change width and
    res@lgPerimOn              = True
    if (dim_MOD .le. 4) then
        res@pmLegendHeightF      = 0.10  ; height of legend.
        res@lgLabelFontHeightF   = .020
        res@pmLegendParallelPosF = 0.10  ; of the legend
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
    res@tiYAxisString          = variable_info@long_name \
        + " (" + variable_info@units + ")"
    res@xyExplicitLegendLabels = lgLabels
    res@xyMonoLineColor        = False

    ;; ************************************************************
    ;; XAxis
    Xarr = moarr

    res@vpHeightF = 0.4                    ; change aspect ratio of plot
    res@vpWidthF = 0.8

    res@trXMaxF = 13
    res@trXMinF = 0

    res@tmXBMode = "Explicit"
    res@tmXBValues = Xarr
    res@tmXBLabels = monames

    ;; ************************************************************
    ;; YAxis
    res@trYMinF = min(mYdata)
    res@trYMaxF = max(mYdata)

    ;; ************************************************************
    ;; PLOT
    res@tfPolyDrawOrder = "Predraw"               ; put line on top

    plot(0) = gsn_csm_xy(wks, Xarr, mYdata, res)   ; cases line plot

    ;; ************************************************************
    ;; plot standard deviation
    if (isatt(diag_script_info, "ml_stdv_plot"))

        ii = -99
        do imod = 0, dim_MOD - 1
            if (mod_names(imod) .eq. diag_script_info@ml_stdv_plot) then
                ii = imod
            end if
        end do

        if(ii .lt. 0) then
            error_msg("f", diag_script, "", "ml_stdv_plot model (" + \
                      diag_script_info@ml_stdv_plot + ") not found")
        end if

        ;;DEFINE A POLYGON BASED on selected 'model' (0 = obs)
        xp = new((/2 * nmos/), float)
        yp = new((/2 * nmos/), float)
        do k = 0, nmos - 1
            yp(k) = mYdata(ii, k) + 2. * mSdata(ii, k)
            xp(k) = mYdata&month(k)
            xp(2 * nmos - 1 - k) = mYdata&month(k)
            yp(2 * nmos - 1 - k) = mYdata(ii, k) - 2. * mSdata(ii, k)
        end do

        gsres = True  ; poly res
        gsres@gsFillColor = "Gray"  ; color chosen
        dummy = gsn_add_polygon(wks, plot(0), xp, yp, gsres)
    end if

    ;; Common resources for panel plots
    pan = True
    pan@gsnMaximize = True
    pan@gsnPaperOrientation = "portrait"
    pan@gsnFrame = False
    gsn_panel(wks, plot, (/2, 1/), pan)  ; panel plot
    frame(wks)

    delete(mYdata)
    delete(mSdata)
    delete(mNdata)
    delete(Ydata)
    delete(Xdata)
    delete(res)
    delete(plot)
    delete(pan)
    info_output("<<<<<<<< Leaving " + diag_script, verbosity, 4)
end
