;;#############################################################################
;; GO_ts_line.ncl
;; Author: Oliver Andrews (UEA, Great Britain)
;;         Clare Enright (UEA, Great Britain)
;; EMBRACE project
;;#############################################################################
;; Description
;;      Evaluation of ocean biogeochemical cycles at global scales, as
;;      simulated by both ocean-only and coupled climate-carbon cycle models. 
;;      This diagnostic (developed from SeaIce_tsline.ncl) produces annual
;;      and seasonal time series plots from monthly ocean data, and includes
;;      reformatting options to convert input irregularly gridded data into a
;;      version that NCL can interpret by introducing intermediate
;;      referencing co-ordinates in dimensioning the variable, and writing
;;      cell areas to the same file
;;
;; Required diag_script_info attributes (diagnostics specific)
;;      ts_lonrange    array of longitudes for time series plots if 2D or 3D
;;      ts_latrange    array of latitudes for time series plots if 2D or 3D
;;      ts_levrange    array of pressure levels for time series plots if 3D
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;      ts_ymax        y axis maximum for time series plots
;;      ts_ymin        y axis minimum for time series plots
;;      ts_xmax_year   start year for time series plots
;;      ts_xmin_year   end year for time series plots
;;      ts_anom        plot anomalies.
;;      ts_smooth      smooth...
;;      ts_trend       overlay trend line
;;
;; Required variable_info attributes (variable specific)
;;    long_name:                    Name to store in derived varible
;;    units:                        Units for derived variable
;;
;; Caveats
;;
;; Modification history
;;    20151217-A_laue_ax: changed to new version of "add_legenda" (defined
;;                        in diag_scripts/lib/ncl/misc_function.ncl) introduced
;;                        with namelist_eyring13jgr
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150626-A_andr_ol: written
;; ############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  ; Needed for Wrap
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
; A-laue_ax+
;load "./diag_scripts/to_be_checked/CCMVal_FUNCTIONS/misc_function.ncl"
load "./diag_scripts/lib/ncl/misc_function.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
; A-laue_ax-
load "diag_scripts/aux/GO/GO_functions/check_params.ncl"
load "diag_scripts/aux/GO/GO_functions/set_ref_mask.ncl"
load "diag_scripts/aux/GO/GO_functions/output_area.ncl"

load "diag_scripts/aux/GO/GO_functions/calc_multimean.ncl"
load "diag_scripts/aux/GO/GO_functions/calc_anom.ncl"
load "diag_scripts/aux/GO/GO_functions/calc_metrics.ncl"
;
; Read data for all models and set root for plot file names
;
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "GO_tsline"
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    ;; The 'models' array has been loaded interface_data/ncl.interface
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    variable = variables(0)
    diag_file = variable + "_tsline"
    field_type0 = field_types(0)

    ;;-----------------------
    ;; References Information
    ;;-----------------------

    write_references(diag_script,                                 \  ; script name
                     (/"A_andr_ol", "A_enri_cl" /),               \  ; authors
                     (/"A_eyri_ve", "A_gett_an", "A_cion_ir"/),   \  ; contributors
                     "",                                          \  ; diag_references
                     "",                                          \  ; obs_references
                     "P_embrace")                                    ; proj_acknowledgements

    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, variable, field_type0)
    end do
end
; Load variable_info
begin
    vardeffile = "interface_data/" + variable + "_info.tmp"
    loadscript(vardeffile)
end
;

begin
    ;;****************************
    ;; initialisation
    ;;****************************
    plot_dir = getenv("ESMValTool_plot_dir")
    plot_script_base = basename(diag_script)
    diri_out = get_output_dir(plot_dir, plot_script_base)
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if
    color_type = "COLOR"
    file_type = str_lower(file_type)
    ;; The 'models' array has been loaded interface_data/ncl.interface
    dim_MOD = dimsizes(models@name)
    info_output(" Got n models: " + dim_MOD, verbosity, 3)
    A = extract_data(0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(A))
    delete(A)
    refModel = "NULL"
    latrange = new(2, float)
    lonrange = new(2, float)
    units = "NULL"
    if(isatt(variable_info, "units")) then
        units = variable_info@units
    end if

    ;; Ignore some warnings from ESMF-regridding
    ignore_warning_file = "interface_data/warnings_to_ignore.txt"
    warnings_to_ignore = \
        (/"warning:ContourPlotInitialize: scalar field is constant"/)
    write_ignore_warnings(warnings_to_ignore, ignore_warning_file)

    anom = False
    smooth = False
    trend = False
    multimean = False
    coverage = False
    yMax_lev = -1.e35
    yMin_lev = 1.e35
    xMax_year = 2999.
    xMin_year = 00.
    check_params(diag_script,\
                 field_types,\
                 diag_script_info,\
                 rank,\
                 refModel,\
                 latrange,\
                 lonrange,\
                 yMax_lev,\
                 yMin_lev,\
                 xMax_year,\
                 xMin_year,\
                 anom, smooth,\
                 trend,\
                 multimean,\
                 coverage)
    info_output(" x limits: " + xMax_year + xMin_year, verbosity, 3)
    ;; ******************************************************************
    ;; seasonal averaging
    ;; ******************************************************************
    seasdt = (/ (/-1, 2./), (/2, 5/), (/5, 8/), (/8, 11./) /) / 12.
    seasons = (/"DJF", "MAM", "JJA", "SON"/)
    ns = dimsizes(seasons)

;same as above, but want y+s(0, s):y+s(1, s)
;where  DJF        MAM     JJA    SON
;0      -11./12.   2/12.   5/12.  8/12
;1      2./12.     5/12.   8/12.  11/12
;only use if there are 3 values in each.

    ;; find time range (if not set?)
    ;; radians for cosine of latitude weights
    rad = 4.0 * atan(1.0) / 180.0

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
        info_output("Model " + models@name(imod), verbosity, 3)
        info_output("A0 dimensions here: " + getvardims(A0), verbosity, 3)
        if (imod .eq. 0 .and. (isatt(diag_script_info, "refModel"))) then
            ref_mask = set_ref_mask(diag_script_info, imod, models, refModel, A0)
        end if
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
        info_output("dimsizes: " + dimsizes(A0), verbosity, 3)

        ;; convert to 1D if necessary
        if(rank .gt. 1) then

            ;; convert to T2Ms
            if (field_type0 .eq. "T2Ds") then
                tmp = day_to_month(A0, "T2Ds")
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if

        ;; further processing depends on whether variable is on a lat-lon grid
        if(all(iscoord(A0, (/"lat", "lon"/)))) then
            info_output("GO_tsline.ncl:  lat - lon grid", verbosity, 3)
            area_temp = map_area(A0&lat, A0&lon)

            ;; make 2D to be compatible with irregular grid handling below
            ;; get 1D arrays (area_lat and area_lon) of lats and lons in A0
            area_lat = conform_dims(dimsizes(area_temp), A0&lat, 0)
            area_lon = conform_dims(dimsizes(area_temp), A0&lon, 1)
            ;; set area fill value to be the same as that for data
            if (isatt(A0, "_FillValue")) then
                if (isatt(area_temp, "_FillValue")) then
                    area_temp = where(area_temp .eq. area_temp@_FillValue,\
                                      A0@_FillValue, area_temp)
                end if
                ;; Set area to _FillValue where data is _FillValue
                area_temp@_FillValue = A0@_FillValue

                ;; re-order dimensions of A0
                tmp = A0(time |:, lat |:, lon |:)
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if
        else  ; irregular grid

            info_output("GO_tsline.ncl:  not on a lat-lon grid", verbosity, 3)

            ;; get area data & lat-lon from another file
            ;; (model specific: reformat folder)
            area_lon = extract_area(imod, data, "lon", variable)
            area_lat = extract_area(imod, data, "lat", variable)
            area_temp = extract_area(imod, data, "areacello", variable)
            info_output("area sizes here: " + dimsizes(area_temp), verbosity, 3)
        end if  ; irregular grid
        if (refModel .ne. "NULL") then
            if (imod .eq. 0) then
                if(all(iscoord(A0, (/"lat", "lon"/)))) then
                    tmp = ref_mask(time |:, lat |:, lon |:)
                else
                    ;; Not working yet; set-ref_mask not for irregular grids
                    ;; set_ref_mask shouldn't let you get here
                    tmp = ref_mask(time |:, i |:, j |:)
                end if
                delete(ref_mask)
                ref_mask = tmp
                delete(tmp)
            else
                A0 = A0 * ref_mask
            end if
        end if

        if (field_type0 .eq. "T3M") then
            print("fatal: Not implemented for field type T3M")
            status_exit(1)
        end if
        if (field_type0 .eq. "T2Ms" .or. field_type0 .eq. "TO2Ms") then
            ;; get all longitudes 0 to 360 degrees
            area_lon = where(area_lon .lt. 0., 360. + area_lon, area_lon)
            area_lon = where(area_lon .gt. 360., area_lon - 360., area_lon)

            ;; Select region; set all area outside region to zero
            area_temp = where(area_lat .ge. latrange(0) .and. area_lat .le. latrange(1), area_temp, 0.0)
            if (lonrange(0) .lt.  lonrange(1)) then
                area_temp = where(area_lon .ge. lonrange(0) .and. area_lon .le. lonrange(1), area_temp, 0.0)
            else
                area_temp = where((area_lon .ge. lonrange(0) .and. area_lon .lt. 360.) .or. area_lon .le.  lonrange(1), area_temp, 0.0)
            end if
            info_output("diri_out " + diri_out, verbosity, 3)
            info_output("imod " + imod, verbosity, 3)

            ;; output area extracted to a file for testing purposes only
            output_area(diri_out, imod, area_lon, area_lat, area_temp)
            info_output("A0 sizes here: " + dimsizes(A0), verbosity, 3)

            fact = conform(A0, area_temp, (/1, 2/))
            A1 = A0 * fact

            total_area = new(dimsizes(A0), typeof(A0), A0@_FillValue)
            total_area = where(A1 .ne. A1@_FillValue, fact, total_area)
            ;; Calculate area of cells filled for each year
            tot_area_yr = dim_sum(dim_sum(total_area))
            tot_area_yr = where(tot_area_yr .ne. 0., tot_area_yr, A1@_FillValue)
            copy_VarMeta(A0, A1)
            dimNames = getvardims(A1)
            tmpln = dim_sum_Wrap(A1)
            tmp = dim_sum_Wrap(tmpln)

            tmp = tmp / tot_area_yr
            ;; Store area covered if model is the reference (observations)
            if (models@name(imod) .eq. refModel) then
                obs_area = tot_area_yr
            end if

            delete(tmpln)
            delete(total_area)
            delete(fact)
            delete(tot_area_yr)
        end if
        if (field_type0 .eq. "T2Mz") then
            print("fatal:Not implemented for field type T2Mz")
            status_exit(1)
        end if
        if (field_type0 .eq. "T1Ms") then
            print("fatal: Not implemented for field type T1Ms")
            status_exit(1)
        end if

        delete(A0)
        delete(A1)
        A0 = tmp
        delete(tmp)
        end if  ; rank = 1 loop

        ;; assign time data (X)
        ;; convert to years, months, etc (0 = years, 1 = months)
        ;; needs units origin attribute e.g. "days since ...."
        ;; note_cd_calendar requires NCL v6.0
        ti = A0&time
        ti2 = cd_calendar(ti, 0)

        ;; first: make array in decimal years...
        tt = ti2(:, 0) + (ti2(:, 1) - 1) / 12.

        Xdata(imod, 0:dimX(imod) - 1) = (/tt/)
        ;; assign data (Y)
        Ydata(imod, 0:dimX(imod) - 1) = (/A0/)
        info_output("Model: " + imod, verbosity, 3)

        info_output("average of data: " + avg(Ydata(imod, :)), verbosity, 3)

        delete(A0)

        delete(area_temp)
        delete(area_lat)
        delete(area_lon)
        delete(ti)
        delete(ti2)
        delete(tt)
    end do  ; imod

    ;; MAKE MULTIMODEL MEAN HERE.... Xdata, Ydata and dim_MOD = dim_MOD + 1
    ;; need to average appropropriate times!  Quick and dirty: take from model 0
    if (multimean .eq. True) then
        ;; Return new arays
        new_arr = calc_multimean(models, Xdata, Ydata,  dim_MOD, nMax, dimX)
        delete(Xdata)
        delete(Ydata)
        Xdata = new_arr[0]
        Ydata = new_arr[1]
        delete(new_arr)
    end if

    ;; average seasonal and annual
    ;; sYdata, aYdata
    aYdata = new((/dim_MOD, nMax/), "float")
    aXdata = new((/dim_MOD, nMax/), "float")
    ;; set number of years available for each model
    nyr = new((/ nMax /), "integer")

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
        ;; ispan returns integers from min to max, stride ( = 1)
        tmp = ispan(floattoint(min(Xdata(imod, :))), \
                    floattoint(max(Xdata(imod, :))), 1) * 1.

        nyr(imod) = dimsizes(tmp)
        info_output("nyr: " + nyr(imod), verbosity, 3)

        aXdata(imod, 0:nyr(imod) - 1) = tmp
        delete(tmp)
        ;; annual average - copy times and data to tmp
        tmp = Ydata(imod, :)
        tmp!0 = "time"
        tmp&time = Xdata(imod, :)
        info_output("Model: " + imod, verbosity, 3)
        info_output("average of tmp: " + avg(tmp), verbosity, 3)
        info_output("nyr: " + nyr(imod), verbosity, 3)
        ;; For reference model calculate average area
        ;; where observations exist
        if (models@name(imod) .eq. refModel) then
            annl_avg_obs_area = new((/nyr(imod)/), "float")
        end if
        do y = 0, nyr(imod) - 1
            ;; should use 'month_to_annual_weighted' function...
            n = 0

            ttt = ind(floattoint(tmp&time) .eq. floattoint(aXdata(imod, y)))
            aa = tmp(ttt)

            n = dimsizes(aa)
            ;; no. of months it is sensible to average  depends on the variable
            ;; and sort of observations - set this but default to 12 for now

            if (n(0) .ge. 12.) then
                if(.not. all(ismissing(aa))) then
                    aYdata(imod, y) = avg(aa)
                    if (models@name(imod) .eq. refModel) then
                        oa = obs_area(ttt)
                        annl_avg_obs_area(y) = avg(oa)
                    end if
                end if
            end if
            delete(aa)
            delete(ttt)
            delete(n)

            ;; loop for seasons
            ;; same as above, but want y+s(0, s):y+s(1, s)
            ;; where  DJF        MAM     JJA    SON
            ;; 0      -11./12.   2/12.   5/12.  8/12
            ;; 1      2./12.     5/12.   8/12.  11/12
            ;; only use if there are 3 values in each.
            do s = 0, ns - 1
                if (aXdata(imod, y) + seasdt(s, 0) .lt. max(aXdata(imod, :) + 1)) then
                    n = 0

                    ti = tmp&time
                    tr = aXdata(imod, y)

                    ttt = ind(ti .ge. tr + seasdt(s, 0) .and. ti .lt. tr + seasdt(s, 1))

                    ssz = dimsizes(ttt)
                    if (ssz(0) .gt. 1) then
                        aa = tmp(ttt)

                        n = dimsizes(aa)

                        N = num(.not. ismissing(aa))

                        if (n(0) .ge. 2) then
                            if (N .gt. 2) then
                                sYdata(imod, y, s) = avg(aa)
                                if (models@name(imod) .eq. refModel) then
                                    seasnl_avg_obs_area = avg(obs_area)
                                end if
                            end if
                        else
                            info_output("Missing months: " + imod + " year:" + y + "season: " + s, verbosity, 4)
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
            calc_anom(imod, aYdata, ns, sYdata)
        end if

        ;; calculate trends here using regline (seasonal and annual).
        ;; save lines for plotting (add line). may also want data.
        ;; annual
        rc = regline(aXdata(imod, :), aYdata(imod, :))
        df = rc@nptxy - 2
        aYreg(imod, :) = rc * aXdata(imod, :) + rc@yintercept
        arc(imod) = rc
        asig(imod) = (1. - betainc(df / (df + rc@tval ^ 2), df / 2.0, 0.5)) * 100.

        delete(rc)
        delete(df)

        ;; seasonal
        do s = 0, ns - 1
            rc = regline(aXdata(imod, :), sYdata(imod, :, s))
            ;; sXdata dowsn't exist so if you wan to different lengths for
            ;; time will need to introduce it.
            df = rc@nptxy - 2
            sYreg(imod, :, s) = rc * aXdata(imod, :) + rc@yintercept
            src(imod, s) = rc
            ssig(imod, s) = (1. - betainc(df / (df + rc@tval ^ 2), df / 2.0, 0.5)) * 100.
            delete(rc)
            delete(df)
            info_output(models@name(imod) + "  SEAS  Trend ("\
                                          + variable_info@units + "/season) = " \
                                          + src(imod, s)\
                                          + "  Significance (%) = " \
                                          + ssig(imod, s), verbosity, 3)
        end do  ; seasonal lines

        ;; write trends and significance
        info_output(models@name(imod) + "  ANN  Trend (" + variable_info@units + "/yr) = " \
                                      + arc(imod) + "  Significance (%) = " \
                                      + asig(imod), verbosity, 3)

    end do     ; model loop

    ;; METRICS
    ;; Make sure time period is the same!
    iref = -99
    if (isatt(diag_script_info, "refModel")) then
        calc_metrics(models, dim_MOD, refModel, iref,  aYdata)
    end if
    ;; end metrics

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

    ;; ******************************************************************
    ;; Use project_style to associate model with color  and dash pattern
    ;; ******************************************************************
    ;; ******************************************************************
    ;; select colors and other plotting attributes
    ;; colours are set in diag_script_info@colors. Dashes used if not available
    ;; styleset EMAC for setting colours by experimnet rathewr than model name
    colors = project_style(diag_script_info, "colors")
    dashs = project_style(diag_script_info, "dashes")

; A_laue_ax+
    LineScatters = new((/dim_MOD/), string)
; A_laue_ax-

    if (multimean .eq. False) then
        lgLabels = models@name + " " + models@experiment
    else
        lgLabels = new((/dim_MOD/), "string")
        info_output("size of names: " + dimsizes(models@name), verbosity, 4)
        info_output("size of experiments: " + dimsizes(models@experiment), verbosity, 4)

        lgLabels(0) = models@name(0)
        lgLabels(1:) = models@name(1:) + " " + models@experiment(0:)
    end if
    TLineThick = new((/dim_MOD/), integer)
    TLineThick(:) = 1

    LineThick = project_style(diag_script_info, "thicks")

    ;; make string for region
    rgstr = variable_info@long_name
    filename_range = ""
    if(isvar("lonrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", lonrange(0)) + "-" \
                             + sprintf("%3.0f", lonrange(1)) + "lon"
        filename_range = filename_range + sprintf("%.1f", lonrange(0)) + "-" \
                                        + sprintf("%.1f", lonrange(1)) + "lon"
    end if
    if(isvar("latrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", latrange(0)) + "-" \
                             + sprintf("%3.0f", latrange(1)) + "lat"
    filename_range = filename_range + sprintf("%.1f", latrange(0)) + "-" \
                                    + sprintf("%.1f", latrange(1)) + "lat"
    end if
    if(isvar("levrange")) then
        rgstr = rgstr + ", " + sprintf("%3.0f", levrange(0)) + "-" \
                             + sprintf("%3.0f", levrange(1)) + "m"
    filename_range = filename_range + sprintf("%.1f", levrange(0)) + "-" \
                                    + sprintf("%.1f", levrange(1)) + "m"
    end if

    if(anom .eq. True) then
        filename_aux = filename_range + "_ANN_anom"
    else
        filename_aux = filename_range + "_ANN_c"
    end if
    info_output(plot_script_base + variable + field_type0 + filename_aux, verbosity, 4)
    output_filename = interface_get_figure_filename(plot_script_base, \
                                                    variable, \
                                                    field_type0,\
                                                    filename_aux,\
                                                    -1)
    output_filename = output_filename + "." + file_type

    wks = gsn_open_wks(file_type, diri_out + "/" + output_filename)
    plot = new(2, "graphic")
    res = True
    res@gsnDraw = False               ; Do not draw plot
    res@gsnFrame = False              ; Do not advance frame

    res@tiXAxisString = "Year"
    res@tiXAxisFontHeightF = 0.020
    res@tiYAxisFontHeightF = 0.020
    res@tiMainFontHeightF = 0.025

    res@vpHeightF = 0.4               ; change aspect ratio of plot
    res@vpWidthF = 0.8

    res@pmLegendSide = "Right"        ; Change location of

    res@pmLegendOrthogonalPosF = 0.01  ; of the legend
    res@pmLegendWidthF = 0.2           ; Change width and
    res@lgPerimOn = True
    if (dim_MOD .le. 4) then
        res@pmLegendHeightF = 0.10     ; height of legend.
        res@lgLabelFontHeightF = .020
        res@pmLegendParallelPosF = 0.10
    else
        if (dim_MOD .le. 8) then
            res@pmLegendHeightF = 0.18
            res@lgLabelFontHeightF = 0.016
            res@pmLegendParallelPosF = 0.18
        else
            res@pmLegendHeightF = 0.5
            res@lgLabelFontHeightF = 0.016
            res@pmLegendParallelPosF = 0.50
        end if
    end if
    res@txFontHeightF = 0.015
    res@xyDashPatterns = dashs
    res@xyLineThicknesses = LineThick
    if (color_type .eq. "COLOR") then
        res@xyLineColors = colors
    end if

    res@tiMainString = "ANN " + rgstr
    res@pmLegendDisplayMode = "Always"
    res@tiYAxisString = variable_info@long_name + " (" + variable_info@units + ")"
    res@xyExplicitLegendLabels = lgLabels
    res@xyMonoLineColor = False

    if(isatt(variable_info, "units") .and. variable_info@units .eq. "hPa") then
        res@trYReverse = True
    else
        res@trYReverse = False
    end if
    ;; XAxis
    if (isvar("xMin_year") .eq. True) then
        res@trXMinF = xMin_year
    else
        res@trXMinF = min(aXdata)
    end if
    if (isvar("xMax_year") .eq. True) then
        res@trXMaxF = xMax_year
    else
        res@trXMaxF = max(aXdata)
    end if

    ;; ************************************************************
    ;; YAxis
    res@trYMinF = min(yMin_lev)

    info_output("ymax: " + yMax_lev, verbosity, 4)
    res@trYMaxF = max(yMax_lev)

    res2 = True
    res2@xyLineThicknesses = 2
    res2@tiYAxisString = "Area covered for reference model km^2"
    res2@tiYAxisFontHeightF = 0.020

    info_output(" coverage: " + coverage, verbosity, 4)
    info_output(" iref " + iref, verbosity, 4)
    if (coverage .eq. True .and.  iref .gt. -99) then
        ;; Plot the reference model with the area coverage;
        res@pmLegendDisplayMode = "Never"
        plot(0) = gsn_csm_xy2(wks, aXdata(iref, 0:nyr(0) - 1), aYdata(iref, 0:nyr(0) - 1), annl_avg_obs_area, res, res2)

        info_output("colors: " + colors, verbosity, 4)
        info_output("dashs: " + dashs, verbosity, 4)

        do imod = 0, dim_MOD - 1
            if (imod .ne. iref) then
                add_line(wks, plot(0), aXdata(imod, 0:nyr(0) - 1), aYdata(imod, 0:nyr(0) - 1), colors(imod), dashs(imod))
            end if
        end do
; A-laue_ax+
        ; change to extended version of "add_legenda" introduced with namelist_eyring13jgr
;        add_legenda(wks, lgLabels, colors, stringtoint(dashs), "bottom")
        add_legenda(wks, lgLabels, colors, stringtoint(dashs), LineScatters, LineThick, "bottom")
; A-laue_ax-
        res@pmLegendDisplayMode = "Always"
    else
        info_output(" max x; " + max(aXdata), verbosity, 4)
        plot(0) = gsn_csm_xy(wks, aXdata, aYdata, res)

        info_output("colors: " + colors, verbosity, 4)
        info_output("dashs: " + dashs, verbosity, 4)
    end if

    ;; cases line plot
    ;; add trend lines
    LineThick(:) = 2
    if (trend .eq. True) then
        res@xyLineThicknesses = LineThick
        do imod = 0, dim_MOD - 1
            add_line(wks,\
                     plot(0),\
                     aXdata(imod, :),\
                     aYreg(imod, :),\
                     colors(imod),\
                     dashs(imod))
        end do
    end if

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
        filename_aux = filename_range + "_SEAS_anom"
    else
        filename_aux = filename_range + "_SEAS_c"
    end if

    output_filename = interface_get_figure_filename(plot_script_base, \
                                                    variable, \
                                                    field_type0,\
                                                    filename_aux, -1)
    output_filename = output_filename + "." + file_type

    wks = gsn_open_wks(file_type, diri_out + "/" + output_filename)

    do s = 0, ns - 1
        ;; ***************************************************************
        ;; --------------------- Open files for plots --------------------
        res@tiMainString = seasons(s) + rgstr
        LineThick(:) = 1
        res@xyLineThicknesses = LineThick

        ;; ************************************************************
        ;; YAxis
        res@trYMinF = min(sYdata(:, :, s))
        res@trYMaxF = max(sYdata(:, :, s))
        info_output("Sy data min: " + res@trYMinF, verbosity, 4)
        if (yMin_lev .lt. res@trYMinF) then
            res@trYMinF = min(yMin_lev)
        end if
        res@trYMaxF = max(sYdata)
        info_output("y data max: " + res@trYMaxF, verbosity, 4)
        if (yMax_lev .gt. res@trYMinF) then
            res@trYMaxF = max(yMax_lev)
        end if
        plot(s) = gsn_csm_xy(wks, aXdata, sYdata(:, :, s), res)

        ;; add trend lines
        if (trend .eq. True) then
            LineThick(:) = 2
            res@xyLineThicknesses = LineThick
            do imod = 0, dim_MOD - 1
                add_line(wks,\
                         plot(s),\
                         aXdata(imod, :),\
                         sYreg(imod, :, s),\
                         colors(imod),\
                         dashs(imod))
            end do
        end if
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
    info_output("<<<<<<<< Leaving " + diag_script + ".ncl", verbosity, 4)
end
