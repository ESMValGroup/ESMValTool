;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# it is based on E06FIG12B.ncl by
;# Irene Cionni (i.cionni@niwa.co.nz) and Veronika Eyring (veronika.eyring@dlr.de)#
;##################################################################################
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package (PI: Andrew Gettelman, co-I: Hisako Shiona)
;
; The CCMVal script has been furthered developed as a part of the
; ESMValTool, an ESM Evaluation Tool within the EMBRACE project. Package
; leaders are Klaus Wyser (SMHI) and Veronika Eyring (DLR).
;
; Description:
; The goal is to take a daily variable of arbitrary dimensions, average it to
; 1D in space and 1D in time, compute its running average and plot it as a line.
;
; info attributes required:
; either read data and sort by variables, or read type...
;
; contact: Martin.Evaldsson@smhi.se
; ts_lonrange    array of pressure levels for time series plots if 2D or 3D
; ts_latrange    array of latitudes for time series plots if 2D or 3D
; ts_levrange    array of pressure levels for time series plots if 3D
; ts_ymax        (optional) y axis maximum for time series plots
; ts_ymin        (optional) y axis minimum for time series plots
; ts_xmax_year   (optional) start year for time series plots
; ts_xmin_year   (optional) end year for time series plots
; ts_anom        (optional) plot anomalies.
; ts_smooth      (optional) smooth...

;Note: cd_calendar requires NCL 6.0
;  A) Remove annual cycle option (anomalies) DONE: needs to be 'right' however
;     small adjustment for variable season/month length (minor)
;
; NOTICE - this routine is _not_ complete yet


;*******************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/time_axis_labels.ncl"

begin

    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot line_runave (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    ;; fields supported
    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ds"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    diri_out = plot_dir + "/line_runave/"
    file_type = "PS"

    dim_MOD = dimsizes(models@name)

    A = extract_data(0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(A))
    delete(A)

    ;;******************************************************************
    ;;check the input parameters
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

    end if  ; rank = 1

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

    ;; radians for cosine of latitude weights
    rad    = 4.0 * atan(1.0) / 180.0

    ;;==================================================================
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

    ;; get data into arrays (X) & (Y)
    ;;loop for models
    imod = 0

    ;;get data
    A0 = extract_data(imod, data, -1, 0, 0)
    uniform_coordinates(A0, field_number)

    ;; convert to 1D if necessary
    if(rank .gt. 1) then

        ;;cosine of latitude for averaging (=weights)
        coslat   = cos(A0&lat * rad)
        coslat!0 = "lat"
        coslat&lat = A0&lat

        ;; average the data over dimensions (may need more details here)
        tmp_runave = runave_n_Wrap(A0, 30, 0, 0)
        time = A0&time
        gregdate = cd_calendar(time, -2)
        yyyyddd = yyyymmdd_to_yyyyddd(gregdate)
        climdayave = clmDayTLL(tmp_runave, yyyyddd)
        tmplt = dim_avg_Wrap(climdayave(year_day | :, {lat | latrange(0):latrange(1)}, \
                               {lon | lonrange(0):lonrange(1)}))
        tmp = dim_avg_Wrap(tmplt(year_day |:, {lat|latrange(0):latrange(1)}))
        delete(tmplt)
        delete(A0)

        A0 = tmp
        delete(climdayave)
        delete(gregdate)
        delete(tmp)

    end if  ; rank =1 loop

    ;;***************************************************************
    ;;--------------------- Open files for plots --------------------
    wks = gsn_open_wks(file_type, diri_out + "line_runave_" + variable + "_ANN_c")

    ;;***************************************************************
    ;; common resources for xy plots
    plot = new(2, "graphic")
    res                        = True
    restick                    = True
    restick                = True
    restick@ttmFormat      = "%c"
    restick@ttmAxis        = "XB"
    restick@ttmMajorStride = 33
    delete(time)
    time = ispan(1, 366, 1)
    time!0 = "time"
    time&time = time
    time@bounds = "time_bnds"
    time@units = "days since 1991-1-1"
    time@calendar = "gregorian"
    time@axis = "T"
    time@long_name = "time"
    time@standard_name = "time"

    time_axis_labels(time, res, restick)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", \
                       "A", "S", "O", "N", "D"/)
    plot(0) = gsn_csm_xy(wks, time, A0, res)   ; cases line plot

    delete(plot)
    delete(Ydata)
    delete(Xdata)
    delete(res)
end
;##############################################################################
