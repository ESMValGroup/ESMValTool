;;#############################################################################
;; PLOTS TIME SERIES
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; ESMVal project
;;#############################################################################
;;
;; Description
;;    Takes an input field of arbitrary dimensions and averages it over the
;;    horizontal and vertical coordinates creating a 1D time-series.
;;    It optionally calculates yearly means.
;;    This script is based on the original CCMVal script tsline.ncl and has
;;    been further developed as part of the ESMValTool.
;;
;; Required variable_info attributes (variable specific)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    ts_minlev: minimum level for area-averaging (hPa), -1 for min level.
;;    ts_maxlev: maximum level for area-averaging (hPa), -1 for max level.
;;    ts_minlat: minimum latitude for area-averaging.
;;    ts_maxlat: maximum latitude for area-averaging.
;;    ts_minlon: minimum longitude for area-averaging (mind the order!!).
;;    ts_maxlon: maximum longitude for area-averaging (mind the order!!).
;;    multi_model_mean: multi-model mean and stddev calculation.
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    time_avg: type of time average (currently only "yearly" is available).
;;    ts_lsmask: extract only selected mask ("land", "ocean" or "obs")
;;    ref_obs: reference model for subsampling if lsmask = "obs"
;;    ts_anomaly: calculates anomalies with respect to the first 10-year average
;;    ref_start: start year of reference period for anomalies
;;    ref_end: end year of reference period for anomalies
;;    area_opper: type of area opperation (currently only "average" (default) 
;;                and "sum" are available).
;;
;; Caveats
;;    Multi-model mean to be implemented.
;;
;; Modification history
;;    20160905-A-Bock_li: implemented multi-model mean and variable refence 
;;                        period for anomaly with refence mean value.
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150622-A_wenz_sa: added optional anomaly calculation, choise of 
;;                        area opperations (sum, average) and detrending of time
;;                        series.
;;    20150420-A_righ_ma: adapted to new structure. Some of the original
;;                        features not ported yet (seasonal average, smoothing,
;;                        multi-model mean, etc.).
;;    20??????-A_cion_ir: written.
;;
;;#############################################################################
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; References
    write_references(diag_script, "A_cion_ir", (/"A_righ_ma","A_wenz_sa"/), "D_0001", \
                     "", "P_esmval")

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

    ;; Check field type
    plot_fields = (/"T3M", "T2Ms", "T2Mz", "T1Ms", "T0M"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field " + field_type0)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; Set output directory and file type
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)
    if (isdefined("output_file_type")) then
        file_type = str_lower(output_file_type)
    else
        file_type = "ps"
    end if
    color_type = "COLOR"

    A = extract_data(0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(A))
    delete(A)

    ;; Check requried diag_script_info attributes
    req_atts1 = (/"multi_model_mean"/)
    if (field_type0.eq."T3M") then
        req_atts2 = (/"ts_minlev", "ts_maxlev", "ts_minlat", "ts_maxlat", \
                      "ts_minlon", "ts_maxlon"/)
    end if
    if (field_type0.eq."T2Ms") then
        req_atts2 = (/"ts_minlat", "ts_maxlat", "ts_minlon", "ts_maxlon"/)
    end if
    if (field_type0.eq."T2Mz") then
        req_atts2 = (/"ts_minlev", "ts_maxlev", "ts_minlat", "ts_maxlat"/)
    end if
    if (field_type0.eq."T1Ms") then
        req_atts2 = (/"ts_minlat", "ts_maxlat"/)
    end if
    if (field_type0.ne."T0M") then
        exit_if_missing_atts(diag_script_info, \
                             array_append_record(req_atts1, req_atts2, 0))
    end if


    ;; Time averaging
    if (isatt(diag_script_info, "time_avg")) then
        time_avg = diag_script_info@time_avg
        if (all(time_avg.ne.(/"yearly","monthly"/))) then
            error_msg("f", diag_script, "", \
                      "time averaging option " + time_avg + " not yet " + \
                      "implemented")
        end if
    else
        time_avg = "monthly"
    end if

    ;; Area opperation
    if (isatt(diag_script_info, "area_opper")) then
        area_opper = diag_script_info@area_opper
        if (all(area_opper.ne.(/"sum","average"/))) then
            error_msg("f", diag_script, "", \
                      "area opperation option " + area_opper + " not yet " + \
                      "implemented")
        end if
    else
        area_opper = "average"
    end if

    ;; Masking
    if (isatt(diag_script_info, "ts_lsmask")) then
        lsmask = diag_script_info@ts_lsmask
        if(lsmask.eq."land" .or. lsmask.eq."ocean") then
            fm = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
        end if 
    else
        lsmask = "nomask"
    end if

    ;; Determine time range
    start_year = min(stringtoint(models@start_year))
    end_year = max(stringtoint(models@end_year))
    all_years = ispan(start_year, end_year, 1)

    ;; Determine time range
    ;; Create time coordinate
    if (time_avg.eq."monthly") then
        ntime = 12 * (end_year - start_year + 1)
        time = new(ntime, integer)
        do yy = start_year, end_year
            do mm = 1, 12
                time(12 * (yy - start_year) + mm - 1) = 100 * yy + mm
            end do
        end do
    else if (time_avg.eq."yearly") then
        ntime = end_year - start_year + 1
        time = new(ntime, integer)
        time = ispan(start_year, end_year, 1)
    end if
    end if

    ;; Anomaly
    if (isatt(diag_script_info, "ts_anomaly")) then
        anom = diag_script_info@ts_anomaly
        if(anom .eq. "anom") then     
            if(isatt(diag_script_info, "ref_start") .and. \
               isatt(diag_script_info, "ref_end")) then
                 ref_start = diag_script_info@ref_start
                 ref_end   = diag_script_info@ref_end
                 if((ref_start.lt.start_year) .or. \
                    (ref_end.gt.end_year)   .or. \
                    (ref_end.lt.ref_start)) then
                      error_msg("f", diag_script, "", \
                          "period for reference years is not properly defined")
                 end if
            else
                error_msg("f", diag_script, "", \
                    "period for reference years is not defined (needed for anomaly)")
            end if
         
            anom_ref = new((/dim_MOD/), float)
            anom_ref!0 = "model"
            anom_ref&model = project_style(diag_script_info, "annots")
        end if
    else
        anom = "noanom"
    end if

    ;; Create model array
    model_arr = new((/dim_MOD, ntime/), float)
    model_arr!0 = "model"
    model_arr!1 = "time"
    model_arr&model = project_style(diag_script_info, "annots")
    model_arr&time = time
    delete(time)
    delete(ntime)

    ;; Select attributes
    avgstd = project_style(diag_script_info, "avgstd")

    ;; Prepare calculation of standard deviation of models
    if (diag_script_info@multi_model_mean .eq. "y") then
        model_arr_stddev = new((/4, end_year - start_year + 1/), float)
    else  ; initialize anyway, because both are needed as arguments of xy_line
        model_arr_stddev = 0
    end if



    ;; Loop over models
    do imod = 0, dim_MOD - 1

        ;; Read data
        A0 = extract_data(imod, data, -1, 0, 0)

        ;; Apply mask
        if (lsmask.eq."ocean") then  ; exclude land, small island, ice shelf
            lsm  = landsea_mask(fm->LSMASK, A0&lat, A0&lon)
            A0 = mask(A0, lsm.eq.1 .or. lsm.eq.3 .or. lsm.eq.4, False)
            delete(lsm)
        else if (lsmask.eq."land") then  ; exclude ocean, lake
            lsm  = landsea_mask(fm->LSMASK, A0&lat, A0&lon)
            A0 = mask(A0, lsm.eq.0 .or. lsm.eq.2, False)
            delete(lsm)
        end if
        end if

        ;; time, lev, lat, lon
        if (field_type0.eq."T3M") then

            ;; Anomaly
            if(anom .eq. "anom") then
                A0_timavg_ref = time_operations(A0,ref_start,ref_end,"average", "annualclim", True)
                do i=0,dimsizes(A0&time)-1
                  A0(i,:,:,:) = A0(i,:,:,:) - A0_timavg_ref(:,:,:)
                end do
            end if
            
            tmp = area_operations(A0, diag_script_info@ts_minlat, \
                                  diag_script_info@ts_maxlat, \
                                  diag_script_info@ts_minlon, \
                                  diag_script_info@ts_maxlon, area_opper, True)
            if (diag_script_info@ts_minlev.lt.0) then
                minlev = min(A0&plev)
            else
                minlev = 100. * diag_script_info@ts_minlev
            end if
            if (diag_script_info@ts_maxlev.lt.0) then
                maxlev = max(A0&plev)
            else
                maxlev = 100. * diag_script_info@ts_maxlev
            end if
            procmod = dim_avg_Wrap(tmp(:, {maxlev:minlev}))
            delete(tmp)

        end if

        ;; time, lat, lon
        if (field_type0.eq."T2Ms") then

            ;; Anomaly
            if(anom .eq. "anom") then
                A0_timavg_ref = time_operations(A0,ref_start,ref_end,"average", "annualclim", True)
                do i=0,dimsizes(A0&time)-1
                  A0(i,:,:) = A0(i,:,:) - A0_timavg_ref(:,:)
                end do
            end if
            
            procmod = area_operations(A0, diag_script_info@ts_minlat, \
                                      diag_script_info@ts_maxlat, \
                                      diag_script_info@ts_minlon, \
                                      diag_script_info@ts_maxlon, area_opper, \
	                              True)

        end if

        ;; time, lev, lat
        if (field_type0.eq."T2Mz") then
            error_msg("f", diag_script, "", \
                      "field type T2Mz not yet implemented")
        end if

        ;; time, lat
        if (field_type0.eq."T1Ms") then
            error_msg("f", diag_script, "", \
                      "field type T1Ms not yet implemented")
        end if

        ;; time
        if (field_type0.eq."T0M") then

            ;; Anomaly
            if(anom .eq. "anom") then
                A0_timavg_ref = time_operations(A0,ref_start,ref_end,"average", "annualclim", True)
                do i=0,dimsizes(A0&time)-1
                  A0(i) = A0(i) - A0_timavg_ref
                end do
            end if
           
            procmod = A0

        end if

        if (.not.isdefined("procmod")) then
            error_msg("f", diag_script_info, "", "cannot process this " + \
                      "field type " + field_type0)
        end if
        delete(A0)

        ;; Detrend
        if (isatt(diag_script_info, "ts_detrend")) then
            detr = diag_script_info@ts_detrend
            tmp = runave_Wrap(procmod,2,0)
            delete(procmod)
            procmod = tmp
            delete(tmp)
        else
            detr = "nodetr"
        end if

        ;; Calcualte time average
        if (time_avg.ne."monthly") then
            A0_timavg = time_operations(procmod, -1, -1, "average", time_avg, True)
        else
            A0_timavg = procmod
        end if
        delete(procmod)
        procmod = A0_timavg

        ;; Anomaly
        if(anom .eq. "anom") then
             anom_ref(imod) = area_operations(A0_timavg_ref, diag_script_info@ts_minlat, \
                                  diag_script_info@ts_maxlat, \
                                  diag_script_info@ts_minlon, \
                                  diag_script_info@ts_maxlon, area_opper, \
                                  True)
             delete(A0_timavg_ref)
        end if

        ;; Match time coordinate
        if (time_avg.eq."monthly") then
            date = cd_calendar(procmod&time, -1)
        else if (time_avg.eq."yearly") then
            date = procmod&year
        end if
        end if
        idx1 = ind(date(0).eq.model_arr&time)
        idx2 = ind(date(dimsizes(date) - 1).eq.model_arr&time)
        model_arr(imod, idx1:idx2) = (/procmod/)
        if (imod.eq.0) then
            copy_VarAtts(procmod, model_arr)
        end if
        delete(procmod)
        delete(date)

    end do

    ;; Calculate multi model mean of models (stddev would be possible)
    if (diag_script_info@multi_model_mean .eq. "y") then
        ;; See style_$project$.ncl for which data shall
        ;; be included in the statistics
        temp = ind(avgstd .eq. 0)
        do imonth = 0, dimsizes(model_arr&time) - 1
            model_arr_stddev(0, imonth) = \
                dim_avg_n_Wrap(model_arr(temp, imonth), (/0/))
            ;no stddev in the plot
            model_arr_stddev(1, imonth) = model_arr_stddev(0, imonth)
            model_arr_stddev(2, imonth) = model_arr_stddev(0, imonth)
            model_arr_stddev(3, imonth) = model_arr_stddev(0, imonth)
            ;stddev for grey area
            ;model_arr_stddev(1, imonth) = \
            ;    dim_stddev_n_Wrap(model_arr(temp, imonth), (/0/))
            ;model_arr_stddev(2, imonth) = \
            ;    model_arr_stddev(0, imonth) - model_arr_stddev(1, imonth)
            ;model_arr_stddev(3, imonth) = \
            ;    model_arr_stddev(0, imonth) + model_arr_stddev(1, imonth)
        end do
        delete(temp)

    end if


    ;; Convert time coordinate to years (required by the plot routine)
    if (time_avg.eq."monthly") then
        year = model_arr&time / 100
        xmin = min(year)
        xmax = max(year)
        month = model_arr&time - 100 * year
        time = todouble(year + month/12.)
        delete(model_arr&time)
        model_arr&time = time
        delete(time)
        delete(year)
        delete(month)
    else if (time_avg.eq."yearly")
        xmin = min(model_arr&time)
        xmax = max(model_arr&time)
        tmp = todouble(model_arr&time)
        delete(model_arr&time)
        model_arr&time = tmp
        delete(tmp)
    end if
    end if

    ;; Convert units if appropriate
    if (isatt(variable_info, "plot_units")) then
        model_arr = convert_units(model_arr, variable_info@plot_units)
    end if

    ;; Optional output to NetCDF
    if (getenv("ESMValTool_write_netcdf").eq."True") then
        out_path = get_output_dir(get_work_dir(), basename(diag_script))
        out_path = out_path + "tsline_" + var0 + "_" + lsmask + "_" + anom + \
            "_" + detr + "_" + models@experiment(dim_MOD-1) + "_" + \
            start_year + "-" + end_year + ".nc"
        model_arr@ncdf = out_path
        model_arr@diag_script = diag_script
        model_arr@var = var0
        ncdf_outfile = ncdf_write(model_arr, out_path)
    end if

    ;; Define workstation
    outfile = var0 + "_" + lsmask + "_" + anom + "_" + detr + "_" + \
              models@experiment(dim_MOD-1) + "_" + start_year + "-" + end_year
    wks = get_wks(0, diag_script, outfile)

    ;; Set resources
    res = True
    res@trXMinF = xmin
    res@trXMaxF = xmax 
    ;res@trXMaxF = xmax + 0.25 * (xmax - xmin)
    res@tmXBMode = "Explicit"
    if (xmax - xmin.gt.20) then
        res@tmXBValues = ispan(xmin, xmax, 10)
        res@tmXBLabels = ispan(xmin, xmax, 10)
        res@tmXBMinorValues =  ispan(xmin, xmax, 5)
    else
        res@tmXBValues = ispan(xmin, xmax, 5)
        res@tmXBLabels = ispan(xmin, xmax, 5)
        res@tmXBMinorValues =  ispan(xmin, xmax, 1)
    end if
    res@tmXBLabelAngleF = 45
    res@tmXBLabelJust     = "CenterRight"

    if isatt(model_arr,"long_name")
        res@tiMainString = model_arr@long_name
    end if
 
    if (isatt(diag_script_info, "ts_anomaly")) then
        res@tiYAxisString = var0 + " Anomaly" + " [" + model_arr@units + "]"
    else
        res@tiYAxisString = var0 + " [" + model_arr@units + "]"
    end if

    res0 = True
    ref_start = diag_script_info@ref_start
    ref_end = diag_script_info@ref_end
    res0@tiYAxisString = tostring(ref_start) + "-" + tostring(ref_end) +" "+ var0 + " Mean" + " [" + model_arr@units + "]"

    if(anom .eq. "anom") then
        xy_line_anom(wks, anom_ref, model_arr, model_arr&time, model_arr_stddev, ref_start, ref_end, res, res0, False)
    else
        xy_line(wks, model_arr, model_arr&time, model_arr_stddev, res, False)
    end if

    info_output(" gv " + wks@fullname + " &", verbosity, 1)

    ;; References
    write_references(diag_script, "A_cion_ir", (/"A_righ_ma","A_wenz_sa"/), "D_0001", \
                     "", "P_esmval")
    info_output("<<<<<<<< Leaving " + diag_script, verbosity, 4)

end
