;;#############################################################################
;; PROCEDURE FOR THE ANNUAL CYCLE PLOT OF THE PERFORMANCE METRICS
;; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################

procedure perfmetrics_ptype_script()
begin

    ;; Define output array
    if (diag_script_info@time_avg.eq."seasonalclim") then
        var_all = new((/dim_MOD, 4, 2/), float)
        var_all!1 = "season"
        var_all&season = (/"DJF", "MAM", "JJA", "SON"/)
    else if (diag_script_info@time_avg.eq."monthlyclim") then
        var_all = new((/dim_MOD, 12, 2/), float)
        var_all!1 = "month"
        var_all&month = (/"J", "F", "M", "A", "M", "J",\
                          "J", "A", "S", "O", "N", "D"/)
    else
        error_msg("f", diag_script, "", "time_avg option " + \
                  diag_script_info@time_avg + \
                  " not compatible with plot_type cycle")
    end if
    end if
    var_all!0 = "model"
    var_all!2 = "statistic"
    var_all&model = modelnames
    var_all&statistic = (/"mean", "stddev"/)

    ;; Attach attributes
    var_all@var = var0
    var_all@diag_script = (/diag_script/)
    copy_VarAtts(diag_script_info, var_all)
    var_all@ref_model = variable_info@reference_model

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        log_debug("Processing " + modelnames(imod))

        ;; Set path for saving processed data ;; FIX-ME add preproc_id
        fullpath = config_user_info@work_dir + "/"
        system("mkdir -p " + fullpath)
        fname = basename(systemfunc("basename " + model_info@filename(imod)))
        fname = fname + "_" + basename(diag_script) + "_cycle"
        fname = fname + "_" + diag_script_info@time_avg
        if (isatt(diag_script_info, "level")) then
            fname = fname + "_" + diag_script_info@level  ; FIX-ME
        end if
        fname = fname + "_" + str_sub_str(diag_script_info@region, " ", "")
        procpath = fullpath + fname + ".nc"

        ;; Determine start/end year
        start_year = get_start_year(imod)
        end_year = get_end_year(imod)

        ;; Read data
        var = read_data(imod, var0, field_type0)
        dnames = getVarDimNames(var)

        ;; Extract region and average over latitude and longitude
        if (any(field_type0.eq.(/"T0M", "T0Ms", "T1M"/))) then
            var_reg = var
        else if (any(dnames.eq."lat") .and. any(dnames.eq."lon")) then
            var_reg = area_operations(var, region(0), region(1), \
                                      region(2), region(3), "average", True)
        else
            error_msg("f", diag_script, "", "dimensionality not " + \
                      "implemented")  ; FIXME
        end if
        end if
        delete(var)

        ;; Calculate time average
        var_avg = time_operations(var_reg, start_year, end_year, "average", \
                                  diag_script_info@time_avg, True)

        ;; Calculate time standard deviation (with lower/upper bounds)
        if (start_year.lt.end_year) then
            var_std = time_operations(var_reg, start_year, end_year, \
                                      "stddev", diag_script_info@time_avg, \
                                      True)
        else
            var_std = 0.
        end if
        delete(var_reg)

        ;; Write output
        var_all(imod, :, 0) = var_avg
        var_all(imod, :, 1) = var_std
        delete(var_avg)
        delete(var_std)
        if (config_user_info@write_netcdf.eq."True") then
            var_all@ncdf = procpath
            ncdf_outfile = ncdf_write(var_all(imod, :, :), procpath)
        end if

    end do

    ;; Plotting
    if (.not.diag_script_info@draw_plots) then
        return
    end if

    ;; Convert units for plotting (if required)
    if (isatt(diag_script_info, "plot_units")) then
        var_all = convert_units(var_all, diag_script_info@plot_units)
    end if

    ;; Annotation and file names
    region_name = ""  ; priority 3
    location = ""  ; priority 3
    if (isatt(diag_script_info, "aux_info")) then
        region_name = "_" + diag_script_info@aux_info  ; priority 2
        location = " - " + diag_script_info@aux_info  ; priority 2
    end if
    if (isatt(diag_script_info, "location")) then
        location = " - " + diag_script_info@location  ; priority 1
    end if
    if (isatt(diag_script_info, "region")) then
        region_name = "_" + region@name  ; priority 1
        location = " - " + diag_script_info@region
    end if
    if (isatt(var_all, "long_name")) then
        varstring = var_all@long_name  ; priority 3
    end if
    if (isatt(var_all, "short_name")) then
        varstring = var_all@short_name  ; priority 2
    end if
    if (isStrSubset(var0, "vmr").and.isStrSubset(var0, "_")) then
        varstring = var0  ; priority 1
    end if
    altitude = ""  ; priority 2
    if (any(field_type0.eq.(/"T0M", "T0Ms"/))) then
        if (isatt(diag_script_info, "altitude")) then
            altitude = " - " + diag_script_info@altitude  ; priority 1
        end if
    else
        if (isatt(diag_script_info, "level")) then  ; FIX-ME
            if (diag_script_info@level.ne."all") then
                altitude = " - " + diag_script_info@level + " hPa"
            end if
        end if
    end if
    outfile = config_user_info@plot_dir + var0 + "_cycle_" + \
        diag_script_info@time_avg + "_" + region_name

    wks = gsn_open_wks(file_type, outfile)
    wks@legendfile = outfile + "_legend"
    var_all@res_tiMainString = varstring
    var_all@res_tiMainString = var_all@res_tiMainString + \
        location + altitude
    plot = cycle_plot(wks, var_all, var0)
    draw(plot)
    frame(wks)

    log_info(" gv " + outfile + "." + file_type)

end
