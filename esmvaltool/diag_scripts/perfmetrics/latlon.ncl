;;#############################################################################
;; PROCEDURE FOR THE LATLON PLOT OF THE PERFORMANCE METRICS
;; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################

procedure perfmetrics_ptype_script()
begin

    ;; Set model indexes, with reference model first
    modidx = ispan(0, dim_MOD - 1, 1)
    modidx := array_append_record(ref_ind, modidx(ind(modidx.ne.ref_ind)), 0)

    ;; Loop over models
    do ii = 0, dimsizes(modidx) - 1
        
        imod = modidx(ii)

        log_debug("Processing " + modelnames(imod))

        ;; Set path for saving processed data
        fullpath = config_user_info@work_dir + "/"
        system("mkdir -p " + fullpath)
        fname = basename(systemfunc("basename " + model_info@filename(imod)))
        fname = fname + "_" + basename(diag_script) + "_latlon"
        fname = fname + "_" + diag_script_info@time_avg
        procpath = fullpath + fname + "_Global.nc"

        ;; Determine start/end year
        start_year = get_start_year(imod)
        end_year = get_end_year(imod)

        ;; Read data
        var = read_data(imod, var0, field_type0)
        dnames = getVarDimNames(var)

        ;; Calculate time average and standard deviation w.r.t.
        ;; interannual variability
        var_avg = time_operations(var, start_year, end_year, \
                                  "average", diag_script_info@time_avg, True)
        var_std = interannual_variability(var, start_year, end_year, \
                                          diag_script_info@time_avg)
        delete(var)

        ;; Write output
        var_all = new(array_append_record(dimsizes(var_avg), 2, 0), float)
        var_all!2 = "statistic"
        var_all&statistic = (/"mean", "stddev"/)
        var_all(:, :, 0) = var_avg
        var_all(:, :, 1) = var_std
        var_all@var = var0
        var_all@diag_script = (/diag_script/)
        copy_VarAtts(diag_script_info, var_all)
        var_all@ref_model = variable_info@reference_model
        if (config_user_info@write_netcdf.eq."True") then
            var_all@ncdf = procpath
            ncdf_outfile = ncdf_write(var_all, procpath)
        end if

        ;; Store reference model
        if (imod.eq.ref_ind) then
            var_ref = var_all
        end if

        ;; Plotting
        if (.not.diag_script_info@draw_plots) then
            continue
        end if

        ;; Convert units for plotting (if required)
        if (isatt(diag_script_info, "plot_units")) then
            var_all = convert_units(var_all, diag_script_info@plot_units)
        end if

        ;; Absolute plot
        outfile = config_user_info@plot_dir + modelnames(imod) + "_" + var0 + \
            "_latlon_" + diag_script_info@time_avg + "_Glob"
        plot_var = var_all(:, :, 0)
        plot_var@projection = diag_script_info@projection
        plot_var@res_cnLinesOn = False
        plot_var@res_cnLevelSelectionMode = "ExplicitLevels"
        plot_var@res_tiMainString = modelnames(imod)
        plot_var@res_cnLevels = diag_script_info@abs_levs
        plot_var@res_gsnLeftString = \
            plot_var@long_name + " [" + format_units(plot_var@units) + "]"
        if (diag_script_info@show_global_avg) then  ;; append global average
            gavg = \
                area_operations(plot_var, -90., 90., 0., 360., "average", True)
            plot_var@res_gsnRightString = \
                sprintf("%5.2f", gavg) + " " + format_units(plot_var@units)
        end if
        wks = gsn_open_wks(file_type, outfile)
        gsn_define_colormap(wks, diag_script_info@latlon_cmap)
        plot = contour_map_ce(wks, plot_var, var0)
        draw(plot)
        frame(wks)
        log_info(" gv " + outfile + "." + file_type)
        delete(plot_var@res_cnLevels)
        delete(wks)
        delete(plot)

        ;; Difference plot to the reference
        if (.not.diag_script_info@plot_diff .or. imod.eq.ref_ind) then
            continue
        end if
        plot_var = var_all(:, :, 0) - var_ref(:, :, 0)
        outfile = config_user_info@plot_dir + modelnames(imod) + "-" + \
            modelnames(ref_ind) + "_" + var0 + "_latlon_" + \
            diag_script_info@time_avg + "_Glob"
        plot_var@res_tiMainString = \
            modelnames(imod) + " - " + modelnames(ref_ind)
        if (isatt(diag_script_info, "t_test")) then
            plot_var@res_gsnLeftString = plot_var@res_gsnLeftString + " - " + \
            sprinti("%2i", toint(100 * diag_script_info@conf_level)) + "% c.l."
        end if
        plot_var@res_cnLevelSelectionMode = "ExplicitLevels"
        plot_var@res_cnLevels = diag_script_info@diff_levs
        plot_var@res_cnMissingValFillColor = "gray70"
        plot_var@res_gsnSpreadColorEnd = -2
        plot_var@projection = diag_script_info@projection
        if (diag_script_info@show_global_avg) then  ;; append global average
            gavg = \
                area_operations(plot_var, -90., 90., 0., 360., "average", True)
            plot_var@res_gsnRightString = \
                sprintf("%5.2f", gavg) + " " + format_units(plot_var@units)
        end if
        wks = gsn_open_wks(file_type, outfile)
        gsn_merge_colormaps(wks, "temp_19lev", "gray70")
        plot = contour_map_ce(wks, plot_var, var0)
        delete(plot_var@res_cnLevels)

        ;; Apply t-test (assuming different population variances)
        if (.not.diag_script_info@t_test) then
            continue
        end if
        x1 = var_all(:, :, 0)
        x2 = var_ref(:, :, 0)
        s1 = var_all(:, :, 1) ^ 2
        s2 = var_ref(:, :, 1) ^ 2
        n1 = get_end_year(imod) - get_start_year(imod) + 1
        n2 = get_end_year(ref_ind) - get_start_year(ref_ind) + 1
        prob = ttest(x1, s1, n1, x2, s2, n2, True, False)
        mask_var = 1. - prob
        copy_VarCoords(plot_var, mask_var)
        delete(prob)

        ;; Mask non-significant values
        res2 = True
        res2@cnLevelSelectionMode = "ExplicitLevels"
        res2@cnLevels = (/0.0, diag_script_info@conf_level, 1.0/)
        res2@cnInfoLabelOn = False
        res2@cnLinesOn = False
        res2@cnLineLabelsOn = False
        res2@gsnDraw = False
        res2@gsnFrame = False
        plot2 = gsn_csm_contour(wks, mask_var, res2)
        opt = True
        if (diag_script_info@stippling) then
            opt@gsnShadeFillType = "pattern"
            opt@gsnShadeLow = 17
        else
            opt@gsnShadeFillType = "color"
            opt@gsnShadeLow = "gray70"
        end if
        plot2 = \
            gsn_contour_shade(plot2, diag_script_info@conf_level, 999., opt)
        overlay(plot, plot2)
        draw(plot)
        frame(wks)
        log_info(" gv " + outfile + "." + file_type)
        delete(wks)
        delete(plot)
        delete(plot2)
        delete(mask_var)
        delete(plot_var)

    end do

end
