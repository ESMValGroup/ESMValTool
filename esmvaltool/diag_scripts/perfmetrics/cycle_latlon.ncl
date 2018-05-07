;;#############################################################################
;; PROCEDURE FOR THE CYCLE-LATLON PLOT OF THE PERFORMANCE METRICS
;; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################

procedure perfmetrics_ptype_script()
begin

    ;; Set model indexes, with reference and alternative models first
    modidx = ispan(0, dim_MOD - 1, 1)
    if (l_altern) then
        modidx := array_append_record((/ref_ind, alt_ind/), \
                                      modidx(ind(modidx.ne.ref_ind .and. \
                                                 modidx.ne.alt_ind)), 0)
    else
        modidx := array_append_record(ref_ind, \
                                      modidx(ind(modidx.ne.ref_ind)), 0)
    end if

    ;; Loop over models
    do ii = 0, dimsizes(modidx) - 1
        
        imod = modidx(ii)

        log_debug("Processing " + modelnames(imod))

        ;; Set path for saving processed data ;; FIX-ME add preproc_id
        fullpath = config_user_info@work_dir + "/"
        system("mkdir -p " + fullpath)
        fname = basename(systemfunc("basename " + model_info@filename(imod)))
        fname = fname + "_" + basename(diag_script) + "_cycle_latlon"
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

        ;; Extract region
        var_reg = area_operations(var, region(0), region(1), \
                                  region(2), region(3), "extract", True)
        delete(var)

        ;; Calculate time average
        var_avg = time_operations(var_reg, start_year, end_year, \
                                  "average", diag_script_info@time_avg, True)
        delete(var_reg)

        ;; Standard deviation calculation for this ptype is not implemented
        ;; yet in statistics.ncl (but is anyway not required here)

        ;; Write output
        var_all = new(array_append_record(dimsizes(var_avg), 2, 0), float)
        var_all!3 = "statistic"
        var_all&statistic = (/"mean", "stddev"/)
        var_all(:, :, :, 0) = var_avg
        var_all@var = var0
        var_all@diag_script = (/diag_script/)
        copy_VarAtts(diag_script_info, var_all)
        var_all@ref_model = variable_info@reference_model
        if (config_user_info@write_netcdf.eq."True") then
            var_all@ncdf = procpath
            ncdf_outfile = ncdf_write(var_all, procpath)
        end if

        ;; Calculate grading
        if (.not.diag_script_info@calc_grading) then
            continue
        end if

        ;; Store reference model
        if (imod.eq.ref_ind) then
            var_ref = var_all
        end if

        ;; Store alternative model
        if (l_altern) then
            if (imod.eq.alt_ind) then
                var_alt = var_all
            end if
        end if

        ;; Loop over requested metrics
        do met = 0, nmetrics - 1

            if (diag_script_info@metric(met).eq."taylor") then
                locidx = ind(modelnames(imod).eq.taylor&models)
                if (ismissing(locidx)) then
                    continue
                end if
                taylor(0, locidx, 0) = \
                    calculate_metric(var_all(:, :, :, 0), \
                                     var_ref(:, :, :, 0), \
                                     "stddev_ratio")
                taylor(0, locidx, 1) = \
                    calculate_metric(var_all(:, :, :, 0), \
                                     var_ref(:, :, :, 0), \
                                     "correlation")
            else
                locidx = ind(modelnames(imod).eq.grading&models)
                if (ismissing(locidx)) then
                    continue
                end if
                grading(met, 0, locidx, 0) = \
                    calculate_metric(var_all(:, :, :, 0), \
                                     var_ref(:, :, :, 0), \
                                     diag_script_info@metric(met))
                if (l_altern) then
                    grading(met, 0, locidx, 1) = \
                        calculate_metric(var_all(:, :, :, 0), \
                                         var_alt(:, :, :, 0), \
                                         diag_script_info@metric(met))
                end if
            end if

        end do

    end do

end
