load "./interface_data/ncl.interface"

undef("check_params")
procedure check_params(plot_type,\
                       field_types,\
                       diag_script_info,\
                       rank:integer,\
                       refModel,\
                       latrange[2]:float,\
                       lonrange[2]:float,\
                       yMax_lev:float,\
                       yMin_lev:float,\
                       xMax_year:float,\
                       xMin_year:float,\
                       anom:logical,\
                       smooth:logical,\
                       trend:logical,\
                       multimean:logical,\
                       coverage:logical)
begin
    ;; check field type is (might eventually be) supported;
    ;;  only TO2Ms at present
    plot_fields = (/ "TO2Ms"/)
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering check_params", verbosity, 4)

    ind_f = ind(plot_fields .eq. field_types)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_types\
                                        + " for plot_type "\
                                        + plot_type)
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)
    ;; ******************************************************************
    ;; check the input parameters
    refModel = "NULL"
    if (isatt(diag_script_info, "refModel")) then
        refModel = diag_script_info@refModel
    end if
    if (rank .gt. 1) then
        if(isatt(diag_script_info, "ts_latrange")) then
            latrange = diag_script_info@ts_latrange
        else
            print("fatal:plot requires variable attribute 'ts_latrange' defined in varDiag.cfg")
            status_exit(1)
        end if

        if(isatt(diag_script_info, "ts_lonrange")) then
            lonrange = diag_script_info@ts_lonrange
        else
            if (isStrSubset(field_type0, "z") .eq. False) then
                print("fatal:plot requires variable attribute 'ts_lonrange'  defined in varDiag.cfg ")
                status_exit(1)
            end if
        end if
    end if

    yMax_lev = -1.e20
    yMin_lev = 1.e20
    ;; optional input parameters
    ;; axes limits
    if(isatt(diag_script_info, "ts_ymax")) then
        yMax_lev = diag_script_info@ts_ymax
    end if
    if(isatt(diag_script_info, "ts_ymin")) then
        yMin_lev = diag_script_info@ts_ymin
    end if
    if(isatt(diag_script_info, "ts_max_year")) then
        xMax_year = diag_script_info@ts_max_year
    end if
    if(isatt(diag_script_info, "ts_min_year")) then
        xMin_year = diag_script_info@ts_min_year
    end if
    info_output("Limits: " + xMin_year + xMax_year, verbosity, 4)

    ;; plot anomaly flag
    if(isatt(diag_script_info, "ts_anom")) then
        anom = diag_script_info@ts_anom
    else
        anom = True
    end if

    ;; plot trend flag
    if(isatt(diag_script_info, "ts_trend")) then
        trend = diag_script_info@ts_trend
    end if

    ;; plot smoothed data flag
    if(isatt(diag_script_info, "ts_smooth")) then
        smooth = diag_script_info@ts_smooth
    end if

    ;; plot mean of models flag
    if(isatt(diag_script_info, "ts_multimean")) then
        multimean = diag_script_info@ts_multimean
    end if

    ;; include area covered by observations in plot
    if(isatt(diag_script_info, "ts_coverage")) then
        coverage = diag_script_info@ts_coverage
    end if
    info_output(">>>>>>>> Leaving check_params", verbosity, 4)
end
