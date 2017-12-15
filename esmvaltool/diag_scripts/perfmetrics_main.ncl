;;#############################################################################
;; MAIN SCRIPT FOR PERFORMANCE METRICS
;; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Calculates and (optionally) plots annual/seasonal cycles, zonal means,
;;    lat-lon fields and time-lat-lon fields from input T3M or T2Ms data.
;;    The calculated fields can be also plotted as difference w.r.t. a given
;;    reference model. It also calculates grading and taylor metrics.
;;    Input data have to be regridded to a common grid in the preprocessor.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    plot_type: cycle (time), zonal (plev, lat), latlon (lat, lon) or
;;               cycle_latlon (time, lat, lon)
;;    time_avg: type of time average (see time_operations in 
;;              diag_scripts/lib/ncl/statistics.ncl)
;;    region: selected region (see select_region in 
;;            diag_scripts/lib/ncl/latlon.ncl)
;;    styleset (for cycle): as in lib/ncl/style.ncl functions
;;    plot_stddev (for cycle): plot standard deviation
;;    legend_outside (for cycle): save legend in a separate file
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    t_test (for zonal and latlon): calculate t-test in difference plots
;;                                   (default: False)
;;    conf_level (for zonal and latlon): confidence level for the t-test
;;                                       (default: False)
;;    range_option: time range selection option (default: 0)
;;    projection: map projection for lat-lon plots (default:
;;                CylindricalEquidistant)
;;    draw_plots: draw plots (default: True)
;;    plot_diff: draw difference plots (default: False)
;;    calc_grading: calculate grading (default: False)
;;    stippling: use stippling to mark stat. significant differences (default:
;;               False = mask out non-significant differences in gray)
;;    show_global_avg: diplay global avaerage as right string on lat-lon plots
;;                     (default: False)
;;    metric: grading metric (if calc_grading is True)
;;    normalization: metric normalization (for RMSD and BIAS metrics)
;;    aux_info: additional information for stations data
;;    location: additional information for stations data (location)
;;    altitude: additional information for stations data (altitude)
;;
;; Required variable_info attributes (variable specific)
;;    map_ref_Levels: contour levels for the latlon plot
;;    map_diff_Levels: contour levels for the latlon difference plot
;;    zonal_ref_Levels: contour levels for the zonal plot
;;    zonal_diff_Levels: contour levels for the zonal difference plot
;;
;; Optional variable_info attributes (variable specific)
;;    zonal_ref_cmap: color table for the zonal plot
;;    zonal_ymin: minimum pressure for the zonal plots
;;    global_avg: attach global field average as right string on the plot
;;    plot_units: plotting units (if different from standard CMOR units)
;;
;; Caveats
;;    The reference model must be specified in the main namelist in the
;;    variables dictionary as ref_model key. An alternative model can be also
;;    given, using a python list.
;;
;; Modification history
;;    20171215-A_righ_ma: merged with perfmetrics_grading/taylor.ncl
;;    20171124-A_righ_ma: completely revised to adapt it to the new backend
;;                        (level selection, regridding and masking now done
;;                        by the python preprocessor)
;;    20161220-A_laue_ax: added option to set map projection for lat-lon plots
;;                        (diag_script_info@projection)
;;                        added option to choose how to plot t-test results:
;;                        stippling or masking out in gray (lat-lon plots only)
;;    20161019-A_laue_ax: changed plotting of t-test results:
;;                        now stippling significant grid cells (old version:
;;                        masking out non-significant values in gray)
;;    20160628-A_righ_ma: moving ref_model specification from cfg- files to
;;                        main namelist
;;    20160628-A_senf_da: added regridding for irregular grids (ESMF_regrid).
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code.
;;    20151013-A_righ_ma: fixed t-test mask in lat-lon difference plots.
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc.).
;;    20150119-A-gott_kl: removed "grid", "region" from req_atts
;;                        (for T2Ms vmrco).
;;    20150113-A_gott_kl: reconciled generalised regridding with T1* & T0*
;;    20140905-A_righ_ma: consistent regridding and missing values mask.
;;    20140701-A_gott_kl: Adapted for T1M.
;;    20140630-A_gott_kl: Adapted for T0Ms.
;;    20131203-A_fran_fr: written.
;;
;;#############################################################################

load "$ESMValTool_interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"

load "./plot_scripts/ncl/xy_line.ncl"
load "./plot_scripts/ncl/zonalmean_profile.ncl"
load "./plot_scripts/ncl/contour_maps.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Write references
    write_references(diag_script, \
                     "A_fran_fr", \
                     (/"A_righ_ma", "A_eyri_ve", "A_gott_kl", "A_senf_da"/), \
                     (/"D_righi15gmd", "D_gleckler08jgr"/), \
                     (/"E_ncep", "E_erainterim", "E_airs", \
                       "E_ceresebaf", "E_srb"/), \
                     (/"P_embrace", "P_esmval"/))
end

begin
    vardeffile = "$ESMValTool_interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin

    ;; ### CHECK CONFIGURATION SETTINGS AND SET DEFAULTS ######################

    ;; Check required diag_script_info attributes
    exit_if_missing_atts(diag_script_info, \
                         (/"plot_type", "time_avg", "region"/))

    ;; Store required attributes
    ptype = diag_script_info@plot_type
    time_avg = diag_script_info@time_avg
    region = select_region(diag_script_info@region)

    if (all(ptype.ne.(/"cycle", "zonal", "latlon", "cycle_latlon"/))) then
        error_msg("f", diag_script, "", "plot_type " + ptype + " is not a " + \
                  "supported plot_type in this diagnostic")
    end if

    ;; Check for plot-type specific settings
    if (ptype.eq."cycle") then
        exit_if_missing_atts(diag_script_info, \
                             (/"legend_outside", "styleset", "plot_stddev"/))
    end if

    ;; Check valid field
    if ((ptype.eq."zonal" .and. \
         all(field_type0.ne.(/"T3M", "T2Mz"/))) .or. \
        (ptype.eq."cycle" .and. \
         all(field_type0.ne.(/"T3M", "T2Ms", "T1M", "T0Ms"/))) .or. \
        (ptype.eq."latlon" .and. \
         all(field_type0.ne.(/"T3M", "T2Ms"/))) .or. \
        (ptype.eq."cycle_latlon" .and. \
         all(field_type0.ne.(/"T3M", "T2Ms"/)))) then
        error_msg("f", diag_script, "", "input field " + field_type0 + \
                  " is not compatible with plot_type " + ptype)
    end if

    ;; Set default values for non-required diag_script_info attributes
    set_default_att(diag_script_info, "range_option", 0)
    set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
    set_default_att(diag_script_info, "draw_plots", True)
    set_default_att(diag_script_info, "plot_diff", False)
    set_default_att(diag_script_info, "calc_grading", False)
    set_default_att(diag_script_info, "stippling", False)
    set_default_att(diag_script_info, "t_test", False)
    set_default_att(diag_script_info, "show_global_avg", False)

    ;; Check consistency of diff plots settings
    if (diag_script_info@t_test .and. .not.diag_script_info@plot_diff) then
        error_msg("f", diag_script, "", \
                  "plot_diff must be True to apply t-test")
    end if
    if (diag_script_info@t_test .and. .not.diag_script_info@conf_level) then
        error_msg("f", diag_script, "", \
                  "conf_level must be specified to apply t-test")
    end if

    ;; Check metric
    if (diag_script_info@calc_grading) then
        exit_if_missing_atts(diag_script_info, (/"metric"/))
        if (diag_script_info@metric.ne."taylor") then
            exit_if_missing_atts(diag_script_info, (/"normalization"/))
        end if
    end if

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Unique names for models
    modelnames = project_style(diag_script_info, "annots")

    ;; Check for reference model definition
    if (var_attr_ref(0).eq."None") then
        error_msg("f", diag_script, "", "no reference model is specified")
    end if

    ;; Set index of the reference (and alternative) model
    ref_inds = get_ref_model_idx(models, var_attr_ref)
    ref_ind = ref_inds(0)
    if (dimsizes(ref_inds).gt.1) then
        alt_ind = ref_inds(1)
        l_altern = True
    end if

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Plot file type
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; ### DEFINE GLOBAL ARRAYS FOR COLLECTING MODEL DATA #####################

    ;; Read common grid from the first model
    ;; (all models are expected to be regridded to a common grid)
    data = read_data(0, var0, field_type0)
    tmp = extract_data_extended(0, var0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(tmp))
    if (rank.eq.4) then
        grid = tmp(0, :, :, :)  ; plev, lat, lon
    end if
    if (rank.eq.3) then
        grid = tmp(0, :, :)  ; lat, lon
    end if
    delete(data)
    delete(tmp)
    delete(rank)

    ;; Case cycle: (model, time, stats)
    if (ptype.eq."cycle") then
        if (time_avg.eq."seasonalclim") then
            var_all = new((/dim_MOD, 4, 2/), float)
            var_all!1 = "season"
            var_all&season = (/"DJF", "MAM", "JJA", "SON"/)
        else if (time_avg.eq."monthlyclim") then
            var_all = new((/dim_MOD, 12, 2/), float)
            var_all!1 = "month"
            var_all&month = (/"J", "F", "M", "A", "M", "J",\
                              "J", "A", "S", "O", "N", "D"/)
        else
            error_msg("f", diag_script, "", "time_avg option " + time_avg + \
                      " not compatible with plot_type " + ptype)
        end if
        end if
        var_all!0 = "model"
        var_all!2 = "statistic"
        var_all&model = modelnames
        var_all&statistic = (/"mean", "stddev"/)
        if (diag_script_infoq@plot_diff) then
            var_diff = var_all(:, :, 0)
        end if
    end if

    ;; Case cycle_latlon: (model, time, lat, lon, stats)
    if (ptype.eq."cycle_latlon") then
        sgrid = area_operations(grid, region(0), region(1), region(2), \
                                region(3), "extract", False)
        if (time_avg.eq."seasonalclim") then
            var_all = new((/dim_MOD, 4, dimsizes(sgrid&lat), \
                          dimsizes(sgrid&lon), 2/), float)
            var_all!1 = "season"
            var_all&season = (/"DJF", "MAM", "JJA", "SON"/)
        else if (time_avg.eq."monthlyclim") then
            var_all = new((/dim_MOD, 12, dimsizes(sgrid&lat), \
                          dimsizes(sgrid&lon), 2/), float)
            var_all!1 = "month"
            var_all&month = (/"J", "F", "M", "A", "M", "J",\
                              "J", "A", "S", "O", "N", "D"/)
        else
            error_msg("f", diag_script, "", "time_avg option " + time_avg + \
                      " not compatible with plot_type " + ptype)
        end if
        end if
        var_all!0 = "model"
        var_all!2 = "lat"
        var_all!3 = "lon"
        var_all!4 = "statistic"
        var_all&model = modelnames
        var_all&lat = sgrid&lat
        var_all&lon = sgrid&lon
        var_all&statistic = (/"mean", "stddev"/)
        if (diag_script_infoq@plot_diff) then
            var_diff = var_all(:, :, :, :, 0)
        end if
        delete(sgrid)
    end if

    ;; Case zonal: (model, level, lat, stats)
    if (ptype.eq."zonal") then
        var_all = new((/dim_MOD, dimsizes(grid&plev), \
                      dimsizes(grid&lat), 2/), float)
        var_all!0 = "model"
        var_all!1 = "plev"
        var_all!2 = "lat"
        var_all!3 = "statistic"
        var_all&model = modelnames
        var_all&plev = grid&plev
        var_all&lat = grid&lat
        var_all&statistic = (/"mean", "stddev"/)
        if (diag_script_infoq@plot_diff) then
            var_diff = var_all(:, :, :, 0)
        end if
    end if

    ;; Case latlon: (model, lat, lon, stats)
    if (ptype.eq."latlon") then
        var_all = new((/dim_MOD, dimsizes(grid&lat), \
                      dimsizes(grid&lon), 2/), float)
        var_all!0 = "model"
        var_all!1 = "lat"
        var_all!2 = "lon"
        var_all!3 = "statistic"
        var_all&model = modelnames
        var_all&lat = grid&lat
        var_all&lon = grid&lon
        var_all&statistic = (/"mean", "stddev"/)
        if (diag_script_infoq@plot_diff) then
            var_diff = var_all(:, :, :, 0)
        end if
    end if

    ; If needed, create array for statistically significant grid cells
    if (diag_script_info@t_test) then
        var_sig = var_diff
    end if

    ;; Variable rank: cycle:3, zonal/latlon:4, cycle_latlon:5
    irank = dimsizes(dimsizes(var_all))

    ;; Attach attributes
    var_all@var = var0
    var_all@diag_script = (/diag_script/)
    copy_VarAtts(diag_script_info, var_all)
    var_all@ref_model = str_get_field(var_attr_ref(0), 1, ",")

    ;; ### MAIN LOOP OVER MODELS AND DATA PROCESSING ##########################

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + modelnames(imod), verbosity, 2)

        ;; Set path for saving processed data ;; FIX-ME add preproc_id
        fullpath = getenv("ESMValTool_work_dir") + "/" + \
            basename(diag_script) + "/"
        if (.not.isfilepresent_esmval(fullpath)) then
            system("mkdir -p " + fullpath)
        end if
        fname = basename(systemfunc("basename " + \
                         interface_get_fullpath(var0, field_type0, imod)))
        fname = fname + "_" + basename(diag_script)
        fname = fname + "_" + diag_script_info@plot_type
        fname = fname + "_" + diag_script_info@time_avg
        if (isatt(diag_script_info, "level")) then
            fname = fname + "_" + diag_script_info@level ;; FIX-ME
        end if
        fname = fname + "_" + str_sub_str(diag_script_info@region, " ", "")
        procpath = fullpath + fname + ".nc"

;;DELETE        ;; Do not process data if already processed
;;DELETE        if (isfilepresent_esmval(procpath)) then
;;DELETE            info_output("No recalculation!", verbosity, 1)
;;DELETE            info_output("Reading in file = " + procpath, verbosity, 1)
;;DELETE            if (irank.eq.5) then
;;DELETE                var_all(imod, :, :, :, :) = ncdf_read(procpath, var0)
;;DELETE            end if
;;DELETE            if (irank.eq.4) then
;;DELETE                var_all(imod, :, :, :) = ncdf_read(procpath, var0)
;;DELETE            end if
;;DELETE            if (irank.eq.3) then
;;DELETE                var_all(imod, :, :) = ncdf_read(procpath, var0)
;;DELETE            end if
;;DELETE            continue
;;DELETE        end if

        ;; Determine start/end year
        start_year = get_start_year(imod)
        end_year = get_end_year(imod)

        ;; Read data
        data = read_data(imod, var0, field_type0)
        var = extract_data_extended(imod, var0, data, -1, 0, 0)
        dnames = getVarDimNames(var)

        ;; CASE: cycle ########################################################
        if (ptype.eq."cycle") then

            ;; Extract region and average over latitude and longitude
            if (any(field_type0.eq.(/"T0M", "T0Ms", "T1M"/))) then
                var_reg = var
            else if (any(dnames.eq."lat") .and. any(dnames.eq."lon")) then
                var_reg = area_operations(var, region(0), region(1), \
                                          region(2), region(3), "average", \
                                          True)
            else
                error_msg("f", diag_script, "", "dimensionality not " + \
                          "implemented")  ; FIXME
            end if
            end if
            delete(var)

            ;; Calculate time average
            var_avg = time_operations(var_reg, start_year, end_year, \
                                      "average", time_avg, True)

            ;; Calculate time standard deviation (with lower/upper bounds)
            if (start_year.lt.end_year) then
                var_std = time_operations(var_reg, start_year, end_year, \
                                          "stddev", time_avg, True)
            else
                var_std = 0.
            end if
            delete(var_reg)

            ;; Assign to global array
            var_all(imod, :, 0) = var_avg
            var_all(imod, :, 1) = var_std
            var = var_all(imod, :, :)
            delete(var_avg)
            delete(var_std)

        end if

        ;; CASE: cycle_latlon #################################################
        if (ptype.eq."cycle_latlon") then

            ;; Extract region
            var_reg = area_operations(var, region(0), region(1), \
                                      region(2), region(3), "extract", True)
            delete(var)

            ;; Calculate time average
            var_avg = time_operations(var_reg, start_year, end_year, \
                                      "average", time_avg, True)
            delete(var_reg)

            ;; Standard deviation calculation for this rank is not implemented
            ;; yet in statistics.ncl (but is anyway not required here)

            ;; Assign to global array
            var_all(imod, :, :, :, 0) = var_avg
            delete(var_avg)
            var = var_all(imod, :, :, :, :)

        end if

        ;; CASE: zonal ########################################################
        if (ptype.eq."zonal") then

            ;; Calculate zonal mean
            var_zon = dim_avg_Wrap(var)
            delete(var)

            ;; Calculate time average and standard deviation w.r.t.
            ;; interannual variability
            var_avg = time_operations(var_zon, start_year, end_year, \
                                      "average", time_avg, True)
            var_std = interannual_variability(var_zon, start_year, end_year, \
                                              time_avg)
            delete(var_zon)

            ;; Assign to global array
            var_all(imod, :, :, 0) = var_avg
            var_all(imod, :, :, 1) = var_std
            var = var_all(imod, :, :, :)
            delete(var_avg)
            delete(var_std)

        end if

        ;; CASE: latlon #######################################################
        if (ptype.eq."latlon") then

            ;; Calculate time average and standard deviation w.r.t.
            ;; interannual variability
            var_avg = time_operations(var, start_year, end_year, \
                                      "average", time_avg, True)
            var_std = interannual_variability(var, start_year, end_year, \
                                              time_avg)
            delete(var)

            ;; Assign to global array
            var_all(imod, :, :, 0) = var_avg
            var_all(imod, :, :, 1) = var_std
            var = var_all(imod, :, :, :)
            delete(var_avg)
            delete(var_std)
        end if

        ;; Optional output to NetCDF
        if (getenv("ESMValTool_write_netcdf").eq."True") then
            var@ncdf = procpath
            ncdf_outfile = ncdf_write(var, procpath)
        end if
        delete(var)
        delete(procpath)

    end do

    ;; ### CALCULATE DIFFERENCES AND APPLY STATISTICAL TESTS ##################

    ;; Calculate differences and apply t-test
    if (diag_script_info@plot_diff) then

        do imod = 0, dim_MOD - 1

            ;; Skip the reference model as result will be zero
            if (imod.eq.ref_ind) then
                continue
            end if

            if (irank.eq.4) then

                ;; Difference
                var_diff(imod, :, :) = var_all(imod, :, :, 0)  ; save metadata
                var_diff(imod, :, :) = var_all(imod, :, :, 0) - \
                    var_all(ref_ind, :, :, 0)

                ;; t-test (assuming different population variances)
                if (diag_script_info@t_test) then
                    x1 = var_all(imod, :, :, 0)
                    x2 = var_all(ref_ind, :, :, 0)
                    s1 = var_all(imod, :, :, 1) ^ 2
                    s2 = var_all(ref_ind, :, :, 1) ^ 2
                    n1 = get_end_year(imod) - get_start_year(imod) + 1
                    n2 = get_end_year(ref_ind) - get_start_year(ref_ind) + 1
                    prob = ttest(x1, s1, n1, x2, s2, n2, True, False)
                    var_sig(imod, :, :) = 1. - prob
                    var_diff@conf_level = diag_script_info@conf_level
                    delete(prob)
                end if

            end if

            if (irank.eq.3) then

                ;; Difference
                var_diff(imod, :) = var_all(imod, :, 0)  ; save metadata
                var_diff(imod, :) = var_all(imod, :, 0) - \
                    var_all(ref_ind, :, 0)

                ;; t-test (assuming different population variances)
                if (diag_script_info@t_test) then
                    x1 = var_all(imod, :, 0)
                    x2 = var_all(ref_ind, :, 0)
                    s1 = var_all(imod, :, 1) ^ 2
                    s2 = var_all(ref_ind, :, 1) ^ 2
                    n1 = get_end_year(imod) - get_start_year(imod) + 1
                    n2 = get_end_year(ref_ind) - get_start_year(ref_ind) + 1
                    prob = ttest(x1, s1, n1, x2, s2, n2, True, False)
                    var_sig(imod, :) = 1. - prob
                    delete(prob)
                end if

            end if

        end do
    end if

    ;; ### GRADING ############################################################

    if (diag_script_info@calc_grading) then ;; FIX-ME loop over metrics?

        ;; Annotation & file names
        region_name = ""  ; priority 3
        location = ""  ; priority 3
        if(isatt(diag_script_info, "aux_info")) then
            region_name = "_" + diag_script_info@aux_info  ; priority 2
            location = " - " + diag_script_info@aux_info  ; priority 2
        end if
        if(isatt(diag_script_info, "location")) then
            location = " - " + diag_script_info@location  ; priority 1
        end if
        if(isatt(diag_script_info, "region")) then
            region_name = "_" + region@name  ; priority 1
            location = " - " + diag_script_info@region
        end if
        altitude = ""  ; priority 2
        if(any(field_type0.eq.(/"T0M", "T0Ms"/))) then
            if(isatt(diag_script_info, "altitude")) then
                altitude = " - " + diag_script_info@altitude  ; priority 1
            end if
        else
            if(isatt(diag_script_info, "level")) then
                altitude = " - " + diag_script_info@level + " hPa"  ; priority 1
            end if
        end if

        ;; Set variable name
        var0_new = var0
        var0_new = var0_new + region_name
        if(isatt(diag_script_info, "level").and. \
            any(field_type0.eq.(/"T3M", "T1M"/))) then
            var0_new = var0_new + "-" + diag_script_info@level
        end if

        ;; Define grading file directory
        work_dir = getenv("ESMValTool_work_dir")
        yml_name = getenv("ESMValTool_yml_name")
        suffix = get_file_suffix(yml_name, 0)
        name = str_sub_str(yml_name, suffix, "")
        ncdf_dir = get_output_dir(work_dir, diag_script_base) + \ ;; FIX-ME ?
            name + "_" + diag_script_info@metric + "_" + var0_new + ".nc"

        ;; Start calculation
;;DELETE        if (isfilepresent_esmval(ncdf_dir)) then
;;DELETE
;;DELETE            info_output("No recalculation!", verbosity, 1)
;;DELETE            info_output("Reading in file = " + ncdf_dir, verbosity, 1)
;;DELETE
;;DELETE        else

;; TO BE REMOVED ONCE AVAILABLE IN BACKEND+++

        ;; Calculate multi-model mean/median
        if (diag_script_info@metric.eq."taylor") then
            id_remove = ref_ind
        else
            id_remove = ref_inds
        end if
        models_only = remove_index(var_all, id_remove)
        super_array = var_all
        if (isatt(diag_script_info, "MultiModelMedian")) then
            if (diag_script_info@MultiModelMean) then
                if (irank.eq.3) then
                    median = models_only(0:0, :, :)
                    median(0, :, :) = dim_median_n(models_only, 0)
                end if
                if (irank.eq.4) then
                    median = models_only(0:0, :, :, :)
                    median(0, :, :, :) = dim_median_n(models_only, 0)
                end if
                if (irank.eq.5) then
                    median = models_only(0:0, :, :, :, :)
                    median(0, :, :, :, :) = dim_median_n(models_only, 0)
                end if
                median&model = "multi-model-median"
                super_array := array_append_record(median, super_array, 0)
                delete(median)
                id_remove = id_remove + 1  ; update ref model index
            end if
        end if

        if (isatt(diag_script_info, "MultiModelMean")) then ;; FIX-ME
            if (diag_script_info@MultiModelMean) then
                if (irank.eq.3) then
                    mean = models_only(0:0, :, :)
                    mean(0, :, :) = dim_avg_n(models_only, 0)
                end if
                if (irank.eq.4) then
                    mean = models_only(0:0, :, :, :)
                    mean(0, :, :, :) = dim_avg_n(models_only, 0)
                end if
                if (irank.eq.5) then
                    mean = models_only(0:0, :, :, :, :)
                    mean(0, :, :, :, :) = dim_avg_n(models_only, 0)
                end if
                mean&model = "multi-model-mean"
                super_array := array_append_record(mean, super_array, 0)
                delete(mean)
                id_remove = id_remove + 1  ; update ref model index
            end if
        end if
        delete(models_only)
;; TO BE REMOVED ONCE AVAILABLE IN BACKEND---
      
        ;; Create separate arrays for models and observations
        models_only = remove_index(super_array, id_remove)
        nmodels = dimsizes(models_only&model)
        if (l_altern .and. diag_script_info@metric.ne."taylor") then
            if (irank.eq.3) then
                obs_only = super_array(id_remove, :, :)
            end if
            if (irank.eq.4) then
                obs_only = super_array(id_remove, :, :, :)
            end if
            if (irank.eq.5) then
                obs_only = super_array(id_remove, :, :, :, :)
            end if
        else  ; avoid dimension reduction
            if (irank.eq.3) then
                obs_only = super_array(id_remove:id_remove, :, :)
            end if
            if (irank.eq.4) then
                obs_only = super_array(id_remove:id_remove, :, :, :)
            end if
            if (irank.eq.5) then
                obs_only = super_array(id_remove:id_remove, :, :, :, :)
            end if
        end if
        nobs = dimsizes(obs_only&model)
        delete(super_array)

        ;; Define result variable
        ;; A dummy coordinate "diagnostics" (size 1) is added to facilitate
        ;; appending different variable in the _collect script
        if (diag_script_info@metric.eq."taylor") then
            grading = new((/dim_VAR, nmodels, 2/), float)
            grading!2 = "statistic"
            grading&statistic = (/"stddev_ratio", "correlation"/)
        else
            grading = new((/dim_VAR, nmodels, nobs/), float)
            grading!2 = "reference"
            grading&reference = obs_only&model
        end if
        grading!0 = "diagnostics"
        grading!1 = "models"
        grading&diagnostics = var0_new
        grading&models = models_only&model

        grading@ncdf_dir = ncdf_dir

        ;; Loop over models, calculate metric
        do imod = 0, nmodels - 1
            do iobs = 0, nobs - 1
                if (irank.eq.3) then
                    model = models_only(imod, :, 0)
                    obs = obs_only(iobs, :, 0)
                end if
                if (irank.eq.4) then
                    model = models_only(imod, :, :, 0)
                    obs = obs_only(iobs, :, :, 0)
                end if
                if (irank.eq.5) then
                    model = models_only(imod, :, :, :, 0)
                    obs = obs_only(iobs, :, :, :, 0)
                end if
                if (diag_script_info@metric.eq."taylor") then
                    grading(0, imod, 0) = \
                        calculate_metric(model, obs, "stddev_ratio")
                    grading(0, imod, 1) = \
                        calculate_metric(model, obs, "correlation")
                else
                    grading(0, imod, iobs) = \
                        calculate_metric(model, obs, diag_script_info@metric)
                end if
            end do
        end do
        
        ;; Apply normalization
        if (diag_script_info@metric.ne."taylor") then
            do iobs = 0, nobs - 1
                grading(:, :, iobs) = \
                    normalize_metric(grading(:, :, iobs), \
                                     diag_script_info@normalization)
            end do
        end if

        ;; Reduce dimensionality if no alternative model
        if (.not.l_altern .and. diag_script_info@metric.ne."taylor") then
            grading := grading(:,:,0)
            delete(grading@reference)
        end if

        ;; Attach attributes to the results
        grading@title = diag_script_info@metric + "metric"
        grading@long_name = \
            "Grading table of metric " + diag_script_info@metric
        grading@metric = diag_script_info@metric
        grading@diag_script = (/diag_script/)
        grading@var = "grade" ;; FIX-ME ?
        grading@region = location
        grading@num_climofiles = dimsizes(models@name)  ; for tagging
        do imod = 0, dimsizes(models@name) - 1
            num_climo = "climofile_"+imod
            climofile = interface_get_inpaths(imod) + "/" + \
                interface_get_infile(variables(0), field_types(0), imod)
                grading@$num_climo$ = climofile
        end do 

        ;; Write NetCDF output
        ncdf_outfile = ncdf_write(grading, grading@ncdf_dir)

        ;; Write results of temporary grading list
        yml_name = getenv("ESMValTool_yml_name")
        suffix = get_file_suffix(yml_name, 0)
        name = str_sub_str(yml_name, suffix, "")
        temp_dir = get_output_dir(getenv("ESMValTool_work_dir"), "temp") + \
            name
        if (diag_script_info@metric.eq."taylor") then
            temp_dir = temp_dir + "_taylor.nc"
        else
            temp_dir = temp_dir + "_grading.nc"
        end if

        if (isfilepresent_esmval(temp_dir)) then
            temp_file = addfile(temp_dir, "r")
            temp_list = temp_file->temp_list
            temp_list_string = tostring(temp_list)
            temp_list_new_string = array_append_record(temp_list_string, \
                                                       ncdf_dir, 0)
            temp_list_new = tochar(temp_list_new_string)
            system("rm -f " + temp_dir)
        else
            ncdf_char = tochar(ncdf_dir)
            temp_list_new = new((/1, dimsizes(ncdf_char)/), character)
            temp_list_new(0, :) = ncdf_char
        end if

        ;; Create new file and add list
        temp = addfile(temp_dir, "c")
        temp->temp_list = temp_list_new

    end if

    ;; ### PLOTTING ###########################################################

    if (diag_script_info@draw_plots) then

        ;; Convert units for plotting (if required)
        if (isatt(variable_info, "plot_units")) then
            var_all = convert_units(var_all, variable_info@plot_units)
            if (diag_script_info@plot_diff) then
                var_diff = convert_units(var_diff, variable_info@plot_units)
            end if
        end if

        ;; Plot output directory
        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, diag_script_base)

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
            if (isatt(diag_script_info, "level")) then ;; FIX-ME
                if (diag_script_info@level.ne."all") then
                    altitude = " - " + diag_script_info@level + " hPa"
                end if
            end if
        end if

        ;; CASE: cycle ########################################################
        if (ptype.eq."cycle") then
            outfile = output_dir + diag_script_base + "_" + var0 + "_" + \
                ptype + "_" + time_avg + "_" + region_name
            wks = gsn_open_wks(file_type, outfile)
            wks@legendfile = outfile + "_legend"
            var_all@res_tiMainString = varstring
            var_all@res_tiMainString = var_all@res_tiMainString + \
                location + altitude
            plot = cycle_plot(wks, var_all, var0)
            draw(plot)
            frame(wks)
            info_output(" gv " + outfile + "." + file_type, verbosity, 1)
        end if

        ;; CASE: zonal ########################################################
        if (ptype.eq."zonal") then
            do imod = 0, dim_MOD - 1
                outfile = output_dir + diag_script_base + "_" + \
                    modelnames(imod) + "_" + var0 + "_" + ptype + "_" + \
                    time_avg + "_" + level + region_name
                var = var_all(imod, :, :, 0)
                var@res_tiMainString = modelnames(imod)
                var@res_cnLevelSelectionMode = "ExplicitLevels"
                var@res_cnLevels = variable_info@zonal_ref_Levels
                if (isatt(variable_info, "zonal_ymin")) then
                    var@res_trYMinF = variable_info@zonal_ymin
                else
                    var@res_trYMinF = 5.
                end if
                var@res_gsnLeftString = var@long_name + " [" + \
                    format_units(var@units) + "]"
                wks = gsn_open_wks(file_type, outfile)
                if (isatt(variable_info, "zonal_ref_cmap")) then
                    gsn_define_colormap(wks, variable_info@zonal_ref_cmap)
                else
                    gsn_define_colormap(wks, "amwg_blueyellowred")
                end if
                plot = zonalmean_profile(wks, var, var0)
                draw(plot)
                frame(wks)
                info_output(" gv " + outfile + "." + file_type, \
                            verbosity, 1)
                delete(var@res_cnLevels)
                delete(wks)
                delete(plot)

                ;; Difference plot to the reference
                if (.not.diag_script_info@plot_diff .or. imod.eq.ref_ind) then
                    continue
                end if
                var = var_diff(imod, :, :)
                outfile = output_dir + diag_script_base + "_" + \
                    modelnames(imod) + "-" + modelnames(ref_ind) + "_" + \
                    var0 + "_" + ptype + "_" + time_avg + "_" + \
                    level + region_name
                var@res_tiMainString = modelnames(imod) + " - " + \
                modelnames(ref_ind)
                if (isatt(diag_script_info, "conf_level")) then
                    var@res_gsnLeftString = var@res_gsnLeftString + " - " + \
                    sprinti("%2i", toint(100 * diag_script_info@conf_level)) \
                            + "% c.l."
                end if
                var@res_cnLevelSelectionMode = "ExplicitLevels"
                var@res_cnLevels = variable_info@zonal_diff_Levels
                var@res_cnMissingValFillColor = "gray70"
                var@res_gsnSpreadColorEnd = -2
                wks = gsn_open_wks(file_type, outfile)
                gsn_merge_colormaps(wks, "temp_19lev", "gray70")
                plot = zonalmean_profile(wks, var, var0)

                ;; Stippling ;; FIX-ME False option not considered (see below)
                var = var_sig(imod, :, :)
                var@res_cnMissingValFillColor = -1
                var@res_gsnSpreadColorEnd = -2
                var@res_gsnDraw = False
                var@res_gsnFrame = False
                var@res_cnLevelSelectionMode = "ExplicitLevels"
                var@res_cnFillColors = (/"transparent", "black"/)
                var@res_cnLevels := diag_script_info@conf_level
                var@res_cnInfoLabelOn = False
                var@res_cnLinesOn = False
                var@res_cnLineLabelsOn = False
                var@res_lbLabelBarOn = False
                var@res_cnFillPattern = 17
                delete(var@long_name)
                plot2 = zonalmean_profile(wks, var, var0)
                overlay(plot, plot2)
                draw(plot)
                frame(wks)
                info_output(" gv " + outfile + "." + file_type, verbosity, 1)
                delete(var@res_cnLevels)
                delete(wks)
                delete(plot)
                delete(plot2)
                delete(var)
            end do
        end if

        ;; CASE: latlon #######################################################
        if (ptype.eq."latlon") then
            do imod = 0, dim_MOD - 1
                outfile = output_dir + diag_script_base + "_" + \
                    modelnames(imod) + "_" + var0 + "_" + ptype + "_" + \
                    time_avg + "_" + level + region_name
                var = var_all(imod, :, :, 0)
                var@projection = diag_script_info@projection
                var@res_cnLinesOn = False
                var@res_cnLevelSelectionMode = "ExplicitLevels"
                var@res_tiMainString = modelnames(imod)
                var@res_cnLevels = variable_info@map_ref_Levels
                var@res_gsnLeftString = var@long_name + " [" + \
                    format_units(var@units) + "]"

                ;; Append global field average
                if (variable_info@show_global_avg) then
                    gavg = area_operations(var, -90., 90., 0., 360., \
                                           "average", True)
                    var@res_gsnRightString = \
                        sprintf("%5.2f", gavg) + " " + format_units(var@units)
                end if
                wks = gsn_open_wks(file_type, outfile)
                if (isatt(variable_info, "map_ref_cmap")) then
                    gsn_define_colormap(wks, variable_info@map_ref_cmap)
                else
                    gsn_define_colormap(wks, "amwg_blueyellowred")
                end if
                plot = contour_map_ce(wks, var, var0)
                draw(plot)
                frame(wks)
                info_output(" gv " + outfile + "." + file_type, \
                            verbosity, 1)
                delete(var@res_cnLevels)
                delete(wks)
                delete(plot)

                ;; Difference plot to the reference
                if (.not.diag_script_info@plot_diff .or. imod.eq.ref_ind) then
                    continue
                end if
                var = var_diff(imod, :, :)
                outfile = output_dir + diag_script_base + "_" + \
                    modelnames(imod) + "-" + modelnames(ref_ind) + "_" + \
                    var0 + "_" + ptype + "_" + time_avg + "_" + level + \
                    region_name
                var@res_tiMainString = modelnames(imod) + " - " + \
                modelnames(ref_ind)
                if (isatt(diag_script_info, "conf_level")) then
                    var@res_gsnLeftString = var@res_gsnLeftString + " - " + \
                    sprinti("%2i", toint(100 * diag_script_info@conf_level)) \
                    + "% c.l."
                end if
                var@res_cnLevelSelectionMode = "ExplicitLevels"
                var@res_cnLevels = variable_info@map_diff_Levels
                var@res_cnMissingValFillColor = "gray70"
                var@res_gsnSpreadColorEnd = -2
                var@projection = diag_script_info@projection

                ;; Append global field averag
                if (variable_info@show_global_avg) then
                    gavg = area_operations(var, -90., 90., 0., 360., \
                                           "average", True)
                    var@res_gsnRightString = \
                        sprintf("%5.2f", gavg) + " " + format_units(var@units)
                end if
                wks = gsn_open_wks(file_type, outfile)
                gsn_merge_colormaps(wks, "temp_19lev", "gray70")
                plot = contour_map_ce(wks, var, var0)

                ;; Stippling/masking
                res2 = True
                res2@cnLevelSelectionMode = "ExplicitLevels"
                res2@cnLevels = (/0.0, diag_script_info@conf_level, 1.0/)
                res2@cnInfoLabelOn = False
                res2@cnLinesOn = False
                res2@cnLineLabelsOn = False
                res2@gsnDraw = False
                res2@gsnFrame = False
                plot2 = gsn_csm_contour(wks, var_sig(imod, :, :), res2)
                opt = True
                if (diag_script_info@stippling) then
                    opt@gsnShadeFillType = "pattern"
                    opt@gsnShadeMid = 17
                else
                    opt@gsnShadeFillType = "color"
                    opt@gsnShadeLow = "gray70"
                end if
                plot2 = gsn_contour_shade(plot2, diag_script_info@conf_level, \
                                          999., opt)
                overlay(plot, plot2)
                draw(plot)
                frame(wks)
                info_output(" gv " + outfile + "." + file_type, verbosity, 1)
                delete(var@res_cnLevels)
                delete(wks)
                delete(plot)
                delete(plot2)
                delete(var)
            end do

        end if
    end if

    leave_msg(diag_script, "", 4)

end
