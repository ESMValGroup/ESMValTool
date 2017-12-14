;;#############################################################################
;; MAIN SCRIPT FOR PERFORMANCE METRICS
;; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Calculates and (optionally) plots annual/seasonal cycles, zonal means,
;;    lat-lon fields and time-lat-lon fields from input T3M or T2Ms data.
;;    The calculated fields can be also plotted as difference w.r.t. a given
;;    reference model. They are also used as input to calculate grading metrics
;;    (see perfmetrics_grading.ncl).
;;    Input data are regridded to a common grid and a consistent mask of
;;    missing values is applied to all models and observations.
;;    The user can select the target grid (the coarsest grid is used as
;;    default). For irregular grids, the 1x1 grid must be selected.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    plot_type: cycle (time), zonal (plev, lat), latlon (lat, lon) or
;;               cycle_latlon (time, lat, lon)
;;    time_avg: type of time average
;;    region: selected region
;;    draw_plots: switch on/off plots drawing
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    legend_outside: save legend in a separate file
;;    styleset: as in lib/ncl/style.ncl functions
;;    plot_stddev: plot standard deviation
;;    t_test: calculate t-test in difference plots
;;    conf_level: confidence level for the t-test
;;    range_option: time range selection option (default 0)
;;    projection: map projection for lat-lon plots (default = cylind. equid.)
;;    stippling: True = use stippling to mark stat. significant differences
;;               False = mask out non-significant differences in gray (default)
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
;;    20171124_A_righ_ma: completely revised to adapt it to the new backend
;;                        (level selection, regridding and masking now done
;;                        by the python preprocessor)
;;    20161220_A_laue_ax: added option to set map projection for lat-lon plots
;;                        (diag_script_info@projection)
;;                        added option to choose how to plot t-test results:
;;                        stippling or masking out in gray (lat-lon plots only)
;;    20161019_A_laue_ax: changed plotting of t-test results:
;;                        now stippling significant grid cells (old version:
;;                        masking out non-significant values in gray)
;;    20160628_A_righ_ma: moving ref_model specification from cfg- files to
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
load "./interface_scripts/interface.ncl"
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

    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"plot_type", "time_avg", "draw_plots"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Set plot type
    ptype = diag_script_info@plot_type
    if (all(ptype.ne.(/"cycle", "zonal", "latlon", "cycle_latlon"/))) then
        error_msg("f", diag_script, "", "plot_type " + ptype + " is not a " + \
                  "valid plot_type")
    end if

    ;; Check for plot-type specific settings
    if (ptype.eq."cycle") then
        exit_if_missing_atts(diag_script_info, \
                             (/"legend_outside", "styleset", "plot_stddev"/))
    end if
    if (ptype.eq."zonal" .or. ptype.eq."lat-lon") then
        exit_if_missing_atts(diag_script_info, (/"t_test", "conf_level"/))
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

    ;; Set stippling or masking of significant/non-significant differences
    if (isatt(diag_script_info, "stippling")) then
        stippling = diag_script_info@stippling
    else
        stippling = False  ; default = mask non-significant values in gray
    end if

    ;; Selected time average
    if (isatt(diag_script_info, "time_avg")) then
        time_avg = diag_script_info@time_avg
    end if

    ;; Selected region
    if (isatt(diag_script_info, "region")) then
        region = select_region(diag_script_info@region)
    end if

    ;; Determine time range (years) of all models
    if (.not.isatt(diag_script_info, "range_option")) then
        diag_script_info@range_option = 0  ; every model uses its own range
    end if

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Unique names for models
    modelnames = project_style(diag_script_info, "annots")

    ;; Check for reference model definition
    if (var_attr_ref(0).eq."None") then
        error_msg("f", diag_script, "", "no reference model is specified")
    end if

    ;; Determine reference (and alternative) model
    ref_ind = get_ref_model_idx(models, var_attr_ref)
    if (dimsizes(ref_ind).gt.1) then
        alt_ind = ref_ind(1)
        ref_ind := ref_ind(0)
    end if

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

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
        var_diff = var_all(:, :, 0)
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
        var_diff = var_all(:, :, :, :, 0)
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
        var_diff = var_all(:, :, :, 0)
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
        var_diff = var_all(:, :, :, 0)
    end if

    ; If needed, create array for statistically significant grid cells
    if (isatt(diag_script_info, "t_test")) then
        if (diag_script_info@t_test) then
            var_sig = var_diff
        end if
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

        ;; Set path for saving processed data
        fullpath = getenv("ESMValTool_work_dir") + "/" + \
            basename(diag_script) + "/"
        if (.not.isfilepresent_esmval(fullpath)) then
            system("mkdir -p " + fullpath)
        end if
        fname = basename(systemfunc("basename " + \
                         interface_get_fullpath(var0, field_type0, imod)))
        fname = fname + "_" + basename(diag_script)
        if (isatt(diag_script_info, "plot_type")) then
            fname = fname + "_" + diag_script_info@plot_type
        end if
        if (isatt(diag_script_info, "time_avg")) then
            fname = fname + "_" + diag_script_info@time_avg
        end if
        if (isatt(diag_script_info, "level")) then
            fname = fname + "_" + diag_script_info@level
        end if
        if (isatt(diag_script_info, "region")) then
            fname = fname + "_" + str_sub_str(diag_script_info@region, " ", "")
        end if
        procpath = fullpath + fname + ".nc"

        ;; Do not process data if already processed
        if (isfilepresent_esmval(procpath)) then
            info_output("No recalculation!", verbosity, 1)
            info_output("Reading in file = " + procpath, verbosity, 1)
            if (irank.eq.5) then
                var_all(imod, :, :, :, :) = ncdf_read(procpath, var0)
            end if
            if (irank.eq.4) then
                var_all(imod, :, :, :) = ncdf_read(procpath, var0)
            end if
            if (irank.eq.3) then
                var_all(imod, :, :) = ncdf_read(procpath, var0)
            end if
            continue
        end if

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
    if (.not.isatt(diag_script_info, "plot_diff")) then
        diag_script_info@plot_diff = False
    end if
    if (.not.isatt(diag_script_info, "t_test")) then
        diag_script_info@t_test = False
    end if
    if (diag_script_info@plot_diff .or. diag_script_info@t_test) then

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
;                    prob = 1. - prob
;                    var_diff(imod, :, :) = \
;                        where(prob.gt.diag_script_info@conf_level, \
;                              var_diff(imod, :, :), var_diff@_FillValue)
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
;                    prob = 1. - prob
;                    var_diff(imod, :) = \
;                        where(prob.gt.diag_script_info@conf_level, \
;                              var_diff(imod, :), var_diff@_FillValue)
;                    var_diff@conf_level = diag_script_info@conf_level
                    var_sig(imod, :) = 1. - prob
                    delete(prob)
                end if

            end if

        end do
    end if

    ;; ### PLOTTING ###########################################################

    if (diag_script_info@draw_plots) then

        ;; Convert units for plotting (if required)
        if (isatt(variable_info, "plot_units")) then
            var_all = convert_units(var_all, variable_info@plot_units)
            if (isdefined("var_diff")) then
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
            if (isatt(diag_script_info, "level")) then
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
                if (imod.ne.ref_ind) then
                    if (diag_script_info@plot_diff) then
                        var = var_diff(imod, :, :)
                        outfile = output_dir + diag_script_base + "_" + \
                            modelnames(imod) + "-" + modelnames(ref_ind) + \
                            "_" + var0 + "_" + ptype + "_" + \
                            time_avg + "_" + level + region_name
                        var@res_tiMainString = modelnames(imod) + " - " + \
                            modelnames(ref_ind)
                        if (isatt(diag_script_info, "conf_level")) then
                            var@res_gsnLeftString = var@res_gsnLeftString + \
                                " - " + \
                                sprinti("%2i", \
                                        toint(100 * \
                                              diag_script_info@conf_level)) + \
                                "% c.l."
                        end if
                        var@res_cnLevelSelectionMode = "ExplicitLevels"
                        var@res_cnLevels = variable_info@zonal_diff_Levels
                        var@res_cnMissingValFillColor = "gray70"
                        var@res_gsnSpreadColorEnd = -2
                        wks = gsn_open_wks(file_type, outfile)
                        gsn_merge_colormaps(wks, "temp_19lev", "gray70")
                        plot = zonalmean_profile(wks, var, var0)

                        ; Stippling
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
                        info_output(" gv " + outfile + "." + file_type, \
                                    verbosity, 1)
                        delete(var@res_cnLevels)
                        delete(wks)
                        delete(plot)
                        delete(plot2)
                    end if
                end if
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
                if (isatt(diag_script_info, "projection")) then
                    var@projection = diag_script_info@projection
                end if
                var@res_cnLinesOn = False
                var@res_cnLevelSelectionMode = "ExplicitLevels"
                var@res_tiMainString = modelnames(imod)
                var@res_cnLevels = variable_info@map_ref_Levels
                var@res_gsnLeftString = var@long_name + " [" + \
                    format_units(var@units) + "]"
                ;; Append global field average
                if (isatt(variable_info, "show_global_avg")) then
                    if (variable_info@show_global_avg) then
                        gavg = area_operations(var, -90., 90., 0., 360., \
                                               "average", True)
                        var@res_gsnRightString = \
                            sprintf("%5.2f", gavg) + " " + \
                            format_units(var@units)
                    end if
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
                if (imod.ne.ref_ind) then
                    if (diag_script_info@plot_diff) then
                        var = var_diff(imod, :, :)
                        outfile = output_dir + diag_script_base + "_" + \
                            modelnames(imod) + "-" + modelnames(ref_ind) + \
                            "_" + var0 + "_" + ptype + "_" + \
                            time_avg + "_" + level + region_name
                        var@res_tiMainString = modelnames(imod) + " - " + \
                            modelnames(ref_ind)
                        if (isatt(diag_script_info, "conf_level")) then
                            var@res_gsnLeftString = var@res_gsnLeftString + \
                                " - " + \
                                sprinti("%2i", \
                                        toint(100 * \
                                              diag_script_info@conf_level)) + \
                                "% c.l."
                        end if
                        var@res_cnLevelSelectionMode = "ExplicitLevels"
                        var@res_cnLevels = variable_info@map_diff_Levels
                        var@res_cnMissingValFillColor = "gray70"
                        var@res_gsnSpreadColorEnd = -2
                        if (isatt(diag_script_info, "projection")) then
                            var@projection = diag_script_info@projection
                        end if

                        ;; Append global field average
                        if (isatt(variable_info, "show_global_avg")) then
                            if (variable_info@show_global_avg) then
                                gavg = area_operations(var, -90., 90., \
                                                       0., 360., \
                                                       "average", True)
                                var@res_gsnRightString = \
                                    sprintf("%5.2f", gavg) + " " + \
                                    format_units(var@units)
                            end if
                        end if
                        wks = gsn_open_wks(file_type, outfile)
                        gsn_merge_colormaps(wks, "temp_19lev", "gray70")
                        plot = contour_map_ce(wks, var, var0)

                        ; Stippling/masking
                        res2 = True
                        res2@cnLevelSelectionMode = "ExplicitLevels"
                        res2@cnLevels = \
                            (/0.0, diag_script_info@conf_level, 1.0/)
                        res2@cnInfoLabelOn = False
                        res2@cnLinesOn = False
                        res2@cnLineLabelsOn = False
                        res2@gsnDraw = False
                        res2@gsnFrame = False
                        plot2 = gsn_csm_contour(wks, var_sig(imod, :, :), res2)
                        opt = True
                        if (stippling) then
                            opt@gsnShadeFillType = "pattern"
                            opt@gsnShadeMid = 17
                        else
                            opt@gsnShadeFillType = "color"
                            opt@gsnShadeLow = "gray70"
                        end if
                        plot2 = gsn_contour_shade(plot2, \
                                                  diag_script_info@conf_level,\
                                                  999., opt)
                        overlay(plot, plot2)
                        draw(plot)
                        frame(wks)
                        info_output(" gv " + outfile + "." + file_type, \
                                    verbosity, 1)
                        delete(var@res_cnLevels)
                        delete(wks)
                        delete(plot)
                        delete(plot2)
                    end if
                end if
                delete(var)
            end do

        end if
    end if

    leave_msg(diag_script, "", 4)

end
