;;#############################################################################
;; MAIN SCRIPT FOR PERFORMANCE METRICS
;; Author: Franziska Frank (DLR, Germany)
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
;;    valid_fraction: exclude gridbox with less than this fraction of non-
;;                    missing values in the time coordinate
;;    draw_plots: switch on/off plots drawing
;;    plot_diff: switch on/off difference plots
;;    dev: plot standard deviation
;;    legend_outside: save legend in a separate file
;;    styleset: as in lib/ncl/style.ncl functions
;;    t_test: calculate t-test in difference plots
;;    conf_level: confidence level for the t-test
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    grid: regridding option
;;    region: selected region
;;    level: pressure level in hPa (for T3M fields only)
;;    range_option: time range selection option (default 0)
;;    projection: map projection for lat-lon plots (default = cylindrical equidistant)
;;    stippling: True = use stippling to mark statistically significant differences
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
;;    The reference model must be specified in the main namelsit in the 
;;    <variable> tag as <variable ref_model="XX">. An alternative model
;;    can be also given, separated by a comme <variable ref_model="XX,"YY">.
;;    Consistent regridding and a uniform missing value mask are applied to all
;;    models before the diagnostics are calculated. Regarding the masking, the
;;    following shall be considered:
;;        - timesteps containg more than 100 * (1 - valid_fraction) % of
;;          missing values are excluded from this procedure, to avoid
;;          propagating missing values due to incomplete data coverage (e.g.,
;;          vmrh2o in HALOE or hus in AIRS). The threshold can be set via
;;          namelist (0.95 seem to work fine for most cases)
;;        - for T3M fields, the operation is performed on a level basis, to
;;          avoid propagating missing levels (from models/observations with
;;          lower vertical extent, e.g., hus in AIRS) to other models.
;;        - in any case, masks resulting only in missing values are excluded
;;          and a warning is issued.
;;
;; Modification history
;;    20170512_A_rett_ri: added tags
;;    20161220_A_laue_ax: added option to set map projection for lat-lon plots
;;                        (diag_script_info@projection)
;;                        added option to choose how to plot t-test results:
;;                        stippling or masking out in gray (lat-lon plots only)
;;    20161019_A_laue_ax: changed plotting of t-test results:
;;                        now stippling significant grid cells
;;                        (old version: masking out non-significant values in gray)
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

load "./interface_data/ncl.interface"
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
load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"
load "./diag_scripts/lib/ncl/meta_data.ncl"
load "./diag_scripts/lib/ncl/meta_data.ncl"

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
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin

    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"plot_type", "time_avg", "valid_fraction", "draw_plots", \
                 "plot_diff", "plot_stddev", "legend_outside", "styleset", \
                 "t_test", "conf_level"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Set plot type
    ptype = diag_script_info@plot_type
    if (all(ptype.ne.(/"cycle", "zonal", "latlon", "cycle_latlon"/))) then
        error_msg("f", diag_script, "", "plot_type " + ptype + " is not a " + \
                  "valid plot_type")
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

    ;; stippling / or masking out of significant/non-significant differences
    if (isatt(diag_script_info, "stippling")) then
        stippling = diag_script_info@stippling
    else
        stippling = False  ; default = mask out non-significant differences in gray
    end if

    ;; Selected time average
    if (isatt(diag_script_info, "time_avg")) then
        time_avg = diag_script_info@time_avg
    end if

    ;; Selected level
    if (isatt(diag_script_info, "level")) then
        level = diag_script_info@level
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

    ;; Check for confidence level
    if (diag_script_info@t_test .and. \
        .not.isatt(diag_script_info, "conf_level")) then
        error_msg("f", diag_script, "", "a confidence level must be " + \
                  "specified for the t-test to be executed")
    end if

    ;; Check for regridding method
    if (.not.(isatt(diag_script_info, "grid"))) then
        error_msg("w", diag_script, "", "regridding method not defined. " + \
                  "Using coarsest")
        selected_grid = "coarsest"
    else
        selected_grid = diag_script_info@grid
    end if

    ;; Check for reference model definition
    if (var_attr_ref(0).eq."None") then
        error_msg("f", diag_script, "", "no reference model is specified")
    end if

    ;; Determine reference (and alternative) model
    ref_ind = get_ref_model_idx(models, var_attr_ref(0))
    if (dimsizes(ref_ind).gt.1) then
        alt_ind = ref_ind(1)
        ref_ind := ref_ind(0)
    end if

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Check for irregular grid
    l_irreg = new(dim_MOD, logical)
    l_irreg = False
    do imod = 0, dim_MOD - 1
        data = read_data(imod, var0, field_type0)
        cnt = interface_get_dictkeys(var0, imod) + "_dimNames"
        dims = data@$cnt$
        ndims = dimsizes(dims)
        if (ndims.lt.2) then
            continue
        end if
        if (ndims.eq.2 .and. dims(0).eq."time") then
            continue
        end if
        if (dims(ndims - 1).ne."lon" .and. dims(ndims - 2).ne."lat") then
            l_irreg(imod) = True
            continue
        end if
    end do
    delete([/data, cnt, dims, ndims/])

    ;; Define target grid for regridding
    if (any(l_irreg) .and. selected_grid.ne."1x1") then
        error_msg("f", diag_script, "", "models are on irregular grid, " + \
                  "diag_script_info@grid = '1x1' must be used")
    end if

    models_ind = ispan(0, dim_MOD - 1, 1)
    info_output("Setting target grid:", verbosity, 2)
    if (selected_grid.eq."ref") then
        grid = find_destination_grid(ref_ind, \
                                     var0, field_type0, "coarsest")
    else if (selected_grid.eq."1x1") then
        grid     = new((/180, 360/), float)
        grid!0   = "lat"
        grid!1   = "lon"
        grid&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
        grid&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
        grid&lat@standard_name = "latitude"
        grid&lon@standard_name = "longitude"
        grid&lat@long_name     = "latitude"
        grid&lon@long_name     = "longitude"
        grid&lat@units         = "degrees_north"
        grid&lon@units         = "degrees_east"
        grid&lat@axis          = "Y"
        grid&lon@axis          = "X"
    else
        grid = find_destination_grid(models_ind, \
                                     var0, field_type0, selected_grid)
    end if
    end if
    delete(models_ind)

    ;; Regrid
    opt = True
    opt@plot_type = "regridded"
    opt@time_avg = selected_grid
    do imod = 0, dim_MOD - 1

        info_output("Regridding " + modelnames(imod), verbosity, 2)
        regr_path = clim_diag_file_path(var0, field_type0, imod, \
                                        diag_script, opt)

        ;; Do not process data if already processed
        if (isfilepresent(regr_path) .and. \
            getenv("ESMValTool_force_calc").eq."False") then

            info_output("No recalculation!", verbosity, 1)
            info_output("Reading in file = " + regr_path, verbosity, 1)
            pre_var = ncdf_read(regr_path, var0)
            dnames = getVarDimNames(pre_var)

        else
            ;; Extract model corresponding data
            data = read_data(imod, var0, field_type0)
            var = extract_data_extended(imod, var0, data, -1, 0, 0)
            dnames = getVarDimNames(var)

            ;; Regular or irregular grid?
            if (l_irreg(imod)) then  ; irregular grid

                ;; Regrid using ESMF_regrid

                ;; ESMF_regrid cannot handle integer coordinates
                do ii = 0, dimsizes(dnames) - 1
                    if (dnames(ii).eq."time") then
                        continue
                    end if
                    if (typeof(var&$dnames(ii)$).eq."integer") then
                        var&$dnames(ii)$ := tofloat(var&$dnames(ii)$)
                    end if
                end do

                ;; Pass information of irregular grid as attributes of var
                ;;  - lat/lon grid information comes from another file
                ;; this is NOT CMOR standard! Reformat below.
                var@lat2d := extract_area(imod, data, "lat", var0)
                var@lon2d := extract_area(imod, data, "lon", var0)

                ;; Naming of intermediate regrid files: will be re-used if
                ;; already present - ensure unambiguous naming!
                ;; Because of the mask of missing values, even different
                ;; time-steps / model levels need their own wgt file
                stem   = interface_get_inpaths(imod) + \
                    "/ESMF_regrid_" + models@project(imod) + "_" + \
                    models@name(imod) + "_" + models@experiment(imod) + "_" + \
                    models@ensemble(imod)
                ending = var0 + "_" + models@start_year(imod) + \
                    "-" + models@end_year(imod) + ".nc"

                ;; Set options
                Opt = True
                Opt@SrcMask2D  = where(ismissing(var(0, :, :)), 0, 1)
                Opt@DstGridLat = grid&lat
                Opt@DstGridLon = grid&lon
                ;Opt@DstGridType    = "rectilinear"  ; problem with that Opt
                Opt@ForceOverwrite = True
                Opt@PrintTimings   = True
                Opt@Debug          = True
                Opt@InterpMethod   = "patch"  ; "bilinear", "patch", "conserve"

                ;; Save regrid files for later use: careful with naming!
                Opt@WgtFileName = stem + "_Wgt_" + var0 + ".nc"
                Opt@SrcFileName = stem + "_Src_" + ending
                Opt@DstFileName = stem + "_Dst_" + ending

                ;; Reuse regrid files if present
                if (isfilepresent(Opt@SrcFileName)) then
                    Opt@SkipSrcGen = True
                end if
                if (isfilepresent(Opt@DstFileName)) then
                    Opt@SkipDstGen = True
                end if
                if (isfilepresent(Opt@WgtFileName)) then
                    Opt@SkipWgtGen = True
                    var_hreg = ESMF_regrid_with_weights(var, \
                                                        Opt@WgtFileName, Opt)
                else
                    var_hreg = ESMF_regrid(var, Opt)
                end if
                delete(Opt)

            else  ; regular grid

                ;; Determine dimensions to be regridded
                xregrid = False  ; lon
                yregrid = False  ; lat
                if (any(dnames.eq."lon")) then
                    xregrid = True
                    if (dimsizes(var&lon).eq.dimsizes(grid&lon)) then
                        if (max(abs(var&lon - grid&lon)).lt.1.d-5) then
                            xregrid = False
                        end if
                    end if
                end if
                if (any(dnames.eq."lat")) then
                    yregrid = True
                    if (dimsizes(var&lat).eq.dimsizes(grid&lat)) then
                        if (max(abs(var&lat - grid&lat)).lt.1.d-5) then
                            yregrid = False
                        end if
                    end if
                end if

                ;; Regrid in the horizontal (if necessary)
                if (xregrid.and.yregrid) then
                    var_hreg = linint2_Wrap(var&lon, var&lat, var, True, \
                                            grid&lon, grid&lat, 0)
                else if (xregrid.and..not.yregrid) then
                    dx = ind(dnames.eq."lon")
                    var_hreg = linint1_n_Wrap(var&lon, var, True, grid&lon, \
                                              0, dx)
                    info_output("Latitude regridding not required", \
                                verbosity, 4)
                else if (.not.xregrid.and.yregrid) then
                    dx = ind(dnames.eq."lat")
                    var_hreg = linint1_n_Wrap(var&lat, var, False, \
                                              grid&lat, 0, dx)
                    info_output("Longitude regridding not required", \
                                verbosity, 4)
                else if (.not.xregrid .and. .not.yregrid) then
                    var_hreg = var
                    info_output("Horizontal regridding not required", \
                                verbosity, 4)
                end if
                end if
                end if
                end if

            end if  ; regular or irregular grid?

            ;; Regrid in the vertical (if necessary)
            zregrid = False  ; plev
            if (any(dnames.eq."plev")) then
                zregrid = True
                if (dimsizes(var&plev).eq.dimsizes(grid&plev)) then
                    if (max(abs(var&plev - grid&plev)).lt.1.d-5) then
                        zregrid = False
                    end if
                end if
            end if
            if (zregrid) then
                var_vreg = int2p_n_Wrap(var_hreg&plev, var_hreg, \
                                        grid&plev, 1, 1)
            else
                var_vreg = var_hreg
                info_output("Vertical regridding not required", verbosity, 4)
            end if
            delete(var_hreg)
            pre_var = var_vreg
            delete(var_vreg)

            ;; Output to NetCDF
            pre_var@var = var0
            pre_var@diag_script = (/diag_script/)
            pre_var@ncdf = regr_path
            ncdf_outfile = ncdf_write(pre_var, regr_path)
            delete([/data, var/])

        end if  ; data already processed?

        ;; Create a missing value mask for each model
        ;; The dim_avg_wgt function with option limit is used on the time
        ;; coordinate: in this way each grid box containing less than
        ;; valid_fraction of missing values along the time series is
        ;; set to missing.
        if (imod.eq.0) then
            global_mask = dim_avg_n_Wrap(pre_var, 0)  ; to get rid of time dim
            global_mask = 0.
        end if

        ww = pre_var&time
        ww = 1.
        limit = toint(diag_script_info@valid_fraction * \
                      dimsizes(pre_var&time))

        if (any(dnames.eq."plev")) then  ; process levels individually

            do kk = 0, dimsizes(pre_var&plev) - 1

                if (dimsizes(dimsizes(global_mask)).eq.3) then
                    lsel = pre_var(:, kk, :, :)
                else if (dimsizes(dimsizes(global_mask)).eq.2) then
                    lsel = pre_var(:, kk, :)
                else if (dimsizes(dimsizes(global_mask)).eq.1) then
                    lsel = pre_var(:, kk)
                else
                    error_msg("f", diag_script, "", "case not covered")
                end if
                end if
                end if

                if (.not.all(ismissing(lsel))) then
                    ;; dim_avg_wgt_n_Wrap produces an invalid result for rank=1
                    if (dimsizes(dimsizes(lsel)).eq.1) then
                        tmp = dim_avg_wgt_Wrap(lsel, ww, limit)
                    else
                        tmp = dim_avg_wgt_n_Wrap(lsel, ww, limit, 0)
                    end if
                    tmp = where(ismissing(tmp), tmp, 0.)
                    if (all(ismissing(tmp))) then
                        error_msg("w", diag_script, "", modelnames(imod) + \
                                  " (level " + kk + ") contains only " + \
                                  "missing values, excluding from global " + \
                                  "masking.")
                    else
                        if (dimsizes(dimsizes(global_mask)).eq.3) then
                            global_mask(kk, :, :) = global_mask(kk, :, :) + tmp
                        else if (dimsizes(dimsizes(global_mask)).eq.2) then
                            global_mask(kk, :) = global_mask(kk, :) + tmp
                        else if (dimsizes(dimsizes(global_mask)).eq.1) then
                            global_mask(kk) = global_mask(kk) + tmp
                        end if
                        end if
                        end if
                    end if
                    delete(tmp)
                end if
                delete(lsel)

            end do  ; levels

        else  ; no plev dimension
            ;; FIXME: Consider distinguishing other dimensions as well?

            ;; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1
            if (dimsizes(dimsizes(pre_var)).eq.1) then
                tmp = dim_avg_wgt_Wrap(pre_var, ww, limit)
            else
                tmp = dim_avg_wgt_n_Wrap(pre_var, ww, limit, 0)
            end if
            tmp = where(ismissing(tmp), tmp, 0.)
            if (all(ismissing(tmp))) then
                error_msg("w", diag_script, "", modelnames(imod) + \
                          "contains only missing values, excluding " + \
                          "from global masking.")
            else
                global_mask = global_mask + tmp
            end if
            delete(tmp)
        end if
        delete(ww)
        delete(pre_var)

    end do  ; imod

    ;; Define global arrays for collecting model data
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

    ; if needed, create array for statistically significant grid cells
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

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + modelnames(imod), verbosity, 2)
        proc_path = clim_diag_file_path(var0, field_type0, imod, \
                                        diag_script, diag_script_info)

        ;; Do not process data if already processed
        if (isfilepresent_esmval(proc_path) .and. \
            getenv("ESMValTool_force_calc").eq."False") then
            info_output("No recalculation!", verbosity, 1)
            info_output("Reading in file = " + proc_path, verbosity, 1)
            if (irank.eq.5) then
                var_all(imod, :, :, :, :) = ncdf_read(proc_path, var0)
            end if
            if (irank.eq.4) then
                var_all(imod, :, :, :) = ncdf_read(proc_path, var0)
            end if
            if (irank.eq.3) then
                var_all(imod, :, :) = ncdf_read(proc_path, var0)
            end if
            continue
        end if

        ;; Apply global mask to conform missing values
        regr_path = clim_diag_file_path(var0, field_type0, imod, \
                                        diag_script, opt)
        pre_var = ncdf_read(regr_path, var0)
        if (sum(dimsizes(global_mask)).eq.1) then  ; scalar
            cdims = -1  ; special treatment, see NCL documentation of conform
        else
            cdims = ispan(1, dimsizes(dimsizes(pre_var)) - 1, 1)
        end if
        pre_var = pre_var + conform(pre_var, global_mask, cdims)
        delete(cdims)

        ;; Determine start/end year
        start_year = get_start_year(imod)
        end_year = get_end_year(imod)

        ;; Process data - Cycle case
        if (ptype.eq."cycle") then

            ;; Extract level
            ;; (level in hPa, pre_var&plev in Pa according to CMOR)
            if (any(dnames.eq."plev")) then

                if ((min(pre_var&plev).ge.100 * diag_script_info@level).or. \
                    (max(pre_var&plev).le.100 * diag_script_info@level)) then
                    info_output("Specified level outside range, " + \
                                "setting to missing values", verbosity, 2)
                    var_lev = dim_avg_n_Wrap(pre_var, 1)  ; get rid of plev
                    var_lev = var_lev@_FillValue
                else
                    tmp = int2p_n_Wrap(pre_var&plev, pre_var, \
                                       100. * diag_script_info@level, 1, 1)
                    if (dimsizes(dimsizes(tmp)).eq.4) then
                        var_lev = tmp(:, 0, :, :)
                    else if (dimsizes(dimsizes(tmp)).eq.3) then
                        var_lev = tmp(:, 0, :)
                    else if (dimsizes(dimsizes(tmp)).eq.2) then
                        var_lev = tmp(:, 0)
                    end if
                    end if
                    end if
                    delete(tmp)
                end if

            else  ; no plev dimension

                var_lev = pre_var

            end if
            delete(pre_var)

            ;; Extract region and average over latitude and longitude
            if (any(field_type0.eq.(/"T0M", "T0Ms", "T1M"/))) then
                var_reg = var_lev
            else if (any(dnames.eq."lat") .and. any(dnames.eq."lon")) then
                var_reg = area_operations(var_lev, region(0), region(1), \
                                          region(2), region(3), "average", \
                                          True)
            else
                error_msg("f", diag_script, "", "dimensionality not " + \
                          "implemented")  ; FIXME
            end if
            end if
            delete(var_lev)

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

        ;; Process data - Cycle lat lon case
        if (ptype.eq."cycle_latlon") then

            ;; Extract level
            ;; (level in hPa, pre_var&plev in Pa according to CMOR)
            if (field_type0.eq."T3M") then
                if ((min(pre_var&plev).ge.100 * diag_script_info@level) .or. \
                    (max(pre_var&plev).le.100 * diag_script_info@level)) then
                    info_output("Specified level outside range, " + \
                                "setting to missing values", verbosity, 2)
                    var_lev = pre_var(:, 0, :, :)
                    var_lev = var_lev@_FillValue
                else
                    tmp = int2p_n_Wrap(pre_var&plev, pre_var, \
                                       100. * diag_script_info@level, 1, 1)
                    var_lev = tmp(:, 0, :, :)
                    delete(tmp)
                end if
            else
                var_lev = pre_var
            end if
            delete(pre_var)

            ;; Extract region
            var = area_operations(var_lev, region(0), region(1), \
                                  region(2), region(3), "extract", True)
            delete(var_lev)

            ;; Calculate time average
            var_avg = time_operations(var, start_year, end_year, \
                                      "average", time_avg, True)

            ;; Standard deviation calculation for this rank is not
            ;; implemented yet in statistics.ncl (but is anyway not
            ;; required here)
            delete(var)

            ;; Assign to global array
            var_all(imod, :, :, :, 0) = var_avg
            delete(var_avg)
            var = var_all(imod, :, :, :, :)

        end if

        ;; Process data - Zonal mean case
        if (ptype.eq."zonal") then

            ;; Calculate zonal mean
            var_zon = dim_avg_Wrap(pre_var)
            delete(pre_var)

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

        ;; Process data - Lat-lon case
        if (ptype.eq."latlon") then

            ;; Calculate time average and standard deviation w.r.t.
            ;; interannual variability
            var_avg = time_operations(pre_var, start_year, end_year, \
                                      "average", time_avg, True)
            var_std = interannual_variability(pre_var, start_year, end_year, \
                                              time_avg)
            delete(pre_var)

            ;; Assign to global array
            var_all(imod, :, :, 0) = var_avg
            var_all(imod, :, :, 1) = var_std
            var = var_all(imod, :, :, :)
            delete(var_avg)
            delete(var_std)
        end if

        ;; Optional output to NetCDF
        if (getenv("ESMValTool_write_netcdf").eq."True") then
            var@ncdf = proc_path
            ncdf_outfile = ncdf_write(var, proc_path)
        end if
        delete(var)
        delete(proc_path)

    end do

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

    ;; Plotting
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

        ;; Annotation & file names
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

        ;; Plot - Cycle case
        if (ptype.eq."cycle") then
            outfile = output_dir + diag_script_base + "_" + var0 + "_" + \
                ptype + "_" + time_avg + "_" + level + region_name
            wks = gsn_open_wks(file_type, outfile)
            wks@legendfile = outfile + "_legend"
            var_all@res_tiMainString = varstring
            var_all@res_tiMainString = var_all@res_tiMainString + \
                location + altitude
            plot = cycle_plot(wks, var_all, var0)
            plot@outfile = outfile_plot + file_type
            draw(plot)
            frame(wks)
            info_output(" gv " + outfile + "." + file_type, verbosity, 1)

            ;; tagging
            ST_tag = (/"ST_mean", "ST_avg"/) 
            if (var_all@metric .eq. "RMSD") then
                ST_tag := (/"ST_mean", "ST_rmsd"/)
            end if 
            plot_fname = wks@fullname
            alltags = array_append_record(tags,\
                      (/region@DM_tag, "PT_times", ST_tag, variables, models@name/), 0)
            caption = "Cycle plot. Models: " + models@name + \
                      ", variable: " + variables(0) + \
                      ", start year: " + start_year + \
                      ", end year: " + end_year
            id = diag_script + "_example_id"   
            do imod = 0, dimsizes(models@name) -1
                climofiles(imod) = interface_get_inpaths(imod) + "/" + \
                    interface_get_infile(variables(0), field_types(0), imod)
            end do
            delete([/alltags, caption, id, climofiles/])
         end if

        ;; Plot - Zonal mean case
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
                plot@outfile = outfile_plot + file_type
                draw(plot)
                frame(wks)
                info_output(" gv " + outfile + "." + file_type, \
                            verbosity, 1)
                            
                ;; tagging
                ST_tag = (/"ST_mean", "ST_avg"/)
                if (var_all@metric .eq. "RMSD") then
                    ST_tag := (/"ST_mean", "ST_rmsd"/)
                end if                 
                plot_fname = wks@fullname
                alltags = array_append_record(tags,\
                          (/region@DM_tag, "PT_zonal", ST_tag, variables, models@name/), 0)
                caption = "Cycle plot. Models: " + models@name + \
                          ", variable: " + variables(0) + \
                          ", start year: " + start_year + \
                          ", end year: " + end_year
                id = diag_script + "_example_id"   
                climofiles(imod) = interface_get_inpaths(imod) + "/" + \
                        interface_get_infile(variables(0), field_types(0), imod)
                ESMValMD(plot@outfile, alltags, caption, id, variables, modesl@name(imod), climofiles)
                delete([/alltags, caption, id, climofiles/])                            

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
                        plot@outfile = outfile_plot + file_type
                
                        ; add stippling for statistically significant grid cells
                        var = var_sig(imod, :, :)
                        var@res_cnMissingValFillColor = -1 ; transparent
                        var@res_gsnSpreadColorEnd   = -2
                        var@res_gsnDraw = False       ; Do not draw plot
                        var@res_gsnFrame= False       ; Do not advance frome
                        var@res_cnLevelSelectionMode = "ExplicitLevels"
                        var@res_cnFillColors        = (/"transparent", "black"/)
                        delete(var@res_cnLevels)
                        var@res_cnLevels            = diag_script_info@conf_level
                        var@res_cnInfoLabelOn       = False       ; turn off info label
                        var@res_cnLinesOn           = False       ; do not draw contour lines
                        var@res_cnLineLabelsOn      = False       ; do not draw contour labels
                        var@res_lbLabelBarOn        = False
                        var@res_cnFillPattern       = 17

                        delete(var@long_name)
                        plot2 = zonalmean_profile(wks, var, var0)

                        overlay (plot, plot2)

                        draw(plot)
                        frame(wks)
                        
                        ;; tagging
                        ST_tag = (/"ST_mean", "ST_avg", "ST_diff"/)  
                        if (var_all@metric .eq. "RMSD") then
                            ST_tag := (/"ST_mean", "ST_diff", "ST_avg", "ST_rmsd"/)
                        end if 
                        plot_fname = wks@fullname
                        alltags = array_append_record(tags,\
                              (/region@DM_tag, "PT_zonal", ST_tag, variables, models@name/), 0)
                        caption = "Cycle plot. Models: " + models@name + \
                              ", variable: " + variables(0) + \
                              ", start year: " + start_year + \
                              ", end year: " + end_year
                        id = diag_script + "_example_id"   
                        climofiles(imod) = interface_get_inpaths(imod) + "/" + \
                               interface_get_infile(variables(0), field_types(0), imod)
                        ESMValMD(plot@outfile, alltags, caption, id, variables, modesl@name(imod), climofiles)
                        delete([/alltags, caption, id, climofiles/])  
                        
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

        ;; Plot - Lat-lon case
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
                plot@outfile = outfile_plot + file_type
                draw(plot)
                frame(wks)
                        
                ;; tagging
                ST_tag = (/"ST_mean", "ST_avg"/)
                if (var_all@metric .eq. "RMSD") then
                    ST_tag := (/"ST_mean", "ST_avg", "ST_rmsd"/)
                end if 
                plot_fname = wks@fullname
                alltags = array_append_record(tags,\
                    (/region@DM_tag, "PT_geo", ST_tag, variables, models@name/), 0)
                caption = "Cycle plot. Models: " + models@name + \
                    ", variable: " + variables(0) + \
                    ", start year: " + start_year + \
                    ", end year: " + end_year
                id = diag_script + "_example_id"   
                climofiles(imod) = interface_get_inpaths(imod) + "/" + \
                     interface_get_infile(variables(0), field_types(0), imod)
                ESMValMD(plot@outfile, alltags, caption, id, variables, modesl@name(imod), climofiles)
                delete([/alltags, caption, id, climofiles/])  
                
                info_output(" gv " + outfile + "." + file_type, \
                            verbosity, 1)                                           
                
                delete(var@res_cnLevels)
                delete(wks)
                delete(plot)

                ;; Difference plot to the reference
                if (imod.ne.ref_ind) then
                    if (diag_script_info@_diff) then
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
                        plot@outfile = outfile_plot + file_type
                        
                        ; add stippling/masking of statistically significant/not-significant grid cells

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
                        if (stippling) then
                            opt@gsnShadeFillType = "pattern"
                            opt@gsnShadeMid = 17
                        else
                            opt@gsnShadeFillType = "color"
                            opt@gsnShadeLow = "gray70"
                        end if

                        plot2 = gsn_contour_shade(plot2, diag_script_info@conf_level, 999., opt)
                        overlay (plot, plot2)

                        draw(plot)
                        frame(wks)
                        
                        ;; tagging
                        ST_tag = (/"ST_mean", "ST_avg", "ST_diff"/)
                        if (var_all@metric .eq. "RMSD") then
                            ST_tag := (/"ST_mean", "ST_avg", "ST_diff", "ST_rmsd"/)
                        end if 
                        plot_fname = wks@fullname
                        alltags = array_append_record(tags,\
                              (/region@DM_tag, "PT_geo", ST_tag, variables, models@name/), 0)
                        caption = "Cycle plot. Models: " + models@name + \
                              ", variable: " + variables(0) + \
                              ", start year: " + start_year + \
                              ", end year: " + end_year
                        id = diag_script + "_example_id"   
                        climofiles(imod) = interface_get_inpaths(imod) + "/" + \
                               interface_get_infile(variables(0), field_types(0), imod)
                        ESMValMD(plot@outfile, alltags, caption, id, variables, modesl@name(imod), climofiles)
                        delete([/alltags, caption, id, climofiles/])  
                        
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
