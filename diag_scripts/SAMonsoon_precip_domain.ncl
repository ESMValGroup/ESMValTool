;;#############################################################################
;; SAMonsoon_precip_intensity.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    Global monsoon intensity and monsoon domain plots.
;;    Intensity is computed as pr_summmer - pr_winter where summer and winter
;;    are assumed opposite on the nothern and southern hemisphere.
;;    Monsoon domains are areas where both,
;;        * pr_summer - pr_winter > mm/day cutoff (default 2.0 mm/day)
;;        * pr_summer > 0.55 * pr_annual
;;
;; Required diag_script_info attributes (diagnostics specific)
;;                      abs_cmap: colormap for domain plots
;;                     diff_cmap: colormap for domain diff plots
;;   high_intensity_cutoff_mmday: cutoff for monsoon domain definition
;;               latrange_global: latitude range to include in the diagnostic
;;                 summer_season: summer season, defaults to MJJAS
;;                 winter_season: winter season, defaults to NDJFM
;;      cn_levels_intensity_mean: contor levels for intensity mean plots
;;      cn_levels_intensity_diff: contor levels for intensity mean diff plots
;;  cn_levels_global_domain_mean: contor levels for domain plots
;;  cn_levels_global_domain_diff: contor levels for domain diff plots
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;
;; Caveats
;;   The settings for GridType and Regional are set explicitly in this script
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150702-A_eval_ma: written.
;;
;;#############################################################################

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip_world.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon_panel_positions.ncl"
load "./diag_scripts/aux/SAMonsoon/exact_panel_positions_generic.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "./plot_scripts/ncl/monsoon_domain_panels.ncl"
;*********************************************************************

begin
    ;; ============================================
    ;; ||                                        ||
    ;; ||         READING CONFIGURATION          ||
    ;; ||                                        ||
    ;; ============================================
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    ;; 'diag_script' if fetched from ncl.interface
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot - " + diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    diag_script_base = basename(diag_script)

    ;; ============================================
    ;; ||                                        ||
    ;; ||          REFERENCES INFO               ||
    ;; ||                                        ||
    ;; ============================================
    write_references(diag_script,            \  ; script name
                     "A_eval_ma",               \  ; authors
                     (/"A_levi_ri", "A_mart_gi"/), \  ; contributors
                     "",                     \  ; diag_references
                     "",                     \  ; obs_references
                     "P_embrace")                    ; proj_acknowledgements

    ;; 'show_debuginfo' if fetched from ncl.interface
    if (.not. isdefined("show_debuginfo")) then
        debuginfo = False
    else
        debuginfo = show_debuginfo
    end if

    ;; Check valid plot fields
    valid_plot_field = (/"T2Ms"/)
    ind_f = inlist(field_type0, valid_plot_field)
    if (any(ind_f .eq. False)) then
        print("fatal:can't use field, " + field_type0 + " for this diag_script")
        status_exit(1)
    end if
    delete(valid_plot_field)
    delete(ind_f)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the SAMonsoon cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"summer_season",\
                            "winter_season",\
                            "latrange_global",\
                            "cn_levels_intensity_mean",\
                            "cn_levels_intensity_diff",\
                            "cn_levels_global_domain_mean",\
                            "cn_levels_global_domain_diff",\
                            "high_intensity_cutoff_mmday",\
                            "abs_cmap",\
                            "diff_cmap"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"season", "latrange"/)
    transforms@latrange = di@latrange_global
    min_lat = di@latrange_global(0)
    max_lat = di@latrange_global(1)

    debuginfo@variable = var0
    debuginfo@latrange = min_lat + "-" + max_lat

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

    min_no_models = 2
    max_no_models = 99
    check_min_max_models(dim_MOD, min_no_models, max_no_models, diag_script)

    ;; Ignore some warnings from ESMF-regridding
    ignore_warning_file = "interface_data/warnings_to_ignore.txt"
    warnings_to_ignore = (/"ESMF_regrid_with_weights: warning: destination grid is not completely"/)
    write_ignore_warnings(warnings_to_ignore, ignore_warning_file)



    ;; ============================================
    ;; ||                                        ||
    ;; ||         COMPUTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================

    ;; ---------------------------------------------
    ;; Data structure for storing the models defined
    ;; in the namelist
    ;; ---------------------------------------------
    nml_data_vault = True

    do idx_mod = 0, dim_MOD - 1
        nml_data_vault = read_data(idx_mod, var0, field_type0)
    end do

    ;; ----------------------------------------------
    ;; Data structure for storing the processed field
    ;; ----------------------------------------------
    storage_vault = True

    ;; ----------------------------------------------
    ;; Diagnostic loop, compute diagnostic and store
    ;;                  in 'storage_vault'
    ;; ----------------------------------------------
    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process data
            ;;------------------------
            A0_tmp = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
            A0_tmp@units = ""
            transforms@season = di@summer_season
            A0_summer_nh = dim_avg_n_Wrap(apply_transforms(transforms, A0_tmp), 0)

            transforms@season = di@winter_season
            A0_winter_nh = dim_avg_n_Wrap(apply_transforms(transforms, A0_tmp), 0)

            ;; Summer on northern hemisphere is winter on sourthern...
            A0_summer = A0_summer_nh
            A0_summer({lat|-90:0}, lon|:) = A0_winter_nh({lat|-90:0}, lon|:)

            ;; Winter on northern hemisphere is summer on sourthern...
            A0_winter = A0_winter_nh
            A0_winter({lat|-90:0}, lon|:) = A0_summer_nh({lat|-90:0}, lon|:)
            delete(A0_summer_nh)
            delete(A0_winter_nh)

            A0_avg = dim_avg_n_Wrap(A0_tmp(time|:, {lat|min_lat:max_lat}, lon|:), 0)
            delete(A0_tmp)

            intensity = A0_summer - A0_winter
            delete(A0_winter)

            ;; Add meta data
            intensity@model_name = models@name(idx_mod)
            intensity@idx = idx_mod
            intensity@GridType = "rectilinear"
            intensity@Regional = False
            copy_VarMeta(A0_summer, intensity)

            storage_variable_name = str_join((/"intensity", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(intensity, storage_vault, storage_variable_name)

            summer_pr_above_55_percent = where(A0_summer .gt. A0_avg * 0.55, True, False)
            intensity_above_cutoff = where(intensity .gt. di@high_intensity_cutoff_mmday, True, False)
            delete(A0_avg)
            delete(A0_summer)

            global_domain = where(intensity_above_cutoff .and. summer_pr_above_55_percent, 0.99, -0.01)
            copy_VarMeta(intensity, global_domain)
            delete(summer_pr_above_55_percent)
            delete(intensity_above_cutoff)

            ;; Add some meta data
            global_domain@model_name = models@name(idx_mod)
            global_domain@idx = idx_mod
            global_domain@GridType = "rectilinear"
            global_domain@Regional = False

            ;; Store in data vault
            storage_variable_name = str_join((/"global", "domain", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(global_domain, storage_vault, storage_variable_name)
            delete(intensity)
            delete(global_domain)
        end do  ; idx_mod=0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"intensity", "global-domain"/)
        do idx_var = 0, dimsizes(stored_variables) - 1
            do idx_mod = 0, dim_MOD - 1
                var_name = str_join((/stored_variables(idx_var), sprinti("%i", idx_mod)/), str_vault_sep)
                storage_vault = netcdf_to_vault(storage_vault,\
                                                var_name,\
                                                diag_script_base,\
                                                var0,\
                                                field_type0,\
                                                idx_mod)
            end do
        end do
    end if




    ;; ============================================
    ;; ||                                        ||
    ;; ||          PLOTTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================
    dictkey_units = interface_get_dictkeys(var0, 0) + "_attributes_units"
    lb_units = nml_data_vault@$dictkey_units$

    ;; -------------------------------
    ;; General resources for this plot
    ;; -------------------------------
    res = True
    res = map_zoom(res, transforms)

    ;; Intensity plot
    plot_settings = True
    plot_settings@lb_title = lb_units
    plot_settings@lb_units = lb_units
    plot_settings@idx_ref = idx_ref
    plot_settings@part_of_header = "Precipitation Intensity"

    valid_statistics = (/"yrs"/)

    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ pr_summmer - pr_winter,"\
        + "~C~                 On the northern hemisphere, "\
        + "~C~                 summer = " + di@summer_season + " and winter = " + di@winter_season\
        + "~C~~F22~ Diag ref:~F21~ DOI 10.1007/s00382-011-1266-z"

    plot_settings@type = "mean"
    precip_domain(storage_vault, di, plot_settings, "intensity", debuginfo, valid_statistics, res)

    ;; Diff intensity plot
    plot_settings@type = "diff"
    plot_settings@part_of_header = "Precipitation Intensity: Model minus Reference"
    precip_domain(storage_vault, di, plot_settings, "intensity", debuginfo, valid_statistics, res)

    ;; Monsoon domain plot
    BlueWhite_domain = read_colormap_file(di@abs_cmap)
    plot_settings@default_palette = BlueWhite_domain
    plot_settings@lb_units = "~F22~Domain def:~F21~ pr~B~summer~N~ - pr~B~winter~N~ > "\
                             + flt2string(di@high_intensity_cutoff_mmday)\
                             + " ~F22~.and.~F21~ pr~B~summer~N~ > 0.55*pr~B~annual~N~"
    plot_settings@lbres_lblabels = (/"outside domain", "inside domain"/)
    plot_settings@lbres_lbLabelAlignment = "BoxCenters"

    debuginfo@description = "~F22~Diag def:~F21~ pr_summer - pr_winter > " + flt2string(di@high_intensity_cutoff_mmday)\
        + "~C~         ~F22~and~F21~ pr_summer > 0.55 * pr_annual" \
        + "~C~              On the northern hemisphere, "\
        + di@summer_season + " and winter = " + di@winter_season\
        + "~C~              Blue indicates areas where the above condition is met"\
        + "~C~~F22~ Diag ref:~F21~ DOI 10.1007/s00382-011-1266-z"

    plot_settings@type = "mean"
    plot_settings@part_of_header = "Global Domain"
    plot_settings@skip_statistics_in_plots = True
    precip_domain(storage_vault, di, plot_settings, "global-domain", debuginfo, valid_statistics, res)

    ;; Diff domain plot
    plot_settings@type = "diff"
    plot_settings@part_of_header = "Global Domain: Model minus Reference"
    BlueWhiteRed_domain = read_colormap_file(di@diff_cmap)
    plot_settings@diff_palette = BlueWhiteRed_domain
    plot_settings@lbres_shift_down_maxcolor = True
    plot_settings@lb_units = "Model domain minus reference domain"
    delete(plot_settings@lbres_lblabels)
    plot_settings@lbres_lblabels = (/"only pr~B~REF_domain~N~", "no net value", "only pr~B~model_domain~N~"/)
    precip_domain(storage_vault, di, plot_settings, "global-domain", debuginfo, valid_statistics, res)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
