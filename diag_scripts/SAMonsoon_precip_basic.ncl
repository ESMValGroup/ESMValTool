;;#############################################################################
;; SAMonsoon_precip_basic.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    This script computes contors a mean and standard deviation for
;;    for the South Asian Monsoon over the indicated season (typcically JJAS)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;                          season: season in letters, e.g., JJAS
;;                  latrange_basic: latitude crop area
;;                  lonrange_basic: longitude crop area
;;            cn_levels_mean_basic: contour levels for mean plots
;;       cn_levels_mean_diff_basic: contour levels for mean difference plots
;;          cn_levels_stddev_basic: contour levels for stddev plots
;;     cn_levels_stddev_diff_basic: contour levels for stddev difference plots
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

exact_panel_position_file = "./diag_scripts/aux/" + diag_script_info@my_region + "Monsoon/exact_panel_positions_wind.ncl"

load "./plot_scripts/ncl/monsoon_panels.ncl"
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
end

; Load variable_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
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
        error_msg("fatal", diag_script, "", "can't use field " + field_type0)
        status_exit(1)
    end if
    delete(valid_plot_field)
    delete(ind_f)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the SAMonsoon cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"season",\
                            "latrange_basic",\
                            "lonrange_basic",\
                            "cn_levels_mean_basic",\
                            "cn_levels_mean_diff_basic",\
                            "cn_levels_stddev_basic",\
                            "cn_levels_stddev_diff_basic"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"season", "latrange", "lonrange"/)
    transforms@season = di@season
    transforms@latrange = di@latrange_basic
    transforms@lonrange = di@lonrange_basic

    debuginfo@latrange = di@latrange_basic(0) + "-" + di@latrange_basic(1)
    debuginfo@lonrange = di@lonrange_basic(0) + "-" + di@lonrange_basic(1)
    debuginfo@variable = var0

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

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

    ;; -----------------------------------------------
    ;; Data structure for storing the processed fields
    ;; -----------------------------------------------
    storage_vault = True

    ;; ----------------------------------------------
    ;; Diagnostic loop, compute diagnostic and store
    ;;         in 'storage_vault'
    ;; ----------------------------------------------
    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process data
            ;;------------------------
            A0_tmp = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
            A0 = apply_transforms(transforms, A0_tmp)
            delete(A0_tmp)

            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------
            var_mean = dim_avg_n_Wrap(A0, 0)
            var_mean@idx = idx_mod
            var_mean@model_name = models@name(idx_mod)
            var_mean@GridType = "rectilinear"
            var_mean@Regional = True
            storage_variable_name = str_join((/"mean", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var_mean, storage_vault, storage_variable_name)

            ;; ------------------------------------------
            ;; Compute and store stddev for model idx_mod
            ;; ------------------------------------------
            var_stddev = dim_stddev_n_Wrap(A0, 0)
            var_stddev@idx = idx_mod
            var_stddev@model_name = models@name(idx_mod)
            var_stddev@GridType = "rectilinear"
            var_stddev@Regional = True
            storage_variable_name = str_join((/"stddev", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var_stddev, storage_vault, storage_variable_name)

            delete(A0)
            delete(var_mean)
            delete(var_stddev)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"mean", "stddev"/)
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

    res@tmXBLabelFontHeightF = 0.008
    res@tmYLLabelFontHeightF = 0.008
    res@mpFillOn = False
    orig_res = res

    ;; Three panel plot (reference, model, and diff)
    plot_settings = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_title = lb_units
    plot_settings@lb_units = lb_units
    plot_settings@part_of_header = "Precipitation"
    plot_settings@type = "Not applicable"
    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ " + di@season + " average across indicated years" +\
                            "~C~~F22~rmse def:~F21~ NCL function wgt_arearmse(...) computed with " +\
                            "~C~                  cosine weights in latitude."+\
                            "~C~~F22~corr def:~F21~ NCL function pattern_cor(...) computed with " +\
                            "~C~                  cosine weights in latitude. Correlation is between" +\
                            "~C~                  reference and model data."

    plot_settings@gsn_right_string = variable_info@long_name

    ;; Three panel plots (ref, model, diff)
    valid_statistics = (/"yrs", "mean", "rmse", "corr"/)

    three_by_one_diff(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)
    plot_settings@prefix_to_header = "Interannual "
    three_by_one_diff(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    res = orig_res
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    plot_settings@type = "mean"
    plot_settings@part_of_header = " of " + variable_info@long_name
    plot_settings@prefix_to_header = ""
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    plot_settings@type = "diff"
    plot_settings@part_of_header = " of Precipitation"
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    plot_settings@type = "mean"
    plot_settings@prefix_to_header = "Interannual "
    plot_settings@part_of_header = " stddev of Precipitation"
    multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    plot_settings@type = "diff"
    plot_settings@prefix_to_header = "Interannual "
    plot_settings@part_of_header = " of stddev of Precipitation"
    multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
