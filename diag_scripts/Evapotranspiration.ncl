;;#############################################################################
;; Evapotranspiration.ncl
;; Author: Martin Evaldsson (SMHI, Sweden)
;;         Edouard Davit (ETHZ, Switzerland)
;;         Boris Orlowski (former ETHZ, Switzerland)
;; EMBRACE project
;;#############################################################################
;; Description
;; Computes and plots contour plots for global monthly climatology
;; for evapotranspiration. An example of this diagnostic on CMIP5 models
;; is presented in Mueller and Seneviratne, 2014,
;;
;;      http://onlinelibrary.wiley.com/doi/10.1002/2013GL058055/abstract
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    et_colormap:                  Colormap to use for the absolute values
;;    cn_levels_mean_0_basic:       Contour levels - January   absolute values
;;    cn_levels_mean_1_basic:       Contour levels - February  absolute values
;;    cn_levels_mean_2_basic:       Contour levels - March     absolute values
;;    cn_levels_mean_3_basic:       Contour levels - April     absolute values
;;    cn_levels_mean_4_basic:       Contour levels - May       absolute values
;;    cn_levels_mean_5_basic:       Contour levels - June      absolute values
;;    cn_levels_mean_6_basic:       Contour levels - July      absolute values
;;    cn_levels_mean_7_basic:       Contour levels - August    absolute values
;;    cn_levels_mean_8_basic:       Contour levels - September absolute values
;;    cn_levels_mean_9_basic:       Contour levels - October   absolute values
;;    cn_levels_mean_10_basic:      Contour levels - November  absolute values
;;    cn_levels_mean_11_basic:      Contour levels - December  absolute values
;;
;;    cn_levels_mean_0_diff_basic:  Contour levels - January   diff values
;;    cn_levels_mean_1_diff_basic:  Contour levels - February  diff values
;;    cn_levels_mean_2_diff_basic:  Contour levels - March     diff values
;;    cn_levels_mean_3_diff_basic:  Contour levels - April     diff values
;;    cn_levels_mean_4_diff_basic:  Contour levels - May       diff values
;;    cn_levels_mean_5_diff_basic:  Contour levels - June      diff values
;;    cn_levels_mean_6_diff_basic:  Contour levels - July      diff values
;;    cn_levels_mean_7_diff_basic:  Contour levels - August    diff values
;;    cn_levels_mean_8_diff_basic:  Contour levels - September diff values
;;    cn_levels_mean_9_diff_basic:  Contour levels - October   diff values
;;    cn_levels_mean_10_diff_basic: Contour levels - November  diff values
;;    cn_levels_mean_11_diff_basic: Contour levels - December  diff values
;;
;;    sftlf_cutoff:                 Cutoff percentage to use in land sea mask
;;                                  (>50 -> grid cell is land)
;;
;; Required variable_info attributes (variable specific)
;;    long_name:                    Name to store in derived varible
;;    units:                        Units for derived variable
;;
;; Caveats
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150520-A_eval_ma: written
;; ############################################################################

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/aux/SAMonsoon/descriptions.ncl"

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip.ncl"

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
    write_references(diag_script,                  \  ; script name
                     (/"A_eval_ma"/),              \  ; authors
                     (/"A_davi_ed", "A_orlo_bo"/), \  ; contributors
                     "",                           \  ; diag_references
                     "",                           \  ; obs_references
                     "P_embrace")                     ; proj_acknowledgements
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
    required_attributes = (/"cn_levels_mean_basic",\
                            "cn_levels_mean_diff_basic"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True

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
    months = (/"Jan", "Feb", "Mar",\
               "Apr", "May", "Jun",\
               "Jul", "Aug", "Sep",\
               "Oct", "Nov", "Dec"/)
    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process data
            ;;------------------------
            A0 = extract_data(idx_mod, nml_data_vault, -1, 0, 0)

            ;; Check if sftlf fx file is defined
            var_fx = models@case_name(idx_mod)
            fx_idx = ind(var_fx .eq. fx_keys)
            if (.not. ismissing(fx_idx)) then
                fx_file = addfile(fx_values(fx_idx), "r")
                sftlf = fx_file->sftlf
            end if


            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------
            var_mean_all = clmMonTLL(A0)
            do idx_month = 0, dimsizes(months) - 1
                var_mean_fullfield = var_mean_all(idx_month, :, :)
                var_mean = var_mean_fullfield

                ;; If sftlf fx file is defined, apply it
                if (.not. ismissing(fx_idx)) then
                    var_mean = where(sftlf .gt. di@sftlf_cutoff,\
                                     var_mean_fullfield,\
                                     var_mean_fullfield@_FillValue)
                end if

                var_mean@idx = idx_mod
                var_mean@model_name = models@name(idx_mod)
                delete(var_mean@month)
                var_mean@month = months(idx_month)
                var_mean@GridType = "rectilinear"
                var_mean@Regional = False

                idx_mon_str = sprinti("%i", idx_month)
                idx_mod_str = sprinti("%i", idx_mod)
                storage_variable_name = str_join((/"mean",\
                                                   idx_mon_str,\
                                                   idx_mod_str/),\
                                                 str_vault_sep)
                store_variable(var_mean, storage_vault, storage_variable_name)

                delete(var_mean)
                delete(var_mean_fullfield)
            end do
            delete(A0)
            delete(var_mean_all)
            if (.not. ismissing(fx_idx)) then
                delete(sftlf)
            end if
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        do idx_mod = 0, dim_MOD - 1
            do idx_month = 0, dimsizes(months) - 1
                idx_mon_str = sprinti("%i", idx_month)
                idx_mod_str = sprinti("%i", idx_mod)
                var_name = str_join((/"mean",\
                                      idx_mon_str,\
                                      idx_mod_str/),\
                                    str_vault_sep)
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
    plot_settings@part_of_header = "Evapotranspiration"
    plot_settings@type = "Not applicable"
    debuginfo@description_ycoord = 0.8
    debuginfo@description = description_for_evapotranspiration(di)

    plot_settings@gsn_right_string = variable_info@long_name

    ;; Colormap settings
    ;; Reverse an existing palette for diff plots
    BlueDarkRed18 = read_colormap_file("BlueDarkRed18")
    BlueDarkRed18 = BlueDarkRed18(::-1, :)
    plot_settings@diff_palette = BlueDarkRed18
    evap_colormap = read_colormap_file(di@et_colormap)
    plot_settings@default_palette = evap_colormap

    ;; Three panel plots (ref, model, diff)
    valid_statistics = (/"yrs", "mean", "rmse", "corr"/)
    var_name_base = "mean"
    var_base_length = strlen(var_name_base)
    plot_settings@var_base_length = var_base_length

    do idx_month = 0, dimsizes(months) - 1
        di@season = months(idx_month)
        debuginfo@month = months(idx_month)
        idx_mon_str = sprinti("%i", idx_month)
        var_name = str_join((/"mean", idx_mon_str/), str_vault_sep)
        three_by_one_diff(storage_vault,\
                          di,\
                          plot_settings,\
                          var_name,\
                          debuginfo,\
                          valid_statistics,\
                          res)
    end do

    res = orig_res
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    do idx_month = 0, dimsizes(months) - 1
        di@season = months(idx_month)
        debuginfo@month = months(idx_month)
        idx_mon_str = sprinti("%i", idx_month)
        var_name = str_join((/"mean", idx_mon_str/), str_vault_sep)

        plot_settings@type = "mean"
        plot_settings@part_of_header = " of " + variable_info@long_name
        multipanel(storage_vault,\
                   di,\
                   plot_settings,\
                   var_name,\
                   debuginfo,\
                   valid_statistics,\
                   res)

        plot_settings@type = "diff"
        plot_settings@part_of_header = " of Evapotranspiration"
        multipanel(storage_vault,\
                   di,\
                   plot_settings,\
                   var_name,\
                   debuginfo,\
                   valid_statistics,\
                   res)
    end do

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
