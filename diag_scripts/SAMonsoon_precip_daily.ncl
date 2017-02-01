;;#############################################################################
;; SAMonsoon_precip_basic_daily.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    JJAS average of standard deviation across indicated years of filtered
;;    daily rainfall. Filtering done by the NCL function band_pass_latlon(...).
;;    Data less than 1 mm/day is masked away.
;;    for the South Asian Monsoon.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;                               season: season in letters, e.g., JJAS
;;                           filter_min: lower limit used for filtering
;;                           filter_max: upper limit used for filtering
;;                       filter_weights: number of filter weights to use
;;                         isv_latrange: latitude crop area
;;                         isv_lonrange: longitude crop area
;;                         bob_latrange: latitude crop area for Bay of Bengal
;;                         bob_lonrange: longitude crop area for Bay of Bengal
;;                         eio_latrange: latitude crop area for eastern
;;                                       equatorial Indian Ocean
;;                         eio_lonrange: longitude crop area for eastern
;;                                       equatorial Indian Ocean
;;          cn_levels_stddev_norm_basic: contour levels for stddev plots
;;     cn_levels_stddev_norm_diff_basic: contour levels for stddev difference
;;                                       plots 
;;
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

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
load "./diag_scripts/aux/SAMonsoon/descriptions.ncl"

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip_isv.ncl"

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
    valid_plot_field = (/"T2Ds"/)
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
    required_attributes = (/"season",\
                            "filter_min",\
                            "filter_max",\
                            "filter_weights",\
                            "isv_latrange",\
                            "isv_lonrange",\
                            "bob_latrange",\
                            "bob_lonrange",\
                            "eio_latrange",\
                            "eio_lonrange",\
                            "cn_levels_stddev_norm_basic",\
                            "cn_levels_stddev_norm_diff_basic"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"latrange",\
                               "lonrange",\
                               "fill_missing_w_poisson"/)

    transforms@season = di@season
    transforms@latrange = di@isv_latrange
    transforms@lonrange = di@isv_lonrange
    transforms@is_cyclic = False

    debuginfo@latrange = di@isv_latrange(0) + "-" + di@isv_latrange(1)
    debuginfo@lonrange = di@isv_lonrange(0) + "-" + di@isv_lonrange(1)
    debuginfo@variable = var0
    debuginfo@filter_min = di@filter_min
    debuginfo@filter_max = di@filter_max
    debuginfo@filter_no_weights = di@filter_weights

    min_lat = di@isv_latrange(0)
    max_lat = di@isv_latrange(1)
    min_lon = di@isv_lonrange(0)
    max_lon = di@isv_lonrange(1)

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

    if (.not. read_from_vault) then
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
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process data
            ;;------------------------
            A0_tmp = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
            A0 = apply_transforms(transforms, A0_tmp)
            delete(A0_tmp)

            precip_jjas_tmp = extract_season(A0, di@season)
            precip_jjas = dim_avg_Wrap(precip_jjas_tmp(lat|:, lon|:, time|:))
            delete(precip_jjas_tmp)

            opt = True
            opt@detrend = True
            srate = 1
            bpf = (/di@filter_min, di@filter_max, di@filter_weights/)
            bpass = band_pass_latlon_time(A0, srate, bpf, opt)

            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------

            season_bpass = extract_season(bpass, di@season)
            var_stddev = dim_stddev_n_Wrap(season_bpass, 0)

            ;; Compute reference statistical values over Bay of Bengal (bob)
            stat_bob_const = coswgt_areaave(var_stddev({lat|di@bob_latrange(0):di@bob_latrange(1)},\
                                                       {lon|di@bob_lonrange(0):di@bob_lonrange(1)}))

            clim_mean = coswgt_areaave(precip_jjas({lat|di@bob_latrange(0):di@bob_latrange(1)},\
                                                   {lon|di@bob_lonrange(0):di@bob_lonrange(1)}))
            stat_bob_normalized = stat_bob_const / clim_mean

            ;; Compute reference statistical values over
            ;; Eastern Equatorial Indian Ocean (eio)
            stat_eio_const = coswgt_areaave(var_stddev({lat|di@eio_latrange(0):di@eio_latrange(1)},\
                                                       {lon|di@eio_lonrange(0):di@eio_lonrange(1)}))

            clim_mean = coswgt_areaave(precip_jjas({lat|di@eio_latrange(0):di@eio_latrange(1)},\
                                                   {lon|di@eio_lonrange(0):di@eio_lonrange(1)}))
            stat_eio_normalized = stat_eio_const / clim_mean

            var_stddev@bob = stat_bob_normalized
            var_stddev@eio = stat_eio_normalized
            var_stddev@idx = idx_mod
            var_stddev@GridType = "rectilinear"
            var_stddev@Regional = True
            var_stddev@model_name = models@name(idx_mod)
            storage_variable_name = str_join((/"stddev", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var_stddev, storage_vault, storage_variable_name)

            ;; Create map normalised by climatological mean, masking out
            ;; points with less than 1mm/day climatological rainfall
            var_stddev_norm = var_stddev  ; Copy meta data
            var_FillValue = getFillValue_with_default(var_stddev_norm)
            precip_jjas_with_fillvalues = where(precip_jjas .lt. 1., var_FillValue, precip_jjas)
            var_stddev_norm = var_stddev / precip_jjas_with_fillvalues

            var_stddev_norm@bob = stat_bob_normalized
            var_stddev_norm@eio = stat_eio_normalized
            var_stddev_norm@idx = idx_mod
            var_stddev_norm@GridType = "rectilinear"
            var_stddev_norm@Regional = True
            var_stddev_norm@model_name = models@name(idx_mod)
            storage_variable_name = str_join((/"stddev_norm", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var_stddev_norm, storage_vault, storage_variable_name)

            delete(A0)
            delete(precip_jjas)
            delete(precip_jjas_with_fillvalues)
            delete(bpass)
            delete(season_bpass)
            delete(var_stddev)
            delete(var_stddev_norm)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"stddev", "stddev-norm"/)
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
    plot_settings@part_of_header = "Precipitation"
    plot_settings@type = "Not applicable"
    debuginfo@description_ycoord = 0.8
    debuginfo@description = SAMonsoon_daily_precip_filter_descr(di)
    filter_period = di@filter_min + "-" + di@filter_max

    plot_settings@gsn_right_string = variable_info@long_name
    plot_settings@prefix_to_header = "Intraseasonal "
    plot_settings@main_title_suffix = " (filtered " + filter_period + " daily rainfall)"
    valid_statistics = (/"yrs", "mean", "rmse", "corr", "bob", "eio"/)

    ;; Three panel plots (ref, model, diff), here for stddev
    three_by_one_diff(storage_vault,\
                      di,\
                      plot_settings,\
                      "stddev",\
                      debuginfo,\
                      valid_statistics,\
                      res)

    ;; Three panel plots (ref, model, diff),
    ;; here for normalized and masked stddev
    plot_settings@main_title_suffix = " (filtered " + filter_period + " daily rainfall)"

    plot_settings@lb_units = "pr/pr~B~clim~N~"
    plot_settings@var_base_length = 6  ; Number of string cols to keep of
                                       ; storage name in header
    debuginfo@description = SAMonsoon_daily_precip_filter_norm_descr(di)
    three_by_one_diff(storage_vault,\
                      di,\
                      plot_settings,\
                      "stddev_norm",\
                      debuginfo,\
                      valid_statistics,\
                      res)

    res = orig_res
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    plot_settings@type = "mean"

    plot_settings@lb_units = lb_units
    plot_settings@part_of_header = " filtered " + filter_period + " daily rainfall"
    plot_settings@var_base_length = 0
    debuginfo@description = SAMonsoon_daily_precip_filter_descr(di)
    multipanel(storage_vault,\
               di,\
               plot_settings,\
               "stddev",\
               debuginfo,\
               valid_statistics,\
               res)

    plot_settings@lb_units = "pr/pr~B~clim~N~ (areas with pr~B~clim~N~<1 mm/day removed)"
    plot_settings@part_of_header = " of filtered " + filter_period + " daily rainfall"
    plot_settings@var_base_length = 6  ; Number of string cols to keep of
                                       ; storage name in header
    debuginfo@description = SAMonsoon_daily_precip_filter_norm_descr(di)
    multipanel(storage_vault,\
               di,\
               plot_settings,\
               "stddev_norm",\
               debuginfo,\
               valid_statistics,\
               res)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
