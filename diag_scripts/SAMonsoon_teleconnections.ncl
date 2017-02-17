;;#############################################################################
;; SAMonsoon_teleconnections.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    This script computes Nino3.4 -- SST (skin temp) correlation
;;
;; Required diag_script_info attributes (diagnostics specific)
;;              monsoon_latrange: latrange crop area for Indian Monsoon
;;              monsoon_lonrange: lonrange crop area for Indian Monsoon
;;                 tele_latrange: latrange crop area for overview area
;;                 tele_lonrange: lonrange crop area for overview area
;;               nino34_latrange: latrange crop area for Nino34 area
;;               nino34_lonrange: lonrange crop area for Nino34 area
;;          cn_levels_mean_basic: contour levels for mean plot
;;     cn_levels_mean_diff_basic: contour levels for mean diff plot
;;        cn_levels_stddev_basic: contour levels for stddev plot
;;   cn_levels_stddev_diff_basic: contour levels for stddev diff plot
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

load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
load "./diag_scripts/aux/SAMonsoon/descriptions.ncl"

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip_tele.ncl"

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
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

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
        print("fatal:can't use field, " + field_type0 + " for this diag_script")
        status_exit(1)
    end if
    delete(valid_plot_field)
    delete(ind_f)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the SAMonsoon cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"tele_latrange",\
                            "tele_lonrange",\
                            "monsoon_latrange",\
                            "monsoon_lonrange",\
                            "nino34_latrange",\
                            "nino34_lonrange",\
                            "cn_levels_mean_basic",\
                            "cn_levels_mean_diff_basic",\
                            "cn_levels_stddev_basic",\
                            "cn_levels_stddev_diff_basic"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;; Get interscection of years
    yrs_intersection = ispan(stringtointeger(models@start_year(0)),\
                             stringtointeger(models@end_year(0)),\
                             1)
    do idx_mod = 1, dim_MOD - 1
        yrs_idx = ispan(stringtointeger(models@start_year(idx_mod)),\
                        stringtointeger(models@end_year(idx_mod)),\
                        1)
        yrs_tmp = intersection(yrs_intersection, yrs_idx)
        delete(yrs_intersection)
        delete(yrs_idx)
        yrs_intersection = yrs_tmp
        delete(yrs_tmp)
    end do
    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"latrange", "lonrange", "years"/)

    transforms@latrange = di@tele_latrange
    transforms@lonrange = di@tele_lonrange
    transforms@years = yrs_intersection

    debuginfo@latrange_monsoon = di@tele_latrange(0) + "-" + di@tele_latrange(1)
    debuginfo@lonrange_monsoon = di@tele_lonrange(0) + "-" + di@tele_lonrange(1)
    debuginfo@latrange_nino34 = di@nino34_latrange(0) + "-" + di@nino34_latrange(1)
    debuginfo@lonrange_nino34 = di@nino34_lonrange(0) + "-" + di@nino34_lonrange(1)
    debuginfo@variable = var0

    min_lat = di@tele_latrange(0)
    max_lat = di@tele_latrange(1)
    min_lon = di@tele_lonrange(0)
    max_lon = di@tele_lonrange(1)

    ;; Reference model index
    idx_ref_pr = get_ref_model_idx(models, var_attr_ref(0))
    idx_ref_ts = get_ref_model_idx(models, var_attr_ref(1))

    ;; Ignore some warnings from ESMF-regridding
    ignore_warning_file = "interface_data/warnings_to_ignore.txt"
    warnings_to_ignore = (/"ESMF_regrid_with_weights: warning: destination grid is not completely"/)
    write_ignore_warnings(warnings_to_ignore, ignore_warning_file)





    ;; ============================================
    ;; ||                                        ||
    ;; ||         COMPUTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================

    ;; -----------------------------------------------
    ;; Data structure for storing the processed fields
    ;; -----------------------------------------------
    storage_vault = True

    ;; -----------------------------------------------
    ;; Reset min_idx/max_idx
    ;; -----------------------------------------------
    min_idx = min((/idx_ref_pr, idx_ref_ts/))
    max_idx = max((/idx_ref_pr, idx_ref_ts/))

    if (.not. read_from_vault) then

        ;; ---------------------------------------------
        ;; Data structure for storing the models defined
        ;; in the namelist
        ;; ---------------------------------------------
        precip_nml_data = True
        ts_nml_data = True

        do idx_mod = 0, dim_MOD - 1
        ;; Skip models<->parameter combinations that are explicitly
        ;; excluded in the xml-file
        curr_var_idx = ind(var0 .eq. variables)
        if (var_attr_exclude(curr_var_idx) .ne. model_attr_id(idx_mod)) then
            precip_nml_data = read_data(idx_mod, var0, field_type0)
        end if
        curr_var_idx = ind(var1 .eq. variables)
        if (var_attr_exclude(curr_var_idx) .ne. model_attr_id(idx_mod)) then
            ts_nml_data = read_data(idx_mod, var1, field_type1)
        end if

        end do

        ;; ----------------------------------------------
        ;; Diagnostic loop, compute diagnostic and store
        ;;         in 'storage_vault'
        ;; ----------------------------------------------
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process data
            ;;------------------------
            ;; Check for reference data sets
            if (any(idx_mod .eq. (/idx_ref_pr, idx_ref_ts/))) then

                ;; Check if there's a single reference data set
                if (idx_ref_pr .eq. idx_ref_ts) then
                    precip_tmp = extract_data_extended(idx_mod, var0, precip_nml_data, -1, 0, 0)
                    ts_tmp = extract_data_extended(idx_mod, var1, ts_nml_data, -1, 0, 0)

                ;; If there are two separate data sets, only use them once
                ;; Skip the second time (max_idx)
                else if (idx_mod .eq. max_idx) then
                    continue
                else  ; Use them the first time (min_idx)
                    precip_tmp = extract_data_extended(idx_ref_pr, var0, precip_nml_data, -1, 0, 0)
                    ts_tmp = extract_data_extended(idx_ref_ts, var1, ts_nml_data, -1, 0, 0)

                    ;; Interpolate onto low res grid
                    precip_grid_size = guestimate_average_grid_area(precip_tmp(0, :, :))
                    ts_grid_size = guestimate_average_grid_area(ts_tmp(0, :, :))

                    weight_file = interface_get_figure_filename(diag_script_base,\
                                                                var0 + var1,\
                                                                field_type0,\
                                                                "weight",\
                                                                idx_mod)
                    source_file = interface_get_figure_filename(diag_script_base,\
                                                                var0 + var1,\
                                                                field_type0,\
                                                                "source",\
                                                                idx_mod)
                    destination_file = interface_get_figure_filename(diag_script_base,\
                                                                     var0 + var1,\
                                                                     field_type0,\
                                                                     "destination",\
                                                                     idx_mod)
                    if (precip_grid_size .gt. ts_grid_size) then
                        local_ts = esmf_conserve_wrapper_time(ts_tmp,\
                                                              precip_tmp,\
                                                              weight_file,\
                                                              source_file,\
                                                              destination_file)
                        local_precip = precip_tmp
                    else
                        local_precip = esmf_conserve_wrapper_time(precip_tmp,\
                                                                  ts_tmp,\
                                                                  weight_file,\
                                                                  source_file,\
                                                                  destination_file)
                        local_ts = ts_tmp
                    end if
                    delete(precip_tmp)
                    delete(ts_tmp)
                    precip_tmp = local_precip
                    ts_tmp = local_ts
                    delete(local_precip)
                    delete(local_ts)
                end if
                end if
            else  ; not a reference data set
                precip_tmp = extract_data_extended(idx_mod, var0, precip_nml_data, -1, 0, 0)
                ts_tmp = extract_data_extended(idx_mod, var1, ts_nml_data, -1, 0, 0)
            end if

            transforms@latrange = di@tele_latrange
            transforms@lonrange = di@tele_lonrange
            precip2 = apply_transforms(transforms, precip_tmp)
            transforms@latrange = di@monsoon_latrange
            transforms@lonrange = di@monsoon_lonrange
            precip1 = apply_transforms(transforms, precip_tmp)

            transforms@latrange = di@nino34_latrange
            transforms@lonrange = di@nino34_lonrange
            precip34 = apply_transforms(transforms, precip_tmp)
            delete(precip_tmp)

            transforms@latrange = di@nino34_latrange
            transforms@lonrange = di@nino34_lonrange
            ts1 = apply_transforms(transforms, ts_tmp)
            delete(ts_tmp)

            precip1_ts = coswgt_areaave(precip1)
            ts1_ts = coswgt_areaave(ts1)
            precip_ts_corr_const = escorc(precip1_ts, ts1_ts)

            precip34_ts = coswgt_areaave(precip34)
            precip34_ts_corr_const = escorc(precip34_ts, ts1_ts)

            precip_ts_corr = escorc(precip2(lat|:, lon|:, time|:), ts1_ts)

            precip_ts_corr!0 = "lat"
            precip_ts_corr!1 = "lon"
            precip_ts_corr&lat = precip2&lat
            precip_ts_corr&lon = precip2&lon

            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------

            precip_ts_corr@idx = idx_mod
            precip_ts_corr@model_name = models@name(idx_mod)
            precip_ts_corr@yrs_intersection = yrs_intersection
            precip_ts_corr@sa34 = precip_ts_corr_const
            precip_ts_corr@en34 = precip34_ts_corr_const
            precip_ts_corr@GridType = "rectilinear"
            precip_ts_corr@Regional = "True"
            precip_ts_corr@GridType = "rectilinear"
            precip_ts_corr@Regional = "True"
            storage_variable_name = str_join((/"mean", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(precip_ts_corr, storage_vault, storage_variable_name)

            delete(precip_ts_corr)
            delete(ts1)
            delete(ts1_ts)
            delete(precip1)
            delete(precip34)
            delete(precip1_ts)
            delete(precip34_ts)
            delete(precip2)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"mean"/)
        do idx_var = 0, dimsizes(stored_variables) - 1
            do idx_mod = 0, dim_MOD - 1
                if (idx_ref_pr .eq. idx_ref_ts) then
                    var_name = str_join((/stored_variables(idx_var), sprinti("%i", idx_mod)/), str_vault_sep)
                else if (idx_mod .eq. max_idx) then
                    continue
                else
                    var_name = str_join((/stored_variables(idx_var), sprinti("%i", idx_mod)/), str_vault_sep)
                end if
                end if
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
    lb_units = "correlation"

    ;; -------------------------------
    ;; General resources for this plot
    ;; -------------------------------
    res = True
    transforms@latrange = di@tele_latrange
    transforms@lonrange = di@tele_lonrange
    res = map_zoom(res, transforms)

    res@tmXBLabelFontHeightF = 0.008
    res@tmYLLabelFontHeightF = 0.008
    res@tiMainFontHeightF = 0.011
    res@mpCenterLonF = 215
    res@mpFillOn = False
    orig_res = res

    ;; Three panel plot (reference, model, and diff)
    plot_settings = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = (/min_idx, max_idx/)
    plot_settings@lb_title = lb_units
    plot_settings@lb_units = lb_units
    plot_settings@part_of_header = "Precipitation"
    plot_settings@type = "Not applicable"
    debuginfo@description_ycoord = 0.8
    debuginfo@description = SAMonsoon_teleconnections_descr(di)

    plot_settings@gsn_right_string = "Precip and SST correlation"
    plot_settings@main_title_suffix = " precip and SST correlation"
    valid_statistics = (/"yrs", "mean", "rmse", "corr", "sa34", "en34"/)

    ;; Set colormaps
    BlueDarkRed18 = read_colormap_file("BlueDarkRed18")
    plot_settings@diff_palette = BlueDarkRed18
    plot_settings@default_palette = "ncl_default"

    ;; Three panel plots (ref, model, diff), here for stddev
    three_by_one_diff(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    res = orig_res
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    plot_settings@type = "mean"
    plot_settings@part_of_header = " precip-Nino3.4 SST teleconnection, mean plots"
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    plot_settings@type = "diff"
    plot_settings@part_of_header = " precip-Nino3.4 SST teleconnection, diff plots"
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
