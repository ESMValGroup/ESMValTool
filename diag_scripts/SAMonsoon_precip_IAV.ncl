;;#############################################################################
;; SAMonsoon_precip_IAV.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    This script computes Standard deviation of anomaly, averaged over
;;    JJAS, across all years providing a measure of the interannual
;;    monsoon "strength".
;;
;; Required diag_script_info attributes (diagnostics specific)
;;                season: Season in letters, e.g., JJAS
;;              styleset: Styleset to use for lines (color/dash/etc..)
;;     latrange_seasonal: Latitude crop area
;;     lonrange_seasonal: Longitude crop area
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;   Only valid for prescribed SST experiments (AMIP)
;;
;; Caveats
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

    ;;-----------------------
    ; References Information
    ;;-----------------------
    write_references(diag_script,            \  ; script name
                     "A_eval_ma",               \  ; authors
                     (/"A_levi_ri", "A_mart_gi"/), \  ; contributors
                     "",                     \  ; diag_references
                     "",                     \  ; obs_references
                     "P_embrace")                    ; proj_acknowledgements

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "PS"
    end if

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
    required_attributes = (/"season",\
                            "styleset",\
                            "latrange_seasonal",\
                            "lonrange_seasonal"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;------------------------------------------
    ;; Define meta data to debuginfo onto the plot
    ;;------------------------------------------
    min_lat = di@latrange_seasonal(0)
    max_lat = di@latrange_seasonal(1)
    debuginfo@latrange = min_lat + "-" + max_lat
    min_lon = di@lonrange_seasonal(0)
    max_lon = di@lonrange_seasonal(1)
    debuginfo@lonrange = min_lon + "-" + max_lon
    latlon_header_string = get_latlon_string(min_lat, max_lat, min_lon, max_lon)

    ;; Read all models and extract number of years
    nml_data_vault = True
    years = (/-999.0/)
    do idx_mod = 0, dim_MOD - 1
        nml_data_vault = read_data(idx_mod, var0, field_type0)
        tmp_model = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
        years_tmp = set_inclusive_OR(years, fetch_years(tmp_model))
        delete(years)
        years = years_tmp
        delete(years_tmp)

        delete(tmp_model)
    end do
    no_years = dimsizes(years)
    di@years = years

    model_years_consistent = True
    first_start_year = models@start_year(0)
    first_end_year = models@end_year(0)





    ;; ============================================
    ;; ||                                        ||
    ;; ||         COMPUTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================
    ;; -----------------------------------------------
    ;; Data structure for storing the processed fields
    ;; -----------------------------------------------
    storage_vault = True

    ;; Diagnostic data arrays
    model_tmp = models@name
    model_tmp@long_name = "Model names"
    month_tmp = ispan(1, 12, 1)
    month_tmp@long_name = "Month number"

    precip_iav   = new((/dim_MOD, no_years/), float)
    precip_iav   = 0.0
    precip_iav!0 = "model"
    precip_iav!1 = "years"
    precip_iav&years = years
    precip_iav&model = model_tmp
    delete(model_tmp)

    ;; Loop through models and compute diagnostics
    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1

            ;; Extract and modify data
            A0_full = extract_data(idx_mod, nml_data_vault, -1, 0, 0)

            season = "JJAS"
            precip_jjas_seasonal = month_to_season_extended(A0_full(time|:,\
                                                                    {lat|min_lat:max_lat},\
                                                                    {lon|min_lon:max_lon}),\
                                                            season)
            precip_jjas_clim = coswgt_areaave(precip_jjas_seasonal)
            clim_years = fetch_first_last_year(precip_jjas_clim)
            precip_iav(idx_mod, {clim_years(0):clim_years(1)}) = (/dim_standardize_Wrap(precip_jjas_clim, 0)/)
            precip_iav@units = precip_jjas_clim@units

            set_missing_years(precip_iav, clim_years, di@years, idx_mod)

            delete(precip_jjas_clim)
            delete(precip_jjas_seasonal)
            delete(A0_full)

            ;; Check if start/end years are consistent across all models
            if ((first_start_year .ne. models@start_year(idx_mod))\
                 .or. (first_end_year .ne. models@end_year(idx_mod))) then
                model_years_consistent = False
            end if
        end do  ; idx_mod = 0, dim_MOD - 1

        storage_variable_name = str_join((/"precip", "iav"/), str_vault_sep)
        store_variable(precip_iav, storage_vault, storage_variable_name)

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"precip-iav"/)
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

    ;; Set years_string depending on whether all models include the same years
    if (model_years_consistent) then
        years_string = "~F22~yrs:~F21~ " + models@start_year(0) + "-" + models@end_year(0)
    else
        years_string = "~F22~yrs:~F21~ inconsistent"
    end if
    di@left_plot_subheader = latlon_header_string + ", " + years_string
    di@years_string = years_string





    ;; ============================================
    ;; ||                                        ||
    ;; ||          PLOTTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================
    plot_settings = True
    plot_settings@yaxis_string = "Standardized anomalies"
    debuginfo@description_ycoord = 0.7
    debuginfo@description = "~F22~Diag def:~F21~ Standardized anomalies of Precipitation, normalized"\
                                + "~C~   with climatology across specified years and averaged over the"\
                                + "~C~   area,"\
                                + "~C~      lat: " + min_lat + "-" + max_lat\
                                + ", lon: " + min_lon + "-" + max_lon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~Caveat:~F21~ Should only be used for comparing AMIP runs with"\
                                + "~C~   the same prescribed SST"
    plot_settings@xaxis_string = "years"
    plot_settings@title_string = "IAV Precipitation (" + di@season + ")"
    plot_settings@type = "iav"
    plot_settings@aux_filename_info = di@season + "-IAV"
    plot_settings@left_plot_subheader = latlon_header_string + ", " + years_string
    ts_line_wrapper(storage_vault, di, plot_settings, "precip-iav", debuginfo)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
