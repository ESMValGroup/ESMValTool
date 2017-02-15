;;#############################################################################
;; WAMonsoon_precip_seasonal.ncl
;; Author: Romain Roehrig (Meteo France, France)
;; EMBRACE project
;;#############################################################################
;; Description
;;    Computes and outputs seasonal Monsoon precipitation figures. Three types
;;    plot that is considered over the WAM region is: standard deviation of
;;    plots are considered, all concerns precip in the WAM region.
;;      i) Monthly climatologies
;;     ii) Monthly anomaly with respect to each models average precip
;;    iii) Standard deviation of anomaly across all years (a measure of
;;         the monsoon "strength"). This plot is averaged over JJAS.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    cn_levels_mean_basic:         Contour levels for mean plot
;;    cn_levels_mean_diff_basic:    Contour levels for diff plot
;;    cn_levels_stddev_basic:       Contour levels for stddev plot
;;    cn_levels_stddev_diff_basic:  Contour levels for diff stddev plot
;;    latrange_seasonal:            Latitude crop range for plot
;;    lonrange_seasonal:            Longitude crop range for plot
;;    season:                       Season, typically "JJAS"
;;    styleset:                     Stylesheet to use
;;    multi_model_mean:             Average and display average model
;;    my_region:                    Switch to indicate Monsoon region
;;    supporting_gridlines:         Display supporting grid line
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;    long_name:                    Name displayed in plot
;;    skip_long_name_in_subheader:  Logical to switch wether 'long_name'
;;                                  is displayed or not
;;    units:                        Displayed units
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20140515-A_roeh_ro: written
;;
;; #############################################################################

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/auxiliary.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"

exact_panel_position_file = "./diag_scripts/aux/" + diag_script_info@my_region + "Monsoon/exact_panel_positions_" + variables(0) + ".ncl"

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
                     "A_roeh_ro",            \  ; authors
                     "",                     \  ; contributors
                     "",                     \  ; diag_references
                     "",                     \  ; obs_references
                     "P_embrace")               ; proj_acknowledgements

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
    month_tmp = ispan(0, 11, 1)
    month_tmp@long_name = "Month number"

    precip_seasonal_cycle = new((/dim_MOD, 12/), float)
    precip_seasonal_cycle!0 = "model"
    precip_seasonal_cycle!1 = "month"
    precip_seasonal_cycle&model = model_tmp
    precip_seasonal_cycle&month = month_tmp

    precip_seasonal_anom  = new((/dim_MOD, 12/), float)
    precip_seasonal_anom!0 = "model"
    precip_seasonal_anom!1 = "month"
    precip_seasonal_anom&model = model_tmp
    precip_seasonal_anom&month = month_tmp
    delete(month_tmp)

    ;; Loop through models and compute diagnostics
    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1

            ;; Extract and modify data
            A0_full = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
            monthly_climatology = clmMonTLL(A0_full(time|:, {lat|min_lat:max_lat}, {lon|min_lon:max_lon}))
            long_name = monthly_climatology@long_name

            ;; From diag_scripts/lib/ncl/statistics.ncl
            precip_clim = coswgt_areaave(monthly_climatology)
            precip_seasonal_cycle(idx_mod, :) = precip_clim
            precip_ave = avg(precip_clim)
            precip_seasonal_anom(idx_mod, :) = (/precip_clim - precip_ave/)

            delete(precip_ave)
            delete(monthly_climatology)
            delete(A0_full)

            ;; Check if start/end years are consistent across all models
            if ((first_start_year .ne. models@start_year(idx_mod))\
                .or. (first_end_year .ne. models@end_year(idx_mod))) then
                model_years_consistent = False
            end if
        end do  ; idx_mod = 0, dim_MOD - 1

        precip_seasonal_anom@units = precip_clim@units

        storage_variable_name = str_join((/"precip", "seasonal", "cycle"/), str_vault_sep)
        store_variable(precip_seasonal_cycle, storage_vault, storage_variable_name)

        storage_variable_name = str_join((/"precip", "seasonal", "anom"/), str_vault_sep)
        store_variable(precip_seasonal_anom, storage_vault, storage_variable_name)

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"precip-seasonal-cycle",\
                             "precip-seasonal-anom"/)
        do idx_var = 0, dimsizes(stored_variables) - 1
            var_name = stored_variables(idx_var)
            storage_vault = netcdf_to_vault(storage_vault,\
                                            var_name,\
                                            diag_script_base,\
                                            var0,\
                                            field_type0,\
                                            -1)
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
    plot_settings@aux_filename_info = di@season + "-seasonal-cycle"
    debuginfo@description_ycoord = 0.7
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology across specified years and"\
                                + "~C~   averaged over the area,"\
                                + "~C~      lat: " + min_lat + "-" + max_lat\
                                + ", lon: " + min_lon + "-" + max_lon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"
    plot_settings@title_string = "Seasonal cycle"
    plot_settings@yaxis_string = "Precipitation"
    plot_settings@type = "monthly_climatology"
    plot_settings@left_plot_subheader = latlon_header_string + ", " + years_string

    xy_plot_wrapper(storage_vault, di, plot_settings, "precip-seasonal-cycle", debuginfo)

    debuginfo@description = "~F22~Diag def:~F21~ Monthly anomalies across specified years and averaged"\
                                + "~C~   over the area,"\
                                + "~C~      lat: " + min_lat + "-" + max_lat\
                                + ", lon: " + min_lon + "-" + max_lon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"
    plot_settings@title_string = "Seasonal anomaly"
    plot_settings@aux_filename_info = di@season + "-seasonal-anom"
    xy_plot_wrapper(storage_vault, di, plot_settings, "precip-seasonal-anom", debuginfo)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
