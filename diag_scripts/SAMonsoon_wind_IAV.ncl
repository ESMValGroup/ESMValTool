;;#############################################################################
;; SAMonsoon_wind_IAV.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;; Computes and outputs the interannual variability of the
;; Goswami, Webster-Yang and Wang-Fan monsoon indicies.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;       season: season in letters, e.g., JJAS
;;     styleset: chose styles set for lines
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;
;; Caveats
;;   Only valid for prescribed SST experiments (AMIP)
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

load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_wind.ncl"

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
    info_output("plot - " + diag_script + " (var: " + variables(0) \
                          + ", " + variables(1) + ")", verbosity, 1)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

;    ;;-----------------------
;    ;; References Information
;    ;;-----------------------
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
    nml_data_vault_1 = True
    nml_data_vault_2 = True
    do idx_mod = 0, dim_MOD - 1
        nml_data_vault_1 = read_data(idx_mod, var0, field_type0)
        nml_data_vault_2 = read_data(idx_mod, var1, field_type1)
    end do

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
                            "styleset"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;; Read all models and extract number of years
    data = True
    years = (/-999.0/)
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
        tmp_model = extract_data(idx_mod, data, -1, 0, 0)

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

    ;;----------------------
    ;; Data arrays for plots
    ;;----------------------
    ;; Interannual variability
    goswami_iav      = new((/dim_MOD, no_years/), float)
    webster_yang_iav = new((/dim_MOD, no_years/), float)
    wang_fan_iav     = new((/dim_MOD, no_years/), float)

    ;; Add dimensions to interannual variability plots
    model_list = NewList("lifo")
    ListPush(model_list, goswami_iav)
    ListPush(model_list, webster_yang_iav)
    ListPush(model_list, wang_fan_iav)

    model_tmp = models@name
    model_tmp@long_name = "Model names"

    do idx = 0, ListCount(model_list) - 1
        model_list[idx] = 0.0
        model_list[idx]!0 = "model"
        model_list[idx]!1 = "years"
        model_list[idx]&model = model_tmp
        model_list[idx]&years = years
    end do
    delete(model_list)

    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and modify data
            ;;------------------------
            A0_ua_full = extract_data(idx_mod, nml_data_vault_1, -1, 0, 0)
            A0_va_full = extract_data(idx_mod, nml_data_vault_2, -1, 0, 0)

            ;; ------------------------------------------------------
            ;; First monsoon index (Monsoon Hadley Index)
            ;; ------------------------------------------------------
            ;; Goswami et al. A broad-scale circulation index for the
            ;; interannual variability of the Indian summer monsoon,
            ;; Quarterly Journal of the Royal Meteorological Society,
            ;; vol. 125, issue 554, pp. 611-633
            ;; ------------------------------------------------------
            gos_minlat = 10
            gos_maxlat = 30
            gos_minlon = 70
            gos_maxlon = 110

            goswami_200_jjas_seasonal = month_to_season_extended(A0_va_full(time|:,\
                                                                            {plev|20000},\
                                                                            {lat|gos_minlat:gos_maxlat},\
                                                                            {lon|gos_minlon:gos_maxlon}),\
                                                                 di@season)
            goswami_200_iav = coswgt_areaave(goswami_200_jjas_seasonal)

            goswami_850_jjas_seasonal = month_to_season_extended(A0_va_full(time|:,\
                                                                            {plev|85000},\
                                                                            {lat|gos_minlat:gos_maxlat},\
                                                                            {lon|gos_minlon:gos_maxlon}),\
                                                                 di@season)
            goswami_850_iav = coswgt_areaave(goswami_850_jjas_seasonal)
            clim_years = fetch_first_last_year(goswami_850_iav)
            goswami_iav(idx_mod, {clim_years(0):clim_years(1)})\
                = dim_standardize(goswami_850_iav - goswami_200_iav, 0)
            goswami_iav@units = goswami_850_iav@units

            set_missing_years(goswami_iav, clim_years, di@years, idx_mod)

            ;; ------------------------------------------------------
            ;; Second monsoon index (Webster-Yang monsoon index)
            ;; ------------------------------------------------------
            ;; Webster, P. J. and S. Yang, 1992: Monsoon and ENSO:
            ;; Selectively interactive systems. Quart. J. Roy. Meteor.
            ;; Soc., 118, 877-926
            ;; ------------------------------------------------------
            wy_minlat = 0
            wy_maxlat = 20
            wy_minlon = 40
            wy_maxlon = 110

            webster_yang_200_seasonal = month_to_season_extended(A0_ua_full(time|:,\
                                                                            {plev|20000},\
                                                                            {lat|wy_minlat:wy_maxlat},\
                                                                            {lon|wy_minlon:wy_maxlon}),\
                                                                 di@season)
            webster_yang_200_iav = coswgt_areaave(webster_yang_200_seasonal)

            webster_yang_850_seasonal = month_to_season_extended(A0_ua_full(time|:,\
                                                                            {plev|85000},\
                                                                            {lat|wy_minlat:wy_maxlat},\
                                                                            {lon|wy_minlon:wy_maxlon}),\
                                                                 di@season)
            webster_yang_850_iav = coswgt_areaave(webster_yang_850_seasonal)
            clim_years = fetch_first_last_year(webster_yang_850_iav)

            webster_yang_iav(idx_mod, {clim_years(0):clim_years(1)})\
                = dim_standardize(webster_yang_850_iav - webster_yang_200_iav, 0)
            webster_yang_iav@units = webster_yang_850_iav@units

            set_missing_years(webster_yang_iav, clim_years, di@years, idx_mod)

            ;; ------------------------------------------------------
            ;; Third  monsoon index (Dynamical Monsoon index)
            ;; ------------------------------------------------------
            ;; Wang, B., and Z. Fan, 1999: Choice of south Asian summer
            ;; monsoon indices.
            ;; Bull. Amer. Meteor. Soc., 80, 629-638.
            ;; ------------------------------------------------------
            wf1_minlat = 5
            wf1_maxlat = 15
            wf1_minlon = 40
            wf1_maxlon = 80

            wf2_minlat = 20
            wf2_maxlat = 30
            wf2_minlon = 70
            wf2_maxlon = 90

            wang_fan_850_1_seasonal = month_to_season_extended(A0_ua_full(time|:,\
                                                                          {plev|85000},\
                                                                          {lat|wf1_minlat:wf1_maxlat},\
                                                                          {lon|wf1_minlon:wf1_maxlon}),\
                                                               di@season)
            wang_fan_850_1_iav = coswgt_areaave(wang_fan_850_1_seasonal)

            wang_fan_850_2_seasonal = month_to_season_extended(A0_ua_full(time|:,\
                                                                          {plev|85000},\
                                                                          {lat|wf2_minlat:wf2_maxlat},\
                                                                          {lon|wf2_minlon:wf2_maxlon}),\
                                                               di@season)
            wang_fan_850_2_iav = coswgt_areaave(wang_fan_850_2_seasonal)
            clim_years = fetch_first_last_year(wang_fan_850_2_iav)
            wang_fan_iav(idx_mod, {clim_years(0):clim_years(1)})\
                = dim_standardize(wang_fan_850_1_iav - wang_fan_850_2_iav, 0)
            wang_fan_iav@units = wang_fan_850_2_iav@units

            set_missing_years(wang_fan_iav, clim_years, di@years, idx_mod)

            delete(A0_ua_full)
            delete(A0_va_full)

            delete(goswami_200_iav)
            delete(goswami_850_iav)
            delete(goswami_200_jjas_seasonal)
            delete(goswami_850_jjas_seasonal)
            delete(wang_fan_850_1_seasonal)
            delete(wang_fan_850_2_seasonal)
            delete(wang_fan_850_1_iav)
            delete(wang_fan_850_2_iav)
            delete(webster_yang_200_seasonal)
            delete(webster_yang_850_seasonal)
            delete(webster_yang_200_iav)
            delete(webster_yang_850_iav)

            ;; Check if start/end years are consistent across all models
            if ((first_start_year .ne. models@start_year(idx_mod))\
                .or. (first_end_year .ne. models@end_year(idx_mod))) then
                model_years_consistent = False
            end if
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Set years_string depending on whether all
        ;; models include the same years
        if (model_years_consistent) then
            years_string = "~F22~yrs:~F21~ " + models@start_year(0) + "-" + models@end_year(0)
        else
            years_string = "~F22~yrs:~F21~ inconsistent"
        end if
        di@years_string = years_string

        storage_variable_name = "goswami-iav"
        store_variable(goswami_iav, storage_vault, storage_variable_name)

        storage_variable_name = "webster-yang-iav"
        store_variable(webster_yang_iav, storage_vault, storage_variable_name)

        storage_variable_name = "wang-fan-iav"
        store_variable(wang_fan_iav, storage_vault, storage_variable_name)

        delete(goswami_iav)
        delete(webster_yang_iav)
        delete(wang_fan_iav)

        ;; Write vault to netCDF
        if (write_netcdf) then
            diag_script_base = basename(diag_script)
            vault_to_netcdf(storage_vault,\
                            diag_script_base,\
                            var0 + var1,\
                            field_type0 + field_type1)
        end if
    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"goswami-iav", "webster-yang-iav", "wang-fan-iav"/)
        do idx_var = 0, dimsizes(stored_variables) - 1
            var_name = str_join((/stored_variables(idx_var), sprinti("%i", idx_mod)/), str_vault_sep)
            storage_vault = netcdf_to_vault(storage_vault,\
                                            var_name,\
                                            diag_script_base,\
                                            var0,\
                                            field_type0,\
                                            -1)
        end do
    end if





    ;; ============================================
    ;; ||                                        ||
    ;; ||          PLOTTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================
    plot_settings = True

    debuginfo@description_ycoord = 0.7
    plot_settings@title_string = "Goswami IAV (" + di@season + ")"
    plot_settings@yaxis_string = "Standardized anomalies"
    plot_settings@xaxis_string = "years"
    plot_settings@type = "iav"
    plot_settings@aux_filename_info = "goswami-" + di@season + "-IAV"
    plot_settings@left_plot_subheader = years_string
    debuginfo@description = "~F22~Diag def:~F21~ Interannual variability of standardized anomalies"\
                                + "~C~   (normalized with climatology) for the Goswami Monsoon Index"\
                                + "~C~   (GMI). GMI is the vertical northward wind shear for the area,"\
                                + "~C~      lat: " + gos_minlat + "-" + gos_maxlat\
                                + ", lon: " + gos_minlon + "-" + gos_maxlon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~GMI ref:~F21~ DOI: 10.1002/qj.49712555412"\
                                + "~C~~F22~Caveat:~F21~ Should only be used for comparing AMIP runs with"\
                                + "~C~   the same prescribed SST"
    ts_line_wrapper(storage_vault, di, plot_settings, "goswami-iav", debuginfo)

; --------------------------

    plot_settings@title_string = "Webster Yang IAV (" + di@season + ")"
    plot_settings@aux_filename_info = "webster-yang-" + di@season + "-IAV"
    debuginfo@description = "~F22~Diag def:~F21~ Interannual variability of standardized anomalies"\
                                + "~C~   (normalized with climatology) for the Webster Yang Monsoon Index"\
                                + "~C~   (WYI). WYI is the vertical eastward wind shear for the area,"\
                                + "~C~      lat: " + wy_minlat + "-" + wy_maxlat\
                                + ", lon: " + wy_minlon + "-" + wy_maxlon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~WYI ref:~F21~ DOI: 10.1002/qj.49711850705"\
                                + "~C~~F22~Caveat:~F21~ Should only be used for comparing AMIP runs with"\
                                + "~C~   the same prescribed SST"
    ts_line_wrapper(storage_vault, di, plot_settings, "webster-yang-iav", debuginfo)

; --------------------------

    plot_settings@title_string = "Wang Fan IAV (" + di@season + ")"
    plot_settings@aux_filename_info = "wang-fan-" + di@season + "-IAV"
    debuginfo@description = "~F22~Diag def:~F21~ Interannual variability of standardized anomalies"\
                                + "~C~   (normalized with climatology) for the Wang Fan Monsoon Index"\
                                + "~C~   (WFI) across specified years. WYI is the eastward wind"\
                                + "~C~   speed difference between the areas,"\
                                + "~C~      lat: " + wf1_minlat + "-" + wf1_maxlat\
                                + ", lon: " + wf1_minlon + "-" + wf1_maxlon + " and"\
                                + "~C~      lat: " + wf2_minlat + "-" + wf2_maxlat\
                                + ", lon: " + wf2_minlon + "-" + wf2_maxlon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~Diag ref:~F21~ DOI: 10.1175/1520-0477(1999)080<0629:COSASM>2.0.CO;2"\
                                + "~C~~F22~Caveat:~F21~ Should only be used for comparing AMIP runs with"\
                                + "~C~   the same prescribed SST"
    ts_line_wrapper(storage_vault, di, plot_settings, "wang-fan-iav", debuginfo)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
