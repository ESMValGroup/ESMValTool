;;#############################################################################
;; SAMonsoon_wind_seasonal.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    Computes the seasonal cycle of the Goswami, Webster-Yang and Wang-Fan
;;    monsoon indicies
;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;; Optional diag_script_info attributes (diagnostics specific)
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

    diag_script_base = basename(diag_script)

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
        error_msg("fatal", diag_script, "", "can't use field " + field_type0)
        status_exit(1)
    end if
    delete(ind_f)
    ind_f = inlist(field_type1, valid_plot_field)
    if (any(ind_f) .eq. False) then
        error_msg("fatal", diag_script, "", "can't use field " + field_type1)
        status_exit(1)
    end if
    delete(ind_f)
    delete(valid_plot_field)

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
    ;; Monthly climatology
    goswami          = new((/dim_MOD, 12/), float)
    goswami_200      = new((/dim_MOD, 12/), float)
    goswami_850      = new((/dim_MOD, 12/), float)
    webster_yang     = new((/dim_MOD, 12/), float)
    webster_yang_200 = new((/dim_MOD, 12/), float)
    webster_yang_850 = new((/dim_MOD, 12/), float)
    wang_fan         = new((/dim_MOD, 12/), float)
    wang_fan_850_1   = new((/dim_MOD, 12/), float)
    wang_fan_850_2   = new((/dim_MOD, 12/), float)

    ;; Add dimensions to monthly climatology arrays
    model_list = NewList("lifo")
    ListPush(model_list, goswami)
    ListPush(model_list, webster_yang)
    ListPush(model_list, wang_fan)
    ListPush(model_list, goswami_200)
    ListPush(model_list, webster_yang_200)
    ListPush(model_list, goswami_850)
    ListPush(model_list, webster_yang_850)
    ListPush(model_list, wang_fan_850_1)
    ListPush(model_list, wang_fan_850_2)

    ;; Model dimensions
    model_tmp = models@name
    model_tmp@long_name = "Model names"
    month_tmp = ispan(0, 11, 1)
    month_tmp@long_name = "Month number"

    do idx = 0, ListCount(model_list) - 1
        model_list[idx]!0 = "model"
        model_list[idx]!1 = "month"
        model_list[idx]&model = model_tmp
        model_list[idx]&month = month_tmp
    end do

    ;; Define regions for the various indicies

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

    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and modify data
            ;;------------------------
            A0_ua_full = extract_data(idx_mod, nml_data_vault_1, -1, 0, 0)
            A0_va_full = extract_data(idx_mod, nml_data_vault_2, -1, 0, 0)
            ua_200 = clmMonTLL(A0_ua_full(time|:, {plev|20000}, lat|:, lon|:))
            ua_850 = clmMonTLL(A0_ua_full(time|:, {plev|85000}, lat|:, lon|:))
            va_200 = clmMonTLL(A0_va_full(time|:, {plev|20000}, lat|:, lon|:))
            va_850 = clmMonTLL(A0_va_full(time|:, {plev|85000}, lat|:, lon|:))
            ua_200&month = ispan(0, 11, 1)
            va_200&month = ispan(0, 11, 1)
            ua_850&month = ispan(0, 11, 1)
            va_850&month = ispan(0, 11, 1)

            ;; ------------------------------------------------------
            ;; First monsoon index (Monsoon Hadley Index)
            ;; ------------------------------------------------------
            ;; Goswami et al. A broad-scale circulation index for the
            ;; interannual variability of the Indian summer monsoon,
            ;; Quarterly Journal of the Royal Meteorological Society,
            ;; vol. 125, issue 554, pp. 611-633
            ;; ------------------------------------------------------
            goswami_200(idx_mod, :) = coswgt_areaave(va_200(month|:,\
                                                            {lat|gos_minlat:gos_maxlat},\
                                                            {lon|gos_minlon:gos_maxlon}))
            goswami_850(idx_mod, :) = coswgt_areaave(va_850(month|:,\
                                                            {lat|gos_minlat:gos_maxlat},\
                                                            {lon|gos_minlon:gos_maxlon}))
            goswami(idx_mod, :) = goswami_850(idx_mod, :) - goswami_200(idx_mod, :)
            goswami@units = goswami_850@units

            ;; ------------------------------------------------------
            ;; Second monsoon index (Webster-Yang monsoon index)
            ;; ------------------------------------------------------
            ;; Webster, P. J. and S. Yang, 1992: Monsoon and ENSO:
            ;; Selectively interactive systems. Quart. J. Roy. Meteor.
            ;; Soc., 118, 877-926
            ;; ------------------------------------------------------
            webster_yang_200(idx_mod, :) = coswgt_areaave(ua_200(month|:,\
                                                                 {lat|wy_minlat:wy_maxlat},\
                                                                 {lon|wy_minlon:wy_maxlon}))
            webster_yang_850(idx_mod, :) = coswgt_areaave(ua_850(month|:,\
                                                                 {lat|wy_minlat:wy_maxlat},\
                                                                 {lon|wy_minlon:wy_maxlon}))
            webster_yang(idx_mod, :) = webster_yang_850(idx_mod, :) - webster_yang_200(idx_mod, :)
            webster_yang@units = webster_yang_850@units

            ;; ------------------------------------------------------
            ;; Third  monsoon index (Dynamical Monsoon index)
            ;; ------------------------------------------------------
            ;; Wang, B., and Z. Fan, 1999: Choice of south Asian summer
            ;; monsoon indices.
            ;; Bull. Amer. Meteor. Soc., 80, 629-638.
            ;; ------------------------------------------------------
            wang_fan_850_1(idx_mod, :)  = coswgt_areaave(ua_850(month|:,\
                                                                {lat|wf1_minlat:wf1_maxlat},\
                                                                {lon|wf1_minlon:wf1_maxlon}))
            wang_fan_850_2(idx_mod, :)  = coswgt_areaave(ua_850(month|:,\
                                                                {lat|wf2_minlat:wf2_maxlat},\
                                                                {lon|wf2_minlon:wf2_maxlon}))
            wang_fan(idx_mod, :) = wang_fan_850_1(idx_mod, :) - wang_fan_850_2(idx_mod, :)
            wang_fan@units = wang_fan_850_2@units

            delete(A0_ua_full)
            delete(A0_va_full)

            delete(ua_200)
            delete(ua_850)
            delete(va_200)
            delete(va_850)

            ;; Check if start/end years are consistent across all models
            if ((first_start_year .ne. models@start_year(idx_mod))\
                .or. (first_end_year .ne. models@end_year(idx_mod))) then
                model_years_consistent = False
            end if
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Store computed arrays in vault
        storage_variable_name = "goswami"
        store_variable(goswami, storage_vault, storage_variable_name)

        storage_variable_name = "goswami-200"
        store_variable(goswami_200, storage_vault, storage_variable_name)

        storage_variable_name = "goswami-850"
        store_variable(goswami_850, storage_vault, storage_variable_name)

        storage_variable_name = "webster-yang"
        store_variable(webster_yang, storage_vault, storage_variable_name)

        storage_variable_name = "webster-yang-200"
        store_variable(webster_yang_200, storage_vault, storage_variable_name)

        storage_variable_name = "webster-yang-850"
        store_variable(webster_yang_850, storage_vault, storage_variable_name)

        storage_variable_name = "wang-fan"
        store_variable(wang_fan, storage_vault, storage_variable_name)

        storage_variable_name = "wang-fan-850-1"
        store_variable(wang_fan_850_1, storage_vault, storage_variable_name)

        storage_variable_name = "wang-fan-850-2"
        store_variable(wang_fan_850_2, storage_vault, storage_variable_name)

        delete(goswami)
        delete(goswami_200)
        delete(goswami_850)
        delete(webster_yang)
        delete(webster_yang_200)
        delete(webster_yang_850)
        delete(wang_fan)
        delete(wang_fan_850_1)
        delete(wang_fan_850_2)

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"goswami",          "goswami-200",\
                             "goswami-850",      "webster-yang",\
                             "webster-yang-200", "webster-yang-850",\
                             "wang-fan",         "wang-fan-850-1",\
                             "wang-fan-850-2"/)
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
    di@years_string = years_string




    ;; ============================================
    ;; ||                                        ||
    ;; ||          PLOTTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================
    plot_settings = True

; --------------------------

    plot_settings@title_string = "Goswami Monsoon Index"
    plot_settings@aux_filename_info = "goswami-seasonal-cycle"
    plot_settings@yaxis_string = "v850-v200"
    plot_settings@type = "monthly_climatology"
    plot_settings@left_plot_subheader = "Seasonal cycle, " + years_string
    debuginfo@description_ycoord = 0.7
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology for the Goswami Monsoon Index (GMI)"\
                                + "~C~   across specified years. GMI is the vertical northward wind shear for"\
                                + "~C~   the area,"\
                                + "~C~      lat: " + gos_minlat + "-" + gos_maxlat\
                                + ", lon: " + gos_minlon + "-" + gos_maxlon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~Diag ref:~F21~ DOI: 10.1002/qj.49712555412"
    xy_plot_wrapper(storage_vault, di, plot_settings, "goswami", debuginfo)

; --------------------------

    plot_settings@title_string = "Goswami Monsoon Wind at v200"
    plot_settings@aux_filename_info = "goswami-seasonal-cycle-v200"
    plot_settings@yaxis_string = "v200"
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology for the northward wind speed. This"\
                                + "~C~   is one of the components of the Goswami Monsoon Index (GMI)"\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~GMI ref:~F21~ DOI: 10.1002/qj.49712555412"
    xy_plot_wrapper(storage_vault, di, plot_settings, "goswami-200", debuginfo)

; --------------------------

    plot_settings@title_string = "Goswami Monsoon Wind at v850"
    plot_settings@aux_filename_info = "goswami-seasonal-cycle-v850"
    plot_settings@yaxis_string = "v850"
    xy_plot_wrapper(storage_vault, di, plot_settings, "goswami-850", debuginfo)

; --------------------------

    plot_settings@title_string = "Webster Yang Monsoon Index"
    plot_settings@aux_filename_info = "webster-yang-seasonal-cycle"
    plot_settings@yaxis_string = "u850-u200"
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology for the Webster Yang Monsoon Index"\
                                + "~C~   (WYI) across specified years. WYI is the eastward vertical wind"\
                                + "~C~   shear average for the area,"\
                                + "~C~      lat: " + wy_minlat + "-" + wy_maxlat\
                                + ", lon: " + wy_minlon + "-" + wy_maxlon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~Diag ref:~F21~ DOI: 10.1002/qj.49711850705"
    xy_plot_wrapper(storage_vault, di, plot_settings, "webster-yang", debuginfo)

; --------------------------

    plot_settings@title_string = "Webster Yang Monsoon Wind at u200"
    plot_settings@aux_filename_info = "webster-yang-seasonal-cycle-u200"
    plot_settings@yaxis_string = "u200"
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology for the eastward wind speed"\
                                + "~C~   average. This is one of the components of the Webster Yang "\
                                + "~C~   Monsoon Index (WYI)."\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~WYI ref:~F21~ DOI: 10.1002/qj.49711850705"
    xy_plot_wrapper(storage_vault, di, plot_settings, "webster-yang-200", debuginfo)

; --------------------------

    plot_settings@title_string = "Webster Yang Monsoon Wind at u850"
    plot_settings@aux_filename_info = "webster-yang-seasonal-cycle-u850"
    plot_settings@yaxis_string = "u850"
    xy_plot_wrapper(storage_vault, di, plot_settings, "webster-yang-850", debuginfo)

; --------------------------

    plot_settings@title_string = "Wang Fan Monsoon Index"
    plot_settings@aux_filename_info = "wang-fan-seasonal-cycle"
    plot_settings@yaxis_string = "u850~B~area_1~E~        - u850~B~area_2~E~"
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology for the Wang Fan Monsoon Index (WFI)"\
                                + "~C~   WYI is the eastward wind speed difference between the areas,"\
                                + "~C~       lat: " + wf1_minlat + "-" + wf1_maxlat\
                                + ", lon: " + wf1_minlon + "-" + wf1_maxlon + " and"\
                                + "~C~       lat: " + wf2_minlat + "-" + wf2_maxlat\
                                + ", lon: " + wf2_minlon + "-" + wf2_maxlon\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~Diag ref:~F21~ DOI: 10.1175/1520-0477(1999)080<0629:COSASM>2.0.CO;2"
    xy_plot_wrapper(storage_vault, di, plot_settings, "wang-fan", debuginfo)

; --------------------------

    plot_settings@title_string = "Wang Fan Monsoon Wind at Area 1"
    plot_settings@aux_filename_info = "wang-fan-seasonal-cycle-area1"
    plot_settings@yaxis_string = "u850~B~area_1~E~"
    debuginfo@description = "~F22~Diag def:~F21~ Monthly climatology for the eastward wind speed average."\
                                + "~C~   This is one of the components of the Wang Fan Monsoon Index (WFI)."\
                                + "~C~ The gray area indicates standard deviation from model mean to"\
                                + "~C~ indicate the spread between models (obs/reanalysis are excluded)"\
                                + "~C~~F22~WFI ref:~F21~ DOI: 10.1175/1520-0477(1999)080<0629:COSASM>2.0.CO;2"
    xy_plot_wrapper(storage_vault, di, plot_settings, "wang-fan-850-1", debuginfo)

; --------------------------

    plot_settings@title_string = "Wang Fan Monsoon Wind at Area 2"
    plot_settings@aux_filename_info = "wang-fan-seasonal-cycle-area2"
    plot_settings@yaxis_string = "u850~B~area_2~E~"
    xy_plot_wrapper(storage_vault, di, plot_settings, "wang-fan-850-2", debuginfo)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
