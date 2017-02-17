;;#############################################################################
;; SAMonsoon_precip_propagation.ncl
;; Author(s): Ricard Levine (UKMET, Great Britain)
;;            Shiyu Wang (SMHI, Sweden)
;;            Martin Evaldsson (SMHI, Sweden)
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    Lagged correlation of precipitation anomalies with respect to an area
;;    avearged region in the Indian Ocean
;;
;; Required diag_script_info attributes (diagnostics specific)
;;                  fca: start frequency to compute Lancoz weights
;;                  fcb: end frequency to compute Lancoz weights
;;                  ihp: lancoz switch for low/high/band pass filter
;;                sigma: lancoz sigma
;;                 nWgt: total number of weights for lancoz
;;             min_days: min days in filtering
;;             max_days: max days in filtering
;;    prop_filter_mxlag: size of lag
;;    prop_isv_latrange: latrange crop area
;;    prop_isv_lonrange: lonrange crop area
;;    prop_lag_latrange: latrange crop area for lag computation along longitude
;;    prop_lag_lonrange: lonrange crop area for lag computation along latitude
;;    prop_ref_latrange: latrange crop area for lag reference area
;;    prop_ref_lonrange: lonrange crop area for lag reference area
;;               season: season in letters, e.g., JJAS
;; cn_levels_mean_basic: contour levels for plot
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

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_propagation.ncl"
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
    diag_script_base = basename(diag_script)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

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
    required_attributes = (/"fca",\
                            "fcb",\
                            "ihp",\
                            "min_days",\
                            "max_days",\
                            "nWgt",\
                            "prop_filter_mxlag",\
                            "prop_isv_latrange",\
                            "prop_isv_lonrange",\
                            "prop_lag_latrange",\
                            "prop_lag_lonrange",\
                            "prop_ref_latrange",\
                            "prop_ref_lonrange",\
                            "season",\
                            "sigma",\
                            "cn_levels_mean_basic"/)
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
    transforms@latrange = di@prop_isv_latrange
    transforms@lonrange = di@prop_isv_lonrange
    transforms@diag_script_base = diag_script_base
    transforms@variable = var0
    transforms@field_type = field_type0
    transforms@is_cyclic = False

    debuginfo@latrange = di@prop_isv_latrange(0) + "-" + di@prop_isv_latrange(1)
    debuginfo@lonrange = di@prop_isv_lonrange(0) + "-" + di@prop_isv_lonrange(1)
    debuginfo@variable = var0
    debuginfo@filter_min = di@max_days
    debuginfo@filter_max = di@min_days
    debuginfo@filter_no_weights = di@nWgt

    debuginfo@ref_lat = write_tuple_with_dash(di@prop_ref_latrange)
    debuginfo@ref_lon = write_tuple_with_dash(di@prop_ref_lonrange)
    debuginfo@lag_lon = write_tuple_with_dash(di@prop_lag_latrange)
    debuginfo@lag_lat = write_tuple_with_dash(di@prop_lag_lonrange)

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

    mxlag = di@prop_filter_mxlag





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

    ;; -----------------------------------------------
    ;; Resources for contour plot
    ;; -----------------------------------------------
    res = True

    ;; Compute lancoz weights
    ihp   = di@ihp     ; bpf=>band pass filter
    nWgt  = di@nWgt
    sigma = di@sigma   ; Lanczos sigma
    fca   = di@fca     ; start freq
    fcb   = di@fcb     ; end freq
    wgt   = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

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
            transforms@idx = idx_mod
            A0 = apply_transforms(transforms, A0_tmp)
            delete(A0_tmp)

            ;; Extract time dim in format yyyyddd to use when
            ;; computing daily anomaly
            ymd  = cd_calendar(A0&time, -2)   ; entire file
            if (isatt(ymd, "calendar")) then
                if (ymd@calendar .eq. "proleptic_gregorian") then
                    ymd@calendar = "gregorian"
                end if
            end if
            yyyyddd = yyyymmdd_to_yyyyddd(ymd)
            prClmDay = clmDayTLL_calendar(A0, yyyyddd)
            no_of_harmonics = 3
            prClmDay_sm = smthClmDayTLL_calendar(prClmDay, no_of_harmonics)
            pr_sm = calcDayAnomTLL_calendar(A0, yyyyddd, prClmDay_sm)

            ;; IO in 'wyIO'/pr_smIO refers to the reference area 'Indian Ocean'
            lat = A0&lat
            wyIO = lat(ind(lat .ge. di@prop_ref_latrange(0) .and.\
                           lat .le. di@prop_ref_latrange(1)))
            wyIO = cos(0.017459 * wyIO)         ; spatial wgt

            pr_smIO    = wgt_areaave_Wrap(pr_sm(time|:,\
                                                {lat|di@prop_ref_latrange(0):di@prop_ref_latrange(1)},\
                                                {lon|di@prop_ref_lonrange(0):di@prop_ref_lonrange(1)}),\
                                                wyIO, 1., 0)  ; (time)
            pr_smIO = dtrend(pr_smIO, False)  ; rmv overall trend
            pr_smIO = wgt_runave_leftdim(pr_smIO, wgt, 0)  ; apply filter

            ;; (time, lon)
            P_timeLon = dim_avg_Wrap(pr_sm(time|:,\
                                           lon|:,\
                                           {lat|di@prop_lag_latrange(0):di@prop_lag_latrange(1)}))
            P_timeLon = dtrend_leftdim(P_timeLon, False)  ; rmv overall trend
            P_timeLon = wgt_runave_leftdim(P_timeLon, wgt, 0)  ; apply filter

            ;; (time, lat)
            P_timeLat = dim_avg_Wrap(pr_sm(time|:,\
                                           lat|:,\
                                           {lon|di@prop_lag_lonrange(0):di@prop_lag_lonrange(1)}))
            P_timeLat = dtrend_leftdim(P_timeLat, False)  ; rmv overall trend
            P_timeLat = wgt_runave_leftdim(P_timeLat, wgt, 0)  ; apply filter

            optXcor = False
            dates = cd_calendar(pr_sm&time, -2)
            ;; JJAS is extracted inside the routines below
            rp_timelon = time_xcor_lag(pr_smIO, P_timeLon, dates, mxlag, optXcor)
            rp_timelat = time_xcor_lag(pr_smIO, P_timeLat, dates, mxlag, optXcor)

            rp_timelon@idx = idx_mod
            rp_timelon@model_name = models@name(idx_mod)

            ; JJAS filtered precip anomaly
            rp_timelon@gsn_right_string = "lag-lon corr"
            rp_timelon@refA = "lat:" + debuginfo@ref_lat + ", lon:" + debuginfo@ref_lon
            rp_timelon@latB = debuginfo@latrange
            storage_variable_name = str_join((/"rp", "timelon", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(rp_timelon, storage_vault, storage_variable_name)

            rp_timelat@idx = idx_mod
            rp_timelat@model_name = models@name(idx_mod)
            rp_timelat@gsn_right_string = "lag-lat corr"
            rp_timelat@refA = "lat:" + debuginfo@ref_lat + ", lon:" + debuginfo@ref_lon
            rp_timelat@lonB = debuginfo@lonrange
            storage_variable_name = str_join((/"rp", "timelat", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(rp_timelat, storage_vault, storage_variable_name)

            delete(A0)
            delete(rp_timelat)
            delete(rp_timelon)
            delete(pr_smIO)
            delete(ymd)
            delete(P_timeLon)
            delete(P_timeLat)
            delete(wyIO)
            delete(dates)
            delete(prClmDay)
            delete(yyyyddd)
            delete(prClmDay_sm)
            delete(lat)
            delete(pr_sm)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"rp-timelon", "rp-timelat"/)
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
    res@tmXBLabelFontHeightF = 0.012
    res@tmYLLabelFontHeightF = 0.012
    res@cnInfoLabelOn = False
    orig_res = res

    ;; Three panel plot (reference, model, and diff)
    plot_settings = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_units = "correlation"
    plot_settings@lb_title = "correlation"
    plot_settings@part_of_header = ", JJAS filtered precip anomaly"
    plot_settings@type = "Not applicable"
    plot_settings@xcor_lagged_plot = True
    plot_settings@gsn_right_string = False

    ;; Update debug info
    delete_VarAtts(debuginfo, (/"latrange", "lonrange"/))
    debuginfo@day_filter_range = di@min_days + "-" + di@max_days
    debuginfo@description_ycoord = 0.8
    debuginfo@description = SAMonsoon_daily_precip_propagation_descr(di)

    plot_settings@main_title_suffix = " (filtered daily propagation)"
    valid_statistics = (/"yrs", "refA", "latB", "lonB"/)

    ;; Three panel plots (ref, model, diff), here for mean
    two_by_one(storage_vault, di, plot_settings, "rp-timelat", "rp-timelon", debuginfo, valid_statistics, res)

    res = orig_res

    ;; Units label,
    lb_units = "~F22~correlation~F21~ for, reference area lat: "\
        + debuginfo@ref_lat\
        + ", lon: " + debuginfo@ref_lon

    ;; Multipanel plots
    plot_settings@type = "mean"

    delete(valid_statistics)
    valid_statistics = (/"yrs"/)
    plot_settings@lb_units = lb_units + ", longitude extent: " + write_tuple_with_dash(di@prop_isv_lonrange)
    plot_settings@gsn_right_string = False
    plot_settings@skip_default_title = False
    plot_settings@part_of_header = " lag-lat propagation of filtered daily Precip anom"
    debuginfo@description = SAMonsoon_daily_precip_propagation_descr(di)
    multipanel(storage_vault, di, plot_settings, "rp-timelat", debuginfo, valid_statistics, res)

    delete(valid_statistics)
    valid_statistics = (/"yrs"/)
    plot_settings@lb_units = lb_units + ", latitude extent: " + write_tuple_with_dash(di@prop_isv_latrange)
    plot_settings@gsn_right_string = False
    plot_settings@skip_default_title = False
    plot_settings@part_of_header = " lag-lon propagation of filtered daily Precip anom"
    debuginfo@description = SAMonsoon_daily_precip_propagation_descr(di)
    multipanel(storage_vault, di, plot_settings, "rp-timelon", debuginfo, valid_statistics, res)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
