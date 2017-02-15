;;#############################################################################
;; WAMonsoon_isv_filtered.ncl
;; Author: Romain Roehrig (Meteo France, France)
;; EMBRACE project
;;#############################################################################
;; Description
;;    Filters a daily variable of contours its standard deviation
;;    over West Africa
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    cn_levels_stddev_basic:            Contour levels for stddev plot
;;    cn_levels_stddev_diff_basic:       Contour levels for diff stddev plot
;;    cn_levels_stddev_norm_basic:       Contour levels for mean plot
;;    cn_levels_stddev_norm_diff_basic:  Contour levels for diff plot
;;    latrange:                          Latitude range for plot
;;    lonrange:                          Longitude range for plot (-10, 10)
;;    my_region:                         Switch to indicate Monsoon region
;;    sahel_latrange:                    Latitude crop range for Sahel region
;;    sahel_lonrange:                    Longitude crop range for Sahel region
;;    season:                            Season, typically "JJAS"
;;    destgrid:                          Destination grid for ESMF regridding
;;    styleset:                          Stylesheet to use
;;    filter_hp:                         Filter settings
;;    filter_max:                        Filter settings
;;    filter_min:                        Filter settings
;;    filter_type:                       Filter settings
;;    filter_weights:                    Filter settings
;;    plot_norm:                         Also plot normalized stddev
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;    long_name:                         Name displayed in plot
;;    skip_long_name_in_subheader:       Logical to switch wether 'long_name' is displayed or not
;;    units:                             Displayed units
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150515-A_roeh_ro: written
;;
;; #############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/data_handling.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
;load "./diag_scripts/aux/WAMonsoon/descriptions.ncl"

exact_panel_position_file = "./diag_scripts/aux/" + diag_script_info@my_region + "Monsoon/exact_panel_positions_pr-mmday.ncl"

load "./plot_scripts/ncl/monsoon_panels.ncl"
;*********************************************************************
undef("high_pass_latlon_time")
function high_pass_latlon_time(x[*][*][*]:numeric,
                               spd[1]:numeric,
                               hpf[2]:numeric,
                               opt:logical)
local nMsg, dimx, ntim, nlat, mlon, dNam, n, tim_taper,\
    hpfc, hpfNwgt, fca, fcb, ihp, sigma,
    hpfWgt, work, WORK, cf, WCF, error_string
begin
    nMsg = num(ismissing(x))
    error_string = "band_pass_latlon_time: currently, missing data not allowed: nMsg=" + nMsg
    if (nMsg .gt. 0) then
        error_msg("fatal",\
                  "WAMonsoon_isv_filtered.ncl",\
                  "high_pass_latlon_time",\
                  error_string)
        status_exit(1)
    end if

    dimx = dimsizes(x)
    ntim = dimx(0)
    nlat = dimx(1)
    mlon = dimx(2)

    hpfc = hpf(0)                     ; days
    hpfNwgt = hpf(1)                     ; effective # weights

    hpfNwgt = (hpfNwgt / 2) * 2 + 1          ; make sure it is odd

    fc = 1.0 / (spd * hpfc)

    dNam = getvardims(x)                 ; get dimension names
    do n = 0, 2                             ; only used if detrending
        if (ismissing(dNam(n))) then      ; or tapering in time
            x!n = "dim" + n                ; assign name
            dNam =  x!n
        end if
    end do

    ihp = 1                               ; hpf=>high pass filter
    sigma = 1.0                             ; Lanczos sigma
    hpfwgt = filwgts_lanczos(hpfNwgt, ihp, fc, -999, sigma)

    if (opt .and. isatt(opt, "fft") .and. opt@fft) then

; As of 6.2.0 there are no ezfft_n / exfftb_n
; until there is an ezfft_n (exfftb_n) the opt@fft=True option
; must be handled in the 'old-fashioned' way' (dimension reordering)

    ;;work = x(lat|:, lon|:, time|:)       ; reorder so time is rightmost
        work = x($dNam(1)$|:, $dNam(2)$|:, $dNam(0)$|:)   ; generic

        if (isatt(work, "_FillValue")) then
            delete(work@_FillValue)          ; avoid annoying warning msg
        end if                               ; when calling ezfftf

                                             ; By default ... no detrend
        if (opt .and.(isatt(opt, "detrend") .and. opt@detrend)) then
            work = dtrend(work, False)         ; detrend reordered array
            work@detrend = "data detrended in time"
        end if
                                             ; By default ... no taper
        if (isatt(opt, "taper")) then
            tim_taper = opt@taper
        else
            tim_taper = 0.10                 ; default taper is 10%
        end if

        work = taper(work, tim_taper, 0)  ; rightmost dimension (time)
        work@taper = "data tapered in time"
                                             ; fft in time
        cf = ezfftf(work)                 ; cf(2, nlat, mlon, ntim/2)
                                             ; map response of digitial
                                             ; filter to fft space
        fcf = fspan(0, 0.5, ntim / 2)         ; fft freq
        wcf = linint1(hpfwgt@freq, hpfwgt@resp, False, fcf, 0)
        WCF = conform(cf(0, :, :, :), wcf, 2)
        delete(wcf)

        cf(0, :, :, :) = cf(0, :, :, :) * WCF        ; apply response coef
        cf(1, :, :, :) = cf(1, :, :, :) * WCF
        delete(WCF)

        work = ezfftb(cf, 0.0)               ; fourier synthesis
        delete(cf)
        work@process =  "FFT with digital respoonse mapped tp FFT space"
    else
        work = x                          ; _n
                                            ; By default ... no detrend
        work = dtrend_n(work, False, 0)   ; detrend
        work@detrend = "data detrended in time"

        work = wgt_runave_n_Wrap(work, hpfwgt, 0, 0)
        work@process =  "wgt_runave"
    end if

    work@high_pass_fc = hpfc
    work@high_pass_Nwgts = hpfNwgt

;;X = work($dNam(0)$|:, $dNam(1)$|:, $dNam(2)$|:)
;;copy_VarMeta (x, X)
;;X@long_name = "High Pass: "
;;if (isatt(x, "long_name")) then
;;    X@long_name = "High Pass: "+x@long_name
;;end if
;;return( X )

    if (isatt(work, "long_name")) then
        work@long_name = "High Pass: " + work@long_name
    else
        work@long_name = "High Pass: "
    end if

    return(work)
end
;*********************************************************************
begin
    ;; ============================================
    ;; ||                                        ||
    ;; ||         READING CONFIGURATION          ||
    ;; ||                                        ||
    ;; ============================================

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
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
                     "A_roeh_ro",            \  ; authors
                     "",                     \  ; contributors
                     "",                     \  ; diag_references
                     "",                     \  ; obs_references
                     "P_embrace")               ; proj_acknowledgements
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
    ;; (from the WAMonsoon cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"season",\
                            "latrange_basic",\
                            "lonrange_basic",\
                            "cn_levels_stddev_norm_basic",\
                            "cn_levels_stddev_norm_diff_basic",\
                            "filter_type"/)
    exit_if_missing_atts(diag_script_info, required_attributes)
    delete(required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    if (di@filter_type .eq. "hp") then
        required_attributes = (/"filter_hp",\
                                "filter_weights"/)
    else if (di@filter_type .eq. "bp") then
        required_attributes = (/"filter_min",\
                                "filter_max",\
                                "filter_weights"/)
    else
        print("fatal:filtering unknown:" + di@filter_type)
        status_exit(1)
    end if
    end if

    do idx_req = 0, dimsizes(required_attributes) - 1
        if (.not. isatt(di, required_attributes(idx_req))) then
            print("fatal:plot requires diag_script_info attribute '" + required_attributes(idx_req) + "'")
            status_exit(1)
        end if
    end do

    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"interp",\
                               "latrange",\
                               "lonrange",\
                               "fill_missing_w_poisson",\
                               "lonflip"/)

    transforms@season = di@season
    transforms@latrange = di@latrange_basic
    transforms@lonrange = di@lonrange_basic
    transforms@is_cyclic = False

    transforms@destgrid = di@destgrid  ; "1x1"

    debuginfo@latrange = di@latrange_basic(0) + "-" + di@latrange_basic(1)
    debuginfo@lonrange = di@lonrange_basic(0) + "-" + di@lonrange_basic(1)
    debuginfo@variable = var0
    if (di@filter_type .eq. "hp") then
        debuginfo@filter_hp = di@filter_hp
        debuginfo@filter_no_weights = di@filter_weights
    end if
    if (di@filter_type .eq. "bp") then
        debuginfo@filter_min = di@filter_min
        debuginfo@filter_max = di@filter_max
        debuginfo@filter_no_weights = di@filter_weights
    end if

    min_lat = di@latrange_basic(0)
    max_lat = di@latrange_basic(1)
    min_lon = di@lonrange_basic(0)
    max_lon = di@lonrange_basic(1)

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

            data_jjas_tmp = extract_season(A0, di@season)
            data_jjas = dim_avg_Wrap(data_jjas_tmp(lat|:, lon|:, time|:))
            delete(data_jjas_tmp)

            opt = True
            opt@detrend = True
            srate = 1
            if (di@filter_type .eq. "hp") then
                hpf = (/di@filter_hp, di@filter_weights/)
                filtered = high_pass_latlon_time(A0, srate, hpf, opt)
                filter = "1-" + di@filter_hp + "-day"
                di@filter_name = "hp" + di@filter_hp
            end if
            if (di@filter_type .eq. "bp") then
                bpf = (/di@filter_min, di@filter_max, di@filter_weights/)
                filtered = band_pass_latlon_time(A0, srate, bpf, opt)
                filter = di@filter_min + "-" + di@filter_max + "-day"
                di@filter_name = "bp" + di@filter_min + "-" + di@filter_max
            end if

            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------

            season_filtered = extract_season(filtered, di@season)
            var_stddev = dim_stddev_n_Wrap(season_filtered, 0)

            ;; Compute reference statistical values over Sahel (sahel)
            stat_sahel_const = coswgt_areaave(var_stddev({lat|di@sahel_latrange(0):di@sahel_latrange(1)},\
                                                         {lon|di@sahel_lonrange(0):di@sahel_lonrange(1)}))

            clim_mean = coswgt_areaave(data_jjas({lat|di@sahel_latrange(0):di@sahel_latrange(1)},\
                                                 {lon|di@sahel_lonrange(0):di@sahel_lonrange(1)}))
            stat_sahel_normalized = stat_sahel_const / clim_mean

            var_stddev@sahel = stat_sahel_const
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
            data_jjas_with_fillvalues = where(data_jjas .lt. 1., var_FillValue, data_jjas)
            var_stddev_norm = var_stddev / data_jjas_with_fillvalues

            var_stddev_norm@sahel = stat_sahel_normalized
            var_stddev_norm@idx = idx_mod
            var_stddev_norm@GridType = "rectilinear"
            var_stddev_norm@Regional = True
            var_stddev_norm@model_name = models@name(idx_mod)
            storage_variable_name = str_join((/"stddev_norm", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var_stddev_norm, storage_vault, storage_variable_name)

            delete(A0)
            delete(data_jjas)
            delete(data_jjas_with_fillvalues)
            delete(filtered)
            delete(season_filtered)
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
    plot_settings@part_of_header = variable_info@long_name
    plot_settings@type = "Not applicable"
    debuginfo@description_ycoord = 0.8
    ;debuginfo@description = WAMonsoon_daily_precip_filter_descr(di)

    plot_settings@gsn_right_string = variable_info@long_name
    plot_settings@main_title_suffix = " (" + filter + " " + variable_info@long_name + ")"
    valid_statistics = (/"yrs", "mean", "rmse", "sahel"/)

    ;; Colormaps
    if (isatt(di, "diff_colormap")) then
        diff_colormap = read_colormap_file(di@diff_colormap)
        plot_settings@diff_palette = diff_colormap
    end if

    ;; Three panel plots (ref, model, diff), here for stddev
    three_by_one_diff(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    if (di@plot_norm) then
        ;; Three panel plots (ref, model, diff), here for normalized and masked stddev
        plot_settings@main_title_suffix = " (" + filter + " " + variable_info@long_name + " normalized)"

        plot_settings@lb_units = "unitless"
        ;debuginfo@description = WAMonsoon_daily_precip_filter_norm_descr(di)
        three_by_one_diff(storage_vault, di, plot_settings, "stddev_norm", debuginfo, valid_statistics, res)
    end if

    res = orig_res
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    plot_settings@type = "mean"

    plot_settings@lb_units = lb_units
    plot_settings@part_of_header = " stddev of " + filter + " " + variable_info@long_name
    ;debuginfo@description = WAMonsoon_daily_precip_filter_descr(di)
    multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    if (di@plot_norm) then
        plot_settings@lb_units = "unitless (normalized with climatology)"
        plot_settings@part_of_header = " stddev of " + filter + " " + variable_info@long_name + " (climatology mask+norm)"
        ;debuginfo@description = WAMonsoon_daily_precip_filter_norm_descr(di)
        multipanel(storage_vault, di, plot_settings, "stddev_norm", debuginfo, valid_statistics, res)
    end if

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
