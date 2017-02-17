;;#############################################################################
;; DiurnalCycle_precip_harmonic.ncl
;; Author(s): Rio Catherine (LMD, France)
;;
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    Season Maps of the amplitude and phase of the first harmonic of the mean
;;    diurnal cycle of precipitation
;;
;; Required diag_script_info attributes (diagnostics specific)
;;                latrange_basic:latitude range of box
;;                lonrange_basic:longitude range of box
;;                        season: season
;;     cn_levels_amplitude_basic: contor levels for amplitude plots
;;        cn_levels_phase_basics: contor levels for phase plots
;;           cn_levels_var_basic: contor levels for pourcentage of variance
;;                                plots
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;
;; Caveats
;;   The settings for GridType and Regional are set explicitly in this script
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;#############################################################################

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip_world.ncl"

load "./plot_scripts/ncl/monsoon_panels.ncl"
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
                     "A_rio_ca",             \  ; authors
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
    valid_plot_field = (/"T2Is"/)
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
                            "cn_levels_amplitude_basic",\
                            "cn_levels_phase_basic",\
                            "cn_levels_var_basic" /)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"interp", "lonflip", "season"/)
    transforms@season = di@season
    transforms@latrange = di@latrange_basic
    transforms@lonrange = di@lonrange_basic
    transforms@destgrid = di@destgrid  ; "2.5x2.5"

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

    ;; ----------------------------------------------
    ;; Data structure for storing the processed field
    ;; ----------------------------------------------
    storage_vault = True

    ;; Diagnostic data arrays
    ;; using CMIP5 3hr data
    nHR = 8
    model_tmp = models@name
    model_tmp@long_name = "Model names"
    hour_tmp = new((/nHR/), float)
    min_tmp = new((/nHR/), float)
    hour_tmp@long_name = "hour of day"
    min_tmp@long_name = "minute of day"
    ;; Interpolated grid dimension: valid for 2.5x2.5 grid
    mlon = toint(360/2.5)  ; number of output longitudes
    mlat = toint(180/2.5)  ; number of output latitudes

    klon0 = ispan(-mlon/2 + 1, mlon/2 - 1, 1) * 2.5
    klat0 = ispan(-mlat/2 + 1, mlat/2 - 1, 1) * 2.5

    klat = tofloat(klat0)
    klon = tofloat(klon0)

    klat@units = "degrees_north"
    klon@units = "degrees_east"

    ;; Mean diurnal cycle
    dc = new((/dimsizes(klat), dimsizes(klon), nHR/), "float")
    dc!0 = "lat"
    dc&lat = klat
    dc!1 = "lon"
    dc&lon = klon

    utc = new((/nHR + 1/), float)

    finfo_dc = new((/3, dimsizes(klat), dimsizes(klon), 1/), "float")
    finfo_dc!1 = "lat"
    finfo_dc&lat = klat
    finfo_dc!2 = "lon"
    finfo_dc&lon = klon

    ;; ----------------------------------------------
    ;; Diagnostic loop, compute diagnostic and store
    ;;                  in 'storage_vault'
    ;; ----------------------------------------------
    if (.not. read_from_vault) then

        do idx_mod = 0, dim_MOD - 1

            ;; Extract and modify data
            A0_full = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
            A0 = apply_transforms(transforms, A0_full)
            delete(A0_full)

            ;; Compute mean for each hour
            timeh = cd_calendar(A0&time, 0)
            hour_tmp = timeh(0:nHR - 1, 3)
            min_tmp = timeh(0:nHR - 1, 4)

            n = -1
            do nh = 0, nHR - 1
                n = n + 1

                ;; rFfind indices for each hour
                indh = ind(timeh(:, 3) .eq. hour_tmp(nh))

                ;; Compute time
                utc(n) = hour_tmp(nh) + min_tmp(nh) / 60.

                ;; Compute mean diurnal cycle
                dc(:, :, n) = dim_avg_n(A0(indh, :, :), 0)
                delete([/indh/])
            end do
            utc(nHR) = utc(nHR - 1) + utc(nHR - 1) - utc(nHR - 2)

            dc@units = A0@units
            dc@_FillValue = A0@_FillValue

            ;;****************************************
            ;; Calculate temporal harmonics: reorder so 'time' is right dimension
            ;;   finfo(0,:, :, nhar) - amplitude of each harmonic
            ;;   finfo(1,:, :, nhar) - phase of each harmonic
            ;;   finfo(2,:, :, nhar) - percent variance explained by each harmonic
            ;;****************************************
            nhar = 1  ; Return 1st 'nhar' harmonics

            ;; Check for missing values
            do i = 0, dimsizes(klat) - 1, 1
                do j = 0, dimsizes(klon) - 1, 1
                    if(any(ismissing(dc(i, j, :)))) then
                        info_output("Your data dc contains some missing values. Beware.", verbosity, 1)
                        finfo_dc(:, i, j, :) = dc@_FillValue
                    else
                        ;; Fourier_info works on right dim
                        finfo_dc(:, i, j, :) = fourier_info(dc(i, j, :), nhar, 1.)
                    end if
                end do
            end do
            copy_VarCoords(dc(:, :, 0), finfo_dc(0, :, :, 0))

            ;; Store amplitude
            amplitude = finfo_dc(0, :, :, 0)
            amplitude@_FillValue = A0@_FillValue
            ;; Store phase
            phase = finfo_dc(1, :, :, 0)
            phase@_FillValue = A0@_FillValue
            ;; Compute phase in utc and local hour
            do i = 0, dimsizes(klat) - 1, 1
                do j = 0, dimsizes(klon) - 1, 1
                    if (.not. ismissing(phase(i, j))) then
                        utc_ceil = utc(floattointeger(ceil(phase(i, j))))
                        utc_floor = utc(floattointeger(floor(phase(i, j))))
                        utc_diff = utc_ceil - utc_floor
                        phase_frac = phase(i, j) - floattointeger(phase(i, j))

                        phase(i, j) = utc(floattointeger(phase(i, j))) + ((phase_frac) * (utc_diff))
                        phase(i, j) = phase(i, j) + klon(j) * 24. / 360.
                        if (phase(i, j) .ge. 24.) then
                            phase(i, j) = phase(i, j) - 24.
                        end if
                        if (phase(i, j) .lt. 0.) then
                            phase(i, j) = phase(i, j) + 24.
                        end if
                    else
                        info_output("Missing values", verbosity, 1)
                    end if
                end do
            end do

            ;; Store pourcentage of variance
            var = finfo_dc(2, :, :, 0)
            var@_FillValue = A0@_FillValue

            delete(A0)
            delete(timeh)

            amplitude@idx = idx_mod
            amplitude@model_name = models@name(idx_mod)
            amplitude@GridType = "rectilinear"
            amplitude@Regional = True
            storage_variable_name = str_join((/"amplitude", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(amplitude, storage_vault, storage_variable_name)

            phase@idx = idx_mod
            phase@model_name = models@name(idx_mod)
            phase@GridType = "rectilinear"
            phase@Regional = True
            storage_variable_name = str_join((/"phase", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(phase, storage_vault, storage_variable_name)

            var@idx = idx_mod
            var@model_name = models@name(idx_mod)
            var@GridType = "rectilinear"
            var@Regional = True
            storage_variable_name = str_join((/"var", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var, storage_vault, storage_variable_name)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"amplitude", "phase", "var"/)
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

    plot_settings = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_title = lb_units
    plot_settings@lb_units = lb_units
    plot_settings@default_palette = "BlAqGrYeOrReVi200"

    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ " + di@season + " average across indicated years" +\
                            "~C~~F22~rmse def:~F21~ NCL function wgt_arearmse(...) computed with " +\
                            "~C~                  cosine weights in latitude." +\
                            "~C~~F22~corr def:~F21~ NCL function pattern_cor(...) computed with " +\
                            "~C~                  cosine weights in latitude. Correlation is between" +\
                            "~C~                  reference and model data."
    valid_statistics = (/"mean", "rmse", "corr"/)
    res = orig_res
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    plot_settings@type = "mean"
    plot_settings@part_of_header = " of the 1st harmonic of the mean diurnal cycle of " + variable_info@long_name
    plot_settings@prefix_to_header = ""
    multipanel(storage_vault, di, plot_settings, "amplitude", debuginfo, valid_statistics, res)

    plot_settings@type = "mean"
    plot_settings@part_of_header = " of the 1st harmonic of the mean diurnal cycle of " + variable_info@long_name
    plot_settings@lb_units = "hour"
    multipanel(storage_vault, di, plot_settings, "phase", debuginfo, valid_statistics, res)

    plot_settings@type = "mean"
    plot_settings@prefix_to_header = " "
    plot_settings@part_of_header = " of the 1st harmonic of the mean diurnal cycle of " + variable_info@long_name
    plot_settings@lb_units = "%"
    multipanel(storage_vault, di, plot_settings, "var", debuginfo, valid_statistics, res)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
