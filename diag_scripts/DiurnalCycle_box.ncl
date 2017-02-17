;;#############################################################################
;; DiurnalCycle_box.ncl
;; Author(s): Rio Catherine (LMD, France)
;;
;; PROJECT-NAME: EMBRACE
;;#############################################################################
;; Description
;;    This script computes and outputs averaged diurnal cycle
;; over specific boxes and seasons of precipitation and radiative fluxes
;;
;; Required diag_script_info attributes (diagnostics specific)
;;   latrange_seasonal: latitude crop area
;;   lonrange_seasonal: lonitude crop area
;;              season: season in letters, e.g., JJAS
;;            styleset: chose styles set for lines
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

    ;;-----------------------
    ; References Information
    ;;-----------------------
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
    ;; (from the DiurnalCycle cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"season",\
                            "styleset",\
                            "latrange",\
                            "lonrange"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

     ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"lonflip", "season", "latrange", "lonrange"/)
    transforms@season = di@season
    transforms@latrange = di@latrange
    transforms@lonrange = di@lonrange

    ;;------------------------------------------
    ;; Define meta data to debuginfo onto the plot
    ;;------------------------------------------
    min_lat = di@latrange(0)
    max_lat = di@latrange(1)
    debuginfo@latrange = min_lat + "-" + max_lat
    min_lon = di@lonrange(0)
    max_lon = di@lonrange(1)
    debuginfo@lonrange = min_lon + "-" + max_lon
    latlon_header_string = get_latlon_string(min_lat, max_lat, min_lon, max_lon)

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

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

    ;; Loop through models and compute diagnostics
    if (.not. read_from_vault) then
        do idx_mod = 0, dim_MOD - 1

            ;; Extract and modify data
            A0_full = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
            A0 = apply_transforms(transforms, A0_full)
            delete(A0_full)

            ;; Diagnostic data arrays
            ;; based on CMIP5 3hr data
            ;; If some models do not have 3hr data, change code to
            ;; extract nHR from input file
            nHR = 8
            model_tmp = models@name
            model_tmp@long_name = "Model names"
            hour_tmp = new((/nHR/), float)
            min_tmp = new((/nHR/), float)
            hour_tmp@long_name = "hour of day"
            min_tmp@long_name = "minute of day"

            ;; Mean diurnal cycle
            dc = new((/nHR/), float)
            dc!0 = "hour"

            ;; Dc sorted by increasing LT
            dc_sort = new((/nHR/), float)

            ;; Diurnalcycle output
            diurnalcycle = new((/nHR + 2/), float)
            diurnalcycle!0 = "hour"

            ;; UTC and LT variables
            utc = new((/nHR/), float)
            LT = new((/nHR/), float)
            LT_sort = new((/nHR/), float)
            localhour = new((/nHR + 2/), float)

            ;; Compute box mean
            var_box = dim_avg_n(A0(:, :, :), (/1, 2/))

            ;; Compute mean for each hour
            timeh = cd_calendar(A0&time, 0)
            hour_tmp = timeh(0:nHR - 1, 3)
            min_tmp = timeh(0:nHR - 1, 4)
            n = -1
            do nh = 0, nHR - 1
                n = n + 1
                ;; Find indices for each hour
                indh = ind(timeh(:, 3) .eq. hour_tmp(nh))

                ;; Compute time
                utc(n) = hour_tmp(nh) + min_tmp(nh) / 60.

                ;; Compute mean diurnal cycle
                dc(n) = dim_avg(var_box(indh))
                delete([/indh/])
            end do

            ;; Compute local hour
            do nh = 0, nHR - 1
                LT(nh) = utc(nh) + (min_lon + max_lon) / 2. * 24. / 360.
                if (LT(nh) .gt. 24.) then
                    LT(nh) = LT(nh) - 24.
                end if
                if (LT(nh) .lt. 0.) then
                    LT(nh) = LT(nh) + 24.
                end if
            end do

            ;; Sort by increasing local time
            ip = dim_pqsort_n(LT, 1, 0)
            LT_sort(:) = LT(ip(:))
            dc_sort(:) = dc(ip(:))

            ;; Duplicate first and last hour values to make plots periodic

            do nh = 1, nHR
                localhour(nh) = LT_sort(nh - 1)
                diurnalcycle(nh) = dc_sort(nh - 1)
            end do
            localhour(0) = LT_sort(0) - (LT_sort(1) - LT_sort(0))
            diurnalcycle(0) = dc_sort(nHR - 1)
            localhour(nHR + 1) = LT_sort(nHR - 1)\
                                 + (LT_sort(nHR - 1) - LT_sort(nHR - 2))
            diurnalcycle(nHR + 1) = dc_sort(0)

            diurnalcycle&hour = localhour
            diurnalcycle@units = A0@units

            delete(A0)
            delete(timeh)
            delete(var_box)
            delete(ip)
            delete(dc)
            delete(dc_sort)
            delete(LT)
            delete(LT_sort)
            delete(utc)
            delete(hour_tmp)
            delete(min_tmp)

            ;; Check if start/end years are consistent across all models
            if ((first_start_year .ne. models@start_year(idx_mod))\
                 .or. (first_end_year .ne. models@end_year(idx_mod))) then
                model_years_consistent = False
            end if

            ;; Store diurnal cycle
            diurnalcycle@idx = idx_mod
            diurnalcycle@model_name = models@name(idx_mod)
            storage_variable_name = str_join((/"diurnalcycle", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(diurnalcycle, storage_vault, storage_variable_name)

            delete(diurnalcycle)
            delete(localhour)
            delete(nHR)

        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/"diurnalcycle"/)
        do idx_var = 0, dimsizes(stored_variables) - 1
            var_name = stored_variables(idx_var)
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

    dictkey_units = interface_get_dictkeys(var0, 0) + "_attributes_units"
    lb_units = nml_data_vault@$dictkey_units$

    res = True

    res@tmXBLabelFontHeightF = 0.008
    res@tmYLLabelFontHeightF = 0.008
    res@cnInfoLabelOn = False
    res@tmYRLabelsOn = False
    res@tmYRMode = "Automatic"
    res@tiMainFontHeightF = 0.015
    res@trXMinF = 0. 
    res@trXMaxF = 24.
    orig_res = res

    plot_settings = True
    plot_settings@x_axis_coord = "hour"
    plot_settings@gsnAddCyclic = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_title = lb_units
    plot_settings@part_of_header = "Precipitation"

    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ "\
                            + di@season\
                            + " average across indicated years"

    plot_settings@gsn_right_string = variable_info@long_name
    plot_settings@plot_pres_hgt = True

    if (model_years_consistent) then
        years_string = "~F22~yrs:~F21~ "\
                       + models@start_year(0)\
                       + "-" + models@end_year(0)
    else
        years_string = "~F22~yrs:~F21~ inconsistent"
    end if

    plot_settings@aux_filename_info = di@box + "-" + di@season + "-diurnal-cycle"
    plot_settings@title_string = di@season + " Diurnal cycle over " + di@box
    plot_settings@yaxis_string = variable_info@long_name
    plot_settings@xaxis_string = "Local hour"
    plot_settings@type = "Not applicable"
    plot_settings@left_plot_subheader = latlon_header_string + ", " + years_string
    di@minmax = variable_info@minmax

    xy_line_overlap(storage_vault, di, plot_settings, "diurnalcycle", debuginfo)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
