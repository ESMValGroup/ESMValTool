;;#############################################################################
;; WAMonsoon_10W10E_3D_basic.ncl
;; Author: Romain Roehrig (Meteo France, France)
;; EMBRACE project
;;#############################################################################
;; Description
;;    Computes an average over 10W-10E of a lev/lat/lon variable
;;    (e.g. ua/va or theta) for the West African Monsoon
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    cn_levels_mean_basic:        Contour levels for mean plot
;;    cn_levels_mean_diff_basic:   Contour levels for diff plot
;;    latrange:                    Latitude crop range for plot
;;    levrange:                    Pressure level crop range for plot
;;    lonrange:                    Longitude crop range for plot (-10, 10)
;;    plottype_lbres:              Handle for labelbar settings
;;    season:                      Season, typically "JJAS"
;;    plot_stddev:                 Switch for plotting stddev
;;    diff_colormap:               Colormap for diff plots
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    cn_levels_stddev_basic:      Contour levels for stddev plot
;;    cn_levels_stddev_diff_basic: Contour levels for diff stddev plot
;;
;; Required variable_info attributes (variable specific)
;;    zonal_ref_cmap:    Color map for zonal plots
;;    zonal_ref_Levels:  Contour levels for zonal plot
;;    zonal_diff_Levels: Contour levels for zonal diff plot
;;    long_name:         Name displayed in plot
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

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"

load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./interface_scripts/auxiliary.ncl"
;load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"
exact_panel_position_file = "./diag_scripts/aux/WAMonsoon/exact_panel_positions_pres.ncl"
load "./plot_scripts/ncl/monsoon_panels.ncl"
load "./diag_scripts/aux/SAMonsoon/exact_panel_positions_generic.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon_panel_positions.ncl"

;*********************************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    ;; 'diag_script' if fetched from ncl.interface
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot - " + diag_script + " (var: " + variables(0) \
                          + ")", verbosity, 1)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    diag_script_base = basename(diag_script)

    ;;-----------------------
    ;; References Information
    ;;-----------------------
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

    data = True
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do

    ;; 'show_debuginfo' if fetched from ncl.interface
    if (.not. isdefined("show_debuginfo")) then
        debuginfo = False
    else
        debuginfo = show_debuginfo
    end if

    ;; Check valid plot fields
    valid_plot_fields = (/"T3M"/)
    ind_f = inlist(field_type0, valid_plot_fields)
    if (any(ind_f .eq. False)) then
        error_msg("fatal", diag_script, "", "can't use field " + field_type0)
        status_exit(1)
    end if
    delete(valid_plot_fields)
    delete(ind_f)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the SAMonsoon cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"season",\
                            "latrange",\
                            "levrange",\
                            "cn_levels_mean_basic",\
                            "cn_levels_mean_diff_basic",\
                            "plot_stddev",\
                            "diff_colormap"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    if (di@plot_stddev) then
      delete(required_attributes)
        required_attributes = (/"cn_levels_stddev_basic",\
                                "cn_levels_stddev_diff_basic"/)
        exit_if_missing_atts(diag_script_info, required_attributes)
    end if
    ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"lonflip", "season", "latrange", "levrange", "lonrange"/)
    transforms@season = di@season
    transforms@latrange = di@latrange
    transforms@lonrange = di@lonrange
    transforms@levrange = di@levrange

    debuginfo@latrange = di@latrange(0) + "-" + di@latrange(1)
    debuginfo@latrange = di@lonrange(0) + "-" + di@lonrange(1)
    debuginfo@levrange = di@levrange(0) + "-" + di@levrange(1)
    debuginfo@variable = var0

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))


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
            A0 = apply_transforms(transforms, A0_tmp)
            delete(A0_tmp)

            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------
            A0_lonavg = dim_avg_n_Wrap(A0, 3)
            var_mean = dim_avg_n_Wrap(A0_lonavg, 0)
            var_mean@units = var_mean@long_name + " (" + var_mean@units + ")"

            var_mean@idx = idx_mod
            var_mean@model_name = models@name(idx_mod)
            storage_variable_name = str_join((/"mean", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(var_mean, storage_vault, storage_variable_name)

            if (di@plot_stddev) then
                ;; ------------------------------------------
                ;; Compute and store stddev for model idx_mod
                ;; ------------------------------------------
                var_stddev = dim_stddev_n_Wrap(A0_lonavg, 0)
                var_stddev@idx = idx_mod
                var_stddev@model_name = models@name(idx_mod)
                storage_variable_name = str_join((/"stddev", sprinti("%i", idx_mod)/), str_vault_sep)
                store_variable(var_stddev, storage_vault, storage_variable_name)
  
                delete(var_stddev)
            end if
            delete(A0)
            delete(A0_lonavg)
            delete(var_mean)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        if (di@plot_stddev) then
            stored_variables = (/"mean", "stddev"/)
        else
            stored_variables = (/"mean"/)
        end if
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

    res@tmXBLabelFontHeightF = 0.008
    res@tmYLLabelFontHeightF = 0.008
    res@cnInfoLabelOn = False
    res@tmYRLabelsOn = False
    res@tmYROn = False
    res@tmYRMode = "Automatic"
    res@tiYAxisString = "pressure (mbar)"
    res@tiMainFontHeightF = 0.015
    orig_res = res

    ;; Three panel plot (reference, model, and diff)
    plot_settings = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_title = lb_units
    plot_settings@lb_units = lb_units
    plot_settings@part_of_header = "Precipitation"
    plot_settings@type = "Not applicable"
    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ " + di@season + " average across indicated years" +\
                            "~C~~F22~rmse def:~F21~ NCL function wgt_arearmse(...) computed with " +\
                            "~C~                  cosine weights in latitude." +\
                            "~C~~F22~corr def:~F21~ NCL function pattern_cor(...) computed with " +\
                            "~C~                  cosine weights in latitude. Correlation is between" +\
                            "~C~                  reference and model data."

    ;; Special cases for wind zonal/meriodnal winds
    if (variable_info@long_name .eq. "Eastward Wind (m s-1)") then
        variable_info@long_name = "Zonal wind"
    else if (variable_info@long_name .eq. "Northward Wind (m s-1)") then
        variable_info@long_name = "Meridional wind"
    end if
    end if

    plot_settings@gsn_right_string = variable_info@long_name
    plot_settings@plot_pres_hgt = True

    ;; Three panel plots (ref, model, diff)
    valid_statistics = (/"yrs"/)

    ;; Colormaps
    diff_colormap = read_colormap_file(di@diff_colormap)
    plot_settings@diff_palette = diff_colormap
    plot_settings@default_palette = "ncl_default"

    three_by_one_diff(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)
    if (di@plot_stddev) then
        three_by_one_diff(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)
    end if

    ;; Multipanel plots
    delete_VarAtts(res, (/"tiMainFontHeightF"/))
    plot_settings@type = "mean"
    plot_settings@part_of_header = " of " + variable_info@long_name
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    plot_settings@type = "diff"
    plot_settings@part_of_header = " of " + variable_info@long_name
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    if (di@plot_stddev) then
        plot_settings@type = "mean"
        plot_settings@part_of_header = " stddev of " + variable_info@long_name
        multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

        plot_settings@type = "diff"
        plot_settings@part_of_header = " of stddev of " + variable_info@long_name
        multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)
    end if

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
