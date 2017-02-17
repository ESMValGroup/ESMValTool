;;#############################################################################
;; GO_comp_map.ncl
;; Author: Oliver Andrews (UEA, Great Britain)
;;         Clare Enright (UEA, Great Britain)
;; EMBRACE project
;;#############################################################################
;; Description
;;      Evaluation of ocean biogeochemical cycles at global scales, as
;;      simulated by both ocean-only and coupled climate-carbon cycle models. 
;;      This diagnostic takes input from up to four models, one of which is
;;      defined as the reference, to produce contour maps of the time-mean
;;      (climatology) and standard deviation (inter-annual or inter-seasonal
;;      [e.g. JJAS] variability) together with the difference between
;;      the reference and each model.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;      GO_season:                    Season
;;      GO_detrend:                   Detrend data
;;      latrange_basic:               Latitude range cropping
;;      lonrange_basic:               Longitude range cropping
;;      cn_levels_mean_basic:         Contour levels for contours
;;      cn_levels_mean_diff_basic:    Contour levels for diffs
;;      cn_levels_stddev_basic:       Contour levels for stddev
;;      cn_levels_stddev_diff_basic:  Contour levels for stddev diffs
;;
;; Required variable_info attributes (variable specific)
;;    long_name:                    Name to store in derived varible
;;    units:                        Units for derived variable
;;
;; Caveats
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150626-A_andr_ol: written
;; ############################################################################

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
;;load "./interface_scripts/read_data.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"

exact_panel_position_file = \
    "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip.ncl"

load "./plot_scripts/ncl/GO_panels.ncl"
load "./diag_scripts/aux/GO/GO_apply_transforms.ncl"
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
    info_output("+++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot - " + diag_script + " (var: " + variables(0) \
                          + ")", verbosity, 1)
    info_output("+++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
    ;; Can only plot 2d field; should have defined a derived var containing
    ;; just one level by this stage
    if (field_type0 .eq. "TO3M") then
        field_type0 = "TO2M"
    end if
    diag_script_base = basename(diag_script)

    ;;-----------------------
    ;; References Information
    ;;-----------------------
    write_references(diag_script,            \  ; script name
                     (/"A_andr_ol", "A_enri_cl" /),  \  ; authors
                     (/"A_eval_ma", "A_levi_ri", "A_mart_gi"/), \  ; contributors
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
    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    ;; 'show_debuginfo' if fetched from ncl.interface
    if (.not. isdefined("show_debuginfo")) then
        debuginfo = False
    else
        debuginfo = show_debuginfo
    end if

    ;; Check valid plot fields
    valid_plot_field = (/"TO2Ms", "TO2M", "TO3M"/)
    ind_f = inlist(field_type0, valid_plot_field)
    if (any(ind_f .eq. False)) then
        error_msg("fatal", diag_script, "", "can't use field " + field_type0)
        status_exit(1)
    end if
    delete(ind_f)
    delete(valid_plot_field)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the GO_cfg/*.ncl) are present
    ;; -------------------------------------------------
    required_attributes = (/"GO_season",\
                            "GO_detrend", \
                            "latrange_basic",\
                            "lonrange_basic",\
                            "cn_levels_mean_basic",\
                            "cn_levels_mean_diff_basic",\
                            "cn_levels_stddev_basic",\
                            "cn_levels_stddev_diff_basic"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;;-----------------------------------------------------
    ;; Define transforms and meta data to apply to the plot
    ;;-----------------------------------------------------
    transforms = True
    if (di@GO_detrend .eq. True) then
        transforms@apply_these = (/"GO_season",\
                                   "GO_detrending",\
                                   "latrange",\
                                   "lonrange"/)
    else
        transforms@apply_these = (/"GO_season", "latrange", "lonrange"/)
    end if
    transforms@GO_season = di@GO_season
    transforms@GO_detrend = di@GO_detrend
    transforms@latrange = di@latrange_basic
    transforms@lonrange = di@lonrange_basic

    debuginfo@GO_season = di@GO_season
    debuginfo@GO_detrend = di@GO_detrend
    debuginfo@latrange = di@latrange_basic(0) + "-" + di@latrange_basic(1)
    debuginfo@lonrange = di@lonrange_basic(0) + "-" + di@lonrange_basic(1)
    debuginfo@variable = var0

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

    ;; Ignore some warnings from ESMF-regridding
    ignore_warning_file = "interface_data/warnings_to_ignore.txt"
    warnings_to_ignore = \
        (/"ESMF_regrid_with_weights: warning: destination grid is not", \
          "warning:Attempt to reference attribute \(region\) which is undefined", \
          "warning:Attempt to reference attribute \(long_name\) which is undefined", \
          "warning:Attempt to reference attribute \(SrcFileName\) which is undefined",\
          "warning:ContourPlotSetValues: Data values out of range of levels set by EXPLICITLEVELS mode",\
          "warning:ContourPlotInitialize: scalar field is constant\; no contour lines will appear\; use cnConstFEnableFill to enable fill"/)
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
    nml_data_vault_var = True
    do idx_mod = 0, dim_MOD - 1
        nml_data_vault_var = read_data(idx_mod, var0, field_type0)
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
            var_tmp = extract_data(idx_mod, nml_data_vault_var, -1, 0, 0)
            var_ref = extract_data(idx_ref, nml_data_vault_var, -1, 0, 0)

            filename = interface_get_fullpath(var0, field_type0, idx_mod)
            fil = addfile(filename, "r")
            varNames = getfilevarnames(fil)
            delete(fil)

            ;; Check for irregular grid
            regrid_projects = di@regrid_for_these_projects
            ind_r = inlist(models@project(idx_mod), regrid_projects)

            if (all(ind_r .eq. False)) then
                info_output("no of times before transforms " + dimsizes(var_tmp), verbosity, 1)
                var1 = apply_transforms(transforms, var_tmp)
                var = GO_apply_transforms(transforms, var1)
                delete(var1)
                delete(var_tmp)
            else
                ;; Regrid; data must already have been extracted for correct level
                ;;         Get 2d lat & lon coordinates
                var = make_latlon2D(var0, idx_mod, nml_data_vault_var, var_tmp)
                delete(var_tmp)
                ;; Regridding
                stem = interface_get_inpaths(idx_mod) + "/ESMF_regrid_"\
                                                      + models@project(idx_mod) + "_"\
                                                      + models@name(idx_mod) + "_"
                ending = var0 + "_" + idx_mod + "_" \
                              + diag_script_info@GO_season + "_"\
                              + models@start_year(idx_mod) + "-"\
                              + models@end_year(idx_mod) + ".nc"

                ;; This could be done with transofrom 2interp" but want more
                ;; control of interpolation at this stage
                Opt = True
                Opt@DstGridLat = var_ref&lat
                Opt@DstGridLon = var_ref&lon
                ;; Set mask for missing values according to data for first year
                ;; this should be ok for model output
                Opt@SrcMask2D = where(ismissing(var(1, :, :)), 0, 1)

                Opt@DstGridType = "rectilinear"      ; Destination grid
                Opt@ForceOverwrite = True
                Opt@PrintTimings = True
                Opt@Debug = False
                Opt@ESMFDebugFileName = stem + "Dbg_" + ending

                ; "bilinear" ; "patch"  ;
                ;  "conserve" fails with a 'concave source cell' for HadGEM models
                ; see https://www.ncl.ucar.edu/Document/Functions/ESMF/ESMF_regrid.sh
                ;     if (ismissing(str_match_ind_ic(models@name, "HadGEM2-CC"))) then
                ;         Opt@InterpMethod = "conserve"
                ;     else
                Opt@InterpMethod = "bilinear"

                ;; Save regrid files for later use: careful with naming!
                Opt@WgtFileName = stem + "Wgt_" + ending
                Opt@SrcFileName = stem + "Src_" + ending
                Opt@DstFileName = stem + "Dst_" + ending

                ;; ESMF_regrid may issue a warning
                ;; (e.g. if Dst is not covered Src grid)
                ;; -> set "exit_on_warning" to False (in main nml) to continue to the
                ;; next diag

                ;; Setting Opt@SrcRegional=True avoids artefacts
                ;; from incomplete Src grids
                ;; Activate if only part of the sphere is given
                Opt@SrcRegional = True

                ;; Reuse regrid files if present
                if (isfilepresent_esmval(Opt@SrcFileName)) then
                    Opt@SkipSrcGen = True
                end if
                if (isfilepresent_esmval(Opt@DstFileName)) then
                    Opt@SkipDstGen = True
                end if

                if (isfilepresent_esmval(Opt@WgtFileName)) then
                    Opt@SkipWgtGen = True
                    var_regrid = ESMF_regrid_with_weights(var, Opt@WgtFileName, Opt)
                else
                    var_regrid = ESMF_regrid(var, Opt)
                end if

                var_regrid@long_name = variable_info@long_name
                ;; Clean up
                ;; Could re-use file sif you knew all input models on same grid
                system("\rm " + Opt@WgtFileName)
                system("\rm " + Opt@SrcFileName)
                system("\rm " + Opt@DstFileName)
                delete([/var, Opt/])
                var1 = apply_transforms(transforms, var_regrid)
                var = GO_apply_transforms(transforms, var1)
                delete(var1)
                delete(var_regrid)
            end if
            delete(varNames)
            mean_var = dim_avg_n_Wrap(var, 0)

            ;; ----------------------------------------
            ;; Compute and store mean for model idx_mod
            ;; ----------------------------------------
            mean_var@idx = idx_mod
            mean_var@model_name = models@name(idx_mod)
            mean_var@Regional = True
            mean_var@GridType = "rectilinear"
            ;; must have variable name somewhere?
            storage_variable_name = str_join((/"mean",  sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(mean_var, storage_vault, storage_variable_name)
            delete(mean_var)

            ;; ------------------------------------------
            ;; Compute and store stddev for model idx_mod
            ;; ------------------------------------------
            stddev_var = dim_stddev_n_Wrap(var, 0)
            stddev_var@idx = idx_mod
            stddev_var@long_name = variable_info@long_name
            stddev_var@model_name = models@name(idx_mod)
            stddev_var@Regional = True
            stddev_var@GridType = "rectilinear"
            storage_variable_name = str_join((/"stddev", sprinti("%i", idx_mod)/), str_vault_sep)
            store_variable(stddev_var, storage_vault, storage_variable_name)
            delete(var)
            delete(stddev_var)
        end do  ; idx_mod = 0, dim_MOD - 1

        ;; Write vault to netCDF
        if (write_netcdf) then
            vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
        end if

    else  ; Retrieve computed diagnostic fields from netCDF
        stored_variables = (/ "mean", "stddev"/)
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
    lb_units = nml_data_vault_var@$dictkey_units$
    ;; Copy season string to di@season for plotting routines to read
    di@season = di@GO_season
    ;; -------------------------------
    ;; General resources for this plot
    ;; -------------------------------
    res = True
    res = map_zoom(res, transforms)
    res@mpCenterLonF = 180.
    res@mpFillOn = False

    orig_res = res

    ;; Three panel plot (reference, model, and diff)
    plot_settings = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_title = lb_units
    plot_settings@lb_units = lb_units
    plot_settings@type = "Not applicable"

    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ " + di@GO_season + " average across indicated years"

    ;; Three panel plots (ref, model, diff)
    valid_statistics = (/"mean"/)
    plot_settings@gsn_right_string = variable_info@long_name
    plot_settings@skip_long_name_in_plots = True
    three_by_one_diff(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    three_by_one_diff(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    res = orig_res
    plot_settings@skip_long_name_in_plots = True
    plot_settings@gsn_right_string = False

    ;; Multipanel plots
    plot_settings@type = "mean"
    plot_settings@part_of_header = " of " + variable_info@long_name
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    plot_settings@type = "diff"
    plot_settings@part_of_header = " of " + variable_info@long_name
    multipanel(storage_vault, di, plot_settings, "mean", debuginfo, valid_statistics, res)

    plot_settings@type = "mean"
    plot_settings@part_of_header = " of stddev of " + variable_info@long_name
    multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    plot_settings@type = "diff"
    plot_settings@part_of_header = " of stddev of " + variable_info@long_name
    multipanel(storage_vault, di, plot_settings, "stddev", debuginfo, valid_statistics, res)

    ;; Reference plot

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
