;;#############################################################################
;; Author: Shiyu Wang (SMHI, Sweden)
;; Contributor: Mark Stevens (NCAR, USA), wrote the original script
;; Contributor: Richard Neale (NCAR, USA), wrote the original script
;; Contributor: Martin Evaldssson (SMHI, Sweden), some ESMValTool adapations
;; EMBRACE project
;;#############################################################################
;; Description:
;;    Creates longitude vs latitude 2D-field contour plots (set5) for the
;;    AMWG diagnostics package (http://www.cgd.ucar.edu/cms/diagnostics/)
;;    The original script has been modified for use within the CCMVal
;;    diagnostic software package. The CCMVal script script has been
;;    furthered developed as a part of the ESMValTool.
;;
;; Required diag_script_info attributes (variable specific):
;;    latrange_basic:    Latitude range for plot
;;    lonrange_basic:    Longituge range for plot
;;    cn_levels_pr:      Contour levels for precip
;;    season_names:      Name of the seasons (for display)
;;    seasons:           Seasons
;;    season_mjo:
;;
;; Optional diag_script_info attributes (variable specific):
;;
;; Required variable_info attributes (variable specific):
;;    long_name:         Variable long name
;;        units:         Variable unit
;;
;; Optional variable_info attributes (variable specific):
;;
;; Caveats:
;;
;; Reference: Adapted from the MJO Clivar scripts available at
;;            http://ncl.ucar.edu/Applications/mjoclivar.shtml
;;
;; Modification history
;;     20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                         of the code
;;     20150828-A_wang_sh: written
;;    * 200801: ERA40 reanalysis
;;    * 200303: ISCCP fields
;;    * 200101: updated with significance tests
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/auxiliary.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"

load "./diag_scripts/aux/mjo/mjo.ncl"
load "./plot_scripts/ncl/mjo_level1.ncl"
load "./plot_scripts/ncl/mjo_level2.ncl"

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

    ;; ============================================
    ;; ||                                        ||
    ;; ||          REFERENCES INFO               ||
    ;; ||                                        ||
    ;; ============================================
    write_references(diag_script,        \  ; script name
                     "A_wang_sh",        \  ; authors
                     "",                 \  ; contributors
                     "",                 \  ; diag_references
                     "",                 \  ; obs_references
                     "P_embrace")           ; proj_acknowledgements
end

; Load var_info
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
    valid_plot_field = (/"T2Ms"/)
    ind_f = inlist(field_type0, valid_plot_field)
    if (any(ind_f .eq. False)) then
        error_msg("fatal", diag_script, "", "can't use field " + field_type0)
        status_exit(1)
    end if
    delete(ind_f)
    ind_f = inlist(field_type1, valid_plot_field)
    if (any(ind_f .eq. False)) then
        error_msg("fatal", diag_script, "", "can't use field " + field_type1)
        status_exit(1)
    end if
    delete(ind_f)
    delete(valid_plot_field)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the mjo cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"season_names", \
                            "seasons", \
                            "season_mjo", \
                            "latrange_basic", \
                            "lonrange_basic", \
                            "cn_levels_pr"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    ;; ---------------------------------------------
    ;; Data structure for storing the models defined
    ;; in the namelist
    ;; ---------------------------------------------
    nml_data_vault_pr = True
    nml_data_vault_ua = True

    do idx_mod = 0, dim_MOD - 1
        ;; Skip models<->parameter combinations that are explicitly
        ;; excluded in the xml-file
        curr_var_idx = ind(var0 .eq. variables)
        if (var_attr_exclude(curr_var_idx) .ne. model_attr_id(idx_mod)) then
            nml_data_vault_pr = read_data(idx_mod, var0, field_type0)
        end if
        curr_var_idx = ind(var1 .eq. variables)
        if (var_attr_exclude(curr_var_idx) .ne. model_attr_id(idx_mod)) then
            nml_data_vault_ua = read_data(idx_mod, var1, field_type1)
        end if
    end do

  ;; Reference model index
    idx_ref_pr = get_ref_model_idx(models, var_attr_ref(0))
    idx_ref_ua = get_ref_model_idx(models, var_attr_ref(1))

    ;; -----------------------------------------------
    ;; Reset min_idx/max_idx
    ;; -----------------------------------------------
    min_idx = min((/idx_ref_pr, idx_ref_ua/))
    max_idx = max((/idx_ref_pr, idx_ref_ua/))

    ;; 'output_file_type' if fetched from ncl.interface
    if (.not. isdefined("output_file_type")) then
        output_file_type = "ps"
    end if

    season_names = di@season_names
    seasons = di@seasons
    i_season = di@season_mjo - 1
    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    optPlot = True
    optPlot@gsnRightString  = "U850 (m/s): contour"

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;;loop through seasons
    do s = 0, dimsizes(seasons) - 1

        ;; ***************************************************************
        ;; ------------------ Open files for plots -----------------------
        aux_title_info = season_names(s)
        output_filename = interface_get_figure_filename(diag_script_base, \
                                                        var0, \
                                                        field_type0, \
                                                        aux_title_info, \
                                                        -1)
        wks = gsn_open_wks(output_file_type, output_dir + output_filename)

        plot = new(dim_MOD, "graphic")
        gsn_define_colormap(wks, "amwg")

        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process data
            ;;------------------------
            ;; Check for reference data sets and plot reference
            if (any(idx_mod .eq. (/idx_ref_pr, idx_ref_ua/))) then

                ;; Check if there's a single reference data set
                if (idx_ref_pr .eq. idx_ref_ua) then
                    idx_mod_pr = idx_mod
                    idx_mod_ua = idx_mod

                ;; If there are two separate data sets, only use them once
                ;; Skip the second time (max_idx)
                else if (idx_mod .eq. max_idx) then
                    continue
                else  ; Use them the first time (min_idx)
                    idx_mod_pr = idx_ref_pr
                    idx_mod_ua = idx_ref_ua

                end if
                end if
            end if
        end do

        pr_ua_mean = calc_pr_ua_mean(nml_data_vault_pr,\
                                     nml_data_vault_ua,\
                                     di,\
                                     idx_mod_pr,\
                                     idx_mod_ua,\
                                     i_season(s, :))
        mean_pr = pr_ua_mean[0]
        mean_ua = pr_ua_mean[1]
        optPlot@gsnLeftString = models@name(idx_ref_pr)\
                              + " precip (mm/day): color"
        plot(0) = pr_u850_mean_plot(wks, mean_pr, mean_ua, di, optPlot)

        delete(pr_ua_mean)
        delete(mean_pr)
        delete(mean_ua)

        plot_n = 1

        do idx_mod = 0, dim_MOD - 1
            ;;------------------------
            ;; Extract and process non-reference data
            ;;------------------------
            ;; Check for reference data sets and plot reference
            if (any(idx_mod .eq. (/idx_ref_pr, idx_ref_ua/))) then
                continue
            else  ; non reference dataset
                idx_mod_pr = idx_mod
                idx_mod_ua = idx_mod

                pr_ua_mean = calc_pr_ua_mean(nml_data_vault_pr,\
                                             nml_data_vault_ua,\
                                             di,\
                                             idx_mod_pr,\
                                             idx_mod_ua,\
                                             i_season(s, :))
                mean_pr = pr_ua_mean[0]
                mean_ua = pr_ua_mean[1]
                optPlot@gsnLeftString = models@name(idx_mod_pr)\
                                      + " precip (mm/day): color"
                plot(plot_n) = pr_u850_mean_plot(wks,\
                                                 mean_pr,\
                                                 mean_ua,\
                                                 di,\
                                                 optPlot)

                plot_n = plot_n + 1

                delete(pr_ua_mean)
                delete(mean_pr)
                delete(mean_ua)
            end if
        end do  ; idx_mod loop

        ;; common resources for panel plots
        pan = True
        pan@gsnMaximize = True
        pan@gsnPaperOrientation = "portrait"
        pan@txString = season_names(s)
        pan@gsnPanelLabelBar = True

        n_p = dimsizes(plot) / 2
        if (mod(dimsizes(plot), 2) .ne. 0) then
            n_p = n_p + 1
        end if
        gsn_panel(wks, plot, (/n_p, 2/), pan)    ; panel plot

        delete(plot)
        delete(pan)
        delete(n_p)

    end do  ; end seasons

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
