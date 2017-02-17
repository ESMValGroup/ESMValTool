;*********************************************************************
; SAMonsoon_precip_propagation.ncl
;*********************************************************************
; NCL script
; SAMonsoon_precip_propagation.ncl, January 2013
; Martin.Evalsson@smhi.se
;*********************************************************************
; PIs: Klaus Wyser (SMHI) and Veronika Eyring (DLR)
; Lagged correlation of precipitation anomalies with respect to area
; avearged region in Indian Ocean
;
; Reference: Adapted from the MJO Clivar scripts available at,
;
; https://www.ncl.ucar.edu/Applications/Scripts/mjoclivar_8.ncl
; and
; https://www.ncl.ucar.edu/Applications/Scripts/mjoclivar_9.ncl
;http://ncl.ucar.edu/Applications/mjoclivar.shtml
;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;*********************************************************************
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

load "./diag_scripts/aux/mjo/mjo.ncl"
load "./plot_scripts/ncl/mjo_level1.ncl"
load "./plot_scripts/ncl/mjo_level2.ncl"


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
    var1 = variables(1)
    var2 = variables(2)
    field_type0 = field_types(0)
    field_type1 = field_types(1)
    field_type2 = field_types(2)

    ;; ============================================
    ;; ||                                        ||
    ;; ||          REFERENCES INFO               ||
    ;; ||                                        ||
    ;; ============================================
    write_references(diag_script,       \  ; script name
                     "A_wang_sh",       \  ; authors
                     "",                \  ; contributors
                     "",                \  ; diag_references
                     "",                \  ; obs_references
                     "P_embrace")          ; proj_acknowledgements
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
    required_attributes = (/"latrange_life",\
                            "ihp",\
                            "nWgt",\
                            "sigma",\
                            "fca",\
                            "fcb",\
                            "season_life_cycle"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info
    di@is_cyclic = False

    ;; ============================================
    ;; ||                                        ||
    ;; ||         COMPUTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================

   ;; Reference model index
    idx_ref_pr = get_ref_model_idx(models, var_attr_ref(0))
    idx_ref_ua = get_ref_model_idx(models, var_attr_ref(1))

    ;; -----------------------------------------------
    ;; Reset min_idx/max_idx
    ;; -----------------------------------------------
    min_idx = min((/idx_ref_pr, idx_ref_ua/))
    max_idx = max((/idx_ref_pr, idx_ref_ua/))

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
            nml_data_vault_va = read_data(idx_mod, var2, field_type2)
        end if
    end do

    optPlot = True

    ;; ***************************************************************
    ;; ------------------ Open files for plots -----------------------

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

    xuv_all = calc_daily_olr_uv_life_cycle(nml_data_vault_pr,\
                                           nml_data_vault_ua,\
                                           nml_data_vault_va,\
                                           di,\
                                           idx_mod_pr,\
                                           idx_mod_ua)
    olr = xuv_all[0]
    u850 = xuv_all[1]
    v850 = xuv_all[2]

    mjo_life_cycle_plot(olr,\
                        u850,\
                        v850,\
                        idx_mod_pr,\
                        di)

    delete(xuv_all)
    delete(olr)
    delete(u850)
    delete(v850)
    delete(idx_mod_pr)
    delete(idx_mod_ua)

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

            xuv_all = calc_daily_olr_uv_life_cycle(nml_data_vault_pr,\
                                                   nml_data_vault_ua,\
                                                   nml_data_vault_va,\
                                                   di,\
                                                   idx_mod_pr,\
                                                   idx_mod_ua)
            olr = xuv_all[0]
            u850 = xuv_all[1]
            v850 = xuv_all[2]
            mjo_life_cycle_plot(olr,\
                                u850,\
                                v850,\
                                idx_mod_pr,\
                                di)

            delete(xuv_all)
            delete(olr)
            delete(u850)
            delete(v850)
            delete(idx_mod_pr)
            delete(idx_mod_ua)
        end if
    end do

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
