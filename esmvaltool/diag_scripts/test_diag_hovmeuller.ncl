;*********************************************************************
; test_diag_hovmeuller.ncl
;*********************************************************************
; NCL script
; test_diag_hovmeuller.ncl, March 2014
; Martin.Evalsson@smhi.se
;*********************************************************************
; Very simple diagnostics used for testing hovmeuller plots
;
; The actual metrics are ported from the UK MetOffice IDL-monsoon
; evaluation scripts.
;
;*********************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"



;*********************************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot - " + diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

    do idx_mod = 0, dim_MOD - 1
        data1 = read_data(idx_mod, var0, field_type0)
        data2 = read_data(idx_mod, var1, field_type1)
    end do

    ;; -------------------------
    ;; Get environment variables
    ;; -------------------------
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    debuginfo = getenv("ESMValTool_debuginfo")
    if(ismissing(debuginfo)) then
        debuginfo = "False"
    end if

    ;; Check valid plot fields
    plot_fields = (/"T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field " + field_type0)
    end if
    delete(ind_f)

    ind_f = inlist(field_type1, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field " + field_type1)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; Output dir
    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; -------------------------------
    ;; General resources for this plot
    ;; -------------------------------
    res = True
    res@cnFillOn = True
    res@cnLinesOn = False

    ;; ------------------------------------
    ;; Loop over models
    ;; ------------------------------------
    do idx_mod = 0, dim_MOD - 1

        ;; --------------------------------------
        ;; Extract and modify data for variable 0
        ;; --------------------------------------
        A0 = extract_data_extended(idx_mod, var0, data1, -1, 0, 0)

        ;; -------------------------------------------------------
        ;; Define output workstation for mean plots for variable 0
        ;; -------------------------------------------------------
        aux_title_info = ""
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0,\
                                                        field_type0,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(file_type, output_file_path)

        map = gsn_csm_hov(wks, dim_avg_n_Wrap(A0, 0), res)

        delete(A0)
        delete(wks)

        ;; --------------------------------------
        ;; Extract and modify data for variable 1
        ;; --------------------------------------
        A0 = extract_data_extended(idx_mod, var1, data2, -1, 0, 0)

        ;; -------------------------------------------------------
        ;; Define output workstation for mean plots for variable 1
        ;; -------------------------------------------------------
        aux_title_info = ""
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var1,\
                                                        field_type0,\
                                                        aux_title_info,\
                                                        idx_mod)
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(file_type, output_file_path)

        map = gsn_csm_hov(wks, dim_avg_n_Wrap(A0, 0), res)

        delete(A0)
        delete(wks)
    end do
    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
