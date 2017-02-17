; Settings used by aerosol_satellite diag_script
diag_script_info = True

; Reference model
diag_script_info@ref_model = "MIMO"

; Style
diag_script_info@styleset = "DEFAULT"

; Regridding method
diag_script_info@regrid_method = "coarsest"

; Time range option (0: model own time range, 1: intersection between models)
diag_script_info@range_option = 1

; Separate files for each plot
diag_script_info@separate = True

; Customize plot
; Plot resources can be specified using the following syntax:
; diag_script_info@res_(abs|dif)_MODEL_RESOURCE = VALUE
;     (abs|dif): use "abs" for the absolute value plots, "dif" for the
;                difference plot(s)
;     MODEL    : use the model names as specified in your main namelist (in the
;                <model> tags)
;     RESOURCE : NCL resource that is applicable to a contour_map_ce plot
;     VALUE    : value to be specified for RESOURCE
diag_script_info@res_abs_34128_tiMainString = \
    "Aerosol optical depth ~F19~@~F21~ EMAC (MADE3)"
diag_script_info@res_abs_34128_gsnLeftString = ""
attname = "res_abs_" + diag_script_info@ref_model + "_tiMainString"
diag_script_info@$attname$ = "Aerosol optical depth ~F19~@~F21~ MISR/MODIS"
attname = "res_abs_" + diag_script_info@ref_model + "_gsnLeftString"
diag_script_info@$attname$ = ""
diag_script_info@res_dif_34128_tiMainString = \
    "Aerosol optical depth ~F19~@~F21~ EMAC (MADE3) minus MISR/MODIS"
diag_script_info@res_dif_34128_gsnLeftString = ""
