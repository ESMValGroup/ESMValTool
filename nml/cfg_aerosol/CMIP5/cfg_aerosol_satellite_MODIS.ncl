; Settings used by aerosol_satellite diag_script
diag_script_info = True

; Reference model
diag_script_info@ref_model = "MODIS"

; Style
diag_script_info@styleset = "DEFAULT"

; Regridding method
diag_script_info@regrid_method = "finest"

; Time range option (0: model own time range, 1: intersection between models)
diag_script_info@range_option = 1

; Multi-model stats
diag_script_info@multi_model_mean = True
diag_script_info@multi_model_median = True
