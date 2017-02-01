diag_script_info = True

; contour level settings for various plot_types
diag_script_info@acp_cnlvl  = (/1.0625, 2.125, 3.1875, 4.25, 5.3125, 6.375, 7.4375, 8.5,\
                              9.5625, 10.625, 11.6875, 12.75, 13.8125, 14.875, 15.9375/)
diag_script_info@acp_dcnlvl = (/-30,-26,-22,-18,-14,-10,-6,-2, 2, 6, 10, 14, 18, 22, 26, 30/) * 1e-6
diag_script_info@rgb_file   = "amwg.rgb"  ; Colour map

; lat-/longranges (implemented in the resp. plot_type)
;diag_script_info@ac_latrange = (/-20, 40/)
;diag_script_info@ac_lonrange = (/77, 100/)

; Reference model (e.g., observations, used by some plot_types)
diag_script_info@refModel = (/"TRMM-3B43-v7-0.25deg"/)
