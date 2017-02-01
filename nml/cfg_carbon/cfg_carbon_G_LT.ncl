diag_script_info = True
diag_script_info@styleset = "CMIP5"

;info for carbon_constraint.ncl
diag_script_info@con_latrange  = (/-30,30/)  ; Latitude
diag_script_info@con_lonrange  = (/0,360/)   ; Longitude 
diag_script_info@con_units     = "GtC/K"
diag_script_info@pcBGC         = "esmFixClim1"
;specy models to be included in emergent constraint ("all";index)
diag_script_info@reg_models    = (/0,1,2,3,4,5,7/);"all";
diag_script_info@scatter_log   = False
diag_script_info@xy_line_legend = False
