diag_script_info = True
;; Used by: all SeaIce routines
;;diag_script_info@region = "Antarctic"       ; "Arctic" or "Antarctic"; entire hemisphere will be evaluated

diag_script_info@month = "A"                ; A = annual mean, 3 = March, 9 = September
diag_script_info@styleset = "EMAC"         ; "CMIP5", "DEFAULT", "EMAC"
diag_script_info@colors = (/"red",    "green"/)
;; Used by: chl_tsline.ncl
;; see style_FUNCTIONS/style_$project$.ncl for which data shall be included in the statistics
diag_script_info@ts_multimean = False    ; "y" = plot multi-model mean & stddev; "n" = don't
;diag_script_info@refModel=(/"ref"/) ; only if models on same grid
diag_script_info@ts_latrange = (/-90, 90/)
diag_script_info@ts_lonrange = (/-180, 180/)
diag_script_info@ts_ymin =  0.0000000
diag_script_info@ts_ymax = .52
diag_script_info@ts_min_year = 1998
diag_script_info@ts_max_year = 2011
diag_script_info@ts_anom = False
diag_script_info@ts_smooth = False
diag_script_info@ts_trend = False
diag_script_info@ts_coverage = True