

diag_script_info = True

diag_script_info@month = "A"                ; A = annual mean, 
diag_script_info@styleset = "EMAC"         ; "CMIP5", "DEFAULT", "GO"

;; Used by: GO_tsline.ncl
;; see style_FUNCTIONS/style_$project$.ncl for which data shall be included in the statistics
;;diag_script_info@ts_multimean = False    ; "y" = plot multi-model mean & stddev; "n" = don't
;;diag_script_info@refModel=(/"ref"/)

;; North Pacific 
diag_script_info@ts_latrange = (/45,90/)
diag_script_info@ts_lonrange = (/140,280/)
;; mid North Hemisphere Pacific 
diag_script_info@ts_latrange = (/15,45/)
diag_script_info@ts_lonrange = (/140,280/)
;; Equatorial Pacific 
diag_script_info@ts_latrange = (/-15,15/)
diag_script_info@ts_lonrange = (/140, 280/)
;; mid South Hemisphere Pacific 
diag_script_info@ts_latrange = (/-45,-15/)
diag_script_info@ts_lonrange = (/140, 280/)
;;  South Pacific 
diag_script_info@ts_latrange = (/-90,-45/)
diag_script_info@ts_lonrange = (/140, 280/)
;;  All Pacific 
diag_script_info@ts_latrange = (/-90,90/)
diag_script_info@ts_lonrange = (/140, 280/)
;; North Atlantic
diag_script_info@ts_latrange = (/45,90/)
diag_script_info@ts_lonrange = (/290, 30/)
;; mid North Hemisphere Atlantic
diag_script_info@ts_latrange = (/15,45/)
diag_script_info@ts_lonrange = (/290, 30/)
;; Equatorial Atlantic
diag_script_info@ts_latrange = (/-15,15/)
diag_script_info@ts_lonrange = (/290, 30/)
;; mid South Hemisphere Atlantic
diag_script_info@ts_latrange = (/-45,-15/)
diag_script_info@ts_lonrange = (/290, 30/)
;; South Atlantic
diag_script_info@ts_latrange = (/-90,-45/)
diag_script_info@ts_lonrange = (/290, 30/)
;;  Atlantic
diag_script_info@ts_latrange = (/-90,90/)
diag_script_info@ts_lonrange = (/290, 30/)
;; South Indian
;;diag_script_info@ts_latrange = (/-90,-45/)
;;diag_script_info@ts_lonrange = (/20, 146/)
;; Mid Indian
;;diag_script_info@ts_latrange = (/-45,-15/)
;;diag_script_info@ts_lonrange = (/20, 146/)
;; North Indian
;;diag_script_info@ts_latrange = (/-15,30/)
;;diag_script_info@ts_lonrange = (/20, 146/)
;; global
diag_script_info@ts_latrange = (/-90,90/)
diag_script_info@ts_lonrange = (/0, 360/)
diag_script_info@ts_ymin =  10.4
diag_script_info@ts_ymax = 11.8
diag_script_info@ts_min_year = 1985
diag_script_info@ts_max_year = 2013
diag_script_info@ts_anom = False
diag_script_info@ts_smooth = False
diag_script_info@ts_trend = True 
diag_script_info@ts_coverage = False
