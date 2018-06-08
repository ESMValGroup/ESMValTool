diag_script_info = True

;; Used by: all tsline routines for temperature amonalies
diag_script_info@styleset = "CMIP5"         ; "CMIP5", "DEFAULT"

;; Used by: Tas_tsline
;; see style_FUNCTIONS/style_$project$.ncl for which data shall be included in the statistics
diag_script_info@multi_model_mean = "y"     ; "y" = plot multi-model mean & stddev; "n" = don't
                                            ; "mean" = plot multi-model mean; "stddev" = plot multi-model mean & stddev; "n" = don't
;diag_script_info@EMs_in_lg = False          ; create legend label for each individual ensemble member

;; Level selection
diag_script_info@ts_minlev = -1.
diag_script_info@ts_maxlev = -1.

;; Region selection
diag_script_info@ts_minlat = -90.
diag_script_info@ts_maxlat = 90.
diag_script_info@ts_minlon = 0.
diag_script_info@ts_maxlon = 360.

diag_script_info@time_avg   = "yearly"
diag_script_info@ts_anomaly = "anom" ; "anom" for anomaly-plot else "noanom"
diag_script_info@ref_start  = 1961 ;start year of reference period for anomalies
diag_script_info@ref_end    = 1990 ;end year of reference period for anomalies

