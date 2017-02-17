diag_script_info = True

;; Used by SouthernOcean_vector_polcon_diff
;; required settings
diag_script_info@styleset     = "CMIP5"           ; plot style, "CMIP5", "DEFAULT"
diag_script_info@month        = "A"               ; A = annual mean, 3 = March, 9 = September
diag_script_info@range_option = 1                 ; 0 = use each model's time range specified in namelist
                                                  ; 1 = use only intersection of all time ranges 
diag_script_info@region       = "Antarctic"       ; region name, must be Antarctic currently
diag_script_info@ref_model    = "ERA-Interim"     ; reference model, as specified in annotations -> if the specified string is
                                                  ;  not found, the routine will print a list of valid strings before crashing

;; optional settings
diag_script_info@max_lat      = -40.0             ; For Arctic, use min_lat, for Antarctic use max_lat
diag_script_info@max_vert     = 2                 ; max allowed number of rows on a panel page (vertical)
diag_script_info@max_hori     = 1                 ; max allowed number of columns on a panel page (horizontal)
diag_script_info@grid_min     = -0.1              ; min for contour plot scale
diag_script_info@grid_max     =  0.1              ; max for contour plot scale
diag_script_info@grid_step    =  0.01             ; step size for contour plot scale
diag_script_info@colormap     = "ncl_default"     ; contour colormap, e.g. WhiteBlueGreenYellowRed, rainbow
diag_script_info@vcmagnitude  = 0.25              ; magnitude of vectors (larger value, shorter arrows)
diag_script_info@vcmindist    = 0.04              ; controls density of vectors (larger value, less arrows)
diag_script_info@vccolormod   = "green4"          ; color of vectors for current model
diag_script_info@vccolorref   = "black"           ; color of vectors for reference model
diag_script_info@showunits    = True              ; display units in figure title
