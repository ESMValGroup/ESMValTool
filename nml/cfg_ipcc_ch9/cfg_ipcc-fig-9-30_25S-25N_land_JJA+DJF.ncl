diag_script_info = True
diag_script_info@box = "25S-25N"

diag_script_info@season = "JJA"
diag_script_info@season2 = "DJF"

diag_script_info@styleset = "CMIP5"

diag_script_info@plotfile_suffix = "Land"

diag_script_info@latrange = (/ -25, 25/)      
diag_script_info@lonrange = (/ -180, 180/)    ; valid range = -180...180

diag_script_info@masking = 1   ; 0 = no masking (default)
                               ; 1 = land points only
                               ; 2 = ocean points only

diag_script_info@minmax = (/-70.0, 150.0/)

diag_script_info@show_anomalies = True

diag_script_info@curve_interpolation = True  ; False = no interpolation of the calculated
                                             ;         diurnal cycle
                                             ; True  = interpolation of the calculated
                                             ;         diurnal cycle (for input data with
                                             ;         a coarse temporal resolution such as
                                             ;         e.g. 3-hourly values)

diag_script_info@multi_model_mean = "y"
diag_script_info@supporting_gridlines = "y"
diag_script_info@x_gridlines = "n"
diag_script_info@y_gridlines = "n"

