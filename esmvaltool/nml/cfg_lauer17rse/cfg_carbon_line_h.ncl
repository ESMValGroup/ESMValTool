diag_script_info = True
diag_script_info@styleset = "CMIP5"

;info for carbon_tsline.ncl
;diag_script_info@ts_minlev = 
;diag_script_info@ts_maxlev = 
diag_script_info@ts_minlat = -90
diag_script_info@ts_maxlat = 90    
diag_script_info@ts_minlon = 0
diag_script_info@ts_maxlon = 360
diag_script_info@multi_model_mean = "True"
;diag_script_info@ts_maxyear = "2014"
;diag_script_info@ts_minyear = "2000"
;diag_script_info@ts_anomaly = "anom"
;diag_script_info@ts_detrend = "detr"
;diag_script_info@EMs_in_lg = False

diag_script_info@dont_add_per_year = True

if (any(variables(0).eq.(/"xco2","xch4","xco2Nobs"/))) then 
    diag_script_info@ts_diff_ref = "ESACCI-GHG"
    ;diag_script_info@time_avg   ="monthly"
    diag_script_info@ts_lsmask = "refmodel"
else if (any(variables(0).eq.(/"tas"/))) then 
    diag_script_info@ts_diff_ref = "ERA-Interim"
    diag_script_info@time_avg   = "yearly"
end if
end if




