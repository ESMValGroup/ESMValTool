diag_script_info = True
diag_script_info@styleset = "CMIP5"

;info for carbon_tsline.ncl
;diag_script_info@ts_minlev = 
;diag_script_info@ts_maxlev = 
if (variables.ne."tas") then
    diag_script_info@ts_minlat = -90
    diag_script_info@ts_maxlat = 90
else
    diag_script_info@ts_minlat = -30
    diag_script_info@ts_maxlat = 30    
end if
diag_script_info@ts_minlon = 0
diag_script_info@ts_maxlon = 360
diag_script_info@multi_model_mean = False
diag_script_info@ts_maxyear = "2005"
diag_script_info@ts_minyear = "1960"
diag_script_info@time_avg   = "yearly"
diag_script_info@ts_anomaly = "anom"
diag_script_info@ts_detrend = "detr"

if (any(variables.eq.(/"nbp","cumnbp","NBP"/))) then
    diag_script_info@area_opper = "sum"
end if
    
;diag_script_info@ts_anom = True
;diag_script_info@ts_int  = True
;diag_script_info@ts_detr = False




